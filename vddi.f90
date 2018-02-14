!
!  The iterative Voronoi deformation density program
!
program vddi
  use accuracy
  use import_gamess
  use gamess_internal
  use molecular_grid
  use atomdens
  !
  character(100)           :: rdm_file='rdm.dat'  ! Name of the file containing rdm orbital coeffs
  character(2),allocatable :: atypes(:)
  integer(ik)              :: i, ib, ipt, iter, ios, npts, pfile
  integer(ik)              :: rdm_count, natom, nbatch, nor, nuniq
  integer(ik)              :: nrad=70, nang=110, nx=100000, niter=50
  integer(ik),allocatable  :: indi(:), iwhr(:), qlist(:)
  real(rk)                 :: tmp_rdm_sv(1000)
  real(rk)                 :: xmin=0., xmax=3000., thrsh=1e-6
  real(rk),pointer         :: xyzw(:,:)           ! Coordinates and weights of the grid points
  real(rk),allocatable     :: xyzq(:,:), x(:), norm(:), normatm(:), charge(:)
  real(ark),allocatable    :: psi2m(:,:), psi2n(:,:), psi2p(:,:), psi2(:,:)
  real(ark),allocatable    :: rho(:), rhoatm(:), rdm_sv(:)
  type(gam_structure)      :: mol                 ! Gamess basis set and orbital data
  type(mol_grid)           :: den_grid

  call accuracyInitialize

  pfile=11
  open(pfile, file='density.dat', action='write', status='new', iostat=ios)

  call gamess_load_rdmsv(trim(rdm_file), tmp_rdm_sv, rdm_count)
  write (out,"( 'Found ',i4,' singular values')") rdm_count
  write (out,"( 'Values are: '/)")
  write (out,"(10(1x,f12.8))") tmp_rdm_sv(:rdm_count)
  !
  allocate (rdm_sv(rdm_count))
  rdm_sv = tmp_rdm_sv(:rdm_count)
  !
  call gamess_load_orbitals(file=trim(rdm_file), structure=mol)

  allocate (atypes(mol%natoms), xyzq(4,mol%natoms), iwhr(mol%natoms))
  allocate (norm(mol%natoms), normatm(mol%natoms), charge(mol%natoms))
  allocate (psi2(mol%natoms,nx))

  call gamess_report_nuclei(natom, xyzq, mol)
  do i = 1,mol%natoms
    atypes(i) = trim(mol%atoms(i)%name)
  enddo
  !
  !  Evaluate the atomic densities for interpolation
  !
  call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
  allocate (qlist(nuniq), x(nx))
  allocate (psi2m(nuniq,nx), psi2n(nuniq,nx), psi2p(nuniq,nx))
  qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
  x = linspace(xmin, xmax, nx)
  !
  setup_atom: do i = 1,nuniq
    nor = get_norb(qlist(i))
    call psisq(x, nx, qlist(i), -1, nor, 100, thrsh, psi2m(i,:))
    call psisq(x, nx, qlist(i), 0, nor, 100, thrsh, psi2n(i,:))
    call psisq(x, nx, qlist(i), 1, nor, 100, thrsh, psi2p(i,:))
  end do setup_atom
  !
  charge(:) = 0
  iterate_chg: do iter = 1,niter
    call GridInitialize(den_grid, nrad, nang, xyzq(1:3,:), atypes)
    !
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    !
    setup_charge: do i = 1,natom
      if (charge(i) < 0) then
        psi2(i,:) = (1+charge(i))*psi2n(iwhr(i),:) - charge(i)*psi2m(iwhr(i),:)
      else if (charge(i) > 0) then
        psi2(i,:) = (1-charge(i))*psi2n(iwhr(i),:) + charge(i)*psi2p(iwhr(i),:)
      else
        psi2(i,:) = psi2n(iwhr(i),:)
      end if
    end do setup_charge
    !
    charge(:)  = 0
    norm(:)    = 0
    normatm(:) = 0
    !
    !  Numerical integration loop
    !
    nullify(xyzw)
    grid_batches: do ib = 1,nbatch
       !
       !  Get grid points
       !
       call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
       !
       !  If the batch size changed, reallocate rho()
       !
       npts = size(xyzw, dim=2)
       if (allocated(rho)) then
         if (size(rho) /= npts) deallocate (rho)
         if (size(rhoatm) /= npts) deallocate (rhoatm)
         if (size(indi) /= npts) deallocate (indi)
       end if
       if (.not. allocated(rho)) allocate (rho(npts))
       if (.not. allocated(rhoatm)) allocate (rhoatm(npts))
       if (.not. allocated(indi)) allocate (indi(npts))
       !
       !  Evaluate density at grid points
       !
       call evaluate_density(rdm_count, npts, mol, rdm_sv, xyzw(1:3,:), rho)
       !
       !  Evaluate total atomic density at grid points
       !
       call evaluate_atomic(psi2, x, nx, xyzq, natom, xyzw(1:3,:), npts, rhoatm)
       !
       !  Determine Voronoi assignments
       !
       indi = assign_voronoi(xyzq(1:3,:), natom, xyzw(1:3,:), npts)
       !
       !  Ready to integrate!
       !
       integrate: do ipt = 1,npts
         norm(indi(ipt))    = norm(indi(ipt))    + xyzw(4,ipt) * rho(ipt)
         normatm(indi(ipt)) = normatm(indi(ipt)) + xyzw(4,ipt) * rhoatm(ipt)
         charge(indi(ipt))  = charge(indi(ipt))  + xyzw(4,ipt) * (rhoatm(ipt)-rho(ipt))
         if (iter == niter) then
           write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhoatm(ipt)-rho(ipt), indi(ipt)
         end if
       end do integrate
    end do grid_batches

    close(pfile)

    call GridDestroy(den_grid)

    print *,'ITER: ', iter
    print *,'TOTAL MOLECULAR DENSITY: ', sum(norm)
    print *,'CONTRIBUTIONS: ', norm
    print *,'TOTAL ATOMIC DENSITY: ', sum(normatm)
    print *,'CONTRIBUTIONS: ', normatm
    print *,''
    print *,'ATOMIC CHARGES: ', charge
    print *,''
  end do iterate_chg

1000 format(e24.8,f16.8,f16.8,f16.8,es18.10,i10)

contains

!
!  Find unique elements of a list and return a mask
!
subroutine unique(arr, narr, mask, nuniq)
    integer(ik),intent(in)  :: narr, arr(narr)
    integer(ik),intent(out) :: nuniq, mask(narr)
    !
    integer(ik) :: i, j, ident(100)

    nuniq = 0
    mask(:) = 0
    parse_list: do i = 1,narr
        check_list: do j = 1,nuniq
            if (arr(i) == ident(j)) then
                mask(i) = j
                exit
            end if
        end do check_list
        if (mask(i) == 0) then
            nuniq = nuniq + 1
            mask(i) = nuniq
            ident(nuniq) = arr(i)
        end if
    end do parse_list
end subroutine unique

!
!  Return a list of unique elements
!
function unique_list(arr, narr, mask, nuniq)
    integer(ik) :: narr, nuniq, arr(narr), mask(narr), unique_list(nuniq)
    !
    integer(ik) :: i, j

    do i = 1,nuniq
        do j = 1,narr
            if (mask(j) == i) then
                unique_list(i) = arr(j)
                exit
            end if
        end do
    end do
end function unique_list

!
!  Generate an equally spaced list of n elements between specified numbers
!
function linspace(first, last, n)
  integer(ik) :: n
  real(rk)    :: first, last, linspace(n)
  !
  integer(ik) :: i
  real(rk)    :: step, list(n)

  list = (/(i, i=0,n-1)/)
  step = (last - first) / (n - 1)
  linspace = list*step + first
end function linspace

!
!  Determine the density at XYZ coordinates
!
subroutine evaluate_density(rdm_count, npt, mol, rdm_sv, xyz, rho)
    integer(ik),intent(in) :: rdm_count, npt
    type(gam_structure), intent(inout)   :: mol
    real(rk), intent(in)   :: rdm_sv(rdm_count)
    real(rk), intent(in)   :: xyz(3,npt) ! XYZ coordinates of the grid points
    real(rk), intent(out)  :: rho(npt)   ! Transition density at these grid points
    !
    integer(ik)            :: ipt, ird, imo, nmo, nbas
    real(ark), allocatable :: basval(:,:,:)
    real(rk), allocatable  :: moval (:,:)

    nmo  = mol%nvectors
    nbas = mol%nbasis
    !
    allocate (basval(1,nbas,npt), moval(nmo,npt))
    !
    !  First, evaluate basis functions
    !
    evaluate_basis_functions: do ipt = 1,npt
      call gamess_evaluate_functions(xyz(:,ipt), basval(:,:,ipt), mol)
    end do evaluate_basis_functions
    !
    !  Transform AOs to the MOs, for all grid points simultaneously
    !
    moval = matmul(transpose(mol%vectors(:,:nmo)), basval(1,:,:))
    !
    !  Finally, evaluate the transition density at grid points
    !
    rho(:) = 0
    evaluate_rdm: do ird = 1,rdm_count
      imo = 2*ird - 1
      rho = rho + rdm_sv(ird) * moval(imo,:) * moval(imo+1,:)
    end do evaluate_rdm
    !
    deallocate (basval,moval)

end subroutine evaluate_density

!
!  Determine the total atomic density at XYZ coordinates
!
subroutine evaluate_atomic(den, r, nr, xyzq, natm, xyz, npt, rhoatm)
    integer(ik) :: nr, natm, npt
    real(rk)    :: r(nr), xyzq(4,natm), xyz(3,npt)
    real(ark)   :: den(natm, nr), rhoatm(npt)
    !
    integer(ik) :: iatm
    real(rk)    :: ratm(npt)
    real(rk),parameter :: fourpi=12.566370614359

    rhoatm(:) = 0
    evaluate_atom: do iatm = 1,natm
        ratm = sqrt((xyz(1,:)-xyzq(1,iatm))**2 + (xyz(2,:)-xyzq(2,iatm))**2 + &
                    (xyz(3,:)-xyzq(3,iatm))**2)
        rhoatm = rhoatm + interp(ratm, npt, r, den(iatm,:), nr)/fourpi
    end do evaluate_atom

end subroutine evaluate_atomic

!
!  Interpolate a set of points on a function f(x)
!
function interp(xpt, npt, x, fx, nx)
    integer(ik) :: npt, nx
    real(rk)    :: xpt(npt), x(nx)
    real(ark)   :: fx(nx), interp(npt)
    !
    integer(ik) :: ipt, ix, i
    real(rk)    :: d

    interp_point: do ipt = 1,npt
        ! find the nearest index
        find_ind: do ix = 2,nx
            if (xpt(ipt) < x(ix) .or. ix == nx) then
                i = ix - 1
                exit
            end if
        end do find_ind
        ! find the distance
        d = (xpt(ipt) - x(i)) / (x(i+1) - x(i))
        ! interpolate to find f(xpt)
        interp(ipt) = fx(i)*(1-d) + fx(i+1)*d
    end do interp_point
end function interp

!
!  Determine the Voronoi cell assignment for a grid
!
function assign_voronoi(xyzatm, natm, xyz, npt)
    integer(ik) :: natm, npt, assign_voronoi(npt)
    real(rk)    :: xyzatm(3,natm), xyz(3,npt)
    !
    integer(ik) :: ipt, iatm
    real(rk)    :: dist(natm)

    assign_point: do ipt = 1,npt
        dist = sqrt((xyz(1,ipt)-xyzatm(1,:))**2 + &
                    (xyz(2,ipt)-xyzatm(2,:))**2 + (xyz(3,ipt)-xyzatm(3,:))**2)
        assign_voronoi(ipt) = 1
        find_min: do iatm = 2,natm
            if (dist(iatm) < dist(assign_voronoi(ipt))) then
                assign_voronoi(ipt) = iatm
            end if
        end do find_min
    end do assign_point
end function assign_voronoi

end program vddi
