program atom_pop
  use accuracy
  use import_gamess
  use gamess_internal
  use molecular_grid
  use atomdens
  !
  integer(ik)            :: i, ib, ipt, ios, npts, pfile, rdm_count, natom, nbatch
  real(ark), allocatable :: rho(:), rhoatm(:), rdm_sv(:)
  type(gam_structure)    :: mol                 ! Gamess basis set and orbital data
  type(mol_grid)         :: den_grid
  character(len=100)     :: rdm_file='rdm.dat'  ! Name of the file containing rdm orbital coeffs
  real(rk)               :: tmp_rdm_sv(1000)
  real(rk),allocatable   :: xyzq(:,:), norm(:), normatm(:), charge(:)
  integer(ik),allocatable :: indi(:)
  integer(ik)            :: nrad=70, nang=110
  character(len=2),allocatable :: atypes(:)
  real(rk), pointer      :: xyzw(:,:)           ! Coordinates and weights of the grid points

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

  allocate (atypes(mol%natoms), xyzq(4,mol%natoms))
  allocate (norm(mol%natoms), normatm(mol%natoms), charge(mol%natoms))

  call gamess_report_nuclei(natom, xyzq, mol)
  do i = 1,mol%natoms
    atypes(i) = trim(mol%atoms(i)%name)
  enddo
  !
  call GridInitialize(den_grid, nrad, nang, xyzq(1:3,:), atypes)

  call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
  !
  !  Numerical integration loop
  !
  charge(:) = 0
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
     call evaluate_atomic(100000, xyzq, natom, xyzw(1:3,:), npts, rhoatm)
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
       write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhoatm(ipt)-rho(ipt), indi(ipt)
     end do integrate
  end do grid_batches

  close(pfile)

  call GridDestroy(den_grid)

  print *,'TOTAL MOLECULAR DENSITY: ', sum(norm)
  print *,'CONTRIBUTIONS: ', norm
  print *,'TOTAL ATOMIC DENSITY: ', sum(normatm)
  print *,'CONTRIBUTIONS: ', normatm
  print *,''
  print *,'ATOMIC CHARGES: ', charge

1000 format(e24.8,f16.8,f16.8,f16.8,es18.10,i10)

contains

!
! Determine the density at XYZ coordinates
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
subroutine evaluate_atomic(nr, xyzq, natm, xyz, npt, rhoatm)
    integer(ik) :: natm, npt
    real(rk)    :: xyzq(4,natm), xyz(3,npt)
    real(ark)   :: rhoatm(npt)
    !
    integer(ik) :: iatm, nr, nor, z
    real(rk)    :: ratm(npt), r(nr)
    real(ark)   :: psi2(nr)
    real(rk),parameter :: rmin=0., rmax=3000., thrsh=1e-6
    real(rk),parameter :: fourpi=12.566370614359

    r = linspace(rmin, rmax, nr)
    rhoatm(:) = 0
    evaluate_atom: do iatm = 1,natm
        ! ideally shouldn't have to evaluate Psi^2 for every atom, just every unique
        z = int(xyzq(4,iatm))
        nor = get_norb(z)
        call psisq(r, nr, z, 0, nor, 100, thrsh, psi2)
        ratm = sqrt((xyz(1,:)-xyzq(1,iatm))**2 + (xyz(2,:)-xyzq(2,iatm))**2 + &
                    (xyz(3,:)-xyzq(3,iatm))**2)
        rhoatm = rhoatm + interp(ratm, npt, r, psi2, nr)/fourpi
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

end program atom_pop
