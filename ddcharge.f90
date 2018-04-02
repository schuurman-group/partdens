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
    character(100)           :: rdm_file='rdm.dat'
    character(2),allocatable :: atypes(:)
    integer(ik)              :: i, j, ib, ipt, iter, ios, npts
    integer(ik)              :: pfile, hfile, cfile, mfile, wfile
    integer(ik)              :: rdm_count, natom, nbatch, nuniq
    integer(ik)              :: nrad=70, nang=110, narad=70, maxiter=1
    integer(ik),allocatable  :: aind(:), iwhr(:), qlist(:)
    real(rk)                 :: tmp_rdm_sv(1000), thrsh=1e-2
    real(rk),pointer         :: xyzw(:,:)
    real(rk),allocatable     :: xyzq(:,:), ax(:,:), awgt(:,:)
    real(rk),allocatable     :: norm(:), normatm(:), charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:)
    real(ark),allocatable    :: rhomol(:), rhoatm(:), rdm_sv(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid

    call accuracyInitialize

    mfile=11
    pfile=12
    wfile=13
    hfile=111
    cfile=116

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
    allocate (norm(mol%natoms), normatm(mol%natoms))
    allocate (charge(mol%natoms), dcharge(mol%natoms))

    call gamess_report_nuclei(natom, xyzq, mol)
    do i = 1,mol%natoms
        atypes(i) = trim(mol%atoms(i)%name)
    enddo
    !
    !  Set up the grid
    !
    call GridInitialize(den_grid, nrad, nang, xyzq(1:3,:), atypes)
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    !
    !  Molecular density numerical integration loop
    !
    open(mfile, file='moldens.dat', action='write', status='new', iostat=ios)
    norm(:) = 0
    nullify(xyzw)
    mol_grid_batches: do ib = 1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
        !
        !  If the batch size changed, reallocate rho()
        !
        npts = size(xyzw, dim=2)
        if (allocated(rhomol)) then
            if (size(rhomol) /= npts) deallocate (rhomol)
        end if
        if (.not. allocated(rhomol)) allocate (rhomol(npts))
        !
        !  Evaluate density at grid points
        !
        call evaluate_density(rdm_count, npts, mol, rdm_sv, xyzw(1:3,:), rhomol)
        !
        !  Integrate and save
        !
        mol_integrate: do ipt = 1,npts
            norm = norm + xyzw(4,ipt) * rhomol(ipt)
            write(mfile,"(e18.10)") rhomol(ipt)
        end do mol_integrate
    end do mol_grid_batches
    close(mfile)
    print *,'TOTAL MOLECULAR DENSITY: ', norm

    !
    !  Atomic density numerical integration loop
    !
    open(mfile, file='moldens.dat', action='read')
    open(pfile, file='density.dat', action='write', status='new', iostat=ios)
    rhoatm(:) = 0
    charge(:) = 0
    iterate_chg: do iter = 1,maxiter
        !
        !  Evaluate the atomic densities for interpolation
        !
        rewind mfile
        ! the two files below are temporary, for tests.
        open(hfile, file='/home/rymac/Projects/PartialCharge/ddcharge/atomlib/hydrog')
        open(cfile, file='/home/rymac/Projects/PartialCharge/ddcharge/atomlib/carbon')
        call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
        allocate (qlist(nuniq))
        qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
        allocate (ax(nuniq,narad), aden(nuniq,narad), awgt(nuniq,narad))
        setup_atom: do i = 1,nuniq
        if (qlist(i) == 1) then
            do j=1,narad
              read(hfile,*) ax(i,j), aden(i,j), awgt(i,j)
            end do
        else if (qlist(i) == 6) then
            do j=1,narad
              read(cfile,*) ax(i,j), aden(i,j), awgt(i,j)
            end do
        end if
        end do setup_atom
        !setup_charge: do i = 1,natom
        !  if (charge(i) < 0) then
        !    psi2(i,:) = (1+charge(i))*psi2n(iwhr(i),:) - charge(i)*psi2m(iwhr(i),:)
        !  else if (charge(i) > 0) then
        !    psi2(i,:) = (1-charge(i))*psi2n(iwhr(i),:) + charge(i)*psi2p(iwhr(i),:)
        !  else
        !    psi2(i,:) = psi2n(iwhr(i),:)
        !  end if
        !end do setup_charge
        !
        call GridPointsBatch(den_grid, 'Reset')
        dcharge(:) = 0
        normatm(:) = 0
        nullify(xyzw)
        grid_batches: do ib = 1,nbatch
            !
            !  Get grid points
            !
            call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
            !
            !  If the batch size changed, reallocate rhoatm() and aind()
            !
            npts = size(xyzw, dim=2)
            if (allocated(rhoatm)) then
                if (size(rhomol) /= npts) deallocate (rhomol)
                if (size(rhoatm) /= npts) deallocate (rhoatm)
                if (size(aind) /= npts) deallocate (aind)
            end if
            if (.not. allocated(rhomol)) allocate (rhomol(npts))
            if (.not. allocated(rhoatm)) allocate (rhoatm(npts))
            if (.not. allocated(aind)) allocate (aind(npts))
            !
            !  Read in molecular densities from file
            !
            do i=1,npts
                read(mfile,*) rhomol(i)
            end do
            !
            !  Evaluate atomic densities at grid points
            !
            call evaluate_atomic(aden, ax, narad, xyzq, iwhr, nuniq, natom, xyzw(1:3,:), npts, rhoatm)
            !
            !  Determine atomic assignments
            !
            aind = assign_voronoi(xyzq(1:3,:), natom, xyzw(1:3,:), npts)
            !
            !  Integrate
            !
            integrate: do ipt = 1,npts
                normatm(aind(ipt)) = normatm(aind(ipt)) + xyzw(4,ipt) * rhoatm(ipt)
                dcharge(aind(ipt)) = dcharge(aind(ipt)) + xyzw(4,ipt) * (rhoatm(ipt)-rhomol(ipt))
                !if (iter == maxiter) then
                write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhoatm(ipt)-rhomol(ipt)
                !end if
            end do integrate
        end do grid_batches
        charge(:) = charge(:) + dcharge(:)

        print *,'ITER: ', iter
        print *,'TOTAL ATOMIC DENSITY: ', sum(normatm)
        print *,'TOTAL CHANGE IN CHARGE: ', sum(dcharge)
        print *,'CONTRIBUTIONS: ', dcharge
        print *,''
    end do iterate_chg
    close(pfile)
    close(mfile)
    call GridDestroy(den_grid)

    print *,'FINAL CHARGES: ', charge

    1000 format(e24.8,f16.8,f16.8,f16.8,es18.10)

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
subroutine evaluate_atomic(den, r, nr, xyzq, uind, nu, natm, xyz, npt, rho)
    integer(ik) :: nr, natm, npt, nu, uind(natm)
    real(rk)    :: r(nu,nr), xyzq(4,natm), xyz(3,npt)
    real(ark)   :: den(nu,nr), rho(npt)
    !
    integer(ik) :: iatm, ipt, ui
    real(rk)    :: ratm

    rho(:) = 0
    evaluate_atom: do iatm = 1,natm
        ui = uind(iatm)
        interp_point: do ipt = 1,npt
            ratm = sqrt((xyz(1,ipt)-xyzq(1,iatm))**2 + &
                        (xyz(2,ipt)-xyzq(2,iatm))**2 + &
                        (xyz(3,ipt)-xyzq(3,iatm))**2)
            ! atomic density = 0 outside the grid
            if (ratm < r(ui,nr)) then
                rho(ipt) = rho(ipt) + interp(ratm, r(ui,:), den(ui,:), nr, "exp", 8)
            end if
        end do interp_point
    end do evaluate_atom

end subroutine evaluate_atomic

!
!  Interpolate a point on a function f(x)
!
function interp(xpt, x, fx, nx, ityp, ordr)
    character(3) :: ityp
    integer(ik)  :: nx, ordr
    real(rk)     :: xpt, x(nx)
    real(ark)    :: fx(nx), interp
    !
    integer(ik)  :: i, j, ix, ind, n(ordr+1)

    ! find the nearest index
    find_ind: do ix = 2,nx
        if (xpt < x(ix) .or. ix == nx) then
            ind = ix - 1
            exit
        end if
    end do find_ind
    ! get the range of indices for interpolation
    if (ind < ordr/2 + 1) then
        n = (/(i, i = 1,ordr+1)/)
    else if (ind+ordr/2+mod(ordr,2) > nx) then
        n = (/(i, i = nx-ordr,nx)/)
    else
        n = (/(i, i = ind-ordr/2,ind+ordr/2+mod(ordr,2))/)
    end if
    select case (ityp)
        case ("pol")
            ! polynomial interpolation
            interp = 0.
            do i = 1,ordr+1
                ! don't bother adding zeros
                if (fx(n(i)) /= 0) then
                    interp = interp + fx(n(i))*x_li(xpt, x, nx, i, n, ordr)
                end if
            end do
        case ("exp")
            ! exponential interpolation
            interp = 1.
            do i = 1,ordr+1
                ! a f(x) value of 0 kills the interpolation
                if (fx(n(i)) == 0) then
                    interp = 0
                    return
                else
                    interp = interp * fx(n(i))**x_li(xpt, x, nx, i, n, ordr)
                end if
            end do
        case default
            write (out,"('interp: Error unrecognized type ',s)") ityp
            stop 'interp - unrecognized type'
    end select
end function interp

!
!  Find the x weighting factor for interpolation
!
function x_li(xpt, x, nx, i, n, ordr)
    integer(ik) :: i, ordr, nx, n(ordr+1)
    real(rk)    :: xpt, x(nx), x_li
    !
    integer(ik) :: j
    x_li = 1.
    do j = 1,ordr+1
        if (j /= i) then
            x_li = x_li * (xpt - x(n(j))) / (x(n(i)) - x(n(j)))
        end if
    end do
end function x_li

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
                    (xyz(2,ipt)-xyzatm(2,:))**2 + &
                    (xyz(3,ipt)-xyzatm(3,:))**2)
        assign_voronoi(ipt) = 1
        find_min: do iatm = 2,natm
            if (dist(iatm) < dist(assign_voronoi(ipt))) then
                assign_voronoi(ipt) = iatm
            end if
        end do find_min
    end do assign_point
end function assign_voronoi

end program vddi
