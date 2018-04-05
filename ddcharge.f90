!
!  The deformation density partial charge program
!
!  Calculates charges from a GAMESS RDM output file using the Voronoi
!  deformation density or Hirshfeld algorithms, or their iterative
!  solutions. Radial atomic densities are from a sperically averaged
!  ab initio library or integrated Slater densities.
!
program vddi
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use atomdens
    !
    character(100)           :: rdm_file="rdm.dat", atyp="slater", alib="scf.cc-pvdz",  wtyp="hirshfeld", ityp="exp"
    character(2),allocatable :: atypes(:)
    integer(ik)              :: pfile, hfile, cfile, mfile, wfile
    integer(ik)              :: i, j, ib, ipt, iter, npts
    integer(ik)              :: rdm_count, natom, nbatch, nuniq, norb
    integer(ik)              :: nrad=70, nang=110, narad=70, maxiter=40, iordr=8
    integer(ik),allocatable  :: iwhr(:), qlist(:)
    real(rk)                 :: norm, normatm, tmp_rdm_sv(1000), thrsh=1e-3
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), ax(:,:), aw(:,:), awgt(:,:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:), acden(:,:,:)
    real(ark),allocatable    :: rhomol(:), rhopro(:), rhoatm(:,:), rdm_sv(:)
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
    open(mfile, file='moldens.dat', action='write', status='new')
    norm = 0
    nullify(xyzw)
    mol_grid_batches: do ib = 1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
        xyz => xyzw(1:3,:)
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
        call evaluate_density(rdm_count, npts, mol, rdm_sv, xyz, rhomol)
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
    !  Import spherically averaged densities
    !
    call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
    allocate (qlist(nuniq), ax(nuniq,narad), aw(nuniq,narad))
    allocate (acden(7,nuniq,narad), aden(natom,narad))
    qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
    select case (atyp)
        case ("abinitio")
            !open(lfile, "/home/rymac/Projects/PartialCharge/ddcharge/atomlib/" // alib)
            ! everything in this section is temporary, need better import
            open(hfile, file='/home/rymac/Projects/PartialCharge/ddcharge/atomlib/hydrog')
            open(cfile, file='/home/rymac/Projects/PartialCharge/ddcharge/atomlib/carbon')
            import_abinit: do i = 1,nuniq
                if (qlist(i) == 1) then
                    do j = 1,narad
                      read(hfile,*) ax(i,j), acden(1,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(hfile,*) ax(i,j), acden(2,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(hfile,*) ax(i,j), acden(3,i,j), aw(i,j)
                    end do
                else if (qlist(i) == 6) then
                    do j = 1,narad
                      read(cfile,*) ax(i,j), acden(1,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(cfile,*) ax(i,j), acden(2,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(cfile,*) ax(i,j), acden(3,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(cfile,*) ax(i,j), acden(4,i,j), aw(i,j)
                    end do
                    do j = 1,narad
                      read(cfile,*) ax(i,j), acden(5,i,j), aw(i,j)
                    end do
                end if
            end do import_abinit
            close(hfile)
            close(cfile)
        case ("slater")
            import_slater: do i = 1,nuniq
                ! these don't work at the moment
                call rlegendre(narad, qlist(i), ax(i,:), aw(i,:))
                norb = get_norb(qlist(i))
                call psisq(ax(i,:), aw(i,:), narad, qlist(i),  0, norb, 50, 1e-6_rk, acden(1,i,:))
                norb = get_norb(qlist(i)+1)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), -1, norb, 50, 1e-6_rk, acden(2,i,:))
                norb = get_norb(qlist(i)-1)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), +1, norb, 50, 1e-6_rk, acden(3,i,:))
                norb = get_norb(qlist(i)+2)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), -2, norb, 50, 1e-6_rk, acden(4,i,:))
                norb = get_norb(qlist(i)-2)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), +2, norb, 50, 1e-6_rk, acden(5,i,:))
                norb = get_norb(qlist(i)+3)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), -3, norb, 50, 1e-6_rk, acden(6,i,:))
                norb = get_norb(qlist(i)-3)
                call psisq(ax(i,:), aw(i,:), narad, qlist(i), +3, norb, 50, 1e-6_rk, acden(7,i,:))
            end do import_slater
        case default
            write(out,"('atomic_dens: Error unrecognized type',s)") atyp
            stop "atomic_dens - unrecognized type"
    end select

    !
    !  Atomic density numerical integration loop
    !
    open(mfile, file='moldens.dat', action='read')
    open(pfile, file='density.dat', action='write', status='new')
    open(wfile, file='weights.dat', action='write', status='new')
    charge(:) = 0
    iterate_chg: do iter = 1,maxiter
        rewind mfile
        !
        !  Set up the radial atomic densities
        !
        setup_atom: do i = 1,natom
            if (abs(charge(i)) > 2) then
                write(out,"('setup_atom: Error |q',d,'| > 2')") i
                stop "setup_atom - charge greater than 2"
            else if (charge(i) < -1) then
                aden(i,:) = (2 + charge(i))*acden(2,iwhr(i),:) - (1 + charge(i))*acden(4,iwhr(i),:)
            else if (charge(i) > +1) then
                aden(i,:) = (2 - charge(i))*acden(3,iwhr(i),:) + (charge(i) - 1)*acden(5,iwhr(i),:)
            else if (charge(i) < 0) then
                aden(i,:) = (1 + charge(i))*acden(1,iwhr(i),:) - charge(i)*acden(2,iwhr(i),:)
            else if (charge(i) > 0) then
                aden(i,:) = (1 - charge(i))*acden(1,iwhr(i),:) + charge(i)*acden(3,iwhr(i),:)
            else
                aden(i,:) = acden(1,iwhr(i),:)
            end if
        end do setup_atom
        call GridPointsBatch(den_grid, 'Reset')
        normatm    = 0
        dcharge(:) = 0
        nullify(xyzw)
        grid_batches: do ib = 1,nbatch
            !
            !  Get grid points
            !
            call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
            !
            !  If the batch size changed, reallocate density and weight arrays
            !
            npts = size(xyzw, dim=2)
            if (allocated(rhopro)) then
                if (size(rhomol) /= npts) deallocate (rhomol)
                if (size(rhopro) /= npts) deallocate (rhopro)
                if (size(rhoatm,2) /= npts) deallocate (rhoatm)
                if (size(awgt,2) /= npts) deallocate (awgt)
            end if
            if (.not. allocated(rhomol)) allocate (rhomol(npts))
            if (.not. allocated(rhopro)) allocate (rhopro(npts))
            if (.not. allocated(rhoatm)) allocate (rhoatm(natom,npts))
            if (.not. allocated(awgt)) allocate (awgt(natom,npts))
            !
            !  Read in molecular densities from file
            !
            read_rho: do ipt = 1,npts
                read(mfile,*) rhomol(ipt)
            end do read_rho
            !
            !  Evaluate atomic densities at grid points
            !
            call evaluate_atomic(aden, ax, narad, xyzq(1:3,:), iwhr, nuniq, natom, ityp, iordr, xyzw(1:3,:), npts, rhoatm)
            rhopro = sum(rhoatm, dim=1)
            !
            !  Determine atomic assignments
            !
            awgt = assign_atom(xyzq(1:3,:), natom, xyzw(1:3,:), npts, rhoatm, wtyp)
            !
            !  Integrate
            !
            integrate: do ipt = 1,npts
                normatm = normatm + xyzw(4,ipt) * rhopro(ipt)
                !if (iter == maxiter) then ! Only want this saved for the last iteration?
                write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhopro(ipt)-rhomol(ipt)
                write(wfile,'(*(e16.8))') awgt(:,ipt)
                !end if
            end do integrate
            !
            !  Find the atomic contributions
            !
            atom_contrib: do i = 1,natom
                dcharge(i) = dcharge(i) + sum(awgt(i,:) * xyzw(4,:) * (rhopro(:) - rhomol(:)))
            end do atom_contrib
        end do grid_batches
        charge(:) = charge(:) + dcharge(:)

        print *,'ITER: ', iter
        print *,'TOTAL ATOMIC DENSITY: ', normatm
        print *,'TOTAL CHANGE IN CHARGE: ', sum(dcharge)
        print *,'CONTRIBUTIONS: ', dcharge
        print *,''

        if (maxval(abs(dcharge)) < thrsh) then
            print *,"CHARGES CONVERGED ON ITER ", iter
            print *,"FINAL CHARGES: ", charge
            exit
        else if (iter == maxiter) then
            print *,"CHARGES NOT CONVERGED FOR MAXITER"
            print *,"UNCONVERGED CHARGES: ", charge
        end if
    end do iterate_chg
    !
    !  Clean up
    !
    close(pfile)
    close(mfile)
    close(wfile)
    call GridDestroy(den_grid)

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
!  Determine the density at XYZ coordinates
!
subroutine evaluate_density(rdm_count, npt, mol, rdm_sv, xyz, rho)
    integer(ik),intent(in) :: rdm_count, npt
    type(gam_structure), intent(inout)   :: mol
    real(rk), intent(in)   :: rdm_sv(rdm_count)
    real(rk), intent(in)   :: xyz(3,npt)
    real(rk), intent(out)  :: rho(npt)
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
!  Determine the atomic density contributions at grid XYZ coordinates
!
subroutine evaluate_atomic(rden, r, nr, xyzc, uind, nu, natm, ityp, iord, xyz, npt, rho)
    character(3) :: ityp
    integer(ik)  :: nr, natm, npt, nu, uind(natm), iord
    real(rk)     :: r(nu,nr), xyzc(3,natm), xyz(3,npt)
    real(ark)    :: rden(natm,nr), rho(natm,npt)
    !
    integer(ik)  :: iatm, ipt, ui
    real(rk)     :: ratm

    rho(:,:) = 0
    evaluate_atom: do iatm = 1,natm
        ui = uind(iatm)
        interp_point: do ipt = 1,npt
            ratm = sqrt((xyz(1,ipt)-xyzc(1,iatm))**2 + &
                        (xyz(2,ipt)-xyzc(2,iatm))**2 + &
                        (xyz(3,ipt)-xyzc(3,iatm))**2)
            ! atomic density = 0 outside the grid
            if (ratm < r(ui,nr)) then
                rho(iatm,ipt) = rho(iatm,ipt) + interp(ratm, r(ui,:), rden(iatm,:), nr, ityp, iord)
            end if
        end do interp_point
    end do evaluate_atom

end subroutine evaluate_atomic

!
!  Determine the atomic weight factors on a grid
!
function assign_atom(xyzatm, natm, xyz, npt, rho, wtyp)
    character(*) :: wtyp
    integer(ik)  :: natm, npt
    real(rk)     :: xyzatm(3,natm), xyz(3,npt), rho(natm,npt)
    real(rk)     :: assign_atom(natm,npt)
    !
    integer(ik)  :: ipt, imin, iatm
    real(rk)     :: dist(natm), rhotot(npt)

    select case (wtyp)
        case ("voronoi")
            assign_atom(:,:) = 0
            do ipt = 1,npt
                dist = sqrt((xyz(1,ipt)-xyzatm(1,:))**2 + &
                            (xyz(2,ipt)-xyzatm(2,:))**2 + &
                            (xyz(3,ipt)-xyzatm(3,:))**2)
                imin = 1
                find_min: do iatm = 2,natm
                    if (dist(iatm) < dist(imin)) then
                        imin = iatm
                    end if
                end do find_min
                assign_atom(imin,ipt) = 1.
            end do
        case ("hirshfeld")
            rhotot = sum(rho, dim=1)
            assign_atom(:,:) = 0
            do ipt = 1,npt
                if (rhotot(ipt) > 0) then
                    do iatm = 1,natm
                        assign_atom(iatm,ipt) = rho(iatm,ipt) / rhotot(ipt)
                    end do
                end if
            end do
        case default
            write(out,"('assign_atom: Error unrecognized type',s)") wtyp
            stop "assign_atom - unrecognized type"
    end select
end function assign_atom

end program vddi
