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
    use atoms
    use atomdens
    !
    character(100)           :: rdm_file="rdm.dat", atyp="abinitio", alib="scf.cc-pvdz",  wtyp="hirshfeld", ityp="exp"
    character(20)            :: pname, wname
    character(2),allocatable :: atypes(:)
    logical,allocatable      :: aload(:,:)
    integer(ik)              :: pfile, mfile, wfile
    integer(ik)              :: i, ib, ipt, iter, npts
    integer(ik)              :: rdm_count, natom, nbatch, nuniq, norb
    integer(ik)              :: nrad=70, nang=110, narad=70, maxiter=12, iordr=8, maxchg=2
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
    open(mfile, file='moldens', form='unformatted', action='write', status='new')
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
            write(mfile) rhomol(ipt)
        end do mol_integrate
    end do mol_grid_batches
    close(mfile)
    print '("Total molecular density: ",f14.8)', norm

    !
    !  Import atomic properties
    !
    call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
    allocate (qlist(nuniq), ax(nuniq,narad), aw(nuniq,narad))
    allocate (aload(2*maxchg+1,nuniq), acden(2*maxchg+1,nuniq,narad), aden(natom,narad))
    qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
    setup_rad: do i = 1,nuniq
        call rlegendre(narad, qlist(i), ax(i,:), aw(i,:))
    end do setup_rad

    !
    !  Atomic density numerical integration loop
    !
    open(mfile, file='moldens', form='unformatted', action='read')
    charge(:) = 0
    aload(:,:) = .false.
    iterate_chg: do iter = 1,maxiter
        print '("ITER = "i4)', iter
        if (iter < 10) then
            write(pname,'("density.",i1)') iter
            write(wname,'("weights.",i1)') iter
        else if (iter < 100) then
            write(pname,'("density.",i2)') iter
            write(wname,'("weights.",i2)') iter
        else
            write(pname,'("density.",i3)') iter
            write(wname,'("weights.",i3)') iter
        end if
        open(pfile, file=trim(pname), action='write', status='new')
        open(wfile, file=trim(wname), action='write', status='new')
        rewind mfile
        !
        !  Set up the radial atomic densities
        !
        call update_atoms(charge, maxchg, qlist, natom, iwhr, nuniq, narad, alib, atyp, aload, ax, aw, acden, aden)
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
                read(mfile) rhomol(ipt)
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
            !  Find the atomic contributions
            !
            normatm = normatm + sum(xyzw(4,:) * rhopro)
            atom_contrib: do i = 1,natom
                dcharge(i) = dcharge(i) + sum(awgt(i,:) * xyzw(4,:) * (rhopro - rhomol))
            end do atom_contrib
            !
            !  Output the densities for each point
            !
            integrate: do ipt = 1,npts
                write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhopro(ipt)-rhomol(ipt)
                write(wfile,'(*(e16.8))') awgt(:,ipt)
            end do integrate
        end do grid_batches
        close(pfile)
        close(wfile)
        charge = charge + dcharge

        do i = 1,natom
            if (charge(i) > qlist(iwhr(i))) then
                write(out,'("iterate_chg: Charge of ",f8.3," for atom ",i3," exceeds nuclear charge")') charge(i), i
                stop "iterate_chg - atomic charge exceeds nuclear charge"
            end if
        end do

        print '("Total atomic density: ",f14.8)', normatm
        print '("Total change in charge: ",f14.8)', sum(dcharge)
        print '("Contributions:")'
        print '("    ",5f14.8)', dcharge
        print '("")'

        if (maxval(abs(dcharge)) < thrsh) then
            print '("Charges converged on ITER = ",i4)', iter
            print '("Final charges:")'
            print '("    ",5f14.8)', charge
            exit
        else if (iter == maxiter) then
            print '("Charges not converged for MAXITER = ",i4)', maxiter
            print '("Unconverged charges:")'
            print '("    ",5f14.8)', charge
        end if
    end do iterate_chg
    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    close(mfile)
    call system("rm moldens")
    ! only keep data from the last iteration
    call system("mv "//pname//" density.dat")
    call system("mv "//wname//" weights.dat")
    call system("rm density.[0-9]* weights.[0-9]*")

    1000 format(e24.8,f16.8,f16.8,f16.8,e20.10)

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
!  Determine the molecular density at XYZ coordinates
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
!  Determine the atomic densities at grid XYZ coordinates
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
!  Update the atomic density functions as required
!
subroutine update_atoms(chg, maxc, ql, natm, uind, nu, narad, alib, atyp, aload, ax, aw, acden, aden)
    character(*) :: alib, atyp
    integer(ik)  :: natm, nu, narad, maxc, ql(nu), uind(natm)
    logical      :: aload(2*maxc+1,nu)
    real(rk)     :: chg(natm), ax(nu,narad), aw(nu,narad)
    real(ark)    :: acden(2*maxc+1,nu,narad), aden(natom,narad)
    !
    integer(ik)  :: ia, il, iu, ui
    real(rk)     :: modc, cthrsh=1e-3

    update_aden: do ia = 1,natm
        if (abs(chg(ia)) > maxc) then
            write(out,'("update_atoms: abs(charge) greater than MAXCHG = ",i3)') maxc
            stop "update_atoms - charge out of bounds"
        end if
        ui = uind(ia)
        ! find the modulus (as it should be defined)
        modc = chg(ia) - floor(chg(ia))
        if (abs(modc) < cthrsh) then
            ! integer charge, only get one contribution
            il = maxc + 1 + nint(chg(ia))
            if (.not. aload(il,ui)) then
                ! import the integer density
                call load_atom(ql(ui), nint(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il,ui,:))
                aload(il,ui) = .true.
            end if
            aden(ia,:) = acden(il,ui,:)
        else
            ! real charge, get ceil(chg) and floor(chg) contributions
            il = maxc + 1 + floor(chg(ia))
            iu = maxc + 1 + ceiling(chg(ia))
            if (.not. aload(il,ui)) then
                ! import the lower integer density
                call load_atom(ql(ui), floor(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il,ui,:))
                aload(il,ui) = .true.
            end if
            if (.not. aload(iu,ui)) then
                ! import the upper integer density
                call load_atom(ql(ui), ceiling(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(iu,ui,:))
                aload(iu,ui) = .true.
            end if
            aden(ia,:) = (1-modc)*acden(il,ui,:) + modc*acden(iu,ui,:)
        end if
    end do update_aden

end subroutine update_atoms

!
!  Load an atomic density
!
subroutine load_atom(q, chg, narad, alib, atyp, ax, aw, aden)
    character(*) :: alib, atyp
    integer(ik)  :: narad, q, chg
    real(rk)     :: ax(narad), aw(narad)
    real(ark)    :: aden(narad)
    !
    character(2) :: elem, fel
    integer(ik)  :: ipt, fch, fnr, ios, afile=111, niter=50
    real(rk)     :: thrsh=1e-6, junk

    select case (atyp)
        case ("abinitio")
            print '("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)', q, chg
            open(afile, file="/home/rymac/Projects/PartialCharge/ddcharge/atomlib/"//alib)
            elem = AtomElementSymbol(1._rk*q)
            read_file: do
                read(afile, *, iostat=ios) fel, fch, fnr
                if (ios /= 0) then
                    write(out,'("load_atom: ab initio density for ",a3,i3,i4," not found")') elem, chg, narad
                    stop "load_atom - ab initio density not found"
                end if
                if ((trim(fel) == trim(elem)) .and. (fch == chg) .and. (fnr == narad)) then
                    do ipt = 1,narad
                        read(afile,*) junk, aden(ipt)
                    end do
                    exit read_file
                end if
            end do read_file
            close(afile)
        case ("slater")
            print '("Loading Slater atomic density for Z = ",i3,", CHG = ",i3)', q, chg
            norb = get_norb(q-chg)
            call psisq(ax, aw, narad, q, chg, norb, niter, thrsh, aden)
        case default
            write(out,'("load_atom: Error unrecognized type",a10)') atyp
            stop "load_atom - unrecognized type"
    end select

end subroutine load_atom

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
            write(out,'("assign_atom: Error unrecognized type",a10)') wtyp
            stop "assign_atom - unrecognized type"
    end select
end function assign_atom

end program vddi
