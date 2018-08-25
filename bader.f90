!
!  The (hopefully temporary) Bader charge algorithm
!
!  Calculates charges from a GAMESS RDM output file using the Bader
!  atoms-in-molecules algorithm which requires density gradients. This
!  can potentially be added to the DDCharge code if there is a good way
!  to cache unassigned points (i.e. save x, y, z, w and rho and propagate
!  the assignment until assigned).
!
program bader
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use atoms
    use atomdens
    !
    character(100),parameter :: infile="ddcharge.inp", outfile="ddcharge.out"
    character(100)           :: rdm_file, vectyp, atyp, alib, wtyp, ityp
    character(20)            :: pname, wname
    character(2),allocatable :: atypes(:)
    logical                  :: first_shell, last, ls
    logical,allocatable      :: aload(:,:), adone(:,:)
    integer(ik)              :: ofile, mfile, pfile, wfile
    integer(ik)              :: i, ib, ipt, iter, npts, nax
    integer(ik)              :: rdm_count, natom, nbatch, nuniq, norb
    integer(ik)              :: nrad, nang, narad, maxiter, iordr, maxchg
    integer(ik),allocatable  :: iwhr(:), qlist(:), nnn(:), nmask(:,:), asgn(:,:)
    real(rk)                 :: norm, normatm, tmp_rdm_sv(1000), thrsh, dist
    real(rk)                 :: avg_now(3), avg_nxt(3)
    real(rk),pointer         :: xyzw_lst(:,:), xyzw_now(:,:), xyzw_nxt(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), ax(:,:), aw(:,:), awgt(:,:), pt_dist(:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:), acden(:,:,:)
    real(ark),allocatable    :: rhomol(:,:), drhomol(:,:,:), rhopro(:), rhoatm(:,:), rdm_sv(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid

    call accuracyInitialize

    ofile=10
    mfile=11
    pfile=12
    wfile=13

    open(ofile, file=outfile)
    call read_input(infile, rdm_file, vectyp, nrad, nang, atyp, alib, narad, ityp, iordr, wtyp, maxchg, maxiter, thrsh)
    call init_output(rdm_file, vectyp, nrad, nang, atyp, alib, narad, ityp, iordr, wtyp, maxchg, maxiter, thrsh)

    write(ofile,'("Loading GAMESS RDM file")')
    call gamess_load_rdmsv(trim(rdm_file), tmp_rdm_sv, rdm_count)
    write(ofile,'("Found ",i4," singular values")') rdm_count
    write(ofile,'("Values are:")')
    write(ofile,'(5(1x,f12.8))') tmp_rdm_sv(:rdm_count)
    write(ofile,'("")')
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
    write(ofile,'("Molecular geometry (in bohr):")')
    do i = 1,mol%natoms
        write(ofile,'("    ",a2,3f14.8)') atypes(i), xyzq(1:3,i)
    end do
    write(ofile,'("")')
    !
    !  Set up the grid
    !
    write(ofile,'("Setting up the molecular grid")')
    write(ofile,'("")')
    call GridInitialize(den_grid, nrad, nang, xyzq(1:3,:), atypes)
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    !!open(mfile, file='moldens', form='unformatted', action='write')
    !
    !  Set initial values
    !
    norm = 0_rk
    charge(:) = 0_rk
    nullify(xyzw_lst, xyzw_now, xyzw_nxt)
    call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_now)
    npts = size(xyzw_now, dim=2)
    allocate (rhomol(2,npts), drhomol(2,3,npts))
    allocate (nnn(npts), nmask(7,npts), pt_dist(npts))
    allocate (asgn(2,npts), adone(2,npts))
    !
    !  Find nearest neighbours
    !
    find_nearest: do ipt = 1, npts
        nax = 0
        do i = 1, 3
            if (abs(xyzw_now(i,ipt)) < 1e-10) nax = nax + 1
        end do
        if (nax == 2) then
            nnn(ipt) = 5
        else
            nnn(ipt) = 7
        end if
        pt_dist = sqrt((xyzw_now(1,:) - xyzw_now(1,ipt))**2 + &
                       (xyzw_now(2,:) - xyzw_now(2,ipt))**2 + &
                       (xyzw_now(3,:) - xyzw_now(3,ipt))**2)
        nmask(:,ipt) = argslow(nnn(ipt), pt_dist, npts)
    end do find_nearest
    !
    !  Molecular density numerical integration loop
    !
    write(ofile,'("Calculating molecular density")')
    asgn(:,:) = 0
    adone(:,:) = .false.
    first_shell = .true.
    mol_grid_batches: do ib = 1, nbatch-1
        if (first_shell) then
            !
            !  Get next grid points
            !
            call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_nxt)
            first_shell = .false.
            last = .false.
            ls = .false.
        else
            if (ib == nbatch-1) then
                last = .true.
                ls = .true.
            else
                call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_nxt)
                !!npts = size(xyzw_now, dim=2)
                avg_now = sum(xyzw_now(1:3,:), dim=2)
                avg_nxt = sum(xyzw_nxt(1:3,:), dim=2)
                dist = sqrt(sum((avg_nxt - avg_now)**2, dim=1))
                if (dist > 0.2) then
                    ls = .true.
                else
                    ls = .false.
                end if
            end if
            xyz => xyzw_now(1:3,:)
            !!!
            !!!  If the batch size changed, reallocate rho() ! These probably should never change...
            !!!
            !!npts = size(xyzw, dim=2)
            !!if (allocated(rhomol)) then
            !!    if (size(rhomol) /= npts) deallocate (rhomol)
            !!    if (size(drhomol) /= npts) deallocate (drhomol)
            !!end if
            !!if (.not. allocated(rhomol)) allocate (rhomol(npts))
            !!if (.not. allocated(drhomol)) allocate (drhomol(npts))
            !
            !  Evaluate density and density gradients at grid points
            !
            call evaluate_density_gradients(vectyp, rdm_count, npts, mol, rdm_sv, xyz, rhomol(2,:), drhomol(2,:,:))
            norm = norm + sum(rhomol(2,:) * xyzw_now(4,:))
            !
            !  Assign the densities for one shell
            !
            rhomol(2,:) = rhomol(2,:) * xyzw_now(4,:)
            call assign_shell(xyzw_lst, xyzw_now, rhomol, drhomol, asgn, adone, npts, xyzq, charge, natom, nmask, nnn, ls)
        end if
        !
        !  Prepare the next shell
        !
        if (.not. last) then
            xyzw_lst => xyzw_now
            xyzw_now => xyzw_nxt
            rhomol(1,:) = rhomol(2,:)
            drhomol(1,:,:) = drhomol(2,:,:)
            asgn(1,:) = asgn(2,:)
            asgn(2,:) = 0
            adone(1,:) = adone(2,:)
            adone(2,:) = .false.
            if (ls) then
                first_shell = .true.
            end if
        end if

        !!!
        !!!  Integrate and save
        !!!
        !!mol_integrate: do ipt = 1,npts
        !!    norm = norm + xyzw(4,ipt) * rhomol(ipt)
        !!    !!write(mfile) rhomol(ipt)
        !!end do mol_integrate
    end do mol_grid_batches
    !!close(mfile)
    write(ofile,'("Total molecular density: ",f14.8)') norm
    write(ofile,'("")')
    charge = xyzq(4,:) - charge
    write(ofile,'("Total charge: ",f14.8)') sum(charge)
    write(ofile,'("Final charges:")')
    write(ofile,'("    ",5f14.8)') charge

    !!!
    !!!  Import atomic properties
    !!!
    !!write(ofile,'("Calculating radial atomic grid")')
    !!write(ofile,'("")')
    !!call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
    !!allocate (qlist(nuniq), ax(nuniq,narad), aw(nuniq,narad))
    !!allocate (aload(2*maxchg+1,nuniq), acden(2*maxchg+1,nuniq,narad), aden(natom,narad))
    !!qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
    !!setup_rad: do i = 1,nuniq
    !!    call rlegendre(narad, qlist(i), ax(i,:), aw(i,:))
    !!end do setup_rad

    !!!
    !!!  Atomic density numerical integration loop
    !!!
    !!write(ofile,'("Starting atomic density evaluation")')
    !!write(ofile,'("")')
    !!open(mfile, file='moldens', form='unformatted', action='read')
    !!charge = 0_rk
    !!aload(:,:) = .false.
    !!iterate_chg: do iter = 1,maxiter
    !!    write(ofile,'("ITER = "i4)') iter
    !!    if (iter < 10) then
    !!        write(pname,'("density.",i1)') iter
    !!        write(wname,'("weights.",i1)') iter
    !!    else if (iter < 100) then
    !!        write(pname,'("density.",i2)') iter
    !!        write(wname,'("weights.",i2)') iter
    !!    else
    !!        write(pname,'("density.",i3)') iter
    !!        write(wname,'("weights.",i3)') iter
    !!    end if
    !!    open(pfile, file=trim(pname), action='write')
    !!    open(wfile, file=trim(wname), action='write')
    !!    rewind mfile
    !!    !
    !!    !  Set up the radial atomic densities
    !!    !
    !!    call update_atoms(charge, maxchg, qlist, natom, iwhr, nuniq, narad, alib, atyp, aload, ax, aw, acden, aden)
    !!    call GridPointsBatch(den_grid, 'Reset')
    !!    normatm = 0_rk
    !!    dcharge = 0_rk
    !!    nullify(xyzw)
    !!    grid_batches: do ib = 1,nbatch
    !!        !
    !!        !  Get grid points
    !!        !
    !!        call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
    !!        !
    !!        !  If the batch size changed, reallocate density and weight arrays
    !!        !
    !!        npts = size(xyzw, dim=2)
    !!        if (allocated(rhopro)) then
    !!            if (size(rhomol) /= npts) deallocate (rhomol)
    !!            if (size(rhopro) /= npts) deallocate (rhopro)
    !!            if (size(rhoatm,2) /= npts) deallocate (rhoatm)
    !!            if (size(awgt,2) /= npts) deallocate (awgt)
    !!        end if
    !!        if (.not. allocated(rhomol)) allocate (rhomol(npts))
    !!        if (.not. allocated(rhopro)) allocate (rhopro(npts))
    !!        if (.not. allocated(rhoatm)) allocate (rhoatm(natom,npts))
    !!        if (.not. allocated(awgt)) allocate (awgt(natom,npts))
    !!        !
    !!        !  Read in molecular densities from file
    !!        !
    !!        read_rho: do ipt = 1,npts
    !!            read(mfile) rhomol(ipt)
    !!        end do read_rho
    !!        !
    !!        !  Evaluate atomic densities at grid points
    !!        !
    !!        call evaluate_atomic(aden, ax, narad, xyzq(1:3,:), iwhr, nuniq, natom, ityp, iordr, xyzw(1:3,:), npts, rhoatm)
    !!        rhopro = sum(rhoatm, dim=1)
    !!        !
    !!        !  Determine atomic assignments
    !!        !
    !!        awgt = assign_atom(xyzq(1:3,:), natom, xyzw(1:3,:), npts, rhoatm, wtyp)
    !!        !
    !!        !  Find the atomic contributions
    !!        !
    !!        normatm = normatm + sum(xyzw(4,:) * rhopro)
    !!        atom_contrib: do i = 1,natom
    !!            dcharge(i) = dcharge(i) + sum(awgt(i,:) * xyzw(4,:) * (rhopro - rhomol))
    !!        end do atom_contrib
    !!        !
    !!        !  Output the densities for each point
    !!        !
    !!        integrate: do ipt = 1,npts
    !!            write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhopro(ipt)-rhomol(ipt)
    !!            write(wfile,'(*(e16.8))') awgt(:,ipt)
    !!        end do integrate
    !!    end do grid_batches
    !!    close(pfile)
    !!    close(wfile)
    !!    charge = charge + dcharge

    !!    do i = 1,natom
    !!        if (charge(i) > qlist(iwhr(i))) then
    !!            write(ofile,'("iterate_chg: Charge of ",f8.3," for atom ",i3," exceeds nuclear charge")') charge(i), i
    !!            stop "iterate_chg - atomic charge exceeds nuclear charge"
    !!        end if
    !!    end do

    !!    write(ofile,'("Total atomic density: ",f14.8)') normatm
    !!    write(ofile,'("Total change in charge: ",f14.8)') sum(dcharge)
    !!    write(ofile,'("Contributions:")')
    !!    write(ofile,'("    ",5f14.8)') dcharge
    !!    write(ofile,'("")')

    !!    if (maxval(abs(dcharge)) < thrsh) then
    !!        write(ofile,'("Charges converged on ITER = ",i4)') iter
    !!        write(ofile,'("Final charges:")')
    !!        write(ofile,'("    ",5f14.8)') charge
    !!        exit
    !!    else if (iter == maxiter) then
    !!        write(ofile,'("Charges not converged for MAXITER = ",i4)') maxiter
    !!        write(ofile,'("Unconverged charges:")')
    !!        write(ofile,'("    ",5f14.8)') charge
    !!    end if
    !!end do iterate_chg
    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    !!close(mfile)
    !!call system("rm moldens")
    !!! only keep data from the last iteration
    !!call system("mv "//pname//" density.dat")
    !!call system("mv "//wname//" weights.dat")
    !!call system("rm density.[0-9]* weights.[0-9]*")
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')

    1000 format(e24.8,f16.8,f16.8,f16.8,e20.10)

contains

!
!  Read a crude input file
!
subroutine read_input(infile, rfile, vtyp, nr, na, atyp, alib, nar, ityp, iord, wtyp, maxc, maxi, thr)
    character(100) :: infile, rfile, vtyp, atyp, alib, wtyp, ityp
    integer(ik)    :: nr, na, nar, maxi, iord, maxc
    real(rk)       :: thr
    !
    character(20)  :: varname, var
    integer(ik)    :: ios, nvar

    ios = 0
    nvar = 0
    open(22, file=infile)
    parse_input: do
        read(22,*,iostat=ios) varname, var
        if (ios /= 0) then
            exit parse_input
        end if
        select case (varname)
            case ("rdm_filename")
                rfile = var
                nvar = nvar + 1
            case ("vec_type")
                vtyp = var
                nvar = nvar + 1
            case ("n_r_grid")
                read(var,*,iostat=ios) nr
                nvar = nvar + 1
            case ("n_ang_grid")
                read(var,*,iostat=ios) na
                nvar = nvar + 1
            case ("atom_type")
                atyp = var
                nvar = nvar + 1
            case ("atom_library")
                alib = var
                nvar = nvar + 1
            case ("n_r_atom")
                read(var,*,iostat=ios) nar
                nvar = nvar + 1
            case ("interp_type")
                ityp = var
                nvar = nvar + 1
            case ("interp_ord")
                read(var,*,iostat=ios) iord
                nvar = nvar + 1
            case ("weight_type")
                wtyp = var
                nvar = nvar + 1
            case ("max_charge")
                read(var,*,iostat=ios) maxc
                nvar = nvar + 1
            case ("max_iter")
                read(var,*,iostat=ios) maxi
                nvar = nvar + 1
            case ("chg_thresh")
                read(var,*,iostat=ios) thr
                nvar = nvar + 1
        end select
        if (ios /= 0) then
            write(ofile,'("read_input: Incorrect type for variable ",a)') varname
            stop "read_input - incorrect variable type"
        end if
    end do parse_input
    if (nvar < 13) then
        write(ofile,'("read_input: Missing variable")')
        stop "read_input - missing variable"
    end if

end subroutine read_input

!
!  Initialize the output file
!
subroutine init_output(rfile, vtyp, nr, na, atyp, alib, nar, ityp, iord, wtyp, maxc, maxi, thr)
    character(100) :: rfile, vtyp, atyp, alib, wtyp, ityp
    integer(ik)    :: nr, na, nar, maxi, iord, maxc
    real(rk)       :: thr

    write(ofile,'("+--------------------------------------------------+")')
    write(ofile,'("|                                                  |")')
    write(ofile,'("|                     DDCharge                     |")')
    write(ofile,'("|                                                  |")')
    write(ofile,'("|    Deformation density partial atomic charges    |")')
    write(ofile,'("|         RJ MacDonell, MS Schuurman 2018          |")')
    write(ofile,'("+--------------------------------------------------+")')
    write(ofile,'("")')
    write(ofile,'("")')
    write(ofile,'("Input summary:")')
    write(ofile,'("    ------- Molecular density --------")')
    write(ofile,'("    rdm_filename   =   ",a15)') trim(rfile)
    write(ofile,'("    vec_type       =   ",a15)') trim(vtyp)
    write(ofile,'("    n_r_grid       =   ",i15)') nr
    write(ofile,'("    n_ang_grid     =   ",i15)') na
    write(ofile,'("")')
    write(ofile,'("    --------- Atomic density ---------")')
    write(ofile,'("    atom_type      =   ",a15)') trim(atyp)
    write(ofile,'("    atom_library   =   ",a15)') trim(alib)
    write(ofile,'("    n_r_atom       =   ",i15)') nar
    write(ofile,'("")')
    write(ofile,'("    --------- Interpolation ----------")')
    write(ofile,'("    interp_type    =   ",a15)') trim(ityp)
    write(ofile,'("    interp_ord     =   ",i15)') iord
    write(ofile,'("")')
    write(ofile,'("    ------------- Charge -------------")')
    write(ofile,'("    weight_type    =   ",a15)') trim(wtyp)
    write(ofile,'("    max_charge     =   ",i15)') maxc
    write(ofile,'("    max_iter       =   ",i15)') maxi
    write(ofile,'("    chg_thresh     =   ",e15.3)') thr
    write(ofile,'("")')
end subroutine init_output

!
!  Find unique elements of a list and return a mask
!
subroutine unique(arr, narr, mask, nuniq)
    integer(ik),intent(in)  :: narr, arr(narr)
    integer(ik),intent(out) :: nuniq, mask(narr)
    !
    integer(ik) :: i, j, ident(100)

    nuniq = 0
    mask = 0
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
!  Find the indices of the smallest n elements of a list the slow way
!
!  For small list sizes (< 10^5) and small nsm, the advantages of
!  quicksort are negligible
!
function argslow(nsm, l, nl) result(sm)
    integer(ik) :: nl, nsm, sm(nsm)
    real(rk)    :: l(nl)
    !
    integer(ik) :: i, j

    sm(:) = 0
    find_smallest: do i = 1, nsm
        l_loop: do j = 1, nl
            if (any(sm == j)) cycle l_loop
            if (sm(i) == 0) then
                sm(i) = j
                cycle l_loop
            end if
            if (l(j) < l(sm(i))) sm(i) = j
        end do l_loop
    end do find_smallest
end function argslow

!
!  Determine the molecular density at XYZ coordinates
!
subroutine evaluate_density(vtyp, rdm_count, npt, mol, rdm_sv, xyz, rho)
    character(100),intent(in) :: vtyp
    integer(ik),intent(in)    :: rdm_count, npt
    type(gam_structure),intent(inout)   :: mol
    real(rk),intent(in)       :: rdm_sv(rdm_count)
    real(rk),intent(in)       :: xyz(3,npt)
    real(ark),intent(out)     :: rho(npt)
    !
    integer(ik)               :: ipt, ird, imo, nmo, nbas
    real(ark),allocatable     :: basval(:,:,:)
    real(rk),allocatable      :: moval(:,:)

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
    rho = 0_ark
    select case (vtyp)
        case ("tr1rdm")
            evaluate_rdm: do ird=1,rdm_count
                imo = 2*ird - 1
                rho = rho + rdm_sv(ird) * moval(imo,:) * moval(imo+1,:)
            end do evaluate_rdm
        case ("natorb")
            evaluate_nat: do ird=1,rdm_count
                rho = rho + rdm_sv(ird) * moval(ird,:)**2
            end do evaluate_nat
        case default
            write(out,'("evaluate_density: Unrecognized VEC type ",a8)') vtyp
            stop "evaluate_density - bad VEC type"
    end select
    !
    deallocate (basval,moval)

end subroutine evaluate_density

!
!  Determine the molecular density and gradients at XYZ coordinates
!
subroutine evaluate_density_gradients(vtyp, rdm_count, npt, mol, rdm_sv, xyz, rho, drho)
    character(100),intent(in) :: vtyp
    integer(ik),intent(in)    :: rdm_count, npt
    type(gam_structure),intent(inout)   :: mol
    real(rk),intent(in)       :: rdm_sv(rdm_count)
    real(rk),intent(in)       :: xyz(3,npt)
    real(ark),intent(out)     :: rho(npt), drho(3,npt)
    !
    integer(ik)               :: ipt, ird, imo, nmo, nbas
    real(ark),allocatable     :: basval(:,:,:)
    real(rk),allocatable      :: moval(:,:), dmoval(:,:,:)

    nmo  = mol%nvectors
    nbas = mol%nbasis
    !
    allocate (basval(4,nbas,npt), moval(nmo,npt), dmoval(3,nmo,npt))
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
    dmoval(1,:,:) = matmul(transpose(mol%vectors(:,:nmo)), basval(2,:,:))
    dmoval(2,:,:) = matmul(transpose(mol%vectors(:,:nmo)), basval(3,:,:))
    dmoval(3,:,:) = matmul(transpose(mol%vectors(:,:nmo)), basval(4,:,:))
    !
    !  Finally, evaluate the transition density and gradients at grid points
    !
    rho = 0_ark
    drho = 0_ark
    select case (vtyp)
        case ("tr1rdm")
            evaluate_rdm: do ird=1,rdm_count
                imo = 2*ird - 1
                rho = rho + rdm_sv(ird) * moval(imo,:) * moval(imo+1,:)
                drho = drho + rdm_sv(ird) * dmoval(:,imo,:) * dmoval(:,imo+1,:)
            end do evaluate_rdm
        case ("natorb")
            evaluate_nat: do ird=1,rdm_count
                rho = rho + rdm_sv(ird) * moval(ird,:)**2
                drho = drho + rdm_sv(ird) * dmoval(:,ird,:)**2
            end do evaluate_nat
        case default
            write(out,'("evaluate_density: Unrecognized VEC type ",a8)') vtyp
            stop "evaluate_density - bad VEC type"
    end select
    !
    deallocate (basval,moval,dmoval)

end subroutine evaluate_density_gradients

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

    rho(:,:) = 0_ark
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
    real(rk)     :: chg(natm), ax(nu,narad), aw(nu,narad)
    real(ark)    :: acden(2*maxc+1,nu,narad), aden(natom,narad)
    logical      :: aload(2*maxc+1,nu)
    !
    integer(ik)  :: ia, il, iu, ui
    real(rk)     :: modc, cthrsh=1e-3

    update_aden: do ia = 1,natm
        if (abs(chg(ia)) > maxc) then
            write(ofile,'("update_atoms: abs(charge) greater than MAXCHG = ",i3)') maxc
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

    if (chg >= q) then
        aden = 0_ark
        return
    end if
    select case (atyp)
        case ("abinitio")
            write(ofile,'("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)') q, chg
            open(afile, file="/home/rymac/Projects/PartialCharge/ddcharge/atomlib/"//alib)
            elem = AtomElementSymbol(1._rk*q)
            read_file: do
                read(afile, *, iostat=ios) fel, fch, fnr
                if (ios /= 0) then
                    write(ofile,'("load_atom: ab initio density for ",a3,i3,i4," not found")') elem, chg, narad
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
            write(ofile,'("Loading Slater atomic density for Z = ",i3,", CHG = ",i3)') q, chg
            norb = get_norb(q-chg)
            call psisq(ax, aw, narad, q, chg, norb, niter, thrsh, aden)
        case default
            write(ofile,'("load_atom: Error unrecognized type",a10)') atyp
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
            assign_atom = 0_rk
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
            assign_atom = 0_rk
            do ipt = 1,npt
                if (rhotot(ipt) > 0) then
                    do iatm = 1,natm
                        assign_atom(iatm,ipt) = rho(iatm,ipt) / rhotot(ipt)
                    end do
                end if
            end do
        case default
            write(ofile,'("assign_atom: Error unrecognized type",a10)') wtyp
            stop "assign_atom - unrecognized type"
    end select
end function assign_atom

!
!  Determine the Bader atomic assignments
!
subroutine assign_shell(xyzw_lst, xyzw_now, wrho, drhomol, asgn, adone, npts, xyzq, atmchg, natm, nrst, nnn, ls)
    integer(ik)           :: npts, natm, nnn(npts), asgn(2,npts), nrst(7,npts)
    real(rk)              :: xyzw_lst(4,npts), xyzw_now(4,npts), xyzq(4,natm), atmchg(natm)
    real(ark)             :: wrho(2,npts), drhomol(2,3,npts)
    logical               :: adone(2,npts), askip(2,npts), ls
    !
    integer(ik)           :: i, j, k, n, iout(2,npts), o_ji, a_ji, o_next, a_next, iatm
    integer(ik),parameter :: maxitr=12
    real(rk)              :: dist(natm), din(3)
    logical               :: shell_done

    iout(:,:) = 1
    askip(:,:) = .false.

    make_assignments: do i = 1, npts
        !if (maxval(drhomol(1,:,i)) < 1e-16) then
        if (wrho(1,i) < 1e-36) then
            adone(1,i) = .true.
        else if (asgn(1,i) == 0) then
            call assign_point(drhomol(1,:,i), xyzw_lst, xyzw_now, npts, nrst(:,i), nnn(i), .true., iout(1,i), asgn(1,i))
        end if

        din = xyzw_lst(1:3,i) - xyzw_now(1:3,i)
        !if (maxval(drhomol(2,:,i)) < 1e-16) then
        if (wrho(2,i) < 1e-36) then
            adone(2,i) = .true.
        else if (unit_dot(din, drhomol(2,:,i), 3) < 0 .and. .not. ls) then
            askip(2,i) = .true.
        else
            call assign_point(drhomol(2,:,i), xyzw_lst, xyzw_now, npts, nrst(:,i), nnn(i), .false., iout(2,i), asgn(2,i))
        end if
    end do make_assignments

    propagate_assignments: do k = 1, maxitr
        shell_done = .true.
        loop_shells: do j = 1, 2
            propagate_shell: do i = 1, npts
                if (adone(j,i) .or. askip(j,i)) cycle propagate_shell
                !if (adone(j,i)) then
                !    print *,npts,k,j,i,"done"
                !    print *,iout(j,i)
                !    print *,asgn(j,i)
                !    cycle propagate_shell
                !else if (askip(j,i)) then
                !    print *,npts,k,j,i,"skip"
                !    print *,iout(j,i)
                !    print *,asgn(j,i)
                !    cycle propagate_shell
                !end if
                o_ji = iout(j,i)
                a_ji = asgn(j,i)

                !print *,npts, k, j, i
                !print *,o_ji
                !print *,a_ji
                o_next = iout(o_ji,a_ji)
                a_next = asgn(o_ji,a_ji)
                if (o_next == j .and. a_next == i) then
                    ! repetition, add both to nearest atom density
                    do n = 1, natm
                        dist(n) = sqrt(sum((xyzw_now(1:3,i) - xyzq(1:3,n))**2))
                    end do
                    iatm = minloc(dist, 1)
                    !if (dist(iatm) > 0.5) then raise error?
                    asgn(j,i) = -iatm
                    asgn(o_ji,a_ji) = -iatm
                    atmchg(iatm) = atmchg(iatm) + wrho(j,i) + wrho(o_ji,a_ji)
                    adone(j,i) = .true.
                    adone(o_ji,a_ji) = .true.
                else if (a_next == 0) then
                    if (j == 1) then
                        ! pointed-to index out but not assigned, add to density (caching would happen here)
                        wrho(o_ji,a_ji) = wrho(o_ji,a_ji) + wrho(j,i)
                        adone(j,i) = .true.
                    else
                        askip(j,i) = .true.
                    end if
                else if (a_next < 0) then
                    ! pointed-to index is an atom, add to atomic density
                    asgn(j,i) = a_next
                    iatm = -a_next
                    atmchg(iatm) = atmchg(iatm) + wrho(j,i)
                    adone(j,i) = .true.
                else
                    ! pointed-to index is assigned, keep propagating
                    asgn(j,i) = a_next
                    iout(j,i) = o_next
                end if

                shell_done = shell_done .and. (adone(j,i) .or. askip(j,i))
            end do propagate_shell
        end do loop_shells

        if (shell_done) then
            ! all points going out or assigned to atoms
            exit propagate_assignments
        else if (k == npts) then
            print *,"Shell assignments not finished in maximum number of iterations"
            stop "propagate_assignments: Maximum iterations exceeded"
        end if
    end do propagate_assignments
    !print *,"Made it through one shell!"
end subroutine assign_shell

!
!  Assign a single point to anoter points based on the steepest ascent
!
subroutine assign_point(drho, xyzw_lst, xyzw_now, npts, inrst, nni, outer, iiout, iasgn)
    integer(ik) :: npts, nni, inrst(7), iiout, iasgn
    real(rk)    :: drho(3), xyzw_lst(4,npts), xyzw_now(4,npts)
    logical     :: outer
    !
    integer(ik) :: i
    real(rk)    :: vdot, vdotmx, xyzi(3), xyz0(3,nni), xyz1(3,nni), dv(3)

    ! find position and neighbours on current shell and other shell
    if (outer) then
        xyzi = xyzw_now(1:3,inrst(1))
        do i = 1, nni
            xyz0(:,i) = xyzw_now(1:3,inrst(i))
            xyz1(:,i) = xyzw_lst(1:3,inrst(i))
        end do
    else
        xyzi = xyzw_lst(1:3,inrst(1))
        do i = 1, nni
            xyz0(:,i) = xyzw_lst(1:3,inrst(i))
            xyz1(:,i) = xyzw_now(1:3,inrst(i))
        end do
    end if

    iiout = 0
    vdotmx = -2

    ! check if in current shell
    current_shell: do i = 2, nni
        dv = xyz0(:,i) - xyzi
        vdot = unit_dot(drho, dv, 3)
        if (vdot > vdotmx) then
            iasgn = i
            vdotmx = vdot
            iiout = 1
        end if
    end do current_shell

    ! check if in other shell
    other_shell: do i = 1, nni
        dv = xyz1(:,i) - xyzi
        vdot = unit_dot(drho, dv, 3)
        if (vdot > vdotmx) then
            iasgn = i
            vdotmx = vdot
            iiout = 2
        end if
    end do other_shell

    if (vdotmx == -2) then
        print *,"Assignment not found for current point"
        stop "assign_point: Assignment not found"
    end if
    if (.not. outer) then
        iiout = 3 - iiout
    end if
end subroutine assign_point

!
!  Evaluate the dot product of unit vectors of u and v
!
function unit_dot(u, v, nd)
    integer(ik) :: nd
    real(rk)    :: u(nd), v(nd), unit_dot
    !
    real(rk)    :: norm_u, norm_v

    norm_u = sqrt(sum(u**2))
    norm_v = sqrt(sum(v**2))
    unit_dot = sum((u * v) / (norm_u * norm_v))
end function unit_dot

end program bader
