!
!  The electron density partial charge program
!
!  Calculates charges from a GAMESS natural orbital data file using
!  grid-based methods (Hirshfeld, VDD, Becke, Bader). Methods can be
!  poulations or deformation densities, with radial atomic densities
!  given in an ab initio library. The maximum number of iterations
!  can be specified for deformation density methods.
!
program gridchg
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use utilities
    use atoms
    use fileio
    !
    character(20)            :: pname, wname
    character(2),allocatable :: atypes(:)
    logical,allocatable      :: dload(:,:)
    integer(ik)              :: ofile, mfile, pfile, qfile, wfile
    integer(ik)              :: i, ib, ipt, iter, npts, iat
    integer(ik)              :: nat_count, natom, nbatch, nuniq
    integer(ik),allocatable  :: iwhr(:), qlist(:), npshell(:)
    real(rk)                 :: norm, normatm, tmp_nat_occ(1000)
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), awgt(:,:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:,:,:), pdens(:,:)
    real(ark),allocatable    :: rhomol(:), rhopro(:), rhoatm(:,:), nat_occ(:)
    type(gam_structure)      :: mol, promol
    type(mol_grid)           :: den_grid
    type(input_data)         :: inp

    call accuracyInitialize

    ofile=10
    mfile=11
    pfile=12
    qfile=13
    wfile=14

    open(ofile, file="gridchg.out")
    call read_input(input, inp)
    call init_gridchg_output(ofile, inp)

    write(ofile,'("Loading GAMESS data file")')
    select case (trim(inp%vec_type))
        case ("tr1rdm")
            call gamess_load_rdmsv(trim(inp%orb_fname), tmp_nat_occ, nat_count)
        case ("natorb")
            call gamess_load_natocc(trim(inp%orb_fname), tmp_nat_occ, nat_count)
        case default
            write(out,'("Unrecognized VEC type ",a8)') trim(inp%vec_type)
            stop "bad VEC type"
    end select
    write(ofile,'("Found ",i4," natural occupations")') nat_count
    write(ofile,'("Values are:")')
    write(ofile,'(5(1x,f12.8))') tmp_nat_occ(:nat_count)
    write(ofile,'("")')
    !
    allocate (nat_occ(nat_count))
    nat_occ = tmp_nat_occ(:nat_count)
    !
    call gamess_load_orbitals(file=trim(inp%orb_fname), structure=mol)

    allocate (atypes(mol%natoms), xyzq(4,mol%natoms), iwhr(mol%natoms))
    allocate (charge(mol%natoms), dcharge(mol%natoms), npshell(mol%natoms))

    call gamess_report_nuclei(natom, xyzq, mol)
    write(ofile,'("Molecular geometry (in bohr):")')
    do i = 1,mol%natoms
        atypes(i) = trim(mol%atoms(i)%name)
        write(ofile,'("    ",a2,3f14.8)') atypes(i), xyzq(1:3,i)
    end do
    write(ofile,'("")')
    !
    !  Import the promolecular basis
    !
    promol = mol
    call gamess_reload_basis(file=inp%extbas, basname=inp%atom_bas, structure=promol)
    do i = 1,promol%natoms
        npshell(i) = promol%atoms(i)%nshell
    end do
    !
    !  Set up the grid
    !
    write(ofile,'("Setting up the molecular grid")')
    write(ofile,'("")')
    call GridInitialize(den_grid, inp%n_rad, inp%n_ang, xyzq(1:3,:), atypes)
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    !
    !  Molecular density numerical integration loop
    !
    nullify(xyzw)
    if (trim(inp%weight_type) /= "qtaim") then
        write(ofile,'("Calculating molecular density")')
        open(mfile, file='moldens', form='unformatted', action='write')
        open(pfile, file='dens.mol', action='write')
        norm = 0_rk
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
            call evaluate_density(inp%vec_type, nat_count, npts, mol, nat_occ, xyz, rhomol)
            !
            !  Integrate and save
            !
            mol_integrate: do ipt = 1,npts
                norm = norm + xyzw(4,ipt) * rhomol(ipt)
                write(mfile) rhomol(ipt)
                write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhomol(ipt)
            end do mol_integrate
        end do mol_grid_batches
        close(mfile)
        close(pfile)
        write(ofile,'("Total molecular density: ",f14.8)') norm
        write(ofile,'("")')
        deallocate (rhomol)
    end if

    !
    !  Import atomic properties
    !
    write(ofile,'("Calculating radial atomic grid")')
    write(ofile,'("")')
    call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
    allocate (qlist(nuniq), pdens(sum(npshell),sum(npshell)))
    allocate (dload(2*inp%max_charge+1,nuniq))
    allocate (aden(maxval(npshell),maxval(npshell),2*inp%max_charge+1,nuniq))
    qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)

    !
    !  Atomic density numerical integration loop
    !
    write(ofile,'("Starting atomic density evaluation")')
    write(ofile,'("")')
    open(mfile, file='moldens', form='unformatted', action='read')
    charge = 0_rk
    dload(:,:) = .false.
    iterate_chg: do iter = 1,inp%max_iter
        if (inp%atom_type == "pro") write(ofile,'("ITER = "i4)') iter
        if (iter < 10) then
            write(pname,'("dens.atm.",i1)') iter
            write(wname,'("wgts.atm.",i1)') iter
        else if (iter < 100) then
            write(pname,'("dens.atm.",i2)') iter
            write(wname,'("wgts.atm.",i2)') iter
        else
            write(pname,'("dens.atm.",i3)') iter
            write(wname,'("wgts.atm.",i3)') iter
        end if
        open(qfile, file=trim(pname), action='write')
        open(wfile, file=trim(wname), action='write')
        rewind mfile
        !
        !  Set up the radial atomic densities (only pro or hirsh?)
        !
        call update_densmat(charge, inp%max_charge, qlist, natom, iwhr, nuniq, aden, dload, npshell, inp%lib_dmat, inp%atom_bas, pdens)
        call GridPointsBatch(den_grid, 'Reset')
        iat = den_grid%next_part
        normatm = 0_rk
        dcharge = 0_rk
        !nullify(xyzw)
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
                if (size(rhomol) /= npts) deallocate (rhomol, rhopro, rhoatm, awgt)
            end if
            if (.not. allocated(rhomol)) allocate (rhomol(npts),rhopro(npts),rhoatm(natom,npts),awgt(natom,npts))
            !
            !  Evaluate atomic densities at grid points if necessary
            !
            call evaluate_atomic(xyzw(1:3,:), npts, promol, pdens, npshell, natom, rhoatm)
            rhopro = sum(rhoatm, dim=1)
            !
            !  Read in molecular densities and determine atomic assignments
            !
            if (trim(inp%weight_type) == "qtaim") then
                awgt = 0_rk
                read_rho_qtaim: do ipt = 1,npts
                    read(mfile) rhomol(ipt), iat
                    awgt(iat,ipt) = 1_rk
                end do read_rho_qtaim
            else
                read_rho: do ipt = 1,npts
                    read(mfile) rhomol(ipt)
                end do read_rho
                awgt = assign_atom(xyzq(1:3,:), natom, xyzw(1:3,:), npts, rhoatm, iat, inp%weight_type)
            end if
            !
            !  Find the atomic contributions and output densities
            !
            normatm = normatm + sum(xyzw(4,:) * rhopro)
            if (inp%atom_type == "pro") then
                atom_contrib_pro: do i = 1,natom
                    dcharge(i) = dcharge(i) + sum(awgt(i,:) * xyzw(4,:) * (rhopro - rhomol))
                end do atom_contrib_pro
                integrate_pro: do ipt = 1,npts
                    write(qfile,'(*(e16.8))') rhoatm(:,ipt)
                    write(wfile,'(*(e16.8))') awgt(:,ipt)
                end do integrate_pro
            else
                atom_contrib_pop: do i = 1,natom
                    dcharge(i) = dcharge(i) - sum(awgt(i,:) * xyzw(4,:) * rhomol)
                end do atom_contrib_pop
                integrate_pop: do ipt = 1,npts
                    write(wfile,'(*(e16.8))') awgt(:,ipt)
                end do integrate_pop
            end if
            ! get iatom for the next shell
            iat = den_grid%next_part
        end do grid_batches
        close(qfile)
        close(wfile)
        charge = charge + dcharge

        write(ofile,'("Total atomic density: ",f14.8)') normatm ! only for pro?
        write(ofile,'("Total change in charge: ",f14.8)') sum(dcharge)
        write(ofile,'("Contributions:")')
        write(ofile,'("    ",5f14.8)') dcharge
        write(ofile,'("")')

        do i = 1,natom
            if (charge(i) > qlist(iwhr(i))) then
                write(ofile,'("iterate_chg: Charge of ",f8.3," for atom ",i3," exceeds nuclear charge")') charge(i), i
                stop "iterate_chg - atomic charge exceeds nuclear charge"
            end if
        end do

        if (inp%atom_type == "pop") then
            write(ofile,'("Final charges:")')
            write(ofile,'("    ",5f14.8)') charge + xyzq(4,:)
            exit
        else if (maxval(abs(dcharge)) < inp%chg_thresh) then
            write(ofile,'("Charges converged on ITER = ",i4)') iter
            write(ofile,'("Final charges:")')
            write(ofile,'("    ",5f14.8)') charge
            exit
        else if (iter == inp%max_iter) then
            write(ofile,'("Charges not converged for MAXITER = ",i4)') inp%max_iter
            write(ofile,'("Unconverged charges:")')
            write(ofile,'("    ",5f14.8)') charge
        end if
    end do iterate_chg
    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    close(mfile)
    call system("rm moldens")
    ! only keep data from the last iteration
    call system("mv "//pname//" dens.atm")
    call system("mv "//wname//" wgts.atm")
    call system("rm -f dens.atm.[0-9]* wgts.atm.[0-9]*")
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')
    close(ofile)

    1000 format(e24.8,f16.8,f16.8,f16.8,e20.10)

contains

!
!  Determine the atomic weight factors on a grid
!
function assign_atom(xyzatm, natm, xyz, npt, rho, iatom, wtyp)
    character(5) :: wtyp
    integer(ik)  :: natm, npt, iatom
    real(rk)     :: xyzatm(3,natm), xyz(3,npt), rho(natm,npt)
    real(rk)     :: assign_atom(natm,npt)
    !
    integer(ik)  :: ipt, imin, iat
    real(rk)     :: dist(natm), rhotot(npt)

    select case (wtyp)
        case ("becke")
            assign_atom = 0_rk
            assign_atom(iatom,:) = 1_rk
        case ("voron")
            assign_atom = 0_rk
            do ipt = 1,npt
                dist = sqrt((xyz(1,ipt)-xyzatm(1,:))**2 + &
                            (xyz(2,ipt)-xyzatm(2,:))**2 + &
                            (xyz(3,ipt)-xyzatm(3,:))**2)
                imin = 1
                find_vmin: do iat = 2,natm
                    if (dist(iat) < dist(imin)) then
                        imin = iat
                    end if
                end do find_vmin
                assign_atom(imin,ipt) = 1_rk
            end do
        case ("hirsh")
            rhotot = sum(rho, dim=1)
            assign_atom = 0_rk
            do ipt = 1,npt
                if (rhotot(ipt) > 0) then
                    do iat = 1,natm
                        assign_atom(iat,ipt) = rho(iat,ipt) / rhotot(ipt)
                    end do
                end if
            end do
        case ("qtaim")
            stop "assign_atom - qtaim not supported yet"
        case default
            write(ofile,'("assign_atom: Error unrecognized type",a10)') wtyp
            stop "assign_atom - unrecognized type"
    end select

end function assign_atom

end program gridchg
