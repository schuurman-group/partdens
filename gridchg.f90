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
    use atoms
    use atomdens
    use fileio
    !
    character(20)            :: pname, wname
    character(100)           :: extbas
    character(2),allocatable :: atypes(:)
    logical,allocatable      :: dload(:,:)
    integer(ik)              :: ofile, mfile, pfile, qfile, wfile
    integer(ik)              :: i, j, ib, ipt, iter, npts, iat
    integer(ik)              :: nat_count, natom, nbatch, nuniq !, norb
    integer(ik),allocatable  :: iwhr(:), qlist(:), npbas(:)
    real(rk)                 :: norm, normatm, tmp_nat_occ(1000)
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), awgt(:,:) !, ax(:,:), aw(:,:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:,:,:), pdens(:,:) !, acden(:,:,:)
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

    extbas = '/globalhome/rymac/Projects/PartialCharge/partdens/atomlib/extbas'

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
    allocate (charge(mol%natoms), dcharge(mol%natoms), npbas(mol%natoms))

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
    call gamess_reload_basis(file=extbas, basname=inp%abas, structure=promol)
    do i = 1,promol%natoms
        npbas(i) = 0
        do j = 1,promol%atoms(i)%nshell
            npbas(i) = npbas(i) + gam_orbcnt(promol%atoms(i)%sh_l(j))
        end do
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
    allocate (qlist(nuniq), pdens(mol%nbasis,mol%nbasis)) !, ax(nuniq,inp%n_rad_atom), aw(nuniq,inp%n_rad_atom))
    allocate (dload(2*inp%max_charge+1,nuniq), aden(maxval(npbas),maxval(npbas),2*inp%max_charge+1,nuniq))
    !allocate (aload(2*inp%max_charge+1,nuniq), acden(2*inp%max_charge+1,nuniq,inp%n_rad_atom), aden(natom,inp%n_rad_atom))
    qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
    !setup_rad: do i = 1,nuniq
    !    call rlegendre(inp%n_rad_atom, qlist(i), ax(i,:), aw(i,:))
    !end do setup_rad

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
        !call update_atoms(charge, inp%max_charge, qlist, natom, iwhr, nuniq, inp%n_rad_atom, inp%atom_lib, inp%atom_type, aload, ax, aw, acden, aden)
        call update_densmat(charge, inp%max_charge, qlist, natom, iwhr, nuniq, aden, dload, npbas, pdens)
        print *,pdens
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
            !call evaluate_atomic(aden, ax, inp%n_rad_atom, xyzq(1:3,:), iwhr, nuniq, natom, inp%interp_type, inp%interp_ord, xyzw(1:3,:), npts, rhoatm)
            call evaluate_atomic_new(xyzw(1:3,:), npts, promol, pdens, npbas, natom, rhoatm)
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

        do i = 1,natom
            if (charge(i) > qlist(iwhr(i))) then
                write(ofile,'("iterate_chg: Charge of ",f8.3," for atom ",i3," exceeds nuclear charge")') charge(i), i
                stop "iterate_chg - atomic charge exceeds nuclear charge"
            end if
        end do

        write(ofile,'("Total atomic density: ",f14.8)') normatm ! only for pro?
        write(ofile,'("Total change in charge: ",f14.8)') sum(dcharge)
        write(ofile,'("Contributions:")')
        write(ofile,'("    ",5f14.8)') dcharge
        write(ofile,'("")')

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
!  Determine the molecular density at XYZ coordinates
!
subroutine evaluate_density(vtyp, nat_count, npt, mol, nat_occ, xyz, rho)
    character(6),intent(in)   :: vtyp
    integer(ik),intent(in)    :: nat_count, npt
    type(gam_structure),intent(inout)   :: mol
    real(rk),intent(in)       :: nat_occ(nat_count)
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
            evaluate_rdm: do ird=1,nat_count
                imo = 2*ird - 1
                rho = rho + nat_occ(ird) * moval(imo,:) * moval(imo+1,:)
            end do evaluate_rdm
        case ("natorb")
            evaluate_nat: do ird=1,nat_count
                rho = rho + nat_occ(ird) * moval(ird,:)**2
            end do evaluate_nat
        case default
            write(out,'("evaluate_density: Unrecognized VEC type ",a8)') vtyp
            stop "evaluate_density - bad VEC type"
    end select
    !
    deallocate (basval,moval)

end subroutine evaluate_density

!!
!!  Determine the atomic densities at grid XYZ coordinates
!!
!subroutine evaluate_atomic(rden, r, nr, xyzc, uind, nu, natm, ityp, iord, xyz, npt, rho)
!    character(3) :: ityp
!    integer(ik)  :: nr, natm, npt, nu, uind(natm), iord
!    real(rk)     :: r(nu,nr), xyzc(3,natm), xyz(3,npt)
!    real(ark)    :: rden(natm,nr), rho(natm,npt)
!    !
!    integer(ik)  :: iat, ipt, ui
!    real(rk)     :: ratm
!
!    rho(:,:) = 0_ark
!    evaluate_atom: do iat = 1,natm
!        ui = uind(iat)
!        interp_point: do ipt = 1,npt
!            ratm = sqrt((xyz(1,ipt)-xyzc(1,iat))**2 + &
!                        (xyz(2,ipt)-xyzc(2,iat))**2 + &
!                        (xyz(3,ipt)-xyzc(3,iat))**2)
!            ! atomic density = 0 outside the grid
!            if (ratm < r(ui,nr)) then
!                rho(iat,ipt) = rho(iat,ipt) + interp(ratm, r(ui,:), rden(iat,:), nr, ityp, iord)
!            end if
!        end do interp_point
!    end do evaluate_atom
!
!end subroutine evaluate_atomic

!
!  Determine the analytic atomic densities at grid XYZ coordinates
!
subroutine evaluate_atomic_new(xyz, npt, pmol, dmat, nbas, natm, rho)
    integer(ik),intent(in)            :: natm, npt
    integer(ik),intent(in)            :: nbas(natm)
    type(gam_structure),intent(inout) :: pmol
    real(rk),intent(in)               :: xyz(3,npt)
    real(ark),intent(in)              :: dmat(sum(nbas),sum(nbas),natm)
    real(ark),intent(out)             :: rho(natm,npt)
    !
    integer(ik)           :: iat, ipt, ibas
    real(ark)             :: basval(1,sum(nbas),npt)
    real(ark),allocatable :: ptval(:), admat(:,:)

    ! evaluate basis functions
    evaluate_basis_functions: do ipt = 1,npt
        call gamess_evaluate_functions(xyz(:,ipt), basval(:,:,ipt), pmol)
    end do evaluate_basis_functions
    rho(:,:) = 0_ark
    ib = 1
    evaluate_atom_new: do iat = 1,natm
        ! divide basval into atomic components
        allocate (ptval(nbas(iat)), admat(nbas(iat),nbas(iat)))
        admat = dmat(ibas:ibas+nbas(iat),ibas:ibas+nbas(iat),iat)
        ! multiply by their respective density matrices
        evaluate_pt_dens: do ipt = 1,npt
            ptval = basval(1,ibas:ibas+nbas(iat),ipt)
            rho(iat,ipt) = dot_product(ptval, matmul(admat, ptval))
        end do evaluate_pt_dens
        deallocate (ptval, admat)
        ibas = ibas + nbas(iat)
    end do evaluate_atom_new

end subroutine evaluate_atomic_new

!!
!!  Update the atomic density functions as required
!!
!subroutine update_atoms(chg, maxc, ql, natm, uind, nu, narad, alib, atyp, aload, ax, aw, acden, aden)
!    character(*) :: alib, atyp
!    integer(ik)  :: natm, nu, narad, maxc, ql(nu), uind(natm)
!    logical      :: aload(2*maxc+1,nu)
!    real(rk)     :: chg(natm), ax(nu,narad), aw(nu,narad)
!    real(ark)    :: acden(2*maxc+1,nu,narad), aden(natom,narad)
!    !
!    integer(ik)  :: ia, il, iu, ui
!    real(rk)     :: modc, cthrsh=1e-3
!
!    update_aden: do ia = 1,natm
!        if (abs(chg(ia)) > maxc) then
!            write(ofile,'("update_atoms: abs(charge) greater than MAXCHG = ",i3)') maxc
!            stop "update_atoms - charge out of bounds"
!        end if
!        ui = uind(ia)
!        ! find the modulus (as it should be defined)
!        modc = chg(ia) - floor(chg(ia))
!        if (abs(modc) < cthrsh) then
!            ! integer charge, only get one contribution
!            il = maxc + 1 + nint(chg(ia))
!            if (.not. aload(il,ui)) then
!                ! import the integer density
!                call load_atom(ql(ui), nint(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il,ui,:))
!                aload(il,ui) = .true.
!            end if
!            aden(ia,:) = acden(il,ui,:)
!        else
!            ! real charge, get ceil(chg) and floor(chg) contributions
!            il = maxc + 1 + floor(chg(ia))
!            iu = maxc + 1 + ceiling(chg(ia))
!            if (.not. aload(il,ui)) then
!                ! import the lower integer density
!                call load_atom(ql(ui), floor(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il,ui,:))
!                aload(il,ui) = .true.
!            end if
!            if (.not. aload(iu,ui)) then
!                ! import the upper integer density
!                call load_atom(ql(ui), ceiling(chg(ia)), narad, alib, atyp, ax(ui,:), aw(ui,:), acden(iu,ui,:))
!                aload(iu,ui) = .true.
!            end if
!            aden(ia,:) = (1-modc)*acden(il,ui,:) + modc*acden(iu,ui,:)
!        end if
!    end do update_aden
!
!end subroutine update_atoms

!
!  Update the set of atomic density matrices based on charges
!
subroutine update_densmat(chg, maxc, ql, natm, uind, nu, atmden, dload, nbas, dens)
    integer(ik),intent(in)  :: natm, nu, maxc
    integer(ik),intent(in)  :: nbas(natm), ql(nu), uind(natm)
    logical,intent(inout)   :: dload(2*maxc+1,nu)
    real(rk),intent(in)     :: chg(natm)
    real(ark),intent(inout) :: atmden(maxval(nbas),maxval(nbas),2*maxc+1,nu)
    real(ark),intent(out)   :: dens(sum(nbas),sum(nbas))
    !
    integer(ik)             :: iat, ibas, il, iu, ui
    real(rk)                :: modc, cthrsh=1e-3

    ibas = 0_ik
    update_dmat: do iat = 1,natm
        if (abs(chg(iat)) > maxc) then
            write(ofile,'("update_densmat: abs(charge) greater than MAXCHG = ",i3)') maxc
            stop "update_densmat - charge out of bounds"
        end if
        ui = uind(iat)
        ! find the modulus (as it should be defined)
        modc = chg(iat) - floor(chg(iat))
        if (abs(modc) < cthrsh) then
            ! integer charge, only get one contribution
            il = maxc + 1 + nint(chg(iat))
            if (.not. dload(il,ui)) then
                ! import the integer density matrix
                call load_atom_dmat(ql(ui), nint(chg(iat)), maxval(nbas), atmden(:,:,il,ui))
                dload(il,ui) = .true.
            end if
            dens(ibas:ibas+nbas(iat),ibas:ibas+nbas(iat)) = atmden(:,:,il,ui)
        else
            ! real charge, get ceil(chg) and floor(chg) contributions
            il = maxc + 1 + floor(chg(iat))
            iu = maxc + 1 + ceiling(chg(iat))
            if (.not. dload(il,ui)) then
                ! import the lower integer density matrix
                call load_atom_dmat(ql(ui), floor(chg(iat)), maxval(nbas), atmden(:,:,il,ui))
                dload(il,ui) = .true.
            end if
            if (.not. dload(iu,ui)) then
                ! import the upper integer density matrix
                call load_atom_dmat(ql(ui), ceiling(chg(iat)), maxval(nbas), atmden(:,:,iu,ui))
                dload(iu,ui) = .true.
            end if
            dens(ibas:ibas+nbas(iat),ibas:ibas+nbas(iat)) = (1-modc)*atmden(:,:,il,ui) + modc*atmden(:,:,iu,ui)
        end if
        ibas = ibas + nbas(iat)
    end do update_dmat

end subroutine update_densmat

!!
!!  Load an atomic density
!!
!subroutine load_atom(q, chg, narad, alib, atyp, ax, aw, aden)
!    character(*) :: alib, atyp
!    integer(ik)  :: narad, q, chg
!    real(rk)     :: ax(narad), aw(narad)
!    real(ark)    :: aden(narad)
!    !
!    character(2) :: elem, fel
!    integer(ik)  :: ipt, fch, fnr, ios, afile=111, niter=50
!    real(rk)     :: thrsh=1e-6, junk
!    logical      :: file_exists
!
!    if (chg >= q) then
!        aden = 0_ark
!        return
!    end if
!    select case (atyp)
!        case ("slater")
!            write(ofile,'("Loading Slater atomic density for Z = ",i3,", CHG = ",i3)') q, chg
!            norb = get_norb(q-chg)
!            call psisq(ax, aw, narad, q, chg, norb, niter, thrsh, aden)
!        case default
!            inquire(file=trim(inp%lib_path)//alib, exist=file_exists)
!            if (.not. file_exists) then
!                write(ofile,'("load_atom: Error unrecognized type",a10)') atyp
!                stop "load_atom - unrecognized type"
!            end if
!            write(ofile,'("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)') q, chg
!            open(afile, file=trim(inp%lib_path)//alib)
!            elem = AtomElementSymbol(1._rk*q)
!            read_file: do
!                read(afile, *, iostat=ios) fel, fch, fnr
!                if (ios /= 0) then
!                    write(ofile,'("load_atom: ab initio density for ",a3,i3,i4," not found")') elem, chg, narad
!                    stop "load_atom - ab initio density not found"
!                end if
!                if ((trim(fel) == trim(elem)) .and. (fch == chg) .and. (fnr == narad)) then
!                    do ipt = 1,narad
!                        read(afile,*) junk, aden(ipt)
!                    end do
!                    exit read_file
!                end if
!            end do read_file
!            close(afile)
!    end select
!
!end subroutine load_atom

!
!  Load the density of an atom/basis/charge combination
!
subroutine load_atom_dmat(q, chg, maxbas, admat)
    integer(ik),intent(in) :: q, chg, maxbas
    real(ark),intent(out)  :: admat(maxbas,maxbas)
    !
    character(20)          :: fbas, sstr
    character(2)           :: elem, fel
    integer(ik)            :: i, j, ios, lstr, fch, fnb, ieq, iao, nao
    integer(ik)            :: row, line, nline, ifield !, afile=111
    integer(ik)            :: chk_row, chk_line, loc_int(5), neq(maxbas)
    logical                :: file_exists, have_dmat
    real(ark)              :: offd(maxbas), adtmp(maxbas,maxbas), loc_val(5)
    real(ark),allocatable  :: tmat(:,:), omat(:,:)

    if (chg >= q) then
        admat = 0_ark
        return
    end if
    inquire(file=trim(inp%lib_path), exist=file_exists)
    if (.not. file_exists) then
        write(ofile,'("load_atom_dmat: Density matrix library not found at ",a)') inp%lib_path
        stop "load_atom_dmat - density matrix file does not exist"
    end if
    write(ofile,'("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)') q, chg
    open (gam_file,file=trim(inp%lib_path),action='read',position='rewind',status='old',iostat=ios)
    elem = AtomElementSymbol(1._rk*q)
    ! find the density matrix of the atom/basis/charge
    if (chg > 0.5) then
        write (sstr,"(a2,' ',a,' ',sp,i3)") adjustl(elem), trim(inp%abas), chg
    else
        write (sstr,"(a2,' ',a,' ',i3)") adjustl(elem), trim(inp%abas), chg
    end if
    lstr = len(trim(sstr))
    have_dmat = .false.
    scan_dmat: do while (.not.have_dmat)
        call gam_readline
        if (gam_line_buf(:lstr)==trim(sstr)) have_dmat = .true.
    end do scan_dmat
    if (.not.have_dmat) then
        write(out,'("load_atom_dmat: Density matrix ",a," not found for atom ",a,"with charge ",sp,i3)') trim(inp%abas), trim(elem), chg
        stop "load_atom_dmat - density matrix not found"
    end if
    ! read the number of blocks, fnb
    read (gam_line_buf,*,iostat=ios) fel, fbas, fch, fnb
    nline = ceiling(fnb / 5._rk)
    ! read the number of each symmetry block
    ieq = 1
    line = 0
    read_nequiv: do
      line = line + 1
      call gam_readline
      read (gam_line_buf,'(i5,5i15)',iostat=ios) chk_line, loc_int
      stuff_ne: do ifield = 1,5
          neq(ieq) = loc_int(ifield)
          ieq = ieq + 1
          if (ieq > fnb) exit stuff_ne
      end do stuff_ne
      if (line == nline) exit read_nequiv
    end do read_nequiv
    ! read the off diagonal elements of symmetry blocks
    ieq = 1
    line = 0
    read_offdiag: do
      line = line + 1
      call gam_readline
      read (gam_line_buf,'(i5,5g15.10)',iostat=ios) chk_line, loc_val
      stuff_od: do ifield = 1,5
          offd(ieq) = loc_val(ifield)
          ieq = ieq + 1
          if (ieq > fnb) exit stuff_od
      end do stuff_od
      if (line == nline) exit read_offdiag
    end do read_offdiag
    ! finally, read the density matrix
    row = 1
    ieq = 1
    line = 0
    read_dmat: do
      line = line + 1
      call gam_readline
      read (gam_line_buf,'(i2,i3,5g15.10)',iostat=ios) chk_row, chk_line, loc_val
      if (ios /= 0) then
          write (out,'("load_atom_dmat: format error ",i8)') ios
          stop 'load_atom_dmat - format error'
      end if
      stuff_d: do ifield = 1,5
          admat(ieq,row) = loc_val(ifield)
          ieq = ieq + 1
          if (ieq > fnb) then
              row = row + 1
              if (row > fnb) exit read_dmat
              line = 0
              ieq = 1
              exit stuff_d
          end if
      end do stuff_d
    end do read_dmat
    close (gam_file,iostat=ios)
    ! construct the full density matrix
    nao = sum(neq(:fnb))
    allocate (tmat(fnb,nao), omat(nao,nao))
    iao = 1
    tmat = 0_ark
    omat = 0_ark
    build_tmp: do ieq = 1,fnb
        do i = 0,neq(ieq)-1
            tmat(ieq,iao+i) = 1_ark
            adtmp(iao+i,:fnb) = admat(ieq,:fnb)
            do j = 0,i-1
                omat(iao+i,iao+j) = offd(ieq) - admat(ieq,ieq)
                omat(iao+j,iao+i) = omat(iao+i,iao+j)
            end do
        end do
        iao = iao + neq(ieq)
    end do build_tmp
    iao = 1
    build_final: do ieq = 1,fnb
        do i = 0,neq(ieq)-1
            admat(:nao,iao+i) = adtmp(:nao,ieq)
        end do
        iao = iao + neq(ieq)
    end do build_final
    admat(:nao,:nao) = admat(:nao,:nao) + omat

end subroutine load_atom_dmat

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
