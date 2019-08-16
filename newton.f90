!
!  The newton fit-density partial charge program
!
!  Calculates charges from a GAMESS natural orbital data file by fitting
!  (charged) atomic densities on a grid to the molecular density.
!
program newton
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use atoms
    use atomdens
    use lapack
    use fileio
    !
    character(20)            :: pname
    character(2),allocatable :: atypes(:)
    logical,allocatable      :: aload(:,:)
    integer(ik)              :: ofile, mfile, pfile, qfile
    integer(ik)              :: i, ib, ipt, iter, npts, iat, jat
    integer(ik)              :: nat_count, natom, nbatch, nuniq, norb
    integer(ik),allocatable  :: iwhr(:), qlist(:)
    real(rk)                 :: norm, normatm, ddsq, qtot, tmp_nat_occ(1000)
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), ax(:,:), aw(:,:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:), agrad(:,:), acden(:,:,:)
    real(ark),allocatable    :: rhomol(:), rhopro(:), rhoatm(:,:), rhoagrad(:,:)
    real(ark),allocatable    :: hess(:,:), rho2rhs(:,:), nat_occ(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid
    type(input_data)         :: inp

    call accuracyInitialize

    ofile=10
    mfile=11
    pfile=12
    qfile=13

    open(ofile, file="newton.out")
    call read_input(input, inp)
    call init_gridchg_output(ofile, inp)
    if (inp%atom_type /= "pro" .or. inp%weight_type /= "newtn") stop "newton.f90 is only for atom_type=pro, weight_type=newtn"

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
    allocate (charge(mol%natoms), dcharge(mol%natoms))

    call gamess_report_nuclei(natom, xyzq, mol)
    write(ofile,'("Molecular geometry (in bohr):")')
    do iat = 1,natom
        atypes(iat) = trim(mol%atoms(iat)%name)
        write(ofile,'("    ",a2,3f14.8)') atypes(iat), xyzq(1:3,iat)
    enddo
    write(ofile,'("")')
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
    write(ofile,'("Calculating molecular density")')
    open(mfile, file='moldens', form='unformatted', action='write')
    open(pfile, file='dens.mol', action='write')
    norm = 0_rk
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
    qtot = sum(xyzq(4,:)) - norm

    !
    !  Import atomic properties
    !
    write(ofile,'("Calculating radial atomic grid")')
    write(ofile,'("")')
    call unique(int(xyzq(4,:)), natom, iwhr, nuniq)
    allocate (qlist(nuniq), ax(nuniq,inp%n_rad_atom), aw(nuniq,inp%n_rad_atom))
    allocate (aload(2*inp%max_charge+1,nuniq), acden(2*inp%max_charge+1,nuniq,inp%n_rad_atom), aden(natom,inp%n_rad_atom))
    allocate (agrad(natom,inp%n_rad_atom), hess(natom+1,natom+1), rho2rhs(natom+1,1))
    qlist = unique_list(int(xyzq(4,:)), natom, iwhr, nuniq)
    setup_rad: do i = 1,nuniq
        call rlegendre(inp%n_rad_atom, qlist(i), ax(i,:), aw(i,:))
    end do setup_rad

    !
    !  Atomic density numerical integration loop
    !
    charge = 0_rk
    aload(:,:) = .false.
    write(ofile,'("Starting atomic density evaluation")')
    write(ofile,'("")')
    open(mfile, file='moldens', form='unformatted', action='read')
    iterate_chg: do iter = 1,inp%max_iter
        if (inp%atom_type == "pro") write(ofile,'("ITER = "i4)') iter
        if (iter < 10) then
            write(pname,'("dens.atm.",i1)') iter
        else if (iter < 100) then
            write(pname,'("dens.atm.",i2)') iter
        else
            write(pname,'("dens.atm.",i3)') iter
        end if
        open(qfile, file=trim(pname), action='write')
        rewind mfile
        call GridPointsBatch(den_grid, 'Reset')
        !
        !  Update the radial atomic densities
        !
        call update_atoms(charge, inp%max_charge, qlist, natom, iwhr, nuniq, inp%n_rad_atom, inp%atom_lib, inp%atom_type, aload, ax, aw, acden, aden)
        update_agrad: do iat = 1,natom
            if (abs(charge(iat) - floor(charge(iat))) < 1e-6) then
                i = inp%max_charge + 1 + nint(charge(iat))
            else
                i = inp%max_charge + 1 + floor(charge(iat))
            end if
            agrad(iat,:) = acden(i+1,iwhr(iat),:) - acden(i,iwhr(iat),:)
        end do update_agrad
        !
        normatm = 0_rk
        ddsq = 0_rk
        hess = 0_ark
        hess(natom+1,:natom) = 1_ark
        hess(:natom,natom+1) = 1_ark
        rho2rhs(:natom,1) = 0_ark
        rho2rhs(natom+1,1) = qtot - sum(charge)
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
            if (allocated(rhomol)) then
                if (size(rhomol) /= npts) deallocate (rhomol, rhopro, rhoatm, rhoagrad)
            end if
            if (.not. allocated(rhomol)) allocate (rhomol(npts),rhopro(npts),rhoatm(natom,npts),rhoagrad(natom,npts))
            !
            !  Evaluate atomic densities at grid points if necessary
            !
            call evaluate_atomic(aden, ax, inp%n_rad_atom, xyzq(1:3,:), iwhr, nuniq, natom, inp%interp_type, inp%interp_ord, xyzw(1:3,:), npts, rhoatm)
            call evaluate_atomic(agrad, ax, inp%n_rad_atom, xyzq(1:3,:), iwhr, nuniq, natom, "pol", 3, xyzw(1:3,:), npts, rhoagrad)
            rhopro = sum(rhoatm, dim=1)
            !
            !  Read in molecular densities
            !
            read_rho: do ipt = 1,npts
                read(mfile) rhomol(ipt)
            end do read_rho
            !
            !  Find the atomic contributions and output densities
            !
            normatm = normatm + sum(xyzw(4,:) * rhopro)
            ddsq = ddsq + sum(xyzw(4,:) * (rhomol - rhopro)**2)
            integrate_pro: do ipt = 1,npts
                write(qfile,'(*(e16.8))') rhoatm(:,ipt)
            end do integrate_pro
            integrate_grad: do iat = 1,natom
                rho2rhs(iat,1) = rho2rhs(iat,1) + 2*sum(xyzw(4,:)*(rhomol - rhopro)*rhoagrad(iat,:))
                do jat = 1,iat
                    hess(iat,jat) = hess(iat,jat) + 2*sum(xyzw(4,:)*rhoagrad(iat,:)*rhoagrad(jat,:))
                end do
            end do integrate_grad
        end do grid_batches
        close(qfile)
        !
        !  Find Hessian and change in charge
        !
        do jat = 1,natom
            do iat = 1,jat-1
                hess(iat,jat) = hess(jat,iat)
            end do
        end do
        !print *,size(hess,dim=1),size(hess,dim=2)
        !print '(5f12.6)',(hess(iat,:), iat=1,natom+1)
        !print *,size(rho2rhs,dim=1),size(rho2rhs,dim=2)
        !print '(5f12.6)',rho2rhs
        !call lapack_dposv(hess, rho2rhs)
        call lapack_dsysv(hess, rho2rhs)
        dcharge = rho2rhs(:natom,1)
        charge = charge + dcharge

        do i = 1,natom
            if (charge(i) > qlist(iwhr(i))) then
                write(ofile,'("iterate_chg: Charge of ",f8.3," for atom ",i3," exceeds nuclear charge")') charge(i), i
                stop "iterate_chg - atomic charge exceeds nuclear charge"
            end if
        end do

        write(ofile,'("Total atomic density: ",f14.8)') normatm
        write(ofile,'("Square of deformation density: ",f14.8)') ddsq
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
    call system("rm -f dens.atm.[0-9]*")
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')

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

!
!  Determine the atomic densities at grid XYZ coordinates
!
subroutine evaluate_atomic(rden, r, nr, xyzc, uind, nu, natm, ityp, iord, xyz, npt, rho)
    character(3) :: ityp
    integer(ik)  :: nr, natm, npt, nu, uind(natm), iord
    real(rk)     :: r(nu,nr), xyzc(3,natm), xyz(3,npt)
    real(ark)    :: rden(natm,nr), rho(natm,npt)
    !
    integer(ik)  :: iat, ipt, ui
    real(rk)     :: ratm

    rho(:,:) = 0_ark
    evaluate_atom: do iat = 1,natm
        ui = uind(iat)
        interp_point: do ipt = 1,npt
            ratm = sqrt((xyz(1,ipt)-xyzc(1,iat))**2 + &
                        (xyz(2,ipt)-xyzc(2,iat))**2 + &
                        (xyz(3,ipt)-xyzc(3,iat))**2)
            ! atomic density = 0 outside the grid
            if (ratm < r(ui,nr)) then
                rho(iat,ipt) = rho(iat,ipt) + interp(ratm, r(ui,:), rden(iat,:), nr, ityp, iord)
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
    integer(ik)  :: ia, il, ui
    real(rk)     :: modc, cthrsh=1e-6

    update_aden: do ia = 1,natm
        if (abs(chg(ia)) > maxc) then
            write(ofile,'("update_atoms: abs(charge) greater than MAXCHG = ",i3)') maxc
            stop "update_atoms - charge out of bounds"
        end if
        ui = uind(ia)
        ! find the modulus (as it should be defined)
        modc = chg(ia) - floor(chg(ia))
        if (abs(modc) < cthrsh) then
            ! integer lower charge
            il = maxc + 1 + nint(chg(ia))
        else
            ! real lower charge
            il = maxc + 1 + floor(chg(ia))
        end if
        if (.not. aload(il,ui)) then
            ! import the lower integer density
            call load_atom(ql(ui), il-maxc-1, narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il,ui,:))
            aload(il,ui) = .true.
        end if
        if (.not. aload(il+1,ui)) then
            ! import the upper integer density
            call load_atom(ql(ui), il-maxc, narad, alib, atyp, ax(ui,:), aw(ui,:), acden(il+1,ui,:))
            aload(il+1,ui) = .true.
        end if
        aden(ia,:) = (1-modc)*acden(il,ui,:) + modc*acden(il+1,ui,:)
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
    logical      :: file_exists

    if (chg >= q) then
        aden = 0_ark
        return
    end if
    select case (atyp)
        case ("slater")
            write(ofile,'("Loading Slater atomic density for Z = ",i3,", CHG = ",i3)') q, chg
            norb = get_norb(q-chg)
            call psisq(ax, aw, narad, q, chg, norb, niter, thrsh, aden)
        case default
            inquire(file=trim(inp%lib_path)//alib, exist=file_exists)
            if (.not. file_exists) then
                write(ofile,'("load_atom: Error unrecognized type",a10)') atyp
                stop "load_atom - unrecognized type"
            end if
            write(ofile,'("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)') q, chg
            open(afile, file=trim(inp%lib_path)//alib)
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
    end select

end subroutine load_atom

end program newton
