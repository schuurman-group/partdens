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
    use utilities
    use atoms
    use lapack
    use fileio
    !
    character(20)            :: pname
    character(2),allocatable :: atypes(:)
    logical,allocatable      :: dload(:,:)
    integer(ik)              :: ofile, mfile, pfile, qfile
    integer(ik)              :: i, ib, ipt, is1, is2, iter, npts, iat, jat
    integer(ik)              :: nat_count, natom, nbatch, nuniq
    integer(ik),allocatable  :: iwhr(:), qlist(:), npshell(:)
    real(rk)                 :: norm, normatm, ddsq, qtot, tmp_nat_occ(1000)
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:)
    real(rk),allocatable     :: charge(:), dcharge(:)
    real(ark),allocatable    :: aden(:,:,:,:), pdens(:,:), pgrad(:,:)
    real(ark),allocatable    :: rhomol(:), rhopro(:), rhoatm(:,:), rhoagrad(:,:)
    real(ark),allocatable    :: hess(:,:), rho2rhs(:,:), nat_occ(:)
    type(gam_structure)      :: mol, promol
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
    allocate (charge(mol%natoms), dcharge(mol%natoms), npshell(mol%natoms))

    call gamess_report_nuclei(natom, xyzq, mol)
    write(ofile,'("Molecular geometry (in bohr):")')
    do iat = 1,natom
        atypes(iat) = trim(mol%atoms(iat)%name)
        write(ofile,'("    ",a2,3f14.8)') atypes(iat), xyzq(1:3,iat)
    enddo
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
    allocate (qlist(nuniq), pdens(sum(npshell),sum(npshell)))
    allocate (dload(2*inp%max_charge+1,nuniq))
    allocate (aden(maxval(npshell),maxval(npshell),2*inp%max_charge+1,nuniq))
    allocate (pgrad(sum(npshell),sum(npshell)))
    allocate (hess(natom+1,natom+1), rho2rhs(natom+1,1))
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
        call update_densmat(charge, inp%max_charge, qlist, natom, iwhr, nuniq, aden, dload, npshell, inp%lib_dmat, inp%atom_bas, pdens)
        is1 = 1
        update_pgrad: do iat = 1,natom
            if (abs(charge(iat) - floor(charge(iat))) < 1e-6) then
                i = inp%max_charge + 1 + nint(charge(iat))
            else
                i = inp%max_charge + 1 + floor(charge(iat))
            end if
            is2 = is1 + npshell(iat) - 1
            pgrad(is1:is2,is1:is1) = aden(:npshell(iat),:npshell(iat),i+1,iwhr(iat)) - aden(:npshell(iat),:npshell(iat),i,iwhr(iat))
            is1 = is1 + npshell(iat)
        end do update_pgrad
        !
        normatm = 0_rk
        ddsq = 0_rk
        hess = 0_ark
        hess(natom+1,:natom) = 1_ark
        hess(:natom,natom+1) = 1_ark
        rho2rhs(:natom,1) = 0_ark
        rho2rhs(natom+1,1) = qtot - sum(charge)
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
            if (allocated(rhomol)) then
                if (size(rhomol) /= npts) deallocate (rhomol, rhopro, rhoatm, rhoagrad)
            end if
            if (.not. allocated(rhomol)) allocate (rhomol(npts),rhopro(npts),rhoatm(natom,npts),rhoagrad(natom,npts))
            !
            !  Evaluate atomic densities at grid points
            !
            call evaluate_atomic(xyzw(1:3,:), npts, promol, pdens, npshell, natom, rhoatm)
            call evaluate_atomic(xyzw(1:3,:), npts, promol, pgrad, npshell, natom, rhoagrad)
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

end program newton
