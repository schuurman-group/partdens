!
!  The noncovalent interaction density calculation
!
!  Calculates densities and reduced density gradients from a GAMESS RDM
!  output file. These can then be plotted
!
program ncicalc
    use accuracy
    use import_gamess
    use gamess_internal
    use utilities
    use atoms
    !
    character(100),parameter :: infile="nci.inp", outfile="nci.out"
    character(100)           :: nat_file, vectyp
    character(2),allocatable :: atypes(:)
    integer(ik)              :: ofile, mfile
    integer(ik)              :: i, ib, ipt, npts, nptmax
    integer(ik)              :: nat_count, natom, nbatch
    integer(ik)              :: nrad, nang
    real(rk)                 :: cmax, cmin, ccons, dxyz
    real(rk)                 :: norm, sconst, covrad, tmp_nat_occ(1000)
    real(rk)                 :: gridbnd(3,2), xyzinit(3)
    real(rk),allocatable     :: xyzq(:,:), xyz(:,:), rbnd(:,:)
    real(ark),allocatable    :: rho(:), sval(:), drho(:,:), nat_occ(:)
    type(gam_structure)      :: mol

    call accuracyInitialize

    ofile=10
    mfile=11
    sconst = 1_ark / (2_ark*(3_ark*pi)**(1_ark/3_ark))

    open(ofile, file=outfile)
    call read_input(infile, nat_file, vectyp, nrad, nang, cmax, cmin, ccons, dxyz)
    call init_output(nat_file, vectyp, nrad, nang, cmax, cmin, ccons, dxyz)

    write(ofile,'("Loading GAMESS data file")')
    select case (vectyp)
        case ("tr1rdm")
            call gamess_load_rdmsv(trim(nat_file), tmp_nat_occ, nat_count)
        case ("natorb")
            call gamess_load_natocc(trim(nat_file), tmp_nat_occ, nat_count)
        case default
            write(out,'("Unrecognized VEC type ",a8)') vectyp
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
    allocate (nat_occ(nat_count))
    nat_occ = tmp_nat_occ(:nat_count)
    !
    call gamess_load_orbitals(file=trim(nat_file), structure=mol)

    allocate (atypes(mol%natoms), xyzq(4,mol%natoms))

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
    nptmax = 1000
    allocate (rbnd(mol%natoms,2), xyz(3,nptmax))
    do i=1,natom
        covrad = AtomCovalentR(atypes(i)) / abohr
        rbnd(i,2) = cmax * covrad + ccons
        rbnd(i,1) = cmin * covrad
    end do
    nbatch = 1
    do i=1,3
        gridbnd(i,1) = minval(xyzq(i,:) - rbnd(:,2))
        gridbnd(i,2) = maxval(xyzq(i,:) + rbnd(:,2))
        nbatch = nbatch * floor((gridbnd(i,2) - gridbnd(i,1)) / dxyz)
    end do
    nbatch = nbatch / nptmax
    xyzinit = gridbnd(:,1)
    !
    !  Molecular density numerical integration loop
    !
    write(ofile,'("Calculating molecular density")')
    open(mfile, file='nci.dat', action='write')
    norm = 0_rk
    mol_grid_batches: do ib = 1,nbatch
        !
        !  Get grid points
        !
        call cartgrid_batch(xyzinit, dxyz, gridbnd, xyzq(1:3,:), rbnd, natom, nptmax, xyz, npts)
        !
        !  If the batch size changed, reallocate rho and drho
        !
        if (allocated(rho)) then
            if (size(rho) /= npts) deallocate (rho)
            if (size(sval) /= npts) deallocate (sval)
        end if
        if (allocated(drho)) then
            if (size(drho, dim=2) /= npts) deallocate (drho)
        end if
        if (.not. allocated(rho)) allocate (rho(npts))
        if (.not. allocated(sval)) allocate (sval(npts))
        if (.not. allocated(drho)) allocate (drho(3,npts))
        !
        !  Evaluate density and gradient at grid points
        !
        call evaluate_density_gradients(vectyp, nat_count, npts, mol, nat_occ, xyz(:,:npts), rho, drho)
        sval = sconst * sqrt(sum(drho**2, dim=1)) / rho**(4_ark/3_ark)
        !
        !  Integrate and save
        !
        mol_integrate: do ipt = 1,npts
            norm = norm + dxyz**3 * rho(ipt)
            write(mfile,1000) xyz(:,ipt), rho(ipt), sval(ipt)
        end do mol_integrate
        if (npts < nptmax) exit mol_grid_batches
        xyzinit = xyz(:,nptmax)
    end do mol_grid_batches
    close(mfile)
    write(ofile,'("Total molecular density: ",f14.8)') norm
    write(ofile,'("")')

    !
    !  Clean up
    !
    close(mfile)
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')

    1000 format(f16.8,f16.8,f16.8,e20.10,e20.10)

contains

!
!  Read a crude input file
!
subroutine read_input(infile, rfile, vtyp, nr, na, fmax, fmin, cnst, dg)
    character(100) :: infile, rfile, vtyp
    integer(ik)    :: nr, na
    real(rk)       :: fmax, fmin, cnst, dg
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
            case ("nat_filename")
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
            case ("covrad_max_scale")
                read(var,*,iostat=ios) fmax
                nvar = nvar + 1
            case ("covrad_min_scale")
                read(var,*,iostat=ios) fmin
                nvar = nvar + 1
            case ("covrad_const")
                read(var,*,iostat=ios) cnst
                nvar = nvar + 1
            case ("dgrid")
                read(var,*,iostat=ios) dg
                nvar = nvar + 1
        end select
        if (ios /= 0) then
            write(ofile,'("read_input: Incorrect type for variable ",a)') varname
            stop "read_input - incorrect variable type"
        end if
    end do parse_input
    if (nvar < 8) then
        write(ofile,'("read_input: Missing variable")')
        stop "read_input - missing variable"
    end if

end subroutine read_input

!
!  Initialize the output file
!
subroutine init_output(rfile, vtyp, nr, na, fmax, fmin, cnst, dg)
    character(100) :: rfile, vtyp
    integer(ik)    :: nr, na
    real(rk)       :: fmax, fmin, cnst, dg

    write(ofile,'("+---------------------------------------------------+")')
    write(ofile,'("|                                                   |")')
    write(ofile,'("|                      NCI                          |")')
    write(ofile,'("|                                                   |")')
    write(ofile,'("|  Noncovalent interaction densities and gradients  |")')
    write(ofile,'("|          RJ MacDonell, MS Schuurman 2018          |")')
    write(ofile,'("+---------------------------------------------------+")')
    write(ofile,'("")')
    write(ofile,'("")')
    write(ofile,'("Input summary:")')
    write(ofile,'("    ------- Molecular density --------")')
    write(ofile,'("    nat_filename   =   ",a15)') trim(rfile)
    write(ofile,'("    vec_type       =   ",a15)') trim(vtyp)
    write(ofile,'("    n_r_grid       =   ",i15)') nr
    write(ofile,'("    n_ang_grid     =   ",i15)') na
    write(ofile,'("")')
    write(ofile,'("    -------------- Grid --------------")')
    write(ofile,'("    covrad_max     =   ",f15.4)') fmax
    write(ofile,'("    covrad_min     =   ",f15.4)') fmin
    write(ofile,'("    covrad_const   =   ",f15.4)') cnst
    write(ofile,'("    dgrid          =   ",f15.4)') dg
end subroutine init_output

!
!  Generate a cartesian grid batch within a certain radius of atoms
!
subroutine cartgrid_batch(xyzi, dx, xyzbnd, xyzatm, rbnd, natm, nmax, xyzgrid, npt)
    integer(ik),intent(in)  :: natm, nmax
    real(rk),intent(in)     :: xyzi(3), xyzbnd(3,2), dx
    real(rk),intent(in)     :: xyzatm(3,natm), rbnd(natm,2)
    integer(ik),intent(out) :: npt
    real(rk),intent(out)    :: xyzgrid(3,nmax)
    !
    integer(ik)             :: i
    real(rk)                :: cxyz(3), dist(natm)

    cxyz = xyzi
    npt = 0
    generate_grid: do while (npt < nmax)
        cxyz(1) = cxyz(1) + dx
        if (cxyz(1) > xyzbnd(1,2)) then
            ! maximum x, go to xmin and next y
            cxyz(1) = xyzbnd(1,1)
            cxyz(2) = cxyz(2) + dx
            if (cxyz(2) > xyzbnd(2,2)) then
                ! maximum y, go to ymin and next z
                cxyz(2) = xyzbnd(2,1)
                cxyz(3) = cxyz(3) + dx
                if (cxyz(3) > xyzbnd(3,2)) then
                    ! maximum z, end of grid
                    return
                end if
            end if
        end if
        ! get the new distance
        do i=1,natm
            dist(i) = sqrt(sum((cxyz - xyzatm(:,i))**2))
        end do
        ! if the point is within the bounds, add it to the grid
        if ((maxval(rbnd(:,2) - dist) > 0) .and. (minval(dist - rbnd(:,1)) > 0)) then
            npt = npt + 1
            xyzgrid(:,npt) = cxyz
        end if
    end do generate_grid

end subroutine cartgrid_batch

end program ncicalc
