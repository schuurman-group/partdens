!
!  The (hopefully temporary) Bader QTAIM charge algorithm
!
!  Calculates charges from a GAMESS data file using the Bader
!  atoms-in-molecules algorithm which requires density gradients. This
!  can potentially be added to the GridCharge code if there is a good way
!  to cache unassigned points (i.e. save x, y, z, w and rho and propagate
!  the assignment until assigned).
!
!  For now, QTAIM deformation densities are calculated by gridchg.x after
!  running bader.x with inp%atom_type = pro. This is done by back-propagating
!  assignments and saving the results in a binary file.
!
!  Rodriguez, J.I.; Kostner, A.M.; Ayers, P.W.; Santos-Valle, A.; Vela, A.;
!  Merino, G. J. Comput. Chem. 2008, 30, 1082-1092.
!
program bader
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use utilities
    use atoms
    use fileio
    !
    character(2),allocatable :: atypes(:)
    logical                  :: ptrust
    integer(ik)              :: ofile, mfile, nfile, pfile, gfile
    integer(ik)              :: i, ib, ipt, npts, iat
    integer(ik)              :: nat_count, natom, ndum, nbatch
    integer(ik),allocatable  :: asgn(:)
    real(rk)                 :: norm, tmp_nat_occ(1000), unchg
    real(rk),pointer         :: xyzw(:,:), xyz(:,:)
    real(rk),allocatable     :: xyzq(:,:), dist(:,:), ed(:,:), er(:,:), cosa(:)
    real(rk),allocatable     :: rtrust(:), charge(:)
    real(ark),allocatable    :: rhomol(:), drho(:,:), nat_occ(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid
    type(input_data)         :: inp

    call accuracyInitialize

    ofile=10
    mfile=11
    nfile=12
    pfile=13
    gfile=14

    open(ofile, file="bader.out")
    call read_input(input, inp)
    call init_gridchg_output(ofile, inp)

    if (inp%weight_type /= "qtaim") then
        write(out,'("Weight type ",a," not supported by bader.x")') trim(inp%weight_type)
        stop "bad weight type"
    end if

    write(ofile,'("Loading GAMESS data file")')
    select case (inp%vec_type)
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

    allocate (atypes(mol%natoms), xyzq(4,mol%natoms))
    allocate (rtrust(mol%natoms), charge(mol%natoms))

    call gamess_report_nuclei(natom, xyzq, mol)
    do i = 1,mol%natoms
        atypes(i) = trim(mol%atoms(i)%name)
    enddo
    write(ofile,'("Molecular geometry (in bohr):")')
    ndum = 0
    do i = 1,mol%natoms
        write(ofile,'("    ",a2,3f14.8)') atypes(i), xyzq(1:3,i)
        ! count dummy atoms (assumes they are at first indices)
        if (xyzq(4,i) < 0.5) ndum = ndum + 1
    end do
    write(ofile,'("")')
    !
    !  Set up the grid
    !
    write(ofile,'("Setting up the molecular grid")')
    if (ndum > 0) write(ofile,'("Removing ",i3," dummy atom(s) from grid")') ndum
    write(ofile,'("")')
    call GridInitialize(den_grid, inp%n_rad, inp%n_ang, xyzq(1:3,ndum+1:), atypes(ndum+1:))
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    open(mfile, file='moldens', form='unformatted', action='write')
    open(pfile, file='dens.mol', action='write')
    open(gfile, file='grad.mol', action='write')
    !
    !  Molecular density numerical integration loop
    !
    write(ofile,'("Calculating molecular density")')
    nullify(xyzw)
    norm = 0_rk
    rtrust = 0_rk
    iat = 1 + ndum
    ptrust = .true.
    mol_grid_batches: do ib = 1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
        xyz => xyzw(1:3,:)
        !
        !  If the batch size changed, reallocate rho and other values
        !
        npts = size(xyzw, dim=2)
        if (allocated(rhomol)) then
            if (size(rhomol) /= npts) deallocate (rhomol, drho, dist, ed, er, cosa)
        end if
        if (.not. allocated(rhomol)) allocate (rhomol(npts), drho(3,npts), dist(3,npts), ed(3,npts), er(3,npts), cosa(npts))
        !
        !  Evaluate density and gradients at grid points
        !
        call evaluate_density_gradients(inp%vec_type, nat_count, npts,  mol, nat_occ, xyz, rhomol, drho)
        !
        !  Propagate the trust sphere, if appropriate
        !
        if (ptrust) then
            ! get the distance vector to the atom centre
            do i = 1,3
                dist(i,:) = xyzq(i,iat)
            end do
            dist = dist - xyz
            ! get unit vectors of distances and gradients
            do ipt = 1,npts
                ed(:,ipt) = drho(:,ipt) / sqrt(sum(drho(:,ipt)**2, dim=1))
                er(:,ipt) = dist(:,ipt) / sqrt(sum(dist(:,ipt)**2, dim=1))
            end do
            ! find cos(ed.er)
            cosa = sum(ed*er, dim=1)
            print *,cosa
            ! check for gradients that point away from the centre
            check_cosa: do ipt = 1,npts
                if (cosa(ipt) < 0.5_rk*sqrt(2.0_rk)) then
                    ptrust = .false.
                    exit check_cosa
                end if
            end do check_cosa
            ! if all pointing to the centre, set rtrust
            if (ptrust) then
                rtrust(iat) = sqrt(sum(dist(:,1)**2, dim=1))
            end if
        end if
        if (den_grid%next_part + ndum /= iat) then
            ! start propagating if the atom centre will change
            ptrust = .true.
            iat = den_grid%next_part + ndum
        end if
        !
        !  Integrate and save
        !
        mol_integrate: do ipt = 1,npts
            norm = norm + xyzw(4,ipt) * rhomol(ipt)
            write(mfile) rhomol(ipt), drho(1:3,ipt)
            write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rhomol(ipt)
            write(gfile,1000) xyzw(4,ipt), drho(1:3,ipt), rhomol(ipt)
        end do mol_integrate
        stop 'Not an error. End of first shell, check dens.mol and grad.mol.'
    end do mol_grid_batches
    close(mfile)
    close(pfile)
    close(gfile)
    write(ofile,'("Total molecular density: ",f14.8)') norm
    write(ofile,'("Trust radii for each atom:")')
    write(ofile,'("    ",5f14.8)') rtrust
    do i = 1,mol%natoms + ndum
        if (rtrust(i) < 1e-8) then
            write(ofile,'("Trust radius of zero for atom ",i4)') i
            stop "bader - zero trust radius"
        end if
    end do
    write(ofile,'("")')
    deallocate (rhomol, drho, dist, ed, er, cosa)
    !
    !  Molecular density numerical integration loop
    !
    write(ofile,'("Starting QTAIM density assignment")')
    write(ofile,'("")')
    open(mfile, file='moldens', form='unformatted', action='read')
    call GridPointsBatch(den_grid, 'Reset')
    charge = 0_rk
    unchg = 0_rk
    grid_batches: do ib = 1,nbatch
        !
        !  Get grid points
        !
        call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
        !
        !  If the batch size changed, reallocate rho and other values
        !
        npts = size(xyzw, dim=2)
        if (allocated(rhomol)) then
            if (size(rhomol) /= npts) deallocate (rhomol, drho, asgn)
        end if
        if (.not. allocated(rhomol)) allocate (rhomol(npts), drho(3,npts), asgn(npts))
        !
        !  Read in molecular densities and gradients
        !
        read_rho: do ipt = 1,npts
            read(mfile) rhomol(ipt), drho(1:3,ipt)
        end do read_rho
        !
        !  Propagate and assign densities
        !
        iat = den_grid%next_part + ndum
        print *,ib
        asgn = assign_atom(xyzq(1:3,:), 0.9*rtrust, natom, xyzw, npts, ndum, rhomol, drho)
        !
        !  Integrate assigned charges
        !
        chg_integrate: do ipt = 1,npts
            iat = asgn(ipt)
            if (iat == 0) then
                unchg = unchg + xyzw(4,ipt) * rhomol(ipt)
            else
                charge(iat) = charge(iat) + xyzw(4,ipt) * rhomol(ipt)
            end if
        end do chg_integrate
    end do grid_batches

    write(ofile,'("Unassigned density: ",f14.8)') unchg
    write(ofile,'("Atomic populations:")')
    write(ofile,'("    ",5f14.8)') charge
    write(ofile,'("")')
    write(ofile,'("Final charges:")')
    write(ofile,'("    ",5f14.8)') xyzq(4,:) - charge
    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    close(mfile)
    call system("rm moldens")
    ! write assignments?
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')
    close(ofile)

    1000 format(e24.8,f16.8,f16.8,f16.8,e20.10)

contains

!
!  Determine the QTAIM atomic assignments for grid points
!
function assign_atom(xyzatm, rad, natm, xyzw, npt, ndum, rho, drho)
    integer(ik) :: natm, npt, ndum
    real(rk)    :: xyzatm(3,natm), rad(natm), xyzw(4,npt)
    real(ark)   :: rho(npt), drho(3,npt)
    integer(ik) :: assign_atom(npt)
    !
    logical     :: adone(npt)
    integer(ik) :: ipt, iat, iter, max_iter=1000
    real(rk)    :: dist(natm), delta, thresh=1e-8
    real(ark)   :: dnorm, rho_copy(npt)

    adone(:) = .false.
    assign_atom(:) = 0
    rho_copy = rho
    propagate_density: do iter = 1,max_iter
        assign_points: do ipt = 1,npt
            if (.not.adone(ipt)) then
                ! assign points within trust spheres
                do iat = 1+ndum,natm
                    dist(iat) = sqrt(sum((xyzw(1:3,ipt) - xyzatm(:,iat))**2))
                    if (dist(iat) <= rad(iat)) then
                        assign_atom(ipt) = iat
                        adone(ipt) = .true.
                        cycle assign_points
                    end if
                end do
                ! screen small density values
                if (rho_copy(ipt)*xyzw(4,ipt) < thresh) then
                    adone(ipt) = .true.
                    cycle assign_points
                end if
                ! find potential non-nuclear attractors (just error for now)
                dnorm = sqrt(sum(drho(:,ipt)**2))
                if ((rho_copy(ipt) > 1e-4).and.(dnorm < thresh)) then
                    write(ofile,'("assign_atom: non-nuclear attractors not supported")')
                    adone(ipt) = .true.
                    !stop "assign_atom - non-nuclear attractors not supported"
                    cycle assign_points
                end if
                ! propagate point by one step
                delta = min(0.99*minval(dist), 0.01)
                !call rkfour(xyzw(1:3,ipt), drho(3,ipt), rho_copy(ipt), delta)
                call euler(xyzw(1:3,ipt), drho(3,ipt), rho_copy(ipt), delta)
            end if
        end do assign_points
        if (all(adone)) then
            return
        else if (iter == max_iter) then
            write(ofile,'("QTAIM trajectories not converged in MAXITER = ",i4)') max_iter
            stop "assign_atom - QTAIM trajectories not converged"
        end if
    end do propagate_density

end function assign_atom

!
!  Propagate grid positions by the density gradients using RK4
!
subroutine rkfour(xyz, dxyz, rho, delta)
    real(rk),intent(in)     :: delta
    real(rk),intent(inout)  :: xyz(3)
    real(ark),intent(inout) :: dxyz(3), rho
    !
    integer(ik)             :: i, it, max_it=100
    real(rk)                :: kxyz(3,1), fxyz(3), k(3,4), wgt(4), coef(3), sc
    real(ark)               :: rho_step(1), drho_step(3,1)

    wgt = (/ 0.5_rk, 1.0_rk, 1.0_rk, 0.5_rk /) / 3.0_rk
    coef = (/ 0.5_rk, 0.5_rk, 1.0_rk /)

    ! maximum step size
    sc = delta / sqrt(sum(dxyz**2))

    try_step: do it = 1,max_it
        fxyz = xyz
        kxyz(:,1) = xyz
        drho_step(:,1) = dxyz
        rk_step: do i = 1,4
            k(:,i) = sc*drho_step(:,1)
            if (i < 4) then
                kxyz(:,1) = kxyz(:,1) + coef(i)*k(:,i)
                call evaluate_density_gradients(inp%vec_type, nat_count, 1,  mol, nat_occ, kxyz, rho_step, drho_step)
            end if
            fxyz = fxyz + wgt(i)*k(:,i)
        end do rk_step
        kxyz(:,1) = fxyz
        !print *,wgt(1)*k(:,1)+wgt(2)*k(:,2)+wgt(3)*k(:,3)+wgt(4)*k(:,4)
        call evaluate_density_gradients(inp%vec_type, nat_count, 1,  mol, nat_occ, kxyz, rho_step, drho_step)
        if (rho_step(1) >= rho) exit try_step ! comment out to test propagation
        ! if density decreases, scale step and try again
        if (it == max_it) then
            write(ofile,'("Density decreasing after maximum number of attemps")')
            !stop "rkfour - density decreasing"
        end if
        sc = 0.25*sc
    end do try_step
    xyz = fxyz
    dxyz = drho_step(:,1)
    rho = rho_step(1)

end subroutine rkfour

!
!  Propagate grid positions by the density gradients using Euler
!
subroutine euler(xyz, dxyz, rho, delta)
    real(rk),intent(in)     :: delta
    real(rk),intent(inout)  :: xyz(3)
    real(ark),intent(inout) :: dxyz(3), rho
    !
    integer(ik)             :: it, max_it=100
    real(rk)                :: fxyz(3,1), sc
    real(ark)               :: rho_step(1), drho_step(3,1)

    ! maximum step size
    sc = delta / sqrt(sum(dxyz**2))

    try_step: do it = 1,max_it
        fxyz(:,1) = xyz + sc*dxyz
        call evaluate_density_gradients(inp%vec_type, nat_count, 1, mol, nat_occ, fxyz, rho_step, drho_step)
        if (rho_step(1) >= rho) exit try_step
        ! if density decreases, scale step and try again
        if (it == max_it) then
            write(ofile,'("Density decreasing after maximum number of attemps")')
            stop "euler - density decreasing"
        end if
        sc = 0.25*sc
    end do try_step
    xyz = fxyz(:,1)
    dxyz = drho_step(:,1)
    rho = rho_step(1)

end subroutine euler

end program bader
