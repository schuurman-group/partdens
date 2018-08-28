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
    character(2),allocatable :: atypes(:), gatyp(:)
    logical                  :: first_shell, last, ls
    logical,allocatable      :: adone(:,:)
    integer(ik)              :: ofile, mfile
    integer(ik)              :: i, ib, ipt, npts, nax
    integer(ik)              :: rdm_count, natom, ndum, nbatch
    integer(ik)              :: nrad, nang, narad, maxiter, iordr, maxchg
    integer(ik),allocatable  :: nnn(:), nmask(:,:), asgn(:,:)
    real(rk)                 :: norm, tmp_rdm_sv(1000), thrsh, dist
    real(rk)                 :: avg_now(3), avg_nxt(3)
    real(rk),pointer         :: xyzw_now(:,:), xyzw_nxt(:,:)
    real(rk),allocatable     :: xyzq(:,:), xyzg(:,:), pt_dist(:), xyz(:,:,:)
    real(rk),allocatable     :: wgt_now(:), charge(:)
    real(ark),allocatable    :: rhomol(:), wrho(:,:), drho(:,:,:), rdm_sv(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid

    call accuracyInitialize

    ofile=10
    mfile=11

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

    allocate (atypes(mol%natoms), xyzq(4,mol%natoms), xyzg(3,mol%natoms))
    allocate (gatyp(mol%natoms), charge(mol%natoms))

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
    ndum = 0
    do i = 1,natom
        if (atypes(i) == "X") then
            write(ofile,'("Skipping dummy atom at index ",i3)') i
        else
            ndum = ndum + 1
            gatyp(ndum) = atypes(i)
            xyzg(:,ndum) = xyzq(1:3,i)
        end if
    end do
    write(ofile,'("")')
    call GridInitialize(den_grid, nrad, nang, xyzg(:,1:ndum), gatyp(1:ndum))
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    !!open(mfile, file='moldens', form='unformatted', action='write')
    !
    !  Set initial values
    !
    nullify(xyzw_now, xyzw_nxt)
    call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_now)
    npts = size(xyzw_now, dim=2)
    allocate (rhomol(npts), wrho(2,npts), drho(2,3,npts))
    allocate (nnn(npts), nmask(7,npts), pt_dist(npts), wgt_now(npts))
    allocate (xyz(2,3,npts), asgn(2,npts), adone(2,npts))
    !
    !  Find nearest neighbours
    !
    write(ofile,'("Finding nearest neighbours in spherical grid")')
    write(ofile,'("")')
    find_nearest: do ipt = 1, npts
        nax = 0
        do i = 1, 3
            if (abs(xyzw_now(i,ipt)) < 1e-10) nax = nax + 1
        end do
        if (nax == 2) then
            ! on-axis points only have 4 nearest neighbours
            nnn(ipt) = 5
        else
            ! all other points have 6 nearest neighbours
            nnn(ipt) = 7
        end if
        pt_dist = sqrt((xyzw_now(1,:) - xyzw_now(1,ipt))**2 + &
                       (xyzw_now(2,:) - xyzw_now(2,ipt))**2 + &
                       (xyzw_now(3,:) - xyzw_now(3,ipt))**2)
        ! get indices of current point and its nearest neighbours
        nmask(:,ipt) = argslow(nnn(ipt), pt_dist, npts)
    end do find_nearest
    !
    !  Molecular density numerical integration loop
    !
    write(ofile,'("Calculating molecular density")')
    ! get density and gradients for the first shell
    call evaluate_density_gradients(vectyp, rdm_count, npts, mol, rdm_sv, xyzw_now(1:3,:), rhomol, drho(2,:,:))
    xyz(2,:,:) = xyzw_now(1:3,:)
    wrho(2,:) = rhomol * xyzw_now(4,:)
    norm = sum(wrho(2,:))
    charge(:) = xyzq(4,:)
    asgn(:,:) = 0
    adone(:,:) = .false.
    first_shell = .true.
    last = .false.
    mol_grid_batches: do ib = 1, nbatch-1
        if (first_shell) then
            ! get next grid points
            call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_nxt)
            first_shell = .false.
            ls = .false.
        else
            if (ib == nbatch-1) then
                ! the next iteration is the last shell
                last = .true.
                ls = .true.
            else
                ! get next grid points
                call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_nxt)

                ! check to see if the next shell is a different atom
                avg_now = sum(xyz(2,:,:), dim=2)
                avg_nxt = sum(xyzw_nxt(1:3,:), dim=2)
                dist = sqrt(sum((avg_nxt - avg_now)**2, dim=1))
                if (dist > 0.2) then
                    ls = .true.
                else
                    ls = .false.
                end if
            end if
            !
            !  Evaluate density and density gradients at grid points
            !
            call evaluate_density_gradients(vectyp, rdm_count, npts, mol, rdm_sv, xyz(2,:,:), rhomol, drho(2,:,:))
            !
            !  Assign the densities for a set of two shells
            !
            wrho(2,:) = rhomol * wgt_now(:)
            norm = norm + sum(wrho(2,:))
            call assign_shell(xyz, wrho, drho, asgn, adone, npts, xyzq, charge, natom, nmask, nnn, ls)
        end if
        !
        !  Prepare the next shell
        !
        if (.not. last) then
            xyz(1,:,:) = xyz(2,:,:)
            xyz(2,:,:) = xyzw_nxt(1:3,:)
            wgt_now(:) = xyzw_nxt(4,:)
            wrho(1,:) = wrho(2,:)
            drho(1,:,:) = drho(2,:,:)
            asgn(1,:) = asgn(2,:)
            asgn(2,:) = 0
            adone(1,:) = adone(2,:)
            adone(2,:) = .false.
            if (ls) then
                first_shell = .true.
            end if
        end if
    end do mol_grid_batches
    !!close(mfile)
    write(ofile,'("Total molecular density: ",f14.8)') norm
    write(ofile,'("")')
    write(ofile,'("Total charge: ",f14.8)') sum(charge)
    write(ofile,'("Final charges:")')
    write(ofile,'("    ",5f14.8)') charge

    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    !!close(mfile)
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
    integer(ik)               :: ipt, ird, imo, ic, nmo, nbas
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
                do ic = 1, 3
                    drho(ic,:) = drho(ic,:) + rdm_sv(ird) * (moval(imo,:) * dmoval(ic,imo+1,:) + moval(imo+1,:) * dmoval(ic,imo,:))
                end do
            end do evaluate_rdm
        case ("natorb")
            evaluate_nat: do ird=1,rdm_count
                rho = rho + rdm_sv(ird) * moval(ird,:)**2
                do ic = 1, 3
                    drho(ic,:) = drho(ic,:) + 2 * rdm_sv(ird) * moval(ird,:) * dmoval(ic,ird,:)
                end do
            end do evaluate_nat
        case default
            write(out,'("evaluate_density: Unrecognized VEC type ",a8)') vtyp
            stop "evaluate_density - bad VEC type"
    end select
    !
    deallocate (basval,moval,dmoval)

end subroutine evaluate_density_gradients

!
!  Determine the Bader atomic assignments
!
subroutine assign_shell(xyz, wrho, drho, asgn, adone, npts, xyzq, atmchg, natm, nrst, nnn, ls)
    integer(ik)           :: npts, natm, nnn(npts), asgn(2,npts), nrst(7,npts)
    real(rk)              :: xyz(2,3,npts), xyzq(4,natm), atmchg(natm)
    real(ark)             :: wrho(2,npts), drho(2,3,npts)
    logical               :: adone(2,npts), ls
    !
    integer(ik)           :: i, j, k, n, iout(2,npts), o_ji, a_ji, o_next, a_next, iatm
    integer(ik),parameter :: maxitr=12
    real(rk)              :: dist(natm), din(3)
    logical               :: shell_done, askip(2,npts)

    iout(:,:) = 1
    askip(:,:) = .false.

    make_assignments: do i = 1, npts
        if (.not. adone(1,i)) then
            if (sqrt(sum(drho(1,:,i)**2)) < 1e-36) then
                ! gradient is zero, skip it
                if (wrho(1,i) > 1e-36) write(out,'("WARNING: skipping zero gradient with density ",e12.4)') wrho(1,i)
                adone(1,i) = .true.
            else
                ! make assignments for inner shell
                call assign_point(drho(1,:,i), xyz, npts, nrst(:,i), nnn(i), .false., iout(1,i), asgn(1,i))
            end if
        end if

        if (.not. adone(2,i)) then
            din = xyz(1,:,i) - xyz(2,:,i)
            if (sqrt(sum(drho(2,:,i)**2)) < 1e-36) then
                ! gradient is zero, skip it
                if (wrho(2,i) > 1e-36) write(out,'("WARNING: skipping zero gradient with density ",e12.4)') wrho(2,i)
                adone(2,i) = .true.
            else if (unit_dot(din, drho(2,:,i), 3) < 0 .and. .not. ls) then
                ! gradient pointing out, come back after next shell is loaded
                askip(2,i) = .true.
            else
                ! make assignments for outer shell
                call assign_point(drho(2,:,i), xyz, npts, nrst(:,i), nnn(i), .true., iout(2,i), asgn(2,i))
            end if
        end if
    end do make_assignments

    propagate_assignments: do k = 1, maxitr
        shell_done = .true.
        loop_shells: do j = 1, 2
            propagate_shell: do i = 1, npts
                o_ji = iout(j,i)
                a_ji = asgn(j,i)
                if (adone(j,i) .or. askip(j,i)) cycle propagate_shell
                o_next = iout(o_ji,a_ji)
                a_next = asgn(o_ji,a_ji)
                if (o_next == j .and. a_next == i) then
                    ! repetition, add both to nearest atom density
                    do n = 1, natm
                        dist(n) = sqrt(sum((xyz(j,:,i) - xyzq(1:3,n))**2))
                    end do
                    iatm = minloc(dist, 1)
                    if (dist(iatm) > 1e-2) write(out,'("WARNING: maximum found ",e12.4," a0 from atom ",i3)') dist(iatm), iatm
                    asgn(j,i) = -iatm
                    asgn(o_ji,a_ji) = -iatm
                    atmchg(iatm) = atmchg(iatm) - wrho(j,i) - wrho(o_ji,a_ji)
                    adone(j,i) = .true.
                    adone(o_ji,a_ji) = .true.
                else if (a_next == 0) then
                    if (j == 1) then
                        ! pointed-to index out but not assigned, add to density (caching would happen here)
                        wrho(o_ji,a_ji) = wrho(o_ji,a_ji) + wrho(j,i)
                        adone(j,i) = .true.
                    else
                        ! outer shell to unassigned outer shell, come back after loading next shell
                        askip(j,i) = .true.
                    end if
                else if (a_next < 0) then
                    ! pointed-to index is an atom, add to atomic density
                    asgn(j,i) = a_next
                    iatm = -a_next
                    atmchg(iatm) = atmchg(iatm) - wrho(j,i)
                    adone(j,i) = .true.
                else
                    ! pointed-to index is assigned, keep propagating
                    asgn(j,i) = a_next
                    iout(j,i) = o_next
                end if

                !shell_done = shell_done .and. (adone(j,i) .or. askip(j,i))
            end do propagate_shell
        end do loop_shells

        shell_done = all(adone .or. askip)
        if (shell_done) then
            ! all points going out or assigned to atoms
            exit propagate_assignments
        else if (k == maxitr) then
            write(out,'("Shell assignments not finished in maximum number of iterations")')
            stop "propagate_assignments: Maximum iterations exceeded"
        end if
    end do propagate_assignments
end subroutine assign_shell

!
!  Assign a single point to anoter points based on the steepest ascent
!
subroutine assign_point(idrho, xyz, npts, inrst, nni, outer, iiout, iasgn)
    integer(ik) :: npts, nni, inrst(7), iiout, iasgn
    real(rk)    :: idrho(3), xyz(2,3,npts)
    logical     :: outer
    !
    integer(ik) :: itr
    real(rk)    :: vdot, vdotmx, xyzi(3), xyz0(3,nni), xyz1(3,nni), dv(3)

    ! find position and neighbours on current shell and other shell
    if (outer) then
        xyzi = xyz(2,:,inrst(1))
        do itr = 1, nni
            xyz0(:,itr) = xyz(2,:,inrst(itr))
            xyz1(:,itr) = xyz(1,:,inrst(itr))
        end do
    else
        xyzi = xyz(1,:,inrst(1))
        do itr = 1, nni
            xyz0(:,itr) = xyz(1,:,inrst(itr))
            xyz1(:,itr) = xyz(2,:,inrst(itr))
        end do
    end if

    iiout = 0
    vdotmx = -2

    ! check if in current shell
    current_shell: do itr = 2, nni
        dv = xyz0(:,itr) - xyzi
        vdot = unit_dot(idrho, dv, 3)
        if (vdot > vdotmx) then
            iasgn = inrst(itr)
            vdotmx = vdot
            iiout = 1
        end if
    end do current_shell

    ! check if in other shell
    other_shell: do itr = 1, nni
        dv = xyz1(:,itr) - xyzi
        vdot = unit_dot(idrho, dv, 3)
        if (vdot > vdotmx) then
            iasgn = inrst(itr)
            vdotmx = vdot
            iiout = 2
        end if
    end do other_shell


    if (vdotmx == -2) then
        write(out,'("Assignment not found for current point")')
        stop "assign_point: Assignment not found"
    end if

    ! switch iout index if current shell is outer shell
    if (outer) then
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
