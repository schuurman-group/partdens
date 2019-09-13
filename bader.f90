!
!  The (hopefully temporary) Bader QTAIM charge algorithm
!
!  Calculates charges from a GAMESS RDM output file using the Bader
!  atoms-in-molecules algorithm which requires density gradients. This
!  can potentially be added to the DDCharge code if there is a good way
!  to cache unassigned points (i.e. save x, y, z, w and rho and propagate
!  the assignment until assigned).
!
!  For now, QTAIM deformation densities are calculated by gridchg.x after
!  running bader.x with inp%atom_type = pro. This is done by back-propagating
!  assignments and saving the results in a binary file.
!
program bader
    use accuracy
    use import_gamess
    use gamess_internal
    use molecular_grid
    use atoms
    use atomdens
    use fileio
    !
    character(2),allocatable :: atypes(:), gatyp(:)
    logical                  :: first_shell, last, ls
    logical,allocatable      :: adone(:,:)
    integer(ik)              :: ofile, mfile, nfile, pfile
    integer(ik)              :: i, ib, ipt, npts, nax, iat
    integer(ik)              :: nat_count, natom, ndum, nbatch
    integer(ik),allocatable  :: nnn(:), nmask(:,:), asgn(:,:)
    real(rk)                 :: norm, tmp_nat_occ(1000)
    real(rk),pointer         :: xyzw_now(:,:), xyzw_nxt(:,:)
    real(rk),allocatable     :: xyzq(:,:), xyzg(:,:), pt_dist(:), xyz(:,:,:)
    real(rk),allocatable     :: wgt_now(:), charge(:)
    real(ark),allocatable    :: rhomol(:,:), wrho(:,:), drho(:,:,:), nat_occ(:)
    type(gam_structure)      :: mol
    type(mol_grid)           :: den_grid
    type(input_data)         :: inp

    call accuracyInitialize

    ofile=10
    mfile=11
    nfile=12
    pfile=13

    open(ofile, file="bader.out")
    call read_input(input, inp)
    call init_gridchg_output(ofile, inp)
    !call read_input(infile, nat_file, vectyp, nrad, nang, atyp, alib, narad, ityp, iordr, wtyp, maxchg, maxiter, thrsh)
    !call init_output(nat_file, vectyp, nrad, nang, atyp, alib, narad, ityp, iordr, wtyp, maxchg, maxiter, thrsh)

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
    call GridInitialize(den_grid, inp%n_rad, inp%n_ang, xyzg(:,1:ndum), gatyp(1:ndum))
    call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
    open(mfile, file='moldens', form='unformatted', action='write')
    open(pfile, file='dens.mol', action='write')
    !
    !  Set initial values
    !
    nullify(xyzw_now, xyzw_nxt)
    call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw_now)
    npts = size(xyzw_now, dim=2)
    allocate (rhomol(2,npts), wrho(2,npts), drho(2,3,npts))
    allocate (nnn(npts), nmask(7,npts), pt_dist(npts), wgt_now(npts))
    allocate (xyz(2,3,npts), asgn(2,npts), adone(2,npts))
    !
    !  Find nearest neighbours
    !
    write(ofile,'("Finding nearest neighbours in spherical grid")')
    write(ofile,'("")')
    find_nearest: do ipt = 1,npts
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
    call evaluate_density_gradients(inp%vec_type, nat_count, npts, mol, nat_occ, xyzw_now(1:3,:), rhomol(2,:), drho(2,:,:))
    xyz(2,:,:) = xyzw_now(1:3,:)
    wrho(2,:) = rhomol(2,:) * xyzw_now(4,:)
    norm = sum(wrho(2,:))
    charge(:) = xyzq(4,:)
    asgn(:,:) = 0
    adone(:,:) = .false.
    first_shell = .true.
    last = .false.
    iat = den_grid%next_part
    mol_grid_batches: do ib = 1,nbatch-1
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
                if (iat /= den_grid%next_part) then
                    ls = .true.
                else
                    ls = .false.
                end if
            end if
            !
            !  Evaluate density and density gradients at grid points
            !
            call evaluate_density_gradients(inp%vec_type, nat_count, npts, mol, nat_occ, xyz(2,:,:), rhomol(2,:), drho(2,:,:))
            !
            !  Assign the densities for a set of two shells
            !
            wrho(2,:) = rhomol(2,:) * wgt_now
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
            rhomol(1,:) = rhomol(2,:)
            wrho(1,:) = wrho(2,:)
            drho(1,:,:) = drho(2,:,:)
            asgn(1,:) = asgn(2,:)
            asgn(2,:) = 0
            adone(1,:) = adone(2,:)
            adone(2,:) = .false.
            iat = den_grid%next_part
            if (ls) then
                first_shell = .true.
            end if
        end if
        !
        !  Integrate and save
        !
        rewind mfile
        mol_integrate: do ipt = 1,npts
            write(mfile) rhomol(1,ipt), asgn(1,ipt)
            write(pfile,1000) wgt_now(ipt), xyz(1,:,ipt), rhomol(1,ipt)
        end do mol_integrate
        call system("cat moldens snedlom > snedlom.tmp")
        call system("mv snedlom.tmp snedlom")
    end do mol_grid_batches
    write(ofile,'("Total molecular density: ",f14.8)') norm
    write(ofile,'("")')
    write(ofile,'("Total charge: ",f14.8)') sum(charge)
    write(ofile,'("Final charges:")')
    write(ofile,'("    ",5f14.8)') charge
    ! write the last shell
    rewind mfile
    mol_integrate_last: do ipt = 1,npts
        if (asgn(2,ipt) == 0) asgn(2,ipt) = -1
        write(mfile) rhomol(2,ipt), -asgn(2,ipt)
        write(pfile,1000) xyzw_nxt(4,ipt), xyz(2,:,ipt), rhomol(2,ipt)
    end do mol_integrate_last
    call system("cp moldens moldens.all")
    !
    !  Propagate assignments backwards
    !
    open(nfile, file='snedlom', form='unformatted', action='read')
    propagate_back: do i = 1,nbatch-1
        ! if any point isn't assigned to an atom, it must point outwards
        read_rho: do ipt = 1,npts
            read(nfile) rhomol(1,ipt), asgn(1,ipt)
        end do read_rho
        rewind mfile
        prop_bk_shell: do ipt = 1,npts
            if (asgn(1,ipt) > 0) then
                asgn(1,ipt) = asgn(2,asgn(1,ipt))
            else if (asgn(1,ipt) == 0) then
                ! this should only happen if weight*density is zero
                asgn(1,ipt) = -1
            end if
            write(mfile) rhomol(1,ipt), -asgn(1,ipt)
        end do prop_bk_shell
        call system("cat moldens moldens.all > moldens.tmp")
        call system("mv moldens.tmp moldens.all")
        rhomol(2,:) = rhomol(1,:)
        asgn(2,:) = asgn(1,:)
    end do propagate_back
    call system("mv moldens.all moldens")
    !
    !  Clean up
    !
    call GridDestroy(den_grid)
    close(mfile)
    close(nfile)
    close(pfile)
    call system("rm snedlom")
    if (trim(inp%atom_type) == "pop") then
        call system("rm moldens")
    end if
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')
    close(ofile)

    1000 format(e24.8,f16.8,f16.8,f16.8,e20.10)

contains

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
subroutine evaluate_density_gradients(vtyp, nat_count, npt, mol, nat_occ, xyz, rho, drho)
    character(6),intent(in) :: vtyp
    integer(ik),intent(in)    :: nat_count, npt
    type(gam_structure),intent(inout)   :: mol
    real(rk),intent(in)       :: nat_occ(nat_count)
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
            evaluate_rdm: do ird=1,nat_count
                imo = 2*ird - 1
                rho = rho + nat_occ(ird) * moval(imo,:) * moval(imo+1,:)
                do ic = 1, 3
                    drho(ic,:) = drho(ic,:) + nat_occ(ird) * (moval(imo,:) * dmoval(ic,imo+1,:) + moval(imo+1,:) * dmoval(ic,imo,:))
                end do
            end do evaluate_rdm
        case ("natorb")
            evaluate_nat: do ird=1,nat_count
                rho = rho + nat_occ(ird) * moval(ird,:)**2
                do ic = 1, 3
                    drho(ic,:) = drho(ic,:) + 2 * nat_occ(ird) * moval(ird,:) * dmoval(ic,ird,:)
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
    integer(ik),parameter :: maxitr=120
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
                if (adone(j,i) .or. askip(j,i)) cycle propagate_shell
                o_ji = iout(j,i)
                a_ji = asgn(j,i)
                o_next = iout(o_ji,a_ji)
                a_next = asgn(o_ji,a_ji)
                if (o_next == j .and. a_next == i) then
                    ! repetition
                    ! try to assign to the same point as a neighbour
                    chk_nbours0: do n = 2, nnn(i)
                        if (asgn(j,nrst(n,i)) < 0) then
                            iatm = -asgn(j,nrst(n,i))
                            asgn(j,i) = -iatm
                            atmchg(iatm) = atmchg(iatm) - wrho(j,i)
                            adone(j,i) = .true.
                            exit chk_nbours0
                        else if (asgn(3-j,nrst(n,i)) < 0) then
                            iatm = -asgn(3-j,nrst(n,i))
                            asgn(j,i) = -iatm
                            atmchg(iatm) = atmchg(iatm) - wrho(j,i)
                            adone(j,i) = .true.
                            exit chk_nbours0
                        end if
                    end do chk_nbours0
                    chk_nbours1: do n = 1, nnn(a_ji)
                        if (asgn(o_ji,nrst(n,a_ji)) < 0) then
                            iatm = -asgn(o_ji,nrst(n,a_ji))
                            asgn(o_ji,a_ji) = -iatm
                            atmchg(iatm) = atmchg(iatm) - wrho(o_ji,a_ji)
                            adone(o_ji,a_ji) = .true.
                            exit chk_nbours1
                        else if (asgn(3-o_ji,nrst(n,a_ji)) < 0) then
                            iatm = -asgn(3-o_ji,nrst(n,a_ji))
                            asgn(o_ji,a_ji) = -iatm
                            atmchg(iatm) = atmchg(iatm) - wrho(o_ji,a_ji)
                            adone(o_ji,a_ji) = .true.
                            exit chk_nbours1
                        end if
                    end do chk_nbours1
                    ! otherwise assign to the nearest nucleus
                    if (.not. adone(j,i)) then
                        do n = 1, natm
                            dist(n) = sqrt(sum((xyz(j,:,i) - xyzq(1:3,n))**2))
                        end do
                        iatm = minloc(dist, 1)
                        asgn(j,i) = -iatm
                        atmchg(iatm) = atmchg(iatm) - wrho(j,i)
                        adone(j,i) = .true.
                    end if
                    if (.not. adone(o_ji,a_ji)) then
                        do n = 1, natm
                            dist(n) = sqrt(sum((xyz(o_ji,:,a_ji) - xyzq(1:3,n))**2))
                        end do
                        iatm = minloc(dist, 1)
                        atmchg(iatm) = atmchg(iatm) - wrho(o_ji,a_ji)
                        asgn(o_ji,a_ji) = -iatm
                        adone(o_ji,a_ji) = .true.
                    end if
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
        print *,"drho = ",idrho
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

!
!  Evaluate the vector v minus its projection on a vector u
!
function projout(u, v, nd)
    integer(ik) :: nd
    real(rk)    :: u(nd), v(nd), projout(nd)
    !
    real(rk)    :: udotv, vdotv

    udotv = sum(u * v)
    vdotv = sum(v**2)
    projout = u - (udotv/vdotv) * v
end function projout

end program bader
