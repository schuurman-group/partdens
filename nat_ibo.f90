!
!  General intrinsic atomic orbitals and intrinsic bond orbitals
!
!  Uses a modified IAO and IBO generation algorithm which extends to
!  natural orbitals.
!
program ibocalc
    use accuracy
    use import_gamess
    use gamess_internal
    use lapack
    use atoms
    use molden
    !
    character(100),parameter :: infile="ibo.inp", outfile="ibo.out"
    character(100)           :: nat_file="nos.dat", extbas
    character(2),allocatable :: atypes(:)
    integer(ik)              :: ofile, i, j, iat, imin, pwr
    integer(ik)              :: n1bas, n2bas, nvec, natom, nsh, nmin
    integer(ik),allocatable  :: mask(:)
    real(rk)                 :: tmp_nat_occ(1000)
    real(rk),allocatable     :: xyzq(:,:), charge(:)
    real(ark),allocatable    :: P12(:,:), P21(:,:), S(:,:), S1(:,:), S2(:,:), Oi(:,:), Ot(:,:)
    real(ark),allocatable    :: Ci(:,:), Ct(:,:), Ciao(:,:), Cibo(:,:), ener(:)
    real(ark),allocatable    :: nat_occ(:), n(:,:), nibo(:,:), nibodiag(:), nibooff(:)
    real(ark),allocatable    :: Iden(:,:), IA(:,:), IB(:,:), IBat(:,:)
    real(ark),allocatable    :: W(:,:,:), Wp(:,:,:), U(:,:), Tt(:,:), nt(:,:)
    logical                  :: localize_virt
    type(gam_structure)      :: mol1, mol2

    call accuracyInitialize

    ofile = 10
    extbas = '/globalhome/rymac/Projects/PartialCharge/partdens/atomlib/extbas'

    localize_virt = .false.
    pwr = 2

    open(ofile, file=outfile)
    !call read_input(infile, rdm_file, vectyp, nrad, nang, cmax, cmin, ccons, dxyz)
    call init_output() !rdm_file, vectyp, nrad, nang, cmax, cmin, ccons, dxyz)
    write(ofile,'("Loading GAMESS data file")')
    call gamess_load_natocc(trim(nat_file), tmp_nat_occ, nvec)
    write(ofile,'("Found ",i4," natural occupations")') nvec
    write(ofile,'("Values are:")')
    write(ofile,'(5(1x,f12.8))') tmp_nat_occ(:nvec)
    write(ofile,'("")')
    !
    allocate (nat_occ(nvec))
    nat_occ = tmp_nat_occ(:nvec)
    !
    call gamess_load_orbitals(file=trim(nat_file), structure=mol1)
    !
    allocate (atypes(mol1%natoms), xyzq(4,mol1%natoms), charge(mol1%natoms))
    !
    call gamess_report_nuclei(natom, xyzq, mol1)
    do i = 1,mol1%natoms
        atypes(i) = trim(mol1%atoms(i)%name)
    enddo
    write(ofile,'("Molecular geometry (in bohr):")')
    do iat = 1,mol1%natoms
        write(ofile,'("    ",a2,3f14.8)') atypes(iat), xyzq(1:3,iat)
    end do
    write(ofile,'("")')
    !
    !  Import the minimal basis
    !
    mol2 = mol1
    call gamess_reload_basis(file=extbas, basname='MINAO', structure=mol2)
    !
    !  Calculate overlaps and projectors
    !
    n1bas = mol1%nbasis
    n2bas = mol2%nbasis
    allocate (P12(n1bas,n2bas), P21(n2bas,n1bas))
    call gamess_1e_integrals('AO OVERLAP', P12, mol1, mol2)
    P21 = transpose(P12)
    allocate (S(n2bas,n2bas), S2(n2bas,n2bas))
    call gamess_1e_integrals('AO OVERLAP', S, mol2, mol2)
    call lapack_dposv(S, P21) ! S2 only used to find P21
    S2 = S
    deallocate (S)
    allocate (S(n1bas,n1bas), S1(n1bas,n1bas))
    call gamess_1e_integrals('AO OVERLAP', S, mol1, mol1)
    S1 = S                    ! S1 used throughout
    call lapack_dposv(S, P12)
    deallocate (S)
    allocate (Oi(n1bas,n1bas), Ot(n2bas,n2bas), IA(n1bas,n2bas), Iden(n1bas,n1bas))
    allocate (Ci(n1bas,nvec), Ct(n2bas,nvec), Ciao(n2bas,nvec), Cibo(n1bas,nvec))
    allocate (n(nvec,nvec), nt(nvec,nvec), nibo(nvec,nvec), nibodiag(nvec), nibooff(nvec), ener(nvec))
    n  = 0_ark
    do i = 1,nvec
        n(i,i)  = tmp_nat_occ(i)
    end do
    Ci = mol1%vectors
    !! alternate (simple) formulation
    !Cii = matmul(S21, Ci)
    !Ct = Cii
    !call lapack_dposv(S2, Ct)
    !St = matmul(transpose(Cii), Ct)
    !Cti = Ct
    !call lapack_dposv(St, transpose(Cti))
    ! ------ error in charge between here ------ !
    !Ct = matmul(matmul(matmul(P12, P21), Ci), sqrt(0.5_ark*n))
    ! mulliken total population of Ct before orthogonalization (test)
    !Ot = matmul(Ct, matmul(n, transpose(Ct))) * S1
    !Ot = matmul(Ct, transpose(Ct)) * S1
    !print *,"Mull(Ct) = ",sum(Ot)," before orth"
    allocate (Tt(nvec,nvec), U(nvec,nvec))
    !call orth(Ct, S1, Tt, n1bas, nvec)
    !Ct = matmul(Ct, Tt)
    !Ct = matmul(matmul(Ct, sqrt(n)), Tt)
    Ct = matmul(P21, Ci)
    call orth(Ct, S2, Tt, n2bas, nvec)
    Ct = matmul(Ct, Tt)
    nt = n
    deallocate (Tt)
    !! decompose the transformation to remove renormalization factors
    !call lapack_dgesvd(Tt, ener, U, VTH)
    !do i = 1,nvec
    !    if (ener(i) < 1e-8) U(:,i) = 0_ark
    !end do
    !U = matmul(U, VTH)
    !deallocate (Tt)
    !nt = matmul(transpose(U), matmul(n, U))
    ! mulliken total population of Ci and Ct (test)
    Oi = matmul(Ci, matmul(n, transpose(Ci))) * S1
    Ot = matmul(Ct, matmul(nt, transpose(Ct))) * S2
    !Ot = matmul(P21, matmul(Oi, P12))
    print *,"Mull(Ci) = ",sum(Oi)
    print *,"Mull(Ct) = ",sum(Ot)
    ! assuming spin-free natural orbitals, need to divide densities by 2
    Oi = 0.5_ark*matmul(matmul(Ci, matmul(n, transpose(Ci))), S1)
    !Ot = 0.5_ark*matmul(matmul(Ct, matmul(nt, transpose(Ct))), S1)
    Ot = 0.5_ark*matmul(matmul(Ct, matmul(nt, transpose(Ct))), S2)
    !
    !  Calculate IAOs
    !
    Iden = 0_ark
    do i = 1,n1bas
        Iden(i,i) = 1_ark
    end do
    !IA = matmul(matmul(Oi, Ot) + matmul(Iden - Oi, Iden - Ot), P12)
    IA = matmul(Oi, P12) + P12 - matmul(P12, Ot)
    allocate (Tt(n2bas,n2bas))
    call orth(IA, S1, Tt, n1bas, n2bas)
    IA = matmul(IA, Tt)
    deallocate (Tt)
    Ciao = matmul(transpose(IA), matmul(S1, Ci))
    ! ---------------- and here ---------------- !
    deallocate (Iden)
    !
    !  Calculate charges and weighted population matrices W_a
    !
    allocate (IB(n2bas,nvec), W(nvec,nvec,natom))
    IB  = matmul(Ciao, sqrt(n))
    W   = 0_ark
    charge = xyzq(4,:)
    imin = 1
    do iat = 1,natom
        nsh = mol2%atoms(iat)%nshell
        ! ignore dummy atoms
        if (nsh > 0) then
            nmin = 0
            do j = 1,nsh
                i = mol2%atoms(iat)%sh_l(j)
                nmin = nmin + ang_loc(i+1) - ang_loc(i)
            end do
            allocate (IBat(nmin,nvec))
            IBat = IB(imin:imin+nmin-1,:)
            W(:,:,iat) = matmul(transpose(IBat), IBat)
            deallocate (IBat)
            imin = imin + nmin
            ! calculate charge from Tr(W_a)
            do j = 1,nvec
                charge(iat) = charge(iat) - W(j,j,iat)
            end do
        end if
    end do
    ! print out IAO charges
    write(ofile,'("Total charge: ",f14.8)') sum(charge)
    write(ofile,'("IAO charges:")')
    write(ofile,'("    ",5f14.8)') charge
    write(ofile,'("")')
    !
    !  Calculate optimal IBO rotation matrix from W
    !
    U = 0_ark
    do i = 1,nvec
        U(i,i) = 1_ark
    end do
    ! localize unoccupied if desired (counterproductive)
    if (localize_virt) then
        write(ofile,'("Starting IBO optimization in unoccupied space")')
        allocate (Wp(nvec,nvec,natom))
        Wp  = 0_ark
        IB  = matmul(Ciao, sqrt(2_ark*U - n))
        imin = 1
        do iat = 1,natom
            nsh = mol2%atoms(iat)%nshell
            ! ignore dummy atoms
            if (nsh > 0) then
                nmin = 0
                do j = 1,nsh
                    i = mol2%atoms(iat)%sh_l(j)
                    nmin = nmin + ang_loc(i+1) - ang_loc(i)
                end do
                allocate (IBat(nmin,nvec))
                IBat = IB(imin:imin+nmin-1,:)
                Wp(:,:,iat) = matmul(transpose(IBat), IBat)
                deallocate (IBat)
                imin = imin + nmin
            end if
        end do
        call localize_mos(U, Wp, pwr, nvec, natom, ofile)
        deallocate (Wp)
        do iat = 1,natom
            W(:,:,iat) = matmul(transpose(U), matmul(W(:,:,iat), U))
        end do
    end if
    deallocate (IB)
    ! then localize occupied
    write(ofile,'("Starting IBO optimization in occupied space")')
    call localize_mos(U, W, pwr, nvec, natom, ofile)
    ! transform NOs and populations to IBOs and IBO populations
    Cibo = matmul(Ci, U)
    nibo = matmul(transpose(U), matmul(n, U))
    ! calculate diagonal and off-diagonal populations
    do i = 1,nvec
        nibodiag(i) = nibo(i,i)
        nibooff(i) = sum(nibo(:,i)) - nibo(i,i)
    end do
    ! sort by diagonal populations
    allocate (mask(nvec))
    mask = argsort(nibodiag, nvec)
    Cibo = Cibo(:,mask)
    nibodiag = nibodiag(mask)
    nibooff = nibooff(mask)
    deallocate (mask)
    ! print IBO populations
    write(ofile,'("IBO diagonal populations:")')
    write(ofile,'(5(1x,f12.8))') nibodiag
    write(ofile,'("")')
    write(ofile,'("Sum of IBO off-diagonal populations:")')
    write(ofile,'(5(1x,f12.8))') nibooff
    write(ofile,'("")')
    ! output IBOs to molden format
    write(ofile,'("Writing IBOs to file ibo.molden")')
    ener = 0_rk
    call write_molden(mol1, 'ibo.molden', n1bas, nvec, Cibo, ener, nibodiag)
    !
    !  Clean up
    !
    write(ofile,'("")')
    write(ofile,'("Exited successfully")')

    1000 format(f16.8,f16.8,f16.8,e20.10,e20.10)

contains

!
!  Return a sorted mask from an array (in reverse order)
!
function argsort(a, na)
    integer(ik),intent(in) :: na
    real(ark),intent(in)   :: a(na)
    integer(ik)            :: argsort(na)
    !
    integer(ik)            :: i, tempi
    real(ark)              :: a2(na), tempr

    a2 = a
    do i = 1,na
        argsort(i) = i
    end do
    do i = 1,na-1
        imin = maxloc(a2(i:),1) + i - 1
        if (imin /= i) then
            tempr = a2(i)
            a2(i) = a2(imin)
            a2(imin) = tempr
            tempi = argsort(i)
            argsort(i) = argsort(imin)
            argsort(imin) = tempi
        end if
    end do

end function argsort

!
!  Symmetrically orthogonalize a matrix (Lowdin)
!
subroutine orth(A, S, U, na, nb)
    integer(ik),intent(in) :: na, nb
    real(ark),intent(in)   :: A(na,nb), S(na,na)
    real(ark),intent(out)  :: U(nb,nb)
    !
    integer(ik)            :: i
    real(ark)              :: V(nb), W(nb,nb)

    U = matmul(transpose(A), matmul(S, A))
    call lapack_dsyev(U, V)
    W = 0_ark
    do i = 1,nb
        if (V(i) < 1e-8) then
            write(out,'("orth - skipping eigenvalue",i3", ",e12.4)') i, V(i)
        else
            W(:,i) = U(:,i) / sqrt(V(i))
        end if
    end do
    U = matmul(W, transpose(U))

end subroutine orth

!
!  Localize molecular orbitals based on a weighted population matrix W
!
subroutine localize_mos(U, W, p, nv, nat, ofile)
    integer(ik),intent(in)  :: nv, nat, p, ofile
    real(ark),intent(inout) :: U(nv,nv), W(nv,nv,nat)
    !
    integer(ik)             :: i, j, itr, maxiter
    real(rk)                :: grad, Loc
    real(ark)               :: Aij, Bij, phi, ss, cs
    real(ark)               :: Wii(nat), Wij(Nat), Wjj(nat)
    real(ark)               :: Wi(nv,nat), Wj(nv,nat), Ui(nv), Uj(nv)

    maxiter = 100
    Loc = 0_rk
    do i = 1,nvec
        Loc = Loc + sum(W(i,i,:)**p)
    end do
    write(ofile,'("           L = ",f14.8)') Loc
    write(ofile,'("")')
    !
    mo_rotations: do itr = 1,maxiter
        write(ofile,'("ITER = ",i4)') itr
        grad = 0_rk
        do j = 2,nvec
            do i = 1,j-1
                Wii = W(i,i,:)
                Wij = W(i,j,:)
                Wjj = W(j,j,:)
                if (all(Wij < 1e-10)) continue
                if (p == 2) then
                    Aij = sum((Wii - Wjj)**2 - 4*Wij**2)
                    Bij = sum(4*Wij*(Wii - Wjj))
                else if (p == 4) then
                    Aij = sum(Wii**4 + Wjj**4 - 6*(Wii**2+Wjj**2)*Wij**2 - Wii**3*Wjj - Wii*Wjj**3)
                    Bij = sum(4*Wij*(Wii**3 - Wjj**3))
                else
                    write(out,'("mo_rotations: localization power ",i3," not supported")')
                    stop "mo_rotations - localization power not supported"
                end if
                phi = 0.25*atan2(Bij,Aij)
                ss = sin(phi)
                cs = cos(phi)
                ! update W
                Wi = W(:,i,:)
                Wj = W(:,j,:)
                W(:,i,:) =  Wi*cs + Wj*ss
                W(:,j,:) = -Wi*ss + Wj*cs
                Wi = W(i,:,:)
                Wj = W(j,:,:)
                W(i,:,:) =  Wi*cs + Wj*ss
                W(j,:,:) = -Wi*ss + Wj*cs
                ! update U
                Ui = U(:,i)
                Uj = U(:,j)
                U(:,i) =  Ui*cs + Uj*ss
                U(:,j) = -Ui*ss + Uj*cs
                ! update the gradient
                grad = grad + Bij**2
            end do
        end do
        Loc = 0_rk
        do i = 1,nvec
            Loc = Loc + sum(W(i,i,:)**p)
        end do
        grad = sqrt(grad)
        write(ofile,'("    gradient = ",f14.8)') grad
        write(ofile,'("           L = ",f14.8)') Loc
        write(ofile,'("")')
        if (grad < 1e-8) then
            write(ofile,'("IBOs converged on ITER = ",i4)') itr
            write(ofile,'("")')
            exit
        else if (itr == maxiter)then
            write(ofile,'("localize_mos: IBOs failed to converge")')
            stop "localize_mos - IBOs failed to converge"
        end if
    end do mo_rotations

end subroutine localize_mos

!
!  Read a crude input file
!
subroutine read_input(infile, nfile, vtyp, nr, na, fmax, fmin, cnst, dg)
    character(100) :: infile, nfile, vtyp
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
                nfile = var
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
subroutine init_output() !nfile, vtyp, nr, na, fmax, fmin, cnst, dg)
    !character(100) :: nfile, vtyp
    !integer(ik)    :: nr, na
    !real(rk)       :: fmax, fmin, cnst, dg

    write(ofile,'("+---------------------------------------------------+")')
    write(ofile,'("|                                                   |")')
    write(ofile,'("|                   IAO-IBO                         |")')
    write(ofile,'("|                                                   |")')
    write(ofile,'("|      Intrinsic atomic/bond orbitals from NOs      |")')
    write(ofile,'("|          RJ MacDonell, MS Schuurman 2019          |")')
    write(ofile,'("+---------------------------------------------------+")')
    write(ofile,'("")')
    write(ofile,'("")')
    !write(ofile,'("Input summary:")')
    !write(ofile,'("    ------- Molecular density --------")')
    !write(ofile,'("    nat_filename   =   ",a15)') trim(nfile)
    !write(ofile,'("    vec_type       =   ",a15)') trim(vtyp)
    !write(ofile,'("    n_r_grid       =   ",i15)') nr
    !write(ofile,'("    n_ang_grid     =   ",i15)') na
    !write(ofile,'("")')
    !write(ofile,'("    covrad_max     =   ",f15.4)') fmax
    !write(ofile,'("    covrad_min     =   ",f15.4)') fmin
    !write(ofile,'("    covrad_const   =   ",f15.4)') cnst
    !write(ofile,'("    dgrid          =   ",f15.4)') dg
end subroutine init_output

end program ibocalc
