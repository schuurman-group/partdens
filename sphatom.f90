!
!  The spherical atomic density program, used to generate
!  atomic density library files
!
program sphatom
    use accuracy
    use import_gamess
    use gamess_internal
    use atomdens
    use lebedev
    !
    character(100),parameter :: rdm_file="rdm.dat", vectyp="natorb"
    integer(ik)              :: ib, ios, achg, pfile, rdm_count, natom
    integer(ik),parameter    :: nrad=70, nang=770
    real(rk)                 :: tmp_rdm_sv(1000), xyzq(4,1)
    real(rk),allocatable     :: r(:), w(:), xyz(:,:)
    real(rk),pointer         :: ang_xyzw(:,:)
    real(ark),allocatable    :: rho(:), rdm_sv(:), rden(:)
    type(gam_structure)      :: mol

    call accuracyInitialize

    pfile=11
    open(pfile,file="raddens.dat",action="write",iostat=ios)

    allocate (r(nrad), w(nrad), rden(nrad), xyz(3,nang), rho(nang))

    call gamess_load_rdmsv(trim(rdm_file),tmp_rdm_sv,rdm_count)
    write(out,'("Found ",i4," singular values")') rdm_count
    write(out,'("Values are: "/)')
    write(out,'(5(1x,f12.8))') tmp_rdm_sv(:rdm_count)
    write(out,'("")')
    !
    allocate (rdm_sv(rdm_count))
    rdm_sv = tmp_rdm_sv(:rdm_count)
    !
    call gamess_load_orbitals(file=trim(rdm_file),structure=mol)
    !
    if (mol%natoms > 1) then
        write(out,'("sphatom: Number of atoms ",i4," is greater than 1.")') mol%natoms
        stop "sphatom - more than one atom"
    end if

    call gamess_report_nuclei(natom,xyzq,mol)
    call rlegendre(nrad, int(xyzq(4,1)), r, w)
    achg = nint(xyzq(4,1) - sum(rdm_sv))
    write(pfile,'(a2,sp,i3,ss,i6)') mol%atoms(1)%name, achg, nrad

    select case (nang)
        case default
            write (out,'("rad_shell: Angular grid ",i6," is not recognized.")') nang
            stop "rad_shell - bad angular grid"
        case (110) ; ang_xyzw => lebedev_gr17
        case (302) ; ang_xyzw => lebedev_gr29
        case (770) ; ang_xyzw => lebedev_gr47
    end select
    !
    !  Numerical integration loop
    !
    rden  = 0_ark
    rad_shell: do ib=1,nrad
        rho = 0_ark
        !
        !  Get grid points
        !
        xyz = spread(xyzq(1:3,1), dim=2, ncopies=nang) + r(ib) * ang_xyzw(1:3,:)
        !
        !  Evaluate density at grid points
        !
        call evaluate_density(vectyp, rdm_count, nang, mol, rdm_sv, xyz, rho)
        !
        !  Integrate the spherical shell
        !
        rden(ib) = sum(rho*ang_xyzw(4,:))
        if (rden(ib) < 1e-99_ark) then
            rden(ib) = 0_ark
        end if
        write(pfile,'(3e20.10)') r(ib), rden(ib), w(ib)
    end do rad_shell

    write(pfile,'("")')
    close(pfile)

    write(out,'("Total atomic density: ",f18.10)') sum(rden*w)
    write(out,'("Charge: ",sp,i3)') achg

contains

!
!  Evaluate the electron density at a set of grid points
!
subroutine evaluate_density(vtyp, rdm_count, npt, mol, rdm_sv, xyz, rho)
    character(100),intent(in)         :: vtyp
    integer(ik),intent(in)            :: rdm_count, npt
    type(gam_structure),intent(inout) :: mol
    real(rk),intent(in)               :: rdm_sv(rdm_count)
    real(rk),intent(in)               :: xyz(3,npt)
    real(ark),intent(out)             :: rho(npt)
    !
    integer(ik)            :: ipt, ird, imo, nmo, nbas
    real(ark), allocatable :: basval(:,:,:)
    real(rk), allocatable  :: moval (:,:)

    nmo  = mol%nvectors
    nbas = mol%nbasis
    !
    allocate (basval(1,nbas,npt),moval(nmo,npt))
    !
    !  First, evaluate basis functions
    !
    evaluate_basis_functions: do ipt=1,npt
        call gamess_evaluate_functions(xyz(:,ipt),basval(:,:,ipt), mol)
    end do evaluate_basis_functions
    !
    !  Transform AOs to the MOs, for all grid points simultaneously
    !
    moval = matmul(transpose(mol%vectors(:,:nmo)),basval(1,:,:))
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

end program sphatom
