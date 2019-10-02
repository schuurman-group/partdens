!
!  Tools for calculating molecular and atomic densities on a grid
!
module utilities
    use accuracy
    use import_gamess
    use gamess_internal
    use atoms

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
!  Determine the analytic atomic densities at grid XYZ coordinates
!
subroutine evaluate_atomic(xyz, npt, pmol, dmat, nsh, natm, rho)
    integer(ik),intent(in)            :: natm, npt
    integer(ik),intent(in)            :: nsh(natm)
    type(gam_structure),intent(inout) :: pmol
    real(rk),intent(in)               :: xyz(3,npt)
    real(ark),intent(in)              :: dmat(sum(nsh),sum(nsh))
    real(ark),intent(out)             :: rho(natm,npt)
    !
    integer(ik)           :: iat, ipt, ish
    real(ark)             :: shlval(1,sum(nsh),npt)
    real(ark),allocatable :: ptval(:), admat(:,:)

    ! evaluate basis function shells
    evaluate_shells: do ipt = 1,npt
        call gamess_evaluate_functions(xyz(:,ipt), shlval(:,:,ipt), pmol, r_only=.true.)
    end do evaluate_shells
    rho(:,:) = 0_ark
    ish = 1
    evaluate_atom: do iat = 1,natm
        ! divide shlval and dmat into atomic components and multiply
        allocate (ptval(nsh(iat)), admat(nsh(iat),nsh(iat)))
        admat = dmat(ish:ish+nsh(iat)-1,ish:ish+nsh(iat)-1)
        evaluate_pt_dens: do ipt = 1,npt
            ptval = shlval(1,ish:ish+nsh(iat)-1,ipt)
            rho(iat,ipt) = dot_product(ptval, matmul(admat, ptval))
        end do evaluate_pt_dens
        deallocate (ptval, admat)
        ish = ish + nsh(iat)
    end do evaluate_atom

end subroutine evaluate_atomic

!
!  Update the set of atomic density matrices based on charges
!
subroutine update_densmat(chg, maxc, ql, natm, uind, nu, atmden, dload, nsh, dlib, abas, dens)
    integer(ik),intent(in)  :: natm, nu, maxc
    integer(ik),intent(in)  :: nsh(natm), ql(nu), uind(natm)
    logical,intent(inout)   :: dload(2*maxc+1,nu)
    real(rk),intent(in)     :: chg(natm)
    real(ark),intent(inout) :: atmden(maxval(nsh),maxval(nsh),2*maxc+1,nu)
    character(100)          :: dlib
    character(8)            :: abas
    real(ark),intent(out)   :: dens(sum(nsh),sum(nsh))
    !
    integer(ik)             :: iat, is1, is2, il, iu, ui
    real(rk)                :: modc, cthrsh=1e-3

    is1 = 1
    update_dmat: do iat = 1,natm
        if (abs(chg(iat)) > maxc + cthrsh) then
            write(out,'("update_densmat: abs(charge) greater than MAXCHG = ",i3)') maxc
            stop "update_densmat - charge out of bounds"
        end if
        ui = uind(iat)
        is2 = is1 + nsh(iat) - 1
        if (abs(chg(iat)-nint(chg(iat))) < cthrsh) then
            ! integer charge, only get one contribution
            il = maxc + 1 + nint(chg(iat))
            if (.not. dload(il,ui)) then
                ! import the integer density matrix
                call load_atom_dmat(ql(ui), nint(chg(iat)), nsh(iat), dlib, abas, atmden(:nsh(iat),:nsh(iat),il,ui))
                dload(il,ui) = .true.
            end if
            dens(is1:is2,is1:is2) = atmden(:nsh(iat),:nsh(iat),il,ui)
        else
            ! find the modulus (as it should be defined)
            modc = chg(iat) - floor(chg(iat))
            ! real charge, get ceil(chg) and floor(chg) contributions
            il = maxc + 1 + floor(chg(iat))
            iu = maxc + 1 + ceiling(chg(iat))
            if (.not. dload(il,ui)) then
                ! import the lower integer density matrix
                call load_atom_dmat(ql(ui), floor(chg(iat)), nsh(iat), dlib, abas, atmden(:nsh(iat),:nsh(iat),il,ui))
                dload(il,ui) = .true.
            end if
            if (.not. dload(iu,ui)) then
                ! import the upper integer density matrix
                call load_atom_dmat(ql(ui), ceiling(chg(iat)), nsh(iat), dlib, abas, atmden(:nsh(iat),:nsh(iat),iu,ui))
                dload(iu,ui) = .true.
            end if
            dens(is1:is2,is1:is2) = (1-modc)*atmden(:nsh(iat),:nsh(iat),il,ui) + modc*atmden(:nsh(iat),:nsh(iat),iu,ui)
        end if
        is1 = is1 + nsh(iat)
    end do update_dmat

end subroutine update_densmat

!
!  Load the radial density matrix of an atom/basis/charge combination
!
subroutine load_atom_dmat(q, chg, nsh, dlib, abas, admat)
    integer(ik),intent(in) :: q, chg, nsh
    character(100)         :: dlib
    character(8)           :: abas
    real(ark),intent(out)  :: admat(nsh,nsh)
    !
    character(20)          :: sstr
    character(2)           :: elem
    integer(ik)            :: ios, ieq
    integer(ik)            :: row, line, ifield
    integer(ik)            :: chk_row, chk_line
    logical                :: file_exists, have_dmat
    real(ark)              :: loc_val(5)

    if (chg == q) then
        admat = 0_ark
        return
    else if (chg > q) then
        write(out,'("load_atom_dmat: Atomic charge greater than nuclear charge")')
        stop "load_atom_dmat - charge out of bounds"
    end if
    inquire(file=trim(dlib), exist=file_exists)
    ! open the radial density matrix library
    if (.not. file_exists) then
        write(out,'("load_atom_dmat: Density matrix library not found at ",a)') dlib
        stop "load_atom_dmat - density matrix file does not exist"
    end if
    write(out,'("Loading ab initio atomic density for Z = ",i3,", CHG = ",i3)') q, chg
    open (gam_file,file=trim(dlib),action='read',position='rewind',status='old',iostat=ios)
    elem = AtomElementSymbol(1._rk*q)
    ! find the radial density matrix of the atom/basis/charge
    if (chg > 0.5) then
        write (sstr,"(a2,' ',a8,sp,i3)") adjustl(elem), adjustl(abas), chg
    else
        write (sstr,"(a2,' ',a8,i3)") adjustl(elem), adjustl(abas), chg
    end if
    have_dmat = .false.
    scan_dmat: do while (.not.have_dmat)
        call gam_readline
        if (trim(gam_line_buf)==trim(sstr)) have_dmat = .true.
    end do scan_dmat
    if (.not.have_dmat) then
        write(out,'("load_atom_dmat: Density matrix ",a," not found for atom ",a,"with charge ",sp,i3)') trim(abas), trim(elem), chg
        stop "load_atom_dmat - density matrix not found"
    end if
    ! read the radial density matrix
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
          if (ieq > nsh) then
              row = row + 1
              if (row > nsh) exit read_dmat
              line = 0
              ieq = 1
              exit stuff_d
          end if
      end do stuff_d
    end do read_dmat
    close (gam_file,iostat=ios)

end subroutine load_atom_dmat

end module utilities
