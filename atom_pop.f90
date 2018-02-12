program atom_pop
  use accuracy
  use import_gamess
  use gamess_internal
  use molecular_grid
  use atomdens
  !
  integer(ik)            :: i, ib, ipt, npts, pfile, rdm_count, natom, nbatch
  real(ark), allocatable :: rho(:), rdm_sv(:)
  type(gam_structure)    :: mol                 ! Gamess basis set and orbital data
  type(mol_grid)         :: den_grid
  character(len=100)     :: rdm_file='rdm.dat'  ! Name of the file containing rdm orbital coeffs
  real(rk)               :: tmp_rdm_sv(1000),overlap
  real(rk),allocatable   :: xyzq(:,:)
  integer(ik)            :: nrad=70, nang=110
  character(len=2),allocatable :: atypes(:)
  real(rk), pointer      :: xyzw(:,:)           ! Coordinates and weights of the grid points

  call accuracyInitialize

  pfile=11
  open(pfile, file='density.dat', action='write', status='new', iostat=ios)

  call gamess_load_rdmsv(trim(rdm_file), tmp_rdm_sv, rdm_count)
  write (out,"( 'Found ',i4,' singular values')") rdm_count
  write (out,"( 'Values are: '/)")
  write (out,"(10(1x,f12.8))") tmp_rdm_sv(:rdm_count)
  !
  allocate (rdm_sv(rdm_count))
  rdm_sv = tmp_rdm_sv(:rdm_count)
  !
  call gamess_load_orbitals(file=trim(rdm_file), structure=mol)

  allocate (atypes(mol%natoms), xyzq(4,mol%natoms))

  call gamess_report_nuclei(natom, xyzq, mol)
  do i = 1,mol%natoms
    atypes(i) = trim(mol%atoms(i)%name)
  enddo
  !
  call GridInitialize(den_grid, nrad, nang, xyzq(1:3,:), atypes)

  call GridPointsBatch(den_grid, 'Batches count', count=nbatch)
  !
  !  Numerical integration loop
  !
  overlap = 0.
  nullify(xyzw)

  grid_batches: do ib = 1,nbatch
     !
     !  Get grid points
     !
     call GridPointsBatch(den_grid, 'Next batch', xyzw=xyzw)
     !
     !  If the batch size changed, reallocate rho()
     !
     npts = size(xyzw, dim=2)
     if (allocated(rho)) then
       if (size(rho) /= npts) deallocate (rho)
     end if
     if (.not. allocated(rho)) allocate (rho(npts))
     !
     !  Evaluate density at grid points
     !
     call evaluate_density(rdm_count, npts, mol, rdm_sv, xyzw(1:3,:), rho)
     !
     !  Evaluate total atomic density at grid points
     !
     !!call evaluate_atomic(xyzq, natoms, xyzw(1:3,:), npts, rhoatm)
     !
     !  Determine Voronoi assignments
     !
     !! indi = assign_voronoi(xyzq(1:3,:), natoms, xyzw(1:3,:), npts)
     !
     !  Ready to integrate!
     !
     integrate: do ipt = 1,npts
       overlap = overlap + xyzw(4,ipt) * rho(ipt)
       write(pfile,1000) xyzw(4,ipt), xyzw(1:3,ipt), rho(ipt)
     end do integrate
  end do grid_batches

  close(pfile)

  call GridDestroy(den_grid)

  print *,'VALUE OF OVERLAP: ', overlap

1000 format(e24.8,f16.8,f16.8,f16.8,es18.10)
end program atom_pop

!
! Determine the density at XYZ coordinates
!
subroutine evaluate_density(rdm_count, npt, mol, rdm_sv, xyz, rho)
    use accuracy
    use import_gamess
    use gamess_internal

    integer(ik),intent(in) :: rdm_count, npt
    type(gam_structure), intent(inout)   :: mol
    real(rk), intent(in)   :: rdm_sv(rdm_count)
    real(rk), intent(in)   :: xyz(3,npt) ! XYZ coordinates of the grid points
    real(rk), intent(out)  :: rho(npt)   ! Transition density at these grid points
    !
    integer(ik)            :: ipt, ird, imo, nmo, nbas
    real(ark), allocatable :: basval(:,:,:)
    real(rk), allocatable  :: moval (:,:)

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
    rho(:) = 0
    evaluate_rdm: do ird = 1,rdm_count
      imo = 2*ird - 1
      rho = rho + rdm_sv(ird) * moval(imo,:) * moval(imo+1,:)
    end do evaluate_rdm
    !
    deallocate (basval,moval)

end subroutine evaluate_density
