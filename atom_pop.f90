program atom_pop
  use accuracy
  use import_gamess

  real(rk)              :: xyz(3,5000)
  real(rk)              :: rho(5000)
  !
  integer(ik)            :: ipt, ird, imo, npt, nmo, nbas, rdm_count, natom
  real(ark), allocatable :: basval(:,:,:)
  real(rk), allocatable  :: moval (:,:), rdm_sv(:)
  type(gam_structure)    :: mol                       ! Gamess basis set and orbital data
  character(len=100)     :: rdm_file = 'rdm.dat'        ! Name of the file containing rdm orbital coefficients
  real(rk)               :: tmp_rdm_sv(1000),xyzq(4,6)


  call accuracyInitialize 

  call read_grid('grid_pts.dat', xyz, npt)

  call gamess_load_rdmsv(trim(rdm_file),tmp_rdm_sv,rdm_count)
  write (out,"( 'Found ',i4,' singular values')") rdm_count
  write (out,"( 'Values are: '/)")
  write (out,"(10(1x,f12.8))") tmp_rdm_sv(:rdm_count)
      !
  allocate (rdm_sv(rdm_count))
  rdm_sv = tmp_rdm_sv(:rdm_count)
  !
  call gamess_load_orbitals(file=trim(rdm_file),structure=mol)

  call gamess_report_nuclei(natom,xyzq,mol)

  nmo  = mol%nvectors
  nbas = mol%nbasis
  !
  allocate (basval(1,nbas,npt),moval(nmo,npt))

  !
  !  First, evaluate basis functions
  !
  evaluate_basis_functions: do ipt=1,npt
    call gamess_evaluate_functions(xyz(:,ipt),basval(:,:,ipt),mol)
  end do evaluate_basis_functions

  !
  !  Transform AOs to the MOs, for all grid points simultaneously
  !
  moval = matmul(transpose(mol%vectors(:,:nmo)),basval(1,:,:))
  !
  !  Finally, evaluate the transition density at grid points
  !
  rho = 0
  evaluate_rdm: do ird=1,rdm_count
    imo = 2*ird - 1
    rho(:npt) = rho(:npt) + rdm_sv(ird) * moval(imo,:npt) * moval(imo+1,:npt)
  end do evaluate_rdm
  !
  deallocate (basval,moval)

  call punch_grid('density.dat', npt, xyz(:,1:npt), rho(1:npt))

end program atom_pop 

!
!
!
subroutine read_grid(grid_file, xyz, npts)
   use accuracy 

   character(len=*), intent(in)   :: grid_file
   real(rk), intent(inout)        :: xyz(3,5000)
   integer(ik), intent(inout)     :: npts
   real(rk)                       :: pt(3)
   integer(ik)                    :: gfile = 10
   integer(ik)                    :: ios, pt_max
   logical                        :: pts_remain

   open(gfile,file=trim(grid_file),action='read',position='rewind',status='old',iostat=ios)

   pts_remain = .true.
   npts = 0
   pt_max = size(xyz,dim=2)
   xyz = -1000.
 
   scan_lines: do while (pts_remain)
     read(gfile,*,iostat=ios)pt
     if(ios/=0) then
         pts_remain = .false.
     else
         npts = npts + 1
         if(npts > pt_max) then
             pts_remain = .false.
         else
             xyz(:,npts) = pt
         endif
     endif
   end do scan_lines
   close(gfile)

   return
end subroutine read_grid

!
!
!
subroutine punch_grid(punch_file, npt, xyz, rho)
   use accuracy 

   character(len=*), intent(in)   :: punch_file
   integer(ik), intent(in)        :: npt
   real(rk), intent(in)           :: xyz(3,npt)
   real(rk), intent(in)           :: rho(npt)
   integer(ik)                    :: pfile = 11
   integer(ik)                    :: ios, ipt

   open(pfile,file=trim(punch_file),action='write',status='new',iostat=ios)

   ipt  = 1

   write_lines: do while (ipt <= npt)
       write(pfile,1000)xyz(:,ipt),rho(ipt)
       ipt = ipt + 1
   end do write_lines
   close(pfile)

   return

1000 format(f16.8,f16.8,f16.8,es18.10)
end subroutine punch_grid
 


