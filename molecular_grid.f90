!
!  Simple non-uniform spatial grid, used for evaluating integrals of bound
!  orbitals and densities numerically. 
!
!  The implementation essentially follows:
!
!   A.D. Becke, J Chem Phys 88, 2547 (1988)
!
!  More sophisticated schemes do exist, especially for larger molecules.
!  However, the original Becke's scheme is simple to implement, efficient
!  enough for small systems, and is accurate enough.
!
!  The routines in this module are intended for "disposable" single-use
!  grids. If the same grid is used multiple times, consider adding 
!  caching of the integration points and weights - perhaps using an
!  external file.
!
!  It is safe and useful to invoke GridPointsBatch('Next') from a parallel
!  region. Each call will receive a unique set of points.
!
  module molecular_grid
    use accuracy
    use atoms
    use timer
    use lebedev
    use math
    !$ use OMP_LIB
    implicit none
    private
    public mol_grid
    public GridInitialize, GridDestroy, GridPointsBatch
    !
    integer(ik), parameter     :: step_order = 3 ! There does not seem to be a good reason to change this one
    !
    !  Radial grid type is used to cache integration points and weights on a
    !  [-1:1] grid. They'll be scaled as appropriate for each atom.
    !
    type radial_grid
      real(rk), allocatable :: rw(:,:)  ! Grid positions (1,:) and weights (2,:)
    end type radial_grid
    !
    !  We may maintain multiple grids; as the result, all grid parameters must be collected in mol_grid type.
    !
    type mol_grid
      private
      !
      !  Parameters defining the grid
      !
      integer(ik)                    :: natom         ! Number of integration centres
      real(rk), allocatable          :: xyz(:,:)      ! Coordinates of the centres
      real(rk), allocatable          :: rat(:)        ! Parameter used for scaling the radial grid - eq. 25
      integer(ik), allocatable       :: nrad(:)       ! Number of radial grid points per centre
      integer(ik), allocatable       :: nang(:)       ! Number of points per angular shell
      !
      !  Position of the next point batch within the grid. We always return an entire angular shell
      !
      integer(ik)                    :: iatom
      integer(ik)                    :: irad
      !$ integer(OMP_LOCK_KIND)      :: update_lock   ! Lock must be acquired before modifying iatom, irad
      !
      type(radial_grid), allocatable :: radial(:)     ! Radial grid cache. Grid of order N is kept
                                                      ! at index N of radial(). Some (most) elements
                                                      ! of the cache will be empty.
    end type mol_grid
    !
    contains
    !
    !  Externally visible interfaces
    !
    subroutine GridInitialize(grid,nrad,nang,xyz,types)
      type(mol_grid), intent(out)   :: grid       ! Grid to initialize
      integer(ik), intent(in)       :: nrad       ! Basic number of radial points; actual number of points 
                                                  ! may depend on this value and atom types
      integer(ik), intent(in)       :: nang       ! Number of angular points; can be 110, 302, or 770,
                                                  ! giving the corresponding Lebedev's grid
      real(rk), intent(in)          :: xyz(:,:)   ! Coordinates of centres, in Bohr
      character(len=*), intent(in)  :: types(:)   ! Types of centres
      !
      integer(ik) :: alloc
      integer(ik) :: natom, iat, jat
      integer(ik) :: irad, rad_max
      real(rk)    :: rat, rij
      !
      call TimerStart('MolGrid Initialize')
      !
      if (size(xyz,dim=1)/=3 .or. size(xyz,dim=2)/=size(types)) then
        stop 'molecular_grid%GridInitialize - bad input array sizes'
      end if
      !
      natom = size(xyz,dim=2)
      grid%natom = natom
      allocate (grid%xyz(3,natom),grid%rat(natom),grid%nrad(natom),grid%nang(natom),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating array for an ',i8,'-atom grid.')") alloc, natom
        stop 'molecular_grid%GridInitialize - out of memory (1)'
      end if
      !
      grid%xyz   = xyz
      grid%nang  = nang  ! Checking will occur later
      grid%nrad  = nrad  ! We'll modify later as appropriate
      grid%rat   = 1._rk ! We'll modify later as appropriate
      grid%iatom = 1_ik  ! Position the grid at the beginning
      grid%irad  = 1_ik
      !$ call omp_init_lock(grid%update_lock)
      !
      !  Fill atom sizes
      !
      scan_atom_types: do iat=1,natom
        rat = AtomCovalentR(types(iat))   ! Returns radius in Angstrom; we use Bohr here
        if (rat<=0) cycle scan_atom_types
        grid%rat(iat) = 0.5_rk * rat / abohr
      end do scan_atom_types
      !
      !  Check for collisions: this will cause problems with integrals!
      !
      collision_test: do iat=2,natom
        do jat=1,iat-1
          rij = sqrt(sum( (xyz(:,iat)-xyz(:,jat))**2 ))
          if (rij>=spacing(1e5_rk)) cycle
          write (out,"(/'molecular_grid%GridInitialize - collision of atomic centres is not allowed')")
          write (out,"(' Atom ',i5,' at ',3g25.16)") iat, xyz(:,iat)
          write (out,"(' Atom ',i5,' at ',3g25.16)") jat, xyz(:,jat)
          stop 'molecular_grid%GridInitialize - atom collision'
        end do
      end do collision_test
      !
      !  Fill radial grid cache - recalculating radial grids on each call may
      !  be quite expensive!
      !
      rad_max = maxval(grid%nrad)
      allocate (grid%radial(rad_max),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' allocating radial grid cache. rad_max = ',i8)") alloc, rad_max
        stop 'molecular_grid%GridInitialize - out of memory (2)'
      end if
      !
      fill_radial_cache: do iat=1,natom
        irad = grid%nrad(iat)
        if (allocated(grid%radial(irad)%rw)) cycle fill_radial_cache
        allocate (grid%radial(irad)%rw(2,irad),stat=alloc)
        if (alloc/=0) then
          write (out,"('Error ',i8,' allocating radial grid ',i8)") alloc, irad
          stop 'molecular_grid%GridInitialize - out of memory (3)'
        end if
        call MathGetQuadrature('Legendre',irad,grid%radial(irad)%rw(1,:),grid%radial(irad)%rw(2,:))
      end do fill_radial_cache
      call TimerStop('MolGrid Initialize')
    end subroutine GridInitialize
    !
    subroutine GridDestroy(grid)
      type(mol_grid), intent(inout) :: grid 
      !
      integer(ik) :: iord
      !
      call TimerStart('MolGrid Destroy')
      free_radial_cache: do iord=1,size(grid%radial)
        if (.not.allocated(grid%radial(iord)%rw)) cycle free_radial_cache
        deallocate (grid%radial(iord)%rw)
      end do free_radial_cache
      !
      deallocate (grid%xyz,grid%rat,grid%nrad,grid%nang)
      call TimerStop('MolGrid Destroy')
    end subroutine GridDestroy
    !
    subroutine GridPointsBatch(grid,action,xyzw,count,done)
      type(mol_grid), intent(inout)                  :: grid       ! Grid
      character(len=*), intent(in)                   :: action     ! What to do
!     real(rk), pointer, intent(inout), optional     :: xyzw(:,:)  ! Next batch of points
!    Declaration above causes Intel compiler 9.1 to barf at the invocation point of fill_grid_shell. Compiler bug?
      real(rk), pointer,                optional     :: xyzw(:,:)  ! Next batch of points
      integer(ik), intent(out), optional             :: count      ! Total number of batches
      logical, intent(out), optional                 :: done       ! Fill be set to .true. if out of points
      !
      integer(ik) :: iatom, irad ! Local copies of iatom and irad; needed to avoid holding
                                 ! a lock for too long.
      integer(ik) :: nang        ! Size of the angular batch
      integer(ik) :: alloc       
      !
      call TimerStart('MolGrid PointsBatch')
      batch_case: select case (action)
        case default
          write (out,"('molecular_grid%GridPointsBatch: unknown action ',a)") trim(action)
          stop 'molecular_grid%GridPointsBatch - Bad action'
        case ('Batches count')
          if (.not.present(count)) stop 'molecular_grid%GridPointsBatch - required argument is missing (1)'
          count = sum(grid%nrad)
        case ('Reset')
          grid%iatom = 1_rk
          grid%irad  = 1_rk
        case ('Next batch')
          if (.not.present(xyzw)) stop 'molecular_grid%GridPointsBatch - required argument is missing (2)'
          !$ call omp_set_lock(grid%update_lock)
          !$omp flush
          iatom = grid%iatom
          irad  = grid%irad
          if (iatom>grid%natom) then
            if (present(done)) then
              done = .true.
            else
              stop 'molecular_grid%GridPointsBatch - run out of points, but "done" is not present!'
            end if
            !$ call omp_unset_lock(grid%update_lock)
            !
            call TimerStop('MolGrid PointsBatch')
            return
          end if
          grid%irad = irad + 1
          if (grid%irad > grid%nrad(iatom)) then
            grid%irad  = 1
            grid%iatom = iatom + 1
          end if
          !$omp flush
          !$ call omp_unset_lock(grid%update_lock)
          !
          !  We are commited to a particular batch of points; the rest of this routine will not
          !  modify anything within the grid structure.
          !
          nang = grid%nang(iatom)
          if (associated(xyzw)) then
            if (size(xyzw,dim=2)/=nang) deallocate (xyzw)
          end if
          if (.not.associated(xyzw)) then
            allocate (xyzw(4,nang),stat=alloc)
            if (alloc/=0) then
              write (out,"('molecular_grid%GridPointsBatch: Error ',i8,' allocating batch buffer for ',i8,' points')") &
                     alloc, nang
              stop 'molecular_grid%GridPointsBatch - out of memory'
            end if
          end if
          !
          !  Done with all preliminaries, fill the grid points
          !
          call fill_grid_shell(grid,iatom,irad,xyzw(:,:))
      end select batch_case
      call TimerStop('MolGrid PointsBatch')
    end subroutine GridPointsBatch
    !
    !  Internal routines beyond this point
    !
    subroutine fill_grid_shell(grid,iatom,irad,xyzw)
      type(mol_grid), intent(in) :: grid      ! Everything we want to know about the grid
      integer(ik), intent(in)    :: iatom     ! Atom index; ignore data within the structure
      integer(ik), intent(in)    :: irad      ! Radial index; ditto
      real(rk), intent(out)      :: xyzw(:,:) ! Grid positions and weights
      !
      !  Step one: fill grid positions and unscaled weights.
      !
      integer(ik)       :: nrad          ! Number of radial points on this atom
      integer(ik)       :: nang          ! Number of angular points on this atom
      integer(ik)       :: iang
      real(rk)          :: rc(3)         ! Center of the atomic sphere
      real(rk)          :: rad           ! Radius of the current sphere
      real(rk)          :: w_rad         ! Integration weight associated with the radial point
      real(rk)          :: r_mid         ! Characteristic radius of the atom
      real(rk), pointer :: ang_xyzw(:,:) ! Our angular grid
      real(rk)          :: p_wgt         ! Becke partitioning weight
      !
      nrad  = grid%nrad(iatom)
      nang  = grid%nang(iatom)
      rc    = grid%xyz(:,iatom)
      r_mid = grid%rat(iatom)
      rad   = grid%radial(nrad)%rw(1,irad) ! Unscaled grid position, [-1:+1] range
      w_rad = grid%radial(nrad)%rw(2,irad)
      !
      !  Scale grid position and weight. We'll use Becke's scaling function:
      !    r_mid * (1+x)/(1-x)
      !  The weight gets scaled by the derivative of the transformation function,
      !  which is:
      !    r_mid * 2/(1-x)**2
      !
      w_rad = w_rad * r_mid * 2._rk / (1-rad)**2
      rad   = r_mid * (1+rad)/(1-rad)
      !
      !  Scale weight by the spherical volume element, since our angular grid 
      !  is normalized to unity
      !
      w_rad = w_rad * 4._rk * pi * rad**2
      !
      !  Choose angular grid
      !
      select case (nang)
        case default
          write (out,"('molecular_grid%fill_grid_shell: Angular grid ',i6,' is not recognized.')") nang
          stop 'molecular_grid%fill_grid_shell - bad angular grid'
        case (110) ; ang_xyzw => lebedev_gr17
        case (302) ; ang_xyzw => lebedev_gr29
        case (770) ; ang_xyzw => lebedev_gr47
      end select
      !
      !  Figure out coordinates and weights of atomic grid points in this shell
      !
      xyzw(1:3,:) = spread(rc,dim=2,ncopies=nang) + rad * ang_xyzw(1:3,:)
      xyzw(  4,:) = w_rad * ang_xyzw(4,:)
      !
      !  Magic sauce: Becke's partitioning weights
      !
      partitioning_weights: do iang=1,nang
        p_wgt        = partitioning_weight(grid,iatom,xyzw(1:3,iang))
        xyzw(4,iang) = p_wgt * xyzw(4,iang)
      end do partitioning_weights
    end subroutine fill_grid_shell
    !
    !  Calculate Becke's partitioning weight for a grid point
    !
    function partitioning_weight(grid,iatom,xyz) result(w)
      type(mol_grid), intent(in) :: grid
      integer(ik), intent(in)    :: iatom  ! Parent atom
      real(rk), intent(in)       :: xyz(:) ! Grid point
      real(rk)                   :: w      ! Partitioning weight
      !
      real(rk)    :: w_tot
      integer(ik) :: jatom
      !
      w     = raw_weight(grid,iatom,xyz)
      w_tot = w
      all_claims: do jatom=1,grid%natom
        if (jatom==iatom) cycle all_claims
        w_tot = w_tot + raw_weight(grid,jatom,xyz)
      end do all_claims
      w = w / w_tot
    end function partitioning_weight
    !
    !  Calculate Becke's primitive cell weight
    !
    function raw_weight(grid,iatom,xyz) result(w)
      type(mol_grid), intent(in) :: grid
      integer(ik), intent(in)    :: iatom  ! Parent atom
      real(rk), intent(in)       :: xyz(:) ! Grid point
      real(rk)                   :: w      ! Partitioning weight
      !
      integer(ik) :: jatom
      real(rk)    :: ri, rj   ! Distance from grid point to atoms I and J
      real(rk)    :: rij      ! Distance between atoms I and J
      real(rk)    :: muij     ! Elliptical distance difference coordinate
      real(rk)    :: nuij     ! "Adjusted" cell boundary
      real(rk)    :: aij      ! Adjustment parameter
      real(rk)    :: chi      ! Ratio of atom sizes
      real(rk)    :: uij      ! 
      real(rk)    :: sf
      !
      ri = sqrt(sum( (xyz-grid%xyz(:,iatom))**2 ))
      w = 1.0_rk
      cell_accumulate: do jatom=1,grid%natom
        if (jatom==iatom) cycle cell_accumulate
        rj   = sqrt(sum( (              xyz-grid%xyz(:,jatom))**2 ))
        rij  = sqrt(sum( (grid%xyz(:,iatom)-grid%xyz(:,jatom))**2 ))
        muij = (ri-rj)/rij
        !
        !  "Adjust" the cell boundary to account for heteronuclear bonding - Becke's appendix A
        !
        chi  = grid%rat(iatom)/grid%rat(jatom)
        uij  = (chi-1._rk)/(chi+1._rk)
        aij  = uij/(uij**2-1._rk)
        aij  = min(0.5_rk,max(-0.5_rk,aij))
        !
        nuij = muij + aij * (1._rk-muij**2)
        sf   = step_function(muij)
        w    = w * sf
      end do cell_accumulate
    end function raw_weight
    !
    function step_function(mu) result(sf)
      real(rk), intent(in)    :: mu
      real(rk)                :: sf
      !
      integer(ik) :: iord
      !
      sf = mu
      recurse_function: do iord=1,step_order
        sf = 1.5_rk * sf - 0.5_rk * sf**3
      end do recurse_function
      sf = 0.5_rk * (1.0_rk - sf)
    end function step_function
  end module molecular_grid
