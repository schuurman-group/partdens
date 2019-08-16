!
!  Routines for handling input and output files
!
module fileio
  use accuracy
  implicit none
  public
  !
  type input_data
    ! input file parameters
    character(100)      :: orb_fname     ! the orbital filename
    character(6)        :: vec_type      ! the $VEC type (natorb or tr1rdm)
    ! grid parameters
    character(5)        :: grid_type     ! the type of grid (becke or carte)
    integer(ik)         :: n_rad         ! the number of radial points on the Becke grid
    integer(ik)         :: n_ang         ! the number of angular points on the Becke grid
    real(rk)            :: covrad_max    ! the maximum covalent radius scale factor for the cartesian grid
    real(rk)            :: covrad_min    ! the minimum covalent radius scale factor for the cartesian grid
    real(rk)            :: covrad_const  ! a constant value added to the maximum covalent radius
    real(rk)            :: dx            ! the cartesian step size
    ! charge method parameters
    character(5)        :: weight_type   ! the weighting method (hirsh, becke, qtaim, voron)
    integer(ik)         :: max_charge    ! the maximum absolute charge for any atom
    ! atomic density parameters
    character(3)        :: atom_type     ! the choice of atomic reference (pop or pro)
    character(100)      :: atom_lib      ! the choice of atomic densities (slater or method.basis)
    character(100)      :: lib_path      ! the path to the library files (default: ./atomlib)
    character(3)        :: interp_type   ! the interpolation for atomic densities (exp or pol)
    integer(ik)         :: n_rad_atom    ! the number of radial points for the atom library
    integer(ik)         :: interp_ord    ! the interpolation order (power)
    integer(ik)         :: max_iter      ! the maximum number of iterations for promolecule methods
    real(rk)            :: chg_thresh    ! the threshold for charge convergence of promolecule methods
  end type input_data

  contains

!
!  Read a crude input file
!
subroutine read_input(inf, struct)
    integer(ik),intent(in)       :: inf
    type(input_data),intent(out) :: struct
    !
    character(20)  :: varname, var
    integer(ik)    :: ios

    ! default values
    struct%orb_fname    = "nos.dat"
    struct%vec_type     = "natorb"
    struct%n_rad        = 70
    struct%n_ang        = 110
    struct%covrad_max   = 1.6
    struct%covrad_min   = 1e-3
    struct%covrad_const = 0.6
    struct%dx           = 0.1
    struct%max_charge   = 3
    struct%atom_lib     = "mp2.aug-cc-pVTZ"
    struct%lib_path     = "/globalhome/rymac/Projects/PartialCharge/partdens/atomlib/"
    struct%interp_type  = "exp"
    struct%n_rad_atom   = 70
    struct%interp_ord   = 8
    struct%max_iter     = 100
    struct%chg_thresh   = 1e-4

    ! loop through each line of the input file
    ios = 0
    parse_input: do
        read(inf,*,iostat=ios) varname, var
        if (ios /= 0) then
            exit parse_input
        end if
        select case (varname)
            case ("orb_fname")
                struct%orb_fname = var
            case ("vec_type")
                struct%vec_type = var
            case ("grid_type")
                struct%grid_type = var
            case ("atom_type")
                struct%atom_type = var
            case ("atom_lib")
                struct%atom_lib = var
            case ("interp_type")
                struct%interp_type = var
            case ("weight_type")
                struct%weight_type = var
            case ("n_rad")
                read(var,*,iostat=ios) struct%n_rad
            case ("n_ang")
                read(var,*,iostat=ios) struct%n_ang
            case ("n_rad_atom")
                read(var,*,iostat=ios) struct%n_rad_atom
            case ("interp_ord")
                read(var,*,iostat=ios) struct%interp_ord
            case ("max_charge")
                read(var,*,iostat=ios) struct%max_charge
            case ("max_iter")
                read(var,*,iostat=ios) struct%max_iter
            case ("chg_thresh")
                read(var,*,iostat=ios) struct%chg_thresh
            case ("covrad_max")
                read(var,*,iostat=ios) struct%covrad_max
            case ("covrad_min")
                read(var,*,iostat=ios) struct%covrad_min
            case ("covrad_const")
                read(var,*,iostat=ios) struct%covrad_const
            case ("dx")
                read(var,*,iostat=ios) struct%dx
        end select
        if (ios /= 0) then
            stop "read_input - incorrect variable type"
        end if
    end do parse_input

end subroutine read_input

!
!  Initialize the output file for gridchg.f90
!
subroutine init_gridchg_output(ouf, struct)
    integer(ik),intent(in)      :: ouf
    type(input_data),intent(in) :: struct

    write(ouf,'("+--------------------------------------------------+")')
    write(ouf,'("|                                                  |")')
    write(ouf,'("|                     GridChg                      |")')
    write(ouf,'("|                                                  |")')
    write(ouf,'("|         Grid-based partial atomic charges        |")')
    write(ouf,'("|       RJ MacDonell, MS Schuurman 2018-2019       |")')
    write(ouf,'("+--------------------------------------------------+")')
    write(ouf,'("")')
    write(ouf,'("")')
    write(ouf,'("Input summary:")')
    write(ouf,'("    --------- Orbital input ----------")')
    write(ouf,'("    orb_fname      =   ",a15)') trim(struct%orb_fname)
    write(ouf,'("    vec_type       =   ",a15)') trim(struct%vec_type)
    write(ouf,'("")')
    write(ouf,'("    -------------- Grid --------------")')
    write(ouf,'("    grid_type      =   ",a15)') struct%grid_type
    if (trim(struct%grid_type) == "becke") then
        write(ouf,'("    n_rad          =   ",i15)') struct%n_rad
        write(ouf,'("    n_ang          =   ",i15)') struct%n_ang
    else if (trim(struct%grid_type) == "carte") then
        write(ouf,'("    covrad_max     =   ",e15.3)') struct%covrad_max
        write(ouf,'("    conrad_min     =   ",e15.3)') struct%covrad_min
        write(ouf,'("    covrad_const   =   ",e15.3)') struct%covrad_const
        write(ouf,'("    dx             =   ",e15.3)') struct%dx
    end if
    write(ouf,'("")')
    write(ouf,'("    ------------- Charge -------------")')
    write(ouf,'("    weight_type    =   ",a15)') trim(struct%weight_type)
    write(ouf,'("    max_charge     =   ",i15)') struct%max_charge
    write(ouf,'("")')
    write(ouf,'("    --------- Atomic density ---------")')
    write(ouf,'("    atom_type      =   ",a15)') trim(struct%atom_type)
    write(ouf,'("    atom_library   =   ",a15)') trim(struct%atom_lib)
    write(ouf,'("    lib_path       =   ",a)') trim(struct%lib_path)
    write(ouf,'("    n_rad_atom     =   ",i15)') struct%n_rad_atom
    write(ouf,'("    max_iter       =   ",i15)') struct%max_iter
    write(ouf,'("    chg_thresh     =   ",e15.3)') struct%chg_thresh
    write(ouf,'("    interp_type    =   ",a15)') trim(struct%interp_type)
    write(ouf,'("    interp_ord     =   ",i15)') struct%interp_ord
    write(ouf,'("")')

end subroutine init_gridchg_output

end module fileio
