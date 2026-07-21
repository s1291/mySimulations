! Updates History:
! Date: 09.09.2025
! Author: S.Ouchene
! Changes description: 
!     - Format the module correctly: spacing, indentations, etc.
!     - Replace idp2 with wp (working precision) that is used commonly in PALM
!
!=======================================================================
module nw_smearing_correction_module 
  !-----------------------------------------------------------------------
  !     Variables used in the new near wake tip correction
  !----------------------------------------------------------------------- 

      
  ! Global precision
  use kinds, ONLY: wp
  
  implicit none
  integer, PARAMETER ::  MyProcNum = 0   
  ! Switch activating the correction
  logical :: smearing_correction_active = .TRUE. !< To disable NWC if set to .FALSE.

  logical :: smearing_correction_needs_init = .TRUE.

  ! At the root and tip the control point lies on a trailed vortex, where the
  ! self-induction is zeroed. When .TRUE. the correction at those two nodes is
  ! set by linear extrapolation from the adjacent interior nodes instead.
  logical :: nwc_extrapolate_endpoint_correction = .TRUE.

  ! Switch invocking the recomputation of constants (large changes
  ! in the blade geometry from deformation could make it necessary)
  logical, dimension(:,:), allocatable :: recompute_smearing_constants

  ! Upper bound on the number of distinct airfoils per turbine.
  ! Used to size the namelist parameter airfoils_filenames in
  ! fastv8_coupler_mod and the matching dummy arguments downstream.
  integer, parameter :: n_max_airfoils = 100

  ! Airfoil polars and blade geometry, loaded once at init from input.f90.
  ! All shared across turbines (today's design has one airfoils_folder_path
  ! and one airfoils_filenames list for the whole farm).
  real(kind=wp), dimension(:),   allocatable :: chord_arr     ! (n_blade_elem)
  real(kind=wp), dimension(:),   allocatable :: twist_arr     ! (n_blade_elem)
  real(kind=wp), dimension(:),   allocatable :: r_cp_arr      ! (n_blade_elem) radial position of CPs
  real(kind=wp), dimension(:),   allocatable :: r_vtp_arr     ! (n_blade_elem + 1) vortex trailing point radii
  real(kind=wp), dimension(:,:), allocatable :: alpha_2d_arr  ! (n_polar_rows, n_blade_elem)
  real(kind=wp), dimension(:,:), allocatable :: cl_2d_arr     ! (n_polar_rows, n_blade_elem)
  real(kind=wp), dimension(:,:), allocatable :: cd_2d_arr     ! (n_polar_rows, n_blade_elem)

  ! Per-element classification used by smooth_trailed_vorticity to localise
  ! the Gamma-smoothing kernel to cyl/airfoil interfaces only. An element is
  ! marked TRUE when its polar table has MAXVAL(|Cl|) above a small cutoff,
  ! i.e. the section actually generates lift. Cylinder sections have Cl
  ! identically zero in their polar and end up FALSE.
  logical, dimension(:), allocatable :: is_lifting_section_arr  ! (n_blade_elem)

  ! Definition of setup: number of rotors, nr. of blades, nr. of 
  ! vortex trailing points (in-between blade sections), nr. of 
  ! sections along blade
  integer :: n_rotors = 0, n_blades = 0, n_v = 0, n_s = 0

  ! Index range of control points (CP): The missing induction is 
  ! only felt in the proximity of the shed vortices, allowing to 
  ! limit the number of CPs to loop over.
  integer, dimension(:,:,:,:), allocatable :: cp_loop_ind 

  ! Pi
  real(kind=wp) :: pi = 4.0_wp * atan( 1.0_wp )

  ! Convergence limit of induced velocity residual (default value)
  real(kind=wp) :: nw_max_resid = 1e-3_wp      !nw_max_resid=1d-6

  ! Constants of the Near-wake model
  real(kind=wp), dimension(4,5) :: nw_N, nw_P

  ! Variables from previous time step: For smooth transition between
  ! steps history needs to be stored
  real(kind=wp), dimension(:,:,:), allocatable:: nw_vn_history

  real(kind=wp), dimension(:,:,:), allocatable:: nw_vt_history

  real(kind=wp), dimension(:,:,:), allocatable:: nw_Gamma_history

  real(kind=wp), dimension(:,:,:), allocatable:: nw_W_history

  real(kind=wp), dimension(:,:,:,:), allocatable:: nw_dyn_comp_history 

  real(kind=wp), dimension(:,:,:,:), allocatable:: nw_X_history

  real(kind=wp), dimension(:,:,:,:), allocatable:: nw_Y_history

  ! Constants that need to be recomputed only after strong changes
  ! in geometry
  real(kind=wp), dimension(:,:,:,:), allocatable:: nw_phi, nw_phi_s, nw_dx, nw_dy, smear_fac

  real(kind=wp), dimension(:,:,:,:,:), allocatable:: nw_a

  ! Per-element diagnostics for wt-output NetCDF. Populated inside
  ! nw_smearing_correction at every call (entry: *_pre; after the
  ! convergence loop: *_post and induced *_nw_*). Snapshot only - no
  ! history is kept. Indexed (turbine, blade, element). Allocated by
  ! f8c_parin under the same gate as the nw_*_history arrays
  ! (use_near_wake_smearing_correction = .TRUE.) and read by
  ! wt_output_mod when writing a NetCDF slice.
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_vn_pre,    nwc_diag_vt_pre
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_vn_post,   nwc_diag_vt_post
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_phi_pre,   nwc_diag_phi_post
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_alpha_pre, nwc_diag_alpha_post
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_cl_pre,    nwc_diag_cl_post
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_cd_pre,    nwc_diag_cd_post
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_nw_vn,     nwc_diag_nw_vt
  real(kind=wp), dimension(:,:,:), allocatable :: nwc_diag_nw_W

end module
