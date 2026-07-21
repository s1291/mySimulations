!> @file fastv8_coupler_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2023 Carl von Ossietzky Universität Oldenburg
! Copyright 2022-2023 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Coupler to multiple FASTv8 wind turbine Simulations using TCP/IP sockets

! Versions and Changelog:
! -----------------------
! TODO: should we keep this here or to a sepate 'changelog' file?
!
! Version 1.0.0:
! Author: Sonja Steinbrück.
! Changes:
!    - Original version of the NWC code.
!
! Version: 1.1.0
! Author: Samir Ouchene
! Changes:
!  - Add 'use_near_wake_smearing_correction' (defaults to .FALSE.) parameter to the fastv8_coupler_parameters
!     namelist and wrap the related variables (e.g., allocations) within IF statements.
!  - Use the _wp suffix for real literals instead of idp2 from NWC code.
!  - Add the missing _wp to 0.00000001_wp used in the expvalue comparison.
!  - Initialize the nw_* arrays with 0.0_wp instead of 0.0d0.
!  - Fix the typo in the initialization of the integer n_dt_fast by changing 0.0_iwp to 0_iwp.
!  - Refactor nw_smearing_correction_module and al_nw_smearing_correction  to convert real literals
!    from 'd' notation to PALM-style _wp literals.
!  - fix bugs in the handling of polars (needed for NWC) in input.f90 file. This file is a temporary solution
!    and will be removed in future versions, where the polar data will be sent directly by OpenFAST instead of
!    being copied manually into files and read through input.f90.
!  - Partially translate the comments from German to English.
!
! Version: 1.2.0
! Author: Samir Ouchene
! Changes:
!  - .....
!  - TODO: For the new code replacing input.f90 make sure the polars are sent for each turbine not only turbine 1.
!  - TODO: Fix the bug in handling the indices of fboxcorners with MPI.
!  - TODO: In the C code from OpenFAST side, the number of blades is incorrectly hardcoded to 3 when allocating arrays,
!    such sa VelMatrix, and when exchanging data between PALM/OpenFAST. As a result, the coupling will not work correctly
!    for two-bladed wind turbines.
!  - NOTE: The blade pitch is now sent from OpenFAST (ED%y%BlPitch, one value per blade) every coupling step
!    and used by the NWC, so pitch-controlled runs work. See blade_pitch and palm_blade_pitch, and the
!    radians-to-degrees conversion in velocity_correction.
!  - TODO: Fix fast_n_blades and fast_n_blade_elem broadcast with wrong MPI datatype. It should be MPI_INTEGER instead of MPI_REAL
!  - TODO: Apparently, the rotor shaft vector (xs_x, xs_y, xs_z) is sent incorrectly: the y and z components are swapped and
!    one component is negated.This bug is 'dormant' when the rotor shaft tilt angle is zero, because the rotor is aligned along
!    the x streamwise direction, so the vector is (1, 0, 0). However, for large tilt values, the sent shaft vector is incorrect.
!  - TODO: x_prime_sec/y_prime_sec/z_prime_sec allocated, filled and deallocated every PALM timestep, but never read or used.
!  - TODO: notice that if palm_dt_3d FALSE (which is not the case during runtime, i.e. set to TRUE) then we advance with 120deg
!    This is confusing. Why do we need that? also this implicitly assumes that the turbine has 3 blades which won't work for
!    2-bladed turbines.
!  - NOTE: In the ASM, it's assumed that the rotor undergoes pure rigid rotation about the shaft axis. But of course if you enable
!    other degrees of freedom in OpenFAST (Tower fore-aft/side-side oscillations between l= 1 and l=N, Blade flapwise/edgewise
!    deflections changes within the PALM step, pitch-angle changes within the PALM step, yaw motion within PALM step, platform
!    surge/sway/heave/roll/ptich/yaw for floating turbines) then all of those are frozen at their l=1 (first FAST sub-step) values
!    and ignored for smearing purposes until the next PALM step. However, these DOF depend on PALM grid spacing and their orders
!    of magnitude. One should assess the validity of ASM, or switch to semi-ALM mode (shrink dt_PALM  towards dt_FAST until
!    dt_PALM/dt_FAST ~ 1 ; assuming these timesteps are sufficient to resolve them)
!--------------------------------------------------------------------------------------------------!
 MODULE fastv8_coupler_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               zu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  pi

    USE control_parameters,                                                                        &
        ONLY:  current_timestep_number,                                                            &
               debug_output,                                                                       &
               dt_3d,                                                                              &
               dz,                                                                                 &
               end_time,                                                                           &
               initializing_actions,                                                               &
               interpolate_to_grid_center,                                                         &
               message_string,                                                                     &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               simulated_time,                                                                     &
               terminate_run,                                                                      &
               simulated_time_at_begin

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               dx,                                                                                 &
               ddy,                                                                                &
               dy


    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  comm2d,                                                                             &
               ierr,                                                                               &
               myid,                                                                               &
               npex,                                                                               &
               npey

! SoS:
        
    USE wt_output_mod,                                                                             &
        ONLY:  write_wt_data, wt_output_meta_t, wt_output_warn_if_slow

    USE nw_smearing_correction_module, only: nw_X_history,nw_Y_history,                &
                                       nw_vn_history, nw_vt_history, nw_Gamma_history, &
                                       nw_dyn_comp_history, nw_W_history, nw_phi,      &
                                       nw_phi_s, nw_dx, nw_dy, nw_a, cp_loop_ind,      &
                                       smear_fac, recompute_smearing_constants,        &
                                       chord_arr, twist_arr, r_cp_arr, r_vtp_arr,      &
                                       alpha_2d_arr, cl_2d_arr, cd_2d_arr,             &
                                       is_lifting_section_arr,                         &
                                       nwc_extrapolate_endpoint_correction,            &
                                       nwc_diag_vn_pre,    nwc_diag_vt_pre,            &
                                       nwc_diag_vn_post,   nwc_diag_vt_post,           &
                                       nwc_diag_phi_pre,   nwc_diag_phi_post,          &
                                       nwc_diag_alpha_pre, nwc_diag_alpha_post,        &
                                       nwc_diag_cl_pre,    nwc_diag_cl_post,           &
                                       nwc_diag_cd_pre,    nwc_diag_cd_post,           &
                                       nwc_diag_nw_vn,     nwc_diag_nw_vt,             &
                                       nwc_diag_nw_W
!s. Ouchene
    USE ISO_C_BINDING, ONLY: C_INT
! SoS                                   

    IMPLICIT NONE


    LOGICAL ::  fastv8_coupler_enabled = .FALSE.  !<
    LOGICAL ::  first_ts_with_wt = .TRUE.         !<
    LOGICAL ::  first_restart_with_wt = .FALSE.    !< SoS: for checking in case of a start of turbine only in restart not at 0s
!-- S. Ouchene: Flag to enable/disable near-wake smearing corrections (Meyer-Forsting model).
!-- Controlled via the fastv8_coupler_parameters namelist. Defaults to .FALSE. (corrections OFF).
!-- When OFF, all associated arrays are not allocated and no extra MPI communication is performed.
    LOGICAL ::  use_near_wake_smearing_correction = .FALSE. !< namelist parameter to enable/disable near-wake velocity corrections
    
!-- Number of (alpha, Cl, Cd) rows in each per-node airfoil polar table. This is
!-- no longer a namelist parameter: the polars now come from OpenFAST (AD15) and
!-- are merged onto a common alpha grid, so the row count is computed by
!-- f8c_build_nwc_arrays (the size of that union grid) at the first coupling step.
    INTEGER(iwp) :: n_airfoil_polar_rows = 0_iwp

!-- S. Ouchene (May 2026): namelist parameter that controls an optional
!-- spanwise Gaussian smoothing of the bound circulation Gamma inside the
!-- near-wake correction. Default 0.0 disables the kernel completely, so
!-- the wake model behaves bit-for-bit identically to earlier versions.
!-- A positive value gives the smoothing length L (in metres of span);
!-- a sensible starting point for the NREL 5 MW is ~0.5 - 1.0 m (about one
!-- chord at the cylinder/airfoil transition). The value is passed down
!-- through velocity_correction -> nw_smearing_correction to where the
!-- kernel is applied (no module-variable state is used).
    REAL(wp) :: nw_trailed_vorticity_smoothing_length = 0.0_wp

!-- S. Ouchene (May 2026): namelist parameters that control the wind-turbine NetCDF output
!-- (see USER_CODE/wt_output_mod.f90 and docs/wt-netcdf-output.md). Default
!-- output_interval = 0 disables the output entirely. A positive value writes every N-th step
!-- in PALM-step units by default; setting output_use_fast_timestep = .TRUE. switches the
!-- counter to FAST-sub-step units instead. Output is rank-0 only and lands in PALM's
!-- per-run working directory.
    INTEGER(iwp) :: output_interval          = 0_iwp
    LOGICAL      :: output_use_fast_timestep = .FALSE.

!-- S. Ouchene (May 2026): control whether the wind-turbine NetCDF file is fsync'd after every
!-- written time slice. Default .TRUE. - crash-safe but slower (every put_var triggers an HDF5
!-- metadata flush, which can be expensive on parallel filesystems). Set to .FALSE. for
!-- performance-sensitive runs: writes are still flushed at file close (wt_output_finalize) and
!-- buffered to OS cache in between, so at most the tail of the most recent buffer is lost on a
!-- mid-run crash.
    LOGICAL      :: output_sync_per_write    = .TRUE.

!-- S. Ouchene: declare the magic numbers for commclient function to use meaningful integer flags instead of 1, 2, 3, or 4
!-- Notice that the convention is specified in fastv8_server.c:commclient function
    INTEGER(C_INT), PARAMETER :: COM_TARGET_OPENFAST   = 0_C_INT
    INTEGER(C_INT), PARAMETER :: COM_INIT_CONNECTION   = 1_C_INT
    INTEGER(C_INT), PARAMETER :: COM_PALM_IS_READY     = 2_C_INT
    INTEGER(C_INT), PARAMETER :: COM_SEND_VELOCITIES   = 3_C_INT
    INTEGER(C_INT), PARAMETER :: COM_RESUME_SIMULATION = 4_C_INT
    
!-- Number of turbines
    INTEGER(iwp), PARAMETER ::  nturbines_max = 300_iwp  !< maximum number of turbines allowed (can be increased if required)
    INTEGER(iwp) ::  nturbines = 0_iwp                   !< namelist parameter
!
!-- Number of time steps in FAST
    INTEGER(iwp) ::  n_dt_fast = 0_iwp  !<

!    
!
!-- Time step size in FAST and PALM derived from the FAST dt
    REAL(wp) ::  dt_fast       !<
    REAL(wp) ::  dt_palm       !<
    REAL(wp) ::  current_time_fast_min  !<
    REAL(wp) ::  palm_time_at_begin  !< SoS: in case the turbine is not started with the first PALM run, which starts at 0s
    REAL(wp) ::  target_time = 0.0_wp  !< s.ouchene: temporary variable to hold the targer time at each palm time step

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  current_time_fast  !<
!
!-- Number of blades (1 value per wind turbine, value is received from FAST)
    INTEGER(iwp), DIMENSION(1:nturbines_max) ::  fast_n_blades = 0_iwp  !<

!-- The ratio of the PALM time step to the OpenFAST time step
    REAL(wp) :: dt_ratio = 0.0_wp   

!-- Trignometric terms used in the rotation matrix
    REAL(wp) :: cos_alpha             = 0.0_wp
    REAL(wp) :: sin_alpha             = 0.0_wp
    REAL(wp) :: one_minus_cos_alpha   = 0.0_wp

!-- Components of the shaft rotation axis used in the rotation matrix
    REAL(wp) :: shaft_x               = 0.0_wp
    REAL(wp) :: shaft_y               = 0.0_wp
    REAL(wp) :: shaft_z               = 0.0_wp

!-- Maximum allowed number of blades per turbine
    INTEGER(iwp) ::  fast_n_blades_max = 0_iwp  !< namelist parameter

!
!-- Number of blade elements per blade in FAST (1 value per wind turbine, value is received from FAST)
    INTEGER(iwp), DIMENSION(1:nturbines_max) ::  fast_n_blade_elem = 0_iwp  !<

!
!-- Maximum allowed number of blade elements per blade
    INTEGER(iwp) ::  fast_n_blade_elem_max = 0_iwp  !< namelist parameter

!
!-- Rotational velocity
    INTEGER(iwp) ::  n_sector_max  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  n_sector  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rotspeed          !< rotational velocity from FAST
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  blade_pitch     !< blade pitch from FAST (radians), indexed (turbine, blade)

!-- S. Ouchene: near-wake-correction (NWC) blade geometry and airfoil polars
!-- received from OpenFAST (AD15) at the first coupling step. They replace the
!-- disk read that input.f90 used to do. These hold the raw received tables at
!-- their native AeroDyn lengths; f8c_build_nwc_arrays turns them into the module
!-- arrays the NWC consumes (chord_arr, twist_arr, alpha_2d_arr, ...). They are
!-- shared across turbines, matching the single blade definition input.f90 read
!-- for the whole farm. Not stored in the restart file: OpenFAST re-sends them on
!-- the first exchange after a restart, before they are used.
    INTEGER(iwp) ::  nwc_recv_nnodes    = 0_iwp     !< blade-node count from OpenFAST
    INTEGER(iwp) ::  nwc_recv_nairfoils = 0_iwp     !< airfoil-table count from OpenFAST
    INTEGER(iwp) ::  nwc_recv_max_nrows = 0_iwp     !< longest airfoil table (column size)
    LOGICAL      ::  nwc_geom_received    = .FALSE.  !< set once the data has arrived
    LOGICAL      ::  nwc_arrays_populated = .FALSE.  !< set once the NWC arrays are built
    REAL(wp),     DIMENSION(:),   ALLOCATABLE ::  nwc_recv_r_cp       !< (nnodes) radial pos (m)
    REAL(wp),     DIMENSION(:),   ALLOCATABLE ::  nwc_recv_twist_deg  !< (nnodes) twist (deg)
    REAL(wp),     DIMENSION(:),   ALLOCATABLE ::  nwc_recv_chord      !< (nnodes) chord (m)
    INTEGER(iwp), DIMENSION(:),   ALLOCATABLE ::  nwc_recv_blafid     !< (nnodes) 1-based airfoil id
    INTEGER(iwp), DIMENSION(:),   ALLOCATABLE ::  nwc_recv_nrows      !< (nairfoils) rows per table
    REAL(wp),     DIMENSION(:,:), ALLOCATABLE ::  nwc_recv_alpha      !< (max_nrows, nairfoils) deg
    REAL(wp),     DIMENSION(:,:), ALLOCATABLE ::  nwc_recv_cl         !< (max_nrows, nairfoils)
    REAL(wp),     DIMENSION(:,:), ALLOCATABLE ::  nwc_recv_cd         !< (max_nrows, nairfoils)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sector_angle      !< angle (in radiance) of sector, the rotor coveres in a PALM timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sector_angle_deg  !< angle (in degree) of sector, the rotor coveres in a PALM timestep
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  alpha_rot         !< angle (in radiance) of sector, the rotor coveres in a FAST timestep

    REAL(wp), DIMENSION(:),     ALLOCATABLE ::  shaft_height_fast  !< shaft height from FAST
    REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::  shaft_coordinates  !< shaft coordinates from FAST
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  r_n                !< rotation matrx

!
!-- Radius of turbines [ToDo]
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_n_radius = 0.0_wp  !<

    CHARACTER (LEN=15)                             ::  fast_host_addr_default = '127.0.0.1'  !<
    CHARACTER (LEN=15), DIMENSION(1:nturbines_max) ::  fast_host_addr = REPEAT(' ', 15)      !<
    CHARACTER (LEN=5),  DIMENSION(1:nturbines_max) ::  fast_host_port = REPEAT(' ', 5)       !<

!
!-- Tower base reference position in FAST (three components)
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_x = 0.0_wp  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_y = 0.0_wp  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  fast_tower_ref_pos_z = 0.0_wp  !<

!
!-- Tower base reference position in PALM (three components, namelist parameter)
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_x = - HUGE(1.0_wp)  !< namelist parameter
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_y = - HUGE(1.0_wp)  !< namelist parameter
    REAL(wp), DIMENSION(1:nturbines_max) ::  palm_tower_ref_pos_z = 0.0_wp          !< namelist parameter

!
!-- Sum of tower base reference position in PALM and FAST
!-- Dimension has to be larger than / equal to nturbines
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_x  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_y  !<
    REAL(wp), DIMENSION(1:nturbines_max) ::  tower_ref_pos_z  !<

!
!-- Hub center position (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_hub_center_pos  !<

!
!-- New 081012:
!-- In the calculations in PALM information on the positions
!-- is required at two different points in time
!-- Therefore, a second array for the positions of the hub center
!-- has to be added
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_hub_center_pos_old  !<

!
!-- Spinner aerodynamic force (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
!    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  fast_spinner_force

!
!-- Hub center velocity (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  palm_hub_center_vel  !<

    !-- Position of rotor blade elements (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_position_blade  !<

!
!-- In the calculations in PALM information on the positions
!-- is required at two different points in time
!-- Therefore, a second array for the positions of the blades
!-- has to be added
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_position_blade_old  !<

!
!-- Aerodynamic forces on blade segments (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_force_blade      !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  fast_force_blade_old  !<

!
!-- Blade velocities (in the course of the simulation)
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions
!-- to be allocated nturbines, fast_n_blades, fast_n_blade_elem, 3
!-- Note: fast_n_blades has to be set to the maximum number of blades
!-- for all simulated wind turbines; the missing value has to be
!-- entered for turbines with less blades than the turbine with the
!-- maximum number of blades
!SoS: palm_vel_blade_b                        
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  palm_vel_blade, palm_vel_blade_b  !<

!
!-- Time at which the wind turbine is switched on in the
!-- simulation (by default it would be switched on at the beginning
!-- of the run, namelist parameter)
    REAL(wp) ::  time_turbine_on = 0.0_wp  !<

!
!-- Simulated time in FAST (about simulated_time-time_turbine_on)
    !REAL(wp)    :: fast_simulated_time

!
!-- 3D arrays of forces caused by wind turbines in the model domain
!-- (forces smeared with a regularization kernel)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_x  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_y  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  force_z  !<

!
!-- Regularization factor for force distribution (epsilon = reg_fac * dx)
    REAL(wp) ::  reg_fac = 2.0_wp  !< namelist parameter

!
!-- Force distribution cut-off factor
!-- Edge length of box around turbine in which forces are distributed (edge length = edgelen_fac * turbine radius)
!-- Requirements for the size of the box:
!-- Around every point of attack there should be the at least the distance of 3 * grid spacing * reg_fac
!-- for force distribution available
    REAL(wp) ::  fbox_fac = 1.5_wp  !< namelist parameter

!-- Force distribution boxes around turbines
!-- has to be allocated in f8c_init
!-- Sizes of the dimensions to be allocated nturbines, 3, 2
!-- number of turbines; x_min, x_max; y_min, y_max; z_min, z_max
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  fboxcorners  !<

!
!-- Field with values of the exponential function
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  expvalue  !<

!
!-- x-, y- and z-components of the turbine induced forces:
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  thrx  !< thrx_av
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tory  !< tory_av
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  torz  !< torz_av


!
! -- Parameters and Variables for the calculation of tower area and thrust applied by the tower
    REAL(wp) ::  thrust_tower_x  !<
    REAL(wp) ::  thrust_tower_y  !<

    REAL(wp), DIMENSION(1:nturbines_max) ::  dtow         = 0.0_wp  !< tower diameter [m]
    REAL(wp), DIMENSION(1:nturbines_max) ::  htow         = 0.0_wp  !< tower height [m]
    REAL(wp), DIMENSION(1:nturbines_max) ::  turb_C_d_tow = 1.2_wp  !< drag coefficient for tower

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tower_area_x  !<
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tower_area_y  !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tower_area_x_4d  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tower_area_y_4d  !<

    SAVE

    PRIVATE
!
!-- Public functions
    PUBLIC                                                                                         &
       f8c_actions,                                                                                &
       f8c_check_data_output,                                                                      &
       f8c_check_parameters,                                                                       &
       f8c_data_output_2d,                                                                         &
       f8c_data_output_3d,                                                                         &
       f8c_define_netcdf_grid,                                                                     &
       f8c_header,                                                                                 &
       f8c_init,                                                                                   &
       f8c_init_arrays,                                                                            &
       f8c_parin
!
!-- Public parameters, constants and initial values
   PUBLIC                                                                                          &
      current_time_fast,                                                                           &
      dt_fast,                                                                                     &
      fast_force_blade,                                                                            &
      fast_position_blade,                                                                         &
      fast_n_blades,                                                                               &
      fast_n_blade_elem,                                                                           &
      fast_n_radius,                                                                               &
      fast_hub_center_pos,                                                                         &
      nturbines,                                                                                   &
      palm_hub_center_vel,                                                                         &
      palm_vel_blade,                                                                              &
      rotspeed,                                                                                    &
      blade_pitch,                                                                                 &
      nwc_recv_nnodes,                                                                             &
      nwc_recv_nairfoils,                                                                          &
      nwc_recv_max_nrows,                                                                          &
      nwc_geom_received,                                                                           &
      nwc_recv_r_cp,                                                                               &
      nwc_recv_twist_deg,                                                                          &
      nwc_recv_chord,                                                                              &
      nwc_recv_blafid,                                                                             &
      nwc_recv_nrows,                                                                              &
      nwc_recv_alpha,                                                                              &
      nwc_recv_cl,                                                                                 &
      nwc_recv_cd,                                                                                 &
      shaft_height_fast,                                                                           &
      fastv8_coupler_enabled,                                                                      &
      shaft_coordinates,                                                                           &
      palm_time_at_begin,                                                                          &        
      first_restart_with_wt


    INTERFACE f8c_parin
       MODULE PROCEDURE f8c_parin
    END INTERFACE f8c_parin

    INTERFACE f8c_check_parameters
       MODULE PROCEDURE f8c_check_parameters
    END INTERFACE f8c_check_parameters

    INTERFACE f8c_check_data_output
       MODULE PROCEDURE f8c_check_data_output
    END INTERFACE f8c_check_data_output

    INTERFACE f8c_define_netcdf_grid
       MODULE PROCEDURE f8c_define_netcdf_grid
    END INTERFACE f8c_define_netcdf_grid

    INTERFACE f8c_init
       MODULE PROCEDURE f8c_init
    END INTERFACE f8c_init

    INTERFACE f8c_init_arrays
       MODULE PROCEDURE f8c_init_arrays
    END INTERFACE f8c_init_arrays

    INTERFACE f8c_header
       MODULE PROCEDURE f8c_header
    END INTERFACE f8c_header

    INTERFACE f8c_actions
       MODULE PROCEDURE f8c_actions
       MODULE PROCEDURE f8c_actions_ij
    END INTERFACE f8c_actions

    INTERFACE f8c_data_output_2d
       MODULE PROCEDURE f8c_data_output_2d
    END INTERFACE f8c_data_output_2d

    INTERFACE f8c_data_output_3d
       MODULE PROCEDURE f8c_data_output_3d
    END INTERFACE f8c_data_output_3d


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &fastv8_coupler_parameters for user module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_parin

    CHARACTER (LEN=80) ::  line  !< string containing the last line read from namelist file

    INTEGER(iwp) ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /fastv8_coupler_parameters/                                                           &
       switch_off_module,                                                                          &
       nturbines,                                                                                  &
       time_turbine_on,                                                                            &
       fast_n_blades_max,                                                                          &
       fast_n_blade_elem_max,                                                                      &
       fast_host_addr_default,                                                                     &
       fast_host_addr,                                                                             &
       fast_host_port,                                                                             &
       palm_tower_ref_pos_x,                                                                       &
       palm_tower_ref_pos_y,                                                                       &
       palm_tower_ref_pos_z,                                                                       &
       reg_fac,                                                                                    &
       fbox_fac,                                                                                   &
       dtow,                                                                                       &
       htow,                                                                                       &
       turb_C_d_tow,                                                                               &
       first_restart_with_wt,                                                                      &
       use_near_wake_smearing_correction,                                                                  &
       nw_trailed_vorticity_smoothing_length,                                                      &
       nwc_extrapolate_endpoint_correction,                                                        &
       output_interval,                                                                            &
       output_use_fast_timestep,                                                                   &
       output_sync_per_write


!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND ( 11 )
    READ( 11, fastv8_coupler_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status.
    IF ( io_status == 0 )  THEN
!
!--    fastv8_coupler_parameters namelist was found and read correctly. Set flag that
!--    fastv8_coupler_mod is switched on.
       IF ( .NOT. switch_off_module )  fastv8_coupler_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    User namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'fastv8_coupler_parameters', line )

    ENDIF
    
        
!-- S. Ouchene: Allocate and initialise near-wake correction arrays only when corrections are enabled.
!-- These arrays store the vortex wake history and smearing parameters used by the Meyer-Forsting model.
!-- Skipping them when use_near_wake_smearing_correction = .FALSE. avoids unnecessary memory consumption.
    IF ( use_near_wake_smearing_correction )  THEN
       allocate(nw_X_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_Y_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_vn_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max))
       allocate(nw_vt_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max))
       allocate(nw_Gamma_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max))
       allocate(nw_dyn_comp_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max,3))
       allocate(nw_W_history(nturbines,fast_n_blades_max,fast_n_blade_elem_max))
       allocate(nw_phi(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_phi_s(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_dx(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_dy(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(nw_a(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max,4))
       allocate(cp_loop_ind(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,2))
       allocate(smear_fac(nturbines,fast_n_blades_max,fast_n_blade_elem_max+1,fast_n_blade_elem_max))
       allocate(recompute_smearing_constants(nturbines,fast_n_blades_max))

       IF ( (TRIM( initializing_actions ) /= 'read_restart_data' ) .OR. first_restart_with_wt)  THEN
          nw_X_history = 0.0_wp
          nw_Y_history = 0.0_wp
          nw_vn_history = 0.0_wp
          nw_vt_history = 0.0_wp
          nw_Gamma_history = 0.0_wp
          nw_dyn_comp_history = 0.0_wp
          nw_W_history =0.0_wp
       ENDIF

       nw_phi = 0.0_wp
       nw_phi_s = 0.0_wp
       nw_dx = 0.0_wp
       nw_dy = 0.0_wp
       nw_a = 0.0_wp
       ! Initialize with first and last index of control point
       cp_loop_ind(:,:,:,1) = 1
       cp_loop_ind(:,:,:,2) = fast_n_blade_elem_max
       smear_fac=0.0_wp
       recompute_smearing_constants=.TRUE.

!--    Pre-/post-correction diagnostics. Snapshot only (no history), populated by
!--    nw_smearing_correction on every call. Read by wt_output_mod for NetCDF output.
       ALLOCATE( nwc_diag_vn_pre   (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_vt_pre   (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_vn_post  (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_vt_post  (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_phi_pre  (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_phi_post (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_alpha_pre(nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_alpha_post(nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_cl_pre   (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_cl_post  (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_cd_pre   (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_cd_post  (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_nw_vn    (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_nw_vt    (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )
       ALLOCATE( nwc_diag_nw_W     (nturbines,fast_n_blades_max,fast_n_blade_elem_max) )

       nwc_diag_vn_pre     = 0.0_wp
       nwc_diag_vt_pre     = 0.0_wp
       nwc_diag_vn_post    = 0.0_wp
       nwc_diag_vt_post    = 0.0_wp
       nwc_diag_phi_pre    = 0.0_wp
       nwc_diag_phi_post   = 0.0_wp
       nwc_diag_alpha_pre  = 0.0_wp
       nwc_diag_alpha_post = 0.0_wp
       nwc_diag_cl_pre     = 0.0_wp
       nwc_diag_cl_post    = 0.0_wp
       nwc_diag_cd_pre     = 0.0_wp
       nwc_diag_cd_post    = 0.0_wp
       nwc_diag_nw_vn      = 0.0_wp
       nwc_diag_nw_vt      = 0.0_wp
       nwc_diag_nw_W       = 0.0_wp
    ENDIF
    
    ALLOCATE( rotspeed(1:nturbines) )
    ALLOCATE( shaft_coordinates(1:nturbines,1:3) )
    ALLOCATE( current_time_fast(1:nturbines) )
    current_time_fast = 0.0_wp

 END SUBROUTINE f8c_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check &fastv8_coupler_parameters control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_check_parameters

    INTEGER(iwp) ::  i  !< loop variable

#ifndef __fastv8
    IF ( fastv8_coupler_enabled )  THEN
       WRITE( message_string, * )  'Using FASTv8 Coupler requires cpp macro "__fastv8" during build'
       CALL message( 'f8c_check_parameters', 'F8C0001', 1, 2, 0, 6, 0 )
    ENDIF
#endif

    IF ( nturbines <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"nturbines" must set to a value greater than 0 ',                                       &
          '(currently it is set to ', nturbines, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0002', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( nturbines > nturbines_max )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter "nturbines" ',    &
          'must be set to a value less or equal to ', nturbines_max,                               &
          '(currently it is set to ', nturbines, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0003', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( fast_n_blades_max <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_n_blades_max" must set to a value greater than 0 ',                               &
          '(currently it is set to ', fast_n_blades_max, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0004', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( fast_n_blade_elem_max <= 0_iwp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_n_blade_elem_max" must set to a value greater than 0 ',                           &
          '(currently it is set to ', fast_n_blade_elem_max, ' )'
       CALL message( 'f8c_check_parameters', 'F8C0005', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM(fast_host_addr_default) == '' )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"fast_host_addr_default" must be set to a valid IP adress e.g 127.0.0.1',               &
          '(currently it is not set)'
       CALL message( 'f8c_check_parameters', 'F8C0006', 1, 2, 0, 6, 0 )
    ENDIF

    DO i = 1, nturbines
       IF ( TRIM(fast_host_port(i)) == '' )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"fast_host_port(', i, ')" must be set to a valid port number ',                      &
             '(currently it is not set)'
          CALL message( 'f8c_check_parameters', 'F8C0007', 1, 2, 0, 6, 0 )
       ENDIF
    END DO

    DO i = 1, nturbines
       IF ( palm_tower_ref_pos_x(i) < 0.0_wp  .OR.  palm_tower_ref_pos_x(i) > (nx+1) * dx )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_x(', i, ')" must be within the model domain. ',                  &
             '(currently it is set to ', palm_tower_ref_pos_x(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0008', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( palm_tower_ref_pos_y(i) < 0.0_wp  .OR.  palm_tower_ref_pos_y(i) > (ny+1) * dy )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_y(', i, ')" must be within the model domain.',                   &
             '(currently it is set to ', palm_tower_ref_pos_y(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0009', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( palm_tower_ref_pos_z(i) < 0.0_wp  .OR.  palm_tower_ref_pos_z(i) > zu(nzt) )  THEN
          WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',             &
             '"palm_tower_ref_pos_z(', i, ')" must be within the model domain.',                   &
             '(currently it is set to ', palm_tower_ref_pos_z(i), ' )'
          CALL message( 'f8c_check_parameters', 'F8C0010', 1, 2, 0, 6, 0 )
       ENDIF
    END DO

!-- The near-wake correction no longer needs airfoils_folder_path,
!-- airfoils_filenames or n_airfoil_polar_rows: the blade geometry and polars are
!-- received from OpenFAST (AD15) at the first coupling step, so there is nothing
!-- to validate here. f8c_build_nwc_arrays checks the received data instead.

!-- S. Ouchene: optional spanwise smoothing of the bound circulation Gamma
!-- inside the near-wake correction. A 3-point Gaussian kernel of width L
!-- (metres of span) is applied to Gamma(k) before dGamma is rebuilt. Default
!-- 0.0 means the kernel is not applied at all, preserving the original
!-- wake-model behaviour. Only non-negative values make physical sense.
    IF ( nw_trailed_vorticity_smoothing_length < 0.0_wp )  THEN
       WRITE( message_string, * )  'fastv8_coupler_parameters namelist parameter ',                &
          '"nw_trailed_vorticity_smoothing_length" must be >= 0.0 ',                               &
          '(currently it is set to ', nw_trailed_vorticity_smoothing_length, ' ). ',               &
          '0.0 disables the kernel; positive values give the spanwise smoothing length in metres.'
       CALL message( 'f8c_check_parameters', 'F8C0024', 1, 2, 0, 6, 0 )
    ENDIF

!    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
!      WRITE( message_string, * )  'Currently initializing_actions = "read_restart_data" is ',      &
!          'not implemented while using the FASTv8 Coupler'
!      CALL message( 'f8c_check_parameters', 'F8C0011', 1, 2, 0, 6, 0 )
!   ENDIF

!-- One-shot startup warning when the wt-output is enabled with the per-write fsync. No-op
!-- in the common case (output disabled, or sync already turned off). Prints once on rank 0.
    CALL wt_output_warn_if_slow( output_interval, output_sync_per_write )

 END SUBROUTINE f8c_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined output quantities. For those variables not recognized by the user,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_check_data_output( variable, unit )

    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<


    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx' )
          unit = 'm/s2'

       CASE ( 'tory' )
          unit = 'm/s2'

       CASE ( 'torz' )
          unit = 'm/s2'

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE f8c_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize user-defined arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_init_arrays

!
!-- Allocation of the arrays required for the output of additional 2D and 3D data.
    ALLOCATE( thrx(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tory(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( torz(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Allocation of position, force and velocity arrays that are
!-- required for the exchange of information between PALM and FAST.
!-- Arrays with dimensions dependent on the number of turbines:
!-- current_time_fast: first dimension number of turbines
!    ALLOCATE( current_time_fast(1:nturbines) )
!
!-- n_sector: first dimension number of turbines
    ALLOCATE( n_sector(1:nturbines) )
!
!-- rotspeed: first dimension number of turbines
!    ALLOCATE( rotspeed(1:nturbines) )
!
!-- sector_angle: first dimension number of turbines
    ALLOCATE( sector_angle(1:nturbines) )
!
!-- sector_angle_deg: first dimension number of turbines
    ALLOCATE( sector_angle_deg(1:nturbines) )
!
!-- alpha_rot: first dimension number of turbines
    ALLOCATE( alpha_rot(1:nturbines) )
!
!-- shaft_height_fast: first dimension number of turbines
    ALLOCATE( shaft_height_fast(1:nturbines) )
!
!-- shaft_coordinates: first dimension number of turbines, second
!-- dimension three directions in space
!    ALLOCATE( shaft_coordinates(1:nturbines,1:3) )
!
!-- r_n: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( r_n(1:nturbines,1:3,1:3) )
!
!-- fast_hub_center_pos: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( fast_hub_center_pos(1:nturbines,1:3) )
!
!-- fast_hub_center_pos_old: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( fast_hub_center_pos_old(1:nturbines,1:3) )
!
!-- palm_hub_center_vel: first dimension number of turbines, second
!-- dimension three directions in space
    ALLOCATE( palm_hub_center_vel(1:nturbines,1:3) )
!
!-- Arrays with dimension dependent on the number of blades
!-- and the number of blade elements per blade
!
!-- fast_position_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_position_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- New: 081012
!-- fast_position_blade_old: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_position_blade_old(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- fast_force_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( fast_force_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
    ALLOCATE( fast_force_blade_old(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!
!-- palm_vel_blade: first dimension number of turbines,
!-- second dimension number of blades, third dimension number
!-- of blade elements, fourth dimension three directions in space
    ALLOCATE( palm_vel_blade(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )
!-- blade_pitch holds the latest blade pitch from FAST (radians) per turbine and
!-- blade. It is filled every coupling step by palm_blade_pitch and read by
!-- velocity_correction. It is not stored in the restart file because FAST
!-- resends it on the first exchange after a restart, before it is used.
    ALLOCATE( blade_pitch(1:nturbines,1:fast_n_blades_max) )
    blade_pitch = 0.0_wp
!-- S. Ouchene: palm_vel_blade_b is the output buffer of velocity_correction(). It holds the
!-- corrected blade velocities before they are copied back into palm_vel_blade. Only needed
!-- when near-wake corrections are active.
    IF ( use_near_wake_smearing_correction )  &
       ALLOCATE( palm_vel_blade_b(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3) )

!-- S. Ouchene: the near-wake-correction blade geometry and airfoil polars used
!-- to be read from disk here by input.f90. They now come from OpenFAST (AD15)
!-- over the TCP coupling. That connection is not open yet at this point in init
!-- (the geometry arrives during the first f8c_actions call), and the polar grid
!-- size is not known until the data is received, so the NWC module arrays
!-- (chord_arr, twist_arr, r_cp_arr, r_vtp_arr, alpha_2d_arr, cl_2d_arr,
!-- cd_2d_arr, is_lifting_section_arr) are allocated and filled lazily by
!-- f8c_build_nwc_arrays at the first coupling step, which also calls
!-- init_nw_smearing_correction. Nothing to do here.
!
!-- Allocation of the arrays of forces (one in x-, y- and z-direction)
    ALLOCATE( force_x(nzb:nzt,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( force_y(nzb:nzt,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( force_z(nzb:nzt,nysg:nyng,nxlg:nxrg) )

!-- fboxcorners: first dimension number of turbines,
!-- second dimension three direction in space, third dimension min and max value
    ALLOCATE( fboxcorners(1:nturbines,1:3,1:2) )

!
!-- Allocate the field with values of the exponential function
    ALLOCATE( expvalue(0:10000000) )

!
!-- First of all the arrays that contain the information on the overlapping
!-- areas for each LES grid box have to be allocated; afterwards they are set to
!-- a basic value of 0 (no overlap between the rotor area and the LES grid box):
    ALLOCATE( tower_area_x(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tower_area_y(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

!
!-- Additionally 4D-arrays with the number of the wind turbines as fourth
!-- dimension are required:
    ALLOCATE( tower_area_x_4d(1:nturbines,nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tower_area_y_4d(1:nturbines,nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

 END SUBROUTINE f8c_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> S. Ouchene: build the near-wake-correction (NWC) blade geometry and airfoil
!> polar arrays from the data received from OpenFAST (AD15). Replaces the disk
!> read that input.f90 used to do in f8c_init_arrays.
!>
!> The airfoil tables arrive at their native AeroDyn lengths, which differ per
!> airfoil. They are merged onto a common alpha grid built as the sorted, de-
!> duplicated union of every airfoil's alpha breakpoints, and each airfoil is
!> linearly interpolated onto that grid. Because the union contains every native
!> breakpoint and the NWC lookup (nw_interp1d) is linear, each airfoil curve is
!> reproduced exactly. The per-node tables are then assembled from the airfoil
!> tables via BlAFID, exactly as input.f90 did. r_vtp follows Convention B.
!>
!> Runs once, on rank 0 (the only rank that runs velocity_correction), at the
!> first coupling step. Also calls init_nw_smearing_correction.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_build_nwc_arrays

    INTEGER(iwp) :: n_total, n_union, ia, ir, j, k, af
    INTEGER(iwp) :: n_airfoils
    REAL(wp), DIMENSION(:),   ALLOCATABLE :: all_alpha, alpha_union
    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: cl_union, cd_union
    REAL(wp), PARAMETER :: dedup_tol = 1.0e-6_wp   !< deg; merges floating-point duplicates only
    REAL(wp) :: tmp

    n_airfoils = nwc_recv_nairfoils

!-- The blade-node count from OpenFAST must match fast_n_blade_elem_max, which
!-- sizes the NWC history arrays. input.f90 enforced the same contract against
!-- Blade.dat; a mismatch would corrupt the per-node lookup, so abort loudly.
    IF ( nwc_recv_nnodes /= fast_n_blade_elem_max )  THEN
       WRITE( message_string, * )  'NWC blade-node count received from OpenFAST (',                 &
          nwc_recv_nnodes, ') does not match fast_n_blade_elem_max (', fast_n_blade_elem_max,       &
          '). Set fast_n_blade_elem_max to the AD15 blade-node count.'
       CALL message( 'f8c_build_nwc_arrays', 'F8C0026', 2, 2, 0, 6, 0 )
    END IF

!-- 1) union alpha grid: gather every airfoil's alpha breakpoints, sort ascending,
!--    then drop floating-point duplicates.
    n_total = 0
    DO ia = 1, n_airfoils
       n_total = n_total + nwc_recv_nrows(ia)
    END DO
    ALLOCATE( all_alpha(n_total) )
    k = 0
    DO ia = 1, n_airfoils
       DO ir = 1, nwc_recv_nrows(ia)
          k = k + 1
          all_alpha(k) = nwc_recv_alpha(ir, ia)
       END DO
    END DO

!-- insertion sort (n_total is a few hundred; this runs once)
    DO j = 2, n_total
       tmp = all_alpha(j)
       k = j - 1
       DO WHILE ( k >= 1 )
          IF ( all_alpha(k) <= tmp )  EXIT
          all_alpha(k+1) = all_alpha(k)
          k = k - 1
       END DO
       all_alpha(k+1) = tmp
    END DO

    ALLOCATE( alpha_union(n_total) )
    n_union = 1
    alpha_union(1) = all_alpha(1)
    DO j = 2, n_total
       IF ( all_alpha(j) - alpha_union(n_union) > dedup_tol )  THEN
          n_union = n_union + 1
          alpha_union(n_union) = all_alpha(j)
       END IF
    END DO

!-- 2) interpolate each airfoil onto the union grid. All NREL 5 MW airfoils span
!--    [-180, 180] deg, so no extrapolation occurs.
    ALLOCATE( cl_union(n_union, n_airfoils), cd_union(n_union, n_airfoils) )
    DO ia = 1, n_airfoils
       DO j = 1, n_union
          cl_union(j, ia) = nwc_interp1d_local( nwc_recv_alpha(1:nwc_recv_nrows(ia), ia),           &
                                                nwc_recv_cl(1:nwc_recv_nrows(ia), ia),              &
                                                nwc_recv_nrows(ia), alpha_union(j) )
          cd_union(j, ia) = nwc_interp1d_local( nwc_recv_alpha(1:nwc_recv_nrows(ia), ia),           &
                                                nwc_recv_cd(1:nwc_recv_nrows(ia), ia),              &
                                                nwc_recv_nrows(ia), alpha_union(j) )
       END DO
    END DO

!-- 3) allocate the NWC module arrays now that the polar grid size is known
    n_airfoil_polar_rows = n_union
    ALLOCATE( chord_arr(fast_n_blade_elem_max) )
    ALLOCATE( twist_arr(fast_n_blade_elem_max) )
    ALLOCATE( r_cp_arr(fast_n_blade_elem_max) )
    ALLOCATE( r_vtp_arr(fast_n_blade_elem_max + 1) )
    ALLOCATE( alpha_2d_arr(n_union, fast_n_blade_elem_max) )
    ALLOCATE( cl_2d_arr(n_union, fast_n_blade_elem_max) )
    ALLOCATE( cd_2d_arr(n_union, fast_n_blade_elem_max) )
    ALLOCATE( is_lifting_section_arr(fast_n_blade_elem_max) )

!-- 4) per-node geometry (twist already in degrees; r_cp from the rotor centre)
    DO j = 1, fast_n_blade_elem_max
       r_cp_arr(j)  = nwc_recv_r_cp(j)
       twist_arr(j) = nwc_recv_twist_deg(j)
       chord_arr(j) = nwc_recv_chord(j)
    END DO

!-- 5) per-node airfoil tables via BlAFID, exactly as input.f90 did
    DO j = 1, fast_n_blade_elem_max
       af = nwc_recv_blafid(j)
       IF ( af < 1  .OR.  af > n_airfoils )  THEN
          WRITE( message_string, * )  'NWC BlAFID out of range at blade node ', j, ': ', af
          CALL message( 'f8c_build_nwc_arrays', 'F8C0027', 2, 2, 0, 6, 0 )
       END IF
       DO k = 1, n_union
          alpha_2d_arr(k, j) = alpha_union(k)
          cl_2d_arr(k, j)    = cl_union(k, af)
          cd_2d_arr(k, j)    = cd_union(k, af)
       END DO
    END DO

!-- 6) trailed-vortex radii, Convention B: root and tip vortices coincide with the
!--    first and last control points, interior ones sit at the CP midpoints.
    r_vtp_arr(1) = r_cp_arr(1)
    DO j = 2, fast_n_blade_elem_max
       r_vtp_arr(j) = 0.5_wp * ( r_cp_arr(j-1) + r_cp_arr(j) )
    END DO
    r_vtp_arr(fast_n_blade_elem_max + 1) = r_cp_arr(fast_n_blade_elem_max)

!-- 7) lifting-section classification (same 0.01 cutoff as before): cylinder
!--    sections have Cl identically zero and end up FALSE.
    is_lifting_section_arr(:) = MAXVAL( ABS( cl_2d_arr(:,:) ), DIM=1 ) > 0.01_wp

!-- 8) initialise the NWC module-level scalars (n_rotors, n_blades, n_s, n_v) and
!--    the Pirrung 2017b constant tables (nw_P, nw_N), once.
    CALL init_nw_smearing_correction( nturbines, fast_n_blades_max, fast_n_blade_elem_max )

    DEALLOCATE( all_alpha, alpha_union, cl_union, cd_union )

    nwc_arrays_populated = .TRUE.

 CONTAINS

!-- Linear interpolation with endpoint clamping (matches the NWC's nw_interp1d).
!-- xs must be sorted ascending, which the native AeroDyn alpha tables are.
    PURE FUNCTION nwc_interp1d_local( xs, ys, n, xq )  RESULT( yq )
       INTEGER(iwp), INTENT(IN) ::  n
       REAL(wp),     INTENT(IN) ::  xs(n), ys(n), xq
       REAL(wp) ::  yq, t
       INTEGER(iwp) ::  lo, hi, mid
       IF ( xq <= xs(1) )  THEN
!--       Extrapolate below the table with the endpoint slope so a re-gridded
!--       value matches the runtime nw_interp1d lookup, which extrapolates rather
!--       than clamps. For airfoils that cover the full union range this branch is
!--       never taken; it only matters when one airfoil's alpha range is shorter.
          IF ( n >= 2 )  THEN
             yq = ys(1) + ( xq - xs(1) ) / ( xs(2) - xs(1) ) * ( ys(2) - ys(1) )
          ELSE
             yq = ys(1)
          END IF
       ELSE IF ( xq >= xs(n) )  THEN
          IF ( n >= 2 )  THEN
             yq = ys(n) + ( xq - xs(n) ) / ( xs(n) - xs(n-1) ) * ( ys(n) - ys(n-1) )
          ELSE
             yq = ys(n)
          END IF
       ELSE
          lo = 1
          hi = n
          DO WHILE ( hi - lo > 1 )
             mid = ( lo + hi ) / 2
             IF ( xs(mid) <= xq )  THEN
                lo = mid
             ELSE
                hi = mid
             END IF
          END DO
          t  = ( xq - xs(lo) ) / ( xs(hi) - xs(lo) )
          yq = ys(lo) + t * ( ys(hi) - ys(lo) )
       END IF
    END FUNCTION nwc_interp1d_local

 END SUBROUTINE f8c_build_nwc_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_init

    USE ISO_C_BINDING,                                                                             &
       ONLY: C_CHAR,                                                                               &
             C_INT,                                                                                &
             C_LOC,                                                                                &
             C_NULL_CHAR,                                                                          &
             C_NULL_PTR,                                                                           &
             C_PTR


    INTEGER(iwp) ::  i  !< loop variable
    INTEGER(iwp) ::  j  !< loop variable
    INTEGER(iwp) ::  k  !< loop variable
    INTEGER(iwp) ::  o  !< loop variable

    INTEGER(iwp) ::  iflag    !<
    INTEGER(iwp) ::  endflag  !<

    INTEGER(iwp) ::  i_hub  !<
    INTEGER(iwp) ::  j_hub  !<
    INTEGER(iwp) ::  k_hub  !<

    INTEGER(iwp) ::  tower_l  !<
    INTEGER(iwp) ::  tower_r  !<
    INTEGER(iwp) ::  tower_n  !<
    INTEGER(iwp) ::  tower_s  !<

    REAL(wp) ::  rcx  !<
    REAL(wp) ::  rcy  !<
    REAL(wp) ::  rcz  !<

    TYPE(C_PTR) ::  fast_host_addr_c(SIZE(fast_host_addr) + 1)  !<
    TYPE(C_PTR) ::  fast_host_port_c(SIZE(fast_host_port) + 1)  !<

    TYPE string
      CHARACTER(LEN=:,KIND=C_CHAR), ALLOCATABLE ::  item  !<
    END TYPE string

    TYPE(string), TARGET :: tmp_addr(SIZE(fast_host_addr))  !<
    TYPE(string), TARGET :: tmp_port(SIZE(fast_host_port))  !<

#if defined( __fastv8 )
    INTERFACE

       FUNCTION init_comm( ifnumturb, host_addr, host_port ) BIND(C,NAME='init_comm' )
          USE ISO_C_BINDING,                                                                       &
              ONLY:  C_INT,                                                                        &
                     C_PTR
          INTEGER(kind=C_INT) ::  init_comm
          INTEGER(kind=C_INT), VALUE ::  ifnumturb
          TYPE(C_PTR), INTENT(IN) ::  host_addr(*)
          TYPE(C_PTR), INTENT(IN) ::  host_port(*)
       END FUNCTION init_comm

       FUNCTION set_nwc_enabled( flag ) BIND(C,NAME='set_nwc_enabled')
          USE ISO_C_BINDING,                                                                       &
              ONLY:  C_INT
          INTEGER(kind=C_INT) ::  set_nwc_enabled
          INTEGER(kind=C_INT), VALUE ::  flag
       END FUNCTION set_nwc_enabled

    END INTERFACE
#endif

!
!-- Compute values of the exponential function.
    expvalue(:) = 0.0_wp
    DO  i = 0, 10000000
       expvalue(i) = EXP(-0.001_wp*i)
       ! s.ouchene add _wp suffix to set the correct REAL kind
       IF ( expvalue(i) < 0.00000001_wp )  THEN
          EXIT
       ENDIF
    ENDDO

    endflag = 0

!
!-- Only do if turbines are supposed to be simulated.
    IF ( nturbines > 0 )  THEN

       IF ( myid == 0 )  THEN

          DO  i = 1, nturbines
             IF ( TRIM( fast_host_addr(i) ) == '' )  THEN
               fast_host_addr(i) = fast_host_addr_default
             ENDIF
          ENDDO

          DO  i = 1, SIZE( fast_host_addr )
!
!--          This may involve kind conversion.
             tmp_addr(i)%item = TRIM( fast_host_addr(i) ) // C_NULL_CHAR
             fast_host_addr_c(i) = C_LOC( tmp_addr(i)%item )
          ENDDO
          fast_host_addr_c(SIZE( fast_host_addr )) = C_NULL_PTR

          DO  i = 1, SIZE( fast_host_port )
!
!--          This may involve kind conversion.
             tmp_port(i)%item = TRIM( fast_host_port(i) ) // C_NULL_CHAR
             fast_host_port_c(i) = C_LOC( tmp_port(i)%item )
          END DO
          fast_host_port_c(SIZE( fast_host_port )) = C_NULL_PTR

!
!--       Initializing the connection with FAST.
#if defined( __fastv8 )
          iflag = init_comm( nturbines, fast_host_addr_c, fast_host_port_c )
#else
          iflag = -1
#endif
          IF ( iflag /= 0 )  THEN
             WRITE( message_string, * )  'Unable to initalize communication with FAST server(s)'
             CALL message( 'f8c_actions', 'F8C0012', 2, 2, 0, 6, 0 )
          ENDIF

!
!--       S. Ouchene: tell the C client whether PALM wants the near-wake-correction
!--       blade geometry and airfoil polars appended to the connection-init
!--       response. This must be set before the COM_INIT_CONNECTION message is
!--       built. With the flag off the connection-init message is unchanged, so
!--       runs without the near-wake correction are byte-for-byte identical.
#if defined( __fastv8 )
          IF ( use_near_wake_smearing_correction )  THEN
             iflag = set_nwc_enabled( 1_C_INT )
          ELSE
             iflag = set_nwc_enabled( 0_C_INT )
          ENDIF
#endif

       ENDIF

    ENDIF !(nturbines > 0)

!    current_time_fast = 0.0_wp

!
!-- Determine the area within each grid cell that overlaps with the area of the nacelle and the
!-- tower (needed for calculation of the forces).
!-- Note: so far this is only a 2D version, in that the mean flow is
!-- perpendicular to the rotor area.
    tower_area_x_4d(:,:,:,:) = 0.0_wp
    tower_area_y_4d(:,:,:,:) = 0.0_wp

    IF ( myid == 0  .AND.  debug_output )  THEN
       WRITE(9,*) "f8c_init - before determination of tower_area"
       FLUSH( 9 )
    ENDIF

!
!-- Loop over all turbines.
    DO  o = 1, nturbines

       rcx = palm_tower_ref_pos_x(o)
       rcy = palm_tower_ref_pos_y(o)
       rcz = palm_tower_ref_pos_z(o) + htow(o)

       tower_area_x(:,:,:) = 0.0_wp
       tower_area_y(:,:,:) = 0.0_wp

!
!--    Determine the indices of the hub height.
       i_hub = INT( rcx / dx )
       j_hub = INT( ( rcy + 0.5_wp * dy ) / dy )
       k_hub = INT( ( rcz + 0.5_wp * dz(1) ) / dz(1) )

!
!--    Determine the indices of the grid boxes containing the left and
!--    the right boundaries of the tower.
       tower_n = INT( ( rcy + 0.5_wp * dtow(o) - 0.5_wp * dy ) / dy )
       tower_s = INT( ( rcy - 0.5_wp * dtow(o) - 0.5_wp * dy ) / dy )
       tower_r = INT( ( rcx + 0.5_wp * dtow(o) - 0.5_wp * dx ) / dx )
       tower_l = INT( ( rcx - 0.5_wp * dtow(o) - 0.5_wp * dx ) / dx )

       IF ( myid == 0  .AND.  debug_output )  THEN
          WRITE(9,*) 'id: ', i, ' tower_n ', tower_n
          WRITE(9,*) 'id: ', i, ' tower_s ', tower_s
          WRITE(9,*) 'id: ', i, ' tower_r ', tower_r
          WRITE(9,*) 'id: ', i, ' tower_l ', tower_l
          FLUSH( 9 )
       ENDIF
!
!--    Determine the fraction of the grid box area overlapping with the tower area.
       IF ( ( nxlg <= i_hub ) .AND. ( nxrg >= i_hub ) .AND.                                        &
            ( nysg <= j_hub ) .AND. ( nyng >= j_hub ) )                                            &
       THEN
!
!--       Loop from the southernmost grid index of the tower to the northernmost.
          DO  k = nzb, k_hub
!
!--          Loop from the southernmost grid index of the tower to the northernmost.
             DO  j = tower_s, tower_n
!
!--             If tower not completely inside one grid box.
                IF ( tower_n - tower_s >= 1 )  THEN

                   IF ( j == tower_s )  THEN
                      tower_area_x(k,j,i_hub) =                                                    &
                                                ! extension in z-direction
                                                MIN(  rcz - ( k * dz(1) - 0.5_wp * dz(1) ),        &
                                                      dz(1) ) *                                    &
                                                ! extension in y-direction
                                                ( ( tower_s + 1 + 0.5_wp ) * dy                    &
                                                - ( rcy - 0.5_wp * dtow(o) ) )
                   ELSEIF ( j == tower_n )  THEN
                      tower_area_x(k,j,i_hub) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ) ,        &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( rcy + 0.5_wp * dtow(o)                           &
                                                - ( tower_n + 0.5_wp ) * dy )
                   ELSE
!
!--                   Grid boxes inbetween (where tower_area = grid box area):
                      tower_area_x(k,j,i_hub) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) * dy
                   ENDIF
                ELSE
!
!--                Tower lies completely within one grid box.
                   tower_area_x(k,j,i_hub) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),            &
                                                  dz(1) ) * dtow(o)
                ENDIF

             ENDDO

!
!--          Loop from the left grid index of the tower to the right index
             DO  i = tower_l, tower_r
!
!--             If tower not completely inside one grid box
                IF ( tower_r - tower_l >= 1 )  THEN
!--                leftmost and rightmost grid box:
                   IF ( i == tower_l )  THEN

                      tower_area_y(k,j_hub,i) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( ( tower_l + 1 ) * dx                             &
                                                - ( rcx - 0.5_wp * dtow(o) ) )
                   ELSEIF ( i == tower_r )  THEN
                      tower_area_y(k,j_hub,i) =                                                    &
                                                ! extension in z-direction
                                                MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) *                                     &
                                                ! extension in y-direction
                                                ( rcx + 0.5 * dtow(o) - tower_r * dx )
                   ELSE
!
!--                   Grid boxes inbetween (where tower_area = grid box area).
                      tower_area_y(k,j_hub,i) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),         &
                                                     dz(1) ) * dx
                   ENDIF
                ELSE
!
!--                Tower lies completely within one grid box.
                   tower_area_y(k,j_hub,i) = MIN( rcz - ( k * dz(1) - 0.5_wp * dz(1) ),            &
                                                  dz(1) ) * dtow(o)
                ENDIF
             ENDDO
          ENDDO
       ENDIF !( ( nxlg <= i_hub ) .AND. ( nxrg >= i_hub ) .AND. ( nysg <= j_hub ) .AND. ( nyng >= j_hub ) )

       CALL exchange_horiz( tower_area_x, nbgp )
       CALL exchange_horiz( tower_area_y, nbgp )
!
!--    Restore tower area field in 4D-field containing each turbine.
       tower_area_x_4d(o,:,:,:) = tower_area_x(:,:,:)
       tower_area_y_4d(o,:,:,:) = tower_area_y(:,:,:)
!
!--    Tabulate the points on the circle that are required in the following for
!--    the calculation of the Riemann integral (node points; they are called
!--    circle_points in the following).

    ENDDO  ! end of loop over turbines

!
!-- Initializing 3D force arrays
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
         DO  k = nzb, nzt
             force_x(k,j,i) = 0.0_wp
             force_y(k,j,i) = 0.0_wp
             force_z(k,j,i) = 0.0_wp
          ENDDO
       ENDDO
   ENDDO

 END SUBROUTINE f8c_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which output quantities are defined.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x  !< x grid of output variable
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y  !< y grid of output variable
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z  !< z grid of output variable
    CHARACTER(LEN=*), INTENT(IN)  ::  var     !< name of output variable

    LOGICAL, INTENT(OUT) ::  found   !< flag if output variable is found

    found  = .TRUE.
!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'thrx', 'thrx_xy', 'thrx_xz', 'thrx_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       u grid
          ELSE
             grid_x = 'xu'
             grid_y = 'y'
             grid_z = 'zu'
          ENDIF

       CASE ( 'tory', 'tory_xy', 'tory_xz', 'tory_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       v grid
          ELSE
             grid_x = 'x'
             grid_y = 'yv'
             grid_z = 'zu'
          ENDIF

       CASE ( 'torz', 'torz_xy', 'torz_xz', 'torz_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       w grid
          ELSE
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zw'
          ENDIF

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE f8c_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with user-defined information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_header( io )

    INTEGER(iwp) ::  io  !<


    WRITE ( io, 1 )
    WRITE ( io, 20 )
    WRITE ( io, 21 ) nturbines

!
!-- Format-descriptors
1   FORMAT ( //' FASTv8 coupler information:'/ ' ------------------------------------------'/ )

20  FORMAT ( '--> Essential parameters:' )
21  FORMAT ( '       number of active turnines    :   nturbines   = ', I3)

 END SUBROUTINE f8c_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_actions( location )

    USE control_parameters,  &
        ONLY: initializing_actions    

    CHARACTER(LEN=*) ::  location  !<
    CHARACTER(LEN=20):: time_str  !!!! Added by Samir for debugging
!
!-- loop variables
    INTEGER(iwp) ::  i    !<
    INTEGER(iwp) ::  j    !<
    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  l    !<
    INTEGER(iwp) ::  o    !<
    INTEGER(iwp) ::  m    !<
    INTEGER(iwp) ::  n    !<
    INTEGER(iwp) ::  n_s  !<
!
!-- grid point numbers
    INTEGER(iwp) ::  ii  !<
    INTEGER(iwp) ::  jj  !<
    INTEGER(iwp) ::  kk  !<
!
!-- run control
    INTEGER(iwp) ::  iflag  !<

    INTEGER(iwp) ::  expargument  !<
    INTEGER(iwp) ::  len_2d       !< airfoil polar table row count, fed to velocity_correction

#if defined( __parallel )
    INTEGER(iwp) ::  size_array_type_1  !<
#endif
!
!-- Velocity interpolation and force distribution
    REAL(wp) ::  x  !<
    REAL(wp) ::  y  !<
    REAL(wp) ::  z  !<

    REAL(wp) ::  aa  !<
    REAL(wp) ::  bb  !<
    REAL(wp) ::  cc  !<
    REAL(wp) ::  dd  !<
    REAL(wp) ::  gg  !<

    REAL(wp) ::  u_int_l  !<
    REAL(wp) ::  u_int_u  !<
    REAL(wp) ::  v_int_l  !<
    REAL(wp) ::  v_int_u  !<
    REAL(wp) ::  w_int_l  !<
    REAL(wp) ::  w_int_u  !<

    REAL(wp) ::  eps_kernel    !<
    REAL(wp) ::  eps_kernel_2  !<
    REAL(wp) ::  nenner        !<
    REAL(wp) ::  kernel_value  !<
    REAL(wp) ::  distance_hub  !<
    REAL(wp) ::  fboxelen      !<

    REAL(wp), DIMENSION(3) ::  tmp_res  !<

    REAL(wp), DIMENSION(:,:),     ALLOCATABLE ::  u_int_h  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  u_int_b  !<

    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  x_prime_sec  !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  y_prime_sec  !<
    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE ::  z_prime_sec  !<

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  x_prime  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  y_prime  !<
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  z_prime  !<

    LOGICAL ::  palm_dt_3d = .FALSE.

    ! SoS >
    REAL(wp)    :: eps
    ! < SoS
#if defined( __fastv8 )
    INTERFACE

       FUNCTION commclient(target, message) BIND(C,NAME='commclient')
          use ISO_C_BINDING
          INTEGER(kind=C_INT) ::  commclient
          INTEGER(kind=C_INT), value ::  target
          INTEGER(kind=C_INT), value ::  message
       END FUNCTION commclient

       FUNCTION time_data(simtime, simdt, simstep, simturbon, simendtime) BIND(C,NAME='time_data')
          use ISO_C_BINDING
          INTEGER(kind=C_INT) ::  time_data
          REAL(kind=c_double), value ::  simtime
          REAL(kind=c_double), value ::  simdt
          INTEGER(kind=C_INT), value ::  simstep
          REAL(kind=c_double), value ::  simturbon
          REAL(kind=c_double), value ::  simendtime
       END FUNCTION time_data

    END INTERFACE
#endif

    CALL cpu_log( log_point(24), 'f8c_actions', 'start' )

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( myid == 0 )  THEN

             IF ( debug_output )  THEN
                WRITE(9,*) 'starting new time step', simulated_time, current_timestep_number
                FLUSH(9)
                PRINT*, 'first_restart_with_wt', first_restart_with_wt
             ENDIF
 
             IF ( first_restart_with_wt ) then
                  palm_time_at_begin = simulated_time_at_begin                
             ENDIF !SoS
             IF ( TRIM( initializing_actions ) /= 'read_restart_data' ) THEN
                palm_time_at_begin = 0.0_wp          
             ENDIF
!
!--          update information between fortran and c code
#if defined( __fastv8 )
             IF ( simulated_time - palm_time_at_begin  < time_turbine_on )  THEN
                iflag = time_data(simulated_time - palm_time_at_begin, dt_3d, current_timestep_number,                  &
                                  time_turbine_on, end_time)
             ELSE
!
!--             Hier wird die Zeit von PALM an die .c Datei geschickt
                iflag = time_data(simulated_time - palm_time_at_begin, dt_fast, current_timestep_number,                &
                                  time_turbine_on, end_time)
             ENDIF !( simulated_time < time_turbine_on )
#else
             iflag = -1
#endif

             IF ( iflag /= 0 )  THEN
                WRITE( message_string, * )  'Unable to update information between ',       &
                   'FORTRAN and C:', simulated_time, ':', current_timestep_number
                CALL message( 'f8c_actions', 'F8C0013', 2, 2, 0, 6, 0 )
             ENDIF !( iflag /= 0 )

          ENDIF !( myid == 0 )

!--       Check if connection to all FAST instances can be established and exchange some data for initialisation.
          IF ( myid == 0 )  THEN
             IF ( current_timestep_number == 0 )  THEN

               CALL location_message( 'Checking connection to FAST server(s)', 'start' )
#if defined( __fastv8 )
                iflag = commclient( COM_TARGET_OPENFAST, COM_INIT_CONNECTION )
#else
                iflag = -1
#endif

                IF ( iflag < 0 )  THEN
                   WRITE( message_string, * )                                                      &
                          'Unable to establish connection with all FAST server(s)'
                   CALL message( 'f8c_actions', 'F8C0014', 2, 2, 0, 6, 0 )
                ELSEIF ( iflag == 1 )  THEN
                   WRITE( message_string, * )  'At least one FAST server finished the simulation'
                   CALL message( 'f8c_actions', 'F8C0015', 2, 2, 0, 6, 0 )
                ELSEIF ( iflag == 0 )  THEN
                   CALL location_message(                                                          &
                           'Connection to FAST server(s) successfully established', '' )
                ENDIF !( iflag < 0 )
                CALL location_message( 'Checking connection to FAST server(s)',    'finished' )
             ENDIF !( current_timestep_number == 0 )

!--          Check if connection to all FAST instances can again be established and receive some data for initialisation
             IF ( TRIM( initializing_actions ) == 'read_restart_data' ) THEN
                IF ( first_ts_with_wt ) THEN
                   CALL location_message(                                                          &
                          'Resume coupled simulation with FAST ...', '' )
#if defined( __fastv8 )
                   IF ( first_restart_with_wt ) THEN  !SoS
                      PRINT*, "First restart with Turbine."
                      iflag = commclient( COM_TARGET_OPENFAST, COM_INIT_CONNECTION)
                   ELSE
                      PRINT*, "Restart with Turbine (not first restart)."
                      iflag = commclient(COM_TARGET_OPENFAST, COM_RESUME_SIMULATION)
                   ENDIF  !SoS
#else
                   iflag = -1
#endif
                   IF (iflag == 0 ) THEN
                      CALL location_message(                                                       &
                           'Connection to FAST server(s) successfully established.', '' )
                   END IF
                END IF
             END IF     


          ENDIF !( myid == 0 )

!--       Errechnen wie viele FAST sub-Zeitschritte es geben muss
          IF ( myid == 0 )  THEN
!--           TODO: s.ouchene for CFL-driven runs adjust dt_3d to be always the smallest integer multiple of dt_fast
!--           using: dt_3d = floot(dt_3d/dt_palm + 4.0_wp*spacing(dt_3d/dt_fast))*dt_fast
!--           I need to figure out that this wouldn't cause any logical errors in the simulation and need to properly
!--           check how to detect if time stepping is driven by CFL reliably. one could compare the current dt_3d
!--           with the previous one to assess that. Need a reliable way.
             dt_palm = dt_3d
             palm_dt_3d = .TRUE.
             dt_ratio = dt_palm / dt_fast    ! no need to broadcast below since this is needed only by rank 0

!--          s.ouchene: Add a tiny correction to FLOOR argument to handle cases where a value that should be an integer is
!--          stored just below due to floating-point rounding, e.g. in double precision, 0.3/0.1 giving
!--          2.9999999999999996 instead of 3.0. SPACING(dt_ratio) makes the correction match the floating-point
!--          gap near dt_ratio. 4.0 is just a safety factor.
             n_dt_fast = FLOOR( dt_ratio + 4.0_wp*SPACING(dt_ratio) ) 
             
          ENDIF !( myid == 0 )

#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_BCAST( dt_fast, 1, MPI_REAL, 0, comm2d, ierr )
          CALL MPI_BCAST( dt_palm, 1, MPI_REAL, 0, comm2d, ierr )
          CALL MPI_BCAST(first_restart_with_wt, 1, MPI_LOGICAL, 0, comm2d, ierr) !SoS
          CALL MPI_BCAST(palm_time_at_begin, 1, MPI_REAL, 0, comm2d, ierr) !SoS
          CALL MPI_BCAST( fast_n_radius, nturbines, MPI_REAL, 0, comm2d, ierr )  ! [todo]
          CALL MPI_BARRIER( comm2d, ierr )
#endif

          IF ( TRIM( initializing_actions ) == 'read_restart_data' ) THEN
             time_turbine_on = 0.0
          END IF     

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN

             dt_3d = dt_palm
!
!--          Besondere Aktionen in Zusammenhang mit der Kopplung von
!--          PALM und FAST sind nur dann durchzufuehren, wenn die
!--          Windenergieanlagen eingeschaltet sind; der Zeitpunkt des
!--          Einschaltens der Anlagen ist ueber den Parameter
!--          time_turbine_on festgelegt, der in der
!--          Namelistparameterdatei unter userpar vorzugeben ist
!
!--          Beim ersten Zeitschritt ist der Ablauf etwas anders
!--          als an den anderen Zeitschritten. Zum ersten Zeitschritt,
!--          anders als zu den anderen Zeitschritten, erhaelt PALM noch
!--          keine Informationen ueber Kraefte aus FAST, sondern
!--          lediglich Informationen ueber die Positionen der
!--          Blattelemente
!
!--          first time step with wind turbine: Sending of signal to FAST to show that PALM is ready
!--          Receive blade element positions from FAST
!--          all other time steps with wind turbine: asking for positions
             IF ( myid == 0 )  THEN

                IF ( first_ts_with_wt )  THEN

                IF ( TRIM( initializing_actions ) /= 'read_restart_data' ) THEN

#if defined( __fastv8 )
                   iflag = commclient(COM_TARGET_OPENFAST, COM_PALM_IS_READY)
#else
                   iflag = -1
#endif

                   IF ( iflag < 0 )  THEN
                      WRITE( message_string, * )  'Unable to send start signale to FAST server(s)'
                      CALL message( 'f8c_actions', 'F8C0016', 2, 2, 0, 6, 0 )
                   ELSEIF ( iflag == 1 )  THEN
                      WRITE( message_string, * )  'At least one FAST server finished the simulation'
                      CALL message( 'f8c_actions', 'F8C0017', 2, 2, 0, 6, 0 )
                   ELSEIF ( iflag == 0 )  THEN
                      CALL location_message(                                                       &
                              'Start signal to FAST servers successfully broadcasted', '' )
                   ENDIF !( iflag < 0 )

                 ELSE

#if defined( __fastv8 )
                   iflag = commclient(COM_TARGET_OPENFAST, COM_PALM_IS_READY)
#else
                   iflag = -1
#endif

                 END IF      
!
!--                Determination of the corners of the boxes around every turbine for force distribution
                   IF ( iflag == 0 )  THEN

                      DO  i = 1, nturbines

                         fboxelen = fast_n_radius(i) * fbox_fac

                         fboxcorners(i, 1, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - fboxelen ) * ddx                          &
                                                     )

                         fboxcorners(i, 1, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) + fboxelen ) * ddx                          &
                                                       )

                         fboxcorners(i, 2, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - fboxelen ) * ddy                          &
                                                     )

                         fboxcorners(i, 2, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) + fboxelen ) * ddy                          &
                                                       )

                         fboxcorners(i, 3, 1) = FLOOR(                                             &
                            ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - fboxelen ) / dz(1)                        &
                                                     )

                         fboxcorners(i, 3, 2) = CEILING(                                           &
                            ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) + fboxelen ) / dz(1)                        &
                                                       )
                      ENDDO

                   ENDIF !( iflag == 0)

                   first_ts_with_wt = .FALSE.

                ENDIF !( first_ts_with_wt )

             ENDIF !( myid == 0 )

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
!
!--          sending of fboxcorners
             CALL MPI_BCAST( fboxcorners(1,1,1), nturbines*3*2, MPI_INTEGER, 0, comm2d, ierr )
!
!--          Senden von fast_n_blades (nturbines) -> senden Anzahl der Blätter für Turbine 1
             CALL MPI_BCAST( fast_n_blades(1), nturbines, MPI_INTEGER, 0, comm2d, ierr )
!
!--          Senden von fast_n_blade_elem (nturbines) -> senden Anzahl der Blattelemente von Turbine 1
             CALL MPI_BCAST( fast_n_blade_elem(1), nturbines, MPI_INTEGER, 0, comm2d, ierr )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
!
!--          Counter for  number of FAST timesteps for a single PALM timestep
             l = 0
             IF ( myid == 0 )  THEN
                current_time_fast_min = MINVAL(current_time_fast(:))

                DO i = 1, nturbines
                   IF ( current_time_fast_min < current_time_fast(i) )  THEN
                      WRITE( message_string, * )  'Different times on different FAST instances: ', &
                         'i, current_time_fast_min, current_time_fast(i) ',                        &
                          i, current_time_fast_min, current_time_fast(i)
                      CALL message( 'f8c_actions', 'F8C0018', 2, 2, 0, 6, 0 )
                   ENDIF
                END DO

             ENDIF

             target_time = simulated_time + dt_palm - palm_time_at_begin
!--          Loop for the section method: while FAST is catching up to PALM,
!--          PALM is frozen while the communication with FAST continues as before.

!--          s.ouchene: due to floating-point arithmetic, the condition is not always satisfied.
!--          A trick is to add half fast timestep and change "<=" to strict "<"
!--          I checked that PALM is accumulating time steps with the simplest approach which would
!--          cause desynchronisation between PALM and FAST for longer runs. PALM should use kahan
!--          summation to fix that. TODO: open a ticket
             Do While ( current_time_fast(1) < target_time + 0.5_wp * dt_fast )

#if defined( __parallel )
!
!--             Kommunizieren der Informationen über Positionen, die
!--             von PE0 erhalten wurden, sind auch an alle anderen PEs
!--             Dies laeuft ab wie in der Kopplung mit FAST
!
!--             Größe der Felder, deren Werte von PE0 allen anderen
!--             Prozessoren mitgeteilt werden
                size_array_type_1 = nturbines * fast_n_blades_max * fast_n_blade_elem_max * 3
                CALL MPI_BARRIER( comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_x(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_y(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_tower_ref_pos_z(1), nturbines,                                &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_hub_center_pos(1,1), nturbines*3,                             &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_position_blade(1,1,1,1), size_array_type_1,                   &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_force_blade(1,1,1,1), size_array_type_1,                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( shaft_height_fast, nturbines,                                      &
                                MPI_REAL, 0, comm2d, ierr )

                CALL MPI_BARRIER( comm2d, ierr )
#endif
!
!--             As in the code for the coupling of PALM and FAST:
!--             Interpolation of the velocity field computed in PALM
!--             onto the blade positions that were provided by FAST
!
!--             Allocation of temporary fields for the interpolated velocities
!--             1. u_int_h: Velocity at the hub position, three components for each turbine
                ALLOCATE(u_int_h(1:nturbines,1:3))
!
!--             2. u_int_b: Velocity at the rotor blade, three components for each blade segment of the n blades of a turbine
                ALLOCATE(u_int_b(1:fast_n_blade_elem_max,1:fast_n_blades_max,1:nturbines,1:3))
!
!--             Initialization of the interpolated fields with zero values; this is actually
!--             not necessary, since these fields will in any case be assigned a value on each
!--             processor during the interpolation procedure (if the point to be interpolated
!--             lies outside the subdomain of a processor, the corresponding field element is
!--             assigned the value 0)
                u_int_h(:,:) = 0.0_wp
                u_int_b(:,:,:,:) = 0.0_wp

!--             count intermediate FAST time step
                l = l + 1

!--             s.ouchene: Wrap hub velocity intepolation inside IF (l == 1) to avoid interpolation every OpenFAST sub-step    
!               l == 1 is the first FAST time step at the beginning of each PALM time step
                IF ( l == 1 ) THEN

                  palm_hub_center_vel(:, :) = 0.0_wp
                  
                  DO  i = 1, nturbines

                   IF ( myid == 0 .AND. debug_output )  THEN
                      WRITE(*,*) "==================================== DEBUG ===================================================="  
                      WRITE(*,*) "Time(s): PALM: ", simulated_time, " ---- ", " FAST:", current_time_fast  
                      WRITE(*,*) "Turbine ID: ", i,  ", Rotspeed: ", rotspeed(i)
                      WRITE(*,*) ""
                      !FLUSH(9)
                   ENDIF
!
!--                Bilineares Interpolationsverfahren:
!--                1. Interpolation des Geschwindigkeitsfeldes auf die Position des Hubs --> palm_hub_center_vel
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von u:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                          fast_hub_center_pos(i,1) ) * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                          fast_hub_center_pos(i,2) - 0.5_wp * dy) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                          fast_hub_center_pos(i,3) + 0.5_wp * dz(1) ) / dz(1) )

#ifdef __VERBOSE_MODE                        
                  print *, myid, 'hub', i, &
                 ' pos=', palm_tower_ref_pos_x(i)+fast_tower_ref_pos_x(i)+fast_hub_center_pos(i,1), &
                  palm_tower_ref_pos_y(i)+fast_tower_ref_pos_y(i)+fast_hub_center_pos(i,2), &
                  palm_tower_ref_pos_z(i)+fast_tower_ref_pos_z(i)+fast_hub_center_pos(i,3), &
                  ' idx=', ii, jj, kk, &
                  ' box=', nxl, nxr, nys, nyn, nzb, nzt
#endif

!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - (ii * dx)
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - 0.5_wp * dy - jj * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1) + 0.5_wp * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         u_int_l = ( ( aa ) * u(kk,jj,ii)     +    &
                                     ( bb ) * u(kk,jj,ii+1)   +    &
                                     ( cc ) * u(kk,jj+1,ii)   +    &
                                     ( dd ) * u(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         u_int_u = ( ( aa ) * u(kk+1,jj,ii)     +    &
                                     ( bb ) * u(kk+1,jj,ii+1)   +    &
                                     ( cc ) * u(kk+1,jj+1,ii)   +    &
                                     ( dd ) * u(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )

                         u_int_h(i,1) = (1/dz(1)) * ((dz(1)-z)*u_int_l + z*u_int_u)

#ifdef __VERBOSE_MODE
print *, 'WEIGHTS i=', i, ' x=', x, ' y=', y, ' z=', z, &
           ' dx=', dx, ' dy=', dy, ' dz(1)=', dz(1)
print *, 'INTERP  i=', i, ' aa=', aa, ' bb=', bb, ' cc=', cc,    &
           ' dd=', dd, ' gg=', gg,                                  &
           ' u_int_l=', u_int_l, ' u_int_u=', u_int_u,              &
           ' formula=', (1/dz(1)) * ((dz(1)-z)*u_int_l + z*u_int_u)


print *, 'HUB_INTERP i=', i,                                                  &
   ' u_stencil_l=', u(kk,jj,ii), u(kk,jj,ii+1), u(kk,jj+1,ii), u(kk,jj+1,ii+1), &
   ' u_stencil_u=', u(kk+1,jj,ii), u(kk+1,jj,ii+1), u(kk+1,jj+1,ii), u(kk+1,jj+1,ii+1), &
   ' u_int_h_1=', u_int_h(i,1)

#endif

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                         u_int_h(i,1) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,1) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von v:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                          fast_hub_center_pos(i,1) - 0.5_wp * dx)  * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                          fast_hub_center_pos(i,2) ) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                          fast_hub_center_pos(i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - ii * dx - 0.5_wp * dx
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - jj * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1) + 0.5_wp * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         v_int_l = ( ( aa ) * v(kk,jj,ii)     +    &
                                     ( bb ) * v(kk,jj,ii+1)   +    &
                                     ( cc ) * v(kk,jj+1,ii)   +    &
                                     ( dd ) * v(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         v_int_u = ( ( aa ) * v(kk+1,jj,ii)     +    &
                                     ( bb ) * v(kk+1,jj,ii+1)   +    &
                                     ( cc ) * v(kk+1,jj+1,ii)   +    &
                                     ( dd ) * v(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )
!
!--                      OLD COMMENT: ToDo check, why the computed values are not used and why u_int_h(i,2) is set here
!                        ! u_int_h(i, 2) = 0.0_wp ! this was set in the original version of the code without a clear explanation
!                        ! s.ouchene change to:
                         u_int_h(i,2) = (1.0_wp / dz(1)) * ((dz(1) - z ) * v_int_l + z * v_int_u)

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                        u_int_h(i,2) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,2) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                Nabenwindgeschwindigkeit an Turbine i, Interpolation von w:
                   ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                &
                        fast_hub_center_pos(i,1) - 0.5_wp * dx) * ddx )
                   jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                &
                        fast_hub_center_pos(i,2) - 0.5_wp * dy) * ddy )
                   kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                &
                        fast_hub_center_pos(i,3) ) / dz(1) )
!
!--                Auf einem Prozessorelement findet nur dann eine Interpolation statt,
!--                wenn alle Stützpunkte der Interpolation auf dem Prozessorelement zur
!--                Verfügung stehen
                   IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                      IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                         x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +                  &
                              fast_hub_center_pos(i,1) - ii * dx - 0.5_wp * dx
                         y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +                  &
                              fast_hub_center_pos(i,2) - jj * dy - 0.5_wp * dy
                         z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                  &
                              fast_hub_center_pos(i,3) - kk * dz(1)

                         aa = abs((dx-x)*(dy-y))
                         bb = abs((x)*(dy-y))
                         cc = abs((dx-x)*(y))
                         dd = x*y
                         gg = dx*dy

                         w_int_l = ( ( aa ) * w(kk,jj,ii)     +    &
                                     ( bb ) * w(kk,jj,ii+1)   +    &
                                     ( cc ) * w(kk,jj+1,ii)   +    &
                                     ( dd ) * w(kk,jj+1,ii+1) ) /  &
                                     ( gg )

                         w_int_u = ( ( aa ) * w(kk+1,jj,ii)     +    &
                                     ( bb ) * w(kk+1,jj,ii+1)   +    &
                                     ( cc ) * w(kk+1,jj+1,ii)   +    &
                                     ( dd ) * w(kk+1,jj+1,ii+1) ) /  &
                                     ( gg )
!
!--                      OLD COMMENT: ToDo check, why the computed values are not used and why u_int_h(i,3) is set here
!                        ! u_int_h(i, 3) = 0.0_wp ! this was set in the original version of the code without a clear explanation
!                        ! s.ouchene change to:
                         u_int_h(i,3) = (1.0_wp / dz(1)) * ((dz(1) - z) * w_int_l + z * w_int_u)

                      ELSE !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                         u_int_h(i,3) = 0.0_wp
                      ENDIF !( ( jj >= nys ) .AND. ( jj <= nyn ) )
                   ELSE !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
                      u_int_h(i,3) = 0.0_wp
                   ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
#ifdef __VERBOSE_MODE
                   print *, 'HUB_INTERP_3C i=', i, &
                  ' u_int_h=', u_int_h(i,1), u_int_h(i,2), u_int_h(i,3)
#endif
                ENDDO !i = 1, nturbines
              END IF ! l == 1 (hub interpolation)
!
!--             2. Interpolation des Geschwindigkeitsfeldes auf die Position der Blattelemente --> palm_vel_blade
! SoS  Added to look at first line of the new sector only ->

              
            IF ( l == 1 ) THEN

! SoS  Added to look at first line of the new sector only <-
                palm_vel_blade(:,:,:,:) = 0.0_wp

                DO  i = 1, nturbines
                   DO  j = 1, fast_n_blades(i)
                      DO  k = 1, fast_n_blade_elem(i)
!
!--                      Velocity at element k of blade j of turbine i, interpolation of u:
                         ii = FLOOR(( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +          &
                                fast_position_blade(k,j,i,1)) * ddx )                              ! wind speed in rotor area
                         jj = FLOOR(( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +          &
                                fast_position_blade(k,j,i,2) - 0.5_wp * dy) * ddy )
                         kk = FLOOR(( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +          &
                                fast_position_blade(k,j,i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                      On a processor element, interpolation only takes place
!--                      if all support points of the interpolation are available
!--                      on that processor element
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) - (ii * dx)                              ! wind speed in rotor area
                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy - 0.5_wp * dy
                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1) + 0.5_wp * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               u_int_l = ( ( aa ) * u(kk,jj,ii)     +    &
                                           ( bb ) * u(kk,jj,ii+1)   +    &
                                           ( cc ) * u(kk,jj+1,ii)   +    &
                                           ( dd ) * u(kk,jj+1,ii+1) ) /  &
                                         ( gg )

                               u_int_u = ( ( aa ) * u(kk+1,jj,ii)     +    &
                                           ( bb ) * u(kk+1,jj,ii+1)   +    &
                                           ( cc ) * u(kk+1,jj+1,ii)   +    &
                                           ( dd ) * u(kk+1,jj+1,ii+1) ) /  &
                                         ( gg )

                               u_int_b(k,j,i,1) = (1/dz(1)) * ((dz(1)-z)*u_int_l + z*u_int_u)

                            ELSE
                               u_int_b(k,j,i,1) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,1) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                      Velocity at element k of blade j of turbine i, interpolation of v:
                         ii = FLOOR( ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +         &
                                       fast_position_blade(k,j,i,1) - 0.5_wp * dx)  * ddx )                              ! wind speed in rotor area

                         jj = FLOOR( ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +         &
                                       fast_position_blade(k,j,i,2) ) * ddy )

                         kk = FLOOR( ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +         &
                                       fast_position_blade(k,j,i,3) + 0.5_wp * dz(1) ) / dz(1) )
!
!--                      Interpolation only takes place on a processor element
!--                      if all interpolation points are available on the processor
!--                      element
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) - ii * dx - 0.5_wp * dx                              ! wind speed in rotor area

                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy

                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1) + 0.5_wp * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               v_int_l = ( ( aa ) * v(kk,jj,ii)     +                              &
                                           ( bb ) * v(kk,jj,ii+1)   +                              &
                                           ( cc ) * v(kk,jj+1,ii)   +                              &
                                           ( dd ) * v(kk,jj+1,ii+1) ) /                            &
                                           ( gg )

                               v_int_u = ( ( aa ) * v(kk+1,jj,ii)     +                            &
                                           ( bb ) * v(kk+1,jj,ii+1)   +                            &
                                           ( cc ) * v(kk+1,jj+1,ii)   +                            &
                                           ( dd ) * v(kk+1,jj+1,ii+1) ) /                          &
                                         ( gg )

                               u_int_b(k,j,i,2) = (1/dz(1)) * ((dz(1)-z)*v_int_l + z*v_int_u)


                            ELSE
                               u_int_b(k,j,i,2) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,2) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )
!
!--                      Geschwindigkeit im Element k des Blattes k der Turbine i,
!--                      Interpolation von w:
                         ii = FLOOR( ( palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +         &
                                       fast_position_blade(k,j,i,1) - 0.5_wp * dx) * ddx )                              ! wind speed in rotor area

                         jj = FLOOR( ( palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +         &
                                       fast_position_blade(k,j,i,2) - 0.5_wp * dy) * ddy )

                         kk = FLOOR( ( palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +         &
                                       fast_position_blade(k,j,i,3) ) / dz(1) )
!
!--                      Auf einem Prozessorelement findet nur dann eine Interpolation
!--                      statt, wenn alle Stützpunkte der Interpolation auf dem
!--                      Prozessorelement zur Verfügung stehen
                         IF ( ( ii >= nxl ) .AND. ( ii <= nxr ) )  THEN
                            IF ( ( jj >= nys ) .AND. ( jj <= nyn ) )  THEN
                               x  = palm_tower_ref_pos_x(i) + fast_tower_ref_pos_x(i) +            &
                                    fast_position_blade(k,j,i,1) - ii * dx - 0.5_wp * dx                              ! wind speed in rotor area

                               y  = palm_tower_ref_pos_y(i) + fast_tower_ref_pos_y(i) +            &
                                    fast_position_blade(k,j,i,2) - jj * dy - 0.5_wp * dy

                               z  = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +            &
                                    fast_position_blade(k,j,i,3) - kk * dz(1)

                               aa = abs((dx-x)*(dy-y))
                               bb = abs((x)*(dy-y))
                               cc = abs((dx-x)*(y))
                               dd = x*y
                               gg = dx*dy

                               w_int_l = ( ( aa ) * w(kk,jj,ii)     +                              &
                                           ( bb ) * w(kk,jj,ii+1)   +                              &
                                           ( cc ) * w(kk,jj+1,ii)   +                              &
                                           ( dd ) * w(kk,jj+1,ii+1) ) /                            &
                                         ( gg )

                               w_int_u = ( ( aa ) * w(kk+1,jj,ii)     +                            &
                                           ( bb ) * w(kk+1,jj,ii+1)   +                            &
                                           ( cc ) * w(kk+1,jj+1,ii)   +                            &
                                           ( dd ) * w(kk+1,jj+1,ii+1) ) /                          &
                                         ( gg )

                               u_int_b(k,j,i,3) = (1/dz(1)) * ((dz(1)-z)*w_int_l + z*w_int_u)


                            ELSE
                               u_int_b(k,j,i,3) = 0.0_wp
                            ENDIF
                         ELSE
                            u_int_b(k,j,i,3) = 0.0_wp
                         ENDIF !( ( ii >= nxl ) .AND. ( ii <= nxr ) )

                      ENDDO !k = 1, fast_n_blade_elem(i)
                   ENDDO !j = 1, fast_n_blades(i)
                ENDDO !i = 1, nturbines


#if defined( __parallel )
!
!--             In preparation for exchanging the information about the velocities
!--             at the positions obtained from FAST, the sizes of the fields
!--             to be sent via MPI command are determined
                size_array_type_1 = nturbines*fast_n_blades_max*fast_n_blade_elem_max*3
!
!--             Exchange of information about the velocities at the positions
!--             received from FAST using MPI commands
!--             It is important that PE0 receives all information, since the
!--             communication with the FAST program later runs through this
!--             processor element
!--             As a precaution, before exchanging information between all
!--             processor elements, all processor elements are synchronized once again
                CALL MPI_BARRIER( comm2d, ierr )
!
!--             Exchange of the velocity components at the rotor hub
!--             The information is then available in the field palm_hub_center_vel
                CALL MPI_ALLREDUCE( u_int_h(1,1), palm_hub_center_vel(1,1), nturbines*3,           &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
!
!--             Exchange of the velocity components at the blade elements
!--             The information is then available in the field palm_vel_blade
                CALL MPI_ALLREDUCE( u_int_b(1,1,1,1), palm_vel_blade(1,1,1,1), size_array_type_1,  &
                                    MPI_REAL, MPI_SUM, comm2d, ierr )
!
!--             Resynchronization of the processor elements
                CALL MPI_BARRIER( comm2d, ierr )
#else
                palm_hub_center_vel = u_int_h
                palm_vel_blade = u_int_b
#endif

#ifdef __VERBOSE_MODE
print *, 'POST_COPY palm_hub_center_vel(1,:)=', palm_hub_center_vel(1,1), &
                                                   palm_hub_center_vel(1,2), &
                                                   palm_hub_center_vel(1,3), &

          ' u_int_h(1,:)=', u_int_h(1,1), u_int_h(1,2), u_int_h(1,3)                                                   

#endif

! DEBUGGING: BY SAMIR: 21:00 07.04.2025. Goal is to output blade velocities for investigation of power fluctuations

!                write(time_str, '(F0.15)') current_time_fast(1)                 
!               open(unit=71, file='blade_01_velocity_' // trim(adjustl(time_str)) // '.txt')
!                open(unit=72, file='blade_02_velocity_' // trim(adjustl(time_str)) // '.txt')
!               open(unit=73, file='blade_03_velocity_' // trim(adjustl(time_str)) // '.txt')

!                do i = 1, 45
!                    write(71, *) palm_vel_blade(i, 1, 1, 1),  palm_vel_blade(i, 1, 1, 2),  palm_vel_blade(i, 1, 1, 3)
!                     write(72, *) palm_vel_blade(i, 2, 1, 1),  palm_vel_blade(i, 2, 1, 2),  palm_vel_blade(i, 2, 1, 3)
!                   write(73, *) palm_vel_blade(i, 3, 1, 1),  palm_vel_blade(i, 3, 1, 2),  palm_vel_blade(i, 3, 1, 3)
!               end do

!                close(71)
!                close(72)
!                close(73)
! END DEBUGGING by Samir

! SoS  Added to look at first line of the new sector only ->
            END IF ! l == 1
! SoS  Added to look at first line of the new sector only <-

! SoS >     

               IF ( l == 1 ) THEN
                  IF ( myid == 0 ) THEN

!--                  S. Ouchene: on the first coupling step, build the NWC blade
!--                  geometry and airfoil polars from the data received from
!--                  OpenFAST. This cannot be done in f8c_init_arrays because the
!--                  TCP connection is not open there; the geometry arrives with
!--                  the connection-init response during this first f8c_actions
!--                  call, before this point. Built once; the routine also calls
!--                  init_nw_smearing_correction.
                     IF ( use_near_wake_smearing_correction  .AND.                                  &
                          .NOT. nwc_arrays_populated )  THEN
                        IF ( .NOT. nwc_geom_received )  THEN
                           message_string = 'use_near_wake_smearing_correction = .TRUE. but the ' //&
                              'blade geometry and airfoil polars were not received from OpenFAST. ' //&
                              'Check that OpenFAST runs AD15 (CompAero = 2) and the coupling is up.'
                           CALL message( 'f8c_actions', 'F8C0025', 2, 2, 0, 6, 0 )
                        END IF
                        CALL f8c_build_nwc_arrays
                     END IF

                     len_2d = n_airfoil_polar_rows
                     eps = reg_fac * dx

                     
!--                     S. Ouchene: Apply near-wake smearing corrections (Meyer-Forsting model) to 
!--                     the blade velocities interpolated from the PALM flow field.
!--                     velocity_correction() writes the corrected velocities into palm_vel_blade_b,
!--                     which is then copied back into palm_vel_blade before the values are sent
!--                     to FAST. Disabled by default (use_near_wake_smearing_correction = .FALSE.).
                        IF ( use_near_wake_smearing_correction )  THEN
                           CALL velocity_correction( len_2d,                                       &
                                                     fast_n_blade_elem_max,                        &
                                                     fast_position_blade,                          &
                                                     shaft_coordinates,                            &
                                                     palm_vel_blade,                               &
                                                     fast_n_blades_max,                            &
                                                     nturbines,                                    &
                                                     rotspeed,                                     &
                                                     blade_pitch,                                  &
                                                     eps,                                          &
                                                     dt_palm,                                      &
                                                     nw_trailed_vorticity_smoothing_length,        &
                                                     palm_vel_blade_b )
 
                          
                           palm_vel_blade = palm_vel_blade_b
                        ENDIF

                  END IF ! IF ( myid == 0 ) THEN
               END IF ! IF ( l == 1 ) THEN

!-- S. Ouchene (May 2026): wind-turbine NetCDF output. Disabled by default
!-- (output_interval = 0). Metadata bundle is populated once (SAVEd inside the helper).
!-- Call placed outside the IF (l == 1) block so FAST-cadence mode can fire on every sub-step;
!-- PALM-cadence mode internally skips when l /= 1.
               CALL save_wt_data( l )

! < SoS
!
!--             Since the temporary fields with the interpolated velocity values are no longer
!--             needed, they can now be deallocated
                DEALLOCATE( u_int_h, u_int_b )
!
!--             PALM-FAST interaction:
!--             The velocities interpolated from the velocity fields of PALM
!--             at the positions obtained from FAST are now transferred from
!--             PALM to FAST using a socket.
!
!--             Essentially two steps:
!--             Step 1: Pack the data into the format expected by FAST
!--             Step 2: Start the PALM client, which connects to the FAST server
!--                     and thus provides the data
!--             Since the exchange with FAST again only takes place via PE0,
!--             a case distinction is made below according to processor elements
                IF ( myid == 0 )  THEN
                   IF ( .NOT. terminate_run )  THEN

                      IF ( debug_output )  THEN
                         WRITE(9,*) 'Sending velocities to FAST server(s)...'
                         FLUSH(9)
                      ENDIF

#if defined( __fastv8 )
                      iflag = commclient(COM_TARGET_OPENFAST, COM_SEND_VELOCITIES)
#else
                      iflag = -1
#endif

                      IF ( iflag < 0 )  THEN
                         WRITE( message_string, * )  'Unable to send velocities to FAST server(s)',&
                                '. Please check if any of the FAST servers finished the simulation.'
                         CALL message( 'f8c_actions', 'F8C0019', 0, 2, 0, 6, 0 )
                         terminate_run = .TRUE.
                      ELSEIF ( iflag == 1 )  THEN
                         WRITE( message_string, * )                                                &
                                'At least one FAST server finished the simulation'
                         CALL message( 'f8c_actions', 'F8C0020', 0, 2, 0, 6, 0 )
                         terminate_run = .TRUE.
                      ELSEIF ( iflag == 0 )  THEN
                         IF ( debug_output )  THEN
                            WRITE(9,*) 'Velocities successfully broadcasted'
                            FLUSH(9)
                         ENDIF
                      ENDIF !( iflag < 0 )

                   ENDIF !(.NOT. terminate_run )
                ENDIF !( myid == 0 )

#if defined( __parallel )
                   CALL MPI_BARRIER( comm2d, ierr )
                   CALL MPI_BCAST( terminate_run, 1, MPI_LOGICAL, 0, comm2d, ierr )
                   CALL MPI_BARRIER( comm2d, ierr )
#endif
                   IF ( terminate_run )  THEN
                      RETURN
                   ENDIF
!
!--             Restoring the positions and forces to have the starting position within
!--             the sector and the forces of the line in the middle of the sector.
                IF ( myid == 0 )  THEN
                   IF ( l == 1 )  THEN

                      DO i = 1, nturbines

                         fast_hub_center_pos_old(i,1) = fast_hub_center_pos(i,1)
                         fast_hub_center_pos_old(i,2) = fast_hub_center_pos(i,2)
                         fast_hub_center_pos_old(i,3) = fast_hub_center_pos(i,3)

                         DO m = 1, fast_n_blades(i)
                            DO n = 1, fast_n_blade_elem(i)

                               fast_position_blade_old(n,m,i,1) = fast_position_blade(n,m,i,1)
                               fast_position_blade_old(n,m,i,2) = fast_position_blade(n,m,i,2)
                               fast_position_blade_old(n,m,i,3) = fast_position_blade(n,m,i,3)

                            ENDDO
                         ENDDO

                      ENDDO

                   ENDIF

                   IF ( l == ceiling( (dt_palm/dt_fast ) / 2 ) )  THEN

                      DO i = 1, nturbines
                         DO m = 1, fast_n_blades(i)
                            DO n = 1, fast_n_blade_elem(i)
                               fast_force_blade_old(n,m,i,1) = fast_force_blade(n,m,i,1)
                               fast_force_blade_old(n,m,i,2) = fast_force_blade(n,m,i,2)
                               fast_force_blade_old(n,m,i,3) = fast_force_blade(n,m,i,3)
                            ENDDO
                         ENDDO
                      ENDDO

                   ENDIF

                ENDIF !(myid == 0)

#if defined( __parallel )
!
!--             Wegen unterschiedlicher Aktivitaeten der einzelnen Prozessorelemente
!--             im Vorfeld werden die Prozessorelemente an dieser Stelle
!--             sicherheitshalber synchronisiert
                CALL MPI_BARRIER( comm2d, ierr )
!
!--             Nach der Kommunikation mit FAST sind Informationen ueber
!--             Kraefte und Positionen zunaechst nur auf PE0 bekannt. Mittels
!--             MPI_BCAST sollen die Informationen auch auf den anderen
!--             Prozessorelementen bekanntgemacht werden
!--             Vorher ist noch die Groesse der mittels MPI_BCAST
!--             auszutauschenden Felder zu bestimmen.
                size_array_type_1 = nturbines * fast_n_blades_max * fast_n_blade_elem_max * 3
!
!--             Broadcasting von fast_hub_center_pos (nturbines,3)
!--             (Position der Rotornabe (x-, y- und z-Komponente))
                CALL MPI_BCAST( fast_hub_center_pos(1,1), nturbines*3,                             &
                                MPI_REAL, 0, comm2d, ierr )
!
!--             Broadcasting von fast_position_blade
!--             (fast_n_blade_elem_max,fast_n_blades_max,nturbines,3)
!--             (Positionen der Blattelemente (x-, y- und z-Komponenten))
                CALL MPI_BCAST( fast_position_blade(1,1,1,1), size_array_type_1,                   &
                                MPI_REAL, 0, comm2d, ierr )
!
!--             Broadcasting von fast_force_blade
!--             (fast_n_blade_elem_max,fast_n_blades_max,nturbines,3)
!--             (Kraefte an den Blattelementen (x-, y- und z-Komponenten))
                CALL MPI_BCAST( fast_force_blade(1,1,1,1), size_array_type_1,                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_position_blade_old(1,1,1,1), size_array_type_1,               &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( fast_force_blade_old(1,1,1,1), size_array_type_1,                  &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BCAST( current_time_fast, nturbines,                                      &
                                MPI_REAL, 0, comm2d, ierr )
                CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDDO ! DO WHILE ( current_time_fast(1) <= simulated_time + dt_palm )
!
!--          Calculating the angle of the sector and the rotation matrix,
!--          to rotate the positions starting from the first line in the sector
             IF ( myid == 0 )  THEN

                DO  i = 1, nturbines

                   IF ( palm_dt_3d )  THEN
                      sector_angle(i) = rotspeed(i) * dt_palm
                      n_sector(i)     =   n_dt_fast
                      alpha_rot(i)    =  sector_angle(i) / n_sector(i)
                   ELSE
                      sector_angle_deg(i) = 120.0_wp
                      n_sector(i)         = MAX(n_dt_fast, 1_iwp)
                      alpha_rot(i)        = ( sector_angle_deg(i) / n_sector(i) ) * pi / 180.0_wp
                   ENDIF
!
!--                Building the rotation matrix:
                   shaft_x = shaft_coordinates(i, 1)
                   shaft_y = shaft_coordinates(i, 2)
                   shaft_z = shaft_coordinates(i, 3)

                   cos_alpha            = COS(alpha_rot(i))
                   sin_alpha            = SIN(alpha_rot(i))
                   one_minus_cos_alpha  = 1.0_wp - cos_alpha

                   r_n(i,1,:) = (/                                                                 &
                      shaft_x**2 * one_minus_cos_alpha + cos_alpha,                                &
                      shaft_x * shaft_y * one_minus_cos_alpha - shaft_z * sin_alpha,                &
                      shaft_x * shaft_z * one_minus_cos_alpha + shaft_y * sin_alpha                &
                                /)

                   r_n(i,2,:) = (/                                                                 &
                      shaft_y * shaft_x * one_minus_cos_alpha + shaft_z * sin_alpha,               &
                      shaft_y**2 * one_minus_cos_alpha + cos_alpha,                                &
                      shaft_y * shaft_z * one_minus_cos_alpha - shaft_x * sin_alpha                &
                                /)

                   r_n(i,3,:) = (/                                                                 &
                      shaft_z * shaft_x * one_minus_cos_alpha - shaft_y * sin_alpha,               &
                      shaft_z * shaft_y * one_minus_cos_alpha + shaft_x * sin_alpha,               &
                      shaft_z**2 * one_minus_cos_alpha + cos_alpha                                 &
                                /)

                ENDDO !DO  i = 1, nturbines

             ENDIF !( myid == 0 )

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
             CALL MPI_BCAST( n_sector, nturbines, MPI_INTEGER, 0, comm2d, ierr )
             CALL MPI_BCAST( r_n, nturbines*9, MPI_REAL, 0, comm2d, ierr )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             n_sector_max = MAXVAL( n_sector(:) )

             ALLOCATE( x_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )
             ALLOCATE( y_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )
             ALLOCATE( z_prime(1:fast_n_blade_elem_max,                                            &
                               1:fast_n_blades_max,                                                &
                               1:nturbines,                                                        &
                               1:n_sector_max) )

             x_prime(:,:,:,:) = 0.0_wp
             y_prime(:,:,:,:) = 0.0_wp
             z_prime(:,:,:,:) = 0.0_wp

             ALLOCATE( x_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )
             ALLOCATE( y_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )
             ALLOCATE( z_prime_sec(1:fast_n_blade_elem_max,                                        &
                                   1:n_sector_max*fast_n_blades_max,                               &
                                   1:nturbines) )

             x_prime_sec(:,:,:) = 0.0_wp
             y_prime_sec(:,:,:) = 0.0_wp
             z_prime_sec(:,:,:) = 0.0_wp
!
!--          Calculating the positions of the points on the line samples distributed over the three blade sectors
             DO  i = 1, nturbines
                DO  j = 1, fast_n_blades(i)
                   DO  k = 1, fast_n_blade_elem(i)

                      x = fast_tower_ref_pos_x(i) +                                                &
                          fast_position_blade_old(k,j,i,1)
                      y = fast_tower_ref_pos_y(i) +                                                &
                          fast_position_blade_old(k,j,i,2)
                      z = palm_tower_ref_pos_z(i) + fast_tower_ref_pos_z(i) +                      &
                          fast_position_blade_old(k,j,i,3) - shaft_height_fast(i)
!
!--                   Multiplication of the points with the rotation matrix -> positions of the lines within the sector
                      x_prime(k,j,i,1) = x
                      y_prime(k,j,i,1) = y
                      z_prime(k,j,i,1) = z
                      DO m = 2, n_sector(i)

                         tmp_res = MATMUL( r_n(i,:,:),                                             &
                            (/ x_prime(k,j,i,m-1), y_prime(k,j,i,m-1), z_prime(k,j,i,m-1) /)       &
                                         )

                         x_prime(k,j,i,m) = tmp_res(1)
                         y_prime(k,j,i,m) = tmp_res(2)
                         z_prime(k,j,i,m) = tmp_res(3)

                      ENDDO !DO m = 2, n_sector(i)

                   ENDDO !DO  k = 1, fast_n_blade_elem(i)
                ENDDO !DO  j = 1, fast_n_blades(i)
!
!--             Positions in the PALM grid
                x_prime(:,:,i,:) = x_prime(:,:,i,:) + palm_tower_ref_pos_x(i)
                y_prime(:,:,i,:) = y_prime(:,:,i,:) + palm_tower_ref_pos_y(i)
                z_prime(:,:,i,:) = z_prime(:,:,i,:) + shaft_height_fast(i)

             ENDDO !DO i = 1, nturbines

             DO  i = 1, nturbines
                DO  j = 1, fast_n_blades(i)
                   DO  k = 1, fast_n_blade_elem(i)
                      DO m = 1, n_sector(i)
                         x_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = x_prime(k,j,i,m)
                         y_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = y_prime(k,j,i,m)
                         z_prime_sec(k,(m-1)*fast_n_blades(i) + j,i) = z_prime(k,j,i,m)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO !DO i = 1, nturbines

             !-Verschmieren und Aufsummieren der Kraefte

             !-Einfuehren eines Regulierungskernels
             eps_kernel = reg_fac * dx
             eps_kernel_2 = -reg_fac * reg_fac * dx * dx

             nenner = (eps_kernel**3.0) * (pi**(3.0/2.0))
!
!--          Informationen ueber die zusaetzlichen durch die Windenergieanlage verursachten
!--          Kraefte werden an jedem Gitterpunkt des Modellgebiets benoetigt
!--          In einer Schleife ueber alle Gitterpunkte des Teilmodellgebiets wird im Folgenden
!--          die Verschmierung der Kraefte durchgefuehrt, so dass jedem Gitterpunkt eine
!--          zusaetzliche Kraft zugeordnet werden kann
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      force_x(k,j,i) = 0.0_wp
                      force_y(k,j,i) = 0.0_wp
                      force_z(k,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDDO

             DO  o = 1, nturbines
!
!--             check if coordinates are even in the limit of the processor
                IF ( .NOT.  ((nxlg >= fboxcorners(o, 1, 2))   .OR.                                  &
                            (nxrg <= fboxcorners(o, 1, 1)))  )  THEN

                   IF ( .NOT. (( nysg >= fboxcorners(o, 2, 2) )  .OR.                               &
                              ( nyng <= fboxcorners(o, 2, 1) )) )  THEN

                      IF ( .NOT. (( nzb >= fboxcorners(o, 3, 2) )  .OR.                             &
                                 ( nzt <= fboxcorners(o, 3, 1) )) )  THEN

                         tower_ref_pos_x(o) = palm_tower_ref_pos_x(o) + fast_tower_ref_pos_x(o)
                         tower_ref_pos_y(o) = palm_tower_ref_pos_y(o) + fast_tower_ref_pos_y(o)
                         tower_ref_pos_z(o) = palm_tower_ref_pos_z(o) + fast_tower_ref_pos_z(o)

                         DO  i = nxlg, nxrg
                            IF ( ( i >= fboxcorners(o, 1, 1) )  .AND.                              &
                                 ( i <= fboxcorners(o, 1, 2) ) )  THEN

                               DO  j = nysg, nyng
                                  IF ( ( j >= fboxcorners(o, 2, 1) )  .AND.                        &
                                       ( j <= fboxcorners(o, 2, 2) ) )  THEN

                                     DO  k = nzb, nzt
                                        IF ( ( k >= fboxcorners(o, 3, 1) )  .AND.                  &
                                             ( k <= fboxcorners(o, 3, 2) ) )  THEN

                                           DO  m = 1, fast_n_blades(o)
                                              DO  n = 1, fast_n_blade_elem(o)
                                                 DO n_s = 1, n_sector(o)
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des u-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy - 0.5_wp * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)            &
                                                                            + 0.5_wp * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )

                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in x-Richtung auf dem u-grid
                                                    IF ( mod(m,3) == 0 )  THEN
                                                       l = 3
                                                    ELSE
                                                       l = mod(m,3)
                                                    ENDIF

                                                    force_x(k,j,i) = force_x(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,1) &
                                                                     * kernel_value
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des v-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx - 0.5_wp * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)            &
                                                                            + 0.5_wp * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )


                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in y-Richtung auf dem v-grid
                                                    force_y(k,j,i) = force_y(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,2) &
                                                                     * kernel_value
!
!--                                                 Bestimmung der Distanz zwischen Blattelement und dem Gitterpunkt des w-grids
                                                    x  = x_prime(n,m,o,n_s) - i * dx - 0.5_wp * dx
                                                    y  = y_prime(n,m,o,n_s) - j * dy - 0.5_wp * dy
                                                    z  = z_prime(n,m,o,n_s) - k * dz(1)

                                                    distance_hub = x**2_iwp + y**2_iwp + z**2_iwp

                                                    expargument = INT(                             &
                                                       - distance_hub / eps_kernel_2 * 1000.0_wp   &
                                                                     )

                                                    IF ( expargument < 10000001_iwp )  THEN
                                                       kernel_value = expvalue(expargument)/nenner
                                                    ELSE
                                                       kernel_value = 0.0_wp
                                                    ENDIF
!
!--                                                 Berechnung der Kraft in z-Richtung auf dem w-grid
                                                    force_z(k,j,i) = force_z(k,j,i) -              &
                                                                     ( 1/float(n_sector(o)) ) *    &
                                                                     fast_force_blade_old(n,m,o,3) &
                                                                     * kernel_value
                                                 ENDDO
                                              ENDDO
                                           ENDDO

!--                                        S. Ouchene: compute signed tower drag components (oppose local velocity component)
!--                                        but also notice that currently the tower is also clipped by fboxcorners which shouldn't
!--                                        be the case
                                           thrust_tower_x = 0.5_wp * turb_C_d_tow(o) *                           &
                                                            tower_area_x_4d(o,k,j,i) / (dx * dy * dz(1)) *       &
                                                            ABS(u(k,j,i)) * u(k,j,i)

                                           thrust_tower_y = 0.5_wp * turb_C_d_tow(o) *                           &
                                                            tower_area_y_4d(o,k,j,i) / (dx * dy * dz(1) ) *      &
                                                            ABS( v(k,j,i) ) * v(k,j,i)

                                           ! subtract tower drag               
                                           force_x(k,j,i) = force_x(k,j,i) - thrust_tower_x
                                           force_y(k,j,i) = force_y(k,j,i) - thrust_tower_y

                                        ENDIF !((k >= fboxcorners(o, 3, 1)) .AND. (k <= fboxcorners(o, 3, 2)))
                                     ENDDO !DO k = nzb, nzt

                                  ENDIF !((j >= fboxcorners(o, 2, 1)) .AND. (j <= fboxcorners(o, 2, 2)))
                               ENDDO !DO j = nysg, nyng

                            ENDIF !((i >= fboxcorners(o, 1, 1)) .AND. (i <= fboxcorners(o, 1, 2)))
                         ENDDO !DO  i = nxlg, nxrg

                      ENDIF !(.NOT.((nzb >= fboxcorners(o, 3, 2)) .OR. (nzt <= fboxcorners(o, 3, 1))))
                   ENDIF !(.NOT.((nysg >= fboxcorners(o, 2, 2)) .OR. (nyng <= fboxcorners(o, 2, 1))))
                ENDIF !(.NOT.((nxlg >= fboxcorners(o, 1, 2)) .OR. (nxrg <= fboxcorners(o, 1, 1))))

             ENDDO !DO o = 1, nturbines

#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif

             DEALLOCATE( x_prime, y_prime, z_prime, stat=ierr )
             DEALLOCATE( x_prime_sec, y_prime_sec, z_prime_sec, stat=ierr )

          ENDIF !( simulated_time >= time_turbine_on )
          
! SoS was in CASE('after_timestep') before... does it work here?
       IF ( myid == 0 ) THEN 
          IF (first_restart_with_wt) THEN
             first_restart_with_wt = .FALSE.
          ENDIF
       ENDIF
! SoS


       CASE ( 'after_integration' )
!
!--       Total thrust and torque components calculated above:
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                   thrx(k,j,i) = force_x(k,j,i)
                   tory(k,j,i) = force_y(k,j,i)
                   torz(k,j,i) = force_z(k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'u-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_x(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'v-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_v_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_y(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'w-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt !nzb_w_inner(j,i)+1, nzt
                      tend(k,j,i) = tend(k,j,i) + force_z(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

    CALL cpu_log( log_point(24), 'f8c_actions', 'stop' )

 END SUBROUTINE f8c_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_actions_ij( i, j, location )


    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_u_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_x(k,j,i)
             ENDDO
          ENDIF

       CASE ( 'v-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_v_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_y(k,j,i)
             ENDDO
          ENDIF

       CASE ( 'w-tendency' )

          IF ( simulated_time - palm_time_at_begin >= time_turbine_on )  THEN
             DO  k = nzb, nzt !nzb_w_inner(j,i)+1, nzt
                tend(k,j,i) = tend(k,j,i) + force_z(k,j,i)
             ENDDO
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE f8c_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) and sets the grid on which it is defined. Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )


    CHARACTER(LEN=*), INTENT(INOUT) ::  grid      !< name of vertical grid
    CHARACTER(LEN=*), INTENT(IN)    ::  mode      !< either 'xy', 'xz' or 'yz'
    CHARACTER(LEN=*), INTENT(IN)    ::  variable  !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av      !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp), INTENT(IN) ::  nzb_do  !< lower limit of the domain (usually nzb)
    INTEGER(iwp), INTENT(IN) ::  nzt_do  !< upper limit of the domain (usually nzt+1)

    INTEGER(iwp) ::  i  !< grid index along x-direction
    INTEGER(iwp) ::  j  !< grid index along y-direction
    INTEGER(iwp) ::  k  !< grid index along z-direction

    LOGICAL, INTENT(INOUT) ::  found  !<
    LOGICAL, INTENT(INOUT) ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp  .OR.  two_d )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx_xy', 'thrx_xz', 'thrx_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = thrx(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'tory_xy', 'tory_xz', 'tory_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = tory(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'torz_xy', 'torz_xz', 'torz_yz' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   local_pf(i,j,k) = torz(k,j,i)
                ENDDO
             ENDDO
          ENDDO
          IF ( mode == 'xy' ) grid = 'zu'


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE f8c_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE f8c_data_output_3d( av, variable, found, local_pf, resorted, flag_nr, nzb_do, nzt_do )


    CHARACTER(LEN=*), INTENT(IN) ::  variable  !<  name of variable

    INTEGER(iwp), INTENT(IN)    ::  av       !< flag for (non-)average output
    INTEGER(iwp), INTENT(INOUT) ::  flag_nr  !< number of masking flag, 0 = scalar, 1 = u, 2 = v, 3 = w
    INTEGER(iwp), INTENT(IN)    ::  nzb_do   !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN)    ::  nzt_do   !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL, INTENT(INOUT) ::  found     !< flag if output variable is found
    LOGICAL, INTENT(INOUT) ::  resorted  !< flag if output is resorted

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf    !< local array
                                                        !< to which output data is resorted to

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<


    found    = .TRUE.
    resorted = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'thrx' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                flag_nr = 0
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( thrx(k,j,i) + thrx(k,j,i+1) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                flag_nr = 1
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = thrx(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF

       CASE ( 'tory' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                flag_nr = 0
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( tory(k,j,i) + tory(k,j+1,i) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                flag_nr = 2
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = tory(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF

       CASE ( 'torz' )
          IF ( av == 0 )  THEN

             IF ( interpolate_to_grid_center )  THEN
                flag_nr = 0
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = 0.5_wp * ( torz(k,j,i) + torz(k-1,j,i) )
                     ENDDO
                   ENDDO
                ENDDO
             ELSE
                flag_nr = 3
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         local_pf(i,j,k) = torz(k,j,i)
                     ENDDO
                   ENDDO
                ENDDO
             ENDIF

          ENDIF

       CASE DEFAULT
          found    = .FALSE.
          resorted = .FALSE.

    END SELECT


 END SUBROUTINE f8c_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Wrapper around write_wt_data that populates the metadata bundle (wt_output_meta_t)
! exactly once - on the first call - using a SAVEd variable, then reuses it on every subsequent
! call. The metadata is read-only after initialisation so the one-shot capture is correct.
!
! Called from inside the FAST sub-step DO WHILE in f8c_actions. Cheap on the non-first calls -
! just a logical check and a forward of the SAVEd meta.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE save_wt_data( fast_sub_iter )

    INTEGER(iwp), INTENT(IN) ::  fast_sub_iter

    TYPE(wt_output_meta_t), SAVE ::  wt_meta
    LOGICAL,                SAVE ::  wt_meta_populated = .FALSE.

    IF ( .NOT. wt_meta_populated )  THEN
       wt_meta%output_interval                       = output_interval
       wt_meta%output_use_fast_timestep              = output_use_fast_timestep
       wt_meta%output_sync_per_write                 = output_sync_per_write
       wt_meta%nturbines                             = nturbines
       wt_meta%fast_n_blades_max                     = fast_n_blades_max
       wt_meta%fast_n_blade_elem_max                 = fast_n_blade_elem_max
       wt_meta%dt_palm                               = dt_palm
       wt_meta%dt_fast                               = dt_fast
       wt_meta%end_simulated_time                    = end_time
       wt_meta%use_near_wake_smearing_correction     = use_near_wake_smearing_correction
       wt_meta%nw_trailed_vorticity_smoothing_length = nw_trailed_vorticity_smoothing_length
       wt_meta%is_restart                            = ( TRIM( initializing_actions )              &
                                                          == 'read_restart_data' )

       wt_meta%dx                                    = dx
       wt_meta%dy                                    = dy
       wt_meta%dz                                    = dz(1)        !< PALM supports stretched grids (dz is rank-1); record the base value
       wt_meta%nx                                    = nx
       wt_meta%ny                                    = ny
       wt_meta%nz                                    = nz
       wt_meta%npex                                  = npex
       wt_meta%npey                                  = npey

       wt_meta%time_turbine_on                       = time_turbine_on
       wt_meta%reg_fac                               = reg_fac
       wt_meta%fbox_fac                              = fbox_fac
       wt_meta%n_airfoil_polar_rows                  = n_airfoil_polar_rows

       wt_meta%palm_tower_ref_pos_x = palm_tower_ref_pos_x(1:nturbines)
       wt_meta%palm_tower_ref_pos_y = palm_tower_ref_pos_y(1:nturbines)
       wt_meta%palm_tower_ref_pos_z = palm_tower_ref_pos_z(1:nturbines)
       wt_meta%dtow                 = dtow(1:nturbines)
       wt_meta%htow                 = htow(1:nturbines)
       wt_meta%turb_C_d_tow         = turb_C_d_tow(1:nturbines)

       wt_meta%fast_host_ports_joined    = join_semicolon( fast_host_port(1:nturbines) )

       wt_meta_populated = .TRUE.
    END IF

    CALL write_wt_data(                                                                    &
            meta                 = wt_meta,                                                        &
            fast_sub_iter        = fast_sub_iter,                                                  &
            current_time_fast    = current_time_fast,                                              &
            fast_position_blade  = fast_position_blade,                                            &
            fast_force_blade     = fast_force_blade,                                               &
            palm_vel_blade       = palm_vel_blade,                                                 &
            rotspeed             = rotspeed,                                                       &
            blade_pitch          = blade_pitch,                                                    &
            shaft_coordinates    = shaft_coordinates )

 END SUBROUTINE save_wt_data


!--------------------------------------------------------------------------------------------------!
! Join an array of CHARACTER values with ';' separators, skipping empty entries. Returned as a
! deferred-length allocatable so the result has exactly the right length.
!--------------------------------------------------------------------------------------------------!
 FUNCTION join_semicolon( strs ) RESULT( joined )
    CHARACTER(LEN=*),    DIMENSION(:), INTENT(IN) ::  strs
    CHARACTER(LEN=:),    ALLOCATABLE              ::  joined
    INTEGER                                       ::  i

    joined = ''
    DO i = 1, SIZE( strs )
       IF ( LEN_TRIM( strs(i) ) == 0 )  CYCLE
       IF ( LEN( joined ) > 0 )  THEN
          joined = joined // ';' // TRIM( strs(i) )
       ELSE
          joined = TRIM( strs(i) )
       END IF
    END DO

 END FUNCTION join_semicolon


 END MODULE fastv8_coupler_mod
