!--------------------------------------------------------------------------------------------------!
! Wind-turbine NetCDF output for the FASTv8-PALM coupling.
!
! Writes blade-element diagnostics (positions, forces, velocities) and per-turbine state
! (rotational speed, shaft coordinates, FAST time) to a self-describing NetCDF file that
! lands in PALM's per-run output directory. Off by default.
!
! Cadence and configuration are bundled in wt_output_meta_t, populated once by the caller
! (from the fastv8_coupler_parameters namelist plus PALM control parameters) and passed in
! on every write. Three knobs:
!   - output_interval         INTEGER  0 = disabled. N > 0 = write every N-th step.
!   - output_use_fast_timestep LOGICAL .FALSE. = step counter ticks once per PALM step (l == 1).
!                                      .TRUE.  = step counter ticks once per FAST sub-step.
!   - output_sync_per_write   LOGICAL .TRUE.  = fsync after every write (crash-safe, slow on
!                                      parallel filesystems).
!                                      .FALSE. = let the OS buffer; flushed on close.
!
! Output file: <run_identifier>_wt_output_<NNNNNNNN>.nc, opened in PALM's working directory.
! NNNNNNNN is INT(simulated_time) at file open, zero-padded to 8 digits, so each restart
! segment lands in its own file with no append/overwrite hazard.
!
! Only rank 0 ever touches the NetCDF library. The wake correction and the data this
! module reads already live on rank 0.
!--------------------------------------------------------------------------------------------------!
MODULE wt_output_mod

   USE kinds,                                                                                      &
       ONLY:  iwp, wp

   USE netcdf

   USE control_parameters,                                                                         &
       ONLY:  run_identifier, simulated_time

   USE pegrid,                                                                                     &
       ONLY:  myid

!-- Per-element pre/post-correction diagnostics produced by nw_smearing_correction.
!-- All 15 arrays share the (turbine, blade, element) layout and are populated on
!-- every nw_smearing_correction call. Only read when meta%use_near_wake_smearing_correction
!-- is .TRUE.; unallocated otherwise.
   USE nw_smearing_correction_module,                                                              &
       ONLY:  nwc_diag_vn_pre,    nwc_diag_vt_pre,                                                 &
              nwc_diag_vn_post,   nwc_diag_vt_post,                                                &
              nwc_diag_phi_pre,   nwc_diag_phi_post,                                               &
              nwc_diag_alpha_pre, nwc_diag_alpha_post,                                             &
              nwc_diag_cl_pre,    nwc_diag_cl_post,                                                &
              nwc_diag_cd_pre,    nwc_diag_cd_post,                                                &
              nwc_diag_nw_vn,     nwc_diag_nw_vt,                                                  &
              nwc_diag_nw_W,                                                                      &
              pi

   IMPLICIT NONE

   PRIVATE


!--------------------------------------------------------------------------------------------------!
! Metadata bundle. Populated once by the caller (or at every call - cheap), passed in on every
! invocation of write_wt_data. Holds everything the module needs to: gate the cadence,
! open the file lazily, and emit the global-attribute block at file-open time. Components are
! grouped by purpose to keep the keyword-argument construction readable.
!--------------------------------------------------------------------------------------------------!
   TYPE, PUBLIC :: wt_output_meta_t
      ! Cadence (from the fastv8_coupler_parameters namelist)
      INTEGER(iwp) :: output_interval                              = 0_iwp
      LOGICAL      :: output_use_fast_timestep                     = .FALSE.
      LOGICAL      :: output_sync_per_write                        = .TRUE.

      ! Dimensions (from the fastv8_coupler_parameters namelist)
      INTEGER(iwp) :: nturbines                                    = 0_iwp
      INTEGER(iwp) :: fast_n_blades_max                            = 0_iwp
      INTEGER(iwp) :: fast_n_blade_elem_max                        = 0_iwp

      ! Simulation context captured for the global-attribute block
      REAL(wp)     :: dt_palm                                      = 0.0_wp
      REAL(wp)     :: dt_fast                                      = 0.0_wp
      REAL(wp)     :: end_simulated_time                           = 0.0_wp
      LOGICAL      :: use_near_wake_smearing_correction            = .FALSE.
      REAL(wp)     :: nw_trailed_vorticity_smoothing_length        = 0.0_wp
      LOGICAL      :: is_restart                                   = .FALSE.

      ! PALM grid (scalars from grid_variables / control_parameters / indices / pegrid)
      REAL(wp)     :: dx                                           = 0.0_wp
      REAL(wp)     :: dy                                           = 0.0_wp
      REAL(wp)     :: dz                                           = 0.0_wp
      INTEGER(iwp) :: nx                                           = 0_iwp
      INTEGER(iwp) :: ny                                           = 0_iwp
      INTEGER(iwp) :: nz                                           = 0_iwp
      INTEGER(iwp) :: npex                                         = 0_iwp
      INTEGER(iwp) :: npey                                         = 0_iwp

      ! Wind-turbine setup scalars (from fastv8_coupler_parameters)
      REAL(wp)     :: time_turbine_on                              = 0.0_wp
      REAL(wp)     :: reg_fac                                      = 0.0_wp
      REAL(wp)     :: fbox_fac                                     = 0.0_wp
      INTEGER(iwp) :: n_airfoil_polar_rows                         = 0_iwp

      ! Per-turbine arrays, length nturbines. Allocated by the caller (auto-allocation on
      ! assignment is fine: `meta%dtow = dtow(1:nturbines)`).
      REAL(wp), ALLOCATABLE :: palm_tower_ref_pos_x(:)
      REAL(wp), ALLOCATABLE :: palm_tower_ref_pos_y(:)
      REAL(wp), ALLOCATABLE :: palm_tower_ref_pos_z(:)
      REAL(wp), ALLOCATABLE :: dtow(:)
      REAL(wp), ALLOCATABLE :: htow(:)
      REAL(wp), ALLOCATABLE :: turb_C_d_tow(:)

      ! Semicolon-joined strings (CHARACTER attributes are scalar, so we flatten string arrays).
      CHARACTER(LEN=:), ALLOCATABLE :: fast_host_ports_joined
   END TYPE wt_output_meta_t


!--------------------------------------------------------------------------------------------------!
! Internal runtime state for the NetCDF file. Opaque to callers.
!--------------------------------------------------------------------------------------------------!
   TYPE :: wt_nc_state_t
      LOGICAL         ::  file_open    = .FALSE.
      INTEGER         ::  ncid         = -1                    !< NetCDF file handle
      INTEGER(iwp)    ::  step_counter = 0_iwp                 !< ticks per eligible step (PALM or FAST)
      INTEGER(iwp)    ::  time_index   = 0_iwp                 !< next slot along the unlimited dim
      INTEGER(KIND=8) ::  wall_open_count      = 0_8           !< SYSTEM_CLOCK count at file open
      INTEGER(KIND=8) ::  wall_open_count_rate = 1_8           !< SYSTEM_CLOCK count_rate at file open
      INTEGER         ::  dim_time, dim_turbine, dim_blade, dim_element, dim_xyz
      INTEGER         ::  var_time, var_time_fast, var_palm_sub_iter
      INTEGER         ::  var_position, var_force, var_palm_vel
      INTEGER         ::  var_rotspeed, var_shaft_coord
      INTEGER         ::  var_blade_pitch   !< per-blade pitch, defined only with the NWC
      ! NWC pre-/post-correction diagnostics (defined only when meta%use_near_wake_smearing_correction)
      INTEGER         ::  var_vn_pre,    var_vt_pre
      INTEGER         ::  var_vn_post,   var_vt_post
      INTEGER         ::  var_phi_pre,   var_phi_post
      INTEGER         ::  var_alpha_pre, var_alpha_post
      INTEGER         ::  var_cl_pre,    var_cl_post
      INTEGER         ::  var_cd_pre,    var_cd_post
      INTEGER         ::  var_nw_vn,     var_nw_vt,    var_nw_W
   END TYPE wt_nc_state_t

   TYPE(wt_nc_state_t), SAVE ::  st


   PUBLIC ::  wt_output_warn_if_slow, write_wt_data, wt_output_finalize


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Print a one-shot warning at startup (on rank 0 only) when wt-output is enabled with the
! crash-safe-but-slow per-write fsync. Called from f8c_check_parameters after the namelist
! has been read, so it fires once per run. Cadence and sync flags are passed as scalars
! because at parin time the dt_palm / dt_fast / end_time values inside wt_output_meta_t are
! not yet populated.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wt_output_warn_if_slow( output_interval, output_sync_per_write )

    INTEGER(iwp), INTENT(IN) ::  output_interval
    LOGICAL,      INTENT(IN) ::  output_sync_per_write

    IF ( myid /= 0 )                    RETURN
    IF ( output_interval <= 0 )         RETURN
    IF ( .NOT. output_sync_per_write )  RETURN

    WRITE(*,'(A)')  ''
    WRITE(*,'(A)')  '+------------------------------------------------------------------+'
    WRITE(*,'(A)')  '| [wt_output_mod] WARNING                                          |'
    WRITE(*,'(A)')  '|                                                                  |'
    WRITE(*,'(A)')  '| Wind-turbine NetCDF output is enabled with                       |'
    WRITE(*,'(A)')  '|     output_sync_per_write = .TRUE.                               |'
    WRITE(*,'(A)')  '|                                                                  |'
    WRITE(*,'(A)')  '| This may slow down the simulation. Set                           |'
    WRITE(*,'(A)')  '|     output_sync_per_write = .FALSE.                              |'
    WRITE(*,'(A)')  '| in fastv8_coupler_parameters if you do not need the per-write    |'
    WRITE(*,'(A)')  '| crash safety (the file is still flushed on close).               |'
    WRITE(*,'(A)')  '+------------------------------------------------------------------+'
    WRITE(*,'(A)')  ''

 END SUBROUTINE wt_output_warn_if_slow


!--------------------------------------------------------------------------------------------------!
! Call site: invoked from fastv8_coupler_mod once per FAST sub-step. Decides whether the
! current step is an output step and, if so, writes one time slice. Lazy-initialises the
! file on the first eligible call.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE write_wt_data( meta, fast_sub_iter,                                            &
                                   current_time_fast, fast_position_blade, fast_force_blade,      &
                                   palm_vel_blade, rotspeed, blade_pitch, shaft_coordinates )

    TYPE(wt_output_meta_t), INTENT(IN) ::  meta
    INTEGER(iwp),           INTENT(IN) ::  fast_sub_iter   !< value of 'l' at the call site

    REAL(wp), DIMENSION(:),         INTENT(IN) ::  current_time_fast, rotspeed
    REAL(wp), DIMENSION(:,:),       INTENT(IN) ::  blade_pitch   !< (turbine, blade), radians
    REAL(wp), DIMENSION(:,:),       INTENT(IN) ::  shaft_coordinates
    REAL(wp), DIMENSION(:,:,:,:),   INTENT(IN) ::  fast_position_blade, fast_force_blade,         &
                                                   palm_vel_blade

    IF ( myid /= 0 )                       RETURN
    IF ( meta%output_interval <= 0 )       RETURN

!-- PALM cadence: only tick on the first FAST sub-step of the PALM step.
    IF ( .NOT. meta%output_use_fast_timestep  .AND.  fast_sub_iter /= 1 )  RETURN

    st%step_counter = st%step_counter + 1
    IF ( MOD( st%step_counter, meta%output_interval ) /= 0 )  RETURN

    IF ( .NOT. st%file_open )  CALL open_file( meta )

    CALL write_one_step( meta, fast_sub_iter, current_time_fast,                                   &
                         fast_position_blade, fast_force_blade, palm_vel_blade,                    &
                         rotspeed, blade_pitch, shaft_coordinates )

 END SUBROUTINE write_wt_data


!--------------------------------------------------------------------------------------------------!
! Closes the file. Called from user_last_actions on every rank; non-rank-0 is a no-op.
!
! Before closing, re-enters define mode and writes two deferred global attributes that can only
! be known at the end of the run:
!   - close_time_iso8601    : wall-clock timestamp at finalize
!   - wall_clock_duration_s : seconds between file open and close, measured with SYSTEM_CLOCK
!
! In netCDF-4 the nf90_redef / nf90_enddef pair is essentially a no-op (HDF5 lets attributes be
! modified at any time), but the calls keep the code valid for classic-format builds too.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wt_output_finalize

    INTEGER(KIND=8) ::  close_count, close_count_rate
    REAL(wp)        ::  duration_s

    IF ( myid /= 0 )            RETURN
    IF ( .NOT. st%file_open )   RETURN

    CALL SYSTEM_CLOCK( count=close_count, count_rate=close_count_rate )
    duration_s = REAL( close_count - st%wall_open_count, KIND=wp ) /                               &
                 REAL( close_count_rate, KIND=wp )

    CALL nc_check( nf90_redef( st%ncid ),  'redef for close-time attributes' )
    CALL put_global_str ( 'close_time_iso8601',        iso8601_now() )
    CALL put_global_real( 'wall_clock_duration_s',     duration_s    )
    CALL put_global_str ( 'wall_clock_duration', seconds_to_human( duration_s ) )
    CALL nc_check( nf90_enddef( st%ncid ), 'enddef after close-time attributes' )

    CALL nc_check( nf90_close( st%ncid ), 'closing wt_output file' )
    st%file_open = .FALSE.
    st%ncid      = -1

 END SUBROUTINE wt_output_finalize


!--------------------------------------------------------------------------------------------------!
! Lazy file-open. Creates dims, vars, global attributes, then leaves define mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE open_file( meta )

    TYPE(wt_output_meta_t), INTENT(IN) ::  meta

    CHARACTER(LEN=256) ::  filename
    CHARACTER(LEN=8)   ::  time_tag

    WRITE( time_tag, '(I8.8)' )  MAX( 0, INT( simulated_time, KIND=iwp ) )
    filename = TRIM( run_identifier ) // '_wt_output_' // time_tag // '.nc'

!-- netCDF-4 (HDF5) format. Matches the PALM run namelist's netcdf_data_format (typically 4
!-- or 5 for PALM 24.04 builds linked against HDF5). Rank-0-only writes - we deliberately do
!-- not request NF90_MPIIO because only this rank holds the data and only this rank writes.
    CALL nc_check( nf90_create( TRIM(filename), IOR( NF90_CLOBBER, NF90_NETCDF4 ), st%ncid ),      &
                   'create ' // TRIM(filename) )

!-- Dimensions
    CALL nc_check( nf90_def_dim( st%ncid, 'time',    NF90_UNLIMITED,                st%dim_time    ), 'def_dim time'    )
    CALL nc_check( nf90_def_dim( st%ncid, 'turbine', meta%nturbines,                st%dim_turbine ), 'def_dim turbine' )
    CALL nc_check( nf90_def_dim( st%ncid, 'blade',   meta%fast_n_blades_max,        st%dim_blade   ), 'def_dim blade'   )
    CALL nc_check( nf90_def_dim( st%ncid, 'element', meta%fast_n_blade_elem_max,    st%dim_element ), 'def_dim element' )
    CALL nc_check( nf90_def_dim( st%ncid, 'xyz',     3,                             st%dim_xyz     ), 'def_dim xyz'     )

!-- Scalar-per-step variables. palm_sub_iter is integer; everything else is double.
    CALL def_var_1d( 'time',           (/ st%dim_time /),                                          &
                     st%var_time,          's',  'PALM simulated time' )
    CALL nc_check( nf90_def_var( st%ncid, 'palm_sub_iter', NF90_INT, (/ st%dim_time /),            &
                                 st%var_palm_sub_iter ),                                           &
                   'def_var palm_sub_iter' )
    CALL nc_check( nf90_put_att( st%ncid, st%var_palm_sub_iter, 'units',     '-' ),                &
                   'attr units palm_sub_iter' )
    CALL nc_check( nf90_put_att( st%ncid, st%var_palm_sub_iter, 'long_name',                       &
                                 'FAST sub-step counter l at write time' ),                       &
                   'attr long_name palm_sub_iter' )

!-- Per-turbine
    CALL def_var_2d( 'time_fast', (/ st%dim_turbine, st%dim_time /),                               &
                     st%var_time_fast, 's', 'FAST internal time per turbine' )
    CALL def_var_2d( 'rotspeed',  (/ st%dim_turbine, st%dim_time /),                               &
                     st%var_rotspeed,  'rad s-1', 'rotor angular speed per turbine' )
    CALL def_var_3d( 'shaft_coord', (/ st%dim_xyz, st%dim_turbine, st%dim_time /),                 &
                     st%var_shaft_coord, 'm', 'shaft coordinates (x,y,z) per turbine' )

!-- Per-(turbine, blade, element). Cartesian xyz on the trailing axis.
    CALL def_var_5d( 'position', (/ st%dim_xyz, st%dim_element, st%dim_blade,                      &
                                    st%dim_turbine, st%dim_time /),                                &
                     st%var_position, 'm',     'control-point position (x,y,z)' )
    CALL def_var_5d( 'force',    (/ st%dim_xyz, st%dim_element, st%dim_blade,                      &
                                    st%dim_turbine, st%dim_time /),                                &
                     st%var_force,    'm4 s-2', 'kinematic aerodynamic force at control point ' // &
                     '(= physical force / rhoairfast); multiply by rhoairfast (=AeroDyn AirDens) to get Newtons' )
    CALL def_var_5d( 'palm_vel', (/ st%dim_xyz, st%dim_element, st%dim_blade,                      &
                                    st%dim_turbine, st%dim_time /),                                &
                     st%var_palm_vel, 'm s-1', 'PALM velocity at control point after wake correction (x,y,z)' )

!-- palm_vel is interpolated and (when enabled) wake-corrected only on the first FAST
!-- sub-step of each PALM time step (l == 1 in fastv8_coupler_mod). For PALM-cadence
!-- output this is invisible because write_wt_data also gates on l == 1. For
!-- FAST-cadence output (output_use_fast_timestep = .TRUE.), the value is repeated
!-- unchanged across all l > 1 sub-steps within a PALM step. Document the cadence
!-- so analysis tools can distinguish per-PALM-step quantities from per-FAST-substep
!-- ones (position, force, time_fast, rotspeed all update per sub-step).
    CALL mark_palm_cadence( (/ st%var_palm_vel /) )

!-- NWC pre/post-correction diagnostics. Only defined when the near-wake smearing
!-- correction is enabled. Useful for investigating whether anomalies (e.g. spikes near
!-- a particular blade element) originate from the raw PALM-interpolated velocity field
!-- (visible in *_pre) or are introduced by the wake-correction model (visible in
!-- *_post but not *_pre, or in the induced nw_* arrays).
    IF ( meta%use_near_wake_smearing_correction )  THEN
       CALL def_var_4d( 'vn_pre',    (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_vn_pre,    'm s-1',                                                  &
                        'shaft-axial velocity at CP, raw (before NWC)' )
       CALL def_var_4d( 'vn_post',   (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_vn_post,   'm s-1',                                                  &
                        'shaft-axial velocity at CP, corrected (after NWC)' )
       CALL def_var_4d( 'vt_pre',    (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_vt_pre,    'm s-1',                                                  &
                        'signed tangential velocity at CP, raw (before NWC)' )
       CALL def_var_4d( 'vt_post',   (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_vt_post,   'm s-1',                                                  &
                        'signed tangential velocity at CP, corrected (after NWC)' )
       CALL def_var_4d( 'phi_pre',   (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_phi_pre,   'rad',                                                    &
                        'flow angle at CP, raw (before NWC, post pi/2 wrap)' )
       CALL def_var_4d( 'phi_post',  (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_phi_post,  'rad',                                                    &
                        'flow angle at CP, corrected (after NWC)' )
       CALL def_var_4d( 'alpha_pre', (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_alpha_pre, 'deg',                                                    &
                        'angle of attack at CP, raw (before NWC)' )
       CALL def_var_4d( 'alpha_post',(/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_alpha_post,'deg',                                                    &
                        'angle of attack at CP, corrected (after NWC)' )
       CALL def_var_4d( 'cl_pre',    (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_cl_pre,    '-',                                                      &
                        'lift coefficient at CP (polar-interpolated at alpha_pre)' )
       CALL def_var_4d( 'cl_post',   (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_cl_post,   '-',                                                      &
                        'lift coefficient at CP (polar-interpolated at alpha_post)' )
       CALL def_var_4d( 'cd_pre',    (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_cd_pre,    '-',                                                      &
                        'drag coefficient at CP (polar-interpolated at alpha_pre)' )
       CALL def_var_4d( 'cd_post',   (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_cd_post,   '-',                                                      &
                        'drag coefficient at CP (polar-interpolated at alpha_post)' )
       CALL def_var_4d( 'nw_vn',     (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_nw_vn,     'm s-1',                                                  &
                        'induced shaft-axial velocity from NWC (= vn_post - vn_pre)' )
       CALL def_var_4d( 'nw_vt',     (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_nw_vt,     'm s-1',                                                  &
                        'induced signed tangential velocity from NWC (= vt_post - vt_pre)' )
       CALL def_var_4d( 'nw_W',      (/ st%dim_element, st%dim_blade, st%dim_turbine, st%dim_time /), &
                        st%var_nw_W,      'm s-1',                                                  &
                        'induced wake velocity magnitude from NWC' )

!--    Blade pitch from FAST, one value per blade. Used by the NWC in the angle of
!--    attack (alpha = phi - twist - pitch). Stored in radians, written in degrees.
       CALL def_var_3d( 'blade_pitch', (/ st%dim_blade, st%dim_turbine, st%dim_time /),             &
                        st%var_blade_pitch, 'deg',                                                  &
                        'blade pitch from FAST, per blade' )

!--    Mark all 15 diag vars with samples_at = palm_step. They are computed inside the
!--    same IF ( l == 1 ) block as palm_vel in fastv8_coupler_mod, so they share its cadence.
       CALL mark_palm_cadence( (/                                                                  &
          st%var_vn_pre,    st%var_vn_post,                                                        &
          st%var_vt_pre,    st%var_vt_post,                                                        &
          st%var_phi_pre,   st%var_phi_post,                                                       &
          st%var_alpha_pre, st%var_alpha_post,                                                     &
          st%var_cl_pre,    st%var_cl_post,                                                        &
          st%var_cd_pre,    st%var_cd_post,                                                        &
          st%var_nw_vn,     st%var_nw_vt,    st%var_nw_W /) )
    ENDIF

!-- Global attributes - all the run / configuration state worth recording, set once. Boolean
!-- values are written as the literal strings "true" / "false" so a reader does not have to
!-- guess what 0 / 1 means.
    CALL put_global_str ( 'simulation_name',                       TRIM( run_identifier ) )
    CALL put_global_str ( 'creation_time_iso8601',                 iso8601_now() )
    CALL put_global_int ( 'nturbines',                             INT( meta%nturbines             ) )
    CALL put_global_int ( 'fast_n_blades_max',                     INT( meta%fast_n_blades_max     ) )
    CALL put_global_int ( 'fast_n_blade_elem_max',                 INT( meta%fast_n_blade_elem_max ) )
    CALL put_global_int ( 'output_interval',                       INT( meta%output_interval       ) )
    CALL put_global_str ( 'output_use_fast_timestep',              bool_str( meta%output_use_fast_timestep ) )
    CALL put_global_str ( 'output_sync_per_write',                 bool_str( meta%output_sync_per_write    ) )
    CALL put_global_str ( 'cadence_unit',                                                          &
                          MERGE( 'fast_sub_step', 'palm_step    ', meta%output_use_fast_timestep ) )
    CALL put_global_real( 'dt_palm',                               meta%dt_palm )
    CALL put_global_real( 'dt_fast',                               meta%dt_fast )
    CALL put_global_real( 'start_simulated_time',                  simulated_time )
    CALL put_global_real( 'end_simulated_time',                    meta%end_simulated_time )
    CALL put_global_str ( 'use_near_wake_smearing_correction',                                     &
                          bool_str( meta%use_near_wake_smearing_correction ) )
    CALL put_global_real( 'nw_trailed_vorticity_smoothing_length',                                 &
                          meta%nw_trailed_vorticity_smoothing_length )
    CALL put_global_str ( 'is_restart',                            bool_str( meta%is_restart ) )

!-- PALM grid and MPI decomposition
    CALL put_global_real( 'dx',                                    meta%dx                       )
    CALL put_global_real( 'dy',                                    meta%dy                       )
    CALL put_global_real( 'dz',                                    meta%dz                       )
    CALL put_global_int ( 'nx',                                    INT( meta%nx   )              )
    CALL put_global_int ( 'ny',                                    INT( meta%ny   )              )
    CALL put_global_int ( 'nz',                                    INT( meta%nz   )              )
    CALL put_global_int ( 'npex',                                  INT( meta%npex )              )
    CALL put_global_int ( 'npey',                                  INT( meta%npey )              )

!-- Wind-turbine setup scalars
    CALL put_global_real( 'time_turbine_on',                       meta%time_turbine_on          )
    CALL put_global_real( 'reg_fac',                               meta%reg_fac                  )
    CALL put_global_real( 'fbox_fac',                              meta%fbox_fac                 )
    CALL put_global_int ( 'n_airfoil_polar_rows',                  INT( meta%n_airfoil_polar_rows ) )

!-- Per-turbine arrays (length nturbines)
    IF ( ALLOCATED( meta%palm_tower_ref_pos_x ) )                                                  &
       CALL put_global_real_array( 'palm_tower_ref_pos_x', meta%palm_tower_ref_pos_x )
    IF ( ALLOCATED( meta%palm_tower_ref_pos_y ) )                                                  &
       CALL put_global_real_array( 'palm_tower_ref_pos_y', meta%palm_tower_ref_pos_y )
    IF ( ALLOCATED( meta%palm_tower_ref_pos_z ) )                                                  &
       CALL put_global_real_array( 'palm_tower_ref_pos_z', meta%palm_tower_ref_pos_z )
    IF ( ALLOCATED( meta%dtow ) )                                                                  &
       CALL put_global_real_array( 'dtow',         meta%dtow         )
    IF ( ALLOCATED( meta%htow ) )                                                                  &
       CALL put_global_real_array( 'htow',         meta%htow         )
    IF ( ALLOCATED( meta%turb_C_d_tow ) )                                                          &
       CALL put_global_real_array( 'turb_C_d_tow', meta%turb_C_d_tow )

!-- Semicolon-joined string arrays (CHARACTER attributes are scalar; the consumer splits on ';')
    IF ( ALLOCATED( meta%fast_host_ports_joined ) )                                                &
       CALL put_global_str( 'fast_host_port',     meta%fast_host_ports_joined )

!-- Free-form notes attribute. Documents where the close-time attributes come from and
!-- the two known reasons they can be absent. Helpful for anyone reading the file without
!-- the PALM log.
    CALL put_global_str( 'notes', build_notes() )

    CALL nc_check( nf90_enddef( st%ncid ), 'enddef' )

!-- Capture wall-clock at file-open. The matching close-time + duration is written by
!-- wt_output_finalize as a deferred set of global attributes (re-enters define mode briefly).
    CALL SYSTEM_CLOCK( count=st%wall_open_count, count_rate=st%wall_open_count_rate )

    st%file_open  = .TRUE.
    st%time_index = 0_iwp

 END SUBROUTINE open_file


!--------------------------------------------------------------------------------------------------!
! Append one time slice to every variable. Indexes are 1-based for NetCDF.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE write_one_step( meta, fast_sub_iter, time_fast_arr, position, force, palm_vel,         &
                            rotspeed_arr, blade_pitch_arr, shaft_arr )

    TYPE(wt_output_meta_t), INTENT(IN) ::  meta
    INTEGER(iwp),           INTENT(IN) ::  fast_sub_iter

    REAL(wp), DIMENSION(:),         INTENT(IN) ::  time_fast_arr, rotspeed_arr
    REAL(wp), DIMENSION(:,:),       INTENT(IN) ::  blade_pitch_arr   !< (turbine, blade), radians
    REAL(wp), DIMENSION(:,:),       INTENT(IN) ::  shaft_arr
    REAL(wp), DIMENSION(:,:,:,:),   INTENT(IN) ::  position, force, palm_vel

    INTEGER(iwp)            ::  i, j, k, c
    INTEGER                 ::  t
    INTEGER(iwp)            ::  n_turb, n_blade, n_elem
    REAL(wp), ALLOCATABLE   ::  pos_5d(:,:,:,:,:), frc_5d(:,:,:,:,:), vel_5d(:,:,:,:,:)
    REAL(wp), ALLOCATABLE   ::  shaft_3d(:,:,:)
    REAL(wp), ALLOCATABLE   ::  pitch_3d(:,:,:)

    n_turb  = meta%nturbines
    n_blade = meta%fast_n_blades_max
    n_elem  = meta%fast_n_blade_elem_max

    st%time_index = st%time_index + 1_iwp
    t = INT( st%time_index )

!-- scalar-per-step
    CALL nc_check( nf90_put_var( st%ncid, st%var_time,          (/ simulated_time /),                       &
                                 start=(/ t /), count=(/ 1 /) ), 'put time' )
    CALL nc_check( nf90_put_var( st%ncid, st%var_palm_sub_iter, (/ INT( fast_sub_iter, KIND=4 ) /),         &
                                 start=(/ t /), count=(/ 1 /) ), 'put palm_sub_iter' )

!-- per-turbine 1D-per-step
    CALL nc_check( nf90_put_var( st%ncid, st%var_time_fast, time_fast_arr,                                  &
                                 start=(/ 1, t /), count=(/ INT(n_turb), 1 /) ), 'put time_fast' )
    CALL nc_check( nf90_put_var( st%ncid, st%var_rotspeed,  rotspeed_arr,                                   &
                                 start=(/ 1, t /), count=(/ INT(n_turb), 1 /) ), 'put rotspeed' )

!-- per-blade pitch, written only when the NWC is active (the variable is defined only
!-- then). Source blade_pitch_arr is (turbine, blade) in radians; reorder to
!-- (blade, turbine, 1) and convert to degrees for the file.
    IF ( meta%use_near_wake_smearing_correction )  THEN
       ALLOCATE( pitch_3d(n_blade, n_turb, 1) )
       DO i = 1, n_turb
          DO j = 1, n_blade
             pitch_3d(j, i, 1) = blade_pitch_arr(i, j) * 180.0_wp / pi
          END DO
       END DO
       CALL nc_check( nf90_put_var( st%ncid, st%var_blade_pitch, pitch_3d,                                  &
                                    start=(/ 1, 1, t /),                                                    &
                                    count=(/ INT(n_blade), INT(n_turb), 1 /) ), 'put blade_pitch' )
       DEALLOCATE( pitch_3d )
    ENDIF

!-- per-turbine xyz: shaft_arr is (nturbines, 3). Reorder to (3, nturbines, 1) for the
!-- (xyz, turbine, time) layout.
    ALLOCATE( shaft_3d(3, n_turb, 1) )
    DO i = 1, n_turb
       DO c = 1, 3
          shaft_3d(c, i, 1) = shaft_arr(i, c)
       END DO
    END DO
    CALL nc_check( nf90_put_var( st%ncid, st%var_shaft_coord, shaft_3d,                                     &
                                 start=(/ 1, 1, t /), count=(/ 3, INT(n_turb), 1 /) ), 'put shaft_coord' )
    DEALLOCATE( shaft_3d )

!-- 4D source arrays in the coupling are indexed (element, blade, turbine, xyz). NetCDF layout
!-- is (xyz, element, blade, turbine, time) so the xyz axis varies fastest on disk - matches the
!-- physical "x,y,z" triplet a reader expects to find adjacent. Reorder here.
    ALLOCATE( pos_5d(3, n_elem, n_blade, n_turb, 1) )
    ALLOCATE( frc_5d(3, n_elem, n_blade, n_turb, 1) )
    ALLOCATE( vel_5d(3, n_elem, n_blade, n_turb, 1) )
    DO i = 1, n_turb
       DO j = 1, n_blade
          DO k = 1, n_elem
             DO c = 1, 3
                pos_5d(c, k, j, i, 1) = position(k, j, i, c)
                frc_5d(c, k, j, i, 1) = force   (k, j, i, c)
                vel_5d(c, k, j, i, 1) = palm_vel(k, j, i, c)
             END DO
          END DO
       END DO
    END DO

    CALL nc_check( nf90_put_var( st%ncid, st%var_position, pos_5d,                                          &
                                 start=(/ 1, 1, 1, 1, t /),                                                 &
                                 count=(/ 3, INT(n_elem), INT(n_blade), INT(n_turb), 1 /) ),                &
                   'put position' )
    CALL nc_check( nf90_put_var( st%ncid, st%var_force,    frc_5d,                                          &
                                 start=(/ 1, 1, 1, 1, t /),                                                 &
                                 count=(/ 3, INT(n_elem), INT(n_blade), INT(n_turb), 1 /) ),                &
                   'put force' )
    CALL nc_check( nf90_put_var( st%ncid, st%var_palm_vel, vel_5d,                                          &
                                 start=(/ 1, 1, 1, 1, t /),                                                 &
                                 count=(/ 3, INT(n_elem), INT(n_blade), INT(n_turb), 1 /) ),                &
                   'put palm_vel' )

    DEALLOCATE( pos_5d, frc_5d, vel_5d )

!-- NWC pre/post diagnostics: source arrays in nw_smearing_correction_module are
!-- (turbine, blade, element). Reorder to (element, blade, turbine, 1) for the
!-- NetCDF (element, blade, turbine, time) layout.
    IF ( meta%use_near_wake_smearing_correction )  CALL write_diag_step( t, n_turb, n_blade, n_elem )

!-- Flush metadata to disk so a killed run still leaves a valid file. On parallel filesystems
!-- this is the dominant per-write cost; meta%output_sync_per_write = .FALSE. skips it.
!-- The file is always synced on close (wt_output_finalize) regardless.
    IF ( meta%output_sync_per_write )  CALL nc_check( nf90_sync( st%ncid ), 'sync' )

 END SUBROUTINE write_one_step


!--------------------------------------------------------------------------------------------------!
! Write the 15 NWC pre/post/induced diagnostic arrays for one time slice. Called only when
! meta%use_near_wake_smearing_correction is set. Pulls the source arrays straight out of
! nw_smearing_correction_module (USEd at module head) and reorders each from the
! (turbine, blade, element) layout used inside the wake-model code into the
! (element, blade, turbine, time) layout used on disk.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE write_diag_step( t, n_turb, n_blade, n_elem )

    INTEGER,      INTENT(IN) ::  t
    INTEGER(iwp), INTENT(IN) ::  n_turb, n_blade, n_elem

    REAL(wp), ALLOCATABLE ::  buf(:,:,:,:)

    ALLOCATE( buf(n_elem, n_blade, n_turb, 1) )

    CALL pack_diag( nwc_diag_vn_pre,     buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_vn_pre,       buf, t, n_turb, n_blade, n_elem, 'vn_pre' )
    CALL pack_diag( nwc_diag_vn_post,    buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_vn_post,      buf, t, n_turb, n_blade, n_elem, 'vn_post' )
    CALL pack_diag( nwc_diag_vt_pre,     buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_vt_pre,       buf, t, n_turb, n_blade, n_elem, 'vt_pre' )
    CALL pack_diag( nwc_diag_vt_post,    buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_vt_post,      buf, t, n_turb, n_blade, n_elem, 'vt_post' )
    CALL pack_diag( nwc_diag_phi_pre,    buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_phi_pre,      buf, t, n_turb, n_blade, n_elem, 'phi_pre' )
    CALL pack_diag( nwc_diag_phi_post,   buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_phi_post,     buf, t, n_turb, n_blade, n_elem, 'phi_post' )
    CALL pack_diag( nwc_diag_alpha_pre,  buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_alpha_pre,    buf, t, n_turb, n_blade, n_elem, 'alpha_pre' )
    CALL pack_diag( nwc_diag_alpha_post, buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_alpha_post,   buf, t, n_turb, n_blade, n_elem, 'alpha_post' )
    CALL pack_diag( nwc_diag_cl_pre,     buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_cl_pre,       buf, t, n_turb, n_blade, n_elem, 'cl_pre' )
    CALL pack_diag( nwc_diag_cl_post,    buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_cl_post,      buf, t, n_turb, n_blade, n_elem, 'cl_post' )
    CALL pack_diag( nwc_diag_cd_pre,     buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_cd_pre,       buf, t, n_turb, n_blade, n_elem, 'cd_pre' )
    CALL pack_diag( nwc_diag_cd_post,    buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_cd_post,      buf, t, n_turb, n_blade, n_elem, 'cd_post' )
    CALL pack_diag( nwc_diag_nw_vn,      buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_nw_vn,        buf, t, n_turb, n_blade, n_elem, 'nw_vn' )
    CALL pack_diag( nwc_diag_nw_vt,      buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_nw_vt,        buf, t, n_turb, n_blade, n_elem, 'nw_vt' )
    CALL pack_diag( nwc_diag_nw_W,       buf, n_turb, n_blade, n_elem )
    CALL put_diag ( st%var_nw_W,         buf, t, n_turb, n_blade, n_elem, 'nw_W' )

    DEALLOCATE( buf )

 END SUBROUTINE write_diag_step


 SUBROUTINE pack_diag( src, dst, n_turb, n_blade, n_elem )
    REAL(wp),     INTENT(IN)  ::  src(:,:,:)              !< (turbine, blade, element)
    REAL(wp),     INTENT(OUT) ::  dst(:,:,:,:)            !< (element, blade, turbine, 1)
    INTEGER(iwp), INTENT(IN)  ::  n_turb, n_blade, n_elem

    INTEGER(iwp) ::  i, j, k

    DO i = 1, n_turb
       DO j = 1, n_blade
          DO k = 1, n_elem
             dst(k, j, i, 1) = src(i, j, k)
          END DO
       END DO
    END DO

 END SUBROUTINE pack_diag


 SUBROUTINE put_diag( varid, buf, t, n_turb, n_blade, n_elem, name )
    INTEGER,          INTENT(IN) ::  varid, t
    REAL(wp),         INTENT(IN) ::  buf(:,:,:,:)
    INTEGER(iwp),     INTENT(IN) ::  n_turb, n_blade, n_elem
    CHARACTER(LEN=*), INTENT(IN) ::  name

    CALL nc_check( nf90_put_var( st%ncid, varid, buf,                                              &
                                 start=(/ 1, 1, 1, t /),                                           &
                                 count=(/ INT(n_elem), INT(n_blade), INT(n_turb), 1 /) ),          &
                   'put ' // name )
 END SUBROUTINE put_diag


!--------------------------------------------------------------------------------------------------!
! Helpers to keep the long def_var calls one-line. All set the variable type to NF90_DOUBLE
! and attach 'units' and 'long_name'.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE def_var_1d( name, dims, varid, units, long_name )
    CHARACTER(LEN=*),       INTENT(IN)  ::  name, units, long_name
    INTEGER, DIMENSION(:),  INTENT(IN)  ::  dims
    INTEGER,                INTENT(OUT) ::  varid

    CALL nc_check( nf90_def_var( st%ncid, name, NF90_DOUBLE, dims, varid ), 'def_var ' // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'units',     units     ), 'attr units '     // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'long_name', long_name ), 'attr long_name ' // name )
 END SUBROUTINE def_var_1d

 SUBROUTINE def_var_2d( name, dims, varid, units, long_name )
    CHARACTER(LEN=*),       INTENT(IN)  ::  name, units, long_name
    INTEGER, DIMENSION(2),  INTENT(IN)  ::  dims
    INTEGER,                INTENT(OUT) ::  varid

    CALL nc_check( nf90_def_var( st%ncid, name, NF90_DOUBLE, dims, varid ), 'def_var ' // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'units',     units     ), 'attr units '     // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'long_name', long_name ), 'attr long_name ' // name )
 END SUBROUTINE def_var_2d

 SUBROUTINE def_var_3d( name, dims, varid, units, long_name )
    CHARACTER(LEN=*),       INTENT(IN)  ::  name, units, long_name
    INTEGER, DIMENSION(3),  INTENT(IN)  ::  dims
    INTEGER,                INTENT(OUT) ::  varid

    CALL nc_check( nf90_def_var( st%ncid, name, NF90_DOUBLE, dims, varid ), 'def_var ' // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'units',     units     ), 'attr units '     // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'long_name', long_name ), 'attr long_name ' // name )
 END SUBROUTINE def_var_3d

 SUBROUTINE def_var_4d( name, dims, varid, units, long_name )
    CHARACTER(LEN=*),       INTENT(IN)  ::  name, units, long_name
    INTEGER, DIMENSION(4),  INTENT(IN)  ::  dims
    INTEGER,                INTENT(OUT) ::  varid

    CALL nc_check( nf90_def_var( st%ncid, name, NF90_DOUBLE, dims, varid ), 'def_var ' // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'units',     units     ), 'attr units '     // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'long_name', long_name ), 'attr long_name ' // name )
 END SUBROUTINE def_var_4d

 SUBROUTINE def_var_5d( name, dims, varid, units, long_name )
    CHARACTER(LEN=*),       INTENT(IN)  ::  name, units, long_name
    INTEGER, DIMENSION(5),  INTENT(IN)  ::  dims
    INTEGER,                INTENT(OUT) ::  varid

    CALL nc_check( nf90_def_var( st%ncid, name, NF90_DOUBLE, dims, varid ), 'def_var ' // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'units',     units     ), 'attr units '     // name )
    CALL nc_check( nf90_put_att( st%ncid, varid, 'long_name', long_name ), 'attr long_name ' // name )
 END SUBROUTINE def_var_5d


!--------------------------------------------------------------------------------------------------!
! Tiny wrappers to keep the global-attribute put calls short.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE put_global_str( name, val )
    CHARACTER(LEN=*), INTENT(IN) ::  name, val
    CALL nc_check( nf90_put_att( st%ncid, NF90_GLOBAL, name, val ), 'attr ' // name )
 END SUBROUTINE put_global_str


!--------------------------------------------------------------------------------------------------!
! Tag a batch of variables with samples_at = "palm_step". Used for quantities that are
! refreshed once per PALM time step (inside the IF ( l == 1 ) block in fastv8_coupler_mod)
! and repeated unchanged across FAST sub-steps in FAST-cadence output. Helps analysis tools
! distinguish them from per-sub-step quantities (position, force, time_fast, rotspeed).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE mark_palm_cadence( varids )
    INTEGER, INTENT(IN) ::  varids(:)
    INTEGER             ::  m

    DO m = 1, SIZE( varids )
       CALL nc_check( nf90_put_att( st%ncid, varids(m), 'samples_at', 'palm_step' ),                &
                      'attr samples_at' )
    END DO
 END SUBROUTINE mark_palm_cadence

 SUBROUTINE put_global_int( name, val )
    CHARACTER(LEN=*), INTENT(IN) ::  name
    INTEGER,          INTENT(IN) ::  val
    CALL nc_check( nf90_put_att( st%ncid, NF90_GLOBAL, name, val ), 'attr ' // name )
 END SUBROUTINE put_global_int

 SUBROUTINE put_global_real( name, val )
    CHARACTER(LEN=*), INTENT(IN) ::  name
    REAL(wp),         INTENT(IN) ::  val
    CALL nc_check( nf90_put_att( st%ncid, NF90_GLOBAL, name, val ), 'attr ' // name )
 END SUBROUTINE put_global_real

 SUBROUTINE put_global_real_array( name, vals )
    CHARACTER(LEN=*),       INTENT(IN) ::  name
    REAL(wp), DIMENSION(:), INTENT(IN) ::  vals
    CALL nc_check( nf90_put_att( st%ncid, NF90_GLOBAL, name, vals ), 'attr ' // name )
 END SUBROUTINE put_global_real_array


!--------------------------------------------------------------------------------------------------!
! Render a Fortran LOGICAL as the literal string "true" or "false". Deferred-length result so
! the returned string carries no trailing pad - what you put is what lands in the NetCDF
! attribute. Used so global-attribute booleans are self-explanatory in the file (no 0/1 to
! guess at).
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
! Build the free-form 'notes' global attribute. Always written. Documents where the close-time
! attributes come from and the two situations in which they can end up absent so a reader of
! the file alone (without the PALM log) understands what to check.
!--------------------------------------------------------------------------------------------------!
 FUNCTION build_notes() RESULT( notes )
    CHARACTER(LEN=:), ALLOCATABLE ::  notes

    notes = 'If close_time_iso8601 and wall_clock_duration_s are missing from this file, '   //    &
            'make sure the _p3d file contains an &user_parameters block and that the '       //    &
            'simulation finished cleanly (no timeout, no abort). '                            //    &
            'Variables tagged samples_at = "palm_step" (palm_vel and the NWC diag arrays) '   //    &
            'are recomputed once per PALM time step (inside IF (l == 1)). With '              //    &
            'output_use_fast_timestep = .TRUE. their values repeat unchanged across FAST '    //    &
            'sub-steps within a PALM step; with the default output_use_fast_timestep '        //    &
            '= .FALSE. only l == 1 is ever written, so the distinction is invisible. '        //    &
            'Other per-element / per-turbine arrays (position, force, time_fast, rotspeed) '  //    &
            'update on every FAST sub-step.'
 END FUNCTION build_notes


 FUNCTION bool_str( b ) RESULT( s )
    LOGICAL, INTENT(IN)           ::  b
    CHARACTER(LEN=:), ALLOCATABLE ::  s
    IF ( b )  THEN
       s = 'true'
    ELSE
       s = 'false'
    END IF
 END FUNCTION bool_str


!--------------------------------------------------------------------------------------------------!
! Format a duration given in seconds as a short human-readable string. Examples:
!   143.2        -> "0 hours 2 minutes 23 seconds"
!   5000.0       -> "1 hours 23 minutes 20 seconds"
!   100000.0     -> "1 days 3 hours 46 minutes 40 seconds"
! Days are emitted only when > 0; the rest are always shown so the unit structure is consistent.
! Returned as a deferred-length allocatable so the result has exactly the right length.
!--------------------------------------------------------------------------------------------------!
 FUNCTION seconds_to_human( secs ) RESULT( s )
    REAL(wp), INTENT(IN)          ::  secs
    CHARACTER(LEN=:), ALLOCATABLE ::  s

    INTEGER(KIND=8)               ::  total, days, hours, mins, sec_i
    CHARACTER(LEN=128)            ::  buf

    total = INT( MAX( secs, 0.0_wp ), KIND=8 )
    sec_i = MOD( total, 60_8 )
    total = total / 60_8
    mins  = MOD( total, 60_8 )
    total = total / 60_8
    hours = MOD( total, 24_8 )
    days  = total / 24_8

    IF ( days > 0_8 )  THEN
       WRITE( buf, '(I0,A,I0,A,I0,A,I0,A)' )                                                       &
          days,  ' days ',  hours, ' hours ', mins, ' minutes ', sec_i, ' seconds'
    ELSE
       WRITE( buf, '(I0,A,I0,A,I0,A)' )                                                            &
          hours, ' hours ', mins, ' minutes ', sec_i, ' seconds'
    END IF

    s = TRIM( buf )
 END FUNCTION seconds_to_human


!--------------------------------------------------------------------------------------------------!
! Return current local wall-clock time as ISO-8601 ("YYYY-MM-DDThh:mm:ss"). Used for the
! creation_time_iso8601 global attribute so a reader can tell when the file was written
! independently of simulated_time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION iso8601_now() RESULT( ts )
    CHARACTER(LEN=19) ::  ts
    CHARACTER(LEN=8)  ::  date
    CHARACTER(LEN=10) ::  time

    CALL DATE_AND_TIME( date=date, time=time )
    ts = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' //                               &
         time(1:2) // ':' // time(3:4) // ':' // time(5:6)
 END FUNCTION iso8601_now


!--------------------------------------------------------------------------------------------------!
! Abort on any NetCDF error. Prints the failing context and the library's error text.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE nc_check( status, context )
    INTEGER,          INTENT(IN) ::  status
    CHARACTER(LEN=*), INTENT(IN) ::  context

    IF ( status /= NF90_NOERR )  THEN
       WRITE(*,*) '[wt_output_mod] NetCDF error during ', TRIM( context ),                          &
                  ' : ', TRIM( nf90_strerror( status ) )
       STOP 'wt_output_mod: NetCDF error'
    END IF
 END SUBROUTINE nc_check


END MODULE wt_output_mod
