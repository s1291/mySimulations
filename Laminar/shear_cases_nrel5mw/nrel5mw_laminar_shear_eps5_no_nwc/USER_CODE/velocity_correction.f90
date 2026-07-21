!------------------------------------------------------------------------------------
! This file should be a linkage between the Palm output of the FAST - PALM coupling
! and the correction model developed by Alexander Meyer-Forsting (DTU) to correct the
! missing induction of trailed vortices of turbine blades.
!------------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------------
! Explanation of parameters used in the correction model

! INTEGER,INTENT(IN)::
! i = ith rotor,  definition in Palm: nturbines
! j = jth blade,  definition in Palm: fast_n_blades 
! len_2d = length of 2D airfoil data. Read from the namelist via n_airfoil_polar_rows

! REAL(kind=idp),INTENT(IN):: 
! r_cp(n_s) = Radial coordinate of control points (CP), definition in the Fast file NRELOffshrBsline5MW_Onshore_AD.ipt   
! r_vtp(n_v) = Radial coordinate of vortex trailing points, defined as the middlepoints between the controlpoints plus tip and rootpoint of the blade
! vn_cfd(n_s) = Normal/axial velocities at CPs from CFD, projected onto the FAST shaft direction
! vt_cfd(n_s) = Tangential velocities at CPs from CFD, projected onto the local positive-rotation direction
! eps_cp = the smearing length epsilon, definition in Palm: epsilon = reg_fac * dx
! omega = Rotational speed in radiens/sec, definition in Palm: RotSpeed
! chord(n_s) = Chord distribution along blade at CPs, definition in the Fast file NRELOffshrBsline5MW_Onshore_AD.ipt
! twist(n_s) = Twist distribution along blade at CPs, definition in the Fast file NRELOffshrBsline5MW_Onshore_AD.ipt
! pitch = Blade pitch in degrees. Taken from FAST (ED%y%BlPitch, radians) through
!         the coupler array blade_pitch and converted to degrees here.
! cl_2d(len_2d,n_s) = Airfoil lift coeff. curve, definition in the Fast files at NREL/AeroData/NREL5M/
! alpha_2d(len_2d,n_s) = Alpha values of lift&drag curves, definition in the Fast files at NREL/AeroData/NREL5M/


!------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------
! Calculation of the correction of the velocities at CPs 

SUBROUTINE velocity_correction( len_2d,                                                            &
                                fast_n_blade_elem_max,                                             &
                                fast_position_blade,                                               &
                                shaft_coordinates,                                                 &
                                palm_vel_blade,                                                    &
                                fast_n_blades_max,                                                 &
                                nturbines,                                                         &
                                RotSpeed,                                                          &
                                blade_pitch,                                                       &
                                eps,                                                               &
                                dt_palm,                                                           &
                                nw_trailed_vorticity_smoothing_length,                             &
                                palm_vel_b )

    use kinds, ONLY: wp, iwp
    USE nw_smearing_correction_module, ONLY: chord_arr, twist_arr, r_cp_arr, r_vtp_arr,            &
                                             alpha_2d_arr, cl_2d_arr, pi

    IMPLICIT NONE

!-- Definition of input- and output-parameter for using the correction model

!-- Input
    INTEGER(iwp), INTENT(IN) ::  fast_n_blade_elem_max
    INTEGER(iwp), INTENT(IN) ::  fast_n_blades_max
    INTEGER(iwp), INTENT(IN) ::  nturbines
    INTEGER(iwp), INTENT(IN) ::  len_2d

    REAL(wp), INTENT(IN) ::  eps
    REAL(wp), INTENT(IN) ::  dt_palm
    REAL(wp), INTENT(IN) ::  nw_trailed_vorticity_smoothing_length  !< namelist parameter from fastv8_coupler_mod
    REAL(wp), DIMENSION(nturbines), INTENT(IN) :: RotSpeed
    REAL(wp), DIMENSION(nturbines, fast_n_blades_max), INTENT(IN) :: blade_pitch  !< blade pitch from FAST (radians)
    REAL(wp), DIMENSION(nturbines, 3), INTENT(IN) ::  shaft_coordinates
    REAL(wp), DIMENSION(fast_n_blade_elem_max, fast_n_blades_max, nturbines, 3), INTENT(IN) ::  palm_vel_blade
    REAL(wp), DIMENSION(fast_n_blade_elem_max, fast_n_blades_max, nturbines, 3), INTENT(IN) ::  fast_position_blade

!-- Output
    REAL(wp), DIMENSION(fast_n_blade_elem_max, fast_n_blades_max, nturbines,3), INTENT(OUT) ::  palm_vel_b  !< Corrected velocities in u,v and w direction

!-- Variables. The airfoil polar tables (alpha_2d_arr, cl_2d_arr) and
!-- blade geometry (chord_arr, twist_arr, r_cp_arr, r_vtp_arr) live in
!-- nw_smearing_correction_module and are loaded once at init by f8c_init_arrays.
    REAL(wp), DIMENSION(fast_n_blade_elem_max,fast_n_blades_max,nturbines,3) ::  palm_vel
    REAL(wp), DIMENSION(fast_n_blade_elem_max) ::  vn_cfd_out
    REAL(wp), DIMENSION(fast_n_blade_elem_max) ::  vt_cfd_out
    REAL(wp), DIMENSION(fast_n_blade_elem_max) ::  vn !< raw shaft-axial velocity; input to nw_smearing_correction
    REAL(wp), DIMENSION(fast_n_blade_elem_max) ::  vt !< raw signed tangential velocity; input to nw_smearing_correction
    REAL(wp) ::  pitch
    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j


!-- palm_vel is used for in- and output in CALL velocity_kart, so palm_vel should be first equal to palm_vel_blade
    palm_vel = palm_vel_blade


!-- Calculate Correction for Cps at blades of the used windturbines

    DO i = 1, nturbines

        DO  j = 1, fast_n_blades_max

!--         Blade pitch from FAST is in radians. The near-wake correction works
!--         in degrees (alpha = phi - twist - pitch, all in degrees), so convert.
            pitch = blade_pitch(i, j) * 180.0_wp / pi

!--         Calculate axial and tangential velocities from cartesian description of velocities given
!--          by PALM 
            CALL velocity_ax_tan( fast_position_blade,                                             &
                                  shaft_coordinates,                                               &
                                  palm_vel_blade,                                                  &
                                  fast_n_blade_elem_max,                                           &
                                  fast_n_blades_max,                                               &
                                  nturbines,                                                       &
                                  i,                                                               &
                                  j,                                                               &
                                  vn,                                                              &
                                  vt )     

!--         Calculate Correction with nw_smearing_correction.f90
            CALL nw_smearing_correction( i,                                                        &
                                         j,                                                        &
                                         r_cp_arr,                                                 &
                                         r_vtp_arr,                                                &
                                         vn,                                                       &
                                         vt,                                                       &
                                         eps,                                                      &
                                         RotSpeed(i),                                              &
                                         chord_arr,                                                &
                                         twist_arr,                                                &
                                         pitch,                                                    &
                                         cl_2d_arr,                                                &
                                         alpha_2d_arr,                                             &
                                         len_2d,                                                   &
                                         dt_palm,                                                  &
                                         nw_trailed_vorticity_smoothing_length,                    &
                                         vn_cfd_out,                                               &
                                         vt_cfd_out )


!--         Calculate velocities in cartesian description from corrected axial and tangential
!--         velocities given by nw_smewaring_correction.f90 
            CALL velocity_kart( fast_position_blade,                                               & 
                                shaft_coordinates,                                                 &
                                fast_n_blade_elem_max,                                             &
                                fast_n_blades_max,                                                 &
                                nturbines,                                                         &
                                i,                                                                 &
                                j,                                                                 &
                                vn_cfd_out,                                                        &
                                vt_cfd_out,                                                        &
                                palm_vel )
    
            palm_vel_b = palm_vel

        END DO ! j = 1, fast_n_blades_max
    END DO ! i = 1, nturbines


END SUBROUTINE velocity_correction
!--------------------------------------------------------------------------------------------------





SUBROUTINE velocity_ax_tan ( fast_position_blade,                                                  &
                             shaft_coordinates,                                                    &
                             palm_vel_blade,                                                       &
                             fast_n_blade_elem_max,                                                &
                             fast_n_blades_max,                                                    &
                             nturbines,                                                            &
                             i,                                                                    &
                             j,                                                                    &
                             vax,                                                                  &
                             vtan )        

    USE kinds, ONLY: wp, iwp
    
    IMPLICIT NONE
!-- Input
    INTEGER(iwp), INTENT(IN) ::  fast_n_blade_elem_max
    INTEGER(iwp), INTENT(IN) ::  fast_n_blades_max
    INTEGER(iwp), INTENT(IN) ::  nturbines
    INTEGER(iwp), INTENT(IN) ::  i
    INTEGER(iwp), INTENT(IN) ::  j
    REAL(wp), DIMENSION(fast_n_blade_elem_max,fast_n_blades_max,nturbines,3), INTENT(IN) ::  palm_vel_blade
    REAL(wp), DIMENSION(fast_n_blade_elem_max,fast_n_blades_max,nturbines,3), INTENT(IN) ::  fast_position_blade
    REAL(wp), DIMENSION(nturbines,3), INTENT(IN) ::  shaft_coordinates

!-- Output
    REAL(wp), DIMENSION(fast_n_blade_elem_max), INTENT(OUT)::  vax
    REAL(wp), DIMENSION(fast_n_blade_elem_max), INTENT(OUT)::  vtan

!-- Variables
    REAL(wp), DIMENSION(3) ::  a
    REAL(wp), DIMENSION(3) ::  e_ax
    REAL(wp), DIMENSION(3) ::  e_rad
    REAL(wp), DIMENSION(3) ::  e_tan
    REAL(wp), DIMENSION(3) ::  u_cfd
    REAL(wp) ::  norm_axis
    REAL(wp) ::  norm_rad
    REAL(wp) ::  norm_tan
    INTEGER(iwp) ::  n

!-- Normalize the FAST shaft axis. For zero shaft tilt this is expected to be
!-- (1,0,0), but normalizing here keeps the projection robust if the coupling
!-- sends a nearly-unit vector.

    norm_axis = sqrt( sum( shaft_coordinates(i,1:3)**2 ) )
    IF ( norm_axis > epsilon( norm_axis ) ) THEN
       e_ax = shaft_coordinates(i,1:3) / norm_axis
    ELSE
       e_ax = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    END IF

!-- Calculate the axial and tangential velocities at CPs in the local rotor
!-- basis. The positive tangential direction is e_ax x e_rad, consistent with
!-- positive FAST RotSpeed and with the old y-z formula when e_ax=(1,0,0).

    DO n = 1, fast_n_blade_elem_max

       u_cfd = palm_vel_blade(n, j, i, 1:3)

       ! Direction vector at blades (Assumption: no bend of blade)
       a = fast_position_blade(3, j, i, 1:3) - fast_position_blade(1, j, i, 1:3)

       ! Remove any component along the shaft so e_rad lies in the rotor plane.
       e_rad = a - dot_product( a, e_ax ) * e_ax
       norm_rad = sqrt( sum( e_rad**2 ) )
       IF ( norm_rad > epsilon( norm_rad ) ) THEN
          e_rad = e_rad / norm_rad
       ELSE
          ! Degenerate blade-axis data: choose a deterministic vector perpendicular
          ! to the shaft, instead of silently using a non-orthogonal fallback.
          IF ( abs( e_ax(1) ) < 0.9_wp ) THEN
             e_rad = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             e_rad = (/ 0.0_wp, 1.0_wp, 0.0_wp /)
          END IF
          e_rad = e_rad - dot_product( e_rad, e_ax ) * e_ax
          e_rad = e_rad / sqrt( sum( e_rad**2 ) )
       END IF

       e_tan(1) = e_ax(2) * e_rad(3) - e_ax(3) * e_rad(2)
       e_tan(2) = e_ax(3) * e_rad(1) - e_ax(1) * e_rad(3)
       e_tan(3) = e_ax(1) * e_rad(2) - e_ax(2) * e_rad(1)
       norm_tan = sqrt( sum( e_tan**2 ) )
       IF ( norm_tan > epsilon( norm_tan ) ) THEN
          e_tan = e_tan / norm_tan
       END IF

       vax(n)  = dot_product( u_cfd, e_ax  )
       vtan(n) = dot_product( u_cfd, e_tan )

    END DO
      
END SUBROUTINE velocity_ax_tan


SUBROUTINE velocity_kart( fast_position_blade, shaft_coordinates,                                  &
                          fast_n_blade_elem_max, fast_n_blades_max, nturbines,                     &
                          i,                                                                       &
                          j,                                                                       &
                          vel_ax,                                                                  &
                          vel_tan,                                                                 &
                          palm_vel_blade )

    USE kinds, ONLY : wp, iwp
    
    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  fast_n_blade_elem_max
    INTEGER(iwp), INTENT(IN) ::  fast_n_blades_max
    INTEGER(iwp), INTENT(IN) ::  nturbines
    INTEGER(iwp), INTENT(IN) ::  i
    INTEGER(iwp), INTENT(IN) ::  j

    REAL(wp), DIMENSION(fast_n_blade_elem_max, fast_n_blades_max, nturbines, 1:3), INTENT(IN) ::  fast_position_blade
    REAL(wp), DIMENSION(nturbines, 3), INTENT(IN) ::  shaft_coordinates

    REAL(wp), DIMENSION(fast_n_blade_elem_max), INTENT(IN) ::  vel_ax
    REAL(wp), DIMENSION(fast_n_blade_elem_max), INTENT(IN) ::  vel_tan

    REAL(wp), DIMENSION(fast_n_blade_elem_max, fast_n_blades_max, nturbines, 1:3), INTENT(INOUT) ::  palm_vel_blade

    REAL(wp), DIMENSION(3) ::  a
    REAL(wp), DIMENSION(3) ::  e_ax
    REAL(wp), DIMENSION(3) ::  e_rad
    REAL(wp), DIMENSION(3) ::  e_tan
    REAL(wp), DIMENSION(3) ::  u_old
    REAL(wp) ::  norm_axis
    REAL(wp) ::  norm_rad
    REAL(wp) ::  norm_tan
    REAL(wp) ::  vel_rad_old

    INTEGER(iwp) ::  n

!-- Normalize the FAST shaft axis. The fallback preserves the old untilted
!-- behavior if an invalid zero vector is ever received.

    norm_axis = sqrt( sum( shaft_coordinates(i,1:3)**2 ) )
    IF ( norm_axis > epsilon( norm_axis ) ) THEN
       e_ax = shaft_coordinates(i,1:3) / norm_axis
    ELSE
       e_ax = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    END IF

!-- Reconstruct Cartesian velocities from the corrected axial/tangential
!-- components while preserving the original radial component in the local
!-- rotor-plane basis.

    DO n = 1, fast_n_blade_elem_max

       u_old = palm_vel_blade(n, j, i, 1:3)

       ! Direction vector of the rotor blade (assumption: no bending of the rotor blade).
       a = fast_position_blade(3, j, i, 1:3) - fast_position_blade(1, j, i, 1:3)

       e_rad = a - dot_product( a, e_ax ) * e_ax
       norm_rad = sqrt( sum( e_rad**2 ) )
       IF ( norm_rad > epsilon( norm_rad ) ) THEN
          e_rad = e_rad / norm_rad
       ELSE
          ! Degenerate blade-axis data: choose a deterministic vector perpendicular
          ! to the shaft, instead of silently using a non-orthogonal fallback.
          IF ( abs( e_ax(1) ) < 0.9_wp ) THEN
             e_rad = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
          ELSE
             e_rad = (/ 0.0_wp, 1.0_wp, 0.0_wp /)
          END IF
          e_rad = e_rad - dot_product( e_rad, e_ax ) * e_ax
          e_rad = e_rad / sqrt( sum( e_rad**2 ) )
       END IF

       e_tan(1) = e_ax(2) * e_rad(3) - e_ax(3) * e_rad(2)
       e_tan(2) = e_ax(3) * e_rad(1) - e_ax(1) * e_rad(3)
       e_tan(3) = e_ax(1) * e_rad(2) - e_ax(2) * e_rad(1)
       norm_tan = sqrt( sum( e_tan**2 ) )
       IF ( norm_tan > epsilon( norm_tan ) ) THEN
          e_tan = e_tan / norm_tan
       END IF

       vel_rad_old = dot_product( u_old, e_rad )
       palm_vel_blade(n, j, i, 1:3) = vel_ax(n)  * e_ax  +                                  &
                                      vel_tan(n) * e_tan +                                  &
                                      vel_rad_old * e_rad

    END DO
      
END SUBROUTINE velocity_kart
