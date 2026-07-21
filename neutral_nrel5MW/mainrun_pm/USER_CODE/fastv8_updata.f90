!> @file fastv8_updata.f90
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
!> Additional routines used in C Codebase of FASTv8 Coupler
!--------------------------------------------------------------------------------------------------!
 MODULE fastv8_updata

    USE fastv8_coupler_mod
    USE kinds
   
    IMPLICIT NONE

    PRIVATE

!
!- Public functions
    PUBLIC                                                                                         &
       palm_bld_data,                                                                              &
       palm_sim_env,                                                                               &
       palm_turb_par,                                                                              &
       palm_blade_pitch,                                                                           &
       palm_force_table_bld,                                                                       &
       palm_position_table,                                                                        &
       palm_vel_value,                                                                             &
       palm_nwc_set_counts,                                                                         &
       palm_nwc_set_node,                                                                           &
       palm_nwc_set_polar_len,                                                                      &
       palm_nwc_set_polar_row,                                                                      &
       palm_nwc_set_done


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates PALM simulation/environment variables based on FAST information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_sim_env(turbine_id, dt_fast_in, shaft_height_fast_in)                             &
            BIND(C,NAME="palm_sim_env")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  dt_fast_in            !< time steps size FAST
    REAL(kind=C_DOUBLE), INTENT(IN) ::  shaft_height_fast_in  !< shaft height from FAST

    dt_fast = dt_fast_in
    shaft_height_fast(turbine_id) = shaft_height_fast_in

 END SUBROUTINE palm_sim_env


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates PALM turbine parameters based on FAST information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_turb_par(turbine_id, simtimef, vel_rot, xs_x, xs_y, xs_z)                         &
            BIND(C,NAME="palm_turb_par")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  simtimef  !< current FAST simulation time
    REAL(kind=C_DOUBLE), INTENT(IN) ::  vel_rot   !< rotational velocity from FAST
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_x      !< shaft coordinate in X direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_y      !< shaft coordinate in Y direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_z      !< shaft coordinate in Z direction

#ifdef __VERBOSE_MODE__
    INTEGER(kind=C_INT), SAVE :: n_calls = 0
#endif   

    rotspeed(turbine_id) = vel_rot
    current_time_fast(turbine_id) = simtimef
    
    shaft_coordinates(turbine_id,1) = xs_x
    shaft_coordinates(turbine_id,2) = xs_y
    shaft_coordinates(turbine_id,3) = xs_z

    !> s.ouchene: add __VERBOSE_MODE for debugging
#ifdef __VERBOSE_MODE__  
    n_calls = n_calls + 1
    print *, "n_calls:", n_calls, "|turbine id:", turbine_id,"|current FAST time ", simtimef,    &
      "|shaft_coordinates (xs_x, xs_y, xs_z) = (", xs_x, xs_y, xs_z, ")", "|rot. speed:", vel_rot
#endif

 END SUBROUTINE palm_turb_par


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Stores the blade pitch sent by FAST for one blade of one turbine.
!> The value is in radians and is read once per coupling step. The near-wake
!> correction converts it to degrees where it is used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_blade_pitch(turbine_id, blade_id, pitch_in)                                       &
            BIND(C,NAME="palm_blade_pitch")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  blade_id    !< blade index
    REAL(kind=C_DOUBLE), INTENT(IN) ::  pitch_in    !< blade pitch from FAST (radians)

    blade_pitch(turbine_id, blade_id) = pitch_in

 END SUBROUTINE palm_blade_pitch


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element and hub positions into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_position_table(turbine_id, max_blade_elem_id, blade_elem_id,            &
                                          pos_x, pos_y, pos_z, error_status)                       &
                      BIND(C,NAME="palm_position_table")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id      !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  max_blade_elem_id  !<
    INTEGER(kind=C_INT), INTENT(IN) ::  blade_elem_id     !<

    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_x  !< input information: positions X
    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_y  !< input information: positions Y
    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_z  !< input information: positions Z

    INTEGER(kind=C_INT), INTENT(OUT) :: error_status  !< determines if an error has occured

    INTEGER ::  cflag
    INTEGER ::  fast_blade_elem_idx
 
    cflag = 0
    error_status = 0
!
!-- check which blade the current blade element is
    IF ( blade_elem_id <= (max_blade_elem_id-1) / 3 )  THEN
       cflag = 1
    ELSEIF ( ( blade_elem_id > (max_blade_elem_id-1) / 3 )  .AND.                                  &
             ( blade_elem_id <= (max_blade_elem_id-1) / 3 * 2 ) )  THEN
       cflag = 2
    ELSEIF ( ( blade_elem_id > (max_blade_elem_id-1) / 3 * 2 )  .AND.                              &
             ( blade_elem_id <= (max_blade_elem_id-1) ) )  THEN
       cflag = 3
!
!-- the last element is the hub position
    ELSEIF ( blade_elem_id == max_blade_elem_id )  THEN
       cflag = 4
    ELSE
       error_status = 1
    END IF
!
!-- update blade element positions
    IF ( (cflag >= 1)  .AND.  (cflag <= 3) )  THEN

       fast_blade_elem_idx = blade_elem_id - (cflag-1)*((max_blade_elem_id-1)/3)

       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 1) = pos_x
       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 2) = pos_y
       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 3) = pos_z
!
!-- update hub position
    ELSEIF ( cflag == 4 )  THEN
       fast_hub_center_pos(turbine_id, 1) = pos_x
       fast_hub_center_pos(turbine_id, 2) = pos_y
       fast_hub_center_pos(turbine_id, 3) = pos_z

    ELSE
       error_status = 1
    END IF

 END SUBROUTINE palm_position_table


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element forces into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_force_table_bld(turbine_id, max_blade_elem_id, blade_elem_id,           &
                                           force_x, force_y, force_z, error_status)                &
                      BIND(C,NAME="palm_force_table_bld")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id      !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  max_blade_elem_id  !< 
    INTEGER(kind=C_INT), INTENT(IN) ::  blade_elem_id     !< blade element ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_x  !< input information: force in X direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_y  !< input information: force in Y direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_z  !< input information: force in Z direction

    INTEGER(kind=C_INT), INTENT(OUT)  :: error_status  !< determines if an error has occured

    INTEGER ::  cflag
    INTEGER ::  fast_blade_elem_idx
 
    cflag = 0
    error_status = 0
!
!-- check which blade the current blade element is
    IF ( blade_elem_id <= (max_blade_elem_id-1) / 3 )  THEN
       cflag = 1
    ELSEIF ( ( blade_elem_id >  (max_blade_elem_id-1) / 3 )  .AND.                                 &
             ( blade_elem_id <= (max_blade_elem_id-1) / 3 * 2 ) )  THEN
       cflag = 2
    ELSEIF ( ( blade_elem_id >  (max_blade_elem_id-1) / 3 * 2 )  .AND.                             &
             ( blade_elem_id <= (max_blade_elem_id-1) ) )  THEN
       cflag = 3
    ELSE
       error_status = 1
    END IF
!
!-- update blade element forces
    IF ( (cflag >= 1)  .AND.  (cflag <= 3) )  THEN

       fast_blade_elem_idx = blade_elem_id - (cflag-1)*((max_blade_elem_id-1)/3)

       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 1) = force_x
       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 2) = force_y
       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 3) = force_z

    ELSE
       error_status = 1
    END IF

 END SUBROUTINE palm_force_table_bld


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element data into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_bld_data(turbine_id, numbld, numbldelem, rotrad)                                  &
            BIND(C,NAME="palm_bld_data")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  numbld      !< number of blades
    INTEGER(kind=C_INT), INTENT(IN) ::  numbldelem  !< number of blade elements

    REAL(kind=C_DOUBLE), INTENT(IN) ::  rotrad      !< rotor radius (length of blades)

    fast_n_blades(turbine_id) =     numbld
    fast_n_blade_elem(turbine_id) = numbldelem
    fast_n_radius(turbine_id) =     rotrad

 END SUBROUTINE palm_bld_data


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the velocities for fast out of the PALM data structures.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_vel_value(turbine_id, bladeid, elemid,                                  &
                                     compu, compv, compw, error_status)                            &
                      BIND(C,NAME="palm_vel_value")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN)  ::  turbine_id  !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN)  ::  bladeid     !< blade ID (index)
    INTEGER(kind=C_INT), INTENT(IN)  ::  elemid      !< blade element ID (index)

    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compu  !< wind velocity component in X direction
    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compv  !< wind velocity component in Y direction
    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compw  !< wind velocity component in Z direction

    INTEGER(kind=C_INT), INTENT(OUT) ::  error_status  !< determines if an error has been encountered

    error_status = 0

#ifdef __VERBOSE_MODE
    IF ( bladeid == 0 .AND. elemid == 0 )  THEN
       print *, 'VEL_VALUE_HUB_ENTRY turb=', turbine_id,                  &
                ' palm_hub_center_vel=',                                   &
                palm_hub_center_vel(turbine_id,1),                         &
                palm_hub_center_vel(turbine_id,2),                         &
                palm_hub_center_vel(turbine_id,3)
    END IF
#endif

    IF ( (turbine_id > 0)  .AND.  (turbine_id <= nturbines) )  THEN
      
       IF ( (bladeid > 0)  .AND.  (bladeid <= fast_n_blades(turbine_id)) )  THEN

          IF ( (elemid > 0)  .AND.  (elemid <= fast_n_blade_elem(turbine_id)) )  THEN 
             compu = palm_vel_blade(elemid, bladeid, turbine_id, 1)
             compv = palm_vel_blade(elemid, bladeid, turbine_id, 2)
             compw = palm_vel_blade(elemid, bladeid, turbine_id, 3)
          ELSE
             error_status = 1
          ENDIF
          
       ELSEIF ( (bladeid == 0)  .AND.  (elemid == 0) )  THEN
          compu = palm_hub_center_vel(turbine_id, 1)
          compv = palm_hub_center_vel(turbine_id, 2)
          compw = palm_hub_center_vel(turbine_id, 3)
#ifdef __VERBOSE_MODE
print *, 'VEL_VALUE_HUB turb=', turbine_id, ' compu,v,w=', compu, compv, compw
#endif
       ELSE
          error_status = 1
       ENDIF

    ELSE
       error_status = 1
    ENDIF

 END SUBROUTINE palm_vel_value


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Near-wake-correction (NWC) blade geometry and airfoil polars received from
!> OpenFAST (AD15) at the first coupling step. These setters store the raw data
!> into the fastv8_coupler_mod module arrays; f8c_build_nwc_arrays later turns it
!> into the arrays the NWC consumes. The data is shared across turbines (one blade
!> definition for the farm), so the turbine_id argument is accepted for interface
!> symmetry but not used to index the arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_nwc_set_counts(turbine_id, nnodes, nairfoils, max_nrows)                          &
            BIND(C,NAME="palm_nwc_set_counts")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID (unused; data is shared)
    INTEGER(kind=C_INT), INTENT(IN) ::  nnodes      !< number of blade nodes
    INTEGER(kind=C_INT), INTENT(IN) ::  nairfoils   !< number of airfoil tables
    INTEGER(kind=C_INT), INTENT(IN) ::  max_nrows   !< longest airfoil table

    nwc_recv_nnodes    = nnodes
    nwc_recv_nairfoils = nairfoils
    nwc_recv_max_nrows = max_nrows

!-- Allocate the raw receive buffers once. Every turbine sends the same blade, so
!-- later connections reuse the existing allocation and overwrite identical data.
    IF ( .NOT. ALLOCATED( nwc_recv_r_cp ) )  THEN
       ALLOCATE( nwc_recv_r_cp(nnodes), nwc_recv_twist_deg(nnodes),                                &
                 nwc_recv_chord(nnodes), nwc_recv_blafid(nnodes) )
       ALLOCATE( nwc_recv_nrows(nairfoils) )
       ALLOCATE( nwc_recv_alpha(max_nrows, nairfoils),                                             &
                 nwc_recv_cl(max_nrows, nairfoils),                                                &
                 nwc_recv_cd(max_nrows, nairfoils) )
       nwc_recv_nrows = 0_iwp
    END IF

 END SUBROUTINE palm_nwc_set_counts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store one blade node's geometry (radial position, twist in degrees, chord,
!> 1-based airfoil id) for the near-wake correction.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_nwc_set_node(turbine_id, node_id, r_cp, twist_deg, chord, blafid)                 &
            BIND(C,NAME="palm_nwc_set_node")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID (unused; data is shared)
    INTEGER(kind=C_INT), INTENT(IN) ::  node_id     !< 1-based blade-node index
    INTEGER(kind=C_INT), INTENT(IN) ::  blafid      !< 1-based airfoil id
    REAL(kind=C_DOUBLE), INTENT(IN) ::  r_cp        !< radial position from rotor centre (m)
    REAL(kind=C_DOUBLE), INTENT(IN) ::  twist_deg   !< aerodynamic twist (deg)
    REAL(kind=C_DOUBLE), INTENT(IN) ::  chord       !< chord (m)

    nwc_recv_r_cp(node_id)      = r_cp
    nwc_recv_twist_deg(node_id) = twist_deg
    nwc_recv_chord(node_id)     = chord
    nwc_recv_blafid(node_id)    = blafid

 END SUBROUTINE palm_nwc_set_node


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store the row count of one airfoil polar table.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_nwc_set_polar_len(turbine_id, airfoil_id, nrows)                                  &
            BIND(C,NAME="palm_nwc_set_polar_len")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID (unused; data is shared)
    INTEGER(kind=C_INT), INTENT(IN) ::  airfoil_id  !< 1-based airfoil index
    INTEGER(kind=C_INT), INTENT(IN) ::  nrows       !< number of rows in this table

    nwc_recv_nrows(airfoil_id) = nrows

 END SUBROUTINE palm_nwc_set_polar_len


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store one row (alpha in degrees, Cl, Cd) of one airfoil polar table.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_nwc_set_polar_row(turbine_id, airfoil_id, row_id, alpha_deg, cl, cd)              &
            BIND(C,NAME="palm_nwc_set_polar_row")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID (unused; data is shared)
    INTEGER(kind=C_INT), INTENT(IN) ::  airfoil_id  !< 1-based airfoil index
    INTEGER(kind=C_INT), INTENT(IN) ::  row_id      !< 1-based row index
    REAL(kind=C_DOUBLE), INTENT(IN) ::  alpha_deg   !< angle of attack (deg)
    REAL(kind=C_DOUBLE), INTENT(IN) ::  cl          !< lift coefficient
    REAL(kind=C_DOUBLE), INTENT(IN) ::  cd          !< drag coefficient

    nwc_recv_alpha(row_id, airfoil_id) = alpha_deg
    nwc_recv_cl(row_id, airfoil_id)    = cl
    nwc_recv_cd(row_id, airfoil_id)    = cd

 END SUBROUTINE palm_nwc_set_polar_row


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Mark the near-wake-correction geometry transfer as complete. The coupler uses
!> this flag to build its NWC arrays at the first coupling step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_nwc_set_done(turbine_id)                                                          &
            BIND(C,NAME="palm_nwc_set_done")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID (unused; data is shared)

    nwc_geom_received = .TRUE.

 END SUBROUTINE palm_nwc_set_done

 END MODULE fastv8_updata
