!-----------------------------------------------------------------------
!            ####  #######  #    #     #     #  #  #    #  ####
!            #   #    #     #    #     #  #  #  #  # #  #  #   #
!            #   #    #     #    #     #  #  #  #  #  # #  #   #
!            ####     #      ####       #   #   #  #    #  ####
!-----------------------------------------------------------------------
!
!     A fast SMEARING CORRECTION for the actuator line 
!
!     Author:        Alexander R. Meyer Forsting, alrf@dtu.dk
!     Contributors:  Georg Pirrung,
!                    Néstor Ramos-García,
!                    Ang Li 
!     Revisions : 
!         date     |      author     |    task
!     - 24-07-2019 | A MeyerForsting | Version 1.0
!     - 09-09-2025 | Samir Ouchene   | Version 1.1
!
!     Changes by S.Ouchene:
!         - Format the code: use correct spacing, indentation
!         - replace the use of d0 (as in 0.0d0) with the suffix '_wp'
!
!=========================================================================
!     AD15 endpoint handling (root and tip control points)
!=========================================================================
!     With the AD15 node layout the first and last control points lie on
!     the root and tip trailed vortices (r_vtp(1) = r_cp(1) and
!     r_vtp(n_v) = r_cp(n_s)), so the distance h = r_vtp(k) - r_cp(l) is
!     zero there and the 1/(|h|*h_r) amplitude in nw_hr_components is
!     singular. The h == 0 branches below zero the kernel terms to avoid
!     Inf/NaN. The endpoint correction velocity is then set by linear
!     extrapolation from the two adjacent interior nodes, controlled by
!     nwc_extrapolate_endpoint_correction.
!=========================================================================
!=========================================================================
!
!     Purpose of file: Correct actuator line forces for the missing
!                      induction, originating from the force smearing at
!                      the control points. This code is an
!                      implementation of the correction presented in
!                      Meyer Forsting et al. (2019).
!
!     Overview of functions contained in this file:
!     - module nw_smearing_correction_module
!     - subroutine nw_smearing_correction()
!     - subroutine smearing_factor()
!     - subroutine nw_hr_components()
!     - subroutine nw_phi_components()
!     - subroutine nw_XY_components()
!     - subroutine nw_circulation()
!     - subroutine nw_calc_relax()
!     - subroutine init_nw_smearing_correction()
!     - subroutine nw_interp1d()
!
!-----------------------------------------------------------------------
!     This file contains all elements necessary for running the smearing
!     correction for the actuator line (AL). The AL code only needs to 
!     invoke (call) nw_smearing_correction() with the given in- and output 
!     variables. Note that it needs to be called for each blade 
!     individually. Call init_* before calling the routine to initialize 
!     all variables. Activate the correction by setting smearing_correction_active=.true. and
!     set a convgerence limit (nw_max_resid, 1d-6 for example). When 
!     using deformation some constants might have to be recomputed which 
!     can be achieved through setting recompute_smearing_constants of the specific rotor
!     and blade to .true.. Variables from previous time steps are 
!     stored for all rotors and blades in single arrays (variables 
!     ending on *_history). 
!
!     This implementation is based on several scientific publications
!     related to the smearing correction and the near-wake model, which
!     are referenced throughout the code:
!
!     Meyer Forsting et al. (2019)
!        Meyer Forsting, A. R., Pirrung, G. R., and Ramos-García, N.: 
!        A vortex-based tip/smearing correction for the actuator line, 
!        Wind Energ. Sci., 4, 369-383, 
!        https://doi.org/10.5194/wes-4-369-2019, 2019. 
!
!     Pirrung et al. (2016)
!        Pirrung, G., Madsen, H. A., Kim, T., and Heinz, J.: 
!        A coupled near and far wake model for wind turbine aerodynamics
!        , Wind Energy, 19, 2053–2069, 
!        https://doi.org/10.1002/we.1969, 2016. 
!
!     Pirrung et al. (2017a)
!        Pirrung, G., Riziotis, V., Madsen, H., Hansen, M., and Kim, T.: 
!        Comparison of a coupled near- and far-wake model with a free-
!        wake vortex code, Wind Energ. Sci., 2, 15–33, 
!        https://doi.org/10.5194/wes-2-15-2017, 2017. 
!
!     Pirrung et al. (2017b)
!     Pirrung, G. R., Madsen, H. A., and Schreck, S.: 
!        Trailed vorticity modeling for aeroelastic wind turbine 
!        simulations in standstill, Wind Energ. Sci., 2, 521–532, 
!        https://doi.org/10.5194/wes-2-521-2017, 2017. 
!
!-----------------------------------------------------------------------


!=======================================================================
!      module nw_smearing_correction_module 
!-----------------------------------------------------------------------
!     Variables used in the new near wake tip correction
!----------------------------------------------------------------------- 

      
      ! Global precision
!      ! Switch activating the correction
!      logical::smearing_correction_active=.true.
!      !TEilers: added smearing_correction_needs_init
!      logical::smearing_correction_needs_init=.true.
!      ! Switch invocking the recomputation of constants (large changes 
!      ! in the blade geometry from deformation could make it necessary)
!      logical,dimension(:,:),allocatable::recompute_smearing_constants
!      ! Definition of setup: number of rotors, nr. of blades, nr. of 
!      ! vortex trailing points (in-between blade sections), nr. of 
!      ! sections along blade
!      integer::n_rotors=0,n_blades=0,n_v=0,n_s=0
!      ! Index range of control points (CP): The missing induction is 
!      ! only felt in the proximity of the shed vortices, allowing to 
!      ! limit the number of CPs to loop over.
!      integer,dimension(:,:,:,:),allocatable::cp_loop_ind 
!      ! Pi
!      real(kind=idp2)::pi=4.d0*atan(1.d0)
!      ! Convergence limit of induced velocity residual (default value)
!      real(kind=idp2)::nw_max_resid=1d-3!nw_max_resid=1d-6
!      ! Constants of the Near-wake model
!      real(kind=idp2),dimension(4,5)::nw_N,nw_P
!      ! Variables from previous time step: For smooth transition between
!      ! steps history needs to be stored
!      real(kind=idp2),dimension(:,:,:),allocatable:: nw_vn_history,nw_vt_history,nw_Gamma_history,nw_W_history          
!      real(kind=idp2),dimension(:,:,:,:),allocatable:: nw_dyn_comp_history,nw_X_history,nw_Y_history
!      ! Constants that need to be recomputed only after strong changes 
!      ! in geometry
!      real(kind=idp2),dimension(:,:,:,:),allocatable:: nw_phi,nw_phi_s,nw_dx,nw_dy,smear_fac



!      real(kind=idp2),dimension(:,:,:,:,:),allocatable:: nw_a   
!      end module 
!=======================================================================
!=======================================================================
SUBROUTINE nw_smearing_correction(i, j,     &
    r_cp, r_vtp, vn_cfd, vt_cfd,            &
    eps_cp, omega, chord, twist, pitch,     &
    cl_2d, alpha_2d, len_2d,                &
    dt_palm,                                &
    nw_trailed_vorticity_smoothing_length,  &
    vn_cfd_out, vt_cfd_out)
  !------------------------------------------------------
  !     The smearing function leads to a loss of vorticity,
  !     which in turn leads to lower induced velocities. The
  !     missing velocity is calculated by a near wake model
  !     and used to correct the angle of attack.
  !-----------------------------------------------------------------------



  USE nw_smearing_correction_module
  USE kinds, ONLY: wp   ! nw_smearing_correction_module is already using wp

  !-----------------------------------------------------------------------

  implicit none !<<< s.ouchene: Force implicit none

  !-----[in]

  ! --------- [In INTEGER variables] -------------
  INTEGER, INTENT(IN) ::         &
    i,                           &   ! ith rotor
    j,                           &   ! jth blade
    len_2d                           ! Length of 2D airfoil data

  !------ [IN REAL variables ] -------------------
  REAL(KIND=wp), INTENT(IN) ::   &
    r_cp(n_s),                   &   ! Radial coordinate of control points (CP)
    r_vtp(n_v),                  &   ! Radial coordinate of vortex trailing points
    vn_cfd(n_s),                 &   ! Normal velocities at CPs from CFD domain
    vt_cfd(n_s),                 &   ! Tangential velocities at CPs from CFD domain
    eps_cp,                      &   ! Epsilon, the smearing length scale
    omega,                       &   ! Rotational speed in radians/sec
    chord(n_s),                  &   ! Chord distribution along blade at CPs
    twist(n_s),                  &   ! Twist distribution along blade at CPs
    pitch,                       &   ! Blade pitch
    cl_2d(len_2d,n_s),           &   ! Airfoil lift coeff. curve at CPs
    alpha_2d(len_2d,n_s),        &   ! Alpha values of lift&drag curves
    dt_palm,                     &   ! PALM time step (interval between successive calls)
    nw_trailed_vorticity_smoothing_length  ! S. Ouchene (May 2026): spanwise smoothing length L (m) for the Gamma kernel. 0.0 = kernel off.

  ! Note that the radial coordiante is defined along the blade and
  ! starts at the rotor centre. The radial coordinate should be 
  ! ascending and r_vtp(1)>0.

  !-----[Out REAL variables] -----------------------
  REAL(KIND=wp), DIMENSION(n_s), INTENT(OUT) :: vn_cfd_out,vt_cfd_out

  !-----[tmp]
  LOGICAL :: nw_recompute_relax,nw_converged
  INTEGER :: nw_sub_iter, k, l
  REAL(KIND=wp) :: dt_global, h, h_r, nw_resmax, nw_relax, &
    nw_relax_safety, core_cutoff, &
    !?rtc,
  cpu_t1, cpu_t2      
  REAL(KIND=wp), DIMENSION(n_s) ::    &
    Gamma_cp, Gamma_cp_steady, phi_cp_cfd, &
    nw_section_res, nw_vn, nw_vt, &
    nw_W, phi_cp, alpha_cp, vrel_cp, nw_W_iter
  
  REAL(KIND=wp), DIMENSION(n_v) ::    &
    phi_vtp, vrel_vtp, dGamma_vtp, nw_dbeta
  !TEilers added nw_phi_b
  
  REAL(KIND=wp) :: nw_dyn_comp(n_s,3),  &
    nw_X(n_v,n_s), nw_Y(n_v,n_s), nw_phi_star(n_v,n_s), nw_phi_b(n_v,n_s), &
    nw_dx_b(n_v,n_s), nw_dy_b(n_v,n_s)

  INTEGER :: count , count2 !!! s.ouchene: count=0, count2= 0 implies SAVE attribute.
  ! Initialization of count and count2 is moved to separate lines

  INTEGER :: counter, count_rate, count_max

  !=======================================================================
  !write(*,*)'Start Smearing Correction'

  ! S.ouchene: Moved count = 0 and count2 = 0 here. Initially they were initialized
  ! upon declaration which is a common pitfall in Fortran( unless intended). In this case
  ! there is an implied SAVE statement which means the local variables will be initialized
  ! only during the first call and they will preserve their values for future calls.
  count = 0
  count2 = 0

  !----[ wt-output diagnostics: pre-correction snapshot ]----
  ! Capture the raw NWC inputs and the polar-interpolated alpha / cl / cd
  ! they imply, before any wake correction is applied. Used by
  ! wt_output_mod (read directly from nw_smearing_correction_module).
  ! Phi is wrapped to <= pi/2 to match what the downstream NWC loop
  ! actually consumes (see the analogous wrap below at line ~347).
  nwc_diag_vn_pre   (i,j,:) = vn_cfd
  nwc_diag_vt_pre   (i,j,:) = vt_cfd
  nwc_diag_phi_pre  (i,j,:) = ATAN2( vn_cfd, omega * r_cp - vt_cfd )
  DO l = 1, n_s
    IF ( nwc_diag_phi_pre(i,j,l) > pi/2.0_wp )                                  &
      nwc_diag_phi_pre(i,j,l) = pi - nwc_diag_phi_pre(i,j,l)
  END DO
  nwc_diag_alpha_pre(i,j,:) = nwc_diag_phi_pre(i,j,:) * 180.0_wp / pi - twist - pitch
  DO l = 1, n_s
    CALL nw_interp1d( alpha_2d(:,l), cl_2d(:,l),     len_2d,                    &
                      nwc_diag_alpha_pre(i,j,l:l), nwc_diag_cl_pre(i,j,l:l), 1 )
    CALL nw_interp1d( alpha_2d(:,l), cd_2d_arr(:,l), len_2d,                    &
                      nwc_diag_alpha_pre(i,j,l:l), nwc_diag_cd_pre(i,j,l:l), 1 )
  END DO

  ! Default post-correction quantities to the pre-correction values and zero
  ! induced velocities. Overwritten at the end of the routine if the wake-model
  ! loop runs to convergence; the early-exit paths below leave these in place.
  nwc_diag_vn_post   (i,j,:) = nwc_diag_vn_pre   (i,j,:)
  nwc_diag_vt_post   (i,j,:) = nwc_diag_vt_pre   (i,j,:)
  nwc_diag_phi_post  (i,j,:) = nwc_diag_phi_pre  (i,j,:)
  nwc_diag_alpha_post(i,j,:) = nwc_diag_alpha_pre(i,j,:)
  nwc_diag_cl_post   (i,j,:) = nwc_diag_cl_pre   (i,j,:)
  nwc_diag_cd_post   (i,j,:) = nwc_diag_cd_pre   (i,j,:)
  nwc_diag_nw_vn     (i,j,:) = 0.0_wp
  nwc_diag_nw_vt     (i,j,:) = 0.0_wp
  nwc_diag_nw_W      (i,j,:) = 0.0_wp


  IF (.NOT. smearing_correction_active) THEN
    ! Correction disabled: pass the input velocities through so callers never
    ! consume undefined output arrays.
    vn_cfd_out = vn_cfd
    vt_cfd_out = vt_cfd
    RETURN
  END IF

  ! Test whether the first coordinate is close to zero. This is a
  ! singularity.
  IF (r_vtp(1) < 0.1_wp ) THEN
     
    PRINT *, "!!!! WARNING smearing correction switched off !!!!"
    PRINT *, "First vortex trailing point is zero, please modify"
   
    smearing_correction_active =.FALSE.

    vn_cfd_out = vn_cfd
    vt_cfd_out = vt_cfd
    return
  END IF
  !#################### INITIALIZE ##################################   

  !========= Starting Preamble =====================================

  nw_relax_safety = 0.2_wp
  nw_recompute_relax = .TRUE.
  nw_converged = .FALSE.
  ! Get the global timestep of the current grid level
  ! dt_global = timestep0(grlev) 

!-- s.ouchene: wake model integration step = intervals between
!-- successive calls. In this coupling that is the PALM time step
!-- the call is gated by l == 1 in fastv8_coupler_mod
  dt_global = dt_palm
  ! Get velocities from the NW model from previous time step
  nw_W = nw_W_history(i,j,:)
  nw_vn =  nw_vn_history(i,j,:)
  nw_vt =  nw_vt_history(i,j,:)

  !+++++++++++++++++++++ ONLY IF REFRESH NEEDED ++++++++++++++++++++

  IF ( recompute_smearing_constants(i,j) ) THEN
    !============= Compute constant nw components ==================
    ! These components only need to be computed once at the very 
    ! first iteration or when the geometry of the blade changes 
    ! drastically.  

    ! Viscous core cut-off at 99th percentile of Gaussian. Outside 
    ! this region the missing induction from the vortex smearing 
    ! is neglegible.
    core_cutoff = 1.82_wp * eps_cp   


    DO k = 1, n_v
    DO l = 1, n_s

    !---------------------------------------------------------
    ! Distance between vortex and control point
    h = r_vtp(k)-r_cp(l)

    ! Root/tip endpoint coincidence: the control point lies on its own
    ! trailed vortex (h = 0). Zero the kernel terms here to avoid the
    ! 1/(|h| * h_r) singularity. The endpoint correction value is set by
    ! extrapolation in the iteration loop below.
    IF (h == 0.0_wp) THEN
      nw_phi   (i,j,k,l)   = 0.0_wp
      nw_phi_s (i,j,k,l)   = 0.0_wp
      nw_dx    (i,j,k,l)   = 0.0_wp
      nw_dy    (i,j,k,l)   = 0.0_wp
      nw_a     (i,j,k,l,:) = 0.0_wp
      smear_fac(i,j,k,l)   = 0.0_wp
      CYCLE
    END IF

    h_r = h/r_vtp(k)

    ! Compute constants of the near-wake model


    call nw_hr_components(h,h_r,nw_phi(i,j,k,l),&
      nw_phi_s(i,j,k,l),nw_dx(i,j,k,l),&
      nw_dy(i,j,k,l),nw_a(i,j,k,l,:))
    ! Compute the smearing factor, which is fixed in this 
    ! implementation i.e. beta=0

    IF (nw_phi(i,j,k,l)==0) THEN

      write(*,*) 'Problem nw_phi' 
      write(*,*) 'h,k,l',h,k,l
      write(*,*) 'h_r',h_r
      write(*,*) 'nw_phi',nw_phi(i,j,k,l)
      !write(*,*) 'nw_phi_s',nw_phi_s(i,j,k,l)         

    END IF


    call smearing_factor(r_vtp(k),h_r,eps_cp,&
      smear_fac(i,j,k,l))
    ! Finally determine which blade sections are influenced
    ! by the vtpex core. For a Gaussian the 99th percentile
    ! is reached at 1.82x the smearing factor. Further from
    ! the influence of the core is assumed negligible.

    IF (h.gt.0.0_wp) THEN
      IF (abs(h).gt.core_cutoff) THEN
        cp_loop_ind(i,j,k,1) = l+1
      END IF
    ELSE
      IF (abs(h).lt.core_cutoff) THEN
        cp_loop_ind(i,j,k,2) = l
      END IF
    END IF
    !--------------------------------------------------------- 
    END DO
    END DO  
    ! Now it's fresh
    recompute_smearing_constants(i,j) =.FALSE.
  END IF ! refresh?    

  !######################### MAIN ##################################


  !?! cpu_t1=rtc() ! Start timing 

  CALL SYSTEM_CLOCK(counter, count_rate, count_max)
  cpu_t1 = REAL( counter, KIND=wp ) / REAL( count_rate, KIND=wp )
  !write(*,*)'CPU Startzeit', cpu_t1
  !===================  Fixed components =======================
  ! The physical flow angle does not change during a CFD
  ! iteration. Therefore there are no sub-iterations needed. 
  !----Flow/helix angle----
  ! Extracted from CFD at control points
  phi_cp_cfd = ATAN2( vn_cfd, omega * r_cp - vt_cfd ) 
  ! Ensure that the angle is below |90| degrees
  DO l = 1, n_s
    IF (phi_cp_cfd(l).gt.pi/2.0_wp) THEN 
      phi_cp_cfd(l)=pi-phi_cp_cfd(l) 
    END IF
  END DO


  ! Trailing points
  ! Interpolate from CPs

  CALL nw_interp1d( r_cp, phi_cp_cfd, n_s, r_vtp(2:n_v-1), phi_vtp(2:n_v-1), n_s-1 ) 

  phi_vtp(1) = phi_cp_cfd(1)
  phi_vtp(n_v) = phi_cp_cfd(n_s)

  !----Near Wake phi components---- 

  DO k = 1, n_v
    DO l =  cp_loop_ind(i,j,k,1), cp_loop_ind(i,j,k,2)
      ! Distance between trailing point and CP
      h = r_vtp(k) - r_cp(l)

      ! Endpoint coincidence (h = 0): nw_phi and nw_phi_s were already
      ! zeroed in the constants loop, so set nw_phi_star = 0 and skip.
      IF (h == 0.0_wp) THEN
        nw_phi_star(k,l) = 0.0_wp
        CYCLE
      END IF

      h_r = h / r_vtp(k)

      CALL nw_phi_components(h_r, nw_a(i,j,k,l,:), phi_vtp(k),    &
      nw_phi(i,j,k,l), nw_phi_s(i,j,k,l),                         &
      nw_phi_star(k,l))
      count = count + 1

    END DO ! l loop

  END DO ! k loop

  !************* Iterative procedure ***************************

  ! Iterate until correctinon velocities converge 
  nw_sub_iter = 0 ! Sub-iteration counter
  DO WHILE (.NOT. nw_converged)

    !============= Velocities & Flow angles  ===================
    !----Control points----
    ! Corrected flow angle
    phi_cp = ATAN2( vn_cfd + nw_vn, omega * r_cp - vt_cfd - nw_vt )


    ! Angle-of-attack at CP
    alpha_cp = phi_cp * 180.0_wp / pi - twist - pitch

    ! Relative wind speed at the control points
    vrel_cp = sqrt( ( vn_cfd + nw_vn )**2 + ( omega * r_cp - vt_cfd - nw_vt )**2 )
    !----Trailing points----
    ! Interpolate the velocity at the trailing points

    CALL nw_interp1d(r_cp, vrel_cp, n_s, r_vtp(2:n_v-1), vrel_vtp(2:n_v-1), n_s-1)
    ! Take velcoities from outer CPs (typical in LL codes)
    vrel_vtp(1) = vrel_cp(1)
    vrel_vtp(n_v) = vrel_cp(n_s)


    !================ Circulation ==============================

    ! Compute circulation at sections
    CALL nw_circulation( dt_global, vrel_cp, alpha_cp,          &
      nw_Gamma_history(i,j,:), nw_dyn_comp_history(i,j,:,:),    &
      alpha_2d, cl_2d, len_2d, chord,                           &
      Gamma_cp_steady, Gamma_cp, dGamma_vtp, nw_dyn_comp )

    ! S. Ouchene (May 2026): optional spanwise smoothing of Gamma to
    ! regularise the wake-model response at sharp Cl(alpha) steps such as
    ! the root cylinder -> airfoil transition. No-op when the namelist
    ! parameter nw_trailed_vorticity_smoothing_length (passed in from
    ! fastv8_coupler_mod via velocity_correction) is 0.0 (default).
    CALL smooth_trailed_vorticity( Gamma_cp, dGamma_vtp, r_cp,  &
                                   nw_trailed_vorticity_smoothing_length )

    !++++++++++++++++++ Velocity history +++++++++++++++++++++++
    ! Local advance in beta (rotation)
    nw_dbeta = ( vrel_vtp * dt_global ) / r_vtp          
    !============ Calculate new velocity =======================

    ! Determine the velocity induced by vorticity shed during  
    ! the current time step. 
    nw_W_iter = 0.0_wp ! Initialize
    nw_X = 0.0_wp 
    nw_Y = 0.0_wp

    DO k = 1, n_v
      DO l =  cp_loop_ind(i,j,k,1), cp_loop_ind(i,j,k,2)
        ! Distance between trailing point and CP
        h = r_vtp(k) - r_cp(l)

        ! Endpoint coincidence (h = 0): skip, otherwise nw_phi_star = 0
        ! would feed 1/nw_phi_star into the terms in nw_XY_components.
        IF (h == 0.0_wp) THEN
          nw_X(k,l) = 0.0_wp
          nw_Y(k,l) = 0.0_wp
          CYCLE
        END IF

        h_r = h / r_vtp(k)

        !----Near Wake velocity components----
        CALL nw_XY_components(nw_phi_star(k,l), nw_dx(i,j,k,l), &
          nw_dy(i,j,k,l), nw_dbeta(k), dGamma_vtp(k),           &
          nw_X_history(i,j,k,l), nw_Y_history(i,j,k,l),         &
          nw_X(k,l), nw_Y(k,l))

        !----Velocity correction----
        nw_W_iter(l) = nw_W_iter(l) + smear_fac(i,j,k,l) *      &
          ( nw_X(k,l) + nw_Y(k,l) )
      END DO
    END DO


    !============== Relaxation factor ==========================


    ! nw_calc_relax declares its nw_dx, nw_dy, nw_phi dummies as 2-D
    ! (n_v, n_s) explicit-shape. The module arrays are 4-D, so passing
    ! them bare would alias (turbine, blade) into the dummy's first
    ! dimension via storage association and feed nw_calc_relax data from
    ! the wrong blade. Pass per-(i, j) slice copies instead.
    nw_phi_b = nw_phi(i,j,:,:)
    nw_dx_b  = nw_dx(i,j,:,:)
    nw_dy_b  = nw_dy(i,j,:,:)


    IF (nw_recompute_relax) THEN
      CALL nw_calc_relax(vrel_cp,alpha_cp,                      &
        (omega * r_cp - vt_cfd - nw_vt), (vn_cfd + nw_vn),      &
        chord, dt_global, nw_dx_b, nw_dy_b, nw_dbeta, nw_phi_b, &
        nw_relax_safety, nw_relax )
      ! After computing it once 
      nw_recompute_relax = .FALSE.
    END IF  

    !============== Advance velocities =========================
    ! Advance velocity correction with relaxation
    nw_W_iter = nw_relax * nw_W + (1.0_wp - nw_relax ) * nw_W_iter

    ! The tip (l=n_s) and root (l=1) self-induction was zeroed above, which
    ! leaves those two nodes with only their inboard contributions. Set the
    ! endpoint correction by linear extrapolation from the two adjacent
    ! interior nodes instead. Done inside the loop so the endpoint feeds the
    ! bound circulation consistently.
    IF ( nwc_extrapolate_endpoint_correction .AND. n_s >= 3 ) THEN
      nw_W_iter(n_s) = 2.0_wp * nw_W_iter(n_s-1) - nw_W_iter(n_s-2)
      nw_W_iter(1)   = 2.0_wp * nw_W_iter(2)     - nw_W_iter(3)
    END IF

    !=================== Residual ==============================
    ! Tip vortex less than 45 degrees away from rotor plane: 
    ! convergence determined by axial induction otherwise 
    ! convergence determined by tangential induction  
    nw_section_res = SQRT( ( nw_W_iter - nw_W)**2 )

    IF ( phi_vtp(n_v) < pi/4.0_wp ) THEN 
      nw_resmax = MAXVAL(nw_section_res * COS(phi_cp))
    ELSE
      nw_resmax = MAXVAL(nw_section_res * SIN(phi_cp))
    END IF

    !bis hier

    !=============== Convergence ===============================
    nw_sub_iter = nw_sub_iter + 1

    !----Check----
    IF ( (nw_resmax < nw_max_resid) .OR. (nw_sub_iter == 500 ) ) THEN
      nw_converged = .TRUE.
    END IF

    !----100 Iterations---- 
    IF ( nw_sub_iter == 100 ) THEN
      nw_recompute_relax = .TRUE.
    END IF

    !----500 Iterations----
    ! No convergence after 500 sub-iterations: fall back to the
    ! previous time step's wake state.
    IF (nw_sub_iter == 500) THEN
      WRITE(*, '(A,X,I0,X,A)') 'NWC: no convergence after', nw_sub_iter, 'iterations.'
      nw_W = nw_W_history(i,j,:)
      nw_vn =  nw_vn_history(i,j,:)
      nw_vt =  nw_vt_history(i,j,:)
      nw_dyn_comp = nw_dyn_comp_history(i,j,:,:)
      Gamma_cp_steady = nw_Gamma_history(i,j,:)
      nw_X = nw_X_history(i,j,:,:) 
      nw_Y = nw_Y_history(i,j,:,:)

    ELSE
      !============== While loop update===========================
      nw_W = nw_W_iter 
      nw_vn = nw_W * COS(phi_cp)
      nw_vt = nw_W * SIN(phi_cp)


      count2 = count2 + 1

    END IF
  END DO 


  !********************* END ITERATING  ************************



  CALL SYSTEM_CLOCK(counter, count_rate, count_max)
  cpu_t2 = REAL( counter, KIND=wp ) / REAL( count_rate, KIND=wp )

  !=================== Output ==================================


  !write(150 + (j-1) * n_rotors + i, '(es15.7,i5,*(es15.7))')        &
  !  time, nw_sub_iter, nw_resmax, alpha_cp_cfd, alpha_cp, phi_vtp,  &
  !  nw_vt, nw_vn, dGamma_vtp, (cpu_t2 - cpu_t1)


  !=================== Saving iteration ========================
  nw_W_history(i,j,:) = nw_W
  nw_dyn_comp_history(i,j,:,:) = nw_dyn_comp
  nw_Gamma_history(i,j,:) = Gamma_cp_steady
  nw_vn_history(i,j,:) = nw_vn
  nw_vt_history(i,j,:) = nw_vt
  nw_X_history(i,j,:,:) = nw_X
  nw_Y_history(i,j,:,:) = nw_Y


  !Extra
  vn_cfd_out = vn_cfd + nw_vn
  vt_cfd_out = vt_cfd + nw_vt

  !----[ wt-output diagnostics: post-correction snapshot ]----
  ! Overwrite the pre-correction defaults set at routine entry with the
  ! actual converged post-correction state (last iteration of the DO WHILE).
  nwc_diag_vn_post   (i,j,:) = vn_cfd_out
  nwc_diag_vt_post   (i,j,:) = vt_cfd_out
  nwc_diag_phi_post  (i,j,:) = phi_cp
  nwc_diag_alpha_post(i,j,:) = alpha_cp
  DO l = 1, n_s
    CALL nw_interp1d( alpha_2d(:,l), cl_2d(:,l),     len_2d,                    &
                      nwc_diag_alpha_post(i,j,l:l), nwc_diag_cl_post(i,j,l:l), 1 )
    CALL nw_interp1d( alpha_2d(:,l), cd_2d_arr(:,l), len_2d,                    &
                      nwc_diag_alpha_post(i,j,l:l), nwc_diag_cd_post(i,j,l:l), 1 )
  END DO
  nwc_diag_nw_vn(i,j,:) = nw_vn
  nwc_diag_nw_vt(i,j,:) = nw_vt
  nwc_diag_nw_W (i,j,:) = nw_W

  !=============================================================
  nw_converged = .FALSE.


END SUBROUTINE nw_smearing_correction
!=======================================================================
!=======================================================================
SUBROUTINE smearing_factor(r, h_r, epsilon, smear_fac)   
  !-----------------------------------------------------------------------
  !     The Gaussian smearing results in a Lamb-Oseen like vortex. The 
  !     smearing in the induced velocity can directly be determined from 
  !     their analytical solution. 
  !     Definition and details in Meyer Forsting et al. (2019).
  !-----------------------------------------------------------------------
  USE kinds, ONLY: wp ! Global precision

  !-----------------------------------------------------------------------      
  IMPLICIT NONE

  !-----[in]        
  REAL(KIND=wp), INTENT(IN)  :: r, h_r, epsilon
  !-----[out]      
  REAL(KIND=wp), INTENT(OUT) :: smear_fac 
  !-----[tmp]      
  REAL(KIND=wp) :: dist_orth

  !=======================================================================      
  ! The definition of the smearing factor is given in Eq.(22) with 
  ! the perpendicular distance between vortex element and control 
  ! point defined in Eq.(21). Setting beta to zero, the expression 
  ! is greatly simplified. 
  ! Distance between CP and vortex trailing point
  dist_orth = (r * h_r)
  ! The smearing factor

  smear_fac = EXP( -dist_orth**2 / epsilon**2)

END SUBROUTINE smearing_factor
!=======================================================================
!=======================================================================
SUBROUTINE nw_hr_components(h, h_r, nw_phi, nw_phi_s, nw_dx, nw_dy, nw_a)
  !-----------------------------------------------------------------------      
  !     Compute the components of the near-wake model that only rely on
  !     geometric definiton of the blade. Unless the blade deforms heavily
  !     these components remain unchaged. Definitions of all equations and
  !     variables follow from Pirrung et al. (2016,2017b).
  !-----------------------------------------------------------------------      
  !      use kinds,only:wp ! Global precision


  USE kinds, ONLY: wp
  use nw_smearing_correction_module, ONLY : pi, nw_P, nw_N ! Constants
  !----------------------------------------------------------------------- 
  IMPLICIT NONE


  !-----[in] Distance between CP and vortex trailing point 
  REAL(KIND=wp), INTENT(IN)  ::  h_r, h
  !-----[out]      
  REAL(KIND=wp), INTENT(OUT) ::  nw_phi,nw_phi_s,nw_a(4),nw_dx,nw_dy
  !-----[tmp]      
  REAL(KIND=wp) ::  dw, root_corr
  INTEGER :: k
  !=======================================================================

  ! Pirrung et al. (2016)
  !------ Wang Coton -----------------
  ! Determine the model phi following Wang and Coton 
  ! Eq.(4), note in the paper the log is erroneously missing
  IF (h_r > 0.0_wp ) THEN
    nw_phi = pi/4.0_wp * ABS((1.0_wp + 0.5_wp * h_r) * LOG(1.0_wp - h_r))
  ELSE
    nw_phi = LOG (1.0_wp - h_r) / (1.5_wp + LOG( 1 - 0.5_wp * h_r) )
  END IF


  !------- Amplitude term -----------------
  ! The start value for the indicial function. Split into slow and 
  ! fast decaying parts. Gamma is added later. 
  ! Eq.(25) 
  dw = 1.0_wp / (4.0_wp * pi * ABS(h) * h_r )
  nw_dx = 1.359_wp * dw
  nw_dy = 0.359_wp / 4.0_wp * dw
  !-------- Root Correction ----------
  ! Implement the root correction that limits the induction to that 
  ! 90 degrees
  ! Eq.(27)
  root_corr = (nw_dx * (1.0_wp - exp( -pi/ (2.0_wp * nw_phi))) - nw_dy * & 
    (1.0_wp - exp(-2.0_wp * pi / nw_phi))) / ( nw_dx - nw_dy)      
  ! Apply root correction 
  nw_dx = root_corr * nw_dx
  nw_dy = root_corr * nw_dy
  nw_phi = root_corr * nw_phi

  ! Pirrung et al. (2017b)
  !----------- Convection correction -----------
  ! The wake can be in-plane solely convecting downstream and 
  ! anything in-between. For this purpose the phi needs to be 
  ! corrected for the convection of the vortex helix.
  ! Phi for straight vortices
  ! Eq.(7)
  IF (h_r < 0.0_wp ) THEN
    nw_phi_s = -0.788_wp * h_r
  ELSE
    nw_phi_s = 0.788_wp * h_r
  END IF 
  ! Determine the mixing ratio of straight and helical vorticity
  ! First compute the constants 
  ! Eq.(10) and (12)
  IF (h_r < 0.0_wp ) THEN
    DO k = 1, 4
    nw_a(k)=nw_N(k,1)+nw_N(k,2)*exp(nw_N(k,3)*h_r)&
      +nw_N(k,4)*exp(nw_N(k,5)*h_r)&
      -nw_N(k,4)-nw_N(k,2)
    END DO
  ELSE
    DO k = 1, 4
    nw_a(k) = nw_P(k,1) + nw_P(k,2) * h_r + nw_P(k,3) * h_r**2  &
      + nw_P(k,4) * h_r**3
    END DO
  END IF


END SUBROUTINE nw_hr_components
!=======================================================================
!=======================================================================
SUBROUTINE nw_phi_components(h_r, nw_a, phi, nw_phi, nw_phi_s,  &
    nw_phi_star)
  !-----------------------------------------------------------------------      
  !     Components of the near-wake model related to the inflow angle
  !     phi. 
  !-----------------------------------------------------------------------  
  USE kinds, ONLY: wp
  USE nw_smearing_correction_module, ONLY : pi 
  !-----------------------------------------------------------------------
  IMPLICIT NONE



  !-----[in]       
  REAL(KIND=wp), INTENT(IN) :: h_r, nw_a(4), phi, nw_phi, nw_phi_s
  !-----[out]     
  REAL(KIND=wp), INTENT(OUT) :: nw_phi_star
  !-----[tmp]      
  REAL(KIND=wp) :: nw_k_phi
  !=======================================================================  

  ! Pirrung et al. (2017b)     
  !----Mixing factor----
  ! Eq. (9) and (11)
  IF (h_r < 0.0_wp) THEN
    nw_k_phi = nw_a(1) + nw_a(2) * exp( nw_a(3) * ( pi / 2.0_wp - phi) ) + &
    nw_a(4) * exp( -8.0_wp * (pi / 2.0_wp - phi) ) - nw_a(2) - nw_a(4)
  ELSE
    nw_k_phi = nw_a(1) + nw_a(2) * phi + nw_a(3) * phi**2 + nw_a(4) * phi**3
  END IF
  ! The factor is only valid for angles below 89 degrees
  nw_k_phi = min(nw_k_phi, 1.0_wp)
  nw_k_phi = max(nw_k_phi, 0.0_wp)
  ! Determine the mix of straight and curved helix
  ! Eq. (8)
  nw_phi_star = nw_k_phi * nw_phi_s + (1.0_wp - nw_k_phi) * nw_phi

END SUBROUTINE nw_phi_components
!=======================================================================
!=======================================================================
SUBROUTINE nw_XY_components(nw_phi_star,nw_dx,nw_dy,dbeta,dGamma, &
    nw_X_history, nw_Y_history, &
    nw_X,nw_Y)
  !-----------------------------------------------------------------------      
  !     Induced velocity components of the near-wake model, with iterative
  !     dependancy through dbeta and dGamma. 
  !-----------------------------------------------------------------------  
  USE kinds, ONLY: wp
  !-----------------------------------------------------------------------        
  IMPLICIT NONE


  !-----[in]       
  REAL(KIND=wp), INTENT(IN) :: nw_phi_star, nw_dx, nw_dy, dbeta, dGamma, &
    nw_X_history, nw_Y_history  
  !-----[out]     
  REAL(KIND=wp), INTENT(OUT) :: nw_X, nw_Y
  !-----[tmp]      
  REAL(KIND=wp) :: nw_X_exp, nw_Y_exp, nw_X_new, nw_Y_new,               &
    nw_X_old, nw_Y_old
  !=======================================================================  


  !----Induced velocity components----
  ! Determine the slow and fast decaying parts 
  ! Pirrung et al. (2016) Eq.(25) with modified phi as described in 
  ! Pirrung et al. (2017b)

  ! Exponential terms
  nw_X_exp=exp(-dbeta/nw_phi_star)
  nw_Y_exp=exp(-4.0_wp*dbeta/nw_phi_star)
  ! New induction from shed vorticity
  nw_X_new = dGamma*nw_dx*nw_phi_star*(1.0_wp-nw_X_exp)
  nw_Y_new = dGamma*nw_dy*nw_phi_star*(1.0_wp-nw_Y_exp)
  ! Contribution from previously shed vortex elements,
  ! which need to be advanced in time
  nw_X_old = nw_X_history*nw_X_exp
  nw_Y_old = nw_Y_history*nw_Y_exp
  ! Total induction from vortex line, fast and slow 
  nw_X = nw_X_new + nw_X_old
  nw_Y = nw_Y_new + nw_Y_old


END SUBROUTINE nw_XY_components
!=======================================================================
!=======================================================================
SUBROUTINE nw_circulation(dt_global, vrel, alpha, Gamma_history,  &
    dyn_comp_history, &
    alpha_2d, cl_2d, len_2d, chord, &     
    Gamma_steady, Gamma, dGamma, dyn_comp)

  !-----------------------------------------------------------------------      
  !     Determine the circulation along the lifting line, including the 
  !     dynamic effects and also returns the delta in circulation. 
  !----------------------------------------------------------------------- 
  USE kinds, ONLY : wp ! Global precision
  USE nw_smearing_correction_module, ONLY: n_s, n_v ! Constants
  !-----------------------------------------------------------------------
  IMPLICIT NONE



  !-----[in]  
  INTEGER, INTENT(IN) :: len_2d
  REAL(KIND=wp), INTENT(IN) :: dt_global, vrel(n_s), alpha(n_s), &
    Gamma_history(n_s), dyn_comp_history(n_s,3), &
    alpha_2d(len_2d,n_s), cl_2d(len_2d,n_s), chord(n_s)      
  !-----[out]  
  REAL(KIND=wp), INTENT(OUT) :: Gamma_steady(n_s), Gamma(n_s),&
    dGamma(n_v), dyn_comp(n_s,3)
  !-----[tmp]        
  INTEGER::l
  REAL(KIND=wp), DIMENSION(n_s)::cl_inter,dt,T0
  !=======================================================================


  !----Circulation at sections----
  ! Interpolate the lift coefficients, output is cl_inter
  DO l = 1,n_s
  CALL nw_interp1d(alpha_2d(:,l), cl_2d(:,l), len_2d, alpha(l), &
    cl_inter(l),1)
  END DO

  ! Calculate the steady circulation at each blade element
  Gamma_steady = 0.5_wp*vrel*chord*cl_inter

  !----Time filter bound circulation----
  ! Madsen and Gaunaa (2004)
  T0 = chord / (2.0_wp * vrel)
  dt = dt_global / T0
  ! Filter functions by Mac Gaunaa
  dyn_comp(:,1) = dyn_comp_history(:,1) &
    * exp(-0.3064_wp * dt) + 0.27735_wp * (Gamma_steady &
    + Gamma_history) * (1.0_wp - exp(-0.3064_wp * dt) )

  dyn_comp(:,2) = dyn_comp_history(:,2)  &
    * exp(-0.0439_wp * dt) + 0.0914_wp * (Gamma_steady &
    + Gamma_history) * (1.0_wp - exp(-0.0439_wp * dt))
  dyn_comp(:,3) = dyn_comp_history(:,3) &
    * exp( -3.227_wp * dt) + 0.13125_wp * (Gamma_steady &
    + Gamma_history) * (1.0_wp - exp( -3.227_wp * dt))

  ! Dynamic filtered circulation
  Gamma = dyn_comp(:,1) + dyn_comp(:,2) + dyn_comp(:,3)
  ! Change in Gamma between sections.At blade ends take all Gamma
  ! from neighboring point. 
  dGamma(2:n_v-1) = Gamma(2:n_s) - Gamma(1:n_s-1)
  dGamma(1) = Gamma(1)
  dGamma(n_v)= -Gamma(n_s);  

END SUBROUTINE nw_circulation
!=======================================================================
!=======================================================================
SUBROUTINE nw_calc_relax(vrel, alpha, vn, vt, chord, dt_global, &
    nw_dx, nw_dy, dbeta, nw_phi, &
    nw_relax_safety, nw_relax)
  !-----------------------------------------------------------------------      
  !     Determine the relaxation factor needed in the near wake model. 
  !     The defintions of Pirrung et al. (2017a) are used and the 
  !     corresponding equations referenced.
  !-----------------------------------------------------------------------       
  USE kinds, ONLY: wp ! Global precision
  use nw_smearing_correction_module, ONLY: pi, n_s,  n_v ! Constants
  !-----------------------------------------------------------------------        
  IMPLICIT NONE



  !-----[in]  
  REAL(KIND=wp), INTENT(IN)::vrel(n_s), alpha(n_s), &
    vt(n_s), &
    dt_global, nw_dx(n_v,n_s), &
    nw_dy(n_v,n_s), &
    dbeta(n_v), nw_phi(n_v,n_s), &
    nw_relax_safety, chord(n_s), &
    vn(n_s)
  !-----[out]  
  REAL(KIND=wp), INTENT(OUT) :: nw_relax
  !-----[tmp]     
  REAL(KIND=wp) ::  A1, A2, gradw, tau(n_s), d(n_s), &
    nw_X(n_v,n_s), nw_Y(n_v,n_s)
  INTEGER :: l, k
  !=======================================================================


  !----Initialize----   
  nw_relax = 0.0_wp

  !----Dynamic to steady circulation ratio----
  ! Eq. (15)

  tau = (2.0_wp * vrel) / chord

  d = 1.0_wp - 0.5547_wp / (tau * dt_global * 0.3064_wp) *         &
    (1.0_wp - exp(-tau * dt_global * 0.3064_wp)) - 0.1828_wp /     &
    (tau * dt_global * 0.0439_wp) * (1.0_wp -exp(-tau * dt_global  &
    * 0.0439_wp)) - 0.2625_wp / (tau * dt_global * 3.2277_wp )     &
    * (1.0_wp - exp(-tau * dt_global * 3.2277_wp ))


  !----Slow and fast decaying components---- 
  ! in-plane only (no phi_star)
  DO l = 1, n_s
    DO k = 1, n_v

      !write(*,*)'In calc relax Term 4 Schleife', nw_phi(k,l)   
      nw_X(k,l) = nw_dx(k,l) * nw_phi(k,l)                         &
        * (1.0_wp - exp( -dbeta(k) / nw_phi(k,l)))
      nw_Y(k,l) = nw_dy(k,l) * nw_phi(k,l) &
        * (1.0_wp - exp(-4.0_wp * dbeta(k) / nw_phi(k,l) ))
    END DO
  END DO



  !----Time derivative of induced velocities----
  ! Eq. (27)
  DO l = 1, n_v -1
    A1 = nw_X(l, l) + nw_Y(l, l)
    A2 = nw_X(l+1, l) + nw_Y(l+1, l)
    gradw = d(l) * pi * chord(l) * (A1-A2) *  &
      ( (alpha(l) * vn(l)) / vrel(l) +  &
      vrel(l) / (ABS(vt(l)) * (vn(l) / ABS(vt(l)))**2 + 1.0_wp) ) 
    ! Update the relaxation factor
    nw_relax = MAX(nw_relax, -(1.0_wp + gradw) / (1.0_wp - gradw))
  END DO


  !----Mix in the safety factor---- 
  nw_relax = nw_relax + nw_relax_safety * (1.0_wp - nw_relax)
  nw_relax = MIN(nw_relax, 0.99_wp)



END SUBROUTINE nw_calc_relax
!=======================================================================
!=======================================================================
SUBROUTINE init_nw_smearing_correction(n_rotors_in, n_blades_in, &
    n_sections_in)
  !-----------------------------------------------------------------------
  !     Initialize the module-level scalars and the Pirrung 2017b
  !     constant tables (nw_P, nw_N) used by the near-wake correction.
  !-----------------------------------------------------------------------
  USE kinds, ONLY: wp
  USE nw_smearing_correction_module ! Assign all
  !-----------------------------------------------------------------------
  IMPLICIT NONE

  !-----[in]
  INTEGER, INTENT(IN) ::  n_rotors_in, n_blades_in, n_sections_in
  !=======================================================================

  IF ( .NOT. smearing_correction_active) RETURN ! Switch

  ! TEilers
  IF ( .NOT. smearing_correction_needs_init) RETURN ! Switch

  !----Copy for module----
  n_rotors = n_rotors_in
  n_blades = n_blades_in

  !----Number of control and vortex trailing points----
  n_s = n_sections_in
  n_v = n_s + 1

  !----Near-wake constants Pirrung et al. (2017b)----
  ! Positve values of h_r given in Eq.(14)
  nw_P(1, :) = [ -1.64636754184988_wp,  8.14821474772595_wp, -12.1784861715161_wp, 5.02652654794305_wp, 21.7713111885100_wp ]
  nw_P(2, :) = [ -0.499006672527007_wp, 6.08465268675599_wp, -15.1712005446736_wp, 14.8254135311781_wp, -2.42318681495122_wp ]
  nw_P(3, :) = [ 3.90835606394205_wp,  -18.7662330478906_wp, 39.1243282362892_wp, -29.4870086106427_wp, -60.2947309276912_wp ]
  nw_P(4, :) = [ -1.60622843516794_wp,  7.42952711674470_wp, -15.8594777074799_wp, 11.6870163960092_wp, -195.060865844227_wp ]

  ! Negative values given in Eq.(13)
  nw_N(1, :) = [ 1.01933206359188_wp,  -0.135668151345576_wp, 0.395524705910666_wp, 0.0801822997046969_wp, 44.8347503471217_wp ]
  nw_N(2, :) = [ 12.9874512870379_wp,   49.9999999999916_wp, 0.00235345218088083_wp, 11.3116074964299_wp, 3935.34323353354_wp ]
  nw_N(3, :) = [ -0.690159136848837_wp, 101.238775997879_wp, -0.00154467589177718_wp, 3.99520376143011_wp, 0.394541495293350_wp ]
  nw_N(4, :) = [ -0.269253409815015_wp, 49.9999999999843_wp, -0.00247975579295973_wp, 0.403642518874381_wp, 1.16610462658971_wp ]

  ! S. Ouchene: the actual nw_* history arrays and cp_loop_ind / smear_fac
  ! / recompute_smearing_constants are allocated in f8c_parin (restart-aware
  ! zeroing there). No allocation lives in this routine.

  ! TEilers
  smearing_correction_needs_init = .FALSE.

END SUBROUTINE init_nw_smearing_correction
!=======================================================================
!=======================================================================
SUBROUTINE nw_interp1d(x, y, n, xin, yin, nin)
  !-----------------------------------------------------------------------      
  !     One-dimensional linear interpolation
  !-----------------------------------------------------------------------      
  USE kinds, ONLY: wp
  !-----------------------------------------------------------------------      
  IMPLICIT NONE

  !-----[in]  
  INTEGER, INTENT(IN) :: n, nin
  REAL(KIND=wp), INTENT(IN) :: x(n), y(n), xin(nin)

  !-----[out]  
  REAL(KIND=wp), INTENT(OUT) :: yin(nin)
  
  !-----[in]        
  INTEGER :: i, j
  REAL(KIND=wp) :: ratio

  !=======================================================================
  DO j = 1, nin
    !------  extrapolate
    IF ( xin(j) > x(n) ) THEN
      ratio  = ( xin(j) - x(n) ) / ( x(n) - x(n-1) )
      yin(j) = y(n) + ratio * ( y(n) - y(n-1) )
    END IF
    !------  extrapolate
    IF ( xin(j) < x(1) ) THEN
      ratio  = ( xin(j) - x(1) ) / ( x(2) - x(1) )
      yin(j) = y(1) + ratio * ( y(2) - y(1) )
    END IF
    !------  interpolate
    DO i = 2, n
      IF ( x(i) >= xin(j) .AND. x(i-1) <= xin(j) ) THEN
        ratio  = (xin(j) - x(i-1) ) / ( x(i) - x(i-1) )
        yin(j) = y(i-1) + ratio * ( y(i) - y(i-1) )
        EXIT
      END IF
    END DO ! i loop
  END DO ! j loop

END SUBROUTINE nw_interp1d
!=======================================================================
SUBROUTINE smooth_trailed_vorticity( Gamma, dGamma, r_cp, smoothing_length )
  !---------------------------------------------------------------------
  ! Optional spanwise Gaussian smoothing of the bound circulation Gamma,
  ! applied LOCALLY at the cyl/airfoil interface only. When
  ! smoothing_length <= 0 (the default behaviour), this routine returns
  ! immediately and Gamma / dGamma are left exactly as nw_circulation
  ! produced them.
  !
  ! Localisation: only elements that sit adjacent to a transition between
  ! a zero-lift (cylinder) and a lifting (airfoil) section get their
  ! Gamma touched by the kernel. The classification is read from
  ! is_lifting_section_arr, populated once by f8c_build_nwc_arrays
  ! from the per-element polar tables. Elements far from any cyl/airfoil
  ! interface keep their raw Gamma(k); dGamma at those boundaries is
  ! identical to the unsmoothed result.
  !
  ! For a blade without any cylinder section (every element lifting) or
  ! every element a cylinder (extreme test case), no transition exists
  ! and the routine produces zero changes regardless of smoothing_length.
  !
  ! When the kernel does fire, a 3-point Gaussian of width L is applied
  ! to Gamma(k); dGamma is then rebuilt from the smoothed Gamma using
  ! the same convention as nw_circulation. The kernel preserves total
  ! bound circulation along the blade, so integrated thrust and torque
  ! are essentially unchanged - only the spanwise concentration of
  ! trailed vorticity at the cyl/airfoil interface is regularised.
  ! Mirror boundary conditions at root and tip prevent endpoint pull.
  !
  ! Motivation: at the root cylinder -> airfoil transition Cl jumps from
  ! ~0 to ~1.3 over a single trailed-vortex boundary. The Pirrung /
  ! Meyer-Forsting wake model then amplifies that delta-function-like
  ! dGamma into a localised alpha spike that does not reflect real
  ! viscous wake-sheet width. A small L (of order one chord) absorbs the
  ! step into a length scale that has physical meaning. Localising the
  ! kernel keeps every other element along the blade untouched.
  !
  ! Added by S. Ouchene, May 2026.
  !---------------------------------------------------------------------
  USE kinds, ONLY: wp
  USE nw_smearing_correction_module, ONLY: n_s, n_v, is_lifting_section_arr
  IMPLICIT NONE

  REAL(KIND=wp), INTENT(INOUT) :: Gamma(n_s)
  REAL(KIND=wp), INTENT(INOUT) :: dGamma(n_v)  ! INOUT so the early-return paths leave the caller's dGamma_vtp as nw_circulation set it
  REAL(KIND=wp), INTENT(IN)    :: r_cp(n_s)
  REAL(KIND=wp), INTENT(IN)    :: smoothing_length

  REAL(KIND=wp) :: Gamma_in(n_s)
  REAL(KIND=wp) :: inv_2L2, w_left, w_right, gamma_left, gamma_right
  LOGICAL :: at_transition(n_s)
  INTEGER :: k

  IF ( smoothing_length <= 0.0_wp ) RETURN

  ! Mark every element whose left- or right-side neighbour belongs to a
  ! different category (lifting vs zero-lift). For the standard NREL 5 MW
  ! blade with Cyl1/Cyl2 root and DU40+ outboard, this flags exactly the
  ! two elements straddling the Cyl2/DU40 boundary.
  at_transition(:) = .FALSE.
  DO k = 2, n_s
    IF ( is_lifting_section_arr(k) .NEQV. is_lifting_section_arr(k-1) ) THEN
      at_transition(k-1) = .TRUE.
      at_transition(k)   = .TRUE.
    END IF
  END DO

  ! If no transition exists, leave Gamma untouched.
  IF ( .NOT. ANY( at_transition ) ) RETURN

  Gamma_in = Gamma
  inv_2L2  = 1.0_wp / ( 2.0_wp * smoothing_length * smoothing_length )

  ! 3-point Gaussian kernel with mirror boundaries, applied only at
  ! flagged elements. Every other element keeps Gamma_in(k) unchanged.
  DO k = 1, n_s
    IF ( .NOT. at_transition(k) ) CYCLE
    IF ( k > 1 ) THEN
      w_left     = EXP( -( r_cp(k)   - r_cp(k-1) )**2 * inv_2L2 )
      gamma_left = Gamma_in(k-1)
    ELSE
      w_left     = EXP( -( r_cp(k+1) - r_cp(k)   )**2 * inv_2L2 )
      gamma_left = Gamma_in(k)
    END IF
    IF ( k < n_s ) THEN
      w_right     = EXP( -( r_cp(k+1) - r_cp(k)   )**2 * inv_2L2 )
      gamma_right = Gamma_in(k+1)
    ELSE
      w_right     = EXP( -( r_cp(k)   - r_cp(k-1) )**2 * inv_2L2 )
      gamma_right = Gamma_in(k)
    END IF
    Gamma(k) = ( w_left * gamma_left + Gamma_in(k) + w_right * gamma_right ) /     &
               ( w_left + 1.0_wp + w_right )
  END DO

  ! Rebuild trailed vorticity at boundaries (same convention as nw_circulation).
  dGamma(2:n_v-1) = Gamma(2:n_s) - Gamma(1:n_s-1)
  dGamma(1)       =  Gamma(1)
  dGamma(n_v)     = -Gamma(n_s)
END SUBROUTINE smooth_trailed_vorticity
!=======================================================================
