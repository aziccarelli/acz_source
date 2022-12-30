c
c           user_routines_umat.f   Distribution version
c
c           Updated: 8/20/2016 rhd
c
c
c           The first routine here, umat_set_features, is called by WARP3D
c           at various times to obtain information about the umat.
c
c           This file is compiled and included in the normal executable
c           build of WARP3D.
c
c           WARP3D calls umat routines in:
c
c              gplns1.f    -- linear stiffness computation
c              rstgp1.f    -- stress update and new [D] computation
c
c           These above two code fies have lots of comments about setting
c           up data arrays, values for the umat from WARP3D data
c           structures.
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine umat_set_features                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/14 rhd                *
c     *                                                              *
c     *  called by warp3d for the umat to obtain statev size and     *
c     *  other characteristic information about the umat             *
c     *                                                              *
c     ****************************************************************
c
      subroutine umat_set_features( info_vector )
      implicit none
      integer info_vector(*)
c
c        set info_data
c
c         1        number of history values per integration
c                  point. Abaqus calls these "statev". Values
c                  double or single precsion based on hardware.
c
c         2        number of values in the symmetric part of the
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
c
c
c        example umat included here is for mises plasticity with
c        bilinear kinematic hardening
c
      info_vector(1) = 27   ! includes back stresses
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 4
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1RPL,DDSDDT,DRPLDE,DRPLDT,
     2STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC,
     5KITER,KOUT,KTHREAD,KNUMTHREADS,NONLOCAL_SHARED,NSHARED)
      ! Note that IMPLICIT definition is active
C      INCLUDE 'ABA_PARAM.INC'
C
      implicit none 
      CHARACTER*80 CMNAME
C
C      DIMENSION STRESS(NTENS),STATEV(NSTATV),
C     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
C     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
C     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
C     4JSTEP(4)
C
      integer :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer, 
     1           kspt, kstep, kinc, kiter, kout, kthread, 
     2           knumthreads, nshared, i, j, n_backstresses, it_num,
     3           converged, ind_alpha
c
      double precision :: stress(ntens), statev(nstatv),
     1 ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens), stran(ntens),
     2 dstran(ntens), predef(1), dpred(1), props(nprops), coords(3),
     3 drot(3,3), dfgrd0(3,3), dfgrd1(3,3), nonlocal_shared(nshared),
     4 jstep(4), olds(ntens)
c
      double precision :: time, dtime, temp, dtemp, sse, spd, scd, rpl,
     1                    drpldt, celent, pnewdt, elastic_modulus, 
     2                    Q_inf, b, D_inf, a, shear_modulus,
     3                    bulk_modulus, poission_ratio, mu2, lame_fist,
     4                    yield_stress, ep_eq, ep_eq_init, a_temp,
     5           hard_iso_Q, hard_iso_D, hard_iso_total, a_dot_n,
     6          plastic_mult, p_mult_numer, p_mult_denom, yield_function,
     7          isotropic_modulus, kin_modulus,
     8          stress_relative_norm, strain_trace, alpha_trace, e_k,
     9          ID2_out_ID2, n_out_n, stress_hydro, sigma_vm,ddczm_on

      double precision :: Lam, n33_check, alpha_out_n, beta, theta_1, 
     2         theta_2, theta_3,srn2, yield_function, isotropic_modulus,
     3         triax_n1
C
c      ! Subroutine control
c      INTEGER :: i, j, n_backstresses, it_num, converged, ind_alpha
c      ! Material properties
c      REAL(8) :: elastic_modulus, Q_inf, b, D_inf, a,
c     1shear_modulus, bulk_modulus, poission_ratio, mu2, lame_first
c      ! Used for intermediate calculations
c      REAL(8) :: yield_stress, ep_eq, ep_eq_init, a_temp,
c     1hard_iso_Q, hard_iso_D, hard_iso_total, a_dot_n,
c     2plastic_mult, p_mult_numer, p_mult_denom, yield_function,
c     3isotropic_modulus, kin_modulus,
c     4stress_relative_norm, strain_trace, alpha_trace, e_k,
c     5ID2_out_ID2, n_out_n, stress_hydro, sigma_vm,
c     6Lam, n33_check, alpha_out_n, beta, theta_1, theta_2, theta_3,srn2
      ! Backstress arrays
      REAL(8), DIMENSION(:, :), ALLOCATABLE :: alpha_k
      REAL(8), DIMENSION(:), ALLOCATABLE :: C_k, gamma_k
      ! Tensors
      REAL(8), DIMENSION(6, 6) :: ID4, c_mat
      REAL(8), DIMENSION(6) :: strain_tens, strain_plastic,
     1yield_normal, alpha, strain_trial, stress_relative,
     2stress_dev, ID2, stress_tens, check, dstran_tens, alpha_diff,
     3alpha_upd, dpe
      ! Parameters
      INTEGER :: N_BASIC_PROPS, TERM_PER_BACK, MAX_ITERATIONS,
     1I_ALPHA, NUM_PROPS
      REAL(8) :: TOL, ONE, TWO, THREE, ZERO, SQRT23
      PARAMETER(TOL=1.0D-10,
     1N_BASIC_PROPS=7, TERM_PER_BACK=2, MAX_ITERATIONS=1000,
     2ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, ZERO=0.D0,
     3SQRT23=SQRT(2.0D0/3.0D0), I_ALPHA=7,NUM_PROPS=1.1D1)
C ----------------------------------------------------------------------C
C
      ! Subroutine start
C
C ----------------------------------------------------------------------C
      ! Get the number of backstresses
      n_backstresses = (NUM_PROPS - N_BASIC_PROPS) / TERM_PER_BACK
C      print*,'num_props:', NUM_PROPS
C      print*,'num_backstresses:',n_backstresses
C      DO i = 1, nprops
C        print*,i
C        print*,props(i)
C      END DO
      IF (n_backstresses .EQ. 0) THEN
        PRINT *, "No backstresses defined, exiting!"
        CALL XIT  ! Exit from analysis command in Abaqus
      END IF
C
      ! Allocate the backstress related arrays
      ALLOCATE(C_k(n_backstresses))
      ALLOCATE(gamma_k(n_backstresses))
      ALLOCATE(alpha_k(n_backstresses, ntens))
C
      ! Initialize
      ddsdde(:, :) = ZERO
      ID4(:, :) = ZERO  ! 4th order symmetric identity tensor
      DO i = 1, ndi
        ID4(i, i) = ONE  
      END DO
      DO i = ndi+1, ntens
        ID4(i, i) = ONE / TWO
      END DO
      ! 2nd order symmetric identity tensor
      ID2(:) = (/ ONE, ONE, ONE, ZERO, ZERO, ZERO /)
C
      ! Read in state variables
      ! 1               = Equivalent plastic strain
      ! 2 - 7           = Plastic strain
      ! 8 - 8 + 6 * N   = Backstress components
      ep_eq = statev(1)
      ep_eq_init = statev(1)
      CALL ROTSIG(statev(2), drot, strain_plastic, 2, ndi, nshr)
      alpha(:) = ZERO
      DO i = 1, n_backstresses
        ind_alpha = I_ALPHA + 1 + (i - 1) * ntens
        CALL ROTSIG(statev(ind_alpha), drot, alpha_k(i, :), 1, ndi,nshr)
        alpha = alpha + alpha_k(i, :)
      END DO
C
      ! Read in the material properties
      elastic_modulus = props(1)
      poission_ratio = props(2)
      yield_stress = props(3)
      q_inf = props(4)
      b = props(5)
      d_inf = props(6)
      a = props(7)
      DO i = 1, n_backstresses  ! First backstress starts at index = 8
        c_k(i) = props((N_BASIC_PROPS - 1) + 2 * i)
        gamma_k(i) = props(N_BASIC_PROPS + 2 * i)
      END DO
      ddczm_on = props(12)
C
C      print*, 'e:',elastic_modulus
C      print*, 'nu:',poission_ratio
C      print*, 'fy:',yield_stress
C      print*, 'qinf:',q_inf
C      print*, 'b:',b
C      print*, 'd_inf:',d_inf
C      print*, 'a:',a
C      print*, 'c1:',props(8)
C      print*, 'gamma1', props(9)
C
      ! Calculate elastic parameters
      shear_modulus = elastic_modulus / (TWO * (ONE + poission_ratio))
      bulk_modulus = elastic_modulus /
     1(THREE * (ONE - TWO * poission_ratio))
      mu2 = TWO * shear_modulus
C
      ! Set-up strain tensor
      ! Tensors are stored: 11, 22, 33, 12, 13, 23
      strain_tens = stran + dstran
C ----------------------------------------------------------------------C
C
      ! Elastic trial step
C
C ----------------------------------------------------------------------C
      ! Tensor of elastic moduli
      DO j = 1, ntens
        DO i = 1, ntens
          ID2_out_ID2 = ID2(i) * ID2(j)
          c_mat(i, j) = ID2_out_ID2 * bulk_modulus +
     1    mu2 * (ID4(i, j) - ONE/THREE * ID2_out_ID2)
        END DO
      END DO
      ! Stress tensor
      olds = stress
      stress_tens = stress + MATMUL(c_mat, dstran)
      !stress_tens = MATMUL(c_mat, (strain_tens - strain_plastic))
C
      stress_hydro = SUM(stress_tens(1:3)) / THREE
      strain_trace = SUM(strain_tens(1:3))
      DO i = 1, ndi
        stress_dev(i) = stress_tens(i) - stress_hydro
        stress_relative(i) = stress_dev(i) - alpha(i)
      END DO
      DO i = ndi+1, ntens
        stress_dev(i) = stress_tens(i)
        stress_relative(i) = stress_dev(i) - alpha(i)
      END DO
      stress_relative_norm = 
     1sqrt(dotprod6(stress_relative, stress_relative))
C
      ! Yield condition
      hard_iso_Q = q_inf * (ONE - EXP(-b * ep_eq))
      hard_iso_D = d_inf * (ONE - EXP(-a * ep_eq))
      hard_iso_total = yield_stress + hard_iso_Q - hard_iso_D
      yield_function = stress_relative_norm - SQRT23 * hard_iso_total
      IF (yield_function .GT. TOL) THEN
        converged = 0
      ELSE
        converged = 1
      END IF
C
      ! Calculate the normal to the yield surface
      yield_normal = stress_relative / (TOL + stress_relative_norm)
C ----------------------------------------------------------------------C
C
      ! Radial return mapping if plastic loading
C
C ----------------------------------------------------------------------C
      ! Calculate the consitency parameter (plastic multiplier)
      plastic_mult = ZERO
      it_num = 0
      DO WHILE ((converged .EQ. 0) .AND. (it_num .LT. MAX_ITERATIONS))
        it_num = it_num + 1
C
        ! Calculate the isotropic hardening parameters
        hard_iso_Q = q_inf * (ONE - EXP(-b * ep_eq))
        hard_iso_D = d_inf * (ONE - EXP(-a * ep_eq))
        hard_iso_total = yield_stress + hard_iso_Q - hard_iso_D
        isotropic_modulus = b * (q_inf - hard_iso_Q) -
     1  a * (d_inf - hard_iso_D)
        ! Calculate the kinematic hardening parameters
        kin_modulus = ZERO
        DO i = 1, n_backstresses
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          kin_modulus = kin_modulus + C_k(i) * e_k 
     1    - SQRT(THREE/TWO)*gamma_k(i)*e_k
     2    * dotprod6(yield_normal, alpha_k(i, :))
        END DO
        a_dot_n = ZERO
        alpha_upd(:) = ZERO
        DO i = 1, n_backstresses
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          alpha_upd = alpha_upd + e_k * alpha_k(i, :)
     1    + SQRT23 * C_k(i) / gamma_k(i) * (ONE - e_k) * yield_normal
        END DO
        a_dot_n = dotprod6(alpha_upd - alpha, yield_normal)  ! n : \Delta \alpha
        
        p_mult_numer = stress_relative_norm -
     1  (a_dot_n + SQRT23 * hard_iso_total + mu2 * plastic_mult)
C
        p_mult_denom = -mu2 *
     1  (ONE + (kin_modulus + isotropic_modulus) /
     2  (THREE * shear_modulus))
C
        ! Update variables
        plastic_mult = plastic_mult - p_mult_numer / p_mult_denom
        ep_eq = ep_eq_init + SQRT23 * plastic_mult
C
        IF (ABS(p_mult_numer) .LT. TOL) THEN
          converged = 1
        END IF
      END DO
C ----------------------------------------------------------------------C
C
      ! Update variables
C
C ----------------------------------------------------------------------C
      IF (it_num .EQ. 0) THEN  ! Elastic loading
        stress = stress_tens
      ELSE  ! Plastic loading
        !strain_plastic = strain_plastic + plastic_mult * yield_normal
        dpe = plastic_mult * yield_normal
        dpe(4:6) = dpe(4:6) + plastic_mult * yield_normal(4:6)
C        strain_plastic(4:6) = strain_plastic(4:6) 
C     1  + plastic_mult * yield_normal(4:6)
        strain_plastic = strain_plastic + dpe
        stress = stress_tens - MATMUL(c_mat, dpe)
        !stress = MATMUL(c_mat, (strain_tens - strain_plastic))
        statev(26) = ep_eq - ep_eq_init
C
        alpha_diff = alpha
        alpha(:) = ZERO
        DO i = 1, n_backstresses  ! Update backstress components
          e_k = EXP(-gamma_k(i) * (ep_eq - ep_eq_init))
          alpha_k(i, :) = e_k * alpha_k(i, :) +
     1    SQRT23 * yield_normal * C_k(i) / gamma_k(i) * (ONE - e_k)
          alpha = alpha + alpha_k(i, :)
        END DO
        alpha_diff = alpha - alpha_diff
      END IF
C
C     Tangent modulus
      IF (it_num .EQ. 0) THEN  ! Elastic loading
        DO j = 1, ntens
          DO i = 1, ntens
            ddsdde(i, j) = c_mat(i, j)
          END DO
        END DO
        DO j = ndi+1, ntens
          ddsdde(j, j) = shear_modulus
        END DO
      ELSE  ! Plastic loading
        beta = ONE +
     1  (kin_modulus + isotropic_modulus) / (THREE * shear_modulus)
        theta_1 = ONE - mu2 * plastic_mult / stress_relative_norm
        theta_3 = ONE / (beta * stress_relative_norm)
        theta_2 = ONE / beta 
     1  + dotprod6(yield_normal, alpha_diff) * theta_3 
     2  - (ONE - theta_1)
        DO j = 1, ntens
          DO i = 1, ntens
            ID2_out_ID2 = ID2(i) * ID2(j)
            n_out_n = yield_normal(i) * yield_normal(j)
            alpha_out_n = alpha_diff(i) * yield_normal(j)
            ddsdde(i, j) = bulk_modulus * ID2_out_ID2
     1      +mu2 * theta_1*(ID4(i, j) - ONE/THREE*ID2_out_ID2)
     2      -mu2 * theta_2 * n_out_n +
     3      +mu2 * theta_3 * alpha_out_n
          END DO
        END DO
        ddsdde = ONE/TWO * (TRANSPOSE(ddsdde) + ddsdde)
      END IF
C ----------------------------------------------------------------------C
C
      ! Update the state variables
C
C ----------------------------------------------------------------------C
      statev(1) = ep_eq
      DO i = 1, ntens
        statev(i + 1) = strain_plastic(i)
      END DO
      DO i = 1, n_backstresses
        DO j = 1, ntens
          statev(I_ALPHA + j + (i-1) * ntens) = alpha_k(i, j)
        END DO
      END DO

      triax_n1 = statev(25)
      if (ddczm_on .lt. one) return
      call umat_compute_ddczm(kinc, kiter, ntens, nstatv, npt,
     &               stress, olds, statev, ep_eq_init, noel, triax_n1)
c
c    Update nonlocal variables
c
      nonlocal_shared(1) = one
      nonlocal_shared(2) = one
      nonlocal_shared(3) = statev(23)
      nonlocal_shared(4) = statev(25)

C ----------------------------------------------------------------------C
C
      ! Reduce time increment if did not converge
C
C ----------------------------------------------------------------------C
      IF (it_num .EQ. MAX_ITERATIONS) THEN
        PRINT *, "WARNING: Return mapping in integration point ", npt,
     1  " of element ", noel, " did not converge."
        PRINT *, "Reducing time increment to 1/4 of current value."
        PNEWDT = 0.25
      END IF
      RETURN
        
      CONTAINS
C ----------------------------------------------------------------------C
C
      ! Define dot product for vectors
C
C ----------------------------------------------------------------------C
      pure function dotprod6(A, B) result(C)
      !! Returns the dot product of two symmetric length 6 vectors
      !! that are reduced from 9 components, with the last 3 symmetric
      REAL(8), intent(in) :: A(6), B(6)
      REAL(8)             :: C
      INTEGER             :: i      
      ! Calculate the dot product
      C = 0.0D0
      DO i = 1, 3
        C = C + A(i) * B(i)
      END DO
      DO i = 4, 6
        C = C + TWO * (A(i) * B(i))
      END DO
      end function
      END


c    ****************************************************************
c    *                                                              *
c    *               subroutine umat_compute_ddczm                  *
c    *                                                              *
c    *         written by : Andy Ziccarelli                         *
c    *      last modified : 05/15/2019 AJZ                          *
c    *                                                              *
c    *     subroutine to compute DDCZM damage parameter.            *
c    *                                                              *
c    ****************************************************************
c
c
      subroutine umat_compute_ddczm(kinc, kiter, ntens, nstatv, npt,
     &                 stress, olds, statev, peeq_n, noel, triax_n1)

      use mod_damage_ddczm
      implicit none
c
      integer, intent(in) :: kinc, kiter, ntens, nstatv, noel, npt
c
      double precision, intent(in) :: 
     &       stress(ntens), olds(ntens)
c
c
      double precision, intent(inout) :: statev(nstatv)
c
      double precision, intent(out) :: triax_n1
c
c            locals
c
      logical :: local_debug
c
      double precision :: 
     &       mises, peeq_n, peeq_n1, lodeang,
     &       peeq_comp_n, peeq_comp_n1,
     &       dmg_intgrnd_n, dmg_intgrnd_n1, 
     &       dmg_intgrl_n, dmg_intgrl_n1, damage,
     &       c_val
c
      local_debug = .false.
c
      c_val = swdm_c
c
c
c            RETRIEVE HISTORY DATA
c
c
      dmg_intgrnd_n = statev(21)
      dmg_intgrl_n = statev(22)
      peeq_comp_n = statev(24)
c
c     retrieve current peeq
      peeq_n1 = statev(1)
c
c
      call compute_swdm_ddczm(stress, mises, triax_n1,
     &               peeq_n, peeq_n1, lodeang, dmg_intgrnd_n,
     &               dmg_intgrnd_n1, dmg_intgrl_n, dmg_intgrl_n1,
     &               peeq_comp_n, peeq_comp_n1, damage,
     &               c_val, swdm_kappa, swdm_lambda, swdm_beta )
c
c
c            UPDATE HISTORY DATA
c
c
      statev(27) = damage-statev(23)
      statev(21) = dmg_intgrnd_n1
      statev(22) = dmg_intgrl_n1
      statev(23) = damage
      statev(24) = peeq_comp_n1
      statev(25) = triax_n1 
c
      if (local_debug) then
           print*, 'mises: ', mises
           print*, 'triax: ', triax_n1
           print*, 'peeq_n1: ', peeq_n1
           print*, 'damage: ', damage
           print*, 'dmg_intgrl_n: ', dmg_intgrl_n
           print*, 'cval: ', c_val
      end if  
c
      return
      end
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine umat_states_labels                    *
c     *                                                              *
c     *   call by WARP3D to get the number of states variables for   *
c     *   output, an 8 character id for each variable and an         *
c     *   descriptor string for each variable                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine umat_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      implicit none
c
c                       parameters
c
      integer :: size_state, num_states, out, max_comment_lines, 
     &           num_comment_lines  
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
      integer :: i
c
      num_states = 4 
      num_comment_lines = 0
c
c      state_labels(1) = "epspls"
c      state_labels(2) = "state"
c      state_labels(3) = "epls_xx"
c      state_labels(4) = "epls_yy"
c      state_labels(5) = "epls_zz"
c      state_labels(6) = "epls_xy"
c      state_labels(7) = "epls_xz"
c      state_labels(8) = "epls_yz"
c      state_labels(9) = "alpha_xx"
c      state_labels(10) = "alpha_yy"
c      state_labels(11) = "alpha_zz"
c      state_labels(12) = "alpha_xy"
c      state_labels(13) = "alpha_xz"   ! abaqus ordering
c      state_labels(14) = "alpha_yz"   ! abaqus ordering
c
c      state_descriptors(1) = "Equivalent plastic strain"
c      state_descriptors(2) = "=1, 2, 3 see states_header"
c      state_descriptors(2:7) = "Plastic strain component"
c      state_descriptors(8:14) = "Backstress component"
c
      state_labels(1) = "peeq"
      state_labels(2) = "SWDM"
      state_labels(3) = "dpeeq"
      state_labels(4) = "dSWDM"
c
c
      state_descriptors(1) = "Equiv Plastic Strain"
      state_descriptors(2) = "Ductile Damage"
      state_descriptors(3) = "Change in plastic strain"
      state_descriptors(4) = "Change in Damage"
c
c
c      comment_lines(1) = "Notes on 'state' quantity"
c      comment_lines(2) = "  = 1 never yielded"
c      comment_lines(3) = "  = 2 active yielding"
c      comment_lines(4) = "  = 3 previous plasticity. " 
c     &                       // " now linear-elastic"
c      comment_lines(5) = "  value is average of int points"
c  
      return
      end
      
c
c *******************************************************************
c *                                                                 *
c *             subroutine umat_output_states                       *
c *                                                                 *
c *     routine called just before output of states. return         *
c *     state values for output at this integration point           *
c *                                                                 *
c *******************************************************************
c
c
      subroutine umat_output_states( statev, output_statev, nvals )
      implicit none
c
c                       parameters
c
      integer :: nvals
      double precision ::
     & statev(*), output_statev(*)
c
c                       locals
c
      integer :: k, k1
c
      nvals = 4
c      
c    First output: PEEQ
c    Second output: SWDM
c
      output_statev(1) = statev(1)
      output_statev(2) = statev(23)
      output_statev(3) = statev(26)
      output_statev(4) = statev(27)
c
      return
      end
c
c *******************************************************************
c *                                                                 *
c *  optional UMAT routine to set additional stress output values   *
c *                                                                 *
c *   set up to 3 material model dependent output values            *
c *                                                                 *
c *******************************************************************
c
c
      subroutine umat_output( stress, mat_vals, statev, nstatv, time,
     &                        cmname, props, nprops, noel,
     &                        npt, kinc, kout )
      implicit double precision (a-h,o-z)
c
c                   parameter declarations
c                   ----------------------
c
      character(len=8) :: cmname
      double precision
     & stress(6), mat_vals(3), statev(nstatv), props(nprops)
c
c               description of parameters
c               -------------------------
c
c     stress            : 6 x 1 stresses for point. Cauchy stresses for
c                         GEONL solution. Abaqus ordering (see below)
c (*) mat_vals          : 3x1. these are the c1, c2, c3 values that may be
c                         set by this routine for including in the stress
c                         output below. Values are floating point. WARP3D
c                         does nothing but print values.
c     statev            : current state vector for point
c     nstatv            : umat decalerd number of state variables/point.
c                         see routine umat_set_features
c     time              : current simultation time (sum of all load step
c                         time increments)
c     cmname            : always "UMAT-WRP"
c     props             : values of umat properties specified by user
c                         in model definition (um_1, um_2, um_3, ...)
c     nprops            : umat spedified number of properties
c                         (see umat_set_features)
c     noel              : current element number in model
c     npt               : current material point number
c     kinc              : current analysis (WARP3D) load step number
c     kout              : Fortran device number for debugging messages
c
c    =>   only mat_vals should be modified by this routine
c
c   stress ordering
c     (1) sig-xx
c     (2) sig-yy
c     (3) sig-zz
c     (4) tau-xy
c     (5) tau-xz
c     (6) tau-yz
c
c            this example implementation supports the umat for
c            bilinear plasticity with kinematic hardening.
c            see umat code for details.
c
c            **** modify code below for your umat -- or delete the   ****
c            **** assignment statements to mat_vals                  ****
c
c            statev(1:6)   = elastic strains
c            statev(7:12)  = plastic strains
c            statev(13:18) = backstress vector
c            statev(19)    = scalar plastic strain (eqpl)
c            statev(20)    = curretn loading state at point
c                            = 0.0 never yielded
c                            = 1.0 actively yielded this step
c                            = 2.0 previosuly yielding, but not this step
c
c
c                   local declarations
c                   ------------------
c
c
      logical debug
c
      debug = .false.
c      debug =  noel .eq. 10 .and. npt .eq. 1
c
      if( debug ) then
          call getnumcpus( numranks )
          call getrank( myrank )
          write(kout,9000) noel, npt, kinc, time, numranks, myrank
          write(kout,9010) stress(1:6)
          write(kout,9020) statev(1:nstatv)
          write(kout,9030) props(1:5)
          write(kout,9040)
      end if
c
      mat_vals(1) = statev(1)
      mat_vals(2) = statev(23)
      mat_vals(3) = statev(1)
c
      return
 9000 format(/,'.... debugging from umat_output ....',
     & /,10x,'noel:                ',i7,
     & /,10x,'npt:                 ',i7,
     & /,10x,'kinc (load step):    ',i7,' simulation time: ',e14.6,
     & /,10x,'numranks:            ',i7,
     & /,10x,'myrank:              ',i7 )
 9010 format(/,'   stresses: ',6e14.6 )
 9020 format(/,'   state variables: ', 30(/10x,6 e14.6) )
 9030 format(/,'   1st 5 umat properties: ',
     & /,10x,5e14.6 )
 9040 format(//)
c
       end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine uhard  (called by UMATs)          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 3/22/12                     *
c     *                                                              *
c     *                                                              *
c     ****************************************************************

      subroutine uhard( syield, hard, eqplas, eqplasrt, time, dtime,
     &     temp, dtemp, noel, npt, layer, kspt, kstep, kinc,
     &     cmname, nstatv, statev, numfieldv,
     &     predef, dpred, nvalue, table, kout )
      implicit double precision (a-h,o-z)
      parameter (nprecd=2)
c
      character(len=80) :: cmname
      dimension hard(3), statev(nstatv), time(*),
     1          predef(numfieldv), dpred(*)
c
      dimension table(2,nvalue)
c
      parameter(zero = 0.d0)
c
c            set yie      ld stress to last value of table,
c            hardening to zero
c
      syield  = table(1,nvalue)
      hard(1) = zero
c
c            if more than one entry, search table
c
      if( nvalue .gt. 1 ) then
         do k1 = 1, nvalue-1
           eqpl1 = table(2,k1+1)
           if( eqplas .lt. eqpl1 ) then
              eqpl0 = table(2,k1)
              if( eqpl1 .le .eqpl0 ) then
                  write(kout,100)
                  call xit
              endif
c
c            current yield stress and hardening
c
              deqpl = eqpl1 - eqpl0
              syiel0 = table(1,k1)
              syiel1 = table(1,k1+1)
              dsyiel = syiel1 - syiel0
              hard(1) = dsyiel / deqpl
              syield = syiel0 + ( eqplas-eqpl0 ) * hard(1)
              goto 10
            endif
           end do
10         continue
      endif
      return
c
 100  format(//,  30x,  '***error - plastic strain must be',
     1         ' entered in ascending order' )
c
      end
