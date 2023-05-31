c     ****************************************************************          
c     *                                                              *          
c     *                   f-90 module damage_data                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *              last modified : 9/19/2021 rhd                   *          
c     *                                                              *          
c     *     define the variables and data structures to support      *          
c     *     crack growth using damage parameters (e.g., the Gurson   *          
c     *     model, ctoa, SMCS, etc.)                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      module damage_data                                                        
c                                                                               
      parameter (mxstp_store = 10)                                              
c                                                                               
c                     arrays                                                    
c                                                                               
      integer, save, allocatable :: dam_ptr(:),
     &          smcs_states_intlst(:), deleted_elist_to_stop(:),
     &          released_nlist_to_stop(:)
c                            
      double precision ::  del_poros(mxstp_store), regular_points(10,2), 
     &                     del_deff(mxstp_store)
c
      integer :: user_kill_list_now(100)      
c                                                                               
c                     scalar double precision/reals                             
c                                                                               
       double precision ::                                                         
     &   porosity_limit, gurson_cell_size,                                      
     &   crack_plane_coord, release_fraction,                                   
     &   critical_angle, release_height,                                        
     &   crack_plane_sign, char_length,                                         
     &   init_crit_ang, 
     &   smcs_gamma, smcs_beta_1, smcs_beta_2, smcs_a_plus,
     &   smcs_a_minus, smcs_kappa, smcs_alpha_1, smcs_alpha_2, 
     &   smcs_cutoff_triaxiality, max_eps_critical, smcs_type_4_A,
     &   smcs_type_4_n, smcs_type_4_c1, smcs_type_4_c2, smcs_type_4_c3,
     &   smcs_type_5_power, smcs_type_5_tp_critical, smcs_alpha,
     &   smcs_adapt_alpha_min, smcs_adapt_alpha_max, 
     &   smcs_beta, control_load_fact, old_load_fact,                                      
     &   min_load_fact, overshoot_limit, CTOA_range,                            
     &   perm_load_fact, max_porosity_change,                                   
     &   max_plast_strain_change,                                               
     &   init_ctoa_dist, ctoa_dist, distortion_plastic_limit,                                            
     &   crkpln_srch_tol, max_deff_change,                                      
     &   critical_cohes_deff_fract, 
     &   ppr_kill_displ_fraction, tolerance_mesh_regularization,
     &   regular_length, regular_up_max, regular_alpha,
     &   regular_Gf, regular_m_power, Oddy_critical_ratio
c                                                                               
c                     scalar integers                                           
c                                                                               
      integer :: crack_growth_type,                                                
     &   num_kill_elem, max_dam_state, csttail, num_print_list,                 
     &   num_kill_order_list, release_type, smcs_type,                                     
     &   crk_pln_normal_idx, num_crack_plane_nodes, crack_front_start,          
     &   crack_front_end, crkfrnt_garbage_start, crkfrnt_garbage_end,           
     &   min_steps_for_release, num_nodes_thick, num_crack_fronts,              
     &   num_nodes_back, num_nodes_grwinc, num_steps_min,                       
     &   num_elements_killed, stop_killed_elist_length,                                                  
     &   num_elements_in_force_release, num_ctoa_released_nodes,
     &   num_user_kill_elems, killed_element_limit, num_top_list,
     &   smcs_allowable_in_release, regular_type, regular_npoints,
     &   stop_released_nlist_length             
c                                                                               
c                     scalar logicals                                           
c                                                                               
      logical :: no_killed_elems, print_status, kill_order,                        
     &  kill_order_now, no_released_nodes, list_crkpln_nodes,                   
     &  list_crkfrnt_nodes, growth_by_kill, growth_by_release,                  
     &  growth_by_cohesive, enforce_node_release,                               
     &  overshoot_control_crk_grth, overshoot_allocated,                        
     &  load_size_control_crk_grth, g_stp_cntrl_allocated,                      
     &  const_front, master_lines_set, load_reduced, all_elems_killed,
     &  print_top_list, smcs_growth, smcs_states, smcs_stream,
     &  smcs_text, smcs_deleted_list_file_flag, use_estiff_at_death,
     &  use_mesh_regularization, use_distortion_metric,
     &  gt_list_file_flag, Oddy_print_initial,
     &  smcs_removed_list_file_flag         
c                                                                               
c                     character                             
c      
      character(len=80) :: smcs_deleted_list_file_name,  
     &                     smcs_removed_list_file_name,
     &                     gt_list_file_name        
c                                                                               
      end module                                                                
c
c     ****************************************************************
c     *                                                              *
c     *                   f-90 module mod_damage_ddczm               *
c     *                                                              *
c     *                       written by : vsp                       *
c     *                                                              *
c     *              last modified : 12/30/17 vsp                    *
c     *                                                              *
c     *     define the variables and data structures to support      *
c     *     crack growth using nonlocal lstar DDCZM                  *
c     *                                                              *
c     ****************************************************************
c
      module mod_damage_ddczm
c
c                     scalar double precision
c
      double precision :: 
c     &   swdm_c, swdm_kappa, swdm_lambda, swdm_beta, vgi_crit
     &   swdfm_c, swdfm_kappa, swdfm_beta, vgi_crit
c
c                     scalar integers
c
      integer :: ddczm_dmg_type
c
c                     scalar logicals
c
      logical :: ddczm_damage_on
c
c                     subroutines
c
      contains ! *** NOTE THIS ***
c     
c
c     ****************************************************************
c     *                                                              *
c     *                 function unsafe_division                     *
c     *                                                              *
c     *                       written by : vsp                       *
c     *                                                              *
c     *                   last modified: 03/02/18 vsp                *
c     *                                                              *
c     *    determines whether dividing two double-precision          *
c     *    floats will cause numerical overflow.                     *
c     *                                                              *
c     *    if function evaluates to .true. overflow will occur       *
c     *                                                              *
c     *    this follows a suggestion by William Long on              *
c     *    comp.lang.fortran google group                            *
c     *                                                              *
c     ****************************************************************
c
      pure logical function unsafe_division( num, den ) result( chk )
!DIR$ ATTRIBUTES INLINE :: unsafe_division
      implicit none
      double precision, intent(in) :: num, den
      chk = ( exponent(num) - exponent(den) >= maxexponent(num) ) .or.
     &      ( den==0 )
      end function unsafe_division
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine getmm_nonlocal_damage_index         *
c    *                                                              *
c    *     given a material model number, return the damage index   *
c    *                                                              *
c    ****************************************************************
c
      subroutine getmm_nonlocal_damage_index( iout, mm, indx )
!DIR$ ATTRIBUTES INLINE :: getmm_nonlocal_damage_index
      implicit none
      integer, intent(in)  :: iout, mm
      integer, intent(out) :: indx

      select case ( mm )
c     ----------------------
c     
      case ( 1, 5, 8 )
c     material model 1
c     material model 5
c     material model 8
      indx = 3
c
      case default 
c     all other material models
      write(iout,9000) mm 
      call die_abort
c
      end select
c
      return
c 
 9000 format("getmm_nonlocal_damage_index :: model ",i2," undefined")
      end subroutine getmm_nonlocal_damage_index
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine getmm_nonlocal_triax_index          *
c    *                                                              *
c    *     given a material model number, return the triax index    *
c    *                                                              *
c    ****************************************************************
c
      subroutine getmm_nonlocal_triax_index( iout, mm, indx )
!DIR$ ATTRIBUTES INLINE :: getmm_nonlocal_triax_index
      implicit none
      integer, intent(in)  :: iout, mm
      integer, intent(out) :: indx

      select case ( mm )
c     ----------------------
c     
      case ( 1, 5, 8 )
c     material model 1 
c     material model 5
c     material model 8
      indx = 4
c
      case default 
c     all other material models
      write(iout,9100) mm 
      call die_abort
c
      end select
c
      return
c 
 9100 format("getmm_nonlocal_triax_index :: model ",i2," undefined")
      end subroutine getmm_nonlocal_triax_index    
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine compute_invariants_ddczm            *
c    *                                                              *
c    *         written by : Vincente Pericoli                       *
c    *      last modified : 04/03/2018 VSP                          *
c    *                                                              *
c    *     subroutine to compute stress invariants                  *
c    *                                                              *
c    ****************************************************************
      subroutine compute_invariants_ddczm( stress, mises, triax,
     &                                     lodeang )
c 
!DIR$ ATTRIBUTES INLINE :: compute_invariants_ddczm
c
      implicit none 
      double precision, intent(in) :: stress(6)
      double precision, intent(out) :: mises, triax, lodeang
c
c               local definitions
c
      double precision :: 
     &    stress_hyd, stress_dev(6), sdev_inv2, sdev_inv3
      double precision, parameter :: 
     &    zero = 0.0d0, one = 1.0d0, two = 2.0d0, three = 3.0d0
c
c               Step 1: Update stress parameters
c
c
c     calculate stress deviator
c
      stress_hyd = (stress(1) + stress(2) + stress(3)) / three
      stress_dev(1) = stress(1) - stress_hyd
      stress_dev(2) = stress(2) - stress_hyd
      stress_dev(3) = stress(3) - stress_hyd
      stress_dev(4) = stress(4)
      stress_dev(5) = stress(5)
      stress_dev(6) = stress(6)
c
c     next, use stress deviator to calculate 
c     the deviatoric invariants:
c        sdev_inv2 = (1/2)*(S_ij * S_ij)
c        sdev_inv3 = (e_ijk)*(S_1i * S_2j * S_3k)
c            mises = sqrt( 3 * sdev_inv2 )
c
      sdev_inv2 = ( stress_dev(1)*stress_dev(1) 
     &            + stress_dev(2)*stress_dev(2) 
     &            + stress_dev(3)*stress_dev(3) ) / two
      sdev_inv2 = sdev_inv2 + stress_dev(4)*stress_dev(4)
     &                      + stress_dev(5)*stress_dev(5) 
     &                      + stress_dev(6)*stress_dev(6)
c
c     assuming: xx, yy, zz, xy, yz, xz (seems correct)
      sdev_inv3 = stress_dev(1)*stress_dev(2)*stress_dev(3) 
     &          + stress_dev(4)*stress_dev(5)*stress_dev(6)*two
     &          - stress_dev(1)*stress_dev(5)*stress_dev(5) 
     &          - stress_dev(2)*stress_dev(6)*stress_dev(6) 
     &          - stress_dev(3)*stress_dev(4)*stress_dev(4)
c
      mises = sqrt( three * sdev_inv2 )
c
c     now update stress triaxiality and normalized lode angle parameter.
c     note that: mises = sqrt( 3 * sdev_inv2 ).
c
      if ( unsafe_division( stress_hyd, mises ) ) then ! avoid divide by zero
         triax   = zero
         lodeang = zero
      else 
         triax   = stress_hyd / mises
         lodeang = 13.5d00 * sdev_inv3 / ( mises**3 ) ! 27/2 = 13.5
      end if
c 
      return 
c 
      end subroutine compute_invariants_ddczm
c
c
c
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine compute_swdm_ddczm                  *
c    *                                                              *
c    *         written by : Vincente Pericoli                       *
c    *      last modified : 08/09/2020 AJZ                          *
c    *                                                              *
c    *     subroutine to compute SWDM damage parameter              *
c    *                                                              *
c    ****************************************************************
c
      subroutine compute_swdm_ddczm( 
     1       stress_n1, mises, triax, peeq_n, peeq_n1, lodeang, 
     2       dmg_intgrnd_n, dmg_intgrnd_n1, dmg_intgrl_n, 
     3       dmg_intgrl_n1, peeq_comp_n, peeq_comp_n1,
c     4       damage, c_prop, kappa, lambda, beta)
     4       damage, c_prop, kappa, beta)
c
      implicit none
      double precision, intent(in) ::
     &     stress_n1(6), peeq_n, peeq_n1, peeq_comp_n, 
c     &     dmg_intgrnd_n, dmg_intgrl_n, c_prop, kappa, lambda, beta
     &     dmg_intgrnd_n, dmg_intgrl_n, c_prop, kappa, beta
      double precision, intent(out) ::
     &     mises, triax, lodeang, dmg_intgrnd_n1, dmg_intgrl_n1, 
     &     damage, peeq_comp_n1
c
c               local definitions
c
      double precision :: dpeeq
      double precision, parameter :: 
     &   zero = 0.0d00, half = 0.5d00, one = 1.0d00, 
     &   triaxfact = 1.3d00, dcrit = 0.6d00, four = 4.0d00,
     &   onepfive = 1.5d00, ten = 1.0d01, tpf=2.5d00
c
c     SUBROUTINE ARGUMENTS:
c     =====================
c     (input)         stress: stresses for solid element, used to compute stress 
c                             invariants. pass in as a vector of 6 components 
c                             xx, yy, zz, xy, yz, xz (any coordinate system).
c    (update)          mises: the 1D mises equivalent stress for current step.
c    (update)          triax: stress triaxiality for current step n+1.
c     (input)         peeq_n: the 1D mises equivalent plastic strain for previous step n.
c     (input)        peeq_n1: the 1D mises equivalent plastic strain for current step n+1.
c    (update)        lodeang: normalized lode angle parameter for current step n+1. 
c                              calculated from stress. see Malvern (1969).
c     (input)  dmg_intgrnd_n: the SWDM damage integrand evaluated at previous step n.
c    (update) dmg_intgrnd_n1: the SWDM damage integrand evaluated at current step n+1.
c     (input)   dmg_intgrl_n: the SWDM damage integral evaluated from step 0 to n.
c    (update)  dmg_intgrl_n1: the SWDM damage integral evaluated from step 0 to n+1.
c     (input)    peeq_comp_n: the total compressive PEEQ at step n 
c    (update)   peeq_comp_n1: the total compressive PEEQ at step n+1
c    (update)         damage: the SWDM damage parameter.
c     (input)         c_prop: material property; ductility parameter.
c     (input)          kappa: material property; lode angle influence parameter.
c     (input)         lambda: material property; stress-independent degredation parameter.
c
c
c               Step 1: compute stress invariants.
c                       also need triax for step n. 
c
      call compute_invariants_ddczm( stress_n1, mises, triax, lodeang )
c
c
c               Step 2: Calculate the total change in PEEQ over step. 
c                       Update compressive PEEQ fraction. Note that
c                       dPEEQ is either fully compressive or fully 
c                       tensile (transition region is elastic). 
c
      dpeeq = (peeq_n1 - peeq_n) ! total change in PEEQ over step
c 
      if( triax .lt. zero ) then 
        peeq_comp_n1 = peeq_comp_n + dpeeq
      else 
        peeq_comp_n1 = peeq_comp_n 
      endif
c
c
c
c
c               Step 3: Calculate SWDM Damage
c
c
c     evaluate the integral using trap rule
c
      if ((triax .gt. -one*four) .and. (triax .lt. four)) then
         dmg_intgrnd_n1 = exp(kappa * (abs(lodeang) - one) ) *
     &                ( exp(triaxfact*triax) - (one/beta)*
     &                            exp(-triaxfact*triax) )
c
         dmg_intgrl_n1 = half * (dmg_intgrnd_n + dmg_intgrnd_n1) * dpeeq
         dmg_intgrl_n1 = dmg_intgrl_n1 + dmg_intgrl_n ! increment + previous
c
      else
         dmg_intgrl_n1 = dmg_intgrl_n
      end if
c
c 
      if( dmg_intgrl_n1 .lt. zero ) dmg_intgrl_n1 = zero 
c 
c      damage = dmg_intgrl_n1 * exp( lambda * peeq_comp_n1 ) * c_prop
      damage = dmg_intgrl_n1 * c_prop
c
c
c
c
      return
      end subroutine compute_swdm_ddczm
c     
c    ****************************************************************
c    *                                                              *
c    *               subroutine compute_vgm_ddczm                   *
c    *                                                              *
c    *         written by : Vincente Pericoli                       *
c    *      last modified : 06/29/2017 VSP                          *
c    *                                                              *
c    *     subroutine to compute monotonic VGM damage parameter     *
c    *                                                              *
c    ****************************************************************
c
      subroutine compute_vgm_ddczm( stress, mises, triax_n1,
     2                             peeq_n, peeq_n1,
     3                             vgi_intgrnd_n, vgi_intgrnd_n1,
     4                             vgi_n, vgi_n1, damage, vgi_crit )
c
      implicit none
      double precision, intent(in) ::
     &     stress(6), peeq_n, peeq_n1, vgi_intgrnd_n,
     &     vgi_n, vgi_crit
      double precision, intent(out) ::
     &     mises, triax_n1, vgi_intgrnd_n1, vgi_n1, damage
c
c               local definitions
c
      double precision :: dpeeq, dumd
      double precision, parameter :: 
     &     zero = 0.0d00, half = 0.5d00, triaxfact = 1.5d00
c
c     SUBROUTINE ARGUMENTS:
c     =====================
c     (input)         stress: stresses for solid element, used to compute stress 
c                             invariants. pass in as a vector of 6 components 
c                             xx, yy, zz, xy, yz, xz (any coordinate system).
c     (input)          mises: the 1D mises equivalent stress for current step.
c    (update)       triax_n1: stress triaxiality for current step n+1.
c     (input)         peeq_n: the 1D mises equivalent plastic strain for previous step n.
c     (input)        peeq_n1: the 1D mises equivalent plastic strain for current step n+1.
c     (input)  vgi_intgrnd_n: the VGI integrand evaluated at previous step n.
c    (update) vgi_intgrnd_n1: the VGI integrand evaluated at current step n+1.
c     (input)          vgi_n: the VGI evaluated from step 0 to n.
c    (update)         vgi_n1: the VGI evaluated from step 0 to n+1.
c     (input)       vgi_crit: the critical VGI (material property)
c    (update)         damage: the VGM damage parameter, vgi_n1 / vgi_crit
c
c               Step 1: calculate stress invariants. Only triax needed. 
c
      call compute_invariants_ddczm( stress, mises, triax_n1, dumd )
c
c               Step 2: evaluate the VGI, use trap rule
c
      if (triax_n1 .lt. zero) then ! monotonic VGI is invalid for compression
         vgi_n1 = vgi_n
      else 
c
         dpeeq = (peeq_n1 - peeq_n)
         vgi_intgrnd_n1 = exp(triaxfact * triax_n1)
c
         vgi_n1 = half * (vgi_intgrnd_n1 + vgi_intgrnd_n) * dpeeq ! increment
         vgi_n1 = vgi_n1 + vgi_n ! increment + previous
      end if
c
c               Step 3: calculate damage
c
      damage = vgi_n1 / vgi_crit
c
      return
      end subroutine compute_vgm_ddczm
c
c
      end module mod_damage_ddczm
