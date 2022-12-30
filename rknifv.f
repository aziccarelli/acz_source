c     ****************************************************************
c     *                                                              *
c     *                      subroutine rknifv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 8/12/2017 rhd              *
c     *                                                              *
c     *     drives comptuation of internal force vectors for a       *
c           block of similar elements. integral trans B * sigma      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknifv( eleifv, updated_element_volumes, props,
     &                   nrow_ifv, local_work )
     &
      use segmental_curves, only : max_seg_points, max_seg_curves
      use constants
      implicit none
      include 'param_def'
      include 'include_sig_up'
c
c                       parameter declarations
c
      integer :: nrow_ifv    ! same as span
      real :: props(mxelpr,*)
      double precision :: eleifv(nrow_ifv,*), updated_element_volumes(*)
c
c                       local variables
c
      integer :: i, span, felem, elem_type, order, ngp, nnode, ndof,
     &           step, iter, mat_type, totdof, surf, gpn, iout
      double precision :: xi, eta, zeta, element_volumes_for_blk(mxvl),
     &                    bar_volumes(mxvl),
     &                    bar_areas_0(mxvl), bar_areas_n1(mxvl)
      logical :: geonl, bbar
      logical, parameter :: local_debug = .false.
c
      span            = local_work%span
      felem           = local_work%felem
      elem_type       = local_work%elem_type
      order           = local_work%int_order
      ngp             = local_work%num_int_points
      nnode           = local_work%num_enodes
      ndof            = local_work%num_enode_dof
      geonl           = local_work%geo_non_flg
      step            = local_work%step
      iter            = local_work%iter
      bbar            = local_work%bbar_flg
      mat_type        = local_work%mat_type
      totdof          = local_work%totdof
      surf            = local_work%surface
      iout            = local_work%iout
c
c                       initiialize ifv's for the block.
c
      eleifv(1:span,1:totdof)         = zero
      element_volumes_for_blk(1:span) = zero
c
c                       handle bar, link and std elements
c                       in service routines
c
      if( local_work%is_bar_elem ) then
        call rknifv_bar
      elseif( local_work%is_link_elem ) then
        call rknifv_link
      else
        call rknifv_std
      end if
c
c                       transform element internal force
c                       vectors to constraint compatible coordinates.
c
      if ( local_work%trn_e_block ) then
         do i = 1, span
           if (  local_work%trn_e_flags(i) ) then
             call trnvecs( eleifv, local_work%trnmte, local_work%trne,
     &                     ndof, nnode, i, 1, span )
           end if
         end do
      end if
c
c                       complete processing updated element volumes.
c                       if element has been killed, don't change its
c                       volume. we freeze the volume as the value
c                       when element is killed.
c
      if ( .not. local_work%is_cohes_elem ) then
        if ( local_debug )
     &    write(iout,*) '.. calling update_element_volumes ..'
       call update_element_volumes( felem, span,
     &                               updated_element_volumes,
     &                               element_volumes_for_blk )
      end if
c
      return
c
      contains
c     ========
c
      subroutine rknifv_bar
      implicit none
c
      integer :: i
      do i = 1, span
        bar_areas_0(i) = props(43,felem+i-1)
      end do
      call gtlsn3_vols( span, mxvl, felem, iout, bar_areas_0(1),
     &                  bar_areas_n1(1),
     &                  local_work%ce_0, local_work%ce_n1,
     &                  bar_volumes(1) )
      updated_element_volumes(1:span) = bar_volumes(1:span)
      call gpifv3( span, mxvl, felem, iout, local_work%elem_type,
     &             geonl, local_work%urcs_blk_n1(1,1,1),
     &             eleifv, local_work%ce_n1, bar_areas_0,
     &             bar_areas_n1 )
c
      return
      end subroutine rknifv_bar
c
c
      subroutine rknifv_link
      implicit none
c
      integer :: i
      if( local_debug ) write(iout,*) '>> eleifv for link2...'
      updated_element_volumes(1:span) = zero
      do i = 1, span
        eleifv(i,1) = -local_work%urcs_blk_n1(i,1,1)
        eleifv(i,2) =  local_work%urcs_blk_n1(i,1,1)
        eleifv(i,3) = -local_work%urcs_blk_n1(i,2,1)
        eleifv(i,4) =  local_work%urcs_blk_n1(i,2,1)
        eleifv(i,5) = -local_work%urcs_blk_n1(i,3,1)
        eleifv(i,6) =  local_work%urcs_blk_n1(i,3,1)
      end do
c
      return
      end subroutine rknifv_link
c
      subroutine rknifv_std
      implicit none
c
c                       compute all the shape function derivates, and
c                       inverse jacobians.  also calculate volume
c                       terms if using bbar. use element shape
c                       at n+1 for geonl.
c
      if( bbar .and. elem_type .eq. 2 ) then
        local_work%vol_block = zero
        local_work%volume_block = zero
      end if
c
c             this subroutine is called only for a block of
c             cohesive elements. here the element coordinates and
c             the global displacements are rotated to a coordinate
c             system in which the normal axis (Z rotated) is
c             perpendicular to the surface ot the cohesive element
c
      if( local_work%is_cohes_elem ) then
           call cohes_rot_mat( span, felem, nnode, elem_type,
     &                         local_work%ce_n1,
     &                         local_work%cohes_rot_block )
        if ( geonl )
     &      call cohes_mirror_refsurf( span, mxvl, totdof, nnode,
     &                                 local_work%ce_n1 )
      end if
c
      do gpn = 1, ngp
        call getgpts( elem_type, order, gpn, xi, eta, zeta,
     &                local_work%weights(gpn) )
        call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),
     &               local_work%neta(1,gpn),local_work%nzeta(1,gpn) )
        call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &    local_work%det_j(1,gpn), local_work%gama(1,1,1,gpn),
     &    local_work%cohes_rot_block,
     &    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &    local_work%nzeta(1,gpn), local_work%ce_n1, nnode )
c
        if ( local_work%is_cohes_elem )
     &     call shapef( elem_type, xi, eta, zeta,
     &                  local_work%shape(1,gpn) )
c
        if ( bbar .and. elem_type .eq. 2 ) then
          call vol_terms( local_work%gama(1,1,1,gpn),
     &                    local_work%det_j(1,gpn),
     &                    local_work%vol_block,
     &                    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                    local_work%nzeta(1,gpn),
     &                    local_work%volume_block, span, mxvl )
        end if
      end do
c
      if ( bbar .and. elem_type .eq. 2 )
     &  call vol_avg( local_work%vol_block, local_work%volume_block,
     &                span, mxvl )
c
c                       compute ifv's, one gauss point at
c                       a time for all elements in this block.
c                       volumetric type elements have their current
c                       volume computed and saved.
c
      do gpn = 1, ngp
         call gpifv1( eleifv, span, gpn, local_work%weights(gpn),
     &                local_work%det_j(1,gpn), local_work%b,
     &                local_work%urcs_blk_n1(1,1,gpn), local_work,
     &                element_volumes_for_blk )
      end do
c
      return
c
 9200 format(/,2x,'... rknifv for block with first element: ',i7)
 9300 format(5x,8e14.6)
c
      end subroutine rknifv_std
c
      end subroutine rknifv

c     ****************************************************************
c     *                                                              *
c     *              subroutine update_element_volumes               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/6/21 rhd                 *
c     *                                                              *
c     *                handle updating of element volumes            *
c     *                                                              *
c     ****************************************************************
c
      subroutine update_element_volumes( felem, span, 
     &                                   updated_element_volumes, 
     &                                   element_volumes_for_blk )
c
      use global_data, only : out
      use elem_extinct_data, only : dam_state, smcs_d_values
      use damage_data, only : dam_ptr, max_dam_state, growth_by_kill,
     &                        use_mesh_regularization
      use constants
c
      implicit none
c
c                       parameter declarations
c
      integer :: felem, span
      double precision :: updated_element_volumes(*), 
     &                    element_volumes_for_blk(*)
c
c                       local variables
c
      integer :: relem, element, elem_ptr
      logical, allocatable :: update_volume(:)
      logical :: simple_death
      logical, parameter :: local_debug = .false.
c
c                       if crack growth by element killing is
c                       *not* being used, all element volumes are
c                       updated.
c
      if( .not. growth_by_kill ) then
        updated_element_volumes(1:span) = 
     &       element_volumes_for_blk(1:span)
        return
      end if
c
c                       solution is using crack growth by element
c                       deletion.  build list of ones to do.
c                       decision based on mesh regularization:
c                          none - skip update if element being or done 
c                                 being killed
c                          yes - skip if damage already exceeds 1.0
c
      allocate( update_volume(1:span) )
      update_volume(1:span) = .false. ! all terms
      simple_death = .not. use_mesh_regularization 
c
      do relem = 1, span
       element     = felem + relem - 1
       elem_ptr    = dam_ptr(element)
       if( elem_ptr == 0 ) then  ! means element is not killable
         update_volume(relem) = .true.
         cycle
       end if
       if( simple_death ) then ! original element death
         if( dam_state(elem_ptr) > 0 ) cycle ! already killed
         update_volume(relem) = .true.
         cycle
       end if
       if( use_mesh_regularization ) then
         if( smcs_d_values(elem_ptr) < one ) 
     &               update_volume(relem) = .true.
         cycle
       end if
       write(out,9900)
       call die_abort
      end do
c
c                       update volume of non-killed elements.
c
      do relem = 1, span
       if( update_volume(relem) ) updated_element_volumes(relem) =
     &                         element_volumes_for_blk(relem)
      end do
c
      return
c
 9900 format(">>>>> FATAL ERROR. invalid condition in routine",
     & " update_element_volumes. Execution terminated.",//)
c
      end
