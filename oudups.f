c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oudups                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 3/22/21   rhd              *          
c     *                                                              *          
c     *     gathers data for a block of elements to support          *          
c     *     generation of a patran, packet or hard copy output.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oudups( span, felem, ngp, geonl, stress, is_cohesive )  
c       
      use global_data, only : mxvl, nstr, nstrs, out
      use main_data, only : elems_to_blocks                                     
      use elem_block_data, only : history_blocks, rot_n1_blocks,                
     &                            eps_n_blocks, urcs_n_blocks,                  
     &                            history_blk_list                              
      use elblk_data, only : elem_hist, blk_size_hist, urcs_blk_n,              
     &                       rot_blk_n1, ddtse, blk_size_gp                     
c                                                                               
      implicit none                                                             
c                                                                               
      integer, intent(in) :: span, felem, ngp                                               
      logical, intent(in) :: geonl, stress, is_cohesive                                     
c                                                                               
c             local declarations                                                
c                                                                               
      integer :: blk, rel_elem, hist_size, hist_offset, rot_offset,             
     &           eps_offset, sig_offset, mxhist, mxngp                          
      logical :: is_solid                                                       
      logical, parameter ::  local_debug = .false.                              
c                                                                               
      blk         = elems_to_blocks(felem,1)                                    
      rel_elem    = elems_to_blocks(felem,2)                                    
      hist_size   = history_blk_list(blk)                                       
      hist_offset = (rel_elem-1)*hist_size*ngp + 1                              
      rot_offset  = (rel_elem-1)*9*ngp + 1                                      
      eps_offset  = (rel_elem-1)*nstr*ngp + 1                                   
      sig_offset  = (rel_elem-1)*nstrs*ngp + 1                                  
      is_solid    = .not. is_cohesive                                           
c                                                                               
      if( local_debug ) then                                                    
         write(out,*) '..... oudups....'                                        
         write(out,*) '....    span, felem, ngp, geonl, stress'                 
         write(out,*) span, felem, ngp, geonl, stress                           
         write(out,*) '.... blk, rel_elm, eps_offset: '                         
         write(out,*) blk, rel_elem, eps_offset                                 
         write(out,*) '.... is_cohesive: ', is_cohesive                         
      end if                                                                    
c                                                                               
c             gather history data. careful: the local                           
c             block size may be larger than stored block size.                  
c             uses non-standard gastr routine !                                 
c                                                                               
c             history data:                                                     
c              o The global blocks are sized(hist_size,ngp,span)                
c              o The local block is sized (mxvl,mxhist,mxngp).                  
c                -> mxhist: for all element blocks, the maximum                 
c                           no. of words of history data per                    
c                           gauss point                                         
c                -> mxngp:  for all elements blocks, the maximum                
c                           no. of integration points for an element            
c                                                                               
c              This makes it possible to pass a 2-D array slice for             
c              all elements of the block for a single gauss point.              
c                                                                               
      mxhist = blk_size_hist                                                    
      mxngp  = blk_size_gp                                                      
      call ou_gastr( elem_hist(1,1,1),                                          
     &               history_blocks(blk)%ptr(hist_offset),                      
     &               ngp, mxhist, mxngp, hist_size, span, mxvl )                
c                                                                               
c             gather stresses. for geonl, gather [Rot], the current             
c             rotation for transforming unrotated cauchy stresses               
c             to cauchy stresses (skip interface-cohesive elements).            
c             if not stresses, gather strain data.                              
c                                                                               
      if( geonl .and. is_solid )                                                
     &   call tanstf_gastr( rot_blk_n1,                                         
     &      rot_n1_blocks(blk)%ptr(rot_offset), ngp, 9, span )                  
c                                                                               
      if( stress ) then                                                         
        call tanstf_gastr( urcs_blk_n,                                          
     &         urcs_n_blocks(blk)%ptr(sig_offset), ngp, nstrs, span )           
      else                                                                      
        call tanstf_gastr( ddtse, eps_n_blocks(blk)%ptr(eps_offset),            
     &         ngp, nstr, span )                                                
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
