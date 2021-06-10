!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Defines tesselation tree for use in local search. ;\n
      if (0+verbose.gt.1) then
         write(6,'(A,F8.4)') ' tesselation_distance_req: ' ,
     $        tesselation_distance_req 
         write(6,'(A)') ' Constructing tesselation tree. '
      end if                    !if (0+verbose.gt.1) then
      
      n_point = nS_0_per;
      call cl1_r8(3*n_point,S_L_)
      do npoint=0,n_point-1
         ns0 = nS_0_sum + npoint
         call get_angle_to_vp_(S_alpha_polar_a_(ns0)
     $        ,S_alpha_azimu_b_(ns0),S_L_(0+3*npoint))
         call normalize_r8(3,S_L_(0 + 3*npoint))
      enddo                     !do npoint=0,n_point-1
      if (verbose.gt.1) then
         do npoint=0,n_point-1
            ns0 = nS_0_sum + npoint
            write(6,'(A,F8.4,F8.4,A,A,I0,1X,3F8.4)') 'S_alpha: ' ,
     $           S_alpha_polar_a_(ns0) , S_alpha_azimu_b_(ns0)
     $           , ' --> ' , ' S_L_: ', npoint , S_L_(0 + 3*npoint),
     $           S_L_(1 + 3 *npoint) , S_L_(2 + 3*npoint)
         enddo                  ! do npoint=0,n_point-1
      end if                    ! if (verbose.gt.1) then
      tradius_min = 1.0d-6
      call tesselation_get_nl_nm_ll(n_point,S_L_,tradius_min,nl_max
     $     ,nm_sum,ll_sum)
      if (verbose.gt.-1) then
         write(6,'(A,F8.6,4(A,I0))') ' tradius_min: ' , tradius_min
     $        ,' n_point: ' , n_point , ' nl_max: ' , nl_max
     $        ,' nm_sum: ' , nm_sum , ' ll_sum: ' , ll_sum
      end if                    !if (verbose.gt.1) then

       call cl1_i4(nm_sum,T_nl_) !level
       call cl1_r8(3*nm_sum,T_vm_) !vertex center
       call cl1_r8(nm_sum,T_tr_) !tradius
       call cl1_i4(nm_sum,T_ll_) !number of points from L_ in T_
       call cl1_l2(nm_sum,T_lf_) !is leaf
       call cl1_i4(nm_sum,T_c0_) !child_0 tesselation_index
       call cl1_i4(nm_sum,T_c1_) !child_1 tesselation_index
       call cl1_i4(nm_sum,T_c2_) !child_2 tesselation_index
       call cl1_i4(nm_sum,T_c3_) !child_3 tesselation_index
       call cl1_i4(nm_sum,T_ls_) !starting index of point_index_list for T_ if leaf (leaves only)
       call cl1_i4(ll_sum,T_LT_) !full point_index_list for all of T_ (leaves only)
       call cl1_i4(8,T_root_base_)

      call tesselation_index_wrapper_0(n_point,S_L_,tradius_min ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)

