!> Doxygen comment: ;\n
!> check memory allocations for test_alpha_update_6.f ;\n
      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_alpha_update_6_checkset]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      call cxs_c16(n_M_sample,C_M_,'C_M_',flag_memory_checkset)
      call cxs_r8(n_azimu_b_polar_a_sum_cur,S_alpha_S_index_all_,
     $     'S_alpha_S_index_all_',flag_memory_checkset)
      call cxs_r8(n_azimu_b_polar_a_sum_cur,S_alpha_polar_a_all_,
     $     'S_alpha_polar_a_all_',flag_memory_checkset)
      call cxs_r8(n_azimu_b_polar_a_sum_cur,S_alpha_azimu_b_all_,
     $     'S_alpha_azimu_b_all_',flag_memory_checkset)
      call cxs_r8(n_k_p_upperbound,grid_cos_polar_a_tmp_,
     $     'grid_cos_polar_a_tmp_',flag_memory_checkset)
      call cxs_r8(n_k_p_upperbound,grid_sin_polar_a_tmp_,
     $     'grid_sin_polar_a_tmp_',flag_memory_checkset)
      call cxs_r8(n_k_p_upperbound,weight_cos_polar_a_tmp_,
     $     'weight_cos_polar_a_tmp_',flag_memory_checkset)
      call cxs_r8(n_k_p_upperbound,azimu_b_step_tmp_,
     $     'azimu_b_step_tmp_',flag_memory_checkset)
      call cxs_i4(n_k_p_upperbound,n_azimu_b_tmp_,'n_azimu_b_tmp_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_S_sample,I_S_sample_,'I_S_sample_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S_sample,S_alpha_S_index_sample_,
     $     'S_alpha_S_index_sample_',flag_memory_checkset)
      call cxs_r8(n_S_sample,S_alpha_polar_a_sample_,
     $     'S_alpha_polar_a_sample_',flag_memory_checkset)
      call cxs_r8(n_S_sample,S_alpha_azimu_b_sample_,
     $     'S_alpha_azimu_b_sample_',flag_memory_checkset)
      if (flag_MS_vs_SM.eqv..false.) then
         call cxs_i4(n_M_sample,n_SM_,'n_SM_',flag_memory_checkset)
         call cxs_r8(n_alpha*n_SM_max*n_M_sample,alpha_SM__,'alpha_SM__'
     $        ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         call cxs_i4(n_S_sample,n_MS_,'n_MS_',flag_memory_checkset)
         call cxs_r8(n_alpha*n_MS_max*n_S_sample,alpha_MS__,'alpha_MS__'
     $        ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      call cxs_r8(n_delta_v_use,delta_x_,'delta_x_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_delta_v_use,delta_y_,'delta_y_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_M_sample,C_M_,'C_M_',flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if !if (flag_memory_checkset.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_alpha_update_6_checkset]'
      end if !if (verbose.gt.1) then
