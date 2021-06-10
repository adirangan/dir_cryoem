      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_ver15_checkset]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      call cxs_i4(n_k_p_max,nlats,'nlats',flag_memory_checkset)
      call cxs_i4(n_k_p_max,n_polar_a_,'n_polar_a_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_max,numonsphere,'numonsphere'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_max,n_azimu_b_polar_a_sum_
     $     ,'n_azimu_b_polar_a_sum_',flag_memory_checkset)
      call cxs_i4(n_k_p_max,nterms_sph,'nterms_sph'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_max,n_Y_l_,'n_Y_l_',flag_memory_checkset)
      call cxs_i4(n_k_p_max,isph_start,'isph_start'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_max,n_Y_lm_sum_,'n_Y_lm_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_lowerbound,ngridc,'ngridc'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_lowerbound,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_p_lowerbound,ngridc_tmp,'ngridc_tmp'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_max,icstart,'icstart',flag_memory_checkset)
      call cxs_i4(n_k_p_max,n_w_sum_,'n_w_sum_',flag_memory_checkset)
      call cxs_i4(n_k_p_lowerbound,ngridps,'ngridps'
     $     ,flag_memory_checkset)
      call cxs_i4(n_k_p_lowerbound,n_azimu_b_,'n_azimu_b_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_k_p_max,grid_k_p_,'grid_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_k_p_max,weight_k_p_,'weight_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,xnodesth,'xnodesth'
     $     ,flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,grid_cos_polar_a_
     $     ,'grid_cos_polar_a_',flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,sthetas,'sthetas'
     $     ,flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,grid_sin_polar_a_
     $     ,'grid_sin_polar_a_',flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,weight_cos_polar_a_
     $     ,'weight_cos_polar_a_',flag_memory_checkset)
      call cxs_r8(n_k_p_lowerbound,azimu_b_step_,'azimu_b_step_'
     $     ,flag_memory_checkset)
      n = 3*n_source_max
      call cxs_r8(n,source_x1_c_x2_c_x3_c_A__
     $     ,'source_x1_c_x2_c_x3_c_A__',flag_memory_checkset)
      if (flag_heterogeneity) then
         n = 3*n_source_max
         call cxs_r8(n,source_x1_c_x2_c_x3_c_B__
     $        ,'source_x1_c_x2_c_x3_c_B__',flag_memory_checkset)
         n = 3*n_source_max
         call cxs_r8(n,source_x1_c_x2_c_x3_c_C__
     $        ,'source_x1_c_x2_c_x3_c_C__',flag_memory_checkset)
      end if !      if (flag_heterogeneity) then
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_r8(n,grid_x1_c___,'grid_x1_c___',flag_memory_checkset)
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_r8(n,grid_x2_c___,'grid_x2_c___',flag_memory_checkset)
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_r8(n,grid_x3_c___,'grid_x3_c___',flag_memory_checkset)
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_c16(n,f_x_c_A___,'f_x_c_A___',flag_memory_checkset)
      if (flag_heterogeneity) then
         n = n_x_c_max*n_x_c_max*n_x_c_max
         call cxs_c16(n,f_x_c_B___,'f_x_c_B___',flag_memory_checkset)
         n = n_x_c_max*n_x_c_max*n_x_c_max
         call cxs_c16(n,f_x_c_C___,'f_x_c_C___',flag_memory_checkset)
      end if !if (flag_heterogeneity) then
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_c16(n,f_x_c_tru___,'f_x_c_tru___'
     $     ,flag_memory_checkset)
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cxs_c16(n,f_x_c_est___,'f_x_c_est___'
     $     ,flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_r8(n,grid_k1_c_,'grid_k1_c_',flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_r8(n,grid_k2_c_,'grid_k2_c_',flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_r8(n,grid_k3_c_,'grid_k3_c_',flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_r8(n,weight_k_c_,'weight_k_c_',flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_c16(n,f_k_c_A_,'f_k_c_A_',flag_memory_checkset)
      if (flag_heterogeneity) then
         n = n_azimu_b_polar_a_k_p_sum
         call cxs_c16(n,f_k_c_B_,'f_k_c_B_',flag_memory_checkset)
         n = n_azimu_b_polar_a_k_p_sum
         call cxs_c16(n,f_k_c_C_,'f_k_c_C_',flag_memory_checkset)
      end if !if (flag_heterogeneity) then
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_c16(n,f_k_c_tru_,'f_k_c_tru_',flag_memory_checkset)
      n = n_azimu_b_polar_a_k_p_sum
      call cxs_c16(n,f_k_c_est_,'f_k_c_est_',flag_memory_checkset)
      n = n_Y_lm_sum_max
      call cxs_c16(n,Y_A_,'Y_A_',flag_memory_checkset)
      if (flag_heterogeneity) then
         n = n_Y_lm_sum_max
         call cxs_c16(n,Y_B_,'Y_B_',flag_memory_checkset)
         n = n_Y_lm_sum_max
         call cxs_c16(n,Y_C_,'Y_C_',flag_memory_checkset)
      end if !if (flag_heterogeneity) then
      n = n_Y_lm_sum_max
      call cxs_c16(n,Y_tru_,'Y_tru_',flag_memory_checkset)
      n = n_Y_lm_sum_max
      call cxs_c16(n,Y_est_,'Y_est_',flag_memory_checkset)
      n = n_alpha*n_M
      call cxs_r8(n,alpha2d_tru__,'alpha2d_tru__'
     $     ,flag_memory_checkset)
      n = n_alpha*n_M
      call cxs_r8(n,alpha2d_est__,'alpha2d_est__'
     $     ,flag_memory_checkset)
      n = ld_M*n_M
      call cxs_c16(n,Y_slice_A__,'Y_slice_A__',flag_memory_checkset)
      if (flag_heterogeneity) then
         n = ld_M*n_M
         call cxs_c16(n,Y_slice_B__,'Y_slice_B__'
     $        ,flag_memory_checkset)
         n = ld_M*n_M
         call cxs_c16(n,Y_slice_B__,'Y_slice_B__'
     $        ,flag_memory_checkset)
      end if !if (flag_heterogeneity) then
      n = ld_M*n_M
      call cxs_c16(n,Y_slice__,'Y_slice__',flag_memory_checkset)
      call cxs_i4(n_M,I_model_,'I_model_',flag_memory_checkset)
      n = n_ctf*ld_CTF
      call cxs_c16(n,CTF_k_p__,'CTF_k_p__',flag_memory_checkset)
      call cxs_i4(n_M,I_ctf_,'I_ctf_',flag_memory_checkset)
      n = ld_S*n_azimu_b_polar_a_k_p_sum
      call cxs_c16(n,S_k_p__,'S_k_p__',flag_memory_checkset)
      n = n_M_sample_max
      call cxs_i4(n,I_M_sample_,'I_M_sample_',flag_memory_checkset)
      n = n_M_sample_max
      call cxs_i4(n,I_model_sample_,'I_model_sample_'
     $     ,flag_memory_checkset)
      n = ld_M*n_M_sample_max
      call cxs_c16(n,M_residual_sample__,'M_residual_sample__'
     $     ,flag_memory_checkset)
      n = n_residual_loading*n_M_sample_max
      call cxs_c16(n,M_residual_loading_,'M_residual_loading_'
     $     ,flag_memory_checkset)
      n = n_alpha*n_M_sample_max
      call cxs_r8(n,alpha_tru__,'alpha_tru__',flag_memory_checkset)
      n = n_alpha*n_M_sample_max
      call cxs_r8(n,alpha_est__,'alpha_est__',flag_memory_checkset)
      n = n_alpha*n_M_sample_max
      call cxs_r8(n,alpha_bkp__,'alpha_bkp__',flag_memory_checkset)

      if (flag_memory_checkset.eqv..true.) then
         write(6,'(A)') '[checkset passed]'
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if !if (flag_memory_checkset.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver15_checkset]'
      end if !if (verbose.gt.1) then
