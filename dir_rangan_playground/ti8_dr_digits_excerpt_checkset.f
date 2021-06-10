      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_dr_digits_checkset]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      call cxs_r8(ceiling(pi*n_delta_x*n_delta_y),delta_x_,'delta_x_'
     $     ,flag_memory_checkset)
      call cxs_r8(ceiling(pi*n_delta_x*n_delta_y),delta_y_,'delta_y_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_r,jacpts_x_,'jacpts_x_',flag_memory_checkset)
      call cxs_r8(quad_n_r,jacpts_w_,'jacpts_w_',flag_memory_checkset)
      call cxs_r8(quad_n_r,quad_grid_k_p_r_,'quad_grid_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_r,quad_weight_k_p_r_,'quad_weight_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_i4(quad_n_r,quad_n_polar_a_,'quad_n_polar_a_'
     $     ,flag_memory_checkset)
      call cxs_i4(quad_n_r,quad_n_w_,'quad_n_w_',flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_grid_k_p_0_,'quad_grid_k_p_0_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_grid_k_p_1_,'quad_grid_k_p_1_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_weight_k_p_,'quad_weight_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)
      call cxs_i4(n_S,I_S_sample_,'I_S_sample_',flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_S_index_,'S_alpha_S_index_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_polar_a_,'S_alpha_polar_a_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_azimu_b_,'S_alpha_azimu_b_'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_S*n_S_max,S_k_p__,'S_k_p__'
     $     ,flag_memory_checkset)
      call cxs_i4(n_M,I_M_sample_,'I_M_sample_',flag_memory_checkset)
      call cxs_r8(n_M,polar_a_est_,'polar_a_est_',flag_memory_checkset)
      call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_',flag_memory_checkset)
      call cxs_r8(n_M,delta_x_est_,'delta_x_est_',flag_memory_checkset)
      call cxs_r8(n_M,delta_y_est_,'delta_y_est_',flag_memory_checkset)
      call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_',flag_memory_checkset)
      call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_',flag_memory_checkset)
      call cxs_r8(n_M,S_index_est_,'S_index_est_',flag_memory_checkset)
      call cxs_c16(n_ctf*ld_CTF,CTF_k_p__,'CTF_k_p__'
     $     ,flag_memory_checkset)
      call cxs_i4(n_M,I_ctf_,'I_ctf_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_M,alpha_est__,'alpha_est__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_M*n_M_max,M_k_p__,'M_k_p__'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S_max,fix_phi_S_,'fix_phi_S_',flag_memory_checkset)
      call cxs_r8(n_S_max,fix_delta_S_x_c_0_,'fix_delta_S_x_c_0_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S_max,fix_delta_S_x_c_1_,'fix_delta_S_x_c_1_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S_max,fix_delta_S_x_p_r_,'fix_delta_S_x_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S_max,fix_delta_S_x_p_w_,'fix_delta_S_x_p_w_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_CTF,fix_phi_CTF_,'fix_phi_CTF_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_M_max,fix_delta_M_x_c_0_,'fix_delta_M_x_c_0_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_M_max,fix_delta_M_x_c_1_,'fix_delta_M_x_c_1_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_M_max,fix_delta_M_x_p_r_,'fix_delta_M_x_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_M_max,fix_delta_M_x_p_w_,'fix_delta_M_x_p_w_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_M,C_M_,'C_M_',flag_memory_checkset)
      if (flag_MS_vs_SM.eqv..true.) then
      call cxs_i4(n_S,n_MS_,'n_MS_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_MS_max*n_S,alpha_MS__,'alpha_MS__'
     $     ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      if (flag_MS_vs_SM.eqv..false.) then
      call cxs_i4(n_M,n_SM_,'n_SM_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_SM_max*n_M,alpha_SM__,'alpha_SM__'
     $     ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         call cxs_r8(n_MS_max*n_S,ferror_SRTRTCM_,'ferror_SRTRTCM_'
     $        ,flag_memory_checkset)
         call cxs_c16(n_MS_max*n_S,ferror_C_Z_opt_,'ferror_C_Z_opt_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_MS_max*n_S,ferror_CTF_R_S_,'ferror_CTF_R_S_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_MS_max*n_S,ferror_l2_norm_,'ferror_l2_norm_'
     $        ,flag_memory_checkset)
         call cxs_c16(n_MS_max*n_S,ferror_C_M_opt_,'ferror_C_M_opt_'
     $        ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      if (flag_MS_vs_SM.eqv..false.) then
         call cxs_r8(n_SM_max*n_M,ferror_SRTRTCM_,'ferror_SRTRTCM_'
     $        ,flag_memory_checkset)
         call cxs_c16(n_SM_max*n_M,ferror_C_Z_opt_,'ferror_C_Z_opt_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_SM_max*n_M,ferror_CTF_R_S_,'ferror_CTF_R_S_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_SM_max*n_M,ferror_l2_norm_,'ferror_l2_norm_'
     $        ,flag_memory_checkset)
         call cxs_c16(n_SM_max*n_M,ferror_C_M_opt_,'ferror_C_M_opt_'
     $        ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..false.) then

      if (flag_memory_checkset.eqv..true.) then
         write(6,'(A)') '[checkset passed]'
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if !if (verbose.gt.1) then
