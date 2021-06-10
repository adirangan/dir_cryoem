      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_dr_alignment_checkset]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      call cxs_i4(n_r,n_polar_a_,'n_polar_a_',flag_memory_checkset)
      call cxs_i4(n_r,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_r,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T0_k_p_,'T0_k_p_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T1_k_p_,'T1_k_p_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T2_k_p_,'T2_k_p_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T3_k_p_,'T3_k_p_',flag_memory_checkset)
      call cxs_i8(n_r,fftw_plan_frwd_,'fftw_plan_frwd_'
     $     ,flag_memory_checkset)
      call cxs_i8(n_r,fftw_plan_back_,'fftw_plan_back_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_w_sum,fftw_0in_,'fftw_0in_',flag_memory_checkset)
      call cxs_c16(n_w_sum,fftw_out_,'fftw_out_',flag_memory_checkset)
      call cxs_r8(2*n_r,grid_x_c_,'grid_x_c_',flag_memory_checkset)
      call cxs_r8(n_r,grid_x_p_,'grid_x_p_',flag_memory_checkset)
      call cxs_r8(2*n_r,grid_k_c_,'grid_k_c_',flag_memory_checkset)
      call cxs_r8(n_r,grid_k_p_r_,'grid_k_p_r_',flag_memory_checkset)
      call cxs_r8(n_r,weight_k_p_r_,'weight_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_w_sum,weight_k_p_,'weight_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(ceiling(pi*n_delta_x*n_delta_y),delta_x_,'delta_x_'
     $     ,flag_memory_checkset)
      call cxs_r8(ceiling(pi*n_delta_x*n_delta_y),delta_y_,'delta_y_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)
      call cxs_i4(n_S,I_S_sample_,'I_S_sample_',flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_S_index_,'S_alpha_S_index_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_polar_a_,'S_alpha_polar_a_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_azimu_b_,'S_alpha_azimu_b_'
     $     ,flag_memory_checkset)
      call cxs_c16(2*n_r*2*n_r*n_S_max,S_x_c__,'S_x_c__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_S*n_S_max,S_x_p__,'S_x_p__'
     $     ,flag_memory_checkset)
      call cxs_c16(2*n_r*2*n_r*n_S_max,S_k_c__,'S_k_c__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_S*n_S_max,S_k_p__,'S_k_p__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_S*n_S_max,S_k_q__,'S_k_q__'
     $     ,flag_memory_checkset)
      call cxs_i4(n_M,I_M_sample_,'I_M_sample_',flag_memory_checkset)
      call cxs_c16(2*n_r*2*n_r*n_M_max,N_x_c__,'N_x_c__'
     $     ,flag_memory_checkset)
      call cxs_c16(2*n_r*2*n_r*n_M_max,N_k_c__,'N_k_c__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_M*n_M_max,N_k_p__,'N_k_p__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_M*n_M_max,N_k_q__,'N_k_q__'
     $     ,flag_memory_checkset)
      call cxs_r8(n_M,polar_a_tru_,'polar_a_tru_',flag_memory_checkset)
      call cxs_r8(n_M,azimu_b_tru_,'azimu_b_tru_',flag_memory_checkset)
      call cxs_r8(n_M,delta_x_tru_,'delta_x_tru_',flag_memory_checkset)
      call cxs_r8(n_M,delta_y_tru_,'delta_y_tru_',flag_memory_checkset)
      call cxs_r8(n_M,gamma_z_tru_,'gamma_z_tru_',flag_memory_checkset)
      call cxs_r8(n_M,l2_norm_tru_,'l2_norm_tru_',flag_memory_checkset)
      call cxs_r8(n_M,S_index_tru_,'S_index_tru_',flag_memory_checkset)
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
      call cxs_r8(n_alpha*n_M,alpha_tru__,'alpha_tru__'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_M,alpha_est__,'alpha_est__'
     $     ,flag_memory_checkset)
      call cxs_c16(ld_M*n_M_max,M_k_p__,'M_k_p__'
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
      call cxs_r8(n_M,polar_a_upd_,'polar_a_upd_',flag_memory_checkset)
      call cxs_r8(n_M,azimu_b_upd_,'azimu_b_upd_',flag_memory_checkset)
      call cxs_r8(n_M,delta_x_upd_,'delta_x_upd_',flag_memory_checkset)
      call cxs_r8(n_M,delta_y_upd_,'delta_y_upd_',flag_memory_checkset)
      call cxs_r8(n_M,gamma_z_upd_,'gamma_z_upd_',flag_memory_checkset)
      call cxs_r8(n_M,l2_norm_upd_,'l2_norm_upd_',flag_memory_checkset)
      call cxs_r8(n_M,S_index_upd_,'S_index_upd_',flag_memory_checkset)
      call cxs_r8(n_M,polar_a_err_,'polar_a_err_',flag_memory_checkset)
      call cxs_r8(n_M,azimu_b_err_,'azimu_b_err_',flag_memory_checkset)
      call cxs_r8(n_M,delta_x_err_,'delta_x_err_',flag_memory_checkset)
      call cxs_r8(n_M,delta_y_err_,'delta_y_err_',flag_memory_checkset)
      call cxs_r8(n_M,gamma_z_err_,'gamma_z_err_',flag_memory_checkset)
      call cxs_r8(n_M,l2_norm_rer_,'l2_norm_rer_',flag_memory_checkset)
      call cxs_r8(n_M,S_index_err_,'S_index_err_',flag_memory_checkset)

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
