      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_bestfitall_excerpt_checkset]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      call cxs_i4(n_r,n_polar_a_,'n_polar_a_',flag_memory_checkset)
      call cxs_i4(n_r,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_r,weight_k_p_r_,'weight_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_w_sum,weight_k_p_,'weight_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_delta_v,delta_x_,'delta_x_',flag_memory_checkset)
      call cxs_r8(n_delta_v,delta_y_,'delta_y_',flag_memory_checkset)
      call cxs_i4(n_S,I_S_sample_,'I_S_sample_',flag_memory_checkset)
      call cxs_r8(n_S,S_alpha_S_index_,'S_alpha_S_index_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_M,I_M_sample_,'I_M_sample_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_M,alpha_est__,'alpha_est__'
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
      call cxs_c16(n_w_sum,S_k_p_,'S_k_p_',flag_memory_checkset)
      call cxs_c16(n_w_sum,CTF_k_p_,'CTF_k_p_',flag_memory_checkset)
      call cxs_c16(n_w_sum,R_S_,'R_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T_R_S_,'T_R_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,R_S_,'T_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,T_R_S_,'R_T_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,CTF_R_S_,'CTF_R_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,CTF_T_R_S_,'CTF_T_R_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,CTF_T_R_S_,'CTF_R_T_S_',flag_memory_checkset)
      call cxs_c16(n_w_sum,M_k_p_,'M_k_p_',flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.0) then
         write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.0) then
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_bestfitall_excerpt_checkset]'
      end if !if (verbose.gt.1) then
