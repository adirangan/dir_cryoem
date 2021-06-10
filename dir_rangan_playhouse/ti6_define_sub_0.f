      n_SM_use_local_(nM_9_sub) = 0
      p_vp_input = nM_9_sub*3
      p_flag_S_use = nM_9_sub*nS_0_per
      p_LT = nM_9_sub*2*nS_0_per
      p_S_alpha_S_index_local = nM_9_sub*nS_1_per_max
      call cl1_r8(nS_1_per_max
     $     ,S_alpha_S_index_local_omp__(p_S_alpha_S_index_local))
      p_S_alpha_polar_a_local = nM_9_sub*nS_1_per_max
      call cl1_r8(nS_1_per_max
     $     ,S_alpha_polar_a_local_omp__(p_S_alpha_polar_a_local))
      p_S_alpha_azimu_b_local = nM_9_sub*nS_1_per_max
      call cl1_r8(nS_1_per_max
     $     ,S_alpha_azimu_b_local_omp__(p_S_alpha_azimu_b_local))
      p_I_S_sample_local = nM_9_sub*nS_1_per_max
      call cl1_i4(nS_1_per_max
     $     ,I_S_sample_local_omp__(p_I_S_sample_local))
      p_CTF_R_S_local = nM_9_sub*n_gamma_z*n_CTF*nS_1_per_max
      call cl1_c16(n_gamma_z*n_CTF*nS_1_per_max
     $     ,CTF_R_S_local_omp__(p_CTF_R_S_local))
      p_O_S_q_local = 0
      p_T_S_q_local = 0
      p_Z_S_q_local = 0
      p_S_T_T_R_CTF_M_q_local = 0
      p_S_Z_T_R_CTF_M_q_local = 0
      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.)) then
         p_O_S_q_local = nM_9_sub*n_r*nS_1_per_max*n_w_max
         call cl1_c16(n_r*nS_1_per_max*n_w_max
     $        ,O_S_q_local_omp__(p_O_S_q_local))
         p_S_T_T_R_CTF_M_q_local = nM_9_sub*n_transf*nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_transf*nS_1_per_max*n_w_max
     $        ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then
      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.)) then
         p_O_S_q_local = nM_9_sub*n_r*nS_1_per_max*n_w_max
         call cl1_c16(n_r*nS_1_per_max*n_w_max
     $        ,O_S_q_local_omp__(p_O_S_q_local))
         p_S_Z_T_R_CTF_M_q_local = nM_9_sub*n_transf*nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_transf*nS_1_per_max*n_w_max
     $        ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local))
         p_S_T_T_R_CTF_M_q_local = nM_9_sub*n_delta_v*nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_delta_v*nS_1_per_max*n_w_max
     $        ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      end if ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then
      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.)) then
         p_T_S_q_local = nM_9_sub*n_r*n_transf*nS_1_per_max*n_w_max
         call cl1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $        ,T_S_q_local_omp__(p_T_S_q_local))
         p_S_T_T_R_CTF_M_q_local = nM_9_sub*n_transf*nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_transf*nS_1_per_max*n_w_max
     $        ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.)) then
         p_Z_S_q_local = nM_9_sub*n_r*n_transf*nS_1_per_max*n_w_max
         call cl1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $        ,Z_S_q_local_omp__(p_Z_S_q_local))
         p_S_Z_T_R_CTF_M_q_local = nM_9_sub*n_transf*nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_transf*nS_1_per_max*n_w_max
     $        ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local))
         p_S_T_T_R_CTF_M_q_local = nM_9_sub*n_delta_v *nS_1_per_max
     $        *n_w_max
         call cl1_c16(n_delta_v*nS_1_per_max*n_w_max
     $        ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' nM_9_sub: ' , nM_9_sub
         write(6,'(A,I0)') ' p_vp_input: ' , p_vp_input
         write(6,'(A,I0)') ' p_LT: ' , p_LT
         write(6,'(A,I0)') ' p_flag_S_use: ' , p_flag_S_use
         write(6,'(A,I0)') ' p_S_alpha_S_index_local: ' ,
     $        p_S_alpha_S_index_local
         write(6,'(A,I0)') ' p_S_alpha_polar_a_local: ' ,
     $        p_S_alpha_polar_a_local
         write(6,'(A,I0)') ' p_S_alpha_azimu_b_local: ' ,
     $        p_S_alpha_azimu_b_local
         write(6,'(A,I0)') ' p_I_S_sample_local: ' , p_I_S_sample_local
         write(6,'(A,I0)') ' p_CTF_R_S_local: ' , p_CTF_R_S_local
         write(6,'(A,I0)') ' p_O_S_q_local: ' , p_O_S_q_local
         write(6,'(A,I0)') ' p_T_S_q_local: ' , p_T_S_q_local
         write(6,'(A,I0)') ' p_Z_S_q_local: ' , p_Z_S_q_local
         write(6,'(A,I0)') ' p_S_T_T_R_CTF_M_q_local: ' ,
     $        p_S_T_T_R_CTF_M_q_local
         write(6,'(A,I0)') ' p_S_Z_T_R_CTF_M_q_local: ' ,
     $        p_S_Z_T_R_CTF_M_q_local
      end if                    !if (verbose.gt.1) then
