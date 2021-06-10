      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if !if (flag_memory_estimate.eqv..true.) then

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti7_checkset_variable_0]'
      end if !if (verbose.gt.1) then
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else !if (flag_tesselation.eqv..false.) then
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if !if (flag_tesselation.eqv..false.) then

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_sum_,'n_w_sum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if !if ((svd_calculation_type.eq.1)) then

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if !if (flag_tesselation.eqv..true.) then

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    !if (flag_tesselation.eqv..true.) then

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      !if (flag_tesselation.eqv..false.) then
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    !if (flag_tesselation.eqv..false.) then
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 !if (flag_tesselation.eqv..false.) then
      end if                    !if ((flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else !if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if !if (flag_memory_checkset.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti7_checkset_variable_0]'
      end if !if (verbose.gt.1) then

      end if !if (flag_memory_estimate.eqv..false.) then
