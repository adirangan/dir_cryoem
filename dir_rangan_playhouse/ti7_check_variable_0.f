      if (flag_memory_estimate.eqv..true.) then
         goto 10
      end if !if (flag_memory_estimate.eqv..true.) then

      write(6,'(A)') '[entering ti7_check_variable_0]'

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else !if (flag_tesselation.eqv..false.) then
         write(6,'(A,I0)') ' n_SM_use_local_ ' , n_M
         call cp1_i4(n_M,n_SM_use_local_,n_SM_use_local_)
         write(6,'(A,I0)') ' S_L_ ' , 3*n_point
         call cp1_r8(3*n_point,S_L_,S_L_)
         write(6,'(A,I0)') ' T_nl_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_nl_,T_nl_)
         write(6,'(A,I0)') ' T_vm_ ' , 3*nm_sum
         call cp1_r8(3*nm_sum,T_vm_,T_vm_)
         write(6,'(A,I0)') ' T_tr_ ' , 1*nm_sum
         call cp1_r8(1*nm_sum,T_tr_,T_tr_)
         write(6,'(A,I0)') ' T_ll_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_ll_,T_ll_)
         write(6,'(A,I0)') ' T_lf_ ' , 1*nm_sum
         call cp1_l2(1*nm_sum,T_lf_,T_lf_)
         write(6,'(A,I0)') ' T_c0_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_c0_,T_c0_)
         write(6,'(A,I0)') ' T_c1_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_c1_,T_c1_)
         write(6,'(A,I0)') ' T_c2_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_c2_,T_c2_)
         write(6,'(A,I0)') ' T_c3_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_c3_,T_c3_)
         write(6,'(A,I0)') ' T_ls_ ' , 1*nm_sum
         call cp1_i4(1*nm_sum,T_ls_,T_ls_)
         write(6,'(A,I0)') ' T_LT_ ' , 1*ll_sum
         call cp1_i4(1*ll_sum,T_LT_,T_LT_)
         write(6,'(A,I0)') ' T_root_base_ ' , 8
         call cp1_i4(8,T_root_base_,T_root_base_)
      end if !if (flag_tesselation.eqv..false.) then

      write(6,'(A,I0)') ' n_w_ ' , n_k_cur
      call cp1_i4(n_k_cur,n_w_,n_w_)
      write(6,'(A,I0)') ' n_w_sum_ ' , n_k_cur
      call cp1_i4(n_k_cur,n_w_sum_,n_w_sum_)

      write(6,'(A,I0)') ' gamma_z_ ' , n_gamma_z
      call cp1_r8(n_gamma_z,gamma_z_,gamma_z_)

      if ((svd_calculation_type.eq.1)) then
         write(6,'(A,I0)') ' svd_r_ ' , n_svd_r
         call cp1_r8(n_svd_r,svd_r_,svd_r_)
         write(6,'(A,I0)') ' svd_d_ ' , n_svd_d
         call cp1_r8(n_svd_d,svd_d_,svd_d_)
         write(6,'(A,I0)') ' svd_l_ ' , n_svd_l
         call cp1_i4(n_svd_l,svd_l_,svd_l_)
         write(6,'(A,I0)') ' svd_U_d_ ' , n_svd_d*n_svd_l
         call cp1_r8(n_svd_d*n_svd_l,svd_U_d_,svd_U_d_)
         write(6,'(A,I0)') ' svd_s_ ' , n_svd_l
         call cp1_r8(n_svd_l,svd_s_,svd_s_)
         write(6,'(A,I0)') ' svd_V_r_ ' , n_svd_r*n_svd_l
         call cp1_r8(n_svd_r*n_svd_l,svd_V_r_,svd_V_r_)
         write(6,'(A,I0)') ' svd_polyval_U_d_ ' , n_svd_l*n_delta_v
         call cp1_r8(n_svd_l*n_delta_v,svd_polyval_U_d_
     $        ,svd_polyval_U_d_)
         write(6,'(A,I0)') ' svd_polyval_V_r_ ' , n_svd_l*n_r
         call cp1_r8(n_svd_l*n_r,svd_polyval_V_r_,svd_polyval_V_r_)
      end if !if ((svd_calculation_type.eq.1)) then

      write(6,'(A,I0)') ' S_p_ ' , n_A
      call cp1_c16(n_A,S_p_,S_p_)
      write(6,'(A,I0)') ' S_q_ ' , n_A
      call cp1_c16(n_A,S_q_,S_q_)
      write(6,'(A,I0)') ' S_p_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,S_p_omp__,S_p_omp__)
      write(6,'(A,I0)') ' S_q_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,S_q_omp__,S_q_omp__)
      write(6,'(A,I0)') ' M_p_ ' , n_A
      call cp1_c16(n_A,M_p_,M_p_)
      write(6,'(A,I0)') ' M_q_ ' , n_A
      call cp1_c16(n_A,M_q_,M_q_)
      write(6,'(A,I0)') ' M_p_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,M_p_omp__,M_p_omp__)
      write(6,'(A,I0)') ' M_q_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,M_q_omp__,M_q_omp__)
      write(6,'(A,I0)') ' CTF_p_ ' , n_A
      call cp1_c16(n_A,CTF_p_,CTF_p_)
      write(6,'(A,I0)') ' CTF_q_ ' , n_A
      call cp1_c16(n_A,CTF_q_,CTF_q_)
      write(6,'(A,I0)') ' CTF_p_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,CTF_p_omp__,CTF_p_omp__)
      write(6,'(A,I0)') ' CTF_q_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,CTF_q_omp__,CTF_q_omp__)
      write(6,'(A,I0)') ' Z_p_ ' , n_A
      call cp1_c16(n_A,Z_p_,Z_p_)
      write(6,'(A,I0)') ' Z_q_ ' , n_A
      call cp1_c16(n_A,Z_q_,Z_q_)
      write(6,'(A,I0)') ' Z_p_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,Z_p_omp__,Z_p_omp__)
      write(6,'(A,I0)') ' Z_q_omp__ ' , n_A*n_omp_sub__in
      call cp1_c16(n_A*n_omp_sub__in,Z_q_omp__,Z_q_omp__)

      write(6,'(A,I0)') ' ZZ__ ' , n_delta_v*n_gamma_z
      call cp1_c16(n_delta_v*n_gamma_z,ZZ__,ZZ__)
      write(6,'(A,I0)') ' ZZ_omp__ ' , n_delta_v*n_gamma_z*n_omp_sub__in
      call cp1_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__,ZZ_omp__)
      write(6,'(A,I0)') ' ZZ_sub_ ' , n_w_max
      call cp1_c16(n_w_max,ZZ_sub_,ZZ_sub_)
      write(6,'(A,I0)') ' ZZ_sub_omp__ ' , n_w_max*n_omp_sub__in
      call cp1_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__,ZZ_sub_omp__)
      write(6,'(A,I0)') ' CTF_R_S_sub__ ' , n_gamma_z
      call cp1_c16(n_gamma_z,CTF_R_S_sub__,CTF_R_S_sub__)
      write(6,'(A,I0)') ' CTF_R_S_sub_omp__ ' , n_gamma_z*n_omp_sub__in
      call cp1_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,CTF_R_S_sub_omp__)

      write(6,'(A,I0)') ' n_S_0_per_ ' , n_S
      call cp1_i4(n_S,n_S_0_per_,n_S_0_per_)
      write(6,'(A,I0)') ' n_S_0_sum_ ' , n_S
      call cp1_i4(n_S,n_S_0_sum_,n_S_0_sum_)

      write(6,'(A,I0)') ' CTF_R_S__ ' , n_gamma_z*n_CTF*nS_0_per_max
      call cp1_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__,CTF_R_S__)

      write(6,'(A,I0)') ' n_S_1_per_ ' , nS_0_per_max
      call cp1_i4(nS_0_per_max,n_S_1_per_,n_S_1_per_)
      write(6,'(A,I0)') ' n_S_1_sum_ ' , nS_0_per_max
      call cp1_i4(nS_0_per_max,n_S_1_sum_,n_S_1_sum_)

      write(6,'(A,I0)') ' C_trn0_ ' , nS_1_per_max*n_transf
      call cp1_c16(nS_1_per_max*n_transf,C_trn0_,C_trn0_)
      write(6,'(A,I0)') ' C_trn0_omp__ ' , nS_1_per_max*n_transf
     $     *n_omp_sub__in
      call cp1_c16(nS_1_per_max*n_transf*n_omp_sub__in,C_trn0_omp__
     $     ,C_trn0_omp__)
      if (flag_tesselation.eqv..true.) then
         write(6,'(A,I0)') ' CTF_R_S_local_omp__ ' , n_gamma_z*n_CTF
     $        *nS_1_per_max*n_omp_sub__in
         call cp1_c16(n_gamma_z*n_CTF*nS_1_per_max*n_omp_sub__in
     $        ,CTF_R_S_local_omp__,CTF_R_S_local_omp__)
      end if !if (flag_tesselation.eqv..true.) then

      if (flag_tesselation.eqv..true.) then
         write(6,'(A,I0)') ' vp_input_omp__ ' , 3*n_omp_sub__in
         call cp1_r8(3*n_omp_sub__in,vp_input_omp__,vp_input_omp__)
         write(6,'(A,I0)') ' flag_S_use_omp__ ' , nS_0_per_max
     $        *n_omp_sub__in
         call cp1_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,flag_S_use_omp__)
         write(6,'(A,I0)') ' LT_omp__ ' , 2*nS_0_per_max*n_omp_sub__in
         call cp1_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__,LT_omp__)
         write(6,'(A,I0)') ' S_alpha_S_index_local_omp__ ' ,
     $        nS_1_per_max*n_omp_sub__in
         call cp1_r8(nS_1_per_max*n_omp_sub__in
     $        ,S_alpha_S_index_local_omp__,S_alpha_S_index_local_omp__)
         write(6,'(A,I0)') ' S_alpha_polar_a_local_omp__ ' ,
     $        nS_1_per_max*n_omp_sub__in
         call cp1_r8(nS_1_per_max*n_omp_sub__in
     $        ,S_alpha_polar_a_local_omp__,S_alpha_polar_a_local_omp__)
         write(6,'(A,I0)') ' S_alpha_azimu_b_local_omp__ ' ,
     $        nS_1_per_max*n_omp_sub__in
         call cp1_r8(nS_1_per_max*n_omp_sub__in
     $        ,S_alpha_azimu_b_local_omp__,S_alpha_azimu_b_local_omp__)
      end if                    !if (flag_tesselation.eqv..true.) then

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         write(6,'(A,I0)') ' O_S_q__ ' , n_r*nS_0_per_max*n_w_max
         call cp1_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,O_S_q__)
         if (flag_tesselation.eqv..true.) then
            write(6,'(A,I0)') ' O_S_q_local_omp__ ' , n_r*nS_1_per_max
     $           *n_w_max*n_omp_sub__in
            call cp1_c16(n_r*nS_1_per_max*n_w_max*n_omp_sub__in
     $           ,O_S_q_local_omp__,O_S_q_local_omp__)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         write(6,'(A,I0)') ' T_S_q__ ' , n_r*n_transf*nS_0_per_max
     $        *n_w_max
         call cp1_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__,T_S_q__)
         if (flag_tesselation.eqv..true.) then
            write(6,'(A,I0)') ' T_S_q_local_omp__ ' , n_r*n_transf
     $           *nS_1_per_max*n_w_max*n_omp_sub__in
            call cp1_c16(n_r*n_transf*nS_1_per_max *n_w_max
     $           *n_omp_sub__in ,T_S_q_local_omp__,T_S_q_local_omp__)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         write(6,'(A,I0)') ' Z_S_q__ ' , n_r*n_transf*nS_0_per_max
     $        *n_w_max
         call cp1_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__,Z_S_q__)
         if (flag_tesselation.eqv..true.) then
            write(6,'(A,I0)') ' Z_S_q_local_omp__ ' , n_r*n_transf
     $           *nS_1_per_max*n_w_max*n_omp_sub__in
            call cp1_c16(n_r*n_transf*nS_1_per_max *n_w_max
     $           *n_omp_sub__in ,Z_S_q_local_omp__,Z_S_q_local_omp__)
         end if !if (flag_tesselation.eqv..true.) then
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      write(6,'(A,I0)') ' n_M_0_per_ ' , n_M
      call cp1_i4(n_M,n_M_0_per_,n_M_0_per_)
      write(6,'(A,I0)') ' n_M_0_sum_ ' , n_M
      call cp1_i4(n_M,n_M_0_sum_,n_M_0_sum_)

      if (flag_tesselation.eqv..false.) then
         write(6,'(A,I0)') ' polar_a_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,polar_a_est_,polar_a_est_)
         write(6,'(A,I0)') ' azimu_b_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,azimu_b_est_,azimu_b_est_)
         write(6,'(A,I0)') ' gamma_z_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,gamma_z_est_,gamma_z_est_)
         write(6,'(A,I0)') ' delta_x_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,delta_x_est_,delta_x_est_)
         write(6,'(A,I0)') ' delta_y_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,delta_y_est_,delta_y_est_)
         write(6,'(A,I0)') ' l2_norm_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,l2_norm_est_,l2_norm_est_)
         write(6,'(A,I0)') ' ctf_ind_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,ctf_ind_est_,ctf_ind_est_)
         write(6,'(A,I0)') ' S_index_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,S_index_est_,S_index_est_)
         write(6,'(A,I0)') ' M_index_est_ ' , nM_0_per_max
         call cp1_r8(nM_0_per_max,M_index_est_,M_index_est_)
      else                      !if (flag_tesselation.eqv..false.) then
         write(6,'(A,I0)') ' polar_a_est_ ' , n_M
         call cp1_r8(n_M,polar_a_est_,polar_a_est_)
         write(6,'(A,I0)') ' azimu_b_est_ ' , n_M
         call cp1_r8(n_M,azimu_b_est_,azimu_b_est_)
         write(6,'(A,I0)') ' gamma_z_est_ ' , n_M
         call cp1_r8(n_M,gamma_z_est_,gamma_z_est_)
         write(6,'(A,I0)') ' delta_x_est_ ' , n_M
         call cp1_r8(n_M,delta_x_est_,delta_x_est_)
         write(6,'(A,I0)') ' delta_y_est_ ' , n_M
         call cp1_r8(n_M,delta_y_est_,delta_y_est_)
         write(6,'(A,I0)') ' l2_norm_est_ ' , n_M
         call cp1_r8(n_M,l2_norm_est_,l2_norm_est_)
         write(6,'(A,I0)') ' ctf_ind_est_ ' , n_M
         call cp1_r8(n_M,ctf_ind_est_,ctf_ind_est_)
         write(6,'(A,I0)') ' S_index_est_ ' , n_M
         call cp1_r8(n_M,S_index_est_,S_index_est_)
         write(6,'(A,I0)') ' M_index_est_ ' , n_M
         call cp1_r8(n_M,M_index_est_,M_index_est_)
      end if                    !if (flag_tesselation.eqv..false.) then
      write(6,'(A,I0)') ' alpha__in_ ' , n_alpha
      call cp1_r8(n_alpha,alpha__in_,alpha__in_)
      write(6,'(A,I0)') ' alpha__in_omp__ ' , n_alpha*n_omp_sub__in
      call cp1_r8(n_alpha*n_omp_sub__in,alpha__in_omp__,alpha__in_omp__)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_r*nM_0_per_max*n_w_max
            write(6,'(A,I0)') ' O_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,O_T_R_CTF_M_q__,O_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_r*n_w_max*n_omp_sub__in
            write(6,'(A,I0)') ' O_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,O_T_R_CTF_M_q_local_omp__
     $           ,O_T_R_CTF_M_q_local_omp__)
         end if                 !if (flag_tesselation.eqv..false.) then
      end if                    !if ((flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_r*n_transf*nM_0_per_max*n_w_max
            write(6,'(A,I0)') ' T_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,T_T_R_CTF_M_q__,T_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_r*n_transf*n_w_max*n_omp_sub__in
            write(6,'(A,I0)') ' T_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,T_T_R_CTF_M_q_local_omp__
     $           ,T_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_r*n_transf*nM_0_per_max*n_w_max
            write(6,'(A,I0)') ' Z_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,Z_T_R_CTF_M_q__,Z_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_r*n_transf*n_w_max*n_omp_sub__in
            write(6,'(A,I0)') ' Z_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,Z_T_R_CTF_M_q_local_omp__
     $           ,Z_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      write(6,'(A,I0)') ' n_M_1_per_ ' , nM_0_per_max
      call cp1_i4(nM_0_per_max,n_M_1_per_,n_M_1_per_)
      write(6,'(A,I0)') ' n_M_1_sum_ ' , nM_0_per_max
      call cp1_i4(nM_0_per_max,n_M_1_sum_,n_M_1_sum_)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_transf*nS_1_per_max*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q_local_omp__
     $           ,S_T_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = nS_1_per_max*n_transf*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q_local_omp__
     $           ,S_T_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         write(6,'(A,I0)') ' Z_S_svdd_ ' , na
         call cp1_c16(na,Z_S_svdd_,Z_S_svdd_)
         if (flag_tesselation.eqv..false.) then
            na = n_transf*nS_1_per_max*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_Z_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_Z_T_R_CTF_M_q__,S_Z_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_Z_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_Z_T_R_CTF_M_q_local_omp__
     $           ,S_Z_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
         if (flag_tesselation.eqv..false.) then
            na = n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q_local_omp__
     $           ,S_T_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         allocate(Z_M_svdd_(0:n_svd_l*n_delta_v-1))
         if (flag_tesselation.eqv..false.) then
            na = nS_1_per_max*n_transf*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_Z_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_Z_T_R_CTF_M_q__,S_Z_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_Z_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_Z_T_R_CTF_M_q_local_omp__
     $           ,S_Z_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
         if (flag_tesselation.eqv..false.) then
            na = n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
         else !if (flag_tesselation.eqv..false.) then
            na = n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $           *n_omp_sub__in
            write(6,'(A,I0)') ' S_T_T_R_CTF_M_q_local_omp__ ' , na
            call cp1_c16(na,S_T_T_R_CTF_M_q_local_omp__
     $           ,S_T_T_R_CTF_M_q_local_omp__)
         end if !if (flag_tesselation.eqv..false.) then
      end if                    ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      write(6,'(A,I0)') ' n_S_9_per_ ' , n_omp_sub__in
      call cp1_i4(n_omp_sub__in,n_S_9_per_,n_S_9_per_)
      write(6,'(A,I0)') ' n_S_9_sum_ ' , n_omp_sub__in
      call cp1_i4(n_omp_sub__in,n_S_9_sum_,n_S_9_sum_)
      write(6,'(A,I0)') ' n_M_9_per_ ' , n_omp_sub__in
      call cp1_i4(n_omp_sub__in,n_M_9_per_,n_M_9_per_)
      write(6,'(A,I0)') ' n_M_9_sum_ ' , n_omp_sub__in
      call cp1_i4(n_omp_sub__in,n_M_9_sum_,n_M_9_sum_)

      write(6,'(A,I0)') ' n_S_use_sum_ ' , n_omp_sub__in
      call cp1_i4(n_omp_sub__in,n_S_use_sum_,n_S_use_sum_)

      write(6,'(A)') '[finished ti7_check_variable_0]'

