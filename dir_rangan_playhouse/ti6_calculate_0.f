      if (verbose.gt.2) then
         call write_all_i4(n_M,n_SM_,12,' n_SM_pre_: ')
      end if !if (verbose.gt.2) then

      nS_1_sum = 0
      do nLT=0,n_LT-1
         ns_local = LT_omp__(nLT+p_LT)
         if (0+verbose.gt.1) then
            write(6,'(4(A,I0))') ' nm0: ' , nm0 , ' nLT: ' , nLT ,
     $           ' ns_local: ' ,ns_local , ' nS_1_sum: ' , nS_1_sum
         end if                 !if (0+verbose.gt.1) then
         S_alpha_S_index_local_omp__(nS_1_sum + p_S_alpha_S_index_local)
     $        = S_alpha_S_index_(ns_local + nS_0_sum)
         S_alpha_polar_a_local_omp__(nS_1_sum + p_S_alpha_polar_a_local)
     $        = S_alpha_polar_a_(ns_local + nS_0_sum)
         S_alpha_azimu_b_local_omp__(nS_1_sum + p_S_alpha_azimu_b_local)
     $        = S_alpha_azimu_b_(ns_local + nS_0_sum)
         I_S_sample_local_omp__(nS_1_sum + p_I_S_sample_local) =
     $        I_S_sample_(ns_local)
         nx0 = n_gamma_z*n_CTF*ns_local
         nx1 = n_gamma_z*n_CTF*nS_0_per_max
         nx2 = n_gamma_z*n_CTF*(nS_1_sum)
         nx3 = n_gamma_z*n_CTF*nS_1_per_max*n_omp_sub__in
         if (nx0+n_gamma_z*n_CTF-1.ge.nx1) then
            write(6,'(A,A,I0,1X,I0,A,I0)') ' Warning: CTF_R_S__: ' ,
     $           ' range: ' ,nx0, nx0+ n_gamma_z*n_CTF-1 , ' max: ' ,
     $           nx1-1
         end if                 !if (nx0+n_gamma_z*n_CTF-1.ge.nx1) then
         if (nx2+n_gamma_z*n_CTF+p_CTF_R_S_local-1.ge.nx3) then
            write(6,'(A,A,I0,1X,I0,A,I0)')
     $           ' Warning: CTF_R_S_local_omp__: ' ,' range: ' ,nx2 +
     $           p_CTF_R_S_local , nx2 +p_CTF_R_S_local +n_gamma_z*n_CTF
     $           -1 , ' max: ' , nx3-1               
         end if                 !if (nx2+n_gamma_z*n_CTF+p_CTF_R_S_local-1.ge.nx3) then
         if (verbose.gt.1) then
            write(6,'(A,A,I0,1X,I0,A,I0)') ' CTF_R_S__: ' , ' range: '
     $           ,nx0, nx0+ n_gamma_z*n_CTF-1 , ' max: ' , nx1-1
            write(6,'(A,A,I0,1X,I0,A,I0)') ' CTF_R_S_local_omp__: ' ,
     $           ' range: ' ,nx2 + p_CTF_R_S_local , nx2 +
     $           p_CTF_R_S_local +n_gamma_z*n_CTF-1 , ' max: ' , nx3-1
         end if                 !if (verbose.gt.1) then
         call cp1_c16(n_gamma_z*n_CTF,CTF_R_S__(nx0)
     $        ,CTF_R_S_local_omp__(nx2 + p_CTF_R_S_local))
         do nw=0,n_w_max-1
            if (flag_RTRT_vs_RTTR.eqv..false.) then
               nx0 = n_r*(ns_local + nS_0_per*nw)
               nx1 = n_r*nS_0_per_max*n_w_max
               nx2 = n_r*(nS_1_sum + nS_1_per_max*nw)
               nx3 = n_r*nS_1_per_max*n_w_max *n_omp_sub__in
               if (nx0+n_r-1.ge.nx1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' Warning: O_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r-1, ' max: ' , n_r
     $                 *nS_0_per_max*n_w_max
               end if !if (nx0+n_r-1.ge.nx1) then
               if (nx2+n_r+p_O_S_q_local-1.ge.nx3) then
                  write(6,'(A,A,I0,1X,I0,A,I0)')
     $                 ' Warning: O_S_q_local_omp__: ',' range: ' ,nx2 +
     $                 p_O_S_q_local , nx2+p_O_S_q_local + n_r-1,
     $                 ' max: ' , n_r*nS_1_per_max*n_w_max
     $                 *n_omp_sub__in
               end if !if (nx2+n_r+p_O_S_q_local-1.ge.nx3) then
               if (verbose.gt.1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' O_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r-1, ' max: ' , n_r
     $                 *nS_0_per_max*n_w_max
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' O_S_q_local_omp__: '
     $                 ,' range: ' ,nx2 + p_O_S_q_local , nx2
     $                 +p_O_S_q_local + n_r-1, ' max: ' , n_r
     $                 *nS_1_per_max*n_w_max *n_omp_sub__in
               end if           !if (verbose.gt.1) then
               call cp1_c16(n_r,O_S_q__(nx0),O_S_q_local_omp__(nx2 +
     $              p_O_S_q_local))
            end if              !if (flag_RTRT_vs_RTTR.eqv..false.) then
            if ((svd_calculation_type.eq.2) .and.
     $           (flag_RTRT_vs_RTTR.eqv..true.)) then
               nx0 = n_r*n_transf*(ns_local + nS_0_per*nw)
               nx1 = n_r*n_transf*nS_0_per_max*n_w_max
               nx2 = n_r*n_transf*(nS_1_sum + nS_1_per_max*nw)
               nx3 = n_r*n_transf*nS_1_per_max *n_w_max*n_omp_sub__in
               if (nx0+n_r*n_transf-1.ge.nx1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' Warning: T_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r*n_transf-1, ' max: ' ,
     $                 n_r*n_transf*nS_0_per_max*n_w_max
               end if !if (nx0+n_r*n_transf-1.ge.nx1) then
               if (nx2+n_r*n_transf+p_T_S_q_local-1.ge.nx3) then
                  write(6,'(A,A,I0,1X,I0,A,I0)')
     $                 ' Warning: T_S_q_local_omp__: ',' range: ' ,nx2 +
     $                 p_T_S_q_local , nx2+p_T_S_q_local + n_r*n_transf
     $                 -1, ' max: ' , n_r*n_transf*nS_1_per_max *n_w_max
     $                 *n_omp_sub__in
               end if !if (nx2+n_r*n_transf+p_T_S_q_local-1.ge.nx3) then
               if (verbose.gt.1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' T_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r*n_transf-1, ' max: ' ,
     $                 n_r*n_transf*nS_0_per_max*n_w_max
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' T_S_q_local_omp__: '
     $                 ,' range: ' ,nx2 + p_T_S_q_local , nx2
     $                 +p_T_S_q_local + n_r*n_transf-1, ' max: ' , n_r
     $                 *n_transf*nS_1_per_max *n_w_max*n_omp_sub__in
               end if           !if (verbose.gt.1) then
               call cp1_c16(n_r*n_transf,T_S_q__(nx0)
     $              ,T_S_q_local_omp__(nx2 + p_T_S_q_local))
            end if              !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
            if ((svd_calculation_type.eq.1) .and.
     $           (flag_RTRT_vs_RTTR.eqv..true.)) then
               nx0 = n_r*n_transf*(ns_local + nS_0_per*nw)
               nx1 = n_r*n_transf*nS_0_per_max*n_w_max
               nx2 = n_r*n_transf*(nS_1_sum + nS_1_per_max*nw)
               nx3 = n_r*n_transf*nS_1_per_max *n_w_max*n_omp_sub__in
               if (nx0+n_r*n_transf-1.ge.nx1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' Warning: Z_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r*n_transf-1, ' max: ' ,
     $                 n_r*n_transf*nS_0_per_max*n_w_max
               end if !if (nx0+n_r*n_transf-1.ge.nx1) then
               if (nx2+n_r*n_transf+p_Z_S_q_local-1.ge.nx3) then
                  write(6,'(A,A,I0,1X,I0,A,I0)')
     $                 ' Warning: Z_S_q_local_omp__: ',' range: ' ,nx2 +
     $                 p_Z_S_q_local , nx2+p_Z_S_q_local + n_r*n_transf
     $                 -1, ' max: ' , n_r*n_transf*nS_1_per_max *n_w_max
     $                 *n_omp_sub__in
               end if !if (nx2+n_r*n_transf+p_Z_S_q_local-1.ge.nx3) then
               if (verbose.gt.1) then
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' Z_S_q__: ' ,
     $                 ' range: ' ,nx0,nx0+ n_r*n_transf-1, ' max: ' ,
     $                 n_r*n_transf*nS_0_per_max*n_w_max
                  write(6,'(A,A,I0,1X,I0,A,I0)') ' Z_S_q_local_omp__: '
     $                 ,' range: ' ,nx2 + p_Z_S_q_local , nx2
     $                 +p_Z_S_q_local + n_r*n_transf-1, ' max: ' , n_r
     $                 *n_transf*nS_1_per_max *n_w_max*n_omp_sub__in
               end if           !if (verbose.gt.1) then
               call cp1_c16(n_r*n_transf,Z_S_q__(nx0)
     $              ,Z_S_q_local_omp__(nx2 + p_Z_S_q_local))
            end if              !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
         enddo !do nw=0,n_w_max-1
         nS_1_sum = nS_1_sum+1
         if ((nS_1_sum.eq.nS_1_per_max) .or. (nLT.eq.n_LT-1)) then
            n_SM_use_local_(nM_9_sub) = n_SM_use_local_(nM_9_sub) +
     $           nS_1_sum
            include 'ti6_STxTRM_0.f'
            include 'ti6_SxTTRM_0.f'
            include 'ti6_SZxTRM_0.f'
            include 'ti6_SxZTRM_0.f'
            include 'ti6_Zstore_1.f'
            if (nS_1_sum.eq.nS_1_per_max) then
               if (verbose.gt.1) then
                  write(6,'(2(A,I0))') ' ns_1_sum: ' , nS_1_sum ,
     $                 ' nS_1_per_max ' , nS_1_per_max
               end if           !if (verbose.gt.1) then
               call cl1_r8(nS_1_per_max ,S_alpha_S_index_local_omp__(
     $              p_S_alpha_S_index_local)) 
               call cl1_r8(nS_1_per_max ,S_alpha_polar_a_local_omp__(
     $              p_S_alpha_polar_a_local)) 
               call cl1_r8(nS_1_per_max ,S_alpha_azimu_b_local_omp__(
     $              p_S_alpha_azimu_b_local)) 
               call cl1_i4(nS_1_per_max
     $              ,I_S_sample_local_omp__(p_I_S_sample_local))
               call cl1_c16(n_gamma_z*n_CTF*nS_1_per_max
     $              ,CTF_R_S_local_omp__(p_CTF_R_S_local))
               if (flag_RTRT_vs_RTTR.eqv..false.) then
                  call cl1_c16(n_r*nS_1_per_max*n_w_max
     $                 ,O_S_q_local_omp__(p_O_S_q_local))
               end if           !if (flag_RTRT_vs_RTTR.eqv..false.) then
               if ((svd_calculation_type.eq.2) .and.
     $              (flag_RTRT_vs_RTTR.eqv..true.)) then
                  call cl1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $                 ,T_S_q_local_omp__(p_T_S_q_local))
               end if           !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
               if ((svd_calculation_type.eq.1) .and.
     $              (flag_RTRT_vs_RTTR.eqv..true.)) then
                  call cl1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $                 ,Z_S_q_local_omp__(p_Z_S_q_local))
               end if           !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
               nS_1_sum = 0
            end if !if (nS_1_sum.eq.nS_1_per_max) then
         end if !if ((nS_1_sum.eq.nS_1_per_max) .or. (nLT.eq.n_LT-1)) then
      enddo !do nLT=0,n_LT-1

      if (verbose.gt.2) then
         call write_all_i4(n_M,n_SM_,12,' n_SM_pos_: ')
      end if !if (verbose.gt.2) then

      if (flag_MS_vs_SM.eqv..false.) then
         if (verbose.gt.1) then
            write(6,'(2(A,I0))') ' nm0: ' , nm0 , ' n_SM: ' , n_SM_(nm0)
            call alpha_SM_write_0(n_SM_max,n_SM_(nm0),alpha_SM__(n_alpha
     $           *n_SM_max*nm0),16,' alpha_SM_pre_: ')
         end if !if (verbose.gt.-1) then
         call alpha_SM_sort_0(n_SM_max,n_SM_(nm0),alpha_SM__(n_alpha
     $        *n_SM_max*nm0))
         if (verbose.gt.1) then
            call alpha_SM_write_0(n_SM_max,n_SM_(nm0),alpha_SM__(n_alpha
     $           *n_SM_max*nm0),16,' alpha_SM_pos_: ')
         end if !if (verbose.gt.-1) then
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         write(6,'(A,A)') ' Warning: flag_MS_vs_SM.eqv..true. ' ,
     $        ' not implemented in test_innerproduct_5_local_0.f'
      end if !if (flag_MS_vs_SM.eqv..true.) then
