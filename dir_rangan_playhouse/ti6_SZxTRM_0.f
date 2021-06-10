      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then

         call cl1_c16(n_transf*nS_1_per_max*n_w_max
     $        ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      do nw=0,n_w_max-1
         nx1 = n_r*n_transf*(0 + nS_1_per_max*nw)
         nx2 = n_r*(0 + 1*nw)
         nx3 = nw*n_transf*nS_1_per_max*1
         nx4 = n_r*n_transf*nS_1_per_max*n_w_max*n_omp_sub__in
         if (nx1+n_r*n_transf*nS_1_per_max+p_Z_S_q_local-1.ge.nx4) then
            write(6,'(A,A,I0,1X,I0,A,I0)')
     $           ' Warning: Z_S_q_local_omp__: ' ,' range: ' ,nx1 +
     $           p_Z_S_q_local , nx1 +p_Z_S_q_local +n_r*n_transf
     $           *nS_1_per_max-1 , ' max: ' , nx4-1
         end if !if (nx1+n_r*n_transf*nS_1_per_max+p_Z_S_q_local-1.ge.nx4) then
         nx4 = n_r*n_w_max*n_omp_sub__in
         if (nx2+n_r*1+p_O_T_R_CTF_M_q_local-1.ge.nx4) then
            write(6,'(A,A,I0,1X,I0,A,I0)')
     $           ' Warning: O_T_R_CTF_M_q_local_omp__: ' ,' range: '
     $           ,nx2 + p_O_T_R_CTF_M_q_local , nx2
     $           +p_O_T_R_CTF_M_q_local +n_r*1-1 ,
     $           ' max: ' , nx4-1
         end if !if (nx2+n_r*1+p_O_T_R_CTF_M_q_local-1.ge.nx4) then
         nx4 = n_transf*nS_1_per_max*n_w_max*n_omp_sub__in
         if (nx3+n_transf*nS_1_per_max+p_S_Z_T_R_CTF_M_q_local-1.ge.nx4)
     $        then
            write(6,'(A,A,I0,1X,I0,A,I0)')
     $           ' Warning: S_Z_T_R_CTF_M_q_local_omp__: ' ,' range: '
     $           ,nx3 + p_S_Z_T_R_CTF_M_q_local , nx3
     $           +p_S_Z_T_R_CTF_M_q_local +n_transf*nS_1_per_max-1 ,
     $           ' max: ' , nx4-1
         end if !if (nx3+n_r*1+p_S_Z_T_R_CTF_M_q_local-1.ge.nx4) then
         call zgemm('T','N',n_transf*nS_1_per_max,1,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_q_local_omp__(nx1+p_Z_S_q_local)
     $        ,n_r,O_T_R_CTF_M_q_local_omp__(nx2+p_O_T_R_CTF_M_q_local)
     $        ,n_r ,0.0d0 *cmplx(1.0d0,0.0d0)
     $        ,S_Z_T_R_CTF_M_q_local_omp__(nx3+p_S_Z_T_R_CTF_M_q_local)
     $        ,n_transf*nS_1_per_max)
      enddo !do nw=0,n_w_max-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_r*1.0d0)*(n_transf*1.0d0)*(nS_1_per_max*1.0d0)
     $     *(n_w_max*1.0d0)
      timing_total_zgemm = timing_total_zgemm + timing_tot
      gnump_total_zgemm = gnump_total_zgemm + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total time: ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total Gflops: ' ,
     $        timing_tmp
      end if !if (verbose_timing.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      n_r_tmp = nM_9_sub*n_r
      n_A_tmp = nM_9_sub*n_A
      n_F_tmp = nM_9_sub*n_fpm
      call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $     ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_9_sub)
     $     ,n_transf*nS_1_per_max ,n_transf *nS_1_per_max
     $     ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local))
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per_max*1.0d0)
      timing_total_fpm_fill = timing_total_fpm_fill + timing_tot
      gnump_total_fpm_fill = gnump_total_fpm_fill + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
         write(6,'(A,F8.4)') ' array fill total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      end if ! if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      n_r_tmp = nM_9_sub*n_r
      n_A_tmp = nM_9_sub*n_A
      n_F_tmp = nM_9_sub*n_fpm
      call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $     ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_9_sub)
     $     ,n_transf*nS_1_per_max ,n_transf *nS_1_per_max
     $     ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local))
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per_max*1.0d0)
      timing_total_fpm_fftw = timing_total_fpm_fftw + timing_tot
      gnump_total_fpm_fftw = gnump_total_fpm_fftw + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
         write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)') ' Now clearing array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_S_svdd_. '
         write(6,'(A)') ' '
      end if
      call cl1_c16(n_delta_v*nS_1_per_max*n_w_max
     $     ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      timing_tic = omp_get_wtime()
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_v*nS_1_per_max*nw
         n_X_tmp = n_transf*nS_1_per_max*nw
         call zgemm('T','N',n_delta_v,nS_1_per_max*1 ,n_svd_l ,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q_local_omp__(n_X_tmp
     $        +p_S_Z_T_R_CTF_M_q_local) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q_local_omp__(n_A_tmp
     $        +p_S_T_T_R_CTF_M_q_local),n_delta_v)
      enddo ! do nw=0,n_w_max-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_svd_l*1.0d0)*(n_delta_v*1.0d0)*(nS_1_per_max
     $     *1.0d0)
      timing_total_zgemm = timing_total_zgemm + timing_tot
      gnump_total_zgemm = gnump_total_zgemm + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_SZxTRM_3 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call test_innerproduct_5_SZxTRM_3(verbose-1,n_svd_r,svd_r_
     $        ,n_svd_l,svd_l_,svd_s_,svd_V_r_,Z_S_svdd_,n_delta_v
     $        ,delta_x_ ,delta_y_,n_gamma_z,gamma_z_ ,delta_x_est_(nm0)
     $        ,delta_y_est_(nm0) ,gamma_z_est_(nm0),ctf_ind_est_(nm0)
     $        ,fftw_plan_frwd_,fftw_plan_back_ ,fftw_in1_ ,fftw_out_
     $        ,fftw_plan_back_last ,fftw_in1_last_ ,fftw_out_last_ ,n_r
     $        ,grid_k_p_,n_w_,n_A ,S_p_ ,S_q_ ,nS_1_per_max
     $        ,I_S_sample_local_omp__(p_I_S_sample_local) ,ld_S,S_p__
     $        ,M_p_,M_q_ ,1 ,I_M_sample_(nm0) ,ld_M ,M_p__,CTF_p_ ,n_CTF
     $        ,ld_CTF ,CTF_p__ ,C_M_(nm0)
     $        ,CTF_R_S_local_omp__(p_CTF_R_S_local)
     $        ,S_Z_T_R_CTF_M_q_local_omp__(p_S_Z_T_R_CTF_M_q_local)
     $        ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local))
      end if !testing S_T_T_R_CTF_M_q__
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
