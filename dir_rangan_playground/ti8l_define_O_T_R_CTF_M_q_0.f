!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calls ti8_build_O_T_R_CTF_M_q_?. ;\n
      if (flag_RTRT_vs_RTTR.eqv..true.) then
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' Calculate bessel coefficients of T(R(CTF.*M)).'
         write(6,'(A)') ' T = translation by -delta_est. '
         write(6,'(A)') ' R = rotation by -gamma_est.'
         write(6,'(A)') ' CTF = dconjg(CTF) in k-space polar coord.'
         write(6,'(A)') ' M = template in k-space polar coord.'
         write(6,'(A)') ' '
         write(6,'(A,A)') ' O_T_R_CTF_M_q__(nr + n_r*(nw ' ,
     $        '+ n_w_max*nM_9_sub)) '
         write(6,'(A)') ' is equal to the bessel-coefficients'
         write(6,'(A)') ' M_q_(ic)) '
         write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
         write(6,'(A)') ' where M_q_ = T(R(CTF.*M)), with '
         write(6,'(A)') ' T <-- -delta_est'
         write(6,'(A)') ' R <-- -gamma_est'
         write(6,'(A)') ' and CTF <-- dconjg(CTF_(nctf))'
         write(6,'(A)') ' and M is derived from '
         write(6,'(A)') ' M_k_p__(I_M_sample_(nM_9_sum + nm)*ld_M).'
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      nm0 = nM_9_sum + nm
      n_r_tmp = nM_9_sub*n_r
      n_A_tmp = nM_9_sub*n_A
      p_O_T_R_CTF_M_q_local = n_r*n_w_max*nM_9_sub
      call cl1_c16(n_r*n_w_max
     $     ,O_T_R_CTF_M_q_local_omp__(p_O_T_R_CTF_M_q_local))
      if (verbose.gt.1) then
         write(6,'(4(A,I0))') ' nM_9_sub: ' , nM_9_sub , ' nm ' , nm ,
     $        ' nm0 ' , nm0, ' p_O_T_R_CTF_M_q_local: '
     $        ,p_O_T_R_CTF_M_q_local
      end if                    !if (verbose.gt.1) then
      call ti8_build_O_T_R_CTF_M_q_3(verbose-1 ,delta_x_est_(nm0)
     $     ,delta_y_est_(nm0) ,gamma_z_est_(nm0) ,ctf_ind_est_(nm0)
     $     ,fftw_plan_frwd__(n_r_tmp) ,fftw_plan_back__(n_r_tmp)
     $     ,fftw_0in__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_
     $     ,weight_k_p_r_,n_w_,n_A ,M_p_omp__(n_A_tmp)
     $     ,M_q_omp__(n_A_tmp),1 ,I_M_sample_(nm0) ,ld_M ,M_k_p__
     $     ,CTF_p_omp__(n_A_tmp) ,n_CTF ,ld_CTF ,CTF_k_p__ ,C_M_(nm0) ,1
     $     ,O_T_R_CTF_M_q_local_omp__(p_O_T_R_CTF_M_q_local))
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (n_A*1.0d0)
      timing_total_O_T_R_CTF_M_q = timing_total_O_T_R_CTF_M_q +
     $     timing_tot
      gnump_total_O_T_R_CTF_M_q = gnump_total_O_T_R_CTF_M_q + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.5)')
     $        ' finished calculating O_T_R_CTF_M_q__ for one M. ' ,
     $        ' total time ' , timing_tot
         write(6,'(A,F8.4)') ' O_T_R_CTF_M_q__ total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose_timing.gt.0) then
      end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
