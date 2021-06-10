!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Uses omp to call blocks of ti8_build_O_T_R_CTF_M_q_?. ;\n
      if (flag_RTRT_vs_RTTR.eqv..true.) then
      call cl1_c16(n_r*nM_0_per*n_w_max,O_T_R_CTF_M_q__)
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' Calculate bessel coefficients of T(R(CTF.*M)).'
         write(6,'(A)') ' T = translation by -delta_est. '
         write(6,'(A)') ' R = rotation by -gamma_est.'
         write(6,'(A)') ' CTF = dconjg(CTF) in k-space polar coord.'
         write(6,'(A)') ' M = template in k-space polar coord.'
         write(6,'(A)') ' '
         write(6,'(A,A)') ' O_T_R_CTF_M_q__(nr ' ,
     $        '+ n_r*(nm + nM_0_per*nw)) '
         write(6,'(A)') ' is equal to the bessel-coefficients'
         write(6,'(A)') ' M_q_(ic)) '
         write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
         write(6,'(A)') ' where M_q_ = T(R(CTF.*M)), with '
         write(6,'(A)') ' T <-- -delta_est'
         write(6,'(A)') ' R <-- -gamma_est'
         write(6,'(A)') ' and CTF <-- dconjg(CTF_(nctf))'
         write(6,'(A)') ' and M is derived from '
         write(6,'(A)') ' M_k_p__(I_M_sample_(nM_0_sum + nm)*ld_M).'
         write(6,'(A)') ' '
      end if
      
      timing_tic = omp_get_wtime()
      if (n_M_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,n_r_tmp,n_A_tmp)
c$OMP DO 
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         call ti8_build_O_T_R_CTF_M_q_3(verbose-1
     $        ,delta_x_est_(nM_9_sum),delta_y_est_(nM_9_sum)
     $        ,gamma_z_est_(nM_9_sum),ctf_ind_est_(nM_9_sum)
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_
     $        ,weight_k_p_r_ ,n_w_,n_A ,M_p_omp__(n_A_tmp)
     $        ,M_q_omp__(n_A_tmp),nM_9_per ,I_M_sample_(nM_0_sum +
     $        nM_9_sum),ld_M ,M_k_p__ ,CTF_p_omp__(n_A_tmp) ,n_CTF
     $        ,ld_CTF ,CTF_k_p__ ,C_M_(nM_0_sum + nM_9_sum) ,nM_0_per
     $        ,O_T_R_CTF_M_q__(n_r *nM_9_sum))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_M_9_sub_use.gt.1) then
         nM_9_sub = 0
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         call ti8_build_O_T_R_CTF_M_q_3(verbose-1
     $        ,delta_x_est_(nM_9_sum),delta_y_est_(nM_9_sum)
     $        ,gamma_z_est_(nM_9_sum),ctf_ind_est_(nM_9_sum)
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_
     $        ,weight_k_p_r_ ,n_w_,n_A ,M_p_omp__(n_A_tmp)
     $        ,M_q_omp__(n_A_tmp),nM_9_per ,I_M_sample_(nM_0_sum +
     $        nM_9_sum),ld_M ,M_k_p__ ,CTF_p_omp__(n_A_tmp) ,n_CTF
     $        ,ld_CTF ,CTF_k_p__ ,C_M_(nM_0_sum + nM_9_sum) ,nM_0_per
     $        ,O_T_R_CTF_M_q__(n_r *nM_9_sum))
      end if !if (n_M_9_sub_use.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (n_A*1.0d0)*(nM_0_per*1.0d0)
      timing_total_O_T_R_CTF_M_q = timing_total_O_T_R_CTF_M_q +
     $     timing_tot
      gnump_total_O_T_R_CTF_M_q = gnump_total_O_T_R_CTF_M_q + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.5)')
     $        ' finished calculating O_T_R_CTF_M_q__ for each M. ' ,
     $        ' total time ' , timing_tot
         write(6,'(A,F8.4)') ' O_T_R_CTF_M_q__ total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose_timing.gt.0) then
         flag_MDA = .false.
         if (flag_MDA.eqv..true.) then
            MDA_n_d = 3
            MDA_d_(0) = n_r
            MDA_d_(1) = nM_0_per
            MDA_d_(2) = n_w_max
            write(MDA_string,'(A,2(A,I0),A)')
     $           './dir_mda/O_T_R_CTF_M_q__', 'nS0_' , nS_0_sub ,
     $           '_nM0_' , nM_0_sub ,'.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,O_T_R_CTF_M_q__
     $           ,MDA_string)
         end if !if (flag_MDA.eqv..true.) then
      end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
