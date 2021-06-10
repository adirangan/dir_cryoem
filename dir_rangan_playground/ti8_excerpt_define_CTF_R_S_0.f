!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Uses omp to call blocks of ti8_build_CTF_R_S_?. ;\n
      if (verbose.gt.1) then
         write(6,'(A)') ' Clearing array CTF_R_S_ to hold'
         write(6,'(A)') ' innerproducts for each CTF-S pair.'
         write(6,'(A)') ' '
      end if
      call cl1_c16(n_gamma_z*n_CTF*nS_0_per,CTF_R_S__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate CTF_R_S_ for each ctf-S pair.'
         write(6,'(A)') ' CTF_R_S = || CTF .* R(S) ||.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)')
     $        ' CTF_R_S_(ngz + n_gamma_z*(nctf + n_CTF*ns)) '
         write(6,'(A)') ' is equal to the l2_norm (not squared)'
         write(6,'(A)') ' of the pointwise product of CTF and R(S),'
         write(6,'(A)') ' where CTF is the nctf-th ctf-function, '
         write(6,'(A)') ' R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' and S is the ns-th template, '
         write(6,'(A)') ' at position '
         write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S). '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      if (n_S_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,
c$OMP& n_r_tmp,n_w_sum_tmp,n_X_tmp)
c$OMP DO 
      do nS_9_sub=0,n_S_9_sub_use-1
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_w_sum_tmp = nS_9_sub*n_w_sum
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_CTF_R_S_3'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_CTF_R_S_quad_3(
     $        verbose-1
     $        ,n_gamma_z 
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp) 
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp) 
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,weight_k_p_r_
     $        ,n_w_
     $        ,n_w_sum 
     $        ,S_p_omp__(n_w_sum_tmp)
     $        ,S_q_omp__(n_w_sum_tmp)
     $        ,nS_9_per 
     $        ,I_S_sample_(nS_0_sum + nS_9_sum)
     $        ,ld_S 
     $        ,S_k_p__ 
     $        ,CTF_p_omp__(n_w_sum_tmp)
     $        ,CTF_q_omp__(n_w_sum_tmp) 
     $        ,n_CTF 
     $        ,ld_CTF 
     $        ,CTF_k_p__
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF) 
     $        ,Z_p_omp__(nS_9_sub*n_Z_x)
     $        ,Z_q_omp__(nS_9_sub*n_Z_x)
     $        )
         call af1_c16(
     $        n_gamma_z*nS_9_per*n_CTF
     $        ,dcmplx(1.0d0/dsqrt(2.0d0),0.0d0)
     $        ,dcmplx(0.0d0,0.0d0)
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        )
      enddo                     !do nS_sub=0,n_S_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_S_9_sub_use.gt.1) then
         nS_9_sub = 0
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_w_sum_tmp = nS_9_sub*n_w_sum
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_CTF_R_S_3'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_CTF_R_S_quad_3(
     $        verbose-1
     $        ,n_gamma_z 
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp) 
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp) 
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,weight_k_p_r_
     $        ,n_w_
     $        ,n_w_sum 
     $        ,S_p_omp__(n_w_sum_tmp)
     $        ,S_q_omp__(n_w_sum_tmp)
     $        ,nS_9_per 
     $        ,I_S_sample_(nS_0_sum + nS_9_sum)
     $        ,ld_S 
     $        ,S_k_p__ 
     $        ,CTF_p_omp__(n_w_sum_tmp)
     $        ,CTF_q_omp__(n_w_sum_tmp) 
     $        ,n_CTF 
     $        ,ld_CTF 
     $        ,CTF_k_p__
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF) 
     $        ,Z_p_omp__(nS_9_sub*n_Z_x)
     $        ,Z_q_omp__(nS_9_sub*n_Z_x)
     $        )
         call af1_c16(
     $        n_gamma_z*nS_9_per*n_CTF
     $        ,dcmplx(1.0d0/dsqrt(2.0d0),0.0d0)
     $        ,dcmplx(0.0d0,0.0d0)
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        )
      end if !if (n_S_9_sub_use.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (n_w_sum*1.0d0)*(n_CTF*1.0d0)*(nS_0_per*1.0d0)
      timing_total_CTF_R_S = timing_total_CTF_R_S + timing_tot
      gnump_total_CTF_R_S = gnump_total_CTF_R_S + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.5)')
     $        ' finished calculating CTF_R_S_ for each ctf-S pair. ' ,
     $        ' total time ' , timing_tot
         write(6,'(A,F8.4)') ' CTF_R_S_ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if                    ! if (verbose_timing.gt.0) then
