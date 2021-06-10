!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Uses omp to call blocks of ti8_build_Z_S_q_?. ;\n
      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.)) then
         call cl1_c16(n_r*n_transf*nS_0_per*n_w_max,Z_S_q__)
         if (verbose.gt.1) then
            write(6,'(A)') ' Calculate bessel coefficients of Z(S).'
            write(6,'(A)') ' Z = svd of translation operator. '
            write(6,'(A)') ' S = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A)') ' Z_S_q__(nr + n_r*(nl + '
     $           ,'n_svd_l*(ns + nS_0_per*nw))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' dconjg(S_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
            write(6,'(A)') ' where S_q_ = Z(S), with '
            write(6,'(A)') ' Z representing the rho-factor '
            write(6,'(A)') ' of the translation-operator'
            write(6,'(A)') ' associated with delta,'
            write(6,'(A)') ' and S is derived from '
            write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S).'
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
c$$$         svd_r_max = 0.0d0
c$$$         svd_r_max = grid_k_p_r_(n_r-1)
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,F8.2)') ' setting svd_r_max: ' , svd_r_max
c$$$         end if                 !if (verbose.gt.0) then
         timing_tic = omp_get_wtime()
         if (n_S_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,n_r_tmp,n_w_sum_tmp)
c$OMP DO 
         do nS_9_sub=0,n_S_9_sub_use-1
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            n_r_tmp = nS_9_sub*n_r
            n_w_sum_tmp = nS_9_sub*n_w_sum
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_build_Z_S_q_quad_5'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_Z_S_q_quad_5(
     $           verbose-1
     $           ,svd_r_max
     $           ,n_svd_r 
     $           ,svd_r_
     $           ,n_svd_l
     $           ,svd_l_
     $           ,svd_s_
     $           ,svd_V_r_
     $           ,svd_polyval_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) 
     $           ,fftw_0in__(n_w_sum_tmp)
     $           ,fftw_out__(n_w_sum_tmp)
     $           ,n_r 
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum 
     $           ,S_p_omp__(n_w_sum_tmp)
     $           ,S_q_omp__(n_w_sum_tmp) 
     $           ,nS_9_per
     $           ,I_S_sample_(nS_0_sum + ns_9_sum)
     $           ,ld_S 
     $           ,S_k_p__
     $           ,nS_0_per 
     $           ,Z_S_q__(n_r*n_svd_l*nS_9_sum)
     $           )
         enddo !do nS_sub=0,n_S_sub_use-1
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
     $              ,' ti8_build_Z_S_q_quad_5'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_Z_S_q_quad_5(
     $           verbose-1
     $           ,svd_r_max
     $           ,n_svd_r 
     $           ,svd_r_
     $           ,n_svd_l
     $           ,svd_l_
     $           ,svd_s_
     $           ,svd_V_r_
     $           ,svd_polyval_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) 
     $           ,fftw_0in__(n_w_sum_tmp)
     $           ,fftw_out__(n_w_sum_tmp)
     $           ,n_r 
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum 
     $           ,S_p_omp__(n_w_sum_tmp)
     $           ,S_q_omp__(n_w_sum_tmp) 
     $           ,nS_9_per
     $           ,I_S_sample_(nS_0_sum + ns_9_sum)
     $           ,ld_S 
     $           ,S_k_p__
     $           ,nS_0_per 
     $           ,Z_S_q__(n_r*n_svd_l*nS_9_sum)
     $           )
         end if !if (n_S_9_sub_use.gt.1) then
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         gnump_tot = (n_w_sum*1.0d0)*(n_svd_l*1.0d0)*(nS_0_per*1.0d0)
         timing_total_Z_S_q = timing_total_Z_S_q + timing_tot
         gnump_total_Z_S_q = gnump_total_Z_S_q + gnump_tot
         timing_tmp = gnump_tot/timing_tot/1e9
         if (verbose_timing.gt.0) then
            write(6,'(A,A,F8.5)')
     $           ' finished calculating Z_S_q__ for each S. ' ,
     $           ' total time ' , timing_tot
            write(6,'(A,F8.4)') ' Z_S_q__ total Gnumps: ' , timing_tmp
            write(6,'(A)') ' '
         end if ! if (verbose_timing.gt.0) then
         flag_MDA = .false.
         if (flag_MDA.eqv..true.) then
            MDA_n_d = 4
            MDA_d_(0) = n_r
            MDA_d_(1) = n_transf
            MDA_d_(2) = nS_0_per
            MDA_d_(3) = n_w_max
            write(MDA_string,'(A,1(A,I0),A)') './dir_mda/Z_S_q__' ,
     $           'nS0_' , nS_0_sub,'.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_q__,MDA_string)
         end if !if (flag_MDA.eqv..true.) then
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then
