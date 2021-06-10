      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         call cl1_c16(nS_1_per*n_transf*nM_1_per*n_w_max
     $        ,S_Z_T_R_CTF_M_q__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      if (n_omp_sub__in.gt.1) then
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = n_r*(nS_1_sum + nS_0_per*nw)
         nx2 = n_r*n_transf*(nM_1_sum + nM_0_per*nw)
         nx3 = nw*nS_1_per*n_transf*nM_1_per
         call zgemm('T','N',nS_1_per,n_transf*nM_1_per,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),O_S_q__(nx1),n_r,Z_T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,nS_1_per)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_omp_sub__in.gt.1) then
      do nw=0,n_w_max-1
         nx1 = n_r*(nS_1_sum + nS_0_per*nw)
         nx2 = n_r*n_transf*(nM_1_sum + nM_0_per*nw)
         nx3 = nw*nS_1_per*n_transf*nM_1_per
         call zgemm('T','N',nS_1_per,n_transf*nM_1_per,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),O_S_q__(nx1),n_r,Z_T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,nS_1_per)
      enddo !do nw=0,n_w_max-1
      end if !if (n_omp_sub__in.gt.1) then
      n_SM_use = nS_1_per*nM_1_per
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_r*1.0d0)*(nS_1_per*1.0d0)*(n_transf*1.0d0)
     $     *(nM_1_per*1.0d0)*(n_w_max*1.0d0)
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
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 4
         MDA_d_(0) = nS_1_per
         MDA_d_(1) = n_transf
         MDA_d_(2) = nM_1_per
         MDA_d_(3) = n_w_max
         write(MDA_string,'(A,4(A,I0),A)')
     $        './dir_mda/S_x_Z_T_R_CTF_M_q__', 'nS0_', nS_0_sub ,
     $'_nS1_' , nS_1_sub , '_nM0_' ,nM_0_sub , '_nM1_' , nM_1_sub ,
     $'.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,S_Z_T_R_CTF_M_q__,MDA_string)
      end if                    !if (flag_MDA.eqv..true.) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      if (n_M_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,n_r_tmp,n_A_tmp,n_F_tmp)
c$OMP DO 
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         n_F_tmp = nM_9_sub*n_fpm
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp)
     $        ,fpm_back_(nM_9_sub) ,nS_1_per*n_transf*nM_9_per ,nS_1_per
     $        *n_transf*nM_1_per,S_Z_T_R_CTF_M_q__(nS_1_per*n_transf
     $        *nM_9_sum))
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_M_9_sub_use.gt.1) then
         nM_9_sub = 0
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         n_F_tmp = nM_9_sub*n_fpm
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp)
     $        ,fpm_back_(nM_9_sub) ,nS_1_per*n_transf*nM_9_per ,nS_1_per
     $        *n_transf*nM_1_per,S_Z_T_R_CTF_M_q__(nS_1_per*n_transf
     $        *nM_9_sum))
      end if !if (n_M_9_sub_use.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(nS_1_per*1.0d0)*(n_transf*1.0d0)
     $     *(nM_1_per*1.0d0)
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
      if (n_M_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,n_r_tmp,n_A_tmp,n_F_tmp)
c$OMP DO 
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         n_F_tmp = nM_9_sub*n_fpm
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp)
     $        ,fpm_back_(nM_9_sub) ,nS_1_per*n_transf*nM_9_per ,nS_1_per
     $        *n_transf*nM_1_per,S_Z_T_R_CTF_M_q__(nS_1_per*n_transf
     $        *nM_9_sum))
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_M_9_sub_use.gt.1) then
         nM_9_sub = 0
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_A_tmp = nM_9_sub*n_A
         n_F_tmp = nM_9_sub*n_fpm
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp)
     $        ,fpm_back_(nM_9_sub) ,nS_1_per*n_transf*nM_9_per ,nS_1_per
     $        *n_transf*nM_1_per,S_Z_T_R_CTF_M_q__(nS_1_per*n_transf
     $        *nM_9_sum))
      end if !if (n_M_9_sub_use.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(nS_1_per*1.0d0)*(n_transf*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_fpm_fftw = timing_total_fpm_fftw + timing_tot
      gnump_total_fpm_fftw = gnump_total_fpm_fftw + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
         write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' transposing dimensions (1,2) of S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      if (n_omp_sub__in.gt.1) then
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,nx0,nx1,nx2)
c$OMP DO 
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         nx2 = nM_9_sub*nS_1_per_max*n_transf
         do nx0=0,nM_9_per*n_w_max-1
            nx1 = (nM_9_sum*n_w_max + nx0)*(nS_1_per*n_transf)
            call trn0_1_c16(nS_1_per,n_transf,S_Z_T_R_CTF_M_q__(nx1)
     $        ,S_Z_T_R_CTF_M_q__(nx1),C_trn0_omp__(nx2))
         enddo !do nx0=0,nM_9_per*n_w_max-1
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_omp_sub__in.gt.1) then
      do nx0=0,nM_1_per*n_w_max-1
         nx1 = nx0*nS_1_per*n_transf
         call trn0_1_c16(nS_1_per,n_transf,S_Z_T_R_CTF_M_q__(nx1)
     $        ,S_Z_T_R_CTF_M_q__(nx1),C_trn0_omp__(0))
      enddo !do nx0=0,nM_1_per*n_w_max-1
      end if !if (n_omp_sub__in.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(nS_1_per*1.0d0)*(n_transf*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_transpose = timing_total_transpose + timing_tot
      gnump_total_transpose = gnump_total_transpose + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,F8.4)') ' finished transposing. ' , timing_tot
         write(6,'(A,F8.4)') ' transpose total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 4
         MDA_d_(0) = n_transf
         MDA_d_(1) = nS_1_per
         MDA_d_(2) = nM_1_per
         MDA_d_(3) = n_w_max
         write(MDA_string,'(A,4(A,I0),A)')
     $        './dir_mda/F_S_Z_T_R_CTF_M_q__', 'nS0_', nS_0_sub ,
     $'_nS1_' , nS_1_sub , '_nM0_' ,nM_0_sub , '_nM1_' , nM_1_sub ,
     $'.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,S_Z_T_R_CTF_M_q__,MDA_string)
      end if                    !if (flag_MDA.eqv..true.) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)') ' Now clearing array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_M_svdd_. '
         write(6,'(A)') ' '
      end if
      call cl1_c16(n_delta_v*nS_1_per*nM_1_per*n_w_max
     $     ,S_T_T_R_CTF_M_q__)
      timing_tic = omp_get_wtime()
      if (n_omp_sub__in.gt.1) then
c$OMP PARALLEL PRIVATE(n_A_tmp,n_X_tmp)
c$OMP DO 
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_v*nS_1_per*nM_1_per*nw
         n_X_tmp = n_transf*nS_1_per*nM_1_per*nw
         call zgemm('T','N',n_delta_v,nS_1_per*nM_1_per*1
     $        ,n_svd_l ,1.0d0*cmplx(1.0d0,0.0d0),Z_M_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q__(n_A_tmp),n_delta_v)
      enddo ! do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      else !if (n_omp_sub__in.gt.1) then
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_v*nS_1_per*nM_1_per*nw
         n_X_tmp = n_transf*nS_1_per*nM_1_per*nw
         call zgemm('T','N',n_delta_v,nS_1_per*nM_1_per*1
     $        ,n_svd_l ,1.0d0*cmplx(1.0d0,0.0d0),Z_M_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q__(n_A_tmp),n_delta_v)
      enddo ! do nw=0,n_w_max-1
      end if !if (n_omp_sub__in.gt.1) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_svd_l*1.0d0)*(n_delta_v*1.0d0)*(nS_1_per
     $     *1.0d0)*(nM_1_per*1.0d0)
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
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 4
         MDA_d_(0) = n_delta_v
         MDA_d_(1) = nS_1_per
         MDA_d_(2) = nM_1_per
         MDA_d_(3) = n_w_max
         write(MDA_string,'(A,4(A,I0),A)')
     $        './dir_mda/Z_M_svdd_F_S_Z_T_R_CTF_M_q__', 'nS0_', nS_0_sub
     $        , '_nS1_' , nS_1_sub , '_nM0_' ,nM_0_sub , '_nM1_' ,
     $        nM_1_sub , '.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,S_T_T_R_CTF_M_q__,MDA_string)
      end if                    !if (flag_MDA.eqv..true.) then
c$$$      %%%%%%%%%%%%%%%%
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling ti7_SxZTRM_check_4 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
c$$$         svd_r_max = 0.0d0
c$$$         svd_r_max = grid_k_p_(n_r-1)
c$$$         svd_d_max = 0.0d0
c$$$         do ndv=0,n_delta_v-1
c$$$            delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
c$$$            if (svd_d_max.lt.delta) then
c$$$               svd_d_max = delta
c$$$            end if              !if (svd_d_max.lt.delta) then
c$$$         enddo                  !do ndv=0,n_delta_v-1
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,F8.6)') ' setting svd_d_max: ' , svd_d_max
c$$$         end if                 !if (verbose.gt.0) then         
         call  ti7_check_SxZTRM_4(verbose-1,svd_d_max,n_svd_d,svd_d_
     $        ,svd_r_max,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_U_d_,svd_s_
     $        ,svd_V_r_ ,Z_M_svdd_,n_delta_v ,delta_x_,delta_y_
     $        ,n_gamma_z ,gamma_z_,delta_x_est_(nM_1_sum)
     $        ,delta_y_est_(nM_1_sum) ,gamma_z_est_(nM_1_sum)
     $        ,ctf_ind_est_(nM_1_sum) ,fftw_plan_frwd_ ,fftw_plan_back_
     $        ,fftw_in1_ ,fftw_out_ ,fftw_plan_back_last,fftw_in1_last_
     $        ,fftw_out_last_ ,n_r ,grid_k_p_,n_w_,n_A,S_p_,S_q_
     $        ,nS_1_per ,I_S_sample_(nS_0_sum + nS_1_sum),ld_S ,S_p__
     $        ,M_p_,M_q_ ,nM_1_per,I_M_sample_(nM_0_sum + nM_1_sum),ld_M
     $        ,M_p__ ,CTF_p_,n_CTF ,ld_CTF ,CTF_p__ ,C_M_(nM_0_sum +
     $        nM_1_sum) ,CTF_R_S__(n_gamma_z*n_CTF*(nS_0_sum +
     $        nS_1_sum)) ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then
