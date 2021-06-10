!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Uses omp to call each part of ti8_build_Zstore_full_0. ;\n

      if (flag_time_Zstore) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_a. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'a'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_a = timing_total_Zstore_a + timing_tot
      gnump_total_Zstore_a = gnump_total_Zstore_a + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_a. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_a total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_b. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'b'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_b = timing_total_Zstore_b + timing_tot
      gnump_total_Zstore_b = gnump_total_Zstore_b + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_b. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_b total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_c. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'c'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_c = timing_total_Zstore_c + timing_tot
      gnump_total_Zstore_c = gnump_total_Zstore_c + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_c. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_c total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_x. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'x'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_x = timing_total_Zstore_x + timing_tot
      gnump_total_Zstore_x = gnump_total_Zstore_x + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_x. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_x total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_d. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'d'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_d = timing_total_Zstore_d + timing_tot
      gnump_total_Zstore_d = gnump_total_Zstore_d + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_d. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_d total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_part_e. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
c$OMP&n_r_tmp,n_w_sum_tmp,n_X_tmp,
c$OMP&nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_r_tmp = nM_9_sub*n_r
         n_w_sum_tmp = nM_9_sub*n_w_sum
         n_X_tmp = nM_9_sub*n_delta_v*n_gamma_z
         nm0 = nM_9_sum
         nm1 = nM_1_sum + nm0
         nm2 = nM_0_sum + nm1
         nm3 = n_alpha*n_SM_max*nm2
         nm4 = n_delta_v*nS_1_per*nm0
         ns0 = 0
         ns1 = nS_1_sum + ns0
         ns2 = nS_0_sum + ns1
         ns3 = n_alpha*n_MS_max*ns2
         ns4 = n_gamma_z*n_CTF*ns1
         call ti8_build_Zstore_full_0(
     $        verbose-1
     $        ,'e'
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp)
     $        ,fftw_0in__(n_w_sum_tmp)
     $        ,fftw_out__(n_w_sum_tmp) 
     $        ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $        ,nufft1df77_lw_(nM_9_sub)
     $        ,nufft1df77_lused_(nM_9_sub)
     $        ,n_r 
     $        ,grid_k_p_r_
     $        ,n_w_
     $        ,n_w_sum
     $        ,nS_1_per
     $        ,nS_1_per
     $        ,ns2
     $        ,S_alpha_S_index_(ns2)
     $        ,S_alpha_polar_a_(ns2)
     $        ,S_alpha_azimu_b_(ns2)
     $        ,nM_9_per
     $        ,polar_a_est_(nm1) 
     $        ,azimu_b_est_(nm1)
     $        ,gamma_z_est_(nm1) 
     $        ,delta_x_est_(nm1)
     $        ,delta_y_est_(nm1) 
     $        ,l2_norm_est_(nm1)
     $        ,ctf_ind_est_(nm1) 
     $        ,S_index_est_(nm1)
     $        ,M_index_est_(nm1)
     $        ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_omp__(ns2 + nM_9_sub*n_S) 
     $        ,alpha_MS_omp__(ns3 + nM_9_sub*n_alpha*n_MS_max*n_S)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_X_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $        ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $        )
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_per*1.0d0)
     $     *(nM_1_per*1.0d0)
      timing_total_Zstore_e = timing_total_Zstore_e + timing_tot
      gnump_total_Zstore_e = gnump_total_Zstore_e + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore_e. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore_e total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose_timing.gt.0) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      end if ! if (flag_time_Zstore) then
