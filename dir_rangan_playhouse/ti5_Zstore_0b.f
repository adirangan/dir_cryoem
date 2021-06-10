      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling test_innerproduct_5_Zstore_3b. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,n_A_tmp,nm0,nm1,nm2,nm3,nm4,
c$OMP&ns0,ns1,ns2,ns3,ns4)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         n_A_tmp = nM_9_sub*n_delta_v*n_gamma_z
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
         call test_innerproduct_5_Zstore_3b(
     $        verbose-1
     $        ,rseed
     $        ,flag_RTRT_vs_RTTR
     $        ,n_delta_v
     $        ,delta_x_
     $        ,delta_y_
     $        ,n_gamma_z
     $        ,gamma_z_
     $        ,n_r 
     $        ,grid_k_p_
     $        ,n_w_
     $        ,n_A
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
     $        ,alpha__in_omp__(nM_9_sub*n_alpha)
     $        ,alpha_update_f 
     $        ,displacement_max
     $        ,flag_MS_vs_SM
     $        ,n_SM_max
     $        ,n_SM_(nm2) 
     $        ,alpha_SM__(nm3)
     $        ,n_MS_max
     $        ,n_MS_(ns2) 
     $        ,alpha_MS__(ns3)
     $        ,C_M_(nm2) 
     $        ,n_CTF
     $        ,CTF_R_S__(ns4)
     $        ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $        ,nM_1_per
     $        ,S_T_T_R_CTF_M_q__(nm4)
     $        ,ZZ_omp__(n_A_tmp)
     $        ,ZZ_sub_omp__(nM_9_sub*n_w_max)
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

