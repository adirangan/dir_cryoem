!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calls ti8_Zstore_? ;\n
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling ti8_build_Zstore_full_0. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      if (flag_MS_vs_SM.eqv..true.) then
         write(6,'(A)') ' Warning, MS not implemented in local search. '
      end if !if (flag_MS_vs_SM.eqv..true.) then
      timing_tic = omp_get_wtime()

      call ti8_build_Zstore_full_0(
     $     verbose-1
     $     ,'0'
     $     ,rseed
     $     ,flag_RTRT_vs_RTTR
     $     ,n_delta_v
     $     ,delta_x_
     $     ,delta_y_
     $     ,n_gamma_z
     $     ,gamma_z_
     $     ,fftw_plan_frwd__(nM_9_sub*n_r)
     $     ,fftw_plan_back__(nM_9_sub*n_r)
     $     ,fftw_0in__(nM_9_sub*n_w_sum)
     $     ,fftw_out__(nM_9_sub*n_w_sum) 
     $     ,nufft1df77_fw__(nM_9_sub*n_gamma_z)
     $     ,nufft1df77_lw_(nM_9_sub)
     $     ,nufft1df77_lused_(nM_9_sub)
     $     ,n_r 
     $     ,grid_k_p_r_
     $     ,n_w_
     $     ,n_w_sum
     $     ,nS_1_per_max
     $     ,nS_1_sum
     $     ,0
     $     ,S_alpha_S_index_local_omp__(p_S_alpha_S_index_local)
     $     ,S_alpha_polar_a_local_omp__(p_S_alpha_polar_a_local)
     $     ,S_alpha_azimu_b_local_omp__(p_S_alpha_azimu_b_local)
     $     ,1
     $     ,polar_a_est_(nm0) 
     $     ,azimu_b_est_(nm0)
     $     ,gamma_z_est_(nm0) 
     $     ,delta_x_est_(nm0)
     $     ,delta_y_est_(nm0) 
     $     ,l2_norm_est_(nm0)
     $     ,ctf_ind_est_(nm0) 
     $     ,S_index_est_(nm0)
     $     ,M_index_est_(nm0)
     $     ,alpha_0in_omp__(nM_9_sub*n_alpha)
     $     ,alpha_update_f 
     $     ,displacement_max
     $     ,flag_MS_vs_SM ! Warning, Assumed .false. ;
     $     ,n_SM_max
     $     ,n_SM_(nm0) 
     $     ,alpha_SM__(n_alpha*n_SM_max*nm0)
     $     ,n_MS_max ! Warning, Not implemented ;
     $     ,n_MS_(0) ! Warning, Not implemented ;
     $     ,alpha_MS__(0) ! Warning, Not implemented ;
     $     ,C_M_(nm0) 
     $     ,n_CTF
     $     ,CTF_R_S_local_omp__(p_CTF_R_S_local)
     $     ,CTF_R_S_sub_omp__(nM_9_sub*n_gamma_z)
     $     ,1
     $     ,S_T_T_R_CTF_M_q_local_omp__(p_S_T_T_R_CTF_M_q_local)
     $     ,ZZ_omp__(nM_9_sub*n_delta_v*n_gamma_z)
     $     ,ZZ_sub_omp__(nM_9_sub*n_w_max)
     $     ,ZZ_pos_omp__(nM_9_sub*n_gamma_z)
     $     )

      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      gnump_tot = (n_w_max*1.0d0)*(n_transf*1.0d0)*(nS_1_sum*1.0d0)
      timing_total_Zstore = timing_total_Zstore + timing_tot
      gnump_total_Zstore = gnump_total_Zstore + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.4)') ' finished Zstore. ' , ' total time ' ,
     $        timing_tot
         write(6,'(A,F8.4)') ' Zstore total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if                    !if (verbose_timing.gt.0) then

