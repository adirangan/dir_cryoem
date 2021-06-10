      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[entering test_innerproduct_5]: '
      end if

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ',verbose
         write(6,'(A,L1)') ' flag_memory_estimate: '
     $        ,flag_memory_estimate
         write(6,'(A,I0)') ' rseed: ',rseed
         write(6,'(A,I0)') ' n_k_cur: ',n_k_cur
         call write_all_i4(n_k_cur,n_polar_a_,13,' n_polar_a_: ')
         call write_all_r8(n_k_cur,grid_k_p_,12,' grid_k_p_: ')
         write(6,'(A,F6.3)') ' half_diameter_x_c: ',half_diameter_x_c
         write(6,'(A,I0)') ' n_S: ',n_S
c$$$         call write_sub_i4(n_S,I_S_sample_,14,' I_S_sample_: ')
         write(6,'(A,I0)') ' ld_S: ',ld_S
c$$$         call write_sub_c16(n_S*ld_S,S_k_p__,11,' S_k_p___: ')
         write(6,'(A,F6.3)') ' tesselation_distance_req: '
     $        ,tesselation_distance_req
c$$$         call write_sub_r8(n_S,S_alpha_S_index_,19
c$$$     $        ,' S_alpha_S_index_: ')
c$$$         call write_sub_r8(n_S,S_alpha_polar_a_,19
c$$$     $        ,' S_alpha_polar_a_: ')
c$$$         call write_sub_r8(n_S,S_alpha_azimu_b_,19
c$$$     $        ,' S_alpha_azimu_b_: ')
         write(6,'(A,I0)') ' n_M: ',n_M
c$$$         call write_sub_i4(n_M,I_M_sample_,14,' I_M_sample_: ')
         write(6,'(A,I0)') ' ld_M: ',ld_M
c$$$         call write_sub_c16(n_M*ld_M,M_k_p__,11,' M_k_p___: ')
         write(6,'(A,I0)') ' n_CTF: ',n_CTF
         write(6,'(A,I0)') ' ld_CTF: ',ld_CTF
c$$$         call write_sub_c16(n_CTF*ld_CTF,CTF_k_p__,13,' CTF_k_p___: ')
         write(6,'(A,I0)') ' n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
c$$$         write(6,'(A)') ' alpha_est__: '
c$$$         write(6,format_string) (alpha_est__(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') ' alpha_update_f: ' , alpha_update_f
         write(6,'(A,L1)') ' flag_MS_vs_SM: ' , flag_MS_vs_SM
c$$$         if (flag_MS_vs_SM.eqv..true.) then
c$$$            write(6,'(A,I0)') ' n_MS_max: ',n_MS_max
c$$$            call write_sub_i4(n_S,n_MS_,8,' n_MS_: ')
c$$$            call write_sub_r8(n_S*n_MS_max*n_alpha,alpha_MS__,13
c$$$     $        ,' alpha_MS__: ')
c$$$         end if !if (flag_MS_vs_SM.eqv..true.) then
c$$$         if (flag_MS_vs_SM.eqv..false.) then
c$$$            write(6,'(A,I0)') ' n_SM_max: ',n_SM_max
c$$$            call write_sub_i4(n_M,n_SM_,8,' n_SM_: ')
c$$$            call write_sub_r8(n_M*n_SM_max*n_alpha,alpha_SM__,13
c$$$     $        ,' alpha_SM__: ')
c$$$         end if !if (flag_MS_vs_SM.eqv..false.) then
         write(6,'(A,F6.3)') ' n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') ' displacement_max: ',displacement_max
         write(6,'(A,I0)') ' n_delta_v: ',n_delta_v
         write(6,'(A,I0)') ' n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') ' svd_calculation_type: '
     $        ,svd_calculation_type
         write(6,'(A,F6.3)') ' eps_svd: ',eps_svd
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' , flag_RTRT_vs_RTTR
         write(6,'(A,I0)') ' fpm_howmany_max: ',fpm_howmany_max
         write(6,'(A,I0)') ' n_omp_sub__in: ',n_omp_sub__in
         write(6,'(A,I0)') ' n_S_0_sub__in: ',n_S_0_sub__in
         write(6,'(A,I0)') ' n_S_1_sub__in: ',n_S_1_sub__in
         write(6,'(A,I0)') ' n_M_0_sub__in: ',n_M_0_sub__in
         write(6,'(A,I0)') ' n_M_1_sub__in: ',n_M_1_sub__in
      end if !if (verbose.gt.1) then

       if (verbose.gt.0) then
          write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
       end if !if (verbose.gt.0) then
       call openblas_set_num_threads(1)
