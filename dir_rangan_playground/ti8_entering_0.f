!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Prints out several parameters at the beginning of test_innerproduct_8. ;\n
      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[entering test_innerproduct_8]: '
      end if

      if (verbose.gt.-1) then
         write(6,'(A,I0)') ' verbose: ',verbose
         write(6,'(A,L1)') ' flag_memory_estimate: '
     $        ,flag_memory_estimate
         write(6,'(A,I0)') ' rseed: ',rseed
         write(6,'(A,I0)') ' n_k_cur: ',n_k_cur
         call print_all_i4(n_k_cur,n_polar_a_,' n_polar_a_: ')
         call print_all_r8(n_k_cur,grid_k_p_,' grid_k_p_: ')
         call print_all_r8(n_k_cur,weight_k_p_r_,' weight_k_p_r_: ')
         call print_sub_r8(ld_S,weight_k_p_,' weight_k_p_: ')
         write(6,'(A,F6.3)') ' half_diameter_x_c: ',half_diameter_x_c
         write(6,'(A,I0)') ' n_S: ',n_S
c$$$         call print_sub_i4(n_S,I_S_sample_,' I_S_sample_: ')
         write(6,'(A,I0)') ' ld_S: ',ld_S
c$$$         call print_sub_c16(n_S*ld_S,S_k_p__,' S_k_p___: ')
         write(6,'(A,F6.3)') ' tesselation_distance_req: '
     $        ,tesselation_distance_req
c$$$         call print_sub_r8(n_S,S_alpha_S_index_
c$$$     $        ,' S_alpha_S_index_: ')
c$$$         call print_sub_r8(n_S,S_alpha_polar_a_
c$$$     $        ,' S_alpha_polar_a_: ')
c$$$         call print_sub_r8(n_S,S_alpha_azimu_b_
c$$$     $        ,' S_alpha_azimu_b_: ')
         write(6,'(A,I0)') ' n_M: ',n_M
c$$$         call print_sub_i4(n_M,I_M_sample_,' I_M_sample_: ')
         write(6,'(A,I0)') ' ld_M: ',ld_M
c$$$         call print_sub_c16(n_M*ld_M,M_k_p__,' M_k_p___: ')
         write(6,'(A,I0)') ' n_CTF: ',n_CTF
         write(6,'(A,I0)') ' ld_CTF: ',ld_CTF
c$$$         call print_sub_c16(n_CTF*ld_CTF,CTF_k_p__,' CTF_k_p___: ')
         write(6,'(A,I0)') ' n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
c$$$         write(6,'(A)') ' alpha_est__: '
c$$$         write(6,format_string) (alpha_est__(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') ' alpha_update_f: ' , alpha_update_f
         write(6,'(A,L1)') ' flag_MS_vs_SM: ' , flag_MS_vs_SM
c$$$         if (flag_MS_vs_SM.eqv..true.) then
c$$$            write(6,'(A,I0)') ' n_MS_max: ',n_MS_max
c$$$            call print_sub_i4(n_S,n_MS_,' n_MS_: ')
c$$$            call print_sub_r8(n_S*n_MS_max*n_alpha,alpha_MS__
c$$$     $        ,' alpha_MS__: ')
c$$$         end if !if (flag_MS_vs_SM.eqv..true.) then
c$$$         if (flag_MS_vs_SM.eqv..false.) then
c$$$            write(6,'(A,I0)') ' n_SM_max: ',n_SM_max
c$$$            call print_sub_i4(n_M,n_SM_,' n_SM_: ')
c$$$            call print_sub_r8(n_M*n_SM_max*n_alpha,alpha_SM__
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
         write(6,'(A,I0)') ' n_omp_sub_0in: ',n_omp_sub_0in
         write(6,'(A,I0)') ' n_S_0_sub_0in: ',n_S_0_sub_0in
         write(6,'(A,I0)') ' n_S_1_sub_0in: ',n_S_1_sub_0in
         write(6,'(A,I0)') ' n_M_0_sub_0in: ',n_M_0_sub_0in
         write(6,'(A,I0)') ' n_M_1_sub_0in: ',n_M_1_sub_0in
      end if !if (verbose.gt.1) then

       if (verbose.gt.0) then
          write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
       end if !if (verbose.gt.0) then
       call openblas_set_num_threads(1)
