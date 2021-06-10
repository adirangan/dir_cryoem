      subroutine test_innerproduct_5(
      include 'ti5_list_var_external_0.f'
     $     )
      include 'ti5_define_var_allocation__on_0.f'
      include 'ti5_entering_0.f'
      include 'ti5_allocate_variable_0.f'
      include 'ti5_define_svdd_1.f'
      include 'ti5_define_fftw_local_0.f'
      include 'ti5_define_fftw_omp_0.f'
      include 'ti5_define_fftw_fpm_0.f'
      do nS_0_sub=0,n_S_0_sub_use-1
         nS_0_per = n_S_0_per_(nS_0_sub)
         nS_0_sum = n_S_0_sum_(nS_0_sub)
         n_S_9_sub_use = min(n_omp_sub__in,nS_0_per)
         call block_0(verbose-2,n_S_9_sub_use,nS_0_per,n_S_9_per_
     $        ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)
         include 'ti5_define_CTF_R_S_0.f'
         include 'ti5_define_O_S_q_0.f'
         include 'ti5_define_T_S_q_0.f'
         include 'ti5_define_Z_S_q_0.f'
         if (tesselation_distance_req.ge.2.0d0) then
            call test_innerproduct_5_global_0(
            include 'ti5_list_var_external_0.f'
            include 'ti5_list_var_internal_0.f'
     $           )
         else !if (tesselation_distance_req.ge.2.0d0) then
            include 'ti5_define_tesselation_0.f'
            call test_innerproduct_5_local_0(
            include 'ti5_list_var_external_0.f'
            include 'ti5_list_var_internal_0.f'
     $           )
         end if !if (tesselation_distance_req.ge.2.0d0) then
      enddo !do nS_0_sub=0,n_S_0_sub_use-1
      include 'ti5_deallocate_0.f'
      include 'ti5_finished_1.f'
