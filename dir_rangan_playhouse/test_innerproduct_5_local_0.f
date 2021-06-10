      subroutine test_innerproduct_5_local_0(
      include 'ti5_list_var_external_0.f'
      include 'ti5_list_var_internal_0.f'
     $ )
      include 'ti5_define_var_allocation_off_0.f'
      include 'ti6_extract_alpha_0.f'
      n_M_9_sub_use = min(n_omp_sub__in,n_M)
      call block_0(verbose-2,n_M_9_sub_use,n_M,n_M_9_per_
     $     ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         include 'ti6_define_sub_0.f'
         do nm=0,nM_9_per-1
            include 'ti6_define_O_T_R_CTF_M_q_0.f'
            include 'ti6_define_T_T_R_CTF_M_q_0.f'
            include 'ti6_define_Z_T_R_CTF_M_q_0.f'
            include 'ti6_search_tesselation_1.f'
         enddo !do nm=0,nM_9_per-1
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
      end     

