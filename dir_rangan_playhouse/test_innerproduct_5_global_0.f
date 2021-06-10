      subroutine test_innerproduct_5_global_0(
      include 'ti5_list_var_external_0.f'
      include 'ti5_list_var_internal_0.f'
     $ )
      include 'ti5_define_var_allocation_off_0.f'
      do nM_0_sub=0,n_M_0_sub_use-1
         nM_0_per = n_M_0_per_(nM_0_sub)
         nM_0_sum = n_M_0_sum_(nM_0_sub)
         n_M_9_sub_use = min(n_omp_sub__in,nM_0_per)
         call block_0(verbose-2,n_M_9_sub_use,nM_0_per,n_M_9_per_
     $        ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)
         include 'ti5_extract_alpha_0.f'
         include 'ti5_define_O_T_R_CTF_M_q_0.f'
         include 'ti5_define_T_T_R_CTF_M_q_0.f'
         include 'ti5_define_Z_T_R_CTF_M_q_0.f'
         n_S_1_sub_use = min(n_S_1_sub__in,nS_0_per)
         call block_0(verbose-2,n_S_1_sub_use,nS_0_per,n_S_1_per_
     $        ,n_S_1_sum_,n_S_1_sub_use,nS_1_per_min,nS_1_per_max)
         n_M_1_sub_use = min(n_M_1_sub__in,nM_0_per)
         call block_0(verbose-2,n_M_1_sub_use,nM_0_per,n_M_1_per_
     $        ,n_M_1_sum_,n_M_1_sub_use,nM_1_per_min,nM_1_per_max)
         do nS_1_sub=0,n_S_1_sub_use-1
            nS_1_per = n_S_1_per_(nS_1_sub)
            nS_1_sum = n_S_1_sum_(nS_1_sub)
            n_S_9_sub_use = min(n_omp_sub__in,nS_1_per)
            call block_0(verbose-2,n_S_9_sub_use,nS_1_per,n_S_9_per_
     $           ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)
            do nM_1_sub=0,n_M_1_sub_use-1
               nM_1_per = n_M_1_per_(nM_1_sub)
               nM_1_sum = n_M_1_sum_(nM_1_sub)
               n_M_9_sub_use = min(n_omp_sub__in,nM_1_per)
               call block_0(verbose-2,n_M_9_sub_use,nM_1_per
     $              ,n_M_9_per_,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min
     $              ,nM_9_per_max)
               include 'ti5_block_display_0.f'
               include 'ti5_STxTRM_0.f'
               include 'ti5_SxTTRM_0.f'
               include 'ti5_SZxTRM_0.f'
               include 'ti5_SxZTRM_0.f'
               if (flag_time_Zstore) then
                  include 'ti5_Zstore_0a.f'
                  include 'ti5_Zstore_0b.f'
                  include 'ti5_Zstore_0c.f'
                  include 'ti5_Zstore_0x.f'
                  include 'ti5_Zstore_0d.f'
                  include 'ti5_Zstore_0e.f'
               end if           !if (flag_time_Zstore) then
               include 'ti5_Zstore_1.f'
            enddo               !do nM_1_sub=0,n_M_1_sub_use-1
         enddo                  !do nS_1_sub=0,n_S_1_sub_use-1
      enddo                     !do nM_0_sub=0,n_M_0_sub_use-1
      end



