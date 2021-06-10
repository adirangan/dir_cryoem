      subroutine test_innerproduct_7_global_0(
      include 'ti7_list_var_external_0.f'
      include 'ti7_list_var_internal_0.f'
     $ )
      include 'ti7_define_var_allocation_off_0.f'
      do nM_0_sub=0,n_M_0_sub_use-1
         nM_0_per = n_M_0_per_(nM_0_sub)
         nM_0_sum = n_M_0_sum_(nM_0_sub)
         n_M_9_sub_use = min(n_omp_sub__in,nM_0_per)
         call block_0(verbose-2,n_M_9_sub_use,nM_0_per,n_M_9_per_
     $        ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_M_9_sub_use ' , n_M_9_sub_use ,
     $           ' nM_0_per ' , nM_0_per , ' n_M_9_sub_use ' ,
     $           n_M_9_sub_use , ' nM_9_per_min ' , nM_9_per_min ,
     $           ' nM_9_per_max ' , nM_9_per_max
            call write_all_i4(n_M_9_sub_use,n_M_9_per_,13
     $           ,' n_M_9_per_: ')
            call write_all_i4(n_M_9_sub_use,n_M_9_sum_,13
     $           ,' n_M_9_sum_: ')
         end if !if (verbose.gt.1) then
      include 'ti7_checkset_variable_0.f'
         include 'ti5_extract_alpha_0.f'
      include 'ti7_checkset_variable_0.f'
         include 'ti7_define_O_T_R_CTF_M_q_0.f'
      include 'ti7_checkset_variable_0.f'
         include 'ti7_define_T_T_R_CTF_M_q_0.f'
      include 'ti7_checkset_variable_0.f'
         include 'ti7_define_Z_T_R_CTF_M_q_0.f'
      include 'ti7_checkset_variable_0.f'
         n_S_1_sub_use = min(n_S_1_sub__in,nS_0_per)
         call block_0(verbose-2,n_S_1_sub_use,nS_0_per,n_S_1_per_
     $        ,n_S_1_sum_,n_S_1_sub_use,nS_1_per_min,nS_1_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_S_1_sub_use ' , n_S_1_sub_use ,
     $           ' nS_0_per ' , nS_0_per , ' n_S_1_sub_use ' ,
     $           n_S_1_sub_use , ' nS_1_per_min ' , nS_1_per_min ,
     $           ' nS_1_per_max ' , nS_1_per_max
            call write_all_i4(n_S_1_sub_use,n_S_1_per_,13
     $           ,' n_S_1_per_: ')
            call write_all_i4(n_S_1_sub_use,n_S_1_sum_,13
     $           ,' n_S_1_sum_: ')
         end if !if (verbose.gt.1) then
         n_M_1_sub_use = min(n_M_1_sub__in,nM_0_per)
         call block_0(verbose-2,n_M_1_sub_use,nM_0_per,n_M_1_per_
     $        ,n_M_1_sum_,n_M_1_sub_use,nM_1_per_min,nM_1_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_M_1_sub_use ' , n_M_1_sub_use ,
     $           ' nM_0_per ' , nM_0_per , ' n_M_1_sub_use ' ,
     $           n_M_1_sub_use , ' nM_1_per_min ' , nM_1_per_min ,
     $           ' nM_1_per_max ' , nM_1_per_max
            call write_all_i4(n_M_1_sub_use,n_M_1_per_,13
     $           ,' n_M_1_per_: ')
            call write_all_i4(n_M_1_sub_use,n_M_1_sum_,13
     $           ,' n_M_1_sum_: ')
         end if !if (verbose.gt.1) then
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
      include 'ti7_checkset_variable_0.f'
               include 'ti7_STxTRM_0.f'
      include 'ti7_checkset_variable_0.f'
               include 'ti7_SxTTRM_0.f'
      include 'ti7_checkset_variable_0.f'
               include 'ti7_SZxTRM_0.f'
      include 'ti7_checkset_variable_0.f'
               include 'ti7_SxZTRM_0.f'
      include 'ti7_checkset_variable_0.f'
               n_S_use_sum = n_S_use_sum + nM_1_per*nS_1_per
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



