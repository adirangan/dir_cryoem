!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Global search proceeds as follows: ;\n
!> For each 0-level S-block: ;\n
!>     For each 0-level M-block: ;\n
!>         Define 1-level S-blocking. ;\n
!>         Define 1-level M-blocking. ;\n
!>         For each 1-level S-block: ;\n
!>             For each 1-level M-block: ;\n
!>                 Calculate innerproducts using omp. ;\n
!>                 (i.e., using 9-level S- and M-blocks). ;\n
!>             end 1-level M-block. ;\n
!>         end 1-level S-block. ;\n
!>     end 0-level M-block. ;\n
!> end 0-level S-block. ;\n
      subroutine ti8_global_0(
      include 'ti8_excerpt_list_var_external_0.f'
      include 'ti8_excerpt_list_var_internal_0.f'
     $ )
      include 'ti8_excerpt_define_var_allocation_off_0.f'
      do nM_0_sub=0,n_M_0_sub_use-1
         nM_0_per = n_M_0_per_(nM_0_sub)
         nM_0_sum = n_M_0_sum_(nM_0_sub)
         n_M_9_sub_use = min(n_omp_sub_0in,nM_0_per)
         call block_0(verbose-2,n_M_9_sub_use,nM_0_per,n_M_9_per_
     $        ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_M_9_sub_use ' , n_M_9_sub_use ,
     $           ' nM_0_per ' , nM_0_per , ' n_M_9_sub_use ' ,
     $           n_M_9_sub_use , ' nM_9_per_min ' , nM_9_per_min ,
     $           ' nM_9_per_max ' , nM_9_per_max
            call print_all_i4(n_M_9_sub_use,n_M_9_per_
     $           ,' n_M_9_per_: ')
            call print_all_i4(n_M_9_sub_use,n_M_9_sum_
     $           ,' n_M_9_sum_: ')
         end if !if (verbose.gt.1) then
      include 'ti8_excerpt_checkset_variable_0.f'
         include 'ti8_excerpt_extract_alpha_sub_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
         include 'ti8_excerpt_define_O_T_R_CTF_M_q_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
         include 'ti8_excerpt_define_T_T_R_CTF_M_q_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
         include 'ti8_excerpt_define_Z_T_R_CTF_M_q_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
         n_S_1_sub_use = min(n_S_1_sub_0in,nS_0_per)
         call block_0(verbose-2,n_S_1_sub_use,nS_0_per,n_S_1_per_
     $        ,n_S_1_sum_,n_S_1_sub_use,nS_1_per_min,nS_1_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_S_1_sub_use ' , n_S_1_sub_use ,
     $           ' nS_0_per ' , nS_0_per , ' n_S_1_sub_use ' ,
     $           n_S_1_sub_use , ' nS_1_per_min ' , nS_1_per_min ,
     $           ' nS_1_per_max ' , nS_1_per_max
            call print_all_i4(n_S_1_sub_use,n_S_1_per_
     $           ,' n_S_1_per_: ')
            call print_all_i4(n_S_1_sub_use,n_S_1_sum_
     $           ,' n_S_1_sum_: ')
         end if !if (verbose.gt.1) then
         n_M_1_sub_use = min(n_M_1_sub_0in,nM_0_per)
         call block_0(verbose-2,n_M_1_sub_use,nM_0_per,n_M_1_per_
     $        ,n_M_1_sum_,n_M_1_sub_use,nM_1_per_min,nM_1_per_max)
         if (verbose.gt.1) then
            write(6,'(5(A,I0))') ' n_M_1_sub_use ' , n_M_1_sub_use ,
     $           ' nM_0_per ' , nM_0_per , ' n_M_1_sub_use ' ,
     $           n_M_1_sub_use , ' nM_1_per_min ' , nM_1_per_min ,
     $           ' nM_1_per_max ' , nM_1_per_max
            call print_all_i4(n_M_1_sub_use,n_M_1_per_
     $           ,' n_M_1_per_: ')
            call print_all_i4(n_M_1_sub_use,n_M_1_sum_
     $           ,' n_M_1_sum_: ')
         end if !if (verbose.gt.1) then
         do nS_1_sub=0,n_S_1_sub_use-1
            nS_1_per = n_S_1_per_(nS_1_sub)
            nS_1_sum = n_S_1_sum_(nS_1_sub)
            n_S_9_sub_use = min(n_omp_sub_0in,nS_1_per)
            call block_0(verbose-2,n_S_9_sub_use,nS_1_per,n_S_9_per_
     $           ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)
            do nM_1_sub=0,n_M_1_sub_use-1
               nM_1_per = n_M_1_per_(nM_1_sub)
               nM_1_sum = n_M_1_sum_(nM_1_sub)
               n_M_9_sub_use = min(n_omp_sub_0in,nM_1_per)
               call block_0(verbose-2,n_M_9_sub_use,nM_1_per
     $              ,n_M_9_per_,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min
     $              ,nM_9_per_max)
               include 'ti8_excerpt_block_display_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
               include 'ti8_excerpt_define_STxTRM_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
               include 'ti8_excerpt_define_SxTTRM_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
               include 'ti8_excerpt_define_SZxTRM_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
               include 'ti8_excerpt_define_SxZTRM_0.f'
      include 'ti8_excerpt_checkset_variable_0.f'
               n_S_use_sum = n_S_use_sum + nM_1_per*nS_1_per
               include 'ti8_excerpt_define_Zstore_part_0.f'
               include 'ti8_excerpt_define_Zstore_full_0.f'
            enddo               !do nM_1_sub=0,n_M_1_sub_use-1
         enddo                  !do nS_1_sub=0,n_S_1_sub_use-1
      enddo                     !do nM_0_sub=0,n_M_0_sub_use-1
      n_S_9_sub_use = min(n_omp_sub_0in,nS_0_per)
      call block_0(verbose-2,n_S_9_sub_use,nS_0_per,n_S_9_per_
     $     ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)
      include 'ti8_excerpt_merge_alpha_MS_0.f'
      end



