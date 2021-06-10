!> Doxygen comment: ;\n
!> This program actually calculates the innerproducts ;\n
!> across n_S templates and n_M images: ;\n
!> ;\n
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
!> ;\n
!> Local search proceeds as follows: ;\n
!> For each 0-level S-block: ;\n
!>     Define 9-level M-blocking. ;\n
!>     Use omp to parallelize across 9-level M-blocks. ;\n
!>     For each 9-level M-block: ;\n
!>         For each image in M-block: ;\n
!>             Search tesselation tree for appropriate templates. ;\n
!>             Calculate innerproducts across those templates. ;\n
!>             Here we use nS_1_per_max as an upper bound for the ;\n
!>             number of templates processed at a single time. ;\n
!>         end image. ;\n
!>     end 9-level M-block. ;\n
!> end 0-level S-block. ;\n
      subroutine test_innerproduct_8(
      include 'ti8_list_var_external_0.f'
     $     )
      include 'ti8_define_var_allocation_0on_0.f'
      include 'ti8_entering_0.f'
      include 'ti8_allocate_variable_0.f'
      include 'ti8_define_svdd_3.f'
      include 'ti8_define_fftw_local_0.f'
      include 'ti8_define_fftw_omp_0.f'
      include 'ti8_define_fftw_fpm_0.f'
      include 'ti8_checkset_variable_0.f'
      do nS_0_sub=0,n_S_0_sub_use-1
         nS_0_per = n_S_0_per_(nS_0_sub)
         nS_0_sum = n_S_0_sum_(nS_0_sub)
         n_S_9_sub_use = min(n_omp_sub_0in,nS_0_per)
      include 'ti8_checkset_variable_0.f'
         call block_0(verbose-2,n_S_9_sub_use,nS_0_per,n_S_9_per_
     $        ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)
      include 'ti8_checkset_variable_0.f'
         include 'ti8_define_CTF_R_S_0.f'
      include 'ti8_checkset_variable_0.f'
         include 'ti8_define_O_S_q_0.f'
      include 'ti8_checkset_variable_0.f'
         include 'ti8_define_T_S_q_0.f'
      include 'ti8_checkset_variable_0.f'
         include 'ti8_define_Z_S_q_0.f'
      include 'ti8_checkset_variable_0.f'
         if (tesselation_distance_req.ge.2.0d0) then
            call ti8_global_0(
            include 'ti8_list_var_external_0.f'
            include 'ti8_list_var_internal_0.f'
     $           )
         else !if (tesselation_distance_req.ge.2.0d0) then
            include 'ti8_define_tesselation_0.f'
            call ti8_local_0(
            include 'ti8_list_var_external_0.f'
            include 'ti8_list_var_internal_0.f'
     $           )
         end if !if (tesselation_distance_req.ge.2.0d0) then
      include 'ti8_checkset_variable_0.f'
      enddo !do nS_0_sub=0,n_S_0_sub_use-1
      include 'ti8_checkset_variable_0.f'
      include 'ti8_deallocate_0.f'
      include 'ti8_finished_1.f'
