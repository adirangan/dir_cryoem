!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
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
      subroutine ti8_local_0(
      include 'ti8_excerpt_list_var_external_0.f'
      include 'ti8_excerpt_list_var_internal_0.f'
     $ )
      include 'ti8_excerpt_define_var_allocation_off_0.f'
      include 'ti8_excerpt_extract_alpha_all_0.f'
      n_M_9_sub_use = min(n_omp_sub_0in,n_M)
      call block_0(verbose-2,n_M_9_sub_use,n_M,n_M_9_per_
     $     ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)

      if (verbose.gt.2) then
         write(6,'(A,I0)') ' n_M_9_sub_use: ' , n_M_9_sub_use      
         write(6,'(A,I0)') ' n_M_9_sum_(0): ' , n_M_9_sum_(0)
         call print_sub_c16(n_w_sum,S_k_p__,' S_k_p__: ')
         call print_sub_c16(n_w_sum,M_k_p__,' M_k_p__: ')
         call print_sub_c16(n_w_sum,CTF_k_p__,' CTF_k_p__: ')
         write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
         call print_sub_r8(n_delta_v,delta_x_,' delta_x_: ');
         call print_sub_r8(n_delta_v,delta_y_,' delta_y_: ');
         write(6,'(A,I0)') ' n_r: ' , n_r
         call print_sub_r8(n_r,grid_k_p_r_,' grid_k_p_r_: ');
         write(6,'(A,I0)') ' n_w_sum: ' , n_w_sum
         write(6,'(A,F8.4)') ' delta_x_est_(0): ' ,
     $        delta_x_est_(0)
         write(6,'(A,F8.4)') ' delta_y_est_(0): ' ,
     $        delta_y_est_(0)
         write(6,'(A,F8.4)') ' gamma_z_est_(0): ' ,
     $        gamma_z_est_(0)
         write(6,'(A,F8.4)') ' ctf_ind_est_(0): ' ,
     $        ctf_ind_est_(0)
         call print_sub_c16(n_w_sum,Z_p_omp__,' Z_p_omp__: ')
         call print_sub_c16(n_w_sum,M_p_omp__,' M_p_omp__: ')
         call print_sub_c16(n_w_sum,Z_q_omp__,' Z_q_omp__: ')
         call print_sub_c16(n_w_sum,M_q_omp__,' M_q_omp__: ')
      end if                    ! if (verbose.gt.2) then

c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
      include 'ti8l_excerpt_private_0.f'
c$OMP&)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         include 'ti8l_excerpt_define_local_pointer_0.f'
         do nm=0,nM_9_per-1
            include 'ti8l_excerpt_define_O_T_R_CTF_M_q_0.f'
            include 'ti8l_excerpt_define_T_T_R_CTF_M_q_0.f'
            include 'ti8l_excerpt_define_Z_T_R_CTF_M_q_0.f'
            include 'ti8l_excerpt_search_tesselation_1.f'
         enddo !do nm=0,nM_9_per-1
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      n_S_use_sum = sum_i4_f(n_M_9_sub_use,n_S_use_sum_)
      end     
