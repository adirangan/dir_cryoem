      subroutine test_innerproduct_7_local_0(
      include 'ti7_list_var_external_0.f'
      include 'ti7_list_var_internal_0.f'
     $ )
      include 'ti7_define_var_allocation_off_0.f'
      include 'ti7_extract_alpha_0.f'
      n_M_9_sub_use = min(n_omp_sub__in,n_M)
      call block_0(verbose-2,n_M_9_sub_use,n_M,n_M_9_per_
     $     ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)

      write(6,'(A,I0)') ' n_M_9_sub_use: ' , n_M_9_sub_use      
      write(6,'(A,I0)') ' n_M_9_sum_(0): ' , n_M_9_sum_(0)
      call block_0(verbose-2,n_M_9_sub_use,n_M,n_M_9_per_
     $     ,n_M_9_sum_,n_M_9_sub_use,nM_9_per_min,nM_9_per_max)
      call print_sub_c16(n_A,S_k_p__,' S_k_p__: ')
      call print_sub_c16(n_A,M_k_p__,' M_k_p__: ')
      call print_sub_c16(n_A,CTF_k_p__,' CTF_k_p__: ')
      write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
      call print_sub_r8(n_delta_v,delta_x_,' delta_x_: ');
      call print_sub_r8(n_delta_v,delta_y_,' delta_y_: ');
      write(6,'(A,I0)') ' n_r: ' , n_r
      call print_sub_r8(n_r,grid_k_p_,' grid_k_p_: ');
      write(6,'(A,I0)') ' n_A: ' , n_A
      write(6,'(A,F8.4)') ' delta_x_est_(0): ' ,
     $     delta_x_est_(0)
      write(6,'(A,F8.4)') ' delta_y_est_(0): ' ,
     $     delta_y_est_(0)
      write(6,'(A,F8.4)') ' gamma_z_est_(0): ' ,
     $     gamma_z_est_(0)
      write(6,'(A,F8.4)') ' ctf_ind_est_(0): ' ,
     $     ctf_ind_est_(0)
      call print_sub_c16(n_A,Z_p_omp__,' Z_p_omp__: ')
      call print_sub_c16(n_A,M_p_omp__,' M_p_omp__: ')
      call print_sub_c16(n_A,Z_q_omp__,' Z_q_omp__: ')
      call print_sub_c16(n_A,M_q_omp__,' M_q_omp__: ')

c$OMP PARALLEL PRIVATE(nM_9_per,nM_9_sum,
      include 'ti7l_private_0.f'
c$OMP&)
c$OMP DO
      do nM_9_sub=0,n_M_9_sub_use-1
         nM_9_per = n_M_9_per_(nM_9_sub)
         nM_9_sum = n_M_9_sum_(nM_9_sub)
         include 'ti7l_define_sub_0.f'
         do nm=0,nM_9_per-1
            include 'ti7l_define_O_T_R_CTF_M_q_0.f'
            include 'ti7l_define_T_T_R_CTF_M_q_0.f'
            include 'ti7l_define_Z_T_R_CTF_M_q_0.f'
            include 'ti7l_search_tesselation_1.f'
         enddo !do nm=0,nM_9_per-1
      enddo !do nM_9_sub=0,n_M_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      n_S_use_sum = sum_i4_f(n_M_9_sub_use,n_S_use_sum_)
      end     

