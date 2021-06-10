      pi = 4*atan(1.0)
      eps_target = eps_svd
      p_S_p__ = loc(S_k_p__(0))
      p_M_p__ = loc(M_k_p__(0))
      p_CTF_p__ = loc(CTF_k_p__(0))
      d_memory_estimate = 0.0d0
      if (flag_memory_estimate.eqv..true.) then
         verbose_mem = 1
      else
         verbose_mem = 0
      end if !if (flag_memory_estimate.eqv..true.) then
      verbose_timing = 0
      timing_total_fftw_plan = 0.0d0
      gnump_total_fftw_plan = 0.0d0
      timing_total_CTF_R_S = 0.0d0
      gnump_total_CTF_R_S = 0.0d0
      timing_total_S_q = 0.0d0
      gnump_total_S_q = 0.0d0
      timing_total_T_S_q = 0.0d0
      gnump_total_T_S_q = 0.0d0
      timing_total_Z_S_q = 0.0d0
      gnump_total_Z_S_q = 0.0d0
      timing_total_T_R_CTF_M_q = 0.0d0
      gnump_total_T_R_CTF_M_q = 0.0d0
      timing_total_T_T_R_CTF_M_q = 0.0d0
      gnump_total_T_T_R_CTF_M_q = 0.0d0
      timing_total_Z_T_R_CTF_M_q = 0.0d0
      gnump_total_Z_T_R_CTF_M_q = 0.0d0
      timing_total_zgemm = 0.0d0
      gnump_total_zgemm = 0.0d0
      timing_total_fpm_fill = 0.0d0
      gnump_total_fpm_fill = 0.0d0
      timing_total_fpm_fftw = 0.0d0
      gnump_total_fpm_fftw = 0.0d0
      timing_total_transpose = 0.0d0
      gnump_total_transpose = 0.0d0
      timing_total_Zstore = 0.0d0
      gnump_total_Zstore = 0.0d0
      flag_time_Zstore = .false.
      timing_total_Zstore_a = 0.0d0
      gnump_total_Zstore_a = 0.0d0
      timing_total_Zstore_b = 0.0d0
      gnump_total_Zstore_b = 0.0d0
      timing_total_Zstore_c = 0.0d0
      gnump_total_Zstore_c = 0.0d0
      timing_total_Zstore_x = 0.0d0
      gnump_total_Zstore_x = 0.0d0
      timing_total_Zstore_d = 0.0d0
      gnump_total_Zstore_d = 0.0d0
      timing_total_Zstore_e = 0.0d0
      gnump_total_Zstore_e = 0.0d0

c$$$      Calculating template size using 'get_template_size'
      allocate(n_w_(0:n_k_cur-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*n_k_cur
      allocate(n_w_sum_(0:n_k_cur-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*n_k_cur
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') ' n_k_cur = ',n_k_cur
     $        ,' calling get_template_size to'
     $        ,' determine n_w_ and n_A'
      end if !if (verbose.gt.1) then
      call get_template_size(n_polar_a_,n_k_cur,n_A,n_w_,n_w_sum_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_A = ',n_A
         call write_all_i4(n_k_cur,n_w_,7,' n_w_: ')
      end if !if (verbose.gt.1) then
      
c$$$      indices
      n_r = n_k_cur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r ' , n_r , ' < 2'
      end if !if (n_r.lt.2) then
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo !do nr=0,n_r-1
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') ' n_w_max ',n_w_max,'; n_A ',n_A
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      allocate(delta_x_(0:n_delta_x-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*n_delta_x
      allocate(delta_y_(0:n_delta_y-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*n_delta_y
      call get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)') ' n_w_max or so.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      allocate(gamma_z_(0:n_gamma_z-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*n_gamma_z
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.1) then
         call write_all_r8(n_delta_x,delta_x_,11,' delta_x_: ')
         call write_all_r8(n_delta_y,delta_y_,11,' delta_y_: ')
         call write_sub_r8(n_gamma_z,gamma_z_,11,' gamma_z_: ')
      end if !if (verbose.gt.1) then

      if ((svd_calculation_type.eq.1)) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if !if (verbose.gt.1) then
      allocate(svd_r_(0:n_svd_max-1))
      allocate(svd_d_(0:n_svd_max-1))
      allocate(svd_l_(0:n_svd_max-1))
      allocate(svd_U_d_(0:n_svd_max*n_svd_max-1))
      allocate(svd_s_(0:n_svd_max-1))
      allocate(svd_V_r_(0:n_svd_max*n_svd_max-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*n_svd_max*(4.0d0 +
     $     2.0d0*n_svd_max)
      call get_svd_0(eps_target,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,grid_k_p_,n_r,n_delta_x,delta_x_,n_delta_y
     $     ,delta_y_,flag_warning,R_max ,K_max,delta_max,n_pixels)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_svd_r: ',n_svd_r
         write(6,'(A,I0)') ' n_svd_d: ',n_svd_d
         write(6,'(A,I0)') ' svd_unitnumber: ',svd_unitnumber
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      if (verbose.gt.0) then
         write(6,'(A,I0,A,A,I0)') ' n_svd_l: ' , n_svd_l , ' versus '
     $        ,' n_delta_x*n_delta_y: ' , n_delta_x*n_delta_y
      end if !if (verbose.gt.0) then      
      end if !if ((svd_calculation_type.eq.1)) then

      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_x*n_delta_y
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0)') ' setting n_transf: ' , n_transf
     $           ,' = n_delta_x*n_delta_y = ' , n_delta_x , '*' ,
     $           n_delta_y
         end if !if (verbose.gt.1) then
      else
         n_transf = n_svd_l
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') ' setting n_transf: ' , n_transf ,
     $           ' = n_svd_l = ' , n_svd_l
         end if !if (verbose.gt.1) then
      end if !if (svd_calculation_type.eq.2) then

      allocate(S_p_(0:n_A-1))
      allocate(S_q_(0:n_A-1))
      allocate(S_p_omp__(0:n_A*n_omp_sub__in-1))
      allocate(S_q_omp__(0:n_A*n_omp_sub__in-1))
      allocate(M_p_(0:n_A-1))
      allocate(M_q_(0:n_A-1))
      allocate(M_p_omp__(0:n_A*n_omp_sub__in-1))
      allocate(M_q_omp__(0:n_A*n_omp_sub__in-1))
      allocate(CTF_p_(0:n_A-1))
      allocate(CTF_q_(0:n_A-1))
      allocate(CTF_p_omp__(0:n_A*n_omp_sub__in-1))
      allocate(CTF_q_omp__(0:n_A*n_omp_sub__in-1))
      allocate(Z_p_(0:n_A-1))
      allocate(Z_q_(0:n_A-1))
      allocate(Z_p_omp__(0:n_A*n_omp_sub__in-1))
      allocate(Z_q_omp__(0:n_A*n_omp_sub__in-1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_A*(8.0d0 + 8.0d0
     $     *n_omp_sub__in)
      allocate(ZZ__(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(ZZ_omp__(0:n_delta_x*n_delta_y*n_gamma_z*n_omp_sub__in
     $     -1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_x*1.0d0
     $     *n_delta_y*1.0d0*n_gamma_z*(1.0d0 + n_omp_sub__in)
      allocate(ZZ_sub_(0:n_w_max-1))
      allocate(ZZ_sub_omp__(0:n_w_max*n_omp_sub__in-1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_w_max*(1.0d0 +
     $     n_omp_sub__in)
      allocate(CTF_R_S_sub__(0:n_gamma_z-1))
      allocate(CTF_R_S_sub_omp__(0:n_gamma_z*n_omp_sub__in-1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_gamma_z*(1.0d0 +
     $     n_omp_sub__in)

      allocate(n_S_0_per_(0:n_S-1))
      allocate(n_S_0_sum_(0:n_S-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*2.0d0*n_S
      n_S_0_sub_use = min(n_S_0_sub__in,n_S)
      call block_0(verbose-2,n_S_0_sub_use,n_S,n_S_0_per_,n_S_0_sum_
     $     ,n_S_0_sub_use,nS_0_per_min,nS_0_per_max)

      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_gamma_z*1.0d0*nS_0_per_max*1.0d0*n_CTF
     $     *16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' CTF_R_S__ requires ' , d_mem*1.0d
     $        -6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(CTF_R_S__(0:n_gamma_z*nS_0_per_max*n_CTF-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_gamma_z*1.0d0
     $     *nS_0_per_max*1.0d0*n_CTF

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each template'
         write(6,'(A)') ' of the form: (S)_q, '
         write(6,'(A)') ' '
      end if                    !if (verbose.gt.1) then
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*nS_0_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_q__ requires ' , d_mem *1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_q__(0:n_r*nS_0_per_max*n_w_max))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $     *nS_0_per_max*1.0d0*n_w_max
      end if !if ((flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array T_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (T_{+delta_upd}(S))_q, '
            write(6,'(A)') ' where T = translation by +delta_upd.'         
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_0_per_max*1.0d0
     $     *n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' T_S_q__ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(T_S_q__(0:n_r*n_transf*nS_0_per_max*n_w_max))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0*n_transf
     $     *1.0d0*nS_0_per_max*1.0d0*n_w_max
      end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (Z_{rho}(S))_q, '
            write(6,'(A)') ' where Z represents the rho-side'
            write(6,'(A)') ' of the svd of the translation operator.'
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_0_per_max*1.0d0
     $     *n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' Z_S_q__ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(Z_S_q__(0:n_r*n_transf*nS_0_per_max*n_w_max))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0*n_transf
     $     *1.0d0*nS_0_per_max*1.0d0*n_w_max
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      allocate(n_M_0_per_(0:n_M-1))
      allocate(n_M_0_sum_(0:n_M-1))
      n_M_0_sub_use = min(n_M_0_sub__in,n_M)
      call block_0(verbose-2,n_M_0_sub_use,n_M,n_M_0_per_,n_M_0_sum_
     $     ,n_M_0_sub_use,nM_0_per_min,nM_0_per_max)

      allocate(polar_a_est_(0:nM_0_per_max-1))
      allocate(azimu_b_est_(0:nM_0_per_max-1))
      allocate(gamma_z_est_(0:nM_0_per_max-1))
      allocate(delta_x_est_(0:nM_0_per_max-1))
      allocate(delta_y_est_(0:nM_0_per_max-1))
      allocate(l2_norm_est_(0:nM_0_per_max-1))
      allocate(ctf_ind_est_(0:nM_0_per_max-1))
      allocate(S_index_est_(0:nM_0_per_max-1))
      allocate(M_index_est_(0:nM_0_per_max-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*9.0d0*nM_0_per_max
      allocate(alpha__in_(0:n_alpha-1))
      allocate(alpha__in_omp__(0:n_alpha*n_omp_sub__in-1))
      d_memory_estimate = d_memory_estimate + 8.0d0*n_alpha*(1.0d0 +
     $     n_omp_sub__in)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array T_R_CTF_M_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each image'
         write(6,'(A)') ' T_{-delta_est}(R_{-gamma_est}(CTF.*M)), '
         write(6,'(A)') ' where T = translation by -delta_est.'
         write(6,'(A)') ' and R = rotation by -gamma_est.'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*nM_0_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' T_R_CTF_M_q_ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(T_R_CTF_M_q__(0:n_r*nM_0_per_max*n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $     *nM_0_per_max*1.0d0*n_w_max
      end if !if ((flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array T_T_R_CTF_M_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each image'
         write(6,'(A,A)') ' T_{-delta_upd}(T_{-delta_est} ' ,
     $        ' (R_{-gamma_est}(CTF.*M))), '
         write(6,'(A)') ' where T = translation '
         write(6,'(A)') ' and R = rotation.'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nM_0_per_max*1.0d0
     $     *n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' T_T_R_CTF_M_q_ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(T_T_R_CTF_M_q__(0:n_r*n_transf*nM_0_per_max*n_w_max
     $        -1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0*n_transf
     $     *1.0d0*nM_0_per_max*1.0d0*n_w_max
      end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array Z_T_R_CTF_M_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each image'
         write(6,'(A)') ' Z(T_{-delta_est}(R_{-gamma_est}(CTF.*M))), '
         write(6,'(A)') ' where Z represents the rho-side'
         write(6,'(A)') ' of the svd of the translation operator.'
         write(6,'(A)') ' and T = translation'
         write(6,'(A)') ' and R = rotation.'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nM_0_per_max*1.0d0
     $     *n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' Z_T_R_CTF_M_q_ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(Z_T_R_CTF_M_q__(0:n_r*n_transf*nM_0_per_max*n_w_max
     $        -1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0*n_transf
     $     *1.0d0*nM_0_per_max*1.0d0*n_w_max
      end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      allocate(n_S_1_per_(0:nS_0_per_max-1))
      allocate(n_S_1_sum_(0:nS_0_per_max-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*2.0d0*nS_0_per_max
      n_S_1_sub_use = min(n_S_1_sub__in,nS_0_per_max)
      call block_0(verbose-2,n_S_1_sub_use,nS_0_per_max,n_S_1_per_
     $     ,n_S_1_sum_,n_S_1_sub_use,nS_1_per_min,nS_1_per_max)
      if (verbose.gt.-1) then
         write(6,'(A,I0)') ' n_S ' , n_S
         write(6,'(3(A,I0))') ' n_S_0_sub_use ' , n_S_0_sub_use ,
     $        ' nS_0_per_min ' , nS_0_per_min ,' nS_0_per_max ' ,
     $        nS_0_per_max
         write(6,'(3(A,I0))') ' n_S_1_sub_use ' , n_S_1_sub_use ,
     $        ' nS_1_per_min ' , nS_1_per_min ,' nS_1_per_max ' ,
     $        nS_1_per_max
      end if !if (verbose.gt.1) then
      allocate(n_M_1_per_(0:nM_0_per_max-1))
      allocate(n_M_1_sum_(0:nM_0_per_max-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*2.0d0*nM_0_per_max
      n_M_1_sub_use = min(n_M_1_sub__in,nM_0_per_max)
      call block_0(verbose-2,n_M_1_sub_use,nM_0_per_max,n_M_1_per_
     $     ,n_M_1_sum_,n_M_1_sub_use,nM_1_per_min,nM_1_per_max)
      if (verbose.gt.-1) then
         write(6,'(A,I0)') ' n_M ' , n_M
         write(6,'(3(A,I0))') ' n_M_0_sub_use ' , n_M_0_sub_use ,
     $        ' nM_0_per_min ' , nM_0_per_min ,' nM_0_per_max ' ,
     $        nM_0_per_max
         write(6,'(3(A,I0))') ' n_M_1_sub_use ' , n_M_1_sub_use ,
     $        ' nM_1_per_min ' , nM_1_per_min ,' nM_1_per_max ' ,
     $        nM_1_per_max
      end if !if (verbose.gt.1) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and T_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_T_T_R_CTF_M_q__(0:n_transf*nS_1_per_max*nM_1_per_max
     $        *n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_transf*1.0d0
     $     *nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
      end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_T_R_CTF_M_q__ and S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_T_T_R_CTF_M_q__(0:nS_1_per_max*n_transf*nM_1_per_max
     $        *n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max*1.0d0
     $     *n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max
      end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array Z_S_svdd_ to store '
         write(6,'(A)') ' displacement-operator for the delta-side '
         write(6,'(A)') ' of the svd-expansion applied to S. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_svd_l*1.0d0*n_delta_x*1.0d0*n_delta_y
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' Z_S_svdd_ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(Z_S_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_svd_l*1.0d0
     $     *n_delta_x*1.0d0*n_delta_y
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and Z_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_Z_T_R_CTF_M_q__(0:n_transf*nS_1_per_max*nM_1_per_max
     $        *n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_transf*1.0d0
     $     *nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
      if (verbose.gt.1) then
         write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_S_svdd_. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_delta_x*n_delta_y*1.0d0*nS_1_per_max*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_T_T_R_CTF_M_q__(0:n_delta_x*n_delta_y*nS_1_per_max
     $        *nM_1_per_max*n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_x*1.0d0
     $     *n_delta_y*1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0
     $     *n_w_max
      end if ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array Z_M_svdd_ to store '
         write(6,'(A)') ' displacement-operator for the delta-side '
         write(6,'(A)') ' of the svd-expansion applied to M. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_svd_l*1.0d0*n_delta_x*1.0d0*n_delta_y
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' Z_M_svdd_ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(Z_M_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      d_memory_estimate = d_memory_estimate + 16.0d0*n_svd_l*1.0d0
     $     *n_delta_x*1.0d0*n_delta_y
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' Z_T_R_CTF_M_q__ and S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_Z_T_R_CTF_M_q__(0:nS_1_per_max*n_transf*nM_1_per_max
     $        *n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max*1.0d0
     $     *n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max
      if (verbose.gt.1) then
         write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_M_svdd_. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_delta_x*n_delta_y*1.0d0*nS_1_per_max*1.0d0
     $     *nM_1_per_max*1.0d0*n_w_max*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      if (flag_memory_estimate.eqv..false.) then
         allocate(S_T_T_R_CTF_M_q__(0:n_delta_x*n_delta_y*nS_1_per_max
     $        *nM_1_per_max*n_w_max-1))
      end if !if (flag_memory_estimate.eqv..false.) then
      d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_x*1.0d0
     $     *n_delta_y*1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0
     $     *n_w_max
      end if ! if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

      allocate(n_S_9_per_(0:n_omp_sub__in-1))
      allocate(n_S_9_sum_(0:n_omp_sub__in-1))
      allocate(n_M_9_per_(0:n_omp_sub__in-1))
      allocate(n_M_9_sum_(0:n_omp_sub__in-1))
      d_memory_estimate = d_memory_estimate + 4.0d0*4.0d0*n_omp_sub__in

      if (flag_memory_estimate.eqv..false.) then
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $           nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $           nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         end if !if (verbose.gt.0 .or. verbose_mem.gt.0) then
      end if !if (flag_memory_estimate.eqv..false.) then

      if (flag_memory_estimate.eqv..true.) then
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $           nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $           nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         end if !if (verbose.gt.1 .or. verbose_mem.gt.0) then
         goto 10
      end if !if (flag_memory_estimate.eqv..true.) then

