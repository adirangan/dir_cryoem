      call get_template_size !Calculating template size using 'get_template_size'. ;
      define n_w_, n_A, n_w_max !Calculate indices. ;
      call get_delta_0 call get_gamma_0 !set array of displacements and rotations. ;
      if (svd_calculation_type.eq.1) call get_svd_0 !select svd library. ;
      if (svd_calculation_type.eq.1) n_transf = n_svd_l !define n_transf. ;
      if (svd_calculation_type.eq.2) n_transf = n_delta_x*n_delta_y !define n_transf. ;
      if (svd_calculation_type.eq.1) 
         call innerproduct_q_k_svdd_1(flag_RTRT_vs_RTTR) !define Z_S_svdd_ or Z_M_svdd_. ;
      call dfftw_plan_dft_1d_ !Generate fftw_plans_for local use (size n_A). ;
      call dfftw_plan_dft_1d_ !Generate fftw_plans for omp sub-blocks (size n_A*n_omp_sub). ;
      call dfftw_plan_many_dft !Generate fftw_plans for fpm (size n_w_max*fpm_howmany). ;
      define n_S_0_use, n_S_0_per, I_S_0_per !level-0 block n_S. ;
      do ns_0=0,n_S_0_use-1 ! run through level-0 blocks for n_S ;
         redefine n_S_9_use, n_S_9_per, I_S_9_per !define level-9 block n_S using nS_0_per. ;
         call test_innerproduct_fast_CTF_R_S_2 !Calculate CTF_R_S_ <-- loop across level-9 blocks for omp ;
         if (flag_RTRT_vs_RTTR.eqv..false.) 
         %   call test_innerproduct_fast_S_q_1 !Calculate S_q__. <-- loop across level-9 blocks for omp ;
         if (flag_RTRT_vs_RTTR.eqv..true .and. svd_calculation_type.eq.1)
         %   call test_innerproduct_fast_Z_S_q_3 !Calculate Z_S_q__. <-- loop across level-9 blocks for omp ;
         if (flag_RTRT_vs_RTTR.eqv..true .and. svd_calculation_type.eq.2)
         %   call test_innerproduct_fast_T_S_q_3 !Calculate T_S_q__. <-- loop across level-9 blocks for omp ;
         redefine n_M_0_use, n_M_0_per, I_M_0_per !level-0 block n_M. ;
         do nm_0=0,n_M_0_use-1 ! run through level-0 blocks for n_M ;
            redefine n_M_9_use, n_M_9_per, I_M_9_per !define level-9 block n_M using nM_0_per. ;
            redefine gamma_z_est, delta_x_est, delta_y_est, etc. from alpha_est.
            if (flag_RTRT_vs_RTTR.eqv..true.)
            %   call test_innerproduct_fast_T_R_CTF_M_q_3 !Calculate T_R_CTF_M_q__. <-- loop level-9 omp ;
            if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.2)
            %   call test_innerproduct_fast_T_T_R_CTF_M_q_1 !Calculate T_T_R_CTF_M_q__. <-- loop level-9 omp ;
            if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.1)
            %   call test_innerproduct_fast_Z_T_R_CTF_M_q_1 !Calculate Z_T_R_CTF_M_q__. <-- loop level-9 omp ;
            redefine n_M_1_use, n_M_1_per, I_M_1_per !define level-1 block n_M using nM_0_per. ;
            do nm_1=0,n_M_1_use-1 ! run through level-1 blocks for n_M ;
               redefine n_S_9_use, n_S_9_per, I_S_9_per !define level-9 block n_S using nS_0_per. ;
               redefine n_M_9_use, n_M_9_per, I_M_9_per !define level-9 block n_M using nM_1_per. ;
               if (flag_RTRT_vs_RTTR.eqv..true. .and. svd_calculation_type.eq.2)
               %   define S_T_T_R_CTF_M_q__ !zgemm T_S_q__ and T_R_CTF_M_q__. <-- loop level-9 omp ;
               %   call dfftw_block_many_1 !fftw.  <-- loop level-9 omp ;
               %   call test_innerproduct_fast_STxTRM_1 !check S_T_T_R_CTF_M_q__. ;
               if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.2)
               %   define S_T_T_R_CTF_M_q__ !zgemm S_q__ and T_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call dfftw_block_many_1 !fftw.  <-- loop level-9 omp ;
               %   call trn0_c16 !transpose dimensions (1,2) of S_T_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call test_innerproduct_fast_SxTTRM_1 !check S_T_T_R_CTF_M_q__. ;
               if (flag_RTRT_vs_RTTR.eqv..true. .and. svd_calculation_type.eq.1)
               %   define S_Z_T_R_CTF_M_q__ !zgemm Z_S_q__ and T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call dfftw_block_many_1 !fftw.  <-- loop level-9 omp ;
               %   define S_T_T_R_CTF_M_q__ !zgemm Z_S_svdd_ and S_Z_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call test_innerproduct_fast_SZxTRM_1 !check S_T_T_R_CTF_M_q__. ;
               if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.1)
               %   define S_Z_T_R_CTF_M_q__ !zgemm S_q__ and Z_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call dfftw_block_many_1 !fftw.  <-- loop level-9 omp ;
               %   call trn0_c16 !transpose dimensions (1,2) of S_Z_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   define S_T_T_R_CTF_M_q__ !zgemm Z_M_svdd_ and S_Z_T_R_CTF_M_q__.  <-- loop level-9 omp ;
               %   call test_innerproduct_fast_SxZTRM_1 !check S_T_T_R_CTF_M_q__. ;
               call test_innerproduct_fast_Zstore_1 !update n_SM_ and alpha_SM__. ;
            enddo !do nm_1=0,n_M_1_use-1 ! run through level-1 blocks for n_M ;
         enddo !do nm_0=0,n_M_0_use-1 ! run through level-0 blocks for n_M ;
      enddo !do ns_0=0,n_S_0_use-1 ! run through level-0 blocks for n_S ;
      deallocate
      
      

      


