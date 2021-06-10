      call get_template_size !Calculating template size using 'get_template_size'. ;
      define n_w_, n_A, n_w_max !Calculate indices. ;
      call get_delta_0 call get_gamma_0 !set array of displacements and rotations. ;
      if (svd_calculation_type.eq.1) call get_svd_0 !select svd library. ;
      if (svd_calculation_type.eq.1) n_transf = n_svd_l !define n_transf. ;
      if (svd_calculation_type.eq.2) n_transf = n_delta_x*n_delta_y !define n_transf. ;
      if (svd_calculation_type.eq.1) 
         call innerproduct_q_k_svdd_1(flag_RTRT_vs_RTTR) !define Z_S_svdd_ or Z_M_svdd_. ;
      call dfftw_plan_dft_1d_ !Generate fftw_plans_for local use (size n_A). ;
      call test_innerproduct_fast_CTF_R_S_1 !Calculate CTF_R_S_ for each ctf-S pair. ;
      define n_S_sub_use, n_S_per, I_S_per !block n_S. ;
      call dfftw_plan_dft_1d_ !Generate fftw_plans for omp sub-blocks (size n_A*n_omp_sub). ;
      if (flag_RTRT_vs_RTTR.eqv..false.) 
         call test_innerproduct_fast_S_q_0 !Calculate S_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..true .and. svd_calculation_type.eq.1)
         call test_innerproduct_fast_Z_S_q_2 !Calculate Z_S_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..true .and. svd_calculation_type.eq.2)
         call test_innerproduct_fast_T_S_q_2 !Calculate T_S_q__. ;
      define n_M_sub_use, n_M_per, I_M_per !block n_M. ;
      define gamma_z_est, delta_x_est, delta_y_est, etc. from alpha_est.
      if (flag_RTRT_vs_RTTR.eqv..true.)
         call test_innerproduct_fast_T_R_CTF_M_q_2 !Calculate T_R_CTF_M_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.2)
         call test_innerproduct_fast_T_T_R_CTF_M_q_0 !Calculate T_T_R_CTF_M_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.1)
         call test_innerproduct_fast_Z_T_R_CTF_M_q_0 !Calculate Z_T_R_CTF_M_q__. ;
      define gamma_z_upd, delta_x_upd, delta_y_upd, etc.
      call dfftw_plan_many_dft !Generate fftw_plans for fpm (size n_w_max*fpm_howmany). ;
      if (flag_RTRT_vs_RTTR.eqv..true. .and. svd_calculation_type.eq.2)
         define S_T_T_R_CTF_M_q__ !zgemm T_S_q__ and T_R_CTF_M_q__. ;
         call dfftw_block_many_1 !fftw. ;
         call test_innerproduct_fast_STxTRM_1 !check S_T_T_R_CTF_M_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.2)
         define S_T_T_R_CTF_M_q__ !zgemm S_q__ and T_T_R_CTF_M_q__. ;
         call dfftw_block_many_1 !fftw. ;
         call trn0_c16 !transpose dimensions (1,2) of S_T_T_R_CTF_M_q__. ;
         call test_innerproduct_fast_SxTTRM_1 !check S_T_T_R_CTF_M_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..true. .and. svd_calculation_type.eq.1)
         define S_Z_T_R_CTF_M_q__ !zgemm Z_S_q__ and T_R_CTF_M_q__. ;
         call dfftw_block_many_1 !fftw. ;
         define S_T_T_R_CTF_M_q__ !zgemm Z_S_svdd_ and S_Z_T_R_CTF_M_q__. ;
         call test_innerproduct_fast_SZxTRM_1 !check S_T_T_R_CTF_M_q__. ;
      if (flag_RTRT_vs_RTTR.eqv..false. .and. svd_calculation_type.eq.1)
         define S_Z_T_R_CTF_M_q__ !zgemm S_q__ and Z_T_R_CTF_M_q__. ;
         call dfftw_block_many_1 !fftw. ;
         call trn0_c16 !transpose dimensions (1,2) of S_Z_T_R_CTF_M_q__. ;
         define S_T_T_R_CTF_M_q__ !zgemm Z_M_svdd_ and S_Z_T_R_CTF_M_q__. ;
         call test_innerproduct_fast_SxZTRM_1 !check S_T_T_R_CTF_M_q__. ;
      call test_innerproduct_fast_Zstore_0 !update n_SM_ and alpha_SM__. ;
      deallocate
      
      

      


