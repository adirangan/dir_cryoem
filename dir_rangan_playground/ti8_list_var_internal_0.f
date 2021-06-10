!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> List of variables used internally. ;\n
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      function outputs
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$     $ ,max_i4_f,sum_i4_f 
c$$$     $ ,max_r8_f
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     $ ,flag_memory_checkset !temporary: used to check memory map. ;
     $ ,flag_tesselation
     $ ,flag_tesselation_ref
     $ ,n_point_init ! number of points used to initialize S_L_. ;
     $ ,n_point,npoint ! number of points used to define tesselation tree. ;
     $ ,S_L_ ! list of n_point points used to define tesselation tree. ;
     $ ,tradius_min ! minimum radius allowed for tesselation element. ;
     $ ,nl_max_init,nm_sum_init,ll_sum_init ! initial variables defining tesselation tree structure. ;
     $ ,nl_max,nm_sum,ll_sum ! variables defining tesselation tree structure. ;
     $ ,n_T_root_base,nT_root_base ! variables defining base of tesselation tree. ;
     $ ,parity_input ! temporary: flag used to determine orientation of tesselation element. ;
c$$$      T_ tesselation tree
     $ ,T_nl_ !level
     $ ,T_vm_ !vertex center
     $ ,T_tr_ !tradius
     $ ,T_ll_ !number of points from S_L_ in T_
     $ ,T_lf_ !is leaf
     $ ,T_c0_ !child_0 tesselation_index
     $ ,T_c1_ !child_1 tesselation_index
     $ ,T_c2_ !child_2 tesselation_index
     $ ,T_c3_ !child_3 tesselation_index
     $ ,T_ls_ !starting index of point_index_list for T_ if leaf (leaves only)
     $ ,T_LT_ !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
     $ ,T_root_base_ !T_ roots
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      bookkeeping tesselation search across threads
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     $ ,tesselation_distance_req_use !actual distance used as lower-bound for size of tesselation elements. ;
     $ ,n_SM_use_local_ ! array storing total number of image-template comparisons made for a particular image. ;
     $ ,n_SM_use ! total number of image-template comparisons made. ;
     $ ,vp_input_omp__ ! array of points on sphere (3-dimensional) (one for each thread). ;
     $ ,p_vp_input ! memory offset for single thread. ;
     $ ,flag_S_use_omp__ ! array of flags indicating which templates have been used (one for each thread). ;
     $ ,p_flag_S_use ! memory offset for single thread. ;
     $ ,n_LT,n_S_use,n_S_use_sum
     $ ,n_S_use_sum_ ! array of n_S_use (one for each thread). ;
     $ ,n_add,n_ref,nref,nLT,ns_local
     $ ,n_pass,n_pass_ref
     $ ,continue_flag
     $ ,LT_omp__ ! array of indices for nearest neighbor list (one for each thread). ;
     $ ,p_LT ! memory offset for single thread. ;
     $ ,S_alpha_S_index_local_omp__ ! template index. (one for each thread). ;
     $ ,p_S_alpha_S_index_local ! memory offset for single thread. ;
     $ ,S_alpha_polar_a_local_omp__ ! polar angle for templates (3-dimensional) (one for each thread). ;
     $ ,p_S_alpha_polar_a_local ! memory offset for single thread. ;
     $ ,S_alpha_azimu_b_local_omp__ ! azimuthal angle for templates (one for each thread). ;
     $ ,p_S_alpha_azimu_b_local ! memory offset for single thread. ;
     $ ,I_S_sample_local_omp__ ! array of template indices for use with S_p__ (one for each thread). ;
     $ ,p_I_S_sample_local ! memory offset for single thread. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
     $ ,eps_target ! copy of eps_svd. ;
     $ ,n_svd_r,n_svd_d,n_svd_l
     $ ,nsvd_r,nsvd_d,nsvd_l
     $ ,svd_r_
     $ ,svd_d_
     $ ,svd_l_
     $ ,svd_U_d_
     $ ,svd_s_
     $ ,svd_V_r_
     $ ,n_svd_max
     $ ,svd_d_max,svd_r_max ! variables used to set ranges for svd-expansion of translation. ;
     $ ,svd_polyval_U_d_ ! array storing polynomial evaluations of svd_U_d_ at the various delta values. ;
     $ ,svd_polyval_V_r_ ! array storing polynomial evaluations of svd_V_r_ at the various grid_k_p_(nr) values. ;
     $ ,flag_warning
     $ ,R_max,K_max,delta,delta_max,n_pixels ! variables used to select svd library. ;
     $ ,Z_S_svdd_ !array storing displacement-operator for the delta-side of the svd-expansion applied to S. ;
     $ ,Z_M_svdd_ !array storing displacement-operator for the delta-side of the svd-expansion applied to M. ;
     $ ,CTF_R_S__ !array to hold innerproducts for each CTF-S pair. ;
     $ ,CTF_R_S_sub__ !array to hold innerproducts for each CTF-S pair. ;
     $ ,CTF_R_S_sub_omp__ !array to hold innerproducts for each CTF-S pair. ;
     $ ,CTF_R_S_local_omp__ !similar to CTF_R_S__, allocated only if we intend to run local search. ;
     $ ,p_CTF_R_S_local ! memory offset for single thread. ;
     $ ,O_S_q__ !array to hold collection of S in bessel-coefficients. ;
     $ ,O_S_q_local_omp__ !similar to O_S_q__, allocated only if we intend to run local search. ;
     $ ,p_O_S_q_local ! memory offset for single thread. ;
     $ ,T_S_q__ !array to hold bessel-coefficients of transformed T_{update}(S), where T represents translation by delta_{update} ;
     $ ,T_S_q_local_omp__ !similar to T_S_q__, allocated only if we intend to run local search. ;
     $ ,p_T_S_q_local ! memory offset for single thread. ;
     $ ,Z_S_q__ !array to hold bessel-coefficients of the transformed Z_{rho}(S), where Z represents the rho-side of the svd-expansion applied to S. ;
     $ ,Z_S_q_local_omp__ !similar to Z_S_q__, allocated only if we intend to run local search. ;
     $ ,p_Z_S_q_local ! memory offset for single thread. ;
c$$$      pointer (p_S_p__,S_p__) !another name for S_k_p__. ;
     $ ,S_p__
c$$$      pointer (p_M_p__,M_p__) !another name for M_k_p__. ;
     $ ,M_p__
c$$$      pointer (p_CTF_p__,CTF_p__) !another name for CTF_k_p__. ;
     $ ,CTF_p__
     $ ,S_p_ !temporary storage for a single template. ;
     $ ,S_q_ !temporary storage for a single template. ;
     $ ,S_p_omp__ !temporary storage for a few templates across multiple threads. ;
     $ ,S_q_omp__ !temporary storage for a few templates across multiple threads. ;
     $ ,M_p_ !temporary storage for a single image. ;
     $ ,M_q_ !temporary storage for a single image. ;
     $ ,M_p_omp__ !temporary storage for a few images across multiple threads. ;
     $ ,M_q_omp__ !temporary storage for a few images across multiple threads. ;
     $ ,CTF_p_ !temporary storage for a single CTF. ;
     $ ,CTF_q_ !temporary storage for a single CTF. ;
     $ ,CTF_p_omp__ !temporary storage for a few CTFs across multiple threads. ;
     $ ,CTF_q_omp__ !temporary storage for a few CTFs across multiple threads. ;
     $ ,Z_p_ !temporary storage for a single array. ;
     $ ,Z_q_ !temporary storage for a single array. ;
     $ ,Z_p_omp__ !temporary storage for a few arrays across multiple threads. ;
     $ ,Z_q_omp__ !temporary storage for a few arrays across multiple threads. ;
     $ ,ZZ__ !temporary storage for a single array. ;
     $ ,ZZ_omp__ !temporary storage for a few arrays across multiple threads. ;
     $ ,ZZ_sub_ !temporary storage for a single array. ;
     $ ,ZZ_sub_omp__ !temporary storage for a few arrays across multiple threads. ;
     $ ,C_trn0_ !temporary storage for a single array. ;
     $ ,C_trn0_omp__ !temporary storage for a few arrays across multiple threads. ;
     $ ,O_T_R_CTF_M_q__ !array to hold bessel-coefficients of the transformed T_est(R_{-gamma_est}(CTF.*M)), where T_est = translation by -delta_est and R = rotation by -gamma_est. ;
     $ ,O_T_R_CTF_M_q_local_omp__ !similar to O_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
     $ ,p_O_T_R_CTF_M_q_local ! memory offset for single thread. ;
     $ ,T_T_R_CTF_M_q__ !array to hold bessel-coefficients of the transformed T_upd(T_est(R_{-gamma_est}(CTF.*M))), where T_upd = translation by -delta_{update} and T_est = translation by -delta_est and R = rotation by -gamma_est. ;
     $ ,T_T_R_CTF_M_q_local_omp__ !similar to T_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
     $ ,p_T_T_R_CTF_M_q_local ! memory offset for single thread. ;
     $ ,Z_T_R_CTF_M_q__ !array to hold bessel-coefficients of the transformed Z_{rho}(T_est(R_{-gamma_est}(CTF.*M))), where Z_{rho} = rho-side of svd-expansion applied to M and T_est = translation by -delta_est and R = rotation by -gamma_est. ;
     $ ,Z_T_R_CTF_M_q_local_omp__ !similar to Z_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
     $ ,p_Z_T_R_CTF_M_q_local ! memory offset for single thread. ;
     $ ,S_T_T_R_CTF_M_q__ !array to hold bessel-coefficients for the product of O_T_R_CTF_M_q__ and T_S_q__ (or T_T_R_CTF_M_q__ and O_S_q__) integrated over k (i.e., nr). ;
     $ ,S_T_T_R_CTF_M_q_local_omp__ !similar to S_T_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
     $ ,p_S_T_T_R_CTF_M_q_local ! memory offset for single thread. ;
     $ ,S_Z_T_R_CTF_M_q__ !array to hold bessel-coefficients for the product of Z_R_CTF_M_q__ and T_S_q__ (or Z_T_R_CTF_M_q__ and O_S_q__) integrated over k (i.e., nr). ;
     $ ,S_Z_T_R_CTF_M_q_local_omp__ !similar to S_Z_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
     $ ,p_S_Z_T_R_CTF_M_q_local ! memory offset for single thread. ;
     $ ,n_C_trn0 !length of corresponding array immediately after allocation. ;
     $ ,n_C_trn0_omp !length of corresponding array immediately after allocation. ;
     $ ,n_CTF_R_S_local_omp !length of corresponding array immediately after allocation. ;
     $ ,n_S_alpha_S_index_local_omp !length of corresponding array immediately after allocation. ;
     $ ,n_O_S_q_local_omp !length of corresponding array immediately after allocation. ;
     $ ,n_T_S_q_local_omp !length of corresponding array immediately after allocation. ;
     $ ,n_Z_S_q_local_omp !length of corresponding array immediately after allocation. ;      
     $ ,n_O_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
     $ ,n_O_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
     $ ,n_T_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
     $ ,n_T_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
     $ ,n_Z_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
     $ ,n_Z_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
     $ ,n_S_T_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
     $ ,n_S_T_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
     $ ,n_S_Z_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
     $ ,n_S_Z_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
c$$$      indices
c$$$     $ ,nr,na,ns,nctf,nm,nw
     $ ,nr,na,ns,nctf,nm,nw
     $ ,n_r !integer: copy of n_k_cur
     $ ,n_w_ !array (length at least n_k_cur) to store number of angles (omega) for each n_k. ;
     $ ,n_w_csum_ !array (length at least n_k_cur) to store cumulative sum of n_w_. ;
     $ ,n_w_max !integer: largest value in n_w_. ;
     $ ,n_A !integer storing the sum of n_w_. ;
c$$$      temporary: image-parameters for each image
     $ ,polar_a_est_
     $ ,azimu_b_est_
     $ ,gamma_z_est_
     $ ,delta_x_est_
     $ ,delta_y_est_
     $ ,l2_norm_est_
     $ ,ctf_ind_est_
     $ ,S_index_est_
     $ ,M_index_est_
     $ ,alpha_0in_
     $ ,alpha_0in_omp__
c$$$      fftw for omp
     $ ,fftw_plan_frwd__
     $ ,fftw_plan_back__
     $ ,fftw_0in__
     $ ,fftw_out__
c$$$      fftw for local use
     $ ,fftw_plan_frwd_
     $ ,fftw_plan_back_
     $ ,fftw_0in_
     $ ,fftw_out_
c$$$      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
c$$$      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
c$$$      pointer (p_fftw_0in_last_,fftw_0in_last_)
c$$$      pointer (p_fftw_out_last_,fftw_out_last_)
     $ ,fftw_plan_frwd_last,fftw_plan_back_last
     $ ,fftw_0in_last_,fftw_out_last_
c$$$      fftw plan many (fpm) for omp
     $ ,n_transf !number of entries in transformed vector. this could be either n_delta_v (if svd_calculation_type.eq.2) or n_svd_l (if svd_calculation_type.eq.1). ;
     $ ,fpm_howmany
     $ ,n_fpm
c$$$     $ ,nfpm
     $ ,nfpm
     $ ,fpm_frwd_
     $ ,fpm_back_
     $ ,fpm_0in__
     $ ,fpm_out__
     $ ,fpm_rank
     $ ,fpm_n__
     $ ,fpm_inembed__
     $ ,fpm_istride
     $ ,fpm_idist
     $ ,fpm_onembed__
     $ ,fpm_ostride
     $ ,fpm_odist
c$$$      array of displacements and rotations to measure
     $ ,ndv,ngz,ndv_optimal,ngz_optimal
     $ ,delta_x,delta_y,gamma_z
     $ ,gamma_z_ !temporary: array of length n_gamma_z holds values for gamma_z_update. ;
c$$$      parameters for blocking S and M
     $ ,nx0,nx1,nx2,nx3,nx4
     $ ,nm0,nm1,nm2,nm3,nm4
     $ ,ns0,ns1,ns2,ns3,ns4
     $ ,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp,n_S_tmp
     $ ,nomp_sub_use,nomp_sub,nomp_per,nomp_sum
     $ ,n_S_9_sub_use
c$$$     $ ,nS_9_sub,nS_9_per,nS_9_sum
     $ ,nS_9_sub,nS_9_per,nS_9_sum
     $ ,n_S_9_per_
     $ ,n_S_9_sum_
     $ ,nS_9_per_max,nS_9_per_min
     $ ,n_M_9_sub_use
c$$$     $ ,nM_9_sub,nM_9_per,nM_9_sum
     $ ,nM_9_sub,nM_9_per,nM_9_sum
     $ ,n_M_9_per_
     $ ,n_M_9_sum_
     $ ,nM_9_per_max,nM_9_per_min
     $ ,n_S_0_sub_use
c$$$      We pass in nS_0_sub, nS_0_per and nS_0_sum
     $ ,nS_0_sub,nS_0_per,nS_0_sum
     $ ,n_S_0_per_
     $ ,n_S_0_sum_
     $ ,nS_0_per_max,nS_0_per_min
     $ ,n_M_0_sub_use
c$$$     $ ,nM_0_sub,nM_0_per,nM_0_sum
     $ ,nM_0_sub,nM_0_per,nM_0_sum
     $ ,n_M_0_per_
     $ ,n_M_0_sum_
     $ ,nM_0_per_max,nM_0_per_min
     $ ,n_S_1_sub_use
c$$$     $ ,nS_1_sub,nS_1_per,nS_1_sum
     $ ,nS_1_sub,nS_1_per,nS_1_sum
     $ ,n_S_1_per_
     $ ,n_S_1_sum_
     $ ,nS_1_per_max,nS_1_per_min
     $ ,n_M_1_sub_use
c$$$     $ ,nM_1_sub,nM_1_per,nM_1_sum
     $ ,nM_1_sub,nM_1_per,nM_1_sum
     $ ,n_M_1_per_
     $ ,n_M_1_sum_
     $ ,nM_1_per_max,nM_1_per_min
c$$$      parameters for memory map
     $ ,d_mem
     $ ,verbose_mem,verbose_timing
c$$$      parameters for timing and testing
c$$$      Note that timing summation is not thread-safe
     $ ,timing_tic,timing_toc,timing_tot,timing_tmp,gnump_tot
     $ ,flag_fill !logical: indicates whether or not to test fill-only fftw fpm. ;
     $ ,flag_test !logical: indicates whether or not to test innerproducts against single image-template calculation. ;
     $ ,timing_total_fftw_plan,gnump_total_fftw_plan
     $ ,timing_total_CTF_R_S,gnump_total_CTF_R_S
     $ ,timing_total_O_S_q,gnump_total_O_S_q
     $ ,timing_total_T_S_q,gnump_total_T_S_q
     $ ,timing_total_Z_S_q,gnump_total_Z_S_q
     $ ,timing_total_O_T_R_CTF_M_q,gnump_total_O_T_R_CTF_M_q
     $ ,timing_total_T_T_R_CTF_M_q,gnump_total_T_T_R_CTF_M_q
     $ ,timing_total_Z_T_R_CTF_M_q,gnump_total_Z_T_R_CTF_M_q
     $ ,timing_total_zgemm,gnump_total_zgemm
     $ ,timing_total_fpm_fill,gnump_total_fpm_fill
     $ ,timing_total_fpm_fftw,gnump_total_fpm_fftw
     $ ,timing_total_transpose,gnump_total_transpose
     $ ,timing_total_Zstore,gnump_total_Zstore
     $ ,flag_time_Zstore
     $ ,timing_total_Zstore_a,gnump_total_Zstore_a
     $ ,timing_total_Zstore_b,gnump_total_Zstore_b
     $ ,timing_total_Zstore_c,gnump_total_Zstore_c
     $ ,timing_total_Zstore_x,gnump_total_Zstore_x
     $ ,timing_total_Zstore_d,gnump_total_Zstore_d
     $ ,timing_total_Zstore_e,gnump_total_Zstore_e
     $ ,timing_string
c$$$c$$$      format_string for printing output
c$$$     $ ,format_string
c$$$     $ ,flag_MDA
c$$$     $ ,MDA_n_d
c$$$     $ ,MDA_d_(0:64-1)
c$$$     $ ,MDA_string
c$$$c$$$      pi
c$$$     $ ,pi
