c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_fast_stage_1(verbose,n_omp_sub
     $     ,svd_calculation_type,eps_svd,n_r,grid_p_,n_w_
     $     ,half_diameter_x_c,n_S ,ld_S,S_p__ ,tesselation_distance_req
     $     ,S_alpha_polar_a_ ,S_alpha_azimu_b_,L_ ,nl_max ,nm_sum,ll_sum
     $     ,T_nl_,T_vm_ ,T_tr_,T_ll_,T_lf_ ,T_c0_ ,T_c1_ ,T_c2_,T_c3_
     $     ,T_ls_,T_LT_ ,n_T_root_base,T_root_base_ ,n_M ,ld_M,M_p__
     $     ,n_CTF ,ld_CTF ,CTF_p__,polar_a_est_ ,azimu_b_est_
     $     ,gamma_z_est_ ,delta_x_est_,delta_y_est_ ,l2_norm_est_
     $     ,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_
     $     ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $     ,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f ,flag_MS_vs_SM
     $     ,N_pixels_in ,displacement_max ,n_delta_x ,n_delta_y
     $     ,n_gamma_z , fpm_howmany_max, C_M_ ,n_SM_use)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer n_omp_sub ! number of batches to use when finding innerproducts
      integer svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,ld_S,n_M,ld_M,n_CTF,ld_CTF
      real *8 grid_p_(0:n_r-1),half_diameter_x_c
      complex *16 S_p__(0:ld_S*n_S-1)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      real *8 L_(0:0)
      real *8 tradius_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
c$$$      T_ structure
      integer *4 T_nl_(0:0) !level
      real *8    T_vm_(0:0) !vertex center
      real *8    T_tr_(0:0) !tradius
      integer *4 T_ll_(0:0) !number of points from L_ in T_
      logical    T_lf_(0:0) !is leaf
      integer *4 T_c0_(0:0) !child_0 tesselation_index
      integer *4 T_c1_(0:0) !child_1 tesselation_index
      integer *4 T_c2_(0:0) !child_2 tesselation_index
      integer *4 T_c3_(0:0) !child_3 tesselation_index
      integer *4 T_ls_(0:0) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4 T_LT_(0:0) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4 T_root_base_(0:0) !T_ roots
      real *8 tesselation_distance_req
c$$$      integer *4, allocatable :: SM_use_sub_(:)
      integer *4 n_SM_use
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      complex *16 M_p__(0:ld_M*n_M-1)
      complex *16 CTF_p__(0:0)
      complex *16, allocatable :: CTF_p_(:)
      complex *16, allocatable :: CTF_q_(:)
      complex *16, allocatable :: CTF_p_omp__(:)
      complex *16, allocatable :: CTF_q_omp__(:)
      real *8 polar_a_est_(0:n_M-1),azimu_b_est_(0:n_M-1)
     $     ,gamma_z_est_(0:n_M-1)
      real *8 delta_x_est_(0:n_M-1),delta_y_est_(0:n_M-1)
     $     ,l2_norm_est_(0:n_M-1)
      real *8 ctf_ind_est_(0:n_M-1)
      real *8 S_index_est_(0:n_M-1)
      real *8 polar_a_upd_(0:n_M-1),azimu_b_upd_(0:n_M-1)
     $     ,gamma_z_upd_(0:n_M-1)
      real *8 delta_x_upd_(0:n_M-1),delta_y_upd_(0:n_M-1)
     $     ,l2_norm_upd_(0:n_M-1)
      real *8 ctf_ind_upd_(0:n_M-1)
      real *8 S_index_upd_(0:n_M-1)
      real *8 alpha_update_f
      logical flag_MS_vs_SM
      real *8 eps_svd,eps_target,N_pixels_in
      real *8 displacement_max
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      complex *16 C_M_(0:n_M-1),C_M
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      integer n_svd_max
      parameter (n_svd_max=512)
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels
      complex *16, allocatable :: Z_svdr_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: CTF_R_S_(:)
      complex *16, allocatable :: T_S_q__(:)
      complex *16, allocatable :: Z_S_q__(:)
      complex *16, allocatable :: S_p_(:)
      complex *16, allocatable :: S_q_(:)
      complex *16, allocatable :: S_p_omp__(:)
      complex *16, allocatable :: S_q_omp__(:)
      complex *16, allocatable :: M_p_(:)
      complex *16, allocatable :: M_q_(:)
      complex *16, allocatable :: M_p_omp__(:)
      complex *16, allocatable :: M_q_omp__(:)
      complex *16, allocatable :: T_R_CTF_M_q__(:)
      complex *16, allocatable :: S_T_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: S_Z_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: ZZ_(:)
      complex *16, allocatable :: ZZ_omp__(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16, allocatable :: C_Z_MS_(:)
      complex *16, allocatable :: C_S_MS_(:)
      real *8, allocatable :: delta_x_MS_(:)
      real *8, allocatable :: delta_y_MS_(:)
      real *8, allocatable :: gamma_z_MS_(:)
      complex *16, allocatable :: C_Z_sort_MS_(:)
      complex *16, allocatable :: C_S_sort_MS_(:)
      real *8, allocatable :: delta_x_sort_MS_(:)
      real *8, allocatable :: delta_y_sort_MS_(:)
      real *8, allocatable :: gamma_z_sort_MS_(:)
      complex *16, allocatable :: I_permute_MS_(:)
      complex *16, allocatable :: I_inverse_MS_(:)
      integer nx1,nx2,nx3
c$$$      fftw for omp
      integer *8, allocatable :: fftw_plan_frwd__(:)
      integer *8, allocatable :: fftw_plan_back__(:)
      complex *16, allocatable :: fftw_in1__(:)
      complex *16, allocatable :: fftw_out__(:)
c$$$      fftw for local use
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*),fftw_out_last_(*)
c$$$      fftw plan many (fpm) for omp
      integer *8, allocatable :: fpm_frwd_(:)
      integer *8, allocatable ::  fpm_back_(:)
      complex *16, allocatable :: fpm_in1__(:)
      complex *16, allocatable :: fpm_out__(:)
      integer fpm_rank
      integer, allocatable :: fpm_n__(:)
      integer fpm_howmany_max
      integer n_transf,fpm_howmany
      integer nfpm
      integer, allocatable :: fpm_inembed__(:)
      integer fpm_istride
      integer fpm_idist
      integer, allocatable :: fpm_onembed__(:)
      integer fpm_ostride
      integer fpm_odist
      integer n_fpm
c$$$      indices
      integer nr,n_w_max,n_A,na,ns,nctf,nm,nw
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      character(len=64) format_string
c$$$      parameters for timing and testing
      logical flag_test
      real *8 timing_tic,timing_toc,timing_tot,timing_tmp
c$$$      parameters for omp
      integer *4 nomp_sub,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp
      integer *4 n_M_sub_use,nM_sub,nM_per,IM_per
      integer *4, allocatable :: n_M_per_(:)
      integer *4, allocatable :: I_M_per_(:)
      integer *4 nM_per_max
      integer *4 n_S_sub_use,nS_sub,nS_per,IS_per
      integer *4, allocatable :: n_S_per_(:)
      integer *4, allocatable :: I_S_per_(:)
      integer *4 nS_per_max
c$$$      pi
      real *8 pi

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_stage_1]'
      end if      

      eps_target = eps_svd

      if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_omp_sub: ',n_omp_sub
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
         write(6,'(A,F6.3)') 'eps_svd: ',eps_svd
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,I0)') 'n_r: ',n_r
         write(format_string,'(A,I0,A)') '(A,',n_r,'(F8.4,1X))'
         write(6,format_string) 'grid_p_: ',(grid_p_(nr),nr=0 ,n_r-1)
         write(format_string,'(A,I0,A)') '(A,',n_r,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0 ,n_r-1)
         write(6,'(A,F6.3)') 'half_diameter_x_c: ',half_diameter_x_c
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,'(A)') 'S_alpha_polar_a_: '
         write(6,format_string) (S_alpha_polar_a_(nr),nr=0,n_S-1)
         write(6,'(A)') 'S_alpha_azimu_b_: '
         write(6,format_string) (S_alpha_azimu_b_(nr),nr=0,n_S-1)
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(6,'(A,I0)') 'n_CTF: ',n_CTF
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         write(6,'(A,F6.3)') 'n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') 'displacement_max: ',displacement_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'fpm_howmany_max: ',fpm_howmany_max
      end if                    !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Generating indices'
         write(6,'(A)') ' '
      end if
      pi = 4.0*atan(1.0)
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_A: ',n_A,'; n_w_max: ',n_w_max
         write(format_string,'(A,I0,A)') '(A,',n_r,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,n_r-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
         write(6,'(A)') ' '
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      call get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)') ' ncur or so.'
         write(6,'(A)') ' '
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(A,',n_delta_x,'(F5.3,1X))'
         write (6,format_string) 'delta_x_: ',(delta_x_(ndx),ndx=0
     $        ,n_delta_x-1)
         write(format_string,'(A,I0,A)') '(A,',n_delta_y,'(F5.3,1X))'
         write (6,format_string) 'delta_y_: ',(delta_y_(ndy),ndy=0
     $        ,n_delta_y-1)
         write(format_string,'(A,I0,A)') '(A,',n_gamma_z,'(F5.3,1X))'
         write (6,format_string) 'gamma_z_: ',(gamma_z_(ngz),ngz=0
     $        ,n_gamma_z-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if
      allocate(svd_r_(0:n_svd_max-1))
      allocate(svd_d_(0:n_svd_max-1))
      allocate(svd_l_(0:n_svd_max-1))
      allocate(svd_U_d_(0:n_svd_max*n_svd_max-1))
      allocate(svd_s_(0:n_svd_max-1))
      allocate(svd_V_r_(0:n_svd_max*n_svd_max-1))
      call get_svd_0(eps_target,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,grid_p_,n_r,n_delta_x,delta_x_,n_delta_y
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

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up displacement-operator Z_svdd_'
         write(6,'(A)') ' associated with svd-expansion.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      timing_tic = omp_get_wtime()
      call innerproduct_q_k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q_k_svdd: total time ',timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y
     $        *1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' innerproduct_q_k_svdd total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for local use'
         write(6,'(A)') ' '
      end if
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      timing_tic = omp_get_wtime()
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nr)
      enddo
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_1d:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array CTF_R_S_ to hold'
         write(6,'(A)') ' innerproducts for each CTF-S pair.'
         write(6,'(A)') ' '
      end if
      allocate(CTF_R_S_(0:n_gamma_z*n_S*n_CTF-1))
      call cl1_c16(n_gamma_z*n_S*n_CTF,CTF_R_S_)
      allocate(S_p_(0:n_A-1))
      allocate(S_q_(0:n_A-1))
      allocate(S_p_omp__(0:n_A*n_omp_sub-1))
      allocate(S_q_omp__(0:n_A*n_omp_sub-1))
      allocate(CTF_p_(0:n_A-1))
      allocate(CTF_q_(0:n_A-1))
      allocate(CTF_p_omp__(0:n_A*n_omp_sub-1))
      allocate(CTF_q_omp__(0:n_A*n_omp_sub-1))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate CTF_R_S_ for each ctf-S pair.'
         write(6,'(A)') ' CTF_R_S = || CTF .* R(S) ||.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)')
     $        ' CTF_R_S_(ngz + ns*n_gamma_z + nctf*n_gamma_z*n_S) '
         write(6,'(A)') ' is equal to the l2_norm (not squared)'
         write(6,'(A)') ' of the pointwise product of CTF '
         write(6,'(A)') ' and R(S), where CTF is the nctf-th '
         write(6,'(A)') ' ctf-function, S is the ns-th template, '
         write(6,'(A)') ' and R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_CTF_R_S_0(verbose-1,n_gamma_z
     $     ,gamma_z_,fftw_plan_frwd_,fftw_in1_,fftw_out_
     $     ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_,n_r
     $     ,grid_p_,n_w_,n_A,S_p_,S_q_,n_S ,ld_S,S_p__,CTF_p_,CTF_q_
     $     ,n_CTF,ld_CTF ,CTF_p__,CTF_R_S_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating CTF_R_S_ for each ctf-S pair. '
     $        ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)*(n_CTF*1.0d0)/timing_tot
     $        /1e9
         write(6,'(A,F8.4)') ' CTF_R_S_ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      allocate(n_S_per_(0:n_S-1))
      allocate(I_S_per_(0:n_S-1))
      n_S_sub_use = min(n_omp_sub,n_S)
      if (n_S_sub_use.eq.1) then
         write(6,'(A,A,I0)')
     $        ' Warning! Turning off omp for n_S ',
     $        ' n_S_sub_use: ',n_S_sub_use
      end if !if (n_S_sub_use.eq.1) then
      if (verbose.gt.0) then
         write (6,'(A,I0,A)') ' Setting up n_S_sub_use = ',n_S_sub_use
     $        ,' sub-blocks for omp parallelization'
      end if ! if (verbose.gt.0) then
      do nS_sub=0,n_S_sub_use-1
         if (nS_sub.eq.0) then
            n_S_per_(0) = n_S/n_S_sub_use
            I_S_per_(0) = 0
         else if (nS_sub.gt.0 .and. nS_sub.lt.n_S_sub_use-1) then
            n_S_per_(nS_sub) = n_S /n_S_sub_use
            I_S_per_(nS_sub) = I_S_per_(nS_sub-1) + n_S_per_(nS_sub-1)
         else if (nS_sub.eq.n_S_sub_use-1) then
            I_S_per_(n_S_sub_use-1) = I_S_per_(n_S_sub_use-2) +
     $           n_S_per_(n_S_sub_use-2)
            n_S_per_(n_S_sub_use-1) = n_S - I_S_per_(n_S_sub_use-1)
         end if ! if nS_sub
      enddo ! do nS_sub=0,n_S_sub_use-1
      nS_per_max = max_i4_f(n_S_sub_use,n_S_per_)
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(A,',n_S_sub_use ,'(I0,1X))'
         write(6,format_string) 'n_S_per_: ' ,(n_S_per_(nS_sub),nS_sub=0
     $        ,n_S_sub_use-1)
         write(6,format_string) 'I_S_per_: ' ,(I_S_per_(nS_sub),nS_sub=0
     $        ,n_S_sub_use-1)
         write(6,'(A,I0)') 'nS_per_max: ' , nS_per_max
      end if ! if (verbose.gt.1) then

      allocate(n_M_per_(0:n_M-1))
      allocate(I_M_per_(0:n_M-1))
      n_M_sub_use = min(n_omp_sub,n_M)
      if (n_M_sub_use.eq.1) then
         write(6,'(A,A,I0)')
     $        ' Warning! Turning off omp for n_M ',
     $        ' n_M_sub_use: ',n_M_sub_use
      end if !if (n_M_sub_use.eq.1) then
      if (verbose.gt.0) then
         write (6,'(A,I0,A)') ' Setting up n_M_sub_use = ',n_M_sub_use
     $        ,' sub-blocks for omp parallelization'
      end if ! if (verbose.gt.0) then
      do nM_sub=0,n_M_sub_use-1
         if (nM_sub.eq.0) then
            n_M_per_(0) = n_M/n_M_sub_use
            I_M_per_(0) = 0
         else if (nM_sub.gt.0 .and. nM_sub.lt.n_M_sub_use-1) then
            n_M_per_(nM_sub) = n_M /n_M_sub_use
            I_M_per_(nM_sub) = I_M_per_(nM_sub-1) + n_M_per_(nM_sub-1)
         else if (nM_sub.eq.n_M_sub_use-1) then
            I_M_per_(n_M_sub_use-1) = I_M_per_(n_M_sub_use-2) +
     $           n_M_per_(n_M_sub_use-2)
            n_M_per_(n_M_sub_use-1) = n_M - I_M_per_(n_M_sub_use-1)
         end if ! if nM_sub
      enddo ! do nM_sub=0,n_M_sub_use-1
      nM_per_max = max_i4_f(n_M_sub_use,n_M_per_)
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(A,',n_M_sub_use ,'(I0,1X))'
         write(6,format_string) 'n_M_per_: ' ,(n_M_per_(nM_sub),nM_sub=0
     $        ,n_M_sub_use-1)
         write(6,format_string) 'I_M_per_: ' ,(I_M_per_(nM_sub),nM_sub=0
     $        ,n_M_sub_use-1)
         write(6,'(A,I0)') 'nM_per_max: ' , nM_per_max
      end if ! if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
      end if ! if (verbose.gt.1) then
      allocate(fftw_plan_frwd__(0:n_r*n_omp_sub-1))
      allocate(fftw_plan_back__(0:n_r*n_omp_sub-1))
      allocate(fftw_in1__(0:n_A*n_omp_sub-1))
      allocate(fftw_out__(0:n_A*n_omp_sub-1))
      do nomp_sub=0,n_omp_sub-1
         n_r_tmp = nomp_sub*n_r
         n_A_tmp = nomp_sub*n_A
         na = 0
         do nr=0,n_r-1
            call dfftw_plan_dft_1d_(fftw_plan_frwd__(n_r_tmp+nr)
     $           ,n_w_(nr),fftw_in1__(n_A_tmp+na),fftw_out__(n_A_tmp
     $           +na),FFTW_FORWARD,FFTW_ESTIMATE) 
            call dfftw_plan_dft_1d_(fftw_plan_back__(n_r_tmp+nr)
     $           ,n_w_(nr),fftw_out__(n_A_tmp+na),fftw_in1__(n_A_tmp
     $           +na),FFTW_BACKWARD,FFTW_ESTIMATE) 
            na = na + n_w_(nr)
         enddo ! do nr=0,n_r-1
      enddo ! do nomp_sub=0,n_omp_sub-1

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_many_plan for omp sub-blocks'
         write(6,'(A)') ' '
      end if
      fpm_rank = 1
      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_x*n_delta_y
      else
         n_transf = n_svd_l
      end if !if (svd_calculation_type.eq.2) then
      fpm_howmany = min(fpm_howmany_max,n_transf*n_S*nM_per_max)
      allocate(fpm_n__(0:fpm_rank*n_omp_sub-1))
      do nomp_sub=0,n_omp_sub-1
         fpm_n__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub-1
      allocate(fpm_inembed__(0:fpm_rank*n_omp_sub-1))
      do nomp_sub=0,n_omp_sub-1
         fpm_inembed__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub-1
      allocate(fpm_onembed__(0:fpm_rank*n_omp_sub-1))
      do nomp_sub=0,n_omp_sub-1
         fpm_onembed__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub-1
      fpm_istride = 1
      fpm_idist = n_w_max
      fpm_ostride = 1
      fpm_odist = n_w_max
      n_fpm = fpm_howmany*fpm_n__(0)
      allocate(fpm_in1__(0:n_fpm*n_omp_sub-1))
      allocate(fpm_out__(0:n_fpm*n_omp_sub-1))     
      allocate(fpm_back_(0:n_omp_sub-1))
      timing_tic = omp_get_wtime()
c$$$c$$$      We do not actually use the forward plan
      do nomp_sub=0,n_omp_sub-1      
         call dfftw_plan_many_dft(fpm_back_(nomp_sub),fpm_rank
     $        ,fpm_n__(nomp_sub),fpm_howmany,fpm_in1__(n_fpm*nomp_sub)
     $        ,fpm_inembed__(nomp_sub) ,fpm_istride,fpm_idist
     $        ,fpm_out__(n_fpm*nomp_sub),fpm_onembed__(nomp_sub)
     $        ,fpm_ostride ,fpm_odist ,FFTW_BACKWARD,FFTW_MEASURE)
      enddo !do nomp_sub=0,n_omp_sub-1
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_many:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then

      if (svd_calculation_type.eq.2) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array T_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (T_{delta}(S))_q, '
            write(6,'(A)') ' where T represents translation by delta.'         
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         allocate(T_S_q__(0:n_r*n_transf*n_S*n_w_max))
         call cl1_c16(n_r*n_transf*n_S*n_w_max,T_S_q__)
         if (verbose.gt.1) then
            write(6,'(A)') ' Calculate bessel coefficients of T(S).'
            write(6,'(A)') ' T = translation by delta_x,delta_y. '
            write(6,'(A)') ' S = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A)') ' T_S_q__(nr + n_r*(ndx + '
     $           ,'n_delta_x*(ndy + n_delta_y*(ns + n_S*nw)))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' conjg(S_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
            write(6,'(A)') ' where S_q_ = T(S), with '
            write(6,'(A)') ' T <-- delta_x_(ndx) and delta_y_(ndy)'
            write(6,'(A)') ' and S <-- S_p__(ns*ld_S).'
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nS_per,IS_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
         do nS_sub=0,n_S_sub_use-1
            nS_per = n_S_per_(nS_sub)
            IS_per = I_S_per_(nS_sub)
            n_r_tmp = nS_sub*n_r
            n_A_tmp = nS_sub*n_A
            n_X_tmp = 0 + n_r*(0 + n_delta_x*(0 + n_delta_y*(IS_per +
     $           n_S*0)))
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' test_innerproduct_fast_T_S_q_1'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call test_innerproduct_fast_T_S_q_1(verbose-1,n_delta_x
     $           ,delta_x_,n_delta_y,delta_y_,fftw_plan_frwd__(n_r_tmp)
     $           ,fftw_in1__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r ,grid_p_
     $           ,n_w_,n_A,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_per
     $           ,ld_S ,S_p__(ld_S*IS_per),n_S ,T_S_q__(n_X_tmp))
         enddo !do nS_sub=0,n_S_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         if (verbose.gt.0) then
            write(6,'(A,F8.5)')
     $           ' finished calculating T_S_q__ for each S. '
     $           ,timing_tot
            timing_tmp = (n_A*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $           *(n_S*1.0d0)/timing_tot/1e9
            write(6,'(A,F8.4)') ' T_S_q__ total Gnumps: ' , timing_tmp
            write(6,'(A)') ' '
         end if ! if (verbose.gt.0) then
      end if !if (svd_calculation_type.eq.2) then

      if (svd_calculation_type.ne.2) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (Z_{delta}(S))_q, '
            write(6,'(A)') ' where Z represents the left hand side'
            write(6,'(A)') ' of the svd of the translation operator.'
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         allocate(Z_S_q__(0:n_r*n_transf*n_S*n_w_max))
         call cl1_c16(n_r*n_transf*n_S*n_w_max,Z_S_q__)
         if (verbose.gt.1) then
            write(6,'(A)') ' Calculate bessel coefficients of Z(S).'
            write(6,'(A)') ' Z = svd of translation operator. '
            write(6,'(A)') ' S = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A)') ' Z_S_q__(nr + n_r*(nl + '
     $           ,'n_svd_l*(ns + n_S*nw))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' conjg(S_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
            write(6,'(A)') ' where S_q_ = Z(S), with '
            write(6,'(A)') ' Z representing the left-hand factor '
            write(6,'(A)') ' of the translation-operator'
            write(6,'(A)') ' associated with delta,'
            write(6,'(A)') ' and S <-- S_p__(ns*ld_S).'
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nS_per,IS_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
         do nS_sub=0,n_S_sub_use-1
            nS_per = n_S_per_(nS_sub)
            IS_per = I_S_per_(nS_sub)
            n_r_tmp = nS_sub*n_r
            n_A_tmp = nS_sub*n_A
            n_X_tmp = 0 + n_r*(0 + n_svd_l*(IS_per + n_S*0))
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' test_innerproduct_fast_T_S_q_1'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call test_innerproduct_fast_Z_S_q_1(verbose-1,n_svd_r,svd_r_
     $           ,n_svd_l,svd_l_,svd_s_,svd_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) ,fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp),n_r ,grid_p_,n_w_,n_A
     $           ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp) ,nS_per,ld_S
     $           ,S_p__(ld_S *IS_per),n_S,Z_S_q__(n_X_tmp))
         enddo !do nS_sub=0,n_S_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         if (verbose.gt.0) then
            write(6,'(A,F8.5)')
     $           ' finished calculating Z_S_q__ for each S. '
     $           ,timing_tot
            timing_tmp = (n_A*1.0d0)*(n_svd_l*1.0d0)*(n_S*1.0d0)
     $           /timing_tot/1e9
            write(6,'(A,F8.4)') ' Z_S_q__ total Gnumps: ' , timing_tmp
            write(6,'(A)') ' '
         end if ! if (verbose.gt.0) then
      end if !if (svd_calculation_type.ne.2) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array T_R_CTF_M_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each image'
         write(6,'(A)') ' T_{-delta_est}(R_{-gamma_est}(CTF.*M)), '
         write(6,'(A)') ' where T = translation by -delta_est.'
         write(6,'(A)') ' and R = rotation by -gamma_est.'
         write(6,'(A)') ' '
      end if
      allocate(T_R_CTF_M_q__(0:n_r*n_M*n_w_max))
      call cl1_c16(n_r*n_M*n_w_max,T_R_CTF_M_q__)
      allocate(M_p_(0:n_A-1))
      allocate(M_q_(0:n_A-1))
      allocate(M_p_omp__(0:n_A*n_omp_sub-1))
      allocate(M_q_omp__(0:n_A*n_omp_sub-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' Calculate bessel coefficients of T(R(CTF.*M)).'
         write(6,'(A)') ' T = translation by -delta_x_est. '
         write(6,'(A)') ' R = rotation by -gamma_est.'
         write(6,'(A)') ' CTF = conjg(CTF) in k-space polar coord.'
         write(6,'(A)') ' M = template in k-space polar coord.'
         write(6,'(A)') ' '
         write(6,'(A,A)') ' T_R_CTF_M_q__(nr + n_r*(nm + n_M*nw)) '
         write(6,'(A)') ' is equal to the bessel-coefficients'
         write(6,'(A)') ' M_q_(ic)) '
         write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
         write(6,'(A)') ' where M_q_ = T(R(CTF.*M)), with '
         write(6,'(A)') ' T <-- -delta_est'
         write(6,'(A)') ' R <-- -gamma_est'
         write(6,'(A)') ' and CTF <-- conjg(CTF_(nctf))'
         write(6,'(A)') ' and M <-- M_p__(nm*ld_M).'
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_X_tmp = 0 + n_r*(IM_per + n_M*0)
         call test_innerproduct_fast_T_R_CTF_M_q_1(verbose-1
     $        ,delta_x_est_(IM_per),delta_y_est_(IM_per)
     $        ,gamma_z_est_(IM_per),ctf_ind_est_(IM_per)
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_p_,n_w_,n_A
     $        ,M_p_omp__(n_A_tmp),M_q_omp__(n_A_tmp),nM_per,ld_M
     $        ,M_p__(ld_M*IM_per) ,CTF_p_omp__(n_A_tmp),n_CTF ,ld_CTF
     $        ,CTF_p__ ,C_M_(IM_per) ,n_M,T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating T_R_CTF_M_q__ for each M. '
     $        ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' T_R_CTF_M_q__ total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      if (svd_calculation_type.eq.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and T_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      allocate(S_T_T_R_CTF_M_q__(0:n_transf*n_S*n_M*n_w_max))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),T_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_T_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_STxTRM_0 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call test_innerproduct_fast_STxTRM_0(verbose-1,n_delta_x
     $        ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $        ,delta_x_est_,delta_y_est_ ,gamma_z_est_,ctf_ind_est_
     $        ,fftw_plan_frwd_,fftw_in1_ ,fftw_out_,fftw_plan_back_last
     $        ,fftw_in1_last_,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_
     $        ,S_q_,n_S,ld_S,S_p__,M_p_,M_q_,n_M ,ld_M ,M_p__,CTF_p_
     $        ,n_CTF ,ld_CTF,CTF_p__,C_M_,CTF_R_S_ ,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.eq.2) then

      if (svd_calculation_type.ne.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and Z_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      allocate(S_Z_T_R_CTF_M_q__(0:n_transf*n_S*n_M*n_w_max))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
      if (verbose.gt.1) then
         write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_svdd_. '
         write(6,'(A)') ' '
      end if
      allocate(S_T_T_R_CTF_M_q__(0:n_delta_x*n_delta_y*n_S*n_M*n_w_max))
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_A_tmp,n_X_tmp)
c$OMP DO 
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_x*n_delta_y*n_S*n_M*nw
         n_X_tmp = n_transf*n_S*n_M*nw
         call zgemm('T','N',n_delta_x*n_delta_y,n_S*n_M*1,n_svd_l ,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q__(n_A_tmp),n_delta_x *n_delta_y)
      enddo ! do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*n_delta_y*1.0d0)*(n_S
     $     *1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_SZxTRM_0 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call  test_innerproduct_fast_SZxTRM_0(verbose-1,n_svd_r,svd_r_
     $     ,n_svd_l,svd_l_,svd_s_,svd_V_r_,Z_svdd_,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,delta_x_est_
     $     ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $     ,fftw_in1_ ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $     ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S
     $     ,S_p__,M_p_,M_q_,n_M ,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF
     $     ,CTF_p__,C_M_,CTF_R_S_ ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.ne.2) then

      if (flag_MS_vs_SM.eqv..false.) then
      allocate(ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(ZZ_omp__(0:n_delta_x*n_delta_y*n_gamma_z*n_omp_sub-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' SM: matching templates to images; '
         write(6,'(A)')
     $        ' calling test_innerproduct_fast_ZZ_1. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_delta_x*n_delta_y*n_gamma_z
         n_X_tmp = 0 + n_delta_x*n_delta_y*(0 + n_S*(IM_per + n_M*0))
         call test_innerproduct_fast_ZZ_1(verbose-1,n_delta_x,delta_x_
     $        ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r ,grid_p_,n_w_
     $        ,n_A ,n_S,S_alpha_polar_a_,S_alpha_azimu_b_,nM_per
     $        ,polar_a_est_(IM_per) ,azimu_b_est_(IM_per)
     $        ,gamma_z_est_(IM_per) ,delta_x_est_(IM_per)
     $        ,delta_y_est_(IM_per) ,l2_norm_est_(IM_per)
     $        ,ctf_ind_est_(IM_per) ,S_index_est_(IM_per)
     $        ,polar_a_upd_(IM_per) ,azimu_b_upd_(IM_per)
     $        ,gamma_z_upd_(IM_per) ,delta_x_upd_(IM_per)
     $        ,delta_y_upd_(IM_per) ,l2_norm_upd_(IM_per)
     $        ,ctf_ind_upd_(IM_per) ,S_index_upd_(IM_per)
     $        ,alpha_update_f ,displacement_max ,C_M_(IM_per),CTF_R_S_
     $        ,n_M,S_T_T_R_CTF_M_q__(n_X_tmp),ZZ_omp__(n_A_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished ZZ. ' , timing_tot
      timing_tmp = (n_gamma_z*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $     *(n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' ZZ total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      deallocate(ZZ_)
      deallocate(ZZ_omp__)
      end if !if (flag_MS_vs_SM.eqv..false.) then

      if (flag_MS_vs_SM.eqv..true.) then
      allocate(ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(ZZ_omp__(0:n_delta_x*n_delta_y*n_gamma_z*n_omp_sub-1))
      allocate(C_Z_MS_(0:n_M*n_S-1))
      allocate(C_S_MS_(0:n_M*n_S-1))
      allocate(delta_x_MS_(0:n_M*n_S-1))
      allocate(delta_y_MS_(0:n_M*n_S-1))
      allocate(gamma_z_MS_(0:n_M*n_S-1))
      allocate(C_Z_sort_MS_(0:n_M*n_S-1))
      allocate(C_S_sort_MS_(0:n_M*n_S-1))
      allocate(delta_x_sort_MS_(0:n_M*n_S-1))
      allocate(delta_y_sort_MS_(0:n_M*n_S-1))
      allocate(gamma_z_sort_MS_(0:n_M*n_S-1))
      allocate(I_permute_MS_(0:n_M*n_S-1))
      allocate(I_inverse_MS_(0:n_M*n_S-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' MS: matching images to templates; '
         write(6,'(A)')
     $        ' calling test_innerproduct_fast_MS_0. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_delta_x*n_delta_y*n_gamma_z
         n_X_tmp = 0 + n_delta_x*n_delta_y*(0 + n_S*(IM_per + n_M*0))
         call test_innerproduct_fast_MS_0(verbose-1,n_delta_x,delta_x_
     $        ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r ,grid_p_,n_w_
     $        ,n_A ,n_S,S_alpha_polar_a_,S_alpha_azimu_b_,nM_per
     $        ,polar_a_est_(IM_per) ,azimu_b_est_(IM_per)
     $        ,gamma_z_est_(IM_per) ,delta_x_est_(IM_per)
     $        ,delta_y_est_(IM_per) ,l2_norm_est_(IM_per)
     $        ,ctf_ind_est_(IM_per) ,S_index_est_(IM_per)
     $        ,delta_x_MS_(IM_per), delta_y_MS_(IM_per),
     $        gamma_z_MS_(IM_per), C_S_MS_(IM_per), C_Z_MS_(IM_per)
     $        ,displacement_max ,C_M_(IM_per),CTF_R_S_ ,n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp),ZZ_omp__(n_A_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_0. ' , timing_tot
      timing_tmp = (n_gamma_z*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $     *(n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_0 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_MS_sort_0(verbose,n_S,n_M
     $     ,delta_x_MS_,delta_y_MS_,gamma_z_MS_,C_S_MS_,C_Z_MS_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_sort_0. ' , timing_tot
      timing_tmp = (n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_sort_0 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_MS_update_0(verbose,n_S,n_M,C_M_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_
     $     ,S_alpha_polar_a_,S_alpha_azimu_b_ , gamma_z_est_
     $     ,delta_x_est_ ,delta_y_est_,polar_a_upd_ ,azimu_b_upd_
     $     ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $     ,S_index_upd_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_update_0. ' , timing_tot
      timing_tmp = (n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_update_0 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      deallocate(ZZ_)
      deallocate(ZZ_omp__)
      deallocate(C_Z_MS_)
      deallocate(C_S_MS_)
      deallocate(delta_x_MS_)
      deallocate(delta_y_MS_)
      deallocate(gamma_z_MS_)
      deallocate(C_Z_sort_MS_)
      deallocate(C_S_sort_MS_)
      deallocate(delta_x_sort_MS_)
      deallocate(delta_y_sort_MS_)
      deallocate(gamma_z_sort_MS_)
      deallocate(I_permute_MS_)
      deallocate(I_inverse_MS_)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      
      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
c$$$      deallocate(SM_use_sub_)
      deallocate(S_T_T_R_CTF_M_q__)
      if (svd_calculation_type.ne.2) then
         deallocate(S_Z_T_R_CTF_M_q__)
      end if ! if (svd_calculation_type.ne.2) then
      deallocate(T_R_CTF_M_q__)
      if (svd_calculation_type.eq.2) then
         deallocate(T_S_q__)
      end if ! if (svd_calculation_type.eq.2) then
      if (svd_calculation_type.ne.2) then
         deallocate(Z_S_q__)
      end if ! if (svd_calculation_type.ne.2) then
      deallocate(S_p_omp__)
      deallocate(S_q_omp__)
      deallocate(M_p_omp__)
      deallocate(M_q_omp__)
      deallocate(CTF_p_omp__)
      deallocate(CTF_q_omp__)
      deallocate(S_p_)
      deallocate(S_q_)
      deallocate(M_p_)
      deallocate(M_q_)
      deallocate(CTF_p_)
      deallocate(CTF_q_)
      deallocate(CTF_R_S_)
      deallocate(Z_svdd_)
      deallocate(gamma_z_)
      deallocate(delta_y_)
      deallocate(delta_x_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan_(fftw_plan_back_(nr))
         call dfftw_destroy_plan_(fftw_plan_frwd_(nr))
      enddo
      deallocate(fftw_out_)
      deallocate(fftw_in1_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_plan_frwd_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_many_plans for omp.'
      end if
      do nomp_sub=0,n_omp_sub-1
         call dfftw_destroy_plan_(fpm_back_(nomp_sub))
c$$$         call dfftw_destroy_plan_(fpm_frwd_(nomp_sub))
      enddo !do nomp_sub=0,n_omp_sub-1
      deallocate(fpm_in1__)
      deallocate(fpm_out__)
      deallocate(fpm_n__)
      deallocate(fpm_inembed__)
      deallocate(fpm_onembed__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for omp.'
      end if
      do nr=0,n_omp_sub*n_r-1
         call dfftw_destroy_plan_(fftw_plan_back__(nr))
         call dfftw_destroy_plan_(fftw_plan_frwd__(nr))
      enddo
      deallocate(fftw_out__)
      deallocate(fftw_in1__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_plan_frwd__)
      deallocate(n_S_per_)
      deallocate(I_S_per_)
      deallocate(n_M_per_)
      deallocate(I_M_per_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating svd arrays.'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_fast_stage_1]'
      end if

      end


