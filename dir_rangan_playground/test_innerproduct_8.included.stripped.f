






























      subroutine test_innerproduct_8(



     $     verbose 
     $     ,flag_memory_estimate 
     $     ,d_memory_estimate 
     $     ,rseed 
     $     ,n_k_cur 
     $     ,n_polar_a_ 
     $     ,grid_k_p_ 
     $     ,half_diameter_x_c 
     $     ,n_S 
     $     ,I_S_sample_ 
     $     ,ld_S 
     $     ,S_k_p__ 
     $     ,tesselation_distance_req 
     $     ,n_LT_add 
     $     ,n_LT_ref 
     $     ,S_alpha_S_index_ 
     $     ,S_alpha_polar_a_ 
     $     ,S_alpha_azimu_b_ 
     $     ,n_M 
     $     ,I_M_sample_ 
     $     ,ld_M 
     $     ,M_k_p__ 
     $     ,n_CTF 
     $     ,ld_CTF 
     $     ,CTF_k_p__ 
     $     ,alpha_est__ 
     $     ,alpha_update_f 
     $     ,flag_MS_vs_SM 
     $     ,n_SM_max 
     $     ,n_SM_ 
     $     ,alpha_SM__ 
     $     ,n_MS_max 
     $     ,n_MS_ 
     $     ,alpha_MS__ 
     $     ,n_pixels_in 
     $     ,displacement_max 
     $     ,n_delta_v 
     $     ,delta_x_ 
     $     ,delta_y_ 
     $     ,n_gamma_z 
     $     ,svd_calculation_type 
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
     $     ,eps_svd 
     $     ,flag_RTRT_vs_RTTR 
     $     ,fpm_howmany_max 
     $     ,n_omp_sub__in 
     $     ,n_S_0_sub__in 
     $     ,n_S_1_sub__in 
     $     ,n_M_0_sub__in 
     $     ,n_M_1_sub__in 
     $     ,C_M_ 
     $     )



      implicit none
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
      include 'omp_lib.h'



      integer *4 n_alpha
      parameter(n_alpha=11)
      integer *4 nalpha_polar_a
      parameter (nalpha_polar_a=0)
      integer *4 nalpha_azimu_b
      parameter (nalpha_azimu_b=1)
      integer *4 nalpha_gamma_z
      parameter (nalpha_gamma_z=2)
      integer *4 nalpha_delta_x
      parameter (nalpha_delta_x=3)
      integer *4 nalpha_delta_y
      parameter (nalpha_delta_y=4)
      integer *4 nalpha_l2_norm
      parameter (nalpha_l2_norm=5)
      integer *4 nalpha_ctf_ind
      parameter (nalpha_ctf_ind=6)
      integer *4 nalpha_S_index
      parameter (nalpha_S_index=7)
      integer *4 nalpha_M_index
      parameter (nalpha_M_index=8)
      integer *4 nalpha_CTF_R_S
      parameter (nalpha_CTF_R_S=9)
      integer *4 nalpha_C_Z_opt
      parameter (nalpha_C_Z_opt=10)
      integer *4 verbose 
      logical flag_memory_estimate 
      logical flag_memory_checkset 
      real *8 d_memory_estimate 
      integer *4 rseed 
      integer *4 n_k_cur 
      integer *4 n_polar_a_(0:0) 
      real *8 grid_k_p_(0:0) 
      real *8 half_diameter_x_c 
      integer *4 n_S 
      integer *4 I_S_sample_(0:0) 
      integer *4 ld_S 
      complex *16 S_k_p__(0:0) 
      real *8 tesselation_distance_req 
      real *8 S_alpha_S_index_(0:0) 
      real *8 S_alpha_polar_a_(0:0) 
      real *8 S_alpha_azimu_b_(0:0) 
      integer *4 n_M 
      integer *4 I_M_sample_(0:0) 
      integer *4 ld_M 
      complex *16 M_k_p__(0:0) 
      integer *4 n_CTF 
      integer *4 ld_CTF 
      complex *16 CTF_k_p__(0:0) 
      real *8 alpha_est__(0:0) 
      real *8 alpha_update_f 
      logical flag_MS_vs_SM 
      integer *4 n_SM_max 
      integer *4 n_SM_(0:0) 
      real *8 alpha_SM__(0:0) 
      integer *4 n_MS_max 
      integer *4 n_MS_(0:0) 
      real *8 alpha_MS__(0:0) 
      real *8 n_pixels_in 
      real *8 displacement_max 
      integer *4 n_delta_v 
      real *8 delta_x_(0:0) 
      real *8 delta_y_(0:0) 
      integer *4 n_gamma_z 
      integer *4 svd_calculation_type 
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
      real *8 eps_svd 
      logical flag_RTRT_vs_RTTR 
      integer *4 fpm_howmany_max 
      integer *4 n_omp_sub__in 
      integer *4 n_S_0_sub__in 
      integer *4 n_S_1_sub__in 
      integer *4 n_M_0_sub__in 
      integer *4 n_M_1_sub__in 
      complex *16 C_M_(0:0) 
c$$$      %%%%%%%%%%%%%%%%
c$$$      function outputs
c$$$      %%%%%%%%%%%%%%%%
      integer *4 min_i4_f,max_i4_f,sum_i4_f 
      real *8 max_r8_f
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%
      logical flag_tesselation 
      logical flag_tesselation_ref 
      integer *4 n_point_init 
      integer *4 n_point,npoint 
      real *8, allocatable :: S_L_(:) 
      real *8 tradius_min 
      integer *4 nl_max_init,nm_sum_init,ll_sum_init 
      integer *4 nl_max,nm_sum,ll_sum 
      integer *4 n_T_root_base,nT_root_base 
      logical parity_input 
c$$$      T_ tesselation tree
      integer *4, allocatable :: T_nl_(:) 
      real *8   , allocatable :: T_vm_(:) 
      real *8   , allocatable :: T_tr_(:) 
      integer *4, allocatable :: T_ll_(:) 
      logical   , allocatable :: T_lf_(:) 
      integer *4, allocatable :: T_c0_(:) 
      integer *4, allocatable :: T_c1_(:) 
      integer *4, allocatable :: T_c2_(:) 
      integer *4, allocatable :: T_c3_(:) 
      integer *4, allocatable :: T_ls_(:) 
      integer *4, allocatable :: T_LT_(:) 
c$$$      T_ roots
      integer *4, allocatable :: T_root_base_(:) 
c$$$      %%%%%%%%%%%%%%%%
c$$$      bookkeeping tesselation search across threads
c$$$      %%%%%%%%%%%%%%%%
      real *8 tesselation_distance_req_use 
      integer *4, allocatable :: n_SM_use_local_(:) 
      integer *4 n_SM_use 
      real *8, allocatable :: vp_input_omp__(:) 
      integer *8 p_vp_input 
      logical, allocatable :: flag_S_use_omp__(:) 
      integer *8 p_flag_S_use 
      integer *4 n_LT,n_S_use,n_S_use_sum
      integer *4 n_LT_add 
      integer *4 n_LT_ref 
      integer *4, allocatable :: n_S_use_sum_(:) 
      integer *4 n_add,n_ref,nref,nLT,ns_local
      integer *4 n_pass,n_pass_ref
      logical continue_flag
      integer *4, allocatable :: LT_omp__(:) 
      integer *8 p_LT 
      real *8, allocatable :: S_alpha_S_index_local_omp__(:) 
      integer *8 p_S_alpha_S_index_local 
      real *8, allocatable :: S_alpha_polar_a_local_omp__(:) 
      integer *8 p_S_alpha_polar_a_local 
      real *8, allocatable :: S_alpha_azimu_b_local_omp__(:) 
      integer *8 p_S_alpha_azimu_b_local 
      integer *4, allocatable :: I_S_sample_local_omp__(:) 
      integer *8 p_I_S_sample_local 
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%      
      real *8 eps_target 
c$$$      declaration of svd-expansion and associated variables
      integer *4 n_svd_r,n_svd_d,n_svd_l 
      integer *4 nsvd_r,nsvd_d,nsvd_l 
      real *8, allocatable :: svd_r_(:) 
      real *8, allocatable :: svd_d_(:) 
      integer *4, allocatable :: svd_l_(:) 
      real *8, allocatable :: svd_U_d_(:) 
      real *8, allocatable :: svd_s_(:) 
      real *8, allocatable :: svd_V_r_(:) 
      integer *4 svd_unitnumber
      parameter (svd_unitnumber=143857)
      character(len=64) svd_fname,svd_format_string
c$$$      include './dir_gen_Jsvd_6/Jsvd_var_allocation_off_0.f'
      integer *4 n_svd_max
      parameter (n_svd_max=512)
      real *8 svd_d_max,svd_r_max 
      real *8 , allocatable :: svd_polyval_U_d_(:) 
      real *8 , allocatable :: svd_polyval_V_r_(:) 
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels 
      complex *16, allocatable :: Z_S_svdd_(:) 
      complex *16, allocatable :: Z_M_svdd_(:) 
      complex *16, allocatable :: CTF_R_S__(:) 
      complex *16, allocatable :: CTF_R_S_sub__(:) 
      complex *16, allocatable :: CTF_R_S_sub_omp__(:) 
      complex *16, allocatable :: CTF_R_S_local_omp__(:) 
      integer *8 p_CTF_R_S_local 
      complex *16, allocatable :: O_S_q__(:) 
      complex *16, allocatable :: O_S_q_local_omp__(:) 
      integer *8 p_O_S_q_local 
      complex *16, allocatable :: T_S_q__(:) 
      complex *16, allocatable :: T_S_q_local_omp__(:) 
      integer *8 p_T_S_q_local 
      complex *16, allocatable :: Z_S_q__(:) 
      complex *16, allocatable :: Z_S_q_local_omp__(:) 
      integer *8 p_Z_S_q_local 
      pointer (p_S_p__,S_p__) 
      complex *16 S_p__(0:0)
      pointer (p_M_p__,M_p__) 
      complex *16 M_p__(0:0)
      pointer (p_CTF_p__,CTF_p__) 
      complex *16 CTF_p__(0:0)
      complex *16, allocatable :: S_p_(:) 
      complex *16, allocatable :: S_q_(:) 
      complex *16, allocatable :: S_p_omp__(:) 
      complex *16, allocatable :: S_q_omp__(:) 
      complex *16, allocatable :: M_p_(:) 
      complex *16, allocatable :: M_q_(:) 
      complex *16, allocatable :: M_p_omp__(:) 
      complex *16, allocatable :: M_q_omp__(:) 
      complex *16, allocatable :: CTF_p_(:) 
      complex *16, allocatable :: CTF_q_(:) 
      complex *16, allocatable :: CTF_p_omp__(:) 
      complex *16, allocatable :: CTF_q_omp__(:) 
      complex *16, allocatable :: Z_p_(:) 
      complex *16, allocatable :: Z_q_(:) 
      complex *16, allocatable :: Z_p_omp__(:) 
      complex *16, allocatable :: Z_q_omp__(:) 
      complex *16, allocatable :: ZZ__(:) 
      complex *16, allocatable :: ZZ_omp__(:) 
      complex *16, allocatable :: ZZ_sub_(:) 
      complex *16, allocatable :: ZZ_sub_omp__(:) 
      complex *16, allocatable :: C_trn0_(:) 
      complex *16, allocatable :: C_trn0_omp__(:) 
c$$$      %%%%%%%%%%%%%%%%
c$$$      precomputed arrays for zgemm
c$$$      %%%%%%%%%%%%%%%%
      complex *16, allocatable :: O_T_R_CTF_M_q__(:) 
      complex *16, allocatable :: O_T_R_CTF_M_q_local_omp__(:) 
      integer *8 p_O_T_R_CTF_M_q_local 
      complex *16, allocatable :: T_T_R_CTF_M_q__(:) 
      complex *16, allocatable :: T_T_R_CTF_M_q_local_omp__(:) 
      integer *8 p_T_T_R_CTF_M_q_local 
      complex *16, allocatable :: Z_T_R_CTF_M_q__(:) 
      complex *16, allocatable :: Z_T_R_CTF_M_q_local_omp__(:) 
      integer *8 p_Z_T_R_CTF_M_q_local 
      complex *16, allocatable :: S_T_T_R_CTF_M_q__(:) 
      complex *16, allocatable :: S_T_T_R_CTF_M_q_local_omp__(:) 
      integer *8 p_S_T_T_R_CTF_M_q_local 
      complex *16, allocatable :: S_Z_T_R_CTF_M_q__(:) 
      complex *16, allocatable :: S_Z_T_R_CTF_M_q_local_omp__(:) 
      integer *8 p_S_Z_T_R_CTF_M_q_local 
      integer *4 n_C_trn0 
      integer *4 n_C_trn0_omp 
      integer *4 n_CTF_R_S_local_omp 
      integer *4 n_S_alpha_S_index_local_omp 
      integer *4 n_O_S_q_local_omp 
      integer *4 n_T_S_q_local_omp 
      integer *4 n_Z_S_q_local_omp 
      integer *4 n_O_T_R_CTF_M_q__ 
      integer *4 n_O_T_R_CTF_M_q_local_omp__ 
      integer *4 n_T_T_R_CTF_M_q__ 
      integer *4 n_T_T_R_CTF_M_q_local_omp__ 
      integer *4 n_Z_T_R_CTF_M_q__ 
      integer *4 n_Z_T_R_CTF_M_q_local_omp__ 
      integer *4 n_S_T_T_R_CTF_M_q__ 
      integer *4 n_S_T_T_R_CTF_M_q_local_omp__ 
      integer *4 n_S_Z_T_R_CTF_M_q__ 
      integer *4 n_S_Z_T_R_CTF_M_q_local_omp__ 
c$$$      indices
      integer *4 nr,na,ns,nctf,nm,nw
      integer *4 n_r 
      integer *4, allocatable :: n_w_(:) 
      integer *4, allocatable :: n_w_csum_(:) 
      integer *4 n_w_max 
      integer *4 n_A 
c$$$      temporary: image-parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: ctf_ind_est_(:)
      real *8, allocatable :: S_index_est_(:)
      real *8, allocatable :: M_index_est_(:)
      real *8, allocatable :: alpha__in_(:)
      real *8, allocatable :: alpha__in_omp__(:)
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
      integer *4 n_transf 
      integer *4 fpm_howmany
      integer *4 n_fpm,nfpm
      integer *8, allocatable :: fpm_frwd_(:)
      integer *8, allocatable ::  fpm_back_(:)
      complex *16, allocatable :: fpm_in1__(:)
      complex *16, allocatable :: fpm_out__(:)
      integer *4 fpm_rank
      integer *4, allocatable :: fpm_n__(:)
      integer *4, allocatable :: fpm_inembed__(:)
      integer *4 fpm_istride
      integer *4 fpm_idist
      integer *4, allocatable :: fpm_onembed__(:)
      integer *4 fpm_ostride
      integer *4 fpm_odist
c$$$      array of displacements and rotations to measure
      integer *4 ndv,ngz,ndv_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: gamma_z_(:) 
c$$$      parameters for blocking S and M
      integer *4 nx0,nx1,nx2,nx3,nx4
      integer *4 nm0,nm1,nm2,nm3,nm4
      integer *4 ns0,ns1,ns2,ns3,ns4
      integer *4 n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp,n_S_tmp
      integer *4 nomp_sub_use,nomp_sub,nomp_per,nomp_sum
      integer *4 n_S_9_sub_use,nS_9_sub,nS_9_per,nS_9_sum
      integer *4, allocatable :: n_S_9_per_(:)
      integer *4, allocatable :: n_S_9_sum_(:)
      integer *4 nS_9_per_max,nS_9_per_min
      integer *4 n_M_9_sub_use,nM_9_sub,nM_9_per,nM_9_sum
      integer *4, allocatable :: n_M_9_per_(:)
      integer *4, allocatable :: n_M_9_sum_(:)
      integer *4 nM_9_per_max,nM_9_per_min
      integer *4 n_S_0_sub_use,nS_0_sub,nS_0_per,nS_0_sum
      integer *4, allocatable :: n_S_0_per_(:)
      integer *4, allocatable :: n_S_0_sum_(:)
      integer *4 nS_0_per_max,nS_0_per_min
      integer *4 n_M_0_sub_use,nM_0_sub,nM_0_per,nM_0_sum
      integer *4, allocatable :: n_M_0_per_(:)
      integer *4, allocatable :: n_M_0_sum_(:)
      integer *4 nM_0_per_max,nM_0_per_min
      integer *4 n_S_1_sub_use,nS_1_sub,nS_1_per,nS_1_sum
      integer *4, allocatable :: n_S_1_per_(:)
      integer *4, allocatable :: n_S_1_sum_(:)
      integer *4 nS_1_per_max,nS_1_per_min
      integer *4 n_M_1_sub_use,nM_1_sub,nM_1_per,nM_1_sum
      integer *4, allocatable :: n_M_1_per_(:)
      integer *4, allocatable :: n_M_1_sum_(:)
      integer *4 nM_1_per_max,nM_1_per_min
c$$$      parameters for memory map
      real *8 d_mem
      integer *4 verbose_mem,verbose_timing
c$$$      parameters for timing and testing
c$$$      Note that timing summation is not thread-safe
      real *8 timing_tic,timing_toc,timing_tot,timing_tmp,gnump_tot
      logical flag_fill 
      parameter(flag_fill=.false.)
      logical flag_test 
      real *8 timing_total_fftw_plan,gnump_total_fftw_plan
      real *8 timing_total_CTF_R_S,gnump_total_CTF_R_S
      real *8 timing_total_O_S_q,gnump_total_O_S_q
      real *8 timing_total_T_S_q,gnump_total_T_S_q
      real *8 timing_total_Z_S_q,gnump_total_Z_S_q
      real *8 timing_total_O_T_R_CTF_M_q,gnump_total_O_T_R_CTF_M_q
      real *8 timing_total_T_T_R_CTF_M_q,gnump_total_T_T_R_CTF_M_q
      real *8 timing_total_Z_T_R_CTF_M_q,gnump_total_Z_T_R_CTF_M_q
      real *8 timing_total_zgemm,gnump_total_zgemm
      real *8 timing_total_fpm_fill,gnump_total_fpm_fill
      real *8 timing_total_fpm_fftw,gnump_total_fpm_fftw
      real *8 timing_total_transpose,gnump_total_transpose
      real *8 timing_total_Zstore,gnump_total_Zstore
      logical flag_time_Zstore
      real *8 timing_total_Zstore_a,gnump_total_Zstore_a
      real *8 timing_total_Zstore_b,gnump_total_Zstore_b
      real *8 timing_total_Zstore_c,gnump_total_Zstore_c
      real *8 timing_total_Zstore_x,gnump_total_Zstore_x
      real *8 timing_total_Zstore_d,gnump_total_Zstore_d
      real *8 timing_total_Zstore_e,gnump_total_Zstore_e
      character(len=1024) timing_string
c$$$      format_string for printing output
      character(len=64) format_string
      logical flag_MDA
      integer *4 MDA_n_d
      integer *4 MDA_d_(0:64-1)
      character(len=1024) MDA_string
c$$$      pi
      real *8 pi
      pi = 4.0d0*atan(1.0d0)



      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[entering test_innerproduct_8]: '
      end if

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ',verbose
         write(6,'(A,L1)') ' flag_memory_estimate: '
     $        ,flag_memory_estimate
         write(6,'(A,I0)') ' rseed: ',rseed
         write(6,'(A,I0)') ' n_k_cur: ',n_k_cur
         call print_all_i4(n_k_cur,n_polar_a_,' n_polar_a_: ')
         call print_all_r8(n_k_cur,grid_k_p_,' grid_k_p_: ')
         write(6,'(A,F6.3)') ' half_diameter_x_c: ',half_diameter_x_c
         write(6,'(A,I0)') ' n_S: ',n_S
c$$$         call print_sub_i4(n_S,I_S_sample_,' I_S_sample_: ')
         write(6,'(A,I0)') ' ld_S: ',ld_S
c$$$         call print_sub_c16(n_S*ld_S,S_k_p__,' S_k_p___: ')
         write(6,'(A,F6.3)') ' tesselation_distance_req: '
     $        ,tesselation_distance_req
c$$$         call print_sub_r8(n_S,S_alpha_S_index_
c$$$     $        ,' S_alpha_S_index_: ')
c$$$         call print_sub_r8(n_S,S_alpha_polar_a_
c$$$     $        ,' S_alpha_polar_a_: ')
c$$$         call print_sub_r8(n_S,S_alpha_azimu_b_
c$$$     $        ,' S_alpha_azimu_b_: ')
         write(6,'(A,I0)') ' n_M: ',n_M
c$$$         call print_sub_i4(n_M,I_M_sample_,' I_M_sample_: ')
         write(6,'(A,I0)') ' ld_M: ',ld_M
c$$$         call print_sub_c16(n_M*ld_M,M_k_p__,' M_k_p___: ')
         write(6,'(A,I0)') ' n_CTF: ',n_CTF
         write(6,'(A,I0)') ' ld_CTF: ',ld_CTF
c$$$         call print_sub_c16(n_CTF*ld_CTF,CTF_k_p__,' CTF_k_p___: ')
         write(6,'(A,I0)') ' n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
c$$$         write(6,'(A)') ' alpha_est__: '
c$$$         write(6,format_string) (alpha_est__(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') ' alpha_update_f: ' , alpha_update_f
         write(6,'(A,L1)') ' flag_MS_vs_SM: ' , flag_MS_vs_SM
c$$$         if (flag_MS_vs_SM.eqv..true.) then
c$$$            write(6,'(A,I0)') ' n_MS_max: ',n_MS_max
c$$$            call print_sub_i4(n_S,n_MS_,' n_MS_: ')
c$$$            call print_sub_r8(n_S*n_MS_max*n_alpha,alpha_MS__
c$$$     $        ,' alpha_MS__: ')
c$$$         end if 
c$$$         if (flag_MS_vs_SM.eqv..false.) then
c$$$            write(6,'(A,I0)') ' n_SM_max: ',n_SM_max
c$$$            call print_sub_i4(n_M,n_SM_,' n_SM_: ')
c$$$            call print_sub_r8(n_M*n_SM_max*n_alpha,alpha_SM__
c$$$     $        ,' alpha_SM__: ')
c$$$         end if 
         write(6,'(A,F6.3)') ' n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') ' displacement_max: ',displacement_max
         write(6,'(A,I0)') ' n_delta_v: ',n_delta_v
         write(6,'(A,I0)') ' n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') ' svd_calculation_type: '
     $        ,svd_calculation_type
         write(6,'(A,F6.3)') ' eps_svd: ',eps_svd
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' , flag_RTRT_vs_RTTR
         write(6,'(A,I0)') ' fpm_howmany_max: ',fpm_howmany_max
         write(6,'(A,I0)') ' n_omp_sub__in: ',n_omp_sub__in
         write(6,'(A,I0)') ' n_S_0_sub__in: ',n_S_0_sub__in
         write(6,'(A,I0)') ' n_S_1_sub__in: ',n_S_1_sub__in
         write(6,'(A,I0)') ' n_M_0_sub__in: ',n_M_0_sub__in
         write(6,'(A,I0)') ' n_M_1_sub__in: ',n_M_1_sub__in
      end if 

       if (verbose.gt.0) then
          write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
       end if 
       call openblas_set_num_threads(1)





      pi = 4*atan(1.0)
      eps_target = eps_svd
      p_S_p__ = loc(S_k_p__(0))
      p_M_p__ = loc(M_k_p__(0))
      p_CTF_p__ = loc(CTF_k_p__(0))
      flag_tesselation = (tesselation_distance_req.lt.2.0d0)
      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         allocate(n_SM_use_local_(0:1+n_M-1))
         call cs1_i4(n_M,n_SM_use_local_)
         d_memory_estimate = d_memory_estimate + 4.0d0*n_M
         if (n_M_0_sub__in.ne.n_M) then
            write(6,'(A,A,I0,A)')
     $           ' Because flag_tesselation.eqv..true.: '
     $           ,' Setting n_M_0_sub__in = n_M: ' , n_M ,
     $           ' and n_M_1_sub__in = 1 '
         end if 
         n_M_0_sub__in = n_M
         n_M_1_sub__in = 1
         if (verbose.gt.1) then
            write(6,'(A)') ''
            write(6,'(A,A)') 'Now we create tesselation'
     $           ,' to access templates efficiently:'
         end if
         n_point = n_S
         n_point_init = n_point
         allocate(S_L_(0:1+3*n_point-1))
         call cs1_r8(3*n_point,S_L_)
         d_memory_estimate = d_memory_estimate + 8.0d0*3*n_point
         do npoint=0,n_point-1
            call get_angle_to_vp_(S_alpha_polar_a_(npoint)
     $           ,S_alpha_azimu_b_(npoint),S_L_(0+3*npoint))
            call normalize_r8(3,S_L_(0 + 3*npoint))
         enddo                  
         if (verbose.gt.1) then
            do npoint=0,n_point-1
               write(6,'(A,F8.4,F8.4,A,A,I0,1X,3F8.4)') 'S_alpha: ' ,
     $              S_alpha_polar_a_(npoint) , S_alpha_azimu_b_(npoint)
     $              , ' --> ' , ' S_L_: ', npoint , S_L_(0 + 3*npoint),
     $              S_L_(1 + 3*npoint) , S_L_(2 + 3*npoint)
            enddo               
         end if                 
         tradius_min = 1.0d-6
         call tesselation_get_nl_nm_ll(n_point,S_L_,tradius_min,nl_max
     $        ,nm_sum,ll_sum)
         nl_max_init = nl_max
         nm_sum_init = nm_sum
         ll_sum_init = ll_sum
         allocate(T_nl_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_nl_)
         d_memory_estimate = d_memory_estimate + 4.0d0*nm_sum
         allocate(T_vm_(0:1+3*nm_sum-1)) 
         call cs1_r8(3*nm_sum,T_vm_)
         allocate(T_tr_(0:1+1*nm_sum-1)) 
         call cs1_r8(1*nm_sum,T_tr_)
         d_memory_estimate = d_memory_estimate + 4.0d0*8.0d0*nm_sum
         allocate(T_ll_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_ll_)
         allocate(T_lf_(0:1+1*nm_sum-1)) 
         call cs1_l2(1*nm_sum,T_lf_)
         allocate(T_c0_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_c0_)
         allocate(T_c1_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_c1_)
         allocate(T_c2_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_c2_)
         allocate(T_c3_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_c3_)
         allocate(T_ls_(0:1+1*nm_sum-1)) 
         call cs1_i4(1*nm_sum,T_ls_)
         d_memory_estimate = d_memory_estimate + 7.0d0*4.0d0*nm_sum
         allocate(T_LT_(0:1+1*ll_sum-1)) 
         call cs1_i4(1*ll_sum,T_LT_)
         d_memory_estimate = d_memory_estimate + 4.0d0*ll_sum
         allocate(T_root_base_(0:1+8-1))
         call cs1_i4(8,T_root_base_)
         d_memory_estimate = d_memory_estimate + 4.0d0*8
      end if 
      d_memory_estimate = 0.0d0
      if (flag_memory_estimate.eqv..true.) then
         verbose_mem = verbose + 1
      else
         verbose_mem = 0
      end if 
      verbose_timing = 0
      timing_total_fftw_plan = 0.0d0
      gnump_total_fftw_plan = 0.0d0
      timing_total_CTF_R_S = 0.0d0
      gnump_total_CTF_R_S = 0.0d0
      timing_total_O_S_q = 0.0d0
      gnump_total_O_S_q = 0.0d0
      timing_total_T_S_q = 0.0d0
      gnump_total_T_S_q = 0.0d0
      timing_total_Z_S_q = 0.0d0
      gnump_total_Z_S_q = 0.0d0
      timing_total_O_T_R_CTF_M_q = 0.0d0
      gnump_total_O_T_R_CTF_M_q = 0.0d0
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
      allocate(n_w_(0:1+n_k_cur-1))
      call cs1_i4(n_k_cur,n_w_)
      d_memory_estimate = d_memory_estimate + 4.0d0*n_k_cur
      allocate(n_w_csum_(0:1+n_k_cur-1))
      call cs1_i4(n_k_cur,n_w_csum_)
      d_memory_estimate = d_memory_estimate + 4.0d0*n_k_cur
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') ' n_k_cur = ',n_k_cur
     $        ,' calling get_template_size to'
     $        ,' determine n_w_ and n_A'
      end if 
      call get_template_size(n_polar_a_,n_k_cur,n_A,n_w_,n_w_csum_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_A = ',n_A
         call print_all_i4(n_k_cur,n_w_,' n_w_: ')
      end if 
      
c$$$      indices
      n_r = n_k_cur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r ' , n_r , ' < 2'
      end if 
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo 
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') ' n_w_max ',n_w_max,'; n_A ',n_A
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') ' We assume that an array of displacements '
         write(6,'(A)') ' has been passed in as input. '
         write(6,'(A)') ' It is usually reasonable for n_delta_v '
         write(6,'(A)') ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal (1+4*N_pixels_in)**2 or so.'
         write(6,'(A)') ' '
      end if 
      if (verbose.gt.1) then
         do ndv=0,n_delta_v-1
            write (6,'(A,I0,A,2(F8.4,1X))') ' ndv: ' , ndv ,
     $           ' (delta_x,delta_y): ' , delta_x_(ndv) , delta_y_(ndv)
         enddo 
      end if 
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)') ' n_w_max or so.'
         write(6,'(A)') ' '
      end if 
      allocate(gamma_z_(0:1+n_gamma_z-1))
      call cs1_r8(n_gamma_z,gamma_z_)
      d_memory_estimate = d_memory_estimate + 8.0d0*n_gamma_z
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.1) then
         call print_sub_r8(n_gamma_z,gamma_z_,' gamma_z_: ')
      end if 

      if ((svd_calculation_type.eq.1)) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if 
      allocate(svd_r_(0:1+n_svd_max-1)) 
      call cs1_r8(n_svd_max,svd_r_)
      allocate(svd_d_(0:1+n_svd_max-1)) 
      call cs1_r8(n_svd_max,svd_d_)
      allocate(svd_l_(0:1+n_svd_max-1))
      call cs1_i4(n_svd_max,svd_l_)
      allocate(svd_U_d_(0:1+n_svd_max*n_svd_max-1))
      call cs1_r8(n_svd_max*n_svd_max,svd_U_d_)
      allocate(svd_s_(0:1+n_svd_max-1))
      call cs1_r8(n_svd_max,svd_s_)
      allocate(svd_V_r_(0:1+n_svd_max*n_svd_max-1))
      call cs1_r8(n_svd_max*n_svd_max,svd_V_r_)
      d_memory_estimate = d_memory_estimate + 8.0d0*n_svd_max*(4.0d0 +
     $     2.0d0*n_svd_max)
      call get_svd_2(eps_target,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,grid_k_p_,n_r,n_delta_v,delta_x_,delta_y_,flag_warning,R_max
     $     ,K_max,delta_max,n_pixels)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_svd_r: ',n_svd_r
         write(6,'(A,I0)') ' n_svd_d: ',n_svd_d
         write(6,'(A,I0)') ' svd_unitnumber: ',svd_unitnumber
         write(6,'(A)') ' '
      end if 
      if (verbose.gt.0) then
         write(6,'(A,I0,A,A,I0)') ' n_svd_l: ' , n_svd_l , ' versus '
     $        ,' n_delta_v: ' , n_delta_v
      end if 
      d_memory_estimate = d_memory_estimate + 8.0*n_svd_l*n_delta_v
      if (flag_memory_estimate.eqv..false.) then
         allocate(svd_polyval_U_d_(0:1+n_svd_l*n_delta_v-1))
         call cs1_r8(n_svd_l*n_delta_v,svd_polyval_U_d_)
      end if 
      d_memory_estimate = d_memory_estimate + 8.0*n_svd_l*n_r
      if (flag_memory_estimate.eqv..false.) then
         allocate(svd_polyval_V_r_(0:1+n_svd_l*n_r-1))
         call cs1_r8(n_svd_l*n_r,svd_polyval_V_r_)
      end if 
      svd_d_max = 0.0d0
      do ndv=0,n_delta_v-1
         delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
         if (svd_d_max.lt.delta) then
            svd_d_max = delta
         end if 
      enddo 
      if (verbose.gt.1) then
         write(6,'(A,F8.6)') ' setting svd_d_max: ' , svd_d_max
      end if 
      if (flag_memory_estimate.eqv..false.) then
      call get_svd_polyval_U_d_(svd_d_max,n_svd_d,svd_d_,n_svd_l
     $     ,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $     ,svd_polyval_U_d_)
      end if 
      svd_r_max = 0.0d0
      svd_r_max = grid_k_p_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,F8.2)') ' setting svd_r_max: ' , svd_r_max
      end if 
      if (flag_memory_estimate.eqv..false.) then
      call get_svd_polyval_V_r_(svd_r_max,n_svd_r,svd_r_,n_svd_l
     $     ,svd_l_,svd_V_r_,n_r,grid_k_p_,svd_polyval_V_r_)
      end if 
      end if 

      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_v
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') ' setting n_transf: ' , n_transf
     $           ,' = n_delta_v = ' , n_delta_v 
         end if 
      else
         n_transf = n_svd_l
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') ' setting n_transf: ' , n_transf ,
     $           ' = n_svd_l = ' , n_svd_l
         end if 
      end if 

      allocate(S_p_(0:1+n_A-1))
      call cs1_c16(n_A,S_p_)
      allocate(S_q_(0:1+n_A-1))
      call cs1_c16(n_A,S_q_)
      allocate(S_p_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,S_p_omp__)
      allocate(S_q_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,S_q_omp__)
      allocate(M_p_(0:1+n_A-1))
      call cs1_c16(n_A,M_p_)
      allocate(M_q_(0:1+n_A-1))
      call cs1_c16(n_A,M_q_)
      allocate(M_p_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,M_p_omp__)
      allocate(M_q_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,M_q_omp__)
      allocate(CTF_p_(0:1+n_A-1))
      call cs1_c16(n_A,CTF_p_)
      allocate(CTF_q_(0:1+n_A-1))
      call cs1_c16(n_A,CTF_q_)
      allocate(CTF_p_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,CTF_p_omp__)
      allocate(CTF_q_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,CTF_q_omp__)
      allocate(Z_p_(0:1+n_A-1))
      call cs1_c16(n_A,Z_p_)
      allocate(Z_q_(0:1+n_A-1))
      call cs1_c16(n_A,Z_q_)
      allocate(Z_p_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,Z_p_omp__)
      allocate(Z_q_omp__(0:1+n_A*n_omp_sub__in-1))
      call cs1_c16(n_A*n_omp_sub__in,Z_q_omp__)
      d_memory_estimate = d_memory_estimate + 16.0d0*n_A*(8.0d0 + 8.0d0
     $     *n_omp_sub__in)
      allocate(ZZ__(0:1+n_delta_v*n_gamma_z-1))
      call cs1_c16(n_delta_v*n_gamma_z,ZZ__)
      allocate(ZZ_omp__(0:1+n_delta_v*n_gamma_z*n_omp_sub__in
     $     -1))
      call cs1_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__)
      d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_v*1.0d0
     $     *n_gamma_z*(1.0d0 + n_omp_sub__in)
      allocate(ZZ_sub_(0:1+n_w_max-1))
      call cs1_c16(n_w_max,ZZ_sub_)
      allocate(ZZ_sub_omp__(0:1+n_w_max*n_omp_sub__in-1))
      call cs1_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__)
      d_memory_estimate = d_memory_estimate + 16.0d0*n_w_max*(1.0d0 +
     $     n_omp_sub__in)
      allocate(CTF_R_S_sub__(0:1+n_gamma_z-1))
      call cs1_c16(n_gamma_z,CTF_R_S_sub__)
      allocate(CTF_R_S_sub_omp__(0:1+n_gamma_z*n_omp_sub__in-1))
      call cs1_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__)
      d_memory_estimate = d_memory_estimate + 16.0d0*n_gamma_z*(1.0d0 +
     $     n_omp_sub__in)

      allocate(n_S_0_per_(0:1+n_S-1))
      call cs1_i4(n_S,n_S_0_per_)
      allocate(n_S_0_sum_(0:1+n_S-1))
      call cs1_i4(n_S,n_S_0_sum_)
      d_memory_estimate = d_memory_estimate + 4.0d0*2.0d0*n_S
      n_S_0_sub_use = min(n_S_0_sub__in,n_S)
      call block_0(verbose-2,n_S_0_sub_use,n_S,n_S_0_per_,n_S_0_sum_
     $     ,n_S_0_sub_use,nS_0_per_min,nS_0_per_max)

      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_gamma_z*1.0d0*n_CTF*1.0d0*nS_0_per_max
     $     *16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' CTF_R_S__ requires ' , d_mem*1.0d
     $        -6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if 
      if (flag_memory_estimate.eqv..false.) then
         allocate(CTF_R_S__(0:1+n_gamma_z*n_CTF*nS_0_per_max-1))
         call cs1_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__)
      end if 
      d_memory_estimate = d_memory_estimate + 16.0d0*n_gamma_z*1.0d0
     $     *n_CTF*1.0d0*nS_0_per_max

      allocate(n_S_1_per_(0:1+nS_0_per_max-1))
      call cs1_i4(nS_0_per_max,n_S_1_per_)
      allocate(n_S_1_sum_(0:1+nS_0_per_max-1))
      call cs1_i4(nS_0_per_max,n_S_1_sum_)
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
      end if 
      if (flag_memory_estimate.eqv..false.) then
         n_C_trn0 = nS_1_per_max*n_transf
         allocate(C_trn0_(0:1+nS_1_per_max*n_transf-1))
         call cs1_c16(nS_1_per_max*n_transf,C_trn0_)
         n_C_trn0_omp = nS_1_per_max*n_transf*n_omp_sub__in
         allocate(C_trn0_omp__(0:1+nS_1_per_max*n_transf*n_omp_sub__in
     $        -1))
         call cs1_c16(nS_1_per_max*n_transf*n_omp_sub__in,C_trn0_omp__)
      end if 
      d_memory_estimate = d_memory_estimate + 16.0d0*(1.0d0
     $     *nS_1_per_max)*(1.0d0*n_transf)*(1.0d0 +n_omp_sub__in)

      if (flag_tesselation.eqv..true.) then
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_gamma_z*1.0d0*n_CTF*1.0d0*nS_1_per_max
     $        *1.0d0*n_omp_sub__in*16.0d0
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' CTF_R_S_local_omp__ requires '
     $           , d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 
         if (flag_memory_estimate.eqv..false.) then
            n_CTF_R_S_local_omp = n_gamma_z*n_CTF
     $           *nS_1_per_max*n_omp_sub__in
            allocate(CTF_R_S_local_omp__(0:1+n_gamma_z*n_CTF
     $           *nS_1_per_max*n_omp_sub__in-1))
            call cs1_c16(n_gamma_z*n_CTF*nS_1_per_max*n_omp_sub__in
     $           ,CTF_R_S_local_omp__)
         end if 
         d_memory_estimate = d_memory_estimate + 16.0d0*(1.0d0
     $        *n_gamma_z)*(1.0d0*n_CTF)*(1.0d0*nS_1_per_max)*(1.0d0
     $        *n_omp_sub__in)
      end if 

      if (flag_tesselation.eqv..true.) then
         if (flag_memory_estimate.eqv..false.) then
            allocate(vp_input_omp__(0:1+3*n_omp_sub__in-1))
            call cs1_r8(3*n_omp_sub__in,vp_input_omp__)
            allocate(flag_S_use_omp__(0:1+nS_0_per_max*n_omp_sub__in-1))
            call cs1_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__)
            allocate(LT_omp__(0:1+2*nS_0_per_max*n_omp_sub__in-1))
            call cs1_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__)
            n_S_alpha_S_index_local_omp = nS_1_per_max
     $           *n_omp_sub__in
            allocate(S_alpha_S_index_local_omp__(0:1+nS_1_per_max
     $           *n_omp_sub__in-1))
            call cs1_r8(nS_1_per_max*n_omp_sub__in
     $           ,S_alpha_S_index_local_omp__)
            allocate(S_alpha_polar_a_local_omp__(0:1+nS_1_per_max
     $           *n_omp_sub__in-1))
            call cs1_r8(nS_1_per_max*n_omp_sub__in
     $           ,S_alpha_polar_a_local_omp__)
            allocate(S_alpha_azimu_b_local_omp__(0:1+nS_1_per_max
     $           *n_omp_sub__in-1))
            call cs1_r8(nS_1_per_max*n_omp_sub__in
     $           ,S_alpha_azimu_b_local_omp__)
            allocate(I_S_sample_local_omp__(0:1+nS_1_per_max
     $           *n_omp_sub__in-1))
            call cs1_i4(nS_1_per_max*n_omp_sub__in
     $           ,I_S_sample_local_omp__)
         end if                 
         d_memory_estimate = d_memory_estimate + 8.0d0*3*n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 2.0d0*nS_0_per_max
     $        *n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 4.0d0*2*nS_0_per_max
     $        *n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 8.0d0*nS_1_per_max
     $        *n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 8.0d0*nS_1_per_max
     $        *n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 8.0d0*nS_1_per_max
     $        *n_omp_sub__in
         d_memory_estimate = d_memory_estimate + 4.0d0*nS_1_per_max
     $        *n_omp_sub__in
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array O_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (S)_q, '
            write(6,'(A)') ' '
         end if                 
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*nS_0_per_max*1.0d0*n_w_max
     $        *16.0d0
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' O_S_q__ requires ' , d_mem
     $           *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 
         if (flag_memory_estimate.eqv..false.) then
            allocate(O_S_q__(0:1+n_r*nS_0_per_max*n_w_max))
            call cs1_c16(n_r*nS_0_per_max*n_w_max,O_S_q__)
         end if                 
         d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $        *nS_0_per_max*1.0d0*n_w_max
         if (flag_tesselation.eqv..true.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*nS_1_per_max*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' O_S_q_local_omp__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               n_O_S_q_local_omp = n_r*nS_1_per_max*n_w_max
     $              *n_omp_sub__in
               allocate(O_S_q_local_omp__(0:1+n_r*nS_1_per_max*n_w_max
     $              *n_omp_sub__in-1))
               call cs1_c16(n_r*nS_1_per_max*n_w_max*n_omp_sub__in
     $              ,O_S_q_local_omp__)
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *nS_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array T_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (T_{+delta_upd}(S))_q, '
            write(6,'(A)') ' where T = translation by +delta_upd.'         
            write(6,'(A)') ' '
         end if 
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_0_per_max
     $        *1.0d0*n_w_max*16.0d0
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' T_S_q__ requires ' , d_mem
     $           *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 
         if (flag_memory_estimate.eqv..false.) then
            allocate(T_S_q__(0:1+n_r*n_transf*nS_0_per_max*n_w_max))
            call cs1_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__)
         end if                 
         d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $          *n_transf*1.0d0*nS_0_per_max*1.0d0*n_w_max
         if (flag_tesselation.eqv..true.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_1_per_max
     $           *1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' T_S_q_local_omp__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               n_T_S_q_local_omp = n_r*n_transf*nS_1_per_max
     $              *n_w_max*n_omp_sub__in
               allocate(T_S_q_local_omp__(0:1+n_r*n_transf*nS_1_per_max
     $              *n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $              *n_omp_sub__in,T_S_q_local_omp__)
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*nS_1_per_max*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (Z_{rho}(S))_q, '
            write(6,'(A)') ' where Z represents the rho-side'
            write(6,'(A)') ' of the svd of the translation operator.'
            write(6,'(A)') ' '
         end if 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_0_per_max
     $           *1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' Z_S_q__ requires ' , d_mem
     $              *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(Z_S_q__(0:1+n_r*n_transf*nS_0_per_max*n_w_max))
               call cs1_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__)
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*nS_0_per_max*1.0d0*n_w_max
         if (flag_tesselation.eqv..true.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nS_1_per_max
     $           *1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' Z_S_q_local_omp__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               n_Z_S_q_local_omp = n_r*n_transf*nS_1_per_max
     $              *n_w_max*n_omp_sub__in
               allocate(Z_S_q_local_omp__(0:1+n_r*n_transf*nS_1_per_max
     $              *n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_r*n_transf*nS_1_per_max*n_w_max
     $              *n_omp_sub__in,Z_S_q_local_omp__)
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*nS_1_per_max*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
      end if 

      allocate(n_M_0_per_(0:1+n_M-1))
      call cs1_i4(n_M,n_M_0_per_)
      allocate(n_M_0_sum_(0:1+n_M-1))
      call cs1_i4(n_M,n_M_0_sum_)
      n_M_0_sub_use = min(n_M_0_sub__in,n_M)
      call block_0(verbose-2,n_M_0_sub_use,n_M,n_M_0_per_,n_M_0_sum_
     $     ,n_M_0_sub_use,nM_0_per_min,nM_0_per_max)

      if (flag_tesselation.eqv..false.) then
         allocate(polar_a_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,polar_a_est_)
         allocate(azimu_b_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,azimu_b_est_)
         allocate(gamma_z_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,gamma_z_est_)
         allocate(delta_x_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,delta_x_est_)
         allocate(delta_y_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,delta_y_est_)
         allocate(l2_norm_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,l2_norm_est_)
         allocate(ctf_ind_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,ctf_ind_est_)
         allocate(S_index_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,S_index_est_)
         allocate(M_index_est_(0:1+nM_0_per_max-1))
         call cs1_r8(nM_0_per_max,M_index_est_)
         d_memory_estimate = d_memory_estimate + 8.0d0*9.0d0
     $        *nM_0_per_max
      else                      
         allocate(polar_a_est_(0:1+n_M-1))
         call cs1_r8(n_M,polar_a_est_)
         allocate(azimu_b_est_(0:1+n_M-1))
         call cs1_r8(n_M,azimu_b_est_)
         allocate(gamma_z_est_(0:1+n_M-1))
         call cs1_r8(n_M,gamma_z_est_)
         allocate(delta_x_est_(0:1+n_M-1))
         call cs1_r8(n_M,delta_x_est_)
         allocate(delta_y_est_(0:1+n_M-1))
         call cs1_r8(n_M,delta_y_est_)
         allocate(l2_norm_est_(0:1+n_M-1))
         call cs1_r8(n_M,l2_norm_est_)
         allocate(ctf_ind_est_(0:1+n_M-1))
         call cs1_r8(n_M,ctf_ind_est_)
         allocate(S_index_est_(0:1+n_M-1))
         call cs1_r8(n_M,S_index_est_)
         allocate(M_index_est_(0:1+n_M-1))
         call cs1_r8(n_M,M_index_est_)
         d_memory_estimate = d_memory_estimate + 8.0d0*9.0d0*n_M
      end if                    
      allocate(alpha__in_(0:1+n_alpha-1))
      call cs1_r8(n_alpha,alpha__in_)
      allocate(alpha__in_omp__(0:1+n_alpha*n_omp_sub__in-1))
      call cs1_r8(n_alpha*n_omp_sub__in,alpha__in_omp__)
      d_memory_estimate = d_memory_estimate + 8.0d0*n_alpha*(1.0d0 +
     $     n_omp_sub__in)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array O_T_R_CTF_M_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each image'
            write(6,'(A)') ' T_{-delta_est}(R_{-gamma_est}(CTF.*M)), '
            write(6,'(A)') ' where T = translation by -delta_est.'
            write(6,'(A)') ' and R = rotation by -gamma_est.'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*nM_0_per_max*1.0d0*n_w_max
     $           *16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' O_T_R_CTF_M_q__ requires ' ,
     $              d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(O_T_R_CTF_M_q__(0:1+n_r*nM_0_per_max*n_w_max-1))
               call cs1_c16(n_r*nM_0_per_max*n_w_max,O_T_R_CTF_M_q__)
               n_O_T_R_CTF_M_q__ = n_r*nM_0_per_max*n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *nM_0_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' O_T_R_CTF_M_q_local_omp__ requires ',d_mem*1.0d-6
     $              ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(O_T_R_CTF_M_q_local_omp__(0:1+n_r*n_w_max
     $              *n_omp_sub__in-1))
               call cs1_c16(n_r*n_w_max*n_omp_sub__in
     $              ,O_T_R_CTF_M_q_local_omp__)
               n_O_T_R_CTF_M_q_local_omp__ = n_r*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_w_max*1.0d0*n_omp_sub__in
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array T_T_R_CTF_M_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each image'
            write(6,'(A,A)') ' T_{-delta_upd}(T_{-delta_est} ' ,
     $           ' (R_{-gamma_est}(CTF.*M))), '
            write(6,'(A)') ' where T = translation '
            write(6,'(A)') ' and R = rotation.'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nM_0_per_max
     $           *1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' T_T_R_CTF_M_q__ requires ' ,
     $              d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(T_T_R_CTF_M_q__(0:1+n_r*n_transf*nM_0_per_max
     $              *n_w_max-1))
               call cs1_c16(n_r*n_transf*nM_0_per_max*n_w_max
     $              ,T_T_R_CTF_M_q__)
               n_T_T_R_CTF_M_q__ = n_r*n_transf*nM_0_per_max*n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*nM_0_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf
     $           *1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' T_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d-6
     $              ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(T_T_R_CTF_M_q_local_omp__(0:1+n_r*n_transf
     $              *n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_r*n_transf*n_w_max*n_omp_sub__in
     $              ,T_T_R_CTF_M_q_local_omp__)
               n_T_T_R_CTF_M_q_local_omp__ = n_r*n_transf*n_w_max
     $              *n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_T_R_CTF_M_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each image'
            write(6,'(A)')
     $           ' Z(T_{-delta_est}(R_{-gamma_est}(CTF.*M))), '
            write(6,'(A)') ' where Z represents the rho-side'
            write(6,'(A)') ' of the svd of the translation operator.'
            write(6,'(A)') ' and T = translation'
            write(6,'(A)') ' and R = rotation.'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*nM_0_per_max
     $           *1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' Z_T_R_CTF_M_q__ requires ' ,
     $              d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(Z_T_R_CTF_M_q__(0:1+n_r*n_transf*nM_0_per_max
     $              *n_w_max-1))
               call cs1_c16(n_r*n_transf*nM_0_per_max*n_w_max
     $              ,Z_T_R_CTF_M_q__)
               n_Z_T_R_CTF_M_q__ = n_r*n_transf*nM_0_per_max*n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*nM_0_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*
     $           *1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' Z_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d-6
     $              ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(Z_T_R_CTF_M_q_local_omp__(0:1+n_r*n_transf
     $              *n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_r*n_transf*n_w_max*n_omp_sub__in
     $              ,Z_T_R_CTF_M_q_local_omp__)
               n_Z_T_R_CTF_M_q_local_omp__ = n_r*n_transf*n_w_max
     $              *n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_r*1.0d0
     $           *n_transf*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
      end if                    

      allocate(n_M_1_per_(0:1+nM_0_per_max-1))
      call cs1_i4(nM_0_per_max,n_M_1_per_)
      allocate(n_M_1_sum_(0:1+nM_0_per_max-1))
      call cs1_i4(nM_0_per_max,n_M_1_sum_)
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
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Allocating array S_T_T_R_CTF_M_q__ to hold '
            write(6,'(A)') ' bessel-coefficients for the product of '
            write(6,'(A)') ' O_T_R_CTF_M_q__ and T_S_q__ '
            write(6,'(A)') ' integrated over k (i.e., nr).'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q__(0:1+n_transf*nS_1_per_max
     $              *nM_1_per_max*n_w_max-1))
               call cs1_c16(n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $              ,S_T_T_R_CTF_M_q__)
               n_S_T_T_R_CTF_M_q__ = n_transf*nS_1_per_max*nM_1_per_max
     $              *n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_transf
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_T_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q_local_omp__(0:1+n_transf
     $              *nS_1_per_max*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_T_T_R_CTF_M_q_local_omp__)
               n_S_T_T_R_CTF_M_q_local_omp__ = n_transf*nS_1_per_max
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_transf
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Allocating array S_T_T_R_CTF_M_q__ to hold '
            write(6,'(A)') ' bessel-coefficients for the product of '
            write(6,'(A)') ' T_T_R_CTF_M_q__ and O_S_q__ '
            write(6,'(A)') ' integrated over k (i.e., nr).'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q__(0:1+nS_1_per_max*n_transf
     $              *nM_1_per_max*n_w_max-1))
               call cs1_c16(nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $              ,S_T_T_R_CTF_M_q__)
               n_S_T_T_R_CTF_M_q__ = nS_1_per_max*n_transf*nM_1_per_max
     $              *n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max
     $           *1.0d0*n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_T_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q_local_omp__(0:1+nS_1_per_max
     $              *n_transf*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_T_T_R_CTF_M_q_local_omp__)
               n_S_T_T_R_CTF_M_q_local_omp__ = nS_1_per_max*n_transf
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max
     $           *1.0d0*n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_S_svdd_ to store '
            write(6,'(A)') ' displacement-operator for the delta-side '
            write(6,'(A)') ' of the svd-expansion applied to S. '
            write(6,'(A)') ' '
         end if
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_svd_l*1.0d0*n_delta_v
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' Z_S_svdd_ requires ' ,
     $           d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 
         allocate(Z_S_svdd_(0:1+n_svd_l*n_delta_v-1))
         call cs1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
         d_memory_estimate = d_memory_estimate + 16.0d0*n_svd_l*1.0d0
     $        *n_delta_v
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
            write(6,'(A)') ' bessel-coefficients for the product of '
            write(6,'(A)') ' O_T_R_CTF_M_q__ and Z_S_q__ '
            write(6,'(A)') ' integrated over k (i.e., nr).'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_Z_T_R_CTF_M_q__(0:1+n_transf*nS_1_per_max
     $              *nM_1_per_max*n_w_max-1))
               call cs1_c16(n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $              ,S_Z_T_R_CTF_M_q__)
               n_S_Z_T_R_CTF_M_q__ = n_transf*nS_1_per_max*nM_1_per_max
     $              *n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_transf
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_transf*1.0d0*nS_1_per_max*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_Z_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_Z_T_R_CTF_M_q_local_omp__(0:1+n_transf
     $              *nS_1_per_max*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_transf*nS_1_per_max*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_Z_T_R_CTF_M_q_local_omp__)
               n_S_Z_T_R_CTF_M_q_local_omp__ = n_transf*nS_1_per_max
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_transf
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in
         end if 
         if (verbose.gt.1) then
            write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
            write(6,'(A)') ' to hold product of '
            write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_S_svdd_. '
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_delta_v*1.0d0*nS_1_per_max
     $           *1.0d0*nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q__(0:1+n_delta_v*nS_1_per_max
     $              *nM_1_per_max*n_w_max-1))
               call cs1_c16(n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $              ,S_T_T_R_CTF_M_q__)
               n_S_T_T_R_CTF_M_q__ = n_delta_v*nS_1_per_max
     $              *nM_1_per_max*n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_v
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_delta_v*1.0d0*nS_1_per_max
     $           *1.0d0*nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in
     $           *16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_T_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q_local_omp__(0:1+n_delta_v
     $              *nS_1_per_max*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_T_T_R_CTF_M_q_local_omp__)
               n_S_T_T_R_CTF_M_q_local_omp__ = n_delta_v*nS_1_per_max
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_v
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_M_svdd_ to store '
            write(6,'(A)') ' displacement-operator for the delta-side '
            write(6,'(A)') ' of the svd-expansion applied to M. '
            write(6,'(A)') ' '
         end if
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_svd_l*1.0d0*n_delta_v
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' Z_M_svdd_ requires ' ,
     $           d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 
         allocate(Z_M_svdd_(0:1+n_svd_l*n_delta_v-1))
         call cs1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
         d_memory_estimate = d_memory_estimate + 16.0d0*n_svd_l*1.0d0
     $        *n_delta_v
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
            write(6,'(A)') ' bessel-coefficients for the product of '
            write(6,'(A)') ' Z_T_R_CTF_M_q__ and O_S_q__ '
            write(6,'(A)') ' integrated over k (i.e., nr).'
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_Z_T_R_CTF_M_q__(0:1+nS_1_per_max*n_transf
     $              *nM_1_per_max*n_w_max-1))
               call cs1_c16(nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $              ,S_Z_T_R_CTF_M_q__)
               n_S_Z_T_R_CTF_M_q__ = nS_1_per_max*n_transf*nM_1_per_max
     $              *n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max
     $           *1.0d0*n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*nS_1_per_max*1.0d0*n_transf*1.0d0
     $           *nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_Z_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_Z_T_R_CTF_M_q_local_omp__(0:1+nS_1_per_max
     $              *n_transf*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(nS_1_per_max*n_transf*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_Z_T_R_CTF_M_q_local_omp__)
               n_S_Z_T_R_CTF_M_q_local_omp__ = nS_1_per_max*n_transf
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*nS_1_per_max
     $           *1.0d0*n_transf*1.0d0*nM_1_per_max*1.0d0*n_w_max*1.0d0
     $           *n_omp_sub__in
         end if 
         if (verbose.gt.1) then
            write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
            write(6,'(A)') ' to hold product of '
            write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_M_svdd_. '
            write(6,'(A)') ' '
         end if
         if (flag_tesselation.eqv..false.) then
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_delta_v*1.0d0*nS_1_per_max
     $           *1.0d0*nM_1_per_max*1.0d0*n_w_max*16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires '
     $              ,d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q__(0:1+n_delta_v
     $              *nS_1_per_max*nM_1_per_max*n_w_max-1))
               call cs1_c16(n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $              ,S_T_T_R_CTF_M_q__)
               n_S_T_T_R_CTF_M_q__ = n_delta_v*nS_1_per_max
     $              *nM_1_per_max*n_w_max
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_v
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
         else 
            d_mem = 0.0d0
            d_mem = d_mem + 1.0d0*n_delta_v*1.0d0*nS_1_per_max
     $           *1.0d0*nM_1_per_max*1.0d0*n_w_max*1.0d0*n_omp_sub__in
     $           *16.0d0
            if (verbose.gt.1 .or. verbose_mem.gt.0) then
               write (6,'(A,2(F16.8,A))')
     $              ' S_T_T_R_CTF_M_q_local_omp__ requires ' ,d_mem*1.0d
     $              -6 ,' MB, ' , d_mem*1.0d-9 , ' GB'
            end if              
            if (flag_memory_estimate.eqv..false.) then
               allocate(S_T_T_R_CTF_M_q_local_omp__(0:1+n_delta_v
     $              *nS_1_per_max*nM_1_per_max*n_w_max*n_omp_sub__in-1))
               call cs1_c16(n_delta_v*nS_1_per_max*nM_1_per_max*n_w_max
     $              *n_omp_sub__in,S_T_T_R_CTF_M_q_local_omp__)
               n_S_T_T_R_CTF_M_q_local_omp__ = n_delta_v*nS_1_per_max
     $              *nM_1_per_max*n_w_max*n_omp_sub__in
            end if              
            d_memory_estimate = d_memory_estimate + 16.0d0*n_delta_v
     $           *1.0d0*nS_1_per_max*1.0d0*nM_1_per_max*1.0d0*n_w_max
     $           *1.0d0*n_omp_sub__in
         end if 
      end if                    

      allocate(n_S_9_per_(0:1+n_omp_sub__in-1))
      call cs1_i4(n_omp_sub__in,n_S_9_per_)
      allocate(n_S_9_sum_(0:1+n_omp_sub__in-1))
      call cs1_i4(n_omp_sub__in,n_S_9_sum_)
      allocate(n_M_9_per_(0:1+n_omp_sub__in-1))
      call cs1_i4(n_omp_sub__in,n_M_9_per_)
      allocate(n_M_9_sum_(0:1+n_omp_sub__in-1))
      call cs1_i4(n_omp_sub__in,n_M_9_sum_)
      d_memory_estimate = d_memory_estimate + 4.0d0*4.0d0*n_omp_sub__in

      n_S_use = 0
      n_S_use_sum = 0
      allocate(n_S_use_sum_(0:1+n_omp_sub__in-1))
      call cs1_i4(n_omp_sub__in,n_S_use_sum_)
      
      if (flag_memory_estimate.eqv..false.) then
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $           nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $           nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         end if 
      end if 

      if (flag_memory_estimate.eqv..true.) then
         if (verbose.gt.1 .or. verbose_mem.gt.0) then
            write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $           nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $           nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         end if 
         goto 10
      end if 




      if (svd_calculation_type.eq.1) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Depending on the value of flag_RTRT_vs_RTTR,'
         write(6,'(A)') ' either set up displacement-operator Z_S_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to S,'
         write(6,'(A)') ' or set up displacement-operator Z_M_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to M.'
         write(6,'(A)') ' '
      end if 
c$$$      svd_d_max = 0.0d0;
c$$$      do ndv=0,n_delta_v-1
c$$$         delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
c$$$         if (svd_d_max.lt.delta) then
c$$$            svd_d_max = delta
c$$$         end if 
c$$$      enddo 
c$$$      if (verbose.gt.0) then
c$$$         write(6,'(A,F8.6)') ' setting svd_d_max: ' , svd_d_max
c$$$      end if 
      timing_tic = omp_get_wtime()
      if (flag_RTRT_vs_RTTR.eqv..true.) then
         call cl1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
c$$$         call innerproduct_q_k_svdd_3(flag_RTRT_vs_RTTR,svd_d_max
c$$$     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_
c$$$     $        ,delta_y_,Z_S_svdd_)
         call innerproduct_q_k_svdd_4(flag_RTRT_vs_RTTR,svd_d_max
     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $        ,n_delta_v,delta_x_ ,delta_y_,Z_S_svdd_)
c$$$      %%%%%%%%
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 2
         MDA_d_(0) = n_svd_l
         MDA_d_(1) = n_delta_v
         write(MDA_string,'(A)') './dir_mda/Z_S_svdd_.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_svdd_,MDA_string)
      end if                    
c$$$      %%%%%%%%
      end if 
      if (flag_RTRT_vs_RTTR.eqv..false.) then
         call cl1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
c$$$         call innerproduct_q_k_svdd_3(flag_RTRT_vs_RTTR,svd_d_max
c$$$     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_
c$$$     $        ,delta_y_,Z_M_svdd_)
         call innerproduct_q_k_svdd_4(flag_RTRT_vs_RTTR,svd_d_max
     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $        ,n_delta_v,delta_x_ ,delta_y_,Z_M_svdd_)
c$$$      %%%%%%%%
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 2
         MDA_d_(0) = n_svd_l
         MDA_d_(1) = n_delta_v
         write(MDA_string,'(A)') './dir_mda/Z_M_svdd_.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,Z_M_svdd_,MDA_string)
      end if                    
c$$$      %%%%%%%%
      end if 
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.1) then
         write(6,'(A,F8.5)')
     $        ' finished ti8_define_svdd_3: total time '
     $        ,timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_v*1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' ti8_define_svdd_3 total Gnumps: '
     $        ,timing_tmp
         write(6,'(A)') ' '
      end if 
c$$$      %%%%%%%%
      end if 





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
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_A)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_1d:'
     $        ,' total_time ',timing_toc-timing_tic
      end if 
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))



      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
      end if 
      allocate(fftw_plan_frwd__(0:n_r*n_omp_sub__in-1))
      allocate(fftw_plan_back__(0:n_r*n_omp_sub__in-1))
      allocate(fftw_in1__(0:n_A*n_omp_sub__in-1))
      allocate(fftw_out__(0:n_A*n_omp_sub__in-1))
      timing_tic = omp_get_wtime()
      do nomp_sub=0,n_omp_sub__in-1
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
         enddo 
      enddo 
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_A)*(1.0d0*n_omp_sub__in)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_omp:'
     $        ,' total_time ',timing_toc-timing_tic
      end if 



      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_many_plan for omp sub-blocks'
         write(6,'(A)') ' '
      end if
      fpm_rank = 1
      fpm_howmany = min(fpm_howmany_max,n_transf*nS_1_per_max
     $     *nM_1_per_max)
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*fpm_rank*1.0d0*n_omp_sub__in*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_n__ requires ' , d_mem*1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if 
      allocate(fpm_n__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_n__(nomp_sub) = n_w_max
      enddo 
      allocate(fpm_inembed__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_inembed__(nomp_sub) = n_w_max
      enddo 
      allocate(fpm_onembed__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_onembed__(nomp_sub) = n_w_max
      enddo 
      fpm_istride = 1
      fpm_idist = n_w_max
      fpm_ostride = 1
      fpm_odist = n_w_max
      n_fpm = fpm_howmany*fpm_n__(0)
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_fpm*1.0d0*n_omp_sub__in*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_in1__ requires ' , d_mem*1.0d
     $        -6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if 
      allocate(fpm_in1__(0:n_fpm*n_omp_sub__in-1))
      allocate(fpm_out__(0:n_fpm*n_omp_sub__in-1))     
      allocate(fpm_back_(0:n_omp_sub__in-1))
      timing_tic = omp_get_wtime()
c$$$      We do not actually use the forward plan
      do nomp_sub=0,n_omp_sub__in-1      
         call dfftw_plan_many_dft(fpm_back_(nomp_sub),fpm_rank
     $        ,fpm_n__(nomp_sub),fpm_howmany,fpm_in1__(n_fpm*nomp_sub)
     $        ,fpm_inembed__(nomp_sub) ,fpm_istride,fpm_idist
     $        ,fpm_out__(n_fpm*nomp_sub),fpm_onembed__(nomp_sub)
     $        ,fpm_ostride ,fpm_odist ,FFTW_BACKWARD,FFTW_MEASURE)
      enddo 
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_fpm)*(1.0d0*n_omp_sub__in)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_many:'
     $        ,' total_time ',timing_toc-timing_tic
      end if 



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 
      do nS_0_sub=0,n_S_0_sub_use-1
         nS_0_per = n_S_0_per_(nS_0_sub)
         nS_0_sum = n_S_0_sum_(nS_0_sub)
         n_S_9_sub_use = min(n_omp_sub__in,nS_0_per)



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 
         call block_0(verbose-2,n_S_9_sub_use,nS_0_per,n_S_9_per_
     $        ,n_S_9_sum_,n_S_9_sub_use,nS_9_per_min,nS_9_per_max)



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 



      if (verbose.gt.1) then
         write(6,'(A)') ' Clearing array CTF_R_S_ to hold'
         write(6,'(A)') ' innerproducts for each CTF-S pair.'
         write(6,'(A)') ' '
      end if
      call cl1_c16(n_gamma_z*n_CTF*nS_0_per,CTF_R_S__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate CTF_R_S_ for each ctf-S pair.'
         write(6,'(A)') ' CTF_R_S = || CTF .* R(S) ||.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)')
     $        ' CTF_R_S_(ngz + n_gamma_z*(nctf + n_CTF*ns)) '
         write(6,'(A)') ' is equal to the l2_norm (not squared)'
         write(6,'(A)') ' of the pointwise product of CTF and R(S),'
         write(6,'(A)') ' where CTF is the nctf-th ctf-function, '
         write(6,'(A)') ' R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' and S is the ns-th template, '
         write(6,'(A)') ' at position '
         write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S). '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      if (n_S_9_sub_use.gt.1) then
c     $OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,n_r_tmp,n_A_tmp,n_X_tmp)
c     $OMP DO 
      do nS_9_sub=0,n_S_9_sub_use-1
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_A_tmp = nS_9_sub*n_A
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_CTF_R_S_3'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_CTF_R_S_3(verbose-1,n_gamma_z
     $        ,gamma_z_,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_,n_w_,n_A
     $        ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_9_per
     $        ,I_S_sample_(nS_0_sum + nS_9_sum),ld_S ,S_k_p__
     $        ,CTF_p_omp__(n_A_tmp) ,CTF_q_omp__(n_A_tmp) ,n_CTF ,ld_CTF
     $        ,CTF_k_p__ ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        ,Z_q_omp__(n_A_tmp))
      enddo                     
c     $OMP END DO
c     $OMP END PARALLEL
      else 
         nS_9_sub = 0
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_A_tmp = nS_9_sub*n_A
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_CTF_R_S_3'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_CTF_R_S_3(verbose-1,n_gamma_z
     $        ,gamma_z_,fftw_plan_frwd__(n_r_tmp)
     $        ,fftw_plan_back__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_,n_w_,n_A
     $        ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_9_per
     $        ,I_S_sample_(nS_0_sum + nS_9_sum),ld_S ,S_k_p__
     $        ,CTF_p_omp__(n_A_tmp) ,CTF_q_omp__(n_A_tmp) ,n_CTF ,ld_CTF
     $        ,CTF_k_p__ ,CTF_R_S__(nS_9_sum*n_gamma_z*n_CTF)
     $        ,Z_q_omp__(n_A_tmp))
      end if 
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (n_A*1.0d0)*(n_CTF*1.0d0)*(nS_0_per*1.0d0)
      timing_total_CTF_R_S = timing_total_CTF_R_S + timing_tot
      gnump_total_CTF_R_S = gnump_total_CTF_R_S + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.5)')
     $        ' finished calculating CTF_R_S_ for each ctf-S pair. ' ,
     $        ' total time ' , timing_tot
         write(6,'(A,F8.4)') ' CTF_R_S_ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if                    



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 



      if (flag_RTRT_vs_RTTR.eqv..false.) then
      call cl1_c16(n_r*nS_0_per*n_w_max,O_S_q__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate bessel coefficients of S.'
         write(6,'(A)') ' S = template in k-space polar coord.'
         write(6,'(A)') ' '
         write(6,'(A,A)') ' O_S_q__(nr + n_r*(ns + nS_0_per*nw)) '
         write(6,'(A)') ' is equal to the bessel-coefficients'
         write(6,'(A)') ' conjg(S_q_(ic)) '
         write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
         write(6,'(A)') ' where S_q_ derived from '
         write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S).'
         write(6,'(A)') ' '
      end if                    
      timing_tic = omp_get_wtime()
      if (n_S_9_sub_use.gt.1) then
c     $OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,n_r_tmp,n_A_tmp,n_X_tmp)
c     $OMP DO 
      do nS_9_sub=0,n_S_9_sub_use-1
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_A_tmp = nS_9_sub*n_A
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_O_S_q_1'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_O_S_q_1(verbose-1
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_,n_w_,n_A
     $        ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_9_per
     $        ,I_S_sample_(nS_0_sum + nS_9_sum) ,ld_S ,S_k_p__,nS_0_per
     $        ,O_S_q__(n_r*nS_9_sum))
      enddo 
c     $OMP END DO
c     $OMP END PARALLEL
      else 
         nS_9_sub = 0
         nS_9_per = n_S_9_per_(nS_9_sub)
         nS_9_sum = n_S_9_sum_(nS_9_sub)
         n_r_tmp = nS_9_sub*n_r
         n_A_tmp = nS_9_sub*n_A
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' ti8_build_O_S_q_1'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call ti8_build_O_S_q_1(verbose-1
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_k_p_,n_w_,n_A
     $        ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_9_per
     $        ,I_S_sample_(nS_0_sum + nS_9_sum) ,ld_S ,S_k_p__,nS_0_per
     $        ,O_S_q__(n_r*nS_9_sum))
      end if 
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (n_A*1.0d0)*(nS_0_per*1.0d0)
      timing_total_O_S_q = timing_total_O_S_q + timing_tot
      gnump_total_O_S_q = gnump_total_O_S_q + gnump_tot
      timing_tmp = gnump_tot/timing_tot/1e9
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.5)')
     $        ' finished calculating O_S_q__ for each S. ' ,
     $        ' total time ' , timing_tot
         write(6,'(A,F8.4)') ' O_S_q__ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if                    
         flag_MDA = .false.
         if (flag_MDA.eqv..true.) then
            MDA_n_d = 3
            MDA_d_(0) = n_r
            MDA_d_(1) = nS_0_per
            MDA_d_(2) = n_w_max
            write(MDA_string,'(A,1(A,I0),A)') './dir_mda/O_S_q__' ,
     $           'nS0_' , nS_0_sub,'.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,O_S_q__,MDA_string)
         end if 
      end if 



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 



      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.)) then
         call cl1_c16(n_r*n_transf*nS_0_per*n_w_max,T_S_q__)
         if (verbose.gt.1) then
            write(6,'(A)') ' Calculate bessel coefficients of T(S).'
            write(6,'(A)') ' T = translation by delta_x,delta_y. '
            write(6,'(A)') ' S = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A)') ' T_S_q__(nr + n_r*(ndv + '
     $           ,'n_delta_v*(ns + nS_0_per*nw))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' conjg(S_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
            write(6,'(A)') ' where S_q_ = T(S), with '
            write(6,'(A)') ' T <-- delta_x_(ndv) and delta_y_(ndv)'
            write(6,'(A)') ' and S derived from '
            write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S).'
            write(6,'(A)') ' '
         end if 
         timing_tic = omp_get_wtime()
         if (n_S_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,n_r_tmp,n_A_tmp)
c$OMP DO 
         do nS_9_sub=0,n_S_9_sub_use-1
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            n_r_tmp = nS_9_sub*n_r
            n_A_tmp = nS_9_sub*n_A
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_build_T_S_q_4'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_T_S_q_4(verbose-1,n_delta_v
     $           ,delta_x_,delta_y_,fftw_plan_frwd__(n_r_tmp)
     $           ,fftw_in1__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r
     $           ,grid_k_p_ ,n_w_,n_A,S_p_omp__(n_A_tmp)
     $           ,S_q_omp__(n_A_tmp),nS_9_per ,I_S_sample_(nS_0_sum +
     $           nS_9_sum) ,ld_S ,S_k_p__,nS_0_per ,T_S_q__(n_r
     $           *n_delta_v*nS_9_sum))
         enddo 
c$OMP END DO
c$OMP END PARALLEL
         else 
            nS_9_sub = 0
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            n_r_tmp = nS_9_sub*n_r
            n_A_tmp = nS_9_sub*n_A
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_build_T_S_q_4'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_T_S_q_4(verbose-1,n_delta_v
     $           ,delta_x_,delta_y_,fftw_plan_frwd__(n_r_tmp)
     $           ,fftw_in1__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r
     $           ,grid_k_p_ ,n_w_,n_A,S_p_omp__(n_A_tmp)
     $           ,S_q_omp__(n_A_tmp),nS_9_per ,I_S_sample_(nS_0_sum +
     $           nS_9_sum) ,ld_S ,S_k_p__,nS_0_per ,T_S_q__(n_r
     $           *n_delta_v*nS_9_sum))
         end if 
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         gnump_tot = (n_A*1.0d0)*(n_delta_v*1.0d0)
     $        *(nS_0_per*1.0d0)
         timing_total_T_S_q = timing_total_T_S_q + timing_tot
         gnump_total_T_S_q = gnump_total_T_S_q + gnump_tot
         timing_tmp = gnump_tot/timing_tot/1e9
         if (verbose_timing.gt.0) then
            write(6,'(A,A,F8.5)')
     $           ' finished calculating T_S_q__ for each S. ' ,
     $           ' total time ' , timing_tot
            write(6,'(A,F8.4)') ' T_S_q__ total Gnumps: ' , timing_tmp
            write(6,'(A)') ' '
         end if 
         flag_MDA = .false.
         if (flag_MDA.eqv..true.) then
            MDA_n_d = 4
            MDA_d_(0) = n_r
            MDA_d_(1) = n_transf
            MDA_d_(2) = nS_0_per
            MDA_d_(3) = n_w_max
            write(MDA_string,'(A,1(A,I0),A)') './dir_mda/T_S_q__' ,
     $           'nS0_' , nS_0_sub,'.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,T_S_q__,MDA_string)
         end if 
      end if 




      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 



      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.)) then
         call cl1_c16(n_r*n_transf*nS_0_per*n_w_max,Z_S_q__)
         if (verbose.gt.1) then
            write(6,'(A)') ' Calculate bessel coefficients of Z(S).'
            write(6,'(A)') ' Z = svd of translation operator. '
            write(6,'(A)') ' S = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A)') ' Z_S_q__(nr + n_r*(nl + '
     $           ,'n_svd_l*(ns + nS_0_per*nw))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' conjg(S_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_csum_(nr),'
            write(6,'(A)') ' where S_q_ = Z(S), with '
            write(6,'(A)') ' Z representing the rho-factor '
            write(6,'(A)') ' of the translation-operator'
            write(6,'(A)') ' associated with delta,'
            write(6,'(A)') ' and S is derived from '
            write(6,'(A)') ' S_k_p__(I_S_sample_(nS_0_sum + ns)*ld_S).'
            write(6,'(A)') ' '
         end if 
c$$$         svd_r_max = 0.0d0
c$$$         svd_r_max = grid_k_p_(n_r-1)
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,F8.2)') ' setting svd_r_max: ' , svd_r_max
c$$$         end if                 
         timing_tic = omp_get_wtime()
         if (n_S_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum,n_r_tmp,n_A_tmp)
c$OMP DO 
         do nS_9_sub=0,n_S_9_sub_use-1
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            n_r_tmp = nS_9_sub*n_r
            n_A_tmp = nS_9_sub*n_A
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_build_Z_S_q_4'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_Z_S_q_5(verbose-1,svd_r_max,n_svd_r
     $           ,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,svd_polyval_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) ,fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp),n_r ,grid_k_p_,n_w_,n_A
     $           ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp) ,nS_9_per
     $           ,I_S_sample_(nS_0_sum + ns_9_sum),ld_S ,S_k_p__
     $           ,nS_0_per ,Z_S_q__(n_r*n_svd_l*nS_9_sum))
         enddo 
c$OMP END DO
c$OMP END PARALLEL
         else 
            nS_9_sub = 0
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            n_r_tmp = nS_9_sub*n_r
            n_A_tmp = nS_9_sub*n_A
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_build_Z_S_q_4'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_build_Z_S_q_5(verbose-1,svd_r_max,n_svd_r
     $           ,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,svd_polyval_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) ,fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp),n_r ,grid_k_p_,n_w_,n_A
     $           ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp) ,nS_9_per
     $           ,I_S_sample_(nS_0_sum + ns_9_sum),ld_S ,S_k_p__
     $           ,nS_0_per ,Z_S_q__(n_r*n_svd_l*nS_9_sum))
         end if 
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         gnump_tot = (n_A*1.0d0)*(n_svd_l*1.0d0)*(nS_0_per*1.0d0)
         timing_total_Z_S_q = timing_total_Z_S_q + timing_tot
         gnump_total_Z_S_q = gnump_total_Z_S_q + gnump_tot
         timing_tmp = gnump_tot/timing_tot/1e9
         if (verbose_timing.gt.0) then
            write(6,'(A,A,F8.5)')
     $           ' finished calculating Z_S_q__ for each S. ' ,
     $           ' total time ' , timing_tot
            write(6,'(A,F8.4)') ' Z_S_q__ total Gnumps: ' , timing_tmp
            write(6,'(A)') ' '
         end if 
         flag_MDA = .false.
         if (flag_MDA.eqv..true.) then
            MDA_n_d = 4
            MDA_d_(0) = n_r
            MDA_d_(1) = n_transf
            MDA_d_(2) = nS_0_per
            MDA_d_(3) = n_w_max
            write(MDA_string,'(A,1(A,I0),A)') './dir_mda/Z_S_q__' ,
     $           'nS0_' , nS_0_sub,'.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_q__,MDA_string)
         end if 
      end if 



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 
         if (tesselation_distance_req.ge.2.0d0) then
            call ti8_global_0(



     $     verbose 
     $     ,flag_memory_estimate 
     $     ,d_memory_estimate 
     $     ,rseed 
     $     ,n_k_cur 
     $     ,n_polar_a_ 
     $     ,grid_k_p_ 
     $     ,half_diameter_x_c 
     $     ,n_S 
     $     ,I_S_sample_ 
     $     ,ld_S 
     $     ,S_k_p__ 
     $     ,tesselation_distance_req 
     $     ,n_LT_add 
     $     ,n_LT_ref 
     $     ,S_alpha_S_index_ 
     $     ,S_alpha_polar_a_ 
     $     ,S_alpha_azimu_b_ 
     $     ,n_M 
     $     ,I_M_sample_ 
     $     ,ld_M 
     $     ,M_k_p__ 
     $     ,n_CTF 
     $     ,ld_CTF 
     $     ,CTF_k_p__ 
     $     ,alpha_est__ 
     $     ,alpha_update_f 
     $     ,flag_MS_vs_SM 
     $     ,n_SM_max 
     $     ,n_SM_ 
     $     ,alpha_SM__ 
     $     ,n_MS_max 
     $     ,n_MS_ 
     $     ,alpha_MS__ 
     $     ,n_pixels_in 
     $     ,displacement_max 
     $     ,n_delta_v 
     $     ,delta_x_ 
     $     ,delta_y_ 
     $     ,n_gamma_z 
     $     ,svd_calculation_type 
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
     $     ,eps_svd 
     $     ,flag_RTRT_vs_RTTR 
     $     ,fpm_howmany_max 
     $     ,n_omp_sub__in 
     $     ,n_S_0_sub__in 
     $     ,n_S_1_sub__in 
     $     ,n_M_0_sub__in 
     $     ,n_M_1_sub__in 
     $     ,C_M_ 



c$$$      %%%%%%%%%%%%%%%%
c$$$      function outputs
c$$$      %%%%%%%%%%%%%%%%
c$$$     $ ,max_i4_f,sum_i4_f 
c$$$     $ ,max_r8_f
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%
     $ ,flag_memory_checkset 
     $ ,flag_tesselation
     $ ,flag_tesselation_ref
     $ ,n_point_init 
     $ ,n_point,npoint 
     $ ,S_L_ 
     $ ,tradius_min 
     $ ,nl_max_init,nm_sum_init,ll_sum_init 
     $ ,nl_max,nm_sum,ll_sum 
     $ ,n_T_root_base,nT_root_base 
     $ ,parity_input 
c$$$      T_ tesselation tree
     $ ,T_nl_ 
     $ ,T_vm_ 
     $ ,T_tr_ 
     $ ,T_ll_ 
     $ ,T_lf_ 
     $ ,T_c0_ 
     $ ,T_c1_ 
     $ ,T_c2_ 
     $ ,T_c3_ 
     $ ,T_ls_ 
     $ ,T_LT_ 
c$$$      T_ roots
     $ ,T_root_base_ 
c$$$      %%%%%%%%%%%%%%%%
c$$$      bookkeeping tesselation search across threads
c$$$      %%%%%%%%%%%%%%%%
     $ ,tesselation_distance_req_use 
     $ ,n_SM_use_local_ 
     $ ,n_SM_use 
     $ ,vp_input_omp__ 
     $ ,p_vp_input 
     $ ,flag_S_use_omp__ 
     $ ,p_flag_S_use 
     $ ,n_LT,n_S_use,n_S_use_sum
     $ ,n_S_use_sum_ 
     $ ,n_add,n_ref,nref,nLT,ns_local
     $ ,n_pass,n_pass_ref
     $ ,continue_flag
     $ ,LT_omp__ 
     $ ,p_LT 
     $ ,S_alpha_S_index_local_omp__ 
     $ ,p_S_alpha_S_index_local 
     $ ,S_alpha_polar_a_local_omp__ 
     $ ,p_S_alpha_polar_a_local 
     $ ,S_alpha_azimu_b_local_omp__ 
     $ ,p_S_alpha_azimu_b_local 
     $ ,I_S_sample_local_omp__ 
     $ ,p_I_S_sample_local 
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%      
     $ ,eps_target 
     $ ,n_svd_r,n_svd_d,n_svd_l
     $ ,nsvd_r,nsvd_d,nsvd_l
     $ ,svd_r_
     $ ,svd_d_
     $ ,svd_l_
     $ ,svd_U_d_
     $ ,svd_s_
     $ ,svd_V_r_
     $ ,n_svd_max
     $ ,svd_d_max,svd_r_max 
     $ ,svd_polyval_U_d_ 
     $ ,svd_polyval_V_r_ 
     $ ,flag_warning
     $ ,R_max,K_max,delta,delta_max,n_pixels 
     $ ,Z_S_svdd_ 
     $ ,Z_M_svdd_ 
     $ ,CTF_R_S__ 
     $ ,CTF_R_S_sub__ 
     $ ,CTF_R_S_sub_omp__ 
     $ ,CTF_R_S_local_omp__ 
     $ ,p_CTF_R_S_local 
     $ ,O_S_q__ 
     $ ,O_S_q_local_omp__ 
     $ ,p_O_S_q_local 
     $ ,T_S_q__ 
     $ ,T_S_q_local_omp__ 
     $ ,p_T_S_q_local 
     $ ,Z_S_q__ 
     $ ,Z_S_q_local_omp__ 
     $ ,p_Z_S_q_local 
c$$$      pointer (p_S_p__,S_p__) 
     $ ,S_p__
c$$$      pointer (p_M_p__,M_p__) 
     $ ,M_p__
c$$$      pointer (p_CTF_p__,CTF_p__) 
     $ ,CTF_p__
     $ ,S_p_ 
     $ ,S_q_ 
     $ ,S_p_omp__ 
     $ ,S_q_omp__ 
     $ ,M_p_ 
     $ ,M_q_ 
     $ ,M_p_omp__ 
     $ ,M_q_omp__ 
     $ ,CTF_p_ 
     $ ,CTF_q_ 
     $ ,CTF_p_omp__ 
     $ ,CTF_q_omp__ 
     $ ,Z_p_ 
     $ ,Z_q_ 
     $ ,Z_p_omp__ 
     $ ,Z_q_omp__ 
     $ ,ZZ__ 
     $ ,ZZ_omp__ 
     $ ,ZZ_sub_ 
     $ ,ZZ_sub_omp__ 
     $ ,C_trn0_ 
     $ ,C_trn0_omp__ 
     $ ,O_T_R_CTF_M_q__ 
     $ ,O_T_R_CTF_M_q_local_omp__ 
     $ ,p_O_T_R_CTF_M_q_local 
     $ ,T_T_R_CTF_M_q__ 
     $ ,T_T_R_CTF_M_q_local_omp__ 
     $ ,p_T_T_R_CTF_M_q_local 
     $ ,Z_T_R_CTF_M_q__ 
     $ ,Z_T_R_CTF_M_q_local_omp__ 
     $ ,p_Z_T_R_CTF_M_q_local 
     $ ,S_T_T_R_CTF_M_q__ 
     $ ,S_T_T_R_CTF_M_q_local_omp__ 
     $ ,p_S_T_T_R_CTF_M_q_local 
     $ ,S_Z_T_R_CTF_M_q__ 
     $ ,S_Z_T_R_CTF_M_q_local_omp__ 
     $ ,p_S_Z_T_R_CTF_M_q_local 
     $ ,n_C_trn0 
     $ ,n_C_trn0_omp 
     $ ,n_CTF_R_S_local_omp 
     $ ,n_S_alpha_S_index_local_omp 
     $ ,n_O_S_q_local_omp 
     $ ,n_T_S_q_local_omp 
     $ ,n_Z_S_q_local_omp 
     $ ,n_O_T_R_CTF_M_q__ 
     $ ,n_O_T_R_CTF_M_q_local_omp__ 
     $ ,n_T_T_R_CTF_M_q__ 
     $ ,n_T_T_R_CTF_M_q_local_omp__ 
     $ ,n_Z_T_R_CTF_M_q__ 
     $ ,n_Z_T_R_CTF_M_q_local_omp__ 
     $ ,n_S_T_T_R_CTF_M_q__ 
     $ ,n_S_T_T_R_CTF_M_q_local_omp__ 
     $ ,n_S_Z_T_R_CTF_M_q__ 
     $ ,n_S_Z_T_R_CTF_M_q_local_omp__ 
c$$$      indices
c$$$     $ ,nr,na,ns,nctf,nm,nw
     $ ,nr,na,ns,nctf,nm,nw
     $ ,n_r 
     $ ,n_w_ 
     $ ,n_w_csum_ 
     $ ,n_w_max 
     $ ,n_A 
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
     $ ,alpha__in_
     $ ,alpha__in_omp__
c$$$      fftw for omp
     $ ,fftw_plan_frwd__
     $ ,fftw_plan_back__
     $ ,fftw_in1__
     $ ,fftw_out__
c$$$      fftw for local use
     $ ,fftw_plan_frwd_
     $ ,fftw_plan_back_
     $ ,fftw_in1_
     $ ,fftw_out_
c$$$      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
c$$$      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
c$$$      pointer (p_fftw_in1_last_,fftw_in1_last_)
c$$$      pointer (p_fftw_out_last_,fftw_out_last_)
     $ ,fftw_plan_frwd_last,fftw_plan_back_last
     $ ,fftw_in1_last_,fftw_out_last_
c$$$      fftw plan many (fpm) for omp
     $ ,n_transf 
     $ ,fpm_howmany
     $ ,n_fpm
c$$$     $ ,nfpm
     $ ,nfpm
     $ ,fpm_frwd_
     $ ,fpm_back_
     $ ,fpm_in1__
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
     $ ,gamma_z_ 
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
     $ ,flag_fill 
     $ ,flag_test 
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
     $           )
         else 



      if (0+verbose.gt.1) then
         write(6,'(A,F8.4)') ' tesselation_distance_req: ' ,
     $        tesselation_distance_req 
         write(6,'(A)') ' Constructing tesselation tree. '
      end if                    
      
      n_point = nS_0_per;
      call cl1_r8(3*n_point,S_L_)
      do npoint=0,n_point-1
         ns0 = nS_0_sum + npoint
         call get_angle_to_vp_(S_alpha_polar_a_(ns0)
     $        ,S_alpha_azimu_b_(ns0),S_L_(0+3*npoint))
         call normalize_r8(3,S_L_(0 + 3*npoint))
      enddo                     
      if (verbose.gt.1) then
         do npoint=0,n_point-1
            ns0 = nS_0_sum + npoint
            write(6,'(A,F8.4,F8.4,A,A,I0,1X,3F8.4)') 'S_alpha: ' ,
     $           S_alpha_polar_a_(ns0) , S_alpha_azimu_b_(ns0)
     $           , ' --> ' , ' S_L_: ', npoint , S_L_(0 + 3*npoint),
     $           S_L_(1 + 3 *npoint) , S_L_(2 + 3*npoint)
         enddo                  
      end if                    
      tradius_min = 1.0d-6
      call tesselation_get_nl_nm_ll(n_point,S_L_,tradius_min,nl_max
     $     ,nm_sum,ll_sum)
      if (verbose.gt.-1) then
         write(6,'(A,F8.6,4(A,I0))') ' tradius_min: ' , tradius_min
     $        ,' n_point: ' , n_point , ' nl_max: ' , nl_max
     $        ,' nm_sum: ' , nm_sum , ' ll_sum: ' , ll_sum
      end if                    

       call cl1_i4(nm_sum,T_nl_) 
       call cl1_r8(3*nm_sum,T_vm_) 
       call cl1_r8(nm_sum,T_tr_) 
       call cl1_i4(nm_sum,T_ll_) 
       call cl1_l2(nm_sum,T_lf_) 
       call cl1_i4(nm_sum,T_c0_) 
       call cl1_i4(nm_sum,T_c1_) 
       call cl1_i4(nm_sum,T_c2_) 
       call cl1_i4(nm_sum,T_c3_) 
       call cl1_i4(nm_sum,T_ls_) 
       call cl1_i4(ll_sum,T_LT_) 
       call cl1_i4(8,T_root_base_)

      call tesselation_index_wrapper_0(n_point,S_L_,tradius_min ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)

            call ti8_local_0(



     $     verbose 
     $     ,flag_memory_estimate 
     $     ,d_memory_estimate 
     $     ,rseed 
     $     ,n_k_cur 
     $     ,n_polar_a_ 
     $     ,grid_k_p_ 
     $     ,half_diameter_x_c 
     $     ,n_S 
     $     ,I_S_sample_ 
     $     ,ld_S 
     $     ,S_k_p__ 
     $     ,tesselation_distance_req 
     $     ,n_LT_add 
     $     ,n_LT_ref 
     $     ,S_alpha_S_index_ 
     $     ,S_alpha_polar_a_ 
     $     ,S_alpha_azimu_b_ 
     $     ,n_M 
     $     ,I_M_sample_ 
     $     ,ld_M 
     $     ,M_k_p__ 
     $     ,n_CTF 
     $     ,ld_CTF 
     $     ,CTF_k_p__ 
     $     ,alpha_est__ 
     $     ,alpha_update_f 
     $     ,flag_MS_vs_SM 
     $     ,n_SM_max 
     $     ,n_SM_ 
     $     ,alpha_SM__ 
     $     ,n_MS_max 
     $     ,n_MS_ 
     $     ,alpha_MS__ 
     $     ,n_pixels_in 
     $     ,displacement_max 
     $     ,n_delta_v 
     $     ,delta_x_ 
     $     ,delta_y_ 
     $     ,n_gamma_z 
     $     ,svd_calculation_type 
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
     $     ,eps_svd 
     $     ,flag_RTRT_vs_RTTR 
     $     ,fpm_howmany_max 
     $     ,n_omp_sub__in 
     $     ,n_S_0_sub__in 
     $     ,n_S_1_sub__in 
     $     ,n_M_0_sub__in 
     $     ,n_M_1_sub__in 
     $     ,C_M_ 



c$$$      %%%%%%%%%%%%%%%%
c$$$      function outputs
c$$$      %%%%%%%%%%%%%%%%
c$$$     $ ,max_i4_f,sum_i4_f 
c$$$     $ ,max_r8_f
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%
     $ ,flag_memory_checkset 
     $ ,flag_tesselation
     $ ,flag_tesselation_ref
     $ ,n_point_init 
     $ ,n_point,npoint 
     $ ,S_L_ 
     $ ,tradius_min 
     $ ,nl_max_init,nm_sum_init,ll_sum_init 
     $ ,nl_max,nm_sum,ll_sum 
     $ ,n_T_root_base,nT_root_base 
     $ ,parity_input 
c$$$      T_ tesselation tree
     $ ,T_nl_ 
     $ ,T_vm_ 
     $ ,T_tr_ 
     $ ,T_ll_ 
     $ ,T_lf_ 
     $ ,T_c0_ 
     $ ,T_c1_ 
     $ ,T_c2_ 
     $ ,T_c3_ 
     $ ,T_ls_ 
     $ ,T_LT_ 
c$$$      T_ roots
     $ ,T_root_base_ 
c$$$      %%%%%%%%%%%%%%%%
c$$$      bookkeeping tesselation search across threads
c$$$      %%%%%%%%%%%%%%%%
     $ ,tesselation_distance_req_use 
     $ ,n_SM_use_local_ 
     $ ,n_SM_use 
     $ ,vp_input_omp__ 
     $ ,p_vp_input 
     $ ,flag_S_use_omp__ 
     $ ,p_flag_S_use 
     $ ,n_LT,n_S_use,n_S_use_sum
     $ ,n_S_use_sum_ 
     $ ,n_add,n_ref,nref,nLT,ns_local
     $ ,n_pass,n_pass_ref
     $ ,continue_flag
     $ ,LT_omp__ 
     $ ,p_LT 
     $ ,S_alpha_S_index_local_omp__ 
     $ ,p_S_alpha_S_index_local 
     $ ,S_alpha_polar_a_local_omp__ 
     $ ,p_S_alpha_polar_a_local 
     $ ,S_alpha_azimu_b_local_omp__ 
     $ ,p_S_alpha_azimu_b_local 
     $ ,I_S_sample_local_omp__ 
     $ ,p_I_S_sample_local 
c$$$      %%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%      
     $ ,eps_target 
     $ ,n_svd_r,n_svd_d,n_svd_l
     $ ,nsvd_r,nsvd_d,nsvd_l
     $ ,svd_r_
     $ ,svd_d_
     $ ,svd_l_
     $ ,svd_U_d_
     $ ,svd_s_
     $ ,svd_V_r_
     $ ,n_svd_max
     $ ,svd_d_max,svd_r_max 
     $ ,svd_polyval_U_d_ 
     $ ,svd_polyval_V_r_ 
     $ ,flag_warning
     $ ,R_max,K_max,delta,delta_max,n_pixels 
     $ ,Z_S_svdd_ 
     $ ,Z_M_svdd_ 
     $ ,CTF_R_S__ 
     $ ,CTF_R_S_sub__ 
     $ ,CTF_R_S_sub_omp__ 
     $ ,CTF_R_S_local_omp__ 
     $ ,p_CTF_R_S_local 
     $ ,O_S_q__ 
     $ ,O_S_q_local_omp__ 
     $ ,p_O_S_q_local 
     $ ,T_S_q__ 
     $ ,T_S_q_local_omp__ 
     $ ,p_T_S_q_local 
     $ ,Z_S_q__ 
     $ ,Z_S_q_local_omp__ 
     $ ,p_Z_S_q_local 
c$$$      pointer (p_S_p__,S_p__) 
     $ ,S_p__
c$$$      pointer (p_M_p__,M_p__) 
     $ ,M_p__
c$$$      pointer (p_CTF_p__,CTF_p__) 
     $ ,CTF_p__
     $ ,S_p_ 
     $ ,S_q_ 
     $ ,S_p_omp__ 
     $ ,S_q_omp__ 
     $ ,M_p_ 
     $ ,M_q_ 
     $ ,M_p_omp__ 
     $ ,M_q_omp__ 
     $ ,CTF_p_ 
     $ ,CTF_q_ 
     $ ,CTF_p_omp__ 
     $ ,CTF_q_omp__ 
     $ ,Z_p_ 
     $ ,Z_q_ 
     $ ,Z_p_omp__ 
     $ ,Z_q_omp__ 
     $ ,ZZ__ 
     $ ,ZZ_omp__ 
     $ ,ZZ_sub_ 
     $ ,ZZ_sub_omp__ 
     $ ,C_trn0_ 
     $ ,C_trn0_omp__ 
     $ ,O_T_R_CTF_M_q__ 
     $ ,O_T_R_CTF_M_q_local_omp__ 
     $ ,p_O_T_R_CTF_M_q_local 
     $ ,T_T_R_CTF_M_q__ 
     $ ,T_T_R_CTF_M_q_local_omp__ 
     $ ,p_T_T_R_CTF_M_q_local 
     $ ,Z_T_R_CTF_M_q__ 
     $ ,Z_T_R_CTF_M_q_local_omp__ 
     $ ,p_Z_T_R_CTF_M_q_local 
     $ ,S_T_T_R_CTF_M_q__ 
     $ ,S_T_T_R_CTF_M_q_local_omp__ 
     $ ,p_S_T_T_R_CTF_M_q_local 
     $ ,S_Z_T_R_CTF_M_q__ 
     $ ,S_Z_T_R_CTF_M_q_local_omp__ 
     $ ,p_S_Z_T_R_CTF_M_q_local 
     $ ,n_C_trn0 
     $ ,n_C_trn0_omp 
     $ ,n_CTF_R_S_local_omp 
     $ ,n_S_alpha_S_index_local_omp 
     $ ,n_O_S_q_local_omp 
     $ ,n_T_S_q_local_omp 
     $ ,n_Z_S_q_local_omp 
     $ ,n_O_T_R_CTF_M_q__ 
     $ ,n_O_T_R_CTF_M_q_local_omp__ 
     $ ,n_T_T_R_CTF_M_q__ 
     $ ,n_T_T_R_CTF_M_q_local_omp__ 
     $ ,n_Z_T_R_CTF_M_q__ 
     $ ,n_Z_T_R_CTF_M_q_local_omp__ 
     $ ,n_S_T_T_R_CTF_M_q__ 
     $ ,n_S_T_T_R_CTF_M_q_local_omp__ 
     $ ,n_S_Z_T_R_CTF_M_q__ 
     $ ,n_S_Z_T_R_CTF_M_q_local_omp__ 
c$$$      indices
c$$$     $ ,nr,na,ns,nctf,nm,nw
     $ ,nr,na,ns,nctf,nm,nw
     $ ,n_r 
     $ ,n_w_ 
     $ ,n_w_csum_ 
     $ ,n_w_max 
     $ ,n_A 
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
     $ ,alpha__in_
     $ ,alpha__in_omp__
c$$$      fftw for omp
     $ ,fftw_plan_frwd__
     $ ,fftw_plan_back__
     $ ,fftw_in1__
     $ ,fftw_out__
c$$$      fftw for local use
     $ ,fftw_plan_frwd_
     $ ,fftw_plan_back_
     $ ,fftw_in1_
     $ ,fftw_out_
c$$$      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
c$$$      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
c$$$      pointer (p_fftw_in1_last_,fftw_in1_last_)
c$$$      pointer (p_fftw_out_last_,fftw_out_last_)
     $ ,fftw_plan_frwd_last,fftw_plan_back_last
     $ ,fftw_in1_last_,fftw_out_last_
c$$$      fftw plan many (fpm) for omp
     $ ,n_transf 
     $ ,fpm_howmany
     $ ,n_fpm
c$$$     $ ,nfpm
     $ ,nfpm
     $ ,fpm_frwd_
     $ ,fpm_back_
     $ ,fpm_in1__
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
     $ ,gamma_z_ 
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
     $ ,flag_fill 
     $ ,flag_test 
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
     $           )
         end if 



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 
      enddo 



      if (flag_memory_estimate.eqv..true.) then
c$$$         do nothing. ;
      end if 

      if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.1) then
         write(6,'(A)') '[entering ti8_checkset_variable_0]'
      end if 
      
      flag_memory_checkset = .true.

      if (flag_tesselation.eqv..false.) then
c$$$         do nothing
      else 
         call cxs_i4(n_M,n_SM_use_local_,'n_SM_use_local_'
     $        ,flag_memory_checkset)
         call cxs_r8(3*n_point_init,S_L_,'S_L_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_nl_,'T_nl_',flag_memory_checkset)
         call cxs_r8(3*nm_sum_init,T_vm_,'T_vm_',flag_memory_checkset)
         call cxs_r8(1*nm_sum_init,T_tr_,'T_tr_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ll_,'T_ll_',flag_memory_checkset)
         call cxs_l2(1*nm_sum_init,T_lf_,'T_lf_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c0_,'T_c0_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c1_,'T_c1_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c2_,'T_c2_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_c3_,'T_c3_',flag_memory_checkset)
         call cxs_i4(1*nm_sum_init,T_ls_,'T_ls_',flag_memory_checkset)
         call cxs_i4(1*ll_sum_init,T_LT_,'T_LT_',flag_memory_checkset)
         call cxs_i4(8,T_root_base_,'T_root_base_'
     $        ,flag_memory_checkset)
      end if 

      call cxs_i4(n_k_cur,n_w_,'n_w_',flag_memory_checkset)
      call cxs_i4(n_k_cur,n_w_csum_,'n_w_csum_',flag_memory_checkset)
      call cxs_r8(n_gamma_z,gamma_z_,'gamma_z_',flag_memory_checkset)

      if ((svd_calculation_type.eq.1)) then
         call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
         call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_max,svd_s_,8,'svd_s_',flag_memory_checkset)
         call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_ 
     $        ,'svd_polyval_U_d_',flag_memory_checkset)
         call cxs_r8(n_svd_l*n_r,svd_polyval_V_r_
     $        ,'svd_polyval_V_r_',flag_memory_checkset)
      end if 

      call cxs_c16(n_A,S_p_,'S_p_',flag_memory_checkset)
      call cxs_c16(n_A,S_q_,'S_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_p_omp__,'S_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,S_q_omp__,'S_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,M_p_,'M_p_',flag_memory_checkset)
      call cxs_c16(n_A,M_q_,'M_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_p_omp__,'M_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,M_q_omp__,'M_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,CTF_p_,'CTF_p_',flag_memory_checkset)
      call cxs_c16(n_A,CTF_q_,'CTF_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_p_omp__,'CTF_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,CTF_q_omp__,'CTF_q_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A,Z_p_,'Z_p_',flag_memory_checkset)
      call cxs_c16(n_A,Z_q_,'Z_q_',flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_p_omp__,'Z_p_omp__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_A*n_omp_sub__in,Z_q_omp__,'Z_q_omp__'
     $     ,flag_memory_checkset)

      call cxs_c16(n_delta_v*n_gamma_z,ZZ__,'ZZ__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_delta_v*n_gamma_z*n_omp_sub__in,ZZ_omp__
     $     ,'ZZ_omp__',flag_memory_checkset)
      call cxs_c16(n_w_max,ZZ_sub_,'ZZ_sub_',flag_memory_checkset)
      call cxs_c16(n_w_max*n_omp_sub__in,ZZ_sub_omp__
     $     ,'ZZ_sub_omp__',flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_sub__,'CTF_R_S_sub__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z*n_omp_sub__in,CTF_R_S_sub_omp__
     $     ,'CTF_R_S_sub_omp__',flag_memory_checkset)

      call cxs_i4(n_S,n_S_0_per_,'n_S_0_per_',flag_memory_checkset)
      call cxs_i4(n_S,n_S_0_sum_,'n_S_0_sum_',flag_memory_checkset)

      call cxs_c16(n_gamma_z*n_CTF*nS_0_per_max,CTF_R_S__
     $     ,'CTF_R_S__',flag_memory_checkset)

      call cxs_i4(nS_0_per_max,n_S_1_per_,'n_S_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nS_0_per_max,n_S_1_sum_,'n_S_1_sum_'
     $     ,flag_memory_checkset)

      call cxs_c16(n_C_trn0,C_trn0_,'C_trn0_',flag_memory_checkset)
      call cxs_c16(n_C_trn0_omp,C_trn0_omp__,'C_trn0_omp__'
     $     ,flag_memory_checkset)
      if (flag_tesselation.eqv..true.) then
         call cxs_c16(n_CTF_R_S_local_omp
     $        ,CTF_R_S_local_omp__,'CTF_R_S_local_omp__'
     $        ,flag_memory_checkset)
      end if 

      if (flag_tesselation.eqv..true.) then
         call cxs_r8(3*n_omp_sub__in,vp_input_omp__
     $        ,'vp_input_omp__',flag_memory_checkset)
         call cxs_l2(nS_0_per_max*n_omp_sub__in,flag_S_use_omp__
     $        ,'flag_S_use_omp__',flag_memory_checkset)
         call cxs_i4(2*nS_0_per_max*n_omp_sub__in,LT_omp__
     $        ,'LT_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_S_index_local_omp__
     $        ,'S_alpha_S_index_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_polar_a_local_omp__
     $        ,'S_alpha_polar_a_local_omp__',flag_memory_checkset)
         call cxs_r8(n_S_alpha_S_index_local_omp
     $        ,S_alpha_azimu_b_local_omp__
     $        ,'S_alpha_azimu_b_local_omp__',flag_memory_checkset)
         call cxs_i4(n_S_alpha_S_index_local_omp
     $        ,I_S_sample_local_omp__
     $        ,'I_S_sample_local_omp__',flag_memory_checkset)
      end if                    

      if ((flag_RTRT_vs_RTTR.eqv..false.)) then         
         call cxs_c16(n_r*nS_0_per_max*n_w_max,O_S_q__,'O_S_q__'
     $        ,flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_O_S_q_local_omp
     $           ,O_S_q_local_omp__,'O_S_q_local_omp__'
     $           ,flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,T_S_q__
     $        ,'T_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_T_S_q_local_omp,T_S_q_local_omp__
     $           ,'T_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         call cxs_c16(n_r*n_transf*nS_0_per_max*n_w_max,Z_S_q__
     $        ,'Z_S_q__',flag_memory_checkset)
         if (flag_tesselation.eqv..true.) then
            call cxs_c16(n_Z_S_q_local_omp,Z_S_q_local_omp__
     $           ,'Z_S_q_local_omp__',flag_memory_checkset)
         end if 
      end if 

      call cxs_i4(n_M,n_M_0_per_,'n_M_0_per_',flag_memory_checkset)
      call cxs_i4(n_M,n_M_0_sum_,'n_M_0_sum_',flag_memory_checkset)

      if (flag_tesselation.eqv..false.) then
         call cxs_r8(nM_0_per_max,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(nM_0_per_max,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      else                      
         call cxs_r8(n_M,polar_a_est_,'polar_a_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,azimu_b_est_,'azimu_b_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,gamma_z_est_,'gamma_z_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_x_est_,'delta_x_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,delta_y_est_,'delta_y_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,l2_norm_est_,'l2_norm_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,ctf_ind_est_,'ctf_ind_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,S_index_est_,'S_index_est_'
     $        ,flag_memory_checkset)
         call cxs_r8(n_M,M_index_est_,'M_index_est_'
     $        ,flag_memory_checkset)
      end if                    
      call cxs_r8(n_alpha,alpha__in_,'alpha__in_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_omp_sub__in,alpha__in_omp__
     $     ,'alpha__in_omp__',flag_memory_checkset)

      if ((flag_RTRT_vs_RTTR.eqv..true.)) then
         if (flag_tesselation.eqv..false.) then
            na = n_O_T_R_CTF_M_q__
            call cxs_c16(na,O_T_R_CTF_M_q__,'O_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_O_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,O_T_R_CTF_M_q_local_omp__ 
     $           ,'O_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if                 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_T_T_R_CTF_M_q__
            call cxs_c16(na,T_T_R_CTF_M_q__,'T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,T_T_R_CTF_M_q_local_omp__ 
     $           ,'T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_Z_T_R_CTF_M_q__
            call cxs_c16(na,Z_T_R_CTF_M_q__,'Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,Z_T_R_CTF_M_q_local_omp__ 
     $           ,'Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(nM_0_per_max,n_M_1_per_,'n_M_1_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(nM_0_per_max,n_M_1_sum_,'n_M_1_sum_'
     $     ,flag_memory_checkset)

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.2) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..true.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_S_svdd_,'Z_S_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      if ((svd_calculation_type.eq.1) .and.
     $     (flag_RTRT_vs_RTTR.eqv..false.))then
         na = n_svd_l*n_delta_v
         call cxs_c16(na,Z_M_svdd_,'Z_M_svdd_',flag_memory_checkset)
         if (flag_tesselation.eqv..false.) then
            na = n_S_Z_T_R_CTF_M_q__
            call cxs_c16(na,S_Z_T_R_CTF_M_q__,'S_Z_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_Z_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_Z_T_R_CTF_M_q_local_omp__ 
     $           ,'S_Z_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
         if (flag_tesselation.eqv..false.) then
            na = n_S_T_T_R_CTF_M_q__
            call cxs_c16(na,S_T_T_R_CTF_M_q__,'S_T_T_R_CTF_M_q__'
     $           ,flag_memory_checkset)
         else 
            na = n_S_T_T_R_CTF_M_q_local_omp__
            call cxs_c16(na,S_T_T_R_CTF_M_q_local_omp__ 
     $           ,'S_T_T_R_CTF_M_q_local_omp__',flag_memory_checkset)
         end if 
      end if                    

      call cxs_i4(n_omp_sub__in,n_S_9_per_,'n_S_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_S_9_sum_,'n_S_9_sum_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_per_,'n_M_9_per_'
     $     ,flag_memory_checkset)
      call cxs_i4(n_omp_sub__in,n_M_9_sum_,'n_M_9_sum_'
     $     ,flag_memory_checkset)

      call cxs_i4(n_omp_sub__in,n_S_use_sum_,'n_S_use_sum_'
     $     ,flag_memory_checkset)

      if ((flag_memory_checkset.eqv..true.) .and. (verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if 
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
      end if 

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_checkset_variable_0]'
      end if 

      end if 



      if (verbose.gt.1) then
         write(6,'(A)')
     $        '[entering ti8_deallocate_0]'
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nr))
         call dfftw_destroy_plan(fftw_plan_back_(nr))
      enddo 
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_in1_)
      deallocate(fftw_out_)

      if (verbose.gt.1) then
         write(6,'(A)') 'destroying fftw_plans for omp sub-blocks.'
      end if 
      do nomp_sub=0,n_omp_sub__in-1
         n_r_tmp = nomp_sub*n_r
         do nr=0,n_r-1
            call dfftw_destroy_plan(fftw_plan_frwd__(n_r_tmp+nr))
            call dfftw_destroy_plan(fftw_plan_back__(n_r_tmp+nr))
         enddo 
      enddo 
      deallocate(fftw_plan_frwd__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_in1__)
      deallocate(fftw_out__)

      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_many_plans for omp.'
      end if
      do nomp_sub=0,n_omp_sub__in-1
         call dfftw_destroy_plan_(fpm_back_(nomp_sub))
c$$$         call dfftw_destroy_plan_(fpm_frwd_(nomp_sub))
      enddo 
      deallocate(fpm_in1__)
      deallocate(fpm_out__)
      deallocate(fpm_n__)
      deallocate(fpm_inembed__)
      deallocate(fpm_onembed__)

      if (verbose.gt.1) then
         write(6,'(A)')
     $        '[finished ti8_deallocate_0]'
      end if



      if (verbose.gt.0) then

         timing_tot = timing_total_fftw_plan
         gnump_tot = gnump_total_fftw_plan
         write(timing_string,'(A)') 'fftw_plan'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 

         timing_tot = timing_total_CTF_R_S
         gnump_tot = gnump_total_CTF_R_S
         write(timing_string,'(A)') 'CTF_R_S'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_O_S_q
         gnump_tot = gnump_total_O_S_q
         write(timing_string,'(A)') 'S_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.))then
         timing_tot = timing_total_T_S_q
         gnump_tot = gnump_total_T_S_q
         write(timing_string,'(A)') 'T_S_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.)) then
         timing_tot = timing_total_Z_S_q
         gnump_tot = gnump_total_Z_S_q
         write(timing_string,'(A)') 'Z_S_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         if (flag_RTRT_vs_RTTR.eqv..true.) then
         timing_tot = timing_total_O_T_R_CTF_M_q
         gnump_tot = gnump_total_O_T_R_CTF_M_q
         write(timing_string,'(A)') 'O_T_R_CTF_M_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_T_T_R_CTF_M_q
         gnump_tot = gnump_total_T_T_R_CTF_M_q
         write(timing_string,'(A)') 'T_T_R_CTF_M_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_Z_T_R_CTF_M_q
         gnump_tot = gnump_total_Z_T_R_CTF_M_q
         write(timing_string,'(A)') 'Z_T_R_CTF_M_q'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         timing_tot = timing_total_zgemm
         gnump_tot = gnump_total_zgemm
         write(timing_string,'(A)') 'zgemm'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 

         if (flag_fill.eqv..true.) then
         timing_tot = timing_total_fpm_fill
         gnump_tot = gnump_total_fpm_fill
         write(timing_string,'(A)') 'fpm_fill'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         timing_tot = timing_total_fpm_fftw
         gnump_tot = gnump_total_fpm_fftw
         write(timing_string,'(A)') 'fpm_fftw'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_transpose
         gnump_tot = gnump_total_transpose
         write(timing_string,'(A)') 'transpose'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

         timing_tot = timing_total_Zstore
         gnump_tot = gnump_total_Zstore
         write(timing_string,'(A)') 'Zstore'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 

         if (flag_time_Zstore) then
         timing_tot = timing_total_Zstore_a
         gnump_tot = gnump_total_Zstore_a
         write(timing_string,'(A)') 'Zstore_a'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         timing_tot = timing_total_Zstore_b
         gnump_tot = gnump_total_Zstore_b
         write(timing_string,'(A)') 'Zstore_b'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         timing_tot = timing_total_Zstore_c
         gnump_tot = gnump_total_Zstore_c
         write(timing_string,'(A)') 'Zstore_c'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         timing_tot = timing_total_Zstore_x
         gnump_tot = gnump_total_Zstore_x
         write(timing_string,'(A)') 'Zstore_x'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         timing_tot = timing_total_Zstore_d
         gnump_tot = gnump_total_Zstore_d
         write(timing_string,'(A)') 'Zstore_d'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         timing_tot = timing_total_Zstore_e
         gnump_tot = gnump_total_Zstore_e
         write(timing_string,'(A)') 'Zstore_e'



      write(6,'(2(A,F8.4),A,A)') ' timing: ' , timing_tot , '; gnump: '
     $     , gnump_tot/max(1.0d0,timing_tot *1.0d9) , '; <-- ' ,
     $     trim(timing_string) 
         end if 

      end if 

 10   continue
      
      if (verbose.gt.0) then
         if (verbose.gt.1) then
         if (flag_tesselation.eqv..true.) then
            call print_all_i4(n_M_9_sub_use,n_S_use_sum_
     $           ,' n_S_use_sum_: ')
            call print_all_i4(n_M,n_SM_use_local_
     $           ,' n_SM_use_local_: ')
         end if 
         end if 
         if (flag_tesselation.eqv..true.) then
            write(6,'(A,I0,A,I0,A)') ' n_SM_use_local_ in range: [' ,
     $           min_i4_f(n_M,n_SM_use_local_) , ',' , max_i4_f(n_M
     $           ,n_SM_use_local_) , ']'
         end if 
         write(6,'(A,I0,A,I0,A,F6.2,A)')
     $        ' [finished test_innerproduct_8]: calculating ' ,
     $        n_S_use_sum , ' out of ' , n_S*n_M ,
     $        ' image-template pairs: ' , 100.0d0*n_S_use_sum / (1.0d0
     $        *n_S*n_M) , '%'
      end if 

      end
