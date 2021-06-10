      subroutine test_innerproduct_fast_wrapper_4(verbose,rseed
     $     ,n_omp_sub,n_cur,n_polar_a_,grid_p_,half_diameter_x_c,n_S
     $     ,I_S_sample_,ld_S ,S_p__,tesselation_distance_req
     $     ,S_alpha_polar_a_ ,S_alpha_azimu_b_,n_M,I_M_sample_,ld_M
     $     ,M_p__,n_CTF ,ld_CTF ,CTF_p__ ,alpha_est_
     $     ,alpha_update_f,flag_MS_vs_SM ,n_SM_max,n_SM_,alpha_SM__
     $     ,n_pixels_in ,displacement_max ,n_delta_x ,n_delta_y
     $     ,n_gamma_z ,svd_calculation_type ,eps_svd,fpm_howmany_max
     $     ,C_M_)
c$$$      Reads in a few templates and images and compares them. 
c$$$      Uses adaptive tesselation to sample from templates.
c$$$      The output (i.e., which innerproducts are highest) is returned.
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer *4 rseed !random seed. ;
      integer n_omp_sub ! number of batches to use when finding innerproducts
c$$$  ld_S is the stride between adjacent templates in S_p__
c$$$  ld_M is the stride between adjacent images in M_p__
c$$$  ld_CTF is the stride between adjacent ctfs in CTF_p__
      integer *4 n_cur,n_polar_a_(0:n_cur-1),n_S,I_S_sample_(0:n_S
     $     -1),ld_S,n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF
      integer *4 n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type
      real *8 grid_p_(0:0),half_diameter_x_c
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      integer *4 n_point,npoint
      real *8, allocatable :: L_(:)
      real *8 tradius_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
c$$$      T_ structure
      integer *4, allocatable :: T_nl_(:) !level
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_tr_(:) !tradius
      integer *4, allocatable :: T_ll_(:) !number of points from L_ in T_
      logical   , allocatable :: T_lf_(:) !is leaf
      integer *4, allocatable :: T_c0_(:) !child_0 tesselation_index
      integer *4, allocatable :: T_c1_(:) !child_1 tesselation_index
      integer *4, allocatable :: T_c2_(:) !child_2 tesselation_index
      integer *4, allocatable :: T_c3_(:) !child_3 tesselation_index
      integer *4, allocatable :: T_ls_(:) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4, allocatable :: T_LT_(:) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4, allocatable :: T_root_base_(:) !T_ roots
      real *8 tesselation_distance_req
c$$$      SM storage
      integer *4 n_SM_use ! temporary: total number of image-template comparisons made. ;
      integer *4 n_SM_max ! total (maximum) number of templates to store per image. ;
      integer *4 n_SM_(0:0) ! array of size n_M indicating the actual number of templates stored per image. ;
      real *8 alpha_SM__(0:0) ! array of size n_alpha*n_SM_max*n_M storing the image-parameters for each stored template-image pair. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      real *8 alpha_est_(0:0)
      real *8 alpha_update_f
      include 'nalpha_define.f'
      logical flag_MS_vs_SM
      real *8 eps_svd,eps_target,n_pixels_in,displacement_max
      complex *16 S_p__(0:0),M_p__(0:0),CTF_p__(0:0)
      complex *16 C_M_(0:n_M-1)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      integer n_svd_max
      parameter (n_svd_max=512)
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels
      complex *16, allocatable :: Z_S_svdd_(:)
      complex *16, allocatable :: Z_M_svdd_(:)
      complex *16, allocatable :: CTF_R_S_(:)
      complex *16, allocatable :: S_q__(:)
      complex *16, allocatable :: T_S_q__(:)
      complex *16, allocatable :: Z_S_q__(:)
      complex *16, allocatable :: S_p_(:)
      complex *16, allocatable :: S_q_(:)
      complex *16, allocatable :: S_p_omp__(:)
      complex *16, allocatable :: S_q_omp__(:)
      complex *16, allocatable :: M_p_pre_(:)
      complex *16, allocatable :: M_p_(:)
      complex *16, allocatable :: M_q_(:)
      complex *16, allocatable :: M_p_pre_omp__(:)
      complex *16, allocatable :: M_p_omp__(:)
      complex *16, allocatable :: M_q_omp__(:)
      complex *16, allocatable :: CTF_p_(:)
      complex *16, allocatable :: CTF_q_(:)
      complex *16, allocatable :: CTF_p_omp__(:)
      complex *16, allocatable :: CTF_q_omp__(:)
      complex *16, allocatable :: T_R_CTF_M_q__(:)
      complex *16, allocatable :: T_T_R_CTF_M_q__(:)
      complex *16, allocatable :: Z_T_R_CTF_M_q__(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,ns,nctf,nm,nw
      integer *4, allocatable :: n_w_(:)
      real *8 pi
      real *8 max_r8_f
      real *8 max_c
c$$$      estimated (current) displacement and rotation parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
c$$$      estimated l2_norm and image-dependent ctf_index
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: ctf_ind_est_(:)
      real *8, allocatable :: S_index_est_(:)
c$$$      updated (future) displacement and rotation parameters for each image
      real *8, allocatable :: polar_a_upd_(:)
      real *8, allocatable :: azimu_b_upd_(:)
      real *8, allocatable :: gamma_z_upd_(:)
      real *8, allocatable :: delta_x_upd_(:)
      real *8, allocatable :: delta_y_upd_(:)
c$$$      updated l2_norm and image-dependent ctf_index
      real *8, allocatable :: l2_norm_upd_(:)
      real *8, allocatable :: ctf_ind_upd_(:)
      real *8, allocatable :: S_index_upd_(:)
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16 C_M,C_Z
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
      integer fpm_howmany_max
      integer n_transf,fpm_howmany
      integer n_fpm,nfpm
      integer *8, allocatable :: fpm_frwd_(:)
      integer *8, allocatable ::  fpm_back_(:)
      complex *16, allocatable :: fpm_in1__(:)
      complex *16, allocatable :: fpm_out__(:)
      integer fpm_rank
      integer, allocatable :: fpm_n__(:)
      integer, allocatable :: fpm_inembed__(:)
      integer fpm_istride
      integer fpm_idist
      integer, allocatable :: fpm_onembed__(:)
      integer fpm_ostride
      integer fpm_odist
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
c$$$      parameters for timing and testing
      logical flag_fill
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
c$$$      parameters for memory map
      real *8 d_mem
      integer verbose_mem
      parameter (verbose_mem=1)
c$$$      format_string for printing output
      character(len=64) format_string

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_fast_wrapper_4]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_cur: ',n_cur
         write(format_string,'(A,I0,A)') '(A,',n_cur,'(I0,1X))'
         write(6,format_string) 'n_polar_a_: ',(n_polar_a_(nr),nr=0
     $        ,n_cur-1)
         write(format_string,'(A,I0,A)') '(A,',n_cur,'(F8.3,1X))'
         write(6,format_string) 'grid_p_: ',(grid_p_(nr),nr=0
     $        ,n_cur-1)
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
         write(6,'(A,I0)') 'n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
         write(6,'(A)') 'alpha_est_: '
         write(6,format_string) (alpha_est_(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') 'n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') 'displacement_max: ',displacement_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
         write(6,'(A,F6.3)') 'eps_svd: ',eps_svd
      end if !if (verbose.gt.1) then

      pi = 4*atan(1.0)
      eps_target = eps_svd

c$$$      Calculating template size using 'get_template_size'
      allocate(ngridc_(n_cur))
      allocate(icstart_(n_cur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'n_cur = ',n_cur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(n_polar_a_,n_cur,ntemplatesize,ngridc_
     $     ,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
         write(format_string,'(A,I0,A)') '(A,',n_cur,'(I0,1X))'
         write(6,format_string) 'ngridc_: ',(ngridc_(nr),nr=1,n_cur)
      end if
      
c$$$      indices
      n_r = n_cur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r',n_r,'<2'
      end if
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_w_max ',n_w_max,'; n_A ',n_A
         write(format_string,'(A,I0,A)') '(A,',n_cur,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,n_cur-1)
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

      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_x*n_delta_y
      else
         n_transf = n_svd_l
      end if !if (svd_calculation_type.eq.2) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up displacement-operator Z_S_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to S.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      allocate(Z_S_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      timing_tic = omp_get_wtime()
      call innerproduct_q_k_svdd_1(.true.,n_svd_d,svd_d_,n_svd_l,svd_l_
     $     ,svd_U_d_,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_S_svdd_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q_k_svdd_1: total time '
     $        ,timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y
     $        *1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' innerproduct_q_k_svdd_1 total Gnumps: '
     $        ,timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up displacement-operator Z_M_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to M.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      allocate(Z_M_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      timing_tic = omp_get_wtime()
      call innerproduct_q_k_svdd_1(.false.,n_svd_d,svd_d_,n_svd_l
     $     ,svd_l_,svd_U_d_,n_delta_x,delta_x_,n_delta_y,delta_y_
     $     ,Z_M_svdd_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q_k_svdd_1: total time '
     $        ,timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y
     $        *1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' innerproduct_q_k_svdd_1 total Gnumps: '
     $        ,timing_tmp
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
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_gamma_z*1.0d0*n_S*1.0d0*n_CTF*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' CTF_R_S_ requires ' , d_mem*1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
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
         write(6,'(A)') ' of the pointwise product of CTF and R(S),'
         write(6,'(A)') ' where CTF is the nctf-th ctf-function, '
         write(6,'(A)') ' R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' and S is the ns-th template, '
         write(6,'(A)') ' at position S_p__(I_S_sample_(ns)*ld_S). '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_CTF_R_S_1(verbose-1,n_gamma_z
     $     ,gamma_z_,fftw_plan_frwd_,fftw_in1_,fftw_out_
     $     ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_,n_r
     $     ,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,I_S_sample_,ld_S
     $     ,S_p__,CTF_p_ ,CTF_q_ ,n_CTF,ld_CTF ,CTF_p__
     $     ,CTF_R_S_)
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
         write(6,'(A)') ' Allocating array S_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each template'
         write(6,'(A)') ' of the form: (S)_q, '
         write(6,'(A)') ' '
      end if                    !if (verbose.gt.1) then
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_S*1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_q__ requires ' , d_mem *1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(S_q__(0:n_r*n_S*n_w_max))
      call cl1_c16(n_r*n_S*n_w_max,S_q__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate bessel coefficients of S.'
         write(6,'(A)') ' S = template in k-space polar coord.'
         write(6,'(A)') ' '
         write(6,'(A,A)') ' S_q__(nr + n_r*(ns + n_S*nw)) '
         write(6,'(A)') ' is equal to the bessel-coefficients'
         write(6,'(A)') ' conjg(S_q_(ic)) '
         write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
         write(6,'(A)') ' where S_q_ = S, with '
         write(6,'(A)') ' S <-- S_p__(I_S_sample_(ns)*ld_S).'
         write(6,'(A)') ' '
      end if                    !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c     $OMP PARALLEL PRIVATE(nS_per,IS_per,n_r_tmp,n_A_tmp,n_X_tmp)
c     $OMP DO 
      do nS_sub=0,n_S_sub_use-1
         nS_per = n_S_per_(nS_sub)
         IS_per = I_S_per_(nS_sub)
         n_r_tmp = nS_sub*n_r
         n_A_tmp = nS_sub*n_A
         n_X_tmp = 0 + n_r*(0 + IS_per + n_S*0)
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' test_innerproduct_fast_S_q_0'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call test_innerproduct_fast_S_q_0(verbose-1
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_p_,n_w_,n_A
     $        ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_per,I_S_sample_
     $        ,IS_per ,ld_S ,S_p__,n_S,S_q__(n_X_tmp))
      enddo                     !do nS_sub=0,n_S_sub_use-1
c     $OMP END DO
c     $OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating S_q__ for each S. '
     $        ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' S_q__ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if                    ! if (verbose.gt.0) then

      if (svd_calculation_type.eq.2) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array T_S_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each template'
            write(6,'(A)') ' of the form: (T_{delta}(S))_q, '
            write(6,'(A)') ' where T represents translation by delta.'         
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*n_S*1.0d0
     $        *n_w_max*16.0d0
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' T_S_q__ requires ' , d_mem
     $           *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 !if (verbose.gt.0) then
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
            write(6,'(A)') ' and S <-- S_p__(I_S_sample_(ns)*ld_S).'
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
     $              ,' test_innerproduct_fast_T_S_q_2'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call test_innerproduct_fast_T_S_q_2(verbose-1,n_delta_x
     $           ,delta_x_,n_delta_y,delta_y_,fftw_plan_frwd__(n_r_tmp)
     $           ,fftw_in1__(n_A_tmp) ,fftw_out__(n_A_tmp) ,n_r ,grid_p_
     $           ,n_w_,n_A,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp),nS_per
     $           ,I_S_sample_,IS_per ,ld_S ,S_p__,n_S
     $           ,T_S_q__(n_X_tmp))
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
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*n_transf*1.0d0*n_S*1.0d0
     $        *n_w_max*16.0d0
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' Z_S_q__ requires ' , d_mem
     $           *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 !if (verbose.gt.0) then
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
            write(6,'(A)') ' and S <-- S_p__(I_S_sample_(ns)*ld_S).'
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
     $              ,' test_innerproduct_fast_T_S_q_2'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call test_innerproduct_fast_Z_S_q_2(verbose-1,n_svd_r,svd_r_
     $           ,n_svd_l,svd_l_,svd_s_,svd_V_r_
     $           ,fftw_plan_frwd__(n_r_tmp) ,fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp),n_r ,grid_p_,n_w_,n_A
     $           ,S_p_omp__(n_A_tmp),S_q_omp__(n_A_tmp) ,nS_per
     $           ,I_S_sample_,IS_per,ld_S ,S_p__,n_S
     $           ,Z_S_q__(n_X_tmp))
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
      
c$$$      reading in estimated translations and rotations
      allocate(polar_a_est_(0:n_M-1))
      allocate(azimu_b_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      allocate(ctf_ind_est_(0:n_M-1))
      allocate(S_index_est_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_est_(nm) = alpha_est_(nalpha_polar_a + nm*n_alpha)
         azimu_b_est_(nm) = alpha_est_(nalpha_azimu_b + nm*n_alpha)
         gamma_z_est_(nm) = alpha_est_(nalpha_gamma_z + nm*n_alpha)
         delta_x_est_(nm) = alpha_est_(nalpha_delta_x + nm*n_alpha)
         delta_y_est_(nm) = alpha_est_(nalpha_delta_y + nm*n_alpha)
         l2_norm_est_(nm) = alpha_est_(nalpha_l2_norm + nm*n_alpha)
         ctf_ind_est_(nm) = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         S_index_est_(nm) = alpha_est_(nalpha_S_index + nm*n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->',polar_a_est_(nm)
     $           ,azimu_b_est_(nm),gamma_z_est_(nm),delta_x_est_(nm)
     $           ,delta_y_est_(nm),l2_norm_est_(nm),ctf_ind_est_(nm)
     $           ,S_index_est_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array T_R_CTF_M_q__ to hold'
         write(6,'(A)') ' bessel-coefficients for each image'
         write(6,'(A)') ' T_{-delta_est}(R_{-gamma_est}(CTF.*M)), '
         write(6,'(A)') ' where T = translation by -delta_est.'
         write(6,'(A)') ' and R = rotation by -gamma_est.'
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_r*1.0d0*n_M*1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' T_R_CTF_M_q_ requires ' , d_mem
     $        *1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(T_R_CTF_M_q__(0:n_r*n_M*n_w_max))
      call cl1_c16(n_r*n_M*n_w_max,T_R_CTF_M_q__)
      allocate(M_p_pre_(0:n_A-1))
      allocate(M_p_(0:n_A-1))
      allocate(M_q_(0:n_A-1))
      allocate(M_p_pre_omp__(0:n_A*n_omp_sub-1))
      allocate(M_p_omp__(0:n_A*n_omp_sub-1))
      allocate(M_q_omp__(0:n_A*n_omp_sub-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' Calculate bessel coefficients of T(R(CTF.*M)).'
         write(6,'(A)') ' T = translation by -delta_est. '
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
         write(6,'(A)') ' and M <-- M_p__(I_M_sample_(nm)*ld_M).'
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
         call test_innerproduct_fast_T_R_CTF_M_q_2(verbose-1
     $        ,delta_x_est_(IM_per),delta_y_est_(IM_per)
     $        ,gamma_z_est_(IM_per),ctf_ind_est_(IM_per)
     $        ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $        ,fftw_out__(n_A_tmp) ,n_r ,grid_p_,n_w_,n_A
     $        ,M_p_omp__(n_A_tmp),M_q_omp__(n_A_tmp),nM_per
     $        ,I_M_sample_ ,IM_per,ld_M ,M_p__,CTF_p_omp__(n_A_tmp)
     $        ,n_CTF ,ld_CTF ,CTF_p__ ,C_M_(IM_per) ,n_M
     $        ,T_R_CTF_M_q__(n_X_tmp))
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
            write(6,'(A)') ' Allocating array T_T_R_CTF_M_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each image'
            write(6,'(A)') ' T(T_est(R_{-gamma_est}(CTF.*M))), '
            write(6,'(A)') ' where T = translation by -delta.'
            write(6,'(A)') ' and T_est = translation by -delta_est.'
            write(6,'(A)') ' and R = rotation by -gamma_est.'
            write(6,'(A)') ' '
         end if
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*n_delta_x*1.0d0*n_delta_y*1.0d0
     $        *n_M*1.0d0*n_w_max*16.0d0
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' T_T_R_CTF_M_q_ requires ' ,
     $           d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 !if (verbose.gt.0) then
         allocate(T_T_R_CTF_M_q__(0:n_r*n_delta_x*n_delta_y*n_M
     $        *n_w_max))
         call cl1_c16(n_r*n_delta_x*n_delta_y*n_M*n_w_max
     $        ,T_T_R_CTF_M_q__)
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Calc bessel coefficients of T(T_est(R(CTF.*M))).'
            write(6,'(A)') ' T = translation by -delta. '
            write(6,'(A)') ' T_est = translation by -delta_est. '
            write(6,'(A)') ' R = rotation by -gamma_est.'
            write(6,'(A)') ' CTF = conjg(CTF) in k-space polar coord.'
            write(6,'(A)') ' M = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A,A)') ' T_T_R_CTF_M_q__(nr + n_r*(ndx + ' ,
     $           'n_delta_x*(ndy + n_delta_y*(nm ' , '+ n_M*nw)))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' M_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
            write(6,'(A)') ' where M_q_ = T(T_est(R(CTF.*M))), with '
            write(6,'(A)') ' T <-- -delta'
            write(6,'(A)') ' T_est <-- -delta_est'
            write(6,'(A)') ' R <-- -gamma_est'
            write(6,'(A)') ' and CTF <-- conjg(CTF_(nctf))'
            write(6,'(A)') ' and M <-- M_p__(I_M_sample_(nm)*ld_M).'
            write(6,'(A)') ' '
         end if
         timing_tic = omp_get_wtime()
c     $OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c     $OMP DO 
         do nM_sub=0,n_M_sub_use-1
            nM_per = n_M_per_(nM_sub)
            IM_per = I_M_per_(nM_sub)
            n_r_tmp = nM_sub*n_r
            n_A_tmp = nM_sub*n_A
            n_X_tmp = 0 + n_r*(0 + n_delta_x*(0 + n_delta_y*(IM_per +
     $           n_M*0)))
            call test_innerproduct_fast_T_T_R_CTF_M_q_0(verbose-1
     $           ,n_delta_x,delta_x_,n_delta_y,delta_y_
     $           ,delta_x_est_(IM_per),delta_y_est_(IM_per)
     $           ,gamma_z_est_(IM_per),ctf_ind_est_(IM_per)
     $           ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp) ,n_r ,grid_p_,n_w_,n_A
     $           ,M_p_pre_omp__(n_A_tmp),M_p_omp__(n_A_tmp)
     $           ,M_q_omp__(n_A_tmp),nM_per ,I_M_sample_ ,IM_per,ld_M
     $           ,M_p__,CTF_p_omp__(n_A_tmp) ,n_CTF ,ld_CTF ,CTF_p__
     $           ,C_M_(IM_per) ,n_M ,T_T_R_CTF_M_q__(n_X_tmp))
         enddo                  !do nM_sub=0,n_M_sub_use-1
c     $OMP END DO
c     $OMP END PARALLEL
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         if (verbose.gt.0) then
            write(6,'(A,F8.5)')
     $           ' finished calculating T_T_R_CTF_M_q__ for each M. '
     $           ,timing_tot
            timing_tmp = (n_A*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
            write(6,'(A,F8.4)') ' T_T_R_CTF_M_q__ total Gnumps: ' ,
     $           timing_tmp
            write(6,'(A)') ' '
         end if                 ! if (verbose.gt.0) then
      end if !if (svd_calculation_type.eq.2) then

      if (svd_calculation_type.ne.2) then
         if (verbose.gt.1) then
            write(6,'(A)') ' Allocating array Z_T_R_CTF_M_q__ to hold'
            write(6,'(A)') ' bessel-coefficients for each image'
            write(6,'(A)') ' Z(T_est(R_{-gamma_est}(CTF.*M))), '
            write(6,'(A)') ' where Z represents the left hand side'
            write(6,'(A)') ' of the svd of the translation operator.'
            write(6,'(A)') ' and T_est = translation by -delta_est.'
            write(6,'(A)') ' and R = rotation by -gamma_est.'
            write(6,'(A)') ' '
         end if
         d_mem = 0.0d0
         d_mem = d_mem + 1.0d0*n_r*1.0d0*n_svd_l*1.0d0
     $        *n_M*1.0d0*n_w_max*16.0d0
         if (verbose.gt.0 .or. verbose_mem.gt.0) then
            write (6,'(A,2(F16.8,A))') ' Z_T_R_CTF_M_q_ requires ' ,
     $           d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
         end if                 !if (verbose.gt.0) then
         allocate(Z_T_R_CTF_M_q__(0:n_r*n_svd_l*n_M
     $        *n_w_max))
         call cl1_c16(n_r*n_svd_l*n_M*n_w_max ,Z_T_R_CTF_M_q__)
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' Calc bessel coefficients of Z(T_est(R(CTF.*M))).'
            write(6,'(A)') ' Z = left hand side of svd-operator.'
            write(6,'(A)') ' T_est = translation by -delta_est.'
            write(6,'(A)') ' R = rotation by -gamma_est.'
            write(6,'(A)') ' CTF = conjg(CTF) in k-space polar coord.'
            write(6,'(A)') ' M = template in k-space polar coord.'
            write(6,'(A)') ' '
            write(6,'(A,A,A)') ' Z_T_R_CTF_M_q__(nr + n_r*(nl + ' ,
     $           'n_svd_l*(nm ' , '+ n_M*nw))) '
            write(6,'(A)') ' is equal to the bessel-coefficients'
            write(6,'(A)') ' M_q_(ic)) '
            write(6,'(A)') ' for ic = nw + n_w_sum_(nr),'
            write(6,'(A)') ' where M_q_ = Z(T_est(R(CTF.*M))), with '
            write(6,'(A)') ' Z <-- svd-operator'
            write(6,'(A)') ' T_est <-- -delta_est'
            write(6,'(A)') ' R <-- -gamma_est'
            write(6,'(A)') ' and CTF <-- conjg(CTF_(nctf))'
            write(6,'(A)') ' and M <-- M_p__(I_M_sample_(nm)*ld_M).'
            write(6,'(A)') ' '
         end if
         timing_tic = omp_get_wtime()
c     $OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c     $OMP DO 
         do nM_sub=0,n_M_sub_use-1
            nM_per = n_M_per_(nM_sub)
            IM_per = I_M_per_(nM_sub)
            n_r_tmp = nM_sub*n_r
            n_A_tmp = nM_sub*n_A
            n_X_tmp = 0 + n_r*(0 + n_svd_l*(IM_per + n_M*0))
            call test_innerproduct_fast_Z_T_R_CTF_M_q_0(verbose-1
     $           ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_
     $           ,delta_x_est_(IM_per),delta_y_est_(IM_per)
     $           ,gamma_z_est_(IM_per),ctf_ind_est_(IM_per)
     $           ,fftw_plan_frwd__(n_r_tmp),fftw_in1__(n_A_tmp)
     $           ,fftw_out__(n_A_tmp) ,n_r ,grid_p_,n_w_,n_A
     $           ,M_p_omp__(n_A_tmp) ,M_q_omp__(n_A_tmp),nM_per
     $           ,I_M_sample_ ,IM_per,ld_M ,M_p__,CTF_p_omp__(n_A_tmp)
     $           ,n_CTF ,ld_CTF ,CTF_p__ ,C_M_(IM_per) ,n_M
     $           ,Z_T_R_CTF_M_q__(n_X_tmp))
         enddo                  !do nM_sub=0,n_M_sub_use-1
c     $OMP END DO
c     $OMP END PARALLEL
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         if (verbose.gt.0) then
            write(6,'(A,F8.5)')
     $           ' finished calculating Z_T_R_CTF_M_q__ for each M. '
     $           ,timing_tot
            timing_tmp = (n_A*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
            write(6,'(A,F8.4)') ' Z_T_R_CTF_M_q__ total Gnumps: ' ,
     $           timing_tmp
            write(6,'(A)') ' '
         end if                 ! if (verbose.gt.0) then
      end if !if (svd_calculation_type.ne.2) then

c$$$      preparing updated translations and rotations
      allocate(polar_a_upd_(0:n_M-1))
      allocate(azimu_b_upd_(0:n_M-1))
      allocate(gamma_z_upd_(0:n_M-1))
      allocate(delta_x_upd_(0:n_M-1))
      allocate(delta_y_upd_(0:n_M-1))
      allocate(l2_norm_upd_(0:n_M-1))
      allocate(ctf_ind_upd_(0:n_M-1))
      allocate(S_index_upd_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_upd_(nm) = alpha_est_(nalpha_polar_a + nm*n_alpha)
         azimu_b_upd_(nm) = alpha_est_(nalpha_azimu_b + nm*n_alpha)
         gamma_z_upd_(nm) = alpha_est_(nalpha_gamma_z + nm*n_alpha)
         delta_x_upd_(nm) = alpha_est_(nalpha_delta_x + nm*n_alpha)
         delta_y_upd_(nm) = alpha_est_(nalpha_delta_y + nm*n_alpha)
         l2_norm_upd_(nm) = alpha_est_(nalpha_l2_norm + nm*n_alpha)
         ctf_ind_upd_(nm) = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         S_index_upd_(nm) = alpha_est_(nalpha_S_index + nm*n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->',polar_a_upd_(nm)
     $           ,azimu_b_upd_(nm),gamma_z_upd_(nm),delta_x_upd_(nm)
     $           ,delta_y_upd_(nm),l2_norm_upd_(nm),ctf_ind_upd_(nm)
     $           ,S_index_upd_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1

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
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*fpm_rank*1.0d0*n_omp_sub*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_n__ requires ' , d_mem*1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
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
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_fpm*1.0d0*n_omp_sub*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_in1__ requires ' , d_mem*1.0d
     $        -6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
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

      if (tesselation_distance_req.ge.2.0d0) then
         if (verbose.gt.1) then
            write(6,'(A)') ''
            write(6,*) 'Now we call test_innerproduct_fast_omp_stage_4'
     $           ,' to calculate all innerproducts:'
         end if                 !if (verbose.gt.1) then
         timing_tic = omp_get_wtime()
         call test_innerproduct_fast_stage_4(verbose,rseed,n_omp_sub
     $        ,svd_calculation_type,eps_svd,n_svd_d,svd_d_,n_svd_r
     $        ,svd_r_ ,n_svd_l ,svd_l_,svd_U_d_,svd_s_,svd_V_r_
     $        ,Z_S_svdd_ ,Z_M_svdd_,n_r,grid_p_ ,n_w_ ,half_diameter_x_c
     $        ,n_S ,I_S_sample_,ld_S,S_p__ ,S_alpha_polar_a_
     $        ,S_alpha_azimu_b_ ,n_M,I_M_sample_,ld_M ,M_p__
     $        ,n_M_sub_use,n_M_per_,I_M_per_,n_CTF ,ld_CTF ,CTF_p__
     $        ,CTF_R_S_ ,S_q__,T_S_q__ ,Z_S_q__ ,T_R_CTF_M_q__
     $        ,T_T_R_CTF_M_q__ ,Z_T_R_CTF_M_q__ ,S_p_,S_q_ ,S_p_omp__
     $        ,S_q_omp__,M_p_ ,M_q_ ,M_p_omp__ ,M_q_omp__ ,CTF_p_,CTF_q_
     $        ,CTF_p_omp__ ,CTF_q_omp__ ,polar_a_est_ ,azimu_b_est_
     $        ,gamma_z_est_ ,delta_x_est_ ,delta_y_est_ ,l2_norm_est_
     $        ,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_
     $        ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $        ,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f ,flag_MS_vs_SM
     $        , n_SM_max , n_SM_ , alpha_SM__ ,N_pixels_in
     $        ,displacement_max ,n_delta_x ,delta_x_ ,n_delta_y
     $        ,delta_y_ ,n_gamma_z ,gamma_z_, fftw_plan_frwd_
     $        ,fftw_plan_back_ ,fftw_in1_,fftw_out_,fpm_howmany,n_fpm
     $        ,fpm_back_ ,fpm_in1__,fpm_out__, C_M_ , n_SM_use)
         timing_toc = omp_get_wtime()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3,A,F8.6,A,F8.7)')
     $           'test_innerproduct_fast_omp_stage_4:',' total_time: '
     $           ,timing_toc-timing_tic,'; time_per S,M: ',(timing_toc
     $           -timing_tic)/max(1,n_SM_use),'; time_per S,M,dx,dy: '
     $           ,(timing_toc -timing_tic)/(n_SM_use*n_delta_x
     $           *n_delta_y)
         end if
         if (verbose.gt.0) then
            write(6,'(A)')
     $           '[finished test_innerproduct_fast_wrapper_4]'
         end if
      end if                    !if (tesselation_distance_req.ge.2.0d0) then

      if (tesselation_distance_req.lt.2.0d0) then
         if (verbose.gt.1) then
            write(6,'(A)') ''
            write(6,*) 'Now we create tesselation'
     $           ,' to access templates efficiently:'
         end if
         timing_tic = omp_get_wtime()
         d_mem = 0.0d0
         n_point = n_S;
         allocate(L_(0:3*n_point-1))
         d_mem = d_mem + 3.0d0*n_point*8.0d0
         do npoint=0,n_point-1
            call get_angle_to_vp_(S_alpha_polar_a_(npoint)
     $           ,S_alpha_azimu_b_(npoint),L_(0+3*npoint))
            call normalize_r8(3,L_(0 + 3*npoint))
         enddo                  !do npoint=0,n_point-1
         if (verbose.gt.1) then
            do npoint=0,n_point-1
               write(6,'(A,F8.4,F8.4,A,A,I0,1X,3F8.4)') 'S_alpha: ' ,
     $              S_alpha_polar_a_(npoint) , S_alpha_azimu_b_(npoint)
     $              , ' --> ' , ' L_: ', npoint , L_(0 + 3*npoint), L_(1
     $              + 3 *npoint) , L_(2 + 3*npoint)
            enddo               ! do npoint=0,n_point-1
         end if                 ! if (verbose.gt.1) then
         call tesselation_get_nl_nm_ll(n_point,L_,tradius_min,nl_max
     $        ,nm_sum,ll_sum)
         d_mem = d_mem + 4.0d0*nm_sum*8.0d0
         d_mem = d_mem + 8.0d0*nm_sum*4.0d0
         d_mem = d_mem + 1.0d0*ll_sum*4.0d0
         if (verbose.gt.0) then
            write (6,'(A,F8.4,A)') ' tesselation requires ' , d_mem*1.0d
     $           -6,' MB'
         end if                 !if (verbose.gt.0) then
         allocate(T_nl_(0:1*nm_sum-1)) !level
         allocate(T_vm_(0:3*nm_sum-1)) !vertex center
         allocate(T_tr_(0:1*nm_sum-1)) !tradius
         allocate(T_ll_(0:1*nm_sum-1)) !number of points from L_ in T_
         allocate(T_lf_(0:1*nm_sum-1)) !is leaf
         allocate(T_c0_(0:1*nm_sum-1)) !child_0 tesselation_index
         allocate(T_c1_(0:1*nm_sum-1)) !child_1 tesselation_index
         allocate(T_c2_(0:1*nm_sum-1)) !child_2 tesselation_index
         allocate(T_c3_(0:1*nm_sum-1)) !child_3 tesselation_index
         allocate(T_ls_(0:1*nm_sum-1)) !starting index of point_index_list for T_ if leaf (leaves only)
         allocate(T_LT_(0:1*ll_sum-1)) !full point_index_list for all of T_ (leaves only)
         allocate(T_root_base_(0:8-1))
         call tesselation_index_wrapper_0(n_point,L_,tradius_min ,nl_max
     $        ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $        ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)
         timing_toc = omp_get_wtime()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'tesselation_index_wrapper_0:'
     $           ,' total_time: ' ,timing_toc-timing_tic
         end if
      end if !if (tesselation_distance_req.lt.2.0d0) then

      do nm=0,n_M-1
         alpha_est_(nalpha_polar_a + nm*n_alpha) = polar_a_upd_(nm)
         alpha_est_(nalpha_azimu_b + nm*n_alpha) = azimu_b_upd_(nm)
         alpha_est_(nalpha_gamma_z + nm*n_alpha) = gamma_z_upd_(nm)
         alpha_est_(nalpha_delta_x + nm*n_alpha) = delta_x_upd_(nm)
         alpha_est_(nalpha_delta_y + nm*n_alpha) = delta_y_upd_(nm)
         alpha_est_(nalpha_l2_norm + nm*n_alpha) = l2_norm_upd_(nm)
         alpha_est_(nalpha_ctf_ind + nm*n_alpha) = ctf_ind_upd_(nm)
         alpha_est_(nalpha_S_index + nm*n_alpha) = S_index_upd_(nm)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->'
     $           ,polar_a_upd_(nm) ,azimu_b_upd_(nm),gamma_z_upd_(nm)
     $           ,delta_x_upd_(nm) ,delta_y_upd_(nm),l2_norm_upd_(nm)
     $           ,ctf_ind_upd_(nm) ,S_index_upd_(nm)
         end if                 !if (verbose.gt.1) then
      enddo                     !do nm=0,n_M-1

c$$$      deallocate(ngridc_)
c$$$      deallocate(icstart_)
c$$$      deallocate(n_w_)

c$$$      deallocate(L_)
c$$$      deallocate(T_nl_) !level
c$$$      deallocate(T_vm_) !vertex center
c$$$      deallocate(T_tr_) !tradius
c$$$      deallocate(T_ll_) !number of points from L_ in T_
c$$$      deallocate(T_lf_) !is leaf
c$$$      deallocate(T_c0_) !child_0 tesselation_index
c$$$      deallocate(T_c1_) !child_1 tesselation_index
c$$$      deallocate(T_c2_) !child_2 tesselation_index
c$$$      deallocate(T_c3_) !child_3 tesselation_index
c$$$      deallocate(T_ls_) !starting index of point_index_list for T_ if leaf (leaves only)
c$$$      deallocate(T_LT_) !full point_index_list for all of T_ (leaves only)
c$$$      deallocate(T_root_base_)

c$$$      deallocate(polar_a_est_)
c$$$      deallocate(azimu_b_est_)
c$$$      deallocate(gamma_z_est_)
c$$$      deallocate(delta_x_est_)
c$$$      deallocate(delta_y_est_)
c$$$      deallocate(l2_norm_est_)
c$$$      deallocate(ctf_ind_est_)
c$$$      deallocate(S_index_est_)
c$$$      deallocate(polar_a_upd_)
c$$$      deallocate(azimu_b_upd_)
c$$$      deallocate(gamma_z_upd_)
c$$$      deallocate(delta_x_upd_)
c$$$      deallocate(delta_y_upd_)
c$$$      deallocate(l2_norm_upd_)
c$$$      deallocate(ctf_ind_upd_)
c$$$      deallocate(S_index_upd_)
c$$$      deallocate(delta_x_)
c$$$      deallocate(delta_y_)
c$$$      deallocate(gamma_z_)


      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nr))
         call dfftw_destroy_plan(fftw_plan_back_(nr))
      enddo !do nr=0,n_r-1
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_in1_)
      deallocate(fftw_out_)

      if (verbose.gt.1) then
         write(6,'(A)') 'destroying fftw_plans for omp sub-blocks.'
      end if ! if (verbose.gt.1) then
      do nomp_sub=0,n_omp_sub-1
         n_r_tmp = nomp_sub*n_r
         do nr=0,n_r-1
            call dfftw_destroy_plan(fftw_plan_frwd__(n_r_tmp+nr))
            call dfftw_destroy_plan(fftw_plan_back__(n_r_tmp+nr))
         enddo ! do nr=0,n_r-1
      enddo ! do nomp_sub=0,n_omp_sub-1
      deallocate(fftw_plan_frwd__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_in1__)
      deallocate(fftw_out__)

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

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_fast_wrapper_4]'
      end if
      
      end
