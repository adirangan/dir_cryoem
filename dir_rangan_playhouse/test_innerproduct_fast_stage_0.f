c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_fast_stage_0(verbose,n_omp_sub
     $     ,svd_calculation_type,eps_svd,n_r,grid_p_,n_w_
     $     ,half_diameter_x_c,n_S ,ld_S,S_p__ ,tesselation_distance_req
     $     ,S_alpha_polar_a_ ,S_alpha_azimu_b_,L_ ,nl_max ,nm_sum,ll_sum
     $     ,T_nl_,T_vm_ ,T_tr_,T_ll_,T_lf_ ,T_c0_ ,T_c1_ ,T_c2_,T_c3_
     $     ,T_ls_,T_LT_ ,n_T_root_base,T_root_base_ ,n_M ,ld_M,M_p__
     $     ,n_CTF ,ld_CTF ,CTF_p__,polar_a_est_ ,azimu_b_est_
     $     ,gamma_z_est_ ,delta_x_est_,delta_y_est_ ,l2_norm_est_
     $     ,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_
     $     ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $     ,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f ,N_pixels_in
     $     ,displacement_max ,n_delta_x ,n_delta_y ,n_gamma_z ,
     $     fpm_howmany_max, C_M_ ,n_SM_use)
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
      complex *16, allocatable :: M_p_(:)
      complex *16, allocatable :: M_q_(:)
      complex *16, allocatable :: T_R_CTF_M_q__(:)
      complex *16, allocatable :: S_T_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: S_Z_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: ZZ_(:)
      complex *16, allocatable :: Z_tmp_(:)
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
c$$$      fftw plan many (fpm)
      integer *8 fpm_frwd
      integer *8 fpm_back
      complex *16, allocatable :: fpm_in1_(:)
      complex *16, allocatable :: fpm_out_(:)
      integer fpm_rank
      integer, allocatable :: fpm_n_(:)
      integer fpm_howmany_max
      integer n_transf,fpm_howmany
      integer nfpm
      integer, allocatable :: fpm_inembed_(:)
      integer fpm_istride
      integer fpm_idist
      integer, allocatable :: fpm_onembed_(:)
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
c$$$      pi
      real *8 pi
c$$$      parameters for timing and testing
      logical flag_test
      real *8 timing_tic,timing_toc,timing_tot,timing_tmp
      integer *4 n_omp_sub_use,nomp_sub,n_r_tmp,n_A_tmp,nM_per,IM_per
      integer *4, allocatable :: n_M_per_(:)
      integer *4, allocatable :: I_M_per_(:)
      integer *4 nM_per_max

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_stage_0]'
      end if      

      eps_target = eps_svd

      if (verbose.gt.1) then
         write(6,'(A,I0,A,F5.3,A,F5.3,A,F3.1,A,I0,A,I0,A,I0)') 'n_r: '
     $        ,n_r,'; eps_svd: ',eps_svd,'; eps_target: ',eps_target
     $        ,'; N_pixels_in: ',N_pixels_in ,'; n_delta_x: ',n_delta_x
     $        ,'; n_delta_y: ',n_delta_y ,'; n_gamma_z: ',n_gamma_z
         write(6,'(A)') ' '
      end if

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
      allocate(CTF_p_(0:n_A-1))
      allocate(CTF_q_(0:n_A-1))
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
      call test_innerproduct_fast_CTF_R_S_0(verbose,n_gamma_z ,gamma_z_
     $     ,fftw_plan_frwd_,fftw_in1_,fftw_out_ ,fftw_plan_back_last
     $     ,fftw_in1_last_,fftw_out_last_,n_r ,grid_p_,n_w_,n_A,S_p_
     $     ,S_q_,n_S ,ld_S,S_p__,CTF_p_,CTF_q_,n_CTF,ld_CTF ,CTF_p__
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

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_many_plan for local use'
         write(6,'(A)') ' '
      end if
      fpm_rank = 1
      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_x*n_delta_y
      else
         n_transf = n_svd_l
      end if !if (svd_calculation_type.eq.2) then
      fpm_howmany = min(fpm_howmany_max,n_transf*n_S*n_M)
      allocate(fpm_n_(0:fpm_rank-1))
      fpm_n_(0) = n_w_max
      allocate(fpm_inembed_(0:fpm_rank-1))
      fpm_inembed_(0) = n_w_max
      allocate(fpm_onembed_(0:fpm_rank-1))
      fpm_onembed_(0) = n_w_max
      fpm_istride = 1
      fpm_idist = n_w_max
      fpm_ostride = 1
      fpm_odist = n_w_max
      n_fpm = fpm_howmany*fpm_n_(0)
      allocate(fpm_in1_(0:n_fpm-1))
      allocate(fpm_out_(0:n_fpm-1))     
      timing_tic = omp_get_wtime()
c$$$c$$$      We do not actually use the forward plan
c$$$      call dfftw_plan_many_dft(fpm_frwd,fpm_rank,fpm_n_,fpm_howmany
c$$$     $     ,fpm_in1_,fpm_inembed_,fpm_istride,fpm_idist,fpm_out_
c$$$     $     ,fpm_onembed_,fpm_ostride,fpm_odist,FFTW_FORWARD
c$$$     $     ,FFTW_MEASURE)
      call dfftw_plan_many_dft(fpm_back,fpm_rank,fpm_n_,fpm_howmany
     $     ,fpm_in1_,fpm_inembed_,fpm_istride,fpm_idist,fpm_out_
     $     ,fpm_onembed_,fpm_ostride,fpm_odist,FFTW_BACKWARD
     $     ,FFTW_MEASURE)
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
         call test_innerproduct_fast_T_S_q_0(verbose,n_delta_x ,delta_x_
     $        ,n_delta_y,delta_y_,fftw_plan_frwd_,fftw_in1_ ,fftw_out_
     $        ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_ ,n_r
     $        ,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S,S_p__,T_S_q__)
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
         call test_innerproduct_fast_Z_S_q_0(verbose,n_svd_r,svd_r_
     $        ,n_svd_l,svd_l_,svd_s_,svd_V_r_,fftw_plan_frwd_,fftw_in1_
     $        ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $        ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S
     $        ,S_p__,Z_S_q__)
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
      call test_innerproduct_fast_T_R_CTF_M_q_0(verbose ,delta_x_est_
     $     ,delta_y_est_,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $     ,fftw_in1_ ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $     ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,M_p_,M_q_,n_M,ld_M
     $     ,M_p__,CTF_p_,n_CTF ,ld_CTF,CTF_p__,C_M_,T_R_CTF_M_q__)
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
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),T_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_T_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
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
      call dfftw_block_many_0(verbose,.false.,fpm_howmany,n_w_max
     $     ,fpm_in1_,fpm_out_,0,n_transf*n_S*n_M,S_T_T_R_CTF_M_q__)
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
      call dfftw_block_many_0(verbose,.true.,fpm_howmany,n_w_max
     $     ,fpm_in1_,fpm_out_,fpm_back,n_transf*n_S*n_M
     $     ,S_T_T_R_CTF_M_q__)
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
         call test_innerproduct_fast_STxTRM_0(verbose,n_delta_x ,delta_x_
     $        ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,delta_x_est_
     $        ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $        ,fftw_in1_ ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $        ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S
     $        ,S_p__,M_p_,M_q_,n_M ,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF
     $        ,CTF_p__,C_M_,CTF_R_S_ ,S_T_T_R_CTF_M_q__)
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
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
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
      call dfftw_block_many_0(verbose,.false.,fpm_howmany,n_w_max
     $     ,fpm_in1_,fpm_out_,fpm_back,n_transf*n_S*n_M
     $     ,S_Z_T_R_CTF_M_q__)
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
      call dfftw_block_many_0(verbose,.true.,fpm_howmany,n_w_max
     $     ,fpm_in1_,fpm_out_,fpm_back,n_transf*n_S*n_M
     $     ,S_Z_T_R_CTF_M_q__)
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
      call zgemm('T','N',n_delta_x*n_delta_y,n_S*n_M*n_w_max,n_svd_l
     $     ,1.0d0*cmplx(1.0d0,0.0d0),Z_svdd_,n_svd_l,S_Z_T_R_CTF_M_q__
     $     ,n_svd_l,0.0d0*cmplx(1.0d0,0.0d0),S_T_T_R_CTF_M_q__,n_delta_x
     $     *n_delta_y)
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
         call  test_innerproduct_fast_SZxTRM_0(verbose,n_svd_r,svd_r_
     $     ,n_svd_l,svd_l_,svd_s_,svd_V_r_,Z_svdd_,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,delta_x_est_
     $     ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $     ,fftw_in1_ ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $     ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S
     $     ,S_p__,M_p_,M_q_,n_M ,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF
     $     ,CTF_p__,C_M_,CTF_R_S_ ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.ne.2) then

      allocate(ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling test_innerproduct_fast_ZZ_0. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_ZZ_0(verbose,n_delta_x,delta_x_
     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r ,grid_p_,n_w_,n_A
     $     ,n_S,S_alpha_polar_a_,S_alpha_azimu_b_,n_M,polar_a_est_
     $     ,azimu_b_est_,gamma_z_est_ ,delta_x_est_ ,delta_y_est_
     $     ,l2_norm_est_,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_
     $     ,azimu_b_upd_,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_
     $     ,l2_norm_upd_,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f
     $     ,displacement_max ,C_M_,CTF_R_S_ ,S_T_T_R_CTF_M_q__,ZZ_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished ZZ. ' , timing_tot
      timing_tmp = (n_gamma_z*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $     *(n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' ZZ total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      
c$$$      allocate(n_M_per_(0:n_M-1))
c$$$      allocate(I_M_per_(0:n_M-1))
c$$$      n_omp_sub_use = min(n_omp_sub,n_M)
c$$$      if (n_omp_sub_use.eq.1) then
c$$$         write(6,'(A,A,I0)')
c$$$     $        ' Warning! Turning off omp for innerproducts ',
c$$$     $        ' n_omp_sub_use: ',n_omp_sub_use
c$$$      end if !if (n_omp_sub_use.eq.1) then
c$$$      if (verbose.gt.0) then
c$$$         write (6,'(A,I0,A)') ' Setting up n_omp_sub_use = '
c$$$     $        ,n_omp_sub_use,' sub-blocks for omp parallelization'
c$$$         write(6,'(A)') ' '
c$$$      end if ! if (verbose.gt.0) then
c$$$      do nomp_sub=0,n_omp_sub_use-1
c$$$         if (nomp_sub.eq.0) then
c$$$            n_M_per_(0) = n_M/n_omp_sub_use
c$$$            I_M_per_(0) = 0
c$$$         else if (nomp_sub.gt.0 .and. nomp_sub.lt.n_omp_sub_use-1) then
c$$$            n_M_per_(nomp_sub) = n_M /n_omp_sub_use
c$$$            I_M_per_(nomp_sub) = I_M_per_(nomp_sub-1) +
c$$$     $           n_M_per_(nomp_sub-1)
c$$$         else if (nomp_sub.eq.n_omp_sub_use-1) then
c$$$            I_M_per_(n_omp_sub_use-1) = I_M_per_(n_omp_sub_use-2) +
c$$$     $           n_M_per_(n_omp_sub_use-2)
c$$$            n_M_per_(n_omp_sub_use-1) = n_M - I_M_per_(n_omp_sub_use-1)
c$$$         end if ! if nomp_sub
c$$$      enddo ! do nomp_sub=0,n_omp_sub_use-1
c$$$      nM_per_max = max_i4_f(n_omp_sub_use,n_M_per_)
c$$$      if (verbose.gt.1) then
c$$$         write(format_string,'(A,I0,A)') '(A,',n_omp_sub_use ,'(I0,1X))'
c$$$         write(6,format_string) 'n_M_per_: ' ,(n_M_per_(nomp_sub)
c$$$     $        ,nomp_sub=0,n_omp_sub_use-1)
c$$$         write(6,format_string) 'I_M_per_: ' ,(I_M_per_(nomp_sub)
c$$$     $        ,nomp_sub=0,n_omp_sub_use-1)
c$$$         write(6,'(A,I0)') 'nM_per_max: ' , nM_per_max
c$$$      end if ! if (verbose.gt.1) then
c$$$
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
c$$$      end if ! if (verbose.gt.1) then
c$$$      allocate(fftw_plan_frwd__(0:n_r*n_omp_sub_use-1))
c$$$      allocate(fftw_plan_back__(0:n_r*n_omp_sub_use-1))
c$$$      allocate(fftw_in1__(0:n_A*n_omp_sub_use-1))
c$$$      allocate(fftw_out__(0:n_A*n_omp_sub_use-1))
c$$$      do nomp_sub=0,n_omp_sub_use-1
c$$$         n_r_tmp = nomp_sub*n_r
c$$$         n_A_tmp = nomp_sub*n_A
c$$$         na = 0
c$$$         do nr=0,n_r-1
c$$$            call dfftw_plan_dft_1d_(fftw_plan_frwd__(n_r_tmp+nr)
c$$$     $           ,n_w_(nr),fftw_in1__(n_A_tmp+na),fftw_out__(n_A_tmp
c$$$     $           +na),FFTW_FORWARD,FFTW_MEASURE) 
c$$$            call dfftw_plan_dft_1d_(fftw_plan_back__(n_r_tmp+nr)
c$$$     $           ,n_w_(nr),fftw_out__(n_A_tmp+na),fftw_in1__(n_A_tmp
c$$$     $           +na),FFTW_BACKWARD,FFTW_MEASURE) 
c$$$            na = na + n_w_(nr)
c$$$         enddo ! do nr=0,n_r-1
c$$$      enddo ! do nomp_sub=0,n_omp_sub_use-1
c$$$
c$$$      allocate(SM_use_sub_(0:n_omp_sub_use-1))
c$$$      timing_tic = omp_get_wtime()
c$$$c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp)
c$$$c$OMP DO 
c$$$      do nomp_sub=0,n_omp_sub_use-1
c$$$         nM_per = n_M_per_(nomp_sub)
c$$$         IM_per = I_M_per_(nomp_sub)
c$$$         n_r_tmp = nomp_sub*n_r
c$$$         n_A_tmp = nomp_sub*n_A
c$$$         if (verbose.gt.1) then
c$$$            write(6,'(A,A,A,I0)') 'calling'
c$$$     $           ,' test_innerproduct_batch_stage_1'
c$$$     $           ,' with thread: ',omp_get_thread_num()
c$$$         end if
c$$$         call test_innerproduct_batch_SM_stage_1a(verbose
c$$$     $        ,svd_calculation_type,eps_svd,n_r,grid_p_,n_w_
c$$$     $        ,half_diameter_x_c ,n_S,ld_S ,S_p__
c$$$     $        ,tesselation_distance_req,S_alpha_polar_a_
c$$$     $        ,S_alpha_azimu_b_,L_ ,nl_max,nm_sum ,ll_sum,T_nl_,T_vm_
c$$$     $        ,T_tr_,T_ll_ ,T_lf_ ,T_c0_,T_c1_,T_c2_ ,T_c3_,T_ls_,T_LT_
c$$$     $        ,n_T_root_base ,T_root_base_,nM_per,ld_M ,M_p__(ld_M
c$$$     $        *IM_per) ,n_CTF,ld_CTF ,CTF_p__ ,polar_a_est_(IM_per)
c$$$     $        ,azimu_b_est_(IM_per) ,gamma_z_est_(IM_per)
c$$$     $        ,delta_x_est_(IM_per) ,delta_y_est_(IM_per)
c$$$     $        ,l2_norm_est_(IM_per) ,ctf_ind_est_(IM_per)
c$$$     $        ,S_index_est_(IM_per) ,polar_a_upd_(IM_per)
c$$$     $        ,azimu_b_upd_(IM_per) ,gamma_z_upd_(IM_per)
c$$$     $        ,delta_x_upd_(IM_per) ,delta_y_upd_(IM_per)
c$$$     $        ,l2_norm_upd_(IM_per) ,ctf_ind_upd_(IM_per)
c$$$     $        ,S_index_upd_(IM_per), alpha_update_f ,N_pixels_in
c$$$     $        ,displacement_max,n_delta_x,n_delta_y ,n_gamma_z ,
c$$$     $        C_M_(IM_per),fftw_plan_frwd__(n_r_tmp)
c$$$     $        ,fftw_plan_back__(n_r_tmp) ,fftw_in1__(n_A_tmp)
c$$$     $        ,fftw_out__(n_A_tmp),n_svd_r,n_svd_d ,n_svd_l,svd_r_
c$$$     $        ,svd_d_ ,svd_l_,svd_U_d_,svd_s_,svd_V_r_ ,Z_svdd_,CTF_R_S_
c$$$     $        ,SM_use_sub_(nomp_sub))
c$$$      enddo !do nomp_sub=0,n_omp_sub_use-1
c$$$c$OMP END DO
c$$$c$OMP END PARALLEL
c$$$      n_SM_use = sum_i4_f(n_omp_sub_use,SM_use_sub_)
c$$$      timing_toc = omp_get_wtime()
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A,A,A,F8.5)') ' Finished',
c$$$     $        ' test_innerproduct_batch_stage_1:', 
c$$$     $        ' total time ',
c$$$     $        timing_toc-timing_tic
c$$$         write(6,'(A)') ' '
c$$$      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
c$$$      deallocate(SM_use_sub_)
      deallocate(ZZ_)
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
      deallocate(S_p_)
      deallocate(S_q_)
      deallocate(M_p_)
      deallocate(M_q_)
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
      call dfftw_destroy_plan_(fpm_back)
c$$$      call dfftw_destroy_plan_(fpm_frwd)
      deallocate(fftw_out_)
      deallocate(fftw_in1_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_plan_frwd_)
      deallocate(fpm_n_)
      deallocate(fpm_inembed_)
      deallocate(fpm_onembed_)
      deallocate(fpm_in1_)
      deallocate(fpm_out_)
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ' Destroying fftw_plans for omp.'
c$$$      end if
c$$$      do nr=0,n_omp_sub_use*n_r-1
c$$$         call dfftw_destroy_plan_(fftw_plan_back__(nr))
c$$$         call dfftw_destroy_plan_(fftw_plan_frwd__(nr))
c$$$      enddo
c$$$      deallocate(fftw_out__)
c$$$      deallocate(fftw_in1__)
c$$$      deallocate(fftw_plan_back__)
c$$$      deallocate(fftw_plan_frwd__)
c$$$      deallocate(n_M_per_)
c$$$      deallocate(I_M_per_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating svd arrays.'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_fast_stage_0]'
      end if

      end


