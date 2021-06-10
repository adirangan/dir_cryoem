c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_batch_SM_stage_1a(verbose
     $     ,svd_calculation_type_in,eps_svd,n_r,grid_p_,n_w_
     $     ,half_diameter_x_c ,n_S,ld_S ,S_p_,tesselation_distance_req
     $     ,S_alpha_polar_a_ ,S_alpha_azimu_b_,L_ ,nl_max,nm_sum ,ll_sum
     $     ,T_nl_,T_vm_,T_tr_ ,T_ll_ ,T_lf_,T_c0_ ,T_c1_,T_c2_ ,T_c3_
     $     ,T_ls_,T_LT_ ,n_T_root_base ,T_root_base_ ,n_M,ld_M,M_p_
     $     ,n_CTF,ld_CTF ,CTF_p_,polar_a_est_ ,azimu_b_est_,gamma_z_est_
     $     ,delta_x_est_ ,delta_y_est_ ,l2_norm_est_,ctf_ind_est_
     $     ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_,gamma_z_upd_
     $     ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_,ctf_ind_upd_
     $     ,S_index_upd_ ,alpha_update_f,n_pixels_in
     $     ,displacement_max ,n_delta_x ,n_delta_y ,n_gamma_z,C_M_ ,
     $     fftw_plan_frwd_ ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_
     $     ,n_svd_r ,n_svd_d ,n_svd_l,svd_r_,svd_d_ ,svd_l_ ,svd_U_d_
     $     ,svd_s_ ,svd_V_r_ ,Z_svdd_,C_S_,n_SM_use)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This function sets up a local-search across image-template pairs, ;
c$$$      calculating innerproducts for a (possibly strict) subset of all ;
c$$$      image-template pairs. ;
c$$$      
c$$$      The governing parameters for this local-search are: ;
c$$$      tesselation_distance_req (real *8, passed in as input), and
c$$$      n_LT_hlf (integer *4, internally defined below). ;
c$$$      
c$$$      The basic steps behind the local-search are: ;
c$$$      (1). For each image: ;
c$$$      (2). Construct a list of all templates within a required distance. ;
c$$$      (3). Add a few randomly distributed templates to this list. ;
c$$$      (3). Calculate innerproducts for all image-template pairs within this list. ;
c$$$      (4). Look at the top few innerproducts (i.e., best templates). ;
c$$$      (5). Check for nearby templates within a required distance of those few best templates. ;
c$$$      (6). Calculate innerproducts for those nearby templates as well. ;
c$$$      (7). Move on to the next image, going back to step (1) above. ;
c$$$      
c$$$      The details are explained in comments below. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      integer svd_calculation_type_in,svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,ld_S,n_M,ld_M,n_CTF,ld_CTF
      real *8 grid_p_(0:n_r-1),half_diameter_x_c
      complex *16 S_p_(0:ld_S*n_S-1)
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
      real *8, allocatable :: vp_input_(:)
      real *8 tesselation_distance_req,tesselation_distance_req_use
      integer *4 n_LT,nLT
      integer *4, allocatable :: LT_(:)
      integer *4 n_SM_use
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for adaptive sampling
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      logical, allocatable :: S_use_(:)
      integer *4 n_S_use,n_add
      integer *4 n_pass
      logical continue_flag
      integer *4 n_LT_tmp,n_LT_srt,n_LT_use
      integer *4 n_LT_min,n_LT_hlf,ns_tmp
      parameter (n_LT_min = 12)
      complex *16, allocatable :: C_Z_srt_(:)
      complex *16 C_Z_min,C_S_min
      integer *4, allocatable :: LT_srt_(:)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      complex *16 M_p_(0:ld_M*n_M-1)
      complex *16 CTF_p_(0:0)
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
      real *8 eps_svd,n_pixels_in
      real *8 displacement_max,delta_x_tmp,delta_y_tmp
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8, allocatable ::  delta_x_SM_(:)
      real *8, allocatable ::  delta_y_SM_(:)
      real *8, allocatable ::  gamma_z_SM_(:)
      complex *16, allocatable ::  C_S_SM_(:)
      complex *16, allocatable ::  C_Z_SM_(:)
      real *8, allocatable ::  delta_x_sort_SM_(:)
      real *8, allocatable ::  delta_y_sort_SM_(:)
      real *8, allocatable ::  gamma_z_sort_SM_(:)
      complex *16, allocatable ::  C_S_sort_SM_(:)
      complex *16, allocatable ::  C_Z_sort_SM_(:)
      integer, allocatable :: I_permute_SM_(:)
      integer, allocatable :: I_inverse_SM_(:)
      complex *16 C_M_(0:n_M-1)
      integer *4 n_svd_r,n_svd_d,n_svd_l
      real *8 svd_r_(0:n_svd_r-1)
      real *8 svd_d_(0:n_svd_d-1)
      integer *4 svd_l_(0:n_svd_l-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      real *8 svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      complex *16 Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1)
      complex *16 C_S_(0:n_gamma_z*n_S*n_CTF-1)
      integer *4 nsvd_r,nsvd_d,nsvd_l
      logical warning_flag
      data warning_flag / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels
      complex *16, allocatable :: Z_svdr_(:)
c$$$      Z_svdd_ should be passed as input
c$$$      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer svd_l
c$$$      fftw plans and workspace should be passed in as input
c$$$      integer *8, allocatable :: fftw_plan_frwd_(:)
c$$$      integer *8, allocatable :: fftw_plan_back_(:)
c$$$      complex *16, allocatable :: fftw_in1_(:)
c$$$      complex *16, allocatable :: fftw_out_(:)
      integer *8 fftw_plan_frwd_(0:0)
      integer *8 fftw_plan_back_(0:0)
      complex *16 fftw_in1_(0:0)
      complex *16 fftw_out_(0:0)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*)
      complex *16 fftw_out_last_(*)
c$$$      indices
      integer ns,ns_use,nr_use,nm,nctf,nr,n_w_max,n_A,na
      complex *16, allocatable :: M_p(:)
      complex *16, allocatable :: S_p(:)
      complex *16, allocatable :: S_q(:)
      complex *16, allocatable :: M_q(:)
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x_optimal,delta_y_optimal,gamma_z_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      character(len=64) format_string
c$$$      array of innerproducts to hold measurements
      complex *16 C_S,C_M,C_Z,C_S_optimal,C_Z_optimal,Z_q
      complex *16, allocatable :: Z_q_(:)
      complex *16 Z_tmp
      real *8 pi
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      real *8 timing_tic_0,timing_toc_0,timing_tot_0
      real *8 timing_tic_1,timing_toc_1,timing_tot_1
      real *8 timing_tic_2,timing_toc_2,timing_tot_2
      real *8 timing_tic_3,timing_toc_3,timing_tot_3
      real *8 timing_tic_4,timing_toc_4,timing_tot_4

      if (0+verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_SM_stage_1a]'
      end if      
      if (0+verbose.gt.1) then
         write(6,'(A,I0,A,F5.3,A,F3.1,A,I0,A,I0,A,I0)') 'n_r: ' ,n_r
     $        ,'; eps_svd: ',eps_svd,'; n_pixels_in: ',n_pixels_in
     $        ,'; n_delta_x: ',n_delta_x ,'; n_delta_y: ',n_delta_y
     $        ,'; n_gamma_z: ',n_gamma_z
      end if

      if (0+verbose.gt.1) then
         write(6,'(A)') 'Generating indices'
      end if
      pi = 4.0*atan(1.0)
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(n_r-1)
      if (0+verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_A: ',n_A,'; n_w_max: ',n_w_max
         write(format_string,'(A,I0,A)') '(A,',n_r,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,n_r-1)
      end if

c$$$      fftw plans and workspace should be passed in as input
c$$$      if (0+verbose.gt.1) then
c$$$         write(6,'(A)') 'Generating fftw_plans'
c$$$      end if
c$$$      allocate(fftw_plan_frwd_(0:n_r-1))
c$$$      allocate(fftw_plan_back_(0:n_r-1))
c$$$      allocate(fftw_in1_(0:n_A-1))
c$$$      allocate(fftw_out_(0:n_A-1))
c$$$      na = 0
c$$$      do nr=0,n_r-1
c$$$         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
c$$$     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_ESTIMATE) 
c$$$         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
c$$$     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_ESTIMATE) 
c$$$         na = na + n_w_(nr)
c$$$      enddo
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*n_pixels_in or so.'
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      do ndx=0,n_delta_x-1
         if (n_delta_x.gt.1) then
            delta_x = (-n_pixels_in + ndx*2*n_pixels_in/(n_delta_x-1))
     $           /n_r*half_diameter_x_c
         else
            delta_x = 0.0d0
         end if
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            if (n_delta_y.gt.1) then
               delta_y = (-n_pixels_in + ndy*2*n_pixels_in/(n_delta_y
     $              -1))/n_r*half_diameter_x_c
            else
               delta_y = 0.0d0
            end if
            delta_y_(ndy) = delta_y
         enddo
      enddo
      if (0+verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)')
     $        ' ncur or so.'
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      do ngz=0,n_gamma_z-1
         if (n_gamma_z.gt.1) then
            gamma_z = (2*pi*ngz)/n_gamma_z
         else
            gamma_z = 0.0d0
         end if
         gamma_z_(ngz) = gamma_z
      enddo
      if (0+verbose.gt.1) then
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

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if
      R_max = 2*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
            delta = dsqrt(delta_x_(ndx)**2 + delta_y_(ndy)**2)
            if (delta.gt.delta_max) then
               delta_max = delta
            end if
         enddo
      enddo
      n_pixels = delta_max/dsqrt(2.0d0)*2*K_max
      if (0+verbose.gt.1) then
         write(6,'(A,F8.3,A,F8.3)') 'R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') 'K_max: ',K_max,'; n_pixels: '
     $        ,n_pixels
      end if
      if (K_max.gt.96 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, K_max ',K_max
     $        ,' too large in test_innerproduct_batch_SM_stage_1a'
      end if
      if (n_pixels.gt.3 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, n_pixels ',n_pixels
     $        ,' too large in test_innerproduct_batch_SM_stage_1a'
      end if
      if (eps_svd.lt.0.1d0 .and. warning_flag) then
         write(6,'(A,F8.5,A)') 'Warning, eps_svd ',eps_svd
     $        ,' too small in test_innerproduct_batch_SM_stage_1a'
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      The svd parameters should be passed in as input
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      if (.false.) then
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdpick.txt'
c$$$      end if
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdload.txt'

      if (0+verbose.gt.1) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
      end if

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Selecting svd_calculation_type: '
      end if
      if (svd_calculation_type_in.lt.0) then
         if (.false.) then
         else if (n_svd_l.gt.n_delta_x*n_delta_y) then
            svd_calculation_type = 0
         else if (n_svd_l.le.n_delta_x*n_delta_y) then
            svd_calculation_type = 1
         else
            svd_calculation_type = 2
         end if !switch
      else if (svd_calculation_type_in.ge.0) then
         svd_calculation_type = max(0,min(3,svd_calculation_type_in))
         if (0+verbose.gt.1) then
            write(6,'(A,I0)') 'forcing: svd_calculation_type = '
     $           ,svd_calculation_type
         end if !verbose
      end if !choose svd_calculation_type
      if (0+verbose.gt.1) then
         if (.false.) then
         else if (svd_calculation_type.eq.0) then
            write(6,'(A)') 'svd_calculation_type 0: multiply then fft.'
         else if (svd_calculation_type.eq.1) then
            write(6,'(A)') 'svd_calculation_type 1: fft then multiply.'
         else if (svd_calculation_type.eq.2) then
            write(6,'(A,A)') 'svd_calculation_type 2: fft for gamma,' ,
     $           ' brute-force displacements.'
         else if (svd_calculation_type.eq.3) then
            write(6,'(A,A,A)') 'svd_calculation_type 2:' ,
     $           ' brute-force for gamma,'
     $           ,' brute-force displacements.'
         end if
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      The operator Z_svdd_ should be passed in as input
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      if (0+verbose.gt.1) then
c$$$         write(6,'(A)') ' Setting up displacement-operator Z_svdd_'
c$$$         write(6,'(A)') ' associated with svd-expansion.'
c$$$      end if
c$$$      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
c$$$      timing_tic = omp_get_wtime()
c$$$      call innerproduct_q_k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
c$$$     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
c$$$      timing_toc = omp_get_wtime()
c$$$      if (0+verbose.gt.1) then
c$$$         write(6,'(A,F8.5)')
c$$$     $        ' finished innerproduct_q_k_svdd: total time ',timing_toc
c$$$     $        -timing_tic
c$$$      end if

      if (0+verbose.gt.1) then
         write(6,'(A,A)') ' Allocating template-dependent operator'
     $        ,' Z_svdr_ used in svd-expansion.'
      end if
      allocate(Z_svdr_(0:n_svd_l*n_w_max-1))
      if (0+verbose.gt.1) then
         write(6,'(A,A)') ' Allocating temporary array Z_tmpC_ '
     $        ,'to hold innerproducts associated with '
         write(6,'(A,A)') ' the displacements across all rotations '
     $        ,'for any single image-template pair'
      end if
      allocate(Z_tmpC_(0:n_w_max*n_delta_x*n_delta_y-1))
      if (0+verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_q to hold bessel-function'
         write(6,'(A)') ' expansion of each template within S_p.'
         write(6,'(A)') ' Note that S_q holds the template alone,'
         write(6,'(A)') ' without accounting for the ctf.'
         write(6,'(A)') ' The ctf will be accounted for in M_q.'
      end if
      allocate(S_p(0:n_A-1))
      allocate(S_q(0:n_A-1))
      if (0+verbose.gt.1) then
         write(6,'(A,A)') ' Allocating arrays M_p'
     $        ,' to hold transformed images.'
      end if
      allocate(M_p(0:n_A-1))

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Allocating array M_q to hold a single'
         write(6,'(A)') ' bessel-function-expansion of an image.'
      end if
      allocate(M_q(0:n_A-1))
      if (0+verbose.gt.1) then
         write(6,'(A)') ' Allocating array Z_q_ to hold the'
         write(6,'(A)') ' innerproducts associated with displacements'
         write(6,'(A)') ' and rotations for any single '
         write(6,'(A)') ' image-template pair. '
      end if
      allocate(Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1))

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Allocating array vp_input_ to hold a point,'
         write(6,'(A)') ' and array LT_ to hold a list of'
         write(6,'(A)') ' unique (template) indices. Both are used by '
         write(6,'(A)') ' tesselation_neighborhood_wrapper_0. '
      end if
      allocate(vp_input_(0:3-1))
      allocate(LT_(0:2*ll_sum-1))
      allocate(LT_srt_(0:n_S-1))
      allocate(C_Z_srt_(0:n_S-1))
      allocate(S_use_(0:n_S-1))
      allocate(C_Z_SM_(0:n_S-1))
      allocate(C_S_SM_(0:n_S-1))
      allocate(delta_x_SM_(0:n_S-1))
      allocate(delta_y_SM_(0:n_S-1))
      allocate(gamma_z_SM_(0:n_S-1))
      allocate(C_Z_sort_SM_(0:n_S-1))
      allocate(C_S_sort_SM_(0:n_S-1))
      allocate(delta_x_sort_SM_(0:n_S-1))
      allocate(delta_y_sort_SM_(0:n_S-1))
      allocate(gamma_z_sort_SM_(0:n_S-1))
      allocate(I_permute_SM_(0:n_S-1))
      allocate(I_inverse_SM_(0:n_S-1))

      if (verbose.gt.2) then
         do ns=0,n_S-1
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_p_(ns*ld_S)
     $           ,S_p_(ns*ld_S),Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_p(',ns,')|: ' , Z_tmp
         enddo                  ! do ns=0,n_S-1
      end if                    !if (verbose.gt.2) then


      if (0+verbose.gt.1) then
         write(6,'(A)') ' Calculating innerproducts...'
      end if
      timing_tot_0 = 0.0d0
      timing_tot_1 = 0.0d0
      timing_tot_2 = 0.0d0
      timing_tot_3 = 0.0d0
      timing_tot_4 = 0.0d0
      n_SM_use = 0
      do nm=0,n_M-1
         do ns=0,n_S-1
            delta_x_SM_(ns) = 0.0d0
            delta_y_SM_(ns) = 0.0d0
            gamma_z_SM_(ns) = 0.0d0
            C_S_SM_(ns) = (1.0d0,0.0d0)
            C_Z_SM_(ns) = (0.0d0,0.0d0)
         enddo !do ns=0,n_S-1
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         if (0+verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0,A,I0)') 'Processing image ',nm,'/'
     $           ,n_M,' associated with ctf ',nctf,'/',n_CTF
            write(6,'(A)') ' First we calculate l2_norm C_M.'
         end if
         timing_tic_0 = omp_get_wtime()
         call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_p_(nm*ld_M),M_p_(nm
     $        *ld_M),C_M)
         C_M = zsqrt(C_M)/(n_r*n_r)
         C_M_(nm) = C_M
         if (verbose.gt.2) then
            write(6,'(A,I0,A,2F16.3)') ' nm: ',nm,'; C_M: ',C_M
         end if
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply ctf-star to M_p '
         end if
         call cp1_c16(n_A,M_p_(nm*ld_M),M_p)
         call xc1_c16(n_A,M_p,CTF_p_(nctf*ld_CTF),M_p)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,M_p,M_q)
         timing_toc_0 = omp_get_wtime()
         timing_tot_0 = timing_tot_0 + (timing_toc_0 - timing_tic_0)

         if (0+verbose.gt.1) then
            write(6,'(A)') ' Resetting S_use_ to false.'
         end if !if (0+verbose.gt.1) then
         do ns=0,n_S-1
            S_use_(ns) = .false.
         enddo !do ns=0,n_S-1
         n_S_use = 0
         C_Z_min = cmplx( 1.0d0 , 0.0d0)
         C_S_min = cmplx( 1.0d0 , 0.0d0)

         if (verbose.gt.2) then
            write(6,'(A)') ' Now we calculate innerproducts.'
         end if

         if (tesselation_distance_req.ge.2.0d0) then
            if (0+verbose.gt.1) then
               write(6,'(A)') ' considering all templates'
            end if !if (0+verbose.gt.1) then
            do ns=0,n_S-1
            ns_use = ns
            S_use_(ns) = .true.
            n_S_use = n_S_use + 1
            call test_innerproduct_batch_SM_stage_2a(verbose ,nm,ns_use
     $           ,n_r,grid_p_,n_w_,n_S,ld_S,n_A,nctf,n_w_max ,S_p_,S_p
     $           ,S_q,M_p,M_q,Z_q_ ,C_M,C_S_,C_S_optimal ,C_Z_optimal
     $           ,svd_calculation_type ,n_svd_r,svd_r_,n_svd_l,svd_l_
     $           ,svd_s_,svd_V_r_ ,Z_svdd_ ,Z_svdr_,Z_tmpC_
     $           ,displacement_max ,delta_x_est,delta_y_est ,gamma_z_est
     $           ,n_delta_x,n_delta_y,n_gamma_z ,delta_x_ ,delta_y_
     $           ,gamma_z_ ,delta_x_est_,delta_y_est_ ,gamma_z_est_
     $           ,ndx_optimal,ndy_optimal,ngz_optimal ,delta_x_optimal
     $           ,delta_y_optimal,gamma_z_optimal ,fftw_plan_frwd_
     $           ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_ ,timing_tic_1
     $           ,timing_toc_1,timing_tot_1 ,timing_tic_2 ,timing_toc_2
     $           ,timing_tot_2 ,timing_tic_3,timing_toc_3 ,timing_tot_3
     $           ,timing_tic_4,timing_toc_4,timing_tot_4)
            delta_x_SM_(ns_use) = delta_x_optimal
            delta_y_SM_(ns_use) = delta_y_optimal
            gamma_z_SM_(ns_use) = gamma_z_optimal
            C_S_SM_(ns_use) = C_S_optimal
            C_Z_SM_(ns_use) = C_Z_optimal
            if (real(C_Z_optimal).lt.real(C_Z_min)) then
               C_Z_min = C_Z_optimal
            end if !if (real(C_Z_optimal).lt.real(C_Z_min)) then
            if (real(C_S_optimal).lt.real(C_S_min)) then
               C_S_min = C_S_optimal
            end if !if (real(C_S_optimal).lt.real(C_S_min)) then
            enddo !do ns=0,n_S-1
            goto 5
         end if !if (tesselation_distance_req.ge.2.0d0) then

c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$         Otherwise, if (tesselation_distance_req.lt.2.0d0) then...
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if (0+verbose.gt.1) then
            write(6,'(A,I0,A)') ' Generating initial list for image ',nm
     $           ,' based on estimated '
            write(6,'(A,F8.4,A,F8.4,A)') ' polar_a ',polar_a_est_(nm)
     $           ,' and azimu_b ',azimu_b_est_(nm), '.'
            write(6,'(A)') ' '
            write(6,'(A)') ' This initial list comprises only the '
            write(6,'(A)') ' templates that have a polar_a and an '
            write(6,'(A)') ' azimu_b that are within distance  '
            write(6,'(A)') ' tesselation_distance_req of the image '
            write(6,'(A)') ' angles (as considered as vectors on '
            write(6,'(A)') ' the sphere). '
            write(6,'(A)') ' '
            write(6,'(A)') ' Note that if there are no templates '
            write(6,'(A)') ' within the requested distance, we  '
            write(6,'(A)') ' broaden our search until there are '
            write(6,'(A)') ' at least two templates in range. '
            write(6,'(A)') ' '
            write(6,'(A)') ' Note also that if you wish to perform a '
            write(6,'(A)') ' global search, simply set the required '
            write(6,'(A)') ' tesselation_distance_req to encompass '
            write(6,'(A)') ' the entire sphere (i.e., .ge. 2.0d0). '
            write(6,'(A)') ' '
            write(6,'(A)') ' Note also that we track which templates ' 
            write(6,'(A)') ' have been used for this particular image. '
            write(6,'(A)') ' (see temporary logical array S_use_, '
            write(6,'(A)') ' as well as temporary indexing array '
            write(6,'(A)') ' LT_). '
            write(6,'(A)') ' '
         end if !if (0+verbose.gt.1) then
         call get_angle_to_vp_(polar_a_est_(nm),azimu_b_est_(nm)
     $        ,vp_input_)
         n_LT = 0
         tesselation_distance_req_use = tesselation_distance_req
         do while (n_S.gt.0 .and. n_LT.lt.min(n_S,2))
            call tesselation_neighborhood_wrapper_0(n_S,L_ ,nl_max
     $           ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_
     $           ,T_c1_ ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base
     $           ,T_root_base_ ,vp_input_,tesselation_distance_req_use
     $           ,n_LT ,LT_(0))
            if (n_LT.lt.min(n_S,2)) then
               if (verbose.gt.1) then
                  write(6,'(A,F8.4)') ' increasing '
     $                 ,tesselation_distance_req_use
               end if !if (verbose.gt.1) then
               tesselation_distance_req_use =
     $              tesselation_distance_req_use + 0.0125
            end if !if (n_LT.lt.min(n_S,2)) then
         enddo !do while (n_S.gt.0 .and. n_LT.lt.min(n_S,2))
         if (0+verbose.gt.1) then
            write(6,'(A,I0,A)') ' found ',n_LT,' templates.'
         end if !if (0+verbose.gt.1) then
         do nLT=0,n_LT-1
            S_use_(LT_(nLT))=.true.
            n_S_use = n_S_use+1
         enddo !do nLT=0,n_LT-1
         n_LT_hlf = n_S_use/2
         n_add = max(n_LT_min,n_LT_hlf)
         n_add = max(0,min(n_S-n_S_use,n_add))
         if (0+verbose.gt.1) then
            write(6,'(A,I0)') ' n_LT: ' , n_LT
            write(6,'(A,I0)') ' n_S_use: ' , n_S_use
            write(6,'(A,I0)') ' n_LT_hlf: ' , n_LT_hlf
            write(6,'(A,I0,A)') ' adding additional n_add ',n_add
     $           ,' templates.'
         end if !if (0+verbose.gt.1) then
         if (0+verbose.gt.1) then
            write(6,'(A)') ' '
            write(6,'(A)') ' This step adds an additional number of'
            write(6,'(A)') ' randomly chosen (unique) templates from '
            write(6,'(A)') ' the list of templates. These additional '
            write(6,'(A)') ' templates will not necessarily  be close '
            write(6,'(A)') ' to the image in terms of polar_a and '
            write(6,'(A)') ' azimu_b. '
            write(6,'(A)') ' '
            write(6,'(A)') ' The number of additional templates '
            write(6,'(A)') ' is half the number originally considered. '
            write(6,'(A)') ' (i.e., half the number within '
            write(6,'(A)') ' tesselation_distance_req of the original '
            write(6,'(A)') ' estimated image-angle). '
            write(6,'(A)') ' If you wish to increase this number, '
            write(6,'(A)') ' modify the variable: n_LT_hlf above. '
            write(6,'(A)') ' '
            write(6,'(A)') ' Note that we always try to add at least '
            write(6,'(A)') ' a dozen templates if possible '
            write(6,'(A)') ' (see parameter n_LT_min above).'
            write(6,'(A)') ' '
            write(6,'(A)') ' Note that these extra templates are also '
            write(6,'(A)') ' accounted for in arrays S_use_ and LT_. '
            write(6,'(A)') ' '
         end if !if (0+verbose.gt.1) then
         call randinclude(n_LT,LT_,n_S,S_use_,n_S_use,n_add)
         if (0+verbose.gt.1) then
            write(6,'(A,I0,A,I0)') ' now using n_S_use: '
     $           ,n_S_use,' n_LT_: ' , n_LT
         end if !if (0+verbose.gt.1) then
         
         if (0+verbose.gt.1) then
            write(6,'(A)') ' '
            write(6,'(A)') ' At this point we step through the '
            write(6,'(A)') ' list of templates that has been '
            write(6,'(A)') ' generated for this particular '
            write(6,'(A)') ' image (i.e., including both nearby '
            write(6,'(A)') ' templates and a few randomly '
            write(6,'(A)') ' chosen templates from elsewhere on'
            write(6,'(A)') ' the sphere). '
            write(6,'(A)') ' '
         end if !if (0+verbose.gt.1) then
         do nLT=0,n_LT-1
            ns_use = LT_(nLT)            
            if (0+verbose.gt.1) then
               write(6,'(A,I0,A,I0)')
     $              '  Calculating innerproduct for image ' , nm ,
     $              ' and template ' , ns_use
            end if !if (0+verbose.gt.1) then
            call test_innerproduct_batch_SM_stage_2a(verbose ,nm,ns_use
     $           ,n_r,grid_p_,n_w_,n_S,ld_S,n_A,nctf,n_w_max ,S_p_,S_p
     $           ,S_q,M_p,M_q,Z_q_ ,C_M,C_S_,C_S_optimal ,C_Z_optimal
     $           ,svd_calculation_type ,n_svd_r,svd_r_,n_svd_l,svd_l_
     $           ,svd_s_,svd_V_r_ ,Z_svdd_ ,Z_svdr_,Z_tmpC_
     $           ,displacement_max ,delta_x_est,delta_y_est ,gamma_z_est
     $           ,n_delta_x,n_delta_y,n_gamma_z ,delta_x_ ,delta_y_
     $           ,gamma_z_ ,delta_x_est_,delta_y_est_ ,gamma_z_est_
     $           ,ndx_optimal,ndy_optimal,ngz_optimal ,delta_x_optimal
     $           ,delta_y_optimal,gamma_z_optimal ,fftw_plan_frwd_
     $           ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_ ,timing_tic_1
     $           ,timing_toc_1,timing_tot_1 ,timing_tic_2 ,timing_toc_2
     $           ,timing_tot_2 ,timing_tic_3,timing_toc_3 ,timing_tot_3
     $           ,timing_tic_4,timing_toc_4,timing_tot_4)
            delta_x_SM_(ns_use) = delta_x_optimal
            delta_y_SM_(ns_use) = delta_y_optimal
            gamma_z_SM_(ns_use) = gamma_z_optimal
            C_S_SM_(ns_use) = C_S_optimal
            C_Z_SM_(ns_use) = C_Z_optimal
         enddo !do nLT=0,n_LT-1

         if (0+verbose.gt.1) then
            write (6,'(A)') ' After calculating the initial list of '
            write (6,'(A)') ' innerproducts for this particular image '
            write (6,'(A)') ' (see above), we pass through this list '
            write (6,'(A)') ' several times: '
            write (6,'(A)') ' '
            write (6,'(A)') ' Step 1: '
            write (6,'(A)') ' Each time we pass through the list we '
            write (6,'(A)') ' find the templates that correspond to '
            write (6,'(A)') ' the largest n_LT_hlf innerproducts. '
            write (6,'(A)') ' Note that this number can certainly be '
            write (6,'(A)') ' increased or decreased (just by '
            write (6,'(A)') ' modifying or replacing n_LT_hlf within '
            write (6,'(A)') ' the following code). '
            write (6,'(A)') ' '
            write (6,'(A)') ' Step 2: '
            write (6,'(A)') ' Once we find these "good" templates, '
            write (6,'(A)') ' we search for other templates nearby. '
            write (6,'(A)') ' As before, we use the distance '
            write (6,'(A)') ' tesselation_distance_req_use, which '
            write (6,'(A)') ' was previously used to generate our '
            write (6,'(A)') ' original list of templates. '
            write (6,'(A)') ' '
            write (6,'(A)') ' Now (hopefully) many of these nearby '
            write (6,'(A)') ' templates close to the "good" templates '
            write (6,'(A)') ' will have already been considered for '
            write (6,'(A)') ' this particular image. (i.e, many '
            write (6,'(A)') ' nearby templates will alread be listed '
            write (6,'(A)') ' in LT_ and tagged in S_use_). '
            write (6,'(A)') ' However, some of these nearby templates '
            write (6,'(A)') ' will be new. '
            write (6,'(A)') ' '
            write (6,'(A)') ' Step 3: '
            write (6,'(A)') ' We add each of these new nearby '
            write (6,'(A)') ' templates to our list of templates '
            write (6,'(A)') ' (updating LT_ and S_use_ as we go) '
            write (6,'(A)') ' and then compute the innerproducts for '
            write (6,'(A)') ' these new templates. '
            write (6,'(A)') ' '
            write (6,'(A)') ' After calculating these new '
            write (6,'(A)') ' innerproducts (for the new templates in '
            write (6,'(A)') ' our now-larger list) we go back to '
            write (6,'(A)') ' Step 1, and once again find the '
            write (6,'(A)') ' templates that correspond to the largest '
            write (6,'(A)') ' innerproducts. '
            write (6,'(A)') ' '
            write (6,'(A)') ' After passing through our slowly growing '
            write (6,'(A)') ' list multiple times, we will eventually '
            write (6,'(A)') ' reach a point where each of the nearby '
            write (6,'(A)') ' templates was already considered. '
            write (6,'(A)') ' (i.e., after Step-2 the list does not '
            write (6,'(A)') ' grow). At this point we consider our '
            write (6,'(A)') ' local-search complete.' 
            write (6,'(A)') ' '
            write (6,'(A)') ' Step 4: '
            write (6,'(A)') ' Now we run through our list of '
            write (6,'(A)') ' innerproducts. The largest few will be '
            write (6,'(A)') ' used to update the alpha1d_ parameters. '
            write (6,'(A)') ' '
         end if !if (0+verbose.gt.1) then
               
         n_pass = 0
         continue_flag = .true.

         do while (continue_flag)

            n_LT_srt = 0
            C_Z_optimal = (1.0d0,0.0d0)
            C_S_optimal = (0.0d0,0.0d0)
            do ns=0,n_S-1
               if (S_use_(ns).eqv..true.) then
                  C_Z_optimal = C_Z_SM_(ns)
                  C_S_optimal = C_S_SM_(ns)
                  C_Z_srt_(n_LT_srt) = C_Z_optimal
                  if (real(C_Z_optimal).lt.real(C_Z_min)) then
                     C_Z_min = C_Z_optimal
                  end if !if (real(C_Z_optimal).lt.real(C_Z_min)) then
                  if (real(C_S_optimal).lt.real(C_S_min)) then
                     C_S_min = C_S_optimal
                  end if !if (real(C_S_optimal).lt.real(C_S_min)) then
                  LT_srt_(n_LT_srt) = ns
                  n_LT_srt = n_LT_srt + 1
               end if           !if (S_use_(ns).eqv..true.) then
            enddo
            
            if (0+verbose.gt.1) then
               write(6,'(A,I0)') 'n_pass: ' , n_pass
               write(6,'(A,I0)') 'found n_LT_srt: ' , n_LT_srt
               write(format_string,'(A,I0,A)') '(A,' , 2*n_LT_srt ,
     $              '(F8.4,1X))'
               write(6,format_string) 'pre: C_Z_srt_: ' , (C_Z_srt_(nLT)
     $              ,nLT=0,n_LT_srt-1)
               write(format_string,'(A,I0,A)') '(A,' , n_LT_srt ,
     $              '(I0,1X))'
               write(6,format_string) 'pre: LT_srt_: ' , (LT_srt_(nLT)
     $              ,nLT=0,n_LT_srt-1)
            end if              !if (0+verbose.gt.1) then
            call quicksort_c16(0,n_LT_srt-1,C_Z_srt_,1,LT_srt_,1
     $           ,quicksort_c16)
            if (0+verbose.gt.1) then
               write(format_string,'(A,I0,A)') '(A,' , 2*n_LT_srt ,
     $              '(F8.4,1X))'
               write(6,format_string) 'pos: C_Z_srt_: ' , (C_Z_srt_(nLT)
     $              ,nLT=0,n_LT_srt-1)
               write(format_string,'(A,I0,A)') '(A,' , n_LT_srt ,
     $              '(I0,1X))'
               write(6,format_string) 'pos: LT_srt_: ' , (LT_srt_(nLT)
     $              ,nLT=0,n_LT_srt-1)
            end if              !if (0+verbose.gt.1) then

            n_LT = 0
            if (0+verbose.gt.1) then
               write(6,'(A,I0)') ' using last entries: ' , min(n_LT_srt
     $              ,n_LT_hlf)
            end if              ! if (0+verbose.gt.1) then
            do nLT=0,min(n_LT_srt,n_LT_hlf)-1
               ns_use = LT_srt_(n_LT_srt-1-nLT)
               call get_angle_to_vp_(S_alpha_polar_a_(ns_use)
     $              ,S_alpha_azimu_b_(ns_use),vp_input_)
               if (0+verbose.gt.1) then
                  write(6,'(A,I0,A,F8.4,F8.4)') ' accessing ns_use ' ,
     $                 ns_use , ' with polar_a and azimu_b: ' ,
     $                 S_alpha_polar_a_(ns_use) ,
     $                 S_alpha_azimu_b_(ns_use)
               end if           !if (0+verbose.gt.1) then
               call tesselation_neighborhood_wrapper_0(n_S,L_ ,nl_max
     $              ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_
     $              ,T_c1_ ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base
     $              ,T_root_base_ ,vp_input_
     $              ,tesselation_distance_req_use,n_LT_tmp ,LT_(n_LT))
               n_LT = n_LT + n_LT_tmp
               if (0+verbose.gt.1) then
                  write(6,'(A,I0,A,I0)') ' pre: n_LT_tmp ' , n_LT_tmp ,
     $                 ' n_LT ', n_LT
               end if           !if (0+verbose.gt.1) then
               call unique_i4(n_LT,LT_,n_LT,LT_)
               if (0+verbose.gt.1) then
                  write(6,'(A,I0)') ' pos: n_LT ' , n_LT
               end if           !if (0+verbose.gt.1) then
            enddo               !do nLT=0,max(n_LT_srt-1,n_LT_hlf-1)

            n_LT_use=0
            do nLT=0,n_LT-1
               ns_use = LT_(nLT)
               if (S_use_(ns_use).eqv..true.) then
                  if (0+verbose.gt.1) then
                     write(6,'(A,I0,A,I0)') ' nLT ' , nLT ,
     $                    ', skipping ns_use ' , ns_use
                  end if        !if (0+verbose.gt.1) then
               end if !if (S_use_(ns_use).eqv..true.) then
               if (S_use_(ns_use).eqv..false.) then

                  if (0+verbose.gt.1) then
                     write(6,'(A,I0,A,I0,A,I0)') ' nLT ' , nLT ,
     $                    ' Calculating innerproduct for image ' , nm ,
     $                    ' and template ' , ns_use
                  end if        !if (0+verbose.gt.1) then
                  call test_innerproduct_batch_SM_stage_2a(verbose ,nm
     $                 ,ns_use ,n_r,grid_p_,n_w_,n_S,ld_S,n_A,nctf
     $                 ,n_w_max ,S_p_,S_p,S_q,M_p,M_q,Z_q_ ,C_M,C_S_
     $                 ,C_S_optimal,C_Z_optimal ,svd_calculation_type
     $                 ,n_svd_r ,svd_r_,n_svd_l ,svd_l_,svd_s_,svd_V_r_
     $                 ,Z_svdd_ ,Z_svdr_,Z_tmpC_ ,displacement_max
     $                 ,delta_x_est ,delta_y_est,gamma_z_est ,n_delta_x
     $                 ,n_delta_y,n_gamma_z ,delta_x_,delta_y_ ,gamma_z_
     $                 ,delta_x_est_ ,delta_y_est_ ,gamma_z_est_
     $                 ,ndx_optimal ,ndy_optimal ,ngz_optimal
     $                 ,delta_x_optimal ,delta_y_optimal
     $                 ,gamma_z_optimal ,fftw_plan_frwd_
     $                 ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_
     $                 ,timing_tic_1,timing_toc_1 ,timing_tot_1
     $                 ,timing_tic_2,timing_toc_2 ,timing_tot_2
     $                 ,timing_tic_3,timing_toc_3 ,timing_tot_3
     $                 ,timing_tic_4,timing_toc_4 ,timing_tot_4)
                  delta_x_SM_(ns_use) = delta_x_optimal
                  delta_y_SM_(ns_use) = delta_y_optimal
                  gamma_z_SM_(ns_use) = gamma_z_optimal
                  C_S_SM_(ns_use) = C_S_optimal
                  C_Z_SM_(ns_use) = C_Z_optimal
                  S_use_(ns_use) = .true.
                  n_S_use = n_S_use+1
                  n_LT_use = n_LT_use + 1
               end if           ! if (S_use_(ns_use).eqv..true.) then
            enddo               !do nLT=0,n_LT-1

            if (0+verbose.gt.1) then
               write(6,'(A,I0)') ' n_LT_use ' , n_LT_use
               write(6,'(A,I0)') ' n_S_use ' , n_S_use
            end if              !if (0+verbose.gt.1) then
            if (n_LT_use.gt.0) then
               if (0+verbose.gt.1) then
                  write(6,'(A)') ' Our list grew: do another pass. '
               end if           !if (0+verbose.gt.1) then
               continue_flag = .true.
            else
               if (0+verbose.gt.1) then
                  write(6,'(A)') ' Our list did not grow: stop. '
               end if           !if (0+verbose.gt.1) then
               continue_flag = .false.
            end if !if (n_LT_use.gt.0) then

            n_pass = n_pass + 1
         enddo !continue_flag

 5       continue

         if (0+verbose.gt.1) then
            write(6,'(A,I0)') ' finished nm ' , nm
            write(6,'(A,I0)') ' n_S_use ' , n_S_use
            write(6,'(A,2F8.4)') ' C_Z_min ' , C_Z_min
            write(6,'(A,2F8.4)') ' C_S_min ' , C_S_min
         end if                 !if (0+verbose.gt.1) then
         n_SM_use = n_SM_use + n_S_use

         if (verbose.gt.1) then
            write(6,'(A)') ' Fill in uncalculated innerproducts.'
         end if !if (verbose.gt.1) then
         do ns=0,n_S-1
            if (S_use_(ns).eqv..false.) then
               C_Z_SM_(ns) = C_Z_min
               C_S_SM_(ns) = C_S_min
               delta_x_SM_(ns) = 0.0d0
               delta_y_SM_(ns) = 0.0d0
               gamma_z_SM_(ns) = 0.0d0
            end if ! if (S_use_(ns).eqv..false.) then
         enddo !do ns=0,n_S-1

         if (0+verbose.gt.1) then
            write(6,'(A)') ' Sort innerproducts.'
         end if !if (verbose.gt.1) then
         if (0+verbose.gt.1) then
            write(6,'(A)') 'delta_x_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_x_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'delta_y_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_y_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'gamma_z_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (gamma_z_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_S_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_S_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_Z_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_Z_SM_(ns),ns=0,n_S-1)
         end if !if (verbose.gt.1) then
         call test_innerproduct_batch_sort_a(n_S,1,delta_x_SM_
     $        ,delta_y_SM_,gamma_z_SM_,C_S_SM_,C_Z_SM_
     $        ,delta_x_sort_SM_ ,delta_y_sort_SM_,gamma_z_sort_SM_
     $        ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_ ,I_inverse_SM_)
         if (0+verbose.gt.1) then
            write(6,'(A)') 'delta_x_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_x_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'delta_y_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_y_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'gamma_z_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (gamma_z_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_S_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_S_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_Z_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_Z_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'I_permute_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
            write(6,format_string) (I_permute_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'I_inverse_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
            write(6,format_string) (I_inverse_SM_(ns),ns=0,n_S-1)
         end if !if (verbose.gt.1) then

         if (verbose.gt.1) then
            write(6,'(A)') ' Update image parameters.'
         end if !if (verbose.gt.1) then
         nr_use = max(0,min(n_S-1,floor((1.0d0-rand()*alpha_update_f)
     $        *n_S)))
         ns_use = I_permute_SM_(nr_use)
         call get_interchange_delta(delta_x_SM_(ns_use)
     $        ,delta_y_SM_(ns_use),gamma_z_est_(nm),delta_x_tmp
     $        ,delta_y_tmp)
         polar_a_upd_(nm) = S_alpha_polar_a_(ns_use)
         azimu_b_upd_(nm) = S_alpha_azimu_b_(ns_use)
         gamma_z_upd_(nm) = gamma_z_est_(nm) + gamma_z_SM_(ns_use)
         delta_x_upd_(nm) = delta_x_est_(nm) + delta_x_tmp
         delta_y_upd_(nm) = delta_y_est_(nm) + delta_y_tmp
         S_index_upd_(nm) = ns_use
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,2F8.3,A,2F8.3,A,2F8.3)') 'nm ' ,nm
     $           ,'; ns_use ',ns_use,'; C_S ' ,C_S_SM_(ns_use),'; C_M '
     $           ,C_M_(nm) ,'; C_Z ', C_Z_SM_(ns_use)
         end if !if (verbose.gt.1) then
         if (zabs(C_S_SM_(ns_use)).le.1.0d-15) then
            l2_norm_upd_(nm) = 1.0d0
         else 
            l2_norm_upd_(nm) = zabs(C_M_(nm))
     $           *zabs(C_Z_SM_(ns_use)) /
     $           zabs(C_S_SM_(ns_use))
         end if ! if (zabs(C_S_SM_(ns_use)).le.1.0d-15) then

      enddo !do nm=0,n_M-1

      if (0+verbose.gt.0) then
         write(6,'(A,I0,A,I0)')
     $        ' Total n_SM_use/(n_S templates*n_M images): ' ,n_SM_use ,
     $        '/' , n_S*n_M
         write(6,'(A,F16.8,A,F16.8)')
     $        ' transforming M_p_ to M_q_: total time ' ,timing_tot_0
     $        ,'; time_per ',timing_tot_0/max(1,n_M)
         write(6,'(A,F16.8,A,F16.8)')
     $        ' transforming S_p_ to S_p_: total time ' ,timing_tot_1
     $        ,'; time_per ',timing_tot_1/max(1,min(n_SM_use,n_S *n_M))
         write(6,'(A,F16.8,A,F16.8)')
     $        '    innerproduct_q_k_svdr: total time ' ,timing_tot_2
     $        ,'; time_per ',timing_tot_2/max(1,min(n_SM_use,n_S *n_M))
         write(6,'(A,I0,A,F16.8,A,F16.8)') '    svd_calculation_type_'
     $        ,svd_calculation_type,': total time ' ,timing_tot_3
     $        ,'; time_per ',timing_tot_3/max(1,min(n_SM_use,n_S *n_M))
         write(6,'(A,F16.8,A,F16.8)')
     $        ' excerpt_5 + excerpt_6: total time ' ,timing_tot_4
     $        ,'; time_per ',timing_tot_4/max(1,min(n_SM_use,n_S *n_M))
      end if !verbose

      if (0+verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
      deallocate(C_Z_sort_SM_)
      deallocate(C_S_sort_SM_)
      deallocate(delta_x_sort_SM_)
      deallocate(delta_y_sort_SM_)
      deallocate(gamma_z_sort_SM_)
      deallocate(C_Z_SM_)
      deallocate(C_S_SM_)
      deallocate(delta_x_SM_)
      deallocate(delta_y_SM_)
      deallocate(gamma_z_SM_)
      deallocate(vp_input_)
      deallocate(LT_)
      deallocate(LT_srt_)
      deallocate(C_Z_srt_)
      deallocate(S_use_)
      deallocate(Z_q_)
      deallocate(M_q)
      deallocate(S_p)
      deallocate(S_q)
      deallocate(M_p)
      deallocate(Z_tmpC_)
      deallocate(Z_svdr_)
c$$$      Z_svdd_ should be passed in as input
c$$$      deallocate(Z_svdd_)
      deallocate(gamma_z_)
      deallocate(delta_y_)
      deallocate(delta_x_)
c$$$      fftw plans and workspace should be passed in as input
c$$$      if (0+verbose.gt.1) then
c$$$         write(6,'(A)') ' Destroying fftw_plans.'
c$$$      end if
c$$$      do nr=0,n_r-1
c$$$         call dfftw_destroy_plan_(fftw_plan_back_(nr))
c$$$         call dfftw_destroy_plan_(fftw_plan_frwd_(nr))
c$$$      enddo
c$$$      deallocate(fftw_out_)
c$$$      deallocate(fftw_in1_)
c$$$      deallocate(fftw_plan_back_)
c$$$      deallocate(fftw_plan_frwd_)

c$$$      svd arrays should be passed in as input
c$$$      if (0+verbose.gt.1) then
c$$$         write(6,'(A)') ' Deallocating svd arrays.'
c$$$      end if
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (0+verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_SM_stage_1a]'
      end if

 10   continue
      end


