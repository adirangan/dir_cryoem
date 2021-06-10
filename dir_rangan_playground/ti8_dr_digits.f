c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program ti8_dr_digits
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      integer nMode ! integer: Mode of operation for the driver. Mode -1 corresponds to the assortment of hardcoded internal parameters below. Nonnegative modes (from 0 to ..) test different combinations of internal parameters. ;
      integer n_Mode ! integer: total number of Modes. ;
      integer tmp_nMode_r, tmp_nMode_q ! temporary: integer: used to factor nMode. ;
      integer n_omp_sub_0in
      integer n_S_0_sub_0in
      integer n_S_1_sub_0in
      integer n_M_0_sub_0in
      integer n_M_1_sub_0in
      logical flag_ctf
      logical flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
c$$$      indices
      integer *4 nr,na,nw
c$$$      grids
      real *8 pi
      real *8 half_diameter_x_c ! real *8: half_diameter of box in x_c coordinates. ;
      real *8 diameter_k_c ! real *8: diameter of image in k_c coordinates. ;
      real *8 half_diameter_k_c ! real *8: half_diameter of image in k_c coordinates. ;
      integer quad_n_r ! integer: number of quadrature nodes in the radial direction. ;
      integer quad_n_w_sum ! integer: number of points in quad_k_p. ;
      integer quad_n_w_max ! integer: number of points at quad_grid_k_p_r_(quad_n_r-1). ;
      real *8, allocatable :: jacpts_x_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: jacpts_w_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: quad_grid_k_p_r_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: quad_weight_k_p_r_(:) !real *8 array, size quad_n_r. ;
      integer *4, allocatable :: quad_n_polar_a_(:) !integer *4 array, size quad_n_r. ;
      integer *4, allocatable :: quad_n_w_(:) !integer *4 array, size quad_n_r. ;
      real *8, allocatable :: quad_grid_k_p_0_(:) !real *8 array, size quad_n_w_sum. ;
      real *8, allocatable :: quad_grid_k_p_1_(:) !real *8 array, size quad_n_w_sum. ;
      real *8, allocatable :: quad_weight_k_p_(:) !real *8 array, size quad_n_w_sum. ;
c$$$      Templates S_
      integer n_S_max,n_S,ld_S,ns,nS_sample
      complex *16, allocatable :: S_k_p__(:)
      integer *4, allocatable :: I_S_sample_(:)
c$$$      Images M_
      integer n_M_max,n_M,ld_M,nm,nM_sample
      complex *16, allocatable :: M_k_p__(:)
      integer *4, allocatable :: I_M_sample_(:)
c$$$      CTF-functions
      integer *4 n_ctf,ld_CTF,nctf
      integer *4, allocatable :: I_ctf_(:)
      complex *16, allocatable :: CTF_k_p__(:)
c$$$      fftw_plan_many 
      integer fpm_howmany_max
      real *8 alpha_update_f
      parameter (alpha_update_f=0.0d0)
c$$$      image parameters 
      include 'excerpt_define_nalpha.f'
      real *8, allocatable :: alpha_est__(:)
      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: S_alpha_S_index_(:)
      real *8, allocatable :: S_alpha_polar_a_(:)
      real *8, allocatable :: S_alpha_azimu_b_(:)
      real *8, allocatable :: fix_phi_S_(:)
      real *8, allocatable :: fix_delta_S_x_c_0_(:)
      real *8, allocatable :: fix_delta_S_x_c_1_(:)
      real *8, allocatable :: fix_delta_S_x_p_r_(:)
      real *8, allocatable :: fix_delta_S_x_p_w_(:)
      real *8, allocatable :: fix_phi_CTF_(:)
      real *8, allocatable :: fix_delta_M_x_c_0_(:)
      real *8, allocatable :: fix_delta_M_x_c_1_(:)
      real *8, allocatable :: fix_delta_M_x_p_r_(:)
      real *8, allocatable :: fix_delta_M_x_p_w_(:)
      real *8 fix_phi_S
      real *8 fix_delta_S_x_c_0, fix_delta_S_x_c_1
      real *8 fix_delta_S_x_p_r, fix_delta_S_x_p_w, fix_nwav_S
      real *8 fix_phi_CTF
      real *8 fix_delta_M_x_c_0, fix_delta_M_x_c_1
      real *8 fix_delta_M_x_p_r, fix_delta_M_x_p_w, fix_nwav_M
c$$$      estimated (current) displacement, rotation and scaling parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: S_index_est_(:)
c$$$      range of displacement and rotation parameters to measure
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 displacement_max
      real *8 tmp_k,tmp_w,tmp_e,tmp_f
      integer n_delta_x,n_delta_y,n_delta_v,ndx,ndy,ndv
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_upd,delta_y_upd,gamma_z_upd,nwav_upd
      real *8 ctf_ind_upd,S_index_upd,M_index_upd
      real *8 l2_norm_upd,CTF_R_S_upd,C_Z_opt_upd
      real *8 SRTRTCM_upd
      real *8 delta_x_est,delta_y_est,gamma_z_est,nwav_est
      real *8 ctf_ind_est,S_index_est,M_index_est
      real *8 l2_norm_est,CTF_R_S_est,C_Z_opt_est
      real *8 form_SRTRTCM,error_SRTRTCM
      real *8, allocatable :: ferror_SRTRTCM_(:)
      complex *16 form_C_Z_opt,error_C_Z_opt
      complex *16, allocatable :: ferror_C_Z_opt_(:)
      complex *16 C_M_opt_upd,form_C_M_opt,error_C_M_opt
      complex *16, allocatable :: ferror_C_M_opt_(:)
      real *8 form_CTF_R_S,error_CTF_R_S
      real *8, allocatable :: ferror_CTF_R_S_(:)
      real *8 form_l2_norm,error_l2_norm
      real *8, allocatable :: ferror_l2_norm_(:)
      real *8 eps_svd,N_pixels_in,delta_max
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      integer n_gamma_z,ngz
      real *8, allocatable :: gamma_z_(:)      
      integer svd_calculation_type
c$$$      tesselation parameters
      real *8 tesselation_distance_req !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_rq determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere).
      integer *4 n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      integer *4 n_LT_ref ! number of image-template pairs to consider when refining local search.
c$$$      SM storage
      integer *4 n_SM_max ! total (maximum) number of templates to store per image. ;
      integer *4, allocatable :: n_SM_(:) ! array of size n_M indicating the actual number of templates stored per image. ;
      real *8, allocatable :: alpha_SM__(:) ! array of size n_alpha*n_SM_max*n_M storing the image-parameters for each stored template-image pair. ;
      integer *4 n_MS_max ! total (maximum) number of images to store per template. ;
      integer *4, allocatable :: n_MS_(:) ! array of size n_S indicating the actual number of images stored per template. ;
      real *8, allocatable :: alpha_MS__(:) ! array of size n_alpha*n_MS_max*n_S storing the image-parameters for each stored image-template pair. ;
      character(len=1024) fname,format_string,prefix_string
      real *8 max_r8_f,avg_r8_f,al2_r8_f !function output. ;
      real *8 maxr_c16_f,fnod_c16_f,al2_c16_f,sum_r8_f !function output. ;
      complex *16 dotw_c16_f !function output. ;
      real *8 I_P_f,I_xP_f,I_PxP_f,I_PxxP_f,l2_xRPx_f !function output. ;
      real *8 I_PxRTRTxP_f !function output. ;
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      random seed
      integer *4 rseed !random seed. ;
c$$$      
      logical flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
      logical flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
      real *8 d_memory_estimate !real *8: estimate of required memory (in bytes). ;
      logical flag_memory_checkset !temporary: logical flag used in ti8_dr_digits_excerpt_checkset.f ;
c$$$      command line input variables
      character(len=64) :: cmd_argstring
      integer ncmd
      real *8 cmd_Mode
c$$$      
      rseed = 1
      pi = 4.0d0*datan(1.0d0)
c$$$
      ncmd = 0
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_Mode
      ncmd = ncmd + 1
      nMode = nint(cmd_Mode)
c$$$      
c$$$      List of parameters: ;
      verbose = 1
      eps_svd = 0.000001d0
      N_pixels_in = 2.5d0*dsqrt(2.0d0) ! Note that this should match one of the svd-libraries. ;
      N_pixels_in = 4.5d0 ! e.g., we have a library for N_pixels = 4.5d0. ;
      n_delta_x = 15
      n_delta_y = 15
c$$$  /*%%%%%%%%*/
c$$$      n_delta_x = 1
c$$$      n_delta_y = 1
c$$$  /*%%%%%%%%*/
      svd_calculation_type = 1
      flag_RTRT_vs_RTTR = .false.
      n_omp_sub_0in = 8
      n_S_0_sub_0in = 2
      n_S_1_sub_0in = 3
      n_M_0_sub_0in = 5
      n_M_1_sub_0in = 2
c$$$ /*%%%%%%%%*/
c$$$      n_omp_sub_0in = 1
c$$$      n_S_0_sub_0in = 1
c$$$      n_S_1_sub_0in = 1
c$$$      n_M_0_sub_0in = 1
c$$$      n_M_1_sub_0in = 1
c$$$ /*%%%%%%%%*/
      flag_ctf = .true.
      n_S_max = 12*3
      n_S = 5*3
      n_M_max = 17*3
      n_M = 15*3
      n_CTF = 10
c$$$  /*%%%%%%%%*/
      n_S = 1
      n_M = 1
      n_CTF = 1
c$$$  /*%%%%%%%%*/
      fpm_howmany_max = 16
      flag_MS_vs_SM = .true.
      n_SM_max = 45
      n_MS_max = 30
      flag_memory_estimate = .false.
      flag_time_Zstore = .false.
c$$$      tesselation_distance_req = 0.25
      tesselation_distance_req = 2.0d0
      n_LT_add = max(1,n_M/10) ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      n_LT_ref = max(1,n_SM_max/2) ! number of image-template pairs to consider when refining local search.

      n_Mode = 2*2*3
      nMode=max(-1,min(nMode,n_Mode-1))

      if (nMode.ge.0) then
         verbose = 0
         tmp_nMode_q = nMode
         tmp_nMode_r = mod(tmp_nMode_q,2)
         tmp_nMode_q = tmp_nMode_q - tmp_nMode_r
         tmp_nMode_q = tmp_nMode_q/2
         if (tmp_nMode_r.eq.0) then
            svd_calculation_type=1
         else
            svd_calculation_type=2
         end if !if (tmp_nMode_r.eq.0) then
         tmp_nMode_r = mod(tmp_nMode_q,2)
         tmp_nMode_q = tmp_nMode_q - tmp_nMode_r
         tmp_nMode_q = tmp_nMode_q/2
         if (tmp_nMode_r.eq.0) then
            flag_RTRT_vs_RTTR = .true.
         else
            flag_RTRT_vs_RTTR = .false.
         end if !if (tmp_nMode_r.eq.0) then
         tmp_nMode_r = mod(tmp_nMode_q,3)
         tmp_nMode_q = tmp_nMode_q - tmp_nMode_r
         tmp_nMode_q = tmp_nMode_q/3
         if (tmp_nMode_r.eq.0) then
            flag_MS_vs_SM = .true.
            tesselation_distance_req = 2.0d0
         end if !if (tmp_nMode_r.eq.0) then
         if (tmp_nMode_r.eq.1) then
            flag_MS_vs_SM = .false.
            tesselation_distance_req = 2.0d0
         end if !if (tmp_nMode_r.eq.1) then
         if (tmp_nMode_r.eq.2) then
            flag_MS_vs_SM = .false.
            tesselation_distance_req = 0.25d0
         end if !if (tmp_nMode_r.eq.2) then
         if (tmp_nMode_q.ne.0) then
            write(6,'(A,I0)') ' Warning, tmp_nMode_q: ' , tmp_nMode_q
         end if !if (tmp_nMode_q.ne.0) then
         if (verbose.gt.-1) then
            write(6,'(A,I0,A,I0,A,L2,A,L2,A,F6.3)')
     $           ' nMode '
     $           ,nMode
     $           ,' svd_calculation_type '
     $           ,svd_calculation_type
     $           ,' flag_RTRT_vs_RTTR '
     $           ,flag_RTRT_vs_RTTR
     $           ,' flag_MS_vs_SM '
     $           ,flag_MS_vs_SM
     $           ,' tesselation_distance_req '
     $           ,tesselation_distance_req
         end if !if (verbose.gt.-1) then
      end if !if (nmode.eq.0) then
      
      fpm_howmany_max = max(1,fpm_howmany_max)
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering ti8_dr_digits]: '
       end if
       if (verbose.gt.0) then
          write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
       end if !if (verbose.gt.0) then
       call openblas_set_num_threads(1)
       if (verbose.gt.0) then
         write(6,'(A,I0)') ' verbose: ',verbose
         write(6,'(A,F24.16)') ' eps_svd: ',eps_svd
         write(6,'(A,F24.16)') ' N_pixels_in: ',N_pixels_in
         write(6,'(A,I0)') ' n_delta_x: ',n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ',n_delta_y
         write(6,'(A,I0)') ' svd_calculation_type: '
     $        ,svd_calculation_type
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' ,flag_RTRT_vs_RTTR
         write(6,'(A,I0)') ' n_omp_sub_0in: ' ,n_omp_sub_0in
         write(6,'(A,I0)') ' n_S_0_sub_0in: ' ,n_S_0_sub_0in
         write(6,'(A,I0)') ' n_S_1_sub_0in: ' ,n_S_1_sub_0in
         write(6,'(A,I0)') ' n_M_0_sub_0in: ' ,n_M_0_sub_0in
         write(6,'(A,I0)') ' n_M_1_sub_0in: ' ,n_M_1_sub_0in
         write(6,'(A,I0)') ' n_S_max: ' ,n_S_max
         write(6,'(A,I0)') ' n_M_max: ' ,n_M_max
         write(6,'(A,I0)') ' fpm_howmany_max: ' ,fpm_howmany_max
         write(6,'(A,L1)') ' flag_ctf: ' ,flag_ctf
      end if
      
      half_diameter_x_c = 1.0d0
      diameter_k_c = 32.0d0
      half_diameter_k_c = diameter_k_c/2.0d0
      delta_max = N_pixels_in*dsqrt(2.0d0)/2.0d0/half_diameter_k_c
      if (verbose.gt.1) then
      write(6,'(A,F21.16)') ' delta_max: ' , delta_max
      end if !if (verbose.gt.1) then
c$$$      %%%%%%%%
      quad_n_r = 48
      allocate(jacpts_x_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,jacpts_x_)
      allocate(jacpts_w_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,jacpts_w_)
      call jacpts_48(quad_n_r,jacpts_x_,jacpts_w_)
      if (verbose.gt.1) then
      call print_sub_r8(quad_n_r,jacpts_x_,' jacpts_x_: ')
      call print_sub_r8(quad_n_r,jacpts_w_,' jacpts_w_: ')
      end if !if (verbose.gt.1) then
      allocate(quad_grid_k_p_r_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,quad_grid_k_p_r_)
      allocate(quad_weight_k_p_r_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,quad_weight_k_p_r_)
      call af1_r8(quad_n_r,0.5d0*half_diameter_k_c,0.5d0
     $     *half_diameter_k_c,jacpts_x_,quad_grid_k_p_r_)
      call af1_r8(quad_n_r,half_diameter_k_c**2 / sum_r8_f(quad_n_r
     $     ,jacpts_w_),0.0d0,jacpts_w_,quad_weight_k_p_r_)
      if (verbose.gt.1) then
      call print_sub_r8(quad_n_r,quad_grid_k_p_r_,' quad_grid_k_p_r_: ')
      call print_sub_r8(quad_n_r,quad_weight_k_p_r_
     $     ,' quad_weight_k_p_r_: ')
      end if !if (verbose.gt.1) then
c$$$      %%%%%%%%
      allocate(quad_n_polar_a_(0:1+quad_n_r-1))
      call cs1_i4(quad_n_r,quad_n_polar_a_)
      do nr=0,quad_n_r-1
         quad_n_polar_a_(nr) = ceiling(pi*(quad_grid_k_p_r_(nr)+1))
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.1) then
      call print_sub_i4(quad_n_r,quad_n_polar_a_,' quad_n_polar_a_: ')
      end if !if (verbose.gt.1) then
      allocate(quad_n_w_(0:1+quad_n_r-1))
      call cs1_i4(quad_n_r,quad_n_w_)
      quad_n_w_sum = 0
      do nr=0,quad_n_r-1
         quad_n_w_(nr) = 2*quad_n_polar_a_(nr)
         quad_n_w_sum = quad_n_w_sum + quad_n_w_(nr)
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.1) then
      call print_sub_i4(quad_n_r,quad_n_w_,' quad_n_w_: ')
      end if !if (verbose.gt.1) then
      quad_n_w_max = quad_n_w_(quad_n_r-1)
      allocate(quad_grid_k_p_0_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_grid_k_p_0_)
      allocate(quad_grid_k_p_1_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_grid_k_p_1_)
      allocate(quad_weight_k_p_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_weight_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = 2.0d0*pi*nw/quad_n_w_(nr)
            quad_grid_k_p_0_(na) = pi*tmp_k/half_diameter_k_c
     $           *dcos(tmp_w)
            quad_grid_k_p_1_(na) = pi*tmp_k/half_diameter_k_c
     $           *dsin(tmp_w)
            quad_weight_k_p_(na) = quad_weight_k_p_r_(nr) * pi /
     $           quad_n_w_(nr)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.1) then
         call print_sub_r8(quad_n_w_sum,quad_grid_k_p_0_
     $        ,' quad_grid_k_p_0_: ')
         call print_sub_r8(quad_n_w_sum,quad_grid_k_p_1_
     $        ,' quad_grid_k_p_1_: ')
         call print_sub_r8(quad_n_w_sum,quad_weight_k_p_
     $        ,' quad_weight_k_p_: ')
      write(6,'(A,F24.16)') ' sum(quad_weight_k_p_): ' ,
     $     sum_r8_f(quad_n_w_sum,quad_weight_k_p_)
      end if !if (verbose.gt.1) then
c$$$      %%%%%%%%
      
c$$$      indices
      ld_S = quad_n_w_sum + 3 ! Note: extending length arbitrarily. ;
      ld_M = quad_n_w_sum + 5 ! Note: extending length arbitrarily. ;
      ld_CTF = quad_n_w_sum + 7 ! Note: extending length arbitrarily. ;

c$$$      defining displacement_max to be 3.0d0 wavelengths (can change later)
      displacement_max = 3.0d0/half_diameter_k_c

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
      end if
      allocate(delta_x_(0:1+ceiling(pi*n_delta_x*n_delta_y)-1))
      call cs1_r8(ceiling(pi*n_delta_x*n_delta_y),delta_x_)
      allocate(delta_y_(0:1+ceiling(pi*n_delta_x*n_delta_y)-1))
      call cs1_r8(ceiling(pi*n_delta_x*n_delta_y),delta_y_)
      call get_delta_1(
     $     N_pixels_in
     $     ,quad_n_r
     $     ,half_diameter_x_c
     $     ,n_delta_x
     $     ,n_delta_y
     $     ,n_delta_v
     $     ,delta_x_
     $     ,delta_y_
     $     )
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)')
     $        ' quad_n_w_max or so.'
      end if
      n_gamma_z = quad_n_w_max
      allocate(gamma_z_(0:1+n_gamma_z-1))
      call cs1_r8(n_gamma_z,gamma_z_)
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.0) then
         call print_sub_r8(n_delta_v,delta_x_,' delta_x_: ')
         call print_sub_r8(n_delta_v,delta_y_,' delta_y_: ')
         call print_sub_r8(n_gamma_z,gamma_z_,' gamma_z_: ')
      end if !if (verbose.gt.1) then

c$$$      generate templates, images and CTFs
      allocate(I_S_sample_(0:1+n_S-1))
      call cs1_i4(n_S,I_S_sample_)
      do ns=0,n_S-1
         I_S_sample_(ns)=max(0,min(n_S_max-1,nint(n_S_max*(1.0d0*ns)
     $        /(1.0d0*n_S))))
      enddo !do ns=0,n_S-1
      if (verbose.gt.1) then
      call print_sub_i4(n_S,I_S_sample_,' I_S_sample_: ')
      end if !if (verbose.gt.1) then
      allocate(S_alpha_S_index_(0:1+n_S-1))
      call cs1_r8(n_S,S_alpha_S_index_)
      allocate(S_alpha_polar_a_(0:1+n_S-1))
      call cs1_r8(n_S,S_alpha_polar_a_)
      allocate(S_alpha_azimu_b_(0:1+n_S-1))
      call cs1_r8(n_S,S_alpha_azimu_b_)
      do ns=0,n_S-1
         S_alpha_S_index_(ns) = ns ! array storing the numbers 0:n_S-1, from 0 up to the number of templates used - 1. ;
         S_alpha_polar_a_(ns) = pi*(1.0d0*ns + 0.5d0)/(1.0d0*n_S)
         S_alpha_azimu_b_(ns) = 2.0d0*pi*(1.0d0*ns + 0.5d0)/(1.0d0*n_S)
      enddo
      allocate(I_M_sample_(0:1+n_M-1))
      call cs1_i4(n_M,I_M_sample_)
      do nm=0,n_M-1
         I_M_sample_(nm)=max(0,min(n_M_max-1,nint(n_M_max*(1.0d0*nm)
     $        /(1.0d0*n_M))))
      enddo !do nm=0,n_M-1
      if (verbose.gt.1) then
      call print_sub_i4(n_M,I_M_sample_,' I_M_sample_: ')
      end if !if (verbose.gt.1) then
      if (flag_memory_estimate.eqv..false.) then
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') ' Here we construct the templates "S_",'
         write(6,'(A)') ' images "M_", and CTFs "CTF_".'
         write(6,'(A,I0,A)') ' We generate n_S = ',n_S
     $        ,' templates in total,'
         write(6,*) 'each using a plane-wave multiplied'
     $        ,' by a planar-function.'
         write(6,'(A,I0,A)') ' We generate n_M = ',n_M
     $        ,' templates in total,'
         write(6,*) 'each using a plane-wave.'
         write(6,'(A,I0,A)') ' We generate n_CTF = ',n_CTF
     $        ,' CTFs in total,'
         write(6,*) 'each using a planar-function.'
         write(6,*) 'For each of these, '
     $        ,' The array ?_k_p__ holds the k-space values on a'
     $        ,' quasi-uniform polar-grid defined by'
     $        ,' half_diameter_k_c and quad_n_w_.'
      end if !if (verbose.gt.1) then
      allocate(S_k_p__(0:1+ld_S*n_S_max-1))
      call cs1_c16(ld_S*n_S_max,S_k_p__)
      allocate(M_k_p__(0:1+ld_M*n_M_max-1))
      call cs1_c16(ld_M*n_M_max,M_k_p__)
      allocate(CTF_k_p__(0:1+ld_CTF*n_CTF-1))
      call cs1_c16(ld_CTF*n_CTF,CTF_k_p__)
      allocate(fix_phi_S_(0:1+n_S_max-1))
      call cs1_r8(n_S_max,fix_phi_S_)
      allocate(fix_delta_S_x_c_0_(0:1+n_S_max-1))
      call cs1_r8(n_S_max,fix_delta_S_x_c_0_)
      allocate(fix_delta_S_x_c_1_(0:1+n_S_max-1))
      call cs1_r8(n_S_max,fix_delta_S_x_c_1_)
      allocate(fix_delta_S_x_p_r_(0:1+n_S_max-1))
      call cs1_r8(n_S_max,fix_delta_S_x_p_r_)
      allocate(fix_delta_S_x_p_w_(0:1+n_S_max-1))
      call cs1_r8(n_S_max,fix_delta_S_x_p_w_)
      allocate(fix_phi_CTF_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,fix_phi_CTF_)
      allocate(fix_delta_M_x_c_0_(0:1+n_M_max-1))
      call cs1_r8(n_M_max,fix_delta_M_x_c_0_)
      allocate(fix_delta_M_x_c_1_(0:1+n_M_max-1))
      call cs1_r8(n_M_max,fix_delta_M_x_c_1_)
      allocate(fix_delta_M_x_p_r_(0:1+n_M_max-1))
      call cs1_r8(n_M_max,fix_delta_M_x_p_r_)
      allocate(fix_delta_M_x_p_w_(0:1+n_M_max-1))
      call cs1_r8(n_M_max,fix_delta_M_x_p_w_)
      do ns=0,n_S_max-1
         tmp_f = 1.0d0*((2.0d0*ns)/(1.0d0*n_S_max)-1.0d0)
         fix_phi_S_(ns) = 3.5d0*pi/12.0d0 + pi*tmp_f
         fix_delta_S_x_c_0_(ns) = +4.0d0*0.75d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         fix_delta_S_x_c_1_(ns) = +4.0d0*0.25d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         fix_delta_S_x_p_r_(ns) = dsqrt(fix_delta_S_x_c_0_(ns)**2 +
     $        fix_delta_S_x_c_1_(ns)**2)
         fix_delta_S_x_p_w_(ns) = datan2(fix_delta_S_x_c_1_(ns)
     $        ,fix_delta_S_x_c_0_(ns))
      enddo !do ns=0,n_S_max-1
      if (verbose.gt.1) then
         call print_all_r8(n_S_max,fix_phi_S_, ' fix_phi_S_: ')
         call print_all_r8(n_S_max,fix_delta_S_x_c_0_,
     $        ' fix_delta_S_x_c_0_: ')
         call print_all_r8(n_S_max,fix_delta_S_x_c_1_,
     $        ' fix_delta_S_x_c_1_: ')
         call print_all_r8(n_S_max,fix_delta_S_x_p_r_,
     $        ' fix_delta_S_x_p_r_: ')
         call print_all_r8(n_S_max,fix_delta_S_x_p_w_,
     $        ' fix_delta_S_x_p_w_: ')
      end if !if (verbose.gt.1) then
      do nm=0,n_M_max-1
         tmp_f = 1.0d0*((2.0d0*nm)/(1.0d0*n_M_max)-1.0d0)
         fix_delta_M_x_c_0_(nm) = -4.0d0*0.65d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         fix_delta_M_x_c_1_(nm) = +4.0d0*0.65d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         fix_delta_M_x_p_r_(nm) = dsqrt(fix_delta_M_x_c_0_(nm)**2 +
     $        fix_delta_M_x_c_1_(nm)**2)
         fix_delta_M_x_p_w_(nm) = datan2(fix_delta_M_x_c_1_(nm)
     $        ,fix_delta_M_x_c_0_(nm))
      enddo !do nm=0,n_M_max-1
      do nctf=0,n_CTF-1
         tmp_f = 1.0d0*((2.0d0*nctf)/(1.0d0*n_CTF)-1.0d0)
         fix_phi_CTF_(nctf) = 8.2d0*pi/12.0d0 + pi*tmp_f
      enddo !do nctf=0,n_CTF-1
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            do ns=0,n_S_max-1
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_S_x_p_r_(ns)*dcos(tmp_w
     $           -fix_delta_S_x_p_w_(ns))
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_S_(ns))
            S_k_p__(na+ns*ld_S) = dcmplx(tmp_f,0.0d0) * dcmplx(
     $           +dcos(tmp_e),+dsin(tmp_e))
            enddo !do ns=0,n_S_max-1
            do nm=0,n_M_max-1
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_M_x_p_r_(nm)*dcos(tmp_w
     $           -fix_delta_M_x_p_w_(nm))
            M_k_p__(na+nm*ld_M) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            enddo !do nm=0,n_M_max-1
            do nctf=0,n_CTF-1
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_CTF_(nctf))
            CTF_k_p__(na+nctf*ld_CTF) = dcmplx(tmp_f,0.0d0)
            enddo !do nctf=0,n_CTF-1
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1      
      end if !if (flag_memory_estimate.eqv..false.) then

      if (flag_memory_estimate.eqv..false.) then
      allocate(polar_a_est_(0:1+n_M-1))
      call cs1_r8(n_M,polar_a_est_)
      allocate(azimu_b_est_(0:1+n_M-1))
      call cs1_r8(n_M,azimu_b_est_)
      allocate(delta_x_est_(0:1+n_M-1))
      call cs1_r8(n_M,delta_x_est_)
      allocate(delta_y_est_(0:1+n_M-1))
      call cs1_r8(n_M,delta_y_est_)
      allocate(gamma_z_est_(0:1+n_M-1))
      call cs1_r8(n_M,gamma_z_est_)
      allocate(l2_norm_est_(0:1+n_M-1))
      call cs1_r8(n_M,l2_norm_est_)
      allocate(S_index_est_(0:1+n_M-1))
      call cs1_r8(n_M,S_index_est_)
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we set up the estimated translations'
     $        ,' and rotations for each image.'
      end if !if (verbose.gt.1) then      
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'For this test we choose these estimated'
     $        ,' translations and rotations to be discrete shifts'
     $        ,' away from the true translations and rotations.'
      end if
      do nm=0,n_M-1
         tmp_f = ((2.0d0*nm)/(1.0d0*n_M_max)-1.0d0)
         polar_a_est_(nm) = pi/2.0d0 + pi/2.0d0*tmp_f
         azimu_b_est_(nm) = pi*tmp_f
         l2_norm_est_(nm) = 1.0d0
         S_index_est_(nm) = 0.0d0
         delta_x_est_(nm) = +1.0d0*0.35d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         delta_y_est_(nm) = -1.0d0*0.95d0/half_diameter_k_c +
     $        tmp_f/half_diameter_k_c
         gamma_z_est_(nm) = 0.25d0*pi + pi*tmp_f
      enddo !do nm=0,n_M-1
      end if !if (flag_memory_estimate.eqv..false.) then

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         Determine which images correspond to which ctf-function.
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_memory_estimate.eqv..false.) then
      allocate(I_ctf_(0:1+n_M-1))
      call cs1_i4(n_M,I_ctf_)
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo
      end if !if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.0) then
         write(6,'(A)') ' Setting displacement_max to 2.0d0 '
      end if ! if (verbose.gt.0) then
      displacement_max = 2.0d0
      if (flag_memory_estimate.eqv..false.) then
      allocate(alpha_est__(0:1+n_alpha*n_M-1))
      call cs1_r8(n_alpha*n_M,alpha_est__)
      call cl1_r8(n_alpha*n_M,alpha_est__)
      do nm=0,n_M-1
         alpha_est__(nalpha_polar_a + n_alpha*nm) = polar_a_est_(nm)
         alpha_est__(nalpha_azimu_b + n_alpha*nm) = azimu_b_est_(nm)
         alpha_est__(nalpha_gamma_z + n_alpha*nm) = gamma_z_est_(nm)
         alpha_est__(nalpha_delta_x + n_alpha*nm) = delta_x_est_(nm)
         alpha_est__(nalpha_delta_y + n_alpha*nm) = delta_y_est_(nm)
         alpha_est__(nalpha_l2_norm + n_alpha*nm) = l2_norm_est_(nm)
         alpha_est__(nalpha_ctf_ind + n_alpha*nm) = 1.0d0*I_ctf_(nm)
         alpha_est__(nalpha_S_index + n_alpha*nm) = S_index_est_(nm)
         alpha_est__(nalpha_M_index + n_alpha*nm) = 1.0d0*nm
      enddo ! do nm=0,n_M-1
      end if ! if (flag_memory_estimate.eqv..false.) then

      if (flag_memory_estimate.eqv..false.) then
      allocate(C_M_(0:1+n_M-1))
      call cs1_c16(n_M,C_M_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  allocate storage for image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..true.) then
      allocate(n_MS_(0:1+n_S-1))
      call cs1_i4(n_S,n_MS_)
      call cl1_i4(n_S,n_MS_)
      allocate(alpha_MS__(0:1+n_alpha*n_MS_max*n_S-1))
      call cs1_r8(n_alpha*n_MS_max*n_S,alpha_MS__)
      call cl1_r8(n_alpha*n_MS_max*n_S,alpha_MS__)
      end if !if (flag_MS_vs_SM.eqv..true.) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  allocate storage for template-image pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..false.) then
      allocate(n_SM_(0:1+n_M-1))
      call cs1_i4(n_M,n_SM_)
      call cl1_i4(n_M,n_SM_)
      allocate(alpha_SM__(0:1+n_alpha*n_SM_max*n_M-1))
      call cs1_r8(n_alpha*n_SM_max*n_M,alpha_SM__)
      call cl1_r8(n_alpha*n_SM_max*n_M,alpha_SM__)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      end if !if (flag_memory_estimate.eqv..false.) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  add unnecessary overfill to check for memory leaks. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A)') ' add unnecessary overfill '
      do ns=0,n_S-1
         do na=quad_n_w_sum,ld_S-1
         S_k_p__(na+ns*ld_S) = dcmplx(1.0d6+1.0d0*ns,2.0d6+1.0d0*ns)
         enddo !do na=quad_n_w_sum,ld_S-1
      enddo !do ns=0,n_S-1
      do nm=0,n_M-1
         do na=quad_n_w_sum,ld_M-1
         M_k_p__(na+nm*ld_M) = dcmplx(1.0d6+1.0d0*nm,2.0d6+1.0d0*nm)
         enddo !do na=quad_n_w_sum,ld_M-1
      enddo !do nm=0,n_M-1
      do nctf=0,n_CTF-1
         do na=quad_n_w_sum,ld_CTF-1
            CTF_k_p__(na+nctf*ld_CTF) = dcmplx(1.0d6+1.0d0*nctf,2.0d6
     $           +1.0d0*nctf)
         enddo !do na=quad_n_w_sum,ld_CTF-1
      enddo !do nctf=0,n_CTF-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $     ,n_M,')-x-(n_S = ',n_S, ') = (', n_M
     $     *n_S,') innerproducts'
      timing_tic = omp_get_wtime()
      call ti8(
     $     verbose !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,quad_n_r !integer *4: length of polar grid (in radial direction). ;
     $     ,quad_n_polar_a_ !integer *4 array (length at least quad_n_r): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $     ,quad_grid_k_p_r_ !real *8 array (length at least quad_n_r): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,quad_weight_k_p_r_ !real *8 array (length at least quad_n_r): weights for radial quadrature. ;
     $     ,quad_weight_k_p_ !real *8 array (length at least quad_n_w_sum_cur): all weights (unrolled) for quadrature, for use with k_p coordinates. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least quad_n_w_sum. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_rq determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_ !real *8 array (length at least n_S): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,n_M !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $     ,I_M_sample_ !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $     ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least quad_n_w_sum. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $     ,n_CTF !integer *4: total number of CTF functions. ;
     $     ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least quad_n_w_sum. ;
     $     ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $     ,alpha_est__ !real *8 array (length at least n_alpha*n_M): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $     ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $     ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
     $     ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_SM_ !integer *4 array (length at least n_M): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
     $     ,alpha_SM__ !real *8 array (length at least n_alpha*n_SM_max*n_M): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
     $     ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_MS_ !integer *4 array (length at least n_S): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
     $     ,alpha_MS__ !real *8 array (length at least n_alpha*n_MS_max*n_S): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
     $     ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'quad_n_r' under consideration (which can change from iteration to iteration). ;
     $     ,displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
     $     ,n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
     $     ,delta_x_ !real *8 array: (length at least n_delta_v): x-coordinates of displacements. ;
     $     ,delta_y_ !real *8 array: (length at least n_delta_v): y-coordinates of displacements. ;
     $     ,n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
     $     ,svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
     $     ,eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
     $     ,flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
     $     ,fpm_howmany_max !integer *4: Maximum number of fftws to call simultaneously within the fftw_plan_many. 
     $     ,n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
     $     ,n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M): actual l2-norm (out to frequency quad_n_r) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
      timing_toc = omp_get_wtime()
      if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
         write(6,'(A,A,F8.3)') 'ti8:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then

      if (flag_memory_estimate.eqv..true.) then
         write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $        nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $        nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         goto 10
      end if !if (flag_memory_estimate.eqv..true.) then

      flag_memory_checkset = .true.
      if (flag_MS_vs_SM.eqv..true.) then
      call cxs_i4(n_S,n_MS_,'n_MS_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_MS_max*n_S,alpha_MS__,'alpha_MS__'
     $     ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      if (flag_MS_vs_SM.eqv..false.) then
      call cxs_i4(n_M,n_SM_,'n_SM_',flag_memory_checkset)
      call cxs_r8(n_alpha*n_SM_max*n_M,alpha_SM__,'alpha_SM__'
     $     ,flag_memory_checkset)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.1) then
         write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.1) then
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

      if (flag_memory_estimate.eqv..false.) then
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (flag_MS_vs_SM.eqv..true.) then
         allocate(ferror_SRTRTCM_(0:1+n_MS_max*n_S-1))
         call cs1_r8(n_MS_max*n_S,ferror_SRTRTCM_)
         allocate(ferror_C_Z_opt_(0:1+n_MS_max*n_S-1))
         call cs1_c16(n_MS_max*n_S,ferror_C_Z_opt_)
         allocate(ferror_CTF_R_S_(0:1+n_MS_max*n_S-1))
         call cs1_r8(n_MS_max*n_S,ferror_CTF_R_S_)
         allocate(ferror_l2_norm_(0:1+n_MS_max*n_S-1))
         call cs1_r8(n_MS_max*n_S,ferror_l2_norm_)
         allocate(ferror_C_M_opt_(0:1+n_MS_max*n_S-1))
         call cs1_c16(n_MS_max*n_S,ferror_C_M_opt_)
         if (verbose.gt.1) then
         write(6,'(1(A,F16.8))')
     $           ' half_diameter_k_c: ' , half_diameter_k_c
         call print_sub_i4(n_S,n_MS_,' n_MS: ')
         end if !if (verbose.gt.1) then
         na=0
         do ns=0,n_S-1
            if (verbose.gt.1) then
               call print_sub_c16(quad_n_w_sum,S_k_p__(I_S_sample_(ns)
     $              *ld_S),' S_k_p_: ')
            end if !if (verbose.gt.1) then
            do nm=0,n_MS_(ns)-1
               if (verbose.gt.1) then
                  write(6,'(A)') ' %%%%%%%%%%%%%%%% '
               end if !if (verbose.gt.1) then
               gamma_z_upd = alpha_MS__(nalpha_gamma_z + (nm + ns
     $              *n_MS_max)*n_alpha)
               delta_x_upd = alpha_MS__(nalpha_delta_x + (nm + ns
     $              *n_MS_max)*n_alpha)
               delta_y_upd = alpha_MS__(nalpha_delta_y + (nm + ns
     $              *n_MS_max)*n_alpha)
               nwav_upd = dsqrt(2.0d0)*half_diameter_k_c
     $              *dsqrt(delta_x_upd**2 +delta_y_upd**2)
               l2_norm_upd = alpha_MS__(nalpha_l2_norm + (nm + ns
     $              *n_MS_max)*n_alpha)
               ctf_ind_upd = alpha_MS__(nalpha_ctf_ind + (nm + ns
     $              *n_MS_max)*n_alpha)
               S_index_upd = alpha_MS__(nalpha_S_index + (nm + ns
     $              *n_MS_max)*n_alpha)
               M_index_upd = alpha_MS__(nalpha_M_index + (nm + ns
     $              *n_MS_max)*n_alpha)
               CTF_R_S_upd = alpha_MS__(nalpha_CTF_R_S + (nm + ns
     $              *n_MS_max)*n_alpha)
               C_Z_opt_upd = alpha_MS__(nalpha_C_Z_opt + (nm + ns
     $              *n_MS_max)*n_alpha)
               call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $              ,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum
     $              ,M_k_p__(I_M_sample_(nint(M_index_upd))*ld_M)
     $              ,M_k_p__(I_M_sample_(nint(M_index_upd)) *ld_M)
     $              ,form_C_M_opt)
               form_C_M_opt = zsqrt(form_C_M_opt)/dsqrt(2.0d0)
               form_C_M_opt = zsqrt(dcmplx(pi*half_diameter_k_c**2
     $              ,0.0d0))
               C_M_opt_upd = C_M_(nint(M_index_upd))
               if (verbose.gt.1) then
                  write(6,'(1(A,2F16.8))') ' C_M_opt_upd: ' ,
     $                 C_M_opt_upd
                  write(6,'(1(A,2F16.8))') ' C_M_(nint(M_index_upd)): '
     $                 , C_M_(nint(M_index_upd))
                  write(6,'(1(A,2F16.8))') ' form_C_M_opt: ' ,
     $                 form_C_M_opt
               end if !if (verbose.gt.1) then
               fix_delta_M_x_c_0 =
     $              fix_delta_M_x_c_0_(I_M_sample_(nint(M_index_upd)))
               fix_delta_M_x_c_1 =
     $              fix_delta_M_x_c_1_(I_M_sample_(nint(M_index_upd)))
               fix_delta_M_x_p_r =
     $              fix_delta_M_x_p_r_(I_M_sample_(nint(M_index_upd)))
               fix_delta_M_x_p_w =
     $              fix_delta_M_x_p_w_(I_M_sample_(nint(M_index_upd)))
               fix_nwav_M = dsqrt(2.0d0)*half_diameter_k_c
     $              *fix_delta_M_x_p_r
               if (verbose.gt.1) then
               write(6,'(2(A,I0),5(A,F16.8))')
     $            ' Image: ' , I_M_sample_(nint(M_index_upd)) 
     $            ,' n_MS: ' , n_MS_(nint(M_index_upd))
     $            ,' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_0 
     $            ,' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_1 
     $            ,' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r 
     $            ,' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w 
     $            ,' fix_nwav_M: ' , fix_nwav_M 
               end if !if (verbose.gt.1) then
               gamma_z_est = gamma_z_est_(nint(M_index_upd))
               delta_x_est = delta_x_est_(nint(M_index_upd))
               delta_y_est = delta_y_est_(nint(M_index_upd))
               nwav_est = dsqrt(2.0d0)*half_diameter_k_c
     $              *dsqrt(delta_x_est**2 +delta_y_est**2)
               if (verbose.gt.1) then
               write(6,'(4(A,F16.8))')
     $           ' gamma_z_est: ' , gamma_z_est
     $           ,' delta_x_est: ' , delta_x_est
     $           ,' delta_y_est: ' , delta_y_est
     $           ,' nwav_est: ' , nwav_est
               end if !if (verbose.gt.1) then
               fix_phi_CTF = fix_phi_CTF_(nint(ctf_ind_upd))
               if (verbose.gt.1) then
               write(6,'(A,F16.8,A,I0,3(A,F16.8))') 
     $              ' M_index: ' , M_index_upd
     $              ,' ns: ' , ns
     $              ,' S_index: ' , S_index_upd
     $              ,' ctf_ind: ' , ctf_ind_upd
     $              ,' fix_phi_CTF: ' , fix_phi_CTF
               end if !if (verbose.gt.1) then
               fix_phi_S =
     $              fix_phi_S_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_c_0 =
     $              fix_delta_S_x_c_0_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_c_1 =
     $              fix_delta_S_x_c_1_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_p_r =
     $              fix_delta_S_x_p_r_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_p_w =
     $              fix_delta_S_x_p_w_(I_S_sample_(nint(S_index_upd)))
               fix_nwav_S = dsqrt(2.0d0)*half_diameter_k_c
     $              *fix_delta_S_x_p_r
               if (verbose.gt.1) then
               write(6,'(1(A,I0),6(A,F16.8))')
     $              ' Template: ' , I_S_sample_(nint(S_index_upd)) 
     $              ,' fix_phi_S: ' , fix_phi_S
     $              ,' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_0 
     $              ,' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_1 
     $              ,' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r 
     $              ,' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w 
     $              ,' fix_nwav_S: ' , fix_nwav_S
               write(6,'(4(A,F16.8))')
     $              ' gamma_z_upd: ' , gamma_z_upd
     $              ,' delta_x_upd: ' , delta_x_upd
     $              ,' delta_y_upd: ' , delta_y_upd
     $              ,' nwav_upd: ' , nwav_upd
               write(6,'(3(A,F16.8))')
     $              ' l2_norm_upd: ' , l2_norm_upd
     $              ,' CTF_R_S_upd: ' , CTF_R_S_upd
     $              ,' C_Z_opt_upd: ' , C_Z_opt_upd
               end if !if (verbose.gt.1) then
               SRTRTCM_upd = C_Z_opt_upd*CTF_R_S_upd*C_M_opt_upd
               if (verbose.gt.1) then
               write(6,'(A,F16.8)') ' SRTRTCM_upd: ' , SRTRTCM_upd
               call print_sub_c16(quad_n_w_sum
     $              ,S_k_p__(I_S_sample_(nint(S_index_upd))*ld_S)
     $              ,' S_k_p_: ')
               call print_sub_c16(quad_n_w_sum
     $              ,CTF_k_p__(nint(ctf_ind_upd)*ld_CTF)
     $              ,' CTF_k_p_: ')
               write(6,'(A,F16.12)') ' gamma_z_est + gamma_z_upd: ' ,
     $              gamma_z_est +gamma_z_upd
               end if !if (verbose.gt.1) then
               form_CTF_R_S = l2_xRPx_f(
     $              half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_phi_S
     $              ,+gamma_z_est + gamma_z_upd
     $              )
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ' form_CTF_R_S: ' , form_CTF_R_S
               end if !if (verbose.gt.1) then
               form_SRTRTCM = I_PxRTRTxP_f(
     $              half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,flag_RTRT_vs_RTTR
     $              ,gamma_z_est
     $              ,delta_x_est
     $              ,delta_y_est
     $              ,gamma_z_upd
     $              ,delta_x_upd
     $              ,delta_y_upd
     $              )
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ' form_SRTRTCM: ' , form_SRTRTCM
               end if !if (verbose.gt.1) then
               form_C_Z_opt = form_SRTRTCM/(form_CTF_R_S*form_C_M_opt)
               if (verbose.gt.1) then
               write(6,'(A,2F16.8)')
     $              ,' form_C_Z_opt: ' , form_C_Z_opt
               end if !if (verbose.gt.1) then
               form_l2_norm = form_SRTRTCM/form_CTF_R_S**2
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ,' form_l2_norm: ' , form_l2_norm
               end if !if (verbose.gt.1) then
               error_l2_norm = dabs(form_l2_norm - l2_norm_upd)
     $              /dabs(form_l2_norm)
               error_C_Z_opt = zabs(form_C_Z_opt - C_Z_opt_upd)
     $              /zabs(form_C_Z_opt)
               error_SRTRTCM = dabs(form_SRTRTCM - SRTRTCM_upd)
     $              /dabs(form_SRTRTCM)
               error_CTF_R_S = dabs(form_CTF_R_S - CTF_R_S_upd)
     $              /dabs(form_CTF_R_S)
               error_C_M_opt = zabs(form_C_M_opt - C_M_opt_upd)
     $              /zabs(form_C_M_opt)
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ,' error_l2_norm: ' , error_l2_norm
               write(6,'(A,2F16.8)')
     $              ,' error_C_Z_opt: ' , error_C_Z_opt
               write(6,'(A,F16.8)')
     $              ,' error_SRTRTCM: ' , error_SRTRTCM
               write(6,'(A,F16.8)')
     $              ,' error_CTF_R_S: ' , error_CTF_R_S
               write(6,'(A,2F16.8)')
     $              ,' error_C_M_opt: ' , error_C_M_opt
               end if !if (verbose.gt.1) then
               ferror_l2_norm_(na) = error_l2_norm
               ferror_C_Z_opt_(na) = error_C_Z_opt
               ferror_SRTRTCM_(na) = error_SRTRTCM
               ferror_CTF_R_S_(na) = error_CTF_R_S
               ferror_C_M_opt_(na) = error_C_M_opt
               if (verbose.gt.1) then
               call ti8_dr_digits_single(
     $              verbose
     $              ,half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,flag_RTRT_vs_RTTR
     $              ,gamma_z_est
     $              ,delta_x_est
     $              ,delta_y_est
     $              ,gamma_z_upd
     $              ,delta_x_upd
     $              ,delta_y_upd
     $              )
               end if !if (verbose.gt.1) then
              na=na+1
            enddo !do nm=0,n_MS_(ns)-1
            if (verbose.gt.1) then
               write(6,'(A)') ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
               write(6,'(A)') ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
            end if !if (verbose.gt.1) then
         enddo !do ns=0,n_M-1
         if (verbose.gt.-1) then
            write(6,'(A,5A)') ' '
     $           ,'  l2  ferror_C_Z_opt'
     $           ,'  l2  ferror_l2_norm'
     $           ,'  l2  ferror_SRTRTCM'
     $           ,'  l2  ferror_CTF_R_S'
     $           ,'  l2  ferror_C_M_opt'
            write(6,'(A,5F20.16)') ' '
     $           ,al2_c16_f(na,ferror_C_Z_opt_)
     $           ,al2_r8_f(na,ferror_l2_norm_)
     $           ,al2_r8_f(na,ferror_SRTRTCM_)
     $           ,al2_r8_f(na,ferror_CTF_R_S_)
     $           ,al2_c16_f(na,ferror_C_M_opt_)
            write(6,'(A,5A)') ' '
     $           ,'  max ferror_C_Z_opt'
     $           ,'  max ferror_l2_norm'
     $           ,'  max ferror_SRTRTCM'
     $           ,'  max ferror_CTF_R_S'
     $           ,'  max ferror_C_M_opt'
            write(6,'(A,5F20.16)') ' '
     $           ,maxr_c16_f(na,ferror_C_Z_opt_)
     $           ,max_r8_f(na,ferror_l2_norm_)
     $           ,max_r8_f(na,ferror_SRTRTCM_)
     $           ,max_r8_f(na,ferror_CTF_R_S_)
     $           ,maxr_c16_f(na,ferror_C_M_opt_)
         end if !if (verbose.gt.-1) then
      end if !if (flag_MS_vs_SM.eqv..true.) then
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (flag_MS_vs_SM.eqv..false.) then
         allocate(ferror_SRTRTCM_(0:1+n_SM_max*n_M-1))
         call cs1_r8(n_SM_max*n_M,ferror_SRTRTCM_)
         allocate(ferror_C_Z_opt_(0:1+n_SM_max*n_M-1))
         call cs1_c16(n_SM_max*n_M,ferror_C_Z_opt_)
         allocate(ferror_CTF_R_S_(0:1+n_SM_max*n_M-1))
         call cs1_r8(n_SM_max*n_M,ferror_CTF_R_S_)
         allocate(ferror_l2_norm_(0:1+n_SM_max*n_M-1))
         call cs1_r8(n_SM_max*n_M,ferror_l2_norm_)
         allocate(ferror_C_M_opt_(0:1+n_SM_max*n_M-1))
         call cs1_c16(n_SM_max*n_M,ferror_C_M_opt_)
         if (verbose.gt.1) then
         write(6,'(1(A,F16.8))')
     $           ' half_diameter_k_c: ' , half_diameter_k_c
         call print_sub_i4(n_M,n_SM_,' n_SM: ')
         end if !if (verbose.gt.1) then
         na=0
         do nm=0,n_M-1
            if (verbose.gt.1) then
               call print_sub_c16(quad_n_w_sum,M_k_p__(I_M_sample_(nm)
     $              *ld_M),' M_k_p_: ')
            end if !if (verbose.gt.1) then
            call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $           ,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum
     $           ,M_k_p__(I_M_sample_(nm)*ld_M),M_k_p__(I_M_sample_(nm)
     $           *ld_M),form_C_M_opt)
            form_C_M_opt = zsqrt(form_C_M_opt)/dsqrt(2.0d0)
            form_C_M_opt = zsqrt(dcmplx(pi*half_diameter_k_c**2,0.0d0))
            C_M_opt_upd = C_M_(nm)
            if (verbose.gt.1) then
            write(6,'(1(A,2F16.8))')
     $           ' C_M_opt_upd: ' , C_M_opt_upd
            write(6,'(1(A,2F16.8))')
     $           ' C_M_(nm): ' , C_M_(nm)
            write(6,'(1(A,2F16.8))')
     $           ' form_C_M_opt: ' , form_C_M_opt
            end if !if (verbose.gt.1) then
            fix_delta_M_x_c_0 = fix_delta_M_x_c_0_(I_M_sample_(nm))
            fix_delta_M_x_c_1 = fix_delta_M_x_c_1_(I_M_sample_(nm))
            fix_delta_M_x_p_r = fix_delta_M_x_p_r_(I_M_sample_(nm))
            fix_delta_M_x_p_w = fix_delta_M_x_p_w_(I_M_sample_(nm))
            fix_nwav_M = dsqrt(2.0d0)*half_diameter_k_c
     $           *fix_delta_M_x_p_r
            if (verbose.gt.1) then
            write(6,'(2(A,I0),5(A,F16.8))')
     $            ' Image: ' , I_M_sample_(nm) , ' n_SM: ' , n_SM_(nm)
     $            ,' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_0 
     $            ,' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_1 
     $            ,' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r 
     $            ,' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w 
     $            ,' fix_nwav_M: ' , fix_nwav_M 
            end if !if (verbose.gt.1) then
            gamma_z_est = gamma_z_est_(nm)
            delta_x_est = delta_x_est_(nm)
            delta_y_est = delta_y_est_(nm)
            nwav_est = dsqrt(2.0d0)*half_diameter_k_c*dsqrt(delta_x_est
     $           **2 +delta_y_est**2)
            if (verbose.gt.1) then
            write(6,'(4(A,F16.8))')
     $           ' gamma_z_est: ' , gamma_z_est
     $           ,' delta_x_est: ' , delta_x_est
     $           ,' delta_y_est: ' , delta_y_est
     $           ,' nwav_est: ' , nwav_est
            end if !if (verbose.gt.1) then
            do ns=0,n_SM_(nm)-1
               if (verbose.gt.1) then
                  write(6,'(A)') ' %%%%%%%%%%%%%%%% '
               end if !if (verbose.gt.1) then
               gamma_z_upd = alpha_SM__(nalpha_gamma_z + (ns + nm
     $              *n_SM_max)*n_alpha)
               delta_x_upd = alpha_SM__(nalpha_delta_x + (ns + nm
     $              *n_SM_max)*n_alpha)
               delta_y_upd = alpha_SM__(nalpha_delta_y + (ns + nm
     $              *n_SM_max)*n_alpha)
               nwav_upd = dsqrt(2.0d0)*half_diameter_k_c
     $              *dsqrt(delta_x_upd**2 +delta_y_upd**2)
               l2_norm_upd = alpha_SM__(nalpha_l2_norm + (ns + nm
     $              *n_SM_max)*n_alpha)
               ctf_ind_upd = alpha_SM__(nalpha_ctf_ind + (ns + nm
     $              *n_SM_max)*n_alpha)
               S_index_upd = alpha_SM__(nalpha_S_index + (ns + nm
     $              *n_SM_max)*n_alpha)
               M_index_upd = alpha_SM__(nalpha_M_index + (ns + nm
     $              *n_SM_max)*n_alpha)
               CTF_R_S_upd = alpha_SM__(nalpha_CTF_R_S + (ns + nm
     $              *n_SM_max)*n_alpha)
               C_Z_opt_upd = alpha_SM__(nalpha_C_Z_opt + (ns + nm
     $              *n_SM_max)*n_alpha)
               fix_phi_CTF = fix_phi_CTF_(nint(ctf_ind_upd))
               if (verbose.gt.1) then
               write(6,'(A,F16.8,A,I0,3(A,F16.8))') 
     $              ' M_index: ' , M_index_upd
     $              ,' ns: ' , ns
     $              ,' S_index: ' , S_index_upd
     $              ,' ctf_ind: ' , ctf_ind_upd
     $              ,' fix_phi_CTF: ' , fix_phi_CTF
               end if !if (verbose.gt.1) then
               fix_phi_S =
     $              fix_phi_S_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_c_0 =
     $              fix_delta_S_x_c_0_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_c_1 =
     $              fix_delta_S_x_c_1_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_p_r =
     $              fix_delta_S_x_p_r_(I_S_sample_(nint(S_index_upd)))
               fix_delta_S_x_p_w =
     $              fix_delta_S_x_p_w_(I_S_sample_(nint(S_index_upd)))
               fix_nwav_S = dsqrt(2.0d0)*half_diameter_k_c
     $              *fix_delta_S_x_p_r
               if (verbose.gt.1) then
               write(6,'(1(A,I0),6(A,F16.8))')
     $              ' Template: ' , I_S_sample_(nint(S_index_upd)) 
     $              ,' fix_phi_S: ' , fix_phi_S
     $              ,' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_0 
     $              ,' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_1 
     $              ,' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r 
     $              ,' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w 
     $              ,' fix_nwav_S: ' , fix_nwav_S
               write(6,'(4(A,F16.8))')
     $              ' gamma_z_upd: ' , gamma_z_upd
     $              ,' delta_x_upd: ' , delta_x_upd
     $              ,' delta_y_upd: ' , delta_y_upd
     $              ,' nwav_upd: ' , nwav_upd
               write(6,'(3(A,F16.8))')
     $              ' l2_norm_upd: ' , l2_norm_upd
     $              ,' CTF_R_S_upd: ' , CTF_R_S_upd
     $              ,' C_Z_opt_upd: ' , C_Z_opt_upd
               end if !if (verbose.gt.1) then
               SRTRTCM_upd = C_Z_opt_upd*CTF_R_S_upd*C_M_opt_upd
               if (verbose.gt.1) then
               write(6,'(A,F16.8)') ' SRTRTCM_upd: ' , SRTRTCM_upd
               call print_sub_c16(quad_n_w_sum
     $              ,S_k_p__(I_S_sample_(nint(S_index_upd))*ld_S)
     $              ,' S_k_p_: ')
               call print_sub_c16(quad_n_w_sum
     $              ,CTF_k_p__(nint(ctf_ind_upd)*ld_CTF)
     $              ,' CTF_k_p_: ')
               write(6,'(A,F16.12)') ' gamma_z_est + gamma_z_upd: ' ,
     $              gamma_z_est +gamma_z_upd
               end if !if (verbose.gt.1) then
               form_CTF_R_S = l2_xRPx_f(
     $              half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_phi_S
     $              ,+gamma_z_est + gamma_z_upd
     $              )
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ' form_CTF_R_S: ' , form_CTF_R_S
               end if !if (verbose.gt.1) then
               form_SRTRTCM = I_PxRTRTxP_f(
     $              half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,flag_RTRT_vs_RTTR
     $              ,gamma_z_est
     $              ,delta_x_est
     $              ,delta_y_est
     $              ,gamma_z_upd
     $              ,delta_x_upd
     $              ,delta_y_upd
     $              )
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ' form_SRTRTCM: ' , form_SRTRTCM
               end if !if (verbose.gt.1) then
               form_C_Z_opt = form_SRTRTCM/(form_CTF_R_S*form_C_M_opt)
               if (verbose.gt.1) then
               write(6,'(A,2F16.8)')
     $              ,' form_C_Z_opt: ' , form_C_Z_opt
               end if !if (verbose.gt.1) then
               form_l2_norm = form_SRTRTCM/form_CTF_R_S**2
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ,' form_l2_norm: ' , form_l2_norm
               end if !if (verbose.gt.1) then
               error_l2_norm = dabs(form_l2_norm - l2_norm_upd)
     $              /dabs(form_l2_norm)
               error_C_Z_opt = zabs(form_C_Z_opt - C_Z_opt_upd)
     $              /zabs(form_C_Z_opt)
               error_SRTRTCM = dabs(form_SRTRTCM - SRTRTCM_upd)
     $              /dabs(form_SRTRTCM)
               error_CTF_R_S = dabs(form_CTF_R_S - CTF_R_S_upd)
     $              /dabs(form_CTF_R_S)
               error_C_M_opt = zabs(form_C_M_opt - C_M_opt_upd)
     $              /zabs(form_C_M_opt)
               if (verbose.gt.1) then
               write(6,'(A,F16.8)')
     $              ,' error_l2_norm: ' , error_l2_norm
               write(6,'(A,2F16.8)')
     $              ,' error_C_Z_opt: ' , error_C_Z_opt
               write(6,'(A,F16.8)')
     $              ,' error_SRTRTCM: ' , error_SRTRTCM
               write(6,'(A,F16.8)')
     $              ,' error_CTF_R_S: ' , error_CTF_R_S
               write(6,'(A,2F16.8)')
     $              ,' error_C_M_opt: ' , error_C_M_opt
               end if !if (verbose.gt.1) then
               ferror_l2_norm_(na) = error_l2_norm
               ferror_C_Z_opt_(na) = error_C_Z_opt
               ferror_SRTRTCM_(na) = error_SRTRTCM
               ferror_CTF_R_S_(na) = error_CTF_R_S
               ferror_C_M_opt_(na) = error_C_M_opt
               if (verbose.gt.1) then
               call ti8_dr_digits_single(
     $              verbose
     $              ,half_diameter_k_c
     $              ,fix_phi_CTF
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,flag_RTRT_vs_RTTR
     $              ,gamma_z_est
     $              ,delta_x_est
     $              ,delta_y_est
     $              ,gamma_z_upd
     $              ,delta_x_upd
     $              ,delta_y_upd
     $              )
               end if !if (verbose.gt.1) then
               na=na+1
            enddo !do ns=0,n_SM_(nm)-1
            if (verbose.gt.1) then
               write(6,'(A)') ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
               write(6,'(A)') ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
            end if !if (verbose.gt.1) then
         enddo !do nm=0,n_M-1
         if (verbose.gt.-1) then
            write(6,'(A,5A)') ' '
     $           ,'  l2  ferror_C_Z_opt'
     $           ,'  l2  ferror_l2_norm'
     $           ,'  l2  ferror_SRTRTCM'
     $           ,'  l2  ferror_CTF_R_S'
     $           ,'  l2  ferror_C_M_opt'
            write(6,'(A,5F20.16)') ' '
     $           ,al2_c16_f(na,ferror_C_Z_opt_)
     $           ,al2_r8_f(na,ferror_l2_norm_)
     $           ,al2_r8_f(na,ferror_SRTRTCM_)
     $           ,al2_r8_f(na,ferror_CTF_R_S_)
     $           ,al2_c16_f(na,ferror_C_M_opt_)
            write(6,'(A,5A)') ' '
     $           ,'  max ferror_C_Z_opt'
     $           ,'  max ferror_l2_norm'
     $           ,'  max ferror_SRTRTCM'
     $           ,'  max ferror_CTF_R_S'
     $           ,'  max ferror_C_M_opt'
            write(6,'(A,5F20.16)') ' '
     $           ,maxr_c16_f(na,ferror_C_Z_opt_)
     $           ,max_r8_f(na,ferror_l2_norm_)
     $           ,max_r8_f(na,ferror_SRTRTCM_)
     $           ,max_r8_f(na,ferror_CTF_R_S_)
     $           ,maxr_c16_f(na,ferror_C_M_opt_)
         end if !if (verbose.gt.-1) then
      end if !if (flag_MS_vs_SM.eqv..false.) then
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end if !if (flag_memory_estimate.eqv..false.) then

      include 'ti8_dr_digits_excerpt_checkset.f'

 10    continue

      if (verbose.gt.1) then
         write(6,'(A)') '[finished ti8_dr_digits]'
      end if

      stop
      end
