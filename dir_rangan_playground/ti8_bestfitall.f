      subroutine ti8_getbestfit(
     $     verbose !integer *4: verbosity level. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*n_M): array of image values in k_p coordinates. also called ximagesk. ;
     $     ,ld_M !integer *4: size of each individual image (i.e., each individual M_k_p_). also called imagesize ;
     $     ,n_M !integer *4: number of images. also called nimages ;
     $     ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_CTF): array of CTF values in k_p_ coordinates. We assume here that ld_CTF.eq.ld_M. also called ctfw ;
     $     ,n_CTF !integer *4: number of individual CTFS. also called nctf ;
     $     ,CTF_1pind_ !integer *4 array (length at least n_M): list of indices (base 1), indicating which of the n_CTF CTFs is associated with each image. also called ctf_ind ;
     $     ,n_w_csum_0in_ !integer *4 array (length at least n_k_p_r_cur); cumulative sum of n_w_. also called icstart ;
     $     ,n_k_p_r_cur !integer *4: current number of shells. also called ncur ;
     $     ,n_w_ !integer *4 array (length at least n_k_p_r_cur); number of points (angles) for each nk. also called ngridc ; Note that we expect each entry of n_w_ to be even (for consistency with get_template_size). ;
     $     ,alpha_first5__ !real *8 array (length at least n_alpha_first5*n_M): estimated image parameters associated with sampled image set stored in 1-dimensional array. The ordering of the image parameters is [polar_a, azimu_b, gamma_z, delta_x, delta_y], which should match the first5 elements in excerpt_define_nalpha.f. also called alphas ;
     $     ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): radial values for each shell. also called xnodesr ;
     $     ,kweight_k_p_r_ !real *8 array (length at least n_k_p_r_cur): radial quadrature weights for each shell. also called wtsr ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*n_S): array of template values in k_p coordinates. also called templates ;
     $     ,ld_S !integer *4: size of each individual template (i.e., each individual S_k_p_). also called itemplatesize ;
     $     ,n_S !integer *4: number of templates. also called ntemplates ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S): polar_a associated with each template S_k_p_. also called thetas ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S): azimu_b associated with each template S_k_p_. also called phis ;
     $     ,n_delta_v !integer *4: number of translations. also called ndelta ;
     $     ,delta__ !real *8 array (length at least 2*n_delta_v): array of translations, such that the delta_x = delta__(0+ndv*2) and delta_y = delta__(1+ndv*2) for translation ndv. also called alltrans ;
     $     )
      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      complex *16 M_k_p__(0:0) !complex *16 array (length at least ld_M*n_M): array of image values in k_p coordinates. also called ximagesk. ;
      integer *4 ld_M !integer *4: size of each individual image (i.e., each individual M_k_p_). also called imagesize ;
      integer *4 n_M !integer *4: number of images. also called nimages ;
      integer *4 n_M_sample !integer *4: number of sampled images. .eq.n_M ;
      complex *16 CTF_k_p__(0:0) !complex *16 array (length at least ld_CTF*n_CTF): array of CTF values in k_p_ coordinates. We assume here that ld_CTF.eq.ld_M. also called ctfw ;
      integer *4 n_CTF !integer *4: number of individual CTFS. also called nctf ;
      integer *4 CTF_1pind_(0:0) !integer *4 array (length at least n_M): list of indices (base 1), indicating which of the n_CTF CTFs is associated with each image. also called ctf_ind ;
      integer *4 n_w_csum_0in_(0:0) !integer *4 array (length at least n_k_p_r_cur); cumulative sum of n_w_. also called icstart ;
      integer *4 n_k_p_r_cur !integer *4: current number of shells. also called ncur ;
      integer *4 n_w_(0:0) !integer *4 array (length at least n_k_p_r_cur); number of points (angles) for each nk. also called ngridc ; Note that we expect each entry of n_w_ to be even (for consistency with get_template_size). ;
      real *8 alpha_first5__(0:0) !real *8 array (length at least n_alpha_first5*n_M): estimated image parameters associated with sampled image set stored in 1-dimensional array. The ordering of the image parameters is [polar_a, azimu_b, gamma_z, delta_x, delta_y], which should match the first5 elements in excerpt_define_nalpha.f. also called alphas ;
      real *8 grid_k_p_r_(0:0) !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
      real *8 kweight_k_p_r_(0:0) !real *8 array (length at least n_k_p_r_cur): weights for radial quadrature. ;
      complex *16 S_k_p__(0:0) !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
      integer *4 ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_w_sum. ;
      integer *4 n_S !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
      integer *4 n_S_sample !integer *4: number of sampled templates. .eq. n_S ;
      real *8 S_alpha_polar_a_(0:0) !real *8 array (length at least n_S): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8 S_alpha_azimu_b_(0:0) !real *8 array (length at least n_S): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      integer *4 n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
      real *8 delta__(0:0) !real *8 array (length at least 2*n_delta_v): array of translations, such that the delta_x = delta__(0+ndv*2) and delta_y = delta__(1+ndv*2) for translation ndv. also called alltrans ;
c$$$  %%%%%%%%
      real *8 pi
c$$$  %%%%%%%%
      integer *4 n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
      parameter(n_omp_sub_0in=1) !This should be set by the user. perhaps 8? ;
      integer *4 n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ; This will be determined internally. ;
      integer *4 n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ; This will be determined internally. ;
      integer *4 n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ; This will be determined internally. ;
      integer *4 n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ; This will be determined internally. ;
      real *8 d_memory_limit !real *8: upper limit on the memory allowed by ti8 (in bytes). ;
      parameter (d_memory_limit=4.0d9) !This should be set by the user. perhaps 4GB = 4.0d9? ;
      real *8 d_memory_estimate !real *8: estimate of required memory (in bytes). ;
      logical flag_continue,flag_proceed
      integer *4 nS_0_sub_tmp,nS_1_sub_tmp,nM_0_sub_tmp,nM_1_sub_tmp
      integer *4 n_S_0_sub_max,n_S_1_sub_max,n_M_0_sub_max,n_M_1_sub_max
c$$$  %%%%%%%%
      logical flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ; This will be set to .false. to assign templates to images (as leslie likes to do). ;
      integer *4 n_w_sum !integer *4: total number of points on all rings. ;
c$$$      indices
      integer *4 n_r,nr,n_w_max,na,nw
      integer n_polar_a,npolar_a,n_azimu_b,n_polar_a_azimu_b_sum
      integer *4, allocatable :: n_polar_a_(:) !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ; This must be half of n_w_. ;
      integer *4, allocatable :: n_w_csum_(:) !integer *4 array (length at least n_k_p_r_cur): cumulative sum of n_w_. ;
      integer *4 max_i4_f,sum_i4_f !function output. ;
      real *8 max_r8_f,avg_r8_f,al2_r8_f !function output. ;
      real *8 half_diameter_x_c !real *8: half the diameter of the box in x_c. ;
      real *8, allocatable :: weight_k_p_r_(:) !real *8 array (length at least n_w_sum): quadrature weights for all points. ;
      real *8, allocatable :: weight_k_p_(:) !real *8 array (length at least n_w_sum): quadrature weights for all points. ;
c$$$      Templates S_
      integer *4 n_S_max !integer *4: number of templates in stack (equal to n_S). ;
      integer *4 ns,nS_sample
      integer *4, allocatable :: I_S_sample_(:) !integer *4 array (length at least n_S): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ; in this case set to 0:n_S-1. ;
c$$$      Images M_
      integer *4 n_M_max !integer *4; number of images in stack (equal to n_M). ;
      integer *4 nm,nM_sample
      integer *4, allocatable :: I_M_sample_(:) !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ; in this case set to 0:n_M-1. ;
c$$$      CTF-functions
      integer *4 ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_w_sum. ; we assume ld_CTF.eq.ld_M ;
      integer *4 nctf
c$$$      fftw_plan_many 
      integer  *4 fpm_howmany_max !integer *4: number of ffts to call simultaneously with fftw. I've tinkered with this number, 16 to 32 seems reasonable. ;
c$$$      image parameters 
      include 'excerpt_define_nalpha.f'
      integer *4 n_alpha_first5 !integer *4: set to 5. ;
      parameter(n_alpha_first5=5)
      real *8, allocatable :: alpha_est__(:) !real *8 array (length at least n_alpha*n_M): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
      complex *16, allocatable :: C_M_(:) !complex *16 array (length at least n_M): actual l2-norm (out to frequency n_k_p_r_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
      real *8, allocatable :: S_alpha_S_index_(:) !real *8 array (length at least n_S): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to equal ns (ranging from 0:n_S-1) which in turn refers to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8 alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
      parameter (alpha_update_f=0.0d0) !when set to 0.0d0, we simply pick the best template for each image. ;
c$$$      range of displacement and rotation parameters to measure
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
      parameter (displacement_max=7.0d0) !as a default we set this to 7.0d0, which should exceed 2.0d0*pi*half_diameter_x_c. This essentially corresponds to no constraint. ;
      real *8 angle_tmp,delta_tmp,dtmp,dr,tmp_k,tmp_w
      integer *4 ndv
      real *8 delta_x,delta_y,gamma_z !temporary. ;
      real *8 delta_0in_(0:1) !temporary: for use with rotate_delta. ;
      real *8 delta_out_(0:1) !temporary: for use with rotate_delta. ;
      integer svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
      real *8 eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
      parameter (eps_svd=0.01d0) !e.g., 2 digits of precision for the svd expansion. ;
      real *8 N_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_p_r_cur' under consideration (which can change from iteration to iteration). ;
      real *8, allocatable :: delta_x_(:) !real *8 array: (length at least n_delta_v): x-coordinates of displacements. ;
      real *8, allocatable :: delta_y_(:) !real *8 array: (length at least n_delta_v): y-coordinates of displacements. ;
      integer *4 n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
      integer *4 ngz
      real *8, allocatable :: gamma_z_(:) !real *8 array: (length at least n_gamma_z): list of angles for in-plane rotation. ;
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
c$$$      parameters for timing
      real *8 timing_tic,timing_toc,timing_innerproduct
c$$$      temporary arrays for testing
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: CTF_k_p_(:)
      complex *16, allocatable :: R_S_(:)
      complex *16, allocatable :: T_R_S_(:)
      complex *16, allocatable :: T_S_(:)
      complex *16, allocatable :: R_T_S_(:)
      complex *16, allocatable :: CTF_R_S_(:)
      complex *16, allocatable :: CTF_T_R_S_(:)
      complex *16, allocatable :: CTF_R_T_S_(:)
      complex *16, allocatable :: M_k_p_(:)
      complex *16 CTF_R_S_norm
      complex *16 M_norm
      complex *16 S_T_R_CTF_M,S_R_T_CTF_M
      complex *16 l2_norm_RTRT,l2_norm_RTTR
      complex *16 C_Z_opt_RTRT,C_Z_opt_RTTR
c$$$      random seed
      integer *4 rseed !random seed. ;
c$$$      
      logical flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
      logical flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
      logical flag_memory_checkset !temporary: logical flag used in ti8_dr_alignment_excerpt_checkset.f ;
      
      ld_CTF = ld_M
      n_M_sample = n_M
      n_S_sample = n_S
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering ti8_bestfitall] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
      end if !if (verbose.gt.0) then
      call openblas_set_num_threads(1)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ',verbose
         call print_sub_c16(ld_M*n_M,M_k_p__,' M_k_p__: ')
         write(6,'(A,I0)') ' ld_M: ',ld_M
         write(6,'(A,I0)') ' n_M: ',n_M
         call print_sub_c16(ld_M*n_CTF,CTF_k_p__,' CTF_k_p__: ')
         write(6,'(A,I0)') ' n_CTF: ',n_CTF
         call print_sub_i4(n_M,CTF_1pind_,' CTF_1pind_: ')
         call print_sub_i4(n_k_p_r_cur,n_w_csum_0in_,' n_w_csum_0in_: ')
         write(6,'(A,I0)') ' n_k_p_r_cur: ',n_k_p_r_cur
         call print_sub_i4(n_k_p_r_cur,n_w_,' n_w_: ')
         call print_sub_r8(n_alpha_first5*n_M,alpha_first5__
     $        ,' alpha_first5__: ')
         call print_sub_r8(n_k_p_r_cur,grid_k_p_r_,' grid_k_p_r_: ')
         call print_sub_r8(n_k_p_r_cur,kweight_k_p_r_
     $        ,' kweight_k_p_r_: ')
         call print_sub_c16(ld_S*n_S,S_k_p__,' S_k_p__: ')
         write(6,'(A,I0)') ' ld_S: ',ld_S
         write(6,'(A,I0)') ' n_S: ',n_S
         call print_sub_r8(n_S,S_alpha_polar_a_,' S_alpha_polar_a_: ')
         call print_sub_r8(n_S,S_alpha_azimu_b_,' S_alpha_azimu_b_: ')
         write(6,'(A,I0)') ' n_delta_v: ',n_delta_v
         call print_sub_r8(2*n_delta_v,delta__,' delta__: ')
      end if !if (verbose.gt.0) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  initialize parameters
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      rseed = 1 ! used for random numbers in ti8. ;
      pi = 4.0d0*datan(1.0d0)
      svd_calculation_type = 2 ! do not use svd expansion. Set to 1 to use svd expansion. if set to 1, then eps_svd and N_pixels_in will become important. ;
      flag_RTRT_vs_RTTR = .false. ! Use this to calculate < T_{\delta_{upd}} R_{\gamma_{upd}} S , CTF.*M > . ;
      n_S_0_sub_0in = 1 ! will be udated automatically later. ;
      n_S_1_sub_0in = 1 ! will be udated automatically later. ;
      n_M_0_sub_0in = 1 ! will be udated automatically later. ;
      n_M_1_sub_0in = 1 ! will be udated automatically later. ;
      n_S_max = n_S ! use all templates. ;
      n_M_max = n_M ! use all images. ;
      fpm_howmany_max = 16 !integer: number of ffts to call simultaneously with fftw. I've tinkered with this number, 16 to 32 seems reasonable. ;
      flag_MS_vs_SM = .false. !logical: set to .false. to store the best templates for each image. ;
      n_SM_max = min(n_S,10) ! store 10 templates per image. ;
      n_MS_max = min(n_M,10) ! store 10 images per template. ; 
      flag_time_Zstore = .false. ! do not time Zstore step in ti8. ;
c$$$      tesselation_distance_req = 0.25
      tesselation_distance_req = 2.0d0 !reduce this to use local search. ;
      n_LT_add = max(1,n_M/10) ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      n_LT_ref = max(1,n_SM_max/2) ! number of image-template pairs to consider when refining local search.
c$$$  %%%%%%%%
      n_r = n_k_p_r_cur
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  define n_polar_a_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_w_max = max_i4_f(n_r,n_w_)
      if (mod(n_w_max,2).ne.0) then
         write(6,'(A,A)') ' Warning, get_template_size implicitly ' ,
     $        'assumes that n_w_max is even '
         stop !exit program due to error. ;
      end if !if (mod(n_w_max,2).ne.0) then
      n_polar_a = n_w_max/2
      allocate(n_polar_a_(0:1+n_r-1))
      call cs1_i4(n_r,n_polar_a_)
      do nr=0,n_r-1
         nw = n_w_(nr)
         if (mod(nw,2).ne.0) then
            write(6,'(A,A)') ' Warning, get_template_size implicitly ' ,
     $           'assumes that n_w is even '
            stop !exit program due to error. ;
         end if !if (mod(nw,2).ne.0) then
         n_polar_a_(nr) = n_w_(nr)/2
      enddo !do nr=0,n_r-1
      if (verbose.gt.1) then
         call print_sub_i4(n_r,n_polar_a_,' n_polar_a_: ')
      end if !if (verbose.gt.1) then
c$$$  %%%%%%%%
      n_w_sum = sum_i4_f(n_r,n_w_)
      n_gamma_z = n_w_max
      fpm_howmany_max = max(1,fpm_howmany_max)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  define n_w_csum_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(n_w_csum_(0:1+n_r-1))
      call cs1_i4(n_r,n_w_csum_)
      nw=1
      do nr=0,n_r-1
         n_w_csum_(nr) = nw
         if (n_w_csum_(nr).ne.n_w_csum_0in_(nr)) then
            write(6,'(3(A,I0))')
     $           ' Warning, at nr ' , nr
     $           , ' n_w_csum ' , n_w_csum_(nr)
     $           , ' n_w_csum_0in ' , n_w_csum_0in_(nr)
            stop !exit program due to error ;
         end if !if (n_w_csum_(nr).ne.n_w_csum_0in_(nr)) then
         nw = nw + n_w_(nr)
      enddo !do nr=0,n_r-1
      if (verbose.gt.1) then
         write(6,'(3(A,I0))') 
     $        '  n_w_max ',n_w_max
     $        ,'  n_w_sum ',n_w_sum
     $        ,'  n_gamma_z ',n_gamma_z
         call print_sub_i4(n_r,n_w_,' n_w_: ')
         call print_sub_i4(n_r,n_w_csum_,' n_w_csum_: ')
      end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  define quadrature weight across all points
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(weight_k_p_r_(0:1+n_r-1))
      call cs1_r8(n_r,weight_k_p_r_)
      allocate(weight_k_p_(0:1+n_w_sum-1))
      call cs1_r8(n_w_sum,weight_k_p_)
      na=0
      do nr=0,n_r-1
         weight_k_p_r_(nr) = kweight_k_p_r_(nr)/grid_k_p_r_(nr)
         do nw=0,n_w_(nr)-1
            weight_k_p_(na) = (2.0d0*pi*weight_k_p_r_(nr))/(1.0d0
     $           *n_w_(nr))
            na = na+1
         enddo !do nw=0,n_w_(nr)-1
      enddo !do nr=0,n_r-1
      if (na.ne.n_w_sum) then
         write(6,'(A,I0)') ' Warning, na: ' , na
         stop !exit program due to error ;
      end if !if (na.ne.n_w_sum) then
      if (verbose.gt.0) then
         call print_sub_r8(n_r,grid_k_p_r_,' grid_k_p_r_: ')
         call print_sub_r8(n_r,kweight_k_p_r_,' kweight_k_p_r_: ')
         call print_sub_r8(n_r,weight_k_p_r_,' weight_k_p_r_: ')
         call print_sub_r8(n_w_sum,weight_k_p_,' weight_k_p_: ')
      end if !if (verbose.gt.0) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  reorganize displacement arrays
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(delta_x_(0:1+n_delta_v-1))
      call cs1_r8(n_delta_v,delta_x_)
      allocate(delta_y_(0:1+n_delta_v-1))
      call cs1_r8(n_delta_v,delta_y_)
      do ndv=0,n_delta_v-1
         delta_out_(0) = +delta__(0+ndv*2)/(2.0d0*pi)
         delta_out_(1) = +delta__(1+ndv*2)/(2.0d0*pi)
         delta_x_(ndv) = delta_out_(0)
         delta_y_(ndv) = delta_out_(1)
      enddo !do ndv=0,n_delta_v-1
      if (verbose.gt.-0) then
         call print_sub_r8(n_delta_v,delta_x_,' delta_x_: ')
         call print_sub_r8(n_delta_v,delta_y_,' delta_y_: ')
      end if !if (verbose.gt.0) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  allocate index arrays for S and M
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_S_sample_(0:1+n_S-1))
      call cs1_i4(n_S,I_S_sample_)
      do ns=0,n_S-1
         I_S_sample_(ns)=1.0d0*ns
      enddo !do ns=0,n_S-1
      allocate(S_alpha_S_index_(0:1+n_S-1))
      call cs1_r8(n_S,S_alpha_S_index_)
      do ns=0,n_S-1
         S_alpha_S_index_(ns) = 1.0d0*ns ! array storing the numbers 0:n_S-1, from 0 up to the number of templates used - 1. ;
      enddo
      allocate(I_M_sample_(0:1+n_M-1))
      call cs1_i4(n_M,I_M_sample_)
      do nm=0,n_M-1
         I_M_sample_(nm)=1.0d0*nm
      enddo !do nm=0,n_M-1
      allocate(C_M_(0:1+n_M-1))
      call cs1_c16(n_M,C_M_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  collect alpha_est__
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(alpha_est__(0:1+n_alpha*n_M-1))
      call cs1_r8(n_alpha*n_M,alpha_est__)
      do nm=0,n_M-1
         alpha_est__(nalpha_polar_a + nm*n_alpha) = 
     $        0.0d0*alpha_first5__(nalpha_polar_a + nm*n_alpha_first5) ! restore for local search. ;
         alpha_est__(nalpha_azimu_b + nm*n_alpha) = 
     $        0.0d0*alpha_first5__(nalpha_azimu_b + nm*n_alpha_first5) ! restore for local search. ;
         alpha_est__(nalpha_gamma_z + nm*n_alpha) = 
     $        0.0d0*alpha_first5__(nalpha_gamma_z + nm*n_alpha_first5) ! restore for successive updates. ;
         alpha_est__(nalpha_delta_x + nm*n_alpha) = 
     $        0.0d0*alpha_first5__(nalpha_delta_x + nm*n_alpha_first5) ! restore for successive updates. ;
         alpha_est__(nalpha_delta_y + nm*n_alpha) = 
     $        0.0d0*alpha_first5__(nalpha_delta_y + nm*n_alpha_first5) ! restore for successive updates. ;
         alpha_est__(nalpha_l2_norm + nm*n_alpha) = 
     $        1.0d0
         alpha_est__(nalpha_ctf_ind + nm*n_alpha) = 
     $        CTF_1pind_(nm)-1 ! shift to base-0 indexing. ;
         alpha_est__(nalpha_S_index + nm*n_alpha) = 
     $        0.0d0
         alpha_est__(nalpha_M_index + nm*n_alpha) = 
     $        1.0d0*nm
         alpha_est__(nalpha_CTF_R_S + nm*n_alpha) = 
     $        0.0d0
         alpha_est__(nalpha_C_Z_opt + nm*n_alpha) = 
     $        0.0d0
      enddo !do nm=0,n_M-1
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
      allocate(S_k_p_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,S_k_p_)
      allocate(CTF_k_p_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,CTF_k_p_)
      allocate(R_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,R_S_)
      allocate(T_R_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,T_R_S_)
      allocate(T_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,T_S_)
      allocate(R_T_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,R_T_S_)
      allocate(CTF_R_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,CTF_R_S_)
      allocate(CTF_T_R_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,CTF_T_R_S_)
      allocate(CTF_R_T_S_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,CTF_R_T_S_)
      allocate(M_k_p_(0:1+n_w_sum-1))
      call cs1_c16(n_w_sum,M_k_p_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'ti8_bestfitall_excerpt_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A,I0)') 'get (n_M = '
     $     ,n_M_sample,')-x-(n_S_sample = ',n_S_sample, ') = (',
     $     n_M_sample *n_S_sample,') innerproducts of length ' ,
     $     n_w_csum_(n_k_p_r_cur -1) + n_w_(n_k_p_r_cur-1)
      timing_tic = omp_get_wtime()
      if (verbose.gt.1) then
         write(6,'(A)') ' Estimating memory requirements: '
      end if !if (verbose.gt.1) then
      n_S_0_sub_max = ceiling(dsqrt(1.0d0*n_S_sample))
      n_S_1_sub_max = ceiling(dsqrt(1.0d0*n_S_sample))
      n_M_0_sub_max = ceiling(dsqrt(1.0d0*n_M_sample))
      n_M_1_sub_max = ceiling(dsqrt(1.0d0*n_M_sample))
      nS_0_sub_tmp = 1
      nS_1_sub_tmp = 1
      nM_0_sub_tmp = 1
      nM_1_sub_tmp = 1
      flag_continue = .true.
      flag_proceed = .true.
      do while (flag_continue.eqv..true.)
         flag_memory_estimate = .true.
         d_memory_estimate = 0.0d0
         if (verbose.gt.1) write(6,'(A)') ' Testing memory: '
         include 'ti8_bestfitall_excerpt_checkset.f'
      call ti8(
     $     verbose-2 !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_r_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $     ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,weight_k_p_r_ !real *8 array (length at least n_k_p_r_cur): weights for radial quadrature. ;
     $     ,weight_k_p_ !real *8 array (length at least n_A_cur): all weights (unrolled) for quadrature, for use with k_p coordinates. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S_sample !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S_sample): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_ !real *8 array (length at least n_S_sample): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S_sample): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S_sample): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $     ,I_M_sample_ !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $     ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $     ,n_CTF !integer *4: total number of CTF functions. ;
     $     ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
     $     ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $     ,alpha_est__ !real *8 array (length at least n_alpha*n_M_sample): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $     ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $     ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
     $     ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_SM_ !integer *4 array (length at least n_M_sample): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
     $     ,alpha_SM__ !real *8 array (length at least n_alpha*n_SM_max*n_M_sample): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
     $     ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_MS_ !integer *4 array (length at least n_S_sample): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
     $     ,alpha_MS__ !real *8 array (length at least n_alpha*n_MS_max*n_S_sample): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
     $     ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_p_r_cur' under consideration (which can change from iteration to iteration). ;
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
     $     ,nS_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,nS_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,nM_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,nM_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M_sample): actual l2-norm (out to frequency n_k_p_r_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
        if (d_memory_estimate.le.d_memory_limit) then
           if (verbose.gt.-1) then
              write(6,'(A,4(A,I0))') ' Memory limit not exceeded: '
     $             , ' nS_0_sub_tmp ' ,nS_0_sub_tmp , ' nS_1_sub_tmp ' ,
     $             nS_1_sub_tmp ,' nM_0_sub_tmp ' , nM_0_sub_tmp ,
     $             ' nM_1_sub_tmp ' , nM_1_sub_tmp 
              write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $             nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $             nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
           end if !if (verbose.gt.1) then
           flag_continue = .false.
           flag_proceed = .true.
        end if !if (d_memory_estimate.le.d_memory_limit) then
        if (d_memory_estimate.gt.d_memory_limit) then
           flag_continue = .true.
           if ((nS_0_sub_tmp.ge.n_S_0_sub_max) .and.
     $          (nS_1_sub_tmp.ge.n_S_1_sub_max) .and.
     $          (nM_0_sub_tmp.ge.n_M_0_sub_max) .and.
     $          (nM_1_sub_tmp.ge.n_M_1_sub_max)) then
              write(6,'(A,4(A,I0))') ' Warning, memory limit exceeded: '
     $             , ' nS_0_sub_tmp ' ,nS_0_sub_tmp , ' nS_1_sub_tmp ' ,
     $             nS_1_sub_tmp ,' nM_0_sub_tmp ' , nM_0_sub_tmp ,
     $             ' nM_1_sub_tmp' , nM_1_sub_tmp 
              write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $             nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $             nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
              flag_continue = .false.
              flag_proceed = .false.
           else
              if (.false.) then
c$$$             do nothing. ;
              else if (nM_0_sub_tmp.lt.n_M_0_sub_max) then
                 nM_0_sub_tmp = nM_0_sub_tmp + 1
              else if (nM_1_sub_tmp.lt.n_M_1_sub_max) then
                 nM_1_sub_tmp = nM_1_sub_tmp + 1
              else if (nS_0_sub_tmp.lt.n_S_0_sub_max) then
                 nS_0_sub_tmp = nS_0_sub_tmp + 1
              else if (nS_1_sub_tmp.lt.n_S_1_sub_max) then
                 nS_1_sub_tmp = nS_1_sub_tmp + 1
              end if ! blocks small. ;
           end if ! blocks too big. ;
        end if !if (d_memory_estimate.le.d_memory_limit) then
      enddo !do while (flag_continue.eqv..true.)
      if (flag_proceed.eqv..true.) then
      flag_memory_estimate = .false.
      if (verbose.gt.1) write(6,'(A)') ' call ti8: '
      include 'ti8_bestfitall_excerpt_checkset.f'
      if (verbose.gt.1) then
         write(6,'(A)') ' first and last images: '
         call print_sub_c16(n_w_sum,M_k_p__(0*ld_M),' first image: ')
         call print_sub_c16(n_w_sum,M_k_p__((n_M-1)*ld_M)
     $        ,'  last image: ')
         write(6,'(A)') ' first and last templates: '
         call print_sub_c16(n_w_sum,S_k_p__(0*ld_S),' first template: ')
         call print_sub_c16(n_w_sum,S_k_p__((n_S-1)*ld_S)
     $        ,'  last template: ')
      end if !if (verbose.gt.1) then
      call ti8(
     $     verbose !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_r_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $     ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,weight_k_p_r_ !real *8 array (length at least n_k_p_r_cur): weights for radial quadrature. ;
     $     ,weight_k_p_ !real *8 array (length at least n_A_cur): all weights (unrolled) for quadrature, for use with k_p coordinates. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S_sample !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S_sample): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_ !real *8 array (length at least n_S_sample): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S_sample): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S_sample): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $     ,I_M_sample_ !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $     ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $     ,n_CTF !integer *4: total number of CTF functions. ;
     $     ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
     $     ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $     ,alpha_est__ !real *8 array (length at least n_alpha*n_M_sample): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $     ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $     ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
     $     ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_SM_ !integer *4 array (length at least n_M_sample): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
     $     ,alpha_SM__ !real *8 array (length at least n_alpha*n_SM_max*n_M_sample): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
     $     ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_MS_ !integer *4 array (length at least n_S_sample): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
     $     ,alpha_MS__ !real *8 array (length at least n_alpha*n_MS_max*n_S_sample): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
     $     ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_p_r_cur' under consideration (which can change from iteration to iteration). ;
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
     $     ,nS_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,nS_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,nM_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,nM_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M_sample): actual l2-norm (out to frequency n_k_p_r_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
      timing_toc = omp_get_wtime()
      timing_innerproduct = timing_toc-timing_tic
      if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
         write(6,'(A,A,F8.3)') 'ti8:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
      end if !if (flag_proceed.eqv..true.) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'ti8_bestfitall_excerpt_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  update alpha_est__
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..false.) then 
         if (verbose.gt.1) then
         write(6,'(A)') ' alpha_SM: '
         do nm=0,n_M_sample-1
            write(prefix_string,'(A,I2,A)') ' alpha_SM: nm:', nm , ' '
            call alpha_SM_write_0(n_SM_max,n_SM_(nm),alpha_SM__(n_alpha
     $           *n_SM_max*nm),len_trim(prefix_string)+1,prefix_string)
         enddo !do nm=0,n_M_sample-1
         end if !if (verbose.gt.1) then
         do nm=0,n_M_sample-1
            call test_alpha_update_SM_3(
     $           verbose
     $           ,rseed
     $           ,n_SM_max
     $           ,n_SM_(nm)
     $           ,alpha_SM__(n_alpha*n_SM_max*nm)
     $           ,flag_RTRT_vs_RTTR
     $           ,alpha_est__(n_alpha*nm)
     $           ,alpha_update_f
     $           ,alpha_est__(n_alpha*nm)
     $           )
         enddo !do nm=0,n_M_sample-1
      end if !if (flag_MS_vs_SM.eqv..false.) then 
      if (flag_MS_vs_SM.eqv..true.) then 
         if (verbose.gt.1) then
         write(6,'(A)') ' alpha_MS: '
         do ns=0,n_S_sample-1
            write(prefix_string,'(A,I2,A)') ' alpha_MS: ns:', ns , ' '
            call alpha_SM_write_0(n_MS_max,n_MS_(ns),alpha_MS__(n_alpha
     $           *n_MS_max*ns),len_trim(prefix_string)+1,prefix_string)
         enddo !do ns=0,n_S_sample-1
         end if !if (verbose.gt.1) then
         call test_alpha_update_MS_3(
     $        verbose
     $        ,rseed
     $        ,n_MS_max
     $        ,n_MS_
     $        ,n_S_sample
     $        ,n_M_sample
     $        ,alpha_MS__
     $        ,flag_RTRT_vs_RTTR
     $        ,alpha_est__
     $        ,alpha_est__
     $        )
      end if !if (flag_MS_vs_SM.eqv..true.) then 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  collect alpha_first5__. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.1) then
      write(6,'(A)') ' alpha_est__ for first 10 images: '
      do nm=0,min(10,n_M-1)
         call print_all_r8(n_alpha,alpha_est__(nm*n_alpha)
     $        ,'alpha_est__: ')
      enddo !do nm=0,min(5,n_M-1)
      end if !if (verbose.gt.1) then
      do nm=0,n_M-1
         gamma_z = alpha_est__(nalpha_gamma_z + nm*n_alpha)
         delta_out_(0) = alpha_est__(nalpha_delta_x + nm
     $        *n_alpha)
         delta_out_(1) = alpha_est__(nalpha_delta_y + nm
     $        *n_alpha)
         if (flag_RTRT_vs_RTTR.eqv..false.) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Rotating delta_upd by +gamma_upd. '
         end if !if (verbose.gt.2) then
         delta_0in_(0) = delta_out_(0)
         delta_0in_(1) = delta_out_(1)
         call rotate_delta(+gamma_z,delta_0in_,delta_out_)
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
         alpha_first5__(nalpha_polar_a + nm*n_alpha_first5) = 
     $        alpha_est__(nalpha_polar_a + nm*n_alpha)
         alpha_first5__(nalpha_azimu_b + nm*n_alpha_first5) = 
     $        alpha_est__(nalpha_azimu_b + nm*n_alpha)
         alpha_first5__(nalpha_gamma_z + nm*n_alpha_first5) = 
     $        2.0d0*pi - gamma_z
         call periodize_r8(
     $        alpha_first5__(nalpha_gamma_z + nm*n_alpha_first5)
     $        ,0.0d0
     $        ,2.0d0*pi
     $        ,alpha_first5__(nalpha_gamma_z + nm*n_alpha_first5)
     $        )
         alpha_first5__(nalpha_delta_x + nm*n_alpha_first5) = 
     $        -delta_out_(0)*(2.0d0*pi)
         alpha_first5__(nalpha_delta_y + nm*n_alpha_first5) = 
     $        -delta_out_(1)*(2.0d0*pi)
         if (verbose.gt.1) then
            write(6,'(3(A,I0),A,F24.16,3(A,F24.16))')
     $           ' checking image: ' , nm
     $           ,' ctf_index: ' , CTF_1pind_(nm)-1
     $           ,' S_index: '
     $           ,nint(alpha_est__(nalpha_S_index+ nm*n_alpha))
     $           ,' C_M: '
     $           ,dreal(C_M_(nm))
     $           ,' gamma_z: ' , gamma_z
     $           ,' delta_x: ' , delta_out_(0)
     $           ,' delta_y: ' , delta_out_(1)
            call cp1_c16(n_w_sum,M_k_p__(nm*ld_M),M_k_p_)
            call cp1_c16(n_w_sum,CTF_k_p__((CTF_1pind_(nm)-1)*ld_CTF)
     $           ,CTF_k_p_)
            call cp1_c16(n_w_sum,S_k_p__(nint(alpha_est__(nalpha_S_index
     $           + nm*n_alpha))*ld_S),S_k_p_)
            call rotate_p_to_p_fftw_plan_0on(
     $           n_r
     $           ,n_w_
     $           ,n_w_sum
     $           ,S_k_p_
     $           ,+gamma_z
     $           ,R_S_
     $           )
            call xc1_c16(n_w_sum,R_S_,CTF_k_p_,CTF_R_S_)
            call transf_p_to_p(
     $           n_r
     $           ,grid_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,R_S_
     $           ,+delta_out_(0)
     $           ,+delta_out_(1)
     $           ,T_R_S_
     $           )
            call xc1_c16(n_w_sum,T_R_S_,CTF_k_p_,CTF_T_R_S_)
            call innerproduct_p_quad(
     $           n_r
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,CTF_R_S_
     $           ,CTF_R_S_
     $           ,CTF_R_S_norm
     $           )
            CTF_R_S_norm = zsqrt(CTF_R_S_norm)/sqrt(2.0d0)
            call innerproduct_p_quad(
     $           n_r
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,M_k_p_
     $           ,M_k_p_
     $           ,M_norm
     $           )
            M_norm = zsqrt(M_norm)/sqrt(2.0d0)
            call innerproduct_p_quad(
     $           n_r
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,CTF_T_R_S_
     $           ,M_k_p_
     $           ,S_R_T_CTF_M
     $           )
            S_R_T_CTF_M = S_R_T_CTF_M/2.0d0
            C_Z_opt_RTTR = S_R_T_CTF_M/(CTF_R_S_norm*M_norm)
            l2_norm_RTTR = S_R_T_CTF_M/CTF_R_S_norm**2
            call transf_p_to_p(
     $           n_r
     $           ,grid_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,S_k_p_
     $           ,+delta_out_(0)
     $           ,+delta_out_(1)
     $           ,T_S_
     $           )
            call rotate_p_to_p_fftw_plan_0on(
     $           n_r
     $           ,n_w_
     $           ,n_w_sum
     $           ,T_S_
     $           ,+gamma_z
     $           ,R_T_S_
     $           )
            call xc1_c16(n_w_sum,R_T_S_,CTF_k_p_,CTF_R_T_S_)
            call innerproduct_p_quad(
     $           n_r
     $           ,grid_k_p_r_
     $           ,weight_k_p_r_
     $           ,n_w_
     $           ,n_w_sum
     $           ,CTF_R_T_S_
     $           ,M_k_p_
     $           ,S_T_R_CTF_M
     $           )
            S_T_R_CTF_M = S_T_R_CTF_M/2.0d0
            C_Z_opt_RTRT = S_T_R_CTF_M/(CTF_R_S_norm*M_norm)
            l2_norm_RTRT = S_T_R_CTF_M/CTF_R_S_norm**2
            write(6,'(2(A,F24.18))')
     $           ' CTF_R_S_norm: ' , dreal(CTF_R_S_norm)
     $           ,' M_norm: ' , dreal(M_norm)
            write(6,'(3(A,F24.18))')
     $           ' S_T_R_CTF_M: ' , dreal(S_T_R_CTF_M)
     $           ,' l2_norm_RTRT: ' , dreal(l2_norm_RTRT)
     $           ,' C_Z_opt_RTRT: ' , dreal(C_Z_opt_RTRT)
            write(6,'(3(A,F24.18))')
     $           ' S_R_T_CTF_M: ' , dreal(S_R_T_CTF_M)
     $           ,' l2_norm_RTTR: ' , dreal(l2_norm_RTTR)
     $           ,' C_Z_opt_RTTR: ' , dreal(C_Z_opt_RTTR)
            if (verbose.gt.2) then
            write(6,'(3(A,F24.16))')
     $           ' gamma_z: ' , gamma_z
     $           ,' delta_x: ' , delta_out_(0)
     $           ,' delta_y: ' , delta_out_(1)
            call print_sub_c16(n_w_sum,T_R_S_,' T_R_S_: ')
            call print_sub_c16(n_w_sum,R_T_S_,' R_T_S_: ')
            call print_sub_c16(n_w_sum,CTF_T_R_S_,' CTF_T_R_S_: ')
            call print_sub_c16(n_w_sum,CTF_R_T_S_,' CTF_R_T_S_: ')
            end if !if (verbose.gt.2) then
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'ti8_bestfitall_excerpt_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  deallocating. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..false.) then
         deallocate(n_SM_)
         deallocate(alpha_SM__)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         deallocate(n_MS_)
         deallocate(alpha_MS__)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      deallocate(C_M_)
      deallocate(alpha_est__)
      deallocate(I_M_sample_)
      deallocate(S_alpha_S_index_)
      deallocate(I_S_sample_)
      deallocate(delta_x_)
      deallocate(delta_y_)
      deallocate(weight_k_p_)
      deallocate(n_w_csum_)
      deallocate(n_polar_a_)

      if (verbose.gt.0) then
         write(6,'(A)') ' [finished ti8_bestfitall] '
      end if !if (verbose.gt.0) then

      end !subroutine
