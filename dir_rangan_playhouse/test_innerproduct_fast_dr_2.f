c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program test_innerproduct_fast_dr_2
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This program is a driver which tests the subroutine
c$$$      test_innerproduct_fast_wrapper_1.f ,
c$$$      which, collectively, calculate innerproducts for a batch of
c$$$      images and templates. These routines account for a collection
c$$$      of ctf-functions, as well as an (unknown) image-normalization 
c$$$      factor.
c$$$      
c$$$      The main steps of this program are:
c$$$      A. First we generate a set of n_S = 5 templates.
c$$$         Each template looks a little different from the next
c$$$      B. Then we generate a set of n_M = 15 images.
c$$$         The first three images are copies of the first template,
c$$$         but with a different alignment (i.e., translated and 
c$$$         rotated slightly).
c$$$         The second three images are copies of the second template,
c$$$         again with a different alignment, and so forth.
c$$$      C. Then we generate n_CTF = 3 ctf-functions.
c$$$         These are deliberately chosen to be very anisotropic.
c$$$         These ctf-functions are assigned arbitrarily to the
c$$$         various images.
c$$$      D. Then we call test_innerproduct_fast_wrapper_1 to calculate
c$$$         the innerproducts between each pair of templates and
c$$$         images. These innerproducts are calculated across a range
c$$$         of displacements delta_x_,delta_y_ as well as a range of 
c$$$         rotations gamma_z_ (all passed in as input).
c$$$      D. The routine test_innerproduct_fast_wrapper_1 is used
c$$$         to determine (among other things) which template correlates
c$$$         most strongly with each image, and what the best
c$$$         alignment (i.e., displacement and rotation) is for each image.
c$$$      E. The true alignment of each image is then compared with
c$$$         the best empirically observed alignment, and the alignment-
c$$$         errors are printed to stdout.
c$$$      F. Finally, if (verbose.gt.1), we generate two figure files:
c$$$         ./dir_jpg/Fig_NMC_[ncur].jpg:
c$$$              R(N_k_p_(nm)) and I(N_k_p_(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid.
c$$$              R(C_k_p_(nm)) and I(C_k_p_(nm)) show real and imaginary
c$$$              parts of ctf-function I_ctf_(nm) on k-space quasi-uniform 
c$$$              polar-grid.
c$$$              R(M_k_p_(nm)) and I(M_k_p_(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid
c$$$              after pointwise multiplication by ctf-function and
c$$$              multiplication by image_specific scaling factor l2_norm.
c$$$         ./dir_jpg/Fig_MTS_[ncur].jpg:
c$$$              Each pair of rows shows the real and imaginary parts of
c$$$              some function on k-space quasi-uniform polar-grid.
c$$$              N_k_p_ is the primitive image
c$$$              M_k_p_ is the image multiplied by ctf and scaled
c$$$              S_k_p_ is the associated best-fit template
c$$$              T_k_p_ is the associated best-fit template aligned to
c$$$                     fit the image (i.e,. translated and rotated).
c$$$              CT_k_p_ is the associated best-fit and aligned template
c$$$                     multiplied by the image-specific ctf-function
c$$$
c$$$      Note regarding alignment:
c$$$      Assume that we are given an image M and template S.
c$$$      The innerproduct between M and S (with no alignment shift) is
c$$$      \[ \int M \cdot \bar{S} \].
c$$$      If we were to consider the alignment shift associated with
c$$$      a displacement (delta_x,delta_y) and an (in-plane) rotation
c$$$      gamma_z, then the innerproduct would be:
c$$$      \[ \int M \cdot \bar{T(S)} \], 
c$$$      where T(S) is the transformed template obtained by *first*
c$$$      translating S by (+delta_x,+delta_y), and *then* rotating S
c$$$      by (+gamma_z).
c$$$      This is equivalent to the innerproduct
c$$$      \[ \int T(M) \cdot \bar{S} \],
c$$$      where T(M) is the transformed image obtained by *first*
c$$$      rotating M by (-gamma_z) and *then* translating M by
c$$$      (-delta_x,-delta_y).
c$$$
c$$$      Note regarding naming-conventions:
c$$$      We use the following naming conventions for templates and images:
c$$$      S_x_c_: template in real-space on a cartesian-grid.
c$$$      S_x_p_: template in real-space on a quasi-uniform polar-grid.
c$$$      S_k_c_: template in fourier-space on a cartesian-grid.
c$$$      S_k_p_: template in fourier-space on a quasi-uniform polar-grid.
c$$$      S_k_q_: template in fourier-space, represented using "hat-hat" 
c$$$              coordinates (i.e., bessel-expansion).
c$$$      N_x_c_ through N_k_q_ are defined similarly for the primitive 
c$$$             images (i.e., before applying ctf or scaling).
c$$$      M_x_c_ through M_k_q_ are defined similarly for the images
c$$$             after multiplication by ctf-function.
c$$$      grid_x_c_: grid-values of real-space cartesian-grid
c$$$      grid_k_c_: grid-values of fourier-space cartesian-grid
c$$$      grid_k_p_: radial-values of rings of fourier-space polar-grid
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$
c$$$      INPUTS: This driver depends on the parameters defined below:
c$$$
c$$$      integer verbose       Determines verbosity level.
c$$$                          - verbose.eq.1 prints out only basic 
c$$$                            information and generates no figures.
c$$$                          - verbose.eq.2 prints out pre- and post-sorted
c$$$                            alignment arrays and generates figures.
c$$$
c$$$      integer ncur          Number of rings to pass into 
c$$$                            get_template_size. This parameter is also
c$$$                            used to determine the number of gridpoints 
c$$$                            (i.e., pixels) on each side of the templates
c$$$                            and images in real space (as measured on a
c$$$                            cartesian-grid). 
c$$$                          - Values of ncur~10 --> low resolution. 
c$$$                          - Values of ncur~50 --> moderate resolution.
c$$$                          - Values of ncur~100 --> high resolution.
c$$$                          - We have not generated libraries for values
c$$$                            of ncur.gt.192 yet.
c$$$
c$$$      real *8 eps_svd       The epsilon used as a cutoff when 
c$$$                            determining which terms to retain within the
c$$$                            svd-expansion of the various bessel-
c$$$                            functions used to approximate the term
c$$$                            exp(2*pi*K*delta*cos(psi-omega)) within
c$$$                            the fourier-space translation operator.
c$$$                          - A value of eps_svd.eq.10.0d0 corresponds
c$$$                            to low accuracy, and will probably not give
c$$$                            exactly the same result as an exact brute-
c$$$                            force calculation 
c$$$                            (see test_innerproduct_svdread_dr.f for a 
c$$$                            direct comparison).
c$$$                          - A value of eps_svd.eq.1.0d0 corresponds
c$$$                            to moderate accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to a few digits, which is 
c$$$                            usually enough for the best alignment to 
c$$$                            match exactly.
c$$$                          - A value of eps_svd.eq.0.1d0 corresponds
c$$$                            to high accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to several digits, which is 
c$$$                            more than enough to match the best alignment
c$$$                            exactly.
c$$$                          - Values of eps_svd.lt.0.1d0 are not
c$$$                            usually necessary, and we have not generated
c$$$                            the appropriate libraries for these 
c$$$                            eps_svd yet.
c$$$                            
c$$$      real *8 N_pixels_in   Number of pixels in each direction 
c$$$                            associated with the largest displacement in 
c$$$                            the array of displacements used when 
c$$$                            calculating innerproducts for
c$$$                            each pair of templates and images.
c$$$                          - Each pixel of displacement is measured in
c$$$                            x-space. Thus, the maximum displacement is:
c$$$                            delta_x = N_pixels_in/n_r * max_x_c
c$$$                            delta_y = N_pixels_in/n_r * max_x_c.
c$$$                          - For example, N_pixels_in.eq.1.0d0 will limit
c$$$                            the array of displacements to a single pixel
c$$$                            in each direction (i.e., a 3x3 grid of 
c$$$                            pixels centered at the origin).
c$$$                          - A value of N_pixels_in.eq.3.0d0 will 
c$$$                            correspond to an array of displacements
c$$$                            spanning up to 3 pixels in each direction
c$$$                            (i.e., a 7x7 grid of pixels centered at the
c$$$                            origin).
c$$$                          - Values of N_pixels_in.gt.3.0 are not usually
c$$$                            necessary, and we have not generated the
c$$$                            appropriate libraries for these N_pixels_in
c$$$                            yet.
c$$$
c$$$      integer n_delta_x     Number of displacements to search for in the 
c$$$                            x-direction. These displacements will be 
c$$$                            linearly spaced from -N_pixels_in to 
c$$$                            +N_pixels_in.
c$$$                          - A value of n_delta_x = 1 + 2*N_pixels_in
c$$$                            implies a displacement-vector associated 
c$$$                            with each pixel-center in the grid 
c$$$                            determined by N_pixels_in.
c$$$                          - A value of n_delta_x = 1 + 4*N_pixels_in
c$$$                            implies a displacement-vector associated 
c$$$                            with each half-pixel (i.e., pixel-center
c$$$                            and pixel-side) in the grid determined by
c$$$                            N_pixels_in.
c$$$                          - Note that the total number of displacement-
c$$$                            vectors is n_delta_x * n_delta_y.
c$$$
c$$$      integer n_delta_y     Number of displacements to search for in the 
c$$$                            y-direction. These displacements will be 
c$$$                            linearly spaced from -N_pixels_in to 
c$$$                            +N_pixels_in.
c$$$                          - Note that the total number of displacement-
c$$$                            vectors is n_delta_x * n_delta_y.
c$$$
c$$$      integer n_gamma_z     Number of rotations to search for. These
c$$$                            rotations will be linearly spaced from 0.0d0
c$$$                            to 2*pi. It is usually reasonable for
c$$$                            n_gamma_z to be on the same order as ncur
c$$$                            (e.g., n_gamma_z = ncur/2 or so).
c$$$
c$$$      integer svd_calculation_type Determines the type of calculation
c$$$                            to use when computing innerproducts.
c$$$                          - A value of svd_calculation_type.eq.1
c$$$                            forces the displacement operation to take
c$$$                            place after the fft. This calculation is
c$$$                            useful when the number of terms in the
c$$$                            svd-expansion is smaller than the number
c$$$                            of distinct displacements.
c$$$                          - A value of svd_calculation_type.eq.2
c$$$                            forces a brute-force calculation of each
c$$$                            displacement separately (using the fft
c$$$                            to handle rotations).
c$$$                          - A value of svd_calculation_type.lt.0
c$$$                            allows test_innerproduct_svdread_batch.f
c$$$                            to automatically choose one of the two
c$$$                            calculation_types based on the other
c$$$                            input parameters (e.g., by comparing
c$$$                            n_svd_l with n_delta_x*n_delta_y).
c$$$
c$$$      integer n_omp_sub     Number of sub-blocks for omp parallelization.
c$$$                            Images are batched into groups of roughly
c$$$                            equal size.
c$$$
c$$$      logical flag_ctf      Whether or not to use anisotropic ctf.
c$$$
c$$$      logical flag_tru      Whether or not to set estimated image parameters
c$$$                            to be equal to the true image parameters.
c$$$
c$$$      integer n_S_max       Maximum number of templates.
c$$$
c$$$      integer n_M_max       Maximum number of images.
c$$$
c$$$      integer fpm_howmany_max Maximum number of fftws to call simultaneously
c$$$                            within the fftw_plan_many. 
c$$$
c$$$      logical flag_MS_vs_SM determines whether to assign images to templates (.true.) 
c$$$                            or templates to images (.false.). ;
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_omp_sub
      logical flag_ctf
      logical flag_tru
      logical flag_fig
      logical flag_MS_vs_SM
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: n_polar_a_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na
      integer n_polar_a,npolar_a,n_azimu_b,n_polar_a_azimu_b_sum
      integer *4, allocatable :: n_w_(:)
c$$$      fftw
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*)
      complex *16 fftw_out_last_(*)
c$$$      grids
      real *8 pi
      real *8 max_r8_f,avg_r8_f,al2_r8_f
      real *8 max_x_c,max_k_c
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
c$$$      Templates S_
      integer n_S_max,n_S,ld_S,ns
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: T_x_c_(:)
      complex *16, allocatable :: S_x_p_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: S_k_q_(:)
c$$$      Images N_ and M_
      integer n_M_max,n_M,ld_M,nm
      complex *16, allocatable :: N_x_c_(:)
      complex *16, allocatable :: N_x_p_(:)
      complex *16, allocatable :: N_k_c_(:)
      complex *16, allocatable :: N_k_p_(:)
      complex *16, allocatable :: N_k_q_(:)
      complex *16, allocatable :: M_k_p_(:)
c$$$      CTF-functions
      integer *4 n_ctf,ld_CTF,nctf
      integer *4, allocatable :: I_ctf_(:)
      complex *16, allocatable :: CTF_k_p_(:)
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4
c$$$      fftw_plan_many 
      integer fpm_howmany_max
c$$$      image parameters 
      integer *4 n_alpha
      parameter (n_alpha=8)
      include 'nalpha_define.f'
      real *8 l2_norm
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: alpha_polar_a_(:)
      real *8, allocatable :: alpha_azimu_b_(:)
      real *8 alpha_update_f
      parameter (alpha_update_f=0.0d0)
c$$$      true displacement, rotation and scaling parameters for each image
      real *8, allocatable :: azimu_b_tru_(:)
      real *8, allocatable :: delta_x_tru_(:)
      real *8, allocatable :: delta_y_tru_(:)
      real *8, allocatable :: gamma_z_tru_(:)
      real *8, allocatable :: l2_norm_tru_(:)
c$$$      estimated (current) displacement, rotation and scaling parameters for each image
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: l2_norm_est_(:)
c$$$      updated displacement, rotation and scaling parameters for each image
      real *8, allocatable :: azimu_b_upd_(:)
      real *8, allocatable :: delta_x_upd_(:)
      real *8, allocatable :: delta_y_upd_(:)
      real *8, allocatable :: gamma_z_upd_(:)
      real *8, allocatable :: l2_norm_upd_(:)
c$$$      errors in displacement, rotation and scaling parameters for each image
      real *8, allocatable :: azimu_b_err_(:)
      real *8, allocatable :: delta_x_err_(:)
      real *8, allocatable :: delta_y_err_(:)
      real *8, allocatable :: gamma_z_err_(:)
      real *8, allocatable :: l2_norm_rer_(:)
c$$$      range of displacement and rotation parameters to measure
      real *8 displacement_max,angle_tmp,delta_tmp,dtmp
      integer n_delta_x,n_delta_y,ndx,ndy
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_pre,delta_y_pre,gamma_z_pre
      real *8 delta_x_pos,delta_y_pos,gamma_z_pos
      real *8 eps_svd,N_pixels_in
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      integer n_gamma_z,ngz
      real *8, allocatable :: gamma_z_(:)      
      integer svd_calculation_type
c$$$      testing variables
      complex *16, allocatable :: Z_p_(:)
      integer ndx_optimal,ndy_optimal,ngz_optimal
      complex *16 C_Z_optimal
c$$$      tesselation parameters
      real *8 tesselation_distance_req
      parameter (tesselation_distance_req=2.0d0)
      integer n_SM_use
c$$$      parameters for figure generation
      character(len=1024) fname,format_string
c$$$      function for template generation
      external get_F2_x_c
      real *8 param_1,param_2
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      parameters 
c$$$      list: verbose ncur eps_svd N_pixels_in n_delta_x n_delta_y n_gamma_z svd_calculation_type n_omp_sub

      pi = 4*atan(1.0)

      verbose = 1
      ncur = 30
      eps_svd = 10.0d0
      N_pixels_in = 3.0d0
      n_delta_x = 7
      n_delta_y = 7
      svd_calculation_type = 1
      n_omp_sub = 8
      flag_ctf = .true.
      flag_tru = .false.
      n_S_max = 5
      n_M_max = 15
      fpm_howmany_max = 16
      flag_MS_vs_SM = .true.

      n_gamma_z = ncur
      n_polar_a = max(6,nint(pi*ncur))
      n_w_max = n_polar_a*2
      n_gamma_z = n_w_max
      fpm_howmany_max = max(1,fpm_howmany_max)
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_fast_dr_2]: '
       end if
       if (verbose.gt.0) then
         write(6,'(A,I0)') ' verbose: ',verbose
         write(6,'(A,I0)') ' ncur: ',ncur
         write(6,'(A,F6.3)') ' eps_svd: ',eps_svd
         write(6,'(A,F6.3)') ' N_pixels_in: ',N_pixels_in
         write(6,'(A,I0)') ' n_delta_x: ',n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ',n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') ' svd_calculation_type: '
     $        ,svd_calculation_type
         write(6,'(A,I0)') ' n_omp_sub: ' ,n_omp_sub
         write(6,'(A,I0)') ' n_S_max: ' ,n_S_max
         write(6,'(A,I0)') ' n_M_max: ' ,n_M_max
         write(6,'(A,I0)') ' fpm_howmany_max: ' ,fpm_howmany_max
         write(6,'(A,L1)') ' flag_ctf: ' ,flag_ctf
         write(6,'(A,L1)') ' flag_tru: ' ,flag_tru
      end if
      
      n_polar_a = max(6,nint(pi*ncur))
      n_w_max = n_polar_a*2
      n_polar_a_azimu_b_sum = 0
      do npolar_a=0,n_polar_a-1
         n_azimu_b = max(6,nint(2.0d0*n_polar_a*sin(pi*(npolar_a+0.5d0)
     $        /n_polar_a)))
         n_polar_a_azimu_b_sum = n_polar_a_azimu_b_sum + n_azimu_b
      enddo !do npolar_a=0,n_polar_a-1
      write(6,'(A,I0,A,I0,A,I0)') ' n_polar_a: ' , n_polar_a ,
     $     ' n_w_max: ' , n_w_max , ' n_polar_a_azimu_b_sum: '
     $     ,n_polar_a_azimu_b_sum
      
      if (verbose.gt.1) then
         flag_fig = .true.
      else
         flag_fig = .false.
      end if

c$$$      Calculating template size using 'get_template_size'
      allocate(n_polar_a_(0:ncur-1));
      do nr=0,ncur-1
c$$$         n_polar_a calculation similar to testrebuild.f:
         n_polar_a_(nr) = nint(pi*(nr+1))
      enddo
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') ' ncur = ',ncur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(n_polar_a_,ncur,ntemplatesize,ngridc_
     $     ,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' ntemplatesize = ',ntemplatesize
      end if
      
c$$$      indices
      n_r = ncur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') ' Error n_r',n_r,'<2'
      end if
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.0) then
         write(6,'(A,I0,A,I0)') '  n_w_max ',n_w_max,'; n_A ',n_A
      end if

c$$$      fftw
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
         write(6,'(A,A,F8.3)') ' fftw_plan:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))
      
c$$$      generating grids for templates and images 
      max_x_c = 1.0d0
      allocate(grid_x_c_(0:n_r-1))
      call linspace(0.0d0,max_x_c,n_r,grid_x_c_)
      allocate(grid_x_p_(0:n_r-1))
      dtmp = max_x_c/2.0/n_r
      call linspace(0.0d0+dtmp,max_x_c/2.0+dtmp,n_r,grid_x_p_)
      max_k_c = (1.0d0*n_r)/max_x_c
      allocate(grid_k_c_(0:n_r-1))
      call linspace(0.0d0,max_k_c,n_r,grid_k_c_)
      allocate(grid_k_p_(0:n_r-1))
      dtmp = max_k_c/2.0/n_r
      call linspace(0.0d0+dtmp,max_k_c/2.0+dtmp,n_r,grid_k_p_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      call get_delta_0(N_pixels_in,n_r,max_x_c,n_delta_x ,delta_x_
     $     ,n_delta_y,delta_y_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)')
     $        ' ncur or so.'
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.0) then
         call write_sub_r8(n_delta_x,delta_x_,11,' delta_x_: ')
         call write_sub_r8(n_delta_y,delta_y_,11,' delta_y_: ')
         call write_sub_r8(n_gamma_z,gamma_z_,11,' gamma_z_: ')
      end if !if (verbose.gt.1) then

      allocate(Z_p_(0:n_delta_x*n_delta_y*n_gamma_z-1))

c$$$      generate templates
      n_S = n_S_max
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') ' Here we construct the templates "S_".'
         write(6,'(A,I0,A)') ' We generate n_S = ',n_S
     $        ,' templates in total,'
         write(6,*) 'each using the function get_F2_x_c'
     $        ,' with different parameters.'
         write(6,*) 'For each template: '
         write(6,*)
     $        ' The array S_x_c_ holds the x-space values on a'
     $        ,' regular cartesian grid defined by grid_x_c_.'
         write(6,*)
     $        ' The array S_k_c_ holds the k-space values on a'
     $        ,' regular cartesian grid defined by grid_k_c_.'
         write(6,*)
     $        ' The array S_k_p_ holds the k-space values on a'
     $        ,' quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
         write(6,*) ' The array S_k_q_ holds the "hat-hat"'
     $        ,' (i.e, bessel-expansion) for the k-space values'
         write(6,*) ' on a quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
      end if !if (verbose.gt.1) then
      allocate(S_x_c_(0:n_r*n_r*n_S-1))
      allocate(T_x_c_(0:n_r*n_r))
      allocate(S_x_p_(0:n_A*n_S-1))
      allocate(S_k_c_(0:n_r*n_r*n_S-1))
      allocate(S_k_p_(0:n_A*n_S-1))
      allocate(S_k_q_(0:n_A*n_S-1))
      do ns=0,n_S-1
         param_1 = ns*dsqrt(2.5d0)
         call periodize_r8(param_1,-1.0d0,1.0d0,param_1)
         param_2 = ns+2
c$$$      Calculating x-space template on regular cartesian-grid
         call get_F2_x_c_(n_r,grid_x_c_,max_x_c,n_r ,grid_x_c_,max_x_c
     $        ,S_x_c_(ns*n_r*n_r),get_F2_x_c,param_1,param_2)
c$$$      Calculating x-space template on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_x_c,n_r,max_x_c ,S_x_c_(ns*n_r*n_r)
     $        ,n_r,grid_x_p_,n_w_,n_A,S_x_p_(ns*n_A))
c$$$      Calculating k-space template on regular cartesian-grid
         call adi_fft2(-1,n_r,n_r,S_x_c_(ns*n_r*n_r),S_k_c_(ns*n_r*n_r))
c$$$      Calculating k-space template on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_k_c,n_r,max_k_c ,S_k_c_(ns*n_r*n_r)
     $        ,n_r,grid_k_p_,n_w_,n_A,S_k_p_(ns*n_A))
c$$$      Calculating J-space template on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,S_k_p_(ns*n_A),S_k_q_(ns*n_A))
      enddo !do ns=0,n_S-1

      n_M = n_M_max
      if (verbose.gt.1) then
         write(6,'(A)') ' '
         write(6,'(A,I0,A)') ' Now we do the same thing for the n_M = '
     $        ,n_M,' images "M_".'
         write(6,*) 'The j^th image will be constructed by taking'
     $        ,' the (j/3)^th template,'
         write(6,*) 'translating it by delta_x_tru_(j) and'
     $        ,' delta_y_tru_(j) and then rotating'
     $        ,' it by gamma_z_tru_(j).'
         write(6,*) 'The template-index associated with each image'
     $        ,' will be stored in the azimu_b variable.'
         write(6,*) 'For this example we pick delta_x_tru_(j)'
     $        ,' to be (j%3)-1,',' and delta_y_tru_(j) to be'
     $        ,' 2*((j%3)-1).'
         write(6,*) 'We also pick gamma_z_tru_(j) to be'
     $        ,' 2*pi/3*(j%3).'
      end if !if (verbose.gt.1) then
      allocate(N_x_c_(0:n_r*n_r*n_M-1))
      allocate(N_k_c_(0:n_r*n_r*n_M-1))
      allocate(N_k_p_(0:n_A*n_M-1))
      allocate(N_k_q_(0:n_A*n_M-1))
c$$$      Calculating x-space image on regular cartesian-grid
      allocate(azimu_b_tru_(0:n_M-1))
      allocate(delta_x_tru_(0:n_M-1))
      allocate(delta_y_tru_(0:n_M-1))
      allocate(gamma_z_tru_(0:n_M-1))
      allocate(l2_norm_tru_(0:n_M-1))
      do nm=0,n_M-1
         ns = max(0,min(n_S-1,nm/3))
         ndx = max(0,min(n_delta_x-1,mod(nm+0,3)))
         ndy = max(0,min(n_delta_y-1,mod(nm+1,5)))
         ngz = max(0,min(n_gamma_z-1,mod(13*nm+2,n_gamma_z)))
         azimu_b_tru_(nm) = 1.0d0*ns
         delta_x_tru_(nm) = delta_x_(ndx)
         delta_y_tru_(nm) = delta_y_(ndy)
         gamma_z_tru_(nm) = gamma_z_(ngz)
         l2_norm_tru_(nm) = 1.0d0 + 1*nm
         if (verbose.gt.1) then
            write(6,'(I0,A,4F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_tru_(nm),delta_y_tru_(nm),gamma_z_tru_(nm)
     $           ,l2_norm_tru_(nm)
         end if !if (verbose.gt.1) then
c$$$         In the absence of interpolation errors,
c$$$         these transformations will produce the N_k_p_ used below.
         call cp1_c16(n_r*n_r,S_x_c_(ns*n_r*n_r),N_x_c_(nm*n_r*n_r))
         call transl_c_to_c(n_r,max_x_c,n_r,max_x_c,N_x_c_(nm*n_r*n_r),
     $        +delta_x_tru_(nm),+delta_y_tru_(nm),N_x_c_(nm*n_r*n_r))
         call rotate_c_to_c(n_r,max_x_c,n_r,max_x_c,N_x_c_(nm*n_r*n_r),
     $        +gamma_z_tru_(nm),N_x_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on regular cartesian-grid      
         call adi_fft2(-1,n_r,n_r,N_x_c_(nm*n_r*n_r),N_k_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_k_c,n_r,max_k_c,N_k_c_(nm*n_r*n_r)
     $        ,n_r,grid_k_p_,n_w_,n_A,N_k_p_(nm*n_A))
c$$$      Calculating J-space image on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,N_k_p_(nm*n_A),N_k_q_(nm*n_A))
c$$$         However, due to interpolation errors, we will generate
c$$$         the N_k_p_ in k-space.
         call cp1_c16(n_A,S_k_p_(ns*n_A),N_k_p_(nm*n_A))
         call transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,N_k_p_(nm*n_A),
     $        +delta_x_tru_(nm),+delta_y_tru_(nm),N_k_p_(nm*n_A))
         call rotate_p2p_fz(n_r,n_w_,n_A,N_k_p_(nm*n_A),
     $        +gamma_z_tru_(nm),N_k_p_(nm*n_A))
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,N_k_p_(nm*n_A),N_k_q_(nm*n_A))         
      enddo !do nm=0,n_M-1

c$$$c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$c$$$      Can use to debug
c$$$c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      ns = 0
c$$$      nm = 0
c$$$      write(6,'(A,I0,1X,F8.4)') 'delta_x_tru ', nm ,
c$$$     $     delta_x_tru_(nm)
c$$$      write(6,'(A,I0,1X,F8.4)') 'delta_y_tru ', nm ,
c$$$     $     delta_y_tru_(nm)
c$$$      write(6,'(A,I0,1X,F8.4)') 'gamma_z_tru ', nm ,
c$$$     $     gamma_z_tru_(nm)
c$$$      call test_innerproduct_batch_excerpt_7(n_delta_x,delta_x_
c$$$     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r,grid_k_p_,n_w_,n_A
c$$$     $     ,S_k_p_(ns*n_A),N_k_p_(nm*n_A),Z_p_)
c$$$      call test_innerproduct_batch_excerpt_8(n_delta_x,delta_x_
c$$$     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,Z_p_,ndx_optimal
c$$$     $     ,ndy_optimal ,ngz_optimal ,C_Z_optimal)
c$$$      write(6,'(A)') ' test_innerproduct_fast_dr_2: '
c$$$      write(6,'(A,I0,1X,F8.4)') 'delta_x_ at ', ndx_optimal ,
c$$$     $     delta_x_(ndx_optimal)
c$$$      write(6,'(A,I0,1X,F8.4)') 'delta_y_ at ', ndy_optimal ,
c$$$     $     delta_y_(ndy_optimal)
c$$$      write(6,'(A,I0,1X,F8.4)') 'gamma_z_ at ', ngz_optimal ,
c$$$     $     gamma_z_(ngz_optimal)
c$$$      write(6,'(A,F8.4,1X,F8.4)') 'C_Z_optimal ', C_Z_optimal

      allocate(azimu_b_est_(0:n_M-1))
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we set up the estimated translations'
     $        ,' and rotations for each image.'
      end if !if (verbose.gt.1) then
      
c$$$      The following code sets up a simple test
c$$$      where the estimated alignments correspond to the true ones
      if (flag_tru) then
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'As a first test we choose these estimated'
     $        ,' translations and rotations to be the same as the'
     $        ,' true ones.'
      end if
      do nm=0,n_M-1
         ns = max(0,min(n_S-1,nm/3))
         azimu_b_est_(nm) = azimu_b_tru_(nm) 
         delta_x_est_(nm) = delta_x_tru_(nm) 
         delta_y_est_(nm) = delta_y_tru_(nm) 
         gamma_z_est_(nm) = gamma_z_tru_(nm) 
         l2_norm_est_(nm) = l2_norm_tru_(nm) 
      enddo
      else !if (flag_tru) then
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ''
c$$$         write(6,*) 'For this test we choose these estimated'
c$$$     $        ,' translations and rotations to be zero.'
c$$$      end if
c$$$      do nm=0,n_M-1
c$$$         azimu_b_est_(nm) = 0.0d0
c$$$         delta_x_est_(nm) = 0.0d0
c$$$         delta_y_est_(nm) = 0.0d0
c$$$         gamma_z_est_(nm) = 0.0d0
c$$$         l2_norm_est_(nm) = 1.0d0
c$$$      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'For this test we choose these estimated'
     $        ,' translations and rotations to be discrete shifts'
     $        ,' away from the true translations and rotations.'
      end if
      do nm=0,n_M-1
         azimu_b_est_(nm) = 0.0d0
         l2_norm_est_(nm) = 1.0d0
         ndx = max(0,min(n_delta_x-1,mod(nm+2,3)))
         ndy = max(0,min(n_delta_y-1,mod(nm+3,5)))
         ngz = max(0,min(n_gamma_z-1,mod(nm+5,7)))
         delta_x = delta_x_(ndx)
         delta_y = delta_y_(ndy)
         gamma_z = gamma_z_(ngz)
         gamma_z_est_(nm) = gamma_z_tru_(nm) - gamma_z
         delta_x_pre = delta_x_tru_(nm) - delta_x
         delta_y_pre = delta_y_tru_(nm) - delta_y
         call get_interchange_delta(delta_x_pre,delta_y_pre
     $        ,-gamma_z,delta_x_pos,delta_y_pos)
         delta_x_est_(nm) = delta_x_pos
         delta_y_est_(nm) = delta_y_pos
      enddo
      end if ! if (flag_tru) then

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') ' Now we generate some ctf-functions.'
         write(6,'(A)')
     $        ' These are deliberately chosen to be anisotropic.'
      end if
      n_ctf = 3
      allocate(CTF_k_p_(0:n_ctf*n_A-1))
c$$$      Here we use an unrealistic ctf-function with a high degree
c$$$      of anisotropy to test our code.
      do nctf=0,n_ctf-1
         ctf_p_0 = max_k_c/3.0
         ctf_p_1 = (2*pi*nctf)/n_ctf
         ctf_p_2 = 2
         if (flag_ctf) then
         call get_ctf_star_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_(nctf
     $        *n_A),ctf_p_0 ,ctf_p_1,ctf_p_2)
         else
         call get_ctf_ones_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_(nctf
     $        *n_A))
         end if ! if (flag_ctf) then
      enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         Determine which images correspond to which ctf-function.
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_ctf_(0:n_M-1))
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo

      displacement_max = 3.0d0/max_k_c
      if (verbose.gt.0) then
         write(6,'(A)') ' Setting displacement_max to 1.0d0 '
      end if ! if (verbose.gt.0) then
      displacement_max = 1.0d0
      allocate(alpha_tru_(0:n_alpha*n_M-1))
      allocate(alpha_est_(0:n_alpha*n_M-1))
      do nm=0,n_M-1
         alpha_tru_(nalpha_polar_a + n_alpha*nm) = 0.0d0
         alpha_tru_(nalpha_azimu_b + n_alpha*nm) = azimu_b_tru_(nm)
         alpha_tru_(nalpha_gamma_z + n_alpha*nm) = gamma_z_tru_(nm)
         alpha_tru_(nalpha_delta_x + n_alpha*nm) = delta_x_tru_(nm)
         alpha_tru_(nalpha_delta_y + n_alpha*nm) = delta_y_tru_(nm)
         alpha_tru_(nalpha_l2_norm + n_alpha*nm) = l2_norm_tru_(nm)
         alpha_tru_(nalpha_ctf_ind + n_alpha*nm) = 1.0d0*I_ctf_(nm)
         alpha_tru_(nalpha_S_index + n_alpha*nm) = 0.0d0
         alpha_est_(nalpha_polar_a + n_alpha*nm) = 0.0d0
         alpha_est_(nalpha_azimu_b + n_alpha*nm) = azimu_b_est_(nm)
         alpha_est_(nalpha_gamma_z + n_alpha*nm) = gamma_z_est_(nm)
         alpha_est_(nalpha_delta_x + n_alpha*nm) = delta_x_est_(nm)
         alpha_est_(nalpha_delta_y + n_alpha*nm) = delta_y_est_(nm)
         alpha_est_(nalpha_l2_norm + n_alpha*nm) = l2_norm_est_(nm)
         alpha_est_(nalpha_ctf_ind + n_alpha*nm) = 1.0d0*I_ctf_(nm)
         alpha_est_(nalpha_S_index + n_alpha*nm) = 0.0d0
      enddo ! do nm=0,n_M-1
      allocate(alpha_polar_a_(0:n_S-1))
      allocate(alpha_azimu_b_(0:n_S-1))
      do ns=0,n_S-1
         alpha_polar_a_(ns) = pi/2
         alpha_azimu_b_(ns) = ns*1.0d0
      enddo

      allocate(M_k_p_(0:n_A*n_M-1))
      allocate(C_M_(0:n_M-1))

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         nctf = I_ctf_(nm)
         call xx1_c16(n_A,N_k_p_(nm*n_A),CTF_k_p_(nctf
     $        *n_A),M_k_p_(nm*n_A))
      enddo                     ! do nm=0,n_M-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         l2_norm = alpha_tru_(nalpha_l2_norm+n_alpha*nm)
         call af1_c16(n_A,1.0d0*cmplx(l2_norm,0.0d0),1.0d0
     $        *cmplx(0.0d0,0.0d0),M_k_p_(nm*n_A)
     $        ,M_k_p_(nm*n_A))
      enddo                     ! do nm=0,n_M-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(6,'(A,I0,A)')
     $        ' Trying to generate Figure dir_jpg/Fig_MNC_c_',n_r,'.jpg'
         write(fname,'(A,I0)') './dir_jpg/Fig_MNC_c_',n_r
         call Fig_gen_ctf_ver2(n_r,n_polar_a_,grid_k_p_,n_M,n_A,N_k_p_
     $        ,M_k_p_,n_A,CTF_k_p_,n_alpha,alpha_est_,min(16,n_M),fname)
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $     ,n_M,')-x-(n_S = ',n_S, ') = (', n_M
     $     *n_S,') innerproducts'
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_wrapper_1(verbose,n_omp_sub,n_r
     $     ,n_polar_a_,grid_k_p_,max_x_c,n_S,n_A ,S_k_p_
     $     ,tesselation_distance_req ,alpha_polar_a_,alpha_azimu_b_ ,n_M
     $     ,n_A,M_k_p_ ,n_CTF,n_A ,CTF_k_p_,n_alpha ,alpha_est_
     $     ,alpha_update_f,flag_MS_vs_SM,N_pixels_in ,displacement_max
     $     ,n_delta_x ,n_delta_y,n_gamma_z, svd_calculation_type
     $     ,eps_svd ,fpm_howmany_max,C_M_ ,n_SM_use)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'test_innerproduct_fast_wrapper_1:'
     $        ,' total_time ',timing_toc-timing_tic
      end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(6,'(A,I0,A)')
     $        ' Trying to generate Figure dir_jpg/Fig_MTS_c_',n_r,'.jpg'
         write(fname,'(A,I0)') './dir_jpg/Fig_MTS_c_',n_r
         call Fig_gen_ver4(n_r,n_polar_a_,grid_k_p_,n_S,n_A ,S_k_p_,n_M
     $        ,n_A,M_k_p_ ,n_A,CTF_k_p_ ,n_alpha,alpha_est_ ,min(16
     $        ,n_M),fname)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we update the estimated displacement'
     $        ,' and rotation parameters.'
      end if
      allocate(azimu_b_upd_(0:n_M-1))
      allocate(delta_x_upd_(0:n_M-1))
      allocate(delta_y_upd_(0:n_M-1))
      allocate(gamma_z_upd_(0:n_M-1))
      allocate(l2_norm_upd_(0:n_M-1))
      do nm=0,n_M-1
         azimu_b_upd_(nm) = alpha_est_(nalpha_azimu_b + n_alpha*nm)
         delta_x_upd_(nm) = alpha_est_(nalpha_delta_x + n_alpha*nm)
         delta_y_upd_(nm) = alpha_est_(nalpha_delta_y + n_alpha*nm)
         gamma_z_upd_(nm) = alpha_est_(nalpha_gamma_z + n_alpha*nm)
         l2_norm_upd_(nm) = alpha_est_(nalpha_l2_norm + n_alpha*nm)
      enddo !do nm=0,n_M-1
      if (verbose.gt.0 .and. n_M.le.15) then
         write(6,'(A)') ''
         write(6,*)
     $        'Now we compare updated estimate to true parameters: '
         write(format_string,'(A,A,A)')
     $        '(A,I2,'
     $        ,'A,I0, 1X,F8.3, 1X,F8.3, 1X,F8.3, 1X,F8.3, '
     $        ,'A,I0, 1X,F8.3, 1X,F8.3, 1X,F8.3, 1X,F8.3)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; upd_:' ,
     $           nint(azimu_b_upd_(nm)) , delta_x_upd_(nm) ,
     $           delta_y_upd_(nm), gamma_z_upd_(nm), l2_norm_upd_(nm) ,
     $           '; tru_:' , nint(azimu_b_tru_(nm)) , delta_x_tru_(nm)
     $           ,delta_y_tru_(nm) , gamma_z_tru_(nm) , l2_norm_tru_(nm)
         enddo !do nm=0,n_M-1
      end if !if (verbose.gt.0 .and. n_M.le.15) then
      allocate(azimu_b_err_(0:n_M-1))
      allocate(delta_x_err_(0:n_M-1))
      allocate(delta_y_err_(0:n_M-1))
      allocate(gamma_z_err_(0:n_M-1))
      allocate(l2_norm_rer_(0:n_M-1))
      do nm=0,n_M-1
         azimu_b_err_(nm) = dabs(azimu_b_tru_(nm)-azimu_b_upd_(nm))
         delta_x_err_(nm) = dabs(delta_x_tru_(nm)-delta_x_upd_(nm))
         delta_y_err_(nm) = dabs(delta_y_tru_(nm)-delta_y_upd_(nm))
         gamma_z_err_(nm) = dabs(gamma_z_tru_(nm)-gamma_z_upd_(nm))
         call periodize_r8(gamma_z_err_(nm),-pi,+pi,gamma_z_err_(nm))
         l2_norm_rer_(nm) = dabs(l2_norm_tru_(nm)-l2_norm_upd_(nm))
     $        /l2_norm_tru_(nm)
      enddo !do nm=0,n_M-1
      if (verbose.gt.0 .and. n_M.le.15) then
         write(6,*) 'Alignment Errors:'
         write(format_string,'(A)')
     $        '(A,I2,A,I0,A,F8.6,A,F8.6,A,F9.4,A,F9.4)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; azimu_b_: ' ,
     $           nint(azimu_b_err_(nm)) , '; delta_x_:'
     $           ,delta_x_err_(nm) , '; delta_y_:' , delta_y_err_(nm)
     $           ,'; gamma_z_:' , gamma_z_err_(nm) , '; l2_norm_: '
     $           ,l2_norm_rer_(nm)
         enddo !do nm=0,n_M-1
      end if !if (verbose.gt.0 .and. n_M.le.15) then
      if (verbose.gt.0) then
         write(6,*) 'Average Alignment Errors:'
         write(format_string,'(A)')
     $        '(A,I0,1X,A,F8.6,1X,A,F8.6,1X,A,F9.4,1X,A,F9.4)'
         write(6,format_string) 
     $    'azimu_b_err_:', nint(avg_r8_f(n_M,azimu_b_err_)), 
     $    'delta_x_err_:', avg_r8_f(n_M,delta_x_err_), 
     $    'delta_y_err_:', avg_r8_f(n_M,delta_y_err_), 
     $    'gamma_z_err_:', avg_r8_f(n_M,gamma_z_err_),
     $    'l2_norm_rer_:', avg_r8_f(n_M,l2_norm_rer_)
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,*) '|Alignment Errors|_{2}:'
         write(format_string,'(A)')
     $        '(A,I0,1X,A,F8.6,1X,A,F8.6,1X,A,F9.4,1X,A,F9.4)'
         write(6,format_string) 
     $    'azimu_b_err_:', nint(al2_r8_f(n_M,azimu_b_err_)), 
     $    'delta_x_err_:', al2_r8_f(n_M,delta_x_err_), 
     $    'delta_y_err_:', al2_r8_f(n_M,delta_y_err_), 
     $    'gamma_z_err_:', al2_r8_f(n_M,gamma_z_err_),
     $    'l2_norm_rer_:', al2_r8_f(n_M,l2_norm_rer_)
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,*) 'Largest Alignment Errors:'
         write(format_string,'(A)')
     $        '(A,I0,1X,A,F8.6,1X,A,F8.6,1X,A,F9.4,1X,A,F9.4)'
         write(6,format_string) 
     $    'azimu_b_err_:', nint(max_r8_f(n_M,azimu_b_err_)), 
     $    'delta_x_err_:', max_r8_f(n_M,delta_x_err_), 
     $    'delta_y_err_:', max_r8_f(n_M,delta_y_err_), 
     $    'gamma_z_err_:', max_r8_f(n_M,gamma_z_err_),
     $    'l2_norm_rer_:', max_r8_f(n_M,l2_norm_rer_)
      end if !if (verbose.gt.0) then

 10    continue

      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_innerproduct_fast_dr_2]'
      end if

      stop
      end
