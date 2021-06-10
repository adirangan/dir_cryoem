c     Compile with: gfortran -fcray-pointer -w -O2 -fPIC --openmp -o test_innerproduct_batch_a_dr.out test_innerproduct_batch_a_dr.f -lfftw3 ; 
c     Run with: ./test_innerproduct_batch_a_dr.out 1 24 0.1 3.0 25 25 48 -1 ; 
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program test_innerproduct_batch_a_dr
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This program is a driver which tests the subroutines
c$$$      test_innerproduct_batch_wrapper_2a.f ,
c$$$      test_innerproduct_batch_wrapper_3a.f ,
c$$$      test_innerproduct_batch_wrapper_4a.f ,
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
c$$$      D. Then we call test_innerproduct_batch_wrapper_2a to calculate
c$$$         the innerproducts between each pair of templates and
c$$$         images. These innerproducts are calculated across a range
c$$$         of displacements delta_x_,delta_y_ as well as a range of 
c$$$         rotations gamma_z_ (all passed in as input).
c$$$      D. The output of test_innerproduct_batch_wrapper_2a is used
c$$$         to determine (among other things) which templates correlate
c$$$         most strongly with each image, and what the best
c$$$         alignment (i.e., displacement and rotation) is for each pair.
c$$$      E. The true alignment of each image is then compared with
c$$$         the best empirically observed alignment, and the alignment-
c$$$         errors are printed to stdout.
c$$$      F. Finally, if (verbose.gt.1), we generate two figure files:
c$$$         ./dir_jpg/Fig_NMC_[ncur].jpg:
c$$$              R(M_k_p_(nm)) and I(M_k_p_(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid.
c$$$              R(C_k_p_(nm)) and I(C_k_p_(nm)) show real and imaginary
c$$$              parts of ctf-function I_ctf_(nm) on k-space quasi-uniform 
c$$$              polar-grid.
c$$$              R(N_k_p_(nm)) and I(N_k_p_(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid
c$$$              after pointwise multiplication by ctf-function and
c$$$              multiplication by image_specific scaling factor l2_norm.
c$$$         ./dir_jpg/Fig_MTS_[ncur].jpg:
c$$$              Each pair of rows shows the real and imaginary parts of
c$$$              some function on k-space quasi-uniform polar-grid.
c$$$              M_k_p_ is the primitive image
c$$$              N_k_p_ is the image multiplied by ctf and scaled
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
c$$$      M_x_c_ through M_k_q_ are defined similarly for the primitive 
c$$$             images (i.e., before applying ctf or scaling).
c$$$      N_x_c_ through N_k_q_ are defined similarly for the images
c$$$             after multipliction by ctf-function.
c$$$      grid_x_c_: grid-values of real-space cartesian-grid
c$$$      grid_k_c_: grid-values of fourier-space cartesian-grid
c$$$      grid_k_p_: radial-values of rings of fourier-space polar-grid
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$
c$$$      INPUTS: This driver requires 8 parameters as command-line inputs.
c$$$              These inputs should be listed in the order below:
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
c$$$      real *8 eps_target    The epsilon used as a cutoff when 
c$$$                            determining which terms to retain within the
c$$$                            svd-expansion of the various bessel-
c$$$                            functions used to approximate the term
c$$$                            exp(2*pi*K*delta*cos(psi-omega)) within
c$$$                            the fourier-space translation operator.
c$$$                          - A value of eps_target.eq.10.0d0 corresponds
c$$$                            to low accuracy, and will probably not give
c$$$                            exactly the same result as an exact brute-
c$$$                            force calculation 
c$$$                            (see test_innerproduct_svdread_dr.f for a 
c$$$                            direct comparison).
c$$$                          - A value of eps_target.eq.1.0d0 corresponds
c$$$                            to moderate accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to a few digits, which is 
c$$$                            usually enough for the best alignment to 
c$$$                            match exactly.
c$$$                          - A value of eps_target.eq.0.1d0 corresponds
c$$$                            to high accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to several digits, which is 
c$$$                            more than enough to match the best alignment
c$$$                            exactly.
c$$$                          - Values of eps_target.lt.0.1d0 are not
c$$$                            usually necessary, and we have not generated
c$$$                            the appropriate libraries for these 
c$$$                            eps_target yet.
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
c$$$                          - A value of svd_calculation_type.eq.0
c$$$                            forces the displacement operation to take
c$$$                            place first, and only later uses the fft
c$$$                            to account for the various rotations.
c$$$                            This calculation is useful when the number
c$$$                            of distinct displacements is smaller than
c$$$                            the number of terms in the svd-expansion.
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
c$$$                            to automatically choose one of the three
c$$$                            calculation_types based on the other
c$$$                            input parameters (e.g., by comparing
c$$$                            n_svd_l with n_delta_x*n_delta_y).
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 2 /
      logical plot_flag
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na
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
      real *8 max_r8_f
      real *8 max_x_c,max_k_c
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
c$$$      Templates S_
      integer n_S,ld_S,ns
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: T_x_c_(:)
      complex *16, allocatable :: S_x_p_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: S_k_q_(:)
c$$$      Images M_ and N_
      integer n_M,ld_M,nm
      complex *16, allocatable :: M_x_c_(:)
      complex *16, allocatable :: M_x_p_(:)
      complex *16, allocatable :: M_k_c_(:)
      complex *16, allocatable :: M_k_p_(:)
      complex *16, allocatable :: M_k_q_(:)
      complex *16, allocatable :: N_k_p_(:)
c$$$      CTF-functions
      integer *4 n_ctf,ld_CTF,nctf
      integer *4, allocatable :: I_ctf_(:)
      complex *16, allocatable :: CTF_k_p_(:)
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4
c$$$      image parameters 
      integer *4 n_alpha
      parameter (n_alpha=7)
      include 'nalpha_define.f'
      real *8 l2_norm
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: alpha_polar_a_(:)
      real *8, allocatable :: alpha_azimu_b_(:)
      real *8, allocatable :: delta_x_max_(:)
      real *8, allocatable :: delta_y_max_(:)
      real *8, allocatable :: gamma_z_max_(:)
      complex *16, allocatable :: C_S_max_(:)
      complex *16, allocatable :: C_Z_max_(:)
      real *8, allocatable :: delta_x_sort_SM_(:)
      real *8, allocatable :: delta_y_sort_SM_(:)
      real *8, allocatable :: gamma_z_sort_SM_(:)
      complex *16, allocatable :: C_S_sort_SM_(:)
      complex *16, allocatable :: C_Z_sort_SM_(:)
      integer *4, allocatable :: I_permute_SM_(:)
      integer *4, allocatable :: I_inverse_SM_(:)
      real *8, allocatable :: delta_x_sort_MS_(:)
      real *8, allocatable :: delta_y_sort_MS_(:)
      real *8, allocatable :: gamma_z_sort_MS_(:)
      complex *16, allocatable :: C_S_sort_MS_(:)
      complex *16, allocatable :: C_Z_sort_MS_(:)
      integer *4, allocatable :: I_permute_MS_(:)
      integer *4, allocatable :: I_inverse_MS_(:)
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
      real *8, allocatable :: l2_norm_err_(:)
c$$$      range of displacement and rotation parameters to measure
      real *8 displacement_max,angle_tmp,delta_tmp,dtmp
      integer n_delta_x,n_delta_y,ndx,ndy
      real *8 eps_target,N_pixels_in,delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      integer n_gamma_z,ngz
      real *8, allocatable :: gamma_z_(:)      
      integer svd_calculation_type
c$$$      parameters for figure generation
      character(len=1024) fname,format_string
c$$$      function for template generation
      external get_F2_x_c
      real *8 param_1,param_2
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      parameters for command line input
c$$$      list: verbose ncur eps_target N_pixels_in n_delta_x n_delta_y n_gamma_z svd_calculation_type
      character(len=8) :: cmd_argstring
      integer cmd_i0,cmd_i1,cmd_i4,cmd_i5,cmd_i6,cmd_i7
      real *8 cmd_d2,cmd_d3
      call get_command_argument(1+0,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i0
      call get_command_argument(1+1,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i1
      call get_command_argument(1+2,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d2
      call get_command_argument(1+3,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d3
      call get_command_argument(1+4,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i4
      call get_command_argument(1+5,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i5
      call get_command_argument(1+6,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i6
      call get_command_argument(1+7,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i7

      verbose = cmd_i0
      ncur = cmd_i1
      eps_target = cmd_d2
      N_pixels_in = cmd_d3
      n_delta_x = cmd_i4
      n_delta_y = cmd_i5
      n_gamma_z = cmd_i6
      svd_calculation_type = cmd_i7
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_batch_a_dr]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,F6.3)') 'N_pixels_in: ',N_pixels_in
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
      end if
      
      if (verbose.ge.2) then
         plot_flag = .true.
      else
         plot_flag = .false.
      end if

      pi = 4*atan(1.0)

c$$$      Calculating template size using 'get_template_size'
      allocate(nlats_(0:ncur-1));
      do nr=0,ncur-1
c$$$         nlats calculation similar to testrebuild.f:
         nlats_(nr) = nint(pi*(nr+1))
      enddo
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'ncur = ',ncur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
      end if
      
c$$$      indices
      n_r = ncur
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
      end if

c$$$      fftw
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_ESTIMATE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_ESTIMATE) 
         na = na + n_w_(nr)
      enddo
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

c$$$      generate templates
      n_S = 5
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Here we construct the templates "S_".'
         write(6,'(A,I0,A)') 'We generate n_S = ',n_S
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
         write(6,*) 'on a quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
      end if
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
      enddo

      n_M = 3*n_S
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A,I0,A)') 'Now we do the same thing for the n_M = '
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
      end if
      allocate(M_x_c_(0:n_r*n_r*n_M-1))
      allocate(M_k_c_(0:n_r*n_r*n_M-1))
      allocate(M_k_p_(0:n_A*n_M-1))
      allocate(M_k_q_(0:n_A*n_M-1))
c$$$      Calculating x-space image on regular cartesian-grid
      allocate(azimu_b_tru_(0:n_M-1))
      allocate(delta_x_tru_(0:n_M-1))
      allocate(delta_y_tru_(0:n_M-1))
      allocate(gamma_z_tru_(0:n_M-1))
      allocate(l2_norm_tru_(0:n_M-1))
      do nm=0,n_M-1
         ns = nm/3
         azimu_b_tru_(nm) = 1.0d0*ns
         delta_x_tru_(nm) = 1.0*(mod(nm,3)-1)/n_r*max_x_c
         delta_y_tru_(nm) = 2.0*(mod(nm,3)-1)/n_r*max_x_c
         gamma_z_tru_(nm) = mod(nm,3)*2*pi/3.0
         l2_norm_tru_(nm) = 1.0d0 + nm
         if (verbose.gt.1) then
            write(6,'(I0,A,4F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_tru_(nm),delta_y_tru_(nm),gamma_z_tru_(nm)
     $           ,l2_norm_tru_(nm)
         end if
         call cp1_c16(n_r*n_r,S_x_c_(ns*n_r*n_r),M_x_c_(nm*n_r*n_r))
         call transl_c_to_c(n_r,max_x_c,n_r,max_x_c,M_x_c_(nm*n_r*n_r),
     $        +delta_x_tru_(nm),+delta_y_tru_(nm),M_x_c_(nm*n_r*n_r))
         call rotate_c_to_c(n_r,max_x_c,n_r,max_x_c,M_x_c_(nm*n_r*n_r),
     $        +gamma_z_tru_(nm),M_x_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on regular cartesian-grid      
         call adi_fft2(-1,n_r,n_r,M_x_c_(nm*n_r*n_r),M_k_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_k_c,n_r,max_k_c,M_k_c_(nm*n_r*n_r)
     $        ,n_r,grid_k_p_,n_w_,n_A,M_k_p_(nm*n_A))
c$$$      Calculating J-space image on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,M_k_p_(nm*n_A),M_k_q_(nm*n_A))
      enddo

c$$$c$$$      The following (commented-out) code sets up a simpler test
c$$$c$$$      where the estimated alignments correspond to the true ones
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ''
c$$$         write(6,*) 'Now we set up the estimated translations'
c$$$     $        ,' and rotations for each image.'
c$$$         write(6,*) 'As a first test we choose these estimated'
c$$$     $        ,' translations and rotations to be the same as the'
c$$$     $        ,' true ones.'
c$$$      end if
c$$$      allocate(delta_x_est_(0:n_M-1))
c$$$      allocate(delta_y_est_(0:n_M-1))
c$$$      allocate(gamma_z_est_(0:n_M-1))
c$$$      do nm=0,n_M-1
c$$$         delta_x_est_(nm) = 1.0*(mod(nm,3)-1)/n_r*max_x_c
c$$$         delta_y_est_(nm) = 2.0*(mod(nm,3)-1)/n_r*max_x_c
c$$$         gamma_z_est_(nm) = mod(nm,3)*2*pi/3.0
c$$$      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we set up the estimated translations'
     $        ,' and rotations for each image.'
         write(6,*) 'For this test we choose these estimated'
     $        ,' translations and rotations to be close to the'
     $        ,' true ones (but not exactly the same).'
      end if
      allocate(azimu_b_est_(0:n_M-1))
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      do nm=0,n_M-1
         azimu_b_est_(nm) = 0.0d0
         delta_x_est_(nm) = delta_x_tru_(nm)*(1 + 0.125*(mod(nm+0,3)-1))
         delta_y_est_(nm) = delta_y_tru_(nm)*(1 + 0.125*(mod(nm+1,3)-1))
         gamma_z_est_(nm) = gamma_z_tru_(nm)*(1 + 0.125*(mod(nm+2,3)-1))
         l2_norm_est_(nm) = 1.0d0
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Now we generate some ctf-functions.'
         write(6,'(A)')
     $        'These are deliberately chosen to be anisotropic.'
      end if
      n_ctf = 3
      allocate(CTF_k_p_(0:n_ctf*n_A-1))
c$$$      Here we use an unrealistic ctf-function with a high degree
c$$$      of anisotropy to test our code.
      do nctf=0,n_ctf-1
         ctf_p_0 = max_k_c/3.0
         ctf_p_1 = (2*pi*nctf)/n_ctf
         ctf_p_2 = 2
         call get_ctf_star_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_(nctf
     $        *n_A),ctf_p_0 ,ctf_p_1,ctf_p_2)
c$$$         call get_ctf_ones_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p_(nctf
c$$$     $        *n_A))
      enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         Determine which images correspond to which ctf-function.
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_ctf_(0:n_M-1))
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo

      displacement_max = 3.0d0/max_k_c
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
         alpha_est_(nalpha_polar_a + n_alpha*nm) = 0.0d0
         alpha_est_(nalpha_azimu_b + n_alpha*nm) = azimu_b_est_(nm)
         alpha_est_(nalpha_gamma_z + n_alpha*nm) = gamma_z_est_(nm)
         alpha_est_(nalpha_delta_x + n_alpha*nm) = delta_x_est_(nm)
         alpha_est_(nalpha_delta_y + n_alpha*nm) = delta_y_est_(nm)
         alpha_est_(nalpha_l2_norm + n_alpha*nm) = l2_norm_est_(nm)
         alpha_est_(nalpha_ctf_ind + n_alpha*nm) = 1.0d0*I_ctf_(nm)
      enddo ! do nm=0,n_M-1
      allocate(alpha_polar_a_(0:n_S-1))
      allocate(alpha_azimu_b_(0:n_S-1))
      do ns=0,n_S-1
         alpha_polar_a_(ns) = 0.0d0
         alpha_azimu_b_(ns) = ns*1.0d0
      enddo

      allocate(N_k_p_(0:n_A*n_M-1))
      allocate(C_M_(0:n_M-1))
      allocate(delta_x_max_(0:n_S*n_M-1))
      allocate(delta_y_max_(0:n_S*n_M-1))
      allocate(gamma_z_max_(0:n_S*n_M-1))
      allocate(C_S_max_(0:n_S*n_M-1))
      allocate(C_Z_max_(0:n_S*n_M-1))
      allocate(delta_x_sort_SM_(0:n_S*n_M-1))
      allocate(delta_y_sort_SM_(0:n_S*n_M-1))
      allocate(gamma_z_sort_SM_(0:n_S*n_M-1))
      allocate(C_S_sort_SM_(0:n_S*n_M-1))
      allocate(C_Z_sort_SM_(0:n_S*n_M-1))
      allocate(I_permute_SM_(0:n_S*n_M-1))
      allocate(I_inverse_SM_(0:n_S*n_M-1))
      allocate(delta_x_sort_MS_(0:n_S*n_M-1))
      allocate(delta_y_sort_MS_(0:n_S*n_M-1))
      allocate(gamma_z_sort_MS_(0:n_S*n_M-1))
      allocate(C_S_sort_MS_(0:n_S*n_M-1))
      allocate(C_Z_sort_MS_(0:n_S*n_M-1))
      allocate(I_permute_MS_(0:n_S*n_M-1))
      allocate(I_inverse_MS_(0:n_S*n_M-1))

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         nctf = I_ctf_(nm)
         call xx1_c16(n_A,M_k_p_(nm*n_A),CTF_k_p_(nctf
     $        *n_A),N_k_p_(nm*n_A))
      enddo                     ! do nm=0,n_M-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         l2_norm = alpha_tru_(nalpha_l2_norm+n_alpha*nm)
         call af1_c16(n_A,1.0d0*cmplx(l2_norm,0.0d0),1.0d0
     $        *cmplx(0.0d0,0.0d0),N_k_p_(nm*n_A)
     $        ,N_k_p_(nm*n_A))
      enddo                     ! do nm=0,n_M-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (plot_flag) then
         write(fname,'(A,I0)') './dir_jpg/Fig_MNC_',n_r
         call Fig_gen_ctf_ver1(n_r,nlats_,grid_k_p_,n_M,n_A,M_k_p_
     $        ,N_k_p_,n_A,CTF_k_p_,n_alpha,alpha_est_,min(16,n_M),fname)
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $     ,n_M,')-x-(n_S = ',n_S, ') = (', n_M
     $     *n_S,') innerproducts'
      timing_tic = omp_get_wtime()
      call test_innerproduct_batch_wrapper_2a(verbose,n_r,nlats_
     $     ,grid_k_p_,max_x_c,n_S,n_A,S_k_p_,n_M ,n_A ,N_k_p_,n_CTF,n_A
     $     ,CTF_k_p_,n_alpha ,alpha_est_,0.1d0 ,N_pixels_in
     $     ,displacement_max,n_delta_x ,n_delta_y,max(32,2 *n_r), -1
     $     ,delta_x_max_,delta_y_max_ ,gamma_z_max_,C_M_ ,C_S_max_
     $     ,C_Z_max_)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'test_innerproduct_batch_wrapper_2a:'
     $        ,' total_time ',timing_toc-timing_tic
      end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  sort results
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      call test_innerproduct_batch_wrapper_3a(verbose,n_S,n_M
     $     ,delta_x_max_,delta_y_max_,gamma_z_max_,C_S_max_,C_Z_max_
     $     ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $     ,C_S_sort_SM_ ,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (plot_flag) then
         write(fname,'(A,I0)') './dir_jpg/Fig_MTS_',n_r
         call Fig_gen_ver1a(n_r,nlats_,grid_k_p_,n_S,n_A ,S_k_p_,n_M,n_A
     $        ,M_k_p_,N_k_p_ ,n_A,CTF_k_p_ ,n_alpha,alpha_est_
     $        ,delta_x_sort_SM_ ,delta_y_sort_SM_ ,gamma_z_sort_SM_
     $        ,I_permute_SM_,min(16 ,n_M),fname)
      end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  update alphas
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      call test_innerproduct_batch_wrapper_4a(verbose,n_S,n_M ,C_M_
     $     ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $     ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $     ,delta_x_sort_MS_,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_,I_inverse_MS_
     $     ,alpha_polar_a_ ,alpha_azimu_b_,n_alpha ,alpha_est_) 

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
      enddo
      if (verbose.gt.-1) then
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
         enddo
      end if
      allocate(azimu_b_err_(0:n_M-1))
      allocate(delta_x_err_(0:n_M-1))
      allocate(delta_y_err_(0:n_M-1))
      allocate(gamma_z_err_(0:n_M-1))
      allocate(l2_norm_err_(0:n_M-1))
      do nm=0,n_M-1
         azimu_b_err_(nm) = dabs(azimu_b_tru_(nm)-azimu_b_upd_(nm))
         delta_x_err_(nm) = dabs(delta_x_tru_(nm)-delta_x_upd_(nm))
         delta_y_err_(nm) = dabs(delta_y_tru_(nm)-delta_y_upd_(nm))
         gamma_z_err_(nm) = dabs(gamma_z_tru_(nm)-gamma_z_upd_(nm))
         call periodize_r8(gamma_z_err_(nm),-pi,+pi,gamma_z_err_(nm))
         l2_norm_err_(nm) = dabs(l2_norm_tru_(nm)-l2_norm_upd_(nm))
      enddo
      if (verbose.gt.-1) then
         write(6,*) 'Alignment Errors:'
         write(format_string,'(A)')
     $        '(A,I2,A,I0,A,F8.6,A,F8.6,A,F8.6,A,F8.6)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; azimu_b_: ' ,
     $           nint(azimu_b_err_(nm)) , '; delta_x_:'
     $           ,delta_x_err_(nm) , '; delta_y_:' , delta_y_err_(nm)
     $           ,'; gamma_z_:' , gamma_z_err_(nm) , '; l2_norm_: '
     $           ,l2_norm_err_(nm)
         enddo
      end if
      if (verbose.gt.0) then
         write(6,*) 'Largest Alignment Errors:'
         write(format_string,'(A)')
     $        '(A,I0,1X,A,F8.6,1X,A,F8.6,1X,A,F8.6,1X,A,F8.6)'
         write(6,format_string) 
     $    'azimu_b_err_:', nint(max_r8_f(n_M,azimu_b_err_)), 
     $    'delta_x_err_:', max_r8_f(n_M,delta_x_err_), 
     $    'delta_y_err_:', max_r8_f(n_M,delta_y_err_), 
     $    'gamma_z_err_:', max_r8_f(n_M,gamma_z_err_),
     $    'l2_norm_err_:', max_r8_f(n_M,l2_norm_err_)
      end if
      
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver5]'
      end if

      stop
      end

      real *8 function max_r8_f(n_x,x_)
      integer n_x,nx
      real *8 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0.0d0
      else
         m=x_(0)
         do nx=0,n_x-1
            if (x_(nx).gt.m) then
               m = dabs(x_(nx))
            end if
         enddo
      end if
      max_r8_f = m
      return
      end

      include 'get_template_size.f'
      include 'get_F2_x_c.f'
      include 'get_F2_x_c_.f'
      include 'get_ctf_star_k_p_.f'
      include 'get_ctf_ones_k_p_.f'
      include 'stdlim_c16.f'
      include 'pearson_c16.f'
      include 'cp1_r8.f'
      include 'cp1_c16.f'
      include 'cps_c16.f'
      include 'af1_c16.f'
      include 'afs_c16.f'
      include 'xx1_c16.f'
      include 'xc1_c16.f'
      include 'trn0_r8.f'
      include 'trn0_c16.f'
      include 'linspace.f'
      include 'interp1_c16.f'
      include 'interp2_c16.f'
      include 'interp_c_to_p.f'
      include 'interp_p_to_q.f'
      include 'interp_p_to_q_fftw.f'
      include 'periodize_r8.f'
      include 'periodize_i.f'
      include 'transl_c_to_c.f'
      include 'transf_c_to_c.f'
      include 'transf_p_to_p.f'
      include 'transf_svd_q_to_q.f'
      include 'rotate_c_to_c.f'
      include 'rotate_p_to_p.f'
      include 'rotate_q_to_q.f'
      include 'innerproduct_c.f'
      include 'innerproduct_p.f'
      include 'innerproduct_q__k_only.f'
      include 'innerproduct_q__k_svds.f'
      include 'innerproduct_q__k_svdd.f'
      include 'test_innerproduct_batch_wrapper_2a.f'
      include 'test_innerproduct_batch_wrapper_3a.f'
      include 'test_innerproduct_batch_wrapper_4a.f'
      include 'test_innerproduct_batch_excerpt_00a.f'
      include 'test_innerproduct_batch_excerpt_0a.f'
      include 'test_innerproduct_batch_excerpt_1a.f'
      include 'test_innerproduct_batch_excerpt_2a.f'
      include 'test_innerproduct_batch_excerpt_4a.f'
      include 'test_innerproduct_batch_stage_0a.f'
      include 'test_innerproduct_batch_stage_1a.f'
      include 'test_innerproduct_batch_sort_a.f'
      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Fig_header.f'
      include 'Fig_text.f'
      include 'recenter_c16.f'
      include 'Fig_c16_carte.f'
      include 'Fig_c16_cart2.f'
      include 'Fig_c16_polar.f'
      include 'Fig_c16_bessl.f'
      include 'Fig_c16_bess2.f'
      include 'Fig_gen_ver1a.f'
      include 'Fig_gen_ctf_ver1.f'
      include 'adi_fft1.f'
      include 'adi_fft2.f'
      include 'quicksort_c16.f'
      include 'quicksort_r8.f'
      include 'quicksort_i4.f'
      include 'randperm.f'
