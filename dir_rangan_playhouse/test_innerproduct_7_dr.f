c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program test_innerproduct_dr_7
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This program is a driver which tests the subroutine
c$$$      test_innerproduct_7.f ,
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
c$$$      D. Then we call test_innerproduct_7 to calculate
c$$$         the innerproducts between each pair of templates and
c$$$         images. These innerproducts are calculated across a range
c$$$         of displacements delta_x_,delta_y_ as well as a range of 
c$$$         rotations gamma_z_ (all passed in as input).
c$$$      D. The routine test_innerproduct_7 is used
c$$$         to determine (among other things) which template correlates
c$$$         most strongly with each image, and what the best
c$$$         alignment (i.e., displacement and rotation) is for each image.
c$$$      E. The true alignment of each image is then compared with
c$$$         the best empirically observed alignment, and the alignment-
c$$$         errors are printed to stdout.
c$$$      F. Finally, if (verbose.gt.1), we generate two figure files:
c$$$         ./dir_jpg/Fig_NMC_[ncur].jpg:
c$$$              R(N_k_p__(nm)) and I(N_k_p__(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid.
c$$$              R(C_k_p_(nm)) and I(C_k_p_(nm)) show real and imaginary
c$$$              parts of ctf-function I_ctf_(nm) on k-space quasi-uniform 
c$$$              polar-grid.
c$$$              R(M_k_p__(nm)) and I(M_k_p__(nm)) show real and imaginary
c$$$              parts of image nm on k-space quasi-uniform polar-grid
c$$$              after pointwise multiplication by ctf-function and
c$$$              multiplication by image_specific scaling factor l2_norm.
c$$$         ./dir_jpg/Fig_MTS_[ncur].jpg:
c$$$              Each pair of rows shows the real and imaginary parts of
c$$$              some function on k-space quasi-uniform polar-grid.
c$$$              N_k_p__ is the primitive image
c$$$              M_k_p__ is the image multiplied by ctf and scaled
c$$$              S_k_p__ is the associated best-fit template
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
c$$$      S_x_c__: template in real-space on a cartesian-grid.
c$$$      S_x_p__: template in real-space on a quasi-uniform polar-grid.
c$$$      S_k_c__: template in fourier-space on a cartesian-grid.
c$$$      S_k_p__: template in fourier-space on a quasi-uniform polar-grid.
c$$$      S_k_q__: template in fourier-space, represented using "hat-hat" 
c$$$              coordinates (i.e., bessel-expansion).
c$$$      N_x_c__ through N_k_q__ are defined similarly for the primitive 
c$$$             images (i.e., before applying ctf or scaling).
c$$$      M_x_c__ through M_k_q__ are defined similarly for the images
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
c$$$                            associated with the largest square box which
c$$$                            can be placed within the array of displacements
c$$$                            used when calculating innerproducts for
c$$$                            each pair of templates and images.
c$$$                          - Each pixel of displacement is measured in
c$$$                            x-space. Thus, the maximum displacement is:
c$$$                            delta_x = N_pixels_in/n_r * max_x_c
c$$$                            delta_y = N_pixels_in/n_r * max_x_c.
c$$$                          - For example, N_pixels_in.eq.1.0d0 will limit
c$$$                            the square of displacements to a single pixel
c$$$                            in each direction (i.e., a 3x3 grid of 
c$$$                            pixels centered at the origin).
c$$$                          - A value of N_pixels_in.eq.3.0d0 will 
c$$$                            correspond to a square of displacements
c$$$                            spanning up to 3 pixels in each direction
c$$$                            (i.e., a 7x7 grid of pixels centered at the
c$$$                            origin).
c$$$                            Note: the actual array of displacements
c$$$                            will include the square mentioned above, as
c$$$                            well as a few other displacements which are
c$$$                            within the same maximum distance.
c$$$                          - Note: Values of N_pixels_in.gt.3.0 are not usually
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
c$$$                            vectors is n_delta_v (below).
c$$$
c$$$      integer n_delta_y     Number of displacements to search for in the 
c$$$                            y-direction. These displacements will be 
c$$$                            linearly spaced from -N_pixels_in to 
c$$$                            +N_pixels_in.
c$$$                          - Note that the total number of displacement-
c$$$                            vectors is n_delta_v (below).
c$$$
c$$$      integer n_delta_v     Total number of displacements in the
c$$$                            displacement array. This is not actually an
c$$$                            input, and will be based on the
c$$$                            n_delta_x and n_delta_y above.
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
c$$$                            applies the displacement operation via
c$$$                            svd expansion. Note that the actual
c$$$                            displacement-operator is only applied
c$$$                            after performing the fft.
c$$$                            This particular ordering of operators is
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
c$$$                            n_svd_l with n_delta_v).
c$$$                            (Note: 0-value not implemented here).
c$$$
c$$$      logical flag_RTRT_vs_RTTR flag indicating whether transformation 
c$$$                            operations are ordered as 
c$$$                            either:
c$$$                            < M.*CTF , R_{+est}(T_{+est}(R_{+upd}(T_{+upd}(S)))) >
c$$$                            =
c$$$                            < T_{-est}(R_{-est}(M.*CTF)) , R_{+upd}(T_{+upd}(S)) >
c$$$                            [with flag_RTRT_vs_RTTR.eqv..true.] 
c$$$                            or:
c$$$                            < M.*CTF , R_{+est}(T_{+est}(T_{+upd}(R_{+upd}(S)))) >
c$$$                            =
c$$$                            < T_{-est}(R_{-est}(M.*CTF)) , T_{+upd}(R_{+upd}(S)) >
c$$$                            =
c$$$                            < R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(M.*CTF)))) , S >
c$$$                            [with flag_RTRT_vs_RTTR.eqv..false.].
c$$$                            Note that the former applies translations T_{+upd}
c$$$                            to the templates S, whereas
c$$$                            the latter applies translations T_{-upd} to the
c$$$                            images M.
c$$$
c$$$      integer n_omp_sub     Number of sub-blocks for omp parallelization.
c$$$                            e.g., the number of processors/threads available (e.g., 8 or 16);     
c$$$                            
c$$$      integer n_S_0_sub     Number of sub-blocks for S at level 0. 
c$$$                            This can often be just 1 if you don't have too many (e.g., 1e7) templates.
c$$$                            Often I would guess 1 if possible, or maybe 2-4 if needed?
c$$$                            I exped that if K<200 and you have lots of memory, n_S_0_sub=1 should be fine.
c$$$                            Used to dived templates for precomputation.
c$$$                            Note that precomputatin is slow, so you don't want too many of these.
c$$$                            (Since each n_S_0_sub will require redoing the precomputation for the n_M_0_sub).
c$$$                            
c$$$      integer n_S_1_sub     Number of sub-blocks for S at level 1.
c$$$                            Within each n_S_0_sub block we break the templates into a further n_S_1_sub smaller blocks.
c$$$                            These smaller blocks will be passed into zgemm, and should be just barely large enough for zgemm to operate efficiently.
c$$$                            Typically zgemm operates efficiently if the 'exterior' dimensions are >= 30.
c$$$                            So that means we should choose n_S_1_sub to be large enough that the number of templates processed per zgemm is ~1e1.5 to 1e2 or so.
c$$$                            For example, with n_S = 10000 and n_S_0_sub = 1, then n_S_1_sub should probably be 100-200 ; 
c$$$                            (although 10-20 is fine too if you have enough memory).
c$$$                            
c$$$      integer n_M_0_sub     Number of sub-blocks for M at level 0.
c$$$                            This is often required to be large, especially when we have many images.
c$$$                            It is usually best to make this as small as you can, often 1 or maybe 4 or 10.
c$$$                            
c$$$      integer n_M_1_sub     Number of sub-blocks for M at level 1.
c$$$                            The same reasoning applies here as to n_S_1_sub. 
c$$$                            E.g., with n_M = 100000, and n_M_0_sub = 50, then n_M/n_M_0_sub = 2000 images will be processed at each time.
c$$$                            If you have many translations, this will be much too large for zgemm to handle, so you will need to break up these image-blocks further.
c$$$                            E.g., n_M_1_sub = 100, implying 20 images will be processed at a time (using zgemm).
c$$$                            
c$$$                            Note that, if we are using flag_RTRT_RTTR=.false., the local translation operator will be applied to the images, 
c$$$                            so the effective exterior-dimension associated with the images (for zgemm) will be either n_svd_l*20 or n_delta_v*20. ;
c$$$                            On the other hand, if we are using flag_RTRT_RTTR=.true., the local translation operator will be applied to the templates,
c$$$                            so the effective exterior-dimension associated with the images (for zgemm) will be only 20, 
c$$$                            but the effective exterior-dimension associated with the templates (for zgemm) will be n_svd_l*(n_S/n_S_0_sub/n_S_1_sub). ;
c$$$
c$$$      logical flag_ctf      Whether or not to use anisotropic ctf.
c$$$
c$$$      logical flag_tru      Whether or not to set estimated image parameters
c$$$                            to be equal to the true image parameters.
c$$$
c$$$      integer n_S_max       Maximum number of templates. only n_S will be used.
c$$$
c$$$      integer n_M_max       Maximum number of images. only n_M will be used.
c$$$
c$$$      integer fpm_howmany_max Maximum number of fftws to call simultaneously
c$$$                            within the fftw_plan_many. 
c$$$
c$$$      logical flag_MS_vs_SM determines whether to assign images to templates (.true.) 
c$$$                            or templates to images (.false.). ;
c$$$                            Note: flag_MS_vs_SM=.false. stores the best templates for each image (this is what leslie wants). ;
c$$$
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$       
      implicit none
      include 'omp_lib.h'
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_omp_sub__in
      integer n_S_0_sub__in
      integer n_S_1_sub__in
      integer n_M_0_sub__in
      integer n_M_1_sub__in
      logical flag_ctf
      logical flag_tru
      logical flag_fig
      logical flag_MS_vs_SM
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: n_polar_a_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_k_cur,n_r,nr,n_w_max,n_A,na
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
      real *8 half_diameter_x_c,half_diameter_k_c
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
c$$$      Templates S_
      integer n_S_max,n_S,ld_S,ns,nS_sample
      complex *16, allocatable :: S_x_c__(:)
      complex *16, allocatable :: S_x_p__(:)
      complex *16, allocatable :: S_k_c__(:)
      complex *16, allocatable :: S_k_p__(:)
      complex *16, allocatable :: S_k_q__(:)
      integer *4, allocatable :: I_S_sample_(:)
c$$$      Images N_ and M_
      integer n_M_max,n_M,ld_M,nm,nM_sample
      complex *16, allocatable :: N_x_c__(:)
      complex *16, allocatable :: N_x_p__(:)
      complex *16, allocatable :: N_k_c__(:)
      complex *16, allocatable :: N_k_p__(:)
      complex *16, allocatable :: N_k_q__(:)
      complex *16, allocatable :: M_k_p__(:)
      integer *4, allocatable :: I_M_sample_(:)
c$$$      Temporary Array T?_
      complex *16, allocatable :: T0_k_p_(:)
      complex *16, allocatable :: T1_k_p_(:)
      complex *16, allocatable :: T2_k_p_(:)
      complex *16, allocatable :: T3_k_p_(:)
c$$$      CTF-functions
      integer *4 n_ctf,ld_CTF,nctf
      integer *4, allocatable :: I_ctf_(:)
      complex *16, allocatable :: CTF_k_p__(:)
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4
c$$$      fftw_plan_many 
      integer fpm_howmany_max
c$$$      image parameters 
      include 'nalpha_define.f'
      real *8 l2_norm
      real *8, allocatable :: alpha_tru__(:)
      real *8, allocatable :: alpha_est__(:)
      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: S_alpha_S_index_(:)
      real *8, allocatable :: S_alpha_polar_a_(:)
      real *8, allocatable :: S_alpha_azimu_b_(:)
      real *8 alpha_update_f
      parameter (alpha_update_f=0.0d0)
c$$$      true displacement, rotation and scaling parameters for each image
      real *8, allocatable :: polar_a_tru_(:)
      real *8, allocatable :: azimu_b_tru_(:)
      real *8, allocatable :: delta_x_tru_(:)
      real *8, allocatable :: delta_y_tru_(:)
      real *8, allocatable :: gamma_z_tru_(:)
      real *8, allocatable :: l2_norm_tru_(:)
      real *8, allocatable :: S_index_tru_(:)
c$$$      estimated (current) displacement, rotation and scaling parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: S_index_est_(:)
c$$$      updated displacement, rotation and scaling parameters for each image
      real *8, allocatable :: polar_a_upd_(:)
      real *8, allocatable :: azimu_b_upd_(:)
      real *8, allocatable :: delta_x_upd_(:)
      real *8, allocatable :: delta_y_upd_(:)
      real *8, allocatable :: gamma_z_upd_(:)
      real *8, allocatable :: l2_norm_upd_(:)
      real *8, allocatable :: S_index_upd_(:)
c$$$      errors in displacement, rotation and scaling parameters for each image
      real *8, allocatable :: polar_a_err_(:)
      real *8, allocatable :: azimu_b_err_(:)
      real *8, allocatable :: delta_x_err_(:)
      real *8, allocatable :: delta_y_err_(:)
      real *8, allocatable :: gamma_z_err_(:)
      real *8, allocatable :: l2_norm_rer_(:)
      real *8, allocatable :: S_index_err_(:)
c$$$      range of displacement and rotation parameters to measure
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 displacement_max,angle_tmp,delta_tmp,dtmp
      integer n_delta_x,n_delta_y,n_delta_v,ndx,ndy,ndv
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_pre,delta_y_pre,gamma_z_pre
      real *8 delta_x_pos,delta_y_pos,gamma_z_pos
      real *8 eps_svd,N_pixels_in
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
      integer *4, allocatable :: n_MS_(:) ! array of size n_S indicating the actual number of images stored per tempalte. ;
      real *8, allocatable :: alpha_MS__(:) ! array of size n_alpha*n_MS_max*n_S storing the image-parameters for each stored image-template pair. ;
      character(len=1024) fname,format_string,prefix_string
c$$$      function for template generation
      external get_F2_x_c
      real *8 param_1,param_2
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      random seed
      integer *4 rseed !random seed. ;
c$$$      
      logical flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
      real *8 d_memory_estimate !real *8: estimate of required memory (in bytes). ;
      complex *16 C_X
      integer *4 MDA_n_d,MDA_d_(0:127)
      character(len=1024) MDA_string
c$$$      
      rseed = 1
      pi = 4*atan(1.0)
c$$$      
c$$$      List of parameters: ;
      verbose = 1
      ncur = 20
      eps_svd = 0.00001d0
      N_pixels_in = 2.5d0*dsqrt(2.0d0) ! Note that this should match one of the svd-libraries. ;
      N_pixels_in = 4.0d0 ! e.g., we have a library for N_pixels = 4.0d0. ;
      n_delta_x = 7
      n_delta_y = 7
      svd_calculation_type = 2
      flag_RTRT_vs_RTTR = .false.
      n_omp_sub__in = 8
      n_S_0_sub__in = 2
      n_S_1_sub__in = 3
      n_M_0_sub__in = 5
      n_M_1_sub__in = 2
c$$$ /*%%%%%%%%*/
      n_omp_sub__in = 1
      n_S_0_sub__in = 1
      n_S_1_sub__in = 1
      n_M_0_sub__in = 1
      n_M_1_sub__in = 1
c$$$ /*%%%%%%%%*/
      flag_ctf = .true.
      flag_tru = .false.
      n_S_max = 12*1
      n_S_max = 1*1
      n_S = 1*1
      n_M_max = 17*1
      n_M_max = 1*1
      n_M = 1*1
      fpm_howmany_max = 16
      flag_MS_vs_SM = .false.
      n_SM_max = 45
      n_MS_max = 30
      flag_memory_estimate = .false.
      tesselation_distance_req = 0.25
c$$$      tesselation_distance_req = 2.0d0
      n_LT_add = max(1,n_M/10) ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      n_LT_ref = max(1,n_SM_max/2) ! number of image-template pairs to consider when refining local search.
      
      n_gamma_z = ncur
      n_polar_a = max(6,nint(pi*ncur))
      n_w_max = n_polar_a*2
      n_gamma_z = n_w_max
      fpm_howmany_max = max(1,fpm_howmany_max)
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_dr_7]: '
       end if
       if (verbose.gt.0) then
          write(6,'(A,I0)') ' openblas_set_num_threads: ' , 1
       end if !if (verbose.gt.0) then
       call openblas_set_num_threads(1)
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
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' ,flag_RTRT_vs_RTTR
         write(6,'(A,I0)') ' n_omp_sub__in: ' ,n_omp_sub__in
         write(6,'(A,I0)') ' n_S_0_sub__in: ' ,n_S_0_sub__in
         write(6,'(A,I0)') ' n_S_1_sub__in: ' ,n_S_1_sub__in
         write(6,'(A,I0)') ' n_M_0_sub__in: ' ,n_M_0_sub__in
         write(6,'(A,I0)') ' n_M_1_sub__in: ' ,n_M_1_sub__in
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
      n_k_cur = n_r
      ld_S = n_A
      ld_M = n_A
      ld_CTF = n_A
      allocate(T0_k_p_(0:n_A-1));
      allocate(T1_k_p_(0:n_A-1));
      allocate(T2_k_p_(0:n_A-1));
      allocate(T3_k_p_(0:n_A-1));

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
      if (verbose.gt.1) then
         write(6,'(A,A,F8.3)') ' fftw_plan:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))
      
c$$$      generating grids for templates and images 
      max_x_c = 1.0d0
      half_diameter_x_c = max_x_c/2.0d0
      allocate(grid_x_c_(0:2*n_r-1))
      dtmp = max_x_c/max(1,2*n_r-1)
      call linspace(0.0d0,max_x_c+dtmp,2*n_r,grid_x_c_)
      allocate(grid_x_p_(0:n_r-1))
      dtmp = half_diameter_x_c/n_r
      call linspace(0.0d0+dtmp,half_diameter_x_c+dtmp,n_r,grid_x_p_)
      max_k_c = (1.0d0*2*n_r)/max_x_c
      half_diameter_k_c = max_k_c/2.0d0
      allocate(grid_k_c_(0:2*n_r-1))
      call linspace(0.0d0,max_k_c,2*n_r,grid_k_c_)
      allocate(grid_k_p_(0:n_r-1))
      dtmp = half_diameter_k_c/n_r
      call linspace(0.0d0+dtmp,half_diameter_k_c+dtmp,n_r,grid_k_p_)
      if (verbose.gt.0) then
         call write_sub_r8(2*n_r,grid_x_c_,12,' grid_x_c_: ')
         call write_sub_r8(1*n_r,grid_x_p_,12,' grid_x_p_: ')
         call write_sub_r8(2*n_r,grid_k_c_,12,' grid_k_c_: ')
         call write_sub_r8(1*n_r,grid_k_p_,12,' grid_k_p_: ')
      end if !if (verbose.gt.0) then
c$$$      defining displacement_max to be 3.0d0 wavelengths (can change later)
      displacement_max = 3.0d0/max_k_c

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
      end if
      allocate(delta_x_(0:ceiling(pi*n_delta_x*n_delta_y)-1))
      allocate(delta_y_(0:ceiling(pi*n_delta_x*n_delta_y)-1))
      call get_delta_1(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,n_delta_y,n_delta_v,delta_x_,delta_y_)
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
         call write_sub_r8(n_delta_v,delta_x_,11,' delta_x_: ')
         call write_sub_r8(n_delta_v,delta_y_,11,' delta_y_: ')
         call write_sub_r8(n_gamma_z,gamma_z_,11,' gamma_z_: ')
      end if !if (verbose.gt.1) then

c$$$      generate templates
      allocate(I_S_sample_(0:n_S-1))
      do ns=0,n_S-1
         I_S_sample_(ns)=max(0,min(n_S_max-1,nint(n_S_max*(1.0d0*ns)
     $        /(1.0d0*n_S))))
      enddo !do ns=0,n_S-1
      allocate(S_alpha_S_index_(0:n_S-1))
      allocate(S_alpha_polar_a_(0:n_S-1))
      allocate(S_alpha_azimu_b_(0:n_S-1))
      do ns=0,n_S-1
         S_alpha_S_index_(ns) = ns ! array storing the numbers 0:n_S-1, from 0 up to the number of templates used - 1. ;
         S_alpha_polar_a_(ns) = pi*(1.0d0*ns + 0.5)/(1.0d0*n_S)
         S_alpha_azimu_b_(ns) = 2*pi*(1.0d0*ns + 0.5)/(1.0d0*n_S)
      enddo
      if (flag_memory_estimate.eqv..false.) then
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') ' Here we construct the templates "S_".'
         write(6,'(A,I0,A)') ' We generate n_S = ',n_S
     $        ,' templates in total,'
         write(6,*) 'each using the function get_F2_x_c'
     $        ,' with different parameters.'
         write(6,*) 'For each template: '
         write(6,*)
     $        ' The array S_x_c__ holds the x-space values on a'
     $        ,' regular cartesian grid defined by grid_x_c_.'
         write(6,*)
     $        ' The array S_k_c__ holds the k-space values on a'
     $        ,' regular cartesian grid defined by grid_k_c_.'
         write(6,*)
     $        ' The array S_k_p__ holds the k-space values on a'
     $        ,' quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
         write(6,*) ' The array S_k_q__ holds the "hat-hat"'
     $        ,' (i.e, bessel-expansion) for the k-space values'
         write(6,*) ' on a quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
      end if !if (verbose.gt.1) then
      allocate(S_x_c__(0:2*n_r*2*n_r*n_S_max-1))
      allocate(S_x_p__(0:n_A*n_S_max-1))
      allocate(S_k_c__(0:2*n_r*2*n_r*n_S_max-1))
      allocate(S_k_p__(0:n_A*n_S_max-1))
      allocate(S_k_q__(0:n_A*n_S_max-1))
      do ns=0,n_S-1
         nS_sample = I_S_sample_(ns)
         param_1 = nS_sample*dsqrt(2.5d0)
         call periodize_r8(param_1,-1.0d0,1.0d0,param_1)
         param_2 = 2*(nS_sample+2)
c$$$      Calculating x-space template on regular cartesian-grid
         call get_F2_x_c_(2*n_r,grid_x_c_,max_x_c,2*n_r,grid_x_c_
     $        ,max_x_c,S_x_c__(nS_sample*2*n_r*2*n_r),get_F2_x_c,param_1
     $        ,param_2)
c$$$      Calculating x-space template on quasi-uniform polar-grid
         call interp_c_to_p(2*n_r,max_x_c,2*n_r,max_x_c
     $        ,S_x_c__(nS_sample*2*n_r*2*n_r),n_r,grid_x_p_,n_w_,n_A
     $        ,S_x_p__(nS_sample*n_A))
c$$$      Calculating k-space template on regular cartesian-grid
         call adi_fft2(-1,2*n_r,2*n_r,S_x_c__(nS_sample*2*n_r*2*n_r)
     $        ,S_k_c__(nS_sample*2*n_r*2*n_r))
c$$$      Calculating k-space template on quasi-uniform polar-grid
         call interp_c_to_p(2*n_r,max_k_c,2*n_r,max_k_c
     $        ,S_k_c__(nS_sample*2*n_r*2*n_r),n_r,grid_k_p_,n_w_,n_A
     $        ,S_k_p__(nS_sample*n_A))
c$$$      Calculating J-space template on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,S_k_p__(nS_sample*n_A),S_k_q__(nS_sample*n_A))
      enddo !do ns=0,n_S-1
      end if !if (flag_memory_estimate.eqv..false.) then

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$% Notes: test rotate_p2p_fx with the following fortran code: ;
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      ns=0
c$$$      nS_sample = I_S_sample_(ns)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call cp1_c16(n_A,S_k_p__(nS_sample*n_A),T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T0_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      gamma_z = pi/3;
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p_to_p(n_r,n_w_,n_A,T0_k_p_,+gamma_z,T1_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T1_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T1_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fz(n_r,n_w_,n_A,T0_k_p_,+gamma_z,T2_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T2_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T2_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
c$$$     $     ,n_A,fftw_in1_,fftw_out_,T0_k_p_,+gamma_z,T3_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_T3_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T3_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p_to_p(n_r,n_w_,n_A,T1_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U1_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fz(n_r,n_w_,n_A,T2_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U2_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
c$$$     $     ,n_A,fftw_in1_,fftw_out_,T3_k_p_,-gamma_z,T0_k_p_)
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_A
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_U3_.mda'
c$$$      call MDA_write_c16(MDA_n_d,MDA_d_,T0_k_p_,MDA_string)
c$$$c$$$      %%%%%%%%%%%%%%%%
c$$$      MDA_n_d = 1
c$$$      MDA_d_(0) = n_r
c$$$      write(MDA_string,'(A)') './dir_mda/tmp_n_w_.mda'
c$$$      call MDA_write_i4(MDA_n_d,MDA_d_,n_w_,MDA_string)
c$$$      goto 10
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$% Notes: then run the following matlab code.; 
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$
c$$$n_w_ = MDA_read_i4('dir_mda/tmp_n_w_.mda'); 
c$$$n_r = length(n_w_); grid_p_ = 1:n_r; n_A = sum(n_w_);
c$$$T0_ = MDA_read_c16('dir_mda/tmp_T0_.mda');
c$$$T1_ = MDA_read_c16('dir_mda/tmp_T1_.mda');
c$$$T2_ = MDA_read_c16('dir_mda/tmp_T2_.mda');
c$$$T3_ = MDA_read_c16('dir_mda/tmp_T3_.mda');
c$$$U1_ = MDA_read_c16('dir_mda/tmp_U1_.mda');
c$$$U2_ = MDA_read_c16('dir_mda/tmp_U2_.mda');
c$$$U3_ = MDA_read_c16('dir_mda/tmp_U3_.mda');
c$$$
c$$$disp_flag=0;
c$$$if disp_flag;
c$$$%%%%%%%%%%%%%%%%;
c$$$figure(1);
c$$$pcols = 4;
c$$$np=1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T0_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T0_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T0_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T0_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(T3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(T3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(T3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(T3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$end;%if disp_flag;
c$$$
c$$$E1_ = T0_ - U1_;
c$$$E2_ = T0_ - U2_;
c$$$E3_ = T0_ - U3_;
c$$$
c$$$disp_flag=1;
c$$$if disp_flag;
c$$$%%%%%%%%%%%%%%%%;
c$$$figure(1);
c$$$pcols = 3;
c$$$np=1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E1_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E1_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E2_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E2_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$subplot(2,pcols,np+0*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,real(E3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('real(E3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$subplot(2,pcols,np+1*pcols); clim = imagesc_p(n_r,grid_p_,n_w_,n_A,imag(E3_),[],colormap(colormap_pm(64)));
c$$$xlim(n_r*[-1,1]);ylim(n_r*[-1,1]);
c$$$title(sprintf('imag(E3_) [%0.5f,%0.5f]',min(clim),max(clim)),'Interpreter','none');
c$$$set(gca,'XTick',[],'YTick',[]); axis square;
c$$$np=np+1;
c$$$%%%%%%%%%%%%%%%%;
c$$$end;%if disp_flag;

      allocate(I_M_sample_(0:n_M-1))
      do nm=0,n_M-1
         I_M_sample_(nm)=max(0,min(n_M_max-1,nint(n_M_max*(1.0d0*nm)
     $        /(1.0d0*n_M))))
      enddo !do nm=0,n_M-1
      if (flag_memory_estimate.eqv..false.) then
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
     $        ,' will be stored in the S_index variable.'
         write(6,*) 'For this example we pick delta_x_tru_(j)'
     $        ,' to be (j%3)-1,',' and delta_y_tru_(j) to be'
     $        ,' 2*((j%3)-1).'
         write(6,*) 'We also pick gamma_z_tru_(j) to be'
     $        ,' 2*pi/3*(j%3).'
      end if !if (verbose.gt.1) then
      allocate(N_x_c__(0:2*n_r*2*n_r*n_M_max-1))
      allocate(N_k_c__(0:2*n_r*2*n_r*n_M_max-1))
      allocate(N_k_p__(0:n_A*n_M_max-1))
      allocate(N_k_q__(0:n_A*n_M_max-1))
c$$$      Calculating x-space image on regular cartesian-grid
      allocate(polar_a_tru_(0:n_M-1))
      allocate(azimu_b_tru_(0:n_M-1))
      allocate(delta_x_tru_(0:n_M-1))
      allocate(delta_y_tru_(0:n_M-1))
      allocate(gamma_z_tru_(0:n_M-1))
      allocate(l2_norm_tru_(0:n_M-1))
      allocate(S_index_tru_(0:n_M-1))
      do nm=0,n_M-1
         nM_sample = I_M_sample_(nm)
         call periodize_i(nm,0,n_S,ns)
         nS_sample = I_S_sample_(ns)
         ndv = max(0,min(n_delta_v-1,mod(nm,n_delta_v)))
         ngz = max(0,min(n_gamma_z-1,mod(13*nm+2,n_gamma_z)))
         polar_a_tru_(nm) = S_alpha_polar_a_(ns)
         azimu_b_tru_(nm) = S_alpha_azimu_b_(ns)
         delta_x_tru_(nm) = delta_x_(ndv)
         delta_y_tru_(nm) = delta_y_(ndv)
         gamma_z_tru_(nm) = gamma_z_(ngz)
         l2_norm_tru_(nm) = 1.0d0 + 1*nm
         S_index_tru_(nm) = S_alpha_S_index_(ns)
         if (verbose.gt.1) then
            write(6,'(I0,A,4F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_tru_(nm),delta_y_tru_(nm),gamma_z_tru_(nm)
     $           ,l2_norm_tru_(nm)
         end if !if (verbose.gt.1) then
c$$$         In the absence of interpolation errors,
c$$$         these transformations will produce the N_k_p__ used below.
         call cp1_c16(2*n_r*2*n_r,S_x_c__(nS_sample*2*n_r*2*n_r)
     $        ,N_x_c__(nM_sample*2*n_r*2*n_r))
         call transl_c_to_c(2*n_r,max_x_c,2*n_r,max_x_c
     $        ,N_x_c__(nM_sample*2*n_r*2*n_r),+delta_x_tru_(nm),
     $        +delta_y_tru_(nm),N_x_c__(nM_sample*2*n_r*2*n_r))
         call rotate_c_to_c(2*n_r,max_x_c,2*n_r,max_x_c
     $        ,N_x_c__(nM_sample*2*n_r*2*n_r),+gamma_z_tru_(nm)
     $        ,N_x_c__(nM_sample*2*n_r*2*n_r))
c$$$      Calculating k-space image on regular cartesian-grid      
         call adi_fft2(-1,2*n_r,2*n_r,N_x_c__(nM_sample*2*n_r*2*n_r)
     $        ,N_k_c__(nM_sample*2*n_r*2*n_r))
c$$$      Calculating k-space image on quasi-uniform polar-grid
         call interp_c_to_p(2*n_r,max_k_c,2*n_r,max_k_c
     $        ,N_k_c__(nM_sample*2*n_r*2*n_r),n_r,grid_k_p_,n_w_,n_A
     $        ,N_k_p__(nM_sample*n_A))
c$$$      Calculating J-space image on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,N_k_p__(nM_sample*n_A),N_k_q__(nM_sample*n_A))
c$$$         However, due to interpolation errors, we will generate
c$$$         the N_k_p__ in k-space.
         call cp1_c16(n_A,S_k_p__(nS_sample*n_A),N_k_p__(nM_sample*n_A))
         call transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,N_k_p__(nM_sample
     $        *n_A),+delta_x_tru_(nm),+delta_y_tru_(nm)
     $        ,N_k_p__(nM_sample*n_A))
         call rotate_p2p_fz(n_r,n_w_,n_A,N_k_p__(nM_sample*n_A),
     $        +gamma_z_tru_(nm),N_k_p__(nM_sample*n_A))
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,N_k_p__(nM_sample*n_A),N_k_q__(nM_sample*n_A))
      enddo !do nm=0,n_M-1
      end if !if (flag_memory_estimate.eqv..false.) then

      if (flag_memory_estimate.eqv..false.) then
      allocate(polar_a_est_(0:n_M-1))
      allocate(azimu_b_est_(0:n_M-1))
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      allocate(S_index_est_(0:n_M-1))
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
         polar_a_est_(nm) = polar_a_tru_(nm) 
         azimu_b_est_(nm) = azimu_b_tru_(nm) 
         delta_x_est_(nm) = delta_x_tru_(nm) 
         delta_y_est_(nm) = delta_y_tru_(nm) 
         gamma_z_est_(nm) = gamma_z_tru_(nm) 
         l2_norm_est_(nm) = l2_norm_tru_(nm) 
         S_index_est_(nm) = S_index_tru_(nm) 
      enddo
      else !if (flag_tru) then
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ''
c$$$         write(6,*) 'For this test we choose these estimated'
c$$$     $        ,' translations and rotations to be zero.'
c$$$      end if
c$$$      do nm=0,n_M-1
c$$$         polar_a_est_(nm) = 0.0d0
c$$$         azimu_b_est_(nm) = 0.0d0
c$$$         delta_x_est_(nm) = 0.0d0
c$$$         delta_y_est_(nm) = 0.0d0
c$$$         gamma_z_est_(nm) = 0.0d0
c$$$         l2_norm_est_(nm) = 1.0d0
c$$$         S_index_est_(nm) = 0.0d0
c$$$      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'For this test we choose these estimated'
     $        ,' translations and rotations to be discrete shifts'
     $        ,' away from the true translations and rotations.'
      end if
      do nm=0,n_M-1
         polar_a_est_(nm) = polar_a_tru_(nm) 
         azimu_b_est_(nm) = azimu_b_tru_(nm) 
         l2_norm_est_(nm) = 1.0d0
         S_index_est_(nm) = 0.0d0
         ndv=  max(0,min(n_delta_v-1,mod(nm,n_delta_v)))
         ngz = max(0,min(n_gamma_z-1,mod(nm+5,7)))
         if (verbose.gt.2) then
            write(6,'(3(A,I0))') ' setting alpha_est: nm ' , nm ,
     $           ' ndv ' ,ndv , ' ngz ' , ngz 
         end if !if (verbose.gt.2) then
         if (flag_RTRT_vs_RTTR.eqv..true.) then
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            gamma_z = gamma_z_(ngz)
            gamma_z_est_(nm) = gamma_z_tru_(nm) - gamma_z
            delta_x_pre = delta_x_tru_(nm) - delta_x
            delta_y_pre = delta_y_tru_(nm) - delta_y
            call get_interchange_delta(delta_x_pre,delta_y_pre
     $           ,-gamma_z,delta_x_pos,delta_y_pos)
            delta_x_est_(nm) = delta_x_pos
            delta_y_est_(nm) = delta_y_pos
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
         if (flag_RTRT_vs_RTTR.eqv..false.) then
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            gamma_z = gamma_z_(ngz)
            gamma_z_est_(nm) = gamma_z_tru_(nm) - gamma_z
            delta_x_pre = delta_x_tru_(nm)
            delta_y_pre = delta_y_tru_(nm)
            call get_interchange_delta(delta_x_pre,delta_y_pre
     $           ,-gamma_z,delta_x_pos,delta_y_pos)
            delta_x_est_(nm) = delta_x_pos - delta_x
            delta_y_est_(nm) = delta_y_pos - delta_y
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
      enddo
      end if ! if (flag_tru) then
      end if !if (flag_memory_estimate.eqv..false.) then

      n_ctf = 3
      if (flag_memory_estimate.eqv..false.) then
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') ' Now we generate some ctf-functions.'
         write(6,'(A)')
     $        ' These are deliberately chosen to be anisotropic.'
      end if
      allocate(CTF_k_p__(0:n_ctf*n_A-1))
c$$$      Here we use an unrealistic ctf-function with a high degree
c$$$      of anisotropy to test our code.
      do nctf=0,n_ctf-1
         ctf_p_0 = max_k_c/6.0
         ctf_p_1 = (2*pi*nctf + 0.25*pi)/n_ctf
         ctf_p_2 = 2
         if (flag_ctf) then
         call get_ctf_star_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p__(nctf
     $        *n_A),ctf_p_0 ,ctf_p_1,ctf_p_2)
         else
         call get_ctf_ones_k_p_(n_r,grid_k_p_,n_w_,n_A,CTF_k_p__(nctf
     $        *n_A))
         end if ! if (flag_ctf) then
      enddo
      if (flag_fig) then
         write(6,'(A,I0,A)')
     $        ' Trying to generate Figure dir_jpg/Fig_CTF_p_',n_r,'.jpg'
         write(fname,'(A,I0)') './dir_jpg/Fig_CTF_p_',n_r
         call Fig_gen_ctf_ver0(n_r,n_polar_a_,grid_k_p_,n_ctf,ld_CTF
     $        ,CTF_k_p__,fname)
      end if
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         Determine which images correspond to which ctf-function.
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_ctf_(0:n_M-1))
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo
      end if !if (flag_memory_estimate.eqv..false.) then

      if (verbose.gt.0) then
         write(6,'(A)') ' Setting displacement_max to max_x_c '
      end if ! if (verbose.gt.0) then
      displacement_max = max_x_c
      if (flag_memory_estimate.eqv..false.) then
      allocate(alpha_tru__(0:n_alpha*n_M-1))
      call cl1_r8(n_alpha*n_M,alpha_tru__)
      allocate(alpha_est__(0:n_alpha*n_M-1))
      call cl1_r8(n_alpha*n_M,alpha_est__)
      do nm=0,n_M-1
         alpha_tru__(nalpha_polar_a + n_alpha*nm) = polar_a_tru_(nm)
         alpha_tru__(nalpha_azimu_b + n_alpha*nm) = azimu_b_tru_(nm)
         alpha_tru__(nalpha_gamma_z + n_alpha*nm) = gamma_z_tru_(nm)
         alpha_tru__(nalpha_delta_x + n_alpha*nm) = delta_x_tru_(nm)
         alpha_tru__(nalpha_delta_y + n_alpha*nm) = delta_y_tru_(nm)
         alpha_tru__(nalpha_l2_norm + n_alpha*nm) = l2_norm_tru_(nm)
         alpha_tru__(nalpha_ctf_ind + n_alpha*nm) = 1.0d0*I_ctf_(nm)
         alpha_tru__(nalpha_S_index + n_alpha*nm) = S_index_tru_(nm)
         alpha_tru__(nalpha_M_index + n_alpha*nm) = 1.0d0*nm
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
      allocate(M_k_p__(0:n_A*n_M_max-1))
      allocate(C_M_(0:n_M-1))
      call cl1_c16(n_M,C_M_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         nM_sample = I_M_sample_(nm)
         nctf = I_ctf_(nm)
         call xx1_c16(n_A,N_k_p__(nM_sample*n_A),CTF_k_p__(nctf *n_A)
     $        ,M_k_p__(nM_sample*n_A))
      enddo                     ! do nm=0,n_M-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         nM_sample = I_M_sample_(nm)
         l2_norm = alpha_tru__(nalpha_l2_norm+n_alpha*nm)
         call af1_c16(n_A,1.0d0*cmplx(l2_norm,0.0d0),1.0d0 *cmplx(0.0d0
     $        ,0.0d0),M_k_p__(nM_sample*n_A) ,M_k_p__(nM_sample*n_A))
      enddo                     ! do nm=0,n_M-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(6,'(A,I0,A)')
     $        ' Trying to generate Figure dir_jpg/Fig_MNC_p_',n_r,'.jpg'
         write(fname,'(A,I0)') './dir_jpg/Fig_MNC_p_',n_r
         call Fig_gen_ctf_ver4(n_r,n_polar_a_,grid_k_p_,n_M,I_M_sample_
     $        ,ld_M,N_k_p__,M_k_p__,ld_CTF,CTF_k_p__,alpha_est__
     $        ,min(6,n_M) ,fname)
      end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  allocate storage for image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..true.) then
      allocate(n_MS_(0:n_S-1))
      call cl1_i4(n_S,n_MS_)
      allocate(alpha_MS__(0:n_alpha*n_MS_max*n_S-1))
      call cl1_r8(n_alpha*n_MS_max*n_S,alpha_MS__)
      end if !if (flag_MS_vs_SM.eqv..true.) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  allocate storage for template-image pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..false.) then
      allocate(n_SM_(0:n_M-1))
      call cl1_i4(n_M,n_SM_)
      allocate(alpha_SM__(0:n_alpha*n_SM_max*n_M-1))
      call cl1_r8(n_alpha*n_SM_max*n_M,alpha_SM__)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      end if !if (flag_memory_estimate.eqv..false.) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts via brute force
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_memory_estimate.eqv..false.) then
      if (verbose.gt.2) then
      do ns=0,n_S-1
         nS_sample = I_S_sample_(ns)
         do nm=0,n_M-1
            nM_sample = I_M_sample_(nm)
            nctf = I_ctf_(nm)
            call test_innerproduct_bruteforce_1(verbose ,n_r ,grid_k_p_
     $           ,n_w_ ,half_diameter_x_c,S_k_p__(n_A*nS_sample)
     $           ,M_k_p__(n_A *nM_sample) ,CTF_k_p__(n_A*nctf) ,
     $           delta_x_est_(nm) ,delta_y_est_(nm) ,gamma_z_est_(nm)
     $           ,n_delta_v,delta_x_,delta_y_ ,n_gamma_z
     $           ,gamma_z_)
         enddo !do nm=0,n_M-1
      enddo !do ns=0,n_S-1
      end if !if (verbose.gt.2) then
      end if !if (flag_memory_estimate.eqv..false.) then

      verbose = 1
      call print_sub_c16(n_A,S_k_p__,' S_k_p__: ')
      call print_sub_c16(n_A,M_k_p__,' M_k_p__: ')
      call print_all_r8(n_alpha,alpha_est__,'alpha_est__: ')

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $     ,n_M,')-x-(n_S = ',n_S, ') = (', n_M
     $     *n_S,') innerproducts'
      timing_tic = omp_get_wtime()
      call test_innerproduct_7(
     $     verbose !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
     $     ,grid_k_p_ !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_rq determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_ !real *8 array (length at least n_S): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,n_M !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $     ,I_M_sample_ !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $     ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y associated with reconstructed molecule. 
     $     ,n_CTF !integer *4: total number of CTF functions. ;
     $     ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
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
     $     ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_cur' under consideration (which can change from iteration to iteration). ;
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
     $     ,n_omp_sub__in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
     $     ,n_S_0_sub__in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,n_S_1_sub__in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,n_M_0_sub__in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,n_M_1_sub__in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M): actual l2-norm (out to frequency n_k_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
      timing_toc = omp_get_wtime()
      if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
         write(6,'(A,A,F8.3)') 'test_innerproduct_7:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then

      if (flag_memory_estimate.eqv..true.) then
         write(6,'(A,2(I0,A))') ' d_memory estimate: ' ,
     $        nint(d_memory_estimate*1.0d-6) , ' (MB); ' , 
     $        nint(d_memory_estimate*1.0d-9) , ' (GB); ' 
         goto 10
      end if !if (flag_memory_estimate.eqv..true.) then


c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  update alpha_est__
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_MS_vs_SM.eqv..false.) then 
         if (verbose.gt.1) then
         write(6,'(A)') ' alpha_SM: '
         do nm=0,n_M-1
            write(prefix_string,'(A,I2,A)') ' alpha_SM: nm:', nm , ' '
            call alpha_SM_write_0(n_SM_max,n_SM_(nm),alpha_SM__(n_alpha
     $           *n_SM_max*nm),len_trim(prefix_string)+1,prefix_string)
         enddo !do nm=0,n_M-1
         end if !if (verbose.gt.1) then
         do nm=0,n_M-1
            call test_alpha_update_SM_3(verbose,rseed,n_SM_max,n_SM_(nm)
     $           ,alpha_SM__(n_alpha*n_SM_max*nm),flag_RTRT_vs_RTTR
     $           ,alpha_est__(n_alpha*nm),alpha_update_f
     $           ,alpha_est__(n_alpha*nm))
         enddo !do nm=0,n_M-1
      end if !if (flag_MS_vs_SM.eqv..false.) then 
      if (flag_MS_vs_SM.eqv..true.) then 
         if (verbose.gt.1) then
         write(6,'(A)') ' alpha_MS: '
         do ns=0,n_S-1
            write(prefix_string,'(A,I2,A)') ' alpha_MS: ns:', ns , ' '
            call alpha_SM_write_0(n_MS_max,n_MS_(ns),alpha_MS__(n_alpha
     $           *n_MS_max*ns),len_trim(prefix_string)+1,prefix_string)
         enddo !do ns=0,n_S-1
         end if !if (verbose.gt.1) then
         call test_alpha_update_MS_3(verbose,rseed,n_MS_max,n_MS_,n_S
     $        ,n_M,alpha_MS__,flag_RTRT_vs_RTTR,alpha_est__,alpha_est__)
      end if !if (flag_MS_vs_SM.eqv..true.) then 

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(6,'(A,I0,A)')
     $        ' Trying to generate Figure dir_jpg/Fig_MTS_c_',n_r,'.jpg'
         write(fname,'(A,I0)') './dir_jpg/Fig_MTS_c_',n_r
         call Fig_gen_ver8(n_r,n_polar_a_,grid_k_p_,n_S,I_S_sample_,n_A
     $        ,S_k_p__,n_M,I_M_sample_,n_A,M_k_p__ ,n_A,CTF_k_p__
     $        ,alpha_est__ ,min(16,n_M),fname)
         goto 10
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we update the estimated displacement'
     $        ,' and rotation parameters.'
      end if
      allocate(polar_a_upd_(0:n_M-1))
      allocate(azimu_b_upd_(0:n_M-1))
      allocate(delta_x_upd_(0:n_M-1))
      allocate(delta_y_upd_(0:n_M-1))
      allocate(gamma_z_upd_(0:n_M-1))
      allocate(l2_norm_upd_(0:n_M-1))
      allocate(S_index_upd_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_upd_(nm) = alpha_est__(nalpha_polar_a + n_alpha*nm)
         azimu_b_upd_(nm) = alpha_est__(nalpha_azimu_b + n_alpha*nm)
         delta_x_upd_(nm) = alpha_est__(nalpha_delta_x + n_alpha*nm)
         delta_y_upd_(nm) = alpha_est__(nalpha_delta_y + n_alpha*nm)
         gamma_z_upd_(nm) = alpha_est__(nalpha_gamma_z + n_alpha*nm)
         l2_norm_upd_(nm) = alpha_est__(nalpha_l2_norm + n_alpha*nm)
         S_index_upd_(nm) = alpha_est__(nalpha_S_index + n_alpha*nm)
      enddo !do nm=0,n_M-1
      if (verbose.gt.0 .and. n_M.le.25) then
         write(6,'(A)') ''
         write(6,*)
     $        'Now we compare updated estimate to true parameters: '
         write(format_string,'(A,A,A)')
     $        '(3(A,I3,1X),6(A,F8.3,1X))'
         do nm=0,n_M-1
            nM_sample = I_M_sample_(nm)
            write(6,format_string) 
     $           ' nm: ' , nm 
     $           , ' nM_sample: ' , nM_sample
     $           , ' S_index_tru_: ' , nint(S_index_tru_(nm))
     $           , ' polar_a_tru_: ' , polar_a_tru_(nm) 
     $           , ' azimu_b_tru_: ' , azimu_b_tru_(nm) 
     $           , ' delta_x_tru_: ' , delta_x_tru_(nm) 
     $           , ' delta_y_tru_: ' , delta_y_tru_(nm) 
     $           , ' gamma_z_tru_: ' , gamma_z_tru_(nm) 
     $           , ' l2_norm_tru_: ' , l2_norm_tru_(nm) 
            write(6,format_string) 
     $           ' nm: ' , nm 
     $           , ' nM_sample: ' , nM_sample
     $           , ' S_index_upd_: ' , nint(S_index_upd_(nm))
     $           , ' polar_a_upd_: ' , polar_a_upd_(nm) 
     $           , ' azimu_b_upd_: ' , azimu_b_upd_(nm) 
     $           , ' delta_x_upd_: ' , delta_x_upd_(nm) 
     $           , ' delta_y_upd_: ' , delta_y_upd_(nm) 
     $           , ' gamma_z_upd_: ' , gamma_z_upd_(nm) 
     $           , ' l2_norm_upd_: ' , l2_norm_upd_(nm) 
         enddo !do nm=0,n_M-1
      end if !if (verbose.gt.0 .and. n_M.le.25) then
      allocate(polar_a_err_(0:n_M-1))
      allocate(azimu_b_err_(0:n_M-1))
      allocate(delta_x_err_(0:n_M-1))
      allocate(delta_y_err_(0:n_M-1))
      allocate(gamma_z_err_(0:n_M-1))
      allocate(l2_norm_rer_(0:n_M-1))
      allocate(S_index_err_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_err_(nm) = dabs(polar_a_tru_(nm)-polar_a_upd_(nm))
         azimu_b_err_(nm) = dabs(azimu_b_tru_(nm)-azimu_b_upd_(nm))
         delta_x_err_(nm) = dabs(delta_x_tru_(nm)-delta_x_upd_(nm))
         delta_y_err_(nm) = dabs(delta_y_tru_(nm)-delta_y_upd_(nm))
         gamma_z_err_(nm) = dabs(gamma_z_tru_(nm)-gamma_z_upd_(nm))
         call periodize_r8(gamma_z_err_(nm),-pi,+pi,gamma_z_err_(nm))
         l2_norm_rer_(nm) = dabs(l2_norm_tru_(nm)-l2_norm_upd_(nm))
     $        /l2_norm_tru_(nm)
         S_index_err_(nm) = dabs(S_index_tru_(nm)-S_index_upd_(nm))
      enddo !do nm=0,n_M-1
      if (verbose.gt.0 .and. n_M.le.25) then
         write(6,*) 'Alignment Errors:'
         write(format_string,'(A)')
     $        '(3(A,I3,1X),6(A,F8.6,1X))'
         do nm=0,n_M-1
            nM_sample = I_M_sample_(nm)
            write(6,format_string) 
     $           ' nm: ' , nm 
     $           , ' nM_sample: ' , nM_sample 
     $           , ' S_index_: ' , nint(S_index_err_(nm))
     $           , ' polar_a_: ' , polar_a_err_(nm)
     $           , ' azimu_b_: ' , azimu_b_err_(nm)
     $           , ' delta_x_: ' , delta_x_err_(nm)
     $           , ' delta_y_: ' , delta_y_err_(nm)
     $           , ' gamma_z_: ' , gamma_z_err_(nm)
     $           , ' l2_norm_: ' , l2_norm_rer_(nm)
         enddo !do nm=0,n_M-1
      end if !if (verbose.gt.0 .and. n_M.le.25) then
      if (verbose.gt.0) then
         write(6,*) 'Average Alignment Errors:'
         write(format_string,'(A)')
     $        '(6(A,F8.6,1X),A,I0)'
         write(6,format_string) 
     $    ' polar_a_err_:', avg_r8_f(n_M,polar_a_err_),
     $    ' azimu_b_err_:', avg_r8_f(n_M,azimu_b_err_),
     $    ' delta_x_err_:', avg_r8_f(n_M,delta_x_err_), 
     $    ' delta_y_err_:', avg_r8_f(n_M,delta_y_err_), 
     $    ' gamma_z_err_:', avg_r8_f(n_M,gamma_z_err_),
     $    ' l2_norm_rer_:', avg_r8_f(n_M,l2_norm_rer_),
     $    ' S_index_err_:', nint(avg_r8_f(n_M,S_index_err_))
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,*) '|Alignment Errors|_{2}:'
         write(format_string,'(A)')
     $        '(6(A,F8.6,1X),A,I0)'
         write(6,format_string) 
     $    ' polar_a_err_:', al2_r8_f(n_M,polar_a_err_),
     $    ' azimu_b_err_:', al2_r8_f(n_M,azimu_b_err_),
     $    ' delta_x_err_:', al2_r8_f(n_M,delta_x_err_), 
     $    ' delta_y_err_:', al2_r8_f(n_M,delta_y_err_), 
     $    ' gamma_z_err_:', al2_r8_f(n_M,gamma_z_err_),
     $    ' l2_norm_rer_:', al2_r8_f(n_M,l2_norm_rer_),
     $    ' S_index_err_:', nint(al2_r8_f(n_M,S_index_err_))
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,*) 'Largest Alignment Errors:'
         write(format_string,'(A)')
     $        '(6(A,F8.6,1X),A,I0)'
         write(6,format_string) 
     $    ' polar_a_err_:', max_r8_f(n_M,polar_a_err_),
     $    ' azimu_b_err_:', max_r8_f(n_M,azimu_b_err_),
     $    ' delta_x_err_:', max_r8_f(n_M,delta_x_err_), 
     $    ' delta_y_err_:', max_r8_f(n_M,delta_y_err_), 
     $    ' gamma_z_err_:', max_r8_f(n_M,gamma_z_err_),
     $    ' l2_norm_rer_:', max_r8_f(n_M,l2_norm_rer_),
     $    ' S_index_err_:', nint(max_r8_f(n_M,S_index_err_))
      end if !if (verbose.gt.0) then

 10    continue

      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_innerproduct_dr_7]'
      end if

      stop
      end
