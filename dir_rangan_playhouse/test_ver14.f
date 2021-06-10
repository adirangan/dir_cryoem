c$$$      Driver for molecular reconstruction. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Generally speaking, this driver first generates two molecules ;
c$$$      (i.e., molecule-A and molecule-B) which are rather similar. ;
c$$$      Using these two molecules, we produce a set of 'slices', some from ;
c$$$      molecule-A, and some from molecule-B. ;
c$$$      We then use a set of (typically 3) anisotropic CTF-functions to ;
c$$$      distort these 'slices', producing 'images'. ;
c$$$      We also distort these images by translating and rotating. ;
c$$$      Using these images, along with (randomly chosen) initial image-parameters, ;
c$$$      we try and reconstruct the molecule out to some maximum radius ;
c$$$      in k-space polar coordinates (denoted by n_k_p_max). ;
c$$$      Along the way we try and estimate the cross-correlations ;
c$$$      between residuals, allowing us to classify the images (i.e., either A or B). ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      implicit none
      include 'omp_lib.h'
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Here are several parameters that ;
c$$$      determine what the driver will do. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer verbose ! verbosity level. ;
      data verbose / 1 /
      integer n_k_p_max ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      parameter (n_k_p_max=40) ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      integer n_x_c_max ! maximum number of x-values (in x-space cartesian coordinates) used to generate initial molecule. sometimes named ngrid. ;
      parameter (n_x_c_max=60) ! maximum number of x-values (in x-space cartesian coordinates) used to generate initial molecule. sometimes named ngrid. ;
      logical flag_displacement ! flag determining whether or not displacements will be considered when matching images to templates (and subsequently used in least-squares calculation). ;
      parameter (flag_displacement=.false.)
      integer n_delta_x_set ! if displacements are considered, this value determines the number of displacements considered in the x-direction (in the plane). ;
      parameter (n_delta_x_set=5)
      integer n_delta_y_set ! if displacements are considered, this value determines the number of displacements considered in the y-direction (in the plane). ;
      parameter (n_delta_y_set=5)
      real *8 n_pixels_in_set ! if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength under consideration (which can change from iteration to iteration). ;
      parameter (n_pixels_in_set=1.5d0)
      real *8 displacement_max_set ! if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
      parameter (displacement_max_set=0.1d0)
      integer svd_calculation_type ! flag determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == -1 --> svd_calculation_type chosen automatically from 0,1,2 based on number of translations and image dimensions. ;
c$$$      svd_calculation_type == 0 --> multiply, then fft (see test_innerproduct_batch_stage_2.f). ;
c$$$      svd_calculation_type == 1 --> fft, then multiply (see test_innerproduct_batch_stage_2.f). ;
c$$$      svd_calculation_type == 2 --> brute force displacements (see test_innerproduct_batch_stage_2.f). ;
      parameter (svd_calculation_type=1)
      real *8 eps_svd ! svd tolerance epsilon, typically 1. ;
      parameter (eps_svd=10.0d0)
      logical flag_normalization ! flag determining whether or not l2_norm will be calculated when matching images to templates (and subsequently used in least-squares calculation). ;
      parameter (flag_normalization=.false.)
      logical flag_ctf ! flag determining whether or not ctf function will be different from unity (e.g., an anisotropic ctf). ;
      parameter (flag_ctf=.false.)
      logical flag_tesselation ! flag determining whether or not to use adaptive sampling of templates when calculating innerproducts for each image. ;
      parameter (flag_tesselation=.false.)
      real *8 tesselation_distance_req !temporary: if adaptively sampling templates, this value determines the maximum distance (on the unit-sphere) allowed between the viewing-angle of any particular template and the estimated-viewing-angle of a given image. If this value is large (e.g., 2.0d0), then every template will be considered for each image. If, on the other hand, this value is small (e.g., <0.125d0), then only a few templates will be considered for each image. When adaptively sampling templates, we might consider this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
      logical flag_heterogeneity !flag determining whether or not to include multiple (i.e., 3) molecules in image pool. ;
      parameter (flag_heterogeneity=.false.)
      real *8 r_heterogeneity_A !fraction of heterogeneous components: molecules A,B,C and D. see command line inputs. ;
      real *8 r_heterogeneity_B !fraction of heterogeneous components: molecules A,B,C and D. see command line inputs. ;
      real *8 r_heterogeneity_C !fraction of heterogeneous components: molecules A,B,C and D. see command line inputs. ;
      real *8 r_heterogeneity_D !fraction of heterogeneous components: molecules A,B,C and D. see command line inputs. ;
      integer n_heterogeneity_A !number of images: molecules A,B,C and D
      integer n_heterogeneity_B !number of images: molecules A,B,C and D
      integer n_heterogeneity_C !number of images: molecules A,B,C and D
      integer n_heterogeneity_D !number of images: molecules A,B,C and D
      logical flag_residual !flag determining whether or not to calculate image residuals (and classify / cluster images). ;
      parameter (flag_residual=.false.)
      integer n_k_1 ! lowest value of nk (k-index in k-space polar coordinates) to use when calculating residuals. ;
      parameter (n_k_1=1)
      integer n_residual_loading ! number of loadings to track when analyzing residuals. ;
      parameter (n_residual_loading=2)
      integer n_residual_iteration ! number of iterations to perform when analyzing residuals. ;
      parameter (n_residual_iteration=7)
      logical flag_fig !flag determining whether or not to generate figures (in dir_base)
      parameter (flag_fig=.false.)
      integer n_M_sample_max ! maximum number of images to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
      parameter (n_M_sample_max=32)
      integer n_S_sample_max    ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
      logical flag_snr !flag determining whether or not to add iid gaussian noise to images. ;
      parameter (flag_snr=.true.)
      real *8 snr !actual snr used. see command line inputs. ;
c$$$      parameter (snr=10.0d0)
      parameter (n_S_sample_max=16)
      integer n_omp_sub ! number of batches to use when finding innerproducts
      parameter (n_omp_sub=8)
      integer fpm_howmany_max ! Maximum number of fftws to call simultaneously within the fftw_plan_many. 
      parameter (fpm_howmany_max=128)
      real *8 alpha_update_f ! fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.).
      parameter (alpha_update_f=0.05d0)
      logical flag_MS_vs_SM !determines whether to assign images to templates (.true.) or templates to images (.false.). ;
      parameter (flag_MS_vs_SM=.true.)
      logical flag_skip !determines whether to start iteration at n_k_low = n_k_max (.true.) or 3 (.false.). ;
      parameter (flag_skip=.false.)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Here are other variables are ;
c$$$      internally generated ;
c$$$      or temporary. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer npolar_a,nazimu_b,nA !temporary indices for polar_a, azimu_b and azimu_b_polar_a. ;
      integer n_k_p_lowerbound ! lower bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;
      integer fft_iflag !temporary: stores direction of fft (+1 or -1). sometimes named iflag. ;
      integer ier !temporary: stores error flag (integer). ;
      integer quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
      integer quadrature_type_radial ! quadrature type in the radial direction. sometimes named ityper. ;
      integer lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
      integer n_image_size ! total number of entries in each image (in k-space polar coordinates). sometimes called ld_M or nimagesize. ;
      integer ld_M ! total number of entries in each image (in k-space polar coordinates). i.e., leading-dimension of M array. sometimes called n_image_size or nimagesize. ;
      integer n_image ! total number of images. sometimes called n_M or nimages. ;
      integer n_M ! total number of images. sometimes called n_image or nimages. ;
      integer n_k_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_max. ;
      integer n_k_bot ! lowest value of nk (k-index in k-space polar coordinates) to use when reconstructing molecule. Note that this runs from 1 to n_k_p_max. ;
      integer n_x1_c_x2_c_x3_c_sum ! total number of points in x-space cartesian coordinates. sometimes called npts. ;
      integer n_azimu_b_polar_a_sum_tmp !temporary: total number of points on sphere at a particular radius grid_k_p_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
      integer n_Y_lm_sum_max ! total number of points used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
      integer ld_S ! total number of entries in each template (in k-space polar coordinates). sometimes called ntemplate_size or ntemplatesize. ;
      integer n_S ! total number of templates. sometimes called ntemplates. ;
      integer n_azimu_b_polar_a_k_p_sum ! total number of points on all spheres (out to radius grid_k_p_(n_k_p_max-1)). sometimes called ntot. ;
      integer n_azimu_b_polar_a_sum_cur !temporary: n_azimu_b_polar_a_sum_(n_k_cur-1). sometimes called numonsphere_cur. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      real *8 pi ! 4.0d0*atan(1.0d0) ;
      real *8 half_diameter_x_c ! half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
      real *8 x_c_step ! spacing of grid in x-space cartesian coordinates. sometimes called h. ;
      real *8 diameter_x_c ! diameter of particle support in x-space cartesian grid. sometimes called boxsize. ;
      real *8 cos_polar_a !temporary: cos(polar angle). ;
      real *8 eps_default !temporary: tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
      real *8 f_k_c_l2_norm !temporary: l2-norm of f in k-space cartesian coordinates. sometimes called ffhatnorm. ;
      real *8 f_x_c_l2_norm !temporary: l2-norm of f in x-space cartesian coordinates. sometimes called fnorm. ;
      real *8 azimu_b !temporary: azimuthal angle. sometimes called phi. ;
      real *8 azimu_b_step !temporary: grid spacing for azimuthal angle. sometimes called phistep. ;
      real *8 k_p_max ! maximum value of k in k-space polar coordinates. sometimes called rmax. ;
      real *8 sin_polar_a !temporary: sin(polar angle). ;
      real *8 lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer n_k_cur ! current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
      integer nk ! current index for k-value grid_k_p_(nk). ;
      parameter (n_k_p_lowerbound=4*n_k_p_max) ! lower bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;
      integer n_source_max ! maximum number of gaussian sources for synthetic molecule. ;
      parameter (n_source_max=128)
      integer n_source_use ! actual number of gaussian sources for synthetic molecule. ;
      integer nsource !temporary; indexes n_source_use. ;
      integer n_source_cut_B ! nsource at which molecule B differs from molecule A. ;
      integer n_source_cut_C ! nsource at which molecule C differs from molecule A. ;
      real *8 source_sigma ! standard-deviation of gaussian sources for synthetic molecule. ;
      real *8 source_r ! radial parameter for gaussian source location for synthetic molecule. (x-space cylindrical coordinates). ;
      real *8 source_w ! angular parameter for gaussian source location for synthetic molecule. (x-space cylindrical coordinates). ;
      real *8 source_z ! height parameter for gaussian source location for synthetic molecule. (x-space cylindrical coordinates). ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer, allocatable :: nlats(:) ! number of latitudes per sphere, as a function of sphere index (i.e., nk). also called n_polar_a_(nk). ;
      integer, allocatable :: n_polar_a_(:) ! number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
      integer, allocatable :: numonsphere(:) ! total number of points on sphere at radius grid_k_p_(nk). also called n_azimu_b_polar_a_sum_(nk). ;
      integer, allocatable :: n_azimu_b_polar_a_sum_(:) ! total number of points on sphere at radius grid_k_p_(nk). also called numonsphere(nk). ;
      integer, allocatable :: nterms_sph(:) ! order of spherical harmonic expansion at each radius grid_k_p_(nk). also called n_Y_l_(nk). ;
      integer, allocatable :: n_Y_l_(:) ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
      integer, allocatable :: isph_start(:) ! cumulative sum of n_Y_lm_. also called n_Y_lm_sum_(nk). ;
      integer, allocatable :: n_Y_lm_sum_(:) ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
      integer, allocatable :: ngridc(:) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called n_w_(nk). ;
      integer, allocatable :: n_w_(:) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer, allocatable :: ngridc_tmp(:) !temporary: similar to ngridc. ;
      integer, allocatable :: icstart(:) ! cumulative sum of n_w_. also called n_w_sum_(nk). ;
      integer, allocatable :: n_w_sum_(:) ! cumulative sum of n_w_. also called icstart(nk). ;
      integer, allocatable :: ngridps(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_(nk). also called n_azimu_b_(np). ;
      integer, allocatable :: n_azimu_b_(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_(nk). also called ngridps(np). ;
c$$$      integer, allocatable :: n_azimu_b_tmp_(:) !temporary: similar to n_azimu_b_. ;
      real *8, allocatable :: grid_k_p_(:) ! values for k on successive shells for k-space polar coordinates. sometimes called xnodesr(nk). ;
      real *8, allocatable :: weight_k_p_(:) ! weight associated with grid_k_p_(nk) in k-space polar coordinates. sometimes called wtsr(nk). ;
      real *8, allocatable :: xnodesth(:) !temporary: values for cos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called grid_cos_polar_a_(nk). ;
      real *8, allocatable :: grid_cos_polar_a_(:) !temporary: values for cos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called xnodesth(nk). ;
      real *8, allocatable :: sthetas(:) !temporary: values for sin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called grid_sin_polar_a_(nk). ;
      real *8, allocatable :: grid_sin_polar_a_(:) !temporary: values for sin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called sthetas(nk). ;
      real *8, allocatable :: weight_cos_polar_a_(:) !temporary: weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called wtsth_(np). ;
      real *8, allocatable :: azimu_b_step_(:) !temporary: grid-spacing for azimu_b. ;
c$$$      real *8, allocatable :: grid_k_p_tmp_(:) ! similar to grid_k_p_. ;
c$$$      real *8, allocatable :: grid_cos_polar_a_tmp_(:) ! similar to grid_cos_polar_a_. 
c$$$      real *8, allocatable :: grid_sin_polar_a_tmp_(:) ! similar to grid_sin_polar_a_. ;
c$$$      real *8, allocatable :: weight_cos_polar_a_tmp_(:) ! similar to weight_cos_polar_a_. ;
c$$$      real *8, allocatable :: azimu_b_step_tmp_(:) !temporary: similar to azimu_b_step_. ;
      real *8, allocatable :: source_x1_c_x2_c_x3_c_A__(:,:) ! array containing (x1,x2,x3) coordinates of gaussian sources for molecule-A in x-space cartesian coordinates. ;
      real *8, allocatable :: source_x1_c_x2_c_x3_c_B__(:,:) ! array containing (x1,x2,x3) coordinates of gaussian sources for molecule-B in x-space cartesian coordinates. ;
      real *8, allocatable :: source_x1_c_x2_c_x3_c_C__(:,:) ! array containing (x1,x2,x3) coordinates of gaussian sources for molecule-C in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x1_c___(:,:,:) ! array containing x-values for points in 3d regular grid in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x2_c___(:,:,:) ! array containing y-values for points in 3d regular grid in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x3_c___(:,:,:) ! array containing z-values for points in 3d regular grid in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_A___(:,:,:) ! function corresponding to molecule-A in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_B___(:,:,:) ! function corresponding to molecule-B in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_C___(:,:,:) ! function corresponding to molecule-C in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_tru___(:,:,:) !temporary: function obtained using true-angles to reconstruct molecule in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_est___(:,:,:) !temporary: function obtained using estimated-angles to reconstruct molecule in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_k1_c_(:) ! values for k1 of (k1,k2,k3) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      real *8, allocatable :: grid_k2_c_(:) ! values for k2 of (k1,k2,k3) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      real *8, allocatable :: grid_k3_c_(:) ! values for k3 of (k1,k2,k3) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      real *8, allocatable :: weight_k_c_(:) ! weight for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      complex *16, allocatable :: f_k_c_A_(:) ! function corresponding to molecule-A in k-space cartesian coordinates evaluted on nested collection of spherical shells described by grid_k?_c_. ;
      complex *16, allocatable :: f_k_c_B_(:) ! function corresponding to molecule-B in k-space cartesian coordinates evaluted on nested collection of spherical shells described by grid_k?_c_. ;
      complex *16, allocatable :: f_k_c_C_(:) ! function corresponding to molecule-C in k-space cartesian coordinates evaluted on nested collection of spherical shells described by grid_k?_c_. ;
      complex *16, allocatable :: f_k_c_tru_(:) !temporary: function obtained using true-angles to reconstruct molecule in k-space cartesian coordinates evaluated on nested collection of spherical shells described by grid_k?_c_. ;
      complex *16, allocatable :: f_k_c_est_(:) !temporary: function obtained using estimated-angles to reconstruct molecule in k-space cartesian coordinates evaluated on nested collection of spherical shells described by grid_k?_c_. ;
      complex *16, allocatable :: Y_A_(:) ! spherical harmonic coefficients for molecule-A. sometimes called modsph_A. ;
      complex *16, allocatable :: Y_B_(:) ! spherical harmonic coefficients for molecule-B. sometimes called modsph_B. ;
      complex *16, allocatable :: Y_C_(:) ! spherical harmonic coefficients for molecule-C. sometimes called modsph_C. ;
      complex *16, allocatable :: Y_tru_(:) ! spherical harmonic coefficients obtained using true-angles to reconstruct molecule. sometimes called modsph_tru. ;
      complex *16, allocatable :: Y_est_(:) ! spherical harmonic coefficients obtained using estimated-angles to reconstruct molecule. sometimes called modsph_est. ;
      complex *16, allocatable :: Y_slice_A__(:,:) ! images associated with molecule-A. These are taken directly from the molecular function itself (with no noise or l2-norm or ctf). sometimes called cslices_A. ;
      complex *16, allocatable :: Y_slice_B__(:,:) ! images associated with molecule-A. These are taken directly from the molecular function itself (with no noise or l2-norm or ctf). sometimes called cslices_A. ;
      complex *16, allocatable :: Y_slice_C__(:,:) ! images associated with molecule-A. These are taken directly from the molecular function itself (with no noise or l2-norm or ctf). sometimes called cslices_A. ;
      complex *16, allocatable :: Y_slice__(:,:) ! mixed stack of images associated with molecules-A,B,C. These are taken directly from the molecular function itself (with no noise or l2-norm or ctf). sometimes called cslices. ; These will later be modified (in place) by multiplying by the CTF and/or modifying the l2-norm (depending on which flags are set above). ;
      integer *4, allocatable :: I_model_(:) ! indexes which images correspond to which molecule (i.e., A, B or C). ;
      integer *4, allocatable :: I_model_sample_(:) ! indexes which of the sampled images correspond to which molecule (i.e., A, B or C). ;
c$$$      complex *16, allocatable :: S__(:,:) !temporary: stack of templates associated with reconstructed molecule. sometimes called templates. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer tmp_d_(0:5)       !temporary: dimensions used for MDA_write. ;
      character(len=1024) str_k,str_displacement,str_svd
     $     ,str_normalization,str_ctf,str_tesselation,str_heterogeneity
     $     ,str_residual,str_sample,str_snr,str_MS_vs_SM,str_rseed
      character(len=1024) dir_base,str_system !fixed name based on driver parameters. ;
      character(len=1024) tmp_dname,tmp_fname,format_string !temporary: strings. ;
      real *8 half_diameter_x_c_tmp ! half-diameter in x-space cartesian grid. sometimes called max_x_c. ;
      parameter(half_diameter_x_c_tmp=1.0d0)
      integer n_gamma_z !temporary: value used for n_gamma_z when calculating innerproducts. ;
      integer n_delta_x,n_delta_y !temporary: values used for n_delta_x and n_delta_y, respectively, when calculating innerproducts involving displacements. ;
      integer *4 n_SM_use !temporary: total number of image-template pairs considered when calculating innerproducts. Might be smaller than the total number when calculating innerproducts adaptively. ;
      parameter(n_SM_use=100)
      real *8 n_pixels_in !temporary: value used for n_pixels_in when calculating innerproducts involving displacements. ;
      real *8 displacement_max !temporary: value used for displacement_max when calculating innerproducts involving displacements. ;
      integer n_M_sample,nM_sample,nm !temporary: values used to track number of sampled images. ;
      integer n_M_A_sample,n_M_B_sample,n_M_C_sample,n_M_D_sample !temporary: values used to track number of sampled images from molecule-A,B,C. ;
      integer n_Y_lm_sum_cur !temporary: total number of points used in spherical harmonic expansion (in k-space polar coordinates) out to a particular radius k. sometimes called n_Y_lm_sum_cur
      integer n_S_sample,nS_sample,ns !temporary: values used to track number of sampled templates. ;
      integer *4, allocatable :: I_M_sample_(:) !temporary: indexing variable used to reference images. ;
      integer *4, allocatable :: I_S_sample_(:) !temporary: indexing variable used to reference templates. ;
c$$$      complex *16, allocatable :: Y_slice_sample_(:) !temporary: sampled images. These are taken directly from the image stack Y_slice__ (with no noise or l2-norm or ctf). ;
c$$$      complex *16, allocatable :: M_sample_(:) !temporary: sampled images. These are modified from Y_slice_sample_ by multiplying by the CTF and/or modifying the l2-norm (depending on which flags are set above). ;
c$$$      complex *16, allocatable :: M_residual_transform__(:) !temporary: residuals = images - slices from current model. ;
c$$$      complex *16, allocatable :: M_residual_slice__(:) !temporary: residuals = images - slices from current model. ;
      complex *16, allocatable :: M_residual_sample__(:) !temporary: residuals = images - slices from current model. ;
      complex *16, allocatable :: M_residual_loading_(:) !temporary: loading vectors associated with residuals. ;
      complex *16, allocatable :: S_k_p__(:) !temporary: stack of templates associated with reconstructed molecule. sometimes called templates. ;
      integer *4 n_alpha ! number of image parameters (i.e., 'alphas') per image. defined in 'nalpha_define.f'. ;
      parameter (n_alpha=8)
      include 'nalpha_define.f'
      real *8 tmp_l2_norm !temporary: estimated l2_norm of an image. ;
      real *8, allocatable :: alpha2d_tru__(:,:) ! true image parameters associated with original image set (e.g., stored in Y_slice_) stored in 2-dimensional array. ;
      real *8, allocatable :: alpha1d_tru_(:) ! true image parameters associated with sampled image set (e.g., could be stored in Y_slice_sample_ or M_sample_) stored in 1-dimensional array. ;
      real *8, allocatable :: alpha2d_est__(:,:) !temporary: estimated image parameters associated with original image set (e.g., stored in Y_slice_) stored in 2-dimensional array. ;
      real *8, allocatable :: alpha1d_est_(:) ! estimated image parameters associated with sampled image set (e.g., could be stored in Y_slice_sample_ or M_sample_) stored in 1-dimensional array. ;
      real *8, allocatable :: alpha1d_bkp_(:) ! estimated image parameters associated with sampled image set from previous iteration stored in 1-dimensional array. ;
      real *8 timing_tic,timing_toc !temporary: timing variables. ;
      external fgauss ! function to describe a gaussian in x-space cartesian coordinates. ;
      integer *4 n_ctf ! stores the number of CTF functions (i.e., 3)
      integer ld_CTF ! total number of entries in each CTF (in k-space polar coordinates). i.e., leading-dimension of CTF array. Typically expected to be at least ld_M. ;
      integer *4 nctf !temporary: indexes ctf. ;
      integer *4, allocatable :: I_ctf_(:) ! indexes which ctf is associated with which image. ;
      complex *16, allocatable :: CTF_k_p_(:) ! stores the n_ctf different CTF functions in k-space polar coordinates. ;
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4 ! parameters associated with CTF generation. ; Note that the CTF we use is strongly anisotropic (cartoonishly so). ;
      integer *4 rseed !random seed. see command line inputs. ;
      real *8 tmp_v,adi_rand_f !temporary: random real number between 0 and 1. ;
      real *8 tmp_al1,al1_c16_f !temporary: l1-norm of slice. ;
      real *8 tmp_al2,tmp_al2s,al2_c16_f !temporary: l2-norm of slice. ;
      real *8 tmp_sigma !temporary: sigma used to add noise. ;
c$$$      command line input variables
      character(len=64) :: cmd_argstring
      integer ncmd
      real *8 cmd_A,cmd_B,cmd_C,cmd_D,cmd_snr
      integer *4 cmd_rseed

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_ver14]'
      end if

      n_delta_x = 11
      n_delta_y = 11
      n_gamma_z = max(32,n_w_(n_k_p_max))
      displacement_max = 0.1d0
      n_pixels_in = 1.5d0
      n_M = 100000
      n_CTF = n_M/500
      call test_innerproduct_memory_estimate_0(verbose+2,n_omp_sub
     $     ,n_k_p_max,n_M,n_CTF,n_alpha,n_pixels_in
     $     ,displacement_max,n_delta_x ,n_delta_y,n_gamma_z
     $     ,svd_calculation_type ,eps_svd ,fpm_howmany_max,n_SM_use)
      goto 10

      ncmd = 0
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_A
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_B
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_C
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_D
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_snr
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_rseed
      ncmd = ncmd + 1
      
      if (verbose.gt.0) then
         write(6,'(A,F8.4)') ' cmd_A: ' , cmd_A
         write(6,'(A,F8.4)') ' cmd_B: ' , cmd_B
         write(6,'(A,F8.4)') ' cmd_C: ' , cmd_C
         write(6,'(A,F8.4)') ' cmd_D: ' , cmd_D
         write(6,'(A,F8.4)') ' cmd_snr: ' , cmd_snr
         write(6,'(A,I0)') ' cmd_rseed: ' , cmd_rseed
      end if !if (verbose.gt.0) then

      r_heterogeneity_A = max(0.0d0,cmd_A)
      r_heterogeneity_B = max(0.0d0,cmd_B)
      r_heterogeneity_C = max(0.0d0,cmd_C)
      r_heterogeneity_D = max(0.0d0,cmd_D)
      if (r_heterogeneity_A.le.0.0d0 .and. r_heterogeneity_B.le.0.0d0
     $     .and. r_heterogeneity_C.le.0.0d0 .and.
     $     r_heterogeneity_D.le.0.0d0) r_heterogeneity_A = 1.0d0
      snr = max(0.0d0,cmd_snr)
      rseed = max(0,cmd_rseed)

      tmp_al1 = abs(r_heterogeneity_A) + abs(r_heterogeneity_B) +
     $     abs(r_heterogeneity_C) + abs(r_heterogeneity_D)
      r_heterogeneity_A = abs(r_heterogeneity_A)/tmp_al1
      r_heterogeneity_B = abs(r_heterogeneity_B)/tmp_al1
      r_heterogeneity_C = abs(r_heterogeneity_C)/tmp_al1
      r_heterogeneity_D = abs(r_heterogeneity_D)/tmp_al1

c$$$      %%%%%%%%%%%%%%%%
      write(str_k,'(A,I0)') '_k' , n_k_p_max
c$$$      %%%%%%%%%%%%%%%%
      if (flag_displacement) then
         write(str_displacement,'(4(A,I0))') '_dTx' , n_delta_x_set ,
     $        'y' ,n_delta_y_set , 'p' , nint(100*n_pixels_in_set) , '
     $        m' ,nint(100*displacement_max_set)
      else
         write(str_displacement,'(A)') '_dF'
      end if !if (flag_displacement) then
c$$$      %%%%%%%%%%%%%%%%
      if (svd_calculation_type.eq.2) then
         write(str_svd,'(A)') '_s2'
      else 
         write(str_svd,'(A)') '_s1'
      end if !if (svd_calculation_type.eq.2) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_normalization) then
         write(str_normalization,'(A)') '_nT'
      else
         write(str_normalization,'(A)') '_nF'
      end if !if (flag_normalization) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_ctf) then
         write(str_ctf,'(A)') '_cT'
      else
         write(str_ctf,'(A)') '_cF'
      end if !if (flag_ctf) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_tesselation) then
         write(str_tesselation,'(A)') '_tT'
      else
         write(str_tesselation,'(A)') '_tF'
      end if !if (flag_tesselation) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_heterogeneity) then
         write(str_heterogeneity,'(4(A,I0))') '_hTa' , nint(100
     $        *r_heterogeneity_A) , 'b' , nint(100*r_heterogeneity_B) ,
     $        'c' , nint(100*r_heterogeneity_C) , 'd' , nint(100
     $        *r_heterogeneity_D)
      else
         write(str_heterogeneity,'(A)') '_hF'
      end if !if (flag_heterogeneity) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_residual) then
         write(str_residual,'(A,I0,A,I0)') '_rTl' , n_residual_loading ,
     $        'i', n_residual_iteration
      else
         write(str_residual,'(A)') '_rF'
      end if !if (flag_residual) then
c$$$      %%%%%%%%%%%%%%%%
      write(str_sample,'(A,I0,A,I0)') '_M' , n_M_sample_max , '_S' ,
     $     n_S_sample_max
c$$$      %%%%%%%%%%%%%%%%
      write(str_snr,'(A,I0)') '_s' , nint(100*snr)
c$$$      %%%%%%%%%%%%%%%%
      if (flag_MS_vs_SM) then
         write(str_MS_vs_SM,'(A)') '_MS'
      else
         write(str_MS_vs_SM,'(A,I0)') '_SMa' , nint(100*alpha_update_f)
      end if !if (flag_MS_vs_SM) then
c$$$      %%%%%%%%%%%%%%%%
      write(str_rseed,'(A,I0)') , '_r' , rseed
c$$$      %%%%%%%%%%%%%%%%
      write(dir_base,'(13(A))') './dir' ,trim(str_k)
     $     ,trim(str_displacement),trim(str_svd),trim(str_normalization)
     $     ,trim(str_ctf),trim(str_tesselation),trim(str_heterogeneity)
     $     ,trim(str_residual),trim(str_sample),trim(str_snr)
     $     ,trim(str_MS_vs_SM),trim(str_rseed)

      write(6,'(A,A)') ' dir_base: ' , trim(dir_base)
      write(str_system,'(A,A)') 'mkdir ' , trim(dir_base)
      call system(str_system)

c-----------------------------------------------------------------------
c     0) allocate some arrays (dynamic allocation required for openmp)
c-----------------------------------------------------------------------
      allocate(nlats(n_k_p_max))
      allocate(n_polar_a_(0:n_k_p_max-1))
      allocate(numonsphere(n_k_p_max))
      allocate(n_azimu_b_polar_a_sum_(0:n_k_p_max-1))
      allocate(nterms_sph(n_k_p_max))
      allocate(n_Y_l_(0:n_k_p_max-1))
      allocate(isph_start(n_k_p_max))
      allocate(n_Y_lm_sum_(0:n_k_p_max-1))
      allocate(ngridc(n_k_p_lowerbound))
      allocate(n_w_(0:n_k_p_lowerbound-1))
      allocate(ngridc_tmp(n_k_p_lowerbound))
      allocate(icstart(n_k_p_max))
      allocate(n_w_sum_(0:n_k_p_max-1))
      allocate(ngridps(n_k_p_lowerbound))
      allocate(n_azimu_b_(0:n_k_p_lowerbound-1))
      allocate(grid_k_p_(0:n_k_p_max-1))
      allocate(weight_k_p_(0:n_k_p_max-1))
      allocate(xnodesth(n_k_p_lowerbound))
      allocate(grid_cos_polar_a_(0:n_k_p_lowerbound-1))
      allocate(sthetas(n_k_p_lowerbound))
      allocate(grid_sin_polar_a_(0:n_k_p_lowerbound-1))
      allocate(weight_cos_polar_a_(0:n_k_p_lowerbound-1))
      allocate(azimu_b_step_(0:n_k_p_lowerbound-1))
c$$$      allocate(n_azimu_b_tmp_(0:n_k_p_lowerbound-1))
c$$$      allocate(grid_cos_polar_a_tmp_(0:n_k_p_lowerbound-1))
c$$$      allocate(grid_sin_polar_a_tmp_(0:n_k_p_lowerbound-1))
c$$$      allocate(weight_cos_polar_a_tmp_(0:n_k_p_lowerbound-1))
c$$$      allocate(azimu_b_step_tmp_(0:n_k_p_lowerbound-1))
      allocate(source_x1_c_x2_c_x3_c_A__(0:3-1,0:n_source_max-1))
      if (flag_heterogeneity) then
         allocate(source_x1_c_x2_c_x3_c_B__(0:3-1,0:n_source_max-1))
         allocate(source_x1_c_x2_c_x3_c_C__(0:3-1,0:n_source_max-1))
      end if !if (flag_heterogeneity) then
      allocate(grid_x1_c___(n_x_c_max,n_x_c_max,n_x_c_max))
      allocate(grid_x2_c___(n_x_c_max,n_x_c_max,n_x_c_max))
      allocate(grid_x3_c___(n_x_c_max,n_x_c_max,n_x_c_max))
      allocate(f_x_c_A___(n_x_c_max,n_x_c_max,n_x_c_max))
      if (flag_heterogeneity) then
         allocate(f_x_c_B___(n_x_c_max,n_x_c_max,n_x_c_max))
         allocate(f_x_c_C___(n_x_c_max,n_x_c_max,n_x_c_max))
      end if !if (flag_heterogeneity) then
      allocate(f_x_c_tru___(n_x_c_max,n_x_c_max,n_x_c_max))
      allocate(f_x_c_est___(n_x_c_max,n_x_c_max,n_x_c_max))

      pi=4.0d0*atan(1.0d0)
      call prini(6,13)

c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
      n_source_use = min(n_source_max,32)
      n_source_cut_B = floor(0.825d0*n_source_use)
      n_source_cut_C = floor(0.750d0*n_source_use)
      source_sigma = 1.0d0/32.0d0 
      do nsource=0,n_source_use-1
         source_r = 0.5d0*(nsource)/(n_source_use-1)
         source_w = 2*pi*2*nsource/(n_source_use)
         source_z = 1.0d0*(nsource)/(n_source_use-1) - 0.50d0
         source_x1_c_x2_c_x3_c_A__(1-1,nsource) = source_r*cos(source_w)
         source_x1_c_x2_c_x3_c_A__(2-1,nsource) = source_r*sin(source_w)
         source_x1_c_x2_c_x3_c_A__(3-1,nsource) = source_z
         if (flag_heterogeneity) then
            if (nsource.lt.n_source_cut_B) then
               source_x1_c_x2_c_x3_c_B__(1-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(1-1,nsource)
               source_x1_c_x2_c_x3_c_B__(2-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(2-1,nsource)
               source_x1_c_x2_c_x3_c_B__(3-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(3-1,nsource)
            else
               source_r = 0.5d0*(n_source_cut_B + 3.5d0*(nsource
     $              -n_source_cut_B))/(n_source_use-1)
               source_w = 2*pi*2*(n_source_cut_B - 0.15d0*(nsource
     $              -n_source_cut_B)*(nsource-n_source_use+1))
     $              /(n_source_use)
               source_z = 1.0d0*(n_source_cut_B - 0.50d0*(nsource
     $              -n_source_cut_B)*(nsource-n_source_use+1))
     $              /(n_source_use -1) - 0.50d0
               source_x1_c_x2_c_x3_c_B__(1-1,nsource) = source_r
     $              *cos(source_w)
               source_x1_c_x2_c_x3_c_B__(2-1,nsource) = source_r
     $              *sin(source_w)
               source_x1_c_x2_c_x3_c_B__(3-1,nsource) = source_z
            end if !if (nsource.lt.n_source_cut_B) then
            if (nsource.lt.n_source_cut_C) then
               source_x1_c_x2_c_x3_c_C__(1-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(1-1,nsource)
               source_x1_c_x2_c_x3_c_C__(2-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(2-1,nsource)
               source_x1_c_x2_c_x3_c_C__(3-1,nsource) =
     $              source_x1_c_x2_c_x3_c_A__(3-1,nsource)
            else
               source_r = 0.5d0*(n_source_cut_C + 1.5d0*(nsource
     $              -n_source_cut_C))/(n_source_use-1)
               source_w = 2*pi*2*(n_source_cut_C - 0.15d0*(nsource
     $              -n_source_cut_C)*(nsource-n_source_use+1))
     $              /(n_source_use)
               source_z = 1.0d0*(n_source_cut_C - 0.25d0*(nsource
     $              -n_source_cut_C)*(nsource-n_source_use+1))
     $              /(n_source_use -1) - 0.50d0
               source_x1_c_x2_c_x3_c_C__(1-1,nsource) = source_r
     $              *cos(source_w)
               source_x1_c_x2_c_x3_c_C__(2-1,nsource) = source_r
     $              *sin(source_w)
               source_x1_c_x2_c_x3_c_C__(3-1,nsource) = source_z
            end if !if (nsource.lt.n_source_cut_C) then
         end if !if (flag_heterogeneity) then
      enddo

      tmp_d_(0) = 3
      tmp_d_(1) = n_source_use
      write(tmp_fname,'(A,A)')
     $     trim(dir_base) , '/source_x1_c_x2_c_x3_c_A___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(2,tmp_d_,source_x1_c_x2_c_x3_c_A__,tmp_fname)
      if (flag_heterogeneity) then
         tmp_d_(0) = 3
         tmp_d_(1) = n_source_use
         write(tmp_fname,'(A,A)')
     $        trim(dir_base) , '/source_x1_c_x2_c_x3_c_B___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,source_x1_c_x2_c_x3_c_B__,tmp_fname)
         tmp_d_(0) = 3
         tmp_d_(1) = n_source_use
         write(tmp_fname,'(A,A)')
     $        trim(dir_base) , '/source_x1_c_x2_c_x3_c_C___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,source_x1_c_x2_c_x3_c_C__,tmp_fname)
      end if !if (flag_heterogeneity) then

c     create physical space grid and sample function_A on it.
      quadrature_type_radial = 0
      quadrature_type_azimu_b = 1
      half_diameter_x_c = 1.0d0
      call createfun_xspace(half_diameter_x_c,n_x_c_max,x_c_step
     $     ,grid_x1_c___,grid_x2_c___,grid_x3_c___
     $     ,source_x1_c_x2_c_x3_c_A__,n_source_use,source_sigma,fgauss
     $     ,f_x_c_A___ ,f_x_c_l2_norm)
      if (verbose.gt.1) write(6,'(A,F8.4)') '  f_x_c_l2_norm is '
     $     ,f_x_c_l2_norm
      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x1_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x1_c___,tmp_fname)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x2_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x2_c___,tmp_fname)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x3_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x3_c___,tmp_fname)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_A___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(3,tmp_d_,f_x_c_A___,tmp_fname)
c     create physical space grid and sample function_B on it.
      if (flag_heterogeneity) then
         quadrature_type_radial = 0
         quadrature_type_azimu_b = 1
         half_diameter_x_c = 1.0d0
         call createfun_xspace(half_diameter_x_c,n_x_c_max,x_c_step
     $        ,grid_x1_c___,grid_x2_c___,grid_x3_c___
     $        ,source_x1_c_x2_c_x3_c_B__,n_source_use,source_sigma
     $        ,fgauss ,f_x_c_B___ ,f_x_c_l2_norm)
         if (verbose.gt.1) write(6,'(A,F8.4)') '  f_x_c_l2_norm is '
     $        ,f_x_c_l2_norm
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_B___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_B___,tmp_fname)
         quadrature_type_radial = 0
         quadrature_type_azimu_b = 1
         half_diameter_x_c = 1.0d0
         call createfun_xspace(half_diameter_x_c,n_x_c_max,x_c_step
     $        ,grid_x1_c___,grid_x2_c___,grid_x3_c___
     $        ,source_x1_c_x2_c_x3_c_C__,n_source_use,source_sigma
     $        ,fgauss ,f_x_c_C___ ,f_x_c_l2_norm)
         if (verbose.gt.1) write(6,'(A,F8.4)') '  f_x_c_l2_norm is '
     $        ,f_x_c_l2_norm
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_C___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_C___,tmp_fname)
      end if !if (flag_heterogeneity) then

c-----------------------------------------------------------------------
c     2) Determine size of spherical grid, allocate space and compute
c     Fourier transform at corresponding points in k-space.
c     Also determine size of spherical harmonic model  (n_Y_lm_sum_max)
c-----------------------------------------------------------------------
      k_p_max = 1.0d0*n_k_p_max
      call getgridr(k_p_max,n_k_p_max,quadrature_type_radial,grid_k_p_
     $     ,weight_k_p_)
      n_Y_lm_sum_max = 0
      n_image_size = 0
      ld_M = 0;
      do nk = 0,n_k_p_max-1
         nlats(nk+1) = nint(pi*grid_k_p_(nk))
         if (nlats(nk+1).lt.6) nlats(nk+1) = 6
         if (mod(nlats(nk),2).ne.0) nlats(nk) = nlats(nk) + 1
         n_polar_a_(nk) = nint(pi*grid_k_p_(nk))
         if (n_polar_a_(nk).lt.6) n_polar_a_(nk) = 6
         if (mod(n_polar_a_(nk),2).ne.0) n_polar_a_(nk) = n_polar_a_(nk)
     $        + 1
         ngridc(nk+1) = nlats(nk+1)*2
         n_w_(nk) = n_polar_a_(nk)*2
         isph_start(nk+1) = n_Y_lm_sum_max
         n_Y_lm_sum_(nk) = n_Y_lm_sum_max
         nterms_sph(nk+1) = nint(grid_k_p_(nk) + 2)
         n_Y_l_(nk) = nint(grid_k_p_(nk) + 2)
         n_Y_lm_sum_max = n_Y_lm_sum_max + (n_Y_l_(nk)+1)**2
         n_image_size = n_image_size + ngridc(nk+1)
         ld_M = ld_M + n_w_(nk)
      enddo
      if (verbose.gt.1) call prinf('  nlats is *',nlats,n_k_p_max)
      if (verbose.gt.0) call prinf('  n_polar_a_ is *',n_polar_a_
     $     ,n_k_p_max)
      if (verbose.gt.1) call prinf('  ngridc is *',ngridc,n_k_p_max)
      if (verbose.gt.0) call prinf('  n_w_ is *',n_w_,n_k_p_max)
      if (verbose.gt.1) call prinf('  n_image_size is *',n_image_size,1)
      if (verbose.gt.0) call prinf('  ld_M is *',ld_M,1)
      ld_S = ld_M
c
      call get_ntot(n_k_p_max,n_polar_a_,quadrature_type_azimu_b
     $     ,n_azimu_b_polar_a_k_p_sum)
      if (verbose.gt.1) write(6,'(A,I0)')
     $     '  n_azimu_b_polar_a_k_p_sum is '
     $     ,n_azimu_b_polar_a_k_p_sum
      if (verbose.gt.1) call prinf(' isph_start *',isph_start,n_k_p_max)
      if (verbose.gt.1) call prinf(' n_Y_lm_sum_ *',n_Y_lm_sum_
     $     ,n_k_p_max)
      if (verbose.gt.1) call prinf(' nterms_sph *',nterms_sph,n_k_p_max)
      if (verbose.gt.1) call prinf(' n_Y_l_ *',n_Y_l_,n_k_p_max)
      if (verbose.gt.1) call prinf(' n_Y_lm_sum_max is *',n_Y_lm_sum_max
     $     ,1)
c
      allocate(grid_k1_c_(0:n_azimu_b_polar_a_k_p_sum-1))
      allocate(grid_k2_c_(0:n_azimu_b_polar_a_k_p_sum-1))
      allocate(grid_k3_c_(0:n_azimu_b_polar_a_k_p_sum-1))
      allocate(weight_k_c_(0:n_azimu_b_polar_a_k_p_sum-1))
      allocate(f_k_c_A_(0:n_azimu_b_polar_a_k_p_sum-1))
      if (flag_heterogeneity) then
         allocate(f_k_c_B_(0:n_azimu_b_polar_a_k_p_sum-1))
         allocate(f_k_c_C_(0:n_azimu_b_polar_a_k_p_sum-1))
      end if !if (flag_heterogeneity) then
      allocate(f_k_c_tru_(0:n_azimu_b_polar_a_k_p_sum-1))
      allocate(f_k_c_est_(0:n_azimu_b_polar_a_k_p_sum-1))
c
      allocate(Y_A_(0:n_Y_lm_sum_max-1))
      if (flag_heterogeneity) then
         allocate(Y_B_(0:n_Y_lm_sum_max-1))
         allocate(Y_C_(0:n_Y_lm_sum_max-1))
      end if !if (flag_heterogeneity) then
      allocate(Y_tru_(0:n_Y_lm_sum_max-1))
      allocate(Y_est_(0:n_Y_lm_sum_max-1))
c
c     create kspace grid data from physical space model A
c
      eps_default = 1.0d-6
      call fun_to_kspacegrid(grid_x1_c___,grid_x2_c___,grid_x3_c___
     $     ,f_x_c_A___,x_c_step,n_x_c_max,eps_default,k_p_max,n_k_p_max
     $     ,quadrature_type_radial,n_polar_a_,quadrature_type_azimu_b
     $     ,n_azimu_b_polar_a_k_p_sum,numonsphere , grid_k1_c_
     $     ,grid_k2_c_ ,grid_k3_c_,weight_k_c_,f_k_c_A_,f_k_c_l2_norm
     $     ,ier)
      call cp1_i4(n_k_p_max,numonsphere,n_azimu_b_polar_a_sum_)
      if (verbose.gt.1) then
         write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
         call write_all_i4(n_k_p_max,numonsphere,13,'numonsphere: ')
         call write_all_i4(n_k_p_max,numonsphere,25
     $        ,' n_azimu_b_polar_a_sum_: ')
      end if
c     create kspace grid data from physical space models B,C
c
      if (flag_heterogeneity) then
         eps_default = 1.0d-6
         call fun_to_kspacegrid(grid_x1_c___,grid_x2_c___,grid_x3_c___
     $        ,f_x_c_B___,x_c_step,n_x_c_max,eps_default,k_p_max
     $        ,n_k_p_max ,quadrature_type_radial,n_polar_a_
     $        ,quadrature_type_azimu_b ,n_azimu_b_polar_a_k_p_sum
     $        ,numonsphere , grid_k1_c_ ,grid_k2_c_ ,grid_k3_c_
     $        ,weight_k_c_,f_k_c_B_,f_k_c_l2_norm ,ier)
         if (verbose.gt.1) then
            write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
            call write_all_i4(n_k_p_max,numonsphere,13,'numonsphere: ')
         end if !if (verbose.gt.1) then
         eps_default = 1.0d-6
         call fun_to_kspacegrid(grid_x1_c___,grid_x2_c___,grid_x3_c___
     $        ,f_x_c_C___,x_c_step,n_x_c_max,eps_default,k_p_max
     $        ,n_k_p_max ,quadrature_type_radial,n_polar_a_
     $        ,quadrature_type_azimu_b ,n_azimu_b_polar_a_k_p_sum
     $        ,numonsphere , grid_k1_c_ ,grid_k2_c_ ,grid_k3_c_
     $        ,weight_k_c_,f_k_c_C_,f_k_c_l2_norm ,ier)
         if (verbose.gt.1) then
            write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
            call write_all_i4(n_k_p_max,numonsphere,13,'numonsphere: ')
         end if !if (verbose.gt.1) then
      end if !if (flag_heterogeneity) then
c
c     convert to model_A: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(f_k_c_A_,k_p_max,n_Y_l_,
     $     quadrature_type_radial,n_k_p_max,n_polar_a_
     $     ,quadrature_type_azimu_b,n_azimu_b_polar_a_sum_,Y_A_
     $     ,n_Y_lm_sum_max)
c     convert to model_B,C: spherical harmonic expansions on successive 
c     spheres
c
      if (flag_heterogeneity) then
         call kspacegrid_to_model(f_k_c_B_,k_p_max,n_Y_l_,
     $        quadrature_type_radial,n_k_p_max,n_polar_a_
     $        ,quadrature_type_azimu_b,n_azimu_b_polar_a_sum_,Y_B_
     $        ,n_Y_lm_sum_max)
         call kspacegrid_to_model(f_k_c_C_,k_p_max,n_Y_l_,
     $        quadrature_type_radial,n_k_p_max,n_polar_a_
     $        ,quadrature_type_azimu_b,n_azimu_b_polar_a_sum_,Y_C_
     $        ,n_Y_lm_sum_max)
      end if !if (flag_heterogeneity) then
c$$$      save both Y_A_ and Y_B_ Y_C_ to disk.
      tmp_d_(0) = n_k_p_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/n_Y_lm_sum__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(tmp_fname)
      call MDA_write_i4(2,tmp_d_,n_Y_lm_sum_,tmp_fname)
      tmp_d_(0) = n_k_p_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/n_Y_l__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(tmp_fname)
      call MDA_write_i4(2,tmp_d_,n_Y_l_,tmp_fname)
      tmp_d_(0) = n_Y_lm_sum_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_A_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(2,tmp_d_,Y_A_,tmp_fname)
      if (flag_heterogeneity) then
         tmp_d_(0) = n_Y_lm_sum_max
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_B_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_B_,tmp_fname)
         tmp_d_(0) = n_Y_lm_sum_max
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_C_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_C_,tmp_fname)
      end if !if (flag_heterogeneity) then

c-----------------------------------------------------------------------
c     3) Assign slices to orientations of spherical grid on highest
c        frequency sphere and define alphas.
c-----------------------------------------------------------------------
      n_image = n_azimu_b_polar_a_sum_(n_k_p_max-1)
      n_M = n_azimu_b_polar_a_sum_(n_k_p_max-1)
      n_S = n_azimu_b_polar_a_sum_(n_k_p_max-1)
      call getspheregrid(n_polar_a_(n_k_p_max-1),quadrature_type_azimu_b
     $     ,grid_cos_polar_a_,grid_sin_polar_a_,weight_cos_polar_a_
     $     ,n_azimu_b_ ,azimu_b_step_ ,n_azimu_b_polar_a_sum_tmp)
      if (verbose.gt.0) write(6,'(A,I0)') '  n_image is ',n_image
      if (verbose.gt.0) write(6,'(A,I0)') '  n_M is ',n_M
      if (verbose.gt.0) then
         call write_all_i4(n_k_p_max,n_w_,7
     $        ,' n_w_: ')
         call write_all_i4(n_k_p_max,n_azimu_b_polar_a_sum_,25
     $        ,' n_azimu_b_polar_a_sum_: ')
      end if! if (verbose.gt.0) then
      if (verbose.gt.0) write(6,'(A,I0)')
     $     '  n_azimu_b_polar_a_k_p_sum is '
     $     ,n_azimu_b_polar_a_k_p_sum
      if (verbose.gt.0) write(6,'(A,I0)')
     $     '  n_azimu_b_polar_a_sum_tmp is ',n_azimu_b_polar_a_sum_tmp
c      
      if (flag_displacement) then
         displacement_max = displacement_max_set
      else
         displacement_max = 0.0d0
      end if
      allocate(alpha2d_tru__(0:n_alpha-1,0:n_M-1))
      allocate(alpha2d_est__(0:n_alpha-1,0:n_M-1))
      nA = 0
      do npolar_a = 0,n_polar_a_(n_k_p_max-1)-1
         cos_polar_a = grid_cos_polar_a_(npolar_a)
         sin_polar_a = grid_sin_polar_a_(npolar_a)
         azimu_b_step = azimu_b_step_(npolar_a)
         do nazimu_b = 0,n_azimu_b_(npolar_a)-1
            azimu_b = nazimu_b*azimu_b_step
            alpha2d_tru__(nalpha_polar_a,nA) = dacos(cos_polar_a)
            alpha2d_tru__(nalpha_azimu_b,nA) = azimu_b
            alpha2d_tru__(nalpha_gamma_z,nA) = pi /32.0d0*(mod(1+nA,32)
     $           -16)
            if (flag_displacement) then
               alpha2d_tru__(nalpha_delta_x,nA) = 0.25d0/k_p_max
     $              *(mod(1+nA,5)-1)
               alpha2d_tru__(nalpha_delta_y,nA) = 0.50d0/k_p_max
     $              *(mod(1+nA,3)-1)
            else
               alpha2d_tru__(nalpha_delta_x,nA) = 0.0d0/k_p_max
               alpha2d_tru__(nalpha_delta_y,nA) = 0.0d0/k_p_max
            end if !flag_displacement
            if (flag_normalization) then
c$$$               Below we draw v from [0.5,1.5]
               alpha2d_tru__(nalpha_l2_norm,nA) = 0.5d0 +
     $              adi_rand_f(rseed)
            else
               alpha2d_tru__(nalpha_l2_norm,nA) = 1.0d0
            end if !flag_normalization
c$$$            to start with somewhat arbitrary alpha2d_est__
            alpha2d_tru__(nalpha_ctf_ind,nA) = 0.0d0 + mod(1+nA,3)
            alpha2d_tru__(nalpha_S_index,nA) = 0.0d0
c$$$            estimated parameters used to initialize least-squares-solve
            alpha2d_est__(nalpha_polar_a,nA) = pi *adi_rand_f(rseed)
            alpha2d_est__(nalpha_azimu_b,nA) = 2 *adi_rand_f(rseed)*pi
            alpha2d_est__(nalpha_gamma_z,nA) = pi /2.0d0*(2*mod(1+nA,2)
     $           -1)
            if (flag_displacement) then
               alpha2d_est__(nalpha_delta_x,nA) = 0.125d0/k_p_max
     $              *(mod(1+nA+2,5)-1)
               alpha2d_est__(nalpha_delta_y,nA) = 0.25d0/k_p_max
     $              *(mod(1+nA+1,3)-1)
            else
               alpha2d_est__(nalpha_delta_x,nA) = 0.0d0
               alpha2d_est__(nalpha_delta_y,nA) = 0.0d0
            end if !flag_displacement
            alpha2d_est__(nalpha_l2_norm,nA) = 1.0d0
            alpha2d_est__(nalpha_ctf_ind,nA) = 0.0d0 + mod(1+nA,3)
            alpha2d_est__(nalpha_S_index,nA) = 0.0d0
            nA = nA+1
         enddo !do nazimu_b = 0,n_azimu_b_(npolar_a)-1
      enddo !do npolar_a = 0,n_polar_a_(n_k_p_max-1)-1
      n_azimu_b_polar_a_sum_tmp = nA
      if (verbose.gt.1) write(6,'(A,I0)')
     $     '  n_azimu_b_polar_a_sum_tmp is ',n_azimu_b_polar_a_sum_tmp
      if (n_azimu_b_polar_a_sum_tmp.ne.n_azimu_b_polar_a_sum_(n_k_p_max
     $     -1))then
         write(6,'(A,I0,A,I0,A,I0,A)')
     $        'Warning, n_azimu_b_polar_a_sum_tmp '
     $        ,n_azimu_b_polar_a_sum_tmp,'.ne.',n_azimu_b_polar_a_sum_
     $        ,'(' ,n_k_p_max-1,')'
      end if
c
      allocate(Y_slice_A__(0:ld_M-1,0:n_M-1))
      if (flag_heterogeneity) then
         allocate(Y_slice_B__(0:ld_M-1,0:n_M-1))
         allocate(Y_slice_C__(0:ld_M-1,0:n_M-1))
      end if !if (flag_heterogeneity) then
      allocate(Y_slice__(0:ld_M-1,0:n_M-1))
c$$$      allocate(S__(0:ld_S-1,0:n_M-1))
c$$$      allocate(S_alpha_polar_a_all_(0:n_S-1))
c$$$      allocate(S_alpha_azimu_b_all_(0:n_S-1))
      call get_template_size(n_polar_a_,n_k_p_max,ld_S,ngridc_tmp,
     $     icstart)
      call cp1_i4(n_k_p_max,icstart,n_w_sum_)
      if (verbose.gt.1) write(6,'(A,I0)') '  ld_S is ',ld_S

c-----------------------------------------------------------------------
c     4) compute Y_slice_A__ directly from physical space function (f_x_c_A___).
c        with orientation vectors defined by alpha2d_tru__.
c-----------------------------------------------------------------------
c
      diameter_x_c = 2.0d0
      call mk_simulated_slices_alpha2d(f_x_c_A___,n_x_c_max,diameter_x_c
     $     ,eps_default,n_k_p_max,quadrature_type_radial,n_w_,k_p_max
     $     ,n_alpha ,alpha2d_tru__,n_M ,ld_M,Y_slice_A__)
c-----------------------------------------------------------------------
c     4) compute Y_slice_B__,C__ directly from physical space function 
c        (f_x_c_B___,C___). with orientation vectors defined by alpha2d_tru__.
c-----------------------------------------------------------------------
c
      if (flag_heterogeneity) then
         diameter_x_c = 2.0d0
         call mk_simulated_slices_alpha2d(f_x_c_B___,n_x_c_max
     $        ,diameter_x_c,eps_default,n_k_p_max,quadrature_type_radial
     $        ,n_w_,k_p_max,n_alpha ,alpha2d_tru__,n_M ,ld_M
     $        ,Y_slice_B__)
         diameter_x_c = 2.0d0
         call mk_simulated_slices_alpha2d(f_x_c_C___,n_x_c_max
     $        ,diameter_x_c,eps_default,n_k_p_max,quadrature_type_radial
     $        ,n_w_,k_p_max,n_alpha ,alpha2d_tru__,n_M ,ld_M
     $        ,Y_slice_C__)
      end if !if (flag_heterogeneity) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Determine which images correspond to which model.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_model_(0:n_M-1))
      do nm=0,n_M-1
         if (flag_heterogeneity) then
            tmp_v = adi_rand_f(rseed)
            if (.false.) then
            else if (tmp_v.le.r_heterogeneity_A) then
               I_model_(nm) = 0
            else if (tmp_v.le.r_heterogeneity_A+r_heterogeneity_B) then
               I_model_(nm) = 1
            else if (tmp_v.le.r_heterogeneity_A+r_heterogeneity_B
     $              +r_heterogeneity_C) then
               I_model_(nm) = 2
            else if (tmp_v.le.r_heterogeneity_A+r_heterogeneity_B
     $              +r_heterogeneity_C+r_heterogeneity_D) then
               I_model_(nm) = 3
            else
               I_model_(nm) = 3
            end if !if (.false.) then
         else
            I_model_(nm) = 0
         end if !if (flag_heterogeneity) then
      enddo !do nm=0,n_M-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use Y_slice_A__,B__,C__ to fill Y_slice__ for later use
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_heterogeneity_A=0
      n_heterogeneity_B=0
      n_heterogeneity_C=0
      n_heterogeneity_D=0
      do nm=0,n_M-1
         if (I_model_(nm).eq.0) then
            n_heterogeneity_A=n_heterogeneity_A+1
            call cp1_c16(ld_M,Y_slice_A__(0,nm),Y_slice__(0,nm))
         end if !if (I_model_(nm).eq.0) then
         if (I_model_(nm).eq.1) then
            n_heterogeneity_B=n_heterogeneity_B+1
            call cp1_c16(ld_M,Y_slice_B__(0,nm),Y_slice__(0,nm))
         end if !if (I_model_(nm).eq.0) then
         if (I_model_(nm).eq.2) then
            n_heterogeneity_C=n_heterogeneity_C+1
            call cp1_c16(ld_M,Y_slice_C__(0,nm),Y_slice__(0,nm))
         end if !if (I_model_(nm).eq.0) then
         if (I_model_(nm).eq.3) then
            n_heterogeneity_D=n_heterogeneity_D+1
            call cl1_c16(ld_M,Y_slice__(0,nm))
         end if !if (I_model_(nm).eq.0) then
      enddo !do nm=0,n_M-1
      if (verbose.gt.0) then
         write(6,'(4(A,I0))')
     %           ' n_heterogeneity_A: ' , n_heterogeneity_A
     $        ,  ' n_heterogeneity_B: ' , n_heterogeneity_B
     $        ,  ' n_heterogeneity_C: ' , n_heterogeneity_C
     $        ,  ' n_heterogeneity_D: ' , n_heterogeneity_D
      end if !if (verbose.gt.0) then
c$$$c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$  if flag_SNR, then add noise to slices. 
c$$$c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_snr.eqv..true. .and. snr.gt.0.0d0) then
         tmp_al1 = 0
         do nm=0,n_M-1
            tmp_al1 = tmp_al1 + al1_c16_f(ld_M,Y_slice__(0,nm))/ld_M
         enddo !do nm=0,n_M-1
         tmp_al1 = tmp_al1/max(1,n_M)
         tmp_al2s = 0
         do nm=0,n_M-1
            tmp_al2 = al2_c16_f(ld_M,Y_slice__(0,nm))
            tmp_al2s = tmp_al2s + tmp_al2*tmp_al2
         enddo !do nm=0,n_M-1
         tmp_al2s = tmp_al2s/max(1,n_M)/max(1,ld_M)
         tmp_sigma = dsqrt(tmp_al2s/snr)
         if (verbose.gt.0) then
            write(6,'(A,F16.8)') 'tmp_al1: ' , tmp_al1
            write(6,'(A,F16.8)') 'tmp_al2s: ' , tmp_al2s
            write(6,'(A,F16.8)') 'tmp_sigma: ' , tmp_sigma
            write(6,'(A,F16.8)') 'snr: ' , tmp_al2s/(tmp_sigma
     $           *tmp_sigma)
         end if !if (verbose.gt.0) then
         do nm=0,n_M-1
            call additive_noise_k_p(rseed,n_k_p_max,grid_k_p_,n_w_,ld_M
     $           ,Y_slice__(0,nm),tmp_sigma)
         enddo !do nm=0,n_M-1
      end if !if (flag_snr.eqv..true.) then

c-----------------------------------------------------------------------
c     5) generate ctf functions to associate with slices.
c-----------------------------------------------------------------------
c

c$$$  Number of different ctfs used is n_ctf.
c$$$  Ultimately, we expect this to be proportional to the number of images,
c$$$  with batches of roughly 30-300 images each having the same ctf.
      n_ctf = 3
      ld_CTF = ld_M
      allocate(CTF_k_p_(0:n_ctf*ld_CTF-1))
c$$$  Here we use an unrealistic ctf-function with a high degree
c$$$  of anisotropy to test our code.
      do nctf=0,n_ctf-1
         if (flag_ctf) then
            ctf_p_0 = 10.0d0 + (5.0d0*nctf)/n_ctf
            ctf_p_1 = (pi*nctf)/n_ctf
            ctf_p_2 = 2
            call get_ctf_star_k_p_(n_k_p_max,grid_k_p_,n_w_,ld_CTF
     $           ,CTF_k_p_(nctf*ld_CTF),ctf_p_0,ctf_p_1,ctf_p_2)
         else
            call get_ctf_ones_k_p_(n_k_p_max,grid_k_p_,n_w_,ld_CTF
     $           ,CTF_k_p_(nctf*ld_CTF))
         end if                 ! if (flag_ctf) then
      enddo
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out the ctf-functions
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
      write(tmp_fname,'(A,A)') trim(dir_base) , '/Fig_ctf_c'
      write(6,'(A)') trim(tmp_fname)
      call Fig_gen_ctf_ver0(n_k_p_max,n_polar_a_,grid_k_p_,n_ctf,ld_CTF
     $     ,CTF_k_p_,tmp_fname)
      end if !if (flag_fig) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Randomly determine which images correspond to which ctf-function.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_ctf_(0:n_M-1))
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo
      do nm=0,n_M-1
         alpha2d_tru__(nalpha_ctf_ind,nm) = 1.0d0*I_ctf_(nm)
         alpha2d_est__(nalpha_ctf_ind,nm) = 1.0d0*I_ctf_(nm)
      enddo
     
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Allocate space for templates, as well as other temporary variables.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(S_k_p__(0:ld_S*n_azimu_b_polar_a_k_p_sum-1))
c$$$  n_M_sample_max = maximum number of images to use at each step.
      allocate(I_M_sample_(0:n_M_sample_max-1))
      allocate(I_model_sample_(0:n_M_sample_max-1))
c$$$      allocate(Y_slice_sample_(0:ld_M*n_M_sample_max-1))
c$$$      allocate(M_sample_(0:ld_M*n_M_sample_max-1))
c$$$      allocate(M_residual_transform__(0:ld_M*n_M_sample_max-1))
c$$$      allocate(M_residual_slice__(0:ld_M*n_M_sample_max-1))
      allocate(M_residual_sample__(0:ld_M*n_M_sample_max-1))
      allocate(M_residual_loading_(0:n_residual_loading*n_M_sample_max
     $     -1))
c$$$      allocate(C_M_(0:n_M_sample_max-1))
      allocate(alpha1d_tru_(0:n_alpha*n_M_sample_max-1))
      allocate(alpha1d_est_(0:n_alpha*n_M_sample_max-1))
      allocate(alpha1d_bkp_(0:n_alpha*n_M_sample_max-1))

c-----------------------------------------------------------------------
c     6) Run the main loop (i.e., try and recover molecule from images).
c-----------------------------------------------------------------------
c     

      if (flag_displacement) then
         n_delta_x = n_delta_x_set
         n_delta_y = n_delta_y_set
         n_pixels_in = n_pixels_in_set
      else
         n_delta_x = 1
         n_delta_y = 1
         n_pixels_in = 0.5d0
      end if !flag_displacement

      if (flag_skip) then
         n_k_low = n_k_p_max
      else
         n_k_low = 3
      end if ! if (flag_skip) then
      if (n_k_low.gt.3) then
         write(6,'(A,I0)') 'Warning! n_k_low starting at: ' , n_k_low
      end if !if (n_k_low.gt.3) then
      n_k_bot = n_k_low
      do n_k_cur=n_k_bot,n_k_p_max
         if (verbose.gt.0) write(6,'(A,I0)') 'n_k_cur: ',n_k_cur
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  select images to use 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         n_M_sample = min(n_M,n_M_sample_max)
         n_M_A_sample = 0
         n_M_B_sample = 0
         n_M_C_sample = 0
         n_M_D_sample = 0
         if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_M_sample
     $        ,' images'
         do nM_sample=0,n_M_sample-1
            if (n_M_sample.lt.n_M) then
               I_M_sample_(nM_sample) = min(n_M-1,floor(n_M
     $              *1.0d0*nM_sample/(n_M_sample-1)))
            else
               I_M_sample_(nM_sample) = nM_sample
            end if
            nm = I_M_sample_(nM_sample)
            I_model_sample_(nM_sample) = I_model_(nm)
c$$$            call cp1_c16(ld_M,Y_slice__(0,nm)
c$$$     $           ,Y_slice_sample_(nM_sample*ld_M))
            if (I_model_(nm).eq.0) then
               n_M_A_sample = n_M_A_sample + 1
            end if !if (I_model_(nm).eq.0) then
            if (I_model_(nm).eq.1) then
               n_M_B_sample = n_M_B_sample + 1
            end if !if (I_model_(nm).eq.1) then
            if (I_model_(nm).eq.2) then
               n_M_C_sample = n_M_C_sample + 1
            end if !if (I_model_(nm).eq.1) then
            if (I_model_(nm).eq.3) then
               n_M_D_sample = n_M_D_sample + 1
            end if !if (I_model_(nm).eq.1) then
         enddo                  ! do nM_sample=0,n_M_sample-1
         if (verbose.gt.1 .or. n_k_cur.eq.n_k_bot) then
            write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'n_M_sample: ' ,
     $           n_M_sample,'; n_M_A_sample: ' , n_M_A_sample ,
     $           '; n_M_B_sample: ', n_M_B_sample , '; n_M_C_sample: ' ,
     $           n_M_C_sample , '; n_M_D_sample: ' , n_M_D_sample
         end if !if (verbose.gt.0) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha2d__ to alpha1d_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'copying alpha2d-->alpha'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha2d_tru__(0,nm),alpha1d_tru_(0
     $           +n_alpha*nM_sample))
            call cp1_r8(n_alpha,alpha2d_est__(0,nm),alpha1d_est_(0
     $           +n_alpha*nM_sample))
            if (n_k_cur.eq.n_k_bot) then
               call cp1_r8(n_alpha,alpha2d_est__(0,nm),alpha1d_bkp_(0
     $              +n_alpha*nM_sample))               
            end if !if (n_k_cur.eq.n_k_bot) then
         enddo
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write alphas to disk for postmortem analysis.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = n_alpha
         tmp_d_(1) = n_M_sample
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/alpha1d_tru_'
     $        ,n_k_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,alpha1d_tru_,tmp_fname)
         tmp_d_(0) = n_alpha
         tmp_d_(1) = n_M_sample
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/alpha1d_est_'
     $        ,n_k_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,alpha1d_est_,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images (in place)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            nctf = I_ctf_(nm)
c$$$            call xx1_c16(ld_M,Y_slice__(0,nm),CTF_k_p_(nctf
c$$$     $           *ld_CTF),M_sample_(nM_sample*ld_M))
            call xx1_c16(ld_M,Y_slice__(0,nm),CTF_k_p_(nctf
     $           *ld_CTF),Y_slice__(0,nm))
         enddo                  ! do nM_sample=0,n_M_sample-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor (in place)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            tmp_l2_norm = alpha1d_tru_(nalpha_l2_norm+n_alpha*nM_sample)
c$$$            call af1_c16(ld_M,1.0d0*cmplx(tmp_l2_norm,0.0d0),1.0d0
c$$$     $           *cmplx(0.0d0,0.0d0),M_sample_(nM_sample*ld_M)
c$$$     $           ,M_sample_(nM_sample*ld_M))
            call af1_c16(ld_M,1.0d0*cmplx(tmp_l2_norm,0.0d0),1.0d0
     $           *cmplx(0.0d0,0.0d0),Y_slice__(0,nm),Y_slice__(0,nm))
         enddo                  ! do nM_sample=0,n_M_sample-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) , '/Fig_MC_c_'
     $           ,n_k_cur
            call Fig_gen_ctf_ver3(n_k_cur,n_polar_a_,grid_k_p_
     $           ,n_M_sample,I_M_sample_,ld_M,Y_slice__,ld_CTF,CTF_k_p_
     $           ,n_alpha,alpha1d_est_,min(16,n_M_sample),tmp_fname)
         end if !if (flag_fig) then
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model based on true alpha.
c$$$  Y_tru_ is the model obtained by using the
c$$$  true alpha-parameters for each image.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         lsq_oversample = 2.0d0
         eps_default = 1.0d-6
         lsq_interpolation_order = 5         
         if (verbose.gt.0) write(6,*) 'rebuilding Y_tru_'
         timing_tic = second()
         call rebuild_model_1(n_M_sample,I_M_sample_,ld_M,Y_slice__
     $        ,n_ctf ,ld_CTF,CTF_k_p_,n_alpha,alpha1d_tru_ ,n_w_sum_
     $        ,n_polar_a_,quadrature_type_azimu_b,grid_k_p_,n_w_,n_k_low
     $        ,n_k_cur,n_Y_lm_sum_,n_Y_l_ ,lsq_oversample
     $        ,lsq_interpolation_order,eps_default ,Y_tru_)
         timing_toc = second()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'rebuild_model_0 (est):'
     $           ,' total_time ',timing_toc-timing_tic
         end if !if (verbose.gt.0) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'Y_tru_ to kspacegrid'
         call model_to_kspacegrid(f_k_c_tru_,k_p_max,n_Y_l_
     $        ,quadrature_type_radial,n_k_p_max,n_polar_a_
     $        ,quadrature_type_azimu_b,n_azimu_b_polar_a_sum_,Y_tru_
     $        ,n_Y_lm_sum_max)
         do nA = 0,n_azimu_b_polar_a_k_p_sum-1
            f_k_c_tru_(nA) = f_k_c_tru_(nA)*weight_k_c_(nA)
         enddo
         fft_iflag = -1
         n_x1_c_x2_c_x3_c_sum = n_x_c_max*n_x_c_max*n_x_c_max
         if (verbose.gt.0) write(6,'(A)') 'finufft3d3_f'
         call  finufft3d3_f(n_azimu_b_polar_a_k_p_sum,grid_k1_c_
     $        ,grid_k2_c_,grid_k3_c_ ,f_k_c_tru_,fft_iflag,eps_default
     $        ,n_x1_c_x2_c_x3_c_sum,grid_x1_c___,grid_x2_c___
     $        ,grid_x3_c___,f_x_c_tru___ ,ier)
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/f_x_c_tru___'
     $        ,n_k_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_tru___,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use Y_tru_ and alpha1d_tru_ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         if (flag_residual) then
            n_Y_lm_sum_cur = n_Y_lm_sum_max
            if (n_k_cur.lt.n_k_p_max) then
               n_Y_lm_sum_cur = n_Y_lm_sum_(n_k_cur+1-1)
            end if !if (n_k_cur.lt.n_k_p_max) then
            call test_residual_3(rseed,n_M_sample,I_M_sample_,ld_M
     $           ,Y_slice__,n_ctf,ld_CTF,CTF_k_p_,n_alpha,alpha1d_tru_
     $           ,n_w_sum_,n_polar_a_,quadrature_type_azimu_b,grid_k_p_
     $           ,n_w_,n_k_1,n_k_cur,n_Y_lm_sum_cur ,n_Y_lm_sum_ ,n_Y_l_
     $           ,lsq_oversample,lsq_interpolation_order,Y_tru_
     $           ,M_residual_sample__,n_residual_loading
     $           ,n_residual_iteration,M_residual_loading_)
            tmp_d_(0) = n_M_sample
            tmp_d_(1) = 1
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/I_model_sample_tru_',n_k_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : '
     $           ,trim(tmp_fname)
            call MDA_write_i4(2,tmp_d_,I_model_sample_,tmp_fname)
            tmp_d_(0) = n_residual_loading
            tmp_d_(1) = n_M_sample
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/M_residual_loading_tru_',n_k_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $           ,trim(tmp_fname)
            call MDA_write_c16(2,tmp_d_,M_residual_loading_,tmp_fname)
         end if !if (flag_residual) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for Y_tru_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_residual .and. flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) ,
     $           '/Fig_residual_tru_',n_k_cur
            call Fig_R_ver1(n_k_cur,n_polar_a_,grid_k_p_,n_M_sample
     $           ,I_M_sample_,ld_M,Y_slice__,M_residual_sample__,n_alpha
     $           ,alpha1d_tru_,min(16 ,n_M_sample) ,tmp_fname)
         end if !if (flag_residual .and. flag_fig) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model based on current alpha.
c$$$  Y_est_ is the model obtained by using the
c$$$  best-guess alpha-parameters for each image (which depends
c$$$  on previous iterations).
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         lsq_oversample = 2.0d0
         eps_default = 1.0d-6
         lsq_interpolation_order = 5         

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  update alpha1d_est_ and Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         write(tmp_dname,'(A)') dir_base
         n_gamma_z = max(32,n_w_(n_k_cur))
         if (flag_tesselation) then
            tesselation_distance_req = 2.0d0 / max(1,n_k_cur)
         else
            tesselation_distance_req = 2.0d0
         end if !flag_tesselation
         call cp1_r8(n_alpha*n_M_sample_max,alpha1d_est_,alpha1d_bkp_)
         call test_alpha_update_4(verbose,rseed,n_omp_sub
     $        ,fpm_howmany_max,half_diameter_x_c_tmp,n_k_p_max
     $        ,n_M_sample,I_M_sample_,ld_M,Y_slice__ ,n_w_sum_
     $        ,n_polar_a_,quadrature_type_azimu_b,grid_k_p_ ,n_w_
     $        ,n_k_low,n_k_cur,n_Y_lm_sum_max,n_Y_lm_sum_,n_Y_l_
     $        ,lsq_oversample,lsq_interpolation_order,eps_default,n_ctf
     $        ,ld_CTF,CTF_k_p_,n_alpha,alpha1d_est_,alpha_update_f
     $        ,flag_MS_vs_SM ,Y_est_ ,n_azimu_b_polar_a_sum_,ld_S
     $        ,n_S_sample_max ,S_k_p__,n_delta_x ,n_delta_y,n_gamma_z
     $        ,svd_calculation_type ,eps_svd ,n_pixels_in
     $        ,displacement_max,flag_tesselation
     $        ,tesselation_distance_req, flag_fig,tmp_dname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'Y_est_ to kspacegrid'
         call model_to_kspacegrid(f_k_c_est_,k_p_max,n_Y_l_
     $        ,quadrature_type_radial,n_k_p_max,n_polar_a_
     $        ,quadrature_type_azimu_b,n_azimu_b_polar_a_sum_,Y_est_
     $        ,n_Y_lm_sum_max)
         do nA = 0,n_azimu_b_polar_a_k_p_sum-1
            f_k_c_est_(nA) = f_k_c_est_(nA)*weight_k_c_(nA)
         enddo
         fft_iflag = -1
         n_x1_c_x2_c_x3_c_sum = n_x_c_max*n_x_c_max*n_x_c_max
         if (verbose.gt.0) write(6,'(A)') 'finufft3d3_f'
         call  finufft3d3_f(n_azimu_b_polar_a_k_p_sum,grid_k1_c_
     $        ,grid_k2_c_,grid_k3_c_ ,f_k_c_est_,fft_iflag,eps_default
     $        ,n_x1_c_x2_c_x3_c_sum,grid_x1_c___,grid_x2_c___
     $        ,grid_x3_c___,f_x_c_est___ ,ier)
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/f_x_c_est___'
     $        ,n_k_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_est___,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use Y_est_ and alpha1d_bkp_ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         if (flag_residual) then
            n_Y_lm_sum_cur = n_Y_lm_sum_max
            if (n_k_cur.lt.n_k_p_max) then
               n_Y_lm_sum_cur = n_Y_lm_sum_(n_k_cur+1-1)
            end if !if (n_k_cur.lt.n_k_p_max) then
            call test_residual_3(rseed,n_M_sample,I_M_sample_,ld_M
     $           ,Y_slice__,n_ctf,ld_CTF,CTF_k_p_,n_alpha,alpha1d_bkp_
     $           ,n_w_sum_,n_polar_a_,quadrature_type_azimu_b,grid_k_p_
     $           ,n_w_,n_k_1,n_k_cur,n_Y_lm_sum_cur ,n_Y_lm_sum_ ,n_Y_l_
     $           ,lsq_oversample,lsq_interpolation_order,Y_est_
     $           ,M_residual_sample__,n_residual_loading
     $           ,n_residual_iteration,M_residual_loading_)
            tmp_d_(0) = n_M_sample
            tmp_d_(1) = 1
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/I_model_sample_est_',n_k_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : '
     $           ,trim(tmp_fname)
            call MDA_write_i4(2,tmp_d_,I_model_sample_,tmp_fname)
            tmp_d_(0) = n_residual_loading
            tmp_d_(1) = n_M_sample
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/M_residual_loading_est_',n_k_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $           ,trim(tmp_fname)
            call MDA_write_c16(2,tmp_d_,M_residual_loading_,tmp_fname)
         end if !if (flag_residual) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_residual .and. flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) ,
     $           '/Fig_residual_est_',n_k_cur
            call Fig_R_ver1(n_k_cur,n_polar_a_,grid_k_p_,n_M_sample
     $           ,I_M_sample_,ld_M,Y_slice__,M_residual_sample__,n_alpha
     $           ,alpha1d_bkp_,min(16 ,n_M_sample) ,tmp_fname)
         end if !if (flag_residual .and. flag_fig) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha to alpha2d
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'copying alpha-->alpha2d'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha1d_est_(0+n_alpha*nM_sample)
     $           ,alpha2d_est__(0,nm))
         enddo
         if (flag_normalization) then
c$$$  do nothing
         else
            do nM_sample=0,n_M_sample-1
               nm = I_M_sample_(nM_sample)
               alpha2d_est__(nalpha_l2_norm,nm) = 1.0d0
            enddo
         end if                 !flag_normalization
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  end loop
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         
      enddo                     ! do n_k_cur = n_k_bot,n_k_p_max

 10   continue
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver14]'
      end if
      stop
      end
