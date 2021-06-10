!> Doxygen comment: ;\n
!>        Driver for molecular reconstruction. ;\n
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
!>        Generally speaking, this driver has two modes. ;\n
!>        Mode_0: The driver first generates 1-to-3 molecules ;\n
!>                (i.e., molecules-A,B,C) which are rather similar. ;\n
!>                Using these three molecules, we produce a set of 'slices', some from ;\n
!>                molecule-A, some from B, some from C. ;\n
!>                We then use a set of (typically 3) anisotropic CTF-functions to ;\n
!>                distort these 'slices', producing 'images'. ;\n
!>                These CTF-functions can be artificial (for debugging), ;\n
!>                or more realistic (from niko_ctf). ;\n
!>                We also distort these images by translating and rotating. ;\n
!>                We also distort these images by adding gaussian noise (iid per mode). ;\n
!>                Using these images, along with (randomly chosen) initial image-parameters, ;\n
!>                we try and reconstruct the molecule out to some maximum radius ;\n
!>                in k-space polar coordinates (denoted by n_k_p_r_max). ;\n
!>                Along the way we try and estimate the cross-correlations ;\n
!>                between residuals, allowing us to classify the images (i.e., either A or B). ;\n
!>        Mode_1: The driver reads 1-to-3 densities (e.g., from dir_trpv1/data_nosym). ;\n
!>                From these density a collection of 'slices' are produced. ;\n
!>                From this point on the driver proceeds as in Mode_0. ;\n
!>        Mode_2: The driver reads a particular density (e.g., from dir_trpv1/data_nosym). ;\n
!>                Additionally, the driver reads images (e.g., from dir_trpv1/data_nosym/images). ;\n
!>                The 'tru' image-parameters for each images are drawn from ;\n
!>                dir_trpv1/data_nosym/ as well, and these image-parameters are ;\n
!>                used to determine the 'Y_tru', etc. in the subsequent ;\n
!>                molecular-reconstruction loop. ;\n
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
!>        Make with: make -f test_ver18.make ; 
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
!>        try out Mode_0: ;\n
!>        ./test_ver18.out 0 12 0.4 0.3 0.2 0.1 100.0d0 1 ; 
!>        ./test_ver18.out 0 12 0.4 0.3 0.2 0.1 10.0d0 1 ; 
!>        ./test_ver18.out 0 12 0.4 0.3 0.2 0.1 1.0d0 1 ; 
!>        ./test_ver18.out 0 12 0.4 0.3 0.2 0.1 0.1d0 1 ; 
!>        ./test_ver18.out 0 12 0.4 0.3 0.2 0.1 0.01d0 1 ; 
!>        command line arguments: A,B,C,D,snr,rseed <-- relative abundance of molecules A,B,C and D, as well as snr (signal-to-noise) and random seed. ;\n
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
!>        try out Mode_1: ;\n
!>        ./test_ver18.out 1 12 1 0 0 0 100.0d0 1 ; 
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
!>        try out Mode_2: ;\n
!>        make -f test_ver18.make ; ./test_ver18.out 2 12 1 0 0 0 0 1 ; 
!>        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;\n
      implicit none
      include 'omp_lib.h'
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Here are several parameters that ;
c$$$      determine what the driver will do. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer verbose ! verbosity level. ;
      data verbose / 1 /
      integer nMode ! Mode of operation for the driver. 0=synthetic molecule + synthetic images. 1=molecule from dir_trpv1/data_nosym/density + synthetic images. 2=molecule from dir_trpv1/data_nosym/density + images from dir_trpv1/data_nosym/images. ;
      integer n_k_p_r_max ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
c$$$      parameter (n_k_p_r_max=16) ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ; read from cmd_K ;
      integer n_x_c_max ! maximum number of x-values (in x-space cartesian coordinates) used to generate initial molecule. sometimes named ngrid. ;
      integer n_x_c_from_file ! maximum number of x-values (in x-space cartesian coordinates) read from file (used in Mode_2). ;
      complex *16, allocatable :: M_x_c_from_file__(:,:) ! image read from file (used in Mode_2). ;
      parameter (n_x_c_max=128) ! maximum number of x-values (in x-space cartesian coordinates) used to generate initial molecule. sometimes named ngrid. ;
      logical flag_displacement ! flag determining whether or not displacements will be considered when matching images to templates (and subsequently used in least-squares calculation). ;
c$$$      parameter (flag_displacement=.false.)
      integer n_delta_x_set ! if displacements are considered, this value determines the number of displacements considered in the x-direction (in the plane). ;
c$$$      parameter (n_delta_x_set=3)
c$$$      integer n_delta_y_set ! if displacements are considered, this value determines the number of displacements considered in the y-direction (in the plane). ;
c$$$      parameter (n_delta_y_set=2) ! This is fixed to equal n_delta_x_set
      integer n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
      real *8 n_pixels_in_set ! if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength under consideration (which can change from iteration to iteration). ;
      parameter (n_pixels_in_set=0.5d0)
      real *8 displacement_max_set ! if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
      parameter (displacement_max_set=0.25d0)
      integer svd_calculation_type_set ! flag determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == -1 --> svd_calculation_type chosen automatically from 1,2 based on number of translations and image dimensions. ;
c$$$      svd_calculation_type == 0 --> multiply, then fft (this should only ever be used for debugging);
c$$$      svd_calculation_type == 1 --> fft, then multiply ;
c$$$      svd_calculation_type == 2 --> brute force displacements ;
      parameter (svd_calculation_type_set=1)
      integer svd_calculation_type ! flag determining how innerproducts are computed across rotations and translations. ;
      real *8 eps_svd ! svd tolerance epsilon, typically 1. ;
      parameter (eps_svd=0.1d0)
      logical flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
      parameter (flag_RTRT_vs_RTTR=.false.)
      logical flag_normalization ! flag determining whether or not l2_norm will be calculated when matching images to templates (and subsequently used in least-squares calculation). ;
      parameter (flag_normalization=.false.)
      integer ctf_type ! flag determining which kind of ctf to use (e.g., an anisotropic ctf). 0=unity, 1=adi_ctf, 2=niko_ctf (from marina). Note that Mode_2 assumes that ctf_type.eq.2. ; 
      parameter (ctf_type=2)
      character fname_n_ctf*(*) !character array: n_ctf for Mode 2. ;
      parameter (fname_n_ctf
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/num_ctf')
      character fname_euler_angle_*(*) !character array: euler_angle_ for Mode 2. ;
      parameter (fname_euler_angle_
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/euler_angles')
      character fname_rotation_matrix_*(*) !character array: rotation_matrix_ for Mode 2. ;
      parameter (fname_rotation_matrix_
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/rot_matrices')
      character fname_I_ctf_*(*) !character array: I_ctf_ for Mode 2. ;
      parameter (fname_I_ctf_
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/ctf_idx')
      character fname_microscope_parameter_*(*) !character array: microscope_parameter_ for Mode 2. ;
      parameter (fname_microscope_parameter_
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/mscope_params'
     $)
      character fname_ctf_parameter_*(*) !character array: ctf_parameter_ for Mode 2. ;
      parameter (fname_ctf_parameter_
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/ctf_params')
      logical flag_tesselation ! flag determining whether or not to use adaptive sampling of templates when calculating innerproducts for each image. ;
      parameter (flag_tesselation=.false.)
      real *8 tesselation_distance_req !temporary: if adaptively sampling templates, this value determines the maximum distance (on the unit-sphere) allowed between the viewing-angle of any particular template and the estimated-viewing-angle of a given image. If this value is large (e.g., 2.0d0), then every template will be considered for each image. If, on the other hand, this value is small (e.g., <0.125d0), then only a few templates will be considered for each image. When adaptively sampling templates, we might consider this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
      integer n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      integer n_LT_ref ! number of image-template pairs to consider when refining local search.
      logical flag_heterogeneity !flag determining whether or not to include multiple (i.e., 3) molecules in image pool. ;
      parameter (flag_heterogeneity=.false.)
      character fname_dim_A*(*) !character array: dimensions of molecule A used for Mode 1,2. ;
      character fname_dim_B*(*) !character array: dimensions of molecule B used for Mode 1,2. ;
      character fname_dim_C*(*) !character array: dimensions of molecule C used for Mode 1,2. ;
      character fname_dim_D*(*) !character array: dimensions of molecule D used for Mode 1,2. ;
      character fname_image_A*(*) !character array: images of molecule A used for Mode 1,2. ;
      character fname_density_A*(*) !character array: density of molecule A used for Mode 1,2. ;
      character fname_density_B*(*) !character array: density of molecule B used for Mode 1,2. ;
      character fname_density_C*(*) !character array: density of molecule C used for Mode 1,2. ;
      character fname_density_D*(*) !character array: density of molecule D used for Mode 1,2. ;
      parameter (fname_dim_A
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/dims')
      parameter (fname_dim_B
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/dims')
      parameter (fname_dim_C
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/dims')
      parameter (fname_dim_D
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/dims')
      parameter (fname_image_A
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/images_bin')
      parameter (fname_density_A
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/density')
      parameter (fname_density_B
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/density')
      parameter (fname_density_C
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/density')
      parameter (fname_density_D
     $     ='/data/rangan/dir_cryoem/dir_trpv1/data_nosym/density')
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
      integer n_M_from_file ! n_M read from file. ;
      integer n_M_sample_max_set ! maximum number of images to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ; A value of 0 implies all images will be used. ;
      parameter (n_M_sample_max_set=4096)
      integer n_M_sample_max ! maximum number of images to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ; A value of 0 implies all images will be used. ;
      integer n_S_sample_max_set ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ; A value of 0 implies all templates will be used. ;
      parameter (n_S_sample_max_set=0)
      integer n_S_sample_max ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ; A value of 0 implies all templates will be used. ;
      logical flag_snr !flag determining whether or not to add iid gaussian noise to images. ;
      parameter (flag_snr=.true.)
      real *8 snr !actual snr used. see command line inputs. ;
c$$$      parameter (snr=10.0d0)
      integer fpm_howmany_max ! Maximum number of fftws to call simultaneously within the fftw_plan_many. 
      parameter (fpm_howmany_max=16)
      integer n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
      parameter (n_omp_sub_0in=8)
      integer n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
c$$$      parameter (n_S_0_sub_0in=2)
      integer n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
c$$$      parameter (n_S_1_sub_0in=3)
      integer n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
c$$$      parameter (n_M_0_sub_0in=3)
      integer n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
c$$$      parameter (n_M_1_sub_0in=4)
      real *8 d_memory_limit !real *8: upper limit on the memory allowed by test_innerproduct_8. ;
      parameter (d_memory_limit=4.0d9)
      real *8 alpha_update_f ! fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.).
      parameter (alpha_update_f=0.00d0)
      logical flag_MS_vs_SM !determines whether to assign images to templates (.true.) or templates to images (.false.). ;
c$$$      parameter (flag_MS_vs_SM=.true.)
      integer n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
      parameter (n_SM_max=16)
      integer n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
c$$$      parameter (n_MS_max=0)
      logical flag_skip_marching !determines whether to start iteration at n_k_p_r_low = n_k_p_r_max (.true.) or 2 (.false.). ;
      parameter (flag_skip_marching=.false.)
      logical flag_perform_estimation !determines whether to perform estimation (.true.) or not (.false.). ;
      parameter (flag_perform_estimation=.true.)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$      Here are other variables are ;
c$$$      internally generated ;
c$$$      or temporary. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer min_i4_f,max_i4_f !function output. ;
      integer n,nx_c,ny_c !temporary index for allocating size. ;
      logical flag_memory_checkset !temporary flag for checking memory setting. ;
      integer npolar_a,nazimu_b,nA !temporary indices for polar_a, azimu_b and azimu_b_polar_a. ;
      integer n_polar_a_lowerbound ! lower bound on the number of k-values (in k-space polar coordinates). sometimes named ntmax. ;
c$$$      parameter (n_polar_a_lowerbound=52)
      parameter (n_polar_a_lowerbound=12)
      integer n_polar_a_upperbound ! upper bound on the number of k-values (in k-space polar coordinates). sometimes named ntmax. ;
      integer fft_iflag !temporary: stores direction of fft (+1 or -1). sometimes named iflag. ;
      integer ier !temporary: stores error flag (integer). ;
      integer quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
      integer quadrature_type_radial ! quadrature type in the radial direction. sometimes named ityper. ;
      integer lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
      external multaha ! least-squares matrix. ;
      integer n_image_size ! total number of entries in each image (in k-space polar coordinates). sometimes called ld_M or nimagesize. ;
      integer ld_M ! total number of entries in each image (in k-space polar coordinates). i.e., leading-dimension of M array. sometimes called n_image_size or nimagesize. ;
      integer n_image ! total number of images. sometimes called n_M or nimages. ;
      integer n_M ! total number of images. sometimes called n_image or nimages. ;
      integer n_k_p_r_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_r_max. ;
      integer n_k_p_r_bot ! lowest value of nk (k-index in k-space polar coordinates) to use when reconstructing molecule. Note that this runs from 1 to n_k_p_r_max. ;
      integer n_x0_c_x1_c_x2_c_sum ! total number of points in x-space cartesian coordinates. sometimes called npts. ;
      integer n_azimu_b__sum_max !temporary: total number of points on sphere at a particular radius grid_k_p_r_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
      integer n_Y_lm_sum ! total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
      integer ld_S ! total number of entries in each template (in k-space polar coordinates). sometimes called ntemplate_size or ntemplatesize. ;
      integer n_S ! total number of templates. sometimes called ntemplates. ;
      integer n_azimu_b__sum_sum ! sum of azimu_b__sum_. total number of points on all spheres (out to radius grid_k_p_r_(n_k_p_r_max-1)). sometimes called ntot. ;
      integer n_azimu_b__sum_cur !temporary: n_azimu_b__sum_(n_k_p_r_cur-1). sometimes called numonsphere_cur. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      real *8 pi ! 4.0d0*datan(1.0d0) ;
      real *8 half_diameter_x_c ! half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
      real *8 x_c_step ! spacing of grid in x-space cartesian coordinates. sometimes called h. ;
      real *8 diameter_x_c ! diameter of particle support in x-space cartesian grid. sometimes called boxsize. ;
      real *8 cos_polar_a !temporary: dcos(polar angle). ;
      real *8 eps_default !temporary: tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
      real *8 f_k_c_l2_norm !temporary: l2-norm of f in k-space cartesian coordinates. sometimes called ffhatnorm. ;
      real *8 f_x_c_l2_norm !temporary: l2-norm of f in x-space cartesian coordinates. sometimes called fnorm. ;
      real *8 azimu_b !temporary: azimuthal angle. sometimes called phi. ;
      real *8 azimu_b_step !temporary: grid spacing for azimuthal angle. sometimes called phistep. ;
      real *8 k_p_r_max ! maximum value of k in k-space polar coordinates. sometimes called rmax. ;
      real *8 sin_polar_a !temporary: dsin(polar angle). ;
      real *8 lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer n_k_p_r_cur ! current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_r_max. ;
      integer nk ! current index for k-value grid_k_p_r_(nk). ;
      integer nw ! current angle within n_w_(nk). ;
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
      integer, allocatable :: n_polar_a_(:) ! number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
      integer n_polar_a_min ! minimum number of polar_a_ values over nk. This should be higher than n_polar_a_lowerbound. ;
      integer n_polar_a_max ! maximum number of polar_a_ values over nk. This should be lower than n_polar_a_upperbound. ;
      integer, allocatable :: n_azimu_b__(:) ! n_azimu_b__(np + nk*n_polar_a_max) = number of azimu_b values on sphere of radius grid_k_p_r_(nk) and polar_a value polar_a_(np). ;
      integer, allocatable :: n_azimu_b__sum_(:) ! n_azimu_b__sum_(nk) = sum of n_azimu_b__ over np. sometimes called numonsphere(nk). ;
      integer, allocatable :: numonsphere(:) ! total number of points on sphere at radius grid_k_p_r_(nk). also called n_azimu_b__sum_(nk). ;
      integer, allocatable :: nterms_sph(:) ! order of spherical harmonic expansion at each radius grid_k_p_r_(nk). also called n_Y_l_(nk). ;
      integer, allocatable :: n_Y_l_(:) ! order of spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). also called nterms_sph(nk). ;
      integer, allocatable :: n_Y_lm_(:) ! number of terms (i.e., (l,m) pairs) of the spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). ;
      integer, allocatable :: isph_start(:) ! cumulative sum of n_Y_lm_. also called n_Y_lm_csum_(nk). ;
      integer, allocatable :: n_Y_lm_csum_(:) ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
      integer, allocatable :: ngridc(:) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_r_(nk) for a template or image in k-space polar coordinates. also called n_w_(nk). ;
      integer, allocatable :: n_w_(:) ! n_w_(nk) is the number of points on the ring of radius grid_k_p_r_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer, allocatable :: n_w_tmp_(:) !temporary: similar to ngridc. ;
      integer, allocatable :: icstart(:) ! cumulative sum of n_w_. also called n_w_csum_(nk). ;
      integer, allocatable :: n_w_csum_(:) ! cumulative sum of n_w_. also called icstart(nk). ;
      integer n_w_sum !integer: sum of n_w_, total number of points in a particular image or template (in k_p_ coordinates). ;
      integer, allocatable :: ngridps(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_r_(nk). also called n_azimu_b_(np). ;
      integer, allocatable :: n_azimu_b_(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_r_(nk). also called ngridps(np). ;
      real *8, allocatable :: grid_k_p_r_(:) ! real *8 array (length at least n_k_p_r_max). Contains values for k on successive shells for k-space polar coordinates. sometimes called xnodesr(nk). Note that, for compatibility with test_innerproduct_8, this driver assumes that k is given in terms of ordinary frequency (not angular frequency). ;
      real *8, allocatable :: weight_k_p_r_(:) ! real *8 array (length at least n_k_p_r_max). Weight associated with grid_k_p_r_(nk) in k-space polar coordinates. sometimes called wtsr(nk). ;
      real *8, allocatable :: weight_k_p_(:) ! real *8 array (length at least n_w_sum). Contains weight for each point in unrolled polar grid. ;
      real *8, allocatable :: xnodesth(:) !temporary: values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called grid_cos_polar_a_(nk). ;
      real *8, allocatable :: grid_cos_polar_a_(:) !temporary: values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called xnodesth(nk). ;
      real *8, allocatable :: sthetas(:) !temporary: values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called grid_sin_polar_a_(nk). ;
      real *8, allocatable :: grid_sin_polar_a_(:) !temporary: values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called sthetas(nk). ;
      real *8, allocatable :: weight_cos_polar_a_(:) !temporary: weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called wtsth_(np). ;
      real *8, allocatable :: azimu_b_step_(:) !temporary: grid-spacing for azimu_b. ;
      real *8, allocatable :: source_x0_c_x1_c_x2_c_A__(:,:) ! array containing (x0,x1,x2) coordinates of gaussian sources for molecule-A in x-space cartesian coordinates. ;
      real *8, allocatable :: source_x0_c_x1_c_x2_c_B__(:,:) ! array containing (x0,x1,x2) coordinates of gaussian sources for molecule-B in x-space cartesian coordinates. ;
      real *8, allocatable :: source_x0_c_x1_c_x2_c_C__(:,:) ! array containing (x0,x1,x2) coordinates of gaussian sources for molecule-C in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x0_c___(:,:,:) ! array containing x-values for points in 3d regular grid in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x1_c___(:,:,:) ! array containing y-values for points in 3d regular grid in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_x2_c___(:,:,:) ! array containing z-values for points in 3d regular grid in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_A___(:,:,:) ! function corresponding to molecule-A in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_B___(:,:,:) ! function corresponding to molecule-B in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_C___(:,:,:) ! function corresponding to molecule-C in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_tru___(:,:,:) !temporary: function obtained using true-angles to reconstruct molecule in x-space cartesian coordinates. ;
      complex *16, allocatable :: f_x_c_est___(:,:,:) !temporary: function obtained using estimated-angles to reconstruct molecule in x-space cartesian coordinates. ;
      real *8, allocatable :: grid_k0_c_(:) ! values for k0 of (k0,k1,k2) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      real *8, allocatable :: grid_k1_c_(:) ! values for k1 of (k0,k1,k2) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
      real *8, allocatable :: grid_k2_c_(:) ! values for k2 of (k0,k1,k2) (i.e., k-space cartesian coordinates) for points distributed on nested collection of spherical shells in k-space cartesian coordinates. ;
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
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer tmp_d_(0:5)       !temporary: dimensions used for MDA_write. ;
      character(len=1024) str_Mode,str_k,str_RTXX,str_displacement
     $     ,str_svd,str_normalization,str_ctf,str_tesselation
     $     ,str_heterogeneity,str_residual,str_sample,str_snr
     $     ,str_MS_vs_SM,str_rseed,str_mem
      character(len=1024) dir_base,str_system !fixed name based on driver parameters. ;
      character(len=1024) tmp_dname,tmp_fname,format_string !temporary: strings. ;
      character(len=1024) tmp_dir,tmp_str !temporary: strings. ;
      integer tmp_tab !temporary: index. ;
      integer n_gamma_z !temporary: value used for n_gamma_z when calculating innerproducts. ;
      integer n_delta_x,n_delta_y !temporary: values used for n_delta_x and n_delta_y, respectively, when calculating innerproducts involving displacements. ;
      integer *4 n_SM_use !temporary: total number of image-template pairs considered when calculating innerproducts. Might be smaller than the total number when calculating innerproducts adaptively. ;
      parameter(n_SM_use=100)
      real *8 n_pixels_in !temporary: value used for n_pixels_in when calculating innerproducts involving displacements. ;
      real *8 displacement_max !temporary: value used for displacement_max when calculating innerproducts involving displacements. ;
      integer n_M_sample,nM_sample,nm !temporary: values used to track number of sampled images. ;
      integer n_M_A_sample,n_M_B_sample,n_M_C_sample,n_M_D_sample !temporary: values used to track number of sampled images from molecule-A,B,C. ;
      integer n_Y_lm_csum_cur !temporary: total number of basis function used in spherical harmonic expansion (in k-space polar coordinates) out to a particular radius k. sometimes called n_Y_lm_csum_cur
      integer n_S_sample,nS_sample,ns !temporary: values used to track number of sampled templates. ;
      integer *4, allocatable :: I_M_sample_(:) !temporary: indexing variable used to reference images. ;
      integer *4, allocatable :: I_S_sample_(:) !temporary: indexing variable used to reference templates. ;
      complex *16, allocatable :: M_residual_sample__(:) !temporary: residuals = images - slices from current model. ;
      complex *16, allocatable :: M_residual_loading_(:) !temporary: loading vectors associated with residuals. ;
      complex *16, allocatable :: S_k_p__(:) !temporary: stack of templates associated with reconstructed molecule. sometimes called templates. ;
      include 'excerpt_define_nalpha.f'
      real *8 tmp_l2_norm !temporary: estimated l2_norm of an image. ;
      real *8, allocatable :: alpha2d_tru__(:,:) ! true image parameters associated with original image set (e.g., stored in Y_slice_) stored in 2-dimensional array. ;
      real *8, allocatable :: alpha_tru__(:) ! true image parameters associated with sampled image set (e.g., could be stored in Y_slice_sample_ or M_sample_) stored in 1-dimensional array. ;
      real *8, allocatable :: alpha2d_est__(:,:) !temporary: estimated image parameters associated with original image set (e.g., stored in Y_slice_) stored in 2-dimensional array. ;
      real *8, allocatable :: alpha_est__(:) ! estimated image parameters associated with sampled image set (e.g., could be stored in Y_slice_sample_ or M_sample_) stored in 1-dimensional array. ;
      real *8, allocatable :: alpha_bkp__(:) ! estimated image parameters associated with sampled image set from previous iteration stored in 1-dimensional array. ;
      external fgauss ! function to describe a gaussian in x-space cartesian coordinates. ;
      integer *4 n_CTF ! stores the number of CTF functions (i.e., 6)
      integer *4 n_CTF_from_file ! number of CTF functions read from file (used in Mode_2). ;
      integer ld_CTF ! total number of entries in each CTF (in k-space polar coordinates). i.e., leading-dimension of CTF array. Typically expected to be at least ld_M. ;
      integer *4 nctf !temporary: indexes ctf. ;
      real *8 I_ctf_from_file !temporary: read from file. Note that this is real *8, and not an integer (even though each real *8 in the file has 0 fractional part). ;
      real *8 tmp_shift_(0:2-1) !temporary: euler-angles from disk. ;
      real *8 tmp_rotation_matrix_(0:9-1) !temporary: rotation matrix used for converting euler-angles from relion convention to our convention. ;
      real *8 tmp_euler_angle_(0:3-1) !temporary: array used for converting euler-angles from relion convention to our convention. ;
      integer *4, allocatable :: I_ctf_(:) ! indexes which ctf is associated with which image. ;
      real *8 CTF_Spherical_Aberration !spherical aberation of the lens ;
      real *8, allocatable :: CTF_Spherical_Aberration_(:)
      real *8 CTF_Voltage_kV,CTF_Voltage_1V !voltage ;
      real *8,allocatable :: CTF_Voltage_kV_(:)
      real *8 CTF_Amplitude_Contrast !amplitude contrast ;
      real *8, allocatable :: CTF_Amplitude_Contrast_(:)
      real *8 CTF_Magnification !magnification of microscope ;
      real *8, allocatable :: CTF_Magnification_(:)
      real *8 CTF_Detector_Pixel_Size !pixel size of scanner in microns ;
      real *8 CTF_B_factor
      real *8, allocatable :: CTF_Detector_Pixel_Size_(:)
      real *8 CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle !defocus (in Angstroms) and angle of astigmatism ;
      real *8, allocatable :: CTF_Defocus_U_(:)
      real *8, allocatable :: CTF_Defocus_V_(:)
      real *8, allocatable :: CTF_Defocus_Angle_(:)
      real *8 tmp_w1,tmp_w2,tmp_k_c_1,tmp_k_c_2,tmp_ctf_value,tmp_theta
      real *8 CTF_lambda,CTF_Object_Pixel_Size,CTF_lambda_per_box
      real *8 CTF_Defocus_max, CTF_Defocus_min
      real *8 CTF_Astigmatism_max, CTF_Astigmatism_min
      real *8 CTF_Astigmatism_value
      complex *16, allocatable :: CTF_k_p__(:) ! stores the n_CTF different CTF functions in k-space polar coordinates. ;
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4 ! parameters associated with adi CTF generation. ; Note that these CTF functions are strongly anisotropic (cartoonishly so). ;
      integer *4 rseed !random seed. see command line inputs. ;
      real *8 timing_tic,timing_toc,timing_tot !temporary: timing variables. ;
      real *8 timing_all_tic,timing_all_toc,timing_all_tot !temporary: timing variables. ;
      real *8 timing_function_create ! stores timing. ;
      real *8 timing_function_sample ! stores timing. ;
      real *8 timing_Y_create ! stores timing. ;
      real *8 timing_alpha_define ! stores timing. ;
      real *8 timing_Y_slice_define ! stores timing. ;
      real *8 timing_CTF_apply ! stores timing. ;
      real *8 timing_frequency_marching ! stores timing. ;
      real *8, allocatable :: timing_rebuild_tru_(:) ! stores timing. ;
      real *8, allocatable :: timing_residual_tru_(:) ! stores timing. ;
      real *8, allocatable :: timing_rebuild_est_(:) ! stores timing. ;
      real *8, allocatable :: timing_residual_est_(:) ! stores timing. ;
      real *8, allocatable :: timing_template_create_(:) ! stores timing. ;
      real *8, allocatable :: timing_innerproduct_(:) ! stores timing. ;
      real *8, allocatable :: timing_alpha_update_(:) ! stores timing. ;
      real *8 tmp_v,adi_rand_f !temporary: random real number between 0 and 1. ;
      real *8 tmp_al1,al1_c16_f !temporary: l1-norm of slice. ;
      real *8 tmp_al2,tmp_al2s,al2_c16_f !temporary: l2-norm of slice. ;
      real *8 tmp_sigma !temporary: sigma used to add noise. ;
c$$$      command line input variables
      character(len=64) :: cmd_argstring
      integer ncmd
      real *8 cmd_Mode, cmd_K,cmd_A,cmd_B,cmd_C,cmd_D,cmd_snr
      integer *4 cmd_rseed
      integer *4 cmd_displacement
      integer *4 cmd_MS_vs_SM
      real *8 dist_eq
      integer *4 n_all_0in
      real *8, allocatable :: polar_a_all_(:)
      real *8, allocatable :: azimu_b_all_(:)
      real *8, allocatable :: weight_all_(:)
      pi=4.0d0*datan(1.0d0)

c$$$  %%%%%%%%;
c$$$  % test legendre discretization of sphere. ;
c$$$  %%%%%%%%;
c$$$      n_all_0in=0
c$$$      dist_eq = 0.125d0
c$$$      call sample_sphere_legendre(1,3.0d0,dist_eq,n_all_0in,polar_a_all_
c$$$     $     ,azimu_b_all_,weight_all_)
c$$$      allocate(polar_a_all_(0:1+n_all_0in-1))
c$$$      allocate(azimu_b_all_(0:1+n_all_0in-1))
c$$$      allocate(weight_all_(0:1+n_all_0in-1))
c$$$      dist_eq = 0.125d0
c$$$      call sample_sphere_legendre(1,3.0d0,dist_eq,n_all_0in,polar_a_all_
c$$$     $     ,azimu_b_all_,weight_all_)
c$$$      dist_eq = 0.250d0
c$$$      call sample_sphere_legendre(1,3.0d0,dist_eq,n_all_0in,polar_a_all_
c$$$     $     ,azimu_b_all_,weight_all_)
c$$$      dist_eq = 0.500d0
c$$$      call sample_sphere_legendre(1,3.0d0,dist_eq,n_all_0in,polar_a_all_
c$$$     $     ,azimu_b_all_,weight_all_)
c$$$      dist_eq = 0.750d0
c$$$      call sample_sphere_legendre(1,3.0d0,dist_eq,n_all_0in,polar_a_all_
c$$$     $     ,azimu_b_all_,weight_all_)
c$$$      stop

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_ver18]'
      end if

      ncmd = 0
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_Mode
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_K
      ncmd = ncmd + 1
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
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_displacement
      ncmd = ncmd + 1
      call get_command_argument(1+ncmd,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_MS_vs_SM
      ncmd = ncmd + 1
      
      if (verbose.gt.0) then
         write(6,'(A,F8.4)') ' cmd_Mode: ' , cmd_Mode
      end if !if (verbose.gt.0) then

      nMode = nint(cmd_Mode)

      if (verbose.gt.0) then
      if (nMode.eq.2) then
      write(6,'(A)') ' reading from files: '
      write(6,'(A)') trim(fname_dim_A)
      write(6,'(A)') trim(fname_density_A)
      write(6,'(A)') trim(fname_n_ctf)
      write(6,'(A)') trim(fname_euler_angle_)
      write(6,'(A)') trim(fname_rotation_matrix_)
      write(6,'(A)') trim(fname_I_ctf_)
      write(6,'(A)') trim(fname_microscope_parameter_)
      write(6,'(A)') trim(fname_ctf_parameter_)
      cmd_A = 1.0d0
      cmd_B = 0.0d0
      cmd_C = 0.0d0
      cmd_D = 0.0d0
      cmd_snr = 0.0d0
      end if !if (nMode.eq.2) then
      end if !if (verbose.gt.0) then

      if (verbose.gt.0) then
         write(6,'(A,F8.4)') ' cmd_K: ' , cmd_K
         write(6,'(A,F8.4)') ' cmd_A: ' , cmd_A
         write(6,'(A,F8.4)') ' cmd_B: ' , cmd_B
         write(6,'(A,F8.4)') ' cmd_C: ' , cmd_C
         write(6,'(A,F8.4)') ' cmd_D: ' , cmd_D
         write(6,'(A,F8.4)') ' cmd_snr: ' , cmd_snr
         write(6,'(A,I0)') ' cmd_rseed: ' , cmd_rseed
         write(6,'(A,I0)') ' cmd_displacement: ' , cmd_displacement
         write(6,'(A,I0)') ' cmd_MS_vs_SM: ' , cmd_MS_vs_SM
      end if !if (verbose.gt.0) then
      if (cmd_displacement.eq.0) then
         flag_displacement = .false.
         n_delta_x_set = 1
      end if !if (cmd_displacement.eq.0) then
      if (cmd_displacement.gt.0) then
         flag_displacement = .true.
         n_delta_x_set = cmd_displacement         
      end if !if (cmd_displacement.gt.0) then
      if (cmd_MS_vs_SM.eq.0) then
         flag_MS_vs_SM = .false.
      end if !if (cmd_MS_vs_SM.eq.0) then
      if (cmd_MS_vs_SM.eq.1) then
         flag_MS_vs_SM = .true.
      end if !if (cmd_MS_vs_SM.eq.1) then
      
      n_k_p_r_max = nint(cmd_K)
      k_p_r_max = 1.0d0*n_k_p_r_max/(2.0d0*pi)
      n_k_p_r_max = nint(cmd_K)
      k_p_r_max = 2.0d0*n_k_p_r_max/(2.0d0*pi)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,F8.4)')
     $        ' n_k_p_r_max: ' , n_k_p_r_max
     $        ,' k_p_r_max: ' , k_p_r_max
      end if !if (verbose.gt.1) then
      n_polar_a_upperbound = max(4*n_k_p_r_max,2*nint(pi*k_p_r_max))
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
      write(str_Mode,'(A,I0)') '_Mode' , nMode
      write(str_k,'(A,I0)') '_k' , n_k_p_r_max
c$$$      %%%%%%%%%%%%%%%%
      if (flag_displacement) then
         if (flag_RTRT_vs_RTTR.eqv..true.) then
            write(str_RTXX,'(A)') 'RTRT'            
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
         if (flag_RTRT_vs_RTTR.eqv..false.) then
            write(str_RTXX,'(A)') 'RTTR'            
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
         write(str_displacement,'(A,A,3(A,I0))') '_d', trim(str_RTXX) ,
     $        'x' , n_delta_x_set ,'p' , nint(100*n_pixels_in_set) ,'m'
     $        ,nint(100*displacement_max_set) 
         svd_calculation_type = svd_calculation_type_set
      else
         write(str_displacement,'(A)') '_dF'
         svd_calculation_type = 2
      end if !if (flag_displacement) then
c$$$      %%%%%%%%%%%%%%%%
      if (svd_calculation_type.eq.2) then
         write(str_svd,'(A)') '_s2'
      else 
         write(str_svd,'(A,I0)') '_s1e' , nint(1000*eps_svd)
      end if !if (svd_calculation_type.eq.2) then
c$$$      %%%%%%%%%%%%%%%%
      if (flag_normalization) then
         write(str_normalization,'(A)') '_nT'
      else
         write(str_normalization,'(A)') '_nF'
      end if !if (flag_normalization) then
c$$$      %%%%%%%%%%%%%%%%
      if (ctf_type.eq.1) then
         write(str_ctf,'(A)') '_cA'
      else if (ctf_type.eq.2) then
         write(str_ctf,'(A)') '_cN'
      else
         write(str_ctf,'(A)') '_cF'
      end if !if (ctf_type.eq.1) then
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
      write(str_sample,'(A,I0,A,I0)') '_M' , n_M_sample_max_set , '_S' ,
     $     n_S_sample_max_set
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
      write(str_mem,'(A,I0,A)') , '_m' , nint(d_memory_limit*1e-9) , 'G'
c$$$      %%%%%%%%%%%%%%%%
      write(dir_base,'(15(A))') './dir' ,trim(str_Mode),trim(str_k)
     $     ,trim(str_displacement),trim(str_svd),trim(str_normalization)
     $     ,trim(str_ctf),trim(str_tesselation),trim(str_heterogeneity)
     $     ,trim(str_residual),trim(str_sample),trim(str_snr)
     $     ,trim(str_MS_vs_SM),trim(str_rseed),trim(str_mem)

      write(6,'(A,A)') ' dir_base: ' , trim(dir_base)
      write(str_system,'(A,A)') 'mkdir ' , trim(dir_base)
      write(6,'(A)') trim(str_system)
      call system(str_system)

c-----------------------------------------------------------------------
c     0) allocate some arrays (dynamic allocation required for openmp)
c-----------------------------------------------------------------------
      allocate(nlats(1+n_k_p_r_max))
      call cs1_i4(n_k_p_r_max,nlats)
      allocate(n_polar_a_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_polar_a_)
      allocate(n_azimu_b__(0:1+n_polar_a_upperbound*n_k_p_r_max-1))
      call cs1_i4(n_polar_a_upperbound*n_k_p_r_max,n_azimu_b__)
      allocate(n_azimu_b__sum_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_azimu_b__sum_)
      allocate(numonsphere(1+n_k_p_r_max))
      call cs1_i4(n_k_p_r_max,numonsphere)
      allocate(nterms_sph(1+n_k_p_r_max))
      call cs1_i4(n_k_p_r_max,nterms_sph)
      allocate(n_Y_l_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_Y_l_)
      allocate(n_Y_lm_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_Y_lm_)
      allocate(isph_start(1+n_k_p_r_max))
      call cs1_i4(n_k_p_r_max,isph_start)
      allocate(n_Y_lm_csum_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_Y_lm_csum_)
      allocate(ngridc(1+n_polar_a_upperbound))
      call cs1_i4(n_polar_a_upperbound,ngridc)
      allocate(n_w_(0:1+n_polar_a_upperbound-1))
      call cs1_i4(n_polar_a_upperbound,n_w_)
      allocate(n_w_tmp_(1+n_polar_a_upperbound))
      call cs1_i4(n_polar_a_upperbound,n_w_tmp_)
      allocate(icstart(1+n_k_p_r_max))
      call cs1_i4(n_k_p_r_max,icstart)
      allocate(n_w_csum_(0:1+n_k_p_r_max-1))
      call cs1_i4(n_k_p_r_max,n_w_csum_)
      allocate(ngridps(1+n_polar_a_upperbound))
      call cs1_i4(n_polar_a_upperbound,ngridps)
      allocate(n_azimu_b_(0:1+n_polar_a_upperbound-1))
      call cs1_i4(n_polar_a_upperbound,n_azimu_b_)
      allocate(grid_k_p_r_(0:1+n_k_p_r_max-1))
      call cs1_r8(n_k_p_r_max,grid_k_p_r_)
      allocate(weight_k_p_r_(0:1+n_k_p_r_max-1))
      call cs1_r8(n_k_p_r_max,weight_k_p_r_)
      allocate(xnodesth(1+n_polar_a_upperbound))
      call cs1_r8(n_polar_a_upperbound,xnodesth)
      allocate(grid_cos_polar_a_(0:1+n_polar_a_upperbound-1))
      call cs1_r8(n_polar_a_upperbound,grid_cos_polar_a_)
      allocate(sthetas(1+n_polar_a_upperbound))
      call cs1_r8(n_polar_a_upperbound,sthetas)
      allocate(grid_sin_polar_a_(0:1+n_polar_a_upperbound-1))
      call cs1_r8(n_polar_a_upperbound,grid_sin_polar_a_)
      allocate(weight_cos_polar_a_(0:1+n_polar_a_upperbound-1))
      call cs1_r8(n_polar_a_upperbound,weight_cos_polar_a_)
      allocate(azimu_b_step_(0:1+n_polar_a_upperbound-1))
      call cs1_r8(n_polar_a_upperbound,azimu_b_step_)
      allocate(source_x0_c_x1_c_x2_c_A__(0:3-1,0:1+n_source_max-1))
      n = 3*n_source_max
      call cs1_r8(n,source_x0_c_x1_c_x2_c_A__)
      if (flag_heterogeneity) then
         allocate(source_x0_c_x1_c_x2_c_B__(0:3-1,0:1+n_source_max-1))
         n = 3*n_source_max
         call cs1_r8(n,source_x0_c_x1_c_x2_c_B__)
         allocate(source_x0_c_x1_c_x2_c_C__(0:3-1,0:1+n_source_max-1))
         n = 3*n_source_max
         call cs1_r8(n,source_x0_c_x1_c_x2_c_C__)
      end if !if (flag_heterogeneity) then
      allocate(grid_x0_c___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_r8(n,grid_x0_c___)
      allocate(grid_x1_c___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_r8(n,grid_x1_c___)
      allocate(grid_x2_c___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_r8(n,grid_x2_c___)
      allocate(f_x_c_A___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_c16(n,f_x_c_A___)
      if (flag_heterogeneity) then
         allocate(f_x_c_B___(n_x_c_max,n_x_c_max,1+n_x_c_max))
         n = n_x_c_max*n_x_c_max*n_x_c_max
         call cs1_c16(n,f_x_c_B___)
         allocate(f_x_c_C___(n_x_c_max,n_x_c_max,1+n_x_c_max))
         n = n_x_c_max*n_x_c_max*n_x_c_max
         call cs1_c16(n,f_x_c_C___)
      end if !if (flag_heterogeneity) then
      allocate(f_x_c_tru___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_c16(n,f_x_c_tru___)
      allocate(f_x_c_est___(n_x_c_max,n_x_c_max,1+n_x_c_max))
      n = n_x_c_max*n_x_c_max*n_x_c_max
      call cs1_c16(n,f_x_c_est___)

      call prini(6,13)

      half_diameter_x_c = 1.0d0
      diameter_x_c = 2.0d0*half_diameter_x_c
      call mkphysgrid(
     $     half_diameter_x_c
     $     ,n_x_c_max
     $     ,x_c_step
     $     ,grid_x0_c___
     $     ,grid_x1_c___
     $     ,grid_x2_c___
     $     )
      if (verbose.gt.1) then
         call print_sub_r8(n_x_c_max**3,grid_x0_c___,' grid_x0_c___: ')
         call print_sub_r8(n_x_c_max**3,grid_x1_c___,' grid_x1_c___: ')
         call print_sub_r8(n_x_c_max**3,grid_x2_c___,' grid_x2_c___: ')
      end if !if (verbose.gt.1) then

c-----------------------------------------------------------------------
c     1) Mode_1: Load (and subsample) functions from disk
c-----------------------------------------------------------------------
c$$$  %%%%%%%%
      if (nMode.eq.1) then
c$$$  %%%%%%%%
      call get_f_x_c_from_file_0(
     $        n_x_c_max
     $        ,fname_dim_A
     $        ,fname_density_A
     $        ,f_x_c_A___
     $        )
      f_x_c_l2_norm = al2_c16_f(n_x_c_max**3,f_x_c_A___)
      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_A___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(3,tmp_d_,f_x_c_A___,tmp_fname)
c$$$  %%%%%%%%
      if (flag_heterogeneity) then
c$$$  %%%%%%%%
      call get_f_x_c_from_file_0(
     $        n_x_c_max
     $        ,fname_dim_B
     $        ,fname_density_B
     $        ,f_x_c_B___
     $        )
      f_x_c_l2_norm = al2_c16_f(n_x_c_max**3,f_x_c_B___)
      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_B___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(3,tmp_d_,f_x_c_B___,tmp_fname)
c$$$  %%%%%%%%
      call get_f_x_c_from_file_0(
     $     n_x_c_max
     $     ,fname_dim_C
     $     ,fname_density_C
     $     ,f_x_c_C___
     $     )
      f_x_c_l2_norm = al2_c16_f(n_x_c_max**3,f_x_c_C___)
      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_C___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(3,tmp_d_,f_x_c_C___,tmp_fname)
      f_x_c_l2_norm = al2_c16_f(n_x_c_max**3,f_x_c_C___)
c$$$  %%%%%%%%
      end if !if (flag_heterogeneity) then
c$$$  %%%%%%%%
      end if !if (nMode.eq.1) then

c-----------------------------------------------------------------------
c     1) Mode_0: Create synthetic function in physical space
c-----------------------------------------------------------------------
      if (nMode.eq.0) then
c$$$  %%%%%%%%
      timing_tic = omp_get_wtime()
      n_source_use = min(n_source_max,32)
      n_source_cut_B = floor(0.825d0*n_source_use)
      n_source_cut_C = floor(0.750d0*n_source_use)
      source_sigma = 1.0d0/32.0d0 
      do nsource=0,n_source_use-1
         source_r = 0.5d0*(nsource)/(n_source_use-1)
         source_w = 2.0d0*pi*2*nsource/(n_source_use)
         source_z = 1.0d0*(nsource)/(n_source_use-1) - 0.50d0
         source_x0_c_x1_c_x2_c_A__(1-1,nsource) = source_r
     $        *dcos(source_w)
         source_x0_c_x1_c_x2_c_A__(2-1,nsource) = source_r
     $        *dsin(source_w)
         source_x0_c_x1_c_x2_c_A__(3-1,nsource) = source_z
         if (flag_heterogeneity) then
            if (nsource.lt.n_source_cut_B) then
               source_x0_c_x1_c_x2_c_B__(1-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(1-1,nsource)
               source_x0_c_x1_c_x2_c_B__(2-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(2-1,nsource)
               source_x0_c_x1_c_x2_c_B__(3-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(3-1,nsource)
            else
               source_r = 0.5d0*(n_source_cut_B + 3.5d0*(nsource
     $              -n_source_cut_B))/(n_source_use-1)
               source_w = 2.0d0*pi*2*(n_source_cut_B - 0.15d0*(nsource
     $              -n_source_cut_B)*(nsource-n_source_use+1))
     $              /(n_source_use)
               source_z = 1.0d0*(n_source_cut_B - 0.50d0*(nsource
     $              -n_source_cut_B)*(nsource-n_source_use+1))
     $              /(n_source_use -1) - 0.50d0
               source_x0_c_x1_c_x2_c_B__(1-1,nsource) = source_r
     $              *dcos(source_w)
               source_x0_c_x1_c_x2_c_B__(2-1,nsource) = source_r
     $              *dsin(source_w)
               source_x0_c_x1_c_x2_c_B__(3-1,nsource) = source_z
            end if !if (nsource.lt.n_source_cut_B) then
            if (nsource.lt.n_source_cut_C) then
               source_x0_c_x1_c_x2_c_C__(1-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(1-1,nsource)
               source_x0_c_x1_c_x2_c_C__(2-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(2-1,nsource)
               source_x0_c_x1_c_x2_c_C__(3-1,nsource) =
     $              source_x0_c_x1_c_x2_c_A__(3-1,nsource)
            else
               source_r = 0.5d0*(n_source_cut_C + 1.5d0*(nsource
     $              -n_source_cut_C))/(n_source_use-1)
               source_w = 2.0d0*pi*2*(n_source_cut_C - 0.15d0*(nsource
     $              -n_source_cut_C)*(nsource-n_source_use+1))
     $              /(n_source_use)
               source_z = 1.0d0*(n_source_cut_C - 0.25d0*(nsource
     $              -n_source_cut_C)*(nsource-n_source_use+1))
     $              /(n_source_use -1) - 0.50d0
               source_x0_c_x1_c_x2_c_C__(1-1,nsource) = source_r
     $              *dcos(source_w)
               source_x0_c_x1_c_x2_c_C__(2-1,nsource) = source_r
     $              *dsin(source_w)
               source_x0_c_x1_c_x2_c_C__(3-1,nsource) = source_z
            end if !if (nsource.lt.n_source_cut_C) then
         end if !if (flag_heterogeneity) then
      enddo
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- function_create '
      timing_function_create = timing_tot
c$$$  %%%%%%%%
      tmp_d_(0) = 3
      tmp_d_(1) = n_source_use
      write(tmp_fname,'(A,A)')
     $     trim(dir_base) , '/source_x0_c_x1_c_x2_c_A___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(2,tmp_d_,source_x0_c_x1_c_x2_c_A__,tmp_fname)
      if (flag_heterogeneity) then
         tmp_d_(0) = 3
         tmp_d_(1) = n_source_use
         write(tmp_fname,'(A,A)')
     $        trim(dir_base) , '/source_x0_c_x1_c_x2_c_B___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,source_x0_c_x1_c_x2_c_B__,tmp_fname)
         tmp_d_(0) = 3
         tmp_d_(1) = n_source_use
         write(tmp_fname,'(A,A)')
     $        trim(dir_base) , '/source_x0_c_x1_c_x2_c_C___.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,source_x0_c_x1_c_x2_c_C__,tmp_fname)
      end if !if (flag_heterogeneity) then
c$$$  %%%%%%%%
c     create physical space grid and sample function_A on it.
      timing_tic = omp_get_wtime()
      quadrature_type_radial = 0
      quadrature_type_azimu_b = 1
      call createfun_xspace(
     $     half_diameter_x_c
     $     ,n_x_c_max
     $     ,x_c_step
     $     ,grid_x0_c___
     $     ,grid_x1_c___
     $     ,grid_x2_c___
     $     ,source_x0_c_x1_c_x2_c_A__
     $     ,n_source_use
     $     ,source_sigma
     $     ,fgauss
     $     ,f_x_c_A___
     $     ,f_x_c_l2_norm
     $     )
      if (verbose.gt.1) write(6,'(A,F8.4)') '  f_x_c_l2_norm is '
     $     ,f_x_c_l2_norm
      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/f_x_c_A___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(3,tmp_d_,f_x_c_A___,tmp_fname)
c     create physical space grid and sample function_B on it.
      if (flag_heterogeneity) then
         quadrature_type_radial = 0
         quadrature_type_azimu_b = 1
         call createfun_xspace(
     $        half_diameter_x_c
     $        ,n_x_c_max
     $        ,x_c_step
     $        ,grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,source_x0_c_x1_c_x2_c_B__
     $        ,n_source_use
     $        ,source_sigma
     $        ,fgauss 
     $        ,f_x_c_B___ 
     $        ,f_x_c_l2_norm
     $        )
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
         call createfun_xspace(
     $        half_diameter_x_c
     $        ,n_x_c_max
     $        ,x_c_step
     $        ,grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,source_x0_c_x1_c_x2_c_C__
     $        ,n_source_use
     $        ,source_sigma
     $        ,fgauss 
     $        ,f_x_c_C___
     $        ,f_x_c_l2_norm
     $        )
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
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- function_sample '
      timing_function_sample = timing_tot
c$$$  %%%%%%%%
      end if !if (nMode.eq.0) then

      tmp_d_(0) = n_x_c_max
      tmp_d_(1) = n_x_c_max
      tmp_d_(2) = n_x_c_max
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x0_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x0_c___,tmp_fname)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x1_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x1_c___,tmp_fname)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_x2_c___.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(tmp_fname)
      call MDA_write_r8(3,tmp_d_,grid_x2_c___,tmp_fname)

c-----------------------------------------------------------------------
c     2) Determine size of spherical grid, allocate space and compute
c     Fourier transform at corresponding points in k-space.
c     Also determine size of spherical harmonic model  (n_Y_lm_sum)
c-----------------------------------------------------------------------
      timing_tic = omp_get_wtime()
      call get_grid_k_p_r_(
     $     k_p_r_max
     $     ,n_k_p_r_max
     $     ,quadrature_type_radial
     $     ,grid_k_p_r_
     $     ,weight_k_p_r_
     $     )
      tmp_d_(0) = n_k_p_r_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/grid_k_p_r_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $     ,trim(tmp_fname)
      call MDA_write_r8(2,tmp_d_,grid_k_p_r_,tmp_fname)
      if (verbose.gt.0) then
         call print_sub_r8(n_k_p_r_max,grid_k_p_r_,' grid_k_p_r_: ')
         call print_sub_r8(n_k_p_r_max,weight_k_p_r_,' weight_k_p_r_: ')
      end if !if (verbose.gt.0) then
      call get_n_polar_a_(
     $     n_k_p_r_max
     $     ,grid_k_p_r_
     $     ,n_polar_a_lowerbound
     $     ,n_polar_a_
     $     ,n_w_
     $     ,n_Y_l_
     $     ,n_Y_lm_
     $     ,n_Y_lm_csum_
     $     ,ld_M
     $     ,n_Y_lm_sum
     $     )
      n_polar_a_min = min_i4_f(n_k_p_r_max,n_polar_a_)
      if (n_polar_a_min.lt.n_polar_a_lowerbound) then
         write(6,'(A,I0,A,I0)') 'Warning, n_polar_a_min ' ,
     $        n_polar_a_min ,' < n_polar_a_lowerbound ' ,
     $        n_polar_a_lowerbound
         stop !exit program due to error.
      end if !if (n_polar_a_min.lt.n_polar_a_lowerbound) then
      n_polar_a_max = max_i4_f(n_k_p_r_max,n_polar_a_)
      if (n_polar_a_max.gt.n_polar_a_upperbound) then
         write(6,'(A,I0,A,I0)') 'Warning, n_polar_a_max ' ,
     $        n_polar_a_max ,' > n_polar_a_upperbound ' ,
     $        n_polar_a_upperbound
         stop !exit program due to error.
      end if !if (n_polar_a_max.gt.n_polar_a_upperbound) then
      tmp_d_(0) = n_k_p_r_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/n_w_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : '
     $     ,trim(tmp_fname)
      call MDA_write_i4(2,tmp_d_,n_w_,tmp_fname)
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_polar_a_
     $     ,' n_polar_a_: ')
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_w_,' n_w_: ')
      if (verbose.gt.0) write(6,'(A,I0)') '  ld_M: ',ld_M
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_Y_l_
     $     ,' n_Y_l_: ')
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_Y_lm_,
     $     ' n_Y_lm_: ')
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_Y_lm_csum_,
     $     ' n_Y_lm_csum_: ')
      ld_S = ld_M
      call get_n_azimu_b__(
     $     n_k_p_r_max
     $     ,n_polar_a_upperbound
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,n_azimu_b__
     $     ,n_azimu_b__sum_
     $     ,n_azimu_b__sum_sum
     $     )
      if (verbose.gt.0) call print_sub_i4(n_k_p_r_max,n_azimu_b__sum_
     $     ,' n_azimu_b__sum_: ')
      if (verbose.gt.0) write(6,'(A,I0)') '  n_azimu_b__sum_sum is '
     $     ,n_azimu_b__sum_sum
c
      allocate(grid_k0_c_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_r8(n,grid_k0_c_)
      allocate(grid_k1_c_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_r8(n,grid_k1_c_)
      allocate(grid_k2_c_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_r8(n,grid_k2_c_)
      allocate(weight_k_c_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_r8(n,weight_k_c_)
      allocate(f_k_c_A_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_c16(n,f_k_c_A_)
      if (flag_heterogeneity) then
         allocate(f_k_c_B_(0:1+n_azimu_b__sum_sum-1))
         n = n_azimu_b__sum_sum
         call cs1_c16(n,f_k_c_B_)
         allocate(f_k_c_C_(0:1+n_azimu_b__sum_sum-1))
         n = n_azimu_b__sum_sum
         call cs1_c16(n,f_k_c_C_)
      end if !if (flag_heterogeneity) then
      allocate(f_k_c_tru_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_c16(n,f_k_c_tru_)
      allocate(f_k_c_est_(0:1+n_azimu_b__sum_sum-1))
      n = n_azimu_b__sum_sum
      call cs1_c16(n,f_k_c_est_)
      allocate(Y_A_(0:1+n_Y_lm_sum-1))
      n = n_Y_lm_sum
      call cs1_c16(n,Y_A_)
      if (flag_heterogeneity) then
         allocate(Y_B_(0:1+n_Y_lm_sum-1))
         n = n_Y_lm_sum
         call cs1_c16(n,Y_B_)
         allocate(Y_C_(0:1+n_Y_lm_sum-1))
         n = n_Y_lm_sum
         call cs1_c16(n,Y_C_)
      end if !if (flag_heterogeneity) then
      allocate(Y_tru_(0:1+n_Y_lm_sum-1))
      n = n_Y_lm_sum
      call cs1_c16(n,Y_tru_)
      allocate(Y_est_(0:1+n_Y_lm_sum-1))
      n = n_Y_lm_sum
      call cs1_c16(n,Y_est_)

c
c     create kspace grid data from physical space model A
c
      eps_default = 1.0d-6
      call fun_to_kspacegrid(
     $     grid_x0_c___
     $     ,grid_x1_c___
     $     ,grid_x2_c___
     $     ,f_x_c_A___
     $     ,x_c_step
     $     ,n_x_c_max
     $     ,eps_default
     $     ,k_p_r_max
     $     ,n_k_p_r_max
     $     ,quadrature_type_radial
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,n_azimu_b__sum_sum
     $     ,numonsphere 
     $     ,grid_k0_c_
     $     ,grid_k1_c_ 
     $     ,grid_k2_c_
     $     ,weight_k_c_
     $     ,f_k_c_A_
     $     ,f_k_c_l2_norm
     $     ,ier
     $     )
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k0_c_
     $     ,grid_k0_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k1_c_
     $     ,grid_k1_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k2_c_
     $     ,grid_k2_c_)
      if (verbose.gt.1) then
         call print_sub_r8(n_x_c_max,grid_x0_c___,' grid_x0_c___: ')
         call print_sub_r8(n_x_c_max,grid_x1_c___,' grid_x1_c___: ')
         call print_sub_r8(n_x_c_max,grid_x2_c___,' grid_x2_c___: ')
         call print_sub_r8(n_azimu_b__sum_sum,grid_k0_c_
     $        ,' grid_k0_c_: ')
         call print_sub_r8(n_azimu_b__sum_sum,grid_k1_c_
     $        ,' grid_k1_c_: ')
         call print_sub_r8(n_azimu_b__sum_sum,grid_k2_c_
     $        ,' grid_k2_c_: ')
      end if !if (verbose.gt.1) then
      call cp1_i4(n_k_p_r_max,numonsphere,n_azimu_b__sum_)
      if (verbose.gt.1) then
         write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
         call print_all_i4(n_k_p_r_max,numonsphere,'numonsphere: ')
         call print_all_i4(n_k_p_r_max,numonsphere
     $        ,' n_azimu_b__sum_: ')
      end if
c     create kspace grid data from physical space models B,C
c
      if (flag_heterogeneity) then
         eps_default = 1.0d-6
         call fun_to_kspacegrid(
     $        grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,f_x_c_B___
     $        ,x_c_step
     $        ,n_x_c_max
     $        ,eps_default
     $        ,k_p_r_max
     $        ,n_k_p_r_max 
     $        ,quadrature_type_radial
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b 
     $        ,n_azimu_b__sum_sum
     $        ,numonsphere 
     $        ,grid_k0_c_ 
     $        ,grid_k1_c_ 
     $        ,grid_k2_c_
     $        ,weight_k_c_
     $        ,f_k_c_B_
     $        ,f_k_c_l2_norm 
     $        ,ier
     $        )
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k0_c_
     $     ,grid_k0_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k1_c_
     $     ,grid_k1_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k2_c_
     $     ,grid_k2_c_)
         if (verbose.gt.1) then
            write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
            call print_all_i4(n_k_p_r_max,numonsphere,'numonsphere: ')
         end if !if (verbose.gt.1) then
         eps_default = 1.0d-6
         call fun_to_kspacegrid(
     $        grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,f_x_c_C___
     $        ,x_c_step
     $        ,n_x_c_max
     $        ,eps_default
     $        ,k_p_r_max
     $        ,n_k_p_r_max 
     $        ,quadrature_type_radial
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b 
     $        ,n_azimu_b__sum_sum
     $        ,numonsphere 
     $        ,grid_k0_c_ 
     $        ,grid_k1_c_ 
     $        ,grid_k2_c_
     $        ,weight_k_c_
     $        ,f_k_c_C_
     $        ,f_k_c_l2_norm 
     $        ,ier
     $        )
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k0_c_
     $     ,grid_k0_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k1_c_
     $     ,grid_k1_c_)
      call af1_r8(n_azimu_b__sum_sum,2.0d0*pi,0.0d0,grid_k2_c_
     $     ,grid_k2_c_)
         if (verbose.gt.1) then
            write(6,'(A,F8.4)') 'f_k_c_l2_norm: ',f_k_c_l2_norm
            call print_all_i4(n_k_p_r_max,numonsphere,'numonsphere: ')
         end if !if (verbose.gt.1) then
      end if !if (flag_heterogeneity) then=
c
c     convert to model_A: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(
     $     f_k_c_A_
     $     ,k_p_r_max
     $     ,n_Y_l_
     $     ,quadrature_type_radial
     $     ,n_k_p_r_max
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,n_azimu_b__sum_
     $     ,Y_A_
     $     ,n_Y_lm_sum
     $     )
c     convert to model_B,C: spherical harmonic expansions on successive 
c     spheres
c
      if (flag_heterogeneity) then
         call kspacegrid_to_model(
     $        f_k_c_B_
     $        ,k_p_r_max
     $        ,n_Y_l_
     $        ,quadrature_type_radial
     $        ,n_k_p_r_max
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b
     $        ,n_azimu_b__sum_
     $        ,Y_B_
     $        ,n_Y_lm_sum
     $        )
         call kspacegrid_to_model(
     $        f_k_c_C_
     $        ,k_p_r_max
     $        ,n_Y_l_
     $        ,quadrature_type_radial
     $        ,n_k_p_r_max
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b
     $        ,n_azimu_b__sum_
     $        ,Y_C_
     $        ,n_Y_lm_sum
     $        )
      end if !if (flag_heterogeneity) then
c$$$      save both Y_A_ and Y_B_ Y_C_ to disk.
      tmp_d_(0) = n_k_p_r_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/n_Y_lm_csum__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(tmp_fname)
      call MDA_write_i4(2,tmp_d_,n_Y_lm_csum_,tmp_fname)
      tmp_d_(0) = n_k_p_r_max
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/n_Y_l__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(tmp_fname)
      call MDA_write_i4(2,tmp_d_,n_Y_l_,tmp_fname)
      tmp_d_(0) = n_Y_lm_sum
      tmp_d_(1) = 1
      write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_A_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(2,tmp_d_,Y_A_,tmp_fname)
      if (flag_heterogeneity) then
         tmp_d_(0) = n_Y_lm_sum
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_B_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_B_,tmp_fname)
         tmp_d_(0) = n_Y_lm_sum
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_C_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_C_,tmp_fname)
      end if !if (flag_heterogeneity) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- Y_create '
      timing_Y_create = timing_tot

c-----------------------------------------------------------------------
c     3) Assign slices to orientations of spherical grid on highest
c        frequency sphere and define alphas.
c-----------------------------------------------------------------------
      timing_tic = omp_get_wtime()
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
      n_M = n_azimu_b__sum_(n_k_p_r_max-1)
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
      if (nMode.eq.2) then
      open(20,FILE=fname_dim_A)
      read(20,*) n_x_c_from_file,n_M_from_file
      close(20)
      n_M = n_M_from_file
      if ((n_M_sample_max_set.gt.0) .and.
     $     (n_M_sample_max_set.lt.n_M_from_file)) then
         write(6,'(A,I0)') ' reducing n_M to : ' , n_M_sample_max_set
         n_M = min(n_M_from_file,n_M_sample_max_set)
      end if !if (n_M_sample_max_set.gt.0) then
      end if !if (nMode.eq.2) then
      n_S = n_azimu_b__sum_(n_k_p_r_max-1)
      call getspheregrid(
     $     n_polar_a_(n_k_p_r_max-1) !integer *4: number of polar_a values on sphere of radius grid_k_p_r_(nk). ;
     $     ,quadrature_type_azimu_b !integer: quadrature type in the azimuthal direction. sometimes named itypep. 0=same number of points for all n_polar_a_ (oversampling poles). 1=adaptive (samples with roughly uniform density). ;
     $     ,grid_cos_polar_a_ !temporary: real *8 array (size at least n_k_p_r_max): values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called xnodesth(nk). ;
     $     ,grid_sin_polar_a_ !temporary: real *8 array (size at least n_k_p_r_max): values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called sthetas(nk). ;
     $     ,weight_cos_polar_a_ !temporary real *8 array (size at least n_k_p_r_max): weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called wtsth_(np). ;
     $     ,n_azimu_b_ !temporary integer *4 array (size at least n_k_p_r_max): number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_r_(nk). also called ngridps(np). ;
     $     ,azimu_b_step_ !temporary real *8 array (size at least n_k_p_r_max): grid-spacing for azimu_b. ;
     $     ,n_azimu_b__sum_max !temporary integer *4: total number of points on sphere at a particular radius grid_k_p_r_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
     $     )
      if (verbose.gt.0) write(6,'(A,I0)') '  n_M is ',n_M
      if (verbose.gt.0) then
         call print_all_i4(n_k_p_r_max,n_w_
     $        ,' n_w_: ')
         call print_all_i4(n_k_p_r_max,n_azimu_b__sum_
     $        ,' n_azimu_b__sum_: ')
      end if! if (verbose.gt.0) then
      if (verbose.gt.0) write(6,'(A,I0)')
     $     '  n_azimu_b__sum_sum is '
     $     ,n_azimu_b__sum_sum
      if (verbose.gt.0) write(6,'(A,I0)')
     $     '  n_azimu_b__sum_max is ',n_azimu_b__sum_max
c      
      if (flag_displacement) then
         displacement_max = displacement_max_set
      else
         displacement_max = 0.0d0
      end if

c
      allocate(alpha2d_tru__(0:n_alpha-1,0:1+n_M-1))
      n = n_alpha*n_M
      call cs1_r8(n,alpha2d_tru__)
      allocate(alpha2d_est__(0:n_alpha-1,0:1+n_M-1))
      n = n_alpha*n_M
      call cs1_r8(n,alpha2d_est__)

c
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
      nA = 0
      do npolar_a = 0,n_polar_a_(n_k_p_r_max-1)-1
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
               alpha2d_tru__(nalpha_delta_x,nA) = 0.25d0
     $              *displacement_max*(mod(1+nA,5)-2)
               alpha2d_tru__(nalpha_delta_y,nA) = 0.50d0
     $              *displacement_max*(mod(1+nA,3)-1)
            else
               alpha2d_tru__(nalpha_delta_x,nA) = 0.0d0*displacement_max
               alpha2d_tru__(nalpha_delta_y,nA) = 0.0d0*displacement_max
            end if !flag_displacement
            if (flag_normalization) then
c$$$               Below we draw v from [0.5,1.5]
               alpha2d_tru__(nalpha_l2_norm,nA) = 0.5d0 +
     $              adi_rand_f(rseed)
            else
               alpha2d_tru__(nalpha_l2_norm,nA) = 1.0d0
            end if !flag_normalization
c$$$            to start with somewhat arbitrary alpha2d_est__
            alpha2d_tru__(nalpha_ctf_ind,nA) = 0.0d0 + mod(0+nA,3)
            alpha2d_tru__(nalpha_S_index,nA) = 0.0d0
            alpha2d_tru__(nalpha_M_index,nA) = 0.0d0
c$$$            estimated parameters used to initialize least-squares-solve
            alpha2d_est__(nalpha_polar_a,nA) = pi *adi_rand_f(rseed)
            alpha2d_est__(nalpha_azimu_b,nA) = 2 *adi_rand_f(rseed)*pi
            alpha2d_est__(nalpha_gamma_z,nA) = pi /2.0d0*(2*mod(1+nA,2)
     $           -1)
            if (flag_displacement) then
               alpha2d_est__(nalpha_delta_x,nA) = 0.125d0
     $              *displacement_max*(mod(1+nA+2,5)-2)
               alpha2d_est__(nalpha_delta_y,nA) = 0.25d0
     $              *displacement_max*(mod(1+nA+1,3)-1)
            else
               alpha2d_est__(nalpha_delta_x,nA) = 0.0d0
               alpha2d_est__(nalpha_delta_y,nA) = 0.0d0
            end if !flag_displacement
            alpha2d_est__(nalpha_l2_norm,nA) = 1.0d0
            alpha2d_est__(nalpha_ctf_ind,nA) = 0.0d0 + mod(0+nA,3)
            alpha2d_est__(nalpha_S_index,nA) = 0.0d0
            alpha2d_est__(nalpha_M_index,nA) = 0.0d0
            nA = nA+1
         enddo !do nazimu_b = 0,n_azimu_b_(npolar_a)-1
      enddo !do npolar_a = 0,n_polar_a_(n_k_p_r_max-1)-1
      n_azimu_b__sum_max = nA
      if (verbose.gt.1) write(6,'(A,I0)')
     $     '  n_azimu_b__sum_max is ',n_azimu_b__sum_max
      if (n_azimu_b__sum_max.ne.n_azimu_b__sum_(n_k_p_r_max
     $     -1))then
         write(6,'(A,I0,A,I0,A,I0,A)')
     $        'Warning, n_azimu_b__sum_max '
     $        ,n_azimu_b__sum_max,'.ne.',n_azimu_b__sum_
     $        ,'(' ,n_k_p_r_max-1,')'
      end if
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

      if (nMode.eq.2) then
      open(12,FILE=fname_euler_angle_)
      do nM=0,n_M-1
         if ((mod(nM,1000).eq.0) .or. (nM.eq.n_M-1)) then
            write(6,'(A,I0,A,I0)') ' nM: ' , nM , '/' , n_M
         end if !if (mod(nM,100).eq.0) then
         read(12,*) (tmp_euler_angle_(n),n=0,2)
         call af1_r8(
     $        3
     $        ,pi/180.0d0
     $        ,0.0d0
     $        ,tmp_euler_angle_
     $        ,tmp_euler_angle_
     $        )
         alpha2d_tru__(nalpha_polar_a,nM) = tmp_euler_angle_(0)
         alpha2d_tru__(nalpha_azimu_b,nM) = tmp_euler_angle_(1)
         alpha2d_tru__(nalpha_gamma_z,nM) = tmp_euler_angle_(2)
         if ((verbose.gt.1) .and. (nM.lt.6)) then
            write(6,'(A,F8.3)') ' polar_a: ' ,
     $           alpha2d_tru__(nalpha_polar_a,nM)
            write(6,'(A,F8.3)') ' azimu_b: ' ,
     $           alpha2d_tru__(nalpha_azimu_b,nM)
            write(6,'(A,F8.3)') ' gamma_z: ' ,
     $           alpha2d_tru__(nalpha_gamma_z,nM)
         end if !if (nM.lt.6) then
c$$$  Now convert the euler-angles from relion-convention to our convention. ;
         call relion_euler_rot_mat(
     $        alpha2d_tru__(nalpha_polar_a,nM)
     $        ,alpha2d_tru__(nalpha_azimu_b,nM)
     $        ,alpha2d_tru__(nalpha_gamma_z,nM)
     $        ,tmp_rotation_matrix_
     $        )
         call rot_mat_to_euler(
     $        tmp_rotation_matrix_
     $        ,tmp_euler_angle_
     $        )
         alpha2d_tru__(nalpha_polar_a,nM) = dacos(tmp_euler_angle_(0))
         alpha2d_tru__(nalpha_azimu_b,nM) = tmp_euler_angle_(1)
         alpha2d_tru__(nalpha_gamma_z,nM) = tmp_euler_angle_(2)
         read(12,*) (tmp_shift_(n),n=0,1)
         alpha2d_tru__(nalpha_delta_x,nM) = 2.0d0*pi*tmp_shift_(0)
         alpha2d_tru__(nalpha_delta_y,nM) = 2.0d0*pi*tmp_shift_(1)
         if ((verbose.gt.1) .and. (nM.lt.6)) then
            write(6,'(A,F8.3)') ' delta_x: ' ,
     $           alpha2d_tru__(nalpha_delta_x,nM)
            write(6,'(A,F8.3)') ' delta_y: ' ,
     $           alpha2d_tru__(nalpha_delta_y,nM)
         end if !if (nM.lt.6) then
c     the translations have to be scaled, because in Brubaker the box size is [-64:64),
c     whereas our box size is [-1:1)
         alpha2d_tru__(nalpha_delta_x,nM) = alpha2d_tru__(nalpha_delta_x
     $        ,nM)*(2.0d0/n_x_c_from_file)
         alpha2d_tru__(nalpha_delta_y,nM) = alpha2d_tru__(nalpha_delta_y
     $        ,nM)*(2.0d0/n_x_c_from_file)
         alpha2d_tru__(nalpha_l2_norm,nM) = 1.0d0
         alpha2d_tru__(nalpha_ctf_ind,nM) = 0.0d0
         alpha2d_tru__(nalpha_S_index,nM) = 0.0d0
         alpha2d_tru__(nalpha_M_index,nM) = 0.0d0
         call convert_gamma_delta_marina_to_adi_0(
     $     alpha2d_tru__(nalpha_gamma_z,nM)
     $     ,alpha2d_tru__(nalpha_delta_x,nM)
     $     ,alpha2d_tru__(nalpha_delta_y,nM)
     $     ,alpha2d_tru__(nalpha_gamma_z,nM)
     $     ,alpha2d_tru__(nalpha_delta_x,nM)
     $     ,alpha2d_tru__(nalpha_delta_y,nM)
     $     )
         if (verbose.gt.2) then
            call print_all_r8(n_alpha,alpha2d_tru__(0,nM)
     $        ,' alpha2d_tru__: ')
         end if !if (verbose.gt.2) then
      enddo !do nM=0,n_M-1
      close(12)
      do nM=0,n_M-1
      alpha2d_est__(nalpha_polar_a,nM) = pi *adi_rand_f(rseed)
      alpha2d_est__(nalpha_azimu_b,nM) = 2 *adi_rand_f(rseed)*pi
      alpha2d_est__(nalpha_gamma_z,nM) = pi /2.0d0*(2*mod(1+nM,2)-1)
      alpha2d_est__(nalpha_delta_x,nM) = 0.0d0
      alpha2d_est__(nalpha_delta_y,nM) = 0.0d0
      alpha2d_est__(nalpha_l2_norm,nM) = 1.0d0
      alpha2d_est__(nalpha_ctf_ind,nM) = 0.0d0
      alpha2d_est__(nalpha_S_index,nM) = 0.0d0
      alpha2d_est__(nalpha_M_index,nM) = 0.0d0
      enddo !do nM=0,n_M-1
      open(20,FILE=fname_I_ctf_);
      do nM=0,n_M-1
         read(20,*) I_ctf_from_file !Note that I_ctf_from_file is a real *8, not an integer *4. ;
         alpha2d_tru__(nalpha_ctf_ind,nM) = 1.0d0*I_ctf_from_file-1.0d0
         alpha2d_est__(nalpha_ctf_ind,nM) = 1.0d0*I_ctf_from_file-1.0d0
      enddo !do nM=0,n_M-1
      close(20)
      end if !if (nMode.eq.2) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- alpha_define '
      timing_alpha_define = timing_tot

c
      allocate(Y_slice__(0:ld_M-1,0:1+n_M-1))
      n = ld_M*n_M
      call cs1_c16(n,Y_slice__)
      call get_template_size(
     $     n_polar_a_
     $     ,n_k_p_r_max
     $     ,ld_S
     $     ,n_w_tmp_
     $     ,n_w_csum_
     $     )
      if (verbose.gt.1) write(6,'(A,I0)') '  ld_S is ',ld_S
c$$$  Checking values of n_w_csum_. ;
      nA = 0
      n_w_sum = 0
      do nk=0,n_k_p_r_max-1
         if (n_w_csum_(nk)-1.ne.nA) then
            write(6,'(A,I0,A,I0,A,I0)') 'Warning, n_w_csum_(',nk
     $           ,')-1 = ',n_w_csum_(nk)-1,' .ne. ',nA
            stop !exit program due to error.
         end if !if (n_w_csum_(nk).ne.nA) then
         nA = nA + n_w_(nk)
         n_w_sum = n_w_sum + n_w_(nk)
      enddo !do nk=0,n_k_p_r_max-1;
      allocate(weight_k_p_(0:1+n_w_sum-1))
      call cs1_r8(n_w_sum,weight_k_p_)
      na=0
      do nk=0,n_k_p_r_max-1
         do nw=0,n_w_(nk)-1
            weight_k_p_(na) = weight_k_p_r_(nk)*2.0d0*pi/n_w_(nk)
            na=na+1
         enddo !do nw=0,n_w_(nk)-1
      enddo !do nk=0,n_k_p_r_max-1
c
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
      allocate(Y_slice_A__(0:ld_M-1,0:1+n_M-1))
      n = ld_M*n_M
      call cs1_c16(n,Y_slice_A__)      
      if (flag_heterogeneity) then
         allocate(Y_slice_B__(0:ld_M-1,0:1+n_M-1))
         n = ld_M*n_M
         call cs1_c16(n,Y_slice_B__)
         allocate(Y_slice_C__(0:ld_M-1,0:1+n_M-1))
         n = ld_M*n_M
         call cs1_c16(n,Y_slice_B__)
      end if !if (flag_heterogeneity) then
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c-----------------------------------------------------------------------
c     4) compute Y_slice_A__ directly from physical space function (f_x_c_A___).
c        with orientation vectors defined by alpha2d_tru__.
c-----------------------------------------------------------------------
c
      timing_tic = omp_get_wtime()
      call mk_simulated_slices_alpha2d(
     $     f_x_c_A___
     $     ,n_x_c_max
     $     ,diameter_x_c
     $     ,eps_default
     $     ,n_k_p_r_max
     $     ,quadrature_type_radial
     $     ,n_w_
     $     ,k_p_r_max
     $     ,alpha2d_tru__
     $     ,n_M 
     $     ,ld_M
     $     ,Y_slice_A__
     $     )
c-----------------------------------------------------------------------
c     4) compute Y_slice_B__,C__ directly from physical space function 
c        (f_x_c_B___,C___). with orientation vectors defined by alpha2d_tru__.
c-----------------------------------------------------------------------
c
      if (flag_heterogeneity) then
         call mk_simulated_slices_alpha2d(
     $        f_x_c_B___
     $        ,n_x_c_max
     $        ,diameter_x_c
     $        ,eps_default
     $        ,n_k_p_r_max
     $        ,quadrature_type_radial 
     $        ,n_w_
     $        ,k_p_r_max
     $        ,alpha2d_tru__
     $        ,n_M
     $        ,ld_M 
     $        ,Y_slice_B__
     $        )
         call mk_simulated_slices_alpha2d(
     $        f_x_c_C___
     $        ,n_x_c_max
     $        ,diameter_x_c
     $        ,eps_default
     $        ,n_k_p_r_max
     $        ,quadrature_type_radial 
     $        ,n_w_
     $        ,k_p_r_max
     $        ,alpha2d_tru__
     $        ,n_M
     $        ,ld_M 
     $        ,Y_slice_C__
     $        )
      end if !if (flag_heterogeneity) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- Y_slice_define '
      timing_Y_slice_define = timing_tot
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

      allocate(I_model_(0:1+n_M-1))
      call cs1_i4(n_M,I_model_)
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Determine which images correspond to which model.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
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
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
      if (nMode.eq.2) then
      do nm=0,n_M-1
         I_model_(nm) = 0
      enddo !do nm=0,n_M-1
      end if !if (nMode.eq.2) then

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
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
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
      if (nMode.eq.2) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Read Y_slice__ from disk 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(M_x_c_from_file__(0:n_x_c_from_file-1,0:1
     $        +n_x_c_from_file-1))
      call cs1_c16(n_x_c_from_file**2,M_x_c_from_file__)
      open(20,FILE=fname_image_A,FORM="unformatted");
      do nM=0,n_M-1
      do ny_c=0,n_x_c_from_file-1
      read(20) (M_x_c_from_file__(nx_c,ny_c),nx_c=0,n_x_c_from_file-1)
      enddo !do ny_c=0,n_x_c-1
      call interp_x_c_to_k_p_finufft(
     $     n_x_c_from_file
     $     ,diameter_x_c
     $     ,n_x_c_from_file
     $     ,diameter_x_c
     $     ,M_x_c_from_file__
     $     ,n_k_p_r_max
     $     ,grid_k_p_r_ 
     $     ,n_w_
     $     ,Y_slice__(0,nM)
     $     )
      enddo !do nM=0,n_M-1
      close(20)
      tmp_d_(0) = n_w_sum
      tmp_d_(1) = min(3,n_M)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/Y_slice__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(2,tmp_d_,Y_slice__,tmp_fname)
      end if !if (nMode.eq.2) then

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$  if flag_SNR, then add noise to slices. 
c$$$c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if ((flag_snr.eqv..true.) .and. snr.gt.0.0d0) then
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
            call additive_noise_k_p(
     $           rseed
     $           ,n_k_p_r_max
     $           ,grid_k_p_r_
     $           ,n_w_
     $           ,ld_M,Y_slice__(0,nm)
     $           ,tmp_sigma
     $           )
         enddo !do nm=0,n_M-1
      end if !if (flag_snr.eqv..true.) then
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

c-----------------------------------------------------------------------
c     5) generate ctf functions to associate with slices.
c-----------------------------------------------------------------------
c

c$$$  Number of different ctfs used is n_CTF.
c$$$  Ultimately, we expect this to be proportional to the number of images,
c$$$  with batches of roughly 30-300 images each having the same ctf.
c$$$  %%%%%%%%;
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
      n_CTF = 6
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  %%%%%%%%;
      if (nMode.eq.2) then
      open(20,FILE=fname_n_ctf);
      read(20,*) n_CTF_from_file
      close(20)
      n_CTF = n_CTF_from_file
      end if !if (nMode.eq.2) then
c$$$  %%%%%%%%;
      ld_CTF = ld_M
      allocate(CTF_Voltage_kV_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Voltage_kV_)
      allocate(CTF_Defocus_U_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Defocus_U_)
      allocate(CTF_Defocus_V_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Defocus_V_)
      allocate(CTF_Defocus_Angle_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Defocus_Angle_)
      allocate(CTF_Spherical_Aberration_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Spherical_Aberration_)
      allocate(CTF_Detector_Pixel_Size_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Detector_Pixel_Size_)
      allocate(CTF_Magnification_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Magnification_)
      allocate(CTF_Amplitude_Contrast_(0:1+n_CTF-1))
      call cs1_r8(n_CTF,CTF_Amplitude_Contrast_)
      allocate(CTF_k_p__(0:1+n_CTF*ld_CTF-1))
      n = n_CTF*ld_CTF
      call cs1_c16(n,CTF_k_p__)
c$$$  %%%%%%%%%%%%%%%%
      if (nMode.eq.2) then
      open(20,FILE=fname_microscope_parameter_);
      read(20,*) CTF_Voltage_kV,CTF_Spherical_Aberration
     $     ,CTF_Detector_Pixel_Size,CTF_Amplitude_Contrast
      close(20)
      CTF_B_factor = 50
      CTF_Magnification = 1.0d0
      open(20,FILE=fname_ctf_parameter_);
      do nctf=0,n_CTF-1
         read(20,*) CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle
         CTF_Voltage_kV_(nctf) = CTF_Voltage_kV
         CTF_Defocus_U_(nctf) = CTF_Defocus_U
         CTF_Defocus_V_(nctf) = CTF_Defocus_V
         CTF_Defocus_Angle_(nctf) = CTF_Defocus_Angle
         CTF_Spherical_Aberration_(nctf) = CTF_Spherical_Aberration
         CTF_Detector_Pixel_Size_(nctf) = CTF_Detector_Pixel_Size
         CTF_Magnification_(nctf) = CTF_Magnification
         CTF_Amplitude_Contrast_(nctf) = CTF_Amplitude_Contrast
      enddo !do nctf=0,n_CTF-1
      close(20)
      end if !if (nMode.eq.2) then
c$$$  %%%%%%%%%%%%%%%%
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  Here we use an unrealistic ctf-function with a high degree
c$$$  of anisotropy to test our code.
      if (ctf_type.eq.0) then
         do nctf=0,n_CTF-1
            call get_ctf_ones_k_p_(n_k_p_r_max,grid_k_p_r_,n_w_,ld_CTF
     $           ,CTF_k_p__(nctf*ld_CTF))
         enddo !do nctf=0,n_CTF-1
      end if !if (ctf_type.eq.0) then
      if (ctf_type.eq.1) then
         do nctf=0,n_CTF-1
            ctf_p_0 = 10.0d0 + (5.0d0*nctf)/n_CTF
            ctf_p_1 = (pi*nctf)/n_CTF
            ctf_p_2 = 2
            call get_ctf_star_k_p_(n_k_p_r_max,grid_k_p_r_,n_w_,ld_CTF
     $           ,CTF_k_p__(nctf*ld_CTF),ctf_p_0,ctf_p_1,ctf_p_2)
         enddo !do nctf=0,n_CTF-1
      end if !if (ctf_type.eq.1) then
      if (ctf_type.eq.2) then
c$$$  the CTF parameters are taken from the rib80s dataset. ;
      CTF_Voltage_kV = 300.0d0
      CTF_Spherical_Aberration = 2.0d0
      CTF_Detector_Pixel_Size = 3.7687499999999998d0
      CTF_Amplitude_Contrast = 0.1d0
      CTF_B_factor = 50
      CTF_Magnification = 1.0d0 !not used yet ;
      CTF_Defocus_max = 26462.15039099d0 !maximum possible value of defocus ;
      CTF_Defocus_min = 8112.450195d0 !minimum possible value of defocus ;
      CTF_Astigmatism_max = 457.03125d0 !max astigmatism (CTF_Defocus_U-CTF_Defocus_V) ;
      CTF_Astigmatism_min = 10.509765d0 !min astigmatism ;
      do nctf=0,n_CTF-1
c$$$  generate a value between the min and max defocus ;
         CTF_Defocus_V = adi_rand_f(rseed)*(CTF_Defocus_max
     $        -CTF_Defocus_min)
         CTF_Defocus_V = CTF_Defocus_V + CTF_Defocus_min
c$$$  generate a value between max and min astigmatism ;
         CTF_Astigmatism_value = adi_rand_f(rseed)*(CTF_Astigmatism_max
     $        -CTF_Astigmatism_min)
         CTF_Astigmatism_value = CTF_Astigmatism_value +
     $        CTF_Astigmatism_min
c$$$  ensure that CTF_Defocus_U > CTF_Defocus_U ;
         CTF_Defocus_U = CTF_Defocus_V+CTF_Astigmatism_value
         CTF_Voltage_kV_(nctf) = CTF_Voltage_kV
         CTF_Defocus_U_(nctf) = CTF_Defocus_U
         CTF_Defocus_V_(nctf) = CTF_Defocus_V
         CTF_Defocus_Angle_(nctf) = CTF_Defocus_Angle
         CTF_Spherical_Aberration_(nctf) = CTF_Spherical_Aberration
         CTF_Detector_Pixel_Size_(nctf) = CTF_Detector_Pixel_Size
         CTF_Magnification_(nctf) = CTF_Magnification
         CTF_Amplitude_Contrast_(nctf) = CTF_Amplitude_Contrast
      enddo !do nctf=0,n_CTF-1
      end if !if (ctf_type.eq.2) then
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  %%%%%%%%%%%%%%%%
      
c$$$  %%%%%%%%%%%%%%%%
      if (ctf_type.eq.2) then
c$$$  Generate CTF functions to associate with slices. ;
      do nctf=0,n_CTF-1
c$$$  spherical aberation of the lens in mm ;
         CTF_Spherical_Aberration = CTF_Spherical_Aberration_(nctf)
c$$$  convert into Angstroms ;
         CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0
     $        **7.0d0)
c$$$  voltage in kVolts ;
         CTF_Voltage_kV = CTF_Voltage_kV_(nctf)
c$$$  convert into Volts ;
         CTF_Voltage_1V=CTF_Voltage_kV*1000.0 
c$$$  electron wavelength in Angstroms ;
         CTF_lambda = 12.2643247/dsqrt(CTF_Voltage_1V+CTF_Voltage_1V**2
     $        *0.978466d-6)
c$$$  defocus values (in Angstroms) ;
         CTF_Defocus_U = CTF_Defocus_U_(nctf)
         CTF_Defocus_V = CTF_Defocus_V_(nctf)
c$$$  angle of astigmatism ;
         CTF_Defocus_Angle = CTF_Defocus_Angle_(nctf)
c$$$  convert into radians ;
         CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0
c$$$  CTF_Amplitude Contrast ;
         CTF_Amplitude_Contrast = CTF_Amplitude_Contrast_(nctf)
c$$$  weights for the amplitude and phase contrasts in CTF ;
         tmp_w1=dsqrt(1.0d0-CTF_Amplitude_Contrast**2)
         tmp_w2=CTF_Amplitude_Contrast
c$$$  pixel size of the scanner in physical space (not magnified) in Angstroms ;
c$$$   CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification
      CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size_(nctf)
c$$$  n_x_c_max*CTF_Object_Pixel_Size is the box size in Angstroms ;
      CTF_lambda_per_box = CTF_lambda/(n_x_c_from_file
     $     *CTF_Object_Pixel_Size)
c$$$   call envelope_fxn(ngridr,xnodesr/pi,D,envelope)
      na=0
      do nk = 0,n_k_p_r_max-1
         do nw = 0,n_w_(nk)-1
            tmp_theta = (2.0d0*pi*nw)/n_w_(nk)
            tmp_k_c_1 = (2.0d0*pi)*grid_k_p_r_(nk)*dcos(tmp_theta)
            tmp_k_c_2 = (2.0d0*pi)*grid_k_p_r_(nk)*dsin(tmp_theta)
            call niko_ctf(
     $           CTF_Spherical_Aberration
     $           ,CTF_lambda
     $           ,tmp_w1
     $           ,tmp_w2
     $           ,CTF_Defocus_U
     $           ,CTF_Defocus_V
     $           ,CTF_Defocus_Angle
     $           ,CTF_lambda_per_box
     $           ,tmp_k_c_1/pi
     $           ,tmp_k_c_2/pi
     $           ,tmp_ctf_value
     $           )
            if (verbose.gt.3) then
            if (nctf.eq.0) then
            write(6,'(2(A,I0),A,F24.4,10(A,F24.16))')
     $           ' nk: ' , nk
     $           ,' nw: ' , nw
     $           ,' CTF_Spherical_Aberration: '
     $           ,CTF_Spherical_Aberration
     $           ,' CTF_lambda: '
     $           ,CTF_lambda
     $           ,' tmp_w1: '
     $           ,tmp_w1
     $           ,' tmp_w2: '
     $           ,tmp_w2
     $           ,' CTF_Defocus_U: '
     $           ,CTF_Defocus_U
     $           ,' CTF_Defocus_V: '
     $           ,CTF_Defocus_V
     $           ,' CTF_Defocus_Angle: '
     $           ,CTF_Defocus_Angle
     $           ,' CTF_lambda_per_box: '
     $           ,CTF_lambda_per_box
     $           ,' tmp_k_c_1/pi: '
     $           ,tmp_k_c_1/pi
     $           ,' tmp_k_c_2/pi: '
     $           ,tmp_k_c_2/pi
     $           ,' tmp_ctf_value: '
     $           ,tmp_ctf_value
            end if !if (nctf.eq.0) then
            end if !if (verbose.gt.3) then
            CTF_k_p__(na + nctf*ld_CTF) = -tmp_ctf_value !Note change in sign
            na = na+1
         enddo !do nw = 0,n_w_(nk)-1
      enddo !do nk = 0,n_k_p_r_max-1
      enddo !do nctf=0,n_CTF-1
      if (verbose.gt.2) then
      call print_sub_c16(n_w_sum,CTF_k_p__(0*ld_CTF)
     $     ,' CTF_k_p__ 0: ')
      call print_sub_c16(n_w_sum,CTF_k_p__(1*ld_CTF)
     $     ,' CTF_k_p__ 1: ')
      call print_sub_c16(n_w_sum,CTF_k_p__(2*ld_CTF)
     $     ,' CTF_k_p__ 2: ')
      end if !if (verbose.gt.2) then
      tmp_d_(0) = n_w_sum
      tmp_d_(1) = min(3,n_CTF)
      write(tmp_fname,'(A,A)') trim(dir_base) , '/CTF_k_p__.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $     ,trim(tmp_fname)
      call MDA_write_c16(2,tmp_d_,CTF_k_p__,tmp_fname)
      end if !if (ctf_type.eq.2) then
c$$$  %%%%%%%%%%%%%%%%

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out the ctf-functions
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
      write(tmp_fname,'(A,A)') trim(dir_base) , '/Fig_ctf_c'
      write(6,'(A)') trim(tmp_fname)
      call Fig_gen_ctf_ver0(n_k_p_r_max,n_polar_a_,grid_k_p_r_,min(6
     $     ,n_CTF),ld_CTF,CTF_k_p__,tmp_fname)
      end if !if (flag_fig) then

      allocate(I_ctf_(0:1+n_M-1))
      call cs1_i4(n_M,I_ctf_)
      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Randomly determine which images correspond to which ctf-function.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,n_M-1
         I_ctf_(nm) = mod(nm,n_CTF)
      enddo !do nm=0,n_M-1
      do nm=0,n_M-1
         alpha2d_tru__(nalpha_ctf_ind,nm) = 1.0d0*I_ctf_(nm)
         alpha2d_est__(nalpha_ctf_ind,nm) = 1.0d0*I_ctf_(nm)
      enddo !do nm=0,n_M-1
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then
     
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Allocate space for templates, as well as other temporary variables.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_S_sample_max = n_S_sample_max_set
      if ((n_S_sample_max.le.0) .or. (n_S_sample_max.ge.n_S)) then
         n_S_sample_max = n_S
      end if !if ((n_S_sample_max.le.0) .or. (n_S_sample_max.ge.n_S)) then
      n_M_sample_max = n_M_sample_max_set
      if ((n_M_sample_max.le.0) .or. (n_M_sample_max.ge.n_M)) then
         n_M_sample_max = n_M
      end if !if ((n_M_sample_max.le.0) .or. (n_M_sample_max.ge.n_M)) then
      allocate(S_k_p__(0:1+ld_S*n_azimu_b__sum_(n_k_p_r_max-1)-1))
      n = ld_S*n_azimu_b__sum_(n_k_p_r_max-1)
      call cs1_c16(n,S_k_p__)
      allocate(I_M_sample_(0:1+n_M_sample_max-1))
      n = n_M_sample_max
      call cs1_i4(n,I_M_sample_)
      allocate(I_model_sample_(0:1+n_M_sample_max-1))
      n = n_M_sample_max
      call cs1_i4(n,I_model_sample_)
      allocate(M_residual_sample__(0:1+ld_M*n_M_sample_max-1))
      n = ld_M*n_M_sample_max
      call cs1_c16(n,M_residual_sample__)
      allocate(M_residual_loading_(0:1+n_residual_loading*n_M_sample_max
     $     -1))
      n = n_residual_loading*n_M_sample_max
      call cs1_c16(n,M_residual_loading_)
      allocate(alpha_tru__(0:1+n_alpha*n_M_sample_max-1))
      n = n_alpha*n_M_sample_max
      call cs1_r8(n,alpha_tru__)
      allocate(alpha_est__(0:1+n_alpha*n_M_sample_max-1))
      n = n_alpha*n_M_sample_max
      call cs1_r8(n,alpha_est__)
      allocate(alpha_bkp__(0:1+n_alpha*n_M_sample_max-1))
      n = n_alpha*n_M_sample_max
      call cs1_r8(n,alpha_bkp__)

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  select images to use 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_M_sample = min(n_M,n_M_sample_max)
      n_M_A_sample = 0
      n_M_B_sample = 0
      n_M_C_sample = 0
      n_M_D_sample = 0
      if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_M_sample
     $     ,' images'
      do nM_sample=0,n_M_sample-1
         if (n_M_sample.lt.n_M) then
            I_M_sample_(nM_sample) = min(n_M-1,floor(n_M
     $           *1.0d0*nM_sample/(n_M_sample-1)))
         else
            I_M_sample_(nM_sample) = nM_sample
         end if
         nm = I_M_sample_(nM_sample)
         I_model_sample_(nM_sample) = I_model_(nm)
         if (I_model_(nm).eq.0) then
            n_M_A_sample = n_M_A_sample + 1
         end if                 !if (I_model_(nm).eq.0) then
         if (I_model_(nm).eq.1) then
            n_M_B_sample = n_M_B_sample + 1
         end if                 !if (I_model_(nm).eq.1) then
         if (I_model_(nm).eq.2) then
            n_M_C_sample = n_M_C_sample + 1
         end if                 !if (I_model_(nm).eq.1) then
         if (I_model_(nm).eq.3) then
            n_M_D_sample = n_M_D_sample + 1
         end if                 !if (I_model_(nm).eq.1) then
      enddo                     ! do nM_sample=0,n_M_sample-1
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'n_M_sample: ' ,
     $        n_M_sample,'; n_M_A_sample: ' , n_M_A_sample ,
     $        '; n_M_B_sample: ', n_M_B_sample , '; n_M_C_sample: ' ,
     $        n_M_C_sample , '; n_M_D_sample: ' , n_M_D_sample
      end if                    !if (verbose.gt.0) then
      end if !      if ((nMode.eq.0) .or. (nMode.eq.1)) then

      if (nMode.eq.2) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  select images to use 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_M_sample = min(n_M,n_M_sample_max)
      if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_M_sample
     $     ,' images'
      do nM_sample=0,n_M_sample-1
         if (n_M_sample.lt.n_M) then
            I_M_sample_(nM_sample) = min(n_M-1,floor(n_M
     $           *1.0d0*nM_sample/(n_M_sample-1)))
         else
            I_M_sample_(nM_sample) = nM_sample
         end if
         nm = I_M_sample_(nM_sample)
         I_model_sample_(nM_sample) = I_model_(nm)
      enddo                     ! do nM_sample=0,n_M_sample-1
      end if !if (nMode.eq.2) then

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images (in place)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      timing_tic = omp_get_wtime()
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_(nM_sample)
         nctf = I_ctf_(nm)
         call xx1_c16(ld_M,Y_slice__(0,nm),CTF_k_p__(nctf
     $        *ld_CTF),Y_slice__(0,nm))
      enddo                     ! do nM_sample=0,n_M_sample-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' , timing_tot ,
     $     ' <-- CTF_apply '
      timing_CTF_apply = timing_tot
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha2d__ to alpha1d_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A)') 'copying alpha2d-->alpha'
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_(nM_sample)
         call cp1_r8(n_alpha,alpha2d_tru__(0,nm),alpha_tru__(0
     $        +n_alpha*nM_sample))
         call cp1_r8(n_alpha,alpha2d_est__(0,nm),alpha_est__(0
     $        +n_alpha*nM_sample))
      enddo !do nM_sample=0,n_M_sample-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Set nalpha_M_index 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A)') 'setting nalpha_M_index'
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_(nM_sample)
         alpha_tru__(nalpha_M_index + n_alpha*nM_sample) = 1.0d0
     $        *nM_sample
         alpha_est__(nalpha_M_index + n_alpha*nM_sample) = 1.0d0
     $        *nM_sample
      enddo !do nM_sample=0,n_M_sample-1

      if ((nMode.eq.0) .or. (nMode.eq.1)) then
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor (in place)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_(nM_sample)
         tmp_l2_norm = alpha_tru__(nalpha_l2_norm+n_alpha*nM_sample)
         call af1_c16(ld_M,1.0d0*dcmplx(tmp_l2_norm,0.0d0),1.0d0
     $        *dcmplx(0.0d0,0.0d0),Y_slice__(0,nm),Y_slice__(0,nm))
      enddo                     ! do nM_sample=0,n_M_sample-1
      end if !if ((nMode.eq.0) .or. (nMode.eq.1)) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) , '/Fig_MC_c_'
     $           ,n_k_p_r_max
            call Fig_gen_ctf_ver3(n_k_p_r_max,n_polar_a_,grid_k_p_r_
     $           ,n_M_sample,I_M_sample_,ld_M,Y_slice__,ld_CTF,CTF_k_p__
     $           ,alpha_est__,min(16,n_M_sample),tmp_fname)
            write(6,'(A)') ' Stopping... '
            stop ! examine figures. ;
         end if !if (flag_fig) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  dump timing data. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_function_create.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_function_create,tmp_fname)
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_function_sample.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_function_sample,tmp_fname)
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_Y_create.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_Y_create,tmp_fname)
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_alpha_define.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_alpha_define,tmp_fname)
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_Y_slice_define.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_Y_slice_define,tmp_fname)
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_CTF_apply.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_CTF_apply,tmp_fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  store timing. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(timing_rebuild_tru_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_rebuild_tru_)
      allocate(timing_residual_tru_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_residual_tru_)
      allocate(timing_rebuild_est_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_rebuild_est_)
      allocate(timing_residual_est_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_residual_est_)
      allocate(timing_template_create_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_template_create_)
      allocate(timing_innerproduct_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_innerproduct_)
      allocate(timing_alpha_update_(0:1+n_k_p_r_max + 4 - 1)) ! stores timing. ;
      n = n_k_p_r_max+4
      call cs1_r8(n,timing_alpha_update_)
      
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  check memory allocations. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_ver18_excerpt_checkset.f'

c-----------------------------------------------------------------------
c     6) Run the main loop (i.e., try and recover molecule from images).
c-----------------------------------------------------------------------
c     

      if (flag_displacement) then
         n_delta_x = n_delta_x_set
         n_delta_y = n_delta_x_set
         n_delta_v = ceiling(pi/2.0d0 * n_delta_x*n_delta_y)
         n_pixels_in = n_pixels_in_set
      else
         n_delta_x = 1
         n_delta_y = 1
         n_delta_v = 1
         n_pixels_in = 0.0d0         
      end if !flag_displacement

      timing_all_tic = omp_get_wtime()
      if (flag_skip_marching) then
         n_k_p_r_low = n_k_p_r_max
      else
         n_k_p_r_low = 2
      end if ! if (flag_skip_marching) then
      if (n_k_p_r_low.gt.3) then
         write(6,'(A,I0)') 'Note: n_k_p_r_low starting at: ' ,
     $        n_k_p_r_low
      end if !if (n_k_p_r_low.gt.3) then
      n_k_p_r_bot = n_k_p_r_low
      do n_k_p_r_cur=n_k_p_r_bot,n_k_p_r_max
         if (verbose.gt.0) write(6,'(A,I0)') 'n_k_p_r_cur: ',n_k_p_r_cur

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha2d__ to alpha1d_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'copying alpha2d-->alpha'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha2d_tru__(0,nm),alpha_tru__(0
     $           +n_alpha*nM_sample))
            call cp1_r8(n_alpha,alpha2d_est__(0,nm),alpha_est__(0
     $           +n_alpha*nM_sample))
            if (n_k_p_r_cur.eq.n_k_p_r_bot) then
               call cp1_r8(n_alpha,alpha2d_est__(0,nm),alpha_bkp__(0
     $              +n_alpha*nM_sample))               
            end if !if (n_k_p_r_cur.eq.n_k_p_r_bot) then
         enddo
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Set nalpha_M_index 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'setting nalpha_M_index'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            alpha_tru__(nalpha_M_index + n_alpha*nM_sample) = 1.0d0
     $           *nM_sample
            alpha_est__(nalpha_M_index + n_alpha*nM_sample) = 1.0d0
     $           *nM_sample
            if (n_k_p_r_cur.eq.n_k_p_r_bot) then
               alpha_bkp__(nalpha_M_index + n_alpha*nM_sample) = 1.0d0
     $              *nM_sample
            end if !if (n_k_p_r_cur.eq.n_k_p_r_bot) then
         enddo
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write alphas to disk for postmortem analysis.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = n_alpha
         tmp_d_(1) = n_M_sample
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/alpha_tru__'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,alpha_tru__,tmp_fname)
         tmp_d_(0) = n_alpha
         tmp_d_(1) = n_M_sample
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/alpha_est__'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : '
     $        ,trim(tmp_fname)
         call MDA_write_r8(2,tmp_d_,alpha_est__,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model based on true alpha.
c$$$  Y_tru_ is the model obtained by using the
c$$$  true alpha-parameters for each image.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         lsq_oversample = 2.0d0
         eps_default = 1.0d-6
         lsq_interpolation_order = 7
         if (verbose.gt.0) write(6,'(A)') 'rebuilding Y_tru_'
         timing_tic = omp_get_wtime()
         call rebuild_model_2(
     $        n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $        ,I_M_sample_ !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $        ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $        ,Y_slice__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $        ,n_CTF !integer *4: total number of CTF functions. ;
     $        ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
     $        ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $        ,alpha_tru__ !real *8 array (length at least n_alpha*n_M_sample): true image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_tru__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $        ,n_w_csum_ !integer *4 array (length at least n_k_p_r_cur): cumulative sum of n_w_. also called icstart(nk). ;
     $        ,n_polar_a_ !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $        ,quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
     $        ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $        ,n_w_ !integer *4 array (length at least n_k_p_r_cur): n_w_(nk) is the number of points on the ring of radius grid_k_p_r_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
     $        ,n_k_p_r_low !integer *4: lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_r_max. ;
     $        ,n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_r_max. ;
     $        ,n_Y_lm_csum_ !integer *4 array (size at least n_k_p_r_cur): cumulative sum of n_Y_lm_. also called isph_start(nk). ;
     $        ,n_Y_l_ !integer *4 array (size at least n_k_p_r_cur): order of spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). also called nterms_sph(nk). ;
     $        ,lsq_oversample !real *8: least-squares-solver oversampling parameter. sometimes called oversamp. ;
     $        ,lsq_interpolation_order !integer *4: least-squares-solver interpolation order. sometimes named kord. ;
     $        ,eps_default !real *8: tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
     $        ,Y_tru_ !complex *16 array (size at least n_Y_lm_sum): spherical harmonic coefficients obtained using true-angles to reconstruct molecule. sometimes called modsph_tru. ;
     $        )
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc-timing_tic
         if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' ,
     $        timing_tot ,' <-- rebuild_tru '
         timing_rebuild_tru_(n_k_p_r_cur) = timing_tot
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'rebuild_model_2 (tru):'
     $           ,' total_time ',timing_toc-timing_tic
         end if !if (verbose.gt.0) then
         tmp_d_(0) = n_Y_lm_sum
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/Y_tru_'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_tru_,tmp_fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'Y_tru_ to kspacegrid'
         if (verbose.gt.2) then
            write(6,'(A,F24.16,4(A,I0))')
     $           ' k_p_r_max '
     $           ,k_p_r_max
     $           ,' n_k_p_r_max '
     $           ,n_k_p_r_max
     $           ,' quadrature_type_radial '
     $           ,quadrature_type_radial
     $           ,' quadrature_type_azimu_b '
     $           ,quadrature_type_azimu_b
     $           ,' n_Y_lm_sum '
     $           ,n_Y_lm_sum
            call print_sub_r8(n_x_c_max,grid_x0_c___,' grid_x0_c___: ')
            call print_sub_r8(n_x_c_max,grid_x1_c___,' grid_x1_c___: ')
            call print_sub_r8(n_x_c_max,grid_x2_c___,' grid_x2_c___: ')
            call print_sub_r8(n_azimu_b__sum_sum,grid_k0_c_
     $           ,' grid_k0_c_: ')
            call print_sub_r8(n_azimu_b__sum_sum,grid_k1_c_
     $           ,' grid_k1_c_: ')
            call print_sub_r8(n_azimu_b__sum_sum,grid_k2_c_
     $           ,' grid_k2_c_: ')
         end if !if (verbose.gt.2) then
         call model_to_kspacegrid(
     $        f_k_c_tru_
     $        ,k_p_r_max
     $        ,n_Y_l_
     $        ,quadrature_type_radial
     $        ,n_k_p_r_max
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b
     $        ,n_azimu_b__sum_
     $        ,Y_tru_
     $        ,n_Y_lm_sum
     $        )
         do nA = 0,n_azimu_b__sum_sum-1
            f_k_c_tru_(nA) = f_k_c_tru_(nA)*weight_k_c_(nA)
         enddo
         fft_iflag = -1
         n_x0_c_x1_c_x2_c_sum = n_x_c_max*n_x_c_max*n_x_c_max
         if (verbose.gt.0) write(6,'(A)') 'finufft3d3_f'
         call finufft3d3_f(
     $        n_azimu_b__sum_sum
     $        ,grid_k0_c_
     $        ,grid_k1_c_
     $        ,grid_k2_c_ 
     $        ,f_k_c_tru_
     $        ,fft_iflag
     $        ,eps_default
     $        ,n_x0_c_x1_c_x2_c_sum
     $        ,grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,f_x_c_tru___
     $        ,ier
     $        )
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/f_x_c_tru___'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_tru___,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use Y_tru_ and alpha_tru__ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         if (flag_residual) then
            timing_tic = omp_get_wtime()
            n_Y_lm_csum_cur = n_Y_lm_sum
            if (n_k_p_r_cur.lt.n_k_p_r_max) then
               n_Y_lm_csum_cur = n_Y_lm_csum_(n_k_p_r_cur+1-1)
            end if !if (n_k_p_r_cur.lt.n_k_p_r_max) then
            call test_residual_5(
     $           rseed
     $           ,n_M_sample
     $           ,I_M_sample_
     $           ,ld_M
     $           ,Y_slice__
     $           ,n_CTF
     $           ,ld_CTF
     $           ,CTF_k_p__
     $           ,alpha_tru__
     $           ,n_w_csum_
     $           ,n_polar_a_
     $           ,quadrature_type_azimu_b
     $           ,grid_k_p_r_
     $           ,n_w_
     $           ,n_k_1
     $           ,n_k_p_r_cur
     $           ,n_Y_lm_csum_cur 
     $           ,n_Y_lm_csum_
     $           ,n_Y_l_ 
     $           ,lsq_oversample
     $           ,lsq_interpolation_order
     $           ,Y_tru_
     $           ,M_residual_sample__
     $           ,n_residual_loading
     $           ,n_residual_iteration
     $           ,M_residual_loading_
     $           )
            tmp_d_(0) = n_M_sample
            tmp_d_(1) = 1
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/I_model_sample_tru_',n_k_p_r_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : '
     $           ,trim(tmp_fname)
            call MDA_write_i4(2,tmp_d_,I_model_sample_,tmp_fname)
            tmp_d_(0) = n_residual_loading
            tmp_d_(1) = n_M_sample
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/M_residual_loading_tru_',n_k_p_r_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $           ,trim(tmp_fname)
            call MDA_write_c16(2,tmp_d_,M_residual_loading_,tmp_fname)
            timing_toc = omp_get_wtime()
            timing_tot = timing_toc-timing_tic
            if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' ,
     $           timing_tot ,' <-- residual_tru '
            timing_residual_tru_(n_k_p_r_cur) = timing_tot
         end if !if (flag_residual) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for Y_tru_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_residual .and. flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) ,
     $           '/Fig_residual_tru_',n_k_p_r_cur
            call Fig_R_ver1(n_k_p_r_cur,n_polar_a_,grid_k_p_r_
     $           ,n_M_sample,I_M_sample_,ld_M,Y_slice__
     $           ,M_residual_sample__,alpha_tru__,min(16 ,n_M_sample)
     $           ,tmp_fname)
         end if !if (flag_residual .and. flag_fig) then

         if (flag_perform_estimation.eqv..true.) then ! perform estimation

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
c$$$  update alpha_est__ and Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         timing_tic = omp_get_wtime()
         write(tmp_dname,'(A)') dir_base
         n_gamma_z = max(32,n_w_(n_k_p_r_cur))
         if (flag_tesselation) then
            tesselation_distance_req = 2.0d0 / max(1,n_k_p_r_cur-16)
            n_LT_add = max(1,n_M_sample/10) ! number of templates to add (randomly) after considering local neighborhood in local search. ;
            n_LT_ref = max(1,n_SM_max/2) ! number of image-template pairs to consider when refining local search.
         else
            tesselation_distance_req = 2.0d0
            n_LT_add = 1
            n_LT_ref = 1
         end if !flag_tesselation
         call cp1_r8(n_alpha*n_M_sample_max,alpha_est__,alpha_bkp__)
         n_MS_max = n_M_sample
         call test_alpha_update_8(
     $     verbose ! verbosity level. ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_p_r_max ! integer *4: maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
     $     ,k_p_r_max ! real *8, maximum value of k in k-space polar coordinates. sometimes called rmax. ;
     $     ,n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_r_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $     ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,weight_k_p_r_  ! real *8 array (length at least n_k_p_r_max). Weight associated with grid_k_p_r_(nk) in k-space polar coordinates. sometimes called wtsr(nk). ;
     $     ,weight_k_p_   ! real *8 array (length at least n_w_sum). Contains weight for each point in unrolled polar grid. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_rq determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $      ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $      ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $      ,n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $      ,I_M_sample_ !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $      ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $      ,Y_slice__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $      ,n_CTF !integer *4: total number of CTF functions. ;
     $      ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
     $      ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_CTF): stack of CTF functions. ;
     $      ,alpha_est__ !real *8 array (length at least n_alpha*n_M): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $      ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $      ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
     $      ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $      ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $      ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_p_r_cur' under consideration (which can change from iteration to iteration). ;
     $      ,displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
     $      ,n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
     $      ,n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
     $      ,svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
     $      ,eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
     $      ,flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
     $      ,fpm_howmany_max !integer *4: Maximum number of fftws to call simultaneously within the fftw_plan_many. 
     $      ,n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
     $      ,n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $      ,n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $      ,n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $      ,n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $      ,d_memory_limit !real *8: upper limit on the memory allowed by test_innerproduct_8. ;
     $      ,n_w_csum_ ! cumulative sum of n_w_. also called icstart(nk). ;
     $      ,quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
     $      ,n_w_ ! n_w_(nk) is the number of points on the ring of radius grid_k_p_r_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
     $      ,n_k_p_r_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_r_max. ;
     $      ,n_Y_lm_sum ! total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
     $      ,n_Y_lm_csum_ ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
     $      ,n_Y_l_ ! order of spherical harmonic expansion on sphere at radius grid_k_p_r_(nk). also called nterms_sph(nk). ;
     $      ,lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
     $      ,lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
     $      ,eps_default ! tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
     $      ,Y_est_ ! spherical harmonic coefficients obtained using estimated-angles to reconstruct molecule. sometimes called modsph_est. ;
     $      ,n_azimu_b__sum_ ! total number of points on sphere at radius grid_k_p_r_(nk). also called numonsphere(nk). ;
     $      ,n_S_sample_max    ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
     $      ,flag_fig  ! flag determining whether or not to dump output file. ;
     $      ,tmp_dname ! directory name to dump output (if flag_fig is set). ;
     $      ,timing_rebuild_est_(n_k_p_r_cur) ! time to rebuild Y_est_. ;
     $      ,timing_template_create_(n_k_p_r_cur) ! time to create templates. ;
     $      ,timing_innerproduct_(n_k_p_r_cur) ! time to calculate innerproducts. ;
     $      )
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc-timing_tic
         if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' ,
     $        timing_tot ,' <-- alpha_update '
         timing_alpha_update_(n_k_p_r_cur) = timing_tot

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = n_Y_lm_sum
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/Y_est_'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,Y_est_,tmp_fname)
         if (verbose.gt.0) write(6,'(A)') 'Y_est_ to kspacegrid'
         call model_to_kspacegrid(
     $        f_k_c_est_
     $        ,k_p_r_max
     $        ,n_Y_l_
     $        ,quadrature_type_radial
     $        ,n_k_p_r_max
     $        ,n_polar_a_
     $        ,quadrature_type_azimu_b
     $        ,n_azimu_b__sum_
     $        ,Y_est_
     $        ,n_Y_lm_sum
     $        )
         do nA = 0,n_azimu_b__sum_sum-1
            f_k_c_est_(nA) = f_k_c_est_(nA)*weight_k_c_(nA)
         enddo
         fft_iflag = -1
         n_x0_c_x1_c_x2_c_sum = n_x_c_max*n_x_c_max*n_x_c_max
         if (verbose.gt.0) write(6,'(A)') 'finufft3d3_f'
         call finufft3d3_f(
     $        n_azimu_b__sum_sum
     $        ,grid_k0_c_
     $        ,grid_k1_c_
     $        ,grid_k2_c_ 
     $        ,f_k_c_est_
     $        ,fft_iflag
     $        ,eps_default
     $        ,n_x0_c_x1_c_x2_c_sum
     $        ,grid_x0_c___
     $        ,grid_x1_c___
     $        ,grid_x2_c___
     $        ,f_x_c_est___
     $        ,ier
     $        )
         tmp_d_(0) = n_x_c_max
         tmp_d_(1) = n_x_c_max
         tmp_d_(2) = n_x_c_max
         write(tmp_fname,'(A,A,I0,A)') trim(dir_base) , '/f_x_c_est___'
     $        ,n_k_p_r_cur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(3,tmp_d_,f_x_c_est___,tmp_fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use Y_est_ and alpha_bkp__ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         if (flag_residual) then
            timing_tic = omp_get_wtime()
            n_Y_lm_csum_cur = n_Y_lm_sum
            if (n_k_p_r_cur.lt.n_k_p_r_max) then
               n_Y_lm_csum_cur = n_Y_lm_csum_(n_k_p_r_cur+1-1)
            end if !if (n_k_p_r_cur.lt.n_k_p_r_max) then
            call test_residual_5(
     $           rseed
     $           ,n_M_sample
     $           ,I_M_sample_
     $           ,ld_M
     $           ,Y_slice__
     $           ,n_CTF
     $           ,ld_CTF
     $           ,CTF_k_p__
     $           ,alpha_bkp__
     $           ,n_w_csum_
     $           ,n_polar_a_
     $           ,quadrature_type_azimu_b
     $           ,grid_k_p_r_
     $           ,n_w_
     $           ,n_k_1
     $           ,n_k_p_r_cur
     $           ,n_Y_lm_csum_cur 
     $           ,n_Y_lm_csum_
     $           ,n_Y_l_ 
     $           ,lsq_oversample
     $           ,lsq_interpolation_order
     $           ,Y_est_
     $           ,M_residual_sample__
     $           ,n_residual_loading
     $           ,n_residual_iteration
     $           ,M_residual_loading_
     $           )
            tmp_d_(0) = n_M_sample
            tmp_d_(1) = 1
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/I_model_sample_est_',n_k_p_r_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : '
     $           ,trim(tmp_fname)
            call MDA_write_i4(2,tmp_d_,I_model_sample_,tmp_fname)
            tmp_d_(0) = n_residual_loading
            tmp_d_(1) = n_M_sample
            write(tmp_fname,'(A,A,I0,A)') trim(dir_base) ,
     $           '/M_residual_loading_est_',n_k_p_r_cur,'_.mda'
            if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $           ,trim(tmp_fname)
            call MDA_write_c16(2,tmp_d_,M_residual_loading_,tmp_fname)
            timing_toc = omp_get_wtime()
            timing_tot = timing_toc-timing_tic
            if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' ,
     $           timing_tot ,' <-- residual_est '
            timing_residual_est_(n_k_p_r_cur) = timing_tot
         end if !if (flag_residual) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (flag_residual .and. flag_fig) then
            write(tmp_fname,'(A,A,I0)') trim(dir_base) ,
     $           '/Fig_residual_est_',n_k_p_r_cur
            call Fig_R_ver1(n_k_p_r_cur,n_polar_a_,grid_k_p_r_
     $           ,n_M_sample,I_M_sample_,ld_M,Y_slice__
     $           ,M_residual_sample__,alpha_bkp__,min(16 ,n_M_sample)
     $           ,tmp_fname)
         end if !if (flag_residual .and. flag_fig) then

         end if !if (flag_perform_estimation.eqv..true.) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha to alpha2d
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A)') 'copying alpha-->alpha2d'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha_est__(0+n_alpha*nM_sample)
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
c$$$  dump timing data. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_rebuild_tru_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_rebuild_tru_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_rebuild_est_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_rebuild_est_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_residual_tru_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_residual_tru_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_residual_est_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_residual_est_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_template_create_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_template_create_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_innerproduct_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_innerproduct_,tmp_fname)
         tmp_d_(0) = n_k_p_r_max+1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_alpha_update_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_alpha_update_,tmp_fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  check memory allocations. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_ver18_excerpt_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  end loop
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
      enddo ! do n_k_p_r_cur = n_k_p_r_bot,n_k_p_r_max
      timing_all_toc = omp_get_wtime()
      timing_all_tot = timing_all_toc-timing_all_tic
      if (verbose.gt.0) write(6,'(A,F8.4,A)') ' timing: ' ,
     $     timing_all_tot ,' <-- function_sample '
      timing_frequency_marching = timing_all_tot

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  dump timing data. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         tmp_d_(0) = 1
         tmp_d_(1) = 1
         write(tmp_fname,'(A,A)') trim(dir_base) ,
     $        '/timing_frequency_marching.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : '
     $        ,trim(tmp_fname)
         call MDA_write_c16(2,tmp_d_,timing_frequency_marching
     $        ,tmp_fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  check memory allocations. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_ver18_excerpt_checkset.f'

 10   continue
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver18]'
      end if
      stop
      end
