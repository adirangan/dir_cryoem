!> Doxygen comment: ;\n
!> This function performs many of the tasks associated with molecular reconstruction: ;\n
!> 1.  First use current stack of image-parameters along with least-squares to find Y_est_. ;\n
!> 2a. Generate templates based on current model (using Y_est_). ;\n
!> 2b. Choose some number of templates to use when angle-fitting. ;\n
!> 3a. Given a memory budget, determine blocking structure for angle-fitting. ;\n
!> 3b. Perform angle-fitting. ;\n
!> 4.  Update stack of image-parameters alpha_est__. ;\n
      subroutine test_alpha_update_6(
     $     verbose ! verbosity level. ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_p_max ! integer *4: maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
     $     ,k_p_max ! real *8, maximum value of k in k-space polar coordinates. sometimes called rmax. ;
     $     ,n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
     $     ,grid_k_p_ !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $      ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $      ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $      ,n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $      ,I_M_sample_ !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $      ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
     $      ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $      ,n_CTF !integer *4: total number of CTF functions. ;
     $      ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
     $      ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $      ,alpha_est__ !real *8 array (length at least n_alpha*n_M_sample): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $      ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $      ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
     $      ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $      ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $      ,n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_cur' under consideration (which can change from iteration to iteration). ;
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
     $      ,n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $      ,n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $      ,n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $      ,n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $      ,d_memory_limit !real *8: upper limit on the memory allowed by test_innerproduct_8. ;
     $      ,n_w_csum_ ! cumulative sum of n_w_. also called icstart(nk). ;
     $      ,quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
     $      ,n_w_ ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
     $      ,n_k_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_max. ;
     $      ,n_Y_lm_sum ! total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
     $      ,n_Y_lm_csum_ ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
     $      ,n_Y_l_ ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
     $      ,lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
     $      ,lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
     $      ,eps_default ! tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
     $      ,Y_est_ ! spherical harmonic coefficients obtained using estimated-angles to reconstruct molecule. sometimes called modsph_est. ;
     $      ,n_azimu_b_polar_a_sum_ ! total number of points on sphere at radius grid_k_p_(nk). also called numonsphere(nk). ;
     $      ,n_S_sample_max    ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
     $      ,flag_fig  ! flag determining whether or not to dump output file. ;
     $      ,dname ! directory name to dump output (if flag_fig is set). ;
     $      ,timing_rebuild_est ! time to rebuild Y_est_. ;
     $      ,timing_template_create ! time to create templates. ;
     $      ,timing_innerproduct ! time to calculate innerproducts. ;
     $      )
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model Y_est_ based on current alpha_est__ ;
c$$$  (i.e., array of estimated image parameters).
c$$$  Then use new Y_est_ to update alpha_est__ (using angle-fitting). ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      implicit none
      include 'omp_lib.h'
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$  Variables used for test_innerproduct_8. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer verbose ! verbosity level. ;
      logical flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
c$$$      parameter(flag_memory_estimate=.false.)
      real *8 d_memory_estimate !real *8: estimate of required memory (in bytes). ;
      integer *4 rseed !integer *4: random seed (used for any random permutations). ;
      integer *4 n_k_p_max ! integer *4: maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      real *8 k_p_max ! real *8, maximum value of k in k-space polar coordinates. sometimes called rmax. ;
      integer *4 n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
      integer *4 n_polar_a_(0:n_k_cur-1) !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
      real *8 grid_k_p_(0:n_k_cur-1) !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
      real *8 half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
      integer *4 n_S_sample !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
      integer *4, allocatable :: I_S_sample_(:) !temporary: integer *4 array (length at least n_S_sample): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
      integer *4 ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
      complex *16 S_k_p__(0:0) !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
      real *8 tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
      integer *4 n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      integer *4 n_LT_ref ! number of image-template pairs to consider when refining local search.
      real *8, allocatable :: S_alpha_S_index_sample_(:) !temporary: real *8 array (length at least n_S_sample): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8, allocatable :: S_alpha_polar_a_sample_(:) !temporary: real *8 array (length at least n_S_sample): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8, allocatable :: S_alpha_azimu_b_sample_(:) !temporary: real *8 array (length at least n_S_sample): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      integer *4 n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
      integer *4 I_M_sample_(0:n_M_sample-1) !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
      integer *4 ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
      complex *16 M_k_p__(0:0) !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
      integer *4 n_CTF !integer *4: total number of CTF functions. ;
      integer *4 ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
      complex *16 CTF_k_p__(0:0) !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
      include 'nalpha_define.f'
      real *8 alpha_est__(0:0) !real *8 array (length at least n_alpha*n_M_sample): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
      real *8 alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
      logical flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
      integer *4 n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
      integer *4, allocatable :: n_SM_(:) !temporary: integer *4 array (length at least n_M_sample): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
      real *8, allocatable :: alpha_SM__(:) !temporary: real *8 array (length at least n_alpha*n_SM_max*n_M_sample): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
      integer *4 n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
      integer *4, allocatable :: n_MS_(:) !temporary: integer *4 array (length at least n_S_sample): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
      real *8, allocatable :: alpha_MS__(:) !temporary: real *8 array (length at least n_alpha*n_MS_max*n_S_sample): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
      real *8 n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_cur' under consideration (which can change from iteration to iteration). ;
      real *8 displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
      integer *4 n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
      integer *4 n_delta_v_use !integer *4: actual n_delta_v used to construct delta_x_ and delta_y_. ;
      real *8, allocatable :: delta_x_(:) !temporary: real *8 array: (length at least n_delta_v_use): x-coordinates of displacements. ;
      real *8, allocatable :: delta_y_(:) !temporary: real *8 array: (length at least n_delta_v_use): y-coordinates of displacements. ;
      integer *4 n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
      integer *4 svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
      real *8 eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
      logical flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
      integer *4 fpm_howmany_max !integer *4: Maximum number of fftws to call simultaneously within the fftw_plan_many. 
      integer *4 n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
      integer *4 n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
      integer *4 n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
      integer *4 n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
      integer *4 n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
      real *8 d_memory_limit !real *8: upper limit on the memory allowed by test_innerproduct_8. ;
      logical flag_continue,flag_proceed
      integer *4 nS_0_sub_tmp,nS_1_sub_tmp,nM_0_sub_tmp,nM_1_sub_tmp
      integer *4 n_S_0_sub_max,n_S_1_sub_max,n_M_0_sub_max,n_M_1_sub_max
      complex *16, allocatable :: C_M_(:) !temporary: complex *16 array (length at least n_M_sample): actual l2-norm (out to frequency n_k_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$  Other variables used to generate templates and perform lsq-solve. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer n_w_csum_(0:n_k_cur-1) ! cumulative sum of n_w_. also called icstart(nk). ;
      integer quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
      integer n_w_(0:n_k_cur-1) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer n_k_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_max. ;
      integer n_Y_lm_sum ! total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
      integer n_Y_lm_csum_(0:n_k_cur-1) ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
      integer n_Y_l_(0:n_k_cur-1) ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
      real *8 lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
      integer lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
      real *8 eps_default ! tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
      complex *16 Y_est_(0:n_Y_lm_sum-1) ! spherical harmonic coefficients obtained using estimated-angles to reconstruct molecule. sometimes called modsph_est. ;
      integer n_azimu_b_polar_a_sum_(0:n_k_cur-1) ! total number of points on sphere at radius grid_k_p_(nk). also called numonsphere(nk). ;
      integer n_S_sample_max    ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
      logical flag_fig  ! flag determining whether or not to dump output file. ;
      character(len=1024) dname ! directory name to dump output (if flag_fig is set). ;
      real *8 timing_rebuild_est ! time to rebuild Y_est_. ;
      real *8 timing_template_create ! time to create templates. ;
      real *8 timing_innerproduct ! time to calculate innerproducts. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$  Other temporary variables. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      logical flag_memory_checkset !temporary: used to check memory. ;
      integer nk !temporary: index for grid_k_p_(nk). ;
      real *8 tesselation_distance_req_use !temporary: actual value of tesselation_distance_req passed to test_innerproduct_fast_wrapper_1. ;
      integer n_azimu_b_polar_a_sum_cur !temporary: n_azimu_b_polar_a_sum_(n_k_cur-1). sometimes called numonsphere_cur. ;
      integer n_k_p_upperbound ! upper bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;      
      real *8, allocatable :: grid_cos_polar_a_tmp_(:) !temporary: values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called xnodesth(nk). ;
      real *8, allocatable :: grid_sin_polar_a_tmp_(:) !temporary: values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called sthetas(nk). ;
      real *8, allocatable :: weight_cos_polar_a_tmp_(:) !temporary: weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called wtsth_(np). ;
      integer, allocatable :: n_azimu_b_tmp_(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_(nk). also called ngridps(np). ;
      real *8, allocatable :: azimu_b_step_tmp_(:) !temporary: grid-spacing for azimu_b. ;
      integer n_azimu_b_polar_a_sum_tmp !temporary: total number of points on sphere at a particular radius grid_k_p_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
      integer nM_sample,nm !temporary: values used to track number of sampled images. ;
      integer npolar_a,nazimu_b,nA !temporary indices for polar_a, azimu_b and azimu_b_polar_a. ;
      real *8 cos_polar_a !temporary: dcos(polar angle). ;
      real *8 sin_polar_a !temporary: dsin(polar angle). ;
      real *8 azimu_b !temporary: azimuthal angle. sometimes called phi. ;
      real *8 azimu_b_step !temporary: grid spacing for azimuthal angle. sometimes called phistep. ;
      integer nS_sample,ns !temporary: values used to track number of sampled templates. ;
      real *8, allocatable :: S_alpha_S_index_all_(:) !temporary: array of S_index associated with templates (used for updating image parameters). ;
      real *8, allocatable :: S_alpha_polar_a_all_(:) !temporary: array of polar_a associated with templates (used for updating image parameters). ;
      real *8, allocatable :: S_alpha_azimu_b_all_(:) !temporary: array of azimu_b associated with templates (used for updating image parameters). ;
      integer *4 n_delta_x,n_delta_y !temporary: used to determine displacements delta_x_ and delta_y_. ;
      real *8 timing_tic,timing_toc !temporary: timing variables. ;
      character(len=1024) tmp_fname,format_string,prefix_string !temporary: strings. ;
      real *8 pi
      pi=4.0d0*datan(1.0d0)
      
      if (verbose.gt.0) write(6,'(A)')
     $     '[entering test_alpha_update_6]'
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose ' , verbose ! verbosity level. ;
         write(6,'(A,I0)') ' rseed ' , rseed !integer *4: random seed (used for any random permutations). ;
         write(6,'(A,I0)') ' n_k_cur ' , n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
         call print_sub_i4(n_k_cur,n_polar_a_,' n_polar_a_: ') !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
         call print_sub_i4(n_k_cur,grid_k_p_,' grid_k_p_: ') !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
         write(6,'(A,F8.4)') ' half_diameter_x_c ' , half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
         write(6,'(A,I0)') ' ld_S ' , ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
         write(6,'(A,F8.4)') ' tesselation_distance_req ' ,
     $        tesselation_distance_req !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
         write(6,'(A,I0)') ' n_LT_add ' , n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
         write(6,'(A,I0)') ' n_LT_ref ' , n_LT_ref ! number of image-template pairs to consider when refining local search.
         write(6,'(A,I0)') ' n_M_sample ' , n_M_sample !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
         call print_sub_i4(n_M_sample,I_M_sample_,' I_M_sample_: ') !integer *4 array (length at least n_M_sample): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
         write(6,'(A,I0)') ' ld_M ' , ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_A. ;
         call print_sub_c16(ld_M*I_M_sample_(n_M_sample-1),M_k_p__,
     $        ' M_k_p__: ') !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
         write(6,'(A,I0)') ' n_CTF ' , n_CTF !integer *4: total number of CTF functions. ;
         write(6,'(A,I0)') ' ld_CTF ' , ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_A. ;
         call print_sub_c16(ld_CTF*n_CTF,CTF_k_p__,' CTF_k_p__: ') !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
         call print_sub_r8(n_alpha*n_M_sample,alpha_est__
     $        ,' alpha_est__: ') !real *8 array (length at least n_alpha*n_M_sample): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
         write(6,'(A,F8.4)') ' alpha_update_f ' , alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
         write(6,'(A,L2)') ' flag_MS_vs_SM ' , flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
         write(6,'(A,I0)') ' n_SM_max ' , n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
         write(6,'(A,I0)') ' n_MS_max ' , n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
         write(6,'(A,F8.4)') ' n_pixels_in ' , n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_cur' under consideration (which can change from iteration to iteration). ;
         write(6,'(A,F8.4)') ' displacement_max ' , displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
         write(6,'(A,I0)') ' n_delta_v ' , n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
         write(6,'(A,I0)') ' n_gamma_z ' , n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
         write(6,'(A,I0)') ' svd_calculation_type ' ,
     $        svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$  svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$  svd_calculation_type == 2 --> account for displacements via brute-force. ;
         write(6,'(A,F8.6)') ' eps_svd ' , eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
         write(6,'(A,l2)') ' flag_RTRT_vs_RTTR ' , flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
         write(6,'(A,I0)') ' fpm_howmany_max ' , fpm_howmany_max !integer *4: Maximum number of fftws to call simultaneously within the fftw_plan_many. 
         write(6,'(A,I0)') ' n_omp_sub_0in ' , n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
         write(6,'(A,I0)') ' n_S_0_sub_0in ' , n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
         write(6,'(A,I0)') ' n_S_1_sub_0in ' , n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
         write(6,'(A,I0)') ' n_M_0_sub_0in ' , n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
         write(6,'(A,I0)') ' n_M_1_sub_0in ' , n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
         call print_sub_i4(n_k_cur,n_w_csum_,' n_w_csum_: ') ! cumulative sum of n_w_. also called icstart(nk). ;
         write(6,'(A,I0)') ' quadrature_type_azimu_b ' ,
     $        quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
         call print_sub_i4(n_k_cur,n_w_,' n_w_: ') ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
         write(6,'(A,I0)') ' n_k_low ' , n_k_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_max. ;
         write(6,'(A,I0)') ' n_k_cur ' , n_k_cur ! current value of nk (k-index in k-space polar coordinates) to use when reconstructing mole
         write(6,'(A,I0)') ' n_Y_lm_sum ' , n_Y_lm_sum ! total number of basis functions used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
         call print_sub_i4(n_k_cur,n_Y_lm_csum_,' n_Y_lm_csum_: ') ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
         call print_sub_i4(n_k_cur,n_Y_l_,' n_Y_l_: ') ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
         write(6,'(A,F8.4)') ' lsq_oversample ' , lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
         write(6,'(A,I0)') ' lsq_interpolation_order ' ,
     $        lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
         write(6,'(A,F8.6)') ' eps_default ' , eps_default ! tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
         call print_sub_i4(n_k_cur,n_azimu_b_polar_a_sum_,
     $        ' n_azimu_b_polar_a_sum_: ') ! total number of points on sphere at radius grid_k_p_(nk). also called numonsphere(nk). ;
         write(6,'(A,I0)') ' n_S_sample_max ' , n_S_sample_max ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
         write(6,'(A,L2)') ' flag_fig ' , flag_fig ! flag determining whether or not to dump output file. ;
         write(6,'(A,A)') ' dname ' , trim(dname) ! directory name to dump output (if flag_fig is set). ;
      end if                    !if (verbose.gt.1) then

      n_k_p_upperbound = max(4*n_k_p_max,2*nint(pi*k_p_max)) ! upper bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;      
      n_azimu_b_polar_a_sum_cur = n_azimu_b_polar_a_sum_(n_k_cur-1)
      allocate(C_M_(0:1+n_M_sample-1))
      call cs1_c16(n_M_sample,C_M_)
      allocate(S_alpha_S_index_all_(0:1+n_azimu_b_polar_a_sum_cur-1))
      call cs1_r8(n_azimu_b_polar_a_sum_cur,S_alpha_S_index_all_)
      allocate(S_alpha_polar_a_all_(0:1+n_azimu_b_polar_a_sum_cur-1))
      call cs1_r8(n_azimu_b_polar_a_sum_cur,S_alpha_polar_a_all_)
      allocate(S_alpha_azimu_b_all_(0:1+n_azimu_b_polar_a_sum_cur-1))
      call cs1_r8(n_azimu_b_polar_a_sum_cur,S_alpha_azimu_b_all_)
      allocate(grid_cos_polar_a_tmp_(0:1+n_k_p_upperbound-1))
      call cs1_r8(n_k_p_upperbound,grid_cos_polar_a_tmp_)
      allocate(grid_sin_polar_a_tmp_(0:1+n_k_p_upperbound-1))
      call cs1_r8(n_k_p_upperbound,grid_sin_polar_a_tmp_)
      allocate(weight_cos_polar_a_tmp_(0:1+n_k_p_upperbound-1))
      call cs1_r8(n_k_p_upperbound,weight_cos_polar_a_tmp_)
      allocate(azimu_b_step_tmp_(0:1+n_k_p_upperbound-1))
      call cs1_r8(n_k_p_upperbound,azimu_b_step_tmp_)
      allocate(n_azimu_b_tmp_(0:1+n_k_p_upperbound-1))
      call cs1_i4(n_k_p_upperbound,n_azimu_b_tmp_)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  First find Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,*) 'rebuilding Y_est_'
      timing_tic = omp_get_wtime()
      call rebuild_model_2(n_M_sample,I_M_sample_,ld_M,M_k_p__,n_ctf
     $     ,ld_CTF,CTF_k_p__,alpha_est__ ,n_w_csum_ ,n_polar_a_
     $     ,quadrature_type_azimu_b,grid_k_p_,n_w_,n_k_low,n_k_cur
     $     ,n_Y_lm_csum_,n_Y_l_ ,lsq_oversample ,lsq_interpolation_order
     $     ,eps_default ,Y_est_)
      timing_toc = omp_get_wtime()
      timing_rebuild_est = timing_toc-timing_tic
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'rebuild_model_2 (est):'
     $        ,' total_time ',timing_toc-timing_tic
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  generate templates based on current model (using Y_est_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_azimu_b_polar_a_sum_cur = n_azimu_b_polar_a_sum_(n_k_cur-1)
      if (verbose.gt.0) write(6,'(A,I0,A)')
     $     'generating n_azimu_b_polar_a_sum_cur = '
     $     ,n_azimu_b_polar_a_sum_cur ,' templates'
      timing_tic = omp_get_wtime()
      call cl1_c16(ld_S*n_azimu_b_polar_a_sum_cur,S_k_p__)
      call template_gen(Y_est_,n_Y_lm_sum,n_Y_lm_csum_, n_Y_l_
     $     ,n_k_cur,n_w_,ld_S,n_w_csum_, n_polar_a_
     $     ,quadrature_type_azimu_b,n_k_cur
     $     ,n_azimu_b_polar_a_sum_cur ,S_k_p__)
      timing_toc = omp_get_wtime()
      timing_template_create = timing_toc-timing_tic
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'template_gen:' ,' total_time '
     $        ,timing_toc-timing_tic
      end if

      call getspheregrid(n_polar_a_(n_k_cur-1)
     $     ,quadrature_type_azimu_b,grid_cos_polar_a_tmp_,
     $     grid_sin_polar_a_tmp_,weight_cos_polar_a_tmp_
     $     ,n_azimu_b_tmp_,azimu_b_step_tmp_
     $     ,n_azimu_b_polar_a_sum_tmp)
      if (n_azimu_b_polar_a_sum_tmp.ne.n_azimu_b_polar_a_sum_cur)
     $     then
         write(6,'(A,I0,A,I0)') 'Warning, n_azimu_b_polar_a_sum '
     $        ,n_azimu_b_polar_a_sum_tmp,'.ne.'
     $        ,n_azimu_b_polar_a_sum_(n_k_cur-1)
      end if
      nA = 0
      do npolar_a = 0,n_polar_a_(n_k_cur-1)-1
         cos_polar_a = grid_cos_polar_a_tmp_(npolar_a)
         sin_polar_a = grid_sin_polar_a_tmp_(npolar_a)
         azimu_b_step = azimu_b_step_tmp_(npolar_a)
         do nazimu_b = 0,n_azimu_b_tmp_(npolar_a)-1
            azimu_b = nazimu_b*azimu_b_step
            S_alpha_S_index_all_(nA) = nA
            S_alpha_polar_a_all_(nA) = dacos(cos_polar_a)
            S_alpha_azimu_b_all_(nA) = azimu_b
            nA = nA+1
         enddo
      enddo

      n_azimu_b_polar_a_sum_tmp = nA
      if (verbose.gt.1) write(6,*) '  n_azimu_b_polar_a_sum is '
     $     ,n_azimu_b_polar_a_sum_tmp
      if (n_azimu_b_polar_a_sum_tmp.ne.n_azimu_b_polar_a_sum_cur)
     $     then
         write(6,'(A,I0,A,I0)') 'Warning, n_azimu_b_polar_a_sum '
     $        ,n_azimu_b_polar_a_sum_tmp,'.ne.'
     $        ,n_azimu_b_polar_a_sum_(n_k_cur-1)
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Choose number of templates to use 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_S_sample = min(n_azimu_b_polar_a_sum_cur,n_S_sample_max)
      if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_S_sample
     $     ,' templates'
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Allocate memory appropriately
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_S_sample_(0:1+n_S_sample-1))
      call cs1_i4(n_S_sample,I_S_sample_)
      allocate(S_alpha_S_index_sample_(0:1+n_S_sample-1))
      call cs1_r8(n_S_sample,S_alpha_S_index_sample_)
      allocate(S_alpha_polar_a_sample_(0:1+n_S_sample-1))
      call cs1_r8(n_S_sample,S_alpha_polar_a_sample_)
      allocate(S_alpha_azimu_b_sample_(0:1+n_S_sample-1))
      call cs1_r8(n_S_sample,S_alpha_azimu_b_sample_)
      n_delta_x = ceiling(dsqrt(n_delta_v*2.0d0/pi))
      n_delta_y = ceiling(dsqrt(n_delta_v*2.0d0/pi))
      n_delta_v_use = 1 + ceiling(1.5d0 * pi/2.0d0 * n_delta_x
     $     *n_delta_y)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_delta_v_use ' , n_delta_v_use
      end if !if (verbose.gt.0) then
      allocate(delta_x_(0:1+n_delta_v_use-1))
      call cs1_r8(n_delta_v_use,delta_x_)
      allocate(delta_y_(0:1+n_delta_v_use-1))
      call cs1_r8(n_delta_v_use,delta_y_)
      if (flag_MS_vs_SM.eqv..false.) then
         allocate(n_SM_(0:1+n_M_sample-1))
         call cs1_i4(n_M_sample,n_SM_)
         allocate(alpha_SM__(0:1+n_alpha*n_SM_max*n_M_sample-1))
         call cs1_r8(n_alpha*n_SM_max*n_M_sample,alpha_SM__)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         allocate(n_MS_(0:1+n_S_sample-1))
         call cs1_i4(n_S_sample,n_MS_)
         allocate(alpha_MS__(0:1+n_alpha*n_MS_max*n_S_sample-1))
         call cs1_r8(n_alpha*n_MS_max*n_S_sample,alpha_MS__)
      end if !if (flag_MS_vs_SM.eqv..true.) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_alpha_update_6_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Extracting sample templates
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nS_sample=0,n_S_sample-1
         if (n_S_sample.lt.n_azimu_b_polar_a_sum_cur) then
            I_S_sample_(nS_sample) = min(n_azimu_b_polar_a_sum_cur-1
     $           ,floor(n_azimu_b_polar_a_sum_cur*1.0d0*nS_sample
     $           /(n_S_sample -1)))
         else
            I_S_sample_(nS_sample) = nS_sample
         end if
         ns = I_S_sample_(nS_sample)
         S_alpha_S_index_sample_(nS_sample) =
     $        S_alpha_S_index_all_(ns)
         S_alpha_polar_a_sample_(nS_sample) =
     $        S_alpha_polar_a_all_(ns)
         S_alpha_azimu_b_sample_(nS_sample) =
     $        S_alpha_azimu_b_all_(ns)
      enddo

      include 'test_alpha_update_6_checkset.f'

      call get_delta_1(n_pixels_in,n_k_cur,half_diameter_x_c
     $     ,n_delta_x,n_delta_y,n_delta_v,delta_x_,delta_y_)
      if (verbose.gt.0) then
         call print_sub_r8(n_delta_v,delta_x_,'delta_x_')
         call print_sub_r8(n_delta_v,delta_y_,'delta_y_')
      end if !if (verbose.gt.0) then

      include 'test_alpha_update_6_checkset.f'

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A,I0)') 'get (n_M = '
     $     ,n_M_sample,')-x-(n_S_sample = ',n_S_sample, ') = (',
     $     n_M_sample *n_S_sample,') innerproducts of length ' ,
     $     n_w_csum_(n_k_cur -1) + n_w_(n_k_cur-1)
      timing_tic = omp_get_wtime()
      tesselation_distance_req_use = tesselation_distance_req
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
         include 'test_alpha_update_6_checkset.f'
      call test_innerproduct_8(
     $     verbose-2 !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
     $     ,grid_k_p_ !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S_sample !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S_sample): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req_use !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_sample_ !real *8 array (length at least n_S_sample): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_sample_ !real *8 array (length at least n_S_sample): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_sample_ !real *8 array (length at least n_S_sample): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
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
     $     ,n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
     $     ,nS_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,nS_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,nM_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,nM_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M_sample): actual l2-norm (out to frequency n_k_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
        if (d_memory_estimate.le.d_memory_limit) then
           if (verbose.gt.-1) then
              write(6,'(A,4(A,I0))') ' Memory limit not exceeded: '
     $             , ' nS_0_sub_tmp ' ,nS_0_sub_tmp , ' nS_1_sub_tmp ' ,
     $             nS_1_sub_tmp ,' nM_0_sub_tmp ' , nM_0_sub_tmp ,
     $             ' nM_1_sub_tmp' , nM_1_sub_tmp 
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
      if (verbose.gt.1) write(6,'(A)') ' call test_innerproduct_8: '
      include 'test_alpha_update_6_checkset.f'
      call test_innerproduct_8(
     $     verbose !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_cur): number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
     $     ,grid_k_p_ !real *8 array (length at least n_k_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S_sample !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S_sample): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_A. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req_use !real *8: !determines whether or not to adaptively sample templates. if tesselation_distance_req.ge.2.0d0, then all templates will be compared to all images. However, if tesselation_distance_req.lt.2.0d0, then only a few templates will be considered for each image. Roughly speaking, the value of tesselation_distance_req determines the neighborhood of viewing angles around each image which will be searched (in terms of distance on the sphere). When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_sample_ !real *8 array (length at least n_S_sample): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_sample_ !real *8 array (length at least n_S_sample): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_sample_ !real *8 array (length at least n_S_sample): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
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
     $     ,n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
     $     ,nS_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_S_sample (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,nS_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_S_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,nM_0_sub_tmp !integer *4: number of requested sub-blocks at level-0 for n_M_sample (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,nM_1_sub_tmp !integer *4: number of requested sub-blocks at level-1 for n_M_sample (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M_sample): actual l2-norm (out to frequency n_k_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
     $     )
      timing_toc = omp_get_wtime()
      timing_innerproduct = timing_toc-timing_tic
      if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
         write(6,'(A,A,F8.3)') 'test_innerproduct_8:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if ((verbose.gt.0) .and. (flag_memory_estimate.eqv..false.)) then
      end if !if (flag_proceed.eqv..true.) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_alpha_update_6_checkset.f'

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
            call test_alpha_update_SM_3(verbose,rseed,n_SM_max,n_SM_(nm)
     $           ,alpha_SM__(n_alpha*n_SM_max*nm),flag_RTRT_vs_RTTR
     $           ,alpha_est__(n_alpha*nm),alpha_update_f
     $           ,alpha_est__(n_alpha*nm))
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
         call test_alpha_update_MS_3(verbose,rseed,n_MS_max,n_MS_
     $        ,n_S_sample,n_M_sample,alpha_MS__,flag_RTRT_vs_RTTR
     $        ,alpha_est__,alpha_est__)
      end if !if (flag_MS_vs_SM.eqv..true.) then 

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(tmp_fname,'(A,A,I0)') trim(dname) , '/Fig_MTS_c_' ,
     $        n_k_cur
         call Fig_gen_ver5(n_k_cur,n_polar_a_,grid_k_p_,n_S_sample
     $        ,I_S_sample_,ld_S,S_k_p__,n_M_sample,I_M_sample_,ld_M
     $        ,M_k_p__ ,ld_CTF,CTF_k_p__,alpha_est__ ,min(16
     $        ,n_M_sample) ,tmp_fname)
      end if !if (flag_fig) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calling checkset. ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      include 'test_alpha_update_6_checkset.f'

      if (flag_MS_vs_SM.eqv..false.) then
         deallocate(n_SM_)
         deallocate(alpha_SM__)
      end if !if (flag_MS_vs_SM.eqv..false.) then
      if (flag_MS_vs_SM.eqv..true.) then
         deallocate(n_MS_)
         deallocate(alpha_MS__)
      end if !if (flag_MS_vs_SM.eqv..true.) then
      deallocate(delta_x_)
      deallocate(delta_y_)
      deallocate(S_alpha_S_index_all_)
      deallocate(S_alpha_polar_a_all_)
      deallocate(S_alpha_azimu_b_all_)
      deallocate(grid_cos_polar_a_tmp_)
      deallocate(grid_sin_polar_a_tmp_)
      deallocate(weight_cos_polar_a_tmp_)
      deallocate(azimu_b_step_tmp_)
      deallocate(n_azimu_b_tmp_)
      deallocate(C_M_)
      deallocate(I_S_sample_)
      deallocate(S_alpha_S_index_sample_)
      deallocate(S_alpha_polar_a_sample_)
      deallocate(S_alpha_azimu_b_sample_)

      if (verbose.gt.0) write(6,'(A)')
     $     '[finished test_alpha_update_6]'
      
      end


