!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Defines variables assuming no preallocation (i.e., with allocating). ;\n
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      include 'excerpt_define_nalpha.f' !defines length n_alpha and indexing of image-parameters. For example, the polar angle associated with image I_M_sample_(nm) would be stored within alpha_est__(nalpha_polar_a + n_alpha*nm). ;
      integer *4 verbose !integer *4: verbosity level. ;
      logical flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
      logical flag_memory_checkset !temporary: used to check memory map. ;
      real *8 d_memory_estimate !real *8: estimate of required memory (in bytes). ;
      integer *4 rseed !integer *4: random seed (used for any random permutations). ;
      integer *4 n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
      integer *4 n_polar_a_(0:0) !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
      real *8 grid_k_p_r_(0:0) !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
      real *8 weight_k_p_r_(0:0) !real *8 array (length at least n_k_p_r_cur): weights for radial quadrature. ;
      real *8 weight_k_p_(0:0) !real *8 array (length at least n_w_sum_cur): all weights (unrolled) for quadrature, for use with k_p coordinates. ;
      real *8 half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
      integer *4 n_S !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
      integer *4 I_S_sample_(0:0) !integer *4 array (length at least n_S): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
      integer *4 ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_w_sum. ;
      complex *16 S_k_p__(0:0) !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
      real *8 tesselation_distance_req !real *8: if adaptively sampling templates, this value determines the maximum distance (on the unit-sphere) allowed between the viewing-angle of any particular template and the estimated-viewing-angle of a given image. If this value is large (e.g., 2.0d0), then every template will be considered for each image. If, on the other hand, this value is small (e.g., <0.125d0), then only a few templates will be considered for each image. When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
      real *8 S_alpha_S_index_(0:0) !real *8 array (length at least n_S): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to equal ns (ranging from 0:n_S-1) which in turn refers to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8 S_alpha_polar_a_(0:0) !real *8 array (length at least n_S): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      real *8 S_alpha_azimu_b_(0:0) !real *8 array (length at least n_S): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
      integer *4 n_M !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
      integer *4 I_M_sample_(0:0) !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
      integer *4 ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_w_sum. ;
      complex *16 M_k_p__(0:0) !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
      integer *4 n_CTF !integer *4: total number of CTF functions. ;
      integer *4 ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_w_sum. ;
      complex *16 CTF_k_p__(0:0) !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
      real *8 alpha_est__(0:0) !real *8 array (length at least n_alpha*n_M): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
      real *8 alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
      logical flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
      integer *4 n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
      integer *4 n_SM_(0:0) !integer *4 array (length at least n_M): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
      real *8 alpha_SM__(0:0) !real *8 array (length at least n_alpha*n_SM_max*n_M): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
      integer *4 n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
      integer *4 n_MS_(0:0) !integer *4 array (length at least n_S): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
      real *8 alpha_MS__(0:0) !real *8 array (length at least n_alpha*n_MS_max*n_S): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
      real *8 n_pixels_in !real *8: if displacements are considered, this value determines the number of pixels (in each direction) to be considered. The number of pixels is related to the x-space cartesian coordinates by the maximum wavelength 'n_k_p_r_cur' under consideration (which can change from iteration to iteration). ;
      real *8 displacement_max !real *8: if displacements are considered, this value determines the maximum displacement (in x-space cartesian coordinates) allowed when assigning parameters to each image. This parameter can be set to mitigate the 'runaway' phenomenon that can occur as the displacements for each image are updated iteratively. ;
      integer *4 n_delta_v !integer *4: if displacements are considered, this value determines the number of displacements considered (in x-space cartesian coordinates). ;
      real *8 delta_x_(0:0) !real *8 array: (length at least n_delta_v): x-coordinates of displacements. ;
      real *8 delta_y_(0:0) !real *8 array: (length at least n_delta_v): y-coordinates of displacements. ;
      integer *4 n_gamma_z !integer *4: determines the number of in-plane rotations gamma_z to consider for each image-template pair. ; 
      integer *4 svd_calculation_type !integer *4: integer determining how innerproducts are computed across rotations and translations. ;
c$$$      svd_calculation_type == 1 --> encode displacements using svd, then account for in-plane rotations using the fft, then multiply to access displacements. ;
c$$$      svd_calculation_type == 2 --> account for displacements via brute-force. ;
      real *8 eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
      logical flag_RTRT_vs_RTTR !logical: determines whether to compute <R_{+upd}(T_{+upd}(Z)),T_{-est}(R_{-est}(CTF.*M))> (if .true.) or <Z,R_{-upd}(T_{-upd}(T_{-est}(R_{-est}(CTF.*M))))> (if .false.). ;
      integer *4 fpm_howmany_max !integer *4: Maximum number of fftws to call simultaneously within the fftw_plan_many. 
      integer *4 n_omp_sub_0in !integer *4: number of omp sub-blocks (e.g., number of available processors). ;
      integer *4 n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
      integer *4 n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
      integer *4 n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
      integer *4 n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
      complex *16 C_M_(0:0) !real *8 array (length at least n_M): actual l2-norm (out to frequency n_k_p_r_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      function outputs
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer *4 min_i4_f,max_i4_f,sum_i4_f 
      real *8 max_r8_f
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      logical flag_tesselation ! determines whether or not to use tesselation. ;
      logical flag_tesselation_ref ! determines whether or not to use flag_S_use_omp__. ;
      integer *4 n_point_init ! number of points used to initialize S_L_. ;
      integer *4 n_point,npoint ! number of points used to define tesselation tree. ;
      real *8, allocatable :: S_L_(:) ! list of n_point points used to define tesselation tree. ;
      real *8 tradius_min ! minimum radius allowed for tesselation element. ;
      integer *4 nl_max_init,nm_sum_init,ll_sum_init ! initial variables defining tesselation tree structure. ;
      integer *4 nl_max,nm_sum,ll_sum ! variables defining tesselation tree structure. ;
      integer *4 n_T_root_base,nT_root_base ! variables defining base of tesselation tree. ;
      logical parity_input ! temporary: flag used to determine orientation of tesselation element. ;
c$$$      T_ tesselation tree
      integer *4, allocatable :: T_nl_(:) !level
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_tr_(:) !tradius
      integer *4, allocatable :: T_ll_(:) !number of points from S_L_ in T_
      logical   , allocatable :: T_lf_(:) !is leaf
      integer *4, allocatable :: T_c0_(:) !child_0 tesselation_index
      integer *4, allocatable :: T_c1_(:) !child_1 tesselation_index
      integer *4, allocatable :: T_c2_(:) !child_2 tesselation_index
      integer *4, allocatable :: T_c3_(:) !child_3 tesselation_index
      integer *4, allocatable :: T_ls_(:) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4, allocatable :: T_LT_(:) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4, allocatable :: T_root_base_(:) !T_ roots
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      bookkeeping tesselation search across threads
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 tesselation_distance_req_use !actual distance used as lower-bound for size of tesselation elements. ;
      integer *4, allocatable :: n_SM_use_local_(:) ! array storing total number of image-template comparisons made for a particular image. ;
      integer *4 n_SM_use ! total number of image-template comparisons made. ;
      real *8, allocatable :: vp_input_omp__(:) ! array of points on sphere (3-dimensional) (one for each thread). ;
      integer *8 p_vp_input ! memory offset for single thread. ;
      logical, allocatable :: flag_S_use_omp__(:) ! array of flags indicating which templates have been used (one for each thread). ;
      integer *8 p_flag_S_use ! memory offset for single thread. ;
      integer *4 n_LT,n_S_use,n_S_use_sum
      integer *4 n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
      integer *4 n_LT_ref ! number of image-template pairs to consider when refining local search.
      integer *4, allocatable :: n_S_use_sum_(:) ! array of n_S_use (one for each thread). ;
      integer *4 n_add,n_ref,nref,nLT,ns_local
      integer *4 n_pass,n_pass_ref
      logical continue_flag
      integer *4, allocatable :: LT_omp__(:) ! array of indices for nearest neighbor list (one for each thread). ;
      integer *8 p_LT ! memory offest for single thread. ;
      real *8, allocatable :: S_alpha_S_index_local_omp__(:) ! template index. (one for each thread). ;
      integer *8 p_S_alpha_S_index_local ! memory offset for single thread. ;
      real *8, allocatable :: S_alpha_polar_a_local_omp__(:) ! polar angle for templates (3-dimensional) (one for each thread). ;
      integer *8 p_S_alpha_polar_a_local ! memory offset for single thread. ;
      real *8, allocatable :: S_alpha_azimu_b_local_omp__(:) ! azimuthal angle for templates (one for each thread). ;
      integer *8 p_S_alpha_azimu_b_local ! memory offset for single thread. ;
      integer *4, allocatable :: I_S_sample_local_omp__(:) ! array of template indices for use with S_p__ (one for each thread). ;
      integer *8 p_I_S_sample_local ! memory offset for single thread. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      real *8 eps_target ! copy of eps_svd. ;
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_6/Jsvd_excerpt_var_allocation_0on_0.f'
c$$$      include './dir_gen_Jsvd_6/Jsvd_excerpt_var_allocation_off_0.f'
      integer *4 n_svd_max
      parameter (n_svd_max=512)
      real *8 svd_d_max,svd_r_max ! variables used to set ranges for svd-expansion of translation. ;
      real *8 , allocatable :: svd_polyval_U_d_(:) ! array storing polynomial evaluations of svd_U_d_ at the various delta values. ;
      real *8 , allocatable :: svd_polyval_V_r_(:) ! array storing polynomial evaluations of svd_V_r_ at the various grid_k_p_r_(nr) values. ;
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels ! variables used to select svd library. ;
      complex *16, allocatable :: Z_S_svdd_(:) !array storing displacement-operator for the delta-side of the svd-expansion applied to S. ;
      complex *16, allocatable :: Z_M_svdd_(:) !array storing displacement-operator for the delta-side of the svd-expansion applied to M. ;
      complex *16, allocatable :: CTF_R_S__(:) !array to hold innerproducts for each CTF-S pair. ;
      complex *16, allocatable :: CTF_R_S_sub__(:) !array to hold innerproducts for each CTF-S pair. ;
      complex *16, allocatable :: CTF_R_S_sub_omp__(:) !array to hold innerproducts for each CTF-S pair. ;
      complex *16, allocatable :: CTF_R_S_local_omp__(:) !similar to CTF_R_S__, allocated only if we intend to run local search. ;
      integer *8 p_CTF_R_S_local ! memory offset for single thread. ;
      complex *16, allocatable :: O_S_q__(:) !array to hold collection of S in bessel-coefficients. ;
      complex *16, allocatable :: O_S_q_local_omp__(:) !similar to O_S_q__, allocated only if we intend to run local search. ;
      integer *8 p_O_S_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: T_S_q__(:) !array to hold bessel-coefficients of transformed T_{update}(S), where T represents translation by delta_{update} ;
      complex *16, allocatable :: T_S_q_local_omp__(:) !similar to T_S_q__, allocated only if we intend to run local search. ;
      integer *8 p_T_S_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: Z_S_q__(:) !array to hold bessel-coefficients of the transformed Z_{rho}(S), where Z represents the rho-side of the svd-expansion applied to S. ;
      complex *16, allocatable :: Z_S_q_local_omp__(:) !similar to Z_S_q__, allocated only if we intend to run local search. ;
      integer *8 p_Z_S_q_local ! memory offset for single thread. ;
      pointer (p_S_p__,S_p__) !another name for S_k_p__. ;
      complex *16 S_p__(0:0)
      pointer (p_M_p__,M_p__) !another name for M_k_p__. ;
      complex *16 M_p__(0:0)
      pointer (p_CTF_p__,CTF_p__) !another name for CTF_k_p__. ;
      complex *16 CTF_p__(0:0)
      complex *16, allocatable :: S_p_(:) !temporary storage for a single template. ;
      complex *16, allocatable :: S_q_(:) !temporary storage for a single template. ;
      complex *16, allocatable :: S_p_omp__(:) !temporary storage for a few templates across multiple threads. ;
      complex *16, allocatable :: S_q_omp__(:) !temporary storage for a few templates across multiple threads. ;
      complex *16, allocatable :: M_p_(:) !temporary storage for a single image. ;
      complex *16, allocatable :: M_q_(:) !temporary storage for a single image. ;
      complex *16, allocatable :: M_p_omp__(:) !temporary storage for a few images across multiple threads. ;
      complex *16, allocatable :: M_q_omp__(:) !temporary storage for a few images across multiple threads. ;
      complex *16, allocatable :: CTF_p_(:) !temporary storage for a single CTF. ;
      complex *16, allocatable :: CTF_q_(:) !temporary storage for a single CTF. ;
      complex *16, allocatable :: CTF_p_omp__(:) !temporary storage for a few CTFs across multiple threads. ;
      complex *16, allocatable :: CTF_q_omp__(:) !temporary storage for a few CTFs across multiple threads. ;
      integer *4 n_Z_x !length of Z_p_, Z_q_. ;
      complex *16, allocatable :: Z_p_(:) !temporary storage for a single array. ;
      complex *16, allocatable :: Z_q_(:) !temporary storage for a single array. ;
      complex *16, allocatable :: Z_p_omp__(:) !temporary storage for a few arrays across multiple threads. ;
      complex *16, allocatable :: Z_q_omp__(:) !temporary storage for a few arrays across multiple threads. ;
      complex *16, allocatable :: ZZ__(:) !temporary storage for a single array. ;
      complex *16, allocatable :: ZZ_omp__(:) !temporary storage for a few arrays across multiple threads. ;
      complex *16, allocatable :: ZZ_sub_(:) !temporary storage for a single array. ;
      complex *16, allocatable :: ZZ_sub_omp__(:) !temporary storage for a few arrays across multiple threads. ;
      complex *16, allocatable :: ZZ_pos_(:) !temporary storage for a single array. ;
      complex *16, allocatable :: ZZ_pos_omp__(:) !temporary storage for a few arrays across multiple threads. ;
      complex *16, allocatable :: C_trn0_(:) !temporary storage for a single array. ;
      complex *16, allocatable :: C_trn0_omp__(:) !temporary storage for a few arrays across multiple threads. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      precomputed arrays for zgemm
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      complex *16, allocatable :: O_T_R_CTF_M_q__(:) !array to hold bessel-coefficients of the transformed T_est(R_{-gamma_est}(CTF.*M)), where T_est = translation by -delta_est and R = rotation by -gamma_est. ;
      complex *16, allocatable :: O_T_R_CTF_M_q_local_omp__(:) !similar to O_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
      integer *8 p_O_T_R_CTF_M_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: T_T_R_CTF_M_q__(:) !array to hold bessel-coefficients of the transformed T_upd(T_est(R_{-gamma_est}(CTF.*M))), where T_upd = translation by -delta_{update} and T_est = translation by -delta_est and R = rotation by -gamma_est. ;
      complex *16, allocatable :: T_T_R_CTF_M_q_local_omp__(:) !similar to T_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
      integer *8 p_T_T_R_CTF_M_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: Z_T_R_CTF_M_q__(:) !array to hold bessel-coefficients of the transformed Z_{rho}(T_est(R_{-gamma_est}(CTF.*M))), where Z_{rho} = rho-side of svd-expansion applied to M and T_est = translation by -delta_est and R = rotation by -gamma_est. ;
      complex *16, allocatable :: Z_T_R_CTF_M_q_local_omp__(:) !similar to Z_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
      integer *8 p_Z_T_R_CTF_M_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: S_T_T_R_CTF_M_q__(:) !array to hold bessel-coefficients for the product of O_T_R_CTF_M_q__ and T_S_q__ (or T_T_R_CTF_M_q__ and O_S_q__) integrated over k (i.e., nr). ;
      complex *16, allocatable :: S_T_T_R_CTF_M_q_local_omp__(:) !similar to S_T_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
      integer *8 p_S_T_T_R_CTF_M_q_local ! memory offset for single thread. ;
      complex *16, allocatable :: S_Z_T_R_CTF_M_q__(:) !array to hold bessel-coefficients for the product of Z_R_CTF_M_q__ and T_S_q__ (or Z_T_R_CTF_M_q__ and O_S_q__) integrated over k (i.e., nr). ;
      complex *16, allocatable :: S_Z_T_R_CTF_M_q_local_omp__(:) !similar to S_Z_T_R_CTF_M_q__, allocated only if we intend to run local search. ;
      integer *8 p_S_Z_T_R_CTF_M_q_local ! memory offset for single thread. ;
      integer *4 n_C_trn0 !length of corresponding array immediately after allocation. ;
      integer *4 n_C_trn0_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_CTF_R_S_local_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_S_alpha_S_index_local_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_O_S_q_local_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_T_S_q_local_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_Z_S_q_local_omp !length of corresponding array immediately after allocation. ;
      integer *4 n_O_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
      integer *4 n_O_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
      integer *4 n_T_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
      integer *4 n_T_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
      integer *4 n_Z_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
      integer *4 n_Z_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
      integer *4 n_S_T_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
      integer *4 n_S_T_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
      integer *4 n_S_Z_T_R_CTF_M_q__ !length of corresponding array immediately after allocation. ;
      integer *4 n_S_Z_T_R_CTF_M_q_local_omp__ !length of corresponding array immediately after allocation. ;
c$$$      indices
      integer *4 nr,na,ns,nctf,nm,nw
      integer *4 n_r !integer: copy of n_k_p_r_cur
      integer *4, allocatable :: n_w_(:) !array (length at least n_k_p_r_cur) to store number of angles (omega) for each nk. ;
      integer *4, allocatable :: n_w_csum_(:) !array (length at least n_k_p_r_cur) to store cumulative sum of n_w_. ;
      integer *4 n_w_max !integer: largest value in n_w_. ;
      integer *4 n_w_sum !integer storing the sum of n_w_. ;
c$$$      temporary: image-parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: ctf_ind_est_(:)
      real *8, allocatable :: S_index_est_(:)
      real *8, allocatable :: M_index_est_(:)
      real *8, allocatable :: alpha_0in_(:)
      real *8, allocatable :: alpha_0in_omp__(:)
      integer *4, allocatable :: n_MS_omp__(:)
      real *8, allocatable :: alpha_MS_omp__(:)
c$$$      fftw for omp
      integer *8, allocatable :: fftw_plan_frwd__(:)
      integer *8, allocatable :: fftw_plan_back__(:)
      complex *16, allocatable :: fftw_0in__(:)
      complex *16, allocatable :: fftw_out__(:)
c$$$      fftw for local use
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_0in_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_0in_last_,fftw_0in_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_0in_last_(*),fftw_out_last_(*)
c$$$      fftw plan many (fpm) for omp
      integer *4 n_transf !number of entries in transformed vector. this could be either n_delta_v (if svd_calculation_type.eq.2) or n_svd_l (if svd_calculation_type.eq.1). ;
      integer *4 fpm_howmany
      integer *4 n_fpm,nfpm
      integer *8, allocatable :: fpm_frwd_(:)
      integer *8, allocatable ::  fpm_back_(:)
      complex *16, allocatable :: fpm_0in__(:)
      complex *16, allocatable :: fpm_out__(:)
      integer *4 fpm_rank
      integer *4, allocatable :: fpm_n__(:)
      integer *4, allocatable :: fpm_inembed__(:)
      integer *4 fpm_istride
      integer *4 fpm_idist
      integer *4, allocatable :: fpm_onembed__(:)
      integer *4 fpm_ostride
      integer *4 fpm_odist
c$$$      nufft workspace for local use
      integer *4 nufft1df77_lw
      integer *4 nufft1df77_lused
      real *8, allocatable :: nufft1df77_fw_(:)
c$$$      nufft workspace for omp
      integer *4, allocatable :: nufft1df77_lw_(:)
      integer *4, allocatable :: nufft1df77_lused_(:)
      real *8, allocatable :: nufft1df77_fw__(:)
c$$$      array of displacements and rotations to measure
      integer *4 ndv,ngz,ndv_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: gamma_z_(:) !temporary: array of length n_gamma_z holds values for gamma_z_update. ;
c$$$      parameters for blocking S and M
      integer *4 nomp_sub_0in
      integer *4 nx0,nx1,nx2,nx3,nx4
      integer *4 nm0,nm1,nm2,nm3,nm4
      integer *4 ns0,ns1,ns2,ns3,ns4
      integer *4 n_r_tmp,n_w_sum_tmp,n_X_tmp,n_Y_tmp,n_F_tmp,n_S_tmp
      integer *4 nomp_sub_use,nomp_sub,nomp_per,nomp_sum
      integer *4 n_S_9_sub_use,nS_9_sub,nS_9_per,nS_9_sum
      integer *4, allocatable :: n_S_9_per_(:)
      integer *4, allocatable :: n_S_9_sum_(:)
      integer *4 nS_9_per_max,nS_9_per_min
      integer *4 n_M_9_sub_use,nM_9_sub,nM_9_per,nM_9_sum
      integer *4, allocatable :: n_M_9_per_(:)
      integer *4, allocatable :: n_M_9_sum_(:)
      integer *4 nM_9_per_max,nM_9_per_min
      integer *4 n_S_0_sub_use,nS_0_sub,nS_0_per,nS_0_sum
      integer *4, allocatable :: n_S_0_per_(:)
      integer *4, allocatable :: n_S_0_sum_(:)
      integer *4 nS_0_per_max,nS_0_per_min
      integer *4 n_M_0_sub_use,nM_0_sub,nM_0_per,nM_0_sum
      integer *4, allocatable :: n_M_0_per_(:)
      integer *4, allocatable :: n_M_0_sum_(:)
      integer *4 nM_0_per_max,nM_0_per_min
      integer *4 n_S_1_sub_use,nS_1_sub,nS_1_per,nS_1_sum
      integer *4, allocatable :: n_S_1_per_(:)
      integer *4, allocatable :: n_S_1_sum_(:)
      integer *4 nS_1_per_max,nS_1_per_min
      integer *4 n_M_1_sub_use,nM_1_sub,nM_1_per,nM_1_sum
      integer *4, allocatable :: n_M_1_per_(:)
      integer *4, allocatable :: n_M_1_sum_(:)
      integer *4 nM_1_per_max,nM_1_per_min
c$$$      parameters for memory map
      real *8 d_mem
      integer *4 verbose_mem,verbose_timing
c$$$      parameters for timing and testing
c$$$      Note that timing summation is not thread-safe
      real *8 timing_tic,timing_toc,timing_tot,timing_tmp,gnump_tot
      logical flag_fill !logical: indicates whether or not to test fill-only fftw fpm. ;
      parameter(flag_fill=.false.)
      logical flag_test !logical: indicates whether or not to test innerproducts against single image-template calculation. ;
      real *8 timing_total_fftw_plan,gnump_total_fftw_plan
      real *8 timing_total_CTF_R_S,gnump_total_CTF_R_S
      real *8 timing_total_O_S_q,gnump_total_O_S_q
      real *8 timing_total_T_S_q,gnump_total_T_S_q
      real *8 timing_total_Z_S_q,gnump_total_Z_S_q
      real *8 timing_total_O_T_R_CTF_M_q,gnump_total_O_T_R_CTF_M_q
      real *8 timing_total_T_T_R_CTF_M_q,gnump_total_T_T_R_CTF_M_q
      real *8 timing_total_Z_T_R_CTF_M_q,gnump_total_Z_T_R_CTF_M_q
      real *8 timing_total_zgemm,gnump_total_zgemm
      real *8 timing_total_fpm_fill,gnump_total_fpm_fill
      real *8 timing_total_fpm_fftw,gnump_total_fpm_fftw
      real *8 timing_total_transpose,gnump_total_transpose
      real *8 timing_total_Zstore,gnump_total_Zstore
      logical flag_time_Zstore
      real *8 timing_total_Zstore_a,gnump_total_Zstore_a
      real *8 timing_total_Zstore_b,gnump_total_Zstore_b
      real *8 timing_total_Zstore_c,gnump_total_Zstore_c
      real *8 timing_total_Zstore_x,gnump_total_Zstore_x
      real *8 timing_total_Zstore_d,gnump_total_Zstore_d
      real *8 timing_total_Zstore_e,gnump_total_Zstore_e
      real *8 timing_total_merge_alpha_MS,gnump_total_merge_alpha_MS
      character(len=1024) timing_string
c$$$      format_string for printing output
      character(len=64) format_string
      logical flag_MDA
      integer *4 MDA_n_d
      integer *4 MDA_d_(0:64-1)
      character(len=1024) MDA_string
c$$$      pi
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
