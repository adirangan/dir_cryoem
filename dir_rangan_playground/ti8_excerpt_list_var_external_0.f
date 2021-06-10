!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> List of variables passed in externally. ;\n
     $     verbose !integer *4: verbosity level. ;
     $     ,flag_memory_estimate !logical: if set to .true. will only estimate total memory requirements. ;
     $     ,flag_time_Zstore !logical: if set to .true. will time sections of Zstore algorithm. ;
     $     ,d_memory_estimate !real *8: estimate of required memory (in bytes). ;
     $     ,rseed !integer *4: random seed (used for any random permutations). ;
     $     ,n_k_p_r_cur !integer *4: current index for current maximum value of n_k. Note that this runs from 1 to n_k_p_max. ;
     $     ,n_polar_a_ !integer *4 array (length at least n_k_p_r_cur): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
     $     ,grid_k_p_r_ !real *8 array (length at least n_k_p_r_cur): values for k on successive shells for k-space polar coordinates. sometimes called grid_p_ or xnodesr. ;
     $     ,weight_k_p_r_ !real *8 array (length at least n_k_p_r_cur): weights for radial quadrature. ;
     $     ,weight_k_p_ !real *8 array (length at least n_w_sum_cur): all weights (unrolled) for quadrature, for use with k_p coordinates. ;
     $     ,half_diameter_x_c !real *8: half diameter of particle support in x-space cartesian coordinates. sometimes called a. ;
     $     ,n_S !integer *4: total number of sampled templates (i.e., length of I_S_sample_). sometimes called ntemplates. ;
     $     ,I_S_sample_ !integer *4 array (length at least n_S): indexing variable used to reference templates. Only templates S_k_p__(I_S_sample_(ns)*ld_S) will be accessed. ;
     $     ,ld_S !integer *4: leading dimension of S_k_p__. Must be at least n_w_sum. ;
     $     ,S_k_p__ !complex *16 array (length at least ld_S*max_i4_f_(I_S_sample_)): stack of templates associated with reconstructed molecule. sometimes called Y_slice__ or cslices or templates. ;
     $     ,tesselation_distance_req !real *8: if adaptively sampling templates, this value determines the maximum distance (on the unit-sphere) allowed between the viewing-angle of any particular template and the estimated-viewing-angle of a given image. If this value is large (e.g., 2.0d0), then a 'global' search will be conducted: every template will be considered for each image. If, on the other hand, this value is small (e.g., <2.0d0), then a 'local' search will be conducted: only a few templates will be considered for each image. When adaptively sampling templates, we might allow this value to shrink with the total number of templates (e.g., to ensure that the number of templates considered per image remains roughly constant). ;
     $     ,n_LT_add ! number of templates to add (randomly) after considering local neighborhood in local search. ;
     $     ,n_LT_ref ! number of image-template pairs to consider when refining local search.
     $     ,S_alpha_S_index_ !real *8 array (length at least n_S): array of S_index associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_S_index_(ns) to equal ns (ranging from 0:n_S-1) which in turn refers to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_polar_a_ !real *8 array (length at least n_S): array of polar_a associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_polar_a_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,S_alpha_azimu_b_ !real *8 array (length at least n_S): array of azimu_b associated with templates indexed within I_S_sample_. ; (used for updating image parameters). Note that we expect S_alpha_azimu_b_(ns) to apply to template S_k_p__(I_S_sample_(ns)*ld_S). ;
     $     ,n_M !integer *4: total number of sampled images (i.e., length of I_M_sample_). sometimes called nimages. ;
     $     ,I_M_sample_ !integer *4 array (length at least n_M): indexing variable used to reference images. Only images M_k_p__(I_M_sample_(nm)*ld_M) will be accessed. ;
     $     ,ld_M !integer *4: leading dimension of M_k_p__. Must be at least n_w_sum. ;
     $     ,M_k_p__ !complex *16 array (length at least ld_M*max_i4_f_(I_M_sample_)): stack of images. sometimes called Y_slice__ associated with reconstructed molecule. 
     $     ,n_CTF !integer *4: total number of CTF functions. ;
     $     ,ld_CTF !integer *4: leading dimension of CTF_k_p__. Must be at least n_w_sum. ;
     $     ,CTF_k_p__ !complex *16 array (length at least ld_CTF*n_ctf): stack of CTF functions. ;
     $     ,alpha_est__ !real *8 array (length at least n_alpha*n_M): ! estimated image parameters associated with sampled image set stored in 1-dimensional array. Note that we expect alpha_est__(n_alpha*nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). ;
     $     ,alpha_update_f !real *8: fraction of 'best' templates to select from when updating image-parameters for any particular image (used when flag_MS_vs_SM.eqv..false.). Also interpreted as the fraction of 'best' images to select from when updating image-parameters when flag_MS_vs_SM.eqv..true. ;
     $     ,flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ; Note that this is required to be .false. for a local search (i.e., when tesselation_distance_req .lt. 2.0d0). ;
     $     ,n_SM_max !integer *4: maximum number of templates-per-image whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_SM_ !integer *4 array (length at least n_M): actual number of templates whose innerproduct-information is stored for a particular image. Note that n_SM_(nm) indicates the number of templates whose information is stored for the image M_k_p__(I_M_sample_(nm)). ;
     $     ,alpha_SM__ !real *8 array (length at least n_alpha*n_SM_max*n_M): actual innerproduct-information stored for various template-image pairs. Note that alpha_SM__((ns + nm*n_SM_max)*n_alpha) (with ns.lt.n_SM_(nm)) stores the innerproduct-information for the image M_k_p__(I_M_sample_(nm)) and the ns-th template whose information is stored for that particular image. ;
     $     ,n_MS_max !integer *4: maximum number of images-per-template whose innerproduct-information will be stored when updating the image-parameters for a particular image. ;
     $     ,n_MS_ !integer *4 array (length at least n_S): actual number of images whose innerproduct-information is stored for a particular template. Note that n_MS_(ns) indicates the number of images whose information is stored for the template S_k_p__(I_S_sample_(ns)). ;
     $     ,alpha_MS__ !real *8 array (length at least n_alpha*n_MS_max*n_S): actual innerproduct-information stored for various image-template pairs. Note that alpha_MS__((nm + ns*n_MS_max)*n_alpha) (with nm.lt.n_MS_(ns)) stores the innerproduct-information for the template S_k_p__(I_S_sample_(ns)) and the nm-th image whose information is stored for that particular template. ;
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
     $     ,n_S_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_S (used for O_S_q__, T_S_q__, Z_S_q__). ;
     $     ,n_S_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_S (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,n_M_0_sub_0in !integer *4: number of requested sub-blocks at level-0 for n_M (used for O_T_R_CTF_M_q__, T_T_R_CTF_M_q__, Z_T_R_CTF_M_q__). ;
     $     ,n_M_1_sub_0in !integer *4: number of requested sub-blocks at level-1 for n_M (used for S_T_T_R_CTF_M_q__, S_Z_T_R_CTF_M_q__). ;
     $     ,C_M_ !complex *16 array (length at least n_M): actual l2-norm (out to frequency n_k_p_r_cur) for each sampled image. Note that we expect C_M_(nm) to apply to image M_k_p__(I_M_sample_(nm)*ld_M). Note also that C_M_(nm) is a raw l2-norm, with no consideration for the CTF. ;
