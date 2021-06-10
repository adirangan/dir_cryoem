      subroutine test_alpha_update_SM_1(verbose,n_omp_sub
     $     ,fpm_howmany_max,half_diameter_x_c_tmp,n_k_p_max,n_M_sample
     $     ,ld_M,M_sample_,n_w_sum_,n_polar_a_,quadrature_type_azimu_b
     $     ,grid_k_p_,n_w_,n_k_low,n_k_cur,n_Y_lm_sum_max,n_Y_lm_sum_
     $     ,n_Y_l_,lsq_oversample,lsq_interpolation_order,eps_default
     $     ,n_ctf,ld_CTF ,CTF_k_p_,n_alpha,alpha1d_est_,alpha_update_f
     $     ,Y_est_,n_azimu_b_polar_a_sum_,ld_S ,n_S_sample_max
     $     ,n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type ,eps_svd
     $     ,n_pixels_in ,displacement_max,flag_tesselation
     $     ,tesselation_distance_req, flag_fig,dname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model Y_est_ based on current alpha1d_est_ ;
c$$$  (i.e., array of estimated image parameters).
c$$$  Then use new Y_est_ to update alpha1d_est_ (using angle-fitting). ;
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      implicit none
      include 'omp_lib.h'
      integer verbose ! verbosity level. ;
      integer n_omp_sub ! number of batches to use when finding innerproducts
      integer fpm_howmany_max ! Maximum number of fftws to call simultaneously within the fftw_plan_many. 
      real *8 half_diameter_x_c_tmp ! half-diameter in x-space cartesian grid. sometimes called max_x_c. ;      
      integer n_k_p_max ! maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      integer n_M_sample ! total number of sampled images (passed in as input). ;
      integer ld_M ! total number of entries in each image (in k-space polar coordinates). i.e., leading-dimension of M array. sometimes called n_image_size or nimagesize. ;
      complex *16 M_sample_(0:ld_M*n_M_sample-1) !temporary: sampled images. These are modified from Y_slice_sample_ by multiplying by the CTF and/or modifying the l2-norm (depending on which flags are set above). ;
      integer n_w_sum_(0:n_k_p_max-1) ! cumulative sum of n_w_. also called icstart(nk). ;
      integer n_polar_a_(0:n_k_p_max-1) ! number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
      integer quadrature_type_azimu_b ! quadrature type in the azimuthal direction. sometimes named itypep. ;
      real *8 grid_k_p_(0:n_k_p_max-1) ! values for k on successive shells for k-space polar coordinates. sometimes called xnodesr(nk). ;      
      integer n_w_(0:n_k_p_max-1) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer n_k_low ! lowest value of nk (k-index in k-space polar coordinates) to use when building models. Note that this runs from 1 to n_k_p_max. ;
      integer n_k_cur ! current value of nk (k-index in k-space polar coordinates) to use when reconstructing mole
      integer n_Y_lm_sum_max ! total number of points used in spherical harmonic expansion (in k-space polar coordinates). sometimes called nsphstore. ;
      integer n_Y_lm_sum_(0:n_k_p_max-1) ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
      integer n_Y_l_(0:n_k_p_max-1) ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
      real *8 lsq_oversample ! least-squares-solver oversampling parameter. sometimes called oversamp. ;
      integer lsq_interpolation_order ! least-squares-solver interpolation order. sometimes named kord. ;
      real *8 eps_default ! tolerance epsilon, typically 1e-6. used in many places (e.g., least-squares-solver). ;
      real *8 eps_svd ! svd tolerance epsilon, typically 1. ;
      integer *4 n_ctf ! stores the number of CTF functions (i.e., 3)
      integer ld_CTF ! total number of entries in each CTF (in k-space polar coordinates). i.e., leading-dimension of CTF array. Typically expected to be at least ld_M. ;
      complex *16 CTF_k_p_(0:n_ctf*ld_M-1) ! stores the n_ctf different CTF functions in k-space polar coordinates. ;
      integer *4 n_alpha ! number of image parameters (i.e., 'alphas') per image. defined in 'nalpha_define.f'. ;
      real *8 alpha1d_est_(0:n_alpha*n_M_sample-1) ! estimated image parameters associated with sampled image set (e.g., stored in Y_slice_sample_ or M_sample_) stored in 1-dimensional array. ;
      real *8 alpha_update_f ! fraction of 'best' templates to select from when updating image-parameters for any particular image.
      complex *16 Y_est_(0:n_Y_lm_sum_max-1) ! spherical harmonic coefficients obtained using estimated-angles to reconstruct molecule. sometimes called modsph_est. ;
      integer n_azimu_b_polar_a_sum_(0:n_k_p_max-1) ! total number of points on sphere at radius grid_k_p_(nk). also called numonsphere(nk). ;
      integer ld_S ! total number of entries in each template (in k-space polar coordinates). sometimes called ntemplate_size or ntemplatesize. ;
      integer n_S_sample_max    ! maximum number of templates to consider (if quadrature_type_azimu_b==1, these will be distributed uniformly across the sphere). ;
      integer n_delta_x,n_delta_y ! values used for n_delta_x and n_delta_y, respectively, when calculating innerproducts involving displacements. ;
      integer n_gamma_z !temporary: value used for n_gamma_z when calculating innerproducts. ;
      integer svd_calculation_type ! flag determining how innerproducts are computed across rotations and translations. ;
      real *8 n_pixels_in ! value used for n_pixels_in when calculating innerproducts involving displacements. ;
      real *8 displacement_max ! value used for displacement_max when calculating innerproducts involving displacements. ;
      logical flag_tesselation ! flag determining whether or not to use adaptive sampling of templates when calculating innerproducts for each image. ;
      real *8 tesselation_distance_req ! if adaptively sampling templates, this value determines the maximum distance (on the unit-sphere) allowed between the viewing-angle of any particular template and the estimated-viewing-angle of a given image. ;
      logical flag_fig  ! flag determining whether or not to dump output file. ;
      character(len=1024) dname ! directory name to dump output (if flag_fig is set). ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
c$$$  Temporary variables. ;
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ;
      integer nk !temporary: index for grid_k_p_(nk). ;
      real *8 tesselation_distance_req_use !temporary: actual value of tesselation_distance_req passed to test_innerproduct_fast_wrapper_0. ;
      integer *4 n_SM_use !temporary: total number of image-template pairs considered when calculating innerproducts. Might be smaller than the total number when calculating innerproducts adaptively. ;
      integer n_azimu_b_polar_a_sum_cur !temporary: n_azimu_b_polar_a_sum_(n_k_cur-1). sometimes called numonsphere_cur. ;
      integer n_k_p_lowerbound ! lower bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;      
      real *8, allocatable :: grid_cos_polar_a_tmp_(:) !temporary: values for cos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called xnodesth(nk). ;
      real *8, allocatable :: grid_sin_polar_a_tmp_(:) !temporary: values for sin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called sthetas(nk). ;
      real *8, allocatable :: weight_cos_polar_a_tmp_(:) !temporary: weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_(nk). sometimes called wtsth_(np). ;
      integer, allocatable :: n_azimu_b_tmp_(:) !temporary: number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_(nk). also called ngridps(np). ;
      real *8, allocatable :: azimu_b_step_tmp_(:) !temporary: grid-spacing for azimu_b. ;
      integer n_azimu_b_polar_a_sum_tmp !temporary: total number of points on sphere at a particular radius grid_k_p_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
      integer nM_sample,nm !temporary: values used to track number of sampled images. ;
      integer npolar_a,nazimu_b,nA !temporary indices for polar_a, azimu_b and azimu_b_polar_a. ;
      real *8 cos_polar_a !temporary: cos(polar angle). ;
      real *8 sin_polar_a !temporary: sin(polar angle). ;
      real *8 azimu_b !temporary: azimuthal angle. sometimes called phi. ;
      real *8 azimu_b_step !temporary: grid spacing for azimuthal angle. sometimes called phistep. ;
      integer n_S_sample,nS_sample,ns !temporary: values used to track number of sampled templates. ;
      integer *4, allocatable :: I_S_sample_(:) !temporary: indexing variable used to reference templates. ;
      complex *16, allocatable :: S__(:,:) !temporary: stack of templates associated with reconstructed molecule. sometimes called templates. ;
      complex *16, allocatable :: S_sample_(:) !temporary: sampled templates. ;
      complex *16, allocatable :: C_M_(:) !temporary: estimated l2-norm associated each sampled image. ;
      real *8, allocatable :: S_alpha_polar_a_all_(:) !temporary: array of polar_a associated with templates (used for updating image parameters). ;
      real *8, allocatable :: S_alpha_azimu_b_all_(:) !temporary: array of azimu_b associated with templates (used for updating image parameters). ;
      real *8, allocatable :: S_alpha_polar_a_sample_(:) !temporary: array of polar_a associated with sampled templates (used for updating image parameters). ;
      real *8, allocatable :: S_alpha_azimu_b_sample_(:) !temporary: array of azimu_b associated with sampled templates (used for updating image parameters). ;
      real *8 timing_tic,timing_toc !temporary: timing variables. ;
      character(len=1024) tmp_fname,format_string !temporary: strings. ;
      
      if (verbose.gt.0) write(6,'(A)')
     $     '[entering test_alpha_update_SM_1]'
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_omp_sub: ',n_omp_sub
         write(6,'(A,F8.4)') 'half_diameter_x_c_tmp: '
     $        ,half_diameter_x_c_tmp
         write(6,'(A,I0)') 'n_k_p_max: ',n_k_p_max
         write(6,'(A,I0)') 'n_M_sample: ',n_M_sample
         write(6,'(A,I0)') 'ld_M: ',ld_M
         call write_sub_c16(n_M_sample*ld_M,M_sample_,11,'M_sample_: ')
         call write_all_i4(n_k_cur,n_w_sum_,10,'n_w_sum_: ')
         call write_all_i4(n_k_cur,n_polar_a_,12,'n_polar_a_: ')
         write(6,'(A,I0)') 'quadrature_type_azimu_b: '
     $        ,quadrature_type_azimu_b
         call write_all_r8(n_k_cur,grid_k_p_,11,'grid_k_p_: ')
         call write_all_i4(n_k_cur,n_w_,6,'n_w_: ')
         write(6,'(A,I0)') 'n_k_low: ',n_k_low
         write(6,'(A,I0)') 'n_k_cur: ',n_k_cur
         write(6,'(A,I0)') 'n_Y_lm_sum_max: ',n_Y_lm_sum_max
         call write_all_i4(n_k_cur,n_Y_lm_sum_,13,'n_Y_lm_sum_: ')
         call write_all_i4(n_k_cur,n_Y_l_,8,'n_Y_l_: ')
         write(6,'(A,F8.4)') 'lsq_oversample: ',lsq_oversample
         write(6,'(A,I0)') 'lsq_interpolation_order: '
     $        ,lsq_interpolation_order
         write(6,'(A,F16.14)') 'eps_default: ',eps_default
         write(6,'(A,I0)') 'n_ctf: ',n_ctf
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         call write_sub_c16(n_ctf*ld_CTF,CTF_k_p_,10,'CTF_k_p_: ')
         write(6,'(A,I0)') 'n_alpha: ',n_alpha
         call write_all_r8(n_alpha,alpha1d_est_(0),17
     $        ,'alpha1d_est_(0): ')
         call write_all_r8(n_alpha,alpha1d_est_((n_M_sample-1)*n_alpha)
     $        ,17,'alpha1d_est_(X): ')
         write(6,'(A,F8.4)') 'alpha_update_f: ',alpha_update_f
         write(6,'(A)') 'Y_est_'
         call write_all_i4(n_k_cur,n_azimu_b_polar_a_sum_,24
     $        ,'n_azimu_b_polar_a_sum_: ')         
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(6,'(A,I0)') 'n_S_sample_max: ',n_S_sample_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
         write(6,'(A,F8.4)') 'eps_svd: ',eps_svd
         write(6,'(A,F8.4)') 'n_pixels_in: ',n_pixels_in
         write(6,'(A,F8.4)') 'displacement_max: ',displacement_max
         write(6,'(A,L1)') ' flag_tesselation: ',flag_tesselation
         write(6,'(A,F8.4)') 'tesselation_distance_req: '
     $        ,tesselation_distance_req
         write(6,'(A,L1)') ' flag_fig: ',flag_fig
         write(6,'(A,A)') ' dname: ' , trim(dname)
      end if !if (verbose.gt.1) then

      n_k_p_lowerbound=4*n_k_p_max ! lower bound on the number of k-values (in k-space polar coordintes). sometimes named ntmax. ;      
      n_azimu_b_polar_a_sum_cur = n_azimu_b_polar_a_sum_(n_k_cur-1)
      allocate(C_M_(0:n_M_sample-1))
      allocate(S__(0:ld_S-1,0:n_azimu_b_polar_a_sum_cur-1))
      allocate(S_alpha_polar_a_all_(0:n_azimu_b_polar_a_sum_cur-1))
      allocate(S_alpha_azimu_b_all_(0:n_azimu_b_polar_a_sum_cur-1))
      allocate(grid_cos_polar_a_tmp_(0:n_k_p_lowerbound-1))
      allocate(grid_sin_polar_a_tmp_(0:n_k_p_lowerbound-1))
      allocate(weight_cos_polar_a_tmp_(0:n_k_p_lowerbound-1))
      allocate(azimu_b_step_tmp_(0:n_k_p_lowerbound-1))
      allocate(n_azimu_b_tmp_(0:n_k_p_lowerbound-1))

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  First find Y_est_
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,*) 'rebuilding Y_est_'
      timing_tic = second()
      call rebuild_model_0(n_M_sample,ld_M,M_sample_,n_ctf
     $     ,ld_CTF,CTF_k_p_,n_alpha,alpha1d_est_ ,n_w_sum_ ,n_polar_a_
     $     ,quadrature_type_azimu_b,grid_k_p_,n_w_,n_k_low,n_k_cur
     $     ,n_Y_lm_sum_,n_Y_l_ ,lsq_oversample
     $     ,lsq_interpolation_order,eps_default ,Y_est_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'rebuild_model_0 (est):'
     $        ,' total_time ',timing_toc-timing_tic
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  generate templates based on current model (using Y_est_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      n_azimu_b_polar_a_sum_cur = n_azimu_b_polar_a_sum_(n_k_cur-1)
      if (verbose.gt.0) write(6,'(A,I0,A)')
     $     'generating n_azimu_b_polar_a_sum_cur = '
     $     ,n_azimu_b_polar_a_sum_cur ,' templates'
      timing_tic = second()
      call template_gen(Y_est_,n_Y_lm_sum_max,n_Y_lm_sum_, n_Y_l_
     $     ,n_k_p_max,n_w_,ld_S,n_w_sum_, n_polar_a_
     $     ,quadrature_type_azimu_b,n_k_cur
     $     ,n_azimu_b_polar_a_sum_cur ,S__)
      timing_toc = second()
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
      allocate(I_S_sample_(0:n_S_sample-1))
      allocate(S_sample_(0:ld_S*n_S_sample-1))
      allocate(S_alpha_polar_a_sample_(0:n_S_sample-1))
      allocate(S_alpha_azimu_b_sample_(0:n_S_sample-1))
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
         call cp1_c16(ld_S,S__(0,ns)
     $        ,S_sample_(nS_sample*ld_S))
         S_alpha_polar_a_sample_(nS_sample) =
     $        S_alpha_polar_a_all_(ns)
         S_alpha_azimu_b_sample_(nS_sample) =
     $        S_alpha_azimu_b_all_(ns)
      enddo

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A,I0)') 'get (n_M = '
     $     ,n_M_sample,')-x-(n_S = ',n_S_sample, ') = (', n_M_sample
     $     *n_S_sample,') innerproducts of length ' , n_w_sum_(n_k_cur
     $     -1) + n_w_(n_k_cur-1)
      timing_tic = omp_get_wtime()
      if (flag_tesselation) then
         tesselation_distance_req_use = tesselation_distance_req
      else
         tesselation_distance_req_use = 2.0d0
      end if                    !flag_tesselation
      call test_innerproduct_fast_wrapper_0(0,n_omp_sub,n_k_cur
     $     ,n_polar_a_,grid_k_p_,half_diameter_x_c_tmp,n_S_sample,ld_S
     $     ,S_sample_,tesselation_distance_req_use
     $     ,S_alpha_polar_a_sample_,S_alpha_azimu_b_sample_ ,n_M_sample
     $     ,ld_M,M_sample_ ,n_ctf,ld_CTF ,CTF_k_p_,n_alpha ,alpha1d_est_
     $     ,alpha_update_f ,n_pixels_in ,displacement_max ,n_delta_x
     $     ,n_delta_y,n_gamma_z, svd_calculation_type ,eps_svd
     $     ,fpm_howmany_max,C_M_ ,n_SM_use)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,A,F8.4)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' tesselation_distance_req_use ' ,
     $        tesselation_distance_req_use
         write(6,'(A,A,A,I0,A,I0,A,I0,A)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' calculated '
     $        , n_SM_use , ' / ' , n_M_sample *n_S_sample , ' = ' ,
     $        floor((n_SM_use*100.0d0) /(1.0d0*n_M_sample
     $        *n_S_sample)) , '% '
         write(6,'(A,A,A,F16.9)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' total_time                    ' ,timing_toc-timing_tic  
         write(6,'(A,A,A,F16.9)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' time for each S,M             ' , (timing_toc
     $        -timing_tic) / max(1,n_SM_use)
         write(6,'(A,A,A,F16.9)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' time for each S,M,dx,dy       ' , (timing_toc
     $        -timing_tic) / max(1,n_SM_use) / (n_delta_x*n_delta_y)
         write(6,'(A,A,A,F16.9)') 'test_alpha_update_SM_1: ',
     $        'test_innerproduct_fast_wrapper_0:',
     $        ' time for each S,M,dx,dy,entry ' , (timing_toc
     $        -timing_tic) / max(1,n_SM_use) / (n_delta_x*n_delta_y) /
     $        (n_w_sum_(n_k_cur -1) + n_w_(n_k_cur-1))
      end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      if (flag_fig) then
         write(tmp_fname,'(A,A,I0)') trim(dname) , '/Fig_MTS_c_' ,
     $        n_k_cur
         call Fig_gen_ver4(n_k_cur,n_polar_a_,grid_k_p_,n_S_sample,ld_S
     $        ,S_sample_,n_M_sample,ld_M,M_sample_ ,ld_CTF,CTF_k_p_
     $        ,n_alpha,alpha1d_est_ ,min(16 ,n_M_sample)
     $        ,tmp_fname)
      end if !if (flag_fig) then

      if (verbose.gt.0) write(6,'(A)')
     $     '[finished test_alpha_update_SM_1]'

      deallocate(S__)
      deallocate(S_alpha_polar_a_all_)
      deallocate(S_alpha_azimu_b_all_)
      deallocate(grid_cos_polar_a_tmp_)
      deallocate(grid_sin_polar_a_tmp_)
      deallocate(weight_cos_polar_a_tmp_)
      deallocate(azimu_b_step_tmp_)
      deallocate(n_azimu_b_tmp_)
      deallocate(I_S_sample_)
      deallocate(S_sample_)
      deallocate(S_alpha_polar_a_sample_)
      deallocate(S_alpha_azimu_b_sample_)
      
      end
