!> Doxygen comment: ;\n
!> Simple driver for test_residual_5.f. ;\n
!> Only checks a single special case. ;\n
c$$$      gfortran -o test_residual_5_dr.out test_residual_5_dr.f --openmp -L/home/adi/finufft-master/lib/ -lfinufft -lfftw3 -lfftw3_threads -fopenmp -lopenblas -lstdc++
      implicit none
      include 'omp_lib.h'
      include 'excerpt_define_nalpha.f'

      integer *4 rseed !random seed. see command line inputs. ;
      integer *4 n_M,ld_M,n_CTF,ld_CTF,quadrature_type_azimu_b
      integer *4, allocatable :: I_M_sample_(:) !temporary: indexing variable used to reference images. ;
      complex *16, allocatable :: M_k_p__(:) ! holds images. ;
      complex *16, allocatable :: CTF_k_p__(:) ! holds images. ;
      real *8, allocatable :: alpha2d__(:,:) ! holds image parameters. ;
      integer, allocatable :: n_w_csum_(:) ! cumulative sum of n_w_. also called icstart(nk). ;
      integer, allocatable :: n_polar_a_(:) ! number of polar_a values on sphere of radius grid_k_p_(nk). also called nlats(nk). ;
      real *8, allocatable :: grid_k_p_(:) ! values for k on successive shells for k-space polar coordinates. sometimes called xnodesr(nk). ;
      integer, allocatable :: n_w_(:) ! ngridc(nk) is the number of points on the ring of radius grid_k_p_(nk) for a template or image in k-space polar coordinates. also called ngridc(nk). ;
      integer *4 n_k_p_r_low,n_k_p_r_cur,n_Y_lm_csum_cur
      integer, allocatable :: n_Y_lm_csum_(:) ! cumulative sum of n_Y_lm_. also called isph_start(nk). ;
      integer, allocatable :: n_Y_l_(:) ! order of spherical harmonic expansion on sphere at radius grid_k_p_(nk). also called nterms_sph(nk). ;
      real *8 lsq_oversample
      integer lsq_interpolation_order
      complex *16, allocatable :: Y_(:)
      complex *16, allocatable :: M_residual__(:)
      integer n_residual_loading
      integer n_residual_iteration
      complex *16, allocatable :: M_loading_(:)
      
      n_M = 134
      ld_M = 44
      n_CTF = 3
      ld_CTF = 44
      allocate(I_M_sample_(0:n_M-1));
      call cl1_i4(n_M,I_M_sample_)
      allocate(M_k_p__(0:n_M*ld_M-1))
      call cl1_c16(n_M*ld_M,M_k_p__)
      allocate(CTF_k_p__(0:n_CTF*ld_CTF-1))
      call cl1_c16(n_CTF*ld_CTF,CTF_k_p__)
      allocate(alpha2d__(0:n_alpha-1,0:n_M*ld_M-1))
      call cl1_r8(n_alpha*n_M*ld_M,alpha2d__)
      n_k_p_r_cur = 3
      allocate(n_w_csum_(0:n_k_p_r_cur-1))
      n_w_csum_(0) = 1
      n_w_csum_(1) = 13
      n_w_csum_(2) = 25
      allocate(n_polar_a_(0:n_k_p_r_cur-1))
      n_polar_a_(0) = 6
      n_polar_a_(1) = 6
      n_polar_a_(2) = 10
      quadrature_type_azimu_b = 1
      allocate(grid_k_p_(0:n_k_p_r_cur-1))
      grid_k_p_(0) = 1.0d0
      grid_k_p_(1) = 2.0d0
      grid_k_p_(2) = 3.0d0
      allocate(n_w_(0:n_k_p_r_cur-1))
      n_w_(0) = 12
      n_w_(1) = 12
      n_w_(2) = 20
      n_k_p_r_low = 1
      n_Y_lm_csum_cur = 77
      allocate(n_Y_lm_csum_(0:n_k_p_r_cur-1))
      n_Y_lm_csum_(0) = 0
      n_Y_lm_csum_(1) = 16
      n_Y_lm_csum_(2) = 41
      allocate(n_Y_l_(0:n_k_p_r_cur-1))
      n_Y_l_(0) = 3
      n_Y_l_(1) = 4
      n_Y_l_(2) = 5
      lsq_oversample = 2.0d0
      lsq_interpolation_order = 5
      n_residual_loading = 2
      n_residual_iteration = 7
      allocate(Y_(0:n_Y_lm_csum_cur-1));
      call cl1_c16(n_Y_lm_csum_cur,Y_)
      allocate(M_residual__(0:n_M*ld_M-1))
      call cl1_c16(n_M*ld_M,M_residual__)
      allocate(M_loading_(0:n_residual_loading*n_M-1))
      call cl1_c16(n_residual_loading*n_M,M_loading_)

      call test_residual_5(
      rseed
     $     ,n_M
     $     ,I_M_sample_
     $     ,ld_M
     $     ,M_k_p__
     $     ,n_CTF
     $     ,ld_CTF
     $     ,CTF_k_p__
     $     ,alpha2d__
     $     ,n_w_csum_
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,grid_k_p_
     $     ,n_w_
     $     ,n_k_p_r_low
     $     ,n_k_p_r_cur
     $     ,n_Y_lm_csum_cur
     $     ,n_Y_lm_csum_
     $     ,n_Y_l_
     $     ,lsq_oversample
     $     ,lsq_interpolation_order
     $     ,Y_
     $     ,M_residual__
     $     ,n_residual_loading
     $     ,n_residual_iteration
     $     ,M_loading_
     $     )

      stop
      end

      include 'test_residual_5.f'
      include 'cl1_i4.f'
      include 'cl1_r8.f'
      include 'cl1_c16.f'
      include 'cs1_i4.f'
      include 'cs1_r8.f'
      include 'cs1_c16.f'
      include 'cxs_i4.f'
      include 'cxs_r8.f'
      include 'cxs_c16.f'
      include 'normalize_c16.f'
      include 'print_xxx_xx.f'
      include 'ac1_c16.f'
      include 'al2_c16_f.f'
      include 'dot_c16.f'
      include 'dot_c16_f.f'
      include 'adi_rand_f.f'
