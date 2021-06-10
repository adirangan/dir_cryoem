      subroutine Y_to_f_x_c(
     $     verbose !integer: verbosity level. ;
     $     ,n_k_p_r !integer *4: number of shells. ;
     $     ,grid_k_p_r_ !real *8 array (length n_k_p_r): radius of each shell. ;
     $     ,n_Y_l_ !integer *4 array (length n_k_p_r): index of spherical harmonic expansion on each shell. ;
     $     ,n_x_c !integer *4: number of x_c points to evaluate f_x_c_. ;
     $     ,grid_x_c_0_ !real *8 array (length n_x_c): x_c_0 positions. ;
     $     ,grid_x_c_1_ !real *8 array (length n_x_c): x_c_1 positions. ;
     $     ,grid_x_c_2_ !real *8 array (length n_x_c): x_c_2 positions. ;
     $     ,f_x_c_ !complex *16 array (length n_x_c): output values of f_x_c_. ;
     $     )
      implicit none
      integer verbose           !integer: verbosity level. ;
      integer *4 n_k_p_r        !integer *4: number of shells. ;
      real *8 grid_k_p_r_(0:0)  !real *8 array (length n_k_p_r): radius of each shell. ;
      real *8 weight_k_p_r_(0:0) !real *8 array (length n_k_p_r): quadrature weight for each shell. ;
      integer *4 n_Y_l_(0:0)    !integer *4 array (length n_k_p_r): index of spherical harmonic expansion on each shell. ;
      integer *4 n_x_c          !integer *4: number of x_c points to evaluate f_x_c_. ;
      real *8 grid_x_c_0_(0:0)  !real *8 array (length n_x_c): x_c_0 positions. ;
      real *8 grid_x_c_1_(0:0)  !real *8 array (length n_x_c): x_c_1 positions. ;
      real *8 grid_x_c_2_(0:0)  !real *8 array (length n_x_c): x_c_2 positions. ;
      complex *16 f_x_c_(0:0)   !complex *16 array (length n_x_c): output values of f_x_c_. ;
      integer *4, allocatable :: n_polar_a_(:) !integer *4 array (length n_k_p_r): number of polar_a for each shell. ;
      integer *4, allocatable :: n_polar_a_csum_(:) !integer *4 array (length n_k_p_r): cumulative sum of n_polar_a_. ;
      integer *4 n_polar_a_sum !integer *4: total number of polar_a summed over all shells. ;
      real *8, allocatable :: polar_a__(:) !real *8 array (length n_polar_a_sum): polar_a across all shells. ;
      real *8, allocatable :: cos_polar_a__(:) !real *8 array (length n_polar_a_sum): cos(polar_a) across all shells. ;
      real *8, allocatable :: weight_polar_a__(:) !real *8 array (length n_polar_a_sum): quadrature weight associated with polar_a across all shells. ;
      integer *4, allocatable :: n_azimu_b__(:) !integer *4 array (length n_polar_a_sum): number of azimu_b for each polar_a. ;
      integer *4, allocatable :: n_azimu_b__sum_(:) !integer *4 array (length n_k_p_r) sum of n_azimu_b__ over each shell (i.e., number of points on each shell). ;
      integer *4 n_azimu_b__sum_sum !integer *4: total number of points summed over all shells. ;
      real *8, allocatable :: grid_k_c_0_ !real *8 array (length n_azimu_b__sum_sum): k_c_0 positions. ;
      real *8, allocatable :: grid_k_c_1_ !real *8 array (length n_azimu_b__sum_sum): k_c_1 positions. ;
      real *8, allocatable :: grid_k_c_2_ !real *8 array (length n_azimu_b__sum_sum): k_c_2 positions. ;
      integer *4 n_polar_a_min !integer *4: parameter: minimum number of polar_a on each shell. ;
      parameter (n_polar_a_min=6)
      integer *4 n_azimu_b_min !integer *4: parameter: minimum number of azimu_b on each ring. ;
      parameter (n_azimu_b_min=3)
      integer *4 tab !temporary: index. ;
      integer *4 nk_p_r !temporary: index. ;
      integer *4 npolar_a !temporary: index. ;
      integer *4 nazimu_b__sum_sum !temporary: index. ;
      real *8 grid_k_p_r !temporary: radius. ;
      real *8 polar_a,azimu_b !temporary: angles. ;
      integer *4 sum_i4_f !function output. ;
      logical flag_memory_checkset
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering Y_to_f_x_c] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A)') ' creating grids. '
      end if !if (verbose.gt.0) then
      allocate(n_polar_a_(0:1+n_k_p_r-1))
      call cs1_i4(n_k_p_r,n_polar_a_)
      allocate(n_polar_a_csum_(0:1+n_k_p_r-1))
      call cs1_i4(n_k_p_r,n_polar_a_csum_)
      do nk_p_r=0,n_k_p_r-1
         n_polar_a_(nk_p_r) = max(n_polar_a_min,nint(2.0d0*pi
     $        *grid_k_p_r_(nk_p_r)))
      enddo !do nk_p_r=0,n_k_p_r-1
      n_polar_a_sum=0
      do nk_p_r=0,n_k_p_r-1
         n_polar_a_csum_(nk_p_r) = n_polar_a_sum
         n_polar_a_sum = n_polar_a_sum + n_polar_a_(nk_p_r)
      enddo !do nk_p_r=0,n_k_p_r-1
      if (n_polar_a_sum.ne.sum_i4_f(n_k_p_r,n_polar_a)) then
         write(6,'(A)') ' Warning, sum incorrect in Y_to_f_x_c'
      end if !if (n_polar_a_sum.ne.sum_i4_f(n_k_p_r,n_polar_a)) then
      allocate(polar_a__(0:1+n_polar_a_sum-1))
      call cs1_r8(n_polar_a_sum,polar_a__)
      allocate(cos_polar_a__(0:1+n_polar_a_sum-1))
      call cs1_r8(n_polar_a_sum,cos_polar_a__)
      do nk_p_r=0,n_k_p_r-1
         call chebexps(
     $        1
     $        ,n_polar_a_(nk_p_r)
     $        ,cos_polar_a__(n_polar_a_csum_(nk_p_r))
     $        ,0.0d0
     $        ,0.0d0
     $        ,weight_polar_a__(n_polar_a_csum_(nk_p_r))
     $        )
         do npolar_a=0,n_polar_a_(nk_p_r)
            tab = npolar_a + n_polar_a_csum_(nk_p_r)
            polar_a__(tab) = dacos(cos_polar_a__(tab))
         enddo !do npolar_a=0,n_polar_a_(nk_p_r)
      enddo !do nk_p_r=0,n_k_p_r-1
      allocate(n_azimu_b__(0:1+n_polar_a_sum-1))
      call cs1_i4(n_polar_a_sum,n_azimu_b__)
      allocate(n_azimu_b__sum_(0:1+n_k_p_r-1))
      call cs1_i4(n_k_p_r,n_azimu_b__sum_)
      do nk_p_r=0,n_k_p_r-1
         n_azimu_b__sum_(nk_p_r) = 0
         do npolar_a=0,n_polar_a_(nk_p_r)
            tab = npolar_a + n_polar_a_csum_(nk_p_r)
            n_azimu_b__(tab) = max(n_azimu_b_min,2*n_polar_a_(nk_p_r) *
     $           cos_polar_a__(tab))
         enddo !do npolar_a=0,n_polar_a_(nk_p_r)
         n_azimu_b__sum_(nk_p_r) = sum_i4_f(n_polar_a_(nk_p_r)
     $        ,n_azimu_b__(n_polar_a_csum_(nk_p_r))
      enddo !do nk_p_r=0,n_k_p_r-1
      n_azimu_b__sum_sum = sum_i4_f(n_k_p_r,n_azimu_b__sum_)
      allocate(grid_k_c_0_(0:1+n_azimu_b__sum_sum-1))
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_0_)
      allocate(grid_k_c_1_(0:1+n_azimu_b__sum_sum-1))
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_1_)
      allocate(grid_k_c_2_(0:1+n_azimu_b__sum_sum-1))
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_2_)
      do nazimu_b__sum_sum=0,n_azimu_b_sum_sum
         polar_a = 
      enddo !do nazimu_b__sum_sum=0,n_azimu_b_sum_sum
      if (verbose.gt.0) then
         call print_sub_i4(n_k_p_r,n_polar_a_,' n_polar_a_: ')
         call print_sub_r8(n_k_p_r,grid_k_p_r_,' grid_k_p_r_: ')
         call print_sub_i4(n_k_p_r,n_polar_a_csum_,' n_polar_a_csum_: ')
         call print_sub_r8(n_polar_a_sum,cos_polar_a__
     $        ,' cos_polar_a__: ')
         call print_sub_r8(n_polar_a_sum,polar_a__,' polar_a__: ')
         call print_sub_i4(n_polar_a_sum,n_azimu_b__,' n_azimu_b__: ')
         call print_sub_i4(n_k_p_r,n_azimu_b__sum_,' n_azimu_b__sum_: ')
      end if !if (verbose.gt.0) then

      if (verbose.gt.0) then
         write(6,'(A)') ' creating grids. '
      end if !if (verbose.gt.0) then

      if (verbose.gt.0) then
         write(6,'(A)') ' creating grids. '
      end if !if (verbose.gt.0) then

      if (verbose.gt.0) then
         write(6,'(A)') ' checking memory allocation '
      end if !if (verbose.gt.0) then
      flag_memory_checkset=.true.
      call cxs_i4(n_k_p_r,n_polar_a_,'n_polar_a_',flag_memory_checkset)
      call cxs_i4(n_k_p_r,n_polar_a_csum_,'n_polar_a_csum_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_polar_a_sum,polar_a__,'polar_a__'
     $     ,flag_memory_checkset)
      call cxs_r8(n_polar_a_sum,cos_polar_a__,'cos_polar_a__'
     $     ,flag_memory_checkset)
      call cxs_i4(n_polar_a_sum,n_azimu_b__,'n_azimu_b__'
     $     ,flag_memory_checkset)
      call cs1_i4(n_k_p_r,n_azimu_b__sum_,'n_azimu_b__sum_'
     $     ,flag_memory_checkset)
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_0_,'grid_k_c_0_'
     $     ,flag_memory_checkset)
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_1_,'grid_k_c_1_'
     $     ,flag_memory_checkset)
      call cs1_r8(n_azimu_b__sum_sum,grid_k_c_2_,'grid_k_c_2_'
     $     ,flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.1) then 
            write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.1) then 
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

      
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished Y_to_f_x_c] '
      end if !if (verbose.gt.0) then
      end !subroutine
      
