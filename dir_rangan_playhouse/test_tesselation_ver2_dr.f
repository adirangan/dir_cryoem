c     gfortran -w --openmp -o test_tesselation_ver2_dr.out test_tesselation_ver2_dr.f ;./test_tesselation_ver2_dr.out -
      program test_tesselation_ver2_dr
      implicit none
      include 'omp_lib.h'
      integer verbose
      data verbose / 1 /
c$$$      indices
      integer *4 n_point,npoint
c$$$      arrays
c$$$      real *8, allocatable :: x_(:)
c$$$      real *8, allocatable :: y_(:)
c$$$      real *8, allocatable :: z_(:)
      real *8, allocatable :: L_(:)
      real *8 pi
      real *8 tradius_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
      real *8, allocatable :: vp_input_(:)
      real *8 distance_req
      integer *4 n_LT_unique,nLT_unique,n_LT_unique_sum
      integer *4, allocatable :: LT_unique_(:)
      integer *4 n_lLT_unique,nlLT_unique,n_lLT_unique_sum
     $     ,n_lLT_unique_min,n_lLT_unique_max
      integer *4, allocatable :: lLT_unique_(:)
      character(len=16) str_tmp
      integer *4 rseed
      real *8 adi_rand_f
      real *8 timing_tic,timing_toc

c$$$      T_ structure
      integer *4, allocatable :: T_nl_(:) !level
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_tr_(:) !tradius
      integer *4, allocatable :: T_ll_(:) !number of points from L_ in T_
      logical   , allocatable :: T_lf_(:) !is leaf
      integer *4, allocatable :: T_c0_(:) !child_0 tesselation_index
      integer *4, allocatable :: T_c1_(:) !child_1 tesselation_index
      integer *4, allocatable :: T_c2_(:) !child_2 tesselation_index
      integer *4, allocatable :: T_c3_(:) !child_3 tesselation_index
      integer *4, allocatable :: T_ls_(:) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4, allocatable :: T_LT_(:) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4, allocatable :: T_root_base_(:) !T_ roots
c$$$      T_ roots
      integer *4 n_leaf,n_leaf_point,n_leaf_point2,n_leaf_point_min
     $     ,n_leaf_point_max
c$$$      memory map
      real *8 d_mem

      rseed = 1
      pi = 4*atan(1.0)
      tradius_min = 0.5d0**12
c$$$      tradius_min = 0.5d0
      n_point = 1024*2

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'building list with n_point ' , n_point
      end if ! if (verbose.gt.0) then
c$$$      allocate(x_(0:n_point-1))
c$$$      allocate(y_(0:n_point-1))
c$$$      allocate(z_(0:n_point-1))
      d_mem = 0.0d0
      d_mem = d_mem + 3.0d0*n_point*8.0d0
      allocate(L_(0:3*n_point-1))
c$$$      call linspace(-2*pi,+2*pi,n_point,x_)
c$$$      call linspace(-2*pi,+2*pi,n_point,y_)
c$$$      call linspace(-2*pi,+2*pi,n_point,z_)
      do npoint=0,n_point-1
         L_(0 + 3*npoint) = adi_rand_f(rseed) - 0.5d0
         L_(1 + 3*npoint) = adi_rand_f(rseed) - 0.5d0
         L_(2 + 3*npoint) = adi_rand_f(rseed) - 0.5d0
c$$$         L_(0 + 3*npoint) = dcos(x_(npoint))
c$$$         L_(1 + 3*npoint) = dsin(y_(npoint))
c$$$         L_(2 + 3*npoint) = z_(npoint)
c$$$         L_(0 + 3*npoint) = rand()-0.5
c$$$         L_(1 + 3*npoint) = rand()-0.5
c$$$         L_(2 + 3*npoint) = rand()-0.5
         call normalize_r8(3,L_(0 + 3*npoint))
      enddo !do npoint=0,n_point-1
      if (verbose.gt.1) then
         do npoint=0,n_point-1
            write(6,'(A,I0,1X,3F8.4)') ' L_ ', npoint , L_(0 + 3*npoint)
     $           , L_(1 + 3*npoint) , L_(2 + 3*npoint)
         enddo ! do npoint=0,n_point-1
      end if! if (verbose.gt.1) then
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'building list:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_get_nl_nm_ll'
      end if ! if (verbose.gt.0) then
      call tesselation_get_nl_nm_ll(n_point,L_,tradius_min,nl_max
     $     ,nm_sum,ll_sum)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,F8.4,3(A,I0))') 'tradius_min: ' , tradius_min ,
     $        ' nl_max ' , nl_max , ' nm_sum ' , nm_sum , ' ll_sum ' ,
     $        ll_sum
      end if
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_get_nl_nm_ll:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A)') 'allocating tesselation structure'
      end if ! if (verbose.gt.0) then
      d_mem = d_mem + 8.0*nm_sum*4.0d0
      d_mem = d_mem + 4.0*nm_sum*8.0d0
      d_mem = d_mem + 1.0*ll_sum*4.0d0
      if (verbose.gt.0) then
         write (6,'(A,F8.4,A)') ' tesselation requires ' , d_mem*1.0d-6
     $        ,' MB'
      end if !if (verbose.gt.0) then
      allocate(T_nl_(0:1*nm_sum-1)) !level
      allocate(T_vm_(0:3*nm_sum-1)) !vertex center
      allocate(T_tr_(0:1*nm_sum-1)) !tradius
      allocate(T_ll_(0:1*nm_sum-1)) !number of points from L_ in T_
      allocate(T_lf_(0:1*nm_sum-1)) !is leaf
      allocate(T_c0_(0:1*nm_sum-1)) !child_0 tesselation_index
      allocate(T_c1_(0:1*nm_sum-1)) !child_1 tesselation_index
      allocate(T_c2_(0:1*nm_sum-1)) !child_2 tesselation_index
      allocate(T_c3_(0:1*nm_sum-1)) !child_3 tesselation_index
      allocate(T_ls_(0:1*nm_sum-1)) !starting index of point_index_list for T_ if leaf (leaves only)
      allocate(T_LT_(0:1*ll_sum-1)) !full point_index_list for all of T_ (leaves only)
      allocate(T_root_base_(0:8-1))
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'allocating tesselation structure:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = omp_get_wtime()

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index_wrapper_0'
      end if ! if (verbose.gt.0) then
      call tesselation_index_wrapper_0(n_point,L_,tradius_min ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_index_wrapper_0:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_leafstats_wrapper_0'
      end if ! if (verbose.gt.0) then
      call tesselation_leafstats_wrapper_0(n_point,L_ ,nl_max ,nm_sum
     $     ,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_,T_c1_ ,T_c2_
     $     ,T_c3_,T_ls_,T_LT_,n_T_root_base ,T_root_base_ ,n_leaf
     $     ,n_leaf_point,n_leaf_point2,n_leaf_point_min
     $     ,n_leaf_point_max)
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_leaf: ' , n_leaf
         write(6,'(A,I0)') 'n_leaf_point: ' , n_leaf_point
         write(6,'(A,I0)') 'n_leaf_point2: ' , n_leaf_point2
         write(6,'(A,I0)') 'n_leaf_point_min: ' , n_leaf_point_min
         write(6,'(A,I0)') 'n_leaf_point_max: ' , n_leaf_point_max
         write(6,'(A,F8.1)') 'n_leaf_point_avg: ' , (n_leaf_point*1.0d0)
     $        /max(1,n_leaf)
         write(6,'(A,F8.1)') 'n_leaf_point_std: ' , dsqrt((n_leaf_point2
     $        *1.0d0)/max(1,n_leaf) - ((n_leaf_point*1.0d0)/max(1
     $        ,n_leaf))**2)
      end if !if (verbose.gt.0) then
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)')
     $        'calling tesselation_leafstats_wrapper_0:',' total_time
     $',timing_toc-timing_tic
         write(6,'(A,F8.4)') 'each: ' , (timing_toc-timing_tic)/n_point
      end if

      if (verbose.gt.0) then
         write(6,'(A)') 'allocating neighborhood list'
      end if ! if (verbose.gt.0) then
      allocate(LT_unique_(0:ll_sum-1))

      timing_tic = omp_get_wtime()
      distance_req = 2.0d0/dsqrt(1.0d0*n_point)
      if (verbose.gt.0) then
         write(6,'(A,A,F8.6)')
     $        'calling tesselation_neighborhood_wrapper_0 ' ,
     $        'with distance_req ' , distance_req
      end if ! if (verbose.gt.0) then
      allocate(vp_input_(0:3-1))
      n_LT_unique_sum = 0
      do npoint=0,n_point-1
         call cp1_r8(3,L_(3*npoint),vp_input_)
         call tesselation_neighborhood_wrapper_0(n_point,L_ ,nl_max
     $        ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_,T_c1_
     $        ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base ,T_root_base_
     $        ,vp_input_ ,distance_req,n_LT_unique,LT_unique_)
         n_LT_unique_sum = n_LT_unique_sum + n_LT_unique
      enddo !do npoint=0,n_point-1
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_LT_unique_sum ' , n_LT_unique_sum
         write(6,'(A,F8.1)') 'n_LT_unique_avg ' , (n_LT_unique_sum
     $        *1.0d0)/(1.0d0*n_point)
      end if
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)')
     $        'calling tesselation_neighborhood_wrapper_0:',' total_time
     $',timing_toc-timing_tic
         write(6,'(A,F8.4)') 'each: ' , (timing_toc-timing_tic)/n_point
      end if

      if (verbose.gt.0) then
         write(6,'(A)') 'allocating leafneighborhood list'
      end if ! if (verbose.gt.0) then
      allocate(lLT_unique_(0:n_leaf-1))

      timing_tic = omp_get_wtime()
      distance_req = 2.0d0/dsqrt(1.0d0*n_point)
      if (verbose.gt.0) then
         write(6,'(A,A,F8.6)')
     $        'calling tesselation_leafneighborhood_wrapper_0 ' ,
     $        'with distance_req ' , distance_req
      end if ! if (verbose.gt.0) then
      n_lLT_unique_sum = 0
      n_lLT_unique_min = n_leaf
      n_lLT_unique_max = 0
      do npoint=0,n_point-1
         call cp1_r8(3,L_(3*npoint),vp_input_)
         call tesselation_leafneighborhood_wrapper_0(n_point,L_ ,nl_max
     $        ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_,T_c1_
     $        ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base ,T_root_base_
     $        ,vp_input_ ,distance_req,n_lLT_unique,lLT_unique_)
         n_lLT_unique_sum = n_lLT_unique_sum + n_lLT_unique
         if (n_lLT_unique.lt.n_lLT_unique_min) n_lLT_unique_min =
     $        n_lLT_unique
         if (n_lLT_unique.gt.n_lLT_unique_max) n_lLT_unique_max =
     $        n_lLT_unique
      enddo !do npoint=0,n_point-1
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_lLT_unique_min ' , n_lLT_unique_min
         write(6,'(A,I0)') 'n_lLT_unique_max ' , n_lLT_unique_max
         write(6,'(A,I0)') 'n_lLT_unique_sum ' , n_lLT_unique_sum
         write(6,'(A,F8.1)') 'n_lLT_unique_avg ' , (n_lLT_unique_sum
     $        *1.0d0)/(1.0d0*n_point)
      end if
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)')
     $        'calling tesselation_leafneighborhood_wrapper_0:'
     $        ,' total_time ',timing_toc-timing_tic
         write(6,'(A,F8.4)') 'each: ' , (timing_toc-timing_tic)/n_point
      end if

 10   stop
      end

      include 'adi_rand_f.f'
      include 'max_i4_f.f'
      include 'max_r8_f.f'
      include 'sum_i4_f.f'
      include 'cp1_i4.f'
      include 'cp1_r8.f'
      include 'cp1_c16.f'
      include 'af1_c16.f'
      include 'quicksort_i4.f'
      include 'unique_i4.f'
      include 'linspace.f'
      include 'cross_r8.f'
      include 'dot_r8.f'
      include 'distance_r8.f'
      include 'normalize_r8.f'
      include 'tesselation_size.f'
      include 'tesselation_index.f'
      include 'tesselation_define_octants.f'
      include 'tesselation_get_nl_nm_ll.f'
      include 'tesselation_index_wrapper_0.f'
      include 'tesselation_neighborhood.f'
      include 'tesselation_neighborhood_wrapper_0.f'
      include 'tesselation_leafstats.f'
      include 'tesselation_leafstats_wrapper_0.f'
      include 'tesselation_leafneighborhood.f'
      include 'tesselation_leafneighborhood_wrapper_0.f'


