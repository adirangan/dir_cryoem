c     gfortran -w -o test_tesselation_ver1_dr.out test_tesselation_ver1_dr.f ;./test_tesselation_ver1_dr.out
      program test_tesselation_ver1_dr
      implicit none
      integer verbose
      data verbose / 1 /
c$$$      indices
      integer *4 n_point,npoint
c$$$      arrays
      real *8, allocatable :: x_(:)
      real *8, allocatable :: y_(:)
      real *8, allocatable :: z_(:)
      real *8, allocatable :: L_(:)
      real *8 pi
      real *8 diameter_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
      real *8, allocatable :: vp_input_(:)
      real *8 distance_req
      integer *4 n_LT_unique,nLT_unique
      integer *4, allocatable :: LT_unique_(:)
      character(len=16) str_tmp
      real *8 timing_tic,timing_toc

c$$$      T_ structure
      integer *4, allocatable :: T_nl_(:) !level
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_dd_(:) !diameter
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

      pi = 4*atan(1.0)
      diameter_min = 0.5**12
      n_point = 128*4

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'building list'
      end if ! if (verbose.gt.0) then
      allocate(x_(0:n_point-1))
      allocate(y_(0:n_point-1))
      allocate(z_(0:n_point-1))
      allocate(L_(0:3*n_point-1))
      call linspace(-2*pi,+2*pi,n_point,x_)
      call linspace(-2*pi,+2*pi,n_point,y_)
      call linspace(-2*pi,+2*pi,n_point,z_)
      do npoint=0,n_point-1
         L_(0 + 3*npoint) = dcos(x_(npoint))
         L_(1 + 3*npoint) = dsin(y_(npoint))
         L_(2 + 3*npoint) = z_(npoint)
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
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'building list:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_get_nl_nm_ll'
      end if ! if (verbose.gt.0) then
      call tesselation_get_nl_nm_ll(n_point,L_,diameter_min,nl_max
     $     ,nm_sum,ll_sum)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_get_nl_nm_ll:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'allocating tesselation structure'
      end if ! if (verbose.gt.0) then
      allocate(T_nl_(0:1*nm_sum-1)) !level
      allocate(T_vm_(0:3*nm_sum-1)) !vertex center
      allocate(T_dd_(0:1*nm_sum-1)) !diameter
      allocate(T_ll_(0:1*nm_sum-1)) !number of points from L_ in T_
      allocate(T_lf_(0:1*nm_sum-1)) !is leaf
      allocate(T_c0_(0:1*nm_sum-1)) !child_0 tesselation_index
      allocate(T_c1_(0:1*nm_sum-1)) !child_1 tesselation_index
      allocate(T_c2_(0:1*nm_sum-1)) !child_2 tesselation_index
      allocate(T_c3_(0:1*nm_sum-1)) !child_3 tesselation_index
      allocate(T_ls_(0:1*nm_sum-1)) !starting index of point_index_list for T_ if leaf (leaves only)
      allocate(T_LT_(0:1*ll_sum-1)) !full point_index_list for all of T_ (leaves only)
      allocate(T_root_base_(0:8-1))
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'allocating tesselation structure:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index_wrapper_0'
      end if ! if (verbose.gt.0) then
      call tesselation_index_wrapper_0(n_point,L_,diameter_min ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_dd_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_index_wrapper_0:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      if (verbose.gt.0) then
         write(6,'(A)') 'allocating neighborhood list'
      end if ! if (verbose.gt.0) then
      allocate(LT_unique_(0:ll_sum-1))

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_neighborhood_wrapper_0'
      end if ! if (verbose.gt.0) then
      distance_req = 0.25d0
      allocate(vp_input_(0:3-1))
      do npoint=0,n_point-1
         call cp1_r8(3,L_(3*npoint),vp_input_)
         call tesselation_neighborhood_wrapper_0(n_point,L_ ,nl_max
     $        ,nm_sum,ll_sum,T_nl_,T_vm_,T_dd_ ,T_ll_,T_lf_,T_c0_,T_c1_
     $        ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base ,T_root_base_
     $        ,vp_input_ ,distance_req,n_LT_unique,LT_unique_)
      enddo !do npoint=0,n_point-1
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)')
     $        'calling tesselation_neighborhood_wrapper_0:',' total_time
     $',timing_toc-timing_tic
         write(6,'(A,F8.4)') 'each: ' , (timing_toc-timing_tic)/n_point
      end if

 10   stop
      end

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
      include 'tesselation_neighborhood.f'
      include 'tesselation_define_octants.f'
      include 'tesselation_get_nl_nm_ll.f'
      include 'tesselation_index_wrapper_0.f'
      include 'tesselation_neighborhood_wrapper_0.f'


