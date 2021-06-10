c     gfortran -w -o test_tesselation_dr.out test_tesselation_dr.f ;./test_tesselation_dr.out
      program test_tesselation_dr
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
      real *8, allocatable :: v_n00_(:)
      real *8, allocatable :: v_0n0_(:)
      real *8, allocatable :: v_00n_(:)
      real *8, allocatable :: v_p00_(:)
      real *8, allocatable :: v_0p0_(:)
      real *8, allocatable :: v_00p_(:)
      integer *4 nl_nnn,nl_nnp,nl_npn,nl_npp
      integer *4 nl_pnn,nl_pnp,nl_ppn,nl_ppp
      integer *4 nl_(0:7),nl_max
      integer *4 nm_nnn,nm_nnp,nm_npn,nm_npp
      integer *4 nm_pnn,nm_pnp,nm_ppn,nm_ppp
      integer *4 nm_(0:7),nm_sum
      integer *4 ll_nnn,ll_nnp,ll_npn,ll_npp
      integer *4 ll_pnn,ll_pnp,ll_ppn,ll_ppp
      integer *4 ll_(0:7),ll_sum
      real *8 pi
      real *8 diameter_min
      integer *4 nl,nm,ll
      integer *4 max_i4_f,sum_i4_f
      integer *4 nl__in
      integer *4 index__input,index_output,index_output_old
      integer *4 ls__input,ls_output
      integer *4 index_parent
      integer *4 n_root_base,nroot_base
      real *8, allocatable :: v0_input_(:)
      real *8, allocatable :: v1_input_(:)
      real *8, allocatable :: v2_input_(:)
      logical parity_input
      real *8, allocatable :: vp_input_(:)
      real *8 distance_req
      integer *4 n_LT_neighborhood,nLT_neighborhood
      integer *4, allocatable :: LT_neighborhood_(:)
      integer *4 n_LT_unique,nLT_unique
      integer *4, allocatable :: LT_unique_(:)
      logical neighborhood_unique
      character(len=16) str_tmp
      real *8 timing_tic,timing_toc

c$$$      T_ structure
      integer *4, allocatable :: LT_(:) ! LT_(j) = j = index of point j within L_
      integer *4, allocatable :: TL_(:) ! TL_(j) = tesselation_index associated with L_(j)
      integer *4, allocatable :: T_nl_(:) !level
      logical   , allocatable :: T_up_(:) !parity
      integer *4, allocatable :: T_id_(:) !self tesselation_index
      integer *4, allocatable :: T_pa_(:) !parent tesselation_index
      real *8   , allocatable :: T_v0_(:) !vertex 0
      real *8   , allocatable :: T_v1_(:) !vertex 1
      real *8   , allocatable :: T_v2_(:) !vertex 2
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_m0_(:) !edge midpoint 0
      real *8   , allocatable :: T_m1_(:) !edge midpoint 1
      real *8   , allocatable :: T_m2_(:) !edge midpoint 2
      real *8   , allocatable :: T_e0_(:) !edge vector 0
      real *8   , allocatable :: T_e1_(:) !edge vector 1
      real *8   , allocatable :: T_e2_(:) !edge vector 2
      real *8   , allocatable :: T_n0_(:) !edge normal 0
      real *8   , allocatable :: T_n1_(:) !edge normal 1
      real *8   , allocatable :: T_n2_(:) !edge normal 2
      real *8   , allocatable :: T_nn_(:) !center normal
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
      integer *4, allocatable :: root_base_(:) !T_ roots

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
      allocate(LT_(0:n_point-1)) ! LT_(j) = j = index of L_(j)
      do npoint=0,n_point-1
         LT_(npoint) = npoint
      enddo !do npoint=0,n_point-1
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'building list:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'defining octants'
      end if ! if (verbose.gt.0) then
      allocate(v_n00_(0:2))
      allocate(v_0n0_(0:2))
      allocate(v_00n_(0:2))
      allocate(v_p00_(0:2))
      allocate(v_0p0_(0:2))
      allocate(v_00p_(0:2))
      v_n00_(0) = -1
      v_n00_(1) = 0
      v_n00_(2) = 0
      v_0n0_(0) = 0
      v_0n0_(1) = -1
      v_0n0_(2) = 0
      v_00n_(0) = 0
      v_00n_(1) = 0
      v_00n_(2) = -1
      v_p00_(0) = +1
      v_p00_(1) = 0
      v_p00_(2) = 0
      v_0p0_(0) = 0
      v_0p0_(1) = +1
      v_0p0_(2) = 0
      v_00p_(0) = 0
      v_00p_(1) = 0
      v_00p_(2) = +1
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'defining octants:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_size'
      end if ! if (verbose.gt.0) then
      call tesselation_size(diameter_min,0,v_n00_,v_0n0_,v_00n_,.false.
     $     ,n_point,L_,LT_,nl_nnn,nm_nnn,ll_nnn)
      call tesselation_size(diameter_min,0,v_n00_,v_0n0_,v_00p_,.true.
     $     ,n_point,L_,LT_,nl_nnp,nm_nnp,ll_nnp)
      call tesselation_size(diameter_min,0,v_n00_,v_0p0_,v_00n_,.true.
     $     ,n_point,L_,LT_,nl_npn,nm_npn,ll_npn)
      call tesselation_size(diameter_min,0,v_n00_,v_0p0_,v_00p_,.false.
     $     ,n_point,L_,LT_,nl_npp,nm_npp,ll_npp)
      call tesselation_size(diameter_min,0,v_p00_,v_0n0_,v_00n_,.true.
     $     ,n_point,L_,LT_,nl_pnn,nm_pnn,ll_pnn)
      call tesselation_size(diameter_min,0,v_p00_,v_0n0_,v_00p_,.false.
     $     ,n_point,L_,LT_,nl_pnp,nm_pnp,ll_pnp)
      call tesselation_size(diameter_min,0,v_p00_,v_0p0_,v_00n_,.false.
     $     ,n_point,L_,LT_,nl_ppn,nm_ppn,ll_ppn)
      call tesselation_size(diameter_min,0,v_p00_,v_0p0_,v_00p_,.true.
     $     ,n_point,L_,LT_,nl_ppp,nm_ppp,ll_ppp)
      nl_(0) = nl_nnn
      nl_(1) = nl_nnp
      nl_(2) = nl_npn
      nl_(3) = nl_npp
      nl_(4) = nl_pnn
      nl_(5) = nl_pnp
      nl_(6) = nl_ppn
      nl_(7) = nl_ppp
      nl_max = max_i4_f(8,nl_)
      nm_(0) = nm_nnn
      nm_(1) = nm_nnp
      nm_(2) = nm_npn
      nm_(3) = nm_npp
      nm_(4) = nm_pnn
      nm_(5) = nm_pnp
      nm_(6) = nm_ppn
      nm_(7) = nm_ppp
      nm_sum = sum_i4_f(8,nm_)
      ll_(0) = ll_nnn
      ll_(1) = ll_nnp
      ll_(2) = ll_npn
      ll_(3) = ll_npp
      ll_(4) = ll_pnn
      ll_(5) = ll_pnp
      ll_(6) = ll_ppn
      ll_(7) = ll_ppp
      ll_sum = sum_i4_f(8,ll_)
      if (verbose.gt.0) then
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' nl_: '
     $        ,(nl_(nl),nl=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' nl_max: ',nl_max
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' nm_: '
     $        ,(nm_(nm),nm=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' nm_sum: ',nm_sum
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' ll_: '
     $        ,(ll_(ll),ll=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' ll_sum: ',ll_sum
      end if                    !verbose
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_size:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'allocating tesselation structure'
      end if ! if (verbose.gt.0) then
      allocate(TL_(0:n_point-1)) ! TL_(j) = tesselation_index associated with L_(j)
      allocate(T_nl_(0:1*nm_sum-1)) !level
      allocate(T_up_(0:1*nm_sum-1)) !parity
      allocate(T_id_(0:1*nm_sum-1)) !self tesselation_index
      allocate(T_pa_(0:1*nm_sum-1)) !parent tesselation_index
      allocate(T_v0_(0:3*nm_sum-1)) !vertex 0
      allocate(T_v1_(0:3*nm_sum-1)) !vertex 1
      allocate(T_v2_(0:3*nm_sum-1)) !vertex 2
      allocate(T_vm_(0:3*nm_sum-1)) !vertex center
      allocate(T_m0_(0:3*nm_sum-1)) !edge midpoint 0
      allocate(T_m1_(0:3*nm_sum-1)) !edge midpoint 1
      allocate(T_m2_(0:3*nm_sum-1)) !edge midpoint 2
      allocate(T_e0_(0:3*nm_sum-1)) !edge vector 0
      allocate(T_e1_(0:3*nm_sum-1)) !edge vector 1
      allocate(T_e2_(0:3*nm_sum-1)) !edge vector 2
      allocate(T_n0_(0:3*nm_sum-1)) !edge normal 0
      allocate(T_n1_(0:3*nm_sum-1)) !edge normal 1
      allocate(T_n2_(0:3*nm_sum-1)) !edge normal 2
      allocate(T_nn_(0:3*nm_sum-1)) !center normal
      allocate(T_dd_(0:1*nm_sum-1)) !diameter
      allocate(T_ll_(0:1*nm_sum-1)) !number of points from L_ in T_
      allocate(T_lf_(0:1*nm_sum-1)) !is leaf
      allocate(T_c0_(0:1*nm_sum-1)) !child_0 tesselation_index
      allocate(T_c1_(0:1*nm_sum-1)) !child_1 tesselation_index
      allocate(T_c2_(0:1*nm_sum-1)) !child_2 tesselation_index
      allocate(T_c3_(0:1*nm_sum-1)) !child_3 tesselation_index
      allocate(T_ls_(0:1*nm_sum-1)) !starting index of point_index_list for T_ if leaf (leaves only)
      allocate(T_LT_(0:1*ll_sum-1)) !full point_index_list for all of T_ (leaves only)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'allocating tesselation structure:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index'
      end if ! if (verbose.gt.0) then

      allocate(v0_input_(0:3-1))
      allocate(v1_input_(0:3-1))
      allocate(v2_input_(0:3-1))
      index_output = 0
      ls_output = 0
      nl = 0
      n_root_base = 0
      allocate(root_base_(0:8-1))

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant nnn'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_n00_,v0_input_)
      call cp1_r8(3,v_0n0_,v1_input_)
      call cp1_r8(3,v_00n_,v2_input_)
      parity_input = .false.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant nnp'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_n00_,v0_input_)
      call cp1_r8(3,v_0n0_,v1_input_)
      call cp1_r8(3,v_00p_,v2_input_)
      parity_input = .true.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant npn'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_n00_,v0_input_)
      call cp1_r8(3,v_0p0_,v1_input_)
      call cp1_r8(3,v_00n_,v2_input_)
      parity_input = .true.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant npp'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_n00_,v0_input_)
      call cp1_r8(3,v_0p0_,v1_input_)
      call cp1_r8(3,v_00p_,v2_input_)
      parity_input = .false.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant pnn'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_p00_,v0_input_)
      call cp1_r8(3,v_0n0_,v1_input_)
      call cp1_r8(3,v_00n_,v2_input_)
      parity_input = .true.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant pnp'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_p00_,v0_input_)
      call cp1_r8(3,v_0n0_,v1_input_)
      call cp1_r8(3,v_00p_,v2_input_)
      parity_input = .false.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant ppn'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_p00_,v0_input_)
      call cp1_r8(3,v_0p0_,v1_input_)
      call cp1_r8(3,v_00n_,v2_input_)
      parity_input = .false.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_index for octant ppp'
      end if ! if (verbose.gt.0) then      
      call cp1_r8(3,v_p00_,v0_input_)
      call cp1_r8(3,v_0p0_,v1_input_)
      call cp1_r8(3,v_00p_,v2_input_)
      parity_input = .true.
      index__input = index_output
      index_output_old = index_output
      ls__input = ls_output
      nl__in = 0
      call tesselation_index(diameter_min,nl_max,nl__in,
     $     v0_input_,v1_input_,v2_input_ , 
     $ parity_input ,
     $ n_point , L_ , LT_ , TL_ , 
     $ index_parent , index__input , ls__input ,
     $ T_nl_ , 
     $ T_up_ , 
     $ T_id_ , 
     $ T_pa_ , 
     $ T_v0_ , 
     $ T_v1_ , 
     $ T_v2_ , 
     $ T_vm_ , 
     $ T_m0_ , 
     $ T_m1_ , 
     $ T_m2_ , 
     $ T_e0_ , 
     $ T_e1_ , 
     $ T_e2_ , 
     $ T_n0_ , 
     $ T_n1_ , 
     $ T_n2_ , 
     $ T_nn_ , 
     $ T_dd_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ index_output ,
     $ ls_output)
      if (index_output.gt.index_output_old) then
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'index_output: ' , index_output_old ,
     $           ' --> ' , index_output
         end if !if (verbose.gt.0) then
         root_base_(n_root_base) = index_output_old
         n_root_base = n_root_base + 1
      end if !if (index_output.gt.index_output_old) then

      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_root_base ' , n_root_base
         write(6,'(A,I0)') 'index_output ' , index_output
         write(6,'(A,I0)') 'ls_output ' , ls_output
      end if ! if (verbose.gt.0) then      

      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_index:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      if (verbose.gt.0) then
         write(6,'(A)') 'allocating neighborhood list'
      end if ! if (verbose.gt.0) then
      allocate(LT_neighborhood_(0:ll_sum-1))
      allocate(LT_unique_(0:ll_sum-1))


      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_neighborhood'
      end if ! if (verbose.gt.0) then
      distance_req = 0.25d0
      allocate(vp_input_(0:3-1))
      do npoint=0,n_point-1
         call cp1_r8(3,L_(3*npoint),vp_input_)
         do nLT_neighborhood=0,ll_sum-1
            LT_neighborhood_(nLT_neighborhood)=-1
         enddo !do nLT_neighborhood=0,ll_sum-1
         n_LT_neighborhood = 0
         nl__in = 0
         do nroot_base=0,n_root_base-1
            index__input = root_base_(nroot_base)
            if (verbose.gt.1) then
               write(6,'(A,I0,A,I0)') 'nroot_base: ' , nroot_base ,
     $              ' index__input: ' , index__input
            end if ! if (verbose.gt.1) then
            call tesselation_neighborhood(
     $           nl_max , nl__in ,
     $           vp_input_ ,
     $           distance_req ,
     $           n_point , L_ , 
     $           index__input , 
     $           T_nl_ , 
     $           T_vm_ , 
     $           T_dd_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           n_LT_neighborhood ,
     $           LT_neighborhood_)
         enddo !do nroot_base=0,n_root_base-1
         call unique_i4(n_LT_neighborhood,LT_neighborhood_,n_LT_unique
     $        ,LT_unique_)
         if (n_LT_unique.eq.n_LT_neighborhood) then
            neighborhood_unique = .true.
            write(str_tmp,'(A)') '  '
         end if !if (n_LT_unique.eq.n_LT_neighborhood) then
         if (n_LT_unique.lt.n_LT_neighborhood) then
            neighborhood_unique = .false.
            write(str_tmp,'(A)') ' !'
         end if !if (n_LT_unique.lt.n_LT_neighborhood) then
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0,A)') 'vp ' , npoint ,
     $           ' n_LT_neighborhood: ' , n_LT_neighborhood , 
     $           ' --> n_LT_unique: ' , n_LT_unique , 
     $           str_tmp
            if (neighborhood_unique.eqv..false.) then
               write(6,*) (LT_neighborhood_(nLT_neighborhood) ,
     $              nLT_neighborhood=0 , n_LT_neighborhood-1)
               write(6,*) (LT_unique_(nLT_unique) ,
     $              nLT_unique=0 , n_LT_unique-1)
            end if !if (neighborhood_unique.eqv..false.) then
         end if ! if (verbose.gt.1) then         
         if (verbose.gt.1) then
            write(6,*) (LT_neighborhood_(nLT_neighborhood) ,
     $           nLT_neighborhood=0 , n_LT_neighborhood-1)
         end if ! if (verbose.gt.1) then         
      enddo !do npoint=0,n_point-1
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_neighborhood:'
     $        ,' total_time ',timing_toc-timing_tic
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


