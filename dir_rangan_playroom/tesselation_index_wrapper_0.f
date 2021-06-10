      subroutine tesselation_index_wrapper_0(n_point,L_,diameter_min
     $     ,nl_max,nm_sum,ll_sum,T_nl_,T_vm_,T_dd_,T_ll_,T_lf_,T_c0_
     $     ,T_c1_,T_c2_,T_c3_,T_ls_,T_LT_,n_root_base,root_base_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_point,npoint
      real *8 L_(0:3*n_point-1)
      integer *4, allocatable :: LT_(:) ! LT_(j) = j = index of point j within L_
      real *8 diameter_min
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
      integer *4 nl,nm,ll
      integer *4 max_i4_f,sum_i4_f
      real *8 timing_tic,timing_toc
      integer *4 nl__in
      integer *4 index__input,index_output,index_output_old
      integer *4 ls__input,ls_output
      integer *4 index_parent
      integer *4 n_root_base,nroot_base
      real *8, allocatable :: v0_input_(:)
      real *8, allocatable :: v1_input_(:)
      real *8, allocatable :: v2_input_(:)
      logical parity_input

c$$$      preallocated T_ structure 
      integer *4 T_nl_(0:0) !level
      real *8    T_vm_(0:0) !vertex center
      real *8    T_dd_(0:0) !diameter
      integer *4 T_ll_(0:0) !number of points from L_ in T_
      logical    T_lf_(0:0) !is leaf
      integer *4 T_c0_(0:0) !child_0 tesselation_index
      integer *4 T_c1_(0:0) !child_1 tesselation_index
      integer *4 T_c2_(0:0) !child_2 tesselation_index
      integer *4 T_c3_(0:0) !child_3 tesselation_index
      integer *4 T_ls_(0:0) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4 T_LT_(0:0) !full point_index_list for all of T_ (leaves only)
c$$$      preallocated T_ roots
      integer *4 root_base_(0:0) !T_ roots

c$$$      remaining T_ structure
      integer *4, allocatable :: TL_(:) ! TL_(j) = tesselation_index associated with L_(j)
      logical   , allocatable :: T_up_(:) !parity
      integer *4, allocatable :: T_id_(:) !self tesselation_index
      integer *4, allocatable :: T_pa_(:) !parent tesselation_index
      real *8   , allocatable :: T_v0_(:) !vertex 0
      real *8   , allocatable :: T_v1_(:) !vertex 1
      real *8   , allocatable :: T_v2_(:) !vertex 2
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

      if (verbose.gt.0) then
         write (6,'(A)') ' [entering tesselation_index_wrapper_0]'
      end if !if (verbose.gt.0) then

      allocate(LT_(0:n_point-1)) ! LT_(j) = j = index of L_(j)
      do npoint=0,n_point-1
         LT_(npoint) = npoint
      enddo !do npoint=0,n_point-1

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
      call tesselation_define_octants(v_n00_,v_0n0_,v_00n_,v_p00_
     $     ,v_0p0_,v_00p_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'defining octants:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'allocating tesselation structure'
      end if ! if (verbose.gt.0) then
      allocate(TL_(0:n_point-1)) ! TL_(j) = tesselation_index associated with L_(j)
      allocate(T_up_(0:1*nm_sum-1)) !parity
      allocate(T_id_(0:1*nm_sum-1)) !self tesselation_index
      allocate(T_pa_(0:1*nm_sum-1)) !parent tesselation_index
      allocate(T_v0_(0:3*nm_sum-1)) !vertex 0
      allocate(T_v1_(0:3*nm_sum-1)) !vertex 1
      allocate(T_v2_(0:3*nm_sum-1)) !vertex 2
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
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'allocating tesselation structure:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

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
         write (6,'(A)') ' [finished tesselation_index_wrapper_0]'
      end if !if (verbose.gt.0) then

      deallocate(T_up_) !parity
      deallocate(T_id_) !self tesselation_index
      deallocate(T_pa_) !parent tesselation_index
      deallocate(T_v0_) !vertex 0
      deallocate(T_v1_) !vertex 1
      deallocate(T_v2_) !vertex 2
      deallocate(T_m0_) !edge midpoint 0
      deallocate(T_m1_) !edge midpoint 1
      deallocate(T_m2_) !edge midpoint 2
      deallocate(T_e0_) !edge vector 0
      deallocate(T_e1_) !edge vector 1
      deallocate(T_e2_) !edge vector 2
      deallocate(T_n0_) !edge normal 0
      deallocate(T_n1_) !edge normal 1
      deallocate(T_n2_) !edge normal 2
      deallocate(T_nn_) !center normal

      deallocate(v_n00_)
      deallocate(v_0n0_)
      deallocate(v_00n_)
      deallocate(v_p00_)
      deallocate(v_0p0_)
      deallocate(v_00p_)

      deallocate(LT_) ! LT_(j) = j = index of L_(j)

      end
