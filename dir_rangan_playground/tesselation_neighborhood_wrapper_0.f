!> Doxygen comment: ;\n
!> For tesselation-tree. ;\n
!> Finds neighborhood of a given tesselation element. ;\n
      subroutine tesselation_neighborhood_wrapper_0(n_point,L_ ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_ ,T_ll_,T_lf_,T_c0_,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_root_base ,root_base_,vp_input_
     $     ,tesselation_distance_req,n_LT_unique,LT_unique_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_point,npoint
      real *8 L_(0:3*n_point-1)
      integer *4, allocatable :: LT_(:) ! LT_(j) = j = index of point j within L_
c$$$      real *8 diamater_min
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
      integer *4 nl_0in
      integer *4 index_0input,index_output,index_output_old
      integer *4 ls_0input,ls_output
      integer *4 index_parent
      integer *4 n_root_base,nroot_base
      real *8, allocatable :: v0_input_(:)
      real *8, allocatable :: v1_input_(:)
      real *8, allocatable :: v2_input_(:)
      logical parity_input
      real *8 vp_input_(0:2),tesselation_distance_req
      integer *4 n_LT_neighborhood,nLT_neighborhood
      integer *4, allocatable :: LT_neighborhood_(:)
      integer *4 n_LT_unique,nLT_unique
      integer *4  LT_unique_(0:0)
      logical neighborhood_unique
      character(len=16) str_tmp

c$$$      preallocated T_ structure 
      integer *4 T_nl_(0:0) !level
      real *8    T_vm_(0:0) !vertex center
      real *8    T_tr_(0:0) !tradius
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

      if (verbose.gt.0) then
         write (6,'(A)')
     $        ' [entering tesselation_neighborhood_wrapper_0]'
      end if !if (verbose.gt.0) then

      allocate(LT_(0:n_point-1)) ! LT_(j) = j = index of L_(j)
      do npoint=0,n_point-1
         LT_(npoint) = npoint
      enddo !do npoint=0,n_point-1
      if (verbose.gt.0) then
         write(6,'(A)') 'allocating neighborhood list'
      end if ! if (verbose.gt.0) then
      allocate(LT_neighborhood_(0:ll_sum-1))

      timing_tic = second()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_neighborhood_wrapper_0'
      end if ! if (verbose.gt.0) then

      do nLT_neighborhood=0,ll_sum-1
         LT_neighborhood_(nLT_neighborhood)=-1
      enddo                     !do nLT_neighborhood=0,ll_sum-1
      n_LT_neighborhood = 0
      nl_0in = 0
      do nroot_base=0,n_root_base-1
         index_0input = root_base_(nroot_base)
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') 'nroot_base: ' , nroot_base ,
     $           ' index_0input: ' , index_0input
         end if                 ! if (verbose.gt.1) then
         call tesselation_neighborhood(
     $        nl_max , nl_0in ,
     $        vp_input_ ,
     $        tesselation_distance_req ,
     $        n_point , L_ , 
     $        index_0input , 
     $        T_nl_ , 
     $        T_vm_ , 
     $        T_tr_ , 
     $        T_ll_ , 
     $        T_lf_ , 
     $        T_c0_ , 
     $        T_c1_ , 
     $        T_c2_ , 
     $        T_c3_ , 
     $        T_ls_ ,
     $        T_LT_ ,
     $        n_LT_neighborhood ,
     $        LT_neighborhood_)
      enddo                     !do nroot_base=0,n_root_base-1
      call unique_i4(n_LT_neighborhood,LT_neighborhood_,n_LT_unique
     $     ,LT_unique_)
      if (n_LT_unique.eq.n_LT_neighborhood) then
         neighborhood_unique = .true.
         write(str_tmp,'(A)') '  '
      end if                    !if (n_LT_unique.eq.n_LT_neighborhood) then
      if (n_LT_unique.lt.n_LT_neighborhood) then
         neighborhood_unique = .false.
         write(str_tmp,'(A)') ' !'
      end if                    !if (n_LT_unique.lt.n_LT_neighborhood) then
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0,A)') ' n_LT_neighborhood: ' ,
     $        n_LT_neighborhood , ' --> n_LT_unique: ' , n_LT_unique ,
     $        str_tmp
         if (neighborhood_unique.eqv..false.) then
            write(6,*) (LT_neighborhood_(nLT_neighborhood) ,
     $           nLT_neighborhood=0 , n_LT_neighborhood-1)
            write(6,*) (LT_unique_(nLT_unique) ,
     $           nLT_unique=0 , n_LT_unique-1)
         end if                 !if (neighborhood_unique.eqv..false.) then
      end if                    ! if (verbose.gt.1) then         
      if (verbose.gt.1) then
         write(6,*) (LT_neighborhood_(nLT_neighborhood) ,
     $        nLT_neighborhood=0 , n_LT_neighborhood-1)
      end if                    ! if (verbose.gt.1) then         

      if (verbose.gt.0) then
         write (6,'(A)')
     $        ' [finished tesselation_neighborhood_wrapper_0]'
      end if !if (verbose.gt.0) then

      deallocate(LT_) ! LT_(j) = j = index of L_(j)
      deallocate(LT_neighborhood_)

      end
