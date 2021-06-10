      recursive subroutine tesselation_index(
     $ tradius_min , nl_max , nl__in ,
     $ v0_i , v1_i , v2_i ,
     $ parity ,
     $ n_L , L_ , LT_ , TL_ , 
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
     $ T_tr_ , 
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
      implicit none
      integer *4 verbose
      data verbose / 0 /
      real *8 tradius_min
      integer *4 nl_max,nl__in,n_L
      integer *4 index_parent
      integer *4 index__input,index_output
      integer *4 ls__input,ls_output      
      real *8 v0_i(0:3-1),v1_i(0:3-1),v2_i(0:3-1)
      logical parity,lf
      real *8 L_(0:3*n_L-1)
      integer *4 LT_(0:0) ! LT_(j) = j = index of L_(j)
      integer *4 TL_(0:0) ! TL_(j) = tesselation_index associated with L_(j)
      integer *4 T_nl_(0:0) !level
      logical    T_up_(0:0) !parity
      integer *4 T_id_(0:0) !self tesselation_index
      integer *4 T_pa_(0:0) !parent tesselation_index
      real *8    T_v0_(0:0) !vertex 0
      real *8    T_v1_(0:0) !vertex 1
      real *8    T_v2_(0:0) !vertex 2
      real *8    T_vm_(0:0) !vertex center
      real *8    T_m0_(0:0) !edge midpoint 0
      real *8    T_m1_(0:0) !edge midpoint 1
      real *8    T_m2_(0:0) !edge midpoint 2
      real *8    T_e0_(0:0) !edge vector 0
      real *8    T_e1_(0:0) !edge vector 1
      real *8    T_e2_(0:0) !edge vector 2
      real *8    T_n0_(0:0) !edge normal 0
      real *8    T_n1_(0:0) !edge normal 1
      real *8    T_n2_(0:0) !edge normal 2
      real *8    T_nn_(0:0) !center normal
      real *8    T_tr_(0:0) !tradius
      integer *4 T_ll_(0:0) !number of points from L_ in T_
      logical    T_lf_(0:0) !is leaf
      integer *4 T_c0_(0:0) !child_0 tesselation_index
      integer *4 T_c1_(0:0) !child_1 tesselation_index
      integer *4 T_c2_(0:0) !child_2 tesselation_index
      integer *4 T_c3_(0:0) !child_3 tesselation_index
      integer *4 T_ls_(0:0) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4 T_LT_(0:0) !full point_index_list for all of T_ (leaves only)
      integer *4 nd
      real *8, allocatable :: v0_j(:)
      real *8, allocatable :: v1_j(:)
      real *8, allocatable :: v2_j(:)
      real *8, allocatable :: v0_k(:)
      real *8, allocatable :: v1_k(:)
      real *8, allocatable :: v2_k(:)
      real *8, allocatable :: vm(:)
      real *8, allocatable :: nn(:)
      real *8, allocatable :: e0(:)
      real *8, allocatable :: e1(:)
      real *8, allocatable :: e2(:)
      real *8, allocatable :: m0(:)
      real *8, allocatable :: m1(:)
      real *8, allocatable :: m2(:)
      real *8, allocatable :: n0(:)
      real *8, allocatable :: n1(:)
      real *8, allocatable :: n2(:)
      real *8, allocatable :: vtA(:)
      real *8, allocatable :: vtB(:)
      real *8, allocatable :: L_sub_(:)
      integer *4, allocatable :: LT_sub_(:)
      integer *4 n_L_sub,nL,nL_sub
      real *8 aa,a0,a1,a2
      logical isin
      integer index_tempA,index_tempB,ls_tempA
      integer *4 max_i4_f,sum_i4_f
      real *8, allocatable :: tradius_(:)
      real *8 tradius,max_r8_f
      character(len=1024) prefix,format_string
      write(format_string,'(A,I0,A)') '(A',nl__in+1,')'
      write(prefix,format_string) '%'      
      if (verbose.gt.0) then
         write(6,'(A,A,I0)') trim(prefix) ,
     $        '[entering tesselation_index]: ',nl__in
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         do nL=0,n_L-1
            write(6,'(A,A,I0,1X,3F8.4)') trim(prefix) , 'L_ ' , nL ,L_(0
     $           +3*nL) , L_(1+3*nL) , L_(2+3*nL)
         enddo !do nL=0,n_L-1
      end if ! if (verbose.gt.0) then
      include 'tesselation_excerpt.f'
      if (n_L_sub.eq.0) then
         if (verbose.gt.0) then
            write(6,'(A,A)') trim(prefix) , 'skipping self'
         end if !if (verbose.gt.0) then
         index_output = index__input
         ls_output = ls__input
      end if !if (n_L_sub.eq.0) then
      if (n_L_sub.ge.1) then
         if (verbose.gt.0) then
            write(6,'(A,A)') trim(prefix) , 'creating self'
         end if !if (verbose.gt.0) then
         T_nl_(1*index__input) = nl__in
         T_up_(1*index__input) = parity
         T_id_(1*index__input) = index__input
         T_pa_(1*index__input) = index_parent
         call cp1_r8(3,v0_k,T_v0_(3*index__input))
         call cp1_r8(3,v1_k,T_v1_(3*index__input))
         call cp1_r8(3,v2_k,T_v2_(3*index__input))
         call cp1_r8(3,vm,T_vm_(3*index__input))
         call cp1_r8(3,m0,T_m0_(3*index__input))
         call cp1_r8(3,m1,T_m1_(3*index__input))
         call cp1_r8(3,m2,T_m2_(3*index__input))
         call cp1_r8(3,e0,T_e0_(3*index__input))
         call cp1_r8(3,e1,T_e1_(3*index__input))
         call cp1_r8(3,e2,T_e2_(3*index__input))
         call cp1_r8(3,n0,T_n0_(3*index__input))
         call cp1_r8(3,n1,T_n1_(3*index__input))
         call cp1_r8(3,n2,T_n2_(3*index__input))
         call cp1_r8(3,nn,T_nn_(3*index__input))
         T_tr_(1*index__input) = tradius
         T_ll_(1*index__input) = n_L_sub
         T_lf_(1*index__input) = lf
         T_c0_(1*index__input) = -1
         T_c1_(1*index__input) = -1
         T_c2_(1*index__input) = -1
         T_c3_(1*index__input) = -1
         T_ls_(1*index__input) = -1
         index_output = index__input + 1
         ls_output = ls__input
         if (lf) then
            if (verbose.gt.0) then
               write(6,'(A,A)') trim(prefix) ,
     $              'leaf, copying LT_sub to T_LT_ '
            end if !if (verbose.gt.0) then
            T_ls_(1*index__input) = ls__input
            do nL_sub=0,n_L_sub-1
               T_LT_(ls__input + nL_sub) = LT_sub_(nL_sub)
               TL_(LT_sub_(nL_sub)) = index__input
            enddo !do nL_sub=0,n_L_sub-1
            ls_output = ls__input + n_L_sub
         end if !if (lf) then
         if (.not. lf) then
            if (verbose.gt.1) then
               write(6,'(A,A)') trim(prefix) , 'creating children'
            end if ! if (verbose.gt.0) then

            index_tempB = index_output
            if (verbose.gt.1) then
               write(6,'(A,A,I0,A,I0)') trim(prefix) ,
     $              'creating child_0 at index ' , index_output , ' ls '
     $              , ls_output
            end if ! if (verbose.gt.0) then
            call tesselation_index(
     $           tradius_min , nl_max , nl__in + 1 ,
     $           v0_k , m2 , m1 ,
     $           parity ,
     $           n_L_sub , L_sub_ , LT_sub_ , TL_ , 
     $           T_id_(index__input) , index_output , ls_output ,
     $           T_nl_ , 
     $           T_up_ , 
     $           T_id_ , 
     $           T_pa_ , 
     $           T_v0_ , 
     $           T_v1_ , 
     $           T_v2_ , 
     $           T_vm_ , 
     $           T_m0_ , 
     $           T_m1_ , 
     $           T_m2_ , 
     $           T_e0_ , 
     $           T_e1_ , 
     $           T_e2_ , 
     $           T_n0_ , 
     $           T_n1_ , 
     $           T_n2_ , 
     $           T_nn_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           index_tempA ,
     $           ls_tempA)
            index_output = index_tempA
            ls_output = ls_tempA
            if (index_output.gt.index_tempB) then
               T_c0_(index__input) = index_tempB
            end if !if (index_output.gt.index_tempB) then

            index_tempB = index_output
            if (verbose.gt.1) then
               write(6,'(A,A,I0,A,I0)') trim(prefix) ,
     $              'creating child_1 at index ' , index_output , ' ls '
     $              , ls_output
            end if ! if (verbose.gt.0) then
            call tesselation_index(
     $           tradius_min , nl_max , nl__in + 1 ,
     $           v1_k , m0 , m2 ,
     $           parity ,
     $           n_L_sub , L_sub_ , LT_sub_ , TL_ , 
     $           T_id_(index__input) , index_output , ls_output ,
     $           T_nl_ , 
     $           T_up_ , 
     $           T_id_ , 
     $           T_pa_ , 
     $           T_v0_ , 
     $           T_v1_ , 
     $           T_v2_ , 
     $           T_vm_ , 
     $           T_m0_ , 
     $           T_m1_ , 
     $           T_m2_ , 
     $           T_e0_ , 
     $           T_e1_ , 
     $           T_e2_ , 
     $           T_n0_ , 
     $           T_n1_ , 
     $           T_n2_ , 
     $           T_nn_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           index_tempA ,
     $           ls_tempA)
            index_output = index_tempA
            ls_output = ls_tempA
            if (index_output.gt.index_tempB) then
               T_c1_(index__input) = index_tempB
            end if !if (index_output.gt.index_tempB) then

            index_tempB = index_output
            if (verbose.gt.1) then
               write(6,'(A,A,I0,A,I0)') trim(prefix) ,
     $              'creating child_2 at index ' , index_output , ' ls '
     $              , ls_output
            end if ! if (verbose.gt.0) then
            call tesselation_index(
     $           tradius_min , nl_max , nl__in + 1 ,
     $           v2_k , m1 , m0 ,
     $           parity ,
     $           n_L_sub , L_sub_ , LT_sub_ , TL_ , 
     $           T_id_(index__input) , index_output , ls_output ,
     $           T_nl_ , 
     $           T_up_ , 
     $           T_id_ , 
     $           T_pa_ , 
     $           T_v0_ , 
     $           T_v1_ , 
     $           T_v2_ , 
     $           T_vm_ , 
     $           T_m0_ , 
     $           T_m1_ , 
     $           T_m2_ , 
     $           T_e0_ , 
     $           T_e1_ , 
     $           T_e2_ , 
     $           T_n0_ , 
     $           T_n1_ , 
     $           T_n2_ , 
     $           T_nn_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           index_tempA ,
     $           ls_tempA)
            index_output = index_tempA
            ls_output = ls_tempA
            if (index_output.gt.index_tempB) then
               T_c2_(index__input) = index_tempB
            end if !if (index_output.gt.index_tempB) then

            index_tempB = index_output
            if (verbose.gt.1) then
               write(6,'(A,A,I0,A,I0)') trim(prefix) ,
     $              'creating child_3 at index ' , index_output , ' ls '
     $              , ls_output
            end if ! if (verbose.gt.0) then
            call tesselation_index(
     $           tradius_min , nl_max , nl__in + 1 ,
     $           m0 , m1 , m2 ,
     $           .not. parity ,
     $           n_L_sub , L_sub_ , LT_sub_ , TL_ , 
     $           T_id_(index__input) , index_output , ls_output ,
     $           T_nl_ , 
     $           T_up_ , 
     $           T_id_ , 
     $           T_pa_ , 
     $           T_v0_ , 
     $           T_v1_ , 
     $           T_v2_ , 
     $           T_vm_ , 
     $           T_m0_ , 
     $           T_m1_ , 
     $           T_m2_ , 
     $           T_e0_ , 
     $           T_e1_ , 
     $           T_e2_ , 
     $           T_n0_ , 
     $           T_n1_ , 
     $           T_n2_ , 
     $           T_nn_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           index_tempA ,
     $           ls_tempA)
            index_output = index_tempA
            ls_output = ls_tempA
            if (index_output.gt.index_tempB) then
               T_c3_(index__input) = index_tempB
            end if !if (index_output.gt.index_tempB) then

         end if !if (.not. lf) then
      end if !if (n_L_sub.ge.1) then
      if (verbose.gt.0) then
         write(6,'(A,A,I0,A,I0)') trim(prefix) , ' final index_output '
     $        , index_output , ' ls_output ' , ls_output
      end if !if (verbose.gt.0) then
      end



