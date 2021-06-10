      recursive subroutine tesselation_leafstats(
     $ nl_max , nl__in ,
     $ n_L , L_ , 
     $ index__input , 
     $ T_nl_ , 
     $ T_vm_ , 
     $ T_tr_ , 
     $ T_ll_ , 
     $ T_lf_ , 
     $ T_c0_ , 
     $ T_c1_ , 
     $ T_c2_ , 
     $ T_c3_ , 
     $ T_ls_ ,
     $ T_LT_ ,
     $ n_leaf ,
     $ n_leaf_point , 
     $ n_leaf_point2 , 
     $ n_leaf_point_min ,
     $ n_leaf_point_max )
      implicit none
      integer *4 verbose
      data verbose / 0 /
      integer *4 nl_max,nl__in,n_L,index__input
      real *8 L_(0:3*n_L-1)
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
      integer *4 n_leaf
      integer *4 n_leaf_point
      integer *4 n_leaf_point2
      integer *4 n_leaf_point_min
      integer *4 n_leaf_point_max
      integer *4 nd,nL
      real *8, allocatable :: vs(:)
      integer *4 n_L_sub,ls_start,nL_sub,LT_tmp
      character(len=1024) prefix,format_string
      write(format_string,'(A,I0,A)') '(A',nl__in+1,')'
      write(prefix,format_string) '%'      
      if (verbose.gt.0) then
         write(6,'(A,A,I0,1X,I0)') trim(prefix) ,
     $        '[entering tesselation_leafstats]: ',nl__in ,
     $        T_nl_(index__input)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         do nL=0,n_L-1
            write(6,'(A,A,I0,1X,3F8.4)') trim(prefix) , 'L_ ' , nL ,L_(0
     $           +3*nL) , L_(1+3*nL) , L_(2+3*nL)
         enddo !do nL=0,n_L-1
      end if ! if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'T_vm_: ' , T_vm_(0 + 3
     $        *index__input) , T_vm_(1 + 3*index__input) ,T_vm_(2 + 3
     $        *index__input)
      end if !if (verbose.gt.0) then
      if (T_lf_(index__input)) then
         if (verbose.gt.0) then
            write(6,'(A,A)') trim(prefix) , 'leaf, checking self'
         end if                 ! if (verbose.gt.0) then
         allocate(vs(0:3-1))
         n_L_sub = T_ll_(index__input)
         ls_start = T_ls_(index__input)
         n_leaf = n_leaf + 1
         n_leaf_point = n_leaf_point + n_L_sub
         n_leaf_point2 = n_leaf_point2 + n_L_sub**2
         if (n_L_sub.lt.n_leaf_point_min) n_leaf_point_min = n_L_sub
         if (n_L_sub.gt.n_leaf_point_max) n_leaf_point_max = n_L_sub
      end if                    !if (T_lf_(index__input)) then
      if (.not. T_lf_(index__input)) then
         if (verbose.gt.0) then
            write(6,'(A,A)') trim(prefix) ,
     $           'branch, checking children'
         end if                 ! if (verbose.gt.0) then
         if (T_c0_(index__input).ge.0) then
            call tesselation_leafstats(
     $           nl_max , nl__in+1 ,
     $           n_L , L_ , 
     $           T_c0_(index__input) , 
     $           T_nl_ , 
     $           T_vm_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           n_leaf ,
     $           n_leaf_point , 
     $           n_leaf_point2 , 
     $           n_leaf_point_min ,
     $           n_leaf_point_max )
         end if                 !if (T_c0_(index__input).ge.0) then
         if (T_c1_(index__input).ge.0) then
            call tesselation_leafstats(
     $           nl_max , nl__in+1 ,
     $           n_L , L_ , 
     $           T_c1_(index__input) , 
     $           T_nl_ , 
     $           T_vm_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           n_leaf ,
     $           n_leaf_point , 
     $           n_leaf_point2 , 
     $           n_leaf_point_min ,
     $           n_leaf_point_max )
         end if                 !if (T_c1_(index__input).ge.0) then
         if (T_c2_(index__input).ge.0) then
            call tesselation_leafstats(
     $           nl_max , nl__in+1 ,
     $           n_L , L_ , 
     $           T_c2_(index__input) , 
     $           T_nl_ , 
     $           T_vm_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           n_leaf ,
     $           n_leaf_point , 
     $           n_leaf_point2 , 
     $           n_leaf_point_min ,
     $           n_leaf_point_max )
         end if                 !if (T_c2_(index__input).ge.0) then
         if (T_c3_(index__input).ge.0) then
            call tesselation_leafstats(
     $           nl_max , nl__in+1 ,
     $           n_L , L_ , 
     $           T_c3_(index__input) , 
     $           T_nl_ , 
     $           T_vm_ , 
     $           T_tr_ , 
     $           T_ll_ , 
     $           T_lf_ , 
     $           T_c0_ , 
     $           T_c1_ , 
     $           T_c2_ , 
     $           T_c3_ , 
     $           T_ls_ ,
     $           T_LT_ ,
     $           n_leaf ,
     $           n_leaf_point , 
     $           n_leaf_point2 , 
     $           n_leaf_point_min ,
     $           n_leaf_point_max )
         end if                 !if (T_c3_(index__input).ge.0) then
      end if                    !if (.not. T_lf_(index__input)) then
      
      end



