      recursive subroutine tesselation_size(diameter_min ,nl__in,v0_i
     $     ,v1_i,v2_i,parity,n_L,L_,LT_,nl_out,nm_out,ll_out)
      implicit none
      integer *4 verbose
      data verbose / 0 /
      real *8 diameter_min
      integer *4 nl__in,n_L,nl_out,nm_out,ll_out
      real *8 v0_i(0:3-1),v1_i(0:3-1),v2_i(0:3-1)
      logical parity,lf
      real *8 L_(0:3*n_L-1)
      integer *4 LT_(0:n_L-1)
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
      integer *4 nl_child_0
      integer *4 nl_child_1
      integer *4 nl_child_2
      integer *4 nl_child_3
      integer *4, allocatable :: nl_child_(:)
      integer *4 nm_child_0
      integer *4 nm_child_1
      integer *4 nm_child_2
      integer *4 nm_child_3
      integer *4, allocatable :: nm_child_(:)
      integer *4 ll_child_0
      integer *4 ll_child_1
      integer *4 ll_child_2
      integer *4 ll_child_3
      integer *4, allocatable :: ll_child_(:)
      integer *4 max_i4_f,sum_i4_f
      real *8, allocatable :: diameter_(:)
      real *8 diameter,max_r8_f
      character(len=64) prefix,format_string
      write(format_string,'(A,I0,A)') '(A',nl__in+1,')'
      write(prefix,format_string) '%'      
      if (verbose.gt.0) then
         write(6,'(A,A,I0)') trim(prefix) ,
     $        '[entering tesselation_size]: ',nl__in
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         do nL=0,n_L-1
            write(6,'(A,A,I0,1X,3F8.4)') trim(prefix) , 'L_ ' , nL ,L_(0
     $           +3*nL) , L_(1+3*nL) , L_(2+3*nL)
         enddo !do nL=0,n_L-1
      end if ! if (verbose.gt.0) then
      nl_out = nl__in
      nm_out = 0
      ll_out = 0
      include 'tesselation_excerpt.f'
      if (n_L_sub.eq.0) then
         nl_out = -1
         nm_out = 0
         ll_out = 0
      end if !if (n_L_sub.eq.0) then
      if (n_L_sub.ge.1) then
         if (lf) then
         nl_out = nl__in
         nm_out = 1
         ll_out = n_L_sub
         end if !if (lf) then
         if (.not. lf) then
            if (verbose.gt.1) then
               write(6,'(A,A)') trim(prefix) , 'calling children'
            end if ! if (verbose.gt.0) then
            call tesselation_size(diameter_min,nl__in+1,v0_k,m2 ,m1
     $           ,parity,n_L_sub,L_sub_,LT_sub_,nl_child_0,nm_child_0
     $           ,ll_child_0)
            call tesselation_size(diameter_min,nl__in+1,v1_k,m0 ,m2
     $           ,parity,n_L_sub,L_sub_,LT_sub_,nl_child_1,nm_child_1
     $           ,ll_child_1)
            call tesselation_size(diameter_min,nl__in+1,v2_k,m1 ,m0
     $           ,parity,n_L_sub,L_sub_,LT_sub_,nl_child_2,nm_child_2
     $           ,ll_child_2)
            call tesselation_size(diameter_min,nl__in+1,m0,m1,m2, .not.
     $           parity,n_L_sub,L_sub_,LT_sub_,nl_child_3,nm_child_3
     $           ,ll_child_3)
            allocate(nl_child_(0:3))
            nl_child_(0) = nl_child_0
            nl_child_(1) = nl_child_1
            nl_child_(2) = nl_child_2
            nl_child_(3) = nl_child_3
            if (verbose.gt.1) then
               write(6,'(A,A,I0,1X,I0,1X,I0,1X,I0)') trim(prefix) ,
     $              'nl_child_: ' , nl_child_0 , nl_child_1 ,
     $              nl_child_2 , nl_child_3
            end if ! if (verbose.gt.1) then
            allocate(nm_child_(0:3))
            nm_child_(0) = nm_child_0
            nm_child_(1) = nm_child_1
            nm_child_(2) = nm_child_2
            nm_child_(3) = nm_child_3
            if (verbose.gt.1) then
               write(6,'(A,A,I0,1X,I0,1X,I0,1X,I0)') trim(prefix) ,
     $              'nm_child_: ' , nm_child_0 , nm_child_1 ,
     $              nm_child_2 , nm_child_3
            end if ! if (verbose.gt.1) then
            allocate(ll_child_(0:3))
            ll_child_(0) = ll_child_0
            ll_child_(1) = ll_child_1
            ll_child_(2) = ll_child_2
            ll_child_(3) = ll_child_3
            if (verbose.gt.1) then
               write(6,'(A,A,I0,1X,I0,1X,I0,1X,I0)') trim(prefix) ,
     $              'll_child_: ' , ll_child_0 , ll_child_1 ,
     $              ll_child_2 , ll_child_3
            end if ! if (verbose.gt.1) then
            nl_out = max_i4_f(4,nl_child_)
            nm_out = 1 + sum_i4_f(4,nm_child_)
            ll_out = sum_i4_f(4,ll_child_)
         end if !if (.not. lf) then
      end if !if (n_L_sub.ge.1) then
      end


