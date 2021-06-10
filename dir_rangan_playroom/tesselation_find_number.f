      recursive subroutine tesselation_find_number(diameter_min
     $     ,nlevel__in,v0_i,v1_i,v2_i,parity,n_L,L_,nnumber_out)
      implicit none
      integer *4 verbose
      data verbose / 0 /
      real *8 diameter_min
      integer *4 nlevel__in,n_L,nnumber_out
      real *8 v0_i(0:3-1),v1_i(0:3-1),v2_i(0:3-1)
      logical parity
      real *8 L_(0:3*n_L-1)
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
      integer *4 n_L_sub,nL
      real *8 aa,a0,a1,a2
      logical isin
      integer *4 nnumber_child_0
      integer *4 nnumber_child_1
      integer *4 nnumber_child_2
      integer *4 nnumber_child_3
      integer *4, allocatable :: nnumber_child_(:)
      integer *4 sum_i4_f
      real *8, allocatable :: diameter_(:)
      real *8 diameter,max_r8_f
      character(len=64) prefix,format_string
      write(format_string,'(A,I0,A)') '(A',nlevel__in+1,')'
      write(prefix,format_string) '%'      
      if (verbose.gt.0) then
         write(6,'(A,A,I0)') trim(prefix) ,
     $        '[entering tesselation_find_number]: ',nlevel__in
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'allocating arrays'
      end if ! if (verbose.gt.0) then
      allocate(v0_j(0:3-1))
      allocate(v1_j(0:3-1))
      allocate(v2_j(0:3-1))
      allocate(v0_k(0:3-1))
      allocate(v1_k(0:3-1))
      allocate(v2_k(0:3-1))
      allocate(vm(0:3-1))
      allocate(nn(0:3-1))
      allocate(e0(0:3-1))
      allocate(e1(0:3-1))
      allocate(e2(0:3-1))
      allocate(m0(0:3-1))
      allocate(m1(0:3-1))
      allocate(m2(0:3-1))
      allocate(n0(0:3-1))
      allocate(n1(0:3-1))
      allocate(n2(0:3-1))
      allocate(vtA(0:3-1))
      allocate(vtB(0:3-1))
      allocate(diameter_(0:3-1));
      nnumber_out = 0
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'copying to vx_j'
      end if ! if (verbose.gt.0) then
      call cp1_r8(3,v0_i,v0_j)
      call cp1_r8(3,v1_i,v1_j)
      call cp1_r8(3,v2_i,v2_j)
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'normalizing'
      end if ! if (verbose.gt.0) then
      call normalize_r8(3,v0_j)
      call normalize_r8(3,v1_j)
      call normalize_r8(3,v2_j)
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'copying to vx_j=k'
      end if ! if (verbose.gt.0) then
      call cp1_r8(3,v0_j,v0_k)
      call cp1_r8(3,v1_j,v1_k)
      call cp1_r8(3,v2_j,v2_k)
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating midpoint vm'
      end if ! if (verbose.gt.0) then
      do nd=0,2
         vm(nd) = (v0_k(nd) + v1_k(nd) + v2_k(nd))/3.0d0
      enddo !do nd=0,2p
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'vm: ',vm(0),vm(1),vm(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating normal nn'
      end if ! if (verbose.gt.0) then
      do nd=0,2
         vtA(nd) = v2_k(nd) - v1_k(nd)
         vtB(nd) = v0_k(nd) - v2_k(nd)
      enddo !do nd=0,2p
      call cross_r8(vtA,vtB,nn)
      call normalize_r8(3,nn)
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'nn: ',nn(0),nn(1),nn(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) ,
     $        'deciding to use (v0,v1,v2) or (v1,v0,v2)'
      end if ! if (verbose.gt.0) then
      call dot_r8(vm,nn,aa)
      if (aa.ge.0.0d0) then
         call cp1_r8(3,v0_j,v0_k)
         call cp1_r8(3,v1_j,v1_k)
         call cp1_r8(3,v2_j,v2_k)
      end if !if (aa.ge.0.0d0) then
      if (aa.lt.0.0d0) then
         call cp1_r8(3,v1_j,v0_k)
         call cp1_r8(3,v0_j,v1_k)
         call cp1_r8(3,v2_j,v2_k)         
      end if !if (aa.lt.0.0d0) then
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'v0_k: ',v0_k(0),v0_k(1)
     $        ,v0_k(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'v1_k: ',v1_k(0),v1_k(1)
     $        ,v1_k(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'v2_k: ',v2_k(0),v2_k(1)
     $        ,v2_k(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating diameter'
      end if ! if (verbose.gt.0) then
      diameter_(0) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v0_k(nd)
         diameter_(0) = diameter_(0) + vtA(nd)**2
      enddo !do nd=0,2
      diameter_(0) = dsqrt(diameter_(0))
      diameter_(1) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v1_k(nd)
         diameter_(1) = diameter_(1) + vtA(nd)**2
      enddo !do nd=0,2
      diameter_(1) = dsqrt(diameter_(1))
      diameter_(2) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v2_k(nd)
         diameter_(2) = diameter_(2) + vtA(nd)**2
      enddo !do nd=0,2
      diameter_(2) = dsqrt(diameter_(2))
      diameter = max_r8_f(3,diameter_)
      if (verbose.gt.1) then
         write(6,'(A,A,F8.3)') trim(prefix) , 'diameter: ' , diameter
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating edges'
      end if ! if (verbose.gt.0) then
      do nd=0,2
         e0(nd) = v2_k(nd) - v1_k(nd)
         e1(nd) = v0_k(nd) - v2_k(nd)
         e2(nd) = v1_k(nd) - v0_k(nd)
      enddo !do nd=0,2p
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'e0: ',e0(0),e0(1),e0(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'e1: ',e1(0),e1(1),e1(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'e2: ',e2(0),e2(1),e2(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating midpoints'
      end if ! if (verbose.gt.0) then
      do nd=0,2
         m0(nd) = v2_k(nd) + v1_k(nd)
         m1(nd) = v0_k(nd) + v2_k(nd)
         m2(nd) = v1_k(nd) + v0_k(nd)
      enddo !do nd=0,2p
      call normalize_r8(3,m0)
      call normalize_r8(3,m1)
      call normalize_r8(3,m2)
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'm0: ',m0(0),m0(1),m0(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'm1: ',m1(0),m1(1),m1(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'm2: ',m2(0),m2(1),m2(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'calculating normals'
      end if ! if (verbose.gt.0) then
      call cross_r8(e0,e1,nn)
      call normalize_r8(3,nn)
      call cross_r8(v2_k,e1,n1)
      call normalize_r8(3,n1)
      call cross_r8(v0_k,e2,n2)
      call normalize_r8(3,n2)
      call cross_r8(v1_k,e0,n0)
      call normalize_r8(3,n0)
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'nn: ',nn(0),nn(1),nn(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'n0: ',n0(0),n0(1),n0(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'n1: ',n1(0),n1(1),n1(2)
         write(6,'(A,A,3F8.4)') trim(prefix) , 'n2: ',n2(0),n2(1),n2(2)
      end if ! if (verbose.gt.0) then
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'creating L_sub_'
      end if ! if (verbose.gt.0) then
      allocate(L_sub_(0:3*n_L-1))
      n_L_sub=0
      if (verbose.gt.1) then
         write(6,'(A,A)') trim(prefix) , 'filling L_sub_'
      end if ! if (verbose.gt.0) then
      do nL=0,n_L-1
         call dot_r8(n0,L_(3*nL),a0)
         call dot_r8(n1,L_(3*nL),a1)
         call dot_r8(n2,L_(3*nL),a2)
         isin = .false.
         if (parity) then
         if (a0.ge.0.0d0 .and. a1.ge.0.0d0 .and. a2.ge.0.0d0) then
            isin = .true.
         end if !if (a0.ge.0.0d0 .and. a1.ge.0.0d0 .and. a2.ge.0.0d0) then
         end if !if (parity) then
         if (.not. parity) then
         if (a0.gt.0.0d0 .and. a1.gt.0.0d0 .and. a2.gt.0.0d0) then
            isin = .true.
         end if !if (a0.gt.0.0d0 .and. a1.gt.0.0d0 .and. a2.gt.0.0d0) then
         end if !if (.not. parity) then
         if (isin) then
            L_sub_(0 + 3*n_L_sub) = L_(0 + 3*nL)
            L_sub_(1 + 3*n_L_sub) = L_(1 + 3*nL)
            L_sub_(2 + 3*n_L_sub) = L_(2 + 3*nL)
            n_L_sub = n_L_sub + 1
         end if !if (isin) then
      enddo !do nL=0,n_L-1
      if (verbose.gt.1) then
         write(6,'(A,A,I0,A,I0)') trim(prefix) , 'n_L: ',n_L
     $        ,'; n_L_sub: ',n_L_sub
         do nL=0,n_L_sub-1
            write(6,'(A,A,I0,3F8.3)') trim(prefix),'L_sub: ',nL,
     $           L_sub_(0+3*nL), L_sub_(1+3*nL), L_sub_(2+3*nL)
         enddo !do nL=0,n_L_sub-1
      end if ! if (verbose.gt.0) then
      if (n_L_sub.eq.0) then
         nnumber_out = 0
      end if !if (n_L_sub.eq.0) then
      if (n_L_sub.eq.1) then
         nnumber_out = 1
      end if !if (n_L_sub.eq.0) then
      if (n_L_sub.gt.1 .and. diameter.gt.diameter_min) then
         if (verbose.gt.1) then
            write(6,'(A,A)') trim(prefix) , 'calling children'
         end if ! if (verbose.gt.0) then
         call tesselation_find_number(diameter_min,nlevel__in+1,v0_k,m2
     $        ,m1,parity,n_L_sub,L_sub_,nnumber_child_0)
         call tesselation_find_number(diameter_min,nlevel__in+1,v1_k,m0
     $        ,m2,parity,n_L_sub,L_sub_,nnumber_child_1)
         call tesselation_find_number(diameter_min,nlevel__in+1,v2_k,m1
     $        ,m0,parity,n_L_sub,L_sub_,nnumber_child_2)
         call tesselation_find_number(diameter_min,nlevel__in+1,m0,m1,m2
     $        , .not.parity,n_L_sub,L_sub_,nnumber_child_3)
         allocate(nnumber_child_(0:3))
         nnumber_child_(0) = nnumber_child_0
         nnumber_child_(1) = nnumber_child_1
         nnumber_child_(2) = nnumber_child_2
         nnumber_child_(3) = nnumber_child_3
         if (verbose.gt.1) then
            write(6,'(A,A,I0,1X,I0,1X,I0,1X,I0)') trim(prefix) ,
     $           'nnumber_child_: ' , nnumber_child_0 , nnumber_child_1
     $           , nnumber_child_2 , nnumber_child_3
         end if ! if (verbose.gt.1) then
         nnumber_out = 1 + sum_i4_f(4,nnumber_child_)
      end if !if (n_L_sub.gt.0) then
      end
