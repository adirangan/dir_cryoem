!> Doxygen comment: ;\n
!> For tesselation-tree. ;\n
!> Defines and orders children of a given tesselation element. ;\n
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
      allocate(tradius_(0:3-1));
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
         write(6,'(A,A)') trim(prefix) , 'preliminary normal nn'
      end if ! if (verbose.gt.0) then
      do nd=0,2
         vtA(nd) = v2_k(nd) - v1_k(nd)
         vtB(nd) = v0_k(nd) - v2_k(nd)
      enddo !do nd=0,2p
      call cross_r8(vtA,vtB,nn)
      call normalize_r8(3,nn)
      if (verbose.gt.1) then
         write(6,'(A,A,3F8.4)') trim(prefix) , 'preliminary nn: ',nn(0)
     $        ,nn(1),nn(2)
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
         write(6,'(A,A)') trim(prefix) , 'calculating tradius'
      end if ! if (verbose.gt.0) then
      tradius_(0) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v0_k(nd)
         tradius_(0) = tradius_(0) + vtA(nd)**2
      enddo !do nd=0,2
      tradius_(0) = dsqrt(tradius_(0))
      tradius_(1) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v1_k(nd)
         tradius_(1) = tradius_(1) + vtA(nd)**2
      enddo !do nd=0,2
      tradius_(1) = dsqrt(tradius_(1))
      tradius_(2) = 0.0d0
      do nd=0,2
         vtA(nd) = vm(nd) - v2_k(nd)
         tradius_(2) = tradius_(2) + vtA(nd)**2
      enddo !do nd=0,2
      tradius_(2) = dsqrt(tradius_(2))
      tradius = max_r8_f(3,tradius_)
      if (verbose.gt.1) then
         write(6,'(A,A,F8.3)') trim(prefix) , 'tradius: ' , tradius
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
      allocate(LT_sub_(0:n_L-1))
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
            LT_sub_(n_L_sub) = LT_(nL)
            L_sub_(0 + 3*n_L_sub) = L_(0 + 3*nL)
            L_sub_(1 + 3*n_L_sub) = L_(1 + 3*nL)
            L_sub_(2 + 3*n_L_sub) = L_(2 + 3*nL)
            n_L_sub = n_L_sub + 1
         end if !if (isin) then
      enddo !do nL=0,n_L-1
      if (verbose.gt.1) then
         write(6,'(A,A,I0,A,I0)') trim(prefix) , 'n_L: ',n_L
     $        ,'; n_L_sub: ',n_L_sub
         do nL_sub=0,n_L_sub-1
            write(6,'(A,A,I0,A,I0,3F8.4)') trim(prefix),'L_sub: ' ,
     $           nL_sub ,'LT_sub: ' , LT_sub_(nL_sub) , L_sub_(0+3
     $           *nL_sub), L_sub_(1+3*nL_sub),L_sub_(2+3*nL_sub)
         enddo !do nL_sub=0,n_L_sub-1
      end if ! if (verbose.gt.0) then
      lf = (n_L_sub.eq.1) .or. ((n_L_sub.gt.1) .and.
     $     (tradius.le.tradius_min))
      if (verbose.gt.0) then
         write(6,'(A,A,I0,1X,L1,A,I0,A,I0,A,I0,A,I0,A,F8.4,A,L1)')
     $        trim(prefix) , 
     $        ' nl_: ' , nl ,
     $        ' parity:' , parity ,
     $        ' n_L:' , n_L ,
     $        ' n_L_sub:' , n_L_sub ,
     $        ' dd:' , tradius ,
     $        ' lf:' , lf 
      end if !if (verbose.gt.0) then
