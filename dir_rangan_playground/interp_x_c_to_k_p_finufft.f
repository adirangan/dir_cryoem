!> Doxygen comment: ;\n
!> Interpolate from cartesian to polar coordinates. ;\n
!> Specifically designed to interpolate ;\n
!> from x_c (real-space cartesin) ;\n
!> to k_p (fourier-space polar). ;\n
!> Uses finufft. ;\n
      subroutine interp_x_c_to_k_p_finufft(
     $     n_x1
     $     ,diameter_x1_c
     $     ,n_x2
     $     ,diameter_x2_c
     $     ,S_x_c_
     $     ,n_r
     $     ,grid_k_p_r_
     $     ,n_w_
     $     ,S_k_p_
     $     )
c$$$      Assuming that S_x_c_ is the original image ordered as:
c$$$      S_x_c_(nx1 + nx2*n_x1) = pixel nx1,nx2
c$$$      we use finufft2d2 to calculate
c$$$      S_k_p_(na) = S_k_p_(nw + n_w_csum_(nr))
c$$$      given by:
c$$$      S_k_p_(na) = 1/Z * \sum_{nx1=0,nx2=0}^{nx1=n_x1-1,nx2=n_x2-1} 
c$$$                   S_x_c_(nx1 + nx2*n_x1) 
c$$$                   * dexp(-i * (x1_(nx1)*k1+x2_(nx2)*k2))
c$$$      where: 
c$$$      Z = dsqrt(n_x1*n_x2)
c$$$      and:
c$$$      k1 = 2.0d0*pi*grid_k_p_r_(nr)*dcos(2.0d0*pi*nw/n_w_(nr))
c$$$      k2 = 2.0d0*pi*grid_k_p_r_(nr)*dsin(2.0d0*pi*nw/n_w_(nr))
c$$$      and
c$$$      x1_(nx1) = half_diameter_x1_c*(-1 + 2*nx1/n_x1)
c$$$      x2_(nx2) = half_diameter_x2_c*(-1 + 2*nx2/n_x2)
c$$$      or equivalently:
c$$$      x1_(nx1) = diameter_x1_c/n_x1*(nx1-n_x1/2)
c$$$      x2_(nx2) = diameter_x2_c/n_x2*(nx2-n_x2/2)
c$$$      where:
c$$$      nx1-n_x1/2 runs from -n_x1/2 up to +n_x1/2-1
c$$$      nx2-n_x2/2 runs from -n_x2/2 up to +n_x2/2-1
c$$$      %%%%%%%%
c$$$      Note that grid_k_p_r_ is assumed to provide the k-magnitude
c$$$      in terms of ordinary frequency (not angular frequency).
c$$$      Hence the 2.0d0*pi in the definitions of k1 and k2. ;
c$$$      %%%%%%%%
c$$$      We use use the type-2 nufft:
c$$$      int finufft2d2(BIGINT nj,FLT* xj,FLT *yj,CPX* cj,int iflag,FLT eps,
c$$$                     BIGINT ms, BIGINT mt, CPX* fk, nufft_opts opts)
c$$$      cj[j] =  SUM   fk[k1,k2] dexp(+/-i (k1 xj[j] + k2 yj[j]))      for j = 0,...,nj-1
c$$$              k1,k2
c$$$      where sum is over -ms/2 <= k1 <= (ms-1)/2, -mt/2 <= k2 <= (mt-1)/2,
c$$$      Inputs:
c$$$      nj     number of targets (int64, aka BIGINT)
c$$$      xj,yj     x,y locations of targets (each a size-nj FLT array) in [-3pi,3pi]
c$$$      fk     FLT complex array of Fourier transform values (size ms*mt,
c$$$             increasing fast in ms then slow in mt, ie Fortran ordering).
c$$$             Along each dimension the ordering is set by opts.modeord.
c$$$      iflag  if >=0, uses + sign in exponential, otherwise - sign (int)
c$$$      eps    precision requested (>1e-16)
c$$$      ms,mt  numbers of Fourier modes given in x and y (int64)
c$$$             each may be even or odd;
c$$$             in either case the mode range is integers lying in [-m/2, (m-1)/2].
c$$$      opts   struct controlling options (see finufft.h)
c$$$      Outputs:
c$$$      cj     size-nj complex FLT array of target values
c$$$             (ie, stored as 2*nj FLTs interleaving Re, Im).
c$$$      returned value - 0 if success, else see ../docs/usage.rst
c$$$      %%%%%%%%
c$$$      Because of the conventions we make regarding the image and frequencies,
c$$$      we call the type-2 nufft using:
c$$$      nj = n_w_sum = n_A
c$$$      xj = k1_*diameter_x1_c/n_x1 = 2.0d0*pi*grid_k_p_r_(nr)*dcos(2.0d0*pi*nw/n_w_(nr))*diameter_x1_c/n_x1
c$$$      yj = k2_*diameter_x2_c/n_x2 = 2.0d0*pi*grid_k_p_r_(nr)*dsin(2.0d0*pi*nw/n_w_(nr))*diameter_x2_c/n_x2
c$$$      cj = S_k_p_
c$$$      iflag = -1
c$$$      eps = 1e-6
c$$$      ms = n_x1
c$$$      mt = n_x2
c$$$      fk = S_x_c_
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_x1 !integer *4: number of pixels in the x1_c-coordinate. ;
      real *8 diameter_x1_c !real *8: diameter in x_c-direction. ;
      integer *4 n_x2 !integer *4: number of pixels in the x2_c-coordinate. ;
      real *8 diameter_x2_c !real *8: diameter in x2_c-direction. ;
      complex *16 S_x_c_(0:0) !complex *16 array (size at least n_x1*n_x2): S_x_c_(nx1 + nx2*n_x1) is the value of the image at pixel nx1,nx2. ;
      integer *4 n_r !integer *4: number of radii in k_p. ;
      real *8 grid_k_p_r_(0:0) !real *8 array (size at least n_r): grid_k_p_r_(nr) is the k-value at radius number nr. Note that, for compatibility with test_innerproduct_8, we assume that k is given in terms of ordinary frequency (not angular frequency). ;
      integer *4 n_w_(0:0) !integer *4 array (size at least n_r): n_w_(nr) is the number of points on the ring of radius grid_k_p_r_(nk). ;
      complex *16 S_k_p_(0:0) !complex *16 array (size at least n_w_sum): S_k_p_(na) is the value of the fourier-transform of the image at a particular radius and angle. ;
      complex *16, allocatable :: tmp_S_k_p_(:) !complex (16 array (size at least n_w_sum): output of type-3 finufft. ;
      real *8 pi
      integer *4 na,n_A,nr,nw !integer *4: here n_A is the cumulative sum of n_w
      real *8, allocatable :: k1_(:) !real *8 array (size at least n_A): values of k1 used in finufft. ;
      real *8, allocatable :: k2_(:) !real *8 array (size at least n_A): values of k2 used in finufft. ;
      integer *4 nx1,nx2 !temporary: integer *4 index. ;
      real *8, allocatable :: x1_(:) !real *8 array (size at least n_x1*n_x2): values of x1 used in finufft. ;
      real *8, allocatable :: x2_(:) !real *8 array (size at least n_x1*n_x2): values of x2 used in finufft. ;
      integer*4 sum_i4_f !function output ;
      real *8 omega !temporary: angle. ;
      real *8 r !temporary: angular frequency. ;
      real *8 dx1,dx2 !temporary: scaling. ;
      integer *4 iflag,ier
      parameter (iflag=-1)
      real *8 eps
      parameter (eps=1.0d-6)
      logical flag_memory_checkset
      complex *16 c16_alpha,c16_beta !temporary: complex placeholders. ;
      logical flag_2_versus_3 !logical: compare type-2 and type-3 results. ;
      parameter (flag_2_versus_3=.false.)
      real *8 fnod_c16_f , al2_c16_f !function output. ;
      pi=4.0d0*datan(1.0d0)
      dx1 = diameter_x1_c/n_x1
      dx2 = diameter_x2_c/n_x2
      n_A = sum_i4_f(n_r,n_w_)
      allocate(k1_(0:1+n_A-1))
      call cs1_r8(n_A,k1_)
      allocate(k2_(0:1+n_A-1))
      call cs1_r8(n_A,k2_)
      if (flag_2_versus_3) then
         allocate(tmp_S_k_p_(0:1+n_A-1))
         call cs1_c16(n_A,tmp_S_k_p_)
         allocate(x1_(0:1+n_x1*n_x2-1))
         call cs1_r8(n_x1*n_x2,x1_)
         allocate(x2_(0:1+n_x1*n_x2-1))
         call cs1_r8(n_x1*n_x2,x2_)
      end if !if (flag_2_versus_3) then
      call cl1_c16(n_A,S_k_p_)
      na=0
      do nr=0,n_r-1
         r = 2.0d0*pi*grid_k_p_r_(nr)
         do nw=0,n_w_(nr)-1
            omega = (2.0d0*pi*nw)/n_w_(nr)
            k1_(na) = r*dcos(omega)*dx1!/(1.0d0*n_x1/2.0d0)
            k2_(na) = r*dsin(omega)*dx2!/(1.0d0*n_x2/2.0d0)
            na = na+1
         enddo !do nw=0,n_w_(nr)-1
      enddo !do nr=0,n_r-1
      if (na.ne.n_A) then
         write(6,'(A,A)') ' Warning, na.ne.n_A in '
     $        ,'interp_x_c_to_k_p_finufft'
      end if !if (na.ne.n_A) then
      call finufft2d2_f(n_A,k1_,k2_,S_k_p_,iflag,eps,n_x1,n_x2,S_x_c_
     $     ,ier)
      if (flag_2_versus_3) then
         do nx1=0,n_x1-1
            do nx2=0,n_x2-1
               x1_(nx1+nx2*n_x1) = nx1-(1.0d0*n_x1)/2.0d0
               x2_(nx1+nx2*n_x1) = nx2-(1.0d0*n_x2)/2.0d0
            enddo               !do nx2=0,n_x2-1
         enddo                  !do nx1=0,n_x1-1
         call finufft2d3_f(n_x1*n_x2,x1_,x2_,S_x_c_,iflag,eps,n_A,k1_
     $        ,k2_,tmp_S_k_p_,ier)
         call print_sub_r8(n_x1*n_x2,x1_,' x1_: ')
         call print_sub_r8(n_x1*n_x2,x2_,' x2_: ')
         call print_sub_c16(n_A,    S_k_p_,' 2: S_k_p_: ')
         call print_sub_c16(n_A,tmp_S_k_p_,' 3: S_k_p_: ')
         write(6,'(A,F16.12)') ' type_2_versus_3: ' ,
     $        fnod_c16_f(n_A,S_k_p_,tmp_S_k_p_)
     $        / max(al2_c16_f(n_A,S_k_p_),al2_c16_f(n_A,tmp_S_k_p_))
      end if !if (flag_2_versus_3) then
      c16_alpha = dcmplx(1.0d0/dsqrt(1.0d0*n_x1*n_x2),0.0d0)
      c16_beta = dcmplx(0.0d0,0.0d0)
      call af1_c16(n_A,c16_alpha,c16_beta,S_k_p_,S_k_p_)
      flag_memory_checkset = .true.
      call cxs_r8(n_A,k1_,'k1_',flag_memory_checkset)
      call cxs_r8(n_A,k2_,'k2_',flag_memory_checkset)
      if (flag_2_versus_3) then
         call cxs_c16(n_A,tmp_S_k_p_,'tmp_S_k_p_',flag_memory_checkset)
         call cxs_r8(n_x1*n_x2,x1_,'x1_',flag_memory_checkset)
         call cxs_r8(n_x1*n_x2,x2_,'x2_',flag_memory_checkset)
      end if !if (flag_2_versus_3) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A,A)')
     $        '[checkset failed in interp_x_c_to_k_p_finufft]'
     $        ,' <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then
      end !subroutine

