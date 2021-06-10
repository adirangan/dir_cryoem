      subroutine innerproduct_q_delta_gamma(n_svd_r,svd_r_,n_svd_d
     $     ,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma,gamma_,n_r,grid_p_,n_w_
     $     ,Wlast_,n_A,S_q_,M_q_,C_Z,Z_q_)
c$$$      Uses svd-expansion defined via n_svd_r, .. , svd_V_r_
      implicit none
      integer n_svd_r,n_svd_d,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_d_(0:n_svd_d-1)
      real *8 svd_s_(0:n_svd_l-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      integer n_delta_x,n_delta_y,n_gamma
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_(0:n_gamma-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 Wlast_(0:4*n_w_(n_r-1)+16-1)
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16 C_Z,Z_q_(0:n_delta_x*n_delta_y*n_gamma-1)
      real *8 pi
      complex *16, allocatable :: Z_svds_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer ndx,ndy,ngamma,nw,n_w_max,l
      real *8 gamma
      integer nC1,nC2,nC3,nC4,nC5
      complex *16 Z_q
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      allocate(Z_svds_(0:n_svd_l*n_w_max-1))
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      allocate(Z_tmpC_(0:n_w_max*n_delta_x*n_delta_y-1))
      call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_
     $     ,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_ ,M_q_,Z_svds_)
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      nC1 = 0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
c$$$            nC4 = ndx*n_w_max + ndy*n_w_max*n_delta_x
            nC4 = nC1
            do nw=0,n_w_max-1
c$$$               nC1 = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
               Z_tmpC_(nC1) = cmplx( 0.0 , 0.0 )
               nC2 = ndx*n_svd_l + ndy*n_svd_l*n_delta_x
               nC3 = nw*n_svd_l
               do l=0,n_svd_l-1
                  Z_tmpC_(nC1) = Z_tmpC_(nC1) + Z_svdd_(nC2)
     $                 *Z_svds_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo
               nC1 = nC1 + 1
            enddo
            call dcfftb(n_w_max,Z_tmpC_(nC4),Wlast_)
c$$$            call adi_fft1(+1,n_w_max,Z_tmpC_,nC4,1,Z_tmpC_,nC4,1)
c$$$            call affine_c16(n_w_max,dcmplx((1.0d0*n_w_max),0.0)
c$$$     $           ,dcmplx(0.0,0.0),n_w_max,Z_tmpC_,nC4,1,n_w_max,Z_tmpC_
c$$$     $           ,nC4,1)
            do ngamma=0,n_gamma-1
               gamma = gamma_(ngamma)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmpC_(nC4),+gamma
     $              ,Z_q)
               nC5 = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               Z_q_(nC5) = Z_q/(n_r**4)/C_Z
c$$$               Z_q_(nC5) = Z_q
            enddo
         enddo
      enddo
      deallocate(Z_svds_)
      deallocate(Z_svdd_)
      deallocate(Z_tmpC_)
      end
