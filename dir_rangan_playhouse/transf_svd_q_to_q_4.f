      subroutine transf_svd_q_to_q_4(svd_r_max,n_svd_r,svd_r_,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,svd_s_,svd_V_r_,svd_polyval_V_r_,n_r ,grid_p_,n_w_,n_A,S_q_
     $     ,delta_x,delta_y ,M_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as S_q_
c$$$      Currently brute force convolution using svd-expansion
c$$$      defined via n_svd_r, .. , svd_V_r_
      implicit none
      logical warning_flag
      data warning_flag / .true. /
      integer n_svd_r,n_svd_d,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_d_(0:n_svd_d-1)
      real *8 svd_s_(0:n_svd_l-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      real *8 svd_polyval_U_d_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      real *8 svd_polyval_V_r_(0:n_svd_l*n_r-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1),delta_x,delta_y
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      real *8 pi
      real *8 svd_d_max,svd_d_m,svd_d_c,svd_d(0:0)
      real *8 svd_r_max,svd_r_m,svd_r_c,svd_r(0:0)
      real *8 D_V_r,D_U_d,D_s
      integer I_l
      real *16 R_q,delta,omega,theta
      integer nl,nr,ic,nw,nwt
      complex *16 C_q,C_w
      complex *16, allocatable :: C_w_(:)
c$$$      real *8, allocatable :: U_d_(:)
c$$$      real *8, allocatable :: V_r_(:)
      pi = 4.0*atan(1.0)
      allocate(Z_q_(0:n_A-1))
c$$$      allocate(U_d_(0:n_svd_l-1))
      svd_d_m = svd_d_max / 2.0
      svd_d_c = svd_d_m
      delta = dsqrt(delta_x**2 + delta_y**2)
      omega = atan2(delta_y,delta_x)
      if (delta.gt.svd_d_max .and. warning_flag) then
         write(6,'(A,F6.3,A,F6.3,A,F6.3,A)') 'Warning, delta '
     $        ,delta,'>',svd_d_max,'; ratio = ',delta/svd_d_max,
     $        ' in innerproduct_q_k_svdd_4'
      end if !if (delta.gt.svd_d_max .and. warning_flag) then
      svd_d(0) = (delta - svd_d_m)/svd_d_c
c$$$      do nl=0,n_svd_l-1
c$$$         call polyval_r8_reverse_0(n_svd_d,svd_U_d_(0+nl*n_svd_d),1
c$$$     $        ,svd_d(0),U_d_(nl))
c$$$      enddo                     !do nl=0,n_svd_l-1
c$$$      allocate(V_r_(0:n_svd_l*n_r-1))
      svd_r_m = svd_r_max / 2.0
      svd_r_c = svd_r_m
      do nr=0,n_r-1
         if (grid_p_(nr).gt.svd_r_max .and. warning_flag) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)')
     $           'Warning, grid_p_(nr) ',grid_p_(nr),'>',svd_r_max
     $           ,'; ratio = ',grid_p_(nr)/svd_r_max
     $           ,' in test_innerproduct_7_Z_S_q_4.f'
         end if
         svd_r(0) = (grid_p_(nr) - svd_r_m)/svd_r_c
c$$$         do nl=0,n_svd_l-1
c$$$            call polyval_r8_reverse_0(n_svd_r,svd_V_r_(0+nl*n_svd_r),1
c$$$     $           ,svd_r(0),V_r_(nl+nr*n_svd_l))
c$$$         enddo !do nl=0,n_svd_l-1         
      enddo !do nr=0,n_r-1
      allocate(C_w_(0:n_svd_l-1))
      ic=0
      do nr=0,n_r-1
         do nl=0,n_svd_l-1
            theta = svd_l_(nl)*(pi/2 - omega)
            C_w = cmplx( +cos(theta) , -sin(theta) )
            D_V_r = svd_polyval_V_r_(nl+nr*n_svd_l);
            D_U_d = svd_polyval_U_d_(nl)
            D_s = svd_s_(nl)
            C_w_(nl) = (D_U_d * D_s * D_V_r) * C_w
         enddo !do nl=0,n_svd_l-1
         do nw=0,n_w_(nr)-1
            Z_q_(ic) = cmplx( 0.0 , 0.0 )
            do nl=0,n_svd_l-1 
               I_l = svd_l_(nl)
               C_q = C_w_(nl)
               call periodize_i(nw+I_l,0,n_w_(nr),nwt)
               Z_q_(ic) = Z_q_(ic) + C_q*S_q_(ic - nw + nwt)
            enddo !do nl=0,n_svd_l-1 
            ic = ic + 1
         enddo !do nw=0,n_w_(nr)-1
      enddo !do nr=0,n_r-1
      call cp1_c16(n_A,Z_q_,M_q_)
      deallocate(C_w_)
c$$$      deallocate(V_r_)
c$$$      deallocate(U_d_)
      deallocate(Z_q_)
      end
