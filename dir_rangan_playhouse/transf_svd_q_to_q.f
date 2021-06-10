      subroutine transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_,n_svd_l
     $     ,svd_l_,svd_U_d_,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_
     $     ,delta_x,delta_y,M_q_)
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
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1),delta_x,delta_y
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      real *8 pi
      real *8 svd_d_max,n_svd_d_v,svd_r_max,n_svd_r_v
      integer n_svd_d_pre,n_svd_d_pos,n_svd_r_pre,n_svd_r_pos
      real *8 d_svd_d_pre,d_svd_d_pos,d_svd_r_pre,d_svd_r_pos
      real *8 alpha_d,beta_d,alpha_r,beta_r
      real *8 D_V_r,D_U_d,D_s
      integer I_l
      real *16 R_q,delta,omega,theta
      integer l,nr,ic,nw,nwt
      complex *16 C_q,C_w
      complex *16, allocatable :: C_w_(:)
      pi = 4.0*atan(1.0)
      allocate(Z_q_(0:n_A-1))
      delta = dsqrt(delta_x**2 + delta_y**2)
      svd_d_max = svd_d_(n_svd_d-1)
      if (delta.gt.svd_d_max*1.0625 .and. warning_flag) then
         write(6,'(A,F6.3,A,F6.3,A)') 'Warning, delta ',delta,'>'
     $        ,svd_d_max,' in transf_svd_q_to_q'
      end if
      svd_r_max = svd_r_(n_svd_r-1)
      R_q = grid_p_(n_r-1)
      if (2*pi*R_q.gt.svd_r_max*1.0625 .and. warning_flag) then
         write(6,'(A,F6.3,A,F6.3,A)') 'Warning, 2*pi*r ',2*pi *R_q,'>'
     $        ,svd_r_max,' in transf_svd_q_to_q'
      end if
      omega = atan2(delta_y,delta_x)
      n_svd_d_v = (n_svd_d-1)*delta/svd_d_max
      n_svd_d_pre = max(0,min(n_svd_d-1,floor(n_svd_d_v)))
      n_svd_d_pos = max(0,min(n_svd_d-1,ceiling(n_svd_d_v)))
      d_svd_d_pre = abs(n_svd_d_v - n_svd_d_pre)
      d_svd_d_pos = abs(n_svd_d_pos - n_svd_d_v)
      if (d_svd_d_pre+d_svd_d_pos.le.0.0d0) then
         alpha_d = 0.0d0
         beta_d = 1.0d0
      else
         alpha_d = d_svd_d_pre / (d_svd_d_pre + d_svd_d_pos)
         beta_d = d_svd_d_pos / (d_svd_d_pre + d_svd_d_pos)
      end if
      allocate(C_w_(0:n_svd_l-1))
      ic=0
      do nr=0,n_r-1
         R_q = 2*pi*grid_p_(nr)
         n_svd_r_v = (n_svd_r-1)*R_q/svd_r_max
         n_svd_r_pre = max(0,min(n_svd_r-1,floor(n_svd_r_v)))
         n_svd_r_pos = max(0,min(n_svd_r-1,ceiling(n_svd_r_v)))
         d_svd_r_pre = abs(n_svd_r_v - n_svd_r_pre)
         d_svd_r_pos = abs(n_svd_r_pos - n_svd_r_v)
         if (d_svd_r_pre+d_svd_r_pos.le.0.0d0) then
            alpha_r = 0.0d0
            beta_r = 1.0d0
         else
            alpha_r = d_svd_r_pre / (d_svd_r_pre + d_svd_r_pos)
            beta_r = d_svd_r_pos / (d_svd_r_pre + d_svd_r_pos)
         end if
         do l=0,n_svd_l-1
            theta = svd_l_(l)*(pi/2 - omega)
            C_w = cmplx( +cos(theta) , -sin(theta) )
            D_V_r = beta_r*svd_V_r_(n_svd_r_pre + l*n_svd_r) + alpha_r
     $           *svd_V_r_(n_svd_r_pos + l*n_svd_r)
            D_U_d = beta_d*svd_U_d_(n_svd_d_pre + l*n_svd_d) + alpha_d
     $           *svd_U_d_(n_svd_d_pos + l*n_svd_d)
            D_s = svd_s_(l)
            C_w_(l) = (D_U_d * D_s * D_V_r) * C_w
         enddo
         do nw=0,n_w_(nr)-1
            Z_q_(ic) = cmplx( 0.0 , 0.0 )
            do l=0,n_svd_l-1
               I_l = svd_l_(l)
               C_q = C_w_(l)
               call periodize_i(nw+I_l,0,n_w_(nr),nwt)
               Z_q_(ic) = Z_q_(ic) + C_q*S_q_(ic - nw + nwt)
            enddo
            ic = ic + 1
         enddo
      enddo
      call cp1_c16(n_A,Z_q_,M_q_)
      deallocate(C_w_)
      deallocate(Z_q_)
      end
