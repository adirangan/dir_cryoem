      subroutine transf_taylor_q_to_q(l_max,n_J,f_,b_,p_,n_r,grid_p_
     $     ,n_w_,n_A,S_q_,delta_x,delta_y,M_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as S_q_
c$$$      Currently brute force convolution using taylor-expansion
c$$$      defined via l_max, n_J, f_, b_, p_
      implicit none
      integer l_max,n_J,b_(0:n_J-1),p_(0:n_J-1)
      complex *16 f_(0:n_J-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1),delta_x,delta_y
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      complex *16, allocatable :: g_(:)
      real *8 pi
      real *16 R_q,z,zl,delta,omega,theta
      integer l,j,jc,nr,ic,nw,nwt
      complex *16 C_q
      if (n_J.ne.(l_max+1)*(l_max+2)/2) then
         write(6,'(A)')
     $        'Warning, n_J does not match l_max in transf_q_to_q'
      end if
      pi = 4.0*atan(1.0)
      delta = dsqrt(delta_x**2 + delta_y**2)
      omega = atan2(delta_y,delta_x)
      allocate(g_(0:n_J-1))
      jc=0
      do l=0,l_max
         do j=0,l
            theta = p_(jc)*omega
            g_(jc) = f_(jc)*b_(jc)*cmplx(cos(theta),-sin(theta))
            jc = jc + 1
         enddo
      enddo
      allocate(Z_q_(0:n_A-1))
      ic=0
      do nr=0,n_r-1
         R_q = grid_p_(nr)
         z = 2*pi*R_q*delta
         do nw=0,n_w_(nr)-1
            Z_q_(ic) = cmplx( 0.0 , 0.0 )
            zl = cmplx( 1.0 , 0.0 )
            jc=0
            do l=0,l_max
               if (l.gt.0) then
                  zl = zl * z
               end if
               do j=0,l
                  C_q = g_(jc)*zl
                  call periodize_i(nw-p_(jc),0,n_w_(nr),nwt)
                  Z_q_(ic) = Z_q_(ic) + C_q*S_q_(ic - nw + nwt)
                  jc = jc + 1
               enddo
            enddo
            ic = ic + 1
         enddo
      enddo
      call cp1_c16(n_A,Z_q_,M_q_)
      deallocate(g_)
      deallocate(Z_q_)
      end
