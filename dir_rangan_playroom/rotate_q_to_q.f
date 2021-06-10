      subroutine rotate_q_to_q(n_r,n_w_,n_A,S_q_,gamma,M_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as S_q_
c$$$      Multiplication performed in place
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma,nt_d
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
c$$$      complex *16, allocatable :: Z_q_(:)
      real *8 pi
      integer nr,ic,nw,nt_pre,nt_pos
      real *8 d_pre,d_pos,A,B,alpha,beta,theta
      complex *16 w1,wq1,wt,wqt,C_q
c$$$      allocate(Z_q_(0:n_A-1))
      pi = 4.0*atan(1.0)
      ic=0
      do nr=0,n_r-1
         nt_d = gamma*n_w_(nr)/(2*pi)
         nt_pre = floor(nt_d)
         nt_pos = ceiling(nt_d)
         d_pre = dabs(nt_d-nt_pre)
         d_pos = dabs(nt_pos-nt_d)
         A = d_pre
         B = d_pos
         if (A+B.le.0.0) then
            alpha = 0.0d0
            beta = 1.0d0
         else
            alpha = A/(A+B)
            beta = B/(A+B)
         end if
         theta = 2*pi/n_w_(nr)
         w1 = cmplx(+cos(theta),-sin(theta))
         wt = cmplx(+cos(theta*nt_pre),-sin(theta*nt_pre))
         wq1 = cmplx(1.0,0.0)
         wqt = cmplx(1.0,0.0)
         do nw=0,n_w_(nr)-1
            C_q = wqt*(beta + alpha*wq1)
c$$$            Z_q_(ic) = S_q_(ic)*C_q
            M_q_(ic) = S_q_(ic)*C_q
            wq1 = wq1 * w1
            wqt = wqt * wt
            ic = ic + 1
         enddo
      enddo
c$$$      call cp1_c16(n_A,Z_q_,M_q_)
c$$$      deallocate(Z_q_)
      end
