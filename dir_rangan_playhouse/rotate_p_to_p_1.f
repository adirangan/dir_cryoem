      subroutine rotate_p_to_p_1(n_w,S_p_1_,gamma,M_p_1_)
c$$$  Assumes that M_p_1_ is the same size and dimensions as S_p_1_
      implicit none
      integer n_w
      real *8 gamma
      complex *16 S_p_1_(0:n_w-1),M_p_1_(0:n_w-1)
      complex *16, allocatable :: Z_p_1_(:)
      real *8 pi
      integer ic,nw,nt_pre,nt_pos
      real *8 W_c,W_pre,W_pos,A,B,alpha,beta
      complex *16 C_p
      allocate(Z_p_1_(0:n_w-1))
      pi = 4.0*atan(1.0)
      ic=0
      do nw=0,n_w-1
         W_c = 0.0 + nw*(2*pi)/(n_w) - gamma
         call periodize_r8(W_c,0.0d0,2.0*pi,W_c)
         nt_pre = floor(W_c*n_w/(2*pi))
         nt_pos = ceiling(W_c*n_w/(2*pi))
         if (nt_pos.ge.n_w) then
            nt_pos = nt_pos - n_w
         end if
         W_pre = 0.0 + nt_pre*(2*pi)/n_w
         W_pos = 0.0 + nt_pos*(2*pi)/n_w
         A = dabs(W_c-W_pre)
         B = dabs(W_pos-W_c)
         if (A+B.le.0.0) then
            C_p = S_p_1_(nt_pre)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            C_p = beta*S_p_1_(nt_pre) + alpha*S_p_1_(nt_pos)
         end if
         Z_p_1_(ic) = C_p
         ic = ic + 1
      enddo
      call cp1_c16(n_w,Z_p_1_,M_p_1_)
      deallocate(Z_p_1_)
      end
