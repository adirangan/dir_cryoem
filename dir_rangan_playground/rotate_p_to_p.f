!> Doxygen comment: ;\n
!> Applies rotation in ?_p (i.e., polar) coordinates. ;\n
!> This uses linear interpolation. ;\n
!> Assumes that M_p_ is the same size and dimensions as S_p_. ;\n
      subroutine rotate_p_to_p(n_r,n_w_,n_A,S_p_,gamma,M_p_)
c$$$      Assumes that M_p_ is the same size and dimensions as S_p_
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma
      complex *16 S_p_(0:n_A-1),M_p_(0:n_A-1)
      complex *16, allocatable :: Z_p_(:)
      real *8 pi
      integer nr,ic,icstart,nw,nt_pre,nt_pos
      real *8 W_c,W_pre,W_pos,A,B,alpha,beta
      complex *16 C_p
      allocate(Z_p_(0:n_A-1))
      pi = 4.0d0*datan(1.0d0)
      ic=0
      icstart=0
      do nr=0,n_r-1         
         do nw=0,n_w_(nr)-1
            W_c = 0.0 + nw*(2.0d0*pi)/(n_w_(nr)) - gamma
            call periodize_r8(W_c,0.0d0,2.0*pi,W_c)
            nt_pre = floor(W_c*n_w_(nr)/(2.0d0*pi))
            nt_pos = ceiling(W_c*n_w_(nr)/(2.0d0*pi))
            if (nt_pos.ge.n_w_(nr)) then
               nt_pos = nt_pos - n_w_(nr)
            end if
            W_pre = 0.0 + nt_pre*(2.0d0*pi)/n_w_(nr)
            W_pos = 0.0 + nt_pos*(2.0d0*pi)/n_w_(nr)
            A = dabs(W_c-W_pre)
            B = dabs(W_pos-W_c)
            if (A+B.le.0.0) then
               C_p = S_p_(icstart + nt_pre)
            else
               alpha = A/(A+B)
               beta = B/(A+B)
               C_p = beta*S_p_(icstart + nt_pre) + alpha*S_p_(icstart +
     $              nt_pos)
            end if
            Z_p_(ic) = C_p
            ic = ic + 1
         enddo
         icstart = icstart + n_w_(nr)
      enddo
      call cp1_c16(n_A,Z_p_,M_p_)
      deallocate(Z_p_)
      end
