!> Doxygen comment: ;\n
!> Applies rotation in ?_p (i.e., polar) coordinates. ;\n
!> This uses lagrange interpolation with n_q=6 nodes. ;\n
!> Assumes that M_p_ is the same size and dimensions as S_p_. ;\n
      subroutine rotate_p2p_l6(n_r,n_w_,n_A,S_p_,gamma,M_p_)
c$$$      Assumes that M_p_ is the same size and dimensions as S_p_
c$$$      uses n_q==6 nodes.
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma
      complex *16 S_p_(0:n_A-1),M_p_(0:n_A-1)
      complex *16, allocatable :: Z_p_(:)
      real *8 pi,pi2i
      integer nr,ic,icstart,nw,nt_pre,nt_pos
      real *8 W_c,W_pre,W_pos,A,B,alpha,beta
      integer n_q,nq
      parameter(n_q=6)
      integer *4 nt_(0:n_q-1)
      complex *16 F_(0:n_q-1)
      real *8 q_(0:n_q-1)
      real *8 x_(0:n_q-1)
      real *8 d_(0:n_q-1)
      real *8 e_(0:n_q-1)
      real *8 x_tmp
      complex *16 C_p
      allocate(Z_p_(0:n_A-1))
      x_(0) = 0.0d0
      x_(1) = 1.0d0
      x_(2) = 2.0d0
      x_(3) = 3.0d0
      x_(4) = 4.0d0
      x_(5) = 5.0d0
      d_(0) = ((x_(0)-x_(1))*(x_(0)-x_(2)) *(x_(0)-x_(3))*(x_(0)-x_(4))
     $     *(x_(0)-x_(5)))
      d_(1) = ((x_(1)-x_(0))*(x_(1)-x_(2)) *(x_(1)-x_(3))*(x_(1)-x_(4))
     $     *(x_(1)-x_(5)))
      d_(2) = ((x_(2)-x_(0))*(x_(2)-x_(1)) *(x_(2)-x_(3))*(x_(2)-x_(4))
     $     *(x_(2)-x_(5)))
      d_(3) = ((x_(3)-x_(0))*(x_(3)-x_(1)) *(x_(3)-x_(2))*(x_(3)-x_(4))
     $     *(x_(3)-x_(5)))
      d_(4) = ((x_(4)-x_(0))*(x_(4)-x_(1)) *(x_(4)-x_(2))*(x_(4)-x_(3))
     $     *(x_(4)-x_(5)))
      d_(5) = ((x_(5)-x_(0))*(x_(5)-x_(1)) *(x_(5)-x_(2))*(x_(5)-x_(3))
     $     *(x_(5)-x_(4)))
      e_(0) = 1.0d0/d_(0)
      e_(1) = 1.0d0/d_(1)
      e_(2) = 1.0d0/d_(2)
      e_(3) = 1.0d0/d_(3)
      e_(4) = 1.0d0/d_(4)
      e_(5) = 1.0d0/d_(5)
      pi = 4.0d0*datan(1.0d0)
      pi2i = 1.0d0/(2.0d0*pi)
      ic=0
      icstart=0
      do nr=0,n_r-1         
         do nw=0,n_w_(nr)-1
            W_c = 0.0 + nw*(2.0d0*pi)/(n_w_(nr)) - gamma
            call periodize_r8(W_c,0.0d0,2.0*pi,W_c)
            nt_pre = floor(W_c*n_w_(nr)*pi2i)
            nt_(0) = nt_pre - ceiling(n_q*0.5d0) + 1 + 0
            nt_(1) = nt_pre - ceiling(n_q*0.5d0) + 1 + 1
            nt_(2) = nt_pre - ceiling(n_q*0.5d0) + 1 + 2
            nt_(3) = nt_pre - ceiling(n_q*0.5d0) + 1 + 3
            nt_(4) = nt_pre - ceiling(n_q*0.5d0) + 1 + 4
            nt_(5) = nt_pre - ceiling(n_q*0.5d0) + 1 + 5
            call periodize_i(nt_(0),0,n_w_(nr))
            call periodize_i(nt_(1),0,n_w_(nr))
            call periodize_i(nt_(2),0,n_w_(nr))
            call periodize_i(nt_(3),0,n_w_(nr))
            call periodize_i(nt_(4),0,n_w_(nr))
            call periodize_i(nt_(5),0,n_w_(nr))
            F_(0) = S_p_(icstart + nt_(0))
            F_(1) = S_p_(icstart + nt_(1))
            F_(2) = S_p_(icstart + nt_(2))
            F_(3) = S_p_(icstart + nt_(3))
            F_(4) = S_p_(icstart + nt_(4))
            F_(5) = S_p_(icstart + nt_(5))
            x_tmp = W_c*n_w_(nr)*pi2i - (nt_pre - ceiling(n_q*0.5d0) + 1
     $           + 0)
            q_(0) = ((x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp-x_(3))*(x_tmp
     $           -x_(4))*(x_tmp-x_(5)))*e_(0)
            q_(1) = ((x_tmp-x_(0))*(x_tmp-x_(2))*(x_tmp-x_(3))*(x_tmp
     $           -x_(4))*(x_tmp-x_(5)))*e_(1)
            q_(2) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(3))*(x_tmp
     $           -x_(4))*(x_tmp-x_(5)))*e_(2)
            q_(3) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $           -x_(4))*(x_tmp-x_(5)))*e_(3)
            q_(4) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $           -x_(3))*(x_tmp-x_(5)))*e_(4)
            q_(5) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $           -x_(3))*(x_tmp-x_(4)))*e_(5)
            Z_p_(ic) = F_(0)*q_(0) + F_(1)*q_(1) + F_(2)*q_(2) + F_(3)
     $           *q_(3) + F_(4)*q_(4) + F_(5)*q_(5)
            ic = ic + 1
         enddo
         icstart = icstart + n_w_(nr)
      enddo
      call cp1_c16(n_A,Z_p_,M_p_)
      deallocate(Z_p_)
      end
