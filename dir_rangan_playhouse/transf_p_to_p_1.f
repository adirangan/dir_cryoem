      subroutine transf_p_to_p_1(grid_p,n_w,S_p_1_,delta_x,delta_y
     $     ,M_p_1_)
c$$$      Assumes that M_p_1_ is the same size and dimensions as S_p_1_
c$$$      Multiplication performed in place
      implicit none
      integer n_w
      real *8 grid_p,delta_x,delta_y
      complex *16 S_p_1_(0:n_w-1),M_p_1_(0:n_w-1)
      real *8 pi
      integer ic,nw,nt_pre,nt_pos
      real *8 R_c,W_c,X_c,Y_c,L_c
      complex *16 C_c
      pi = 4.0*atan(1.0)
      ic=0
      R_c = grid_p
      do nw=0,n_w-1
         W_c = 0.0 + nw*(2*pi)/(n_w)
         X_c = R_c*cos(W_c)
         Y_c = R_c*sin(W_c)
         L_c = (X_c * delta_x) + (Y_c * delta_y)
         C_c = cmplx(+cos(2*pi*L_c),-sin(2*pi*L_c))
         M_p_1_(ic) = C_c*S_p_1_(ic)
         ic = ic + 1
      enddo
      end
