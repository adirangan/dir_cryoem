      subroutine transf_p_to_p(n_r,grid_p_,n_w_,n_A,S_p_,delta_x,delta_y
     $     ,M_p_)
c$$$      Assumes that M_p_ is the same size and dimensions as S_p_
c$$$      Multiplication performed in place
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1),delta_x,delta_y
      complex *16 S_p_(0:n_A-1),M_p_(0:n_A-1)
c$$$      complex *16, allocatable :: Z_p_(:)
      real *8 pi
      integer nr,ic,nw,nt_pre,nt_pos
      real *8 R_c,W_c,X_c,Y_c,L_c
      complex *16 C_c
c$$$      allocate(Z_p_(0:n_A-1))
      pi = 4.0*atan(1.0)
      ic=0
      do nr=0,n_r-1
         R_c = grid_p_(nr)
         do nw=0,n_w_(nr)-1
            W_c = 0.0 + nw*(2*pi)/(n_w_(nr))
            X_c = R_c*cos(W_c)
            Y_c = R_c*sin(W_c)
            L_c = (X_c * delta_x) + (Y_c * delta_y)
            C_c = cmplx(+cos(2*pi*L_c),-sin(2*pi*L_c))
c$$$            Z_p_(ic) = C_c*S_p_(ic)
            M_p_(ic) = C_c*S_p_(ic)
            ic = ic + 1
         enddo
      enddo
c$$$      call copy_c16(n_A,n_A,Z_p_,0,1,n_A,M_p_,0,1)
c$$$      deallocate(Z_p_)
      end
