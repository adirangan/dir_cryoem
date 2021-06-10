!> Doxygen comment: ;\n
!> Applies translation in k_p (i.e., fourier-space polar) coordinates. ;\n
!> This amounts to a multiplication by a plane-wave. ;\n
!> Assumes that M_p_ is the same size and dimensions as S_p_. ;\n
!> This only applies to a single ring of S_p_. ;\n
      subroutine transf_p_to_p_single(grid_p,n_w,S_p_1_,delta_x,delta_y
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
      pi = 4.0d0*datan(1.0d0)
      ic=0
      R_c = grid_p
      do nw=0,n_w-1
         W_c = 0.0 + nw*(2.0d0*pi)/(n_w)
         X_c = R_c*dcos(W_c)
         Y_c = R_c*dsin(W_c)
         L_c = (X_c * delta_x) + (Y_c * delta_y)
         C_c = dcmplx(+dcos(2.0d0*pi*L_c),-dsin(2.0d0*pi*L_c))
         M_p_1_(ic) = C_c*S_p_1_(ic)
         ic = ic + 1
      enddo
      end
