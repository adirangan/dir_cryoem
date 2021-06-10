!> Doxygen comment: ;\n
!> Applies translation in k_p (i.e., fourier-space polar) coordinates. ;\n
!> This amounts to a multiplication by a plane-wave. ;\n
!> Assumes that M_k_p_ is the same size and dimensions as S_k_p_. ;\n
      subroutine transf_p_to_p(
     $     n_k_p_r
     $     ,grid_k_p_r_
     $     ,n_w_
     $     ,n_w_sum
     $     ,S_k_p_
     $     ,delta_x
     $     ,delta_y
     $     ,M_k_p_
     $     )
c$$$      Assumes that M_k_p_ is the same size and dimensions as S_k_p_
c$$$      Multiplication performed in place
      implicit none
      integer n_k_p_r,n_w_(0:n_k_p_r-1),n_w_sum
      real *8 grid_k_p_r_(0:n_k_p_r-1),delta_x,delta_y
      complex *16 S_k_p_(0:n_w_sum-1),M_k_p_(0:n_w_sum-1)
c$$$      complex *16, allocatable :: Z_p_(:)
      real *8 pi
      integer nr,ic,nw,nt_pre,nt_pos
      real *8 R_c,W_c,X_c,Y_c,L_c
      complex *16 C_c
c$$$      allocate(Z_p_(0:n_w_sum-1))
      pi = 4.0d0*datan(1.0d0)
      ic=0
      do nr=0,n_k_p_r-1
         R_c = grid_k_p_r_(nr)
         do nw=0,n_w_(nr)-1
            W_c = 0.0d0 + nw*(2.0d0*pi)/(n_w_(nr))
            X_c = R_c*dcos(W_c)
            Y_c = R_c*dsin(W_c)
            L_c = (X_c * delta_x) + (Y_c * delta_y)
            C_c = dcmplx(+dcos(2.0d0*pi*L_c),-dsin(2.0d0*pi*L_c))
c$$$            Z_p_(ic) = C_c*S_k_p_(ic)
            M_k_p_(ic) = C_c*S_k_p_(ic)
            ic = ic + 1
         enddo
      enddo
c$$$      call copy_c16(n_w_sum,n_w_sum,Z_p_,0,1,n_w_sum,M_k_p_,0,1)
c$$$      deallocate(Z_p_)
      end
