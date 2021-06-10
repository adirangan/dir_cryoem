!> Doxygen comment: ;\n
!> calculates innerproduct in ?_p (i.e., polar) coordinates. ;\n
!> assumes that M_k_p_ is the same size and dimensions as T_k_p_. ;\n
!> assumes quasi-uniform polar-grid. ;\n
!> Upgraded to use radial quadrature weight. ;\n
      subroutine innerproduct_p_quad(
     $     n_k_p_r
     $     ,grid_k_p_r_
     $     ,weight_k_p_r_
     $     ,n_w_
     $     ,n_w_sum
     $     ,T_k_p_
     $     ,M_k_p_
     $     ,C
     $     )
c$$$      assumes that M_k_p_ is the same size and dimensions as T_k_p_
c$$$      assumes quasi-uniform polar-grid
c$$$      Upgraded to use radial quadrature weight. ;\n
      implicit none
      integer n_k_p_r,n_w_(0:n_k_p_r-1),n_w_sum
      real *8 grid_k_p_r_(0:n_k_p_r-1)
      real *8 weight_k_p_r_(0:n_k_p_r-1)
      complex *16 T_k_p_(0:n_w_sum-1),M_k_p_(0:n_w_sum-1)
      complex *16 C,C_tmp
      real *8 pi
      integer ic,nr,nw
      real *8 dw,dA,dAn
      pi = 4.0d0*datan(1.0d0)
      C = dcmplx(0.0d0,0.0d0)
      ic = 0
      do nr=0,n_k_p_r-1
         dw = 2.0d0*pi/(1.0d0*max(1,n_w_(nr)))
         dA = weight_k_p_r_(nr)
         dAn = dA*dw
         C_tmp = dcmplx(0.0d0,0.0d0)
         do nw=0,n_w_(nr)-1
            C_tmp = C_tmp + dconjg(T_k_p_(ic))*M_k_p_(ic)
            ic = ic + 1
         enddo
         C_tmp = C_tmp * dAn
         C = C + C_tmp
      enddo
      end
