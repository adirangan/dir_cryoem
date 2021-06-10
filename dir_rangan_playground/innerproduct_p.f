!> Doxygen comment: ;\n
!> calculates innerproduct in ?_p (i.e., polar) coordinates. ;\n
!> assumes that M_p_ is the same size and dimensions as T_p_. ;\n
!> assumes quasi-uniform polar-grid. ;\n
      subroutine innerproduct_p(n_r,grid_p_,n_w_,n_A,T_p_,M_p_,C_p)
c$$$      assumes that M_p_ is the same size and dimensions as T_p_
c$$$      assumes quasi-uniform polar-grid
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 T_p_(0:n_A-1),M_p_(0:n_A-1),C_p,C_tmp
      real *8 pi
      integer ic,nr,nw
      real *8 R_pos,R_pre,dr,dw,dA
      pi = 4.0d0*datan(1.0d0)
      C_p = dcmplx(0.0d0,0.0d0)
      ic = 0
      do nr=0,n_r-1
         if (nr.gt.0) then
            R_pre = 0.5*(grid_p_(nr-1) + grid_p_(nr))
         else
            R_pre = grid_p_(0)
         end if
         if (nr.lt.n_r-1) then
            R_pos = 0.5*(grid_p_(nr+1) + grid_p_(nr))
         else
            R_pos = grid_p_(n_r-1)
         end if
         dr = R_pos - R_pre
c$$$         We set the zero-mode to zero
         if (grid_p_(nr).le.0.0d0) then
            dr = 0.0d0
         end if
         dw = 2.0d0*pi/(1.0d0*max(1,n_w_(nr)))
         dA = (R_pre*dr + (dr**2)/2)*dw
         C_tmp = dcmplx(0.0d0,0.0d0)
         do nw=0,n_w_(nr)-1
            C_tmp = C_tmp + dconjg(T_p_(ic))*M_p_(ic)
            ic = ic + 1
         enddo
         C_tmp = C_tmp * dA
         C_p = C_p + C_tmp
      enddo
      end
