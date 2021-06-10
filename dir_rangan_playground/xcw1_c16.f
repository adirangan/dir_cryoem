!> Doxygen comment: ;\n
!> conjugate multiply two complex *16 arrays, alongside a real weight. ;\n
!> C = A*dconjg(B) ;\n
      subroutine xcw1_c16(n_x,A_,B_,w_,C_)
      implicit none
      integer n_x
      complex *16 A_(0:0),B_(0:0),C_(0:0)
      real *8 w_(0:0)
      integer nx,na,nb,nc,nw
      na = 0
      nb = 0
      nc = 0
      nw = 0
      do nx=0,n_x-1
         C_(nc) = A_(na)*dconjg(B_(nb))*w_(nw)
         na = na + 1
         nb = nb + 1
         nc = nc + 1
         nw = nw + 1
      enddo
      end
