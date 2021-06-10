!> Doxygen comment: ;\n
!> conjugate multiply two complex *16 arrays ;\n
!> C = A*dconjg(B) ;\n
      subroutine xc1_c16(n_x,A_,B_,C_)
      implicit none
      integer n_x
      complex *16 A_(0:0),B_(0:0),C_(0:0)
      integer nx,na,nb,nc
      na = 0
      nb = 0
      nc = 0
      do nx=0,n_x-1
         C_(nc) = A_(na)*dconjg(B_(nb))
         na = na + 1
         nb = nb + 1
         nc = nc + 1
      enddo
      end
