!> Doxygen comment: ;\n
!> copies real *8 array ; \n
      subroutine cp1_r8(n_x,A_,B_)
      implicit none
      integer n_x
      real *8 A_(0:0),B_(0:0)
      integer nx,na,nb
      na = 0
      nb = 0
      do nx=0,n_x-1
         B_(nb) = A_(na)
         na = na + 1
         nb = nb + 1
      enddo
      end
