!> Doxygen comment: ;\n
!> clears real *8 array ;\n
      subroutine cl1_r8(n_x,A_)
      implicit none
      integer n_x
      real *8 A_(0:0)
      integer nx
      nx = 0
      do nx=0,n_x-1
         A_(nx) = 0.0d0
      enddo
      end
