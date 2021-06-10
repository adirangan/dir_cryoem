!> Doxygen comment: ;\n
!> clears integer *4 array ;\n
      subroutine cl1_i4(n_x,A_)
      implicit none
      integer n_x
      integer *4 A_(0:0)
      integer nx
      nx = 0
      do nx=0,n_x-1
         A_(nx) = 0
      enddo
      end
