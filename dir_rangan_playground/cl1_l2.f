!> Doxygen comment: ;\n
!> clears logical array ;\n
      subroutine cl1_l2(n_x,A_)
      implicit none
      integer n_x
      logical A_(0:0)
      integer nx
      nx = 0
      do nx=0,n_x-1
         A_(nx) = .false.
      enddo
      end
