!> Doxygen comment: ;\n
!> clears integer *8 array ;\n
!> sets final+1 entry for memory checking ;\n
      subroutine cs1_i8(n_x,A_)
      implicit none
      integer n_x
      integer *8 A_(0:0)
      integer nx
      nx = 0
      do nx=0,n_x-1
         A_(nx) = 0
      enddo
      A_(n_x) = -1
      end
