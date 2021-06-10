!> Doxygen comment: ;\n
!> simple print function: ;\n
      subroutine ping(n)
      implicit none
      integer n
      write(6,'(A,I0)') ' % ping: ' , n
      end
