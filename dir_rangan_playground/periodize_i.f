!> Doxygen comment: ;\n
!> periodize integer A within the limits [min,max) ;\n
!> Note that if input=max, the output will equal min. ;\n
      subroutine periodize_i(A,min,max,B)
      implicit none
      integer A,min,max,B
      B = A
      do while (B.lt.min)
         B = B + (max-min)
      end do
      do while (B.ge.max)
         B = B - (max-min)
      end do
      end
