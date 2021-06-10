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
