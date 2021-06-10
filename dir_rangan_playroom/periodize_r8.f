      subroutine periodize_r8(A,min,max,B)
      implicit none
      real *8 A,min,max,B
      B = A
      do while (B.lt.min)
         B = B + (max-min)
      end do
      do while (B.ge.max)
         B = B - (max-min)
      end do
      end
