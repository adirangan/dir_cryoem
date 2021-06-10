      subroutine periodize_r8_(n_A,A_,min,max,B_)
      implicit none
      integer *4 n_A
      real *8 A_(0:0),min,max,B_(0:0)
      integer *4 na
      do na=0,n_A-1
         call periodize_r8(A_(na),min,max,B_(na))
      enddo !do na=0,n_A-1
      end !subroutine
