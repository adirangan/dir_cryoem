      subroutine cl1_c16(n_x,A_)
      implicit none
      integer n_x
      complex *16 A_(0:0)
      integer nx
      nx = 0
      do nx=0,n_x-1
         A_(nx) = cmplx(0.0d0,0.0d0)
      enddo
      end
