      subroutine transpose_c16(n_r,n_c,A_,C_)
      implicit none
      integer n_r,n_c
      complex *16 A_(0:0)
      complex *16 C_(0:0)
      integer nr,nc
      complex *16, allocatable :: B_(:)
      allocate(B_(0:n_r*n_c-1))
      do nc=0,n_c-1
      do nr=0,n_r-1
         B_(nc+nr*n_c) = A_(nr+nc*n_r)
      enddo !do nr=0,n_r-1
      enddo !do nc=0,n_c-1
      call cp1_c16(n_r*n_c,B_,C_)
      deallocate(B_)
      end !subroutine
