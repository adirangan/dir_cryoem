      subroutine trn0_r8(n_r,n_c,A_,B_)
      implicit none
      integer n_r,n_c
      real *8 A_(0:n_r*n_c-1),B_(0:n_r*n_c-1)
      integer nr,nc
      real *8 C_(0:n_r*n_c-1)
      do nc=0,n_c-1
         do nr=0,n_r-1
            C_(nc+nr*n_c) = A_(nr+nc*n_r)
         enddo
      enddo
      call cp1_r8(n_r*n_c,C_,B_)
      end
