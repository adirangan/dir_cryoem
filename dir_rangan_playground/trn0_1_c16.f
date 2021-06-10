!> Doxygen comment: ;\n
!> transpose complex *16 array ;\n
!> B = transpose(A) ;\n
!> assumes workspace C is preallocated ;\n
      subroutine trn0_1_c16(n_r,n_c,A_,B_,C_)
c$$$      C_ is workspace of same size as A_ and B_. ;
      implicit none
      integer n_r,n_c
      complex *16 A_(0:n_r*n_c-1),B_(0:n_r*n_c-1),C_(0:n_r*n_c-1)
      integer nr,nc
      do nc=0,n_c-1
         do nr=0,n_r-1
            C_(nc+nr*n_c) = A_(nr+nc*n_r)
         enddo
      enddo
      call cp1_c16(n_r*n_c,C_,B_)
      end
