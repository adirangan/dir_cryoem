!> Doxygen comment: ;\n
!> Simple gram-schmidt orthonormalization of columns of complex *16 array A_. ;\n
      subroutine gramschmidt_c16(n_r,n_c,A_)
      implicit none
      integer n_r,n_c
      complex *16 A_(0:n_r*n_c-1)
      complex *16 aa
      integer nr,nc,nd
      do nc=0,n_c-1
         call normalize_c16(n_r,A_(0+nc*n_r))
         do nd=nc+1,n_c-1
            call dot_c16(n_r,A_(0+nc*n_r),A_(0+nd*n_r),aa)
            do nr=0,n_r-1
               A_(nr+nd*n_r) = A_(nr+nd*n_r) - aa*A_(nr+nc*n_r)
            enddo !do nr=0,n_r-1
         enddo !do nd=nc+1,n_c-1
      enddo !do nc=0,n_c-1
      end
