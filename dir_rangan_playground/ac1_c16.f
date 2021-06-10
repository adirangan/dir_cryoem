!> Doxygen comment: ;\n
!> affine addition: i.e., B = B + A*alpha + beta ;\n
      subroutine ac1_c16(n_x,alpha,beta,A_,B_)
      implicit none
      integer n_x
      complex *16 alpha,beta
      complex *16 A_(0:0),B_(0:0)
      integer nx,na,nb
      na = 0
      nb = 0
      do nx=0,n_x-1
         B_(nb) = B_(nb) + A_(na)*alpha + beta
         na = na + 1
         nb = nb + 1
      enddo
      end
