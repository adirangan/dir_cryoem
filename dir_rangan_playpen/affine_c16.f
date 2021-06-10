      subroutine affine_c16(n_x,alpha,beta,n_A,A_,A_start,A_stride,n_B
     $     ,B_,B_start,B_stride)
      implicit none
      integer n_x,n_A,n_B
      complex *16 alpha,beta
      complex *16 A_(0:n_A-1),B_(0:n_B-1)
      integer A_start,A_stride,B_start,B_stride
      integer nx,na,nb
      na = A_start
      nb = B_start
      do nx=0,n_x-1
         B_(nb) = A_(na)*alpha + beta
         na = na + A_stride
         nb = nb + B_stride
      enddo
      end
