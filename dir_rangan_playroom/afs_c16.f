      subroutine afs_c16(n_x,alpha,beta,A_,A_stride,B_,B_stride)
      implicit none
      integer n_x
      complex *16 alpha,beta
      complex *16 A_(0:0),B_(0:0)
      integer A_stride,B_stride
      integer nx,na,nb
      na = 0
      nb = 0
      do nx=0,n_x-1
         B_(nb) = A_(na)*alpha + beta
         na = na + A_stride
         nb = nb + B_stride
      enddo
      end
