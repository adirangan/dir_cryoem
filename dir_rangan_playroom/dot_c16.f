      subroutine dot_c16(n_v,v1_,v2_,aa)
      implicit none
      integer n_v
      complex *16 v1_(0:n_v-1)
      complex *16 v2_(0:n_v-1)
      complex *16 aa
      integer nv
      aa = cmplx(0.0d0,0.0d0)
      do nv=0,n_v-1
         aa = aa + dconjg(v1_(nv))*v2_(nv)
      enddo !do nv=0,n_v-1
      end
