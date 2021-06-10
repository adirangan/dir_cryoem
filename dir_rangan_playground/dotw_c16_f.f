!> Doxygen comment: ;\n
!> calculates dotw(v1,v2) with real weight w_ ;\n
      complex *16 function dotw_c16_f(n_v,v1_,v2_,w_)
      implicit none
      integer n_v
      complex *16 v1_(0:n_v-1)
      complex *16 v2_(0:n_v-1)
      real *8 w_(0:n_v-1)
      complex *16 aa
      integer nv
      aa = dcmplx(0.0d0,0.0d0)
      do nv=0,n_v-1
         aa = aa + dconjg(v1_(nv))*v2_(nv)*w_(nv)
      enddo !do nv=0,n_v-1
      dotw_c16_f = aa
      return
      end
