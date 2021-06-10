!> Doxygen comment: ;\n
!> calculates dotw(v1,v2) with real weight w_ ;\n
      real *8 function dotw_r8_f(n_v,v1_,v2_,w_)
      implicit none
      integer n_v
      real *8 v1_(0:n_v-1)
      real *8 v2_(0:n_v-1)
      real *8 w_(0:n_v-1)
      real *8 aa
      integer nv
      aa = 0.0d0
      do nv=0,n_v-1
         aa = aa + v1_(nv)*v2_(nv)*w_(nv)
      enddo !do nv=0,n_v-1
      dotw_r8_f = aa
      return
      end
