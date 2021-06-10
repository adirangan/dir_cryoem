      subroutine distance_r8(n_v,vA_,vB_,dd)
      integer *4 n_v
      real *8 vA_(0:n_v-1)
      real *8 vB_(0:n_v-1)
      real *8 dd
      integer *4 nv
      dd = 0.0d0
      do nv=0,n_v-1
         dd = dd + (vA_(nv) - vB_(nv))**2
      enddo !do nv=0,n_v-1
      dd = dsqrt(dd)
      end
