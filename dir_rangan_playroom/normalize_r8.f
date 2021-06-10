      subroutine normalize_r8(n_v,v_)
      implicit none
      integer n_v,nv
      real *8 v_(0:0)
      real *8 vn
      vn=0.0d0
      do nv=0,n_v-1
         vn = vn + v_(nv)**2
      enddo !do nv=0,n_v-1
      vn = dsqrt(vn)
      do nv=0,n_v-1
         v_(nv) = v_(nv)/vn;
      enddo !do nv=0,n_v-1
      end
      
