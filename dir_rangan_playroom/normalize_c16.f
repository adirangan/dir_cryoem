      subroutine normalize_c16(n_v,v_)
      implicit none
      integer n_v,nv
      complex *16 v_(0:0)
      real *8 vn
      vn=0.0d0
      do nv=0,n_v-1
         vn = vn + dreal(dconjg(v_(nv))*v_(nv))
      enddo !do nv=0,n_v-1
      vn = dsqrt(vn)
      vn = max(1.0d-15,vn)
      do nv=0,n_v-1
         v_(nv) = v_(nv)/vn;
      enddo !do nv=0,n_v-1
      end
      
