      subroutine cross_r8(v1,v2,v3)
      real *8 v1(0:3-1)
      real *8 v2(0:3-1)
      real *8 v3(0:3-1)
      v3(0) = v1(1)*v2(2) - v2(1)*v1(2)
      v3(1) = v1(2)*v2(0) - v2(2)*v1(0)
      v3(2) = v1(0)*v2(1) - v2(0)*v1(1)
      end
