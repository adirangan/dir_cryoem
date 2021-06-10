!> Doxygen comment: ;\n
!> Calculates dot-product aa of real *8 vectors v1_ and v2_ ;\n
      subroutine dot_r8(v1,v2,aa)
      real *8 v1(0:3-1)
      real *8 v2(0:3-1)
      real *8 aa
      aa = v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2)
      end
