!> Doxygen comment: ;\n
!> sum of real *8 array. ;\n
      real *8 function sum_r8_f(n_x,x_)
      integer *4 n_x,nx
      real *8 x_(0:n_x-1),s
      if (n_x.le.0) then
         s=0
      else
         s=0
         do nx=0,n_x-1
            s = s + x_(nx)
         enddo
      end if
      sum_r8_f = s
      return
      end
