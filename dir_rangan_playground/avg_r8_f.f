!> Doxygen comment: ;\n
!> calculates average of real *8 array x_ ;\n
      real *8 function avg_r8_f(n_x,x_)
      integer n_x,nx
      real *8 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0.0d0
      else
         m=x_(0)
         do nx=1,n_x-1
            m = m + x_(nx)
         enddo
         m = m/n_x
      end if
      avg_r8_f = m
      return
      end
