!> Doxygen comment: ;\n
!> sum of complex *16 array. ;\n
      complex *16 function sum_c16_f(n_x,x_)
      integer *4 n_x,nx
      complex *16 x_(0:n_x-1),s
      if (n_x.le.0) then
         s=0
      else
         s=0
         do nx=0,n_x-1
            s = s + x_(nx)
         enddo
      end if
      sum_c16_f = s
      return
      end
