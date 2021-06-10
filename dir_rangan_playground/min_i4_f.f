!> Doxygen comment: ;\n
!> minimum of integer *4 array. ;\n
      integer *4 function min_i4_f(n_x,x_)
      integer *4 n_x,nx
      integer *4 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0
      else
         m=x_(0)
         do nx=0,n_x-1
            if (x_(nx).lt.m) then
               m = x_(nx)
            end if
         enddo
      end if
      min_i4_f = m
      return
      end
