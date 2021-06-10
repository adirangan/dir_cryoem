      real *8 function max_r8_f(n_x,x_)
      integer n_x,nx
      real *8 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0.0d0
      else
         m=x_(0)
         do nx=0,n_x-1
            if (x_(nx).gt.m) then
               m = dabs(x_(nx))
            end if
         enddo
      end if
      max_r8_f = m
      return
      end
