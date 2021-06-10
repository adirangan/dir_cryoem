      subroutine max_i4(n_x,x_,m)
      integer *4 n_x,nx
      integer *4 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0
      else
         m=x_(0)
         do nx=0,n_x-1
            if (x_(nx).gt.m) then
               m = x_(nx)
            end if
         enddo
      end if
      end
