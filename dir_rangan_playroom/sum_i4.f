      subroutine sum_i4(n_x,x_,s)
      integer *4 n_x,nx
      integer *4 x_(0:n_x-1),s
      if (n_x.le.0) then
         s=0
      else
         s=0
         do nx=0,n_x-1
            s = s + x_(nx)
         enddo
      end if
      end
