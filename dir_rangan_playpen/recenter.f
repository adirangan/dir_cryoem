      subroutine recenter(n_x,n_y,a_,b_)
      integer n_x,n_y,nx,ny,nxt,nyt
      real *8 a_(0:n_x*n_y-1),b_(0:n_x*n_y-1)
      do ny=0,n_y-1
         if (ny.lt.n_y/2) then
            nyt = ny + n_y/2
         else
            nyt = ny - n_y/2
         end if
         do nx=0,n_x-1
            if (nx.lt.n_x/2) then
               nxt = nx + n_x/2
            else
               nxt = nx - n_x/2
            end if
            b_(nxt+nyt*n_x) = a_(nx+ny*n_x)
         enddo
      enddo
      end      
