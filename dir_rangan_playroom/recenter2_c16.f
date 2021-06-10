      subroutine recenter2_c16(n_x,n_y,a_,b_)
      integer *4 n_x,n_y
      complex *16 a_(0:n_x*n_y-1),b_(0:n_x*n_y-1)
      integer *4 fx,fy,nxt,nyt
      integer *4, allocatable :: x_(:)
      integer *4, allocatable :: y_(:)
      complex *16, allocatable :: c_(:)
      fx = floor(n_x/2.0d0)
      fy = floor(n_y/2.0d0)
      allocate(x_(0:n_x-1))
      allocate(y_(0:n_y-1))
      nxt = fx
      do nx=0,n_x-1
         x_(nx) = nxt
         nxt = nxt+1
         if (nxt.ge.n_x) nxt=0
      enddo !do nx=0,n_x-1
      nyt = fy
      do ny=0,n_y-1
         y_(ny) = nyt
         nyt = nyt+1
         if (nyt.ge.n_y) nyt=0
      enddo !do ny=0,n_y-1
      allocate(c_(0:n_x*n_y-1))
      do ny=0,n_y-1
         nyt = y_(ny)
         do nx=0,n_x-1
            nxt = x_(nx)
            c_(nx+ny*n_x) = a_(nxt+nyt*n_x)
         enddo
      enddo
      do ny=0,n_y-1
         do nx=0,n_x-1
            b_(nx+ny*n_x) = c_(nx+ny*n_x)
         enddo
      enddo
      deallocate(c_)
      end      
