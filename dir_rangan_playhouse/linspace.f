      subroutine linspace(d_start,d_endat,n_steps,grid_)
      implicit none
      real *8 d_start,d_endat
      integer n_steps
      real *8 grid_(0:n_steps-1)
      integer ns
      do ns=0,n_steps-1
         grid_(ns) = d_start + ns*(d_endat-d_start)/(n_steps)
      enddo
      end
