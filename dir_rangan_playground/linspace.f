!> Doxygen comment: ;\n
!> generate linearly spaced array. ;\n
!> note that this does not include the upper endpoint! ;\n
!> this is specifically for use with periodic array. ;\n
!> for example: linspace(0,2.0d0*pi,16,grid_) generates ;\n
!> 16 equispaced points from 0 to (15/16)*2.0d0*pi. ;\n
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
