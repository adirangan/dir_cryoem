      subroutine polyval_r8_reverse_0(n_p,p_,n_x,x_,v_)
c$$$      Evaluates a polynomial given by: ;
c$$$      v_(j) = p_(0) + p_(1)*x_(j) + p_(2)*x_(j)**2 + .. + p_(n_p-1)*x_(j)**(n_p-1) ;
      implicit none
      integer n_p,n_x
      real *8 p_(0:n_p-1)
      real *8 x_(0:n_x-1)
      real *8 v_(0:n_x-1)
      integer np,nx,nj
      real *8 p,x
      do nx=0,n_x-1
         x = x_(nx)
         p = 0.0d0
         nj=n_p-1
         do np=0,n_p-1
            p = p_(nj) + p*x
            nj = nj-1
         enddo                  !do np=0,n_p-1
         v_(nx) = p
      enddo                     !do nx=0,n_x-1
      end
