      subroutine innerproduct_c(n_x,grid_x_c_,n_y,grid_y_c_,T_c_,M_c_
     $     ,C_c)
c$$$      assumes that M_c_ is the same size and dimensions as T_c_
c$$$      assumes uniform grid and periodic boundary conditions
      implicit none
      integer n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      complex *16 T_c_(0:n_x*n_y-1),M_c_(0:n_x*n_y-1),C_c,C_T,C_M
      integer nx,ny,nc
      real *8 dx,dy
      C_c = cmplx( 0.0 , 0.0 )
c$$$      C_T = cmplx( 0.0 , 0.0 )
c$$$      C_M = cmplx( 0.0 , 0.0 )
      dx = (grid_x_c_(n_x-1)-grid_x_c_(0))/max(1,n_x-1)
      dy = (grid_y_c_(n_y-1)-grid_y_c_(0))/max(1,n_y-1)
      nc = 0
      do ny=0,n_y-1
         do nx=0,n_x-1
c$$$            nc = nx + ny*n_x
            C_c = C_c + conjg(T_c_(nc))*M_c_(nc)
c$$$            C_T = C_T + conjg(T_c_(nc))*T_c_(nc)
c$$$            C_M = C_M + conjg(M_c_(nc))*M_c_(nc)
            nc = nc+1  
         enddo
      enddo
c$$$      write(6,'(A,2F8.3)') 'C_T: ',C_T*dx*dy
c$$$      write(6,'(A,2F8.3)') 'C_M: ',C_M*dx*dy
c$$$      write(6,'(A,2F8.3)') 'C_c: ',C_c*dx*dy
      C_c = C_c*dx*dy
c$$$      C_c = C_c / (zsqrt(C_T) * zsqrt(C_M))
c$$$      write(6,'(A,2F8.3)') 'C_c: ',C_c
      end      
