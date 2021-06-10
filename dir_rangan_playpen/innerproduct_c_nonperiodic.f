      subroutine innerproduct_c_nonperiodic(n_x,grid_x_c_,n_y,grid_y_c_
     $     ,T_c_,M_c_,C_c)
c$$$      assumes that M_c_ is the same size and dimensions as T_c_
c$$$      assumes non-periodic boundary conditions
      implicit none
      integer n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      complex *16 T_c_(0:n_x*n_y-1),M_c_(0:n_x*n_y-1),C_c
      integer nx,ny,nx_pre,nx_pos,ny_pre,ny_pos,nc
      real *8 dx,dy
      C_c = cmplx( 0.0 , 0.0 )
      nc = 0
      do ny=0,n_y-1
         if (ny.gt.0) then
            ny_pre = ny-1
         else
            ny_pre = ny
         end if
         if (ny.lt.n_y-1) then
            ny_pos = ny+1
         else
            ny_pos = ny
         end if
         dy = 0.5*(grid_y_c_(ny_pos) + grid_y_c(ny)) - 0.5
     $        *(grid_y_c_(ny) + grid_y_c(ny_pre))
         do nx=0,n_x-1
            if (nx.gt.0) then
               nx_pre = nx-1
            else
               nx_pre = nx
            end if
            if (nx.lt.n_x-1) then
               nx_pos = nx+1
            else
               nx_pos = nx
            end if
            dx = 0.5*(grid_x_c_(nx_pos) + grid_x_c(nx)) - 0.5
     $           *(grid_x_c_(nx) + grid_x_c(nx_pre))
c$$$            nc = nx + ny*n_x
            C_c = C_c + conjg(T_c_(nc))*M_c_(nc)*dx*dy
            nc = nc+1  
         enddo
      enddo
      end
