      subroutine get_F3_x_c_(n_x,grid_x_c_,max_x_c,n_y,grid_y_c_,max_y_c
     $     ,S_x_c_,get_F3_x_c,param_1,param_2,param_3)
      implicit none
      integer n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      real *8 max_x_c,max_y_c,param_1,param_2,param_3
      complex *16 S_x_c_(0:n_x*n_y-1)
      integer ny,nx
      real *8 X_x_c,Y_x_c
      complex *16 F_x_c
      external get_F3_x_c
      do ny=0,n_y-1
         if (ny.lt.n_y/2) then
            Y_x_c = grid_x_c_(ny)
         else
            Y_x_c = grid_x_c_(ny) - max_x_c
         end if
         do nx=0,n_x-1
            if (nx.lt.n_x/2) then
               X_x_c = grid_x_c_(nx)
            else
               X_x_c = grid_x_c_(nx) - max_x_c
            end if
            call get_F3_x_c(max_x_c,X_x_c,Y_x_c,F_x_c,param_1,param_2
     $           ,param_3)
            S_x_c_(nx+ny*n_x) = F_x_c
         enddo
      enddo
      end
