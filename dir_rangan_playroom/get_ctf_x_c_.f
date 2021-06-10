      subroutine get_ctf_x_c_(n_x,grid_x_c_,max_x_c,n_y,grid_y_c_
     $     ,max_y_c,S_x_c_,get_ctf_x_c,param_0,param_1,param_2,param_3
     $     ,param_4)
      implicit none
      integer n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      real *8 max_x_c,max_y_c
      real *8 param_0,param_1,param_2,param_3,param_4
      complex *16 S_x_c_(0:n_x*n_y-1)
      integer ny,nx
      real *8 X_x_c,Y_x_c
      complex *16 C_x_c
      external get_ctf_x_c
      do ny=0,n_y-1
         if (ny.lt.n_y/2) then
            Y_x_c = grid_y_c_(ny)
         else
            Y_x_c = grid_y_c_(ny) - max_y_c
         end if
c$$$         Y_x_c = grid_y_c_(ny)
         do nx=0,n_x-1
            if (nx.lt.n_x/2) then
               X_x_c = grid_x_c_(nx)
            else
               X_x_c = grid_x_c_(nx) - max_x_c
            end if
c$$$            X_x_c = grid_x_c_(nx)
            call get_ctf_x_c(max_x_c,X_x_c,Y_x_c,C_x_c,param_0,param_1
     $           ,param_2,param_3,param_4)
            S_x_c_(nx+ny*n_x) = C_x_c
         enddo
      enddo
      end
