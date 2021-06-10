      subroutine Fig_carte(unitnumber,n_x,grid_x_c_,n_y,grid_y_c_,F_
     $     ,Fmin,Fmax,x_offset,y_offset,x_side,y_side)
      integer *4 unitnumber
      integer *4 n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      real *8 F_(0:n_x*n_y-1),Fmin,Fmax
      real *8 x_offset,y_offset,x_side,y_side
      real *8 x_pos,y_pos,x_pre,y_pre
      integer *4 diameter,nx,ny
      integer *4 xyval_(0:9)
      integer *4 color_index
      character(len=20) format_string
      write(format_string,'(A,I0,A)') '(A,1X,',10,'(I0,1X))'
      diameter = 1000
      do ny=0,n_y-1
         if (ny.gt.0) then
            y_pre = 0.5*(grid_y_c_(ny-1) + grid_y_c_(ny))
         else
            y_pre = grid_y_c_(0)
         end if
         if (ny.lt.n_y-1) then
            y_pos = 0.5*(grid_y_c_(ny) + grid_y_c_(ny+1))
         else
            y_pos = grid_y_c_(n_y-1)
         end if
         do nx=0,n_x-1
            if (nx.gt.0) then
               x_pre = 0.5*(grid_x_c_(nx-1) + grid_x_c_(nx))
            else
               x_pre = grid_x_c_(0)
            end if
            if (nx.lt.n_x-1) then
               x_pos = 0.5*(grid_x_c_(nx) + grid_x_c_(nx+1))
            else
               x_pos = grid_x_c_(n_x-1)
            end if
            xyval_(0 + 2*0) = +int(diameter*(x_pre*x_side + x_offset))
            xyval_(1 + 2*0) = -int(diameter*(y_pre*y_side + y_offset))
            xyval_(0 + 2*1) = +int(diameter*(x_pos*x_side + x_offset))
            xyval_(1 + 2*1) = -int(diameter*(y_pre*y_side + y_offset))
            xyval_(0 + 2*2) = +int(diameter*(x_pos*x_side + x_offset))
            xyval_(1 + 2*2) = -int(diameter*(y_pos*y_side + y_offset))
            xyval_(0 + 2*3) = +int(diameter*(x_pre*x_side + x_offset))
            xyval_(1 + 2*3) = -int(diameter*(y_pos*y_side + y_offset))
            xyval_(0 + 2*4) = +int(diameter*(x_pre*x_side + x_offset))
            xyval_(1 + 2*4) = -int(diameter*(y_pre*y_side + y_offset))
            color_index = 32 + int(dmax1(0.0,dmin1(511.0,511.0*(F_(nx+ny
     $           *n_x)-Fmin)/(Fmax-Fmin))))
            write(unitnumber,'(A,I0,A,I0)') '2 3 0 0 0 ',color_index
     $           ,' 50 -1 20 0.000 0 0 -1 0 0 ',5
            write(unitnumber,format_string) char(9),(xyval_(j),j=0,9)
         enddo
      enddo
      end
