      subroutine Fig_c16_carte(unitnumber,n_x,grid_x_c_,n_y,grid_y_c_
     $     ,real_or_imag,F_,Fmin,Fmax,x_offset,y_offset,x_side,y_side)
      implicit none
      integer *4 unitnumber
      integer *4 n_x,n_y
      real *8 grid_x_c_(0:n_x-1),grid_y_c_(0:n_y-1)
      logical real_or_imag
      complex *16 F_(0:n_x*n_y-1)
      real *8 Fmin,Fmax
      real *8 x_offset,y_offset,x_side,y_side
      real *8 F_avg,F_std,Fmin_use,Fmax_use
      complex *16, allocatable :: T_(:)
      real *8 x_pos,y_pos,x_pre,y_pre
      integer *4 diameter,nx,ny,nd
      integer *4 xyval_(0:9)
      integer *4 color_index
      character(len=20) format_string
      allocate(T_(0:n_x*n_y-1))
      call recenter_c16(n_x,n_y,F_,T_)
      write(format_string,'(A,I0,A)') '(A,1X,',10,'(I0,1X))'
      Fmin_use = Fmin
      Fmax_use = Fmax
      if (Fmin_use.ge.Fmax_use) then
         call stdlim_c16(n_x*n_y,real_or_imag,F_,F_avg,F_std)
         Fmin_use = F_avg - 2.0*F_std
         Fmax_use = F_avg + 2.0*F_std
      end if
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
            if (real_or_imag.eqv..true.) then
               color_index = 32 + int(dmax1(0.0,dmin1(511.0,511.0
     $              *(real(T_(nx+ny*n_x))-Fmin_use)/(Fmax_use
     $              -Fmin_use))))
            else
               color_index = 32 + int(dmax1(0.0,dmin1(511.0,511.0
     $              *(aimag(T_(nx+ny*n_x))-Fmin_use)/(Fmax_use
     $              -Fmin_use))))
            end if
            write(unitnumber,'(A,I0,A,I0)') '2 3 0 0 0 ',color_index
     $           ,' 50 -1 20 0.000 0 0 -1 0 0 ',5
            write(unitnumber,format_string) char(9),(xyval_(nd),nd=0,9)
         enddo
      enddo
      deallocate(T_)
      end
