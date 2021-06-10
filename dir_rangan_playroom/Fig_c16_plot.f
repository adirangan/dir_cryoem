      subroutine Fig_c16_plot(unitnumber,line_style,line_thickness,cval
     $     ,dash_length,axis_type,n_x,grid_x_c_,xmin,xmax,real_or_imag
     $     ,F_,Fmin ,Fmax,x_offset ,y_offset,x_side,y_side)
      implicit none
      integer *4 unitnumber
      integer *4 line_style,line_thickness,axis_type
      real *8 cval,dash_length
      integer *4 n_x
      real *8 grid_x_c_(0:n_x-1)
      real *8 xmin,xmax
      logical real_or_imag
      complex *16 F_(0:n_x-1)
      real *8 Fmin,Fmax
      real *8 x_offset,y_offset,x_side,y_side
      real *8 xmin_use,xmax_use
      real *8 F_avg,F_std,Fmin_use,Fmax_use
      real *8 x_pos,y_pos
      integer *4 diameter,nx,ny,nd
      integer *4, allocatable :: xyval_(:)
      integer *4, allocatable :: axval_(:)
      integer *4 color_index
      character(len=20) format_string
      xmin_use = xmin
      xmax_use = xmax
      if (xmin_use.ge.xmax_use) then
         xmin_use = grid_x_c_(0)
         xmax_use = grid_x_c_(n_x-1)
      end if
      Fmin_use = Fmin
      Fmax_use = Fmax
      if (Fmin_use.ge.Fmax_use) then
         call stdlim_c16(n_x,real_or_imag,F_,F_avg,F_std)
         Fmin_use = F_avg - 2.0*F_std
         Fmax_use = F_avg + 2.0*F_std
      end if
      allocate(xyval_(0:2*n_x-1))
      diameter = 1000
      x_pos = 0.0d0
      y_pos = 0.0d0
      do nx=0,n_x-1
         x_pos = (grid_x_c_(nx)-xmin_use)/(xmax_use-xmin_use)
         if (real_or_imag.eqv..true.) then
            y_pos = (real(F_(nx))-Fmin_use)/(Fmax_use -Fmin_use)
         else
            y_pos = (aimag(F_(nx))-Fmin_use)/(Fmax_use -Fmin_use)
         end if
         xyval_(0 + 2*nx) = +int(diameter*(x_pos*x_side + x_offset))
         xyval_(1 + 2*nx) = -int(diameter*(y_pos*y_side + y_offset))
      enddo
      color_index = 32 + int(dmax1(0.0,dmin1(511.0,511.0*cval)))
      write(format_string,'(A,I0,A)') '(A,1X,',2*n_x,'(I0,1X))'
      write(unitnumber,'(A,I0,A,I0,A,I0,A,F8.3,A,I0)') '2 1 ',line_style
     $     ,' ',line_thickness,' ',color_index,' 7 50 -1 -1 '
     $     ,dash_length,' 0 0 -1 0 0 ',n_x
      write(unitnumber,format_string) char(9),(xyval_(nd),nd=0,2*n_x-1)
      allocate(axval_(0:2*5-1))
      if (axis_type.eq.1) then
         do nx=0,5-1
            x_pos = mod((0+nx)/2,2)
            y_pos = mod((1+nx)/2,2)
            axval_(0 + 2*nx) = +int(diameter*(x_pos*x_side + x_offset))
            axval_(1 + 2*nx) = -int(diameter*(y_pos*y_side + y_offset))
         enddo
         write(format_string,'(A,I0,A)') '(A,1X,',10,'(I0,1X))'
         write(unitnumber,'(A)') '2 1 0 2 0 7 51 -1 -1 0.0 0 0 -1 0 0 5'
         write(unitnumber,format_string) char(9),(axval_(nd),nd=0,9)
      end if                    !axis_type.eq.1

      deallocate(axval_)
      deallocate(xyval_)
      end
