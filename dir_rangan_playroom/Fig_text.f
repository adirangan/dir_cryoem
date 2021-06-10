      subroutine Fig_text(unitnumber,n_t,text,color_index,font_type
     $     ,font_size,x_loc,y_loc,x_offset,y_offset,x_side,y_side)
      integer *4 unitnumber,n_t
      real *8 x_loc,y_loc
      character text(0:n_t-1)
      integer *4 color_index,font_type,font_size
      real *8 x_offset,y_offset,x_side,y_side
      integer *4 diameter
      integer *4 xyval_(0:3)
      character(len=48) format_string
      write(format_string,'(A,I0,A)') '(A,I0,A,I0,1X,I0,A,(4(I0,1X)),('
     $     ,n_t,'(A)),A)'
      diameter = 1000
      xyval_(0) = 0
      xyval_(1) = 0
      xyval_(2) = +int(diameter*(x_loc*x_side + x_offset))
      xyval_(3) = -int(diameter*(y_loc*y_side + y_offset))
      write(unitnumber,format_string) '4 0 ' ,color_index,' 50 -1 '
     $     ,font_type,font_size,' 0.0000 4 ' ,(xyval_(j),j=0,3)
     $     ,text(0:n_t-1),'\001'
      end
