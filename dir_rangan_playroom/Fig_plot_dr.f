c     Compile with: gfortran -O2 -o Fig_plot_dr.out Fig_plot_dr.f ; 
c     Run with: ./Fig_plot_dr.out ;
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program Fig_plot_dr
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This program is a driver which tests Fig_plot.
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      integer verbose
      data verbose / 2 /
c$$$      grids
      real *8 pi
      integer n_x,nx,n_w,nw
      real *8 min_x_c,max_x_c
      real *8, allocatable :: grid_x_c_(:)
c$$$      data
      complex *16, allocatable :: F_x_c_(:)
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber
      integer *4 line_style,line_thickness,axis_type
      real *8 cval,dash_length
      integer *4 text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering Fig_plot_dr]: '
       end if

      pi = 4*atan(1.0)

      n_x = 128
      min_x_c = -pi
      max_x_c = +pi

      allocate(grid_x_c_(0:n_x-1))
      call linspace(min_x_c,max_x_c,n_x,grid_x_c_)
      allocate(F_x_c_(0:n_x-1))

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We initialize the figure file.' 
      end if

      unitnumber = 7
      write(fname_pre,'(A)') './dir_jpg/Fig_plot_dr'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)')
     $        'We plot the data to the figure file.'
      end if
      
      n_w = 6
      do nw=0,n_w-1
         do nx=0,n_x-1
            F_x_c_(nx) = cmplx(cos((1+nw)*grid_x_c_(nx)) , sin((1+nw)
     $           *grid_x_c_(nx)))
         enddo !nx
c$$$  call stdlim_c16(n_x,.true.,F_x_c_,F_avg,F_std)
c$$$  F_max = F_avg + 2.5*F_std
c$$$  F_min = F_avg - 2.5*F_std
         F_max = +1.0d0
         F_min = -1.0d0
         x_offset = 0.0d0 + 11.0d0*(nw/2)
         y_offset = 0.0d0 + 11.0d0*(mod(nw,2))
         x_side = 10.0d0
         y_side = 10.0d0
         line_style = 0
         line_thickness = 5
         cval = 0.75d0
         dash_length = 5.0d0
         axis_type = 1
         call Fig_c16_plot(unitnumber,line_style,line_thickness,cval
     $        ,dash_length,axis_type,n_x,grid_x_c_,min_x_c,max_x_c
     $        ,.true.,F_x_c_ ,F_min ,F_max,x_offset ,y_offset,x_side
     $        ,y_side)
         line_style = 1
         line_thickness = 4
         cval = 0.25d0
         dash_length = 10.0d0
         call Fig_c16_plot(unitnumber,line_style,line_thickness,cval
     $        ,dash_length,axis_type,n_x,grid_x_c_,min_x_c,max_x_c
     $        ,.false.,F_x_c_ ,F_min ,F_max,x_offset ,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0)') 'F_',(1+nw)
         text_length = 3
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = 0.0d0
         y_loc = +0.90d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
      enddo !nw
      
      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: "gthumb '
     $        ,'Fig_plot_dr.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished Fig_plot_dr]'
      end if

      stop
      end
      
      include 'stdlim_c16.f'
      include 'cp1_c16.f'
      include 'cps_c16.f'
      include 'af1_c16.f'
      include 'afs_c16.f'
      include 'linspace.f'
      include 'interp1_c16.f'
      include 'periodize_r8.f'
      include 'periodize_i.f'
      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Fig_header.f'
      include 'Fig_text.f'
      include 'Fig_c16_plot.f'

