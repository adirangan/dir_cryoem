c     gfortran -o adi_fft.out adi_fft.f ; ./adi_fft.out ; 
      program adi_fft
      integer *4 n,j,k,nx,ny
      real *8 pi
      complex *16 ci
      complex *16, allocatable ::  a_(:)
      complex *16, allocatable ::  b_(:)
      real *8, allocatable :: grid_kx_(:)
      real *8, allocatable :: grid_ky_(:)
      character(len=64) fig_title
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      pi = 4.0*atan(1.0)
      ci = ( 0.0 , 1.0 )
      n = 8
      write(6,'(A)') ' %% testing adi_fft1 forwards '
      allocate(a_(0:n-1))
      allocate(b_(0:n-1))
      do j=0,n-1
c$$$         a_(j) = exp(-2*pi*ci*1.0*j/n)
         a_(j) = sin(-2*pi*2.0*j/n)
      enddo
      call adi_fft1(-1,n,a_,0,1,b_,0,1)
      do j=0,n-1
         write(6,'(I0,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3)') j,real(a_(j))
     $        ,aimag(a_(j)),real(b_(j)),aimag(b_(j))
      enddo
      deallocate(a_)
      deallocate(b_)
      write(6,'(A)') ' %% testing adi_fft1 backwards '
      allocate(a_(0:n-1))
      allocate(b_(0:n-1))
      do j=0,n-1
c$$$         a_(j) = exp(-2*pi*ci*1.0*j/n)
         a_(j) = sin(-2*pi*2.0*j/n)
      enddo
      call adi_fft1(+1,n,a_,0,1,b_,0,1)
      do j=0,n-1
         write(6,'(I0,1X,F6.3,1X,F6.3,1X,F6.3,1X,F6.3)') j,real(a_(j))
     $        ,aimag(a_(j)),real(b_(j)),aimag(b_(j))
      enddo
      deallocate(a_)
      deallocate(b_)

      write(6,'(A)') ' %% testing adi_fft2 forwards '
      nx = 2*n
      ny = 1*n
      allocate(grid_ky_(0:ny-1))
      allocate(grid_kx_(0:nx-1))
      allocate(a_(0:nx*ny-1))
      allocate(b_(0:nx*ny-1))
      do j=0,nx-1
         grid_kx_(j) = (1.0d0*j)/(nx)
      enddo
      do k=0,ny-1
         grid_ky_(k) = (1.0d0*k)/(ny)
      enddo
      do k=0,ny-1
         do j=0,nx-1
            a_(j+k*(nx)) = +0.5*exp(+2*pi*ci*5.0*j/(nx))*exp(+2*pi*ci
     $           *3.0 *k/(ny)) - 0.5*exp(+2*pi*ci*7.0*j/(nx))*exp(+2*pi
     $           *ci*2.0*k/(ny))
         enddo
      enddo
      call adi_fft2(-1,nx,ny,a_,b_)
      unitnumber = 7
      open(unitnumber,file='adi_fft2_frwd.fig',status='replace',form
     $     ='formatted')
      call Figheader(unitnumber);
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 0.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.true.,a_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'real(a_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 10.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.false.,a_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'imag(a_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0*nx*ny
      F_min = -1.0d0*nx*ny
      x_offset = 0.0d0
      y_offset = 10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.true.,b_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'real(b_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0*nx*ny
      F_min = -1.0d0*nx*ny
      x_offset = 10.0d0
      y_offset = 10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.false.,b_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'imag(b_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A)')
     $     'fig2dev -Ljpeg -q 10 adi_fft2_frwd.fig adi_fft2_frwd.jpg;'
      system_error = system(fig2dev_system_call)
      deallocate(grid_kx_)
      deallocate(grid_ky_)
      deallocate(a_)
      deallocate(b_)

      write(6,'(A)') ' %% testing adi_fft2 backwards '
      nx = 2*n
      ny = 1*n
      allocate(grid_ky_(0:ny-1))
      allocate(grid_kx_(0:nx-1))
      allocate(a_(0:nx*ny-1))
      allocate(b_(0:nx*ny-1))
      do j=0,nx-1
         grid_kx_(j) = (1.0d0*j)/(nx)
      enddo
      do k=0,ny-1
         grid_ky_(k) = (1.0d0*k)/(ny)
      enddo
      do k=0,ny-1
         do j=0,nx-1
            a_(j+k*(nx)) = +0.5*exp(+2*pi*ci*5.0*j/(nx))*exp(+2*pi*ci
     $           *3.0 *k/(ny)) - 0.5*exp(+2*pi*ci*7.0*j/(nx))*exp(+2*pi
     $           *ci*2.0*k/(ny))
         enddo
      enddo
      call adi_fft2(+1,nx,ny,a_,b_)
      unitnumber = 7
      open(unitnumber,file='adi_fft2_back.fig',status='replace',form
     $     ='formatted')
      call Figheader(unitnumber);
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 0.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.true.,a_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'real(a_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 10.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.false.,a_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'imag(a_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 0.0d0
      y_offset = 10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.true.,b_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'real(b_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +1.0d0
      F_min = -1.0d0
      x_offset = 10.0d0
      y_offset = 10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,nx,grid_kx_,ny,grid_ky_,.false.,b_
     $     ,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'imag(b_)'
      text_length = 8
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A)')
     $     'fig2dev -Ljpeg -q 10 adi_fft2_back.fig adi_fft2_back.jpg;'
      system_error = system(fig2dev_system_call)
      deallocate(grid_kx_)
      deallocate(grid_ky_)
      deallocate(a_)
      deallocate(b_)
      stop
      end

      include 'adi_fft1.f'
      include 'adi_fft2.f'
      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Figheader.f'
      include 'Fig_text.f'
      include 'Fig_c16_carte.f'

