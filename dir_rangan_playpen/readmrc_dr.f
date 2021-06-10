c$$$      gfortran -o readmrc_dr.out readmrc_dr.f ; ./readmrc_dr ;
      program readmrc_dr
c$$$      This function tests readmrc.f using the file "normst_66033.mrc"
c$$$      found at : 
c$$$      https://www.ebi.ac.uk/pdbe/emdb/empiar/entry/10065/
c$$$      with a (local) name of :
c$$$      '/data/rangan/dir_cryoem/archive/10065/data/normst_66033.mrc'
c$$$      You should change this local filename (below).
c$$$      This function generates an image named:
c$$$      readmrc_dr.jpg
c$$$      which illustrates the first 9 images in the mrc file.
      implicit none
      integer verbose
      data verbose / 2 /
      character(len=1024) filename
      integer n_x,n_y,n_M_all,n_M,nm,n_r
      complex *16, allocatable :: M_(:)
      complex *16, allocatable :: T_(:)
      real *8 max_x_c,F_avg,F_std
      real *8, allocatable :: grid_x_c_(:)
      character(len=1024) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      write(filename,'(A,A)') '/data/rangan/dir_cryoem/'
     $     ,'archive/10065/data/normst_66033.mrc'
      if (verbose.gt.0) then
         write(6,'(A,A)') 'opening: ',trim(filename)
      end if
      call readmrc_size(filename,n_x,n_y,n_M_all)
      n_r = n_x
      max_x_c = 1.0d0
      allocate(grid_x_c_(0:n_r-1))
      call linspace(0.0d0,max_x_c,n_r,grid_x_c_)
      n_M = min(9,n_M_all)
      allocate(M_(n_x*n_y*n_M))
      allocate(T_(n_x*n_y))
      call readmrc_data(filename,n_M,M_)
      do nm=0,n_M-1
         call recenter_c16(n_r,n_r,M_(nm*n_r*n_r),T_)
         call cp1_c16(n_r*n_r,T_,M_(nm*n_r*n_r))
      enddo
      unitnumber = 7
      write(fname_pre,'(A)') 'readmrc_dr'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);
      do nm=0,n_M-1
         call cp1_c16(n_r*n_r,M_(nm*n_r*n_r),T_)
         call stdlim_c16(n_r*n_r,.true.,T_,F_avg,F_std)
         F_max = F_avg + 3.5*F_std
         F_min = F_avg - 3.5*F_std
         x_offset = 0.0d0 + 10.0d0*mod(nm,3)
         y_offset = 0.0d0 + 10.0d0*(nm/3)
         x_side = 10.0d0
         y_side = 10.0d0
         call Fig_c16_carte(unitnumber,n_r,grid_x_c_,n_r,grid_x_c_
     $        ,.true.,T_,F_min,F_max,x_offset,y_offset
     $        ,x_side,y_side)
         write(fig_title,'(A,I0,A)') 'M_(',nm,')'
         text_length = 5
         color_index = 18
         font_type = 0
         font_size = 72
         x_loc = +0.0d0
         y_loc = +0.0d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
      enddo
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)
      
      stop
      end

      include "readmrc.f"
      include 'stdlim_c16.f'
      include 'pearson_c16.f'
      include 'cp1_c16.f'
      include 'cps_c16.f'
      include 'af1_c16.f'
      include 'afs_c16.f'
      include 'linspace.f'
      include 'periodize_r8.f'
      include 'periodize_i.f'
      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Fig_header.f'
      include 'Fig_text.f'
      include 'recenter_c16.f'
      include 'Fig_c16_carte.f'
