      subroutine Fig_gen_ver0(ncur,nlats_,grid_k_p_,n_S,ld_S,S_k_p_,n_M
     $     ,ld_M,M_k_p_)
      implicit none
      integer verbose
      data verbose / 2 /
      integer *4 ncur,nlats_(0:ncur-1),n_S,ld_S,n_M,ld_M
      real *8 grid_k_p_(0:0)
      complex *16 S_k_p_(0:0),M_k_p_(0:0)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,ns,nm
      integer *4, allocatable :: n_w_(:)
c$$$      grids
      real *8 pi
      real *8 max_k_c
      real *8 Flim_x_c,Flim_x_p,Flim_k_c,Flim_k_p,Flim_k_q
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering Fig_gen_ver0]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
      end if

      pi = 4*atan(1.0)

c$$$      Calculating template size using 'get_template_size'
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'ncur = ',ncur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
      end if
      
c$$$      indices
      n_r = ncur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r',n_r,'<2'
      end if
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_w_max ',n_w_max,'; n_A ',n_A
      end if

c$$$      generating grids for templates and images 
      max_k_c = 2*grid_k_p_(n_r-1)

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We initialize the first figure file.' 
      end if

      unitnumber = 7
      write(fname_pre,'(A)') 'Fig_gen_ver0'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)')
     $        'We display the templates and images to the figure file.'
      end if
c$$$      printing templates and images
      do ns=0,n_S-1
         call stdlim_c16(n_A,.true.,S_k_p_(ns*ld_S),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*ns
         y_offset = -5.0d0 - 0.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,S_k_p_(ns*ld_S),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'R(S_k_p(',ns,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 18
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = -0.5d0*(max_k_c)
         y_loc = -0.5d0*(max_k_c)
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         call stdlim_c16(n_A,.false.,S_k_p_(ns*ld_S),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*ns
         y_offset = -5.0d0 - 10.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,S_k_p_(ns*ld_S),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'I(S_k_p(',ns,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 18
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = -0.5d0*(max_k_c)
         y_loc = -0.5d0*(max_k_c)
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
      enddo
      do nm=0,n_M-1
         call stdlim_c16(n_A,.true.,M_k_p_(nm*ld_M),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nm
         y_offset = +17.0d0 - 0.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,M_k_p_(nm*ld_M),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'R(M_k_p(',nm,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 18
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = -0.5d0*(max_k_c)
         y_loc = -0.5d0*(max_k_c)
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         call stdlim_c16(n_A,.false.,M_k_p_(nm*ld_M),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nm
         y_offset = +17.0d0 - 10.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,M_k_p_(nm*ld_M),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'I(M_k_p(',nm,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 18
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = -0.5d0*(max_k_c)
         y_loc = -0.5d0*(max_k_c)
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: "gthumb '
     $        ,'Fig_gen_ver0_pre.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

 10   continue ! if plot_flag.eqv..false.
      if (verbose.gt.0) then
         write(6,'(A)') '[finished Fig_gen_ver0]'
      end if

      end
