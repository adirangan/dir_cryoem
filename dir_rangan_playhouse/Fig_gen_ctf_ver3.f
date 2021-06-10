      subroutine Fig_gen_ctf_ver3(ncur,nlats_,grid_k_p_,n_M,I_M_sample_
     $     ,ld_M,M_k_p_,ld_CTF,CTF_k_p_,alpha_est_
     $     ,n_M_sample,prefix)
c$$$      Reads in a few images and plots them after
c$$$      convolution with their ctf-function
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 2 /
      logical plot_flag
      integer *4 ncur,nlats_(0:ncur-1),n_M,I_M_sample_(0:n_M-1),ld_M
     $     ,ld_CTF,n_M_sample
      real *8 grid_k_p_(0:0)
      real *8 alpha_est_(0:0)
      include 'nalpha_define.f'
      complex *16 M_k_p_(0:0),CTF_k_p_(0:0)
      character(len=1024) prefix
      integer *4 I_M_sample_tmp_(0:n_M_sample)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,nctf,nm,nM_sample
      integer *4, allocatable :: n_w_(:)
c$$$      temporary storage for ctf
      complex *16, allocatable :: T_k_p_(:)
      real *8 pi
      real *8 max_r8_f
      real *8 max_k_c
      character(len=64) format_string
c$$$      parameters for figure generation
      character(len=1024) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering Fig_gen_ctf_ver3]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'nlats_: ',(nlats_(nr),nr=0,ncur-1)
         write(format_string,'(A,I0,A)') '(A,',ncur,'(F8.3,1X))'
         write(6,format_string) 'grid_k_p_: ',(grid_k_p_(nr),nr=0,ncur
     $        -1)
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         write(6,'(A,I0)') 'n_alpha: ',n_alpha
         write(6,'(A,I0)') 'n_M_sample: ',n_M_sample
         write(6,'(A)') trim(prefix)
      end if

      do nM_sample=0,n_M_sample-1
         if (n_M_sample.lt.n_M) then
            I_M_sample_tmp_(nM_sample) = min(n_M-1,floor(n_M
     $           *1.0d0*nM_sample/(n_M_sample-1)))
         else
            I_M_sample_tmp_(nM_sample) = nM_sample
         end if
      enddo
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(',n_M_sample,'(I0,1X))'
         write(6,format_string) (I_M_sample_tmp_(nM_sample),nM_sample=0
     $        ,n_M_sample-1)
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
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'ngridc_: ',(ngridc_(nr),nr=1,ncur)
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
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,ncur-1)
      end if
      
c$$$      generating grids for templates and images 
      max_k_c = 2*grid_k_p_(n_r-1)

c$$$      allocating temporary storage for images and templates 
      allocate(T_k_p_(0:n_A-1))

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We initialize the first figure file.' 
      end if

      unitnumber = 7
      write(fname_pre,'(A)') trim(prefix)
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)')
     $        'We display the images on a k-space polar grid'
     $        ,' along with their values after ctf-convolution.'
      end if
c$$$      printing images
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_tmp_(nM_sample)
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') 'sample ',nM_sample,' is image ',nm
         end if
         call stdlim_c16(n_A,.true.,M_k_p_(I_M_sample_(nm)*ld_M),F_avg
     $        ,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,M_k_p_(I_M_sample_(nm)*ld_M),F_min,F_max,x_offset
     $        ,y_offset,x_side ,y_side)
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
         call stdlim_c16(n_A,.false.,M_k_p_(I_M_sample_(nm)*ld_M),F_avg
     $        ,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 10.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,M_k_p_(I_M_sample_(nm)*ld_M),F_min,F_max,x_offset
     $        ,y_offset,x_side ,y_side)
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
c$$$         Print ctf-function associated with image
         nctf = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') 'ctf for image ',nm,' is ',nctf
         end if
c$$$         Display ctf-function
         call cp1_c16(n_A,CTF_k_p_(nctf*ld_CTF),T_k_p_)
         call stdlim_c16(n_A,.true.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0 + 22.5d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'R(C_k_p(',nctf,')) '
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
         call stdlim_c16(n_A,.false.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 10.0d0 + 22.5d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'I(C_k_p(',nctf,')) '
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
      enddo !do nM_sample=0,n_M_sample-1

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: '
     $        ,'"gthumb ?.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

      deallocate(T_k_p_)
      deallocate(n_w_)
      deallocate(ngridc_)
      deallocate(icstart_)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished Fig_gen_ctf_ver3]'
      end if

      end
