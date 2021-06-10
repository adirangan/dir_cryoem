      subroutine Fig_gen_ver1c(ncur,nlats_,grid_k_p_,n_S,ld_S,S_k_p_,n_M
     $     ,ld_M,M_k_p_,N_k_p_,ld_CTF,CTF_k_p_,n_alpha,alpha_est_
     $     ,delta_x_sort_ ,delta_y_sort_ ,gamma_z_sort_,I_permute_
     $     ,n_M_sample,prefix)
c$$$      Reads in a few templates and images and plots a comparison
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 0 /
      logical plot_flag
      integer *4 ncur,nlats_(0:ncur-1)
      integer *4 n_S,ld_S,n_M,ld_M,ld_CTF,n_alpha,n_M_sample
      real *8 grid_k_p_(0:0)
      real *8 alpha_est_(0:0)
      include 'nalpha_define.f'
      real *8 delta_x_sort_(0:n_S*n_M-1)
      real *8 delta_y_sort_(0:n_S*n_M-1)
      real *8 gamma_z_sort_(0:n_S*n_M-1)
      integer *4 I_permute_(0:n_S*n_M-1)
      complex *16 S_k_p_(0:0),M_k_p_(0:0),N_k_p_(0:0),CTF_k_p_(0:0)
      character(len=1024) prefix
      integer *4 I_M_sample_(0:n_M_sample)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,ns,nm,nctf,nM_sample
      integer *4, allocatable :: n_w_(:)
c$$$      temporary storage for templates and images
      complex *16, allocatable :: T_k_p_(:)
      real *8 pi
      real *8 max_r8_f
      real *8 max_k_c
      real *8 delta_x_tmp,delta_y_tmp
      real *8 delta_x,delta_y,gamma_z
      character(len=64) format_string
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering Fig_gen_ver1c]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'nlats_: ',(nlats_(nr),nr=0,ncur-1)
         write(format_string,'(A,I0,A)') '(A,',ncur,'(F8.3,1X))'
         write(6,format_string) 'grid_k_p_: ',(grid_k_p_(nr),nr=0,ncur
     $        -1)
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         if (verbose.gt.2) then
            write(6,'(A)') 'delta_x_sort^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_x_sort_(na),na=0,n_S*n_M-1)
            write(6,'(A)') 'delta_y_sort^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_y_sort_(na),na=0,n_S*n_M-1)
            write(6,'(A)') 'gamma_z_sort^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (gamma_z_sort_(na),na=0,n_S*n_M-1)
            write(6,'(A)') 'I_permute_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(I0,1X))'
            write(6,format_string) (I_permute_(nr),nr=0,n_S*n_M-1)
         end if
         write(6,'(A,I0)') 'n_M_sample: ',n_M_sample
         write(6,'(A)') trim(prefix)
      end if

      do nM_sample=0,n_M_sample-1
         if (n_M_sample.lt.n_M) then
            I_M_sample_(nM_sample) = min(n_M-1,floor(n_M
     $           *1.0d0*nM_sample/(n_M_sample-1)))
         else
            I_M_sample_(nM_sample) = nM_sample
         end if
      enddo
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(',n_M_sample,'(I0,1X))'
         write(6,format_string) (I_M_sample_(nM_sample),nM_sample=0
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
     $        ,' along with their best matching (transformed) templates'
     $        ,' also on a k-space polar grid.'
      end if
c$$$      printing templates and images
      do nM_sample=0,n_M_sample-1
         nm = I_M_sample_(nM_sample)
         nctf = nint(alpha_est_(nalpha_ctf_ind + nm*n_alpha))
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') 'sample ',nM_sample,' is image ',nm
            write(6,'(A,I0,A,I0)') 'sample ',nM_sample,' has ctf ',nctf
         end if
c$$$         Display image without ctf
         call stdlim_c16(n_A,.true.,M_k_p_(nm*ld_M),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0 - 22.5d0
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
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 10.0d0 - 22.5d0
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
c$$$         Display image with ctf
         call stdlim_c16(n_A,.true.,N_k_p_(nm*ld_M),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,N_k_p_(nm*ld_M),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'R(N_k_p(',nm,')) '
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
         call stdlim_c16(n_A,.false.,N_k_p_(nm*ld_M),F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 10.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,N_k_p_(nm*ld_M),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'I(N_k_p(',nm,')) '
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
c$$$         Pick best matching template
         ns = I_permute_(n_S-1 + nm*n_S)
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0)') 'Best match to image ',nm,' is ',ns
         end if
c$$$         Display template without perturbation
         ns = I_permute_(n_S-1 + nm*n_S)
         call cp1_c16(n_A,S_k_p_(ns*ld_S),T_k_p_)
         call stdlim_c16(n_A,.true.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0 + 67.5d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
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
         call stdlim_c16(n_A,.false.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 10.0d0 + 67.5d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
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
c$$$         Display template with perturbation but no ctf
         ns = I_permute_(n_S-1 + nm*n_S)
         call get_interchange_delta(delta_x_sort_(n_S-1 + nm*n_S)
     $        ,delta_y_sort_(n_S-1 + nm*n_S),alpha_est_(nalpha_gamma_z +
     $        nm*n_alpha),delta_x_tmp,delta_y_tmp)
         delta_x = alpha_est_(nalpha_delta_x + nm*n_alpha) +
     $        delta_x_tmp
         delta_y = alpha_est_(nalpha_delta_y + nm*n_alpha) +
     $        delta_y_tmp
         gamma_z = alpha_est_(nalpha_gamma_z + nm*n_alpha) +
     $        gamma_z_sort_(n_S-1 + nm*n_S)
         call cp1_c16(n_A,S_k_p_(ns*ld_S),T_k_p_)
         call transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,+delta_x,
     $        +delta_y,T_k_p_)
         call rotate_p_to_p(n_r,n_w_,n_A,T_k_p_,+gamma_z,T_k_p_)
         call stdlim_c16(n_A,.true.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nM_sample
         y_offset = +17.0d0 - 0.0d0 + 45.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'R(T_k_p(',ns,')) '
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
         y_offset = +17.0d0 - 10.0d0 + 45.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,T_k_p_,F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A,I4,A)') 'I(T_k_p(',ns,')) '
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
c$$$         Display template with perturbation and ctf
         ns = I_permute_(n_S-1 + nm*n_S)
         call get_interchange_delta(delta_x_sort_(n_S-1 + nm*n_S)
     $        ,delta_y_sort_(n_S-1 + nm*n_S),alpha_est_(nalpha_gamma_z +
     $        nm*n_alpha),delta_x_tmp,delta_y_tmp)
         delta_x = alpha_est_(nalpha_delta_x + nm*n_alpha) +
     $        delta_x_tmp
         delta_y = alpha_est_(nalpha_delta_y + nm*n_alpha) +
     $        delta_y_tmp
         gamma_z = alpha_est_(nalpha_gamma_z + nm*n_alpha) +
     $        gamma_z_sort_(n_S-1 + nm*n_S)
         call cp1_c16(n_A,S_k_p_(ns*ld_S),T_k_p_)
         call transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,+delta_x,
     $        +delta_y,T_k_p_)
         call rotate_p_to_p(n_r,n_w_,n_A,T_k_p_,+gamma_z,T_k_p_)
         call xx1_c16(n_A,T_k_p_,CTF_k_p_(nctf*ld_CTF),T_k_p_)
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
         write(fig_title,'(A,I0,A,I4,A)') 'R(CT_k_p(',ns,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 19
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
         write(fig_title,'(A,I0,A,I4,A)') 'I(CT_k_p(',ns,')) '
     $        ,nint(dlog(F_std)),'    '
         text_length = 19
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
         write(6,'(A)') '[finished Fig_gen_ver1c]'
      end if

      end
