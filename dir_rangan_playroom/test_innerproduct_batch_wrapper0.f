      subroutine test_innerproduct_batch_wrapper_0(ncur,nlats_,grid_k_p_
     $     ,max_x_c,n_S,ld_S,S_k_p_,n_M,ld_M,M_k_p_,alpha_tru_
     $     ,alpha_est_,eps_target,N_pixels_in,n_delta_x,n_delta_y
     $     ,n_gamma_z,svd_calculation_type)
c$$$      Reads in a few templates and images and compares them, plotting 
c$$$      output to figure files. The output of this analysis (i.e., which
c$$$      template-image-pairs are best correlated) is not returned.
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 2 /
      logical plot_flag
      integer *4 ncur,nlats_(0:ncur-1),n_S,ld_S,n_M,ld_M
      integer *4 n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type
      real *8 grid_k_p_(0:0),max_x_c,alpha_tru_(0:0),alpha_est_(0:0)
      real *8 eps_target,N_pixels_in
      complex *16 S_k_p_(0:0),M_k_p_(0:0)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,ns,nm
      integer *4, allocatable :: n_w_(:)
c$$$      temporary storage for templates and images
      complex *16, allocatable :: T_k_p_(:)
      real *8 pi
      real *8 max_r8_f
      real *8 max_k_c
c$$$      true displacement and rotation parameters for each image
      real *8, allocatable :: delta_x_tru_(:)
      real *8, allocatable :: delta_y_tru_(:)
      real *8, allocatable :: gamma_z_tru_(:)
c$$$      estimated (current) displacement and rotation parameters for each image
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
c$$$      updated displacement and rotation parameters for each image
      real *8, allocatable :: delta_x_upd_(:)
      real *8, allocatable :: delta_y_upd_(:)
      real *8, allocatable :: gamma_z_upd_(:)
c$$$      errors in displacement and rotation parameters for each image
      real *8, allocatable :: delta_x_err_(:)
      real *8, allocatable :: delta_y_err_(:)
      real *8, allocatable :: gamma_z_err_(:)
c$$$      range of displacement and rotation parameters to measure
      integer ndx,ndy,ngz
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)      
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16 C_S,C_M,C_Z
      real *8, allocatable :: delta_x_max_(:)
      real *8, allocatable :: delta_y_max_(:)
      real *8, allocatable :: gamma_z_max_(:)
      complex *16, allocatable :: C_Z_max_(:)
      real *8, allocatable :: delta_x_sort_(:)
      real *8, allocatable :: delta_y_sort_(:)
      real *8, allocatable :: gamma_z_sort_(:)
      complex *16, allocatable :: C_Z_sort_(:)
      integer *4, allocatable :: I_permute_(:)
      integer *4, allocatable :: I_inverse_(:)
      character(len=64) format_string
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
c$$$      function for template generation
      external get_F2_x_c
      real *8 param_1,param_2
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_batch_wrapper_0]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'nlats_: ',(nlats_(nr),nr=0,ncur-1)
         write(format_string,'(A,I0,A)') '(A,',ncur,'(F8.3,1X))'
         write(6,format_string) 'grid_k_p_: ',(grid_k_p_(nr),nr=0,ncur
     $        -1)
         write(6,'(A,F6.3)') 'max_x_c: ',max_x_c
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(format_string,'(A,I0,A)') '(',5,'(F8.3,1X))'
         write(6,'(A)') 'alpha_tru_: '
         write(6,format_string) (alpha_tru_(nr),nr=0,5*n_M-1)
         write(6,'(A)') 'alpha_est_: '
         write(6,format_string) (alpha_est_(nr),nr=0,5*n_M-1)
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,F6.3)') 'N_pixels_in: ',N_pixels_in
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
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

c$$$      reading in true translations and rotations
      allocate(delta_x_tru_(0:n_M-1))
      allocate(delta_y_tru_(0:n_M-1))
      allocate(gamma_z_tru_(0:n_M-1))
      do nm=0,n_M-1
         gamma_z_tru_(nm) = alpha_tru_(2 + nm*5)
         delta_x_tru_(nm) = alpha_tru_(3 + nm*5)
         delta_y_tru_(nm) = alpha_tru_(4 + nm*5)
         if (verbose.gt.1) then
            write(6,'(I0,A,3F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_tru_(nm),delta_y_tru_(nm),gamma_z_tru_(nm)
         end if
      enddo

c$$$      reading in estimated translations and rotations
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      do nm=0,n_M-1
         gamma_z_est_(nm) = alpha_est_(2 + nm*5)
         delta_x_est_(nm) = alpha_est_(3 + nm*5)
         delta_y_est_(nm) = alpha_est_(4 + nm*5)
         if (verbose.gt.1) then
            write(6,'(I0,A,3F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_est_(nm),delta_y_est_(nm),gamma_z_est_(nm)
         end if
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_svdread_batch'
     $        ,' to calculate all innerproducts:'
      end if
c$$$      ld_S is the stride between adjacent templates in S_k_q_
c$$$      ld_M is the stride between adjacent images in M_k_q_
c$$$      calculating all innerproducts
      allocate(delta_x_max_(0:n_S*n_M-1))
      allocate(delta_y_max_(0:n_S*n_M-1))
      allocate(gamma_z_max_(0:n_S*n_M-1))
      allocate(C_Z_max_(0:n_S*n_M-1))
      timing_tic = second()
      call test_innerproduct_svdread_batch(verbose,svd_calculation_type
     $     ,n_r,grid_k_p_,n_w_,max_x_c,n_S,ld_S,S_k_p_,n_M,ld_M,M_k_p_
     $     ,delta_x_est_,delta_y_est_,gamma_z_est_,eps_target
     $     ,N_pixels_in,n_delta_x,n_delta_y,n_gamma_z,delta_x_max_
     $     ,delta_y_max_,gamma_z_max_,C_Z_max_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') 'test_innerproduct_svdread_batch:'
     $        ,' total_time ',timing_toc-timing_tic
      end if
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_max_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_batch_sort'
     $        ,' to sort the results:'
      end if
      allocate(delta_x_sort_(0:n_S*n_M-1))
      allocate(delta_y_sort_(0:n_S*n_M-1))
      allocate(gamma_z_sort_(0:n_S*n_M-1))
      allocate(C_Z_sort_(0:n_S*n_M-1))
      allocate(I_permute_(0:n_S*n_M-1))
      allocate(I_inverse_(0:n_S*n_M-1))
      call test_innerproduct_batch_sort(n_S,n_M,delta_x_max_
     $     ,delta_y_max_,gamma_z_max_,C_Z_max_,delta_x_sort_
     $     ,delta_y_sort_,gamma_z_sort_,C_Z_sort_,I_permute_,I_inverse_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_permute_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_permute_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_inverse_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_inverse_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we update the estimated displacement'
     $        ,' and rotation parameters.'
      end if
      allocate(delta_x_upd_(0:n_M-1))
      allocate(delta_y_upd_(0:n_M-1))
      allocate(gamma_z_upd_(0:n_M-1))
      do nm=0,n_M-1
         delta_x_upd_(nm) = delta_x_est_(nm) + delta_x_sort_(n_S-1 + nm
     $        *n_S)
         delta_y_upd_(nm) = delta_y_est_(nm) + delta_y_sort_(n_S-1 + nm
     $        *n_S)
         gamma_z_upd_(nm) = gamma_z_est_(nm) + gamma_z_sort_(n_S-1 + nm
     $        *n_S)
      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*)
     $        'Now we compare updated estimate to true parameters: '
         write(format_string,'(A)')
     $        '(A,I2,A,F8.3,1X,F8.3,1X,F8.3,A,F8.3,1X,F8.3,1X,F8.3)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; upd_:' ,
     $           delta_x_upd_(nm) , delta_y_upd_(nm), gamma_z_ upd_(nm)
     $           , '; tru_:' , delta_x_tru_(nm) ,delta_y_tru_(nm) ,
     $           gamma_z_tru_(nm)
         enddo
      end if
      allocate(delta_x_err_(0:n_M-1))
      allocate(delta_y_err_(0:n_M-1))
      allocate(gamma_z_err_(0:n_M-1))
      do nm=0,n_M-1
         delta_x_err_(nm) = dabs(delta_x_tru_(nm)-delta_x_upd_(nm))
         delta_y_err_(nm) = dabs(delta_y_tru_(nm)-delta_y_upd_(nm))
         gamma_z_err_(nm) = dabs(gamma_z_tru_(nm)-gamma_z_upd_(nm))
         call periodize_r8(gamma_z_err_(nm),-pi,+pi,gamma_z_err_(nm))
      enddo
      if (verbose.gt.1) then
         write(6,*) 'Alignment Errors:'
         write(format_string,'(A)') '(A,I2,A,F8.4,A,F8.4,A,F8.4)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; delta_x_:' ,
     $           delta_x_err_(nm) , '; delta_y_:' , delta_y_err_(nm),
     $           '; gamma_z_:' , gamma_z_err_(nm)
         enddo
      end if
      if (verbose.gt.0) then
         write(6,*) 'Largest Alignment Errors:'
         write(format_string,'(A)') '(A,F8.4,1X,A,F8.4,1X,A,F8.4)'
         write(6,format_string) 
     $    'delta_x_err_:', max_r8_f(n_M,delta_x_err_), 
     $    'delta_y_err_:', max_r8_f(n_M,delta_y_err_), 
     $    'gamma_z_err_:', max_r8_f(n_M,gamma_z_err_)
      end if

      plot_flag = .false.
      if (verbose.gt.1) then
         plot_flag = .true.
      end if
      if (plot_flag) then
         continue
      else
         goto 10
      end if

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We initialize the first figure file.' 
      end if

      unitnumber = 7
      write(fname_pre,'(A)') 'test_innerproduct_batch_wrapper_0_pre'
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
     $        ,'test_innerproduct_batch_wrapper_0_pre.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We initialize the second figure file.' 
      end if
c$$$      Initializing figure file for image
      unitnumber = 7
      write(fname_pre,'(A)') 'test_innerproduct_batch_wrapper_0_pos'
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
         ns = I_permute_(n_S-1 + nm*n_S)
         delta_x = delta_x_upd_(nm)
         delta_y = delta_y_upd_(nm)
         gamma_z = gamma_z_upd_(nm)
         call cp1_c16(n_A,S_k_p_(ns*ld_S),T_k_p_)
         call transf_p_to_p(n_r,grid_k_p_,n_w_,n_A,T_k_p_,+delta_x,
     $        +delta_y,T_k_p_)
         call rotate_p_to_p(n_r,n_w_,n_A,T_k_p_,+gamma_z,T_k_p_)
         call stdlim_c16(n_A,.true.,T_k_p_,F_avg,F_std)
         F_max = F_avg + 2.5*F_std
         F_min = F_avg - 2.5*F_std
         x_offset = +5.0d0 + 10.0d0*nm
         y_offset = +17.0d0 - 0.0d0 + 22.5d0
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
         x_offset = +5.0d0 + 10.0d0*nm
         y_offset = +17.0d0 - 10.0d0 + 22.5d0
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
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: '
     $        ,'"gthumb test_innerproduct_batch_wrapper_0_pos.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

 10   continue ! if plot_flag.eqv..false.
      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_0]'
      end if

      end
