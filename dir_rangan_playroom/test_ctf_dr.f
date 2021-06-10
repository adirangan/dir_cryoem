c     gfortran -w -o test_ctf_dr.out test_ctf_dr.f -lfftw3; ./test_ctf_dr.out 0.25 0.75 10.0 10.0 5.0;
c     gthumb test_ctf_S.jpg; 
      program test_ctf_dr
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 0 /
c$$$      cryoem parameters
      integer ngridr,ityper
      real *8 rmax
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na
      integer *4, allocatable :: n_w_(:)
c$$$      grids 
      integer *4 n_gridpoints,ng
      real *8 pi
      real *8 max_x_c,max_k_c,dtmp
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
c$$$      Template S
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: S_x_p_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: S_k_q_(:)
c$$$      fftw
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
c$$$      innerproducts
      complex *16 C_S_x_c,C_S_x_p,C_S_k_c,C_S_k_p,C_S_k_q
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_avg,F_std,F_max,F_min
      real *8 x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      external get_ctf_x_c
      character(len=8) :: cmd_argstring
      real *8 cmd_d0,cmd_d1,cmd_d2,cmd_d3,cmd_d4
      call get_command_argument(1+0,cmd_argstring)
      read (cmd_argstring, '(F8.5)') cmd_d0
      write(6,'(A,F8.5)') 'param_0: ',cmd_d0
      call get_command_argument(1+1,cmd_argstring)
      read (cmd_argstring, '(F8.5)') cmd_d1
      write(6,'(A,F8.5)') 'param_1: ',cmd_d1
      call get_command_argument(1+2,cmd_argstring)
      read (cmd_argstring, '(F8.5)') cmd_d2
      write(6,'(A,F8.5)') 'param_2: ',cmd_d2
      call get_command_argument(1+3,cmd_argstring)
      read (cmd_argstring, '(F8.5)') cmd_d3
      write(6,'(A,F8.5)') 'param_3: ',cmd_d3
      call get_command_argument(1+4,cmd_argstring)
      read (cmd_argstring, '(F8.5)') cmd_d4
      write(6,'(A,F8.5)') 'param_4: ',cmd_d4

c$$$         We define xnodesr via 'getgridr'
      ngridr = 16 
      rmax = 1.0d0*ngridr
      ityper = 0
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)

c$$$      Calculating template size using 'get_template_size'
      pi = 4*atan(1.0)
      ncur = ngridr
      allocate(nlats_(0:ncur-1));
      do ng=0,ncur-1
         nlats_(ng) = nint(pi*xnodesr(1+ng))
         if (nlats_(ng).lt.6) nlats_(ng) = 6
      enddo
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
      end if
      n_gridpoints = ncur
      if (n_gridpoints.lt.2) then
         write(6,'(A,I0,A)') 'Error n_gridpoints',n_gridpoints,'<2'
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

c$$$      fftw
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_ESTIMATE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_ESTIMATE) 
         na = na + n_w_(nr)
      enddo

c$$$      Initializing figure file for template
      unitnumber = 7
      write(fname_pre,'(A)') './dir_jpg/test_ctf_S'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

c$$$      Calculating x-space template on regular cartesian-grid
      max_x_c = 1.0d0
      allocate(grid_x_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c,n_gridpoints,grid_x_c_)
      allocate(S_x_c_(0:n_gridpoints*n_gridpoints-1))
      call get_ctf_x_c_(n_gridpoints,grid_x_c_,max_x_c,n_gridpoints
     $     ,grid_x_c_,max_x_c,S_x_c_,get_ctf_x_c,cmd_d0,cmd_d1,cmd_d2
     $     ,cmd_d3,cmd_d4)
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,S_x_c_,S_x_c_,C_S_x_c)
      C_S_x_c = zsqrt(C_S_x_c)
      write(6,'(A,2F16.3)') ' C_S_x_c: ',C_S_x_c
      call stdlim_c16(n_gridpoints*n_gridpoints,.true.,S_x_c_,F_avg
     $     ,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = 0.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,S_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      write(fig_title,'(A,I4)') 'S_x_c_ ',nint(dlog(F_std))
      text_length = 11
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating x-space template on quasi-uniform polar-grid
      allocate(S_x_p_(0:ntemplatesize-1))
      allocate(grid_x_p_(0:n_gridpoints-1))
c$$$      instead of:
c$$$      call linspace(0.0d0,max_x_c/2.0,n_gridpoints,grid_x_p_)
c$$$      we use:
      dtmp = max_x_c/2.0/n_gridpoints
      call linspace(dtmp,max_x_c/2.0 + dtmp,n_gridpoints,grid_x_p_)
      call interp_c_to_p(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,S_x_c_,n_gridpoints,grid_x_p_,ngridc_,ntemplatesize,S_x_p_)
      call innerproduct_p(n_gridpoints,grid_x_p_,ngridc_ ,ntemplatesize
     $     ,S_x_p_,S_x_p_,C_S_x_p)
      C_S_x_p = zsqrt(C_S_x_p)
      write(6,'(A,2F16.3)') ' C_S_x_p: ',C_S_x_p
      call stdlim_c16(n_A,.true.,S_x_p_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +5.0d0
      y_offset = -5.5d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_x_p_,ngridc_
     $     ,.true.,S_x_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      write(fig_title,'(A,I4)') 'S_x_p_ ',nint(dlog(F_std))
      text_length = 11
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0
      y_loc = -0.5d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      
c$$$      Calculating k-space template on regular cartesian-grid      
      max_k_c = (1.0d0*n_gridpoints)/max_x_c
      allocate(grid_k_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c,n_gridpoints,grid_k_c_)
      allocate(S_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,S_x_c_,S_k_c_)
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,S_k_c_,S_k_c_,C_S_k_c)
      C_S_k_c = zsqrt(C_S_k_c)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_S_k_c: ',C_S_k_c
      call stdlim_c16(n_gridpoints*n_gridpoints,.true.,S_k_c_,F_avg
     $     ,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +10.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,S_k_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      write(fig_title,'(A,I4)') 'R(S_k_c_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      call stdlim_c16(n_gridpoints*n_gridpoints,.false.,S_k_c_,F_avg
     $     ,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +20.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.false.,S_k_c_,F_min,F_max,x_offset,y_offset
     $     ,x_side ,y_side)
      write(fig_title,'(A,I4)') 'I(S_k_c_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating k-space template on quasi-uniform polar-grid
      allocate(S_k_p_(0:ntemplatesize-1))
      allocate(grid_k_p_(0:n_gridpoints-1))
c$$$      instead of:
c$$$      call linspace(0.0d0,max_k_c/2.0,n_gridpoints,grid_k_p_)
c$$$      we use:
      dtmp = max_k_c/2.0/n_gridpoints
      call linspace(dtmp,max_k_c/2.0+dtmp,n_gridpoints,grid_k_p_)
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,S_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,S_k_p_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_p_,S_k_p_,C_S_k_p)
      C_S_k_p = zsqrt(C_S_k_p)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S_k_p: ',C_S_k_p
      call stdlim_c16(n_A,.true.,S_k_p_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +15.0d0
      y_offset = -5.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      write(fig_title,'(A,I4)') 'R(S_k_p_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      call stdlim_c16(n_A,.false.,S_k_p_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +25.0d0
      y_offset = -5.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      write(fig_title,'(A,I4)') 'I(S_k_p_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating J-space template on quasi-uniform polar-grid
      allocate(S_k_q_(0:ntemplatesize-1))
      call interp_p_to_q_fftw(n_gridpoints,fftw_plan_frwd_,n_w_,n_A
     $     ,fftw_in1_,fftw_out_,S_k_p_,S_k_q_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S_k_q)
      C_S_k_q = zsqrt(C_S_k_q)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S_k_q: ',C_S_k_q
      call stdlim_c16(n_A,.true.,S_k_q_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +30.0d0
      y_offset = -0.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bessl(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      write(fig_title,'(A,I4)') 'R(S_k_q_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      call stdlim_c16(n_A,.false.,S_k_q_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +30.0d0
      y_offset = -10.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bessl(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      write(fig_title,'(A,I4)') 'I(S_k_q_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      call stdlim_c16(n_A,.true.,S_k_q_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +40.0d0
      y_offset = -0.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      write(fig_title,'(A,I4)') 'R(S_k_q_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      call stdlim_c16(n_A,.false.,S_k_q_,F_avg,F_std)
      F_max = F_avg + 2.5*F_std
      F_min = F_avg - 2.5*F_std
      x_offset = +40.0d0
      y_offset = -10.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      write(fig_title,'(A,I4)') 'I(S_k_q_) ',nint(dlog(F_std))
      text_length = 14
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      
c$$$      Printing figure file for template
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)


 10   stop
      end

      include 'dfftpack.f'
      include 'getgridr.f'
      include 'get_template_size.f'
      include 'get_ctf_x_c.f'
      include 'get_ctf_x_c_.f'
      include 'stdlim_c16.f'
      include 'cp1_r8.f'
      include 'cp1_c16.f'
      include 'af1_c16.f'
      include 'linspace.f'
      include 'interp1_c16.f'
      include 'interp2_c16.f'
      include 'interp_c_to_p.f'
      include 'interp_p_to_q.f'
      include 'interp_p_to_q_fftw.f'
      include 'periodize_r8.f'
      include 'periodize_i.f'
      include 'transl_c_to_c.f'
      include 'transf_c_to_c.f'
      include 'transf_p_to_p.f'
      include 'transf_svd_q_to_q.f'
      include 'rotate_c_to_c.f'
      include 'rotate_p_to_p.f'
      include 'rotate_q_to_q.f'
      include 'innerproduct_c.f'
      include 'innerproduct_p.f'
      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Fig_header.f'
      include 'Fig_text.f'
      include 'recenter_c16.f'
      include 'Fig_c16_carte.f'
      include 'Fig_c16_cart2.f'
      include 'Fig_c16_polar.f'
      include 'Fig_c16_bessl.f'
      include 'Fig_c16_bess2.f'
      include 'adi_fft1.f'
      include 'adi_fft2.f'
