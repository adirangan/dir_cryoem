c     gfortran -w -o test_transforms_dr.out test_transforms_dr.f ; ./test_transforms_dr.out ; gthumb test_transforms_M.jpg; 
c     svn commit dir_rangan_playpen/ --username adirangan --password githubpass0 -m "testing transformations"
      program test_transforms_dr
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      dcfft
      integer *4 n_r,nr,n_w_max,n_A,n_A_W,ic_W
      integer *4, allocatable :: n_w_(:)
      integer *4, allocatable :: n_w_W_(:)
      complex *16, allocatable :: Wsave_(:)
c$$$      complex *16, allocatable :: Wdata_(:)
c$$$      grids 
      integer *4 n_gridpoints,ng
      real *8 pi
      real *8 max_x_c,max_k_c
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
      real *8 Flim_x_c,Flim_x_p,Flim_k_c,Flim_k_p,Flim_k_q
c$$$      Template S
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: S_x_p_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: S_k_q_(:)
c$$$      displacement and rotation parameters
      real *8 delta_x,delta_y,gamma
c$$$      Image M
      complex *16, allocatable :: M_x_c_(:)
      complex *16, allocatable :: M_x_p_(:)
      complex *16, allocatable :: M_k_c_(:)
      complex *16, allocatable :: M_k_p_(:)
      complex *16, allocatable :: M_k_q_(:)
c$$$      Transformed copy of Template T for comparison
      complex *16, allocatable :: T_x_c_(:)
      complex *16, allocatable :: T_x_p_(:)
      complex *16, allocatable :: T_k_c_(:)
      complex *16, allocatable :: T_k_p_(:)
      complex *16, allocatable :: T_k_q_(:)
c$$$      arrays hold Jtaylor expansion
      integer l_max,n_J
      complex *16, allocatable :: f_(:)
      integer, allocatable :: b_(:)
      integer, allocatable :: p_(:)
c$$$      arrays hold svd-expansion
      integer n_svd_r,n_svd_d,n_svd_l
      data n_svd_r / 512 /
      data n_svd_d / 256 /
      data n_svd_l / 34 /
      include 'svd_r_R075z088e1000'
      include 'svd_d_R075z088e1000'
      include 'svd_l_R075z088e1000'
      include 'svd_U_d_R075z088e1000'
      include 'svd_s_R075z088e1000'
      include 'svd_V_r_R075z088e1000'
c$$$      innerproducts
      complex *16 C_S_x_c,C_S_x_p,C_S_k_c,C_S_k_p,C_S_k_q
      complex *16 C_M_x_c,C_M_x_p,C_M_k_c,C_M_k_p,C_M_k_q
      complex *16 C_T_x_c,C_T_x_p,C_T_k_c,C_T_k_p,C_T_k_q
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      external get_F_x_c

c$$$      Calculating template size using 'get_template_size'
      pi = 4*atan(1.0)
      ncur = 24
      allocate(nlats_(0:ncur-1));
      do ng=0,ncur-1
         nlats_(ng) = 3 + 2*ng
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

c$$$      dcfft
      n_r = n_gridpoints
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         if (verbose.gt.3) then
            write(6,*) ' n_w_(nr): ',n_w_(nr)
         end if
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.3) then
         write(6,*) 'n_A ',n_A,'; n_w_max ',n_w_max
      end if
      allocate(n_w_W_(0:n_r-1))
      n_A_W = 0
      do nr=0,n_r-1
         n_w_W_(nr) = 4*n_w_(nr)+16
         n_A_W = n_A_W + n_w_W_(nr)
      enddo
      if (verbose.gt.3) then
         write(6,*) 'n_A_W ',n_A_W
      end if
      allocate(Wsave_(0:n_A_W-1))
      ic_w = 0
      do nr=0,n_r-1
         if (n_w_(nr).gt.0) then
            if (verbose.gt.3) then
               write(6,*) 'Calling dcffti with ', n_w_(nr)
            end if
            call dcffti(n_w_(nr),Wsave_(ic_w))
         end if
            ic_w = ic_w + n_w_W_(nr)
      enddo
      if (verbose.gt.3) then
         ic_w = 0
         do nr=0,n_r-1
            ic_w = ic_w + n_w_W_(nr)
            write(6,*) 'n_w_(',nr,')=',n_w_(nr),' --> n_w_W(',nr,')='
     $           ,n_w_W_(nr),'; ic_w = ',ic_w
         enddo
      end if
c$$$      allocate(Wdata_(0:n_w_max-1))
c$$$      call test_dcfft(n_r,n_w_,n_w_W_,n_A_W,Wsave_,n_A,Wdata_)

c$$$      Calculating Jtaylor expansion for use in transf_q_to_q
      l_max = 4
      n_J = ((l_max+1) * (l_max+2)) / 2
      allocate(f_(0:n_J-1))
      allocate(b_(0:n_J-1))
      allocate(p_(0:n_J-1))
      call get_Jtaylor(l_max,n_J,f_,b_,p_)

c$$$      Initializing figure file for template
      unitnumber = 7
      write(fname_pre,'(A)') 'test_transforms_S'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

c$$$      Calculating x-space template on regular cartesian-grid
      max_x_c = 1.0d0
      Flim_x_c = max_x_c*4.0
      allocate(grid_x_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c,n_gridpoints,grid_x_c_)
      allocate(S_x_c_(0:n_gridpoints*n_gridpoints-1))
      call get_F_x_c_(n_gridpoints,grid_x_c_,max_x_c,n_gridpoints
     $     ,grid_x_c_,max_x_c,S_x_c_,get_F_x_c)
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,S_x_c_,S_x_c_,C_S_x_c)
      C_S_x_c = zsqrt(C_S_x_c)
      write(6,'(A,2F16.3)') ' C_S_x_c: ',C_S_x_c
      F_max = Flim_x_c
      F_min = 0.0d0
      x_offset = 0.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,S_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'S_x_c'
      text_length = 5
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating x-space template on quasi-uniform polar-grid
      Flim_x_p = Flim_x_c
      allocate(S_x_p_(0:ntemplatesize-1))
      allocate(grid_x_p_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c/2.0,n_gridpoints,grid_x_p_)
      call interp_c_to_p(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,S_x_c_,n_gridpoints,grid_x_p_,ngridc_,ntemplatesize,S_x_p_)
      call innerproduct_p(n_gridpoints,grid_x_p_,ngridc_ ,ntemplatesize
     $     ,S_x_p_,S_x_p_,C_S_x_p)
      C_S_x_p = zsqrt(C_S_x_p)
      write(6,'(A,2F16.3)') ' C_S_x_p: ',C_S_x_p
      F_max = Flim_x_p
      F_min = 0.0d0
      x_offset = +5.0d0
      y_offset = -5.5d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_x_p_,ngridc_
     $     ,.true.,S_x_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'S_x_p'
      text_length = 5
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0
      y_loc = -0.5d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      
c$$$      Calculating k-space template on regular cartesian-grid      
      Flim_k_c = 1.0d0*n_gridpoints*n_gridpoints/8.0
      max_k_c = (1.0d0*n_gridpoints)/max_x_c
      allocate(grid_k_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c,n_gridpoints,grid_k_c_)
      allocate(S_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,S_x_c_,S_k_c_)
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,S_k_c_,S_k_c_,C_S_k_c)
      C_S_k_c = zsqrt(C_S_k_c)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_S_k_c: ',C_S_k_c
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +10.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,S_k_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'R(S_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +20.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.false.,S_k_c_,F_min,F_max,x_offset,y_offset
     $     ,x_side ,y_side)
      fig_title = 'I(S_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating k-space template on quasi-uniform polar-grid
      Flim_k_p = 1.0d0*n_gridpoints*n_gridpoints/8.0
      allocate(S_k_p_(0:ntemplatesize-1))
      allocate(grid_k_p_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c/2.0,n_gridpoints,grid_k_p_)
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,S_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,S_k_p_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_p_,S_k_p_,C_S_k_p)
      C_S_k_p = zsqrt(C_S_k_p)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S_k_p: ',C_S_k_p
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +15.0d0
      y_offset = -5.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'R(S_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +25.0d0
      y_offset = -5.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(S_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating J-space template on quasi-uniform polar-grid
      Flim_k_q = 1.0d0*n_gridpoints*n_gridpoints/8.0
      allocate(S_k_q_(0:ntemplatesize-1))
c$$$      call interp_p_to_q(n_gridpoints,ngridc_,ntemplatesize,S_k_p_
c$$$     $     ,S_k_q_)
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,S_k_p_,S_k_q_)
c$$$      call interp_p_to_q_dcfft2(n_gridpoints,ngridc_,ntemplatesize
c$$$     $     ,S_k_p_,S_k_q_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S_k_q)
      C_S_k_q = zsqrt(C_S_k_q)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S_k_q: ',C_S_k_q

      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +30.0d0
      y_offset = -0.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bessl(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      fig_title = 'R(S_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +30.0d0
      y_offset = -10.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bessl(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(S_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +40.0d0
      y_offset = -0.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      fig_title = 'R(S_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.0d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +40.0d0
      y_offset = -10.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,S_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(S_k_q)'
      text_length = 8
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

c$$$      Initializing figure file for image
      unitnumber = 7
      write(fname_pre,'(A)') 'test_transforms_M'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

c$$$      Calculating x-space image on regular cartesian-grid
      delta_x = +2.0/n_gridpoints*max_x_c
      delta_y = -1.0/n_gridpoints*max_x_c
      gamma = 3.0/n_gridpoints*2.0*pi
c$$$      delta_x = +0.050*max_x_c
c$$$      delta_y = -0.025*max_x_c
c$$$      gamma = 0.10*2.0*pi
      allocate(M_x_c_(0:n_gridpoints*n_gridpoints-1))
      call cp1_c16(n_gridpoints*n_gridpoints,S_x_c_,M_x_c_)
      call rotate_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+gamma,M_x_c_)
      call transl_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+delta_x,+delta_y,M_x_c_)
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,M_x_c_,M_x_c_,C_M_x_c)
      C_M_x_c = zsqrt(C_M_x_c)
      write(6,'(A,2F16.3)') ' C_M_x_c: ',C_M_x_c
      F_max = Flim_x_c
      F_min = 0.0d0
      x_offset = 0.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,M_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'M_x_c'
      text_length = 5
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = Flim_x_c
      F_min = 0.0d0
      x_offset = 10.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,M_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'M_x_c'
      text_length = 5
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      
c$$$      Comparing with transformed copy of template
      allocate(T_x_c_(0:n_gridpoints*n_gridpoints-1))
      call cp1_c16(n_gridpoints*n_gridpoints,S_x_c_,T_x_c_)
      call rotate_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,T_x_c_,+gamma,T_x_c_)
      call transl_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,T_x_c_,+delta_x,+delta_y,T_x_c_)
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,T_x_c_,T_x_c_,C_T_x_c)
      C_T_x_c = zsqrt(C_T_x_c)
      write(6,'(A,2F16.3)') ' C_T_x_c: ',C_T_x_c
      F_max = Flim_x_c
      F_min = 0.0d0
      x_offset = 0.0d0
      y_offset = -10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,T_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'T_x_c'
      text_length = 5
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = Flim_x_c
      F_min = 0.0d0
      x_offset = 10.0d0
      y_offset = -10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,S_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'S_x_c'
      text_length = 5
      color_index = 7
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating k-space image on regular cartesian-grid      
      allocate(M_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,M_x_c_,M_k_c_)
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,M_k_c_,M_k_c_,C_M_k_c)
      C_M_k_c = zsqrt(C_M_k_c)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_M_k_c: ',C_M_k_c
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +20.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,M_k_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'R(M_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +30.0d0
      y_offset = 0.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.false.,M_k_c_,F_min,F_max,x_offset,y_offset
     $     ,x_side ,y_side)
      fig_title = 'I(M_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Comparing with transformed k-space copy on regular cartesian-grid
      allocate(T_k_c_(0:n_gridpoints*n_gridpoints-1))
      call cp1_c16(n_gridpoints*n_gridpoints,S_k_c_,T_k_c_)
      call rotate_c_to_c(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,T_k_c_,+gamma,T_k_c_)
      call transf_c_to_c(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,T_k_c_,+delta_x,+delta_y,T_k_c_)
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,T_k_c_,T_k_c_,C_T_k_c)
      C_T_k_c = zsqrt(C_T_k_c)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_T_k_c: ',C_T_k_c
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +20.0d0
      y_offset = -10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.true.,T_k_c_,F_min,F_max,x_offset,y_offset,x_side
     $     ,y_side)
      fig_title = 'R(T_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_c
      F_min = -Flim_k_c
      x_offset = +30.0d0
      y_offset = -10.0d0
      x_side = 10.0d0
      y_side = 10.0d0
      call Fig_c16_carte(unitnumber,n_gridpoints,grid_x_c_,n_gridpoints
     $     ,grid_x_c_,.false.,T_k_c_,F_min,F_max,x_offset,y_offset
     $     ,x_side ,y_side)
      fig_title = 'I(T_k_c)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = +0.0d0
      y_loc = +0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating k-space image on quasi-uniform polar-grid
      allocate(M_k_p_(0:ntemplatesize-1))
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,M_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,M_k_p_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_p_,M_k_p_,C_M_k_p)
      C_M_k_p = zsqrt(C_M_k_p)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M_k_p: ',C_M_k_p
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +25.0d0
      y_offset = -15.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,M_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'R(M_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +35.0d0
      y_offset = -15.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,M_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(M_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Comparing with transformed k-space copy on quasi-uniform polar-grid
      allocate(T_k_p_(0:ntemplatesize-1))
      call cp1_c16(ntemplatesize,S_k_p_,T_k_p_)
      call rotate_p_to_p(n_gridpoints,ngridc_,ntemplatesize,T_k_p_,
     $     +gamma,T_k_p_)
      call transf_p_to_p(n_gridpoints,grid_k_p_,ngridc_,ntemplatesize
     $     ,T_k_p_,+delta_x,+delta_y,T_k_p_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,T_k_p_,T_k_p_,C_T_k_p)
      C_T_k_p = zsqrt(C_T_k_p)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_T_k_p: ',C_T_k_p
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +25.0d0
      y_offset = -25.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,T_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'R(T_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_p
      F_min = -Flim_k_p
      x_offset = +35.0d0
      y_offset = -25.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0/(max_k_c)
      call Fig_c16_polar(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,T_k_p_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(T_k_p)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = -0.5d0*(max_k_c)
      y_loc = -0.5d0*(max_k_c)
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Calculating J-space image on quasi-uniform polar-grid
      allocate(M_k_q_(0:ntemplatesize-1))
c$$$      call interp_p_to_q(n_gridpoints,ngridc_,ntemplatesize,M_k_p_
c$$$     $     ,M_k_q_)
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,M_k_p_,M_k_q_)
c$$$      call interp_p_to_q_dcfft2(n_gridpoints,ngridc_,ntemplatesize
c$$$     $     ,M_k_p_,M_k_q_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M_k_q)
      C_M_k_q = zsqrt(C_M_k_q)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M_k_q: ',C_M_k_q
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +0.0d0
      y_offset = -20.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,M_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      fig_title = 'R(M_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.5d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +10.0d0
      y_offset = -20.5d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,M_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(M_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.5d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Comparing with transformed J-space template on quasi-uniform polar-grid
      allocate(T_k_q_(0:ntemplatesize-1))
      call cp1_c16(ntemplatesize,S_k_q_,T_k_q_)
      call rotate_q_to_q(n_gridpoints,ngridc_,ntemplatesize,T_k_q_,
     $     +gamma,T_k_q_)
c$$$      call transf_taylor_q_to_q(l_max,n_J,f_,b_,p_,n_gridpoints
c$$$     $     ,grid_k_p_,ngridc_,ntemplatesize,T_k_q_,delta_x,delta_y
c$$$     $     ,T_k_q_)
      call transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_,n_svd_l
     $     ,svd_l_,svd_U_d_,svd_s_,svd_V_r_,n_gridpoints ,grid_k_p_
     $     ,ngridc_,ntemplatesize,T_k_q_,delta_x,delta_y ,T_k_q_)
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,T_k_q_,T_k_q_,C_T_k_q)
      C_T_k_q = zsqrt(C_T_k_q)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_T_k_q: ',C_T_k_q
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +0.0d0
      y_offset = -31.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.true.,T_k_q_,F_min,F_max,x_offset,y_offset,x_side ,y_side)
      fig_title = 'R(T_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.5d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)
      F_max = +Flim_k_q
      F_min = -Flim_k_q
      x_offset = +10.0d0
      y_offset = -31.0d0
      x_side = 10.0d0/(max_k_c)
      y_side = 10.0d0
      call Fig_c16_bess2(unitnumber,n_gridpoints,grid_k_p_,ngridc_
     $     ,.false.,T_k_q_,F_min,F_max,x_offset,y_offset,x_side,y_side)
      fig_title = 'I(T_k_q)'
      text_length = 8
      color_index = 0
      font_type = 0
      font_size = 72
      x_loc = 0.5d0*(max_k_c)
      y_loc = 0.0d0
      call Fig_text(unitnumber,text_length,fig_title,color_index
     $     ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $     ,y_side)

c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

 10   stop
      end

      include 'dfftpack.f'
      include 'test_dcfft.f'
      include 'get_template_size.f'
      include 'get_Jtaylor.f'
      include 'get_F_x_c.f'
      include 'get_F_x_c_.f'
      include 'stdlim_c16.f'
      include 'cp1_c16.f'
      include 'af1_c16.f'
      include 'linspace.f'
      include 'interp1_c16.f'
      include 'interp2_c16.f'
      include 'interp_c_to_p.f'
      include 'interp_p_to_q.f'
      include 'interp_p_to_q_dcfft.f'
      include 'interp_p_to_q_dcfft2.f'
      include 'periodize_r8.f'
      include 'periodize_i.f'
      include 'transl_c_to_c.f'
      include 'transf_c_to_c.f'
      include 'transf_p_to_p.f'
      include 'transf_taylor_q_to_q.f'
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
      include 'Fig_c16_polar.f'
      include 'Fig_c16_bessl.f'
      include 'Fig_c16_bess2.f'
      include 'adi_fft1.f'
      include 'adi_fft2.f'
