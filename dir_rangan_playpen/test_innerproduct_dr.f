c     gfortran -w -o test_innerproduct_dr.out test_innerproduct_dr.f ; ./test_innerproduct_dr.out ; gthumb test_innerproduct_C_x_c.jpg
c     svn commit dir_rangan_playpen/ --username adirangan --password githubpass0 -m "testing innerproduct"
      program test_innerproduct_dr
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
c$$$      displacement and rotation parameters for image
      real *8 delta_x_M,delta_y_M,gamma_M
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
c$$$      array of displacements to measure
      integer n_delta_x,n_delta_y,ndx,ndy
      real *8 delta_x,delta_y,delta_r,delta_w
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: delta_r_(:)
      real *8, allocatable :: delta_w_(:)
c$$$      array of angles to measure
      integer n_gamma,ngamma
      real *8 gamma
      real *8, allocatable :: gamma_(:)
c$$$      array of innerproducts to hold measurements
      integer n_C,nC,nC2,nC3,nC4,nC5,nw,l
      real *8 Clim
      complex *16 C_S,C_M,C_Z,C_x_c,C_k_c,C_k_p,C_k_q,Y_k_q,Z_k_q
      complex *16, allocatable :: C_x_c_(:)
      complex *16, allocatable :: C_k_c_(:)
      complex *16, allocatable :: C_k_p_(:)
      complex *16, allocatable :: C_k_q_(:)
      complex *16, allocatable :: Y_k_q_(:)
      complex *16, allocatable :: Y_tmp_(:)
      complex *16, allocatable :: Z_k_q_(:)
      complex *16, allocatable :: Z_svds_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmp_(:)
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      external get_F_x_c
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

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
      n_w_max = ngridc_(ncur)
      write(6,'(A,I0,A,I0)') 'n_w_max = ',n_w_max,'; ntemplatesize = '
     $     ,ntemplatesize
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
      if (verbose.gt.-1) then
         write(6,*) 'n_w_max ',n_w_max,'; n_A ',n_A
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
            write(6,*) 'n_w_(',nr,')=',n_w_(nr),' --> n_w_W_(',nr,')='
     $           ,n_w_W_(nr),'; ic_w = ',ic_w
         enddo
      end if
      
c$$$      Calculating Jtaylor expansion for use in transf_q_to_q
      l_max = 4
      n_J = ((l_max+1) * (l_max+2)) / 2
      allocate(f_(0:n_J-1))
      allocate(b_(0:n_J-1))
      allocate(p_(0:n_J-1))
      call get_Jtaylor(l_max,n_J,f_,b_,p_)

c$$$      Calculating x-space template on regular cartesian-grid
      max_x_c = 1.0d0
      Flim_x_c = max_x_c*4.0
      allocate(grid_x_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c,n_gridpoints,grid_x_c_)
      allocate(S_x_c_(0:n_gridpoints*n_gridpoints-1))
      call get_F_x_c_(n_gridpoints,grid_x_c_,max_x_c,n_gridpoints
     $     ,grid_x_c_,max_x_c,S_x_c_,get_F_x_c)
c$$$      Calculating x-space template on quasi-uniform polar-grid
      Flim_x_p = Flim_x_c
      allocate(S_x_p_(0:ntemplatesize-1))
      allocate(grid_x_p_(0:n_gridpoints-1))
      call linspace(0.0d0,max_x_c/2.0,n_gridpoints,grid_x_p_)
      call interp_c_to_p(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,S_x_c_,n_gridpoints,grid_x_p_,ngridc_,ntemplatesize,S_x_p_)      
c$$$      Calculating k-space template on regular cartesian-grid      
      Flim_k_c = 1.0d0*n_gridpoints*n_gridpoints/8.0
      max_k_c = (1.0d0*n_gridpoints)/max_x_c
      allocate(grid_k_c_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c,n_gridpoints,grid_k_c_)
      allocate(S_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,S_x_c_,S_k_c_)
c$$$      Calculating k-space template on quasi-uniform polar-grid
      Flim_k_p = 1.0d0*n_gridpoints*n_gridpoints/8.0
      allocate(S_k_p_(0:ntemplatesize-1))
      allocate(grid_k_p_(0:n_gridpoints-1))
      call linspace(0.0d0,max_k_c/2.0,n_gridpoints,grid_k_p_)
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,S_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,S_k_p_)
c$$$      Calculating J-space template on quasi-uniform polar-grid
      Flim_k_q = 1.0d0*n_gridpoints*n_gridpoints/8.0
      allocate(S_k_q_(0:ntemplatesize-1))
c$$$      call interp_p_to_q(n_gridpoints,ngridc_,ntemplatesize,S_k_p_
c$$$     $     ,S_k_q_)
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,S_k_p_,S_k_q_)

c$$$      Calculating x-space image on regular cartesian-grid
      delta_x_M = -1.5/n_gridpoints*max_x_c
      delta_y_M = +1.0/n_gridpoints*max_x_c
      gamma_M = -6.0/n_gridpoints*2.0*pi
      allocate(M_x_c_(0:n_gridpoints*n_gridpoints-1))
      call cp1_c16(n_gridpoints*n_gridpoints,S_x_c_,M_x_c_)
      call transl_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+delta_x_M,+delta_y_M,M_x_c_)
      call rotate_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+gamma_M,M_x_c_)
c$$$      Calculating k-space image on regular cartesian-grid      
      allocate(M_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,M_x_c_,M_k_c_)
c$$$      Calculating k-space image on quasi-uniform polar-grid
      allocate(M_k_p_(0:ntemplatesize-1))
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,M_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,M_k_p_)
c$$$      Calculating J-space image on quasi-uniform polar-grid
      allocate(M_k_q_(0:ntemplatesize-1))
c$$$      call interp_p_to_q(n_gridpoints,ngridc_,ntemplatesize,M_k_p_
c$$$     $     ,M_k_q_)
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,M_k_p_,M_k_q_)

c$$$      setting up array of displacements to measure
      n_delta_x = 7
      n_delta_y = 7
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      allocate(delta_r_(0:n_delta_x*n_delta_y-1))
      allocate(delta_w_(0:n_delta_x*n_delta_y-1))
      do ndx=0,n_delta_x-1
         delta_x = (-1.0 + ndx*2.0/(n_delta_x-1))/n_gridpoints*max_x_c
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            delta_y = (-1.0 + ndy*2.0/(n_delta_y-1))/n_gridpoints
     $           *max_x_c
            delta_y_(ndy) = delta_y
            delta_r_(ndx + ndy*n_delta_x) = dsqrt(delta_x**2+delta_y**2)
            delta_w_(ndx + ndy*n_delta_x) = atan2(delta_y,delta_x)
         enddo
      enddo
c$$$      setting up array of angles to measure
      n_gamma = n_gridpoints/3
      allocate(gamma_(0:n_gamma-1))
      do ngamma=0,n_gamma-1
         gamma = (2*pi*ngamma)/n_gamma
         gamma_(ngamma) = gamma
      enddo
c$$$      setting up arrays of innerproducts to hold measurements
      Clim = 1.0d0
      allocate(C_x_c_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(C_k_c_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(C_k_p_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(C_k_q_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(Y_k_q_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(Y_tmp_(0:n_w_max-1))
      allocate(Z_k_q_(0:n_delta_x*n_delta_y*n_gamma-1))
      allocate(Z_svds_(0:n_svd_l*n_w_max-1))
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      allocate(Z_tmp_(0:n_w_max*n_delta_x*n_delta_y-1))
      allocate(T_x_c_(0:n_gridpoints*n_gridpoints-1))
      allocate(T_k_c_(0:n_gridpoints*n_gridpoints-1))
      allocate(T_k_p_(0:ntemplatesize-1))
      allocate(T_k_q_(0:ntemplatesize-1))

c$$$      Initializing figure file for image
      unitnumber = 7
      write(fname_pre,'(A)') 'test_innerproduct_C_x_c'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

c$$$      calculate innerproducts in x-space on regular cartesian grid
      write(6,'(A)') ' entering calculating C_x_c_ '
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,S_x_c_,S_x_c_,C_S)
      C_S = zsqrt(C_S)
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints ,grid_x_c_
     $     ,M_x_c_,M_x_c_,C_M)
      C_M = zsqrt(C_M)
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      nC = 0
      do ngamma=0,n_gamma-1
         gamma = gamma_(ngamma)
         do ndy=0,n_delta_y-1
            delta_y = delta_y_(ndy)
            do ndx=0,n_delta_x-1
               delta_x = delta_x_(ndx)
               call cp1_c16(n_gridpoints*n_gridpoints,S_x_c_,T_x_c_)
               call transl_c_to_c(n_gridpoints,max_x_c,n_gridpoints
     $              ,max_x_c,T_x_c_,+delta_x,+delta_y,T_x_c_)
               call rotate_c_to_c(n_gridpoints,max_x_c,n_gridpoints
     $              ,max_x_c,T_x_c_,+gamma,T_x_c_)
               call innerproduct_c(n_gridpoints,grid_x_c_,n_gridpoints
     $              ,grid_x_c_,T_x_c_,M_x_c_,C_x_c)
c$$$               nC = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               C_x_c_(nC) = C_x_c/C_Z
               nC = nC + 1
           enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating C_x_c_, total time '
     $     ,timing_toc-timing_tic

c$$$      calculate innerproducts in k-space on regular cartesian grid
      write(6,'(A)') ' entering calculating C_k_c_ '
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,S_k_c_,S_k_c_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints ,grid_k_c_
     $     ,M_k_c_,M_k_c_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      nC = 0
      do ngamma=0,n_gamma-1
         gamma = gamma_(ngamma)
         do ndy=0,n_delta_y-1
            delta_y = delta_y_(ndy)
            do ndx=0,n_delta_x-1
               delta_x = delta_x_(ndx)
               call cp1_c16(n_gridpoints*n_gridpoints,S_k_c_,T_k_c_)
               call transf_c_to_c(n_gridpoints,max_k_c,n_gridpoints
     $              ,max_k_c,T_k_c_,+delta_x,+delta_y,T_k_c_)
               call rotate_c_to_c(n_gridpoints,max_k_c,n_gridpoints
     $              ,max_k_c,T_k_c_,+gamma,T_k_c_)
               call innerproduct_c(n_gridpoints,grid_k_c_,n_gridpoints
     $              ,grid_k_c_,T_k_c_,M_k_c_,C_k_c)
c$$$               nC = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               C_k_c_(nC) = C_k_c/(n_gridpoints**4)/C_Z
               nC = nC + 1
           enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating C_k_c_, total time '
     $     ,timing_toc-timing_tic

c$$$      calculate innerproducts in k-space on quasi-uniform polar-grid
      write(6,'(A)') ' entering calculating C_k_p_ '
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_p_,S_k_p_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_p_,M_k_p_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      nC = 0
      do ngamma=0,n_gamma-1
         gamma = gamma_(ngamma)
         do ndy=0,n_delta_y-1
            delta_y = delta_y_(ndy)
            do ndx=0,n_delta_x-1
               delta_x = delta_x_(ndx)
               call cp1_c16(ntemplatesize,S_k_p_,T_k_p_)
               call transf_p_to_p(n_gridpoints,grid_k_p_,ngridc_
     $              ,ntemplatesize,T_k_p_,+delta_x,+delta_y,T_k_p_)
               call rotate_p_to_p(n_gridpoints,ngridc_,ntemplatesize
     $              ,T_k_p_,+gamma,T_k_p_)
               call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_
     $              ,ntemplatesize,T_k_p_,M_k_p_,C_k_p)
c$$$               nC = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               C_k_p_(nC) = C_k_p/(n_gridpoints**4)/C_Z
               nC = nC + 1
           enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating C_k_p_, total time '
     $     ,timing_toc-timing_tic

c$$$      calculate innerproducts in J-space on quasi-uniform polar-grid
      write(6,'(A)') ' entering calculating C_k_q_ '
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      nC = 0
      do ngamma=0,n_gamma-1
         gamma = gamma_(ngamma)
         do ndy=0,n_delta_y-1
            delta_y = delta_y_(ndy)
            do ndx=0,n_delta_x-1
               delta_x = delta_x_(ndx)
               call cp1_c16(ntemplatesize,S_k_q_,T_k_q_)
c$$$               call transf_taylor_q_to_q(l_max,n_J,f_,b_,p_,n_gridpoints
c$$$     $              ,grid_k_p_,ngridc_,ntemplatesize,T_k_q_,delta_x
c$$$     $              ,delta_y ,T_k_q_)
               call transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_
     $              ,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_
     $              ,n_gridpoints ,grid_k_p_,ngridc_,ntemplatesize
     $              ,T_k_q_,+delta_x,+delta_y ,T_k_q_)
               call rotate_q_to_q(n_gridpoints,ngridc_,ntemplatesize
     $              ,T_k_q_,+gamma,T_k_q_)
               call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_
     $              ,ntemplatesize,T_k_q_,M_k_q_,C_k_q)
c$$$               nC = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               C_k_q_(nC) = C_k_q/(n_gridpoints**4)/C_Z
               nC = nC + 1
           enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating C_k_q_, total time '
     $     ,timing_toc-timing_tic

c$$$      calculate innerproducts in J-space on quasi-uniform polar-grid
c$$$      once again, this time using k_only integral to facilitate
c$$$      calculation across rotations.
      write(6,'(A)') ' entering calculating Y_k_q_ '
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      do ndy=0,n_delta_y-1
         delta_y = delta_y_(ndy)
         do ndx=0,n_delta_x-1
            delta_x = delta_x_(ndx)
c$$$            At this point we can either translate S_k_q_ directly, via:
c$$$            call cp1_c16(ntemplatesize,S_k_q_,T_k_q_)
c$$$            call transf_svd_q_to_q(n_svd_r,svd_r_,n_svd_d,svd_d_
c$$$     $           ,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_
c$$$     $           ,n_gridpoints ,grid_k_p_,ngridc_,ntemplatesize
c$$$     $           ,T_k_q_,+delta_x,+delta_y ,T_k_q_)
c$$$            Or we can translate S_k_p_ and then interpolate to form T_k_q_:
            call cp1_c16(ntemplatesize,S_k_p_,T_k_p_)
            call transf_p_to_p(n_gridpoints,grid_k_p_,ngridc_
     $           ,ntemplatesize,T_k_p_,+delta_x,+delta_y,T_k_p_)
c$$$            call interp_p_to_q(n_gridpoints,ngridc_,ntemplatesize,T_k_p_
c$$$     $           ,T_k_q_)
            call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W
     $           ,Wsave_,ntemplatesize,T_k_p_,T_k_q_)
c$$$            Either way, we deal with rotations using fft:
            call innerproduct_q__k_only(n_gridpoints,grid_k_p_,ngridc_
     $           ,ntemplatesize,T_k_q_,M_k_q_,Y_tmp_)
c$$$            call adi_fft1(+1,n_w_max,Y_tmp_,0,1,Y_tmp_,0,1)
c$$$            call affine_c16(n_w_max,dcmplx((1.0d0*n_w_max),0.0)
c$$$     $           ,dcmplx(0.0,0.0),n_w_max,Y_tmp_,0,1,n_w_max,Y_tmp_,0,1)
            call dcfftb(n_w_max,Y_tmp_,Wsave_(n_A_W-n_w_W_(n_r-1)))
            do ngamma=0,n_gamma-1
               gamma = gamma_(ngamma)
               call interp1_c16(n_w_max,0.0d0,2*pi,Y_tmp_,+gamma,Y_k_q)
               nC = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               Y_k_q_(nC) = Y_k_q/(n_gridpoints**4)/C_Z
            enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating Y_k_q_, total time '
     $     ,timing_toc-timing_tic

c$$$      calculate innerproducts in J-space on quasi-uniform polar-grid
c$$$      once again, this time using svd-expansion coupled with k_only
c$$$      integral to facilitate calculation across rotations and
c$$$      translations simultaneously.
      write(6,'(A)') ' entering calculating Z_k_q_ '
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_S: ',C_S
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      write(6,'(A,2F16.3)') ' C_M: ',C_M
      C_Z = C_S*C_M
      timing_tic = second()
      call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_
     $     ,svd_V_r_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,S_k_q_
     $     ,M_k_q_,Z_svds_)
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      nC = 0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
c$$$            nC4 = ndx*n_w_max + ndy*n_w_max*n_delta_x
            nC4 = nC
            do nw=0,n_w_max-1
c$$$               nC = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
               Z_tmp_(nC) = cmplx( 0.0 , 0.0 )
               nC2 = ndx*n_svd_l + ndy*n_svd_l*n_delta_x
               nC3 = nw*n_svd_l
               do l=0,n_svd_l-1
                  Z_tmp_(nC) = Z_tmp_(nC) + Z_svdd_(nC2)*Z_svds_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo
               nC = nC + 1
            enddo
c$$$            call adi_fft1(+1,n_w_max,Z_tmp_,nC4,1,Z_tmp_,nC4 ,1)
c$$$            call affine_c16(n_w_max,dcmplx((1.0d0*n_w_max),0.0)
c$$$     $           ,dcmplx(0.0,0.0),n_w_max,Z_tmp_,nC4,1,n_w_max,Z_tmp_
c$$$     $           ,nC4,1)
            call dcfftb(n_w_max,Z_tmp_(nC4),Wsave_(n_A_W-n_w_W_(n_r-1)))
            do ngamma=0,n_gamma-1
               gamma = gamma_(ngamma)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_(nC4),+gamma
     $              ,Z_k_q)
               nC5 = ndx + ndy*n_delta_x + ngamma*n_delta_x*n_delta_y
               Z_k_q_(nC5) = Z_k_q/(n_gridpoints**4)/C_Z
            enddo
         enddo
      enddo
      timing_toc = second()
      write(6,'(A,F6.3)') ' finished calculating Z_k_q_, total time '
     $     ,timing_toc-timing_tic

c$$$      printing innerproducts
      do ngamma=0,n_gamma-1
         gamma = gamma_(ngamma)
         F_max = +1.0*Clim
         F_min = -0.0*Clim
         x_side = 10.0d0
         y_side = 10.0d0
         x_offset = +5.0d0*cos(gamma) - x_side/20.0
         y_offset = +5.0d0*sin(gamma) + 3*y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,C_x_c_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +5.0d0*cos(gamma) - x_side/20.0
         y_offset = +5.0d0*sin(gamma) + y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,C_k_c_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +5.0d0*cos(gamma) - x_side/20.0
         y_offset = +5.0d0*sin(gamma) - y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,C_k_p_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +5.0d0*cos(gamma) + x_side/20.0
         y_offset = +5.0d0*sin(gamma) - y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,C_k_q_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +5.0d0*cos(gamma) + x_side/20.0
         y_offset = +5.0d0*sin(gamma) + y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,Y_k_q_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +5.0d0*cos(gamma) + x_side/20.0
         y_offset = +5.0d0*sin(gamma) + 3*y_side/20.0
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,Z_k_q_(ngamma*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
      enddo

c$$$      adding labels
         text_length = 5
         color_index = 0
         font_type = 0
         font_size = 36
         x_loc = +0.0d0
         y_loc = +0.0d0
         x_offset = 3.0*(- 0.3 )
         y_offset = 3.0*(+ 3*0.15 )
         fig_title = 'x_c'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(- 0.3 )
         y_offset = 3.0*(+ 0.15 )
         fig_title = 'k_c'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(- 0.3 )
         y_offset = 3.0*(- 0.15 )
         fig_title = 'k_p'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(+ 0.3 )
         y_offset = 3.0*(- 0.15 )
         fig_title = 'k_q'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(+ 0.3 )
         y_offset = 3.0*(+ 0.15 )
         fig_title = 'Y_q'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(+ 0.3 )
         y_offset = 3.0*(+ 3*0.15 )
         fig_title = 'Z_q'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)

c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

      stop
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
      include 'innerproduct_q__k_only.f'
      include 'innerproduct_q__k_svds.f'
      include 'innerproduct_q__k_svdd.f'
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
