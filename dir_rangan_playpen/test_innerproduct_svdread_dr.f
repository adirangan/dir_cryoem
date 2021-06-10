c     gfortran -w -O2 -o test_innerproduct_svdread_dr.out test_innerproduct_svdread_dr.f -I ./dir_gen_Jsvd -lfftw3 ; ./test_innerproduct_svdread_dr.out 36 10.0 3.0 13 13 12 ;
c     gfortran -w -O2 -o test_innerproduct_svdread_dr.out test_innerproduct_svdread_dr.f -I ./dir_gen_Jsvd -lfftw3 ;
c     ./test_innerproduct_svdread_dr.out 36 10.0 3.0 13 13 12 ;
c     svn commit dir_rangan_playpen/ --username adirangan --password githubpass0 -m "testing innerproduct"
      program test_innerproduct_svdread_dr
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 1 /
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      dcfft
      integer *4 n_r,nr,n_w_max,n_A,ic,n_A_W,ic_W
      integer *4, allocatable :: n_w_(:)
      integer *4, allocatable :: n_w_W_(:)
      complex *16, allocatable :: Wsave_(:)
      complex *16, allocatable :: Wlast_(:)
c$$$      fftw
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16, allocatable :: fftw_in1_last_(:)
      complex *16, allocatable :: fftw_out_last_(:)
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
      real *8 delta_x_M,delta_y_M,gamma_z_M
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
c$$$      array of displacements to measure
      integer n_delta_x,n_delta_y,ndx,ndy
      real *8 eps_target,N_pixels,delta_x,delta_y,delta_r,delta_w
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: delta_r_(:)
      real *8, allocatable :: delta_w_(:)
      character(len=64) format_string
c$$$      array of angles to measure
      integer n_gamma_z,ngz
      real *8 gamma_z
      real *8, allocatable :: gamma_z_(:)
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16 C_S,C_M,C_Z,Y_k_q,Z_k_q
      complex *16, allocatable :: Y_k_q_(:)
      complex *16, allocatable :: Y_tmp_(:)
      complex *16, allocatable :: Z_k_q_(:)
      real *8 Y_k_q_x,Y_k_q_y
      real *8 Z_k_q_x,Z_k_q_y
      integer n_rho,nrho1,nrho2
      real *8, allocatable :: rho_(:)
      complex *16, allocatable :: rho_C_(:)
      complex *16 rho_C
      real *8 rho_D
c$$$      parameters for figure generation
      character(len=64) fig_title,fname_pre,fname_fig,fname_jpg
      integer unitnumber,text_length,color_index,font_type,font_size
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
      external get_F_x_c
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      parameters for command line input
c$$$      list: ncur N_pixels n_delta_x n_delta_y n_gamma_z
      character(len=8) :: cmd_argstring
      integer cmd_i1,cmd_i4,cmd_i5,cmd_i6
      real *8 cmd_d2,cmd_d3
      call get_command_argument(1,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i1
      call get_command_argument(2,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d2
      call get_command_argument(3,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d3
      call get_command_argument(4,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i4
      call get_command_argument(5,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i5
      call get_command_argument(6,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i6

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_svdread_dr]'
      end if

c$$$      Calculating template size using 'get_template_size'
      pi = 4*atan(1.0)
      ncur = cmd_i1
      allocate(nlats_(0:ncur-1));
      do ng=0,ncur-1
         nlats_(ng) = 3 + 2*ng
      enddo
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'ncur = ',ncur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      n_w_max = ngridc_(ncur)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
      end if
      n_gridpoints = ncur
      if (n_gridpoints.lt.2) then
         write(6,'(A,I0,A)') 'Error n_gridpoints',n_gridpoints,'<2'
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Here we construct the Wsave_ array.'
         write(6,*) 'This array holds the various workspaces for each'
     $        ,' call to dcfft (one for each nr).'
         write(6,'(A)') ''
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
      if (verbose.gt.1) then
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
      allocate(Wlast_(0:n_w_W_(n_r-1)-1))
      call dcffti(n_w_max,Wlast_)
      if (verbose.gt.3) then
         ic_w = 0
         do nr=0,n_r-1
            ic_w = ic_w + n_w_W_(nr)
            write(6,*) 'n_w_(',nr,')=',n_w_(nr),' --> n_w_W_(',nr,')='
     $           ,n_w_W_(nr),'; ic_w = ',ic_w
         enddo
      end if
c$$$      fftw
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      ic = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_in1_(ic),fftw_out_(ic),FFTW_FORWARD,FFTW_ESTIMATE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(ic),fftw_in1_(ic),FFTW_BACKWARD,FFTW_ESTIMATE) 
         ic = ic + n_w_(nr)
      enddo
      allocate(fftw_in1_last_(0:n_w_max-1))
      allocate(fftw_out_last_(0:n_w_max-1))
      call dfftw_plan_dft_1d_(fftw_plan_frwd_last,n_w_max
     $     ,fftw_in1_last_,fftw_out_last_,FFTW_FORWARD,FFTW_ESTIMATE) 
      call dfftw_plan_dft_1d_(fftw_plan_back_last,n_w_max
     $     ,fftw_out_last_,fftw_in1_last_,FFTW_BACKWARD,FFTW_ESTIMATE) 
      
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Here we construct the template "S".'
         write(6,*)
     $        'The array S_x_c_ holds the x-space values on a'
     $        ,' regular cartesian grid defined by grid_x_c_.'
         write(6,*)
     $        'The array S_k_c_ holds the k-space values on a'
     $        ,' regular cartesian grid defined by grid_k_c_.'
         write(6,*)
     $        'The array S_k_p_ holds the k-space values on a'
     $        ,' quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
         write(6,*) 'The array S_k_q_ holds the "hat-hat"'
     $        ,' (i.e, bessel-expansion) for the k-space values'
         write(6,*) ,'on a quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
         write(6,'(A)') ''
      end if
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
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,S_k_p_,S_k_q_)

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Now we do the same thing for the image "M".'
         write(6,*) 'This image will be constructed by taking'
     $        ,' the template,'
         write(6,*) 'translating it by delta_x_M and'
     $        ,' delta_y_M and then rotating' ,' it by gamma_z_M.'
         write(6,'(A)') ''
      end if
c$$$      Calculating x-space image on regular cartesian-grid
      delta_x_M = -1.5/n_gridpoints*max_x_c
      delta_y_M = +1.0/n_gridpoints*max_x_c
c$$$      gamma_z_M = +3.0/n_gridpoints*2.0*pi
      gamma_z_M = 1.0*pi/4.0;
      allocate(M_x_c_(0:n_gridpoints*n_gridpoints-1))
      call cp1_c16(n_gridpoints*n_gridpoints,S_x_c_,M_x_c_)
      call transl_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+delta_x_M,+delta_y_M,M_x_c_)
      call rotate_c_to_c(n_gridpoints,max_x_c,n_gridpoints,max_x_c
     $     ,M_x_c_,+gamma_z_M,M_x_c_)
c$$$      Calculating k-space image on regular cartesian-grid      
      allocate(M_k_c_(0:n_gridpoints*n_gridpoints-1))
      call adi_fft2(-1,n_gridpoints,n_gridpoints,M_x_c_,M_k_c_)
c$$$      Calculating k-space image on quasi-uniform polar-grid
      allocate(M_k_p_(0:ntemplatesize-1))
      call interp_c_to_p(n_gridpoints,max_k_c,n_gridpoints,max_k_c
     $     ,M_k_c_,n_gridpoints,grid_k_p_,ngridc_,ntemplatesize,M_k_p_)
c$$$      Calculating J-space image on quasi-uniform polar-grid
      allocate(M_k_q_(0:ntemplatesize-1))
      call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W,Wsave_
     $     ,ntemplatesize,M_k_p_,M_k_q_)

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*)
     $        'Now we set up an array of translations and rotations'
     $        ,' to serve as inputs to the calculation of'
     $        ,' innerproducts.'
         write(6,*) 'The translation array comprises'
     $        ,' n_delta_x by n_delta_y'
     $        ,' translations,'
         write(6,*) 'spanning +/- 2 pixels in each'
     $        ,' direction (i.e., a total area of 5x5 pixels).'
         write(6,*) 'The rotation array comprises'
     $        ,' n_gamma_z rotations uniformly distributed'
     $        ,' on [0,2*pi].'         
         write(6,'(A)') ''
      end if
c$$$      setting up array of displacements to measure
      eps_target = dmax1(0.0d0,cmd_d2)
      N_pixels = dmax1(0.0d0,cmd_d3)
      n_delta_x = max(1,cmd_i4)
      n_delta_y = max(1,cmd_i5)
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      allocate(delta_r_(0:n_delta_x*n_delta_y-1))
      allocate(delta_w_(0:n_delta_x*n_delta_y-1))
      do ndx=0,n_delta_x-1
         if (n_delta_x.gt.1) then
            delta_x = (-N_pixels + ndx*2*N_pixels/(n_delta_x-1))
     $           /n_gridpoints*max_x_c
         else
            delta_x = 0.0d0
         end if
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            if (n_delta_y.gt.1) then
               delta_y = (-N_pixels + ndy*2*N_pixels/(n_delta_y-1))
     $              /n_gridpoints*max_x_c
            else
               delta_y = 0.0d0
            end if
            delta_y_(ndy) = delta_y
         enddo
      enddo
c$$$      setting up array of angles to measure
      n_gamma_z = max(1,cmd_i6)
      allocate(gamma_z_(0:n_gamma_z-1))
      do ngz=0,n_gamma_z-1
         if (n_gamma_z.gt.1) then
            gamma_z = (2*pi*ngz)/n_gamma_z
         else
            gamma_z = 0.0d0
         end if
         gamma_z_(ngz) = gamma_z
      enddo
      if (verbose.gt.0) then
         write(6,'(A,I0,A,F5.2,A,F3.1,A,I0,A,I0,A,I0)') 'n_gridpoints: '
     $        ,n_gridpoints,'; eps_target: ',eps_target,'; N_pixels: '
     $        ,N_pixels,'; n_delta_x: ' ,n_delta_x,'; n_delta_y: '
     $        ,n_delta_y,'; n_gamma_z: ',n_gamma_z
      end if
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_delta_x,'(F5.3,1X))'
         write (6,format_string) 'delta_x_: ',(delta_x_(ndx),ndx=0
     $        ,n_delta_x-1)
         write(format_string,'(A,I0,A)') '(A,',n_delta_y,'(F5.3,1X))'
         write (6,format_string) 'delta_y_: ',(delta_y_(ndy),ndy=0
     $        ,n_delta_y-1)
         write(format_string,'(A,I0,A)') '(A,',n_gamma_z,'(F5.3,1X))'
         write (6,format_string) 'gamma_z_: ',(gamma_z_(ngz),ngz=0
     $        ,n_gamma_z-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we initialize the various arrays of '
     $        ,' innerproducts that we will calculate.'
         write(6,*) '(T_?_?_ serve as temporary storage'
     $        ,' for the translated/rotated templates)'
         write(6,*) 'The innerproduct Y_k_q_ also holds'
     $        ,' innerproducts calculated using the bessel-expansion,'
         write(6,*)
     $        ,' with the various rotations again handled using the'
     $        ,' fft.'
         write(6,*)
     $        '(translations are carried out on in polar-coordinates'
     $        ,' on S_k_p_ and M_k_p_ prior to calculation of'
     $        ,' bessel-expansion).'
         write(6,*) 'The innerproduct Z_k_q_ also holds'
     $        ,' innerproducts calculated using the bessel-expansion,'
         write(6,*)
     $        ,' with the various rotations again handled using the'
     $        ,' fft.'
         write(6,*)
     $        'This time, however, the translations are handled using'
     $        ,' an svd-expansion'
         write(6,'(A)') 'to separate each bessel-function.'
         write(6,'(A)') ''
      end if
c$$$      setting up arrays of innerproducts to hold measurements
      Clim = 1.0d0
      allocate(Y_k_q_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(Y_tmp_(0:n_w_max-1))
      allocate(Z_k_q_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(T_x_c_(0:n_gridpoints*n_gridpoints-1))
      allocate(T_k_c_(0:n_gridpoints*n_gridpoints-1))
      allocate(T_k_p_(0:ntemplatesize-1))
      allocate(T_k_q_(0:ntemplatesize-1))

      if (verbose.gt.1) then
         write(6,'(A)') 'We initialize the figure file.' 
      end if
c$$$      Initializing figure file for image
      unitnumber = 7
      write(fname_pre,'(A)') 'test_innerproduct_svdread'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

c$$$      Calculate innerproducts in J-space on quasi-uniform polar-grid
c$$$      once again, this time using k_only integral to facilitate
c$$$      calculation across rotations.  For this calculation we
c$$$      translate S_k_p (in k-space on polar-grid) before converting
c$$$      to J-space via interpolation. We deal with rotations using fft.
      if (verbose.gt.1) then
         write(6,'(A)') ' entering calculation of Y_k_q_ '
      end if
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      if (verbose.gt.1) then
         write(6,'(A,2F16.3)') ' C_S: ',C_S
      end if
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      if (verbose.gt.1) then
         write(6,'(A,2F16.3)') ' C_M: ',C_M
      end if
      C_Z = C_S*C_M
      timing_tic = second()
      do ndy=0,n_delta_y-1
         delta_y = delta_y_(ndy)
         do ndx=0,n_delta_x-1
            delta_x = delta_x_(ndx)
            if (verbose.gt.2) then
               write (6,'(A,F5.3,1X,F5.3,1X,F5.3)')
     $              'delta_x delta_y gamma_z: ',delta_x,delta_y,gamma_z
            end if
            call cp1_c16(ntemplatesize,S_k_p_,T_k_p_)
            call transf_p_to_p(n_gridpoints,grid_k_p_,ngridc_
     $           ,ntemplatesize,T_k_p_,+delta_x,+delta_y,T_k_p_)
c$$$            call interp_p_to_q_dcfft(n_gridpoints,ngridc_,n_w_W_,n_A_W
c$$$     $           ,Wsave_,ntemplatesize,T_k_p_,T_k_q_)
            call interp_p_to_q_fftw(n_gridpoints,fftw_plan_frwd_,ngridc_
     $           ,ntemplatesize,fftw_in1_,fftw_out_,T_k_p_,T_k_q_)
            call innerproduct_q__k_only(n_gridpoints,grid_k_p_,ngridc_
     $           ,ntemplatesize,T_k_q_,M_k_q_,Y_tmp_)
c$$$            call dcfftb(n_w_max,Y_tmp_,Wlast_)
            call cp1_c16(n_w_max,Y_tmp_,fftw_out_last_)
            call dfftw_execute_(fftw_plan_back_last)
            call cp1_c16(n_w_max,fftw_in1_last_,Y_tmp_)
            do ngz=0,n_gamma_z-1
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Y_tmp_,+gamma_z
     $              ,Y_k_q)
               nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               Y_k_q_(nC) = Y_k_q/(n_gridpoints**4)/C_Z
            enddo
         enddo
      enddo
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculation of Y_k_q_, total time ',timing_toc
     $        -timing_tic
      end if
      if (verbose.gt.1) then
         write(6,'(A)') 'Y_k_q_: '
         write(6,'(4(2F8.3,1X))') (Y_k_q_(nC),nC=0,n_delta_x*n_delta_y
     $        *n_gamma_z-1)
         write(6,'(A)') ' '
      end if

c$$$      Calculate innerproducts in J-space on quasi-uniform polar-grid
c$$$      once again, this time using svd-expansion coupled with k_only
c$$$      integral to facilitate calculation across rotations and
c$$$      translations simultaneously.
      if (verbose.gt.1) then
         write(6,'(A)') ' entering calculation of Z_k_q_ '
      end if
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,S_k_q_,S_k_q_,C_S)
      C_S = zsqrt(C_S)/(n_gridpoints*n_gridpoints)      
      if (verbose.gt.1) then
         write(6,'(A,2F16.3)') ' C_S: ',C_S
      end if
      call innerproduct_p(n_gridpoints,grid_k_p_,ngridc_ ,ntemplatesize
     $     ,M_k_q_,M_k_q_,C_M)
      C_M = zsqrt(C_M)/(n_gridpoints*n_gridpoints)      
      if (verbose.gt.1) then
         write(6,'(A,2F16.3)') ' C_M: ',C_M
      end if
      C_Z = C_S*C_M
      timing_tic = second()
      call innerproduct_q_delta_gamma_svdread(0,eps_target,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_gridpoints
     $     ,grid_k_p_,ngridc_,Wlast_,ntemplatesize,S_k_q_,M_k_q_,C_Z
     $     ,Z_k_q_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculation of Z_k_q_ (1), total time '
     $        ,timing_toc-timing_tic
      end if
      timing_tic = second()
      call innerproduct_q_delta_gamma_svdread(1,eps_target,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_gridpoints
     $     ,grid_k_p_,ngridc_,Wlast_,ntemplatesize,S_k_q_,M_k_q_,C_Z
     $     ,Z_k_q_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculation of Z_k_q_ (2), total time '
     $        ,timing_toc-timing_tic
      end if
      if (verbose.gt.1) then
         write(6,'(A)') 'Z_k_q_: '
         write(6,'(4(2F8.3,1X))') (Z_k_q_(nC),nC=0,n_delta_x*n_delta_y
     $        *n_gamma_z-1)
         write(6,'(A)') ' '
      end if

c$$$      Displaying pearson-correlations between various innerproducts.
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'We calculate the pearson-correlation'
     $        ,'between the various innerproduct arrays.'
         write(6,'(A)') ''
      end if
      n_A = n_delta_x*n_delta_y*n_gamma_z
      n_rho = 2
      allocate(rho_(0:n_rho*n_rho-1))
      allocate(rho_C_(0:n_A*n_rho-1))
      call cp1_c16(n_A,Y_k_q_,rho_C_(0*n_A))
      call cp1_c16(n_A,Z_k_q_,rho_C_(1*n_A))
      do nrho2=0,n_rho-1
         do nrho1=0,n_rho-1
            call pearson_c16(n_A,rho_C_,nrho1*n_A,1,rho_C_,nrho2*n_A,1
     $           ,rho_C)
            rho_(nrho1 + nrho2*n_rho) = zabs(rho_C)
         enddo
      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ' % Pearson-correlation: '
         write(6,'(A,A)') ' % % % %  Y_k_q_, Z_k_q_ '
         write(6,'(A,7(F8.3))') 'Y_k_q_ ',(rho_(nrho1),nrho1=0*n_rho,1
     $        *n_rho-1)
         write(6,'(A,7(F8.3))') 'Z_k_q_ ',(rho_(nrho1),nrho1=1*n_rho,2
     $        *n_rho-1)
      else if (verbose.gt.0) then
         write(6,'(A,F8.3)') ' Pearson-correlation: ',rho_(1)
      end if

c$$$      Arranging figure locations
      Y_k_q_x = - 1.0
      Y_k_q_y = - 0.0
      Z_k_q_x = + 1.0
      Z_k_q_y = - 0.0

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)')
     $        'We write the innerproduct arrays to figure file.' 
         write(6,*)
     $        'The outer angle indicates rotation angle gamma_z,'
     $        ,' each patch illustrates the translations'
     $        ,' delta_x and delta_y.'
      end if
c$$$      printing innerproducts
      do ngz=0,n_gamma_z-1
         gamma_z = gamma_z_(ngz)
         F_max = +1.0*Clim
         F_min = -0.0*Clim
         x_side = 10.0d0
         y_side = 10.0d0
         x_offset = +7.5d0*cos(gamma_z) + Y_k_q_x*x_side/12.5
         y_offset = +7.5d0*sin(gamma_z) + Y_k_q_y*y_side/12.5
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,Y_k_q_(ngz*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
         x_offset = +7.5d0*cos(gamma_z) + Z_k_q_x*x_side/12.5
         y_offset = +7.5d0*sin(gamma_z) + Z_k_q_y*y_side/12.5
         call Fig_c16_cart2(unitnumber,n_delta_x,delta_x_ ,n_delta_y
     $        ,delta_y_,.true.,Z_k_q_(ngz*n_delta_x*n_delta_y),F_min
     $        ,F_max ,x_offset,y_offset ,x_side,y_side)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 'We write the labels to the figure file.' 
         write(6,*) 'The various innerproduct arrays are indexed'
     $        ,' according to their location.'
      end if
c$$$      adding labels
         text_length = 5
         color_index = 0
         font_type = 0
         font_size = 36
         x_loc = +0.0d0
         y_loc = +0.0d0
         x_offset = 3.0*(+ Y_k_q_x*0.30 ) - 0.5
         y_offset = 3.0*(+ Y_k_q_y*0.15 )
         fig_title = 'Y_k_q'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         x_offset = 3.0*(+ Z_k_q_x*0.30 ) - 0.5
         y_offset = 3.0*(+ Z_k_q_y*0.15 )
         fig_title = 'Z_k_q'
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)

      if (verbose.gt.1) then
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A)') 'try:
     $"gthumb test_innerproduct_svdread.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_svdread_dr]'
      end if

      stop
      end

      include 'dfftpack.f'
      include 'test_dcfft.f'
      include 'get_template_size.f'
      include 'get_F_x_c.f'
      include 'get_F_x_c_.f'
      include 'stdlim_c16.f'
      include 'pearson_c16.f'
      include 'cp1_c16.f'
      include 'cps_c16.f'
      include 'af1_c16.f'
      include 'afs_c16.f'
      include 'linspace.f'
      include 'interp1_c16.f'
      include 'interp2_c16.f'
      include 'interp_c_to_p.f'
      include 'interp_p_to_q.f'
      include 'interp_p_to_q_dcfft.f'
      include 'interp_p_to_q_dcfft2.f'
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
      include 'innerproduct_q__k_only.f'
      include 'innerproduct_q__k_svds.f'
      include 'innerproduct_q__k_svdd.f'
      include 'innerproduct_q_delta_gamma.f'
      include 'innerproduct_q_delta_gamma_svdread.f'
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
