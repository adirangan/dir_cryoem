c     Compile with: gfortran -fcray-pointer -w -O2 -o test_innerproduct_svdread_batch_dr.out test_innerproduct_svdread_batch_dr.f -lfftw3 ; 
c     Run with: ./test_innerproduct_svdread_batch_dr.out 1 36 1.0 3.0 13 13 12 -1 ; 
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      program test_innerproduct_svdread_batch_dr
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      This program is a driver which tests the subroutine
c$$$      test_innerproduct_svdread_batch_dr.
c$$$      The main steps of this program are:
c$$$      A. First we generate a set of n_S = 5 templates.
c$$$         Each template looks a little different from the next
c$$$         (see test_innerproduct_svdread_batch_pre.jpg).
c$$$      B. Then we generate a set of n_M = 15 images.
c$$$         The first three images are copies of the first template,
c$$$         but with a different alignment (i.e., translated and 
c$$$         rotated slightly).
c$$$         The second three images are copies of the second template,
c$$$         again with a different alignment, and so forth.
c$$$      C. Then we call test_innerproduct_svdread_batch to calculate
c$$$         the innerproducts between each pair of templates and
c$$$         images. These innerproducts are calculated across as range
c$$$         of displacements delta_x_,delta_y_ as well as a range of 
c$$$         rotations gamma_z_ (all passed in as input).
c$$$      D. The output of test_innerproduct_svdread_batch is used
c$$$         to determine (among other things) which templates correlate
c$$$         most strongly with each image, and what the best
c$$$         alignment (i.e., displacement and rotation) is for each pair
c$$$         (see test_innerproduct_svdread_batch_pos.jpg).
c$$$      E. The true alignment of each image is then compared with
c$$$         the best empirically observed alignment, and the alignment-
c$$$         errors are printed to stdout.
c$$$      F. Finally, if (verbose.gt.1), we generate two figure files:
c$$$         test_innerproduct_svdread_batch_pre.jpg:
c$$$              S_x_c_(ns) shows template ns in x-space cartesian grid.
c$$$              R(S_k_p_(ns)) shows real part of fourier-transform on
c$$$                quasi-uniform grid in k-space.
c$$$              I(S_k_p_(ns)) shows imag part of fourier-transform on
c$$$                quasi-uniform grid in k-space.
c$$$              M_x_c_(nm) shows image nm in x-space cartesian grid.
c$$$         test_innerproduct_svdread_batch_pre.jpg:
c$$$              M_x_c_(nm) shows image nm in x-space cartesian grid.
c$$$              Above each M_x_c_(nm) we display T(S_x_c_(ns)),
c$$$              where ns is the template-index which corresponds most
c$$$              closely to image nm. In each case the template S_x_c_(ns)
c$$$              has been transformed (i.e., producing T(S_x_c_(ns))) 
c$$$              using the best fitting alignment (i.e., delta_x,delta_y
c$$$              and gamma_z) associated with that image-template pair.
c$$$
c$$$      Note regarding alignment:
c$$$      Assume that we are given an image M and template S.
c$$$      The innerproduct between M and S (with no alignment shift) is
c$$$      \[ \int M \cdot \bar{S} \].
c$$$      If we were to consider the alignment shift associated with
c$$$      a displacement (delta_x,delta_y) and an (in-plane) rotation
c$$$      gamma_z, then the innerproduct would be:
c$$$      \[ \int M \cdot \bar{T(S)} \], 
c$$$      where T(S) is the transformed template obtained by *first*
c$$$      translating S by (+delta_x,+delta_y), and *then* rotating S
c$$$      by (+gamma_z).
c$$$      This is equivalent to the innerproduct
c$$$      \[ \int T(M) \cdot \bar{S} \],
c$$$      where T(M) is the transformed image obtained by *first*
c$$$      rotating M by (-gamma_z) and *then* translating M by
c$$$      (-delta_x,-delta_y).
c$$$
c$$$      Note regarding naming-conventions:
c$$$      We use the following naming conventions for templates and images:
c$$$      S_x_c_: template in real-space on a cartesian-grid.
c$$$      S_x_p_: template in real-space on a quasi-uniform polar-grid.
c$$$      S_k_c_: template in fourier-space on a cartesian-grid.
c$$$      S_k_p_: template in fourier-space on a quasi-uniform polar-grid.
c$$$      S_k_q_: template in fourier-space, represented using "hat-hat" 
c$$$              coordinates (i.e., bessel-expansion).
c$$$      M_x_c_ through M_k_q_ are defined similarly.
c$$$      grid_x_c_: grid-values of real-space cartesian-grid
c$$$      grid_k_c_: grid-values of fourier-space cartesian-grid
c$$$      grid_k_p_: radial-values of rings of fourier-space polar-grid
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$
c$$$      INPUTS: This driver requires 8 parameters as command-line inputs.
c$$$              These inputs should be listed in the order below:
c$$$
c$$$      integer verbose       Determines verbosity level.
c$$$                          - verbose.eq.1 prints out only basic 
c$$$                            information and generates no figures.
c$$$                          - verbose.eq.2 prints out pre- and post-sorted
c$$$                            alignment arrays and generates figures.
c$$$
c$$$      integer ncur          Number of rings to pass into 
c$$$                            get_template_size. This parameter is also
c$$$                            used to determine the number of gridpoints 
c$$$                            (i.e., pixels) on each side of the templates
c$$$                            and images in real space (as measured on a
c$$$                            cartesian-grid). 
c$$$                          - Values of ncur~10 --> low resolution. 
c$$$                          - Values of ncur~50 --> moderate resolution.
c$$$                          - Values of ncur~100 --> high resolution.
c$$$                          - We have not generated libraries for values
c$$$                            of ncur.gt.192 yet.
c$$$
c$$$      real *8 eps_target    The epsilon used as a cutoff when 
c$$$                            determining which terms to retain within the
c$$$                            svd-expansion of the various bessel-
c$$$                            functions used to approximate the term
c$$$                            exp(2*pi*K*delta*cos(psi-omega)) within
c$$$                            the fourier-space translation operator.
c$$$                          - A value of eps_target.eq.10.0d0 corresponds
c$$$                            to low accuracy, and will probably not give
c$$$                            exactly the same result as an exact brute-
c$$$                            force calculation 
c$$$                            (see test_innerproduct_svdread_dr.f for a 
c$$$                            direct comparison).
c$$$                          - A value of eps_target.eq.1.0d0 corresponds
c$$$                            to moderate accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to a few digits, which is 
c$$$                            usually enough for the best alignment to 
c$$$                            match exactly.
c$$$                          - A value of eps_target.eq.0.1d0 corresponds
c$$$                            to high accuracy, and will usually match
c$$$                            the exact brute-force calculation of each
c$$$                            innerproduct to several digits, which is 
c$$$                            more than enough to match the best alignment
c$$$                            exactly.
c$$$                          - Values of eps_target.lt.0.1d0 are not
c$$$                            usually necessary, and we have not generated
c$$$                            the appropriate libraries for these 
c$$$                            eps_target yet.
c$$$                            
c$$$      real *8 N_pixels_in   Number of pixels in each direction 
c$$$                            associated with the largest displacement in 
c$$$                            the array of displacements used when 
c$$$                            calculating innerproducts for
c$$$                            each pair of templates and images.
c$$$                          - Each pixel of displacement is measured in
c$$$                            x-space. Thus, the maximum displacement is:
c$$$                            delta_x = N_pixels_in/n_r * max_x_c
c$$$                            delta_y = N_pixels_in/n_r * max_x_c.
c$$$                          - For example, N_pixels_in.eq.1.0d0 will limit
c$$$                            the array of displacements to a single pixel
c$$$                            in each direction (i.e., a 3x3 grid of 
c$$$                            pixels centered at the origin).
c$$$                          - A value of N_pixels_in.eq.3.0d0 will 
c$$$                            correspond to an array of displacements
c$$$                            spanning up to 3 pixels in each direction
c$$$                            (i.e., a 7x7 grid of pixels centered at the
c$$$                            origin).
c$$$                          - Values of N_pixels_in.gt.3.0 are not usually
c$$$                            necessary, and we have not generated the
c$$$                            appropriate libraries for these N_pixels_in
c$$$                            yet.
c$$$
c$$$      integer n_delta_x     Number of displacements to search for in the 
c$$$                            x-direction. These displacements will be 
c$$$                            linearly spaced from -N_pixels_in to 
c$$$                            +N_pixels_in.
c$$$                          - A value of n_delta_x = 1 + 2*N_pixels_in
c$$$                            implies a displacement-vector associated 
c$$$                            with each pixel-center in the grid 
c$$$                            determined by N_pixels_in.
c$$$                          - A value of n_delta_x = 1 + 4*N_pixels_in
c$$$                            implies a displacement-vector associated 
c$$$                            with each half-pixel (i.e., pixel-center
c$$$                            and pixel-side) in the grid determined by
c$$$                            N_pixels_in.
c$$$                          - Note that the total number of displacement-
c$$$                            vectors is n_delta_x * n_delta_y.
c$$$
c$$$      integer n_delta_y     Number of displacements to search for in the 
c$$$                            y-direction. These displacements will be 
c$$$                            linearly spaced from -N_pixels_in to 
c$$$                            +N_pixels_in.
c$$$                          - Note that the total number of displacement-
c$$$                            vectors is n_delta_x * n_delta_y.
c$$$
c$$$      integer n_gamma_z     Number of rotations to search for. These
c$$$                            rotations will be linearly spaced from 0.0d0
c$$$                            to 2*pi. It is usually reasonable for
c$$$                            n_gamma_z to be on the same order as ncur
c$$$                            (e.g., n_gamma_z = ncur/2 or so).
c$$$
c$$$      integer svd_calculation_type Determines the type of calculation
c$$$                            to use when computing innerproducts.
c$$$                          - A value of svd_calculation_type.eq.0
c$$$                            forces the displacement operation to take
c$$$                            place first, and only later uses the fft
c$$$                            to account for the various rotations.
c$$$                            This calculation is useful when the number
c$$$                            of distinct displacements is smaller than
c$$$                            the number of terms in the svd-expansion.
c$$$                          - A value of svd_calculation_type.eq.1
c$$$                            forces the displacement operation to take
c$$$                            place after the fft. This calculation is
c$$$                            useful when the number of terms in the
c$$$                            svd-expansion is smaller than the number
c$$$                            of distinct displacements.
c$$$                          - A value of svd_calculation_type.eq.2
c$$$                            forces a brute-force calculation of each
c$$$                            displacement separately (using the fft
c$$$                            to handle rotations).
c$$$                          - A value of svd_calculation_type.lt.0
c$$$                            allows test_innerproduct_svdread_batch.f
c$$$                            to automatically choose one of the three
c$$$                            calculation_types based on the other
c$$$                            input parameters (e.g., by comparing
c$$$                            n_svd_l with n_delta_x*n_delta_y).
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      data verbose / 2 /
      logical plot_flag
      integer *4 ntemplatesize,ncur
      integer *4, allocatable :: nlats_(:)
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na
      integer *4, allocatable :: n_w_(:)
c$$$      fftw
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*)
      complex *16 fftw_out_last_(*)
c$$$      grids
      real *8 pi
      real *8 max_r8
      real *8 max_x_c,max_k_c
      real *8, allocatable :: grid_x_c_(:)
      real *8, allocatable :: grid_x_p_(:)
      real *8, allocatable :: grid_k_c_(:)
      real *8, allocatable :: grid_k_p_(:)
      real *8 Flim_x_c,Flim_x_p,Flim_k_c,Flim_k_p,Flim_k_q
c$$$      Templates S_
      integer n_S,ld_S,ns
      complex *16, allocatable :: S_x_c_(:)
      complex *16, allocatable :: T_x_c_(:)
      complex *16, allocatable :: S_x_p_(:)
      complex *16, allocatable :: S_k_c_(:)
      complex *16, allocatable :: S_k_p_(:)
      complex *16, allocatable :: S_k_q_(:)
c$$$      Images M_
      integer n_M,ld_M,nm
      complex *16, allocatable :: M_x_c_(:)
      complex *16, allocatable :: M_x_p_(:)
      complex *16, allocatable :: M_k_c_(:)
      complex *16, allocatable :: M_k_p_(:)
      complex *16, allocatable :: M_k_q_(:)
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
      integer n_delta_x,n_delta_y,ndx,ndy
      real *8 eps_target,N_pixels_in,delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      integer n_gamma_z,ngz
      real *8, allocatable :: gamma_z_(:)      
      integer svd_calculation_type
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
      real *8 F_max,F_min,x_loc,y_loc,x_side,y_side,x_offset,y_offset
      character(len=1024) fig2dev_system_call
      integer system,system_error
c$$$      function for template generation
      external get_F2_x_c
      real *8 param_1,param_2
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      parameters for command line input
c$$$      list: verbose ncur eps_target N_pixels_in n_delta_x n_delta_y n_gamma_z svd_calculation_type
      character(len=8) :: cmd_argstring
      integer cmd_i0,cmd_i1,cmd_i4,cmd_i5,cmd_i6,cmd_i7
      real *8 cmd_d2,cmd_d3
      call get_command_argument(1+0,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i0
      call get_command_argument(1+1,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i1
      call get_command_argument(1+2,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d2
      call get_command_argument(1+3,cmd_argstring)
      read (cmd_argstring, '(F8.0)') cmd_d3
      call get_command_argument(1+4,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i4
      call get_command_argument(1+5,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i5
      call get_command_argument(1+6,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i6
      call get_command_argument(1+7,cmd_argstring)
      read (cmd_argstring, '(I10)') cmd_i7

      verbose = cmd_i0
      ncur = cmd_i1
      eps_target = cmd_d2
      N_pixels_in = cmd_d3
      n_delta_x = cmd_i4
      n_delta_y = cmd_i5
      n_gamma_z = cmd_i6
      svd_calculation_type = cmd_i7
      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_svdread_batch_dr]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,F6.3)') 'N_pixels_in: ',N_pixels_in
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
      end if

      pi = 4*atan(1.0)

c$$$      Calculating template size using 'get_template_size'
      allocate(nlats_(0:ncur-1));
      do nr=0,ncur-1
c$$$         nlats_(nr) = 3 + 2*nr
c$$$         052017: updating with nlats calculation from testrebuild.f:
         nlats_(nr) = nint(pi*(nr+1))
      enddo
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
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))
      
c$$$      generating grids for templates and images 
      max_x_c = 1.0d0
      Flim_x_c = max_x_c*4.0
      allocate(grid_x_c_(0:n_r-1))
      call linspace(0.0d0,max_x_c,n_r,grid_x_c_)
      Flim_x_p = Flim_x_c
      allocate(grid_x_p_(0:n_r-1))
      call linspace(0.0d0,max_x_c/2.0,n_r,grid_x_p_)
      max_k_c = (1.0d0*n_r)/max_x_c
      Flim_k_c = 1.0d0*n_r*n_r/8.0
      allocate(grid_k_c_(0:n_r-1))
      call linspace(0.0d0,max_k_c,n_r,grid_k_c_)
      Flim_k_p = 1.0d0*n_r*n_r/8.0
      allocate(grid_k_p_(0:n_r-1))
      call linspace(0.0d0,max_k_c/2.0,n_r,grid_k_p_)
      Flim_k_q = 1.0d0*n_r*n_r/8.0

c$$$      generate templates
      n_S = 5
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)') 'Here we construct the templates "S_".'
         write(6,'(A,I0,A)') 'We generate n_S = ',n_S
     $        ,' templates in total,'
         write(6,*) 'each using the function get_F2_x_c'
     $        ,' with different parameters.'
         write(6,*) 'For each template: '
         write(6,*)
     $        ' The array S_x_c_ holds the x-space values on a'
     $        ,' regular cartesian grid defined by grid_x_c_.'
         write(6,*)
     $        ' The array S_k_c_ holds the k-space values on a'
     $        ,' regular cartesian grid defined by grid_k_c_.'
         write(6,*)
     $        ' The array S_k_p_ holds the k-space values on a'
     $        ,' quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
         write(6,*) ' The array S_k_q_ holds the "hat-hat"'
     $        ,' (i.e, bessel-expansion) for the k-space values'
         write(6,*) 'on a quasi-uniform polar-grid defined by'
     $        ,' ncur and get_template_size.'
      end if
      allocate(S_x_c_(0:n_r*n_r*n_S-1))
      allocate(T_x_c_(0:n_r*n_r))
      allocate(S_x_p_(0:n_A*n_S-1))
      allocate(S_k_c_(0:n_r*n_r*n_S-1))
      allocate(S_k_p_(0:n_A*n_S-1))
      allocate(S_k_q_(0:n_A*n_S-1))
      do ns=0,n_S-1
         param_1 = ns*dsqrt(2.5d0)
         call periodize_r8(param_1,-1.0d0,1.0d0,param_1)
         param_2 = ns+2
c$$$      Calculating x-space template on regular cartesian-grid
         call get_F2_x_c_(n_r,grid_x_c_,max_x_c,n_r ,grid_x_c_,max_x_c
     $        ,S_x_c_(ns*n_r*n_r),get_F2_x_c,param_1,param_2)
c$$$      Calculating x-space template on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_x_c,n_r,max_x_c ,S_x_c_(ns*n_r*n_r)
     $        ,n_r,grid_x_p_,n_w_,n_A,S_x_p_(ns*n_A))
c$$$      Calculating k-space template on regular cartesian-grid
         call adi_fft2(-1,n_r,n_r,S_x_c_(ns*n_r*n_r),S_k_c_(ns*n_r*n_r))
c$$$      Calculating k-space template on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_k_c,n_r,max_k_c ,S_k_c_(ns*n_r*n_r)
     $        ,n_r,grid_k_p_,n_w_,n_A,S_k_p_(ns*n_A))
c$$$      Calculating J-space template on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,S_k_p_(ns*n_A),S_k_q_(ns*n_A))
      enddo

      n_M = 3*n_S
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A,I0,A)') 'Now we do the same thing for the n_M = '
     $        ,n_M,' images "M_".'
         write(6,*) 'The j^th image will be constructed by taking'
     $        ,' the (j/3)^th template,'
         write(6,*) 'translating it by delta_x_tru_(j) and'
     $        ,' delta_y_tru_(j) and then rotating'
     $        ,' it by gamma_z_tru_(j).'
         write(6,*) 'For this example we pick delta_x_tru_(j)'
     $        ,' to be (j%3)-1,',' and delta_y_tru_(j) to be'
     $        ,' 2*((j%3)-1).'
         write(6,*) 'We also pick gamma_z_tru_(j) to be'
     $        ,' 2*pi/3*(j%3).'
      end if
      allocate(M_x_c_(0:n_r*n_r*n_M-1))
      allocate(M_k_c_(0:n_r*n_r*n_M-1))
      allocate(M_k_p_(0:n_A*n_M-1))
      allocate(M_k_q_(0:n_A*n_M-1))
c$$$      Calculating x-space image on regular cartesian-grid
      allocate(delta_x_tru_(0:n_M-1))
      allocate(delta_y_tru_(0:n_M-1))
      allocate(gamma_z_tru_(0:n_M-1))
      do nm=0,n_M-1
         ns = nm/3
         delta_x_tru_(nm) = 1.0*(mod(nm,3)-1)/n_r*max_x_c
         delta_y_tru_(nm) = 2.0*(mod(nm,3)-1)/n_r*max_x_c
         gamma_z_tru_(nm) = mod(nm,3)*2*pi/3.0
         if (verbose.gt.1) then
            write(6,'(I0,A,3F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_tru_(nm),delta_y_tru_(nm),gamma_z_tru_(nm)
         end if
         call cp1_c16(n_r*n_r,S_x_c_(ns*n_r*n_r),M_x_c_(nm*n_r*n_r))
         call transl_c_to_c(n_r,max_x_c,n_r,max_x_c,M_x_c_(nm*n_r*n_r),
     $        +delta_x_tru_(nm),+delta_y_tru_(nm),M_x_c_(nm*n_r*n_r))
         call rotate_c_to_c(n_r,max_x_c,n_r,max_x_c,M_x_c_(nm*n_r*n_r),
     $        +gamma_z_tru_(nm),M_x_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on regular cartesian-grid      
         call adi_fft2(-1,n_r,n_r,M_x_c_(nm*n_r*n_r),M_k_c_(nm*n_r*n_r))
c$$$      Calculating k-space image on quasi-uniform polar-grid
         call interp_c_to_p(n_r,max_k_c,n_r,max_k_c,M_k_c_(nm*n_r*n_r)
     $        ,n_r,grid_k_p_,n_w_,n_A,M_k_p_(nm*n_A))
c$$$      Calculating J-space image on quasi-uniform polar-grid
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,M_k_p_(nm*n_A),M_k_q_(nm*n_A))
      enddo

c$$$c$$$      The following (commented-out) code sets up a simpler test
c$$$c$$$      where the estimated alignments correspond to the true ones
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ''
c$$$         write(6,*) 'Now we set up the estimated translations'
c$$$     $        ,' and rotations for each image.'
c$$$         write(6,*) 'As a first test we choose these estimated'
c$$$     $        ,' translations and rotations to be the same as the'
c$$$     $        ,' true ones.'
c$$$      end if
c$$$      allocate(delta_x_est_(0:n_M-1))
c$$$      allocate(delta_y_est_(0:n_M-1))
c$$$      allocate(gamma_z_est_(0:n_M-1))
c$$$      do nm=0,n_M-1
c$$$         delta_x_est_(nm) = 1.0*(mod(nm,3)-1)/n_r*max_x_c
c$$$         delta_y_est_(nm) = 2.0*(mod(nm,3)-1)/n_r*max_x_c
c$$$         gamma_z_est_(nm) = mod(nm,3)*2*pi/3.0
c$$$      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we set up the estimated translations'
     $        ,' and rotations for each image.'
         write(6,*) 'As a second test we choose these estimated'
     $        ,' translations and rotations to be close to the'
     $        ,' true ones (but not exactly the same).'
      end if
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      do nm=0,n_M-1
         delta_x_est_(nm) = delta_x_tru_(nm)*(1 + 0.125*(mod(nm+0,3)-1))
         delta_y_est_(nm) = delta_y_tru_(nm)*(1 + 0.125*(mod(nm+1,3)-1))
         gamma_z_est_(nm) = gamma_z_tru_(nm)*(1 + 0.125*(mod(nm+2,3)-1))
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_svdread_batch'
     $        ,' to calculate all innerproducts:'
      end if
c$$$      ld_S is the stride between adjacent templates in S_k_q_
      ld_S = n_A
c$$$      ld_M is the stride between adjacent images in M_k_q_
      ld_M = n_A
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
         write(format_string,'(A)') '(A,I2,A,F8.6,A,F8.6,A,F8.6)'
         do nm=0,n_M-1
            write(6,format_string) 'nm: ' , nm , '; delta_x_:' ,
     $           delta_x_err_(nm) , '; delta_y_:' , delta_y_err_(nm),
     $           '; gamma_z_:' , gamma_z_err_(nm)
         enddo
      end if
      if (verbose.gt.0) then
         write(6,*) 'Largest Alignment Errors:'
         write(format_string,'(A)') '(A,F8.6,1X,A,F8.6,1X,A,F8.6)'
         write(6,format_string) 
     $    'delta_x_err_:', max_r8(n_M,delta_x_err_), 
     $    'delta_y_err_:', max_r8(n_M,delta_y_err_), 
     $    'gamma_z_err_:', max_r8(n_M,gamma_z_err_)
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
      write(fname_pre,'(A)') 'test_innerproduct_svdread_batch_pre'
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
         F_max = Flim_x_c
         F_min = 0.0d0
         x_offset = 0.0d0 + 10.0d0*ns
         y_offset = 0.0d0
         x_side = 10.0d0
         y_side = 10.0d0
         call Fig_c16_carte(unitnumber,n_r,grid_x_c_,n_r,grid_x_c_
     $        ,.true.,S_x_c_(ns*n_r*n_r),F_min,F_max,x_offset,y_offset
     $        ,x_side,y_side)
         write(fig_title,'(A,I0,A)') 'S_x_c(',ns,')'
         text_length = 8
         color_index = 7
         font_type = 0
         font_size = 72
         x_loc = +0.0d0
         y_loc = +0.0d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         F_max = +Flim_k_p
         F_min = -Flim_k_p
         x_offset = +5.0d0 + 10.0d0*ns
         y_offset = -5.0d0 - 0.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
     $        ,S_k_p_(ns*n_A),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A)') 'R(S_k_p(',ns,'))'
         text_length = 11
         color_index = 0
         font_type = 0
         font_size = 72
         x_loc = -0.5d0*(max_k_c)
         y_loc = -0.5d0*(max_k_c)
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         F_max = +Flim_k_p
         F_min = -Flim_k_p
         x_offset = +5.0d0 + 10.0d0*ns
         y_offset = -5.0d0 - 10.0d0
         x_side = 10.0d0/(max_k_c)
         y_side = 10.0d0/(max_k_c)
         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
     $        ,S_k_p_(ns*n_A),F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A)') 'I(S_k_p(',ns,'))'
         text_length = 11
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
         F_max = Flim_x_c
         F_min = 0.0d0
         x_offset = 0.0d0 + 10.0d0*(nm/3)
         y_offset = 0.0d0 + 12.0d0*(1 + mod(nm,3))
         x_side = 10.0d0
         y_side = 10.0d0
         call Fig_c16_carte(unitnumber,n_r,grid_x_c_,n_r,grid_x_c_
     $        ,.true.,M_x_c_(nm*n_r*n_r),F_min,F_max,x_offset,y_offset
     $        ,x_side,y_side)
         write(fig_title,'(A,I2,A)') 'M_x_c(',nm,')'
         text_length = 9
         color_index = 7
         font_type = 0
         font_size = 72
         x_loc = +0.0d0
         y_loc = +0.0d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
c$$$c$$$         The following (commented-out) code adds extra subplots 
c$$$c$$$         showing the real and imaginary parts of M_k_p_
c$$$         F_max = +Flim_k_p
c$$$         F_min = -Flim_k_p
c$$$         x_offset = +5.0d0 + 10.0d0*(nm/3)
c$$$         y_offset = +27.0d0 + 10.0d0*(1 + mod(nm,3))
c$$$         x_side = 10.0d0/(max_k_c)
c$$$         y_side = 10.0d0/(max_k_c)
c$$$         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.true.
c$$$     $        ,M_k_p_(nm*n_A),F_min,F_max,x_offset,y_offset,x_side
c$$$     $        ,y_side)
c$$$         write(fig_title,'(A,I2,A)') 'R(M_k_p(',nm,'))'
c$$$         text_length = 12
c$$$         color_index = 0
c$$$         font_type = 0
c$$$         font_size = 72
c$$$         x_loc = -0.5d0*(max_k_c)
c$$$         y_loc = -0.5d0*(max_k_c)
c$$$         call Fig_text(unitnumber,text_length,fig_title,color_index
c$$$     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
c$$$     $        ,y_side)
c$$$         F_max = +Flim_k_p
c$$$         F_min = -Flim_k_p
c$$$         x_offset = +5.0d0 + 10.0d0*(nm/3)
c$$$         y_offset = +27.0d0 + 0.0d0*(1 + mod(nm,3))
c$$$         x_side = 10.0d0/(max_k_c)
c$$$         y_side = 10.0d0/(max_k_c)
c$$$         call Fig_c16_polar(unitnumber,n_r,grid_k_p_,n_w_,.false.
c$$$     $        ,M_k_p_(nm*n_A),F_min,F_max,x_offset,y_offset,x_side
c$$$     $        ,y_side)
c$$$         write(fig_title,'(A,I2,A)') 'I(M_k_p(',nm,'))'
c$$$         text_length = 12
c$$$         color_index = 0
c$$$         font_type = 0
c$$$         font_size = 72
c$$$         x_loc = -0.5d0*(max_k_c)
c$$$         y_loc = -0.5d0*(max_k_c)
c$$$         call Fig_text(unitnumber,text_length,fig_title,color_index
c$$$     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
c$$$     $        ,y_side)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: "gthumb '
     $        ,'test_innerproduct_svdread_batch_pre.jpg"'
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
      write(fname_pre,'(A)') 'test_innerproduct_svdread_batch_pos'
      write(fname_fig,'(A,A)') trim(fname_pre),'.fig'
      write(fname_jpg,'(A,A)') trim(fname_pre),'.jpg'
      open(unitnumber,file=fname_fig,status='replace' ,form='formatted')
      call Fig_header(unitnumber);

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,'(A)')
     $        'We display the images on an x-space cartesian grid'
     $        ,' along with their best matching (transformed) templates'
     $        ,' also on an x-space cartesian grid.'
      end if
c$$$      printing templates and images
      do nm=0,n_M-1
         F_max = Flim_x_c
         F_min = 0.0d0
         x_offset = 0.0d0 + 10.0d0*(nm/3)
         y_offset = 0.0d0 + 22.0d0*(mod(nm,3))
         x_side = 10.0d0
         y_side = 10.0d0
         call Fig_c16_carte(unitnumber,n_r,grid_x_c_,n_r,grid_x_c_
     $        ,.true.,M_x_c_(nm*n_r*n_r),F_min,F_max,x_offset,y_offset
     $        ,x_side,y_side)
         write(fig_title,'(A,I0,A)') 'M_x_c(',nm,')'
         text_length = 8
         color_index = 7
         font_type = 0
         font_size = 72
         x_loc = +0.0d0
         y_loc = +0.0d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
         ns = I_permute_(n_S-1 + nm*n_S)
         delta_x = delta_x_upd_(nm)
         delta_y = delta_y_upd_(nm)
         gamma_z = gamma_z_upd_(nm)
         call cp1_c16(n_r*n_r,S_x_c_(ns*n_r*n_r),T_x_c_)
         call transl_c_to_c(n_r,max_x_c,n_r,max_x_c,T_x_c_,+delta_x,
     $        +delta_y,T_x_c_)
         call rotate_c_to_c(n_r,max_x_c,n_r,max_x_c,T_x_c_,+gamma_z
     $        ,T_x_c_)
         F_max = Flim_x_c
         F_min = 0.0d0
         x_offset = 0.0d0 + 10.0d0*(nm/3)
         y_offset = 0.0d0 + 22.0d0*(mod(nm,3)) + 11.0d0
         x_side = 10.0d0
         y_side = 10.0d0
         call Fig_c16_carte(unitnumber,n_r,grid_x_c_,n_r,grid_x_c_
     $        ,.true.,T_x_c_,F_min,F_max,x_offset,y_offset,x_side
     $        ,y_side)
         write(fig_title,'(A,I0,A)') 'T(S_x_c(',ns,'))'
         text_length = 11
         color_index = 7
         font_type = 0
         font_size = 72
         x_loc = +0.0d0
         y_loc = +0.0d0
         call Fig_text(unitnumber,text_length,fig_title,color_index
     $        ,font_type,font_size,x_loc,y_loc,x_offset,y_offset,x_side
     $        ,y_side)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 
         write(6,'(A)') 'We write the figure file to disk.'
         write(6,'(A,A)') 'try: '
     $        ,'"gthumb test_innerproduct_svdread_batch_pos.jpg"'
      end if
c$$$      Printing figure file for image
      close(unitnumber,status='keep')
      write(fig2dev_system_call,'(A,A,1X,A)') 'fig2dev -Ljpeg -q 10 '
     $     ,trim(fname_fig),trim(fname_jpg)
      system_error = system(fig2dev_system_call)

 10   continue ! if plot_flag.eqv..false.
      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_svdread_batch_dr]'
      end if

      stop
      end
      
      real *8 function max_r8(n_x,x_)
      integer n_x,nx
      real *8 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=0.0d0
      else
         m=x_(0)
         do nx=0,n_x-1
            if (x_(nx).gt.m) then
               m = dabs(x_(nx))
            end if
         enddo
      end if
      max_r8 = m
      return
      end

      include 'get_template_size.f'
      include 'get_F2_x_c.f'
      include 'get_F2_x_c_.f'
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
      include 'test_innerproduct_svdread_batch.f'
      include 'test_innerproduct_batch_excerpt_0.f'
      include 'test_innerproduct_batch_excerpt_1.f'
      include 'test_innerproduct_batch_excerpt_2.f'
      include 'test_innerproduct_batch_excerpt_3.f'
      include 'test_innerproduct_batch_sort.f'
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
      include 'quicksort_c16.f'
      include 'quicksort_i4.f'
