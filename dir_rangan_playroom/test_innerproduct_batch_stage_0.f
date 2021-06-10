c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_batch_stage_0(verbose
     $    ,svd_calculation_type,n_r,grid_p_,n_w_,max_x_c,n_S,ld_S
     $    ,S_p_,n_M,ld_M,M_p_,delta_x_est_,delta_y_est_,gamma_z_est_
     $    ,eps_target,N_pixels_in,displacement_max,n_delta_x,n_delta_y
     $    ,n_gamma_z,delta_x_max_,delta_y_max_,gamma_z_max_,C_Z_max_)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      Summary:
c$$$      This subroutine reads in an array of templates and an array
c$$$      of images, along with parameters that specify an array of
c$$$      displacements and rotations. Then, for each pair of template
c$$$      and image, an array of innerproducts is calculated (spanning
c$$$      all the specified displacements and rotations). For each
c$$$      template-image pair this array of innerproducts is scanned
c$$$      to reveal the alignment (i.e., displacement+rotation) that
c$$$      corresponds to the highest innerproduct. This best alignment
c$$$      is stored for each template-image pair. This stored output
c$$$      is not sorted (although it can be sorted later by calling
c$$$      test_innerproduct_batch_sort.f).
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
c$$$      In this subroutine we often suppress the _k_, referring to S_k_p_ 
c$$$      as S_p_, and to grid_k_p_ as grid_p_, etc.
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$      INPUTS:
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
c$$$                            allows test_innerproduct_batch.f
c$$$                            to automatically choose one of the three
c$$$                            calculation_types based on the other
c$$$                            input parameters (e.g., by comparing
c$$$                            n_svd_l with n_delta_x*n_delta_y).
c$$$
c$$$      integer *4 n_r        Number of rings in quasi-uniform polar-grid.
c$$$                            Equivalent to ncur 
c$$$                            (passed to get_template_size).
c$$$                            Note: We assume (strongly) that n_r.gt.1
c$$$
c$$$      real *8 grid_p_       Array of size n_r holding the k-values
c$$$                            associated with the various rings in the
c$$$                            polar-grid.
c$$$                            (assumed to be a linearly spaced array from
c$$$                            0 to max_k_c/2, where max_k_c = n_r/max_x_c).
c$$$
c$$$      integer *4 n_w_       Array of size n_r holding the number of
c$$$                            points on each ring in the polar-grid.
c$$$                            Equivalent to ngridc
c$$$                            (obtained from get_template_size).
c$$$                          - Note that the sum of n_w_ is referred to
c$$$                            as n_A, and is equivalent to
c$$$                            ntemplatesize.
c$$$
c$$$      real *8 max_x_c       Side-length of templates and images in 
c$$$                            real-space (typically 1.0d0).
c$$$
c$$$      integer *4 n_S        Number of templates stored within S_p_
c$$$
c$$$      integer *4 ld_S       Leading dimension of the array storing the
c$$$                            templates. This is often the size of a
c$$$                            single template (assuming the templates are
c$$$                            stored sequentially in memory), but can be
c$$$                            larger if desired.
c$$$
c$$$      complex *16 S_p_      Array of size n_S*ld_S storing n_S templates.
c$$$                            Each template should be represented in
c$$$                            fourier space on a quasi-uniform polar-grid
c$$$                            (i.e., of form S_k_p_)
c$$$                          - Each template is of size
c$$$                            n_A = n_w_(0) + .. + n_w_(n_r-1).
c$$$                          - This subroutine assumes that templates are
c$$$                            stored sequentially, using a stride of ld_S.
c$$$                          - More specifically, we assume that S_p_ is 
c$$$                            indexed from 0 to n_S*ld_S-1, and that
c$$$                            the jth point in template ns is at 
c$$$                            S_p_(j + ns*ld_S).
c$$$
c$$$      integer *4 n_M        Number of images stored within M_p_
c$$$
c$$$      integer *4 ld_M       Leading dimension of the array storing the
c$$$                            images. This is often the size of a
c$$$                            single image (assuming the images are
c$$$                            stored sequentially in memory), but can be
c$$$                            larger if desired.
c$$$
c$$$      complex *16 M_p_      Array of size n_M*ld_M storing n_M images.
c$$$                            Each image should be represented in
c$$$                            fourier space on a quasi-uniform polar-grid
c$$$                            (i.e., of form M_k_p_)
c$$$                          - Each image is of size
c$$$                            n_A = n_w_(0) + .. + n_w_(n_r-1).
c$$$                          - This subroutine assumes that images are
c$$$                            stored sequentially, using a stride of ld_M.
c$$$                          - More specifically, we assume that M_p_ is 
c$$$                            indexed from 0 to n_M*ld_M-1, and that
c$$$                            the jth point in image nm is at 
c$$$                            M_p_(j + nm*ld_M).
c$$$
c$$$      real *8 delta_x_est_  Array of size n_M storing the current
c$$$                            estimates of the x-displacements associated
c$$$                            with each of the n_M images.
c$$$
c$$$      real *8 delta_y_est_  Array of size n_M storing the current
c$$$                            estimates of the y-displacements associated
c$$$                            with each of the n_M images.
c$$$
c$$$      real *8 gamma_z_est_  Array of size n_M storing the current
c$$$                            estimates of the in-plane rotations
c$$$                            associated with each of the n_M images.
c$$$
c$$$                          - Note: The estimated alignments are accounted
c$$$                            for by *first* rotating each image by
c$$$                            gamma_z = -gamma_z_est_(nm), 
c$$$                            and *then* translating the result by 
c$$$                            delta_x = -delta_x_est_(nm),
c$$$                            delta_y = -delta_y_est_(nm).
c$$$                          - Note: These estimated alignments are
c$$$                            accounted for on the k-space polar-grid
c$$$                            (i.e., in terms of M_p_) and not in terms
c$$$                            of the bessel-expansion (i.e., not in terms
c$$$                            of M_q_). Consequently, the estimated 
c$$$                            alignments can exceed the 3.0d0 pixel bound
c$$$                            expected within the array of displacements.
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
c$$$      real *8 displacement_max Maximum displacement (in the same units 
c$$$                            used for delta_x and delta_y). For each
c$$$                            image-template pair we will search through 
c$$$                            all n_delta_x-by-n_delta_y-by-n_gamma_z
c$$$                            parameters for the largest possible overlap.
c$$$                          - This search will be limited to only those
c$$$                            displacements (i.e., delta_x , delta_y) for
c$$$                            which:
c$$$                            \vec{delta} = [delta_x_est,delta_y_est] + 
c$$$                                          [delta_x , delta_y]
c$$$                            has a magnitude less than
c$$$                            displacement_max.
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
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$
c$$$      OUTPUTS:
c$$$
c$$$      real *8 delta_x_max_  Array of size n_S*n_M.
c$$$                          - Assuming this array is indexed from 0
c$$$                            to n_S*n_M-1, delta_x_max_(ns + nm*n_S) 
c$$$                            holds the delta_x value associated with the 
c$$$                            best alignment between template ns and image nm.
c$$$                          - Note that this "best alignment" is
c$$$                            restricted to have a displacement vector
c$$$                            with a magnitude less than displacement_max.
c$$$
c$$$      real *8 delta_y_max_  Array of size n_S*n_M.
c$$$                          - Assuming this array is indexed from 0
c$$$                            to n_S*n_M-1, delta_y_max_(ns + nm*n_S) 
c$$$                            holds the delta_y value associated with the 
c$$$                            best alignment between template ns and image nm.
c$$$                          - Note that this "best alignment" is
c$$$                            restricted to have a displacement vector
c$$$                            with a magnitude less than displacement_max.
c$$$
c$$$      real *8 gamma_z_max_  Array of size n_S*n_M.
c$$$                          - Assuming this array is indexed from 0
c$$$                            to n_S*n_M-1, gamma_z_max_(ns + nm*n_S) 
c$$$                            holds the gamma_z value associated with the 
c$$$                            best alignment between template ns and image nm.
c$$$                          - Note that this "best alignment" is
c$$$                            restricted to have a displacement vector
c$$$                            with a magnitude less than displacement_max.
c$$$
c$$$      complex *16 C_Z_max_  Array of size n_S*n_M.
c$$$                          - Assuming this array is indexed from 0
c$$$                            to n_S*n_M-1, C_Z_max_(ns + nm*n_S) 
c$$$                            holds the (complex) normalized innerproduct
c$$$                            associated with the best alignment between 
c$$$                            template ns and image nm.
c$$$                          - Note that this "best alignment" is
c$$$                            restricted to have a displacement vector
c$$$                            with a magnitude less than displacement_max.
c$$$                          - Note that this is complex. While the 
c$$$                            imaginary part should be 0.0d0, the 
c$$$                            approximations used will incur a small 
c$$$                            error; consequently, the imaginary part will
c$$$                            in general be nonzero.
c$$$
c$$$                          - Note: For any image-template pair, the best
c$$$                            alignment determined by:
c$$$                            delta_x = delta_x_max_(ns + nm*n_S), 
c$$$                            delta_y = delta_y_max_(ns + nm*n_S), 
c$$$                            gamma_z = gamma_z_max_(ns + nm*n_S) 
c$$$                            can be used to update the estimated 
c$$$                            alignment via:
c$$$                            delta_x_upd_(nm) = delta_x_est_(nm) 
c$$$                                               + delta_x
c$$$                            delta_y_upd_(nm) = delta_y_est_(nm) 
c$$$                                               + delta_y
c$$$                            gamma_z_upd_(nm) = gamma_z_est_(nm) 
c$$$                                               + gamma_z
c$$$                          - Note that the sum delta_x_upd and 
c$$$                            delta_y_upd should (together) yield a
c$$$                            displacement vector with a magnitude less
c$$$                            than displacement_max
c$$$
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,ld_S,n_M,ld_M
      real *8 grid_p_(0:n_r-1),max_x_c
      complex *16 S_p_(0:ld_S*n_S-1)
      complex *16 M_p_(0:ld_M*n_M-1)
      real *8 delta_x_est_(0:n_M-1),delta_y_est_(0:n_M-1)
      real *8 gamma_z_est_(0:n_M-1),eps_target,N_pixels_in
      real *8 displacement_max
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8 delta_x_max_(0:n_S*n_M-1)
      real *8 delta_y_max_(0:n_S*n_M-1)
      real *8 gamma_z_max_(0:n_S*n_M-1)
      complex *16 C_Z_max_(0:n_S*n_M-1)
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      logical warning_flag
      data warning_flag / .true. /
      real *8 R_max,K_max,delta,delta_max,N_pixels
      complex *16, allocatable :: Z_svds_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer l
c$$$      fftw
      integer *8, allocatable :: fftw_plan_frwd__(:)
      integer *8, allocatable :: fftw_plan_back__(:)
      complex *16, allocatable :: fftw_in1__(:)
      complex *16, allocatable :: fftw_out__(:)
c$$$      indices
      integer nr,n_w_max,n_A,na
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_max,ndy_max,ngz_max
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      character(len=64) format_string
c$$$      array of innerproducts to hold measurements
      real *8 pi
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      integer *4 n_M_sub,nM_sub,r_tmp,a_tmp,n_tmp,i_tmp
      integer *4, allocatable :: n_M_sub_(:)
      integer *4, allocatable :: I_M_sub_(:)

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_stage_0]'
      end if      
      if (verbose.gt.0) then
         write(6,'(A,I0,A,F5.3,A,F3.1,A,I0,A,I0,A,I0)') 'n_r: ' ,n_r
     $        ,'; eps_target: ',eps_target,'; N_pixels_in: ',N_pixels_in
     $        ,'; n_delta_x: ',n_delta_x ,'; n_delta_y: ',n_delta_y
     $        ,'; n_gamma_z: ',n_gamma_z
      end if

      do na=0,n_S*n_M-1
         delta_x_max_(na) = 0.0d0
         delta_y_max_(na) = 0.0d0
         gamma_z_max_(na) = 0.0d0
         C_Z_max_(na) = (0.0d0,0.0d0)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating indices'
      end if
      pi = 4.0*atan(1.0)
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_A: ',n_A,'; n_w_max: ',n_w_max
         write(format_string,'(A,I0,A)') '(A,',n_r,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,n_r-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      do ndx=0,n_delta_x-1
         if (n_delta_x.gt.1) then
            delta_x = (-N_pixels_in + ndx*2*N_pixels_in/(n_delta_x-1))
     $           /n_r*max_x_c
         else
            delta_x = 0.0d0
         end if
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            if (n_delta_y.gt.1) then
               delta_y = (-N_pixels_in + ndy*2*N_pixels_in/(n_delta_y
     $              -1))/n_r*max_x_c
            else
               delta_y = 0.0d0
            end if
            delta_y_(ndy) = delta_y
         enddo
      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)')
     $        ' ncur or so.'
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      do ngz=0,n_gamma_z-1
         if (n_gamma_z.gt.1) then
            gamma_z = (2*pi*ngz)/n_gamma_z
         else
            gamma_z = 0.0d0
         end if
         gamma_z_(ngz) = gamma_z
      enddo
      if (verbose.gt.1) then
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
         write(6,'(A)') ' Selecting svd library to use.'
      end if
      R_max = 2*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
            delta = dsqrt(delta_x_(ndx)**2 + delta_y_(ndy)**2)
            if (delta.gt.delta_max) then
               delta_max = delta
            end if
         enddo
      enddo
      N_pixels = delta_max/dsqrt(2.0d0)*2*K_max
      if (verbose.gt.0) then
         write(6,'(A,F8.3,A,F8.3)') 'R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') 'K_max: ',K_max,'; N_pixels: '
     $        ,N_pixels
      end if
      if (K_max.gt.96 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, K_max ',K_max
     $        ,' too large in test_innerproduct_batch_stage_0'
      end if
      if (N_pixels.gt.3 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, N_pixels ',N_pixels
     $        ,' too large in test_innerproduct_batch_stage_0'
      end if
      if (eps_target.lt.0.1d0 .and. warning_flag) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in test_innerproduct_batch_stage_0'
      end if
      if (.false.) then
      include './dir_gen_Jsvd_1/gen_Jsvd_svdpick.txt'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdload.txt'
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
         write(6,'(A,I0)') 'svd_unitnumber: ',svd_unitnumber
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up displacement-operator Z_svdd_'
         write(6,'(A)') ' associated with svd-expansion.'
      end if
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      timing_tic = second()
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svdd: total time ',timing_toc
     $        -timing_tic
      end if

      allocate(n_M_sub_(0:n_M-1))
      allocate(I_M_sub_(0:n_M-1))
      n_M_sub = min(8,n_M)
      do nM_sub=0,n_M_sub-1
         if (nM_sub.eq.0) then
            n_M_sub_(0) = n_M/n_M_sub
            I_M_sub_(0) = 0
         else if (nM_sub.gt.0 .and. nM_sub.lt.n_M_sub-1) then
            n_M_sub_(nM_sub) = n_M /n_M_sub
            I_M_sub_(nM_sub) = I_M_sub_(nM_sub-1) + n_M_sub_(nM_sub-1)
         else if (nM_sub.eq.n_M_sub-1) then
            I_M_sub_(n_M_sub-1) = I_M_sub_(n_M_sub-2) + n_M_sub_(n_M_sub
     $           -2)
            n_M_sub_(n_M_sub-1) = n_M - I_M_sub_(n_M_sub-1)
         end if
      enddo
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_M_sub
     $        ,'(I0,1X))'
         write(6,format_string) 'n_M_sub_: '
     $        ,(n_M_sub_(nM_sub),nM_sub=0
     $        ,n_M_sub-1)
         write(6,format_string) 'I_M_sub_: '
     $        ,(I_M_sub_(nM_sub),nM_sub=0
     $        ,n_M_sub-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans'
      end if
      allocate(fftw_plan_frwd__(0:n_r*n_M_sub-1))
      allocate(fftw_plan_back__(0:n_r*n_M_sub-1))
      allocate(fftw_in1__(0:n_A*n_M_sub-1))
      allocate(fftw_out__(0:n_A*n_M_sub-1))
      do nM_sub=0,n_M_sub-1
         r_tmp = nM_sub*n_r
         a_tmp = nM_sub*n_A
         na = 0
         do nr=0,n_r-1
            call dfftw_plan_dft_1d_(fftw_plan_frwd__(r_tmp+nr),n_w_(nr)
     $          ,fftw_in1__(a_tmp+na),fftw_out__(a_tmp +na)
     $          ,FFTW_FORWARD,FFTW_ESTIMATE) 
            call dfftw_plan_dft_1d_(fftw_plan_back__(r_tmp+nr),n_w_(nr)
     $          ,fftw_out__(a_tmp+na),fftw_in1__(a_tmp +na)
     $          ,FFTW_BACKWARD,FFTW_ESTIMATE) 
            na = na + n_w_(nr)
         enddo
      enddo

      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_tmp,i_tmp,r_tmp,a_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub-1
         n_tmp = n_M_sub_(nM_sub)
         i_tmp = I_M_sub_(nM_sub)
         r_tmp = nM_sub*n_r
         a_tmp = nM_sub*n_A
         if (verbose.gt.0) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' test_innerproduct_batch_stage_1'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call test_innerproduct_batch_stage_1(verbose
     $        ,svd_calculation_type,n_r,grid_p_,n_w_,max_x_c,n_S,ld_S
     $        ,S_p_,n_tmp,ld_M,M_p_(ld_M*i_tmp),delta_x_est_(i_tmp)
     $        ,delta_y_est_(i_tmp),gamma_z_est_(i_tmp),eps_target
     $        ,N_pixels_in,displacement_max,n_delta_x,n_delta_y
     $        ,n_gamma_z,delta_x_max_(n_S*i_tmp),delta_y_max_(n_S
     $        *i_tmp),gamma_z_max_(n_S*i_tmp),C_Z_max_(n_S*i_tmp)
     $        ,fftw_plan_frwd__(r_tmp),fftw_plan_back__(r_tmp)
     $        ,fftw_in1__(a_tmp),fftw_out__(a_tmp),n_svd_r,n_svd_d
     $        ,n_svd_l,svd_r_,svd_d_,svd_l_,svd_U_d_,svd_s_,svd_V_r_
     $        ,Z_svdd_)
      enddo !do nM_sub=0,n_M_sub-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      if (verbose.gt.1) then
         write(6,'(A,A,A,F8.5)') ' Finished',
     $        ' test_innerproduct_batch_stage_1:', 
     $        ' total time ',
     $        timing_toc-timing_tic
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
      deallocate(Z_svdd_)
      deallocate(gamma_z_)
      deallocate(delta_y_)
      deallocate(delta_x_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans.'
      end if
      do nr=0,n_M_sub*n_r-1
         call dfftw_destroy_plan_(fftw_plan_back__(nr))
         call dfftw_destroy_plan_(fftw_plan_frwd__(nr))
      enddo
      deallocate(fftw_out__)
      deallocate(fftw_in1__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_plan_frwd__)

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating svd arrays.'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_stage_0]'
      end if

      end


