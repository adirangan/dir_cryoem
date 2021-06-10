c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_batch_stage_1(verbose
     $     ,svd_calculation_type_in,n_r,grid_p_,n_w_,max_x_c,n_S,ld_S
     $     ,S_p_,n_M,ld_M,M_p_,delta_x_est_,delta_y_est_,gamma_z_est_
     $     ,eps_target,N_pixels_in,displacement_max,n_delta_x,n_delta_y
     $     ,n_gamma_z,delta_x_max_,delta_y_max_,gamma_z_max_,C_Z_max_
     $     ,fftw_plan_frwd_,fftw_plan_back_,fftw_in1_,fftw_out_,n_svd_r
     $     ,n_svd_d,n_svd_l,svd_r_,svd_d_,svd_l_,svd_U_d_,svd_s_
     $     ,svd_V_r_,Z_svdd_)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer svd_calculation_type_in,svd_calculation_type
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
      integer *4 n_svd_r,n_svd_d,n_svd_l
      real *8 svd_r_(0:n_svd_r-1)
      real *8 svd_d_(0:n_svd_d-1)
      integer *4 svd_l_(0:n_svd_l-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      real *8 svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      complex *16 Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1)
      integer *4 nsvd_r,nsvd_d,nsvd_l
      logical warning_flag
      data warning_flag / .true. /
      real *8 R_max,K_max,delta,delta_max,N_pixels
      complex *16, allocatable :: Z_svds_(:)
c$$$      Z_svdd_ should be passed as input
c$$$      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer l
c$$$      fftw plans and workspace should be passed in as input
c$$$      integer *8, allocatable :: fftw_plan_frwd_(:)
c$$$      integer *8, allocatable :: fftw_plan_back_(:)
c$$$      complex *16, allocatable :: fftw_in1_(:)
c$$$      complex *16, allocatable :: fftw_out_(:)
      integer *8 fftw_plan_frwd_(0:0)
      integer *8 fftw_plan_back_(0:0)
      complex *16 fftw_in1_(0:0)
      complex *16 fftw_out_(0:0)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*)
      complex *16 fftw_out_last_(*)
c$$$      indices
      integer ns,nm,nr,n_w_max,n_A,na
      complex *16, allocatable :: S_q_(:)
      complex *16, allocatable :: M_q_(:)
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_max,ndy_max,ngz_max
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      character(len=64) format_string
c$$$      array of innerproducts to hold measurements
      complex *16 C_S_(0:n_S-1),C_M,C_Z,C_Z_max,Z_q
      complex *16, allocatable :: Z_q_(:)
      real *8 pi
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      real *8 timing_tic_1,timing_toc_1,timing_tot_1
      real *8 timing_tic_2,timing_toc_2,timing_tot_2
      real *8 timing_tic_3,timing_toc_3,timing_tot_3
      real *8 timing_tic_4,timing_toc_4,timing_tot_4

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_stage_1]'
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

c$$$      fftw plans and workspace should be passed in as input
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') 'Generating fftw_plans'
c$$$      end if
c$$$      allocate(fftw_plan_frwd_(0:n_r-1))
c$$$      allocate(fftw_plan_back_(0:n_r-1))
c$$$      allocate(fftw_in1_(0:n_A-1))
c$$$      allocate(fftw_out_(0:n_A-1))
c$$$      na = 0
c$$$      do nr=0,n_r-1
c$$$         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
c$$$     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_ESTIMATE) 
c$$$         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
c$$$     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_ESTIMATE) 
c$$$         na = na + n_w_(nr)
c$$$      enddo
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

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
     $        ,' too large in test_innerproduct_batch_stage_1'
      end if
      if (N_pixels.gt.3 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, N_pixels ',N_pixels
     $        ,' too large in test_innerproduct_batch_stage_1'
      end if
      if (eps_target.lt.0.1d0 .and. warning_flag) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in test_innerproduct_batch_stage_1'
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      The svd parameters should be passed in as input
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      if (.false.) then
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdpick.txt'
c$$$      end if
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdload.txt'

      if (verbose.gt.1) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd_calculation_type: '
      end if
      if (svd_calculation_type_in.lt.0) then
         if (.false.) then
         else if (n_svd_l.gt.n_delta_x*n_delta_y) then
            svd_calculation_type = 0
         else if (n_svd_l.le.n_delta_x*n_delta_y) then
            svd_calculation_type = 1
         else
            svd_calculation_type = 2
         end if !switch
      else if (svd_calculation_type_in.ge.0) then
         svd_calculation_type = max(0,min(2,svd_calculation_type_in))
         if (verbose.gt.0) then
            write(6,'(A,I0)') 'forcing: svd_calculation_type = '
     $           ,svd_calculation_type
         end if !verbose
      end if !choose svd_calculation_type
      if (verbose.gt.0) then
         if (.false.) then
         else if (svd_calculation_type.eq.0) then
            write(6,'(A)') 'svd_calculation_type 0: multiply then fft.'
         else if (svd_calculation_type.eq.1) then
            write(6,'(A)') 'svd_calculation_type 1: fft then multiply.'
         else if (svd_calculation_type.eq.2) then
            write(6,'(A)')
     $           'svd_calculation_type 2: brute-force displacements.'
         end if
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      The operator Z_svdd_ should be passed in as input
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ' Setting up displacement-operator Z_svdd_'
c$$$         write(6,'(A)') ' associated with svd-expansion.'
c$$$      end if
c$$$      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
c$$$      timing_tic = second()
c$$$      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
c$$$     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
c$$$      timing_toc = second()
c$$$      if (verbose.gt.0) then
c$$$         write(6,'(A,F8.5)')
c$$$     $        ' finished innerproduct_q__k_svdd: total time ',timing_toc
c$$$     $        -timing_tic
c$$$      end if

      if (verbose.gt.1) then
         write(6,'(A,A)') ' Allocating template-dependent operator'
     $        ,' Z_svd_ used in svd-expansion.'
      end if
      allocate(Z_svds_(0:n_svd_l*n_w_max-1))
      if (verbose.gt.1) then
         write(6,'(A,A)') ' Allocating temporary array Z_tmpC_ '
     $        ,'to hold innerproducts associated with '
         write(6,'(A,A)') ' the displacements across all rotations '
     $        ,'for any single image-template pair'
      end if
      allocate(Z_tmpC_(0:n_w_max*n_delta_x*n_delta_y-1))
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_q_ to hold bessel-function'
         write(6,'(A)') ' expansion of each template within S_p_.'
      end if
      allocate(S_q_(0:n_A*n_S))
      if (verbose.gt.1) then
         write(6,'(A)') ' Converting S_p_ to S_q_.'
      end if
      timing_tic = second()
      do ns=0,n_S-1
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,S_p_(ns*ld_S),S_q_(ns*n_A))
         call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_q_(ns*n_A),S_q_(ns
     $        *n_A),C_S_(ns))
         C_S_(ns) = zsqrt(C_S_(ns))/(n_r*n_r)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,2F16.3)') 'ns: ',ns,'; C_S: ',C_S_(ns)
         end if
      enddo
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished converting S_p_ to S_q_: total time '
     $        ,timing_toc-timing_tic
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array M_q_ to hold a single'
         write(6,'(A)') ' bessel-function-expansion of an image.'
      end if
      allocate(M_q_(0:n_A-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' Allocating array Z_q_ to hold the innerproucts'
         write(6,'(A)') ' associated with displacements and rotations '
         write(6,'(A)') ' for any single image-template pair. '
      end if
      allocate(Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1))

      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating innerproducts...'
      end if
      timing_tot_1 = 0.0d0
      timing_tot_2 = 0.0d0
      timing_tot_3 = 0.0d0
      timing_tot_4 = 0.0d0
      do nm=0,n_M-1
         timing_tic_1 = second()
         delta_x = -delta_x_est_(nm)
         delta_y = -delta_y_est_(nm)
         gamma_z = -gamma_z_est_(nm)
         call rotate_p_to_p(n_r,n_w_,n_A,M_p_(nm*ld_M),+gamma_z,M_q_)
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_q_,+delta_x,+delta_y
     $        ,M_q_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,M_q_,M_q_)
         call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_q_,M_q_,C_M)
         C_M = zsqrt(C_M)/(n_r*n_r)      
         if (verbose.gt.2) then
            write(6,'(A,I0,A,2F16.3)') ' nm: ',nm,'; C_M: ',C_M
         end if
         timing_toc_1 = second()
         timing_tot_1 = timing_tot_1 + (timing_toc_1 - timing_tic_1)
         do ns=0,n_S-1
            C_Z = C_S_(ns)*C_M
            if (.false.) then
            else if (svd_calculation_type.eq.0) then
               timing_tic_2 = second();
               call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_
     $              ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_(ns*n_A)
     $              ,M_q_ ,Z_svds_)
               timing_toc_2 = second()
               timing_tot_2 = timing_tot_2 + (timing_toc_2 -
     $              timing_tic_2)
               timing_tic_3 = second()
               call test_innerproduct_batch_excerpt_0(n_delta_x
     $              ,n_delta_y,n_gamma_z,gamma_z_,n_w_max,n_svd_l
     $              ,Z_svdd_,Z_svds_,Z_tmpC_,fftw_plan_back_last
     $              ,fftw_in1_last_,fftw_out_last_,n_r,C_Z,Z_q_)
               timing_toc_3 = second()
               timing_tot_3 = timing_tot_3 + (timing_toc_3 -
     $              timing_tic_3)
            else if (svd_calculation_type.eq.1) then
               timing_tic_2 = second();
               call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_
     $              ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_(ns*n_A)
     $              ,M_q_ ,Z_svds_)
               do l=0,n_svd_l-1
                  call cps_c16(n_w_max,Z_svds_(l),n_svd_l,fftw_out_last_
     $                 ,1)
                  call dfftw_execute_(fftw_plan_back_last)
                  call cps_c16(n_w_max,fftw_in1_last_,1,Z_svds_(l)
     $                 ,n_svd_l)
               enddo !do l=0,n_svd_l-1
               timing_toc_2 = second()
               timing_tot_2 = timing_tot_2 + (timing_toc_2 -
     $              timing_tic_2)
               timing_tic_3 = second()
               call test_innerproduct_batch_excerpt_1(n_delta_x
     $              ,n_delta_y,n_gamma_z,gamma_z_,n_w_max,n_svd_l
     $              ,Z_svdd_,Z_svds_,Z_tmpC_,n_r,C_Z,Z_q_)
               timing_toc_3 = second()
               timing_tot_3 = timing_tot_3 + (timing_toc_3 -
     $              timing_tic_3)
            else if (svd_calculation_type.eq.2) then
               timing_tic_3 = second()
               call test_innerproduct_batch_excerpt_2(n_delta_x,delta_x_
     $              ,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $              ,fftw_plan_frwd_,fftw_in1_,fftw_out_
     $              ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_
     $              ,n_r,grid_p_,n_w_,n_A,C_Z,S_p_(ns*ld_S),M_q_,Z_q_)
               timing_toc_3 = second()
               timing_tot_3 = timing_tot_3 + (timing_toc_3 -
     $              timing_tic_3)
            end if !svd_calculation_type
            timing_tic_4 = second()
            call test_innerproduct_batch_excerpt_4(delta_x_est_(nm)
     $           ,delta_y_est_(nm),delta_x_,delta_y_,displacement_max
     $           ,n_delta_x,n_delta_y,n_gamma_z,Z_q_,ndx_max,ndy_max
     $           ,ngz_max,C_Z_max)
            delta_x_max_(ns + nm*n_S) = delta_x_(ndx_max)
            delta_y_max_(ns + nm*n_S) = delta_y_(ndy_max)
            gamma_z_max_(ns + nm*n_S) = gamma_z_(ngz_max)
            C_Z_max_(ns + nm*n_S) = C_Z_max
            timing_toc_4 = second()
            timing_tot_4 = timing_tot_4 + (timing_toc_4 - timing_tic_4)
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I2,1X,I2,1X,I2,A,2F8.3)')
     $              ' Match to S_q_(',ns
     $              ,') at ndx_max,ndy_max,ngz_max: ' ,ndx_max,ndy_max
     $              ,ngz_max,'; C_Z_max: ',C_Z_max
            end if !verbose
         enddo !do ns=0,n_S-1
      enddo !do nm=0,n_M-1
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'Total (n_S templates)*(n_M images): ',n_S
     $        *n_M
         write(6,'(A,F8.5,A,F8.5)')
     $        'transforming M_p_ to M_q_: total time ' ,timing_tot_1
     $        ,'; time_per ',timing_tot_1/max(1,n_S*n_M)
         write(6,'(A,F8.5,A,F8.5)')
     $        '   innerproduct_q__k_svds: total time ' ,timing_tot_2
     $        ,'; time_per ',timing_tot_2/max(1,n_S*n_M)
         write(6,'(A,I0,A,F8.5,A,F8.5)') '                excerpt_'
     $        ,svd_calculation_type,': total time ' ,timing_tot_3
     $        ,'; time_per ',timing_tot_3/max(1,n_S*n_M)
         write(6,'(A,F8.5,A,F8.5)')
     $        '                excerpt_4: total time ' ,timing_tot_4
     $        ,'; time_per ',timing_tot_4/max(1,n_S*n_M)
      end if !verbose

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
      deallocate(Z_q_)
      deallocate(M_q_)
      deallocate(S_q_)
      deallocate(Z_tmpC_)
      deallocate(Z_svds_)
c$$$      Z_svdd_ should be passed in as input
c$$$      deallocate(Z_svdd_)
      deallocate(gamma_z_)
      deallocate(delta_y_)
      deallocate(delta_x_)
c$$$      fftw plans and workspace should be passed in as input
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ' Destroying fftw_plans.'
c$$$      end if
c$$$      do nr=0,n_r-1
c$$$         call dfftw_destroy_plan_(fftw_plan_back_(nr))
c$$$         call dfftw_destroy_plan_(fftw_plan_frwd_(nr))
c$$$      enddo
c$$$      deallocate(fftw_out_)
c$$$      deallocate(fftw_in1_)
c$$$      deallocate(fftw_plan_back_)
c$$$      deallocate(fftw_plan_frwd_)

c$$$      svd arrays should be passed in as input
c$$$      if (verbose.gt.1) then
c$$$         write(6,'(A)') ' Deallocating svd arrays.'
c$$$      end if
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_stage_1]'
      end if

      end


