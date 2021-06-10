      subroutine ti8_dr_digits_single(
     $     verbose
     $     ,half_diameter_k_c
     $     ,fix_phi_M
     $     ,fix_delta_M_x_p_r
     $     ,fix_delta_M_x_p_w
     $     ,fix_phi_S
     $     ,fix_delta_S_x_p_r
     $     ,fix_delta_S_x_p_w
     $     ,flag_RTRT_vs_RTTR
     $     ,gamma_z_est
     $     ,delta_est_x_c_0
     $     ,delta_est_x_c_1
     $     ,gamma_z_upd_0in
     $     ,delta_upd_x_c_0
     $     ,delta_upd_x_c_1
     $     )
c$$$      checks ti8_dr_digits for a single image-template pair. ;
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose !verbosity level. ;
      real *8 diameter_x_c_0 ! real *8: diameter of image in x_c_0 coordinates. ;
      real *8 diameter_x_c_1 ! real *8: diameter of image in x_c_1 coordinates. ;
      real *8 half_diameter_x_c_0 ! real *8: half_diameter of image in x_c_0 coordinates. ;
      real *8 half_diameter_x_c_1 ! real *8: half_diameter of image in x_c_1 coordinates. ;
      real *8 dx_c_0 !real *8: stepsize in x_c_0. ;
      real *8 dx_c_1 !real *8: stepsize in x_c_1. ;
      real *8 x_c_0_max !real *8: max in x_c_0. ;
      real *8 x_c_1_max !real *8: max in x_c_1. ;
      real *8 diameter_k_c ! real *8: diameter of image in k_c coordinates. ;
      real *8 half_diameter_k_c ! real *8: half_diameter of image in k_c coordinates. ;
      integer quad_n_r ! integer: number of quadrature nodes in the radial direction. ;
      integer quad_n_w_sum ! integer: number of points in quad_k_p. ;
      integer quad_n_w_max ! integer: number of points at quad_grid_k_p_r_(quad_n_r-1). ;
      real *8, allocatable :: jacpts_x_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: jacpts_w_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: quad_grid_k_p_r_(:) !real *8 array, size quad_n_r. ;
      real *8, allocatable :: quad_weight_k_p_r_(:) !real *8 array, size quad_n_r. ;
      integer *4, allocatable :: quad_n_polar_a_(:) !integer *4 array, size quad_n_r. ;
      integer *4, allocatable :: quad_n_w_(:) !integer *4 array, size quad_n_r. ;
      real *8, allocatable :: quad_grid_k_p_0_(:) !real *8 array, size quad_n_w_sum. ;
      real *8, allocatable :: quad_grid_k_p_1_(:) !real *8 array, size quad_n_w_sum. ;
      real *8, allocatable :: quad_weight_k_p_(:) !real *8 array, size quad_n_w_sum. ;
      integer n_delta_v ! number of displacements to analyze; only 1 here. ;
      real *8 fix_phi,fix_phi_S,fix_phi_M !real *8 parameters for planar function. ;
      real *8 fix_delta_x_p_r,fix_delta_x_p_w !real *8, parameters for fix_delta_x_c_. ;
      real *8, allocatable :: fix_delta_x_c_(:) !real *8 array, size 2. plane-wave displacement. ;
      real *8 tmp_j0,tmp_j2 !temporary: real *8, bessel evaluation. ;      
      real *8 form_I !real *8, integral of plane-wave. ;
      real *8 tmp_p,tmp_a !temporary: real *8, parameters for polynomial function. ;
      real *8 tmp_f,tmp_e !temporary: real *8, exponent and planar function. ;
      real *8 tmp_phi_M,tmp_phi_S,tmp_phi_n !temporary: real *8, phase. ;
      complex *16, allocatable :: quad_P_k_p_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_x_c_ ;
      complex *16 quad_I !complex *16, quadrature approximation to integral of plane-wave. ;
      real *8 fix_delta_S_x_p_r,fix_delta_S_x_p_w !real *8, parameters for fix_delta_S_x_c_. ;
      real *8, allocatable :: fix_delta_S_x_c_(:) !real *8 array, size 2. plane-wave displacement. ;
      real *8 fix_delta_M_x_p_r,fix_delta_M_x_p_w !real *8, parameters for fix_delta_M_x_c_. ;
      real *8, allocatable :: fix_delta_M_x_c_(:) !real *8 array, size 2. plane-wave displacement. ;
      real *8 fix_delta_C_x_p_r,fix_delta_C_x_p_w !real *8, parameters for fix_delta_C_x_c_. ;
      real *8, allocatable :: fix_delta_C_x_c_(:) !real *8 array, size 2. plane-wave displacement. ;
      real *8, allocatable :: tmp_delta_M_x_c_(:) !temporary: real *8 array, size 2. plane-wave displacement. ;
      real *8 tmp_delta_M_x_p_r,tmp_delta_M_x_p_w !temporary: real *8, parameters for tmp_delta_M_x_c_. ;
      real *8, allocatable :: tmp_delta_S_x_c_(:) !temporary: real *8 array, size 2. plane-wave displacement. ;
      real *8 tmp_delta_S_x_p_r,tmp_delta_S_x_p_w !temporary: real *8, parameters for tmp_delta_S_x_c_. ;
      real *8, allocatable :: tmp_delta_Y_x_c_(:) !temporary: real *8 array, size 2. plane-wave displacement. ;
      integer n_gamma_z_upd,n_gamma_z_nuf !integer: size of arrays for uniform and nonuniform CTF_R_S_. ;
      real *8, allocatable :: gamma_z_upd_(:) !real *8 array, size n_gamma_z_upd. uniform rotations for CTF_R_S_. ;
      real *8, allocatable :: gamma_z_nuf_(:) !real *8 array, size n_gamma_z_nuf. nonuniform rotations for CTF_R_S_. ;
      real *8 form_I_SM !real *8, <P_S,P_M>. ;
      real *8 form_I_SCM !real *8, <P_S,conj(P_C)*P_M>. ;
      complex *16, allocatable :: quad_P_S_k_p_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_S_x_c_ ;
      complex *16, allocatable :: quad_P_S_k_q_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_S_x_c_ ;
      complex *16, allocatable :: quad_P_M_k_p_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_M_x_c_ ;
      complex *16, allocatable :: quad_P_M_k_q_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_M_x_c_ ;
      complex *16, allocatable :: quad_P_C_k_p_(:) !complex *16 array, size quad_n_w_sum. plane-wave with displacement fix_delta_C_x_c_ ;
      complex *16, allocatable :: quad_P_T_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_P_T_k_q_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_P_R_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_P_R_k_q_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_P_N_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_P_N_k_q_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_Y_S_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_Y_C_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_Y_T_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_Y_R_k_p_(:) !temporary: complex *16 array, size quad_n_w_sum. ;
      integer *8, allocatable :: quad_fftw_plan_frwd_(:) !integer *8 array, size quad_n_r. ;
      integer *8, allocatable :: quad_fftw_plan_back_(:) !integer *8 array, size quad_n_r. ;
      complex *16, allocatable :: quad_fftw_0in_(:) !complex *16 array, size quad_n_w_sum. ;
      complex *16, allocatable :: quad_fftw_out_(:) !complex *16 array, size quad_n_w_sum. ;
      real *8 delta_est_x_c_0,delta_est_x_c_1,gamma_z_est !parameter estimate. ;
      real *8 gamma_z_upd_0in !parameter update. ;
      integer ngamma_z_upd_0in !parameter update. ;
      real *8 delta_upd_x_c_0,delta_upd_x_c_1,gamma_z_upd !parameter update. ;
      complex *16 quad_I_SM !complex *16, quadrature approximation to <P_S,P_M>. ;
      complex *16 quad_I_SCM !complex *16, quadrature approximation to <P_S,conj(P_C)*P_M>. ;
      complex *16, allocatable :: form_X0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: form_SX0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: form_MX0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_X0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_X1_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_X2_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_X3_(:) !complex *16 array, size quad_n_w_max. ;
      real *8 eps_svd !real *8: svd tolerance epsilon, typically 0.1d0, 0.01d0 or 0.001d0. ;
      integer *4 n_eps_svd,neps_svd !integer *4: number of svd_eps. ;
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_6/Jsvd_excerpt_var_allocation_0on_0.f'
      integer *4 n_svd_max
      parameter (n_svd_max=512)
      real *8 svd_d_max,svd_r_max ! variables used to set ranges for svd-expansion of translation. ;
      real *8 , allocatable :: svd_polyval_U_d_(:) ! array storing polynomial evaluations of svd_U_d_ at the various delta values. ;
      real *8 , allocatable :: svd_polyval_V_r_(:) ! array storing polynomial evaluations of svd_V_r_ at the various grid_k_p_(nr) values. ;
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels ! variables used to select svd library. ;
      real *8 tmp_delta_x_c_0_(0:0),tmp_delta_x_c_1_(0:0) !temporary: real *8 arrays of size 1 for get_svd_2.f ;
      complex *16, allocatable :: quad_SX2_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_MX2_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_SX1_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_MX1_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_SX0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: quad_MX0_(:) !complex *16 array, size quad_n_w_max. ;
      complex *16, allocatable :: Z_S_svdd_(:) !array storing displacement-operator for the delta-side of the svd-expansion applied to S. ;
      complex *16, allocatable :: Z_M_svdd_(:) !array storing displacement-operator for the delta-side of the svd-expansion applied to M. ;
      complex *16, allocatable :: Z_S_svdr_(:) !array storing displacement-operator for the k-side of the svd-expansion applied to S. ;
      complex *16, allocatable :: Z_M_svdr_(:) !array storing displacement-operator for the k-side of the svd-expansion applied to M. ;
      complex *16, allocatable :: tmp_Z_q_(:) !complex *16 array, size max(n_gamma_z_upd,n_gamma_z_nuf). ;
      complex *16, allocatable :: tmp_Z_p_(:) !complex *16 array, size max(n_gamma_z_upd,n_gamma_z_nuf). ;
      complex *16, allocatable :: form_CTF_R_S_0__(:) !complex *16 array, size n_gamma_z_upd. ;
      complex *16, allocatable :: quad_CTF_R_S_0__(:) !complex *16 array, size n_gamma_z_upd. ;
      complex *16, allocatable :: quad_Y_CTF_R_S_build_uniform__(:) !complex *16 array, size n_gamma_z_upd. ;
      complex *16, allocatable :: quad_Y_CTF_R_S_brute_uniform__(:) !complex *16 array, size n_gamma_z_upd. ;
      complex *16, allocatable :: quad_Y_CTF_R_S_build_nonunif__(:) !complex *16 array, size n_gamma_z_nuf. ;
      complex *16, allocatable :: quad_Y_CTF_R_S_brute_nonunif__(:) !complex *16 array, size n_gamma_z_nuf. ;
      logical flag_RTRT_vs_RTTR
      logical flag_memory_checkset
      real *8 pi
      complex *16 c16_alpha,c16_beta
      real *8 dr,dtmp
      real *8 finufft_eps
      integer finufft_iflag,finufft_ier
      integer na,nx_c_0,nx_c_1,nw,nr,ngz
      real *8 tmp_k,tmp_w
      real *8 tmp_gamma_z
      real *8 fnod_c16_f,al2_c16_f,sum_r8_f !function output. ;
      complex *16 dotw_c16_f !function output. ;
      real *8 I_P_f,I_xP_f,I_PxP_f,I_PxxP_f !function output. ;
      real *8 I_PxRTRTxP_f !function output. ;

      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_dr_digits_single]'
      end if !if (verbose.gt.0) then
      if (verbose.gt.2) then
         write(6,'(A,F24.16)') ' half_diameter_k_c: ' ,
     $        half_diameter_k_c
         write(6,'(A,F24.16)') ' fix_phi_M: ' , fix_phi_M
         write(6,'(A,F24.16)') ' fix_delta_M_x_p_r: ' ,
     $        fix_delta_M_x_p_r
         write(6,'(A,F24.16)') ' fix_delta_M_x_p_w: ' ,
     $        fix_delta_M_x_p_w
         write(6,'(A,F24.16)') ' fix_phi_S: ' , fix_phi_S
         write(6,'(A,F24.16)') ' fix_delta_S_x_p_r: ' ,
     $        fix_delta_S_x_p_r
         write(6,'(A,F24.16)') ' fix_delta_S_x_p_w: ' ,
     $        fix_delta_S_x_p_w
         write(6,'(A,L2)') ' flag_RTRT_vs_RTTR: ' ,
     $        flag_RTRT_vs_RTTR
         write(6,'(A,F24.16)') ' gamma_z_est: ' , gamma_z_est
         write(6,'(A,F24.16)') ' delta_est_x_c_0: ' , delta_est_x_c_0
         write(6,'(A,F24.16)') ' delta_est_x_c_1: ' , delta_est_x_c_1
         write(6,'(A,F24.16)') ' gamma_z_upd_0in: ' , gamma_z_upd_0in
         write(6,'(A,F24.16)') ' delta_upd_x_c_0: ' , delta_upd_x_c_0
         write(6,'(A,F24.16)') ' delta_upd_x_c_1: ' , delta_upd_x_c_1
      end if !if (verbose.gt.2) then
      
      pi = 4.0d0*datan(1.0d0)
      diameter_x_c_0 = 2.0d0
      diameter_x_c_1 = 2.0d0
      half_diameter_x_c_0 = diameter_x_c_0/2.0d0
      half_diameter_x_c_1 = diameter_x_c_1/2.0d0
      diameter_k_c = 2.0d0*half_diameter_k_c
      quad_n_r = 48
      allocate(jacpts_x_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,jacpts_x_)
      allocate(jacpts_w_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,jacpts_w_)
      call jacpts_48(quad_n_r,jacpts_x_,jacpts_w_)
      if (verbose.gt.2) then
      call print_sub_r8(quad_n_r,jacpts_x_,' jacpts_x_: ')
      call print_sub_r8(quad_n_r,jacpts_w_,' jacpts_w_: ')
      end if !if (verbose.gt.2) then
      allocate(quad_grid_k_p_r_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,quad_grid_k_p_r_)
      allocate(quad_weight_k_p_r_(0:1+quad_n_r-1))
      call cs1_r8(quad_n_r,quad_weight_k_p_r_)
      call af1_r8(quad_n_r,0.5d0*half_diameter_k_c,0.5d0
     $     *half_diameter_k_c,jacpts_x_,quad_grid_k_p_r_)
      call af1_r8(quad_n_r,half_diameter_k_c**2 / sum_r8_f(quad_n_r
     $     ,jacpts_w_),0.0d0,jacpts_w_,quad_weight_k_p_r_)
      if (verbose.gt.2) then
      call print_sub_r8(quad_n_r,quad_grid_k_p_r_,' quad_grid_k_p_r_: ')
      call print_sub_r8(quad_n_r,quad_weight_k_p_r_
     $     ,' quad_weight_k_p_r_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_n_polar_a_(0:1+quad_n_r-1))
      call cs1_i4(quad_n_r,quad_n_polar_a_)
      do nr=0,quad_n_r-1
         quad_n_polar_a_(nr) = ceiling(pi*(quad_grid_k_p_r_(nr)+1))
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_i4(quad_n_r,quad_n_polar_a_,' quad_n_polar_a_: ')
      end if !if (verbose.gt.2) then
      allocate(quad_n_w_(0:1+quad_n_r-1))
      call cs1_i4(quad_n_r,quad_n_w_)
      quad_n_w_sum = 0
      do nr=0,quad_n_r-1
         quad_n_w_(nr) = 2*quad_n_polar_a_(nr)
         quad_n_w_sum = quad_n_w_sum + quad_n_w_(nr)
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_i4(quad_n_r,quad_n_w_,' quad_n_w_: ')
      end if !if (verbose.gt.2) then
      quad_n_w_max = quad_n_w_(quad_n_r-1)
      ngamma_z_upd_0in = nint(gamma_z_upd_0in*(1.0d0*quad_n_w_max)
     $     /(2.0d0*pi))
      if (verbose.gt.2) then
         write(6,'(A,I0)') ' ngamma_z_upd_0in: ' , ngamma_z_upd_0in
         write(6,'(A,F24.16)') ' gamma_z_upd_0in: ' , (2.0d0*pi
     $        *ngamma_z_upd_0in)/(1.0d0*quad_n_w_max)
      end if !if (verbose.gt.2) then
      allocate(quad_grid_k_p_0_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_grid_k_p_0_)
      allocate(quad_grid_k_p_1_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_grid_k_p_1_)
      allocate(quad_weight_k_p_(0:1+quad_n_w_sum-1))
      call cs1_r8(quad_n_w_sum,quad_weight_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = 2.0d0*pi*nw/quad_n_w_(nr)
            quad_grid_k_p_0_(na) = pi*tmp_k/half_diameter_k_c
     $           *dcos(tmp_w)
            quad_grid_k_p_1_(na) = pi*tmp_k/half_diameter_k_c
     $           *dsin(tmp_w)
            quad_weight_k_p_(na) = quad_weight_k_p_r_(nr) * pi /
     $           quad_n_w_(nr)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
         call print_sub_r8(quad_n_w_sum,quad_grid_k_p_0_
     $        ,' quad_grid_k_p_0_: ')
      call print_sub_r8(quad_n_w_sum,quad_grid_k_p_1_
     $     ,' quad_grid_k_p_1_: ')
      call print_sub_r8(quad_n_w_sum,quad_weight_k_p_
     $     ,' quad_weight_k_p_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(A)') ' Testing integration of polynomial '
      end if !if (verbose.gt.2) then
      tmp_p = 3.5d0
      tmp_a = 1.2d0
      form_I = 2*pi*tmp_a*half_diameter_k_c**(tmp_p+2)/(tmp_p+2.0d0)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_I: ' , form_I
      end if !if (verbose.gt.2) then
      allocate(quad_P_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = tmp_k**tmp_p * ( tmp_a + 0.5d0*cos(1.0d0*tmp_w) +
     $           1.5*cos(2.0d0*tmp_w) + 2.5d0*cos(-3.0d0*tmp_w) )
            quad_P_k_p_(na) = dcmplx(tmp_e,0.0d0)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_k_p_,' quad_P_k_p_: ')
      end if !if (verbose.gt.2) then
      quad_I = dcmplx(0.0d0,0.0d0)
      do na=0,quad_n_w_sum-1
         quad_I = quad_I + quad_P_k_p_(na)*quad_weight_k_p_(na)
      enddo !do na=0,quad_n_w_sum-1
      if (verbose.gt.2) then
         write(6,'(A,F21.16,A,F21.16,A)') ' quad_I error: ( '
     $        ,real((quad_I-form_I)/form_I),' , ' , aimag((quad_I
     $        -form_I)/form_I) , ' )'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(A)') ' Testing integration of 1 plane-wave '
      end if !if (verbose.gt.2) then
      allocate(fix_delta_x_c_(0:1+2-1))
      call cs1_r8(2,fix_delta_x_c_)
      fix_delta_x_c_(0) = 4.0d0*0.85d0/half_diameter_k_c
      fix_delta_x_c_(1) = 4.0d0*0.15d0/half_diameter_k_c
      fix_delta_x_p_r = dsqrt(fix_delta_x_c_(0)**2 + fix_delta_x_c_(1)
     $     **2)
      fix_delta_x_p_w = datan2(fix_delta_x_c_(1),fix_delta_x_c_(0))
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_x_p_r: ' , fix_delta_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_x_p_w: ' , fix_delta_x_p_w
      end if !if (verbose.gt.2) then
      form_I = I_P_f(half_diameter_k_c,fix_delta_x_p_r)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_I: ' , form_I
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      call cl1_c16(quad_n_w_sum,quad_P_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_x_p_r*dcos(tmp_w
     $           -fix_delta_x_p_w)
            quad_P_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_k_p_,' quad_P_k_p_: ')
      end if !if (verbose.gt.2) then
      quad_I = dcmplx(0.0d0,0.0d0)
      do na=0,quad_n_w_sum-1
         quad_I = quad_I + quad_P_k_p_(na)*quad_weight_k_p_(na)
      enddo !do na=0,quad_n_w_sum-1
      write(6,'(A,F21.16,A,F21.16,A)') ' quad_I error: ( ',real((quad_I
     $     -form_I)/form_I),' , ' , aimag((quad_I-form_I)/form_I) , ' )'
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(A)') ' Testing integration of 2 plane-waves: '
      end if !if (verbose.gt.2) then
      allocate(fix_delta_S_x_c_(0:1+2-1))
      call cs1_r8(2,fix_delta_S_x_c_)
      allocate(fix_delta_M_x_c_(0:1+2-1))
      call cs1_r8(2,fix_delta_M_x_c_)
      fix_delta_S_x_c_(0) = fix_delta_S_x_p_r * cos(fix_delta_S_x_p_w)
      fix_delta_S_x_c_(1) = fix_delta_S_x_p_r * sin(fix_delta_S_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_M_x_c_(0) = fix_delta_M_x_p_r * cos(fix_delta_M_x_p_w)
      fix_delta_M_x_c_(1) = fix_delta_M_x_p_r * sin(fix_delta_M_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_x_p_r = dsqrt((fix_delta_S_x_c_(0)-fix_delta_M_x_c_(0))
     $     **2 +(fix_delta_S_x_c_(1)-fix_delta_M_x_c_(1))**2)
      form_I_SM = I_P_f(half_diameter_k_c,fix_delta_x_p_r)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_I_SM: ' , form_I_SM
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_P_S_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_S_k_p_)
      allocate(quad_P_S_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_S_k_q_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_S_x_p_r*dcos(tmp_w
     $           -fix_delta_S_x_p_w)
            quad_P_S_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_S_k_p_,' quad_P_S_k_p_: ')
      end if !if (verbose.gt.2) then
      allocate(quad_P_M_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_M_k_p_)
      allocate(quad_P_M_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_M_k_q_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_M_x_p_r*dcos(tmp_w
     $           -fix_delta_M_x_p_w)
            quad_P_M_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_M_k_p_,' quad_P_M_k_p_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      quad_I_SM = dotw_c16_f(quad_n_w_sum,quad_P_S_k_p_,quad_P_M_k_p_
     $     ,quad_weight_k_p_)
      if (verbose.gt.2) then
      write(6,'(A,F21.16,A,F21.16,A)') ' quad_I_SM error: ( '
     $     ,real((quad_I_SM-form_I_SM)/form_I_SM),' , ' ,
     $     aimag((quad_I_SM -form_I_SM)/form_I_SM) , ' )'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(A)') ' Testing integration of 3 plane-waves: '
      end if !if (verbose.gt.2) then
      allocate(fix_delta_C_x_c_(0:1+2-1))
      call cs1_r8(2,fix_delta_C_x_c_)
      fix_delta_C_x_c_(0) = +4.0d0*0.25d0/half_diameter_k_c
      fix_delta_C_x_c_(1) = -4.0d0*0.15d0/half_diameter_k_c
      fix_delta_C_x_p_r = dsqrt(fix_delta_C_x_c_(0)**2 +
     $     fix_delta_C_x_c_(1)**2)
      fix_delta_C_x_p_w = datan2(fix_delta_C_x_c_(1)
     $     ,fix_delta_C_x_c_(0))
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_C_x_p_r: ' , fix_delta_C_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_C_x_p_w: ' , fix_delta_C_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_x_p_r = dsqrt((fix_delta_S_x_c_(0)+fix_delta_C_x_c_(0)
     $     -fix_delta_M_x_c_(0))**2 +(fix_delta_S_x_c_(1)
     $     +fix_delta_C_x_c_(1)-fix_delta_M_x_c_(1))**2)
      form_I_SCM = I_P_f(half_diameter_k_c,fix_delta_x_p_r)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_I_SCM: ' , form_I_SCM
      end if !if (verbose.gt.2) then
      allocate(quad_P_C_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_C_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_C_x_p_r*dcos(tmp_w
     $           -fix_delta_C_x_p_w)
            quad_P_C_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_C_k_p_,' quad_P_C_k_p_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_P_N_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_N_k_p_)
      allocate(quad_P_T_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_T_k_p_)
      allocate(quad_P_T_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_T_k_q_)
      allocate(quad_P_R_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_R_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      quad_I_SCM = dotw_c16_f(quad_n_w_sum,quad_P_S_k_p_,quad_P_N_k_p_
     $     ,quad_weight_k_p_)
      if (verbose.gt.2) then
      write(6,'(A,F21.16,A,F21.16,A)') ' quad_I_SCM error: ( '
     $     ,real((quad_I_SCM-form_I_SCM)/form_I_SCM),' , ' ,
     $     aimag((quad_I_SCM -form_I_SCM)/form_I_SCM) , ' )'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(2A)') ' Testing integration of 1 plane-wave'
     $     ,' multiplied by a planar function: '
      end if !if (verbose.gt.2) then
      fix_phi = fix_phi_M
      call cl1_r8(2,fix_delta_x_c_)
      fix_delta_x_c_(0) = 4.0d0*0.85d0/half_diameter_k_c
      fix_delta_x_c_(1) = 4.0d0*0.15d0/half_diameter_k_c
      fix_delta_x_p_r = dsqrt(fix_delta_x_c_(0)**2 + fix_delta_x_c_(1)
     $     **2)
      fix_delta_x_p_w = datan2(fix_delta_x_c_(1),fix_delta_x_c_(0))
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_x_p_r: ' , fix_delta_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_x_p_w: ' , fix_delta_x_p_w
      end if !if (verbose.gt.2) then
      form_I = I_xP_f(half_diameter_k_c,fix_phi,fix_delta_x_p_r
     $     ,fix_delta_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_I: ' , form_I
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      call cl1_c16(quad_n_w_sum,quad_P_k_p_)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_x_p_r*dcos(tmp_w
     $           -fix_delta_x_p_w)
            tmp_f = -2*tmp_k*dcos(tmp_w - fix_phi)
            quad_P_k_p_(na) = dcmplx(0.0d0,tmp_f) * dcmplx(+dcos(tmp_e),
     $           +dsin(tmp_e))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_k_p_,' quad_P_k_p_: ')
      end if !if (verbose.gt.2) then
      quad_I = dcmplx(0.0d0,0.0d0)
      do na=0,quad_n_w_sum-1
         quad_I = quad_I + quad_P_k_p_(na)*quad_weight_k_p_(na)
      enddo !do na=0,quad_n_w_sum-1
      if (verbose.gt.2) then
         write(6,'(A,F21.16,A,F21.16,A)') ' quad_I error: ( '
     $        ,real((quad_I-form_I)/form_I),' , ' , aimag((quad_I
     $        -form_I)/form_I) , ' )'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(2A)') ' Testing integration of 2 plane-waves'
     $     ,' multiplied by a planar function: '
      end if !if (verbose.gt.2) then
      fix_phi = fix_phi_M
      fix_delta_S_x_c_(0) = fix_delta_S_x_p_r * cos(fix_delta_S_x_p_w)
      fix_delta_S_x_c_(1) = fix_delta_S_x_p_r * sin(fix_delta_S_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_M_x_c_(0) = fix_delta_M_x_p_r * cos(fix_delta_M_x_p_w)
      fix_delta_M_x_c_(1) = fix_delta_M_x_p_r * sin(fix_delta_M_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w
      end if !if (verbose.gt.2) then
      form_I_SCM = I_PxP_f(half_diameter_k_c,fix_phi,fix_delta_M_x_p_r
     $     ,fix_delta_M_x_p_w,fix_delta_S_x_p_r,fix_delta_S_x_p_w)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_S_x_p_r*dcos(tmp_w
     $           -fix_delta_S_x_p_w)
            quad_P_S_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_M_x_p_r*dcos(tmp_w
     $           -fix_delta_M_x_p_w)
            quad_P_M_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            tmp_f = -2*tmp_k*dcos(tmp_w - fix_phi)
            quad_P_C_k_p_(na) = dcmplx(0.0d0,tmp_f)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_S_k_p_,' quad_P_S_k_p_: ')
      call print_sub_c16(quad_n_w_sum,quad_P_M_k_p_,' quad_P_M_k_p_: ')
      call print_sub_c16(quad_n_w_sum,quad_P_C_k_p_,' quad_P_C_k_p_: ')
      end if !if (verbose.gt.2) then
      call xc1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      quad_I_SCM = dotw_c16_f(quad_n_w_sum,quad_P_S_k_p_,quad_P_N_k_p_
     $     ,quad_weight_k_p_)
      if (verbose.gt.2) then
      write(6,'(A,F21.16,A,F21.16,A)') ' quad_I_SCM error: ( '
     $     ,real((quad_I_SCM-form_I_SCM)/form_I_SCM),' , ' ,
     $     aimag((quad_I_SCM -form_I_SCM)/form_I_SCM) , ' )'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(2A)') ' Testing integration of 2 plane-waves'
     $     ,' each multiplied by a planar function: '
      end if !if (verbose.gt.2) then
c$$$      fix_phi_S read as input. ;
      fix_delta_S_x_c_(0) = fix_delta_S_x_p_r * cos(fix_delta_S_x_p_w)
      fix_delta_S_x_c_(1) = fix_delta_S_x_p_r * sin(fix_delta_S_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_M_x_c_(0) = fix_delta_M_x_p_r * cos(fix_delta_M_x_p_w)
      fix_delta_M_x_c_(1) = fix_delta_M_x_p_r * sin(fix_delta_M_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w
      end if !if (verbose.gt.2) then
c$$$      fix_phi_M read from input. ;
      form_I_SCM = I_PxxP_f(half_diameter_k_c,fix_phi_M
     $     ,fix_delta_M_x_p_r,fix_delta_M_x_p_w,fix_phi_S
     $     ,fix_delta_S_x_p_r,fix_delta_S_x_p_w)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_S_x_p_r*dcos(tmp_w
     $           -fix_delta_S_x_p_w)
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_S)
            quad_P_S_k_p_(na) = dcmplx(tmp_f,0.0d0) * dcmplx(
     $           +dcos(tmp_e),+dsin(tmp_e))
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_M_x_p_r*dcos(tmp_w
     $           -fix_delta_M_x_p_w)
            quad_P_M_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_M)
            quad_P_C_k_p_(na) = dcmplx(tmp_f,0.0d0)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_P_S_k_p_,' quad_P_S_k_p_: ')
      call print_sub_c16(quad_n_w_sum,quad_P_M_k_p_,' quad_P_M_k_p_: ')
      call print_sub_c16(quad_n_w_sum,quad_P_C_k_p_,' quad_P_C_k_p_: ')
      end if !if (verbose.gt.2) then
      call xc1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      quad_I_SCM = dotw_c16_f(quad_n_w_sum,quad_P_S_k_p_,quad_P_N_k_p_
     $     ,quad_weight_k_p_)
      if (verbose.gt.2) then
      write(6,'(A,F21.16,A,F21.16,A)') ' quad_I_SCM error: ( '
     $     ,real((quad_I_SCM-form_I_SCM)/form_I_SCM),' , ' ,
     $     aimag((quad_I_SCM -form_I_SCM)/form_I_SCM) , ' )'
      end if !if (verbose.gt.2) then
c$$$  %%%%%%%%
      allocate(quad_fftw_plan_frwd_(0:1+quad_n_r-1))
      call cs1_i8(quad_n_r,quad_fftw_plan_frwd_)
      allocate(quad_fftw_plan_back_(0:1+quad_n_r-1))
      call cs1_i8(quad_n_r,quad_fftw_plan_back_)
      allocate(quad_fftw_0in_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_fftw_0in_)
      allocate(quad_fftw_out_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_fftw_out_)
      if (verbose.gt.2) then
      write(6,'(A)') ' entering quad_fft'
      end if !if (verbose.gt.2) then
      na=0
      do nr=0,quad_n_r-1
         call dfftw_plan_dft_1d_(quad_fftw_plan_frwd_(nr),quad_n_w_(nr)
     $        ,quad_fftw_0in_(na),quad_fftw_out_(na),FFTW_FORWARD
     $        ,FFTW_ESTIMATE)
         call dfftw_plan_dft_1d_(quad_fftw_plan_back_(nr),quad_n_w_(nr)
     $        ,quad_fftw_out_(na),quad_fftw_0in_(na),FFTW_BACKWARD
     $        ,FFTW_ESTIMATE)
         na = na + quad_n_w_(nr)
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      write(6,'(A)') ' finished quad_fft'
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      n_delta_v = 1
c$$$      delta_est_x_c_0 read from input. ;
c$$$      delta_est_x_c_1 read from input. ;
c$$$      gamma_z_est read from input. ;
c$$$      delta_upd_x_c_0 read from input. ;
c$$$      delta_upd_x_c_1 read from input. ;
c$$$      %%%%%%%%
      write(6,'(A)') ' First: testing SZxTRM: '
c$$$      fix_phi_S read as input. ;
      fix_delta_S_x_c_(0) = fix_delta_S_x_p_r * cos(fix_delta_S_x_p_w)
      fix_delta_S_x_c_(1) = fix_delta_S_x_p_r * sin(fix_delta_S_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_0: ' , fix_delta_S_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_S_x_c_1: ' , fix_delta_S_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_r: ' , fix_delta_S_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_S_x_p_w: ' , fix_delta_S_x_p_w
      end if !if (verbose.gt.2) then
      fix_delta_M_x_c_(0) = fix_delta_M_x_p_r * cos(fix_delta_M_x_p_w)
      fix_delta_M_x_c_(1) = fix_delta_M_x_p_r * sin(fix_delta_M_x_p_w)
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_0: ' , fix_delta_M_x_c_(0)
      write(6,'(A,F21.16)') ' fix_delta_M_x_c_1: ' , fix_delta_M_x_c_(1)
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_r: ' , fix_delta_M_x_p_r
      write(6,'(A,F21.16)') ' fix_delta_M_x_p_w: ' , fix_delta_M_x_p_w
      end if !if (verbose.gt.2) then
c$$$      fix_phi_M read from input. ;
      form_I_SCM = I_PxxP_f(half_diameter_k_c,fix_phi_M
     $     ,fix_delta_M_x_p_r,fix_delta_M_x_p_w,fix_phi_S
     $     ,fix_delta_S_x_p_r,fix_delta_S_x_p_w)
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_S_x_p_r*dcos(tmp_w
     $           -fix_delta_S_x_p_w)
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_S)
            quad_P_S_k_p_(na) = dcmplx(tmp_f,0.0d0) * dcmplx(
     $           +dcos(tmp_e),+dsin(tmp_e))
            tmp_e = 2.0d0*pi*tmp_k*fix_delta_M_x_p_r*dcos(tmp_w
     $           -fix_delta_M_x_p_w)
            quad_P_M_k_p_(na) = dcmplx(+dcos(tmp_e),+dsin(tmp_e))
            tmp_f = 2*tmp_k*dcos(tmp_w - fix_phi_M)
            quad_P_C_k_p_(na) = dcmplx(tmp_f,0.0d0)
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      allocate(quad_X1_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_X1_)
      allocate(quad_P_N_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_N_k_q_)
      allocate(quad_P_R_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_P_R_k_q_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_S_k_p_,+delta_upd_x_c_0,+delta_upd_x_c_1
     $     ,quad_P_T_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_T_k_p_
     $     ,quad_P_R_k_q_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_N_k_p_
     $     ,quad_P_N_k_q_)
      call innerproduct_q_k_stretch_quad_0(quad_n_r,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum,quad_P_R_k_q_
     $     ,quad_P_N_k_q_,quad_X1_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0) 
      c16_beta = dcmplx(0.0d0,0.0d0)
      call af1_c16(quad_n_w_max,c16_alpha,c16_beta
     $     ,quad_X1_,quad_X1_)
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_X1_,' quad_X1_: ')
      end if !if (verbose.gt.2) then
      call cp1_c16(quad_n_w_max,quad_X1_,quad_fftw_out_(quad_n_w_sum
     $     -quad_n_w_max))
      call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
      c16_alpha = dcmplx(1.0d0,0.0d0)
      c16_beta = dcmplx(0.0d0,0.0d0)
      call af1_c16(quad_n_w_max,c16_alpha,c16_beta
     $     ,quad_fftw_0in_(quad_n_w_sum-quad_n_w_max),quad_X1_)
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_X1_,' quad_X1_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_X0_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_X0_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_S_k_p_,+delta_upd_x_c_0,+delta_upd_x_c_1
     $     ,quad_P_T_k_p_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      do nw=0,quad_n_w_max-1
         gamma_z_upd = (2.0d0*pi*nw)/(1.0d0*quad_n_w_max)
         call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $        ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum
     $        ,quad_fftw_0in_ ,quad_fftw_out_,quad_P_T_k_p_,
     $        +gamma_z_upd,quad_P_R_k_p_)
         quad_X0_(nw) = dotw_c16_f(quad_n_w_sum,quad_P_R_k_p_
     $        ,quad_P_N_k_p_,quad_weight_k_p_)
      enddo !do nw=0,quad_n_w_max-1
      if (verbose.gt.2) then
      write(6,'(A,2F24.16)') ' quad_X0_(ngamma_z_upd_0in): ' ,
     $     quad_X0_(ngamma_z_upd_0in)
      call print_sub_c16(quad_n_w_max,quad_X0_,' quad_X0_: ')
      end if !if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' quad_X1 vs quad_X0_: '
     $     ,fnod_c16_f(quad_n_w_max,quad_X1_,quad_X0_)
     $     /al2_c16_f(quad_n_w_max,quad_X0_)
c$$$      %%%%%%%%
      allocate(tmp_delta_S_x_c_(0:1+2-1))
      call cs1_r8(2,tmp_delta_S_x_c_)
      allocate(tmp_delta_M_x_c_(0:1+2-1))
      call cs1_r8(2,tmp_delta_M_x_c_)
      allocate(tmp_delta_Y_x_c_(0:1+2-1))
      call cs1_r8(2,tmp_delta_Y_x_c_)
      allocate(form_X0_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,form_X0_)
      allocate(form_SX0_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,form_SX0_)
      allocate(form_MX0_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,form_MX0_)
      do nw=0,quad_n_w_max-1
         gamma_z_upd = (2.0d0*pi*nw)/(1.0d0*quad_n_w_max)
         tmp_delta_S_x_c_(0) = fix_delta_S_x_c_(0) - delta_upd_x_c_0
         tmp_delta_S_x_c_(1) = fix_delta_S_x_c_(1) - delta_upd_x_c_1
         call rotate_delta(+gamma_z_upd,tmp_delta_S_x_c_
     $        ,tmp_delta_S_x_c_)
         tmp_delta_S_x_p_r = dsqrt(tmp_delta_S_x_c_(0)**2 +
     $        tmp_delta_S_x_c_(1)**2)
         tmp_delta_S_x_p_w = datan2(tmp_delta_S_x_c_(1)
     $        ,tmp_delta_S_x_c_(0))
         tmp_delta_M_x_c_(0) = fix_delta_M_x_c_(0)
         tmp_delta_M_x_c_(1) = fix_delta_M_x_c_(1)
         call rotate_delta(-gamma_z_est,tmp_delta_M_x_c_
     $        ,tmp_delta_M_x_c_)
         tmp_delta_M_x_c_(0) = tmp_delta_M_x_c_(0) + delta_est_x_c_0
         tmp_delta_M_x_c_(1) = tmp_delta_M_x_c_(1) + delta_est_x_c_1
         tmp_delta_M_x_p_r = dsqrt(tmp_delta_M_x_c_(0)**2 +
     $        tmp_delta_M_x_c_(1)**2)
         tmp_delta_M_x_p_w = datan2(tmp_delta_M_x_c_(1)
     $        ,tmp_delta_M_x_c_(0))
         tmp_phi_M = fix_phi_M - gamma_z_est
         tmp_phi_S = fix_phi_S + gamma_z_upd
         form_X0_(nw) = I_PxxP_f(half_diameter_k_c,tmp_phi_M
     $        ,tmp_delta_M_x_p_r,tmp_delta_M_x_p_w,tmp_phi_S
     $        ,tmp_delta_S_x_p_r ,tmp_delta_S_x_p_w)
         form_SX0_(nw) = form_X0_(nw)
         form_SX0_(nw) = I_PxRTRTxP_f(
     $              half_diameter_k_c
     $              ,fix_phi_M
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,.true.
     $              ,gamma_z_est
     $              ,delta_est_x_c_0
     $              ,delta_est_x_c_1
     $              ,gamma_z_upd
     $              ,delta_upd_x_c_0
     $              ,delta_upd_x_c_1
     $              )
         tmp_delta_S_x_c_(0) = fix_delta_S_x_c_(0)
         tmp_delta_S_x_c_(1) = fix_delta_S_x_c_(1)
         call rotate_delta(+gamma_z_upd,tmp_delta_S_x_c_
     $        ,tmp_delta_S_x_c_)
         tmp_delta_S_x_c_(0) = tmp_delta_S_x_c_(0) - delta_upd_x_c_0
         tmp_delta_S_x_c_(1) = tmp_delta_S_x_c_(1) - delta_upd_x_c_1
         tmp_delta_S_x_p_r = dsqrt(tmp_delta_S_x_c_(0)**2 +
     $        tmp_delta_S_x_c_(1)**2)
         tmp_delta_S_x_p_w = datan2(tmp_delta_S_x_c_(1)
     $        ,tmp_delta_S_x_c_(0))
         form_MX0_(nw) = I_PxxP_f(half_diameter_k_c,tmp_phi_M
     $        ,tmp_delta_M_x_p_r,tmp_delta_M_x_p_w,tmp_phi_S
     $        ,tmp_delta_S_x_p_r ,tmp_delta_S_x_p_w)
         form_MX0_(nw) = I_PxRTRTxP_f(
     $              half_diameter_k_c
     $              ,fix_phi_M
     $              ,fix_delta_M_x_p_r
     $              ,fix_delta_M_x_p_w
     $              ,fix_phi_S
     $              ,fix_delta_S_x_p_r
     $              ,fix_delta_S_x_p_w
     $              ,.false.
     $              ,gamma_z_est
     $              ,delta_est_x_c_0
     $              ,delta_est_x_c_1
     $              ,gamma_z_upd
     $              ,delta_upd_x_c_0
     $              ,delta_upd_x_c_1
     $              )
      enddo !do nw=0,quad_n_w_max-1
      if (verbose.gt.1) then
      write(6,'(A,2F24.16)') ' form_SX0_(ngamma_z_upd_0in): ' ,
     $     form_SX0_(ngamma_z_upd_0in)
      write(6,'(A,2F24.16)') ' form_MX0_(ngamma_z_upd_0in): ' ,
     $     form_MX0_(ngamma_z_upd_0in)
      end if !if (verbose.gt.1) then
      if (verbose.gt.2) then
      write(6,'(A,2F24.16)') ' form_X0_(ngamma_z_upd_0in): ' ,
     $     form_X0_(ngamma_z_upd_0in)
      call print_sub_c16(quad_n_w_max,form_X0_,' form_X0_: ')
      call print_sub_c16(quad_n_w_max,form_SX0_,' form_SX0_: ')
      call print_sub_c16(quad_n_w_max,form_MX0_,' form_MX0_: ')
      end if !if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' form_X0 vs quad_X0_: '
     $     ,fnod_c16_f(quad_n_w_max,form_X0_,quad_X0_)
     $     /al2_c16_f(quad_n_w_max,form_X0_)
c$$$      %%%%%%%%
      eps_svd = 1.0d-6
      R_max = 2.0d0*pi*quad_grid_k_p_r_(quad_n_r-1)
      K_max = quad_grid_k_p_r_(quad_n_r-1)
      n_pixels = 4.50d0
      delta_max = n_pixels*dsqrt(2.0d0)/2.0d0/K_max
      tmp_delta_x_c_0_(0) = delta_max
      tmp_delta_x_c_1_(0) = 0.0d0
      allocate(svd_r_(0:1+n_svd_max-1)) !holds jacobi nodes for k ;
      call cs1_r8(n_svd_max,svd_r_)
      allocate(svd_d_(0:1+n_svd_max-1)) !holds jacobi nodes for delta ;
      call cs1_r8(n_svd_max,svd_d_)
      allocate(svd_l_(0:1+n_svd_max-1))
      call cs1_i4(n_svd_max,svd_l_)
      allocate(svd_U_d_(0:1+n_svd_max*n_svd_max-1))
      call cs1_r8(n_svd_max*n_svd_max,svd_U_d_)
      allocate(svd_s_(0:1+n_svd_max-1))
      call cs1_r8(n_svd_max,svd_s_)
      allocate(svd_V_r_(0:1+n_svd_max*n_svd_max-1))
      call cs1_r8(n_svd_max*n_svd_max,svd_V_r_)
      call get_svd_2(eps_svd,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,quad_grid_k_p_r_,quad_n_r,1,tmp_delta_x_c_0_
     $     ,tmp_delta_x_c_1_,flag_warning ,R_max ,K_max,delta_max
     $     ,n_pixels)
      if (verbose.gt.2) then
      write(6,'(A,I0)') ' n_svd_d: ' , n_svd_d
      write(6,'(A,I0)') ' n_svd_r: ' , n_svd_r
      write(6,'(A,I0)') ' n_svd_l: ' , n_svd_l
      if (verbose.gt.2) then
      call print_sub_i4(n_svd_l,svd_l_,' svd_l_: ')
      call print_sub_r8(n_svd_l*n_svd_d,svd_U_d_,' svd_U_d_: ')
      call print_sub_r8(n_svd_l,svd_s_,' svd_s_: ')
      call print_sub_r8(n_svd_l*n_svd_r,svd_V_r_,' svd_V_r_: ')
      end if !if (verbose.gt.2) then
      end if !if (verbose.gt.2) then
      allocate(svd_polyval_U_d_(0:1+n_svd_l*n_delta_v-1))
      call cs1_r8(n_svd_l*n_delta_v,svd_polyval_U_d_)
      allocate(svd_polyval_V_r_(0:1+n_svd_l*quad_n_r-1))
      call cs1_r8(n_svd_l*quad_n_r,svd_polyval_V_r_)
      svd_d_max = delta_max
      call get_svd_polyval_U_d_(svd_d_max,n_svd_d,svd_d_,n_svd_l ,svd_l_
     $     ,svd_U_d_,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1
     $     ,svd_polyval_U_d_)
      svd_r_max = quad_grid_k_p_r_(quad_n_r-1)
      call get_svd_polyval_V_r_(svd_r_max,n_svd_r,svd_r_,n_svd_l ,svd_l_
     $     ,svd_V_r_,quad_n_r,quad_grid_k_p_r_ ,svd_polyval_V_r_)
      if (verbose.gt.2) then
      call print_sub_r8(n_svd_l*n_delta_v,svd_polyval_U_d_
     $     ,' svd_polyval_U_d_: ')
      call print_sub_r8(n_svd_l*quad_n_r,svd_polyval_V_r_
     $     ,' svd_polyval_V_r_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_X2_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_X2_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_S_k_p_
     $     ,quad_P_S_k_q_)
      call transf_svd_q_to_q_polyval_0on_0(svd_r_max,n_svd_r,svd_r_
     $     ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_s_
     $     ,svd_V_r_,quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,+delta_upd_x_c_0,+delta_upd_x_c_1
     $     ,quad_P_R_k_q_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_N_k_p_
     $     ,quad_P_N_k_q_)
      call innerproduct_q_k_stretch_quad_0(quad_n_r,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_,quad_n_w_,quad_n_w_sum,quad_P_R_k_q_
     $     ,quad_P_N_k_q_,quad_X2_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(quad_n_w_max,c16_alpha,c16_beta
     $     ,quad_X2_,quad_X2_)
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_X2_,' quad_X2_: ')
      end if !if (verbose.gt.2) then
      call cp1_c16(quad_n_w_max,quad_X2_,quad_fftw_out_(quad_n_w_sum
     $     -quad_n_w_max))
      call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
      c16_alpha = dcmplx(1.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(quad_n_w_max,c16_alpha,c16_beta
     $     ,quad_fftw_0in_(quad_n_w_sum-quad_n_w_max),quad_X2_)
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_X2_,' quad_X2_: ')
      end if !if (verbose.gt.2) then
      if (verbose.gt.2) then
      write(6,'(A,F8.6,A,F21.16)') ' eps_svd: ' , eps_svd ,
     $     '; form_X0 vs quad_X2_: ',fnod_c16_f(quad_n_w_max,form_X0_
     $     ,quad_X2_)/al2_c16_f(quad_n_w_max,form_X0_)
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      allocate(quad_X3_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_X3_)
      allocate(Z_S_svdd_(0:1+n_svd_l*n_delta_v-1))
      call cs1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
      allocate(Z_M_svdd_(0:1+n_svd_l*n_delta_v-1))
      call cs1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
      allocate(Z_S_svdr_(0:1+n_svd_l*quad_n_w_max-1))
      call cs1_c16(n_svd_l*quad_n_w_max,Z_S_svdr_)
      allocate(Z_M_svdr_(0:1+n_svd_l*quad_n_w_max-1))
      call cs1_c16(n_svd_l*quad_n_w_max,Z_M_svdr_)
      call innerproduct_q_k_svdd_polyval_off_0(.true.,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1,Z_S_svdd_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*n_delta_v,Z_S_svdd_,' Z_S_svdd_: ')
      end if !if (verbose.gt.2) then
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_S_k_p_
     $     ,quad_P_S_k_q_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_N_k_p_
     $     ,quad_P_N_k_q_)
      call innerproduct_q_k_svdr_quad_polyval_0on_0(.true.,svd_r_max
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,quad_n_r
     $     ,quad_grid_k_p_r_,quad_weight_k_p_r_ ,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,quad_P_N_k_q_,Z_S_svdr_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $     ,Z_S_svdr_,Z_S_svdr_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,' Z_S_svdr_: ')
      end if !if (verbose.gt.2) then
      do nsvd_l=0,n_svd_l-1
         call cps_c16(quad_n_w_max,Z_S_svdr_(nsvd_l),n_svd_l
     $        ,quad_fftw_out_(quad_n_w_sum-quad_n_w_max),1)
         call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
         call cps_c16(quad_n_w_max,quad_fftw_0in_(quad_n_w_sum
     $        -quad_n_w_max),1,Z_S_svdr_(nsvd_l),n_svd_l)
      enddo !do nsvd_l=0,n_svd_l-1
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,' Z_S_svdr_: ')
      end if !if (verbose.gt.2) then
      call zgemm('T','N',1,quad_n_w_max,n_svd_l,1.0d0
     $     *dcmplx(1.0d0,0.0d0),Z_S_svdd_
     $     ,n_svd_l,Z_S_svdr_
     $     ,n_svd_l,0.0d0 *dcmplx(1.0d0,0.0d0)
     $     ,quad_X3_
     $     ,1)
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_X3_,' quad_X3_: ')
      end if !if (verbose.gt.2) then
      write(6,'(A,F8.6,A,F21.16)') ' eps_svd: ' , eps_svd ,
     $     '; form_X0 vs quad_X3_: ',fnod_c16_f(quad_n_w_max,form_X0_
     $     ,quad_X3_)/al2_c16_f(quad_n_w_max,form_X0_)
      write(6,'(A,F8.6,A,F21.16)') ' eps_svd: ' , eps_svd ,
     $     '; quad_X2 vs quad_X3_: ',fnod_c16_f(quad_n_w_max,quad_X2_
     $     ,quad_X3_)/al2_c16_f(quad_n_w_max,quad_X2_)
c$$$      %%%%%%%%
      write(6,'(A)') ' Second: testing SxZTRM and SZxTRM: '
      eps_svd = 1.0d-6
      call cl1_r8(n_svd_max,svd_r_)
      call cl1_r8(n_svd_max,svd_d_)
      call cl1_i4(n_svd_max,svd_l_)
      call cl1_r8(n_svd_max*n_svd_max,svd_U_d_)
      call cl1_r8(n_svd_max,svd_s_)
      call cl1_r8(n_svd_max*n_svd_max,svd_V_r_)
      call get_svd_2(eps_svd,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,quad_grid_k_p_r_,quad_n_r,1,tmp_delta_x_c_0_
     $     ,tmp_delta_x_c_1_,flag_warning ,R_max ,K_max,delta_max
     $     ,n_pixels)
      call cl1_r8(n_svd_l*n_delta_v,svd_polyval_U_d_)
      call cl1_r8(n_svd_l*quad_n_r,svd_polyval_V_r_)
      call get_svd_polyval_U_d_(svd_d_max,n_svd_d,svd_d_,n_svd_l ,svd_l_
     $     ,svd_U_d_,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1
     $     ,svd_polyval_U_d_)
      call get_svd_polyval_V_r_(svd_r_max,n_svd_r,svd_r_,n_svd_l ,svd_l_
     $     ,svd_V_r_,quad_n_r,quad_grid_k_p_r_ ,svd_polyval_V_r_)
      call cl1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
      call cl1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
      call cl1_c16(n_svd_l*quad_n_w_max,Z_S_svdr_)
      call cl1_c16(n_svd_l*quad_n_w_max,Z_M_svdr_)
      call innerproduct_q_k_svdd_polyval_off_0(.true.,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1,Z_S_svdd_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*n_delta_v,Z_S_svdd_,' Z_S_svdd_: ')
      end if !if (verbose.gt.2) then
      call innerproduct_q_k_svdd_polyval_off_0(.false.,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1,Z_M_svdd_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*n_delta_v,Z_M_svdd_,' Z_M_svdd_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_S_k_p_
     $     ,quad_P_S_k_q_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_N_k_p_
     $     ,quad_P_N_k_q_)
      call innerproduct_q_k_svdr_quad_polyval_0on_0(.true.,svd_r_max
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,quad_n_r
     $     ,quad_grid_k_p_r_,quad_weight_k_p_r_ ,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,quad_P_N_k_q_,Z_S_svdr_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $     ,Z_S_svdr_,Z_S_svdr_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,' Z_S_svdr_: ')
      end if !if (verbose.gt.2) then
      do nsvd_l=0,n_svd_l-1
         call cps_c16(quad_n_w_max,Z_S_svdr_(nsvd_l),n_svd_l
     $        ,quad_fftw_out_(quad_n_w_sum-quad_n_w_max),1)
         call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
         call cps_c16(quad_n_w_max,quad_fftw_0in_(quad_n_w_sum
     $        -quad_n_w_max),1,Z_S_svdr_(nsvd_l),n_svd_l)
      enddo !do nsvd_l=0,n_svd_l-1
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,' Z_S_svdr_: ')
      end if !if (verbose.gt.2) then
      call innerproduct_q_k_svdr_quad_polyval_0on_0(.false.,svd_r_max
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,quad_n_r
     $     ,quad_grid_k_p_r_,quad_weight_k_p_r_ ,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,quad_P_N_k_q_,Z_M_svdr_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $     ,Z_M_svdr_,Z_M_svdr_)
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_M_svdr_,' Z_M_svdr_: ')
      end if !if (verbose.gt.2) then
      do nsvd_l=0,n_svd_l-1
         call cps_c16(quad_n_w_max,Z_M_svdr_(nsvd_l),n_svd_l
     $        ,quad_fftw_out_(quad_n_w_sum-quad_n_w_max),1)
         call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
         call cps_c16(quad_n_w_max,quad_fftw_0in_(quad_n_w_sum
     $        -quad_n_w_max),1,Z_M_svdr_(nsvd_l),n_svd_l)
      enddo !do nsvd_l=0,n_svd_l-1
      if (verbose.gt.2) then
      call print_sub_c16(n_svd_l*quad_n_w_max,Z_M_svdr_,' Z_M_svdr_: ')
      end if !if (verbose.gt.2) then
      allocate(quad_SX2_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_SX2_)
      allocate(quad_MX2_(0:1+quad_n_w_max-1))     
      call cs1_c16(quad_n_w_max,quad_MX2_)
      call zgemm('T','N',1,quad_n_w_max,n_svd_l,1.0d0
     $     *dcmplx(1.0d0,0.0d0),Z_S_svdd_
     $     ,n_svd_l,Z_S_svdr_
     $     ,n_svd_l,0.0d0 *dcmplx(1.0d0,0.0d0)
     $     ,quad_SX2_
     $     ,1)
      call zgemm('T','N',1,quad_n_w_max,n_svd_l,1.0d0
     $     *dcmplx(1.0d0,0.0d0),Z_M_svdd_
     $     ,n_svd_l,Z_M_svdr_
     $     ,n_svd_l,0.0d0 *dcmplx(1.0d0,0.0d0)
     $     ,quad_MX2_
     $     ,1)
      if (verbose.gt.1) then
      write(6,'(A,2F24.16)') ' quad_SX2_(ngamma_z_upd_0in): ' ,
     $     quad_SX2_(ngamma_z_upd_0in)
      write(6,'(A,2F24.16)') ' quad_MX2_(ngamma_z_upd_0in): ' ,
     $     quad_MX2_(ngamma_z_upd_0in)
      end if !if (verbose.gt.1) then
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_SX2_,' quad_SX2_: ')
      call print_sub_c16(quad_n_w_max,quad_MX2_,' quad_MX2_: ')
      end if !if (verbose.gt.2) then
      allocate(quad_SX1_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,quad_SX1_)
      allocate(quad_MX1_(0:1+quad_n_w_max-1))     
      call cs1_c16(quad_n_w_max,quad_MX1_)
      do nw=0,quad_n_w_max-1
         gamma_z_upd = (2.0d0*pi*nw)/(1.0d0*quad_n_w_max)
         call transf_svd_q_to_q_polyval_0on_0(svd_r_max,n_svd_r,svd_r_
     $        ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_ ,svd_s_
     $        ,svd_V_r_,quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum
     $        ,quad_P_S_k_q_,+delta_upd_x_c_0,+delta_upd_x_c_1
     $        ,quad_P_R_k_q_)
         call rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum
     $        ,quad_P_R_k_q_,+gamma_z_upd,quad_P_R_k_q_)
         quad_SX1_(nw) = dotw_c16_f(quad_n_w_sum,quad_P_R_k_q_
     $        ,quad_P_N_k_q_,quad_weight_k_p_)
         call rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum
     $        ,quad_P_S_k_q_,+gamma_z_upd,quad_P_R_k_q_)
         call transf_svd_q_to_q_polyval_0on_0(svd_r_max,n_svd_r,svd_r_
     $        ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_ ,svd_s_
     $        ,svd_V_r_,quad_n_r,quad_grid_k_p_r_,quad_n_w_,quad_n_w_sum
     $        ,quad_P_R_k_q_,+delta_upd_x_c_0,+delta_upd_x_c_1
     $        ,quad_P_R_k_q_)
         quad_MX1_(nw) = dotw_c16_f(quad_n_w_sum,quad_P_R_k_q_
     $        ,quad_P_N_k_q_,quad_weight_k_p_)
      enddo                     !do nw=0,quad_n_w_max-1
      if (verbose.gt.1) then
      write(6,'(A,2F24.16)') ' quad_SX1_(ngamma_z_upd_0in): ' ,
     $     quad_SX1_(ngamma_z_upd_0in)
      write(6,'(A,2F24.16)') ' quad_MX1_(ngamma_z_upd_0in): ' ,
     $     quad_MX1_(ngamma_z_upd_0in)
      end if !if (verbose.gt.1) then
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_max,quad_SX1_,' quad_SX1_: ')
      call print_sub_c16(quad_n_w_max,quad_MX1_,' quad_MX1_: ')
      end if !if (verbose.gt.2) then
      write(6,'(A,F8.6,A,F21.16)') ' eps_svd: ' , eps_svd ,
     $     '; quad_SX1 vs quad_SX2_: ',fnod_c16_f(quad_n_w_max
     $     ,quad_SX1 _,quad_SX2_)/al2_c16_f(quad_n_w_max
     $     ,quad_SX1_)
      write(6,'(A,F8.6,A,F21.16)') ' eps_svd: ' , eps_svd ,
     $     '; quad_MX1 vs quad_MX2_: ',fnod_c16_f(quad_n_w_max
     $     ,quad_MX1 _,quad_MX2_)/al2_c16_f(quad_n_w_max
     $     ,quad_MX1_)
      write(6,'(A,F21.16)') ' form_SX0 vs quad_SX2_: '
     $     ,fnod_c16_f(quad_n_w_max,form_SX0_,quad_SX2_)
     $     /al2_c16_f(quad_n_w_max,form_SX0_)
      write(6,'(A,F21.16)') ' form_MX0 vs quad_MX2_: '
     $     ,fnod_c16_f(quad_n_w_max,form_MX0_,quad_MX2_)
     $     /al2_c16_f(quad_n_w_max,form_MX0_)
c$$$  %%%%%%%%
      write(6,'(A)') ' Now reducing error of svd-expansion: '
      n_eps_svd = 11
      do neps_svd=0,n_eps_svd-1
      eps_svd = 0.1d0**(1.0d0 + (5.0d0*neps_svd)/(n_eps_svd-1))
      call cl1_r8(n_svd_max,svd_r_)
      call cl1_r8(n_svd_max,svd_d_)
      call cl1_i4(n_svd_max,svd_l_)
      call cl1_r8(n_svd_max*n_svd_max,svd_U_d_)
      call cl1_r8(n_svd_max,svd_s_)
      call cl1_r8(n_svd_max*n_svd_max,svd_V_r_)
      call get_svd_2(eps_svd,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,quad_grid_k_p_r_,quad_n_r,1,tmp_delta_x_c_0_
     $     ,tmp_delta_x_c_1_,flag_warning ,R_max ,K_max,delta_max
     $     ,n_pixels)
      call cl1_r8(n_svd_l*n_delta_v,svd_polyval_U_d_)
      call cl1_r8(n_svd_l*quad_n_r,svd_polyval_V_r_)
      call get_svd_polyval_U_d_(svd_d_max,n_svd_d,svd_d_,n_svd_l ,svd_l_
     $     ,svd_U_d_,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1
     $     ,svd_polyval_U_d_)
      call get_svd_polyval_V_r_(svd_r_max,n_svd_r,svd_r_,n_svd_l ,svd_l_
     $     ,svd_V_r_,quad_n_r,quad_grid_k_p_r_ ,svd_polyval_V_r_)
      call cl1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
      call cl1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
      call cl1_c16(n_svd_l*quad_n_w_max,Z_S_svdr_)
      call cl1_c16(n_svd_l*quad_n_w_max,Z_M_svdr_)
      call innerproduct_q_k_svdd_polyval_off_0(.true.,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1,Z_S_svdd_)
      call innerproduct_q_k_svdd_polyval_off_0(.false.,svd_d_max
     $     ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,svd_polyval_U_d_
     $     ,n_delta_v,delta_upd_x_c_0,delta_upd_x_c_1,Z_M_svdd_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_S_k_p_
     $     ,quad_P_S_k_q_)
      call cp1_c16(quad_n_w_sum,quad_P_M_k_p_,quad_P_N_k_p_)
      call xc1_c16(quad_n_w_sum,quad_P_N_k_p_,quad_P_C_k_p_
     $     ,quad_P_N_k_p_)
      call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum,quad_fftw_0in_
     $     ,quad_fftw_out_,quad_P_N_k_p_,-gamma_z_est,quad_P_N_k_p_)
      call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $     ,quad_n_w_sum,quad_P_N_k_p_,-delta_est_x_c_0,-delta_est_x_c_1
     $     ,quad_P_N_k_p_)
      call interp_p_to_q_fftw(quad_n_r,quad_fftw_plan_frwd_,quad_n_w_
     $     ,quad_n_w_sum,quad_fftw_0in_,quad_fftw_out_,quad_P_N_k_p_
     $     ,quad_P_N_k_q_)
      call innerproduct_q_k_svdr_quad_polyval_0on_0(.true.,svd_r_max
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,quad_n_r
     $     ,quad_grid_k_p_r_,quad_weight_k_p_r_ ,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,quad_P_N_k_q_,Z_S_svdr_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $     ,Z_S_svdr_,Z_S_svdr_)
      do nsvd_l=0,n_svd_l-1
         call cps_c16(quad_n_w_max,Z_S_svdr_(nsvd_l),n_svd_l
     $        ,quad_fftw_out_(quad_n_w_sum-quad_n_w_max),1)
         call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
         call cps_c16(quad_n_w_max,quad_fftw_0in_(quad_n_w_sum
     $        -quad_n_w_max),1,Z_S_svdr_(nsvd_l),n_svd_l)
      enddo !do nsvd_l=0,n_svd_l-1
      call innerproduct_q_k_svdr_quad_polyval_0on_0(.false.,svd_r_max
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,quad_n_r
     $     ,quad_grid_k_p_r_,quad_weight_k_p_r_ ,quad_n_w_,quad_n_w_sum
     $     ,quad_P_S_k_q_,quad_P_N_k_q_,Z_M_svdr_)
      c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $     ,Z_M_svdr_,Z_M_svdr_)
      do nsvd_l=0,n_svd_l-1
         call cps_c16(quad_n_w_max,Z_M_svdr_(nsvd_l),n_svd_l
     $        ,quad_fftw_out_(quad_n_w_sum-quad_n_w_max),1)
         call dfftw_execute_(quad_fftw_plan_back_(quad_n_r-1))
         call cps_c16(quad_n_w_max,quad_fftw_0in_(quad_n_w_sum
     $        -quad_n_w_max),1,Z_M_svdr_(nsvd_l),n_svd_l)
      enddo !do nsvd_l=0,n_svd_l-1
      call cl1_c16(quad_n_w_max,quad_SX2_)
      call cl1_c16(quad_n_w_max,quad_MX2_)
      call zgemm('T','N',1,quad_n_w_max,n_svd_l,1.0d0
     $     *dcmplx(1.0d0,0.0d0),Z_S_svdd_
     $     ,n_svd_l,Z_S_svdr_
     $     ,n_svd_l,0.0d0 *dcmplx(1.0d0,0.0d0)
     $     ,quad_SX2_
     $     ,1)
      call zgemm('T','N',1,quad_n_w_max,n_svd_l,1.0d0
     $     *dcmplx(1.0d0,0.0d0),Z_M_svdd_
     $     ,n_svd_l,Z_M_svdr_
     $     ,n_svd_l,0.0d0 *dcmplx(1.0d0,0.0d0)
     $     ,quad_MX2_
     $     ,1)
      write(6,'(A,F8.6,A,I0,A,F21.16,A,F21.16)')
     $     ' eps_svd: ' , eps_svd
     $     ,'; n_svd_l: ' , n_svd_l
     $     , '; form_SX0 vs quad_SX2_: '
     $     ,fnod_c16_f(quad_n_w_max,form_SX0_,quad_SX2_)
     $     /al2_c16_f(quad_n_w_max,form_SX0_)         
     $     , '; form_MX0 vs quad_MX2_: '
     $     ,fnod_c16_f(quad_n_w_max,form_MX0_,quad_MX2_)
     $     /al2_c16_f(quad_n_w_max,form_MX0_)         
      enddo !do neps_svd=0,n_eps_svd-1
c$$$      %%%%%%%%
      if (verbose.gt.2) then
      write(6,'(A)') ' Third: testing CTF_R_S: '
      end if !if (verbose.gt.2) then
      n_gamma_z_upd = quad_n_w_max
      n_gamma_z_nuf = 2*n_gamma_z_upd
      allocate(gamma_z_upd_(0:1+n_gamma_z_upd-1))
      call cs1_r8(n_gamma_z_upd,gamma_z_upd_)
      allocate(gamma_z_nuf_(0:1+n_gamma_z_nuf-1))
      call cs1_r8(n_gamma_z_nuf,gamma_z_nuf_)
      do ngz=0,n_gamma_z_upd-1
         gamma_z_upd_(ngz) = (2.0d0*pi*ngz)/(1.0d0*n_gamma_z_upd)
      enddo !do ngz=0,n_gamma_z_upd-1
      do ngz=0,n_gamma_z_nuf-1
         gamma_z_nuf_(ngz) = (2.0d0*pi)*((1.0d0*ngz)/(1.0d0
     $        *n_gamma_z_nuf))**2
      enddo !do ngz=0,n_gamma_z_nuf-1
      allocate(quad_Y_S_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_Y_S_k_p_)
      allocate(quad_Y_C_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_Y_C_k_p_)
      allocate(quad_Y_T_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_Y_T_k_p_)
      allocate(quad_Y_R_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_Y_R_k_p_)
      allocate(tmp_Z_p_(0:1+max(n_gamma_z_upd,n_gamma_z_nuf)-1))
      call cs1_c16(max(n_gamma_z_upd,n_gamma_z_nuf),tmp_Z_p_)
      allocate(tmp_Z_q_(0:1+max(n_gamma_z_upd,n_gamma_z_nuf)-1))
      call cs1_c16(max(n_gamma_z_upd,n_gamma_z_nuf),tmp_Z_q_)
      allocate(form_CTF_R_S_0__(0:1+n_gamma_z_upd-1))
      call cs1_c16(n_gamma_z_upd,form_CTF_R_S_0__)
      allocate(quad_CTF_R_S_0__(0:1+n_gamma_z_upd-1))
      call cs1_c16(n_gamma_z_upd,quad_CTF_R_S_0__)
      allocate(quad_Y_CTF_R_S_build_uniform__(0:1+n_gamma_z_upd-1))
      call cs1_c16(n_gamma_z_upd,quad_Y_CTF_R_S_build_uniform__)
      allocate(quad_Y_CTF_R_S_brute_uniform__(0:1+n_gamma_z_upd-1))
      call cs1_c16(n_gamma_z_upd,quad_Y_CTF_R_S_brute_uniform__)
      allocate(quad_Y_CTF_R_S_build_nonunif__(0:1+n_gamma_z_nuf-1))
      call cs1_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_build_nonunif__)
      allocate(quad_Y_CTF_R_S_brute_nonunif__(0:1+n_gamma_z_nuf-1))
      call cs1_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_brute_nonunif__)
c$$$      %%%%%%%%
      call ti8_build_CTF_R_S_quad_3(
     $     0
     $     ,n_gamma_z_upd
     $     ,gamma_z_upd_
     $     ,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_
     $     ,quad_fftw_0in_
     $     ,quad_fftw_out_
     $     ,quad_n_r
     $     ,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_
     $     ,quad_n_w_
     $     ,quad_n_w_sum
     $     ,quad_P_R_k_p_
     $     ,quad_P_R_k_q_
     $     ,1
     $     ,0
     $     ,quad_n_w_sum 
     $     ,quad_P_S_k_p_
     $     ,quad_P_T_k_p_
     $     ,quad_P_T_k_q_
     $     ,1
     $     ,quad_n_w_sum 
     $     ,quad_P_C_k_p_
     $     ,quad_CTF_R_S_0__
     $     ,tmp_Z_q_
     $     ,tmp_Z_p_
     $     )
      c16_alpha = dcmplx(1.0d0/dsqrt(2.0d0)
     $     ,0.0d0)
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_gamma_z_upd,c16_alpha,c16_beta
     $     ,quad_CTF_R_S_0__,quad_CTF_R_S_0__)
      do ngz=0,n_gamma_z_upd-1
         gamma_z_upd = gamma_z_upd_(ngz)
         tmp_phi_n = fix_phi_S + gamma_z_upd - fix_phi_M
         form_CTF_R_S_0__(ngz) = dsqrt(2*pi*half_diameter_k_c**6/3.0d0 *
     $        (1.0d0 + 2.0d0*dcos(tmp_phi_n)*dcos(tmp_phi_n)))
      enddo !do ngz=0,n_gamma_z_upd-1
      if (verbose.gt.2) then
         write(6,'(3(A,F16.8))')
     $   ' half_diameter_k_c: ' , half_diameter_k_c
     $   ,' fix_phi_M: ' , fix_phi_M
     $   ,' fix_phi_S: ' , fix_phi_S
         call print_sub_r8(n_gamma_z_upd,gamma_z_upd_,
     $        ' gamma_z_upd_: ')
         call print_sub_c16(quad_n_w_sum,quad_P_S_k_p_,
     $        ' quad_P_S_k_p_: ')
         call print_sub_c16(quad_n_w_sum,quad_P_C_k_p_,
     $        ' quad_P_C_k_p_: ')
         call print_sub_r8(quad_n_r,quad_weight_k_p_r_,
     $        ' quad_weight_k_p_r_: ')
         call print_sub_c16(n_gamma_z_upd,quad_CTF_R_S_0__,
     $        ' quad_CTF_R_S_0__: ')
         call print_sub_c16(n_gamma_z_upd,form_CTF_R_S_0__,
     $        ' form_CTF_R_S_0__: ')
      end if !if (verbose.gt.2) then
      if (verbose.gt.2) then
      write(6,'(A,F21.16)') ' quad_CTF_R_S_0__ vs form_CTF_R_S_0__: '
     $     ,fnod_c16_f(n_gamma_z_upd,quad_CTF_R_S_0__,form_CTF_R_S_0__)
     $     /al2_c16_f(n_gamma_z_upd,form_CTF_R_S_0__)
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      na=0
      do nr=0,quad_n_r-1
         tmp_k = quad_grid_k_p_r_(nr)
         do nw=0,quad_n_w_(nr)-1
            tmp_w = (2.0d0*pi*nw)/quad_n_w_(nr)
            quad_Y_S_k_p_(na) = tmp_k*dexp(-tmp_k**2/(2.0d0
     $           *half_diameter_k_c**2))*(dcos(tmp_w-pi/3.0d0) - 0.5d0
     $           *dcos(3.0d0 *(tmp_w+pi/2.0d0)) + 0.25d0*dcos(7.0d0
     $           *(tmp_w-pi/4.0d0)) - 0.125d0 *dcos(13.0d0*(tmp_w+pi
     $           /6.0d0)))
            quad_Y_C_k_p_(na) = tmp_k*dexp(-tmp_k**2/(2.0d0
     $           *half_diameter_k_c**2))*(dcos(2.0d0*(tmp_w-0.0d0)) -
     $           0.25d0*cos(5.0d0 *(tmp_w+pi/12.0d0)) - 0.35d0
     $           *dcos(4.0d0*(tmp_w-pi/2.0d0)) + 0.125d0*dcos(9.0d0
     $           *(tmp_w+pi/6.0d0)) + 0.125d0*dcos(17.0d0*(tmp_w-pi
     $           /5.0d0)))
            na = na+1
         enddo !do nw=0,quad_n_w_(nr)-1
      enddo !do nr=0,quad_n_r-1
      if (verbose.gt.2) then
      call print_sub_c16(quad_n_w_sum,quad_Y_S_k_p_,' quad_Y_S_k_p_: ')
      call print_sub_c16(quad_n_w_sum,quad_Y_C_k_p_,' quad_Y_C_k_p_: ')
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      call ti8_build_CTF_R_S_quad_3(
     $     0
     $     ,n_gamma_z_upd
     $     ,gamma_z_upd_
     $     ,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_
     $     ,quad_fftw_0in_
     $     ,quad_fftw_out_
     $     ,quad_n_r
     $     ,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_
     $     ,quad_n_w_
     $     ,quad_n_w_sum
     $     ,quad_P_R_k_p_ 
     $     ,quad_P_R_k_q_
     $     ,1
     $     ,0
     $     ,quad_n_w_sum 
     $     ,quad_Y_S_k_p_
     $     ,quad_P_T_k_p_ 
     $     ,quad_P_T_k_q_
     $     ,1
     $     ,quad_n_w_sum 
     $     ,quad_Y_C_k_p_
     $     ,quad_Y_CTF_R_S_build_uniform__
     $     ,tmp_Z_q_
     $     ,tmp_Z_p_
     $     )
      c16_alpha = dcmplx(1.0d0/dsqrt(2.0d0)
     $     ,0.0d0)
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_gamma_z_upd,c16_alpha,c16_beta
     $     ,quad_Y_CTF_R_S_build_uniform__
     $     ,quad_Y_CTF_R_S_build_uniform__)
      if (verbose.gt.2) then
         call print_sub_c16(n_gamma_z_upd,quad_Y_CTF_R_S_build_uniform__
     $        ,' quad_Y_CTF_R_S_build_uniform__: ')
      end if !if (verbose.gt.2) then
      do ngz=0,n_gamma_z_upd-1
         gamma_z_upd = gamma_z_upd_(ngz)
         call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $        ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum
     $        ,quad_fftw_0in_ ,quad_fftw_out_,quad_Y_S_k_p_,+gamma_z_upd
     $        ,quad_Y_T_k_p_)
         call xc1_c16(quad_n_w_sum,quad_Y_T_k_p_,quad_Y_C_k_p_
     $        ,quad_Y_R_k_p_)
         quad_Y_CTF_R_S_brute_uniform__(ngz) =
     $        zsqrt(dotw_c16_f(quad_n_w_sum,quad_Y_R_k_p_,quad_Y_R_k_p_
     $        ,quad_weight_k_p_))
      enddo !do ngz=0,n_gamma_z_upd-1
      if (verbose.gt.2) then
         call print_sub_c16(n_gamma_z_upd,quad_Y_CTF_R_S_brute_uniform__
     $        ,' quad_Y_CTF_R_S_brute_uniform__: ')
      end if !if (verbose.gt.2) then
      if (verbose.gt.2) then
      write(6,'(A,A,F21.16)')
     $     ' fftw: quad_Y_CTF_R_S_brute_uniform__'
     $     ,' vs quad_Y_CTF_R_S_build_uniform__: '
     $     ,fnod_c16_f(n_gamma_z_upd,quad_Y_CTF_R_S_brute_uniform__
     $     ,quad_Y_CTF_R_S_build_uniform__) /al2_c16_f(n_gamma_z_upd
     $     ,quad_Y_CTF_R_S_build_uniform__)
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%
      call ti8_build_CTF_R_S_quad_3(
     $     0
     $     ,n_gamma_z_nuf
     $     ,gamma_z_nuf_
     $     ,quad_fftw_plan_frwd_
     $     ,quad_fftw_plan_back_
     $     ,quad_fftw_0in_
     $     ,quad_fftw_out_
     $     ,quad_n_r
     $     ,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_
     $     ,quad_n_w_
     $     ,quad_n_w_sum
     $     ,quad_P_R_k_p_ 
     $     ,quad_P_R_k_q_
     $     ,1
     $     ,0
     $     ,quad_n_w_sum 
     $     ,quad_Y_S_k_p_
     $     ,quad_P_T_k_p_ 
     $     ,quad_P_T_k_q_
     $     ,1
     $     ,quad_n_w_sum 
     $     ,quad_Y_C_k_p_
     $     ,quad_Y_CTF_R_S_build_nonunif__
     $     ,tmp_Z_q_
     $     ,tmp_Z_p_
     $     )
      c16_alpha = dcmplx(1.0d0/dsqrt(2.0d0)
     $     ,0.0d0)
      c16_beta = dcmplx(0.0d0,0.0d0); 
      call af1_c16(n_gamma_z_nuf,c16_alpha,c16_beta
     $     ,quad_Y_CTF_R_S_build_nonunif__
     $     ,quad_Y_CTF_R_S_build_nonunif__)
      if (verbose.gt.2) then
         call print_sub_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_build_nonunif__
     $        ,' quad_Y_CTF_R_S_build_nonunif__: ')
      end if !if (verbose.gt.2) then
      do ngz=0,n_gamma_z_nuf-1
         gamma_z_upd = gamma_z_nuf_(ngz)
         call rotate_p_to_p_fftw(quad_n_r,quad_fftw_plan_frwd_
     $        ,quad_fftw_plan_back_,quad_n_w_,quad_n_w_sum
     $        ,quad_fftw_0in_ ,quad_fftw_out_,quad_Y_S_k_p_,+gamma_z_upd
     $        ,quad_Y_T_k_p_)
         call xc1_c16(quad_n_w_sum,quad_Y_T_k_p_,quad_Y_C_k_p_
     $        ,quad_Y_R_k_p_)
         quad_Y_CTF_R_S_brute_nonunif__(ngz) =
     $        zsqrt(dotw_c16_f(quad_n_w_sum,quad_Y_R_k_p_,quad_Y_R_k_p_
     $        ,quad_weight_k_p_))
      enddo !do ngz=0,n_gamma_z_nuf-1
      if (verbose.gt.2) then
         call print_sub_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_brute_nonunif__
     $        ,' quad_Y_CTF_R_S_brute_nonunif__: ')
      end if !if (verbose.gt.2) then
      if (verbose.gt.2) then
      write(6,'(A,A,F21.16)')
     $     ' finufft: quad_Y_CTF_R_S_brute_nonunif__'
     $     ,' vs quad_Y_CTF_R_S_build_nonunif__: '
     $     ,fnod_c16_f(n_gamma_z_nuf,quad_Y_CTF_R_S_brute_nonunif__
     $     ,quad_Y_CTF_R_S_build_nonunif__) /al2_c16_f(n_gamma_z_nuf
     $     ,quad_Y_CTF_R_S_build_nonunif__)
      end if !if (verbose.gt.2) then
c$$$      %%%%%%%%

c$$$      %%%%%%%%
      flag_memory_checkset = .true.
      call cxs_r8(quad_n_r,jacpts_x_,'jacpts_x_',flag_memory_checkset)
      call cxs_r8(quad_n_r,jacpts_w_,'jacpts_w_',flag_memory_checkset)
      call cxs_r8(quad_n_r,quad_grid_k_p_r_,'quad_grid_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_r,quad_weight_k_p_r_,'quad_weight_k_p_r_'
     $     ,flag_memory_checkset)
      call cxs_i4(quad_n_r,quad_n_polar_a_,'quad_n_polar_a_'
     $     ,flag_memory_checkset)
      call cxs_i4(quad_n_r,quad_n_w_,'quad_n_w_',flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_grid_k_p_0_,'quad_grid_k_p_0_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_grid_k_p_1_,'quad_grid_k_p_1_'
     $     ,flag_memory_checkset)
      call cxs_r8(quad_n_w_sum,quad_weight_k_p_,'quad_weight_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,fix_delta_x_c_,'fix_delta_x_c_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_k_p_,'quad_P_k_p_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,fix_delta_S_x_c_,'fix_delta_S_x_c_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_S_k_p_,'quad_P_S_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_S_k_q_,'quad_P_S_k_q_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,fix_delta_M_x_c_,'fix_delta_M_x_c_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_M_k_p_,'quad_P_M_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_M_k_q_,'quad_P_M_k_q_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_C_k_p_,'quad_P_C_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_N_k_p_,'quad_P_N_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_N_k_q_,'quad_P_N_k_q_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_T_k_p_,'quad_P_T_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_T_k_q_,'quad_P_T_k_q_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_R_k_p_,'quad_P_R_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_P_R_k_q_,'quad_P_R_k_q_'
     $     ,flag_memory_checkset)
      call cxs_i8(quad_n_r,quad_fftw_plan_frwd_,',quad_fftw_plan_frwd_'
     $     ,flag_memory_checkset)
      call cxs_i8(quad_n_r,quad_fftw_plan_back_,',quad_fftw_plan_back_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_fftw_0in_,',quad_fftw_0in_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_fftw_out_,',quad_fftw_out_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_X0_,',quad_X0_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,tmp_delta_S_x_c_,'tmp_delta_S_x_c_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,tmp_delta_M_x_c_,'tmp_delta_M_x_c_'
     $     ,flag_memory_checkset)
      call cxs_r8(2,tmp_delta_Y_x_c_,'tmp_delta_Y_x_c_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,form_X0_,',form_X0_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,form_SX0_,',form_SX0_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,form_MX0_,',form_MX0_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_X1_,',quad_X1_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_svd_max,svd_r_,'svd_r_',flag_memory_checkset)
      call cxs_r8(n_svd_max,svd_d_,'svd_d_',flag_memory_checkset)
      call cxs_i4(n_svd_max,svd_l_,'svd_l_',flag_memory_checkset)
      call cxs_r8(n_svd_max*n_svd_max,svd_U_d_,'svd_U_d_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_svd_max,svd_s_,'svd_s_',flag_memory_checkset)
      call cxs_r8(n_svd_max*n_svd_max,svd_V_r_,'svd_V_r_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_svd_l*n_delta_v,svd_polyval_U_d_,'svd_polyval_U_d_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_svd_l*quad_n_r,svd_polyval_V_r_,'svd_polyval_V_r_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_X2_,',quad_X2_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_X3_,',quad_X3_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*n_delta_v,Z_S_svdd_,'Z_S_svdd_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*n_delta_v,Z_M_svdd_,'Z_M_svdd_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,'Z_S_svdr_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*quad_n_w_max,Z_M_svdr_,'Z_M_svdr_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_SX2_,'quad_SX2_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_MX2_,'quad_MX2_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_SX1_,'quad_SX1_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,quad_MX1_,'quad_MX1_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_gamma_z_upd,gamma_z_upd_,'gamma_z_upd_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_gamma_z_nuf,gamma_z_nuf_,'gamma_z_nuf_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_Y_S_k_p_,'quad_Y_S_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_Y_C_k_p_,'quad_Y_C_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_Y_T_k_p_,'quad_Y_T_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_Y_R_k_p_,'quad_Y_R_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(max(n_gamma_z_upd,n_gamma_z_nuf),tmp_Z_p_,'tmp_Z_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(max(n_gamma_z_upd,n_gamma_z_nuf),tmp_Z_q_,'tmp_Z_q_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z_upd,form_CTF_R_S_0__,'form_CTF_R_S_0__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z_upd,quad_CTF_R_S_0__,'quad_CTF_R_S_0__'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z_upd,quad_Y_CTF_R_S_build_uniform__
     $     ,'quad_Y_CTF_R_S_build_uniform__',flag_memory_checkset)
      call cxs_c16(n_gamma_z_upd,quad_Y_CTF_R_S_brute_uniform__
     $     ,'quad_Y_CTF_R_S_brute_uniform__',flag_memory_checkset)
      call cxs_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_build_nonunif__
     $     ,'quad_Y_CTF_R_S_build_nonunif__',flag_memory_checkset)
      call cxs_c16(n_gamma_z_nuf,quad_Y_CTF_R_S_brute_nonunif__
     $     ,'quad_Y_CTF_R_S_brute_nonunif__',flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.2) then
         write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.2) then
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

 10   continue
      if (verbose.gt.0) then
         write(6,'(A)') '[finished ti8_dr_digits_single]'
      end if !if (verbose.gt.0) then
      end !driver. ;

