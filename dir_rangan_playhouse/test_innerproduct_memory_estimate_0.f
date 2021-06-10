      subroutine test_innerproduct_memory_estimate_0(verbose ,n_omp_sub
     $     ,n_k_p_max,n_M,n_CTF,n_alpha ,n_pixels_in
     $     ,displacement_max,n_delta_x ,n_delta_y,n_gamma_z
     $     ,svd_calculation_type ,eps_svd ,fpm_howmany_max,n_SM)
      implicit none
      integer verbose
      integer n_omp_sub,n_k_p_max,n_M,n_CTF,n_alpha
     $     ,n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type,n_SM
      integer n_r,nr,n_A,n_w_max,n_S
      integer fpm_howmany_max,nM_per_max
      real *8 fpm_howmany
      real *8 eps_svd,eps_target,n_pixels_in,displacement_max
      real *8 half_diameter_x_c
      real *8, allocatable :: grid_p_(:)
      real *8, allocatable :: weight_p_(:)
      integer, allocatable :: n_w_(:)
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      integer n_svd_max
      parameter (n_svd_max=512)
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      real *8 pi
      integer n_transf
      integer *4, allocatable :: n_Y_lm_sum_max_(:)
      integer *4, allocatable :: ld_S_(:)
      integer *4, allocatable :: n_polar_a_azimu_b_sum_(:)
      integer *4, allocatable :: n_Y_l_(:)
      integer *4, allocatable :: n_Y_lm_sum_(:)
      integer *4, allocatable :: n_polar_a_(:)
      integer *4 nk,nl,n_Y_lm_sum_max,ld_S,ld_M
      integer *4 n_polar_a,npolar_a
      integer *4 n_azimu_b,n_polar_a_azimu_b_sum
      real *8 dpolar_a
      real *8 phi_over
      parameter(phi_over=2.0d0)
      real *8, allocatable :: polar_a_(:)
      real *8, allocatable :: sin_polar_a_(:)
      real *8 mem_S_p_k__
      real *8 mem_M_p_k__
      real *8 mem_fpm_in1__
      real *8 mem_Z_svdd_
      real *8 mem_CTF_R_S_
      real *8 mem_S_q__
      real *8 mem_T_S_q__
      real *8 mem_Z_S_q__
      real *8 mem_T_R_CTF_M_q__
      real *8 mem_T_T_R_CTF_M_q__
      real *8 mem_Z_T_R_CTF_M_q__
      real *8 mem_S_T_T_R_CTF_M_q__
      real *8 mem_S_Z_T_R_CTF_M_q__
      real *8 mem_ZZ_sub__
      real *8 mem_ZZ_all__
      real *8 mem_output_sub__
      real *8 mem_output_all__
      real *8 mem_per_GB
      parameter (mem_per_GB=1.0737d9)
      real *8 mem_r8
      parameter (mem_r8=8.0d0)
      real *8 mem_c16
      parameter (mem_c16=16.0d0)
            
      if (verbose.gt.1) then
         write(6,'(A)')
     $        '[entering test_innerproduct_memory_estimate_0]'
      end if !if (verbose.gt.1) then

      n_r = n_k_p_max
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r',n_r,'<2'
      end if
      pi = 4.0d0*atan(1.0d0)
      eps_target = eps_svd

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_omp_sub: ' , n_omp_sub
         write(6,'(A,I0)') ' n_k_p_max: ' , n_k_p_max
         write(6,'(A,I0)') ' n_M: ' , n_M
         write(6,'(A,I0)') ' n_CTF: ' , n_CTF
         write(6,'(A,I0)') ' n_alpha: ' , n_alpha
         write(6,'(A,F8.4)') ' n_pixels_in: ' , n_pixels_in
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
         write(6,'(A,I0)') ' n_delta_x: ' , n_delta_x
         write(6,'(A,I0)') ' n_delta_y: ' , n_delta_y
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         write(6,'(A,I0)') ' svd_calculation_type: ' ,
     $        svd_calculation_type
         write(6,'(A,F8.4)') ' eps_svd: ' , eps_svd
         write(6,'(A,I0)') ' fpm_howmany_max: ' , fpm_howmany_max
         write(6,'(A,I0)') ' n_SM: ' , n_SM
      end if !if (verbose.gt.1) then
      
c$$$      Calculating grid_p_ <-- grid_k_p_ using 'getgridr'
      allocate(grid_p_(0:n_r-1))
      allocate(weight_p_(0:n_r-1))
      call getgridr(1.0d0*n_r,n_r,0,grid_p_,weight_p_)
      if (verbose.gt.0) then
         call write_sub_r8(n_r,grid_p_,10,' grid_p_: ')
      end if !if (verbose.gt.0) then


      if (verbose.gt.1) then
         write(6,'(A)') ' determining array sizes'
      end if !if (verbose.gt.1) then

      allocate(n_Y_lm_sum_max_(0:n_r-1))
      allocate(ld_S_(0:n_r-1))
      allocate(n_Y_l_(0:n_r-1))
      allocate(n_Y_lm_sum_(0:n_r-1))
      allocate(n_polar_a_azimu_b_sum_(0:n_r-1))
      allocate(n_polar_a_(0:n_r-1))
      allocate(n_w_(0:n_r-1))
      allocate(polar_a_(0:12*max(1,n_r)-1))
      allocate(sin_polar_a_(0:12*max(1,n_r)-1))

      n_Y_lm_sum_max = 0
      ld_S = 0
      do nr=0,n_r-1
         n_polar_a = 2*ceiling(max(6,nint(pi*grid_p_(nr)))
     $        /2.0d0)
         n_polar_a_(nr) = n_polar_a
         n_w_(nr) = n_polar_a*2
         dpolar_a = 1.0d0/(2.0d0*n_polar_a)
         do nl=0,n_polar_a-1
            polar_a_(nl) = (dpolar_a + 2.0d0*dpolar_a*nl)*pi
            sin_polar_a_(nl) = sin(polar_a_(nl))
         enddo !do nl=0,n_polar_a-1
         n_polar_a_azimu_b_sum = 0
         do npolar_a=0,n_polar_a-1
            n_azimu_b = max(6,nint(n_polar_a*phi_over
     $           *sin_polar_a_(npolar_a)))
            n_polar_a_azimu_b_sum = n_polar_a_azimu_b_sum + n_azimu_b
         enddo !do npolar_a=0,n_polar_a-1
         n_polar_a_azimu_b_sum_(nr) = n_polar_a_azimu_b_sum
         if (n_w_(nr).ne.n_polar_a_(nr)*2) write(6,'(A)')
     $        'Warning, n_w_(nr) .ne. n_polar_a_(nr)'
         n_Y_lm_sum_(nr) = n_Y_lm_sum_max
         n_Y_l_(nr) = nint(grid_p_(nr)+2)
         n_Y_lm_sum_max = n_Y_lm_sum_max + (n_Y_l_(nr)+1)**2
         ld_S = ld_S + n_w_(nr)
         n_Y_lm_sum_max_(nr) = n_Y_lm_sum_max
         ld_S_(nr) = ld_S
      enddo !do nr=0,n_r-1

      call write_sub_i4(n_r,ld_S_,7,' ld_S: ')
      call write_sub_i4(n_r,n_Y_lm_sum_max_,18,' n_Y_lm_sum_max_: ')
      call write_sub_i4(n_r,n_polar_a_azimu_b_sum_,25
     $     ,' n_polar_a_azimu_b_sum_: ')

      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') ' n_w_max ',n_w_max,'; n_A ',n_A
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of displacements to measure'
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      half_diameter_x_c = 1.0d0
      call get_delta_0(N_pixels_in,n_r,half_diameter_x_c,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      call get_gamma_0(n_gamma_z,gamma_z_)
      if (verbose.gt.1) then
         call write_sub_r8(n_delta_x,delta_x_,11,' delta_x_: ')
         call write_sub_r8(n_delta_y,delta_y_,11,' delta_y_: ')
         call write_sub_r8(n_gamma_z,gamma_z_,11,' gamma_z_: ')
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if !if (verbose.gt.1) then
      allocate(svd_r_(0:n_svd_max-1))
      allocate(svd_d_(0:n_svd_max-1))
      allocate(svd_l_(0:n_svd_max-1))
      allocate(svd_U_d_(0:n_svd_max*n_svd_max-1))
      allocate(svd_s_(0:n_svd_max-1))
      allocate(svd_V_r_(0:n_svd_max*n_svd_max-1))
      call get_svd_0(eps_target,n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $     ,svd_l_,svd_U_d_ ,svd_s_,svd_V_r_,svd_unitnumber,svd_fname
     $     ,grid_p_,n_r,n_delta_x,delta_x_,n_delta_y
     $     ,delta_y_,flag_warning,R_max ,K_max,delta_max,n_pixels)
      if (verbose.gt.2) then
         write(6,'(A,I0)') ' n_svd_r: ',n_svd_r
         write(6,'(A,I0)') ' n_svd_d: ',n_svd_d
         write(6,'(A,I0)') ' svd_unitnumber: ',svd_unitnumber
      end if !if (verbose.gt.1) then
      if (verbose.gt.0) then
         write(6,'(A,I0,A,A,I0)') ' n_svd_l: ' , n_svd_l , ' versus '
     $        ,' n_delta_x*n_delta_y: ' , n_delta_x*n_delta_y
      end if !if (verbose.gt.0) then            

      n_transf = min(n_svd_l,n_delta_x*n_delta_y)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_transf: ' , n_transf
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') ' memory requirements: '
      end if !if (verbose.gt.1) then

      n_S = n_polar_a_azimu_b_sum_(n_r-1)
      ld_S = ld_S_(n_r-1)
      ld_M = ld_S_(n_r-1)
      nM_per_max = ceiling((1.0d0*n_M)/(1.0d0*max(1,n_omp_sub)))
      fpm_howmany = min(1.0d0*fpm_howmany_max,1.0d0*n_transf*n_S
     $     *nM_per_max)
      mem_S_p_k__ = (mem_c16*ld_S*n_S)/mem_per_GB
      mem_M_p_k__ = (mem_c16*ld_M*n_M)/mem_per_GB
      mem_fpm_in1__ = (mem_c16*fpm_howmany*n_w_max*n_omp_sub)/mem_per_GB
      mem_Z_svdd_ = (mem_c16*n_svd_l*n_delta_x*n_delta_y)/mem_per_GB
      mem_CTF_R_S_ = (mem_c16*n_gamma_z*n_S*n_CTF)/mem_per_GB
      mem_T_S_q__ = (mem_c16*n_r*n_delta_x*n_delta_y*n_S*n_w_max)
     $     /mem_per_GB
      mem_S_q__ = (mem_c16*n_r*n_S*n_w_max)/mem_per_GB
      mem_Z_S_q__ = (mem_c16*n_r*n_svd_l*n_S*n_w_max)/mem_per_GB
      mem_T_R_CTF_M_q__ = (mem_c16*n_r*n_M*n_w_max)/mem_per_GB
      mem_T_T_R_CTF_M_q__ = (mem_c16*n_r*n_delta_x*n_delta_y*n_M
     $     *n_w_max)/mem_per_GB
      mem_Z_T_R_CTF_M_q__ = (mem_c16*n_r*n_svd_l*n_M*n_w_max)/mem_per_GB
      mem_S_T_T_R_CTF_M_q__ = (mem_c16*n_delta_x*n_delta_y*n_S*n_M
     $     *n_w_max)/mem_per_GB
      mem_S_Z_T_R_CTF_M_q__ = (mem_c16*n_svd_l*n_S*n_M*n_w_max)
     $     /mem_per_GB
      mem_ZZ_sub__ = (mem_c16*n_delta_x*n_delta_y*n_gamma_z*n_omp_sub)
     $     /mem_per_GB
      mem_ZZ_all__ = (mem_c16*n_delta_x*n_delta_y*n_S*n_M)/mem_per_GB
      mem_output_sub__ = (mem_r8*n_alpha*n_SM*n_M)/mem_per_GB
      mem_output_all__ = (mem_r8*n_alpha*n_S*n_M)/mem_per_GB

      write(6,'(A,F16.2,A)') ' mem_S_p_k__: ' , mem_S_p_k__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_M_p_k__: ' , mem_M_p_k__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_fpm_in1__: ' , mem_fpm_in1__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_Z_svdd_: ' , mem_Z_svdd_ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_CTF_R_S_: ' , mem_CTF_R_S_ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_S_q__: ' , mem_S_q__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_T_S_q__: ' , mem_T_S_q__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_Z_S_q__: ' , mem_Z_S_q__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_T_R_CTF_M_q__: ' , mem_T_R_CTF_M_q__
     $     , ' GB '
      write(6,'(A,F16.2,A)') ' mem_T_T_R_CTF_M_q__: ' ,
     $     mem_T_T_R_CTF_M_q__, ' GB '
      write(6,'(A,F16.2,A)') ' mem_Z_T_R_CTF_M_q__: ' ,
     $     mem_Z_T_R_CTF_M_q__, ' GB '
      write(6,'(A,F16.2,A)') ' mem_S_T_T_R_CTF_M_q__: ' ,
     $     mem_S_T_T_R_CTF_M_q__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_S_Z_T_R_CTF_M_q__: ' ,
     $     mem_S_Z_T_R_CTF_M_q__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_ZZ_sub__: ' , mem_ZZ_sub__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_ZZ_all__: ' , mem_ZZ_all__ , ' GB '
      write(6,'(A,F16.2,A)') ' mem_output_sub__: ' , mem_output_sub__ ,
     $     ' GB '
      write(6,'(A,F16.2,A)') ' mem_output_all__: ' , mem_output_all__ ,
     $     ' GB '

      if (verbose.gt.1) then
         write(6,'(A)')
     $        '[finished test_innerproduct_memory_estimate_0]'
      end if !if (verbose.gt.1) then
      
      end
