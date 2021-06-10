!> Doxygen comment: ;\n
!> This program calculates innerproducts via brute force. ;\n
      subroutine test_innerproduct_bruteforce_1(verbose ,n_r ,grid_p_
     $     ,n_w_ ,half_diameter_x_c ,S_p_0in_,M_p_0in_ ,CTF_p_0in_
     $     ,delta_x_est ,delta_y_est ,gamma_z_est ,n_delta_v ,delta_x_
     $     ,delta_y_ ,n_gamma_z ,gamma_z_)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer *4 n_r,n_w_(0:n_r-1)
      real *8 grid_p_(0:n_r-1),half_diameter_x_c
      complex *16 S_p_0in_(0:0)
      complex *16 M_p_0in_(0:0)
      complex *16 CTF_p_0in_(0:0)
c$$$      array of displacements and rotations to measure
      real *8 delta_x_est,delta_y_est,gamma_z_est
      integer n_delta_v,n_gamma_z
      real *8 delta_x_(0:0)
      real *8 delta_y_(0:0)
      real *8 gamma_z_(0:0)
c$$$      fftw for local use
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_0in_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_0in_last_,fftw_0in_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_0in_last_(*),fftw_out_last_(*)
c$$$      temporary arrays for templates and images
      complex *16, allocatable :: S_p_(:)
      complex *16, allocatable :: M_p_(:)
      complex *16, allocatable :: CTF_p_(:)
      complex *16, allocatable :: T_p_(:)
      complex *16, allocatable :: S_q_(:)
      complex *16, allocatable :: M_q_(:)
      complex *16, allocatable :: CTF_q_(:)
      complex *16, allocatable :: T_q_(:)
c$$$      temporary arrays for innerproducts
      integer *4, allocatable :: I_S_sample_(:)
      integer *4, allocatable :: I_M_sample_(:)
      complex *16, allocatable :: CTF_R_S_xx_(:)
      complex *16, allocatable :: CTF_R_S_bf_(:)
      complex *16, allocatable :: S_T_R_T_R_M_xx_(:)
      complex *16, allocatable :: S_T_R_T_R_M_bf_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16 C_p,C_q
      real *8 costheta_c16_f
c$$$      indices
      integer n_S,n_M,n_CTF,ld_S,ld_M,ld_CTF,nr,n_w_max,n_A,na,nC
      integer ndv,ngz,ndv_optimal,ngz_optimal,n_X
      real *8 delta_x,delta_y,gamma_z
      character(len=64) format_string
c$$$      pi
      real *8 pi
c$$$      parameters for timing
      real *8 timing_tic,timing_toc,timing_tmp,timing_tot

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_bruteforce_1]'
      end if      

      if (verbose.gt.1) then
         write(6,'(A)') ' Generating indices'
         write(6,'(A)') ' '
      end if
      pi = 4.0d0*datan(1.0d0)
      n_A = 0
      do nr=0,n_r-1
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_A: ',n_A,'; n_w_max: ',n_w_max
         call print_all_i4(n_r,n_w_,' n_w_: ')
      end if
      ld_S = n_A
      ld_M = n_A
      ld_CTF = n_A
      n_S = 1
      n_M = 1
      n_CTF = 1

      if (verbose.gt.1) then
         write(6,'(A,F8.4)') ' delta_x_est: ' , delta_x_est
         write(6,'(A,F8.4)') ' delta_y_est: ' , delta_y_est
         write(6,'(A,F8.4)') ' gamma_z_est: ' , gamma_z_est
         call print_sub_c16(n_A,S_p_0in_,' S_p_0in_: ')
         call print_sub_c16(n_A,M_p_0in_,' M_p_0in_: ')
         call print_sub_c16(n_A,CTF_p_0in_,' CTF_p_0in_: ')
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for local use'
         write(6,'(A)') ' '
      end if
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_0in_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      timing_tic = omp_get_wtime()
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_0in_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_0in_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nr)
      enddo
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_1d:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_0in_last_ = loc(fftw_0in_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

      allocate(S_p_(0:n_A-1))
      allocate(M_p_(0:n_A-1))
      allocate(CTF_p_(0:n_A-1))
      allocate(T_p_(0:n_A-1))
      allocate(S_q_(0:n_A-1))
      allocate(M_q_(0:n_A-1))
      allocate(CTF_q_(0:n_A-1))
      allocate(T_q_(0:n_A-1))
      allocate(I_S_sample_(0:0))
      I_S_sample_(0)=0
      allocate(I_M_sample_(0:0))
      I_M_sample_(0)=0

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array CTF_R_S_xx_ to hold'
         write(6,'(A)') ' innerproducts for each CTF-S pair.'
         write(6,'(A)') ' '
      end if
      allocate(CTF_R_S_xx_(0:n_gamma_z))
      call cl1_c16(n_gamma_z,CTF_R_S_xx_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate CTF_R_S_xx_.'
         write(6,'(A)') ' CTF_R_S = || CTF .* R(S) ||.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)') ' CTF_R_S_xx_(ngz) '
         write(6,'(A)') ' is equal to the l2_norm (not squared)'
         write(6,'(A)') ' of the pointwise product of CTF and R(S),'
         write(6,'(A)') ' where CTF is the ctf-function, '
         write(6,'(A)') ' R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' and S is the template, '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      call test_innerproduct_bruteforce_CTF_R_S_1(verbose-1,n_gamma_z
     $     ,gamma_z_,fftw_plan_frwd_,fftw_0in_,fftw_out_
     $     ,fftw_plan_back_last,fftw_0in_last_,fftw_out_last_,n_r
     $     ,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,I_S_sample_,ld_S ,S_p_0in_
     $     ,CTF_p_ ,CTF_q_ ,n_CTF,ld_CTF ,CTF_p_0in_ ,CTF_R_S_xx_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating CTF_R_S_xx_ for each ctf-S pair. '
     $        ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)*(n_CTF*1.0d0)/timing_tot
     $        /1e9
         write(6,'(A,F8.4)') ' CTF_R_S_xx_ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array CTF_R_S_bf_ to hold'
         write(6,'(A)') ' innerproducts for each CTF-S pair.'
         write(6,'(A)') ' '
      end if
      allocate(CTF_R_S_bf_(0:n_gamma_z))
      call cl1_c16(n_gamma_z,CTF_R_S_bf_)
      timing_tic = omp_get_wtime()
      do ngz=0,n_gamma_z
         gamma_z = gamma_z_(ngz)
         call rotate_p_to_p_fftw(n_r,fftw_plan_frwd_,fftw_plan_back_
     $        ,n_w_,n_A,fftw_0in_,fftw_out_,S_p_0in_,+gamma_z,S_p_)
         call xx1_c16(n_A,S_p_,CTF_p_0in_,S_p_)
         call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_p_,S_p_,C_p)
         CTF_R_S_bf_(ngz) = zsqrt(C_p)
      enddo !do ngz=0,n_gamma_z
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating CTF_R_S_bf_ for each ctf-S pair. '
     $        ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)*(n_CTF*1.0d0)/timing_tot
     $        /1e9
         write(6,'(A,F8.4)') ' CTF_R_S_bf_ total Gnumps: ' , timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      if (verbose.gt.1) then
         call print_all_c16(n_gamma_z,CTF_R_S_xx_,' CTF_R_S_xx_: ')
         call print_all_c16(n_gamma_z,CTF_R_S_bf_,' CTF_R_S_bf_: ')
         write(6,'(A,F8.4)') ' costheta: ' , costheta_c16_f(n_gamma_z
     $        ,CTF_R_S_xx_,CTF_R_S_bf_)
      end if !if (verbose.gt.1) then

      allocate(Z_tmp_(0:n_w_max-1))
      call cl1_c16(n_w_max,Z_tmp_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_R_T_R_M_xx_ to hold'
         write(6,'(A)') ' innerproducts for each S-M pair.'
         write(6,'(A)') ' '
      end if
      allocate(S_T_R_T_R_M_xx_(0:n_delta_v*n_gamma_z))
      call cl1_c16(n_delta_v*n_gamma_z,S_T_R_T_R_M_xx_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate S_T_R_T_R_M_xx_.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)')
     $        ' S_T_R_T_R_M_xx_(ngz + n_gamma_z*ndv) '
         write(6,'(A)') ' is equal to the innerproduct between '
         write(6,'(A)') ' R_{est}T_{est}R_{upd}T_{upd}(S) and'
         write(6,'(A)')
     $        ' dconjg(CTF).*M, where CTF is the ctf-function '
         write(6,'(A)') ' R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' and S is the template, '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
      call xc1_c16(n_A,M_p_0in_,CTF_p_0in_,M_p_)
      call rotate_p_to_p_fftw(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
     $     ,n_A,fftw_0in_,fftw_out_,M_p_,-gamma_z_est,M_p_)
      call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $     -delta_y_est,M_p_)
c$$$      call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x,
c$$$     $     -delta_y,M_p_)
      call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $     ,fftw_0in_,fftw_out_,M_p_,M_q_)
      do ndv=0,n_delta_v-1
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         call cp1_c16(n_A,S_p_0in_,S_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,S_p_,delta_x,delta_y
     $        ,S_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $        ,fftw_0in_,fftw_out_,S_p_,S_q_)
         call innerproduct_q_k_stretch_0(n_r,grid_p_,n_w_
     $        ,n_A,S_q_,M_q_,Z_tmp_)
         call cp1_c16(n_w_max,Z_tmp_,fftw_out_last_)
         call dfftw_execute_(fftw_plan_back_last)
         call cp1_c16(n_w_max,fftw_0in_last_,Z_tmp_)
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' ndv ' , ndv
            call print_sub_c16(n_w_max,Z_tmp_,' Z_tmp: ')
         end if !if (verbose.gt.1) then
         do ngz=0,n_gamma_z-1
            gamma_z = gamma_z_(ngz)
            call interp1_c16(n_w_max,0.0d0,2.0d0*pi,Z_tmp_,+gamma_z
     $           ,C_q)
            nC = ngz + n_gamma_z*ndv
            S_T_R_T_R_M_xx_(nC) = C_q
            enddo !do ngz=0,n_gamma_z-1
      enddo !do ndv=0,n_delta_v-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,A,F8.5)') ' finished calculating S_T_R_T_R_M_xx_ '
     $        , ' for each S-M pair. ' ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)*(n_M*1.0d0)/timing_tot
     $        /1e9
         write(6,'(A,F8.4)') ' S_T_R_T_R_M_xx total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_R_T_R_M_bf_ to hold'
         write(6,'(A)') ' innerproducts for each S-M pair.'
         write(6,'(A)') ' '
      end if
      allocate(S_T_R_T_R_M_bf_(0:n_delta_v*n_gamma_z))
      call cl1_c16(n_delta_v*n_gamma_z,S_T_R_T_R_M_bf_)
      timing_tic = omp_get_wtime()
      call xc1_c16(n_A,M_p_0in_,CTF_p_0in_,M_p_)
      call rotate_p_to_p_fftw(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_
     $     ,n_A,fftw_0in_,fftw_out_,M_p_,-gamma_z_est,M_p_)
      call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $     -delta_y_est,M_p_)
      do ndv=0,n_delta_v-1
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         do ngz=0,n_gamma_z-1
            gamma_z = gamma_z_(ngz)
            call cp1_c16(n_A,S_p_0in_,S_p_)
            call transf_p_to_p(n_r,grid_p_,n_w_,n_A,S_p_,delta_x
     $           ,delta_y,S_p_)
            call rotate_p_to_p_fftw(n_r,fftw_plan_frwd_,fftw_plan_back_
     $           ,n_w_,n_A,fftw_0in_,fftw_out_,S_p_,gamma_z,S_p_)
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_p_,M_p_,C_p)
            nC = ngz + n_gamma_z*ndv
            S_T_R_T_R_M_bf_(nC) = C_p
         enddo                  !do ngz=0,n_gamma_z-1
      enddo !do ndv=0,n_delta_v-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,A,F8.5)') ' finished calculating S_T_R_T_R_M_bf_ '
     $        , ' for each S-M pair. ' ,timing_tot
         timing_tmp = (n_A*1.0d0)*(n_S*1.0d0)*(n_M*1.0d0)/timing_tot
     $        /1e9
         write(6,'(A,F8.4)') ' S_T_R_T_R_M_bf total Gnumps: ' ,
     $        timing_tmp
         write(6,'(A)') ' '
      end if ! if (verbose.gt.0) then

      if (verbose.gt.0) then
         do ndv=0,n_delta_v-1
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            write(6,'(A,I0)') ' ndv ' , ndv
            n_X = n_gamma_z*ndv
            call print_sub_c16(n_gamma_z ,S_T_R_T_R_M_xx_(n_X)
     $           ,' S_T_R_T_R_M_xx_: ')
            call print_sub_c16(n_gamma_z ,S_T_R_T_R_M_bf_(n_X)
     $           ,' S_T_R_T_R_M_bf_: ')
            write(6,'(A,F8.4)') ' costheta: ' ,
     $           costheta_c16_f(n_gamma_z,S_T_R_T_R_M_xx_(n_X)
     $           ,S_T_R_T_R_M_bf_(n_X))
         enddo                  !do ndv=0,n_delta_v-1
      end if                    !if (verbose.gt.0) then

      deallocate(S_p_)
      deallocate(M_p_)
      deallocate(CTF_p_)
      deallocate(T_p_)
      deallocate(S_q_)
      deallocate(M_q_)
      deallocate(CTF_q_)
      deallocate(T_q_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nr))
         call dfftw_destroy_plan(fftw_plan_back_(nr))
      enddo !do nr=0,n_r-1
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_0in_)
      deallocate(fftw_out_)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_bruteforce_1]'
      end if      

      end
