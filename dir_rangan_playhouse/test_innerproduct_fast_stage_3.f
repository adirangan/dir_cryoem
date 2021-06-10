c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_fast_stage_3(verbose,rseed,n_omp_sub
     $     ,svd_calculation_type,eps_svd,n_svd_d,svd_d_,n_svd_r,svd_r_
     $     ,n_svd_l,svd_l_ ,svd_U_d_,svd_s_,svd_V_r_,Z_S_svdd_,Z_M_svdd_
     $     ,n_r ,grid_p_ ,n_w_ ,half_diameter_x_c ,n_S,I_S_sample_,ld_S
     $     ,S_p__ ,S_alpha_polar_a_ ,S_alpha_azimu_b_ ,n_M,I_M_sample_
     $     ,ld_M ,M_p__ ,n_M_sub_use,n_M_per_,I_M_per_,n_CTF ,ld_CTF
     $     ,CTF_p__ ,CTF_R_S_ ,S_q__,T_S_q__,Z_S_q__ ,T_R_CTF_M_q__
     $     ,T_T_R_CTF_M_q__ ,Z_T_R_CTF_M_q__,S_p_ ,S_q_ ,S_p_omp__
     $     ,S_q_omp__ ,M_p_,M_q_ ,M_p_omp__,M_q_omp__ ,CTF_p_ ,CTF_q_
     $     ,CTF_p_omp__ ,CTF_q_omp__ ,polar_a_est_ ,azimu_b_est_
     $     ,gamma_z_est_ ,delta_x_est_ ,delta_y_est_ ,l2_norm_est_
     $     ,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_
     $     ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $     ,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f ,flag_MS_vs_SM
     $     ,N_pixels_in ,displacement_max ,n_delta_x ,delta_x_,n_delta_y
     $     ,delta_y_ ,n_gamma_z ,gamma_z_, fftw_plan_frwd_
     $     ,fftw_plan_back_ ,fftw_in1_,fftw_out_,fpm_howmany,n_fpm
     $     ,fpm_back_ ,fpm_in1__ ,fpm_out__, C_M_ ,n_SM_use)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer *4 rseed !random seed. ;
      integer n_omp_sub ! number of batches to use when finding innerproducts
      integer svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,I_S_sample_(0:n_S-1),ld_S,n_M
     $     ,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF
      real *8 grid_p_(0:n_r-1),half_diameter_x_c
      complex *16 S_p__(0:0)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      viewing angles for templates
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      integer *4 max_i4_f,sum_i4_f
c$$$      integer *4, allocatable :: SM_use_sub_(:)
      integer *4 n_SM_use
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      complex *16 M_p__(0:0)
      complex *16 CTF_p__(0:0)
      real *8 polar_a_est_(0:n_M-1),azimu_b_est_(0:n_M-1)
     $     ,gamma_z_est_(0:n_M-1)
      real *8 delta_x_est_(0:n_M-1),delta_y_est_(0:n_M-1)
     $     ,l2_norm_est_(0:n_M-1)
      real *8 ctf_ind_est_(0:n_M-1)
      real *8 S_index_est_(0:n_M-1)
      real *8 polar_a_upd_(0:n_M-1),azimu_b_upd_(0:n_M-1)
     $     ,gamma_z_upd_(0:n_M-1)
      real *8 delta_x_upd_(0:n_M-1),delta_y_upd_(0:n_M-1)
     $     ,l2_norm_upd_(0:n_M-1)
      real *8 ctf_ind_upd_(0:n_M-1)
      real *8 S_index_upd_(0:n_M-1)
      real *8 alpha_update_f
      logical flag_MS_vs_SM
      real *8 eps_svd,eps_target,N_pixels_in
      real *8 displacement_max
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      complex *16 C_M_(0:n_M-1),C_M
c$$$      declaration of svd-expansion and associated variables
c$$$      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      integer *4 n_svd_d,n_svd_r,n_svd_l
      real *8 svd_d_(0:0)
      real *8 svd_r_(0:0)
      integer *4 svd_l_(0:0)
      real *8 svd_U_d_(0:0)
      real *8 svd_s_(0:0)
      real *8 svd_V_r_(0:0)
      integer n_svd_max
      parameter (n_svd_max=512)
      logical flag_warning
      data flag_warning / .true. /
      real *8 R_max,K_max,delta,delta_max,n_pixels
      complex *16 Z_S_svdd_(0:0)
      complex *16 Z_M_svdd_(0:0)
c$$$      declaration of innerproduct arrays
      complex *16 CTF_R_S_(0:0)
      complex *16 S_q__(0:0)
      complex *16 T_S_q__(0:0)
      complex *16 Z_S_q__(0:0)
      complex *16 S_p_(0:0)
      complex *16 S_q_(0:0)
      complex *16 S_p_omp__(0:0)
      complex *16 S_q_omp__(0:0)
      complex *16 M_p_(0:0)
      complex *16 M_q_(0:0)
      complex *16 M_p_omp__(0:0)
      complex *16 M_q_omp__(0:0)
      complex *16 CTF_p_(0:0)
      complex *16 CTF_q_(0:0)
      complex *16 CTF_p_omp__(0:0)
      complex *16 CTF_q_omp__(0:0)
      complex *16 T_R_CTF_M_q__(0:0)
      complex *16 T_T_R_CTF_M_q__(0:0)
      complex *16 Z_T_R_CTF_M_q__(0:0)
      complex *16, allocatable :: S_T_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: S_Z_T_R_CTF_M_q__(:)      
      complex *16, allocatable :: ZZ_(:)
      complex *16, allocatable :: ZZ_omp__(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16, allocatable :: C_Z_MS_(:)
      complex *16, allocatable :: C_S_MS_(:)
      real *8, allocatable :: delta_x_MS_(:)
      real *8, allocatable :: delta_y_MS_(:)
      real *8, allocatable :: gamma_z_MS_(:)
      complex *16, allocatable :: C_Z_sort_MS_(:)
      complex *16, allocatable :: C_S_sort_MS_(:)
      real *8, allocatable :: delta_x_sort_MS_(:)
      real *8, allocatable :: delta_y_sort_MS_(:)
      real *8, allocatable :: gamma_z_sort_MS_(:)
      complex *16, allocatable :: I_permute_MS_(:)
      complex *16, allocatable :: I_inverse_MS_(:)
      integer nx0,nx1,nx2,nx3
c$$$      fftw for omp not used
c$$$      fftw for local use
      integer *8 fftw_plan_frwd_(0:0)
      integer *8 fftw_plan_back_(0:0)
      complex *16 fftw_in1_(0:0)
      complex *16 fftw_out_(0:0)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*),fftw_out_last_(*)
c$$$      fftw plan many (fpm) for omp
      integer fpm_howmany,n_fpm
      integer *8 fpm_back_(0:0)
      complex *16 fpm_in1__(0:0)
      complex *16 fpm_out__(0:0)
c$$$      indices
      integer n_transf,nr,n_w_max,n_A,na,ns,nctf,nm,nw
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_(0:0)
      real *8 delta_y_(0:0)
      real *8 gamma_z_(0:0)
      character(len=64) format_string
c$$$      parameters for timing and testing
      logical flag_fill
      parameter (flag_fill=.false.)
      logical flag_test
      real *8 timing_tic,timing_toc,timing_tot,timing_tmp
c$$$      parameters for omp
      integer *4 nomp_sub,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp
      integer *4 n_M_sub_use,nM_sub,nM_per,IM_per
      integer *4 n_M_per_(0:0)
      integer *4 I_M_per_(0:0)
      integer *4 nM_per_max
c$$$      pi
      real *8 pi
c$$$      parameters for memory map
      real *8 d_mem
      integer verbose_mem
      parameter (verbose_mem=1)

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_stage_3]'
      end if      

      eps_target = eps_svd
      if (svd_calculation_type.eq.2) then
         n_transf = n_delta_x*n_delta_y
      else
         n_transf = n_svd_l
      end if !if (svd_calculation_type.eq.2) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_omp_sub: ',n_omp_sub
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
         write(6,'(A,F6.3)') 'eps_svd: ',eps_svd
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,I0)') 'n_r: ',n_r
         write(format_string,'(A,I0,A)') '(A,',n_r,'(F8.4,1X))'
         write(6,format_string) 'grid_p_: ',(grid_p_(nr),nr=0 ,n_r-1)
         write(format_string,'(A,I0,A)') '(A,',n_r,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0 ,n_r-1)
         write(6,'(A,F6.3)') 'half_diameter_x_c: ',half_diameter_x_c
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,'(A)') 'S_alpha_polar_a_: '
         write(6,format_string) (S_alpha_polar_a_(nr),nr=0,n_S-1)
         write(6,'(A)') 'S_alpha_azimu_b_: '
         write(6,format_string) (S_alpha_azimu_b_(nr),nr=0,n_S-1)
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(6,'(A,I0)') 'n_CTF: ',n_CTF
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         write(6,'(A,F6.3)') 'n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') 'displacement_max: ',displacement_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'fpm_howmany: ',fpm_howmany
         write(6,'(A,I0)') 'n_fpm: ',n_fpm
         write(6,'(A,I0)') 'n_transf: ',n_transf
      end if                    !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A)') ' Generating indices'
         write(6,'(A)') ' '
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
         write(6,'(A)') 'linking fftw_plans for local use'
         write(6,'(A)') ' '
      end if
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      svd_calculation_type.eq.2: STxTRM
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (svd_calculation_type.eq.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_T_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and T_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      flag_RTRT_vs_RTTR = .true.
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_transf*1.0d0*n_S*1.0d0
     $     *n_M*1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(S_T_T_R_CTF_M_q__(0:n_transf*n_S*n_M*n_w_max))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),T_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_T_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      end if ! if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_STxTRM_1 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call test_innerproduct_fast_STxTRM_1(verbose-1,n_delta_x
     $        ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $        ,delta_x_est_,delta_y_est_ ,gamma_z_est_,ctf_ind_est_
     $        ,fftw_plan_frwd_,fftw_in1_ ,fftw_out_,fftw_plan_back_last
     $        ,fftw_in1_last_,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_
     $        ,S_q_,n_S,I_S_sample_,ld_S,S_p__,M_p_,M_q_,n_M
     $        ,I_M_sample_,ld_M ,M_p__,CTF_p_ ,n_CTF ,ld_CTF,CTF_p__
     $        ,C_M_,CTF_R_S_ ,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.eq.2) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      svd_calculation_type.eq.2: SxTTRM
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (svd_calculation_type.eq.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Redefining array S_T_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_T_R_CTF_M_q__ and S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      flag_RTRT_vs_RTTR = .true.
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_transf*1.0d0*n_S*1.0d0
     $     *n_M*1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      call cl1_c16(n_transf*n_S*n_M*n_w_max,S_T_T_R_CTF_M_q__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Recalculating S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_S
         nx2 = nw*n_r*n_transf*n_M
         nx3 = nw*n_S*n_transf*n_M
         call zgemm('T','N',n_S,n_transf*n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),S_q__(nx1),n_r,T_T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_T_T_R_CTF_M_q__(nx3)
     $        ,n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_S*(0 + n_transf*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_S*n_transf*nM_per,n_S*n_transf*n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      end if ! if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_S*(0 + n_transf*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_S*n_transf*nM_per,n_S*n_transf*n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' transposing dimensions (1,2) of S_T_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      do nx0=0,n_M*n_w_max-1
         nx1 = nx0*n_S*n_transf
         call trn0_c16(n_S,n_transf,S_T_T_R_CTF_M_q__(nx1)
     $        ,S_T_T_R_CTF_M_q__(nx1))
      enddo !do nx0=0,n_M*n_w_max-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished transposing. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' transpose total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_STxTRM_1 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call test_innerproduct_fast_SxTTRM_1(verbose-1,n_delta_x
     $        ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $        ,delta_x_est_,delta_y_est_ ,gamma_z_est_,ctf_ind_est_
     $        ,fftw_plan_frwd_,fftw_in1_ ,fftw_out_,fftw_plan_back_last
     $        ,fftw_in1_last_,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_
     $        ,S_q_,n_S,I_S_sample_,ld_S,S_p__,M_p_,M_q_,n_M
     $        ,I_M_sample_,ld_M ,M_p__,CTF_p_ ,n_CTF ,ld_CTF,CTF_p__
     $        ,C_M_,CTF_R_S_ ,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.eq.2) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      svd_calculation_type.eq.1: SZxTRM
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (svd_calculation_type.ne.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Allocating array S_Z_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' T_R_CTF_M_q__ and Z_S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      flag_RTRT_vs_RTTR = .true.
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_transf*1.0d0*n_S*1.0d0 *n_M*1.0d0*n_w_max
     $     *16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(S_Z_T_R_CTF_M_q__(0:n_transf*n_S*n_M*n_w_max))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_transf*n_S
         nx2 = nw*n_r*n_M
         nx3 = nw*n_transf*n_S*n_M
         call zgemm('T','N',n_transf*n_S,n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_q__(nx1),n_r,T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,n_transf*n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      end if !if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_transf*(0 + n_S*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_transf *n_S *nM_per,n_transf *n_S *n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_transf*1.0d0)*(n_S*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
      if (verbose.gt.1) then
         write(6,'(A)') ' Now allocating array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_S_svdd_. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_delta_x*n_delta_y*1.0d0*n_S*1.0d0 *n_M
     $     *1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(S_T_T_R_CTF_M_q__(0:n_delta_x*n_delta_y*n_S*n_M*n_w_max))
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_A_tmp,n_X_tmp)
c$OMP DO 
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_x*n_delta_y*n_S*n_M*nw
         n_X_tmp = n_transf*n_S*n_M*nw
         call zgemm('T','N',n_delta_x*n_delta_y,n_S*n_M*1,n_svd_l ,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_S_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q__(n_A_tmp),n_delta_x *n_delta_y)
      enddo ! do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*n_delta_y*1.0d0)*(n_S
     $     *1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_SZxTRM_1 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call  test_innerproduct_fast_SZxTRM_1(verbose-1,n_svd_r,svd_r_
     $        ,n_svd_l,svd_l_,svd_s_,svd_V_r_,Z_S_svdd_,n_delta_x
     $        ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $        ,delta_x_est_ ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_
     $        ,fftw_plan_frwd_ ,fftw_in1_ ,fftw_out_,fftw_plan_back_last
     $        ,fftw_in1_last_ ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_
     $        ,S_q_,n_S,I_S_sample_,ld_S ,S_p__,M_p_,M_q_,n_M
     $        ,I_M_sample_,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF ,CTF_p__
     $        ,C_M_,CTF_R_S_ ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.ne.2) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      svd_calculation_type.eq.1: SxZTRM
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (svd_calculation_type.ne.2) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Reallocating array S_Z_T_R_CTF_M_q__ to hold '
         write(6,'(A)') ' bessel-coefficients for the product of '
         write(6,'(A)') ' Z_T_R_CTF_M_q__ and S_q__ '
         write(6,'(A)') ' integrated over k (i.e., nr).'
         write(6,'(A)') ' '
      end if
      flag_RTRT_vs_RTTR = .false.
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_S*1.0d0*n_transf*1.0d0*n_M*1.0d0*n_w_max
     $     *16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_Z_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      call cl1_c16(n_S*n_transf*n_M*n_w_max,S_Z_T_R_CTF_M_q__)
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculating S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nx1,nx2,nx3)
c$OMP DO 
      do nw=0,n_w_max-1
         nx1 = nw*n_r*n_S
         nx2 = nw*n_r*n_transf*n_M
         nx3 = nw*n_S*n_transf*n_M
         call zgemm('T','N',n_S,n_transf*n_M,n_r,1.0d0
     $        *cmplx(1.0d0,0.0d0),S_q__(nx1),n_r,Z_T_R_CTF_M_q__(nx2)
     $        ,n_r ,0.0d0*cmplx(1.0d0,0.0d0),S_Z_T_R_CTF_M_q__(nx3)
     $        ,n_S)
      enddo !do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      n_SM_use = n_S*n_M
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_r*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)*(n_M
     $     *1.0d0)*(n_w_max*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_Z_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
c$$$      %%%%%%%%%%%%%%%%
      if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling array fill for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_S*(0 + n_transf*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.false.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_S*n_transf*nM_per,n_S*n_transf*n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished array fill. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' array fill total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
c$$$      %%%%%%%%%%%%%%%%
      end if !if (flag_fill) then
c$$$      %%%%%%%%%%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' calling fftw_plan_many for S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,n_F_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_A
         n_F_tmp = nM_sub*n_fpm
         n_X_tmp = 0 + n_S*(0 + n_transf*(IM_per + n_M*0))
         call dfftw_block_many_1(verbose-1,.true.,fpm_howmany,n_w_max
     $        ,fpm_in1__(n_F_tmp),fpm_out__(n_F_tmp) ,fpm_back_(nM_sub)
     $        ,n_S*n_transf*nM_per,n_S*n_transf*n_M
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished fftw_plan_many. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' fftw_plan_many total Gnumps: ' , timing_tmp
      write(6,'(A)') ' '
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' transposing dimensions (1,2) of S_Z_T_R_CTF_M_q__. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
      do nx0=0,n_M*n_w_max-1
         nx1 = nx0*n_S*n_transf
         call trn0_c16(n_S,n_transf,S_Z_T_R_CTF_M_q__(nx1)
     $        ,S_Z_T_R_CTF_M_q__(nx1))
      enddo !do nx0=0,n_M*n_w_max-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished transposing. ' , timing_tot
      timing_tmp = (n_w_max*1.0d0)*(n_S*1.0d0)*(n_transf*1.0d0)
     $     *(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' transpose total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      if (verbose.gt.1) then
         write(6,'(A)') ' Now reallocating array S_T_T_R_CTF_M_q__ '
         write(6,'(A)') ' to hold product of '
         write(6,'(A)') ' S_Z_T_R_CTF_M_q__ and Z_M_svdd_. '
         write(6,'(A)') ' '
      end if
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_delta_x*n_delta_y*1.0d0*n_S*1.0d0 *n_M
     $     *1.0d0*n_w_max*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' S_T_T_R_CTF_M_q__ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      call cl1_c16(n_delta_x*n_delta_y*n_S*n_M*n_w_max
     $     ,S_T_T_R_CTF_M_q__)
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_A_tmp,n_X_tmp)
c$OMP DO 
      do nw=0,n_w_max-1
         n_A_tmp = n_delta_x*n_delta_y*n_S*n_M*nw
         n_X_tmp = n_transf*n_S*n_M*nw
         call zgemm('T','N',n_delta_x*n_delta_y,n_S*n_M*1,n_svd_l ,1.0d0
     $        *cmplx(1.0d0,0.0d0),Z_M_svdd_,n_svd_l
     $        ,S_Z_T_R_CTF_M_q__(n_X_tmp) ,n_svd_l,0.0d0*cmplx(1.0d0
     $        ,0.0d0) ,S_T_T_R_CTF_M_q__(n_A_tmp),n_delta_x *n_delta_y)
      enddo ! do nw=0,n_w_max-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total time: ' , timing_tot
      timing_tmp = (n_svd_l*1.0d0)*(n_delta_x*n_delta_y*1.0d0)*(n_S
     $     *1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' S_T_T_R_CTF_M_q__ total Gflops: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      flag_test = .false.
      if (flag_test) then
         if (verbose.gt.1) then
            write(6,'(A)')
     $           ' calling test_innerproduct_fast_SxZTRM_1 to test. '
            write(6,'(A)') ' '
         end if !if (verbose.gt.1) then
         call  test_innerproduct_fast_SxZTRM_1(verbose-1,n_svd_d,svd_d_
     $        ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_U_d_,svd_s_,svd_V_r_
     $        ,Z_S_svdd_ ,Z_M_svdd_,n_delta_x ,delta_x_,n_delta_y
     $        ,delta_y_ ,n_gamma_z,gamma_z_,delta_x_est_ ,delta_y_est_
     $        ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_ ,fftw_in1_
     $        ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $        ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S
     $        ,I_S_sample_,ld_S ,S_p__,M_p_,M_q_,n_M,I_M_sample_,ld_M
     $        ,M_p__,CTF_p_,n_CTF ,ld_CTF ,CTF_p__,C_M_,CTF_R_S_
     $        ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
      end if !testing S_T_T_R_CTF_M_q__
      end if !if (svd_calculation_type.ne.2) then
      
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      Now update alpha: SM 
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (flag_MS_vs_SM.eqv..false.) then
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_delta_x*1.0d0*n_delta_y*1.0d0*n_gamma_z
     $     *1.0d0*n_omp_sub*16.0d0
      if (verbose.gt.0 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' ZZ_omp_ requires ' ,
     $        d_mem*1.0d-6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if                    !if (verbose.gt.0) then
      allocate(ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(ZZ_omp__(0:n_delta_x*n_delta_y*n_gamma_z*n_omp_sub-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' SM: matching templates to images; '
         write(6,'(A)')
     $        ' calling test_innerproduct_fast_ZZ_1. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp,rseed)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_delta_x*n_delta_y*n_gamma_z
         n_X_tmp = 0 + n_delta_x*n_delta_y*(0 + n_S*(IM_per + n_M*0))
         call test_innerproduct_fast_ZZ_3(verbose-1,rseed
     $        ,flag_RTRT_vs_RTTR,n_delta_x,delta_x_,n_delta_y,delta_y_
     $        ,n_gamma_z,gamma_z_,n_r,grid_p_,n_w_,n_A ,n_S
     $        ,S_alpha_polar_a_,S_alpha_azimu_b_,nM_per
     $        ,polar_a_est_(IM_per) ,azimu_b_est_(IM_per)
     $        ,gamma_z_est_(IM_per) ,delta_x_est_(IM_per)
     $        ,delta_y_est_(IM_per) ,l2_norm_est_(IM_per)
     $        ,ctf_ind_est_(IM_per) ,S_index_est_(IM_per)
     $        ,polar_a_upd_(IM_per) ,azimu_b_upd_(IM_per)
     $        ,gamma_z_upd_(IM_per) ,delta_x_upd_(IM_per)
     $        ,delta_y_upd_(IM_per) ,l2_norm_upd_(IM_per)
     $        ,ctf_ind_upd_(IM_per) ,S_index_upd_(IM_per),alpha_update_f
     $        ,displacement_max ,C_M_(IM_per),CTF_R_S_,n_M
     $        ,S_T_T_R_CTF_M_q__(n_X_tmp),ZZ_omp__(n_A_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished ZZ. ' , timing_tot
      timing_tmp = (n_gamma_z*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $     *(n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' ZZ total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      deallocate(ZZ_)
      deallocate(ZZ_omp__)
      end if !if (flag_MS_vs_SM.eqv..false.) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      Now update alpha: MS 
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (flag_MS_vs_SM.eqv..true.) then
      allocate(ZZ_(0:n_delta_x*n_delta_y*n_gamma_z-1))
      allocate(ZZ_omp__(0:n_delta_x*n_delta_y*n_gamma_z*n_omp_sub-1))
      allocate(C_Z_MS_(0:n_M*n_S-1))
      allocate(C_S_MS_(0:n_M*n_S-1))
      allocate(delta_x_MS_(0:n_M*n_S-1))
      allocate(delta_y_MS_(0:n_M*n_S-1))
      allocate(gamma_z_MS_(0:n_M*n_S-1))
      allocate(C_Z_sort_MS_(0:n_M*n_S-1))
      allocate(C_S_sort_MS_(0:n_M*n_S-1))
      allocate(delta_x_sort_MS_(0:n_M*n_S-1))
      allocate(delta_y_sort_MS_(0:n_M*n_S-1))
      allocate(gamma_z_sort_MS_(0:n_M*n_S-1))
      allocate(I_permute_MS_(0:n_M*n_S-1))
      allocate(I_inverse_MS_(0:n_M*n_S-1))
      if (verbose.gt.1) then
         write(6,'(A)')
     $        ' MS: matching images to templates; '
         write(6,'(A)')
     $        ' calling test_innerproduct_fast_MS_1. '
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(nM_per,IM_per,n_r_tmp,n_A_tmp,n_X_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub_use-1
         nM_per = n_M_per_(nM_sub)
         IM_per = I_M_per_(nM_sub)
         n_r_tmp = nM_sub*n_r
         n_A_tmp = nM_sub*n_delta_x*n_delta_y*n_gamma_z
         n_X_tmp = 0 + n_delta_x*n_delta_y*(0 + n_S*(IM_per + n_M*0))
         call test_innerproduct_fast_MS_1(verbose-1,flag_RTRT_vs_RTTR
     $        ,n_delta_x,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_
     $        ,n_r ,grid_p_,n_w_,n_A ,n_S,S_alpha_polar_a_
     $        ,S_alpha_azimu_b_,nM_per,polar_a_est_(IM_per)
     $        ,azimu_b_est_(IM_per),gamma_z_est_(IM_per)
     $        ,delta_x_est_(IM_per),delta_y_est_(IM_per)
     $        ,l2_norm_est_(IM_per),ctf_ind_est_(IM_per)
     $        ,S_index_est_(IM_per),delta_x_MS_(IM_per),
     $        delta_y_MS_(IM_per),gamma_z_MS_(IM_per), C_S_MS_(IM_per),
     $        C_Z_MS_(IM_per),displacement_max ,C_M_(IM_per),CTF_R_S_
     $        ,n_M,S_T_T_R_CTF_M_q__(n_X_tmp),ZZ_omp__(n_A_tmp))
      enddo !do nM_sub=0,n_M_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_1. ' , timing_tot
      timing_tmp = (n_gamma_z*1.0d0)*(n_delta_x*1.0d0)*(n_delta_y*1.0d0)
     $     *(n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_1 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_MS_sort_0(verbose,n_S,n_M
     $     ,delta_x_MS_,delta_y_MS_,gamma_z_MS_,C_S_MS_,C_Z_MS_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_sort_0. ' , timing_tot
      timing_tmp = (n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_sort_0 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_MS_update_2(verbose,rseed
     $     ,flag_RTRT_vs_RTTR,n_S,n_M,C_M_,delta_x_sort_MS_
     $     ,delta_y_sort_MS_,gamma_z_sort_MS_,C_S_sort_MS_,C_Z_sort_MS_
     $     ,I_permute_MS_,I_inverse_MS_,S_alpha_polar_a_
     $     ,S_alpha_azimu_b_ ,gamma_z_est_,delta_x_est_ ,delta_y_est_
     $     ,polar_a_upd_,azimu_b_upd_,gamma_z_upd_ ,delta_x_upd_
     $     ,delta_y_upd_,l2_norm_upd_,S_index_upd_)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') ' finished MS_update_2. ' , timing_tot
      timing_tmp = (n_S*1.0d0)*(n_M*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') ' MS_update_2 total Gnumps: ' ,
     $     timing_tmp
      write(6,'(A)') ' '
      deallocate(ZZ_)
      deallocate(ZZ_omp__)
      deallocate(C_Z_MS_)
      deallocate(C_S_MS_)
      deallocate(delta_x_MS_)
      deallocate(delta_y_MS_)
      deallocate(gamma_z_MS_)
      deallocate(C_Z_sort_MS_)
      deallocate(C_S_sort_MS_)
      deallocate(delta_x_sort_MS_)
      deallocate(delta_y_sort_MS_)
      deallocate(gamma_z_sort_MS_)
      deallocate(I_permute_MS_)
      deallocate(I_inverse_MS_)
      end if !if (flag_MS_vs_SM.eqv..true.) then

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      finished update alpha 
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if

      deallocate(S_T_T_R_CTF_M_q__)
      if (svd_calculation_type.ne.2) then
         deallocate(S_Z_T_R_CTF_M_q__)
      end if ! if (svd_calculation_type.ne.2) then

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_fast_stage_3]'
      end if

      end


