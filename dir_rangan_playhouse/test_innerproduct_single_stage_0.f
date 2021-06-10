      subroutine test_innerproduct_single_stage_0(
     $      verbose
     $     ,nm,ns
     $     ,n_r,grid_p_,n_w_,n_S,ld_S,n_A,nctf,n_w_max
     $     ,S_p_,S_p,S_q,M_q_,Z_q_
     $     ,C_M,C_S_,C_S_optimal,C_Z_optimal
     $     ,svd_calculation_type
     $     ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_
     $     ,Z_svdd_,Z_svdr_,Z_tmpC_
     $     ,displacement_max 
     $     ,delta_x,delta_y,gamma_z
     $     ,n_delta_x,n_delta_y,n_gamma_z
     $     ,delta_x_,delta_y_,gamma_z_
     $     ,delta_x_est_,delta_y_est_,gamma_z_est_
     $     ,ndx_optimal,ndy_optimal,ngz_optimal 
     $     ,delta_x_optimal,delta_y_optimal,gamma_z_optimal
     $     ,fftw_plan_frwd_
     $     ,fftw_plan_back_
     $     ,fftw_in1_
     $     ,fftw_out_
     $     ,timing_tic_1,timing_toc_1,timing_tot_1
     $     ,timing_tic_2,timing_toc_2,timing_tot_2
     $     ,timing_tic_3,timing_toc_3,timing_tot_3
     $     ,timing_tic_4,timing_toc_4,timing_tot_4
     $     )
      implicit none
      integer verbose
      integer svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,ld_S
      real *8 grid_p_(0:n_r-1)
      complex *16 S_p_(0:0)
      complex *16 M_p_(0:0)
      real *8 gamma_z_est_(0:0)
      real *8 delta_x_est_(0:0),delta_y_est_(0:0)
      real *8 displacement_max
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8 delta_x_optimal
      real *8 delta_y_optimal
      real *8 gamma_z_optimal
      complex *16 C_M_(0:0)
      complex *16 C_S_optimal
      complex *16 C_Z_optimal
      integer *4 n_svd_r,n_svd_d,n_svd_l
      real *8 svd_r_(0:0)
      real *8 svd_d_(0:0)
      integer *4 svd_l_(0:0)
      real *8 svd_s_(0:0)
      real *8 svd_V_r_(0:0)
      complex *16 Z_svdd_(0:0)
      complex *16 C_S_(0:0)
      complex *16 Z_svdr_(0:0)
      complex *16 Z_tmpC_(0:0)
      integer svd_l
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
      integer ns,nm,nctf,n_w_max,n_A
      complex *16 S_p(0:0)
      complex *16 S_q(0:0)
      complex *16 M_q_(0:0)
c$$$      array of displacements and rotations to measure
      integer ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_(0:0)
      real *8 delta_y_(0:0)
      real *8 gamma_z_(0:0)
c$$$      array of innerproducts to hold measurements
      complex *16 C_S,C_M,C_Z,Z_q
      complex *16 Z_q_(0:0)
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      real *8 timing_tic_1,timing_toc_1,timing_tot_1
      real *8 timing_tic_2,timing_toc_2,timing_tot_2
      real *8 timing_tic_3,timing_toc_3,timing_tot_3
      real *8 timing_tic_4,timing_toc_4,timing_tot_4

      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

      timing_tic_1 = second()
      if (verbose.gt.1) then
         write(6,'(A)') ' Transforming S_p.'
      end if
      call cp1_c16(n_A,S_p_(ns*ld_S),S_p)
      call transf_p_to_p(n_r,grid_p_,n_w_,n_A,S_p,+delta_x, +delta_y
     $     ,S_p)
      call rotate_p2p_fz(n_r,n_w_,n_A,S_p,+gamma_z,S_p)            
      timing_toc_1 = second()
      timing_tot_1 = timing_tot_1 + (timing_toc_1 - timing_tic_1)
      if (.false.) then
      else if (svd_calculation_type.eq.0) then
         timing_tic_2 = second();
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A ,fftw_in1_
     $        ,fftw_out_,S_p,S_q)
         call innerproduct_q_k_svdr(n_svd_r,svd_r_,n_svd_l,svd_l_
     $        ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q ,M_q_ ,Z_svdr_)
         timing_toc_2 = second()
         timing_tot_2 = timing_tot_2 + (timing_toc_2 - timing_tic_2)
         timing_tic_3 = second()
         call test_innerproduct_batch_excerpt_1(n_delta_x ,n_delta_y
     $        ,n_gamma_z,gamma_z_,n_w_max,n_svd_l ,Z_svdd_,Z_svdr_
     $        ,Z_tmpC_,fftw_plan_back_last ,fftw_in1_last_
     $        ,fftw_out_last_,n_r,C_M,C_S_(0 + ns *n_gamma_z + nctf
     $        *n_gamma_z*n_S),+gamma_z,Z_q_)
         timing_toc_3 = second()
         timing_tot_3 = timing_tot_3 + (timing_toc_3 - timing_tic_3)
      else if (svd_calculation_type.eq.1) then
         timing_tic_2 = second();
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A ,fftw_in1_
     $        ,fftw_out_,S_p,S_q)
         call innerproduct_q_k_svdr(n_svd_r,svd_r_,n_svd_l,svd_l_
     $        ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q ,M_q_ ,Z_svdr_)
         do svd_l=0,n_svd_l-1
            call cps_c16(n_w_max,Z_svdr_(svd_l),n_svd_l ,fftw_out_last_
     $           ,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(n_w_max,fftw_in1_last_,1,Z_svdr_(svd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
         timing_toc_2 = second()
         timing_tot_2 = timing_tot_2 + (timing_toc_2 - timing_tic_2)
         timing_tic_3 = second()
         call test_innerproduct_batch_excerpt_2(n_delta_x ,n_delta_y
     $        ,n_gamma_z,gamma_z_,n_w_max,n_svd_l ,Z_svdd_,Z_svdr_
     $        ,Z_tmpC_,n_r,C_M,C_S_(0 + ns *n_gamma_z + nctf*n_gamma_z
     $        *n_S),+gamma_z,Z_q_)
         timing_toc_3 = second()
         timing_tot_3 = timing_tot_3 + (timing_toc_3 - timing_tic_3)
      else if (svd_calculation_type.eq.2) then
         timing_tic_3 = second()
         call test_innerproduct_batch_excerpt_3(n_delta_x ,delta_x_
     $        ,n_delta_y,delta_y_,n_gamma_z,gamma_z_ ,fftw_plan_frwd_
     $        ,fftw_in1_,fftw_out_ ,fftw_plan_back_last,fftw_in1_last_
     $        ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,C_M,C_S_(0 + ns
     $        *n_gamma_z + nctf*n_gamma_z*n_S),+gamma_z,S_p,M_q_,Z_q_)
         timing_toc_3 = second()
         timing_tot_3 = timing_tot_3 + (timing_toc_3 -
     $        timing_tic_3)
      end if                    !svd_calculation_type
      timing_tic_4 = second()
      call test_innerproduct_batch_excerpt_4(delta_x_est_(nm)
     $     ,delta_y_est_(nm),gamma_z_est_(nm),delta_x_,delta_y_
     $     ,displacement_max ,n_delta_x,n_delta_y,n_gamma_z,Z_q_
     $     ,ndx_optimal,ndy_optimal ,ngz_optimal ,C_Z_optimal)
      delta_x_optimal = delta_x_(ndx_optimal)
      delta_y_optimal = delta_y_(ndy_optimal)
      gamma_z_optimal = gamma_z_(ngz_optimal)
      call test_innerproduct_batch_excerpt_5(n_gamma_z,gamma_z_ ,C_S_(0
     $     +ns*n_gamma_z + nctf*n_gamma_z*n_S),+gamma_z ,ngz_optimal
     $     ,C_S_optimal)
      timing_toc_4 = second()
      timing_tot_4 = timing_tot_4 + (timing_toc_4 - timing_tic_4)
      if (verbose.gt.2) then
         write(6,'(A,I0,A,I2,1X,I2,1X,I2,A,2F8.3,A,2F8.3)')
     $        ' Match to S_q_(',ns
     $        ,') at ndx_optimal,ndy_optimal,ngz_optimal: ' ,ndx_optimal
     $        ,ndy_optimal ,ngz_optimal ,'; C_S_optimal: ',C_S_optimal
     $        ,'; C_Z_optimal: ' ,C_Z_optimal
      end if                    !verbose
      end
