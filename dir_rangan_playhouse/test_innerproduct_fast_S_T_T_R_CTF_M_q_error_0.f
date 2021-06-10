      subroutine test_innerproduct_fast_STRM_0(verbose,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,,delta_x_est_,delta_y_est_
     $     ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_ ,fftw_in1_
     $     ,fftw_out_,fftw_plan_back_last,fftw_in1_last_ ,fftw_out_last_
     $     ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S,S_p__,M_p_,M_q_,n_M
     $     ,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF,CTF_p__,CTF_M_,C_q_)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_delta_x,n_delta_y
      integer n_r,n_w_(0:n_r-1),n_A
      integer n_S,ld_S
      integer n_M,ld_M
      integer n_CTF,ld_CTF
      real *8 delta_x_(0:0)
      real *8 delta_x_(0:0)
      real *8 delta_x_est_(0:0)
      real *8 delta_y_est_(0:0)
      real *8 gamma_z_est_(0:0)
      real *8 ctf_ind_est_(0:0)
      real *8 grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1),fftw_plan_back_last
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 fftw_in1_last_(0:0),fftw_out_last_(0:0)
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 M_p_(0:0),M_q_(0:0),M_p__(0:0)
      complex *16 CTF_p_(0:0),CTF_q_(0:0),CTF_p__(0:0)
      complex *16, allocatable :: C_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16 CTF_M
      complex *16 Z_tmp
      integer n_w_max,nr,nw,nm,nctf
      integer ndx,ndy,nx,ld_X,nx1,nx2,nx3
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_STRM_0]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)

      allocate(C_q_(0:n_w_max-1))
      call cl1_c16(n_w_max,C_q_)
      allocate(Z_tmp_(0:n_w_max-1))
      call cl1_c16(n_w_max,Z_tmp_)
      do nx1=0,n_delta_x*n_delta_y*n_S*n_M-1
         nx = nx1
         ndx = mod(nx,n_delta_x)
         nx = (nx - ndx)/n_delta_x
         ndy = mod(nx,n_delta_y)
         nx = (nx - ndy)/n_delta_y
         ns = mod(nx,n_S)
         nx = (nx - ns)/n_S
         nm = mod(nx,n_M)
         nx = (nx - nm)/n_M
         if (nx.ne.0) then
            write(6,'(A,I0)') 'Warning! incorrectly unpacked nx1: '
     $           , nx1
         end if                 ! if (nx.ne.0) then
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0)') , ' nx1: '
     $           , nx1 ,' ndx: ' , ndx ,  ' ndy: ' , ndy ,  ' ns: ' ,
     $           ns ,' nm: ' , nm 
         end if                 !if (verbose.gt.2) then
         delta_x = delta_x_(ndx)
         delta_y = delta_y_(ndy)
         call cp1_c16(n_A,S_p__(ns*ld_S),S_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_ ,n_A,S_p_,+delta_x,
     $        +delta_y,S_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $        ,fftw_in1_,fftw_out_,S_p_,S_q_)
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         call cp1_c16(n_A,M_p__(nm*ld_M),M_p_)
         call xc1_c16(n_A,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
         call rotate_p2p_fz(n_r,n_w_,n_A,M_p_,-gamma_z_est,M_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,M_q_)
         call innerproduct_q_k_only(n_r,grid_p_,n_w_,n_A,S_q_,M_q_
     $        ,C_q_)
         call cp1_c16(n_w_max,C_q_,fftw_out_last_)
         call dfftw_execute_(fftw_plan_back_last)
         call cp1_c16(n_w_max,fftw_in1_last_,C_q_)
         C_tmp = (0.0d0,0.0d0)
         Z_tmp = (0.0d0,0.0d0)
         do nw=0,n_w_max-1
            C_tmp = C_tmp + zabs(C_q_(nw))**2
            nx2 = ndx + n_delta_x*(ndy + n_delta_y*(nS + n_S*nm))
            nx3 = nw*n_delta_x*n_delta_y*n_S*n_M
            Z_tmp = Z_tmp + zabs(C_q_(nw) - S_T_T_R_CTF_M_q__(nx2 +
     $           nx3))**2
         enddo                  !do nw=0,n_w_max-1
         C_tmp = zsqrt(C_tmp)
         Z_tmp = zsqrt(Z_tmp)
         if (verbose.gt.-2) then
            write(6,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,F16.8,A,F16.8)') ,
     $           ' nx1: ' , nx1 ,' ndx: ' , ndx ,  ' ndy: ' , ndy ,
     $           ' ns: ' , ns ,' nm: ' , nm , ' absolute error: ' ,
     $           real(zabs(Z_tmp)) , ' relative error: ' ,
     $           real(zabs(Z_tmp)/zabs(C_tmp))
         end if                 !if (verbose.gt.2) then            
      enddo                     !do nx1=0,n_delta_x*n_delta_y*n_S*n_M-1

      deallocate(C_q_)
      deallocate(Z_tmp_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_fast_STRM_0]'
      end if !if (verbose.gt.0) then
      end
