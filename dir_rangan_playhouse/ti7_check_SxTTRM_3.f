      subroutine ti7_check_SxTTRM_3(verbose,n_delta_v
     $     ,delta_x_,delta_y_,n_gamma_z,gamma_z_,delta_x_est_
     $     ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $     ,fftw_plan_back_ ,fftw_in1_ ,fftw_out_,fftw_plan_back_last
     $     ,fftw_in1_last_ ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_
     $     ,S_q_,n_S ,I_S_sample_,ld_S ,S_p__,M_p_,M_q_,n_M ,I_M_sample_
     $     ,ld_M ,M_p__,CTF_p_ ,n_CTF ,ld_CTF ,CTF_p__,C_M_,CTF_R_S_
     $     ,S_T_T_R_CTF_M_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_delta_v,n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      integer n_S,I_S_sample_(0:n_S-1),ld_S
      integer n_M,I_M_sample_(0:n_M-1),ld_M
      integer n_CTF,ld_CTF
      real *8 delta_x_(0:0)
      real *8 delta_y_(0:0)
      real *8 gamma_z_(0:0)
      real *8 delta_x_est_(0:0)
      real *8 delta_y_est_(0:0)
      real *8 gamma_z_est_(0:0)
      real *8 ctf_ind_est_(0:0)
      real *8 grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      integer *8 fftw_plan_back_(0:n_r-1)
      integer *8 fftw_plan_back_last
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 fftw_in1_last_(0:0),fftw_out_last_(0:0)
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 M_p_(0:0),M_q_(0:0),M_p__(0:0)
      complex *16 CTF_p_(0:0),CTF_q_(0:0),CTF_p__(0:0)
      complex *16 C_M_(0:0),C_M
      complex *16 CTF_R_S_(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      complex *16, allocatable :: CTF_R_S_use_(:)
      complex *16, allocatable :: C_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16, allocatable :: ZZ1_(:)
      complex *16, allocatable :: ZZ2_(:)
      complex *16 C_tmp,Z_tmp,CTF_R_S_use,C_Z_use
      integer n_w_max,nr,nw,ns,nm,nctf,nC
      integer ndv,ngz,nx,ld_X,nx1,nx2,nx3
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti7_check_SxTTRM_3]'
      end if !if (verbose.gt.0) then

      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)

      if (verbose.gt.-1) then
         write(6,'(A)') 'calling test_innerproduct_bruteforce: '
         ns=0
         nm=0
         nctf = nint(ctf_ind_est_(nm))
         write(6,'(A,F8.4)') ' delta_x_est: ' , delta_x_est_(nm)
         write(6,'(A,F8.4)') ' delta_y_est: ' , delta_y_est_(nm)
         write(6,'(A,F8.4)') ' gamma_z_est: ' , gamma_z_est_(nm)
         call write_sub_c16(n_A,S_p__(ld_S*ns),11,' S_p__in_: ')
         call write_sub_c16(n_A,M_p__(ld_M*nm),11,' M_p__in_: ')
         call write_sub_c16(n_A,CTF_p__(ld_CTF*nctf),13,' CTF_p__in_: ')
      end if                    !if (verbose.gt.1) then
      call test_innerproduct_bruteforce_1(2,n_r ,grid_p_ ,n_w_
     $     ,0.5d0 ,S_p__(ld_S*ns),M_p__(ld_M*nm)
     $     ,CTF_p__(ld_CTF*nctf) ,delta_x_est_(nm) ,delta_y_est_(nm)
     $     ,gamma_z_est_(nm) ,n_delta_v ,delta_x_ ,delta_y_
     $     ,n_gamma_z ,gamma_z_)

      allocate(CTF_R_S_use_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,CTF_R_S_use_)
      allocate(ZZ1_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,ZZ1_)
      allocate(ZZ2_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,ZZ2_)

      allocate(C_q_(0:n_w_max-1))
      call cl1_c16(n_w_max,C_q_)
      allocate(Z_tmp_(0:n_w_max-1))
      call cl1_c16(n_w_max,Z_tmp_)
      do nx1=0,n_delta_v*n_S*n_M-1
         nx = nx1
         ndv = mod(nx,n_delta_v)
         nx = (nx - ndv)/n_delta_v
         ns = mod(nx,n_S)
         nx = (nx - ns)/n_S
         nm = mod(nx,n_M)
         nx = (nx - nm)/n_M
         if (nx.ne.0) then
            write(6,'(A,I0)') 'Warning! incorrectly unpacked nx1: '
     $           , nx1
         end if                 ! if (nx.ne.0) then
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0,A,I0)') , ' nx1: '
     $           , nx1 ,' ndv: ' , ndv , ' ns: ' ,
     $           ns ,' nm: ' , nm 
         end if                 !if (verbose.gt.2) then
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_)
c$$$         call write_sub_c16(n_A,S_p_,7,' S_p_: ')
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $        ,fftw_in1_,fftw_out_,S_p_,S_q_)
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         call get_CTF_R_S_use_(gamma_z_est,n_gamma_z,CTF_R_S_(0 + ns
     $        *n_gamma_z + nctf*n_gamma_z*n_S),CTF_R_S_use_)
         call cp1_c16(n_A,M_p__(I_M_sample_(nm)*ld_M),M_p_)
c$$$         call write_sub_c16(n_A,M_p_,7,' M_p_: ')
c$$$         call write_sub_c16(n_A,CTF_p__(nctf*ld_CTF),9,' CTF_p_: ')
         call xc1_c16(n_A,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
         call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,-gamma_z_est,M_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_ ,n_A,M_p_,-delta_x,
     $        -delta_y,M_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,M_q_)
         call innerproduct_q_k_only(n_r,grid_p_,n_w_,n_A,S_q_,M_q_
     $        ,C_q_)
         call cp1_c16(n_w_max,C_q_,fftw_out_last_)
         call dfftw_execute_(fftw_plan_back_last)
         call cp1_c16(n_w_max,fftw_in1_last_,C_q_)
c$$$         write(6,'(1(A,I0))') ' ndv ' , ndv
c$$$         call write_sub_c16(n_w_max,C_q_,6,' C_q: ')
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         C_tmp = (0.0d0,0.0d0)
         Z_tmp = (0.0d0,0.0d0)
         do nw=0,n_w_max-1
            C_tmp = C_tmp + zabs(C_q_(nw))**2
            nx2 = ndv + n_delta_v*(ns + n_S*nm)
            nx3 = nw*n_delta_v*n_S*n_M
            Z_tmp = Z_tmp + zabs(C_q_(nw) - S_T_T_R_CTF_M_q__(nx2 +
     $           nx3))**2
         enddo                  !do nw=0,n_w_max-1
         C_tmp = zsqrt(C_tmp)
         Z_tmp = zsqrt(Z_tmp)
         if (verbose.gt.-2) then
            write(6,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A,F16.8)') ,
     $           ' C_q_ vs S_T_T_R_CTF_M_q__: ' ,
     $           ' nx1: ' , nx1 ,' ndv: ' , ndv ,
     $           ' ns: ' , ns ,' nm: ' , nm , ' absolute error: ' ,
     $           real(zabs(Z_tmp)) , ' relative error: ' ,
     $           real(zabs(Z_tmp)/zabs(C_tmp))
         end if                 !if (verbose.gt.2) then            
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         call cp1_c16(n_w_max,C_q_,Z_tmp_)
         do ngz=0,n_gamma_z-1
            CTF_R_S_use = CTF_R_S_use_(ngz)
            gamma_z = gamma_z_(ngz)
            call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $           ,Z_tmp)
c$$$            nC = ndv + n_delta_v*ngz
            nC = ngz
            if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
               C_Z_use = 1.0d0
            else
               C_Z_use = C_M*CTF_R_S_use
            end if
            ZZ1_(nC) = Z_tmp/(n_r**4)/C_Z_use
         enddo !do ngz=0,n_gamma_z-1
c$$$         call write_sub_c16(n_gamma_z,ZZ1_,6,'ZZ1_: ')
         do nw=0,n_w_max-1
            nx2 = ndv + n_delta_v*(ns + n_S*nm)
            nx3 = nw*n_delta_v*n_S*n_M
            Z_tmp_(nw) = S_T_T_R_CTF_M_q__(nx2 + nx3)
         enddo !do nw=0,n_w_max-1
         do ngz=0,n_gamma_z-1
            CTF_R_S_use = CTF_R_S_use_(ngz)
            gamma_z = gamma_z_(ngz)
            call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $           ,Z_tmp)
c$$$            nC = ndv + n_delta_v*ngz
            nC = ngz
            if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
               C_Z_use = 1.0d0
            else
               C_Z_use = C_M*CTF_R_S_use
            end if
            ZZ2_(nC) = Z_tmp/(n_r**4)/C_Z_use
         enddo !do ngz=0,n_gamma_z-1
c$$$         call write_sub_c16(n_gamma_z,ZZ2_,6,'ZZ2_: ')
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         C_tmp = (0.0d0,0.0d0)
         Z_tmp = (0.0d0,0.0d0)
         do ngz=0,n_gamma_z-1
            C_tmp = C_tmp + zabs(ZZ1_(ngz))**2
            Z_tmp = Z_tmp + zabs(ZZ1_(ngz) - ZZ2_(ngz))**2
         enddo                  !do ngz=0,n_gamma_z-1
         C_tmp = zsqrt(C_tmp)
         Z_tmp = zsqrt(Z_tmp)
         if (verbose.gt.-2) then
            write(6, '(A,A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A,F16.8)')
     $           ' Correlations(C_q_) vs'
     $           ,' Correlations(S_T_T_R_CTF_M_q__): ' , ' nx1: ' , nx1
     $           ,' ndv: ' , ndv , ' ns: ' , ns
     $           ,' nm: ' , nm , ' absolute error: ' , real(zabs(Z_tmp))
     $           , ' relative error: ' , real(zabs(Z_tmp)/zabs(C_tmp))
         end if                 !if (verbose.gt.2) then
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$         write(6,'(1(A,I0))') ' ndv ' , ndv 
c$$$         call write_all_c16(n_gamma_z ,ZZ1_,7 ,' ZZ1_: ')
c$$$         call write_all_c16(n_gamma_z ,ZZ2_,7 ,' ZZ2_: ')         
      enddo                     !do nx1=0,n_delta_v*n_S*n_M-1
      
      deallocate(CTF_R_S_use_)
      deallocate(ZZ1_)
      deallocate(ZZ2_)
      deallocate(C_q_)
      deallocate(Z_tmp_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished ti7_check_SxTTRM_3]'
      end if !if (verbose.gt.0) then
      end
