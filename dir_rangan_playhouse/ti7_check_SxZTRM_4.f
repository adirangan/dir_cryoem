      subroutine ti7_check_SxZTRM_4(verbose,svd_d_max,n_svd_d,svd_d_
     $     ,svd_r_max,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_U_d_,svd_s_
     $     ,svd_V_r_ ,Z_M_svdd_,n_delta_v,delta_x_ ,delta_y_ ,n_gamma_z
     $     ,gamma_z_ ,delta_x_est_,delta_y_est_ ,gamma_z_est_
     $     ,ctf_ind_est_ ,fftw_plan_frwd_,fftw_plan_back_,fftw_in1_
     $     ,fftw_out_ ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_
     $     ,n_r ,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,I_S_sample_,ld_S ,S_p__
     $     ,M_p_ ,M_q_,n_M ,I_M_sample_,ld_M,M_p__,CTF_p_ ,n_CTF ,ld_CTF
     $     ,CTF_p__,C_M_,CTF_R_S_,S_Z_T_R_CTF_M_q__ ,S_T_T_R_CTF_M_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_svd_d,n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_d_(0:n_svd_d-1),svd_r_(0:n_svd_r-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      real *8 svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      complex *16, allocatable :: Z_S_svdd_(:)
      complex *16 Z_M_svdd_(0:0)
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
      complex *16 S_Z_T_R_CTF_M_q__(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      real *8 svd_d_max,svd_r_max
      complex *16, allocatable :: CTF_R_S_use_(:)
      complex *16, allocatable :: ZZ1_(:)
      complex *16, allocatable :: ZZ2_(:)
      complex *16, allocatable :: Z_S_svdr_(:)
      complex *16, allocatable :: Z_M_svdr_(:)
      complex *16, allocatable :: C_S_q_(:)
      complex *16, allocatable :: C_M_q_(:)
      complex *16, allocatable :: S_R_T_x_T_R_CTF_M_q__(:)
      complex *16, allocatable :: S_T_R_x_T_R_CTF_M_q__(:)
      complex *16, allocatable :: T_p_(:)
      complex *16, allocatable :: T_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16 C_tmp,Z_tmp,C_Z_use,CTF_R_S_use
      complex *16 C_S_tmp,C_M_tmp,Z_S_tmp,Z_M_tmp,Z_X_tmp
      complex *16 costheta_c16_f
      integer n_w_max,nr,nw,ns,nm,nctf,nC
      integer ndv,ngz,nx,ld_X,nx1,nx2,nx3
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      integer nsvd_l,nC1,nC2,nC3,nC4
      real *8 pi
      integer *4 MDA_n_d,MDA_d_(0:127)
      character(len=1024) MDA_string
      integer max_i4_f
      integer ns_dump
      integer nm_dump
      integer ndelta_dump
      integer ngamma_dump
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti7_check_SxZTRM_4]'
      end if !if (verbose.gt.0) then

      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      
c$$$      dump output to dir_tmp. ;
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/svd_d_max.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_d_max,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_svd_d.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_svd_d,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/svd_r_max.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_r_max,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_svd_r.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_svd_r,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_svd_l.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_svd_l,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_l
      write(MDA_string,'(A)') './dir_tmp/svd_l_.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,svd_l_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_d
      write(MDA_string,'(A)') './dir_tmp/svd_d_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_d_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_r
      write(MDA_string,'(A)') './dir_tmp/svd_r_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_r_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_d*n_svd_l
      write(MDA_string,'(A)') './dir_tmp/svd_U_d_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_U_d_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_l
      write(MDA_string,'(A)') './dir_tmp/svd_s_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_s_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_r*n_svd_l
      write(MDA_string,'(A)') './dir_tmp/svd_V_r_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,svd_V_r_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_delta_v.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_delta_v,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_gamma_z.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_gamma_z,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_r.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_r,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_r
      write(MDA_string,'(A)') './dir_tmp/n_w_.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_w_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_A.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_A,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_S.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_S,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_S
      write(MDA_string,'(A)') './dir_tmp/I_S_sample_.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,I_S_sample_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ld_S.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ld_S,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_M.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_M,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_M
      write(MDA_string,'(A)') './dir_tmp/I_M_sample_.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,I_M_sample_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ld_M.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ld_M,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_CTF.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_CTF,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ld_CTF.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ld_CTF,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_delta_v
      write(MDA_string,'(A)') './dir_tmp/delta_x_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_x_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_delta_v
      write(MDA_string,'(A)') './dir_tmp/delta_y_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_y_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_gamma_z
      write(MDA_string,'(A)') './dir_tmp/gamma_z_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,gamma_z_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_M
      write(MDA_string,'(A)') './dir_tmp/delta_x_est_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_x_est_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_M
      write(MDA_string,'(A)') './dir_tmp/delta_y_est_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_y_est_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_M
      write(MDA_string,'(A)') './dir_tmp/gamma_z_est_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,gamma_z_est_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_M
      write(MDA_string,'(A)') './dir_tmp/ctf_ind_est_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,ctf_ind_est_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_r
      write(MDA_string,'(A)') './dir_tmp/grid_p_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,grid_p_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = ld_S*(1+max_i4_f(n_S,I_S_sample_))
      write(MDA_string,'(A)') './dir_tmp/S_p__.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,S_p__,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = ld_M*(1+max_i4_f(n_M,I_M_sample_))
      write(MDA_string,'(A)') './dir_tmp/M_p__.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,M_p__,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = ld_CTF*n_CTF
      write(MDA_string,'(A)') './dir_tmp/CTF_p__.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,CTF_p__,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_w_max*n_svd_l*n_S*n_M
      write(MDA_string,'(A)') './dir_tmp/S_Z_T_R_CTF_M_q__.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,S_Z_T_R_CTF_M_q__,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_w_max*n_delta_v*n_S*n_M
      write(MDA_string,'(A)') './dir_tmp/S_T_T_R_CTF_M_q__.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,S_T_T_R_CTF_M_q__,MDA_string)
c$$$      %%%%%%%% ;
      ns_dump = 0
      nm_dump = 0
      ndelta_dump = 0
      ngamma_dump = 0
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ns_dump.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ns_dump,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/nm_dump.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,nm_dump,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ndelta_dump.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ndelta_dump,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/ngamma_dump.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,ngamma_dump,MDA_string)
c$$$      %%%%%%%% ;

      allocate(Z_S_svdd_(0:n_svd_l*n_delta_v - 1))
      call innerproduct_q_k_svdd_3(.true.,svd_d_max,n_svd_d,svd_d_
     $     ,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $     ,Z_S_svdd_)

c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_l*n_delta_v
      write(MDA_string,'(A)') './dir_tmp/Z_S_svdd_.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_svdd_,MDA_string)
c$$$      %%%%%%%% ;
      MDA_n_d = 1
      MDA_d_(0) = n_svd_l*n_delta_v
      write(MDA_string,'(A)') './dir_tmp/Z_M_svdd_.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,Z_M_svdd_,MDA_string)
c$$$      %%%%%%%% ;

      allocate(Z_S_svdr_(0:n_svd_l*n_w_max - 1))
      call cl1_c16(n_svd_l*n_w_max,Z_S_svdr_)
      allocate(Z_M_svdr_(0:n_svd_l*n_w_max - 1))
      call cl1_c16(n_svd_l*n_w_max,Z_M_svdr_)
      allocate(C_S_q_(0:n_w_max*n_delta_v-1))
      call cl1_c16(n_w_max*n_delta_v,C_S_q_)
      allocate(C_M_q_(0:n_w_max*n_delta_v-1))
      call cl1_c16(n_w_max*n_delta_v,C_M_q_)
      allocate(S_R_T_x_T_R_CTF_M_q__(0:n_w_max-1))
      call cl1_c16(n_w_max,S_R_T_x_T_R_CTF_M_q__)
      allocate(S_T_R_x_T_R_CTF_M_q__(0:n_w_max-1))
      call cl1_c16(n_w_max,S_T_R_x_T_R_CTF_M_q__)
      allocate(T_p_(0:n_A-1))
      call cl1_c16(n_A,T_p_)
      allocate(T_q_(0:n_A-1))
      call cl1_c16(n_A,T_q_)

      allocate(CTF_R_S_use_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,CTF_R_S_use_)
      allocate(ZZ1_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,ZZ1_)
      allocate(ZZ2_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,ZZ2_)
      allocate(Z_tmp_(0:n_w_max-1))
      call cl1_c16(n_w_max,Z_tmp_)
      do nx1=0,n_S*n_M-1
         nx = nx1
         ns = mod(nx,n_S)
         nx = (nx-ns)/n_S
         nm = mod(nx,n_M)
         nx = (nx-nm)/n_M
         if (nx.ne.0) then
            write(6,'(A,I0)') 'Warning! incorrectly unpacked nx1: '
     $           , nx1
         end if                 ! if (nx.ne.0) then
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') , ' nx1: ' , nx1 , ' ns: ' , ns
     $           ,' nm: ' , nm 
         end if                 !if (verbose.gt.2) then
c$$$         Define CTF_R_S_use_
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         call get_CTF_R_S_use_(gamma_z_est,n_gamma_z,CTF_R_S_(0 + ns
     $        *n_gamma_z + nctf*n_gamma_z*n_S),CTF_R_S_use_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_gamma_z
            write(MDA_string,'(A)') './dir_tmp/CTF_R_S_use_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,CTF_R_S_use_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         Define S_q
         call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/S_p_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,S_p_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $        ,fftw_in1_,fftw_out_,S_p_,S_q_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/S_q_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,S_q_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = 1
            write(MDA_string,'(A)') './dir_tmp/delta_x_est.mda'
            call MDA_write_r8(MDA_n_d,MDA_d_,delta_x_est,MDA_string)            
            write(MDA_string,'(A)') './dir_tmp/delta_y_est.mda'
            call MDA_write_r8(MDA_n_d,MDA_d_,delta_y_est,MDA_string)            
            write(MDA_string,'(A)') './dir_tmp/gamma_z_est.mda'
            call MDA_write_r8(MDA_n_d,MDA_d_,gamma_z_est,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         Define M_q
         nctf = nint(ctf_ind_est_(nm))
         call cp1_c16(n_A,M_p__(I_M_sample_(nm)*ld_M),M_p_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/M_p_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,M_p_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call xc1_c16(n_A,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/M_p_1_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,M_p_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call rotate_p2p_fx(n_r,fftw_plan_frwd_,fftw_plan_back_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,-gamma_z_est,M_p_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/M_p_2_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,M_p_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/M_p_3_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,M_p_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,M_q_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_A
            write(MDA_string,'(A)') './dir_tmp/M_q_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,M_q_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         Combine S_q and M_q into Z_S_svdr_
         call innerproduct_q_k_svdr_5(.true.,svd_r_max,n_svd_r,svd_r_
     $        ,n_svd_l,svd_l_,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_
     $        ,M_q_,Z_S_svdr_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_svd_l*n_w_max
            write(MDA_string,'(A)') './dir_tmp/Z_S_svdr_0_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_svdr_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         fftw each term in Z_S_svdr_
         do nsvd_l=0,n_svd_l-1
            call cps_c16(n_w_max,Z_S_svdr_(nsvd_l),n_svd_l
     $           ,fftw_out_last_,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(n_w_max,fftw_in1_last_,1,Z_S_svdr_(nsvd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
c$$$         Combine S_q and M_q into Z_M_svdr_
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_svd_l*n_w_max
            write(MDA_string,'(A)') './dir_tmp/Z_S_svdr_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_svdr_,MDA_string)
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         call innerproduct_q_k_svdr_5(.false.,svd_r_max,n_svd_r,svd_r_
     $        ,n_svd_l,svd_l_,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_
     $        ,M_q_,Z_M_svdr_)
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_svd_l*n_w_max
            write(MDA_string,'(A)') './dir_tmp/Z_M_svdr_0_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_M_svdr_,MDA_string)            
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         fftw each term in Z_M_svdr_
         do nsvd_l=0,n_svd_l-1
            call cps_c16(n_w_max,Z_M_svdr_(nsvd_l),n_svd_l
     $           ,fftw_out_last_,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(n_w_max,fftw_in1_last_,1,Z_M_svdr_(nsvd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_svd_l*n_w_max
            write(MDA_string,'(A)') './dir_tmp/Z_M_svdr_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,Z_M_svdr_,MDA_string)
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         do nsvd_l=0,n_svd_l-1
            C_M_tmp = (0.0d0,0.0d0)
            Z_M_tmp = (0.0d0,0.0d0)
            do nw=0,n_w_max-1
               C_M_tmp = C_M_tmp + zabs(Z_M_svdr_(nsvd_l + nw*n_svd_l))
     $              **2
               nx2 = nsvd_l + n_svd_l*(ns + n_S*nm)
               nx3 = nw*n_svd_l*n_S*n_M
               Z_M_tmp = Z_M_tmp + zabs(Z_M_svdr_(nsvd_l + nw*n_svd_l) -
     $              S_Z_T_R_CTF_M_q__(nx2+nx3))**2
            enddo               !do nw=0,n_w_max-1
            C_M_tmp = zsqrt(C_M_tmp)
            Z_M_tmp = zsqrt(Z_M_tmp)
            if (verbose.gt.-2) then
               write(6 ,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $              ,' Z_M_svdr_ vs S_Z_T_R_CTF_M_q__: ' ,' nx1: ' , nx1
     $              ,' nsvd_l: ' , nsvd_l , ' ns: ' , ns ,' nm: ' , nm ,
     $              ' absolute error: ' ,real(zabs(Z_M_tmp)) ,
     $              ' relative error: ' ,real(zabs(Z_M_tmp)
     $              /zabs(C_M_tmp))
            end if              !if (verbose.gt.2) then               
         enddo                  !do nsvd_l=0,n_svd_l-1
c$$$         Combining terms to recover delta_x and delta_y, ;
c$$$         but not yet interpolating gamma. ;
         nC1 = 0
         do ndv=0,n_delta_v-1
            do nw=0,n_w_max-1
c$$$  nC1 = nw + ndv*n_w_max
               C_S_q_(nC1) = cmplx( 0.0 , 0.0 )
               C_M_q_(nC1) = cmplx( 0.0 , 0.0 )
               nC2 = ndv*n_svd_l
               nC3 = nw*n_svd_l
               do nsvd_l=0,n_svd_l-1
                  C_S_q_(nC1) = C_S_q_(nC1) + Z_S_svdd_(nC2)
     $                 *Z_S_svdr_(nC3)
                  C_M_q_(nC1) = C_M_q_(nC1) + Z_M_svdd_(nC2)
     $                 *Z_M_svdr_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo            !do nsvd_l=0,n_svd_l-1
               nC1 = nC1 + 1
            enddo               !do nw=0,n_w_max-1
         enddo                  !do ndv=0,n_delta_v-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (verbose.gt.3) then
         do ndv=0,n_delta_v-1
            write(6,'(A,I0)') ' ndv:' , ndv 
            call print_sub_c16(n_w_max,C_S_q_(n_w_max*ndv),' C_S_q_: ')
            call print_sub_c16(n_w_max,C_M_q_(n_w_max*ndv),' C_M_q_: ')
         enddo !do ndv=0,n_delta_v-1
         end if !if (verbose.gt.3) then
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_w_max*n_delta_v
            write(MDA_string,'(A)') './dir_tmp/C_S_q_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,C_S_q_,MDA_string)
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
c$$$         %%%%%%%%
         if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_w_max*n_delta_v
            write(MDA_string,'(A)') './dir_tmp/C_M_q_.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,C_M_q_,MDA_string)
         end if !if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump)) then
c$$$         %%%%%%%%
         do ndv=0,n_delta_v-1
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            C_S_tmp = (0.0d0,0.0d0)
            Z_S_tmp = (0.0d0,0.0d0)
            C_M_tmp = (0.0d0,0.0d0)
            Z_M_tmp = (0.0d0,0.0d0)
            Z_X_tmp = (0.0d0,0.0d0)
            call cl1_c16(n_w_max,S_R_T_x_T_R_CTF_M_q__)
            call cl1_c16(n_w_max,S_T_R_x_T_R_CTF_M_q__)
            do nw=0,n_w_max-1
               gamma_z = (2.0d0*pi*nw)/(1.0d0*n_w_max)
               call cp1_c16(n_A,S_q_,T_q_)
               call transf_svd_q_to_q_5(svd_r_max,n_svd_r,svd_r_
     $              ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $              ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,T_q_,delta_x
     $              ,delta_y,T_q_)
               call rotate_q_to_q_0(n_r,n_w_,n_A,T_q_,gamma_z,T_q_)
               call innerproduct_p(n_r,grid_p_,n_w_,n_A,T_q_,M_q_
     $              ,S_T_R_x_T_R_CTF_M_q__(nw))
               call cp1_c16(n_A,S_q_,T_q_)
               call rotate_q_to_q_0(n_r,n_w_,n_A,T_q_,gamma_z,T_q_)
               call transf_svd_q_to_q_5(svd_r_max,n_svd_r,svd_r_
     $              ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $              ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,T_q_,delta_x
     $              ,delta_y,T_q_)
               call innerproduct_p(n_r,grid_p_,n_w_,n_A,T_q_,M_q_
     $              ,S_R_T_x_T_R_CTF_M_q__(nw))
               nC4 = nw + n_w_max*ndv
c$$$  C_S_tmp = C_S_tmp + zabs(C_S_q_(nC4))**2
c$$$  C_M_tmp = C_M_tmp + zabs(C_M_q_(nC4))**2
               C_S_tmp = C_S_tmp + zabs(S_T_R_x_T_R_CTF_M_q__(nw))**2
               C_M_tmp = C_M_tmp + zabs(S_R_T_x_T_R_CTF_M_q__(nw))**2
               nx2 = ndv + n_delta_v*(ns + n_S*nm)
               nx3 = nw*n_delta_v*n_S*n_M
               Z_S_tmp = Z_S_tmp + zabs(C_S_q_(nC4) -
     $              S_T_R_x_T_R_CTF_M_q__(nw))**2
               Z_M_tmp = Z_M_tmp + zabs(C_M_q_(nC4) -
     $              S_R_T_x_T_R_CTF_M_q__(nw))**2
               Z_X_tmp = Z_X_tmp + zabs(C_M_q_(nC4) -
     $              S_T_T_R_CTF_M_q__(nx2+nx3))**2
            enddo               !do nw=0,n_w_max-1
c$$$         %%%%%%%%
            if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump) .and.
     $           (ndv.eq.ndelta_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_w_max
            write(MDA_string,'(A)')
     $           './dir_tmp/S_R_T_x_T_R_CTF_M_q__.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,S_R_T_x_T_R_CTF_M_q__
     $           ,MDA_string)
         end if !if dump
c$$$         %%%%%%%%
c$$$         %%%%%%%%
            if ((ns.eq.ns_dump) .and. (nm.eq.nm_dump) .and.
     $           (ndv.eq.ndelta_dump)) then
            MDA_n_d = 1
            MDA_d_(0) = n_w_max
            write(MDA_string,'(A)')
     $           './dir_tmp/S_T_R_x_T_R_CTF_M_q__.mda'
            call MDA_write_c16(MDA_n_d,MDA_d_,S_T_R_x_T_R_CTF_M_q__
     $           ,MDA_string)
         end if !if dump
c$$$         %%%%%%%%

            C_S_tmp = zsqrt(C_S_tmp)
            Z_S_tmp = zsqrt(Z_S_tmp)
            C_M_tmp = zsqrt(C_M_tmp)
            Z_M_tmp = zsqrt(Z_M_tmp)
            Z_X_tmp = zsqrt(Z_X_tmp)
            if (verbose.gt.-2) then
               write(6
     $              ,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $              ,' C_S_q_ vs S_T_R_x_T_R_CTF_M_q__: ' ,' nx1: ' ,
     $              nx1 ,' ndv: ' , ndv , ' ns: ' ,
     $              ns ,' nm: ' , nm , ' absolute error: '
     $              ,real(zabs(Z_S_tmp)) , ' relative error: '
     $              ,real(zabs(Z_S_tmp) /zabs(C_S_tmp))
            end if              !if (verbose.gt.2) then               
            if (verbose.gt.-2) then
               write(6
     $              ,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $              ,' C_M_q_ vs S_R_T_x_T_R_CTF_M_q__: ' ,' nx1: ' ,
     $              nx1 ,' ndv: ' , ndv ,' ns: ' ,
     $              ns ,' nm: ' , nm , ' absolute error: '
     $              ,real(zabs(Z_M_tmp)) , ' relative error: '
     $              ,real(zabs(Z_M_tmp) /zabs(C_M_tmp))
            end if              !if (verbose.gt.2) then               
            if (verbose.gt.-2) then
               write(6
     $              ,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $              ,' C_M_q_ vs S_T_T_R_CTF_M_q__: ' ,' nx1: ' , nx1
     $              ,' ndv: ' , ndv ,' ns: ' , ns
     $              ,' nm: ' , nm , ' absolute error: '
     $              ,real(zabs(Z_X_tmp)) , ' relative error: '
     $              ,real(zabs(Z_X_tmp) /zabs(C_M_tmp))
            end if              !if (verbose.gt.2) then               
         enddo !do ndv=0,n_delta_v-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do ndv=0,n_delta_v-1
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call cp1_c16(n_w_max,C_M_q_(0+n_w_max*ndv),Z_tmp_)
            do ngz=0,n_gamma_z-1
               CTF_R_S_use = CTF_R_S_use_(ngz)
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $              ,Z_tmp)
c$$$  nC = ndv + n_delta_v*ngz
               nC = ngz
               if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
                  C_Z_use = 1.0d0
               else
                  C_Z_use = C_M*CTF_R_S_use
               end if
               ZZ1_(nC) = Z_tmp/(n_r**4)/C_Z_use
            enddo               !do ngz=0,n_gamma_z-1
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            do nw=0,n_w_max-1
               nx2 = ndv + n_delta_v*(ns + n_S*nm)
               nx3 = nw*n_delta_v*n_S*n_M
               Z_tmp_(nw) = S_T_T_R_CTF_M_q__(nx2 + nx3)
            enddo               !do nw=0,n_w_max-1
            do ngz=0,n_gamma_z-1
               CTF_R_S_use = CTF_R_S_use_(ngz)
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $              ,Z_tmp)
               nC = ngz
               if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
                  C_Z_use = 1.0d0
               else
                  C_Z_use = C_M*CTF_R_S_use
               end if
               ZZ2_(nC) = Z_tmp/(n_r**4)/C_Z_use
            enddo               !do ngz=0,n_gamma_z-1
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            C_tmp = (0.0d0,0.0d0)
            Z_tmp = (0.0d0,0.0d0)
            do ngz=0,n_gamma_z-1
               C_tmp = C_tmp + zabs(ZZ1_(ngz))**2
               Z_tmp = Z_tmp + zabs(ZZ1_(ngz) - ZZ2_(ngz))**2
            enddo               !do ngz=0,n_gamma_z-1
            C_tmp = zsqrt(C_tmp)
            Z_tmp = zsqrt(Z_tmp)
            if (verbose.gt.-2) then
               write(6 ,
     $              '(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A,F16.8)'
     $              )' ZZ1 vs ZZ2: ' ,' nx1: ' , nx1 ,' ndv: ' , ndv
     $              , ' ns: ' , ns ,' nm: ' , nm ,
     $              ' absolute error: ' ,real(zabs(Z_tmp)) ,
     $              ' relative error: ' ,real(zabs(Z_tmp)
     $              /zabs(C_tmp))
            end if              !if (verbose.gt.2) then            
c$$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         enddo                  !do ndv=0,n_delta_v-1

      enddo ! do nx1=0,n_S*n_M-1

      deallocate(S_R_T_x_T_R_CTF_M_q__)
      deallocate(S_T_R_x_T_R_CTF_M_q__)
      deallocate(T_p_)
      deallocate(T_q_)
      deallocate(Z_S_svdr_)
      deallocate(Z_M_svdr_)
      deallocate(C_S_q_)
      deallocate(C_M_q_)
      deallocate(CTF_R_S_use_)
      deallocate(ZZ1_)
      deallocate(ZZ2_)
      deallocate(Z_tmp_)
      deallocate(Z_S_svdd_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished ti7_check_SxZTRM_4]'
      end if !if (verbose.gt.0) then
      end
