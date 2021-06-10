      subroutine test_innerproduct_fast_SZxTRM_0(verbose,n_svd_r,svd_r_
     $     ,n_svd_l,svd_l_,svd_s_,svd_V_r_,Z_svdd_,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,delta_x_est_
     $     ,delta_y_est_ ,gamma_z_est_,ctf_ind_est_,fftw_plan_frwd_
     $     ,fftw_in1_ ,fftw_out_,fftw_plan_back_last,fftw_in1_last_
     $     ,fftw_out_last_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S,ld_S
     $     ,S_p__,M_p_,M_q_,n_M ,ld_M ,M_p__,CTF_p_,n_CTF ,ld_CTF
     $     ,CTF_p__,C_M_,CTF_R_S_ ,S_Z_T_R_CTF_M_q__,S_T_T_R_CTF_M_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      complex *16 Z_svdd_(0:0)
      integer n_delta_x,n_delta_y,n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      integer n_S,ld_S
      integer n_M,ld_M
      integer n_CTF,ld_CTF
      real *8 delta_x_(0:0)
      real *8 delta_y_(0:0)
      real *8 gamma_z_(0:0)
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
      complex *16 C_M_(0:0),C_M
      complex *16 CTF_R_S_(0:0)
      complex *16 S_Z_T_R_CTF_M_q__(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      complex *16, allocatable :: CTF_R_S_use_(:)
      complex *16, allocatable :: ZZ1_(:)
      complex *16, allocatable :: ZZ2_(:)
      complex *16, allocatable :: Z_svdr_(:)
      complex *16, allocatable :: C_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16 C_tmp,Z_tmp,C_Z_use,CTF_R_S_use
      integer n_w_max,nr,nw,ns,nm,nctf,nC
      integer ndx,ndy,ngz,nx,ld_X,nx1,nx2,nx3
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      integer nsvd_l,nC1,nC2,nC3,nC4
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_SZxTRM_0]'
      end if !if (verbose.gt.0) then

      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)

      allocate(Z_svdr_(0:n_svd_l*n_w_max - 1))
      call cl1_c16(n_svd_l*n_w_max,Z_svdr_)
      allocate(C_q_(0:n_w_max*n_delta_x*n_delta_y-1))
      call cl1_c16(n_w_max*n_delta_x*n_delta_y,C_q_)

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
c$$$         Define S_q
         call cp1_c16(n_A,S_p__(ns*ld_S),S_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $        ,fftw_in1_,fftw_out_,S_p_,S_q_)
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
c$$$         Define M_q
         nctf = nint(ctf_ind_est_(nm))
         call cp1_c16(n_A,M_p__(nm*ld_M),M_p_)
         call xc1_c16(n_A,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
         call rotate_p2p_fz(n_r,n_w_,n_A,M_p_,-gamma_z_est,M_p_)
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A
     $        ,fftw_in1_,fftw_out_,M_p_,M_q_)
c$$$         Combine S_q and M_q into Z_svdr_
         call innerproduct_q_k_svdr(n_svd_r,svd_r_,n_svd_l,svd_l_
     $        ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_,M_q_ ,Z_svdr_)
c$$$         fftw each term in Z_svdr_
         do nsvd_l=0,n_svd_l-1
            call cps_c16(n_w_max,Z_svdr_(nsvd_l),n_svd_l ,fftw_out_last_
     $           ,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(n_w_max,fftw_in1_last_,1,Z_svdr_(nsvd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
         do nsvd_l=0,n_svd_l-1
            C_tmp = (0.0d0,0.0d0)
            Z_tmp = (0.0d0,0.0d0)
            do nw=0,n_w_max-1
               C_tmp = C_tmp + zabs(Z_svdr_(nsvd_l + nw*n_svd_l))**2
               nx2 = nsvd_l + n_svd_l*(ns + n_S*nm)
               nx3 = nw*n_svd_l*n_S*n_M
               Z_tmp = Z_tmp + zabs(Z_svdr_(nsvd_l + nw*n_svd_l) -
     $              S_Z_T_R_CTF_M_q__(nx2+nx3))**2
            enddo               !do nw=0,n_w_max-1
            C_tmp = zsqrt(C_tmp)
            Z_tmp = zsqrt(Z_tmp)
            if (verbose.gt.-2) then
               write(6 ,'(A,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $              ,' Z_svdr_ vs S_Z_T_R_CTF_M_q__: ' ,' nx1: ' , nx1
     $              ,' nsvd_l: ' , nsvd_l , ' ns: ' , ns ,' nm: ' , nm ,
     $              ' absolute error: ' ,real(zabs(Z_tmp)) ,
     $              ' relative error: ' ,real(zabs(Z_tmp) /zabs(C_tmp))
            end if              !if (verbose.gt.2) then               
         enddo                  !do nsvd_l=0,n_svd_l-1
c$$$         Combining terms to recover delta_x and delta_y, ;
c$$$         but not yet interpolating gamma. ;
         nC1 = 0
         do ndy=0,n_delta_y-1
            do ndx=0,n_delta_x-1
               do nw=0,n_w_max-1
c$$$              nC1 = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
                  C_q_(nC1) = cmplx( 0.0 , 0.0 )
                  nC2 = ndx*n_svd_l + ndy*n_svd_l*n_delta_x
                  nC3 = nw*n_svd_l
                  do nsvd_l=0,n_svd_l-1
                     C_q_(nC1) = C_q_(nC1) + Z_svdd_(nC2)
     $                    *Z_svdr_(nC3)
                     nC2 = nC2 + 1
                     nC3 = nC3 + 1
                  enddo !do nsvd_l=0,n_svd_l-1
                  nC1 = nC1 + 1
               enddo            !do nw=0,n_w_max-1
            enddo               !ndx=0,n_delta_x-1
         enddo                  !do ndy=0,n_delta_y-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do ndx=0,n_delta_x-1
            do ndy=0,n_delta_y-1
               C_tmp = (0.0d0,0.0d0)
               Z_tmp = (0.0d0,0.0d0)
               do nw=0,n_w_max-1
                  nC4 = nw + n_w_max*(ndx + n_delta_x*ndy);
                  C_tmp = C_tmp + zabs(C_q_(nC4))**2
                  nx2 = ndx + n_delta_x*(ndy + n_delta_y*(ns + n_S*nm))
                  nx3 = nw*n_delta_x*n_delta_y*n_S*n_M
                  Z_tmp = Z_tmp + zabs(C_q_(nC4) -
     $                 S_T_T_R_CTF_M_q__(nx2+nx3))**2
               enddo !do nw=0,n_w_max-1
               C_tmp = zsqrt(C_tmp)
               Z_tmp = zsqrt(Z_tmp)
               if (verbose.gt.-2) then
                  write(6
     $                 ,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,F16.8,A ,F16.8)')
     $                 ,' C_q_ vs S_T_T_R_CTF_M_Q__: ' ,' nx1: ' , nx1
     $                 ,' ndx: ' , ndx , ' ndy: ' , ndy ,' ns: ' , ns
     $                 ,' nm: ' , nm , ' absolute error: '
     $                 ,real(zabs(Z_tmp)) , ' relative error: '
     $                 ,real(zabs(Z_tmp) /zabs(C_tmp))
               end if !if (verbose.gt.2) then               
            enddo !do ndy=0,n_delta_y-1
         enddo !do ndx=0,n_delta_x-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do ndx=0,n_delta_x-1
            do ndy=0,n_delta_y-1
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               call cp1_c16(n_w_max,C_q_(0+n_w_max*(ndx+n_delta_x *ndy))
     $              ,Z_tmp_)
               do ngz=0,n_gamma_z-1
                  CTF_R_S_use = CTF_R_S_use_(ngz)
                  gamma_z = gamma_z_(ngz)
                  call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $                 ,Z_tmp)
c$$$  nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
                  nC = ngz
                  if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
                     C_Z_use = 1.0d0
                  else
                     C_Z_use = C_M*CTF_R_S_use
                  end if
                  ZZ1_(nC) = Z_tmp/(n_r**4)/C_Z_use
               enddo            !do ngz=0,n_gamma_z-1
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               do nw=0,n_w_max-1
                  nx2 = ndx + n_delta_x*(ndy + n_delta_y*(ns + n_S*nm))
                  nx3 = nw*n_delta_x*n_delta_y*n_S*n_M
                  Z_tmp_(nw) = S_T_T_R_CTF_M_q__(nx2 + nx3)
               enddo            !do nw=0,n_w_max-1
               do ngz=0,n_gamma_z-1
                  CTF_R_S_use = CTF_R_S_use_(ngz)
                  gamma_z = gamma_z_(ngz)
                  call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $                 ,Z_tmp)
                  nC = ngz
                  if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
                     C_Z_use = 1.0d0
                  else
                     C_Z_use = C_M*CTF_R_S_use
                  end if
                  ZZ2_(nC) = Z_tmp/(n_r**4)/C_Z_use
               enddo            !do ngz=0,n_gamma_z-1
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               C_tmp = (0.0d0,0.0d0)
               Z_tmp = (0.0d0,0.0d0)
               do ngz=0,n_gamma_z-1
                  C_tmp = C_tmp + zabs(ZZ1_(ngz))**2
                  Z_tmp = Z_tmp + zabs(ZZ1_(ngz) - ZZ2_(ngz))**2
               enddo            !do ngz=0,n_gamma_z-1
               C_tmp = zsqrt(C_tmp)
               Z_tmp = zsqrt(Z_tmp)
               if (verbose.gt.-2) then
                  write(6 ,
     $                 '(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,F16.8,A,F16.8)'
     $                 )' ZZ1 vs ZZ2: ' ,' nx1: ' , nx1 ,' ndx: ' , ndx
     $                 ,  ' ndy: ' , ndy ,' ns: ' , ns ,' nm: ' , nm ,
     $                 ' absolute error: ' ,real(zabs(Z_tmp)) ,
     $                 ' relative error: ' ,real(zabs(Z_tmp)
     $                 /zabs(C_tmp))
               end if           !if (verbose.gt.2) then            
c$$$               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
            enddo               !do ndy=0,n_delta_y-1
         enddo                  !do ndx=0,n_delta_x-1

      enddo ! do nx1=0,n_S*n_M-1

      deallocate(Z_svdr_)
      deallocate(C_q_)
      deallocate(CTF_R_S_use_)
      deallocate(ZZ1_)
      deallocate(ZZ2_)
      deallocate(Z_tmp_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_fast_SZxTRM_0]'
      end if !if (verbose.gt.0) then
      end
