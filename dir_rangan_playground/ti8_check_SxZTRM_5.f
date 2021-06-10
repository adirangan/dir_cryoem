!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Tests S_T_T_R_CTF_M_q_ and S_Z_T_R_CTF_M_q_ calculation directly. ;\n
      subroutine ti8_check_SxZTRM_5(
     $     verbose
     $     ,svd_d_max
     $     ,n_svd_d
     $     ,svd_d_
     $     ,svd_r_max
     $     ,n_svd_r
     $     ,svd_r_
     $     ,n_svd_l
     $     ,svd_l_
     $     ,svd_U_d_
     $     ,svd_s_
     $     ,svd_V_r_ 
     $     ,Z_M_svdd_
     $     ,n_delta_v
     $     ,delta_x_ 
     $     ,delta_y_ 
     $     ,n_gamma_z
     $     ,gamma_z_ 
     $     ,delta_x_est_
     $     ,delta_y_est_ 
     $     ,gamma_z_est_
     $     ,ctf_ind_est_ 
     $     ,fftw_plan_frwd_
     $     ,fftw_plan_back_
     $     ,fftw_0in_
     $     ,fftw_out_ 
     $     ,fftw_plan_frwd_last
     $     ,fftw_plan_back_last
     $     ,fftw_0in_last_
     $     ,fftw_out_last_
     $     ,quad_n_r 
     $     ,quad_grid_k_p_r_
     $     ,quad_weight_k_p_r_
     $     ,quad_n_w_
     $     ,quad_n_w_sum
     $     ,quad_S_k_p_
     $     ,quad_S_k_q_
     $     ,n_S
     $     ,I_S_sample_
     $     ,ld_S 
     $     ,quad_S_k_p__
     $     ,quad_M_k_p_ 
     $     ,quad_M_k_q_
     $     ,n_M 
     $     ,I_M_sample_
     $     ,ld_M
     $     ,quad_M_k_p__
     $     ,quad_CTF_k_p_ 
     $     ,n_CTF 
     $     ,ld_CTF
     $     ,quad_CTF_k_p__
     $     ,C_M_
     $     ,CTF_R_S_
     $     ,S_Z_T_R_CTF_M_q__ 
     $     ,S_T_T_R_CTF_M_q__
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_ and S_Z_T_R_CTF_M_q_.
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
      integer quad_n_r,quad_n_w_(0:quad_n_r-1),quad_n_w_sum
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
      real *8 quad_grid_k_p_r_(0:quad_n_r-1)
      real *8 quad_weight_k_p_r_(0:quad_n_r-1)
      logical flag_nufft
      integer *8 fftw_plan_frwd_(0:quad_n_r-1)
      integer *8 fftw_plan_back_(0:quad_n_r-1)
      integer *8 fftw_plan_frwd_last
      integer *8 fftw_plan_back_last
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      complex *16 fftw_0in_last_(0:0),fftw_out_last_(0:0)
      complex *16 quad_S_k_p_(0:0),quad_S_k_q_(0:0),quad_S_k_p__(0:0)
      complex *16 quad_M_k_p_(0:0),quad_M_k_q_(0:0),quad_M_k_p__(0:0)
      complex *16 quad_CTF_k_p_(0:0),CTF_q_(0:0),quad_CTF_k_p__(0:0)
      complex *16 C_M_(0:0),C_M
      complex *16 CTF_R_S_(0:0)
      complex *16 S_Z_T_R_CTF_M_q__(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      real *8 svd_d_max,svd_r_max
      complex *16, allocatable :: CTF_R_S_use_(:)
      complex *16, allocatable :: Z_S_svdr_(:)
      complex *16, allocatable :: Z_M_svdr_(:)
      complex *16, allocatable :: tmp_Z_M_svdr_(:)
      complex *16, allocatable :: quad_SX1_(:)
      complex *16, allocatable :: quad_SX2_(:)
      complex *16, allocatable :: quad_MX1_(:)
      complex *16, allocatable :: quad_MX2_(:)
      complex *16, allocatable :: S_R_T_x_T_R_CTF_M_q_svd__(:)
      complex *16, allocatable :: S_T_R_x_T_R_CTF_M_q_svd__(:)
      complex *16, allocatable :: S_R_T_x_T_R_CTF_M_q_ori__(:)
      complex *16, allocatable :: S_T_R_x_T_R_CTF_M_q_ori__(:)
      complex *16, allocatable :: tmp_S_T_T_R_CTF_M_q__(:)
      complex *16, allocatable :: quad_T_k_p_(:)
      complex *16, allocatable :: quad_T_k_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      complex *16 C_tmp,Z_tmp,C_Z_use,CTF_R_S_use
      complex *16 C_S_tmp,C_M_tmp,Z_S_tmp,Z_M_tmp,Z_X_tmp
      complex *16 costheta_c16_f
      complex *16 c16_alpha,c16_beta
      integer quad_n_w_max,nr,nw,ns,nm,nctf,nC
      integer ndv,ngz,nx,ld_X,nx1,nx2,nx3,tab
      real *8 delta_x,delta_y,gamma_z
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 tmp_gamma_z
      integer nsvd_l,nC1,nC2,nC3,nC4
      real *8 pi
      integer max_i4_f !function output. ;
      real *8 al2_c16_f, fnod_c16_f !function output. ;
      integer *4 MDA_n_d,MDA_d_(0:127)
      character(len=1024) MDA_string
      logical flag_memory_checkset
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_check_SxZTRM_5]'
      end if !if (verbose.gt.0) then

      pi = 4.0d0*datan(1.0d0)
      quad_n_w_max = quad_n_w_(quad_n_r-1)
      flag_nufft = .false.
      if (n_gamma_z.ne.quad_n_w_max) then
         flag_nufft = .true.
      end if !if (n_gamma_z.ne.quad_n_w_max) then
      if (n_gamma_z.eq.quad_n_w_max) then
         do nw=0,quad_n_w_max-1
            tmp_gamma_z = (2.0d0*pi*nw)/(1.0d0*n_gamma_z)
            ngz = nw
            gamma_z = gamma_z_(ngz)
            if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
               flag_nufft = .true.
            end if !if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
         enddo !do nw=0,quad_n_w_max-1
      end if !if (n_gamma_z.eq.quad_n_w_max) then
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' quad_n_w_max: ' , quad_n_w_max
         write(6,'(A,L2)') ' flag_nufft: ' , flag_nufft
      end if !if (verbose.gt.1) then      

      allocate(Z_S_svdd_(0:1+n_svd_l*n_delta_v - 1))
      call cs1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
      call innerproduct_q_k_svdd_polyval_0on_0(.true.,svd_d_max,n_svd_d
     $     ,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $     ,Z_S_svdd_)
      allocate(Z_S_svdr_(0:1+n_svd_l*quad_n_w_max - 1))
      call cs1_c16(n_svd_l*quad_n_w_max,Z_S_svdr_)
      allocate(Z_M_svdr_(0:1+n_svd_l*quad_n_w_max - 1))
      call cs1_c16(n_svd_l*quad_n_w_max,Z_M_svdr_)
      allocate(tmp_Z_M_svdr_(0:1+n_svd_l*quad_n_w_max - 1))
      call cs1_c16(n_svd_l*quad_n_w_max,tmp_Z_M_svdr_)
      allocate(quad_SX1_(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,quad_SX1_)
      allocate(quad_SX2_(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,quad_SX2_)
      allocate(quad_MX1_(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,quad_MX1_)
      allocate(quad_MX2_(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,quad_MX2_)
      allocate(S_R_T_x_T_R_CTF_M_q_svd__(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_svd__)
      allocate(S_T_R_x_T_R_CTF_M_q_svd__(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_svd__)
      allocate(S_R_T_x_T_R_CTF_M_q_ori__(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_ori__)
      allocate(S_T_R_x_T_R_CTF_M_q_ori__(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_ori__)
      allocate(tmp_S_T_T_R_CTF_M_q__(0:1+quad_n_w_max*n_delta_v-1))
      call cs1_c16(quad_n_w_max*n_delta_v,tmp_S_T_T_R_CTF_M_q__)
      allocate(quad_T_k_p_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_T_k_p_)
      allocate(quad_T_k_q_(0:1+quad_n_w_sum-1))
      call cs1_c16(quad_n_w_sum,quad_T_k_q_)
      allocate(CTF_R_S_use_(0:1+n_gamma_z-1))
      call cs1_c16(n_gamma_z,CTF_R_S_use_)
      allocate(Z_tmp_(0:1+quad_n_w_max-1))
      call cs1_c16(quad_n_w_max,Z_tmp_)
      do nx1=0,n_S*n_M-1
c$$$         %%%%%%%%
c$$$         % unpack nx1
c$$$         %%%%%%%%
         nx = nx1
         ns = mod(nx,n_S)
         nx = (nx-ns)/n_S
         nm = mod(nx,n_M)
         nx = (nx-nm)/n_M
         if (nx.ne.0) then
            write(6,'(A,I0)') 'Warning, incorrectly unpacked nx1: '
     $           , nx1
         end if                 ! if (nx.ne.0) then
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') , ' nx1: ' , nx1 , ' ns: ' , ns
     $           ,' nm: ' , nm 
         end if                 !if (verbose.gt.2) then
c$$$         %%%%%%%%
c$$$         % Define CTF_R_S_use_
c$$$         %%%%%%%%
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         if (flag_nufft.eqv..false.) then
            call rotate_p_to_p_fftw_single(
     $           fftw_plan_frwd_last
     $           ,fftw_plan_back_last
     $           ,quad_n_w_max
     $           ,fftw_0in_last_
     $           ,fftw_out_last_
     $           ,CTF_R_S_(0+n_gamma_z*(nctf + n_CTF*ns))
     $           ,-gamma_z_est
     $           ,CTF_R_S_use_
     $           )
         end if !if (flag_nufft.eqv..false.) then
         if (flag_nufft.eqv..true.) then
            call interp1_lagrange_c16(
     $           n_gamma_z
     $           ,0.0d0
     $           ,2.0d0*pi
     $           ,CTF_R_S_(0+n_gamma_z*(nctf + n_CTF*ns))
     $           ,0
     $           ,-gamma_z_est
     $           ,CTF_R_S_use_
     $           )
         end if !if (flag_nufft.eqv..true.) then
         if (verbose.gt.1) then
            call print_sub_c16(n_gamma_z,CTF_R_S_use_,' CTF_R_S_use_: ')
         end if !if (verbose.gt.1) then
c$$$         %%%%%%%%
c$$$         % Define S_k_q
c$$$         %%%%%%%%
         call cp1_c16(quad_n_w_sum,quad_S_k_p__(I_S_sample_(ns)*ld_S)
     $        ,quad_S_k_p_)
         call interp_p_to_q_fftw(quad_n_r,fftw_plan_frwd_,quad_n_w_
     $        ,quad_n_w_sum,fftw_0in_,fftw_out_,quad_S_k_p_,quad_S_k_q_)
c$$$         %%%%%%%%
c$$$         % extract parameters
c$$$         %%%%%%%%
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
c$$$         %%%%%%%%
c$$$         % Define M_q
c$$$         %%%%%%%%
         nctf = nint(ctf_ind_est_(nm))
         call cp1_c16(quad_n_w_sum,quad_M_k_p__(I_M_sample_(nm)*ld_M)
     $        ,quad_M_k_p_)
         call xc1_c16(quad_n_w_sum,quad_M_k_p_,quad_CTF_k_p__(nctf
     $        *ld_CTF),quad_M_k_p_)
         call rotate_p_to_p_fftw(quad_n_r,fftw_plan_frwd_
     $        ,fftw_plan_back_,quad_n_w_,quad_n_w_sum,fftw_0in_
     $        ,fftw_out_ ,quad_M_k_p_,-gamma_z_est,quad_M_k_p_)
         call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $        ,quad_n_w_sum,quad_M_k_p_,-delta_x_est,-delta_y_est
     $        ,quad_M_k_p_)
         call interp_p_to_q_fftw(quad_n_r,fftw_plan_frwd_,quad_n_w_
     $        ,quad_n_w_sum,fftw_0in_,fftw_out_,quad_M_k_p_,quad_M_k_q_)
c$$$         %%%%%%%%
c$$$         % Combine S_q and M_q into Z_S_svdr_
c$$$         %%%%%%%%
         call innerproduct_q_k_svdr_quad_polyval_0on_0(.true. ,svd_r_max
     $        ,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_ ,svd_V_r_,quad_n_r
     $        ,quad_grid_k_p_r_,quad_weight_k_p_r_,quad_n_w_
     $        ,quad_n_w_sum ,quad_S_k_q_,quad_M_k_q_ ,Z_S_svdr_)
         c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
         c16_beta = dcmplx(0.0d0,0.0d0); 
         call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $        ,Z_S_svdr_,Z_S_svdr_)
c$$$         %%%%%%%%
c$$$         % fftw each term in Z_S_svdr_
c$$$         %%%%%%%%
         do nsvd_l=0,n_svd_l-1
            call cps_c16(quad_n_w_max,Z_S_svdr_(nsvd_l),n_svd_l
     $           ,fftw_out_last_,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(quad_n_w_max,fftw_0in_last_,1,Z_S_svdr_(nsvd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
c$$$         %%%%%%%%
c$$$         % Combine S_q and M_q into Z_M_svdr_
c$$$         %%%%%%%%
         call innerproduct_q_k_svdr_quad_polyval_0on_0(.false.
     $        ,svd_r_max,n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_ ,svd_V_r_
     $        ,quad_n_r,quad_grid_k_p_r_,quad_weight_k_p_r_,quad_n_w_
     $        ,quad_n_w_sum,quad_S_k_q_,quad_M_k_q_ ,Z_M_svdr_)
         c16_alpha = dcmplx(1.0d0/2.0d0,0.0d0); 
         c16_beta = dcmplx(0.0d0,0.0d0); 
         call af1_c16(n_svd_l*quad_n_w_max,c16_alpha,c16_beta
     $        ,Z_M_svdr_,Z_M_svdr_)
c$$$         %%%%%%%%
c$$$         fftw each term in Z_M_svdr_
c$$$         %%%%%%%%
         do nsvd_l=0,n_svd_l-1
            call cps_c16(quad_n_w_max,Z_M_svdr_(nsvd_l),n_svd_l
     $           ,fftw_out_last_,1)
            call dfftw_execute_(fftw_plan_back_last)
            call cps_c16(quad_n_w_max,fftw_0in_last_,1,Z_M_svdr_(nsvd_l)
     $           ,n_svd_l)
         enddo                  !do svd_l=0,n_svd_l-1
c$$$         %%%%%%%%
         call cl1_c16(quad_n_w_max*n_svd_l,tmp_Z_M_svdr_)
         do nsvd_l=0,n_svd_l-1
            do nw=0,quad_n_w_max-1
               nx2 = nsvd_l + n_svd_l*(ns + n_S*nm)
               nx3 = nw*n_svd_l*n_S*n_M
               tmp_Z_M_svdr_(nsvd_l + nw*n_svd_l) =
     $              0.50d0*S_Z_T_R_CTF_M_q__(nx2+nx3)
            enddo               !do nw=0,quad_n_w_max-1
         enddo                  !do nsvd_l=0,n_svd_l-1
         write(6,'(A,F24.16)') ' Z_M_svdr_ vs S_Z_T_R_CTF_M_q_: ' ,
     $        fnod_c16_f(quad_n_w_max*n_svd_l,Z_M_svdr_ ,tmp_Z_M_svdr_)
     $        /al2_c16_f(quad_n_w_max*n_svd_l,Z_M_svdr_)
c$$$         Combining terms to recover delta_x and delta_y, ;
c$$$         but not yet interpolating gamma. ;
         call zgemm(
     $        'T','N'
     $        ,n_delta_v
     $        ,quad_n_w_max
     $        ,n_svd_l
     $        ,1.0d0*dcmplx(1.0d0,0.0d0)
     $        ,Z_S_svdd_
     $        ,n_svd_l
     $        ,Z_S_svdr_
     $        ,n_svd_l
     $        ,0.0d0*dcmplx(1.0d0,0.0d0)
     $        ,quad_SX1_
     $        ,n_delta_v
     $        )
         call trn0_c16(n_delta_v,quad_n_w_max,quad_SX1_,quad_SX1_)
         call zgemm('T','N'
     $        ,n_delta_v
     $        ,quad_n_w_max
     $        ,n_svd_l
     $        ,1.0d0*dcmplx(1.0d0,0.0d0)
     $        ,Z_M_svdd_
     $        ,n_svd_l
     $        ,Z_M_svdr_
     $        ,n_svd_l
     $        ,0.0d0*dcmplx(1.0d0,0.0d0)
     $        ,quad_MX1_
     $        ,n_delta_v
     $        )
         call trn0_c16(n_delta_v,quad_n_w_max,quad_MX1_,quad_MX1_)
         nC1 = 0
         do ndv=0,n_delta_v-1
            do nw=0,quad_n_w_max-1
c$$$  nC1 = nw + ndv*quad_n_w_max
               quad_SX2_(nC1) = dcmplx(0.0d0,0.0d0)
               quad_MX2_(nC1) = dcmplx(0.0d0,0.0d0)
               nC2 = ndv*n_svd_l
               nC3 = nw*n_svd_l
               do nsvd_l=0,n_svd_l-1
                  quad_SX2_(nC1) = quad_SX2_(nC1) + Z_S_svdd_(nC2)
     $                 *Z_S_svdr_(nC3)
                  quad_MX2_(nC1) = quad_MX2_(nC1) + Z_M_svdd_(nC2)
     $                 *Z_M_svdr_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo            !do nsvd_l=0,n_svd_l-1
               nC1 = nC1 + 1
            enddo               !do nw=0,quad_n_w_max-1
         enddo                  !do ndv=0,n_delta_v-1
         write(6,'(A,F24.16)') ' quad_SX1_ vs quad_SX2_: ' ,
     $        fnod_c16_f(quad_n_w_max*n_delta_v,quad_SX1_,quad_SX2_)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_SX2_)
         write(6,'(A,F24.16)') ' quad_MX1_ vs quad_MX2_: ' ,
     $        fnod_c16_f(quad_n_w_max*n_delta_v,quad_MX1_,quad_MX2_)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_MX2_)
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (verbose.gt.3) then
         do ndv=0,n_delta_v-1
            write(6,'(A,I0)') ' ndv:' , ndv 
            call print_sub_c16(quad_n_w_max,quad_SX2_(quad_n_w_max*ndv)
     $           ,' quad_SX2_: ')
            call print_sub_c16(quad_n_w_max,quad_MX2_(quad_n_w_max*ndv)
     $           ,' quad_MX2_: ')
         enddo !do ndv=0,n_delta_v-1
         end if !if (verbose.gt.3) then
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         call cl1_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_svd__)
         call cl1_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_svd__)
         call cl1_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_ori__)
         call cl1_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_ori__)
         call cl1_c16(quad_n_w_max*n_delta_v,tmp_S_T_T_R_CTF_M_q__)
         do ndv=0,n_delta_v-1
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            Z_X_tmp = (0.0d0,0.0d0)
            do nw=0,quad_n_w_max-1
               tab = nw + ndv*quad_n_w_max
               gamma_z = (2.0d0*pi*nw)/(1.0d0*quad_n_w_max)
c$$$           %%%%%%%%
               call cp1_c16(quad_n_w_sum,quad_S_k_q_,quad_T_k_q_)
               call transf_svd_q_to_q_polyval_0on_0(svd_r_max,n_svd_r
     $              ,svd_r_,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_
     $              ,svd_U_d_,svd_s_,svd_V_r_,quad_n_r,quad_grid_k_p_r_
     $              ,quad_n_w_,quad_n_w_sum,quad_T_k_q_ ,delta_x,delta_y
     $              ,quad_T_k_q_)
               call rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum
     $              ,quad_T_k_q_,gamma_z,quad_T_k_q_)
               call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $              ,quad_weight_k_p_r_,quad_n_w_ ,quad_n_w_sum
     $              ,quad_T_k_q_ ,quad_M_k_q_
     $              ,S_T_R_x_T_R_CTF_M_q_svd__(tab))
               S_T_R_x_T_R_CTF_M_q_svd__(tab) = 0.5d0
     $              *S_T_R_x_T_R_CTF_M_q_svd__(tab)
c$$$           %%%%%%%%
               call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $              ,quad_n_w_sum,quad_S_k_p_,delta_x,delta_y
     $              ,quad_T_k_p_)
               call interp_p_to_q_fftw(quad_n_r,fftw_plan_frwd_
     $              ,quad_n_w_,quad_n_w_sum,fftw_0in_,fftw_out_
     $              ,quad_T_k_p_ ,quad_T_k_q_)
               call rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum
     $              ,quad_T_k_q_,gamma_z,quad_T_k_q_)
               call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $              ,quad_weight_k_p_r_,quad_n_w_ ,quad_n_w_sum
     $              ,quad_T_k_q_ ,quad_M_k_q_
     $              ,S_T_R_x_T_R_CTF_M_q_ori__(tab))
               S_T_R_x_T_R_CTF_M_q_ori__(tab) = 0.5d0
     $              *S_T_R_x_T_R_CTF_M_q_ori__(tab)
c$$$           %%%%%%%%
               call cp1_c16(quad_n_w_sum,quad_S_k_q_,quad_T_k_q_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_q_
     $              ,' 1 quad_T_k_q_: ')
               call rotate_q_to_q(quad_n_r,quad_n_w_,quad_n_w_sum
     $              ,quad_T_k_q_,gamma_z,quad_T_k_q_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_q_
     $              ,' 1 quad_T_k_q_: ')
               call transf_svd_q_to_q_polyval_0on_0(svd_r_max,n_svd_r
     $              ,svd_r_,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_
     $              ,svd_U_d_,svd_s_,svd_V_r_,quad_n_r,quad_grid_k_p_r_
     $              ,quad_n_w_,quad_n_w_sum,quad_T_k_q_ ,delta_x,delta_y
     $              ,quad_T_k_q_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_q_
     $              ,' 1 quad_T_k_q_: ')
               call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $              ,quad_weight_k_p_r_,quad_n_w_ ,quad_n_w_sum
     $              ,quad_T_k_q_ ,quad_M_k_q_
     $              ,S_R_T_x_T_R_CTF_M_q_svd__(tab))
               S_R_T_x_T_R_CTF_M_q_svd__(tab) = 0.5d0
     $              *S_R_T_x_T_R_CTF_M_q_svd__(tab)
c$$$           %%%%%%%%
               call print_sub_c16(quad_n_w_sum,quad_S_k_p_
     $              ,' 2 quad_S_k_p_: ')
               call rotate_p_to_p_fftw(quad_n_r,fftw_plan_frwd_
     $              ,fftw_plan_back_,quad_n_w_,quad_n_w_sum,fftw_0in_
     $              ,fftw_out_,quad_S_k_p_,gamma_z,quad_T_k_p_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_p_
     $              ,' 2 quad_T_k_p_: ')
               call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $              ,quad_n_w_sum,quad_T_k_p_,delta_x,delta_y
     $              ,quad_T_k_p_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_p_
     $              ,' 2 quad_T_k_p_: ')
               call interp_p_to_q_fftw(quad_n_r,fftw_plan_frwd_
     $              ,quad_n_w_,quad_n_w_sum,fftw_0in_,fftw_out_
     $              ,quad_T_k_p_ ,quad_T_k_q_)
               call print_sub_c16(quad_n_w_sum,quad_T_k_q_
     $              ,' 2 quad_T_k_q_: ')
               call print_sub_c16(quad_n_w_sum,quad_M_k_q_
     $              ,' 2 quad_M_k_q_: ')
               call innerproduct_p_quad(quad_n_r,quad_grid_k_p_r_
     $              ,quad_weight_k_p_r_,quad_n_w_ ,quad_n_w_sum
     $              ,quad_T_k_q_ ,quad_M_k_q_
     $              ,S_R_T_x_T_R_CTF_M_q_ori__(tab))
               S_R_T_x_T_R_CTF_M_q_ori__(tab) = 0.5d0
     $              *S_R_T_x_T_R_CTF_M_q_ori__(tab)
c$$$           %%%%%%%%
               nx2 = ndv + n_delta_v*(ns + n_S*nm)
               nx3 = nw*n_delta_v*n_S*n_M
               Z_X_tmp = S_T_T_R_CTF_M_q__(nx2+nx3)
               tmp_S_T_T_R_CTF_M_q__(tab) = Z_X_tmp
c$$$           %%%%%%%%
               if (verbose.gt.-2) then
                  write(6,'(A,I0,2(A,F24.16),A,I0,A,F24.16)')
     $                 ' ndv: ' , ndv 
     $                 ,' delta_x: ' , delta_x
     $                 ,' delta_y: ' , delta_y
     $                 , ' nw: ' , nw
     $                 , ' gamma_z: ' , gamma_z
                  write(6,'(2(A,2F24.16),A,2F24.16)') ' quad_SX2_ ' ,
     $                 quad_SX2_(tab) ,' S_T_R_x_T_R_CTF_M_q_svd__ ' ,
     $                 S_T_R_x_T_R_CTF_M_q_svd__(tab) , ' error ' ,
     $                 zabs(quad_SX2_(tab) -
     $                 S_T_R_x_T_R_CTF_M_q_svd__(tab))
     $                 /zabs(quad_SX2_(tab))
                  write(6,'(2(A,2F24.16),A,2F24.16)') ' quad_SX2_ ' ,
     $                 quad_SX2_(tab) ,' S_T_R_x_T_R_CTF_M_q_ori__ ' ,
     $                 S_T_R_x_T_R_CTF_M_q_ori__(tab) , ' error ' ,
     $                 zabs(quad_SX2_(tab) -
     $                 S_T_R_x_T_R_CTF_M_q_ori__(tab))
     $                 /zabs(quad_SX2_(tab))
                  write(6,'(2(A,2F24.16),A,2F24.16)') ' quad_MX2_ ' ,
     $                 quad_MX2_(tab) ,' S_R_T_x_T_R_CTF_M_q_svd__ ' ,
     $                 S_R_T_x_T_R_CTF_M_q_svd__(tab) , ' error ' ,
     $                 zabs(quad_MX2_(tab) -
     $                 S_R_T_x_T_R_CTF_M_q_svd__(tab))
     $                 /zabs(quad_MX2_(tab))
                  write(6,'(2(A,2F24.16),A,2F24.16)') ' quad_MX2_ ' ,
     $                 quad_MX2_(tab) ,' S_R_T_x_T_R_CTF_M_q_ori__ ' ,
     $                 S_R_T_x_T_R_CTF_M_q_ori__(tab) , ' error ' ,
     $                 zabs(quad_MX2_(tab) -
     $                 S_R_T_x_T_R_CTF_M_q_ori__(tab))
     $                 /zabs(quad_MX2_(tab))
               end if !if (verbose.gt.-2) then

      if ((ndv.eq.0) .and. (nw.eq.84)) then
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/n_delta_v.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,n_delta_v,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = n_delta_v
      write(MDA_string,'(A)') './dir_tmp/delta_x_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_x_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = n_delta_v
      write(MDA_string,'(A)') './dir_tmp/delta_y_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_y_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/quad_n_r.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,quad_n_r,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = quad_n_r
      write(MDA_string,'(A)') './dir_tmp/quad_grid_k_p_r_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,quad_grid_k_p_r_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = quad_n_r
      write(MDA_string,'(A)') './dir_tmp/quad_weight_k_p_r_.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,quad_weight_k_p_r_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = quad_n_r
      write(MDA_string,'(A)') './dir_tmp/quad_n_w_.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,quad_n_w_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/quad_n_w_sum.mda'
      call MDA_write_i4(MDA_n_d,MDA_d_,quad_n_w_sum,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/gamma_z.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,gamma_z,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/delta_x.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_x,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = 1
      write(MDA_string,'(A)') './dir_tmp/delta_y.mda'
      call MDA_write_r8(MDA_n_d,MDA_d_,delta_y,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = quad_n_w_sum
      write(MDA_string,'(A)') './dir_tmp/quad_S_k_p_.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,quad_S_k_p_,MDA_string)
      MDA_n_d = 1 ; MDA_d_(0) = quad_n_w_sum
      write(MDA_string,'(A)') './dir_tmp/quad_M_k_q_.mda'
      call MDA_write_c16(MDA_n_d,MDA_d_,quad_M_k_q_,MDA_string)
      stop !check for errors. ;
      end if !if ((ndv.eq.0) .and. (nw.eq.84)) then

c$$$           %%%%%%%%               
            enddo               !do nw=0,quad_n_w_max-1
         enddo !do ndv=0,n_delta_v-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         write(6,'(A,A,F24.16)') ' quad_SX2_ vs ' ,
     $        'S_T_R_x_T_R_CTF_M_q_svd__: ' , fnod_c16_f(quad_n_w_max
     $        *n_delta_v,quad_SX2_,S_T_R_x_T_R_CTF_M_q_svd__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_SX2_)
         write(6,'(A,A,F24.16)') ' quad_MX2_ vs ' ,
     $        'S_R_T_x_T_R_CTF_M_q_svd__: ' , fnod_c16_f(quad_n_w_max
     $        *n_delta_v,quad_MX2_,S_R_T_x_T_R_CTF_M_q_svd__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_MX2_)
         write(6,'(A,A,F24.16)') ' quad_SX2_ vs ' ,
     $        'S_T_R_x_T_R_CTF_M_q_ori__: ' , fnod_c16_f(quad_n_w_max
     $        *n_delta_v,quad_SX2_,S_T_R_x_T_R_CTF_M_q_ori__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_SX2_)
         write(6,'(A,A,F24.16)') ' quad_MX2_ vs ' ,
     $        'S_R_T_x_T_R_CTF_M_q_ori__: ' , fnod_c16_f(quad_n_w_max
     $        *n_delta_v,quad_MX2_,S_R_T_x_T_R_CTF_M_q_ori__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_MX2_)
         write(6,'(A,A,F24.16)') ' quad_MX2_ vs ' ,
     $        'tmp_S_T_T_R_CTF_M_q__: ' , fnod_c16_f(quad_n_w_max
     $        *n_delta_v,quad_MX2_,tmp_S_T_T_R_CTF_M_q__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v,quad_MX2_)

      enddo ! do nx1=0,n_S*n_M-1

      flag_memory_checkset = .true.
      call cxs_c16(n_svd_l*n_delta_v,Z_S_svdd_,'Z_S_svdd_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*quad_n_w_max,Z_S_svdr_,'Z_S_svdr_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*quad_n_w_max,Z_M_svdr_,'Z_M_svdr_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_svd_l*quad_n_w_max,tmp_Z_M_svdr_,'tmp_Z_M_svdr_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,quad_SX1_,'quad_SX1_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,quad_SX2_,'quad_SX2_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,quad_MX1_,'quad_MX1_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,quad_MX2_,'quad_MX2_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_svd__
     $     ,'S_R_T_x_T_R_CTF_M_q_svd__',flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_svd__
     $     ,'S_T_R_x_T_R_CTF_M_q_svd__',flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,S_R_T_x_T_R_CTF_M_q_ori__
     $     ,'S_R_T_x_T_R_CTF_M_q_ori__',flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,S_T_R_x_T_R_CTF_M_q_ori__
     $     ,'S_T_R_x_T_R_CTF_M_q_ori__',flag_memory_checkset)
      call cxs_c16(quad_n_w_max*n_delta_v,tmp_S_T_T_R_CTF_M_q__
     $     ,'tmp_S_T_T_R_CTF_M_q__',flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_T_k_p_,'quad_T_k_p_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_sum,quad_T_k_q_,'quad_T_k_q_'
     $     ,flag_memory_checkset)
      call cxs_c16(n_gamma_z,CTF_R_S_use_,'CTF_R_S_use_'
     $     ,flag_memory_checkset)
      call cxs_c16(quad_n_w_max,Z_tmp_,'Z_tmp_',flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.2) then
         write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.2) then
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then

      deallocate(S_R_T_x_T_R_CTF_M_q_svd__)
      deallocate(S_T_R_x_T_R_CTF_M_q_svd__)
      deallocate(S_R_T_x_T_R_CTF_M_q_ori__)
      deallocate(S_T_R_x_T_R_CTF_M_q_ori__)
      deallocate(quad_T_k_p_)
      deallocate(quad_T_k_q_)
      deallocate(Z_S_svdr_)
      deallocate(Z_M_svdr_)
      deallocate(tmp_Z_M_svdr_)
      deallocate(quad_SX2_)
      deallocate(quad_MX2_)
      deallocate(CTF_R_S_use_)
      deallocate(Z_tmp_)
      deallocate(Z_S_svdd_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished ti8_check_SxZTRM_5]'
      end if !if (verbose.gt.0) then
      end
