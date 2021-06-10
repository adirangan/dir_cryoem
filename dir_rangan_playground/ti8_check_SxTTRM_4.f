!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Tests S_T_T_R_CTF_M_q_ calculation directly. ;\n
      subroutine ti8_check_SxTTRM_4(
     $     verbose
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
     $     ,S_T_T_R_CTF_M_q__
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      test S_T_T_R_CTF_M_q_.
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
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
      complex *16, allocatable :: CTF_R_S_use_(:)
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
      real *8 pi
      integer max_i4_f !function output. ;
      real *8 al2_c16_f, fnod_c16_f !function output. ;
      integer *4 MDA_n_d,MDA_d_(0:127)
      character(len=1024) MDA_string
      logical flag_memory_checkset
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_check_SxTTRM_4]'
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
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
               call rotate_p_to_p_fftw(quad_n_r,fftw_plan_frwd_
     $              ,fftw_plan_back_,quad_n_w_,quad_n_w_sum,fftw_0in_
     $              ,fftw_out_,quad_S_k_p_,gamma_z,quad_T_k_p_)
               call transf_p_to_p(quad_n_r,quad_grid_k_p_r_,quad_n_w_
     $              ,quad_n_w_sum,quad_T_k_p_,delta_x,delta_y
     $              ,quad_T_k_p_)
               call interp_p_to_q_fftw(quad_n_r,fftw_plan_frwd_
     $              ,quad_n_w_,quad_n_w_sum,fftw_0in_,fftw_out_
     $              ,quad_T_k_p_ ,quad_T_k_q_)
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
               if (verbose.gt.2) then
                  write(6,'(A,I0,2(A,F24.16),A,I0,A,F24.16)')
     $                 ' ndv: ' , ndv 
     $                 ,' delta_x: ' , delta_x
     $                 ,' delta_y: ' , delta_y
     $                 , ' nw: ' , nw
     $                 , ' gamma_z: ' , gamma_z
                  write(6,'(2(A,2F24.16),A,2F24.16)')
     $                 ' tmp_S_T_T_R_CTF_M_q__ ' 
     $                 ,tmp_S_T_T_R_CTF_M_q__(tab)
     $                 ,' S_T_R_x_T_R_CTF_M_q_ori__ ' 
     $                 ,S_T_R_x_T_R_CTF_M_q_ori__(tab)
     $                 , ' error '
     $                 ,zabs(tmp_S_T_T_R_CTF_M_q__(tab) -
     $                 S_T_R_x_T_R_CTF_M_q_ori__(tab))
     $                 /zabs(tmp_S_T_T_R_CTF_M_q__(tab))
                  write(6,'(2(A,2F24.16),A,2F24.16)')
     $                 ' tmp_S_T_T_R_CTF_M_q__ ' 
     $                 ,tmp_S_T_T_R_CTF_M_q__(tab)
     $                 ,' S_R_T_x_T_R_CTF_M_q_ori__ ' 
     $                 ,S_R_T_x_T_R_CTF_M_q_ori__(tab)
     $                 , ' error '
     $                 ,zabs(tmp_S_T_T_R_CTF_M_q__(tab) -
     $                 S_R_T_x_T_R_CTF_M_q_ori__(tab))
     $                 /zabs(tmp_S_T_T_R_CTF_M_q__(tab))
               end if !if (verbose.gt.2) then
c$$$           %%%%%%%%               
            enddo !do nw=0,quad_n_w_max-1
         enddo !do ndv=0,n_delta_v-1
c$$$         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (verbose.gt.1) then
         write(6,'(A,A,F24.16)')
     $        ' tmp_S_T_T_R_CTF_M_q__ vs ' 
     $        ,' S_T_R_x_T_R_CTF_M_q_ori__ ' 
     $        ,fnod_c16_f(quad_n_w_max*n_delta_v
     $        ,tmp_S_T_T_R_CTF_M_q__
     $        ,S_T_R_x_T_R_CTF_M_q_ori__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v
     $        ,tmp_S_T_T_R_CTF_M_q__)
         write(6,'(A,A,F24.16)')
     $        ' tmp_S_T_T_R_CTF_M_q__ vs ' 
     $        ,' S_R_T_x_T_R_CTF_M_q_ori__ ' 
     $        ,fnod_c16_f(quad_n_w_max*n_delta_v
     $        ,tmp_S_T_T_R_CTF_M_q__
     $        ,S_R_T_x_T_R_CTF_M_q_ori__)
     $        /al2_c16_f(quad_n_w_max*n_delta_v
     $        ,tmp_S_T_T_R_CTF_M_q__)
         end if !if (verbose.gt.1) then

      enddo ! do nx1=0,n_S*n_M-1

      flag_memory_checkset = .true.
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

      deallocate(S_R_T_x_T_R_CTF_M_q_ori__)
      deallocate(S_T_R_x_T_R_CTF_M_q_ori__)
      deallocate(quad_T_k_p_)
      deallocate(quad_T_k_q_)
      deallocate(CTF_R_S_use_)
      deallocate(Z_tmp_)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished ti8_check_SxTTRM_4]'
      end if !if (verbose.gt.0) then
      end

