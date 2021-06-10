!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Performs search-and-store routine. ;\n
      subroutine ti8_build_Zstore_full_0(
     $     verbose
     $     ,stage_type
     $     ,rseed
     $     ,flag_RTRT_vs_RTTR
     $     ,n_delta_v
     $     ,delta_x_
     $     ,delta_y_
     $     ,n_gamma_z
     $     ,gamma_z_
     $     ,fftw_plan_frwd_
     $     ,fftw_plan_back_ 
     $     ,fftw_0in_
     $     ,fftw_out_
     $     ,nufft1df77_fw_
     $     ,nufft1df77_lw
     $     ,nufft1df77_lused
     $     ,n_r 
     $     ,grid_p_
     $     ,n_w_
     $     ,n_w_sum
     $     ,n_S_sub
     $     ,n_S_sum
     $     ,n_S_tmp
     $     ,S_alpha_S_index_
     $     ,S_alpha_polar_a_
     $     ,S_alpha_azimu_b_
     $     ,n_M_sub
     $     ,polar_a_est_
     $     ,azimu_b_est_
     $     ,gamma_z_est_ 
     $     ,delta_x_est_
     $     ,delta_y_est_
     $     ,l2_norm_est_
     $     ,ctf_ind_est_ 
     $     ,S_index_est_
     $     ,M_index_est_
     $     ,alpha_0in_
     $     ,alpha_update_f
     $     ,displacement_max
     $     ,flag_MS_vs_SM
     $     ,n_SM_max
     $     ,n_SM_
     $     ,alpha_SM__
     $     ,n_MS_max
     $     ,n_MS_
     $     ,alpha_MS__
     $     ,C_M_
     $     ,n_CTF
     $     ,CTF_R_S_ 
     $     ,CTF_R_S_use_ 
     $     ,n_M_tot
     $     ,S_T_T_R_CTF_M_q__ 
     $     ,ZZ_
     $     ,ZZ_sub_
     $     ,ZZ_pos_
     $     )
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      character stage_type ! character: indicates what part of the calculation to carry out. if stage_type.eq.'0' then the entire calculation will be carried out. otherwise, only part of the subroutine will be executed. used for timing. ;
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      integer *4 rseed !random seed. ;
      real *8 adi_rand_f
      integer n_delta_v,n_gamma_z,n_CTF
      integer n_r,n_w_(0:n_r-1),n_w_sum
      real *8 S_alpha_S_index_(0:0)
      real *8 S_alpha_polar_a_(0:0)
      real *8 S_alpha_azimu_b_(0:0)
      real *8 delta_x_(0:n_delta_v-1),delta_y_(0:n_delta_v-1)
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      integer *8 fftw_plan_back_(0:n_r-1)
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_0in_last_,fftw_0in_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last
      integer *8 fftw_plan_back_last
      complex *16 fftw_0in_last_(0:0),fftw_out_last_(0:0)
      real *8 nufft1df77_fw_(0:0)
      integer *4 nufft1df77_lw,nufft1df77_lused
      logical flag_nufft
      real *8 polar_a_est_(0:n_M_sub-1),azimu_b_est_(0:n_M_sub-1)
     $     ,gamma_z_est_(0:n_M_sub-1)
      real *8 delta_x_est_(0:n_M_sub-1),delta_y_est_(0:n_M_sub-1)
     $     ,l2_norm_est_(0:n_M_sub-1)
      real *8 ctf_ind_est_(0:n_M_sub-1)
      real *8 S_index_est_(0:n_M_sub-1)
      real *8 M_index_est_(0:n_M_sub-1)
      real *8 displacement_max
      logical flag_MS_vs_SM !logical: determines whether to assign images to templates (.true.) or templates to images (.false.). ;
c$$$      SM storage
      integer *4 n_SM_max ! total (maximum) number of templates to store per image. ;
      integer *4 n_SM_(0:0) ! array of size n_M indicating the actual number of templates stored per image. ;
      real *8 alpha_SM__(0:0) ! array of size n_alpha*n_SM_max*n_M storing the image-parameters for each stored template-image pair. ;
c$$$      MS storage
      integer *4 n_MS_max ! total (maximum) number of images to store per template. ;
      integer *4 n_MS_(0:0) ! array of size n_S indicating the actual number of images stored per template. ;
      real *8 alpha_MS__(0:0) ! array of size n_alpha*n_MS_max*n_S storing the image-parameters for each stored iamge-template pair. ;
      real *8 alpha_update_f
      include 'excerpt_define_nalpha.f'
      real *8 alpha_0in_(0:n_alpha-1)
      integer n_S_sub,n_S_sum,n_S_tmp,n_M_sub,n_M_tot
      complex *16 C_M_(0:0),C_M
      complex *16 CTF_R_S_(0:0)
      complex *16 CTF_R_S_use_(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      complex *16 ZZ_(0:0)
      complex *16 ZZ_sub_(0:0)
      complex *16 ZZ_pos_(0:0)
      complex *16 ZZ,ZZ_pos
      complex *16 CTF_R_S_use,C_Z_use
      integer n_w_max,nr,nw,nC,nC_stride
      integer ndv,ngz
      integer ns,nm,nctf,ns_use,nr_use
      integer nx2,nx3,nx3_stride
      real *8 delta_x,delta_y,gamma_z,gamma_z_est,tmp_gamma_z
      integer ndv_optimal,ngz_optimal
      real *8 delta_x_optimal,delta_y_optimal,gamma_z_optimal
      real *8 delta_x_tmp,delta_y_tmp,gamma_z_tmp
      complex *16 C_Z_optimal,CTF_R_S_optimal,l2_norm_optimal
      real *8 pi
      character(len=1024) format_string

      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_build_Zstore_full_0]'
      end if !if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ' , verbose
         write(6,'(A,I0)') ' rseed: ' , rseed
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' , flag_RTRT_vs_RTTR
         write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
         call print_sub_r8(n_delta_v,delta_x_,' delta_x_: ')
         write(6,'(A,I0)') ' n_delta_v: ' , n_delta_v
         call print_sub_r8(n_delta_v,delta_y_,' delta_y_: ')
         write(6,'(A,I0)') ' n_gamma_z: ' , n_gamma_z
         call print_sub_r8(n_gamma_z,gamma_z_,' gamma_z_: ')
         write(6,'(A,I0)') ' n_r : ' , n_r 
         call print_sub_r8(n_r,grid_p_,' grid_p_: ')
         call print_sub_i4(n_r,n_w_,' n_w_: ')
         write(6,'(A,I0)') ' n_w_sum: ' , n_w_sum
         write(6,'(A,I0)') ' n_S_sub: ' , n_S_sub
         write(6,'(A,I0)') ' n_S_sum: ' , n_S_sum
         write(6,'(A,I0)') ' n_S_tmp: ' , n_S_tmp
         call print_sub_r8(n_S_sub,S_alpha_S_index_
     $        ,' S_alpha_S_index_: ')
         call print_sub_r8(n_S_sub,S_alpha_polar_a_
     $        ,' S_alpha_polar_a_: ')
         call print_sub_r8(n_S_sub,S_alpha_azimu_b_
     $        ,' S_alpha_azimu_b_: ')
         write(6,'(A,I0)') ' n_M_sub: ' , n_M_sub
         call print_sub_r8(n_M_sub,polar_a_est_,' polar_a_est_: ')
         call print_sub_r8(n_M_sub,azimu_b_est_,' azimu_b_est_: ')
         call print_sub_r8(n_M_sub,gamma_z_est_,' gamma_z_est_: ')
         call print_sub_r8(n_M_sub,delta_x_est_,' delta_x_est_: ')
         call print_sub_r8(n_M_sub,delta_y_est_,' delta_y_est_: ')
         call print_sub_r8(n_M_sub,l2_norm_est_,' l2_norm_est_: ')
         call print_sub_r8(n_M_sub,ctf_ind_est_,' ctf_ind_est_: ')
         call print_sub_r8(n_M_sub,S_index_est_,' S_index_est_: ')
         call print_sub_r8(n_M_sub,M_index_est_,' M_index_est_: ')
         write(6,'(A,F8.4)') ' alpha_update_f: ' , alpha_update_f
         write(6,'(A,F8.4)') ' displacement_max: ' , displacement_max
         write(6,'(A,L1)') ' flag_MS_vs_SM: ' , flag_MS_vs_SM
         if (flag_MS_vs_SM.eqv..false.) then
         write(6,'(A,I0)') ' n_SM_max: ' , n_SM_max
         call print_sub_i4(n_M_sub,n_SM_,' n_SM_: ')
         call print_sub_r8(n_alpha*n_SM_max*n_M_sub,alpha_SM__
     $        ,' alpha_SM__: ')
         end if !if (flag_MS_vs_SM.eqv..false.) then
         if (flag_MS_vs_SM.eqv..true.) then
         write(6,'(A,I0)') ' n_MS_max: ' , n_MS_max
         call print_sub_i4(n_S_sub,n_MS_,' n_MS_: ')
         call print_sub_r8(n_alpha*n_MS_max*n_S_sub,alpha_MS__
     $        ,' alpha_MS__: ')
         end if !if (flag_MS_vs_SM.eqv..true.) then
         call print_sub_c16(n_M_sub,C_M_,' C_M_: ')
         write(6,'(A,I0)') ' n_CTF: ' , n_CTF
         call print_sub_c16(n_gamma_z*n_CTF*n_S_sub,
     $        CTF_R_S_,' CTF_R_S_: ')
         write(6,'(A,I0)') ' n_M_tot: ' , n_M_tot
         call print_sub_c16(n_delta_v*n_S_sub*n_M_sub
     $        ,S_T_T_R_CTF_M_q__,' S_T_T_R_CTF_M_q__: ')
         call print_sub_c16(n_delta_v*n_gamma_z
     $        ,ZZ_,' ZZ_: ')
      end if !if (verbose.gt.1) then
      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)
      flag_nufft = .false.
      if (n_gamma_z.ne.n_w_max) then
         flag_nufft = .true.
      end if !if (n_gamma_z.ne.n_w_max) then
      if (n_gamma_z.eq.n_w_max) then
         do nw=0,n_w_max-1
            tmp_gamma_z = (2.0d0*pi*nw)/(1.0d0*n_gamma_z)
            ngz = nw
            gamma_z = gamma_z_(ngz)
            if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
               flag_nufft = .true.
            end if !if (dabs(tmp_gamma_z-gamma_z).gt.1.0d-12) then
         enddo !do nw=0,n_w_max-1
      end if !if (n_gamma_z.eq.n_w_max) then
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_w_max: ' , n_w_max
         write(6,'(A,L2)') ' flag_nufft: ' , flag_nufft
         write(6,'(A)') ' linking fftw_plan_frwd_last, etc. '
      end if !if (verbose.gt.1) then      
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_0in_last_ = loc(fftw_0in_(n_w_sum-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_w_sum-n_w_max))

      call cl1_r8(n_alpha,alpha_0in_)
      call cl1_c16(n_gamma_z,CTF_R_S_use_)
      call cl1_c16(n_w_max,ZZ_sub_)
      call cl1_c16(n_gamma_z,ZZ_pos_)

c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do nm=0,n_M_sub-1
         C_M = C_M_(nm)
         nctf = nint(ctf_ind_est_(nm))
         gamma_z_est = gamma_z_est_(nm)
         do ns=0,n_S_sum-1
            if ((stage_type.eq.'0') .or. (stage_type.eq.'a')) then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage a'
            end if !verbose
            if (flag_nufft.eqv..false.) then
            call rotate_p_to_p_fftw_single(
     $           fftw_plan_frwd_last
     $           ,fftw_plan_back_last
     $           ,n_w_max
     $           ,fftw_0in_last_
     $           ,fftw_out_last_
     $           ,CTF_R_S_(0+n_gamma_z*(nctf + n_CTF*ns))
     $           ,-gamma_z_est
     $           ,CTF_R_S_use_
     $           )
            if (verbose.gt.2) then
            write(6,*) ' ns: ' , ns , ' nctf: ' , nctf , ' n_gamma_z: '
     $           , n_gamma_z , ' tab: ' , n_gamma_z*(nctf + n_CTF*ns)
            call print_sub_c16(n_w_max,CTF_R_S_(0+n_gamma_z*(nctf +
     $           n_CTF*ns)),' CTF_R_S_: ')
            call print_sub_c16(n_w_max,CTF_R_S_use_,' CTF_R_S_use_: ')
            end if !if (verbose.gt.2) then
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
            end if !if ((stage_type.eq.'0') .or. (stage_type.eq.'a')) then
            if ((stage_type.eq.'0') .or. (stage_type.eq.'b')) then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage b'
            end if !verbose
            do ngz=0,n_gamma_z-1
               if (zabs(C_M*CTF_R_S_use_(ngz)).le.1.0d-15) then
                  CTF_R_S_use_(ngz) = 1.0d0
               else
                  CTF_R_S_use_(ngz) = 1.0d0/(C_M*CTF_R_S_use_(ngz))
               end if !if (zabs(C_M*CTF_R_S_use_(ngz)).le.1.0d-15) then               
            enddo !do ngz=0,n_gamma_z-1
            call cl1_c16(n_delta_v*n_gamma_z,ZZ_)
            if (flag_nufft.eqv..false.) then
            nx2 = n_delta_v*(ns + n_S_sub*nm)
            nx3_stride = n_delta_v*n_S_sub*n_M_tot
            nC_stride = n_delta_v
            do ndv=0,n_delta_v-1
               nx3 = 0
               nC = ndv
               do nw=0,n_w_max-1
                  ngz = nw
                  CTF_R_S_use = CTF_R_S_use_(ngz)
                  ZZ_(nC) = S_T_T_R_CTF_M_q__(nx2 + nx3)*CTF_R_S_use
                  nx3 = nx3 + nx3_stride
                  nC = nC + nC_stride
               enddo            !do nw
               nx2 = nx2 + 1
            enddo ! do ndv
            else
            if (verbose.gt.1) then
               write(6,'(A)') ' calling interp1_lagrange_c16'
            end if !if (verbose.gt.1) then
            do ndv=0,n_delta_v-1
               delta_y = delta_y_(ndv)
               delta_x = delta_x_(ndv)
               call cl1_c16(n_w_max,ZZ_sub_)
               do nw=0,n_w_max-1
                  nx2 = ndv + n_delta_v*(ns + n_S_sub*nm)
                  nx3 = nw*n_delta_v*n_S_sub*n_M_tot
                  ZZ_sub_(nw) = S_T_T_R_CTF_M_q__(nx2 + nx3)
               enddo            !do nw
c$$$               Unfortunately, I could not get either finufft1d2 or nufft1d2 to work across multiple threads. ;
c$$$               call interp1_nufft1df77_c16(
c$$$     $              n_w_max
c$$$     $              ,0.0d0
c$$$     $              ,2.0d0*pi
c$$$     $              ,ZZ_sub_
c$$$     $              ,fftw_plan_frwd_last
c$$$     $              ,fftw_plan_back_last
c$$$     $              ,fftw_0in_last_
c$$$     $              ,fftw_out_last_
c$$$     $              ,nufft1df77_fw_
c$$$     $              ,nufft1df77_lw
c$$$     $              ,nufft1df77_lused
c$$$     $              ,n_gamma_z
c$$$     $              ,gamma_z_
c$$$     $              ,ZZ_pos_
c$$$     $              )
               call interp1_lagrange_c16(
     $              n_w_max
     $              ,0.0d0
     $              ,2.0d0*pi
     $              ,ZZ_sub_
     $              ,n_gamma_z
     $              ,gamma_z_
     $              ,ZZ_pos_
     $              )
               do ngz=0,n_gamma_z-1
                  CTF_R_S_use = CTF_R_S_use_(ngz)
                  gamma_z = gamma_z_(ngz)
                  ZZ_pos = ZZ_pos_(ngz);
                  nC = ndv + n_delta_v*ngz
                  ZZ_(nC) = ZZ_pos*CTF_R_S_use
               enddo            !do ngz
            enddo ! do ndv
            end if !if (flag_nufft.eqv..false.) then
            end if !if ((stage_type.eq.'0') .or. (stage_type.eq.'b')) then
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ((stage_type.eq.'0') .or. (stage_type.eq.'c')) then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage c'
            end if !verbose
            call ti8_ZZ_scan_4(flag_RTRT_vs_RTTR
     $           ,delta_x_est_(nm),delta_y_est_(nm),gamma_z_est_(nm)
     $           ,displacement_max,n_delta_v ,delta_x_ ,delta_y_
     $           ,n_gamma_z,gamma_z_,ZZ_ ,ndv_optimal ,ngz_optimal
     $           ,C_Z_optimal)
            delta_x_optimal = delta_x_(ndv_optimal)
            delta_y_optimal = delta_y_(ndv_optimal)
            gamma_z_optimal = gamma_z_(ngz_optimal)
            end if !if ((stage_type.eq.'0') .or. (stage_type.eq.'c')) then
            if (stage_type.eq.'x') then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage x'
            end if !verbose
            call ti8_ZZ_scan_4x(flag_RTRT_vs_RTTR
     $           ,delta_x_est_(nm),delta_y_est_(nm),gamma_z_est_(nm)
     $           ,displacement_max,n_delta_v ,delta_x_ ,delta_y_
     $           ,n_gamma_z,gamma_z_,ZZ_ ,ndv_optimal ,ngz_optimal
     $           ,C_Z_optimal)
            delta_x_optimal = delta_x_(ndv_optimal)
            delta_y_optimal = delta_y_(ndv_optimal)
            gamma_z_optimal = gamma_z_(ngz_optimal)
            end if !if (stage_type.eq.'x') then
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ((stage_type.eq.'0') .or. (stage_type.eq.'d')) then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage d'
            end if !verbose
            if (flag_nufft.eqv..false.) then
            call rotate_p_to_p_fftw_single(
     $           fftw_plan_frwd_last
     $           ,fftw_plan_back_last
     $           ,n_w_max
     $           ,fftw_0in_last_
     $           ,fftw_out_last_
     $           ,CTF_R_S_(0+n_gamma_z*(nctf + n_CTF*ns))
     $           ,-gamma_z_est_(nm)-gamma_z_optimal
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
     $           ,-gamma_z_est_(nm)-gamma_z_optimal
     $           ,CTF_R_S_use_
     $           )
            end if !if (flag_nufft.eqv..true.) then
            CTF_R_S_optimal = CTF_R_S_use_(0)
            if (verbose.gt.2) then
               write(6,'(2(A,I0),A,I4,1X,I4,3(A,2F8.3))')
     $              ' Match between M_q_(' , nm , ') and S_q_(',ns
     $              ,') at ndv_optimal,ngz_optimal: '
     $              ,ndv_optimal ,ngz_optimal
     $              ,'; C_M: ',C_M_(nm) 
     $              ,'; C_Z_optimal: ' ,C_Z_optimal
     $              ,'; CTF_R_S_optimal: ',CTF_R_S_optimal
            end if !verbose            
            if (zabs(CTF_R_S_optimal).le.1.0d-15) then
               l2_norm_optimal = 1.0d0
            else 
               l2_norm_optimal = zabs(C_M_(nm)) *zabs(C_Z_optimal)
     $              /zabs(CTF_R_S_optimal)
            end if ! if (zabs(CTF_R_S_optimal).le.1.0d-15) then
            end if !if ((stage_type.eq.'0') .or. (stage_type.eq.'d')) then
            if ((stage_type.eq.'0') .or. (stage_type.eq.'e')) then
            if (verbose.gt.1) then 
               write(6,'(A)') 'stage e'
            end if !verbose
            call cl1_r8(n_alpha,alpha_0in_)
            alpha_0in_(nalpha_polar_a) = S_alpha_polar_a_(ns)
            alpha_0in_(nalpha_azimu_b) = S_alpha_azimu_b_(ns)
            alpha_0in_(nalpha_gamma_z) = gamma_z_optimal
            alpha_0in_(nalpha_delta_x) = delta_x_optimal
            alpha_0in_(nalpha_delta_y) = delta_y_optimal
            alpha_0in_(nalpha_l2_norm) = l2_norm_optimal
            alpha_0in_(nalpha_ctf_ind) = nctf
            alpha_0in_(nalpha_S_index) = S_alpha_S_index_(ns)
            alpha_0in_(nalpha_M_index) = M_index_est_(nm)
            alpha_0in_(nalpha_CTF_R_S) = CTF_R_S_optimal
            alpha_0in_(nalpha_C_Z_opt) = C_Z_optimal
            if (verbose.gt.2) then
               write(6,'(A)') ' calling alpha_SM_update_1 with '
               write(6,'(5(A,I0))') ' n_alpha ' , n_alpha , ' n_SM_max '
     $              , n_SM_max , ' nm ' , nm ,' nm*n_alpha*n_SM_max ' ,
     $              nm*n_alpha*n_SM_max , ' n_SM_(nm) ' , n_SM_(nm) 
               write(6,'(A,I0)') ' calling alpha_SM_update_1 with '
               call print_all_r8(n_alpha,alpha_0in_
     $              ,' alpha_0in_: ')
            end if !if (verbose.gt.2) then
            if (flag_MS_vs_SM.eqv..true.) then
c$$$        Note: we use alpha_SM_update to update MS as well.
            call alpha_SM_update_1(n_MS_max,n_MS_(ns)
     $           ,alpha_MS__(ns*n_alpha*n_MS_max),alpha_0in_)
            end if !if (flag_MS_vs_SM.eqv..true.) then
            if (flag_MS_vs_SM.eqv..false.) then
            call alpha_SM_update_1(n_SM_max,n_SM_(nm)
     $           ,alpha_SM__(nm*n_alpha*n_SM_max),alpha_0in_)
            end if !if (flag_MS_vs_SM.eqv..false.) then
            end if !if ((stage_type.eq.'0') .or. (stage_type.eq.'e')) then
         enddo !do ns=0,n_S_sum-1
      enddo !do nm=0,n_M_sub-1
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (verbose.gt.0) then
         write(6,'(A)') '[finished ti8_build_Zstore_full_0]'
      end if !if (verbose.gt.0) then

      end
