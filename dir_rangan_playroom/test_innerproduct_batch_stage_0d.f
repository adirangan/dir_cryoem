c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine test_innerproduct_batch_stage_0d(verbose
     $     ,svd_calculation_type,n_r,grid_p_,n_w_,max_x_c,n_S,ld_S,S_p_
     $     ,distance_req,S_alpha_polar_a_,S_alpha_azimu_b_,L_,nl_max
     $     ,nm_sum,ll_sum ,T_nl_,T_vm_,T_dd_,T_ll_,T_lf_ ,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_ ,T_ls_,T_LT_,n_T_root_base,T_root_base_ ,n_M
     $     ,ld_M,M_p_,n_CTF ,ld_CTF,CTF_p_,polar_a_est_,azimu_b_est_
     $     ,gamma_z_est_ ,delta_x_est_,delta_y_est_,l2_norm_est_
     $     ,ctf_ind_est_ ,eps_target,N_pixels_in,displacement_max
     $     ,n_delta_x,n_delta_y ,n_gamma_z,delta_x_max_,delta_y_max_
     $     ,gamma_z_max_,C_M_ ,C_S_max_,C_Z_max_,S_used_all)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer svd_calculation_type
      integer *4 n_r,n_w_(0:n_r-1),n_S,ld_S,n_M,ld_M,n_CTF,ld_CTF
      real *8 grid_p_(0:n_r-1),max_x_c
      complex *16 S_p_(0:ld_S*n_S-1)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      real *8 L_(0:0)
      real *8 diameter_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
c$$$      T_ structure
      integer *4 T_nl_(0:0) !level
      real *8    T_vm_(0:0) !vertex center
      real *8    T_dd_(0:0) !diameter
      integer *4 T_ll_(0:0) !number of points from L_ in T_
      logical    T_lf_(0:0) !is leaf
      integer *4 T_c0_(0:0) !child_0 tesselation_index
      integer *4 T_c1_(0:0) !child_1 tesselation_index
      integer *4 T_c2_(0:0) !child_2 tesselation_index
      integer *4 T_c3_(0:0) !child_3 tesselation_index
      integer *4 T_ls_(0:0) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4 T_LT_(0:0) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4 T_root_base_(0:0) !T_ roots
      real *8 distance_req
      integer *4 S_used_all
      integer *4, allocatable :: S_used_sub_(:)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      complex *16 M_p_(0:ld_M*n_M-1)
      complex *16 CTF_p_(0:0)
      real *8 polar_a_est_(0:n_M-1),azimu_b_est_(0:n_M-1)
     $     ,gamma_z_est_(0:n_M-1)
      real *8 delta_x_est_(0:n_M-1),delta_y_est_(0:n_M-1)
     $     ,l2_norm_est_(0:n_M-1)
      real *8 ctf_ind_est_(0:n_M-1)
      real *8 eps_target,N_pixels_in
      real *8 displacement_max
      integer n_delta_x,n_delta_y,n_gamma_z
      real *8 delta_x_max_(0:n_S*n_M-1)
      real *8 delta_y_max_(0:n_S*n_M-1)
      real *8 gamma_z_max_(0:n_S*n_M-1)
      complex *16 C_M_(0:n_M-1)
      complex *16 C_S_max_(0:n_S*n_M-1)
      complex *16 C_Z_max_(0:n_S*n_M-1)
c$$$      declaration of svd-expansion and associated variables
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      logical warning_flag
      data warning_flag / .true. /
      real *8 R_max,K_max,delta,delta_max,N_pixels
      complex *16, allocatable :: Z_svds_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: C_S_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer l
c$$$      fftw for omp
      integer *8, allocatable :: fftw_plan_frwd__(:)
      integer *8, allocatable :: fftw_plan_back__(:)
      complex *16, allocatable :: fftw_in1__(:)
      complex *16, allocatable :: fftw_out__(:)
c$$$      fftw for local use
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      pointer (p_fftw_plan_frwd_last,fftw_plan_frwd_last)
      pointer (p_fftw_plan_back_last,fftw_plan_back_last)
      pointer (p_fftw_in1_last_,fftw_in1_last_)
      pointer (p_fftw_out_last_,fftw_out_last_)
      integer *8 fftw_plan_frwd_last,fftw_plan_back_last
      complex *16 fftw_in1_last_(*)
      complex *16 fftw_out_last_(*)
c$$$      indices
      integer nr,n_w_max,n_A,na,ns,nctf
c$$$      array of displacements and rotations to measure
      integer ndx,ndy,ngz,ndx_max,ndy_max,ngz_max
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)
      character(len=64) format_string
c$$$      array of innerproducts to hold measurements
      real *8 pi
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      integer *4 n_M_sub,nM_sub,r_tmp,a_tmp,n_tmp,i_tmp
      integer *4, allocatable :: n_M_sub_(:)
      integer *4, allocatable :: I_M_sub_(:)

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_stage_0d]'
      end if      
      if (verbose.gt.0) then
         write(6,'(A,I0,A,F5.3,A,F3.1,A,I0,A,I0,A,I0)') 'n_r: ' ,n_r
     $        ,'; eps_target: ',eps_target,'; N_pixels_in: ',N_pixels_in
     $        ,'; n_delta_x: ',n_delta_x ,'; n_delta_y: ',n_delta_y
     $        ,'; n_gamma_z: ',n_gamma_z
      end if

      do na=0,n_S*n_M-1
         delta_x_max_(na) = 0.0d0
         delta_y_max_(na) = 0.0d0
         gamma_z_max_(na) = 0.0d0
         C_S_max_(na) = (1.0d0,0.0d0)
         C_Z_max_(na) = (0.0d0,0.0d0)
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating indices'
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
         write(6,'(A)') ' Setting up array of displacements to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_delta_x and n_delta_y '
         write(6,'(A)')
     $        ' (i.e., displacement array dimension) '
         write(6,'(A)') ' to equal 1+4*N_pixels_in or so.'
      end if
      allocate(delta_x_(0:n_delta_x-1))
      allocate(delta_y_(0:n_delta_y-1))
      do ndx=0,n_delta_x-1
         if (n_delta_x.gt.1) then
            delta_x = (-N_pixels_in + ndx*2*N_pixels_in/(n_delta_x-1))
     $           /n_r*max_x_c
         else
            delta_x = 0.0d0
         end if
         delta_x_(ndx) = delta_x
         do ndy=0,n_delta_y-1
            if (n_delta_y.gt.1) then
               delta_y = (-N_pixels_in + ndy*2*N_pixels_in/(n_delta_y
     $              -1))/n_r*max_x_c
            else
               delta_y = 0.0d0
            end if
            delta_y_(ndy) = delta_y
         enddo
      enddo
      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up array of rotations to measure'
         write(6,'(A)')
     $        ' It is usually reasonable for n_gamma_z (i.e,. the '
         write(6,'(A)') ' dimensions of the rotation array) to equal '
         write(6,'(A)')
     $        ' ncur or so.'
      end if
      allocate(gamma_z_(0:n_gamma_z-1))
      do ngz=0,n_gamma_z-1
         if (n_gamma_z.gt.1) then
            gamma_z = (2*pi*ngz)/n_gamma_z
         else
            gamma_z = 0.0d0
         end if
         gamma_z_(ngz) = gamma_z
      enddo
      if (verbose.gt.1) then
         write(format_string,'(A,I0,A)') '(A,',n_delta_x,'(F5.3,1X))'
         write (6,format_string) 'delta_x_: ',(delta_x_(ndx),ndx=0
     $        ,n_delta_x-1)
         write(format_string,'(A,I0,A)') '(A,',n_delta_y,'(F5.3,1X))'
         write (6,format_string) 'delta_y_: ',(delta_y_(ndy),ndy=0
     $        ,n_delta_y-1)
         write(format_string,'(A,I0,A)') '(A,',n_gamma_z,'(F5.3,1X))'
         write (6,format_string) 'gamma_z_: ',(gamma_z_(ngz),ngz=0
     $        ,n_gamma_z-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Selecting svd library to use.'
      end if
      R_max = 2*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
            delta = dsqrt(delta_x_(ndx)**2 + delta_y_(ndy)**2)
            if (delta.gt.delta_max) then
               delta_max = delta
            end if
         enddo
      enddo
      N_pixels = delta_max/dsqrt(2.0d0)*2*K_max
      if (verbose.gt.0) then
         write(6,'(A,F8.3,A,F8.3)') 'R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') 'K_max: ',K_max,'; N_pixels: '
     $        ,N_pixels
      end if
      if (K_max.gt.96 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, K_max ',K_max
     $        ,' too large in test_innerproduct_batch_stage_0d'
      end if
      if (N_pixels.gt.3 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, N_pixels ',N_pixels
     $        ,' too large in test_innerproduct_batch_stage_0d'
      end if
      if (eps_target.lt.0.1d0 .and. warning_flag) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in test_innerproduct_batch_stage_0d'
      end if
      if (.false.) then
      include './dir_gen_Jsvd_1/gen_Jsvd_svdpick.txt'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdload.txt'
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
         write(6,'(A,I0)') 'svd_unitnumber: ',svd_unitnumber
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Setting up displacement-operator Z_svdd_'
         write(6,'(A)') ' associated with svd-expansion.'
      end if
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      timing_tic = second()
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svdd: total time ',timing_toc
     $        -timing_tic
      end if

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for local use'
      end if
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_ESTIMATE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_ESTIMATE) 
         na = na + n_w_(nr)
      enddo
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_in1_last_ = loc(fftw_in1_(n_A-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_A-n_w_max))

      if (verbose.gt.1) then
         write(6,'(A,A)') ' Allocating array C_S_ to hold'
     $        ,' innerproducts for each CTF-S pair.'
      end if
      allocate(C_S_(0:n_gamma_z*n_S*n_CTF-1))
      if (verbose.gt.1) then
         write(6,'(A)') ' Calculate C_S_ for each ctf-S pair.'
         write(6,'(A)') ' C_S = || CTF .* R(S) ||.'
         write(6,'(A)') ' More specifically, the value: '
         write(6,'(A)')
     $        ' C_S_(ngz + ns*n_gamma_z + nctf*n_gamma_z*n_S) '
         write(6,'(A)') ' is equal to the l2_norm (not squared)'
         write(6,'(A)') ' of the pointwise product of CTF '
         write(6,'(A)') ' and R(S), where CTF is the nctf-th '
         write(6,'(A)') ' ctf-function, S is the ns-th template, '
         write(6,'(A)') ' and R is rotation by +gamma_z_(ngz).'
         write(6,'(A)') ' (see test_innerproduct_timing_dr.f) '
      end if
      timing_tic = second()
      do ns=0,n_S-1
         do nctf=0,n_CTF-1
            call test_innerproduct_batch_excerpt_00a(n_gamma_z,gamma_z_
     $           ,fftw_plan_frwd_,fftw_in1_,fftw_out_
     $           ,fftw_plan_back_last ,fftw_in1_last_,fftw_out_last_,n_r
     $           ,grid_p_,n_w_,n_A,S_p_(ns*ld_S),CTF_p_(nctf*ld_CTF)
     $           ,C_S_(0 + ns*n_gamma_z + nctf*n_gamma_z*n_S))
            if (verbose.gt.2) then
               write(format_string,'(A,I0,A)') '(A,',n_gamma_z
     $              ,'(F5.3,F5.3,1X))'
               write(6,format_string) 'C_S_: ',(C_S_(ngz),ngz= 0 + ns
     $              *n_gamma_z + nctf*n_gamma_z*n_S , n_gamma_z - 1 + ns
     $              *n_gamma_z + nctf*n_gamma_z*n_S)
            end if ! if (verbose.gt.2) then
         enddo ! do nctf=0,n_CTF-1
      enddo ! do ns=0,n_S-1
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculating C_S_ for each ctf-S pair.'
     $        ,timing_toc-timing_tic
      end if

      allocate(n_M_sub_(0:n_M-1))
      allocate(I_M_sub_(0:n_M-1))
      n_M_sub = min(8,n_M)
c$$$      n_M_sub = 1
c$$$      write(6,'(A,A,I0)') ' Warning! Turning off omp for innerproducts '
c$$$     $     , ' Testing: n_M_sub: ',n_M_sub
      if (verbose.gt.0) then
         write (6,'(A,I0,A)') 'Setting up n_M_sub = ',n_M_sub
     $        ,' sub-blocks for omp parallelization'
      end if
      do nM_sub=0,n_M_sub-1
         if (nM_sub.eq.0) then
            n_M_sub_(0) = n_M/n_M_sub
            I_M_sub_(0) = 0
         else if (nM_sub.gt.0 .and. nM_sub.lt.n_M_sub-1) then
            n_M_sub_(nM_sub) = n_M /n_M_sub
            I_M_sub_(nM_sub) = I_M_sub_(nM_sub-1) + n_M_sub_(nM_sub-1)
         else if (nM_sub.eq.n_M_sub-1) then
            I_M_sub_(n_M_sub-1) = I_M_sub_(n_M_sub-2) + n_M_sub_(n_M_sub
     $           -2)
            n_M_sub_(n_M_sub-1) = n_M - I_M_sub_(n_M_sub-1)
         end if
      enddo
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_M_sub
     $        ,'(I0,1X))'
         write(6,format_string) 'n_M_sub_: '
     $        ,(n_M_sub_(nM_sub),nM_sub=0
     $        ,n_M_sub-1)
         write(6,format_string) 'I_M_sub_: '
     $        ,(I_M_sub_(nM_sub),nM_sub=0
     $        ,n_M_sub-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
      end if
      allocate(fftw_plan_frwd__(0:n_r*n_M_sub-1))
      allocate(fftw_plan_back__(0:n_r*n_M_sub-1))
      allocate(fftw_in1__(0:n_A*n_M_sub-1))
      allocate(fftw_out__(0:n_A*n_M_sub-1))
      do nM_sub=0,n_M_sub-1
         r_tmp = nM_sub*n_r
         a_tmp = nM_sub*n_A
         na = 0
         do nr=0,n_r-1
            call dfftw_plan_dft_1d_(fftw_plan_frwd__(r_tmp+nr),n_w_(nr)
     $          ,fftw_in1__(a_tmp+na),fftw_out__(a_tmp +na)
     $          ,FFTW_FORWARD,FFTW_ESTIMATE) 
            call dfftw_plan_dft_1d_(fftw_plan_back__(r_tmp+nr),n_w_(nr)
     $          ,fftw_out__(a_tmp+na),fftw_in1__(a_tmp +na)
     $          ,FFTW_BACKWARD,FFTW_ESTIMATE) 
            na = na + n_w_(nr)
         enddo
      enddo

      allocate(S_used_sub_(0:n_M-1))
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_tmp,i_tmp,r_tmp,a_tmp)
c$OMP DO 
      do nM_sub=0,n_M_sub-1
         n_tmp = n_M_sub_(nM_sub)
         i_tmp = I_M_sub_(nM_sub)
         r_tmp = nM_sub*n_r
         a_tmp = nM_sub*n_A
         if (verbose.gt.0) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' test_innerproduct_batch_stage_1'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call test_innerproduct_batch_stage_1d(verbose
     $        ,svd_calculation_type,n_r,grid_p_,n_w_,max_x_c,n_S,ld_S
     $        ,S_p_,distance_req,S_alpha_polar_a_,S_alpha_azimu_b_,L_
     $        ,nl_max,nm_sum ,ll_sum,T_nl_,T_vm_,T_dd_,T_ll_ ,T_lf_
     $        ,T_c0_,T_c1_,T_c2_ ,T_c3_,T_ls_,T_LT_,n_T_root_base
     $        ,T_root_base_,n_tmp,ld_M ,M_p_(ld_M*i_tmp),n_CTF,ld_CTF
     $        ,CTF_p_ ,polar_a_est_(i_tmp) ,azimu_b_est_(i_tmp)
     $        ,gamma_z_est_(i_tmp) ,delta_x_est_(i_tmp)
     $        ,delta_y_est_(i_tmp) ,l2_norm_est_(i_tmp)
     $        ,ctf_ind_est_(i_tmp),eps_target ,N_pixels_in
     $        ,displacement_max,n_delta_x,n_delta_y ,n_gamma_z
     $        ,delta_x_max_(n_S*i_tmp),delta_y_max_(n_S *i_tmp)
     $        ,gamma_z_max_(n_S*i_tmp),C_M_(i_tmp),C_S_max_(n_S *i_tmp)
     $        ,C_Z_max_(n_S *i_tmp),fftw_plan_frwd__(r_tmp)
     $        ,fftw_plan_back__(r_tmp) ,fftw_in1__(a_tmp)
     $        ,fftw_out__(a_tmp),n_svd_r,n_svd_d ,n_svd_l,svd_r_,svd_d_
     $        ,svd_l_,svd_U_d_,svd_s_,svd_V_r_ ,Z_svdd_,C_S_
     $        ,S_used_sub_(nM_sub))
      enddo !do nM_sub=0,n_M_sub-1
c$OMP END DO
c$OMP END PARALLEL
      S_used_all = sum_i4_f(n_M_sub,S_used_sub_)
      timing_toc = omp_get_wtime()
      if (verbose.gt.1) then
         write(6,'(A,A,A,F8.5)') ' Finished',
     $        ' test_innerproduct_batch_stage_1:', 
     $        ' total time ',
     $        timing_toc-timing_tic
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating temporary arrays.'
      end if
      deallocate(C_S_)
      deallocate(Z_svdd_)
      deallocate(gamma_z_)
      deallocate(delta_y_)
      deallocate(delta_x_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan_(fftw_plan_back_(nr))
         call dfftw_destroy_plan_(fftw_plan_frwd_(nr))
      enddo
      deallocate(fftw_out_)
      deallocate(fftw_in1_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_plan_frwd_)
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for omp.'
      end if
      do nr=0,n_M_sub*n_r-1
         call dfftw_destroy_plan_(fftw_plan_back__(nr))
         call dfftw_destroy_plan_(fftw_plan_frwd__(nr))
      enddo
      deallocate(fftw_out__)
      deallocate(fftw_in1__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_plan_frwd__)
      deallocate(n_M_sub_)
      deallocate(I_M_sub_)

      if (verbose.gt.1) then
         write(6,'(A)') ' Deallocating svd arrays.'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdfree.txt'

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_stage_0d]'
      end if

      end


