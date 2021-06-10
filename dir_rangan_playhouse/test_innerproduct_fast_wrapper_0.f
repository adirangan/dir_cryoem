      subroutine test_innerproduct_fast_wrapper_0 (verbose,n_omp_sub
     $     ,n_k_cur,n_polar_a_,grid_k_p_,half_diameter_x_c,n_S,ld_S
     $     ,S_k_p_,tesselation_distance_req ,S_alpha_polar_a_
     $     ,S_alpha_azimu_b_,n_M,ld_M,M_k_p_,n_CTF ,ld_CTF ,CTF_k_p_
     $     ,n_alpha,alpha_est_,alpha_update_f,n_pixels_in
     $     ,displacement_max ,n_delta_x ,n_delta_y,n_gamma_z
     $     ,svd_calculation_type ,eps_svd,fpm_howmany_max,C_M_
     $     ,n_SM_use)
c$$$      Reads in a few templates and images and compares them. 
c$$$      Uses adaptive tesselation to sample from templates.
c$$$      The output (i.e., which innerproducts are highest) is returned.
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer n_omp_sub ! number of batches to use when finding innerproducts
      integer *4 n_k_cur,n_polar_a_(0:n_k_cur-1),n_S,ld_S,n_M,ld_M,n_CTF
     $     ,ld_CTF,n_alpha
      integer *4 n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type
      real *8 grid_k_p_(0:0),half_diameter_x_c
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for tesselation
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      integer *4 n_point,npoint
      real *8, allocatable :: L_(:)
      real *8 tradius_min
      integer *4 nl_max,nm_sum,ll_sum
      integer *4 max_i4_f,sum_i4_f
      integer *4 n_T_root_base,nT_root_base
      logical parity_input
c$$$      T_ structure
      integer *4, allocatable :: T_nl_(:) !level
      real *8   , allocatable :: T_vm_(:) !vertex center
      real *8   , allocatable :: T_tr_(:) !tradius
      integer *4, allocatable :: T_ll_(:) !number of points from L_ in T_
      logical   , allocatable :: T_lf_(:) !is leaf
      integer *4, allocatable :: T_c0_(:) !child_0 tesselation_index
      integer *4, allocatable :: T_c1_(:) !child_1 tesselation_index
      integer *4, allocatable :: T_c2_(:) !child_2 tesselation_index
      integer *4, allocatable :: T_c3_(:) !child_3 tesselation_index
      integer *4, allocatable :: T_ls_(:) !starting index of point_index_list for T_ if leaf (leaves only)
      integer *4, allocatable :: T_LT_(:) !full point_index_list for all of T_ (leaves only)
c$$$      T_ roots
      integer *4, allocatable :: T_root_base_(:) !T_ roots
      real *8 tesselation_distance_req
      integer *4 n_SM_use
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      begin variables for innerproducts
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      real *8 alpha_est_(0:0)
      real *8 alpha_update_f
      include 'nalpha_define.f'
      real *8 eps_svd,n_pixels_in,displacement_max
      complex *16 S_k_p_(0:0),M_k_p_(0:0),CTF_k_p_(0:0)
      complex *16 C_M_(0:n_M-1)
      integer *4 ntemplatesize
      integer *4, allocatable :: ngridc_(:)
      integer *4, allocatable :: icstart_(:)
c$$$      indices
      integer *4 n_r,nr,n_w_max,n_A,na,ns,nm
      integer *4, allocatable :: n_w_(:)
      real *8 pi
      real *8 max_r8_f
      real *8 max_k_c
c$$$      estimated (current) displacement and rotation parameters for each image
      real *8, allocatable :: polar_a_est_(:)
      real *8, allocatable :: azimu_b_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
c$$$      estimated l2_norm and image-dependent ctf_index
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: ctf_ind_est_(:)
      real *8, allocatable :: S_index_est_(:)
c$$$      updated (future) displacement and rotation parameters for each image
      real *8, allocatable :: polar_a_upd_(:)
      real *8, allocatable :: azimu_b_upd_(:)
      real *8, allocatable :: gamma_z_upd_(:)
      real *8, allocatable :: delta_x_upd_(:)
      real *8, allocatable :: delta_y_upd_(:)
c$$$      updated l2_norm and image-dependent ctf_index
      real *8, allocatable :: l2_norm_upd_(:)
      real *8, allocatable :: ctf_ind_upd_(:)
      real *8, allocatable :: S_index_upd_(:)
c$$$      range of displacement and rotation parameters to measure
      integer ndx,ndy,ngz
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)      
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16 C_M,C_Z
      character(len=64) format_string
c$$$      fftw_plan_many
      integer fpm_howmany_max
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_fast_wrapper_0]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_k_cur: ',n_k_cur
         write(format_string,'(A,I0,A)') '(A,',n_k_cur,'(I0,1X))'
         write(6,format_string) 'n_polar_a_: ',(n_polar_a_(nr),nr=0
     $        ,n_k_cur-1)
         write(format_string,'(A,I0,A)') '(A,',n_k_cur,'(F8.3,1X))'
         write(6,format_string) 'grid_k_p_: ',(grid_k_p_(nr),nr=0
     $        ,n_k_cur-1)
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
         write(6,'(A,I0)') 'n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
         write(6,'(A)') 'alpha_est_: '
         write(6,format_string) (alpha_est_(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') 'n_pixels_in: ',n_pixels_in
         write(6,'(A,F6.3)') 'displacement_max: ',displacement_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
         write(6,'(A,F6.3)') 'eps_svd: ',eps_svd
      end if !if (verbose.gt.1) then

      pi = 4*atan(1.0)

c$$$      Calculating template size using 'get_template_size'
      allocate(ngridc_(n_k_cur))
      allocate(icstart_(n_k_cur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'n_k_cur = ',n_k_cur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(n_polar_a_,n_k_cur,ntemplatesize,ngridc_
     $     ,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
         write(format_string,'(A,I0,A)') '(A,',n_k_cur,'(I0,1X))'
         write(6,format_string) 'ngridc_: ',(ngridc_(nr),nr=1,n_k_cur)
      end if
      
c$$$      indices
      n_r = n_k_cur
      if (n_r.lt.2) then
         write(6,'(A,I0,A)') 'Error n_r',n_r,'<2'
      end if
      allocate(n_w_(0:n_r-1))
      n_A = 0
      do nr=0,n_r-1
         n_w_(nr) = ngridc_(1+nr)
         n_A = n_A + n_w_(nr)
      enddo
      n_w_max = n_w_(nr-1)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0)') 'n_w_max ',n_w_max,'; n_A ',n_A
         write(format_string,'(A,I0,A)') '(A,',n_k_cur,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,n_k_cur-1)
      end if
      
c$$$      generating grids for templates and images 
      max_k_c = 2*grid_k_p_(n_r-1)

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we create tesselation'
     $        ,' to access templates efficiently:'
      end if
      n_point = n_S;
      allocate(L_(0:3*n_point-1))
      do npoint=0,n_point-1
         call get_angle_to_vp_(S_alpha_polar_a_(npoint)
     $        ,S_alpha_azimu_b_(npoint),L_(0+3*npoint))
         call normalize_r8(3,L_(0 + 3*npoint))
      enddo !do npoint=0,n_point-1
      if (verbose.gt.1) then
         do npoint=0,n_point-1
            write(6,'(A,F8.4,F8.4,A,A,I0,1X,3F8.4)') 'S_alpha: ' ,
     $           S_alpha_polar_a_(npoint) , S_alpha_azimu_b_(npoint) ,
     $           ' --> ' , ' L_: ', npoint , L_(0 + 3*npoint), L_(1 + 3
     $           *npoint) , L_(2 + 3*npoint)
         enddo ! do npoint=0,n_point-1
      end if! if (verbose.gt.1) then
      call tesselation_get_nl_nm_ll(n_point,L_,tradius_min,nl_max
     $     ,nm_sum,ll_sum)
      allocate(T_nl_(0:1*nm_sum-1)) !level
      allocate(T_vm_(0:3*nm_sum-1)) !vertex center
      allocate(T_tr_(0:1*nm_sum-1)) !tradius
      allocate(T_ll_(0:1*nm_sum-1)) !number of points from L_ in T_
      allocate(T_lf_(0:1*nm_sum-1)) !is leaf
      allocate(T_c0_(0:1*nm_sum-1)) !child_0 tesselation_index
      allocate(T_c1_(0:1*nm_sum-1)) !child_1 tesselation_index
      allocate(T_c2_(0:1*nm_sum-1)) !child_2 tesselation_index
      allocate(T_c3_(0:1*nm_sum-1)) !child_3 tesselation_index
      allocate(T_ls_(0:1*nm_sum-1)) !starting index of point_index_list for T_ if leaf (leaves only)
      allocate(T_LT_(0:1*ll_sum-1)) !full point_index_list for all of T_ (leaves only)
      allocate(T_root_base_(0:8-1))
      call tesselation_index_wrapper_0(n_point,L_,tradius_min ,nl_max
     $     ,nm_sum,ll_sum,T_nl_,T_vm_,T_tr_,T_ll_,T_lf_,T_c0_ ,T_c1_
     $     ,T_c2_,T_c3_,T_ls_,T_LT_,n_T_root_base,T_root_base_)

c$$$      reading in estimated translations and rotations
      allocate(polar_a_est_(0:n_M-1))
      allocate(azimu_b_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      allocate(ctf_ind_est_(0:n_M-1))
      allocate(S_index_est_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_est_(nm) = alpha_est_(nalpha_polar_a + nm*n_alpha)
         azimu_b_est_(nm) = alpha_est_(nalpha_azimu_b + nm*n_alpha)
         gamma_z_est_(nm) = alpha_est_(nalpha_gamma_z + nm*n_alpha)
         delta_x_est_(nm) = alpha_est_(nalpha_delta_x + nm*n_alpha)
         delta_y_est_(nm) = alpha_est_(nalpha_delta_y + nm*n_alpha)
         l2_norm_est_(nm) = alpha_est_(nalpha_l2_norm + nm*n_alpha)
         ctf_ind_est_(nm) = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         S_index_est_(nm) = alpha_est_(nalpha_S_index + nm*n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->',polar_a_est_(nm)
     $           ,azimu_b_est_(nm),gamma_z_est_(nm),delta_x_est_(nm)
     $           ,delta_y_est_(nm),l2_norm_est_(nm),ctf_ind_est_(nm)
     $           ,S_index_est_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1

c$$$      preparing updated translations and rotations
      allocate(polar_a_upd_(0:n_M-1))
      allocate(azimu_b_upd_(0:n_M-1))
      allocate(gamma_z_upd_(0:n_M-1))
      allocate(delta_x_upd_(0:n_M-1))
      allocate(delta_y_upd_(0:n_M-1))
      allocate(l2_norm_upd_(0:n_M-1))
      allocate(ctf_ind_upd_(0:n_M-1))
      allocate(S_index_upd_(0:n_M-1))
      do nm=0,n_M-1
         polar_a_upd_(nm) = alpha_est_(nalpha_polar_a + nm*n_alpha)
         azimu_b_upd_(nm) = alpha_est_(nalpha_azimu_b + nm*n_alpha)
         gamma_z_upd_(nm) = alpha_est_(nalpha_gamma_z + nm*n_alpha)
         delta_x_upd_(nm) = alpha_est_(nalpha_delta_x + nm*n_alpha)
         delta_y_upd_(nm) = alpha_est_(nalpha_delta_y + nm*n_alpha)
         l2_norm_upd_(nm) = alpha_est_(nalpha_l2_norm + nm*n_alpha)
         ctf_ind_upd_(nm) = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         S_index_upd_(nm) = alpha_est_(nalpha_S_index + nm*n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->',polar_a_upd_(nm)
     $           ,azimu_b_upd_(nm),gamma_z_upd_(nm),delta_x_upd_(nm)
     $           ,delta_y_upd_(nm),l2_norm_upd_(nm),ctf_ind_upd_(nm)
     $           ,S_index_upd_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_fast_omp_stage_0'
     $        ,' to calculate all innerproducts:'
      end if !if (verbose.gt.1) then
c$$$      ld_S is the stride between adjacent templates in S_k_q_
c$$$      ld_M is the stride between adjacent images in M_k_q_
c$$$      ld_CTF is the stride between adjacent ctfs in CTF_k_q_
c$$$      calculating all innerproducts
      timing_tic = omp_get_wtime()
      call test_innerproduct_fast_omp_stage_0(verbose,n_omp_sub
     $     ,svd_calculation_type,eps_svd,n_r,grid_k_p_,n_w_
     $     ,half_diameter_x_c ,n_S,ld_S ,S_k_p_,tesselation_distance_req
     $     ,S_alpha_polar_a_ ,S_alpha_azimu_b_,L_ ,nl_max,nm_sum ,ll_sum
     $     ,T_nl_,T_vm_,T_tr_ ,T_ll_ ,T_lf_,T_c0_ ,T_c1_,T_c2_ ,T_c3_
     $     ,T_ls_,T_LT_ ,n_T_root_base ,T_root_base_ ,n_M,ld_M ,M_k_p_
     $     ,n_CTF,ld_CTF ,CTF_k_p_ ,polar_a_est_ ,azimu_b_est_
     $     ,gamma_z_est_ ,delta_x_est_ ,delta_y_est_ ,l2_norm_est_
     $     ,ctf_ind_est_ ,S_index_est_ ,polar_a_upd_ ,azimu_b_upd_
     $     ,gamma_z_upd_ ,delta_x_upd_ ,delta_y_upd_ ,l2_norm_upd_
     $     ,ctf_ind_upd_ ,S_index_upd_ ,alpha_update_f ,n_pixels_in
     $     ,displacement_max ,n_delta_x ,n_delta_y ,n_gamma_z
     $     ,fpm_howmany_max,C_M_ ,n_SM_use)
      do nm=0,n_M-1
         alpha_est_(nalpha_polar_a + nm*n_alpha) = polar_a_upd_(nm)
         alpha_est_(nalpha_azimu_b + nm*n_alpha) = azimu_b_upd_(nm)
         alpha_est_(nalpha_gamma_z + nm*n_alpha) = gamma_z_upd_(nm)
         alpha_est_(nalpha_delta_x + nm*n_alpha) = delta_x_upd_(nm)
         alpha_est_(nalpha_delta_y + nm*n_alpha) = delta_y_upd_(nm)
         alpha_est_(nalpha_l2_norm + nm*n_alpha) = l2_norm_upd_(nm)
         alpha_est_(nalpha_ctf_ind + nm*n_alpha) = ctf_ind_upd_(nm)
         alpha_est_(nalpha_S_index + nm*n_alpha) = S_index_upd_(nm)
         if (verbose.gt.2) then
            write(6,'(I2,A,8F8.3)') nm
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns -->',polar_a_upd_(nm)
     $           ,azimu_b_upd_(nm),gamma_z_upd_(nm),delta_x_upd_(nm)
     $           ,delta_y_upd_(nm),l2_norm_upd_(nm),ctf_ind_upd_(nm)
     $           ,S_index_upd_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,n_M-1
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3,A,F8.6,A,F8.7)')
     $        'test_innerproduct_fast_omp_stage_0:',' total_time: '
     $        ,timing_toc-timing_tic,'; time_per S,M: ',(timing_toc
     $        -timing_tic)/max(1,n_SM_use),'; time_per S,M,dx,dy: '
     $        ,(timing_toc -timing_tic)/(n_SM_use*n_delta_x *n_delta_y)
      end if

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_fast_wrapper_0]'
      end if

      deallocate(ngridc_)
      deallocate(icstart_)
      deallocate(n_w_)

      deallocate(L_)
      deallocate(T_nl_) !level
      deallocate(T_vm_) !vertex center
      deallocate(T_tr_) !tradius
      deallocate(T_ll_) !number of points from L_ in T_
      deallocate(T_lf_) !is leaf
      deallocate(T_c0_) !child_0 tesselation_index
      deallocate(T_c1_) !child_1 tesselation_index
      deallocate(T_c2_) !child_2 tesselation_index
      deallocate(T_c3_) !child_3 tesselation_index
      deallocate(T_ls_) !starting index of point_index_list for T_ if leaf (leaves only)
      deallocate(T_LT_) !full point_index_list for all of T_ (leaves only)
      deallocate(T_root_base_)

      deallocate(polar_a_est_)
      deallocate(azimu_b_est_)
      deallocate(gamma_z_est_)
      deallocate(delta_x_est_)
      deallocate(delta_y_est_)
      deallocate(l2_norm_est_)
      deallocate(ctf_ind_est_)
      deallocate(S_index_est_)
      deallocate(polar_a_upd_)
      deallocate(azimu_b_upd_)
      deallocate(gamma_z_upd_)
      deallocate(delta_x_upd_)
      deallocate(delta_y_upd_)
      deallocate(l2_norm_upd_)
      deallocate(ctf_ind_upd_)
      deallocate(S_index_upd_)

      end
