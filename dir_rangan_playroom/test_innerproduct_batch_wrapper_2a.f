      subroutine test_innerproduct_batch_wrapper_2a(verbose,ncur,nlats_
     $     ,grid_k_p_,max_x_c,n_S,ld_S,S_k_p_,n_M,ld_M,M_k_p_,n_CTF
     $     ,ld_CTF ,CTF_k_p_,n_alpha,alpha_est_,eps_target,N_pixels_in
     $     ,displacement_max,n_delta_x,n_delta_y,n_gamma_z
     $     ,svd_calculation_type,delta_x_max_,delta_y_max_,gamma_z_max_
     $     ,C_M_,C_S_max_,C_Z_max_)
c$$$      Reads in a few templates and images and compares them. The
c$$$      output (i.e., which innerproducts are highest) is returned.
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer *4 ncur,nlats_(0:ncur-1),n_S,ld_S,n_M,ld_M,n_CTF,ld_CTF
     $     ,n_alpha
      integer *4 n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type
      real *8 grid_k_p_(0:0),max_x_c,alpha_est_(0:0)
      include 'nalpha_define.f'
      real *8 eps_target,N_pixels_in,displacement_max
      complex *16 S_k_p_(0:0),M_k_p_(0:0),CTF_k_p_(0:0)
      real *8 delta_x_max_(0:n_S*n_M-1)
      real *8 delta_y_max_(0:n_S*n_M-1)
      real *8 gamma_z_max_(0:n_S*n_M-1)
      complex *16 C_M_(0:n_M-1)
      complex *16 C_S_max_(0:n_S*n_M-1)
      complex *16 C_Z_max_(0:n_S*n_M-1)
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
      real *8, allocatable :: delta_x_est_(:)
      real *8, allocatable :: delta_y_est_(:)
      real *8, allocatable :: gamma_z_est_(:)
c$$$      estimated l2_norm and image-dependent ctf_index
      real *8, allocatable :: l2_norm_est_(:)
      real *8, allocatable :: ctf_ind_est_(:)
c$$$      range of displacement and rotation parameters to measure
      integer ndx,ndy,ngz
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)      
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16, allocatable :: C_S_(:)
      complex *16 C_M,C_Z
      character(len=64) format_string
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_batch_wrapper_2a]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'ncur: ',ncur
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'nlats_: ',(nlats_(nr),nr=0,ncur-1)
         write(format_string,'(A,I0,A)') '(A,',ncur,'(F8.3,1X))'
         write(6,format_string) 'grid_k_p_: ',(grid_k_p_(nr),nr=0,ncur
     $        -1)
         write(6,'(A,F6.3)') 'max_x_c: ',max_x_c
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'ld_S: ',ld_S
         write(6,'(A,I0)') 'n_M: ',n_M
         write(6,'(A,I0)') 'ld_M: ',ld_M
         write(6,'(A,I0)') 'n_CTF: ',n_CTF
         write(6,'(A,I0)') 'ld_CTF: ',ld_CTF
         write(6,'(A,I0)') 'n_alpha: ',n_alpha
         write(format_string,'(A,I0,A)') '(',n_alpha,'(F8.3,1X))'
         write(6,'(A)') 'alpha_est_: '
         write(6,format_string) (alpha_est_(nr),nr=0 ,n_alpha*n_M-1)
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,F6.3)') 'N_pixels_in: ',N_pixels_in
         write(6,'(A,F6.3)') 'displacement_max: ',displacement_max
         write(6,'(A,I0)') 'n_delta_x: ',n_delta_x
         write(6,'(A,I0)') 'n_delta_y: ',n_delta_y
         write(6,'(A,I0)') 'n_gamma_z: ',n_gamma_z
         write(6,'(A,I0)') 'svd_calculation_type: ',svd_calculation_type
      end if

      pi = 4*atan(1.0)

c$$$      Calculating template size using 'get_template_size'
      allocate(ngridc_(ncur))
      allocate(icstart_(ncur))
      if (verbose.gt.1) then
         write(6,'(A,I0,A,A)') 'ncur = ',ncur
     $        ,'; calling get_template_size to'
     $        ,' determine ngridc_ and ntemplatesize'
      end if
      call get_template_size(nlats_,ncur,ntemplatesize,ngridc_,icstart_)
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'ntemplatesize = ',ntemplatesize
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'ngridc_: ',(ngridc_(nr),nr=1,ncur)
      end if
      
c$$$      indices
      n_r = ncur
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
         write(format_string,'(A,I0,A)') '(A,',ncur,'(I0,1X))'
         write(6,format_string) 'n_w_: ',(n_w_(nr),nr=0,ncur-1)
      end if
      
c$$$      generating grids for templates and images 
      max_k_c = 2*grid_k_p_(n_r-1)

c$$$      reading in estimated translations and rotations
      allocate(delta_x_est_(0:n_M-1))
      allocate(delta_y_est_(0:n_M-1))
      allocate(gamma_z_est_(0:n_M-1))
      allocate(l2_norm_est_(0:n_M-1))
      allocate(ctf_ind_est_(0:n_M-1))
      do nm=0,n_M-1
         gamma_z_est_(nm) = alpha_est_(nalpha_gamma_z + nm*n_alpha)
         delta_x_est_(nm) = alpha_est_(nalpha_delta_x + nm*n_alpha)
         delta_y_est_(nm) = alpha_est_(nalpha_delta_y + nm*n_alpha)
         l2_norm_est_(nm) = alpha_est_(nalpha_l2_norm + nm*n_alpha)
         ctf_ind_est_(nm) = alpha_est_(nalpha_ctf_ind + nm*n_alpha)
         if (verbose.gt.1) then
            write(6,'(I0,A,5F8.3)') nm,'<-- nm: dx,dy,gz,l2,ci-->'
     $           ,delta_x_est_(nm),delta_y_est_(nm),gamma_z_est_(nm)
     $           ,l2_norm_est_(nm),ctf_ind_est_(nm)
         end if
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_svdread_batch'
     $        ,' to calculate all innerproducts:'
      end if
c$$$      ld_S is the stride between adjacent templates in S_k_q_
c$$$      ld_M is the stride between adjacent images in M_k_q_
c$$$      ld_CTF is the stride between adjacent ctfs in CTF_k_q_
c$$$      calculating all innerproducts
      timing_tic = second()
      call test_innerproduct_batch_stage_0a(verbose
     $     ,svd_calculation_type,n_r,grid_k_p_,n_w_,max_x_c,n_S,ld_S
     $     ,S_k_p_,n_M,ld_M,M_k_p_,n_CTF,ld_CTF,CTF_k_p_,delta_x_est_
     $     ,delta_y_est_,gamma_z_est_,l2_norm_est_,ctf_ind_est_
     $     ,eps_target,N_pixels_in,displacement_max,n_delta_x,n_delta_y
     $     ,n_gamma_z,delta_x_max_,delta_y_max_,gamma_z_max_,C_M_
     $     ,C_S_max_ ,C_Z_max_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3,A,F8.6,A,F8.7)')
     $        'test_innerproduct_batch_stage_0:',' total_time: '
     $        ,timing_toc-timing_tic,'; time_per S,M: ',(timing_toc
     $        -timing_tic)/(n_S*n_M),'; time_per S,M,dx,dy: '
     $        ,(timing_toc -timing_tic)/(n_S*n_M*n_delta_x*n_delta_y)
      end if
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_max_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_max_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_max_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_S_max_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_S_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_max_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_max_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_2a]'
      end if

      deallocate(ngridc_)
      deallocate(icstart_)
      deallocate(n_w_)
      deallocate(delta_x_est_)
      deallocate(delta_y_est_)
      deallocate(gamma_z_est_)
      deallocate(l2_norm_est_)
      deallocate(ctf_ind_est_)

      end
