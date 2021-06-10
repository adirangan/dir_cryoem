      subroutine test_innerproduct_batch_wrapper_1(verbose,ncur,nlats_
     $     ,grid_k_p_,max_x_c,n_S,ld_S,S_k_p_,n_M,ld_M,M_k_p_,alpha_est_
     $     ,eps_target,N_pixels_in,n_delta_x,n_delta_y,n_gamma_z
     $     ,svd_calculation_type,delta_x_sort_,delta_y_sort_
     $     ,gamma_z_sort_,C_Z_sort_,I_permute_,I_inverse_)
c$$$      Reads in a few templates and images and compares them. The
c$$$      output (i.e., which innerproducts are highest) is returned.
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer *4 ncur,nlats_(0:ncur-1),n_S,ld_S,n_M,ld_M
      integer *4 n_delta_x,n_delta_y,n_gamma_z,svd_calculation_type
      real *8 grid_k_p_(0:0),max_x_c,alpha_est_(0:0)
      real *8 eps_target,N_pixels_in
      complex *16 S_k_p_(0:0),M_k_p_(0:0)
      real *8 delta_x_sort_(0:n_S*n_M-1)
      real *8 delta_y_sort_(0:n_S*n_M-1)
      real *8 gamma_z_sort_(0:n_S*n_M-1)
      complex *16 C_Z_sort_(0:n_S*n_M-1)
      integer I_permute_(0:n_S*n_M-1)
      integer I_inverse_(0:n_S*n_M-1)
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
c$$$      range of displacement and rotation parameters to measure
      integer ndx,ndy,ngz
      real *8 delta_x,delta_y,gamma_z
      real *8, allocatable :: delta_x_(:)
      real *8, allocatable :: delta_y_(:)
      real *8, allocatable :: gamma_z_(:)      
c$$$      array of innerproducts to hold measurements
      integer n_C,nC
      real *8 Clim
      complex *16 C_S,C_M,C_Z
      real *8, allocatable :: delta_x_max_(:)
      real *8, allocatable :: delta_y_max_(:)
      real *8, allocatable :: gamma_z_max_(:)
      complex *16, allocatable :: C_Z_max_(:)
      character(len=64) format_string
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_batch_wrapper_1]: '
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
         write(format_string,'(A,I0,A)') '(',5,'(F8.3,1X))'
         write(6,'(A)') 'alpha_est_: '
         write(6,format_string) (alpha_est_(nr),nr=0,5*n_M-1)
         write(6,'(A,F6.3)') 'eps_target: ',eps_target
         write(6,'(A,F6.3)') 'N_pixels_in: ',N_pixels_in
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
      do nm=0,n_M-1
         gamma_z_est_(nm) = alpha_est_(2 + nm*5)
         delta_x_est_(nm) = alpha_est_(3 + nm*5)
         delta_y_est_(nm) = alpha_est_(4 + nm*5)
         if (verbose.gt.1) then
            write(6,'(I0,A,3F8.3)') nm,'<-- nm: dx,dy,gz-->'
     $           ,delta_x_est_(nm),delta_y_est_(nm),gamma_z_est_(nm)
         end if
      enddo

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_svdread_batch'
     $        ,' to calculate all innerproducts:'
      end if
c$$$      ld_S is the stride between adjacent templates in S_k_q_
c$$$      ld_M is the stride between adjacent images in M_k_q_
c$$$      calculating all innerproducts
      allocate(delta_x_max_(0:n_S*n_M-1))
      allocate(delta_y_max_(0:n_S*n_M-1))
      allocate(gamma_z_max_(0:n_S*n_M-1))
      allocate(C_Z_max_(0:n_S*n_M-1))
      timing_tic = second()
      call test_innerproduct_svdread_batch(verbose,svd_calculation_type
     $     ,n_r,grid_k_p_,n_w_,max_x_c,n_S,ld_S,S_k_p_,n_M,ld_M,M_k_p_
     $     ,delta_x_est_,delta_y_est_,gamma_z_est_,eps_target
     $     ,N_pixels_in,n_delta_x,n_delta_y,n_gamma_z,delta_x_max_
     $     ,delta_y_max_,gamma_z_max_,C_Z_max_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.3,A,F8.6,A,F8.7)')
     $        'test_innerproduct_svdread_batch:',' total_time: '
     $        ,timing_toc-timing_tic,'; time_per S,M: ',(timing_toc
     $        -timing_tic)/(n_S*n_M),'; time_per S,M,dx,dy: '
     $        ,(timing_toc -timing_tic)/(n_S*n_M*n_delta_x*n_delta_y)
      end if
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_max_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_max^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_max_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_batch_sort'
     $        ,' to sort the results:'
      end if
      call test_innerproduct_batch_sort(n_S,n_M,delta_x_max_
     $     ,delta_y_max_,gamma_z_max_,C_Z_max_,delta_x_sort_
     $     ,delta_y_sort_,gamma_z_sort_,C_Z_sort_,I_permute_,I_inverse_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_sort_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_permute_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_permute_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_inverse_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_inverse_(na),na=0,n_S*n_M-1)
      end if


      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_1]'
      end if

      deallocate(ngridc_)
      deallocate(icstart_)
      deallocate(n_w_)
      deallocate(delta_x_est_)
      deallocate(delta_y_est_)
      deallocate(gamma_z_est_)
      deallocate(delta_x_max_)
      deallocate(delta_y_max_)
      deallocate(gamma_z_max_)
      deallocate(C_Z_max_)

      end
