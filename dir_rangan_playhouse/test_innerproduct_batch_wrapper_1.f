      subroutine test_innerproduct_batch_wrapper_1(verbose,n_S,n_M
     $     ,delta_x_max_,delta_y_max_,gamma_z_max_,C_S_max_,C_Z_max_
     $     ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $     ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
c$$$      sorts output of test_innerproduct_batch.
      implicit none
      integer verbose
      integer *4 n_S,n_M
      real *8 delta_x_max_(0:n_S*n_M-1)
      real *8 delta_y_max_(0:n_S*n_M-1)
      real *8 gamma_z_max_(0:n_S*n_M-1)
      complex *16 C_S_max_(0:n_S*n_M-1)
      complex *16 C_Z_max_(0:n_S*n_M-1)
      real *8 delta_x_sort_SM_(0:n_S*n_M-1)
      real *8 delta_y_sort_SM_(0:n_S*n_M-1)
      real *8 gamma_z_sort_SM_(0:n_S*n_M-1)
      complex *16 C_S_sort_SM_(0:n_S*n_M-1)
      complex *16 C_Z_sort_SM_(0:n_S*n_M-1)
      integer I_permute_SM_(0:n_S*n_M-1)
      integer I_inverse_SM_(0:n_S*n_M-1)
      real *8 delta_x_sort_MS_(0:n_S*n_M-1)
      real *8 delta_y_sort_MS_(0:n_S*n_M-1)
      real *8 gamma_z_sort_MS_(0:n_S*n_M-1)
      complex *16 C_S_sort_MS_(0:n_S*n_M-1)
      complex *16 C_Z_sort_MS_(0:n_S*n_M-1)
      integer I_permute_MS_(0:n_S*n_M-1)
      integer I_inverse_MS_(0:n_S*n_M-1)
c$$$      temporary arrays
      real *8, allocatable :: delta_x_max_SM_(:)
      real *8, allocatable :: delta_y_max_SM_(:)
      real *8, allocatable :: gamma_z_max_SM_(:)
      complex *16, allocatable :: C_S_max_SM_(:)
      complex *16, allocatable :: C_Z_max_SM_(:)
      real *8, allocatable :: delta_x_max_MS_(:)
      real *8, allocatable :: delta_y_max_MS_(:)
      real *8, allocatable :: gamma_z_max_MS_(:)
      complex *16, allocatable :: C_S_max_MS_(:)
      complex *16, allocatable :: C_Z_max_MS_(:)
      character(len=64) format_string
      integer na
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
          write(6,'(A)')
     $        '[entering test_innerproduct_batch_wrapper_1]: '
       end if
       if (verbose.gt.1) then
         write(6,'(A,I0)') 'verbose: ',verbose
         write(6,'(A,I0)') 'n_S: ',n_S
         write(6,'(A,I0)') 'n_M: ',n_M
      end if

      allocate(delta_x_max_SM_(0:n_S*n_M-1))
      allocate(delta_y_max_SM_(0:n_S*n_M-1))
      allocate(gamma_z_max_SM_(0:n_S*n_M-1))
      allocate(C_S_max_SM_(0:n_S*n_M-1))
      allocate(C_Z_max_SM_(0:n_S*n_M-1))
      allocate(delta_x_max_MS_(0:n_S*n_M-1))
      allocate(delta_y_max_MS_(0:n_S*n_M-1))
      allocate(gamma_z_max_MS_(0:n_S*n_M-1))
      allocate(C_S_max_MS_(0:n_S*n_M-1))
      allocate(C_Z_max_MS_(0:n_S*n_M-1))
      call cp1_r8(n_S*n_M,delta_x_max_,delta_x_max_SM_)
      call cp1_r8(n_S*n_M,delta_y_max_,delta_y_max_SM_)
      call cp1_r8(n_S*n_M,gamma_z_max_,gamma_z_max_SM_)
      call cp1_c16(n_S*n_M,C_S_max_,C_S_max_SM_)
      call cp1_c16(n_S*n_M,C_Z_max_,C_Z_max_SM_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_max_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_max_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_max_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_max_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_max_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_max_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_S_max_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_S_max_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_max_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_max_SM_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call trn0_r8 to transpose results:'
      end if
      call trn0_r8(n_S,n_M,delta_x_max_SM_,delta_x_max_MS_)
      call trn0_r8(n_S,n_M,delta_y_max_SM_,delta_y_max_MS_)
      call trn0_r8(n_S,n_M,gamma_z_max_SM_,gamma_z_max_MS_)
      call trn0_c16(n_S,n_M,C_S_max_SM_,C_S_max_MS_)
      call trn0_c16(n_S,n_M,C_Z_max_SM_,C_Z_max_MS_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_max_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (delta_x_max_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_max_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (delta_y_max_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_max_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (gamma_z_max_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_S_max_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M*2,'(F8.3,1X))'
         write(6,format_string) (C_S_max_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_max_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_max_MS_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') ''
         write(6,*) 'Now we call test_innerproduct_batch_sort'
     $        ,' to sort the results:'
      end if
      call test_innerproduct_batch_sort_a(n_S,n_M,delta_x_max_SM_
     $     ,delta_y_max_SM_,gamma_z_max_SM_,C_S_max_SM_,C_Z_max_SM_
     $     ,delta_x_sort_SM_ ,delta_y_sort_SM_,gamma_z_sort_SM_
     $     ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_ ,I_inverse_SM_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_sort_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_x_sort_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_sort_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (delta_y_sort_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_sort_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
         write(6,format_string) (gamma_z_sort_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_S_sort_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_S_sort_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_sort_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_sort_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_permute_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_permute_SM_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_inverse_SM_^T: '
         write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
         write(6,format_string) (I_inverse_SM_(na),na=0,n_S*n_M-1)
      end if
      call test_innerproduct_batch_sort_a(n_M,n_S,delta_x_max_MS_
     $     ,delta_y_max_MS_,gamma_z_max_MS_,C_S_max_MS_,C_Z_max_MS_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_sort_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (delta_x_sort_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'delta_y_sort_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (delta_y_sort_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'gamma_z_sort_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(F8.3,1X))'
         write(6,format_string) (gamma_z_sort_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_S_sort_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M*2,'(F8.3,1X))'
         write(6,format_string) (C_S_sort_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'C_Z_sort_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_sort_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_permute_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(I8,1X))'
         write(6,format_string) (I_permute_MS_(na),na=0,n_S*n_M-1)
         write(6,'(A)') 'I_inverse_MS_^T: '
         write(format_string,'(A,I0,A)') '(',n_M,'(I8,1X))'
         write(6,format_string) (I_inverse_MS_(na),na=0,n_S*n_M-1)
      end if

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_1]'
      end if

      deallocate(delta_x_max_SM_)
      deallocate(delta_y_max_SM_)
      deallocate(gamma_z_max_SM_)
      deallocate(C_S_max_SM_)
      deallocate(C_Z_max_SM_)
      deallocate(delta_x_max_MS_)
      deallocate(delta_y_max_MS_)
      deallocate(gamma_z_max_MS_)
      deallocate(C_S_max_MS_)
      deallocate(C_Z_max_MS_)

      end
