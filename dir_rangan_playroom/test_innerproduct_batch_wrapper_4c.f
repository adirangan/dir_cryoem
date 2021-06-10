      subroutine test_innerproduct_batch_wrapper_4c(verbose,n_S,n_M,C_M_
     $     ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $     ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $     ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $     ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_
     $     ,alpha_polar_a_,alpha_azimu_b_ ,n_alpha,alpha_est_)
c$$$      Update alphas
      implicit none
      integer verbose
      integer *4 n_S,n_M,n_alpha
      include 'nalpha_define.f'
      complex *16 C_M_(0:n_M-1)
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
      real *8 alpha_polar_a_(0:n_S-1),alpha_azimu_b_(0:n_S-1)
      real *8 alpha_est_(0:n_alpha*n_M-1)
      real *8 delta_x_tmp,delta_y_tmp
c$$$      temporary arrays
      integer *4, allocatable :: I_S_(:)
      logical, allocatable :: B_M_(:)
      logical bm_use,continue_flag
      integer ns,ns_tmp,nm,nm_sum,nm_use
      character(len=64) format_string
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_wrapper_4c]'
      end if
      if (verbose.gt.0) then
         write(6,'(A)') 'update alphas'
     $        ,' using best image for each template'
     $        ,' assuming n_S<n_M'
      end if

      allocate(I_S_(0:n_S-1))
      allocate(B_M_(0:n_M-1))
      call randperm(n_S,I_S_)
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_S,'(I0,1X))'
         write(6,format_string) 'I_S_: ',(I_S_(ns),ns=0,n_S-1)
      end if
      do nm=0,n_M-1
         B_M_(nm) = .false.
      enddo
      nm_sum = 0
      ns_tmp = 0
      do while (nm_sum.le.n_M-1)
         ns = I_S_(ns_tmp)
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0,A,I0)') 'nm_sum ',nm_sum,'; ns_tmp ',
     $           ns_tmp,'; ns ',ns
         end if
         nm_use = n_M-1
         continue_flag=.true.
         do while (continue_flag)
            nm = I_permute_MS_(nm_use + ns*n_M)
            bm_use = B_M_(nm)
            if (bm_use) then
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,I0,A)') ' nm_use ',nm_use,'; nm ',nm
     $                 ,'; skipping'
               end if
               nm_use = nm_use - 1
            else ! if unused
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,I0,A)') ' nm_use ',nm_use,'; nm ',nm
     $                 ,'; updating'
               end if
               bm_use = .true.
               B_M_(nm) = .true.
               continue_flag = .false.
               nm_sum = nm_sum + 1
               call get_interchange_delta(delta_x_sort_MS_(nm_use + ns
     $              *n_M),delta_y_sort_MS_(nm_use + ns*n_M)
     $              ,alpha_est_(nalpha_gamma_z + nm*n_alpha),delta_x_tmp
     $              ,delta_y_tmp)
               alpha_est_(nalpha_polar_a + nm*n_alpha) =
     $              alpha_polar_a_(ns)
               alpha_est_(nalpha_azimu_b + nm*n_alpha) =
     $              alpha_azimu_b_(ns)
               alpha_est_(nalpha_gamma_z + nm*n_alpha) =
     $              alpha_est_(nalpha_gamma_z + nm*n_alpha)
     $              +gamma_z_sort_MS_(nm_use + ns*n_M)
               alpha_est_(nalpha_delta_x + nm*n_alpha) =
     $              alpha_est_(nalpha_delta_x + nm*n_alpha)
     $              +delta_x_tmp
               alpha_est_(nalpha_delta_y + nm*n_alpha) =
     $              alpha_est_(nalpha_delta_y + nm*n_alpha)
     $              +delta_y_tmp
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,I0,A,2F8.3,A,2F8.3,A,2F8.3)') 'nm '
     $                 ,nm,'; nm_use ',nm_use,'; C_S '
     $                 ,C_S_sort_MS_(nm_use + ns*n_M),'; C_M ',C_M_(nm)
     $                 ,'; C_Z ', C_Z_sort_MS_(nm_use + ns *n_M)
               end if
               if (zabs(C_S_sort_MS_(nm_use + ns*n_M)).le.1.0d-15) then
                  alpha_est_(nalpha_l2_norm + nm*n_alpha) = 1.0d0
               else 
                  alpha_est_(nalpha_l2_norm + nm*n_alpha) =
     $                 zabs(C_M_(nm)) * zabs(C_Z_sort_MS_(nm_use + ns
     $                 *n_M)) / zabs(C_S_sort_MS_(nm_use + ns*n_M))
               end if
            end if ! check bm_used
         end do ! continue_flag
         ns_tmp = ns_tmp + 1
         if (ns_tmp.ge.n_S) then
            ns_tmp = 0
         end if
      end do ! while (nm_sum.lt.n_M-1)
      ns = I_S_(ns_tmp)
      if (verbose.gt.0) then
         write(6,'(A,I0,A,I0,A,I0,A)') 'nm_sum ',nm_sum,'; ns_tmp '
     $        ,ns_tmp,'; ns ',ns,'; finished'
      end if

      deallocate(I_S_)
      deallocate(B_M_)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_4c]'
      end if
      end
