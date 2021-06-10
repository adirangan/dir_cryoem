      subroutine test_innerproduct_batch_wrapper_4(verbose,n_S,n_M
     $    ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $    ,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_,delta_x_sort_MS_
     $    ,delta_y_sort_MS_,gamma_z_sort_MS_,C_Z_sort_MS_,I_permute_MS_
     $    ,I_inverse_MS_,alpha_0_,alpha_1_,alpha_est_)
c$$$      Update alphas
      implicit none
      integer verbose
      integer *4 n_S,n_M
      real *8 delta_x_sort_SM_(0:n_S*n_M-1)
      real *8 delta_y_sort_SM_(0:n_S*n_M-1)
      real *8 gamma_z_sort_SM_(0:n_S*n_M-1)
      complex *16 C_Z_sort_SM_(0:n_S*n_M-1)
      integer I_permute_SM_(0:n_S*n_M-1)
      integer I_inverse_SM_(0:n_S*n_M-1)
      real *8 delta_x_sort_MS_(0:n_S*n_M-1)
      real *8 delta_y_sort_MS_(0:n_S*n_M-1)
      real *8 gamma_z_sort_MS_(0:n_S*n_M-1)
      complex *16 C_Z_sort_MS_(0:n_S*n_M-1)
      integer I_permute_MS_(0:n_S*n_M-1)
      integer I_inverse_MS_(0:n_S*n_M-1)
      real *8 alpha_0_(0:n_S-1),alpha_1_(0:n_S-1)
      real *8 alpha_est_(0:5*n_M-1)
c$$$      temporary arrays
      integer *4 I_S_(0:n_S-1)
      logical B_M_(0:n_M-1),bm_use
      integer ns,ns_tmp,nm,nm_sum,nm_use
      character(len=64) format_string
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_batch_wrapper_4]'
      end if
      if (verbose.gt.0) then
         write(6,'(A)') 'update alphas'
     $        ,' using best image for each template'
     $        ,' assuming n_S<n_M'
      end if
      call randperm(n_S,I_S_)
      do nm=0,n_M-1
         B_M_(nm) = .false.
      enddo
      nm_sum = 0
      ns_tmp = 0
      do while (nm_sum.lt.n_M-1)
         ns = I_S_(ns_tmp)
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0)') 'nm_sum ',nm_sum,'; ns_tmp ',ns_tmp
     $           ,'; ns ',ns
         end if
         nm_use = n_M-1
         continue_flag=.true.
         do while (continue_flag)
            nm = I_permute_MS_(nm_use + ns*n_M)
            bm_use = B_M_(nm)
            if (bm_use) then
               if (verbose.gt.0) then
                  write(6,'(A,I0,A,I0,A)') ' nm_use ',nm_use,'; nm ',nm
     $                 ,'; skipping'
               end if
               nm_use = nm_use - 1
            else ! if unused
               if (verbose.gt.0) then
                  write(6,'(A,I0,A,I0,A)') ' nm_use ',nm_use,'; nm ',nm
     $                 ,'; updating'
               end if
               bm_use = .true.
               B_M_(nm) = .true.
               continue_flag = .false.
               nm_sum = nm_sum + 1
               alpha_est_(0 + nm*5) = alpha_0_(ns)
               alpha_est_(1 + nm*5) = alpha_1_(ns)
               alpha_est_(2 + nm*5) = alpha_est_(2 + nm*5) +
     $              gamma_z_sort_MS_(nm_use + ns*n_M)
               alpha_est_(3 + nm*5) = alpha_est_(3 + nm*5) +
     $              delta_x_sort_MS_(nm_use + ns*n_M)
               alpha_est_(4 + nm*5) = alpha_est_(4 + nm*5) +
     $              delta_y_sort_MS_(nm_use + ns*n_M)
            end if ! check bm_used
         end do ! continue_flag
         ns_tmp = ns_tmp + 1
         if (ns_tmp.ge.n_S) then
            ns_tmp = 0
         end if
      end do

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_batch_wrapper_4]'
      end if
      end
