      subroutine test_innerproduct_fast_MS_update_2(verbose,rseed
     $     ,flag_RTRT_vs_RTTR,n_S,n_M,C_M_,delta_x_sort_MS_
     $     ,delta_y_sort_MS_,gamma_z_sort_MS_,C_S_sort_MS_,C_Z_sort_MS_
     $     ,I_permute_MS_ ,I_inverse_MS_,S_alpha_polar_a_
     $     ,S_alpha_azimu_b_ ,gamma_z_est_,delta_x_est_ ,delta_y_est_
     $     ,polar_a_upd_ ,azimu_b_upd_,gamma_z_upd_ ,delta_x_upd_
     $     ,delta_y_upd_ ,l2_norm_upd_,S_index_upd_)
c$$$      Update alphas
      implicit none
      integer verbose
      integer *4 rseed
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      integer *4 n_S,n_M
      include 'nalpha_define.f'
      complex *16 C_M_(0:n_M-1)
      real *8 delta_x_sort_MS_(0:n_M*n_S-1)
      real *8 delta_y_sort_MS_(0:n_M*n_S-1)
      real *8 gamma_z_sort_MS_(0:n_M*n_S-1)
      complex *16 C_S_sort_MS_(0:n_M*n_S-1)
      complex *16 C_Z_sort_MS_(0:n_M*n_S-1)
      integer I_permute_MS_(0:n_M*n_S-1)
      integer I_inverse_MS_(0:n_M*n_S-1)
      real *8 S_alpha_polar_a_(0:n_S-1),S_alpha_azimu_b_(0:n_S-1)
      real *8 gamma_z_est_(0:0)
      real *8 delta_x_est_(0:0)
      real *8 delta_y_est_(0:0)
      real *8 polar_a_upd_(0:0)
      real *8 azimu_b_upd_(0:0)
      real *8 gamma_z_upd_(0:0)
      real *8 delta_x_upd_(0:0)
      real *8 delta_y_upd_(0:0)
      real *8 l2_norm_upd_(0:0)
      real *8 S_index_upd_(0:0)
      real *8 delta_x_tmp,delta_y_tmp
c$$$      temporary arrays
      integer *4, allocatable :: I_S_(:)
      logical, allocatable :: B_M_(:)
      logical bm_use,continue_flag
      integer ns,ns_tmp,nm,nm_sum,nm_use,ns_use,na_use
      character(len=64) format_string
c$$$      parameters for timing
      real *8 timing_tic,timing_toc

      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_MS_update_2]'
      end if
      if (verbose.gt.0) then
         write(6,'(A)') 'update alphas'
     $        ,' using best image for each template'
     $        ,' assuming n_S<n_M'
      end if

      allocate(I_S_(0:n_S-1))
      allocate(B_M_(0:n_M-1))
      call adi_randperm(rseed,n_S,I_S_)
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
         ns_use = ns
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
c$$$               %%%%%%%%%%%%%%%%
               na_use = nm_use + ns_use*n_M
               polar_a_upd_(nm) = S_alpha_polar_a_(ns_use)
               azimu_b_upd_(nm) = S_alpha_azimu_b_(ns_use)
               S_index_upd_(nm) = ns_use
               if (flag_RTRT_vs_RTTR.eqv..true.) then
                  call get_interchange_delta(delta_x_est_(nm)
     $                 ,delta_y_est_(nm),gamma_z_sort_MS_(na_use)
     $                 ,delta_x_tmp ,delta_y_tmp)
                  gamma_z_upd_(nm) = gamma_z_sort_MS_(na_use) +
     $                 gamma_z_est_(nm)
                  delta_x_upd_(nm) = delta_x_sort_MS_(na_use) +
     $                 delta_x_tmp
                  delta_y_upd_(nm) = delta_y_sort_MS_(na_use) +
     $                 delta_y_tmp
               end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
               if (flag_RTRT_vs_RTTR.eqv..false.) then
                  delta_x_tmp = delta_x_est_(nm) +
     $                 delta_x_sort_MS_(na_use)
                  delta_y_tmp = delta_y_est_(nm) +
     $                 delta_y_sort_MS_(na_use)
                  call get_interchange_delta(delta_x_tmp ,delta_y_tmp
     $                 ,gamma_z_sort_MS_(na_use),delta_x_upd_(nm)
     $                 ,delta_y_upd_(nm))
                  gamma_z_upd_(nm) = gamma_z_sort_MS_(na_use) +
     $                 gamma_z_est_(nm)
               end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,I0,A,2F8.3,A,2F8.3,A,2F8.3)')
     $                 'nm ',nm,'; ns_use ',ns_use,'; C_S '
     $                 ,C_S_sort_MS_(na_use),'; C_M ',C_M_(nm) , 
     $                 '; C_Z ', C_Z_sort_MS_(na_use)
               end if !if (verbose.gt.1) then
               if (zabs(C_S_sort_MS_(na_use)).le.1.0d-15) then
                  l2_norm_upd_(nm) = 1.0d0
               else !if (zabs(C_S_sort_MS_(na_use)).le.1.0d-15) then
                  l2_norm_upd_(nm) = zabs(C_M_(nm)) *
     $                 zabs(C_Z_sort_MS_(na_use))
     $                 /zabs(C_S_sort_MS_(na_use))
               end if ! if (zabs(C_S_sort_MS_(na_use)).le.1.0d-15) then
c$$$               %%%%%%%%%%%%%%%%
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
         write(6,'(A)') '[finished test_innerproduct_fast_MS_update_2]'
      end if
      end

