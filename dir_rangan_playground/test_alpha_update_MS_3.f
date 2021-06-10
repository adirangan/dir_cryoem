!> Doxygen comment: ;\n
!> This function updates a stack of image-parameters by matching images to templates. ;\n
!> (This enforces a uniform distribution of viewing angles, which is not what angle-fitting typically does). ;\n
      subroutine test_alpha_update_MS_3(verbose,rseed,n_MS_max,n_MS_,n_S
     $     ,n_M,alpha_MS__,flag_RTRT_vs_RTTR,alpha_est__,alpha_upd__)
c$$$      Update alphas
      implicit none
      integer verbose
      integer *4 rseed
      integer *4 n_MS_max,n_MS_(0:0),n_S,n_M
      real *8 alpha_MS__(0:0)
      logical flag_RTRT_vs_RTTR ! flag indicating whether transformation operations are ordered as: R_{est}T_{est}R_{upd}T_{upd}(S) [flag_RTRT_vs_RTTR.eqv..true.] or R_{est}T_{est}T_{upd}R_{upd}(S) [flag_RTRT_vs_RTTR.eqv..false.] ;
      real *8 alpha_est__(0:0)
      real *8 alpha_upd__(0:0)
      include 'excerpt_define_nalpha.f'
      real *8 delta_x_upd,delta_y_upd,gamma_z_upd
      real *8 delta_x_est,delta_y_est,gamma_z_est
c$$$      temporary arrays
      integer *4, allocatable :: I_S_(:)
      logical, allocatable :: B_M_(:)
      integer *4, allocatable :: n_MS_tmp_(:)
      logical bm,continue_flag
      integer sum_i4_f
      integer ns,ns_tmp,nm,nm_sum,n_MS_sum,nm_srt,na,nb
      character(len=64) format_string

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_alpha_update_MS_3]'
      end if

      if (verbose.gt.2) then
         write(6,'(A,I0)') ' verbose: ' , verbose
         write(6,'(A,I0)') ' rseed: ' , rseed
         write(6,'(A,I0)') ' n_MS_max: ' , n_MS_max
         call print_sub_i4(n_S,n_MS_,' n_MS_: ')
         write(6,'(A,I0)') ' n_S: ' , n_S
         write(6,'(A,I0)') ' n_M: ' , n_M
         call print_sub_r8(n_alpha*n_MS_max*n_S,alpha_MS__
     $        ,' alpha_MS__: ')
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' , flag_RTRT_vs_RTTR
         do nM=0,n_M-1
            call print_all_r8(n_alpha,alpha_est__(n_alpha*nM)
     $           ,' alpha_est__: ')
            call print_all_r8(n_alpha,alpha_upd__(n_alpha*nM)
     $           ,' alpha_upd__: ')
         enddo ! do nM=0,n_M-1
      end if !if (verbose.gt.2) then

      do ns=0,n_S-1
      if (n_MS_(ns).lt.n_MS_max) then
         call alpha_SM_sort_0(n_MS_max,n_MS_(ns),alpha_MS__(n_alpha
     $        *n_MS_max*ns))
      end if !if (n_MS_(ns).lt.n_MS_max) then
      enddo !do ns=0,n_S-1

      allocate(I_S_(0:n_S-1))
      allocate(B_M_(0:n_M-1))
      allocate(n_MS_tmp_(0:n_S-1))
      call cp1_i4(n_S,n_MS_,n_MS_tmp_)
      n_MS_sum = sum_i4_f(n_s,n_MS_)

      call adi_randperm(rseed,n_S,I_S_)
      if (verbose.gt.1) then
         call print_all_i4(n_S,I_S_,' I_S_: ')
         call print_all_i4(n_S,n_MS_,' n_MS_: ')
      end if !if (verbose.gt.1) then

      do nm=0,n_M-1
         B_M_(nm) = .false.
      enddo
      nm_sum = 0
      ns_tmp = 0
      do while (nm_sum.le.n_M-1 .and. nm_sum.le.n_MS_sum-1)
         ns = I_S_(ns_tmp)
         if (verbose.gt.2) then
            call print_all_i4(n_S,n_MS_tmp_,' n_MS_tmp_: ')
            write(6,'(4(A,I0))') 'nm_sum ',nm_sum,'; ns_tmp ', ns_tmp
     $           ,'; ns ',ns , '; n_MS_tmp_: ' , n_MS_tmp_(ns)
         end if !if (verbose.gt.2) then
         nm_srt = n_MS_tmp_(ns)-1
         continue_flag = .true.
         do while (continue_flag)
            continue_flag = .false.
            na = n_alpha*(nm_srt + n_MS_max*ns)
            nm = alpha_MS__(nalpha_M_index + na)
            bm = B_M_(nm)
            if (bm) then
               if (verbose.gt.2) then
                  write(6,'(A,I0,A,I0,A)') ' nm_srt ',nm_srt,'; nm ',nm
     $                 ,'; skipping'
               end if
               n_MS_tmp_(ns) = n_MS_tmp_(ns) - 1
               nm_srt = nm_srt - 1
               if (nm_srt.ge.0) then
                  continue_flag = .true.
               end if !if (nm_srt.ge.0) then
               if (nm_srt.lt.0) then
                  continue_flag = .false.
               end if !if (nm_srt.lt.0) then
            else ! if unused
               if (verbose.gt.2) then
                  write(6,'(A,I0,A,I0,A)') ' nm_srt ',nm_srt,'; nm ',nm
     $                 ,'; updating'
               end if
               bm = .true.
               B_M_(nm) = .true.
               n_MS_tmp_(ns) = n_MS_tmp_(ns) - 1
               continue_flag = .false.
               nm_sum = nm_sum + 1
c$$$           %%%%%%%%%%%%%%%%
               nb = n_alpha*nm
               delta_x_upd = alpha_MS__(nalpha_delta_x + na)
               delta_y_upd = alpha_MS__(nalpha_delta_y + na)
               gamma_z_upd = alpha_MS__(nalpha_gamma_z + na)
               delta_x_est = alpha_est__(nalpha_delta_x + nb)
               delta_y_est = alpha_est__(nalpha_delta_y + nb)
               gamma_z_est = alpha_est__(nalpha_gamma_z + nb)
               call get_interchange_delta_RTRT_vs_RTTR(flag_RTRT_vs_RTTR
     $              ,delta_x_est,delta_y_est,gamma_z_est,delta_x_upd
     $              ,delta_y_upd ,gamma_z_upd)
               call cp1_r8(n_alpha,alpha_MS__(na),alpha_upd__(nb))
               alpha_upd__(nalpha_delta_x + nb) = delta_x_upd
               alpha_upd__(nalpha_delta_y + nb) = delta_y_upd
               alpha_upd__(nalpha_gamma_z + nb) = gamma_z_upd
c$$$           %%%%%%%%%%%%%%%%
            end if ! check bm
         end do ! continue_flag
         ns_tmp = ns_tmp + 1
         if (ns_tmp.ge.n_S) then
            ns_tmp = 0
         end if
      end do ! while (nm_sum.le.n_M-1 .and. nm_sum.le.n_MS_sum-1)
      ns = I_S_(ns_tmp)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0,A,I0,A)') 'nm_sum ',nm_sum,'; ns_tmp '
     $        ,ns_tmp,'; ns ',ns,'; finished'
      end if !if (verbose.gt.1) then

      deallocate(I_S_)
      deallocate(B_M_)
      deallocate(n_MS_tmp_)

      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_alpha_update_MS_3]'
      end if
      end

