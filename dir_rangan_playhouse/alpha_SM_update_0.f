      subroutine alpha_SM_update_0(n_SM_max,n_SM,alpha_SM_,alpha__in)
      implicit none
      integer verbose
      data verbose / 0 /
c$$$      SM storage
      integer *4 n_SM_max ! total (maximum) number of templates to store per image. ;
      integer *4 n_SM ! the actual number of templates stored for this particular image. ;
      real *8 alpha_SM_(0:0) ! array of size n_alpha*n_SM_max storing the image-parameters for each stored template-image pair for this particular image. ;
      integer nSM,nSM_sub
      logical flag_continue
      real *8, allocatable :: alpha_tmp_(:)
      real *8, allocatable :: C_Z_opt_(:)
      integer *4, allocatable :: I_permute_(:)
      external quicksort_c16
      include 'nalpha_define.f'
      real *8 alpha__in(0:n_alpha-1)
      character(len=1024) format_string

      if (verbose.gt.1) then
         write(6,'(A)') '[entering alpha_SM_update_0]'
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_SM_max: ' , n_SM_max
         write(6,'(A,I0)') ' n_SM: ' , n_SM
         write(6,'(A,I0)') ' n_alpha: ' , n_alpha
         call write_sub_r8(n_alpha*n_SM_max,alpha_SM_,12
     $        ,' alpha_SM_: ')
         call write_all_r8(n_alpha,alpha__in,13
     $        ,' alpha__in_: ')
      end if !if (verbose.gt.1) then

      if (n_SM.lt.n_SM_max) then
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' adding at n_SM: ' , n_SM
         end if !if (verbose.gt.1) then
         call cp1_r8(n_alpha,alpha__in,alpha_SM_(n_alpha*n_SM))
         n_SM = n_SM+1
         if (n_SM.eq.n_SM_max) then
         if (verbose.gt.1) then
            write(6,'(A,I0,A)') ' n_SM: ' , n_SM , '; sorting. '
         end if !if (verbose.gt.1) then
         if (verbose.gt.2) then
            write(6,'(A)') ' alpha_SM_ pre: '
            write(format_string,'(A,I0,A)') '(' , n_alpha , '(F8.4,1X))'
            do nSM=0,n_SM_max-1
               write(6,format_string) (alpha_SM_(nSM_sub),nSM_sub
     $              =n_alpha*nSM,n_alpha*nSM+n_alpha-1)
            enddo !do nSM=0,n_SM_max-1
         end if !if (verbose.gt.2) then
         allocate(alpha_tmp_(0:n_alpha*n_SM_max-1))
         allocate(C_Z_opt_(0:n_SM_max-1))
         allocate(I_permute_(0:n_SM_max-1))
         do nSM=0,n_SM_max-1
            C_Z_opt_(nSM) = alpha_SM_(nalpha_C_Z_opt + nSM*n_alpha)
            I_permute_(nSM) = nSM
         enddo !do nSM=0,n_SM_max-1
         call quicksort_r8(0,n_SM_max-1,C_Z_opt_,1,I_permute_,1
     $        ,quicksort_r8)
         do nSM=0,n_SM_max-1
            call cp1_r8(n_alpha,alpha_SM_(n_alpha*I_permute_(nSM))
     $           ,alpha_tmp_(n_alpha*nSM))
         enddo !do nSM=0,n_SM_max-1
         call cp1_r8(n_alpha*n_SM_max,alpha_tmp_,alpha_SM_)
         deallocate(alpha_tmp_)
         deallocate(C_Z_opt_)
         deallocate(I_permute_)
         if (verbose.gt.2) then
            write(6,'(A)') ' alpha_SM_ pos: '
            write(format_string,'(A,I0,A)') '(' , n_alpha , '(F8.4,1X))'
            do nSM=0,n_SM_max-1
               write(6,format_string) (alpha_SM_(nSM_sub),nSM_sub
     $              =n_alpha*nSM,n_alpha*nSM+n_alpha-1)
            enddo !do nSM=0,n_SM_max-1
         end if !if (verbose.gt.2) then
         end if !if (n_SM.eq.n_SM_max-1) then
         goto 10
      end if !if (n_SM.lt.n_SM_max) then
      if (n_SM.eq.n_SM_max) then
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' alpha_SM full; n_SM: ' , n_SM
         end if !if (verbose.gt.1) then
         flag_continue=.true.
         nSM=0
         do while (flag_continue)
            if (alpha__in(nalpha_C_Z_opt).lt.alpha_SM_(nalpha_C_Z_opt +
     $           nSM*n_alpha)) then
               flag_continue = .false.
            else
               nSM = nSM+1
               flag_continue = .true.
               if (nSM.ge.n_SM_max) then
                  flag_continue = .false.
               end if !if (nSM.ge.n_SM_max-1) then
            end if !if (alpha__in(nalpha_C_Z_opt).lt.alpha_SM_(nalpha_C_Z_opt + nSM*n_alpha)) then
         enddo !do while (flag_continue)
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' alpha_SM full; insert at nSM-1 : ' , nSM
     $           -1
         end if !if (verbose.gt.1) then
         if (nSM.gt.0) then
            do nSM_sub=0,nSM-2
               if (verbose.gt.1) then
                  write(6,'(A,I0,A,I0)') ' replacing nSM ' , nSM_sub ,
     $                 ' with ' , nSM_sub+1
               end if !if (verbose.gt.1) then
               call cp1_r8(n_alpha,alpha_SM_(n_alpha*(nSM_sub+1))
     $              ,alpha_SM_(n_alpha*(nSM_sub)))
            enddo !do nSM_sub=0,nSM-2
            if (verbose.gt.1) then
               write(6,'(A,I0,A)') ' replacing nSM ' , nSM-1 ,
     $              ' with new data '
            end if !if (verbose.gt.1) then
            call cp1_r8(n_alpha,alpha__in,alpha_SM_(n_alpha*(nSM-1)))
         end if !if (nSM.gt.0) then         
         if (verbose.gt.2) then
            write(6,'(A)') ' alpha_SM_: '
            write(format_string,'(A,I0,A)') '(' , n_alpha , '(F8.4,1X))'
            do nSM=0,n_SM_max-1
               write(6,format_string) (alpha_SM_(nSM_sub),nSM_sub
     $              =n_alpha*nSM,n_alpha*nSM+n_alpha-1)
            enddo !do nSM=0,n_SM_max-1
         end if !if (verbose.gt.2) then
         goto 10
      end if !if (n_SM.eq.n_SM_max) then

 10   continue

      if (verbose.gt.1) then
         write(6,'(A)') '[finished alpha_SM_update_0]'
      end if !if (verbose.gt.1) then
      end
