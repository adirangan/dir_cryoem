!> Doxygen comment: ;\n
!> Maintains a stack of n_SM .le. n_SM_max image parameters called alpha_SM_. ;\n
!> Each time a new set of image parameters (alpha_0in) is about to be added, ;\n
!> we check the 'best' element of alpha_SM_ (ranked by label C_Z_opt) ;\n
!> to decide whether or not to actually add the new set of image parameters. ;\n
!> In this way we can ensure that the stack never references more than n_SM_max images. ;\n
      subroutine alpha_SM_update_1(n_SM_max,n_SM,alpha_SM_,alpha_0in)
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
      include 'excerpt_define_nalpha.f'
      real *8 alpha_0in(0:n_alpha-1)
      character(len=1024) format_string

      if (verbose.gt.1) then
         write(6,'(A)') '[entering alpha_SM_update_1]'
      end if !if (verbose.gt.1) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' n_SM_max: ' , n_SM_max
         write(6,'(A,I0)') ' n_SM: ' , n_SM
         write(6,'(A,I0)') ' n_alpha: ' , n_alpha
         call print_sub_r8(n_alpha*n_SM_max,alpha_SM_
     $        ,' alpha_SM_: ')
         call print_all_r8(n_alpha,alpha_0in
     $        ,' alpha_0in_: ')
      end if !if (verbose.gt.1) then

      if (n_SM.lt.n_SM_max) then
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' adding at n_SM: ' , n_SM
         end if !if (verbose.gt.1) then
         call cp1_r8(n_alpha,alpha_0in,alpha_SM_(n_alpha*n_SM))
         n_SM = n_SM+1
         if (n_SM.eq.n_SM_max) then
         if (verbose.gt.1) then
            write(6,'(A,I0,A)') ' n_SM: ' , n_SM , '; sorting. '
         end if !if (verbose.gt.1) then
         call alpha_SM_sort_0(n_SM_max,n_SM,alpha_SM_)
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
            if (alpha_0in(nalpha_C_Z_opt).lt.alpha_SM_(nalpha_C_Z_opt +
     $           nSM*n_alpha)) then
               flag_continue = .false.
            else
               nSM = nSM+1
               flag_continue = .true.
               if (nSM.ge.n_SM_max) then
                  flag_continue = .false.
               end if !if (nSM.ge.n_SM_max-1) then
            end if !if (alpha_0in(nalpha_C_Z_opt).lt.alpha_SM_(nalpha_C_Z_opt + nSM*n_alpha)) then
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
            call cp1_r8(n_alpha,alpha_0in,alpha_SM_(n_alpha*(nSM-1)))
         end if !if (nSM.gt.0) then         
         if (verbose.gt.2) then
            call alpha_SM_write_0(n_SM_max,n_SM ,alpha_SM_,13
     $           ,' alpha_SM_:  ')
         end if !if (verbose.gt.2) then
         goto 10
      end if !if (n_SM.eq.n_SM_max) then

 10   continue

      if (verbose.gt.1) then
         write(6,'(A)') '[finished alpha_SM_update_1]'
      end if !if (verbose.gt.1) then
      end
