!> Doxygen comment: ;\n
!> Extracts sub-stack of n_SM image parameters from the full stack alpha_SM_, ;\n
!> Then sorts the stack by them by label C_Z_opt and puts them back in. ;\n
      subroutine alpha_SM_sort_0(n_SM_max,n_SM,alpha_SM_)
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
         write(6,'(A,I0,A)') ' n_SM: ' , n_SM , '; sorting. '
      end if                    !if (verbose.gt.1) then

      if (verbose.gt.2) then
         call alpha_SM_write_0(n_SM_max,n_SM,alpha_SM_,16
     $        ,' alpha_SM_ pre: ')
      end if                    !if (verbose.gt.2) then

      allocate(alpha_tmp_(0:n_alpha*n_SM_max-1))
      allocate(C_Z_opt_(0:n_SM_max-1))
      allocate(I_permute_(0:n_SM_max-1))

      do nSM=0,n_SM_max-1
         C_Z_opt_(nSM) = alpha_SM_(nalpha_C_Z_opt + nSM*n_alpha)
         I_permute_(nSM) = nSM
      enddo                     !do nSM=0,n_SM_max-1

      do nSM=0,n_SM-1
         I_permute_(nSM) = nSM
      enddo                     !do nSM=0,n_SM-1

      call quicksort_r8(0,n_SM-1,C_Z_opt_,1,I_permute_,1
     $     ,quicksort_r8)

      do nSM=0,n_SM-1
         call cp1_r8(n_alpha,alpha_SM_(n_alpha*I_permute_(nSM))
     $        ,alpha_tmp_(n_alpha*nSM))
      enddo                     !do nSM=0,n_SM-1

      call cp1_r8(n_alpha*n_SM,alpha_tmp_,alpha_SM_)

      deallocate(alpha_tmp_)
      deallocate(C_Z_opt_)
      deallocate(I_permute_)

      if (verbose.gt.2) then
         call alpha_SM_write_0(n_SM_max,n_SM,alpha_SM_,16
     $        ,' alpha_SM_ pos: ')
      end if                    !if (verbose.gt.2) then

      end
