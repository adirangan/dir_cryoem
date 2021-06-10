!> Doxygen comment: ;\n
!> Generate a random permutation of n_A elements. ;\n
!> Assumes that the integer *4 output array I_ is already preallocated. ;\n
      subroutine randperm(n_A,I_)
      implicit none
      integer n_A,na
      integer *4 I_(0:n_A-1)
      real *8, allocatable :: f_(:)
      allocate(f_(0:n_A-1))
      do na=0,n_A-1
         f_(na) = rand()
         I_(na) = na
      enddo
      call quicksort_r8(0,n_A-1,f_,1,I_,1,quicksort_r8)
      deallocate(f_)
      end
