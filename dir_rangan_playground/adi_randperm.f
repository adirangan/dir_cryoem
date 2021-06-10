!> Doxygen comment: ;\n
!> generate pseudo-random permutation of size n_A ;\n
      subroutine adi_randperm(rseed,n_A,I_)
      implicit none
      integer *4 rseed
      real *8 adi_rand_f
      integer n_A,na
      integer *4 I_(0:n_A-1)
      real *8, allocatable :: f_(:)
      allocate(f_(0:n_A-1))
      do na=0,n_A-1
         f_(na) = adi_rand_f(rseed)
         I_(na) = na
      enddo
      call quicksort_r8(0,n_A-1,f_,1,I_,1,quicksort_r8)
      deallocate(f_)
      end
