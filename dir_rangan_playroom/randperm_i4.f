      subroutine randperm_i4(n_A,I_)
      implicit none
      integer n_A,na
      integer *4 I_(0:n_A-1)
      integer *4, allocatable :: f_(:)
      allocate(f_(0:n_A-1))
      do na=0,n_A-1
         f_(na) = nint(1000*n_A*rand())
         I_(na) = na
      enddo
      call quicksort_i4(0,n_A-1,f_,1,I_,1,quicksort_i4)
      deallocate(f_)
      end
