c     gfortran -w -o test_randperm_dr.out test_randperm_dr.f ;./test_randperm_dr.out
      program test_randperm_dr
      implicit none
      integer n_A,na
      integer, allocatable :: I_(:)
      character(len=64) format_string
      n_A = 32
      write(format_string,'(A,I0,A)') '(',n_A,'(I0,1X))'
      allocate(I_(0:n_A-1))
      call randperm(n_A,I_)
      write(6,format_string) (I_(na),na=0,n_A-1)
      call randperm(n_A,I_)
      write(6,format_string) (I_(na),na=0,n_A-1)
      call randperm(n_A,I_)
      write(6,format_string) (I_(na),na=0,n_A-1)
      call randperm(n_A,I_)
      write(6,format_string) (I_(na),na=0,n_A-1)
      call randperm(n_A,I_)
      write(6,format_string) (I_(na),na=0,n_A-1)
      stop
      end

      include 'quicksort_r8.f'
      include 'randperm.f'
