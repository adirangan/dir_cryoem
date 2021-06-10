c     gfortran -w -o test_unique_i4_dr.out test_unique_i4_dr.f ;./test_unique_i4_dr.out
      program test_unique_i4_dr
      implicit none
      integer verbose
      data verbose / 1 /
c$$$      indices
      integer *4 n_L,nL,n_I,nI
c$$$      arrays
      integer *4, allocatable :: L_(:)
      integer *4, allocatable :: I_(:)
      character(len=64) format_string

      n_L = 32
      allocate(L_(0:n_L-1))
      do nL=0,n_L-1
         L_(nL) = nint(256*rand())
      enddo !do nL=0,n_L-1
      write(format_string,'(A,I0,A)') '(',n_L,'(I0,1X))'
      write(6,format_string) (L_(nL),nL=0,n_L-1)
      allocate(I_(0:n_L-1));
      n_I = 0;
      write(6,'(A)') 'calling unique_i4: '
      call unique_i4(n_L,L_,n_I,I_)
      write(6,'(A,I0)') 'n_I: ',n_I
      write(format_string,'(A,I0,A)') '(',n_I,'(I0,1X))'
      write(6,format_string) (I_(nI),nI=0,n_I-1)
      write(6,'(A)') 'calling unique_i4 in place: '
      call unique_i4(n_L,L_,n_L,L_)
      write(6,'(A,I0)') 'n_L: ',n_L
      write(format_string,'(A,I0,A)') '(',n_L,'(I0,1X))'
      write(6,format_string) (L_(nL),nL=0,n_L-1)

 10   stop
      end

      include 'cp1_i4.f'
      include 'quicksort_i4.f'
      include 'unique_i4.f'


