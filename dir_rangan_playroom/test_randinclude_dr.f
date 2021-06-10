c     gfortran -w -o test_randinclude_dr.out test_randinclude_dr.f ;./test_randinclude_dr.out
      program test_randinclude_dr
      implicit none
      integer *4 n_LT,n_S,S_used_total,n_add
      integer *4, allocatable :: LT_(:)
      logical, allocatable :: S_used_(:)
      character(len=64) format_string_S
      character(len=64) format_string_LT
      integer *4 nLT,nS
      n_S = 32
      write(format_string_S,'(A,I0,A)') '(',n_S,'(L1,1X))'
      n_LT = 0
      allocate(LT_(0:n_S-1))
      allocate(S_used_(0:n_S-1))
      do nS=0,n_S-1
         S_used_(ns) = .false.
      enddo
      S_used_total=0
      S_used_(3) = .true.
      S_used_(5) = .true.
      S_used_(9) = .true.
      S_used_total=3
      n_add = 3
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      n_add = 5
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      n_add = 4
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      n_add = 10
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      n_add = 15
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      n_add = 15
      write(6,'(A,I0)') 'calling randinclude with n_add: ',n_add
      call randinclude(n_LT,LT_,n_S,S_used_,S_used_total,n_add)
      write(format_string_LT,'(A,I0,A)') '(',n_LT,'(I0,1X))'
      write(6,format_string_LT) (LT_(nLT),nLT=0,n_LT-1)
      write(6,format_string_S) (S_used_(nS),nS=0,n_S-1)
      write(6,'(A,I0)') 'finished, now n_LT: ',n_LT
      write(6,'(A,I0)') 'finished, now S_used_total: ',S_used_total
      stop
      end

      include 'quicksort_r8.f'
      include 'quicksort_i4.f'
      include 'randperm.f'
      include 'randinclude.f'
