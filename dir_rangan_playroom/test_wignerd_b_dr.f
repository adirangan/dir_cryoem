c     gfortran -w -o test_wignerd_b_dr.out test_wignerd_b_dr.f ; ./test_wignerd_b_dr.out
      program test_wignerd_b_dr
      implicit none
      integer verbose
      data verbose / 2 /
      integer *4 n_l,nl,n_W_sum,n_W,nw
      real *8 pi,beta
      real *8, allocatable :: W_(:)
      integer *4, allocatable :: n_W_(:)
      character(len=64) format_string

      pi = 4.0*atan(1.0)
      n_l = 5
      beta = -2*pi/5
      if (verbose.gt.0) then
         write(6,'(A,I0,A,F10.6)') 'n_l = ',n_l, ' beta = ',beta
      end if
      
      n_W_sum = 0
      do nl=0,n_l
         n_W = (1+2*nl)*(1+2*nl)
         n_W_sum = n_W_sum + n_W
      enddo !do nl=0,n_l
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_W_sum = ',n_W_sum
      end if
      allocate(W_(0:n_W_sum-1))
      allocate(n_W_(0:n_l))
      
      call wignerd_b(n_l,beta,W_,n_W_)

      if (verbose.gt.0) then
      do nl=0,n_l
         n_W = (1+2*nl)*(1+2*nl)
         write(6,'(A,I0,A,I0,A)') 'nl = ' , nl , ' n_W = ' , n_W ,
     $        ' W_^{t}: '
         write(format_string,'(A,I0,A)') '(' , 1+2*nl , 'F10.6)'
         write(6,format_string) (W_(nw),nw=n_W_(nl),n_W_(nl) + n_W-1)
      enddo !for nl=0,n_l
      end if !verbose

 10   stop
      end

      include 'wignerd_b.f'
