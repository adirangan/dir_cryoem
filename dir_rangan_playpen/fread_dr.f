      program fread_dr
c     gfortran -o fread_dr.out fread_dr.f ; ./fread_dr.out ; 
      implicit none
      character(len=64) fname,format_string
      integer *4 n_A,na
      real *8, allocatable :: A_(:)
      write(fname,'(A)') 'data.test'
      open(7,file=fname,status='old',form='formatted',action='read')
      read(7,*) n_A
      write(6,'(A,I0)') 'n_A: ',n_A
      allocate(A_(0:n_A-1))
      read(7,*) (A_(na),na=0,n_A-1)
      write(format_string,'(A,I0,A)') '(A,',n_A,'(F6.3,1X))'
      write(6,format_string) 'A_:',(A_(na),na=0,n_A-1)
      close(7,status='keep')
      deallocate(A_)
      stop
      end
