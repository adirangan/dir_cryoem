      subroutine MDA_write_i4(n_d,d_,A_,fname)
      implicit none
      integer n_d,nd,n_A,na
      integer d_(0:n_d-1)
      integer *4 A_(0:0)
      character(len=1024) fname
      n_A=1
      do nd=0,n_d-1
         n_A = n_A*d_(nd)
      enddo
      open(7,file=fname,status='replace',form='unformatted')
      write(7) n_d
      write(7) (d_(nd),nd=0,n_d-1)
      write(7) (A_(na),na=0,n_A-1)
      close(7,status='keep')
      end
