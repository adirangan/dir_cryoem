      subroutine MDA_write_c16(n_d,d_,A_,fname)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_d,nd,n_A,na
      integer d_(0:n_d-1)
      complex *16 A_(0:0)
      character(len=1024) fname
      n_A=1
      do nd=0,n_d-1
         n_A = n_A*d_(nd)
      enddo
      open(7,file=fname,status='replace',form='unformatted')
      if (verbose.gt.0) then
         write(6,*) 'Writing: ',n_d
      end if
      write(7) n_d
      if (verbose.gt.0) then
      write(6,*) 'Writing: ',(d_(nd),nd=0,n_d-1)
      end if
      write(7) (d_(nd),nd=0,n_d-1)
      if (verbose.gt.0) then
      write(6,*) 'Writing: ',(A_(na),na=0,n_A-1)
      end if
      write(7) (A_(na),na=0,n_A-1)
      close(7,status='keep')
      end
