c     gfortran -o helloworld.out helloworld.f ; ./helloworld.out ; cat helloworld.txt ;
      program helloworld
      real *8 pi
      integer *4 i,j
      real *8, allocatable :: x_(:,:)
      i = 3
      j = 7
      pi = 4*atan(1.0)
      allocate(x_(0:i-1,0:j-1))
      do i=0,2
         do j=0,6
            x_(i,j) = i + j*3
         enddo
      enddo
      open(7,file='helloworld.txt',status='unknown',form='formatted')
      write(7,100) 'hello','world'
 100  format(A,1X,A)
      write(7,200) ((x_(i,j),j=0,6),i=0,2)
 200  format(7F4.0)
      close(7,status='keep')
      stop
      end
