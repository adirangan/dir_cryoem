c     gfortran -w -o test_recenter2_c16_dr.out test_recenter2_c16_dr.f ; ./test_recenter2_c16_dr.out
      program test_recenter2_c16_dr
      implicit none
      integer verbose
      data verbose / 2 /
      integer *4 n_x,n_y,nx,ny,na
      complex *16, allocatable :: a_(:)
      character(len=64) format_string

      n_x = 7
      n_y = 9
      allocate(a_(0:n_x*n_y-1))
      do ny=0,n_y-1
         do nx=0,n_x-1
            a_(nx+ny*n_x) = cmplx( 1.0d0*nx , 1.0d0*ny )
         enddo
      enddo
      
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(' , n_x , '(2F8.4,4X))'
         write(6,'(A)') 'original a_: '
         write(6,format_string) (a_(na),na=0,n_x*n_y-1)
      end if !verbose

      call recenter2_c16(n_x,n_y,a_,a_)
      
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(' , n_x , '(2F8.4,4X))'
         write(6,'(A)') 'recenter2ed a_: '
         write(6,format_string) (a_(na),na=0,n_x*n_y-1)
      end if !verbose

      call recenter2_c16(n_x,n_y,a_,a_)
      
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(' , n_x , '(2F8.4,4X))'
         write(6,'(A)') 'recenter2ed again a_: '
         write(6,format_string) (a_(na),na=0,n_x*n_y-1)
      end if !verbose

 10   stop
      end

      include 'recenter2_c16.f'
