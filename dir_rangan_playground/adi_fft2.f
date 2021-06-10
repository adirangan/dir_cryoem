!> Doxygen comment: ;\n
!> Simple brute force fourier transform for 2 dimensions. ;\n
!> Only use for testing. ;\n
      subroutine adi_fft2(sgn,m,n,a_,c_)
      integer sgn,m,n,j
      complex *16 a_(0:m*n-1),c_(0:m*n-1)
      complex *16, allocatable :: b_(:)
      allocate(b_(0:m*n-1))
      do k=0,n-1
         call adi_fft1(sgn,m,a_,k*m,1,b_,k*m,1)
      enddo
      do j=0,m-1
         call adi_fft1(sgn,n,b_,j,m,c_,j,m)
      enddo
      deallocate(b_)
      end
