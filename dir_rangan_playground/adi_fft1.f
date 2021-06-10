!> Doxygen comment: ;\n
!> Simple brute force fourier transform for 1 dimension. ;\n
!> Only use for testing. ;\n
      subroutine adi_fft1(sgn,n,a_,a_start,a_stride,c_,c_start,c_stride)
      integer sgn
      integer n,j,k,a_start,a_stride,c_start,c_stride
      complex *16 a_(0:(a_start + n*a_stride)-1)
      complex *16 c_(0:(c_start + n*c_stride)-1)
      complex *16, allocatable :: b_(:)
      complex *16 wk,wkj
      real *8 pi,theta
      integer kb,ja
      pi = 4.0d0*datan(1.0d0)
      allocate(b_(0:n-1))
      do k=0,n-1
         theta = 2.0d0*pi*k/n
         wk = dcmplx(+dcos(theta),sgn*dsin(theta))
         wkj = ( 1.0 , 0.0 )
         b_(k)=0
         ja = a_start
         do j=0,n-1
            b_(k) = b_(k) + wkj*a_(ja)
            wkj = wkj*wk
            ja = ja + a_stride
         enddo
         if (sgn.eq.+1) then
            b_(k) = b_(k) / n
         end if
      enddo
      do k=0,n-1
         c_(c_start + k*c_stride) = b_(k)
      enddo
      deallocate(b_)
      end
