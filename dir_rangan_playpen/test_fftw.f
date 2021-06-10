      program test_fftw
c$$$      gfortran -w -o test_fftw.out test_fftw.f -I/usr/include/ -L/usr/lib -lfftw3 ; ./test_fftw.out ;
      implicit none
      include '/usr/include/fftw3.f'
      integer n_fftw
      parameter (n_fftw=8)
      integer na
      double complex fftw_in1_(0:n_fftw-1)
      double complex fftw_in2_(0:n_fftw-1)
      double complex fftw_out_(0:n_fftw-1)
      real *8 pi,theta
      integer*8 plan_backward
      integer*8 plan_forward      
      pi = 4.0*atan(1.0)
      do na = 0,n_fftw-1
         theta = 1*na*2*pi/n_fftw
         fftw_in1_(na) = dcmplx(cos(theta),sin(theta))
      end do
      write (6,'(A)') 'fftw_in1_:'
      write (6,'(2F8.3)') (fftw_in1_(na),na=0,n_fftw-1)
      call dfftw_plan_dft_1d_(plan_forward,n_fftw,fftw_in1_,fftw_out_,
     $     FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_(plan_forward)
      write (6,'(A)') 'fftw_out_:'
      write (6,'(2F8.3)') (fftw_out_(na),na=0,n_fftw-1)
      call dfftw_plan_dft_1d_(plan_backward,n_fftw,fftw_out_,fftw_in2_,
     $     FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute_(plan_backward)
      write (6,'(A,I0)') 'fftw_in2_ / ',n_fftw
      write (6,'(2F8.3)') (fftw_in2_(na)/n_fftw,na=0,n_fftw-1)
      call dfftw_destroy_plan_(plan_forward)
      call dfftw_destroy_plan_(plan_backward)
      stop
      end

      include 'dfftpack.f'
