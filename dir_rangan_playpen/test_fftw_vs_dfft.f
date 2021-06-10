      program test_fftw_vs_dfft
c$$$      gfortran -w -O2 -o test_fftw_vs_dfft.out test_fftw_vs_dfft.f -I/usr/include/ -L/usr/lib -lfftw3 ; ./test_fftw_vs_dfft.out ;
      implicit none
      include '/usr/include/fftw3.f'
      integer n_iter,ni,n_c_base,n_c,nc
      parameter (n_iter = 10)
      integer n_fftw
      parameter (n_fftw=2048)
      double complex fftw_in1_(0:n_fftw*n_iter-1)
      double complex fftw_in2_(0:n_fftw*n_iter-1)
      double complex fftw_out_(0:n_fftw*n_iter-1)
      complex *16 Wdata_(0:n_fftw-1)
      complex *16 Wsave_(0:4*n_fftw+16-1)
      real *8 pi,theta
      integer*8 fftw_plan_back_(0:n_iter-1)
      integer*8 fftw_plan_frwd_(0:n_iter-1)
c$$$      parameters for timing
      integer n_t,nt
      real *8 timing_tic,timing_toc
      pi = 4.0*atan(1.0)

      n_c_base = 2
      do ni=0,n_iter-1
         n_c_base = n_c_base*2
         n_c = n_c_base - 1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(ni),n_c,fftw_in1_(ni
     $        *n_fftw),fftw_out_(ni*n_fftw),FFTW_FORWARD,FFTW_ESTIMATE)
         call dfftw_plan_dft_1d_(fftw_plan_back_(ni),n_c,fftw_out_(ni
     $        *n_fftw),fftw_in2_(ni*n_fftw),FFTW_BACKWARD,FFTW_ESTIMATE)
      enddo

      n_c_base = 2
      n_t = 1024*1024*4
      do ni=0,n_iter-1
         n_c_base = n_c_base*2
         n_c = n_c_base - 1
         n_t = n_t/2
         write(6,'(A,I0,A,I0)') 'n_c: ',n_c,'; n_t: ',n_t
         call dcffti(n_c,Wsave_)
         do nc = 0,n_c-1
            theta = 1*nc*2*pi/n_c
            Wdata_(nc) = dcmplx(cos(theta),sin(theta))
            fftw_in1_(nc+ni*n_fftw) = dcmplx(cos(theta),sin(theta))
            fftw_out_(nc+ni*n_fftw) = dcmplx(cos(theta),sin(theta))
            fftw_in2_(nc+ni*n_fftw) = dcmplx(cos(theta),sin(theta))
         end do
         timing_tic = second()
         do nt=0,n_t-1
            call dcfftf(n_c,Wdata_,Wsave_)
         enddo
         timing_toc = second()
         write(6,'(A,F8.5)') ' dcfftf         total time: ',timing_toc
     $        -timing_tic
         timing_tic = second()
         do nt=0,n_t-1
            call dcfftb(n_c,Wdata_,Wsave_)
         enddo
         timing_toc = second()
         write(6,'(A,F8.5)') ' dcfftb         total time: ',timing_toc
     $        -timing_tic
         timing_tic = second()
         do nt=0,n_t-1
            call cp1_c16(n_c,Wdata_,fftw_in1_(ni*n_iter))
            call dfftw_execute_(fftw_plan_frwd_(ni))
            call cp1_c16(n_c,fftw_out_(ni*n_fftw),Wdata_)
         enddo
         timing_toc = second()
         write(6,'(A,F8.5)') ' fftw_plan_frwd total time: ',timing_toc
     $        -timing_tic
         timing_tic = second()
         do nt=0,n_t-1
            call cp1_c16(n_c,Wdata_,fftw_out_(ni*n_fftw))
            call dfftw_execute_(fftw_plan_back_(ni))
            call cp1_c16(n_c,fftw_in2_(ni*n_fftw),Wdata_)
         enddo
         timing_toc = second()
         write(6,'(A,F8.5)') ' fftw_plan_back total time: ',timing_toc
     $        -timing_tic
       enddo

      do ni=0,n_iter-1
         call dfftw_destroy_plan_(fftw_plan_frwd_(ni))
         call dfftw_destroy_plan_(fftw_plan_back_(ni))
      enddo

      stop
      end

      include 'dfftpack.f'
      include 'cp1_c16.f'
