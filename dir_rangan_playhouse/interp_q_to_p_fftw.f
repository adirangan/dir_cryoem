      subroutine interp_q_to_p_fftw(n_r,fftw_plan_back_,n_w_,n_A
     $     ,fftw_in1_,fftw_out_,S_q_,S_p_)
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      integer *8 fftw_plan_back_(0:n_r-1)
      complex *16 fftw_in1_(0:n_A-1),fftw_out_(0:n_A-1)
      complex *16 S_q_(0:n_A-1),S_p_(0:n_A-1)
      integer n_w_max,ic,nr
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write(6,*) '[entering interp_q_to_p_fftw]'
         timing_tic = omp_get_wtime()
      end if
      if (verbose.gt.0) then
         write(6,*) '% n_r: ',n_r
         write(6,*) '% n_w_: ',(n_w_(nr),nr=0,n_r-1)
         write(6,*) '% n_A: ',n_A
      end if
      n_w_max = n_w_(n_r-1)
      ic=0
      do nr=0,n_r-1
         if (n_w_(nr).gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if
            call cp1_c16(n_w_(nr),S_q_(ic),fftw_out_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; executing fftw_plan: '
            end if
            call dfftw_execute_(fftw_plan_back_(nr))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling affine_c16: '
            end if
            call af1_c16(n_w_(nr),1.0d0*cmplx(1.0d0*dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*cmplx(0.0d0),fftw_in1_(ic),S_p_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % setting ic ',ic,' to ic ',ic + n_w_(nr)
            end if
            ic = ic + n_w_(nr)
         end if
      enddo
      if (verbose.gt.0) then
         timing_toc = omp_get_wtime()
         write(6,'(A,F6.3)')
     $        ' [finished interp_q_to_p_fftw]: total time: '
     $        ,timing_toc-timing_tic
      end if
      end
