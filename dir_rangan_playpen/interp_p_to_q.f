      subroutine interp_p_to_q(n_r,n_w_,n_A,S_p_,S_q_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      complex *16 S_p_(0:n_A-1),S_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      integer ic,nr
      integer nt,n_t
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write(6,*) '[entering interp_p_to_q]'
         timing_tic = second()
      end if
      if (verbose.gt.0) then
         n_t = 10000
         write(6,*) 'running n_t times with n_t: ',n_t
      else
         n_t = 1
      end if
      allocate(Z_q_(0:n_A-1));
      do nt=0,n_t-1
      ic=0
      do nr=0,n_r-1
         call adi_fft1(-1,n_w_(nr),S_p_,ic,1,Z_q_,ic,1)
         call af1_c16(n_w_(nr),1.0d0*cmplx(1.0d0/dsqrt(1.0d0 *n_w_(nr)))
     $        ,1.0d0*cmplx(0.0d0),Z_q_(ic),Z_q_(ic))
         ic = ic + n_w_(nr)
      enddo
      enddo
      call cp1_c16(n_A,Z_q_,S_q_);
      deallocate(Z_q_)
      if (verbose.gt.0) then
         timing_toc = second()
         write(6,'(A,F6.3)') ' [finished interp_p_to_q]: total time: '
     $        ,timing_toc-timing_tic
      end if
      end
