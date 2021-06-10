      subroutine interp_p_to_q_dcfft(n_r,n_w_,n_w_W_,n_A_W,Wsave_,n_A
     $     ,S_p_,S_q_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_w_W_(0:n_r-1),n_A,n_A_W
      complex *16 Wsave_(0:n_A_W-1)
      complex *16 S_p_(0:n_A-1),S_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      integer ic,nr,ic_W,n_w_max
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write(6,*) '[entering interp_p_to_q_dcfft]'
         timing_tic = second()
      end if
      if (verbose.gt.0) then
         write(6,*) '% n_r: ',n_r
         write(6,*) '% n_w_: ',(n_w_(nr),nr=0,n_r-1)
         write(6,*) '% n_w_W_: ',(n_w_W_(nr),nr=0,n_r-1)
         write(6,*) '% n_A_W: ',n_A_W
         write(6,*) '% n_A: ',n_A
      end if
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,*) '% Allocating temporary storage'
      end if
      allocate(Z_q_(0:n_w_max-1))
      ic=0
      ic_W=0
      do nr=0,n_r-1
         if (n_w_(nr).gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling cp1_c16: '
            end if
            call cp1_c16(n_w_(nr),S_p_(ic),Z_q_)
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling dcfftf: '
            end if
            call dcfftf(n_w_(nr),Z_q_,Wsave_(ic_W))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling af1_c16: '
            end if
            call af1_c16(n_w_(nr),1.0d0*cmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*cmplx(0.0d0),Z_q_,S_q_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % setting ic ',ic,' to ic ',ic + n_w_(nr)
            end if
            ic = ic + n_w_(nr)
         end if
         if (verbose.gt.2) then
            write(6,*) '% % setting ic_W ',ic_W,' to ic_W ',ic_W +
     $           n_w_W_(nr)
         end if
         ic_W = ic_W + n_w_W_(nr)
      enddo
      if (verbose.gt.1) then
         write(6,*) ' % deallocating temporary storage'
      end if
      deallocate(Z_q_)
      if (verbose.gt.0) then
         timing_toc = second()
         write(6,'(A,F6.3)')
     $        ' [finished interp_p_to_q_dcfft]: total time: '
     $        ,timing_toc-timing_tic
      end if
      end
