      subroutine interp_p_to_q_dcfft2(n_r,n_w_,n_A,S_p_,S_q_)
c$$$      allocates within each loop
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      complex *16 S_p_(0:n_A-1),S_q_(0:n_A-1)
      complex *16, allocatable :: Z_q_(:)
      complex *16, allocatable :: Wsave2_(:)
      integer ic,nr,n_w_max
c$$$      parameters for timing
      integer nt,n_t
      real *8 timing_tic,timing_toc      
      if (verbose.gt.0) then
         write(6,*) '[entering interp_p_to_q_dcfft2]'
         timing_tic = second()
      end if
      if (verbose.gt.0) then
         n_t = 10000
         write(6,*) 'running n_t times with n_t: ',n_t
      else
         n_t = 1
      end if
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,*) '% Allocating temporary storage'
      end if
      allocate(Z_q_(0:n_w_max-1))
      do nt=0,n_t-1
      ic=0
      do nr=0,n_r-1
         if (n_w_(nr).gt.0) then
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling copy_c16: '
            end if
            call cp1_c16(n_w_(nr),S_p_(ic),Z_q_)
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; allocating Wsave2_ '
            end if
            allocate(Wsave2_(0:4*n_w_(nr)+16-1))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling dcffti '
            end if
            call dcffti(n_w_(nr),Wsave2_)
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling dcfftf '
            end if
            call dcfftf(n_w_(nr),Z_q_,Wsave2_)
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; calling af1_c16 '
            end if
            call af1_c16(n_w_(nr),1.0d0*cmplx(1.0d0/dsqrt(1.0d0
     $           *n_w_(nr))),1.0d0*cmplx(0.0d0),Z_q_,S_q_(ic))
            if (verbose.gt.2) then
               write(6,*) '% % nr: ',nr,'; freeing Wsave2_ '
            end if
            deallocate(Wsave2_)
            ic = ic + n_w_(nr)
         end if
      enddo
      enddo
      if (verbose.gt.1) then
         write(6,*) ' % deallocating temporary storage'
      end if
      deallocate(Z_q_)
      if (verbose.gt.0) then
         timing_toc = second()
         write(6,'(A,F6.3)')
     $        ' [finished interp_p_to_q_dcfft2]: total time: '
     $        ,timing_toc-timing_tic
      end if
      end
