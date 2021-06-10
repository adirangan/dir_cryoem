      subroutine dfftw_block_many_0(verbose,flag_fftw,fpm_howmany
     $     ,n_w_max,fpm_in1_,fpm_out_,fpm_back,n_X,A_)
c$$$      This subroutine passes blocks of the input array ;
c$$$      A_ to the plan stored in fpm_back, which is generated ;
c$$$      by fftw_plan_many. ;
c$$$      We assume that A_ is stored with ordering ;
c$$$      A_(nx + nw*n_X) ;
c$$$      with the goal of performing the fft on the nw-dimension ;
c$$$      (i.e., the most slowly varying index). ;
c$$$      We also assume that the plan fpm_back has been ;
c$$$      generated with the ordering ;
c$$$      fpm_out_(nw + nx*n_w_max) ;
c$$$      where the nw-dimension is the most quickly varying ;
c$$$      index. ;
c$$$      Due to these assumptions, we rearrange blocks of A_ ;
c$$$      when we copy them into fpm_out_, and we rearrange them ;
c$$$      once again when we copy them back out from fpm_in1_. ;
c$$$       ;
c$$$      timing: ;
c$$$      If flag_fftw.eqv..false., then we do not perform the fftw at all, ;
c$$$      and instead simply copy the array A_ into and out of ;
c$$$      fpm_out_ and fpm_in1_ respectively. ;
      
      implicit none
      integer verbose
      logical flag_fftw
      integer fpm_howmany,n_w_max,n_X
      complex *16 fpm_in1_(0:0),fpm_out_(0:0)
      integer *8 fpm_back
      complex *16 A_(0:0)
      integer n_block,nblock,nx_bgn,nx_end,nx_tot
      integer nx,nw,na,nf
      if (verbose.gt.0) then
         write(6,'(A)') '[entering dfftw_block_many_0]'
      end if !if (verbose.gt.0) then
      n_block = n_X/fpm_howmany
      if (mod(n_X,fpm_howmany).gt.0) then
         n_block = n_block + 1
      end if !if (mod(n_X,fpm_howmany).gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' fpm_howmany: ' , fpm_howmany
         write(6,'(A,I0)') ' n_w_max: ' , n_w_max
         write(6,'(A,I0)') ' n_X: ' , n_X
         write(6,'(A,I0)') ' n_block: ' , n_block
      end if !if (verbose.gt.1) then
      do nblock=0,n_block-1
         nx_bgn = nblock*fpm_howmany
         nx_end = min(n_X,(nblock+1)*fpm_howmany)
         nx_tot = nx_end-nx_bgn
         if (verbose.gt.2) then
            write(6,'(A,I0)') ' nblock: ' , nblock
            write(6,'(A,I0)') ' nx_bgn: ' , nx_bgn
            write(6,'(A,I0)') ' nx_end: ' , nx_end
            write(6,'(A,I0)') ' nx_tot: ' , nx_tot
         end if !if (verbose.gt.2) then
         if (nx_tot.gt.0) then
            nf = 0;
            do nx=nx_bgn,nx_end-1
               na = nx
               do nw=0,n_w_max-1
c$$$                  nf = nw + nx*n_w_max
c$$$                  na = nx + nw*n_X
                  fpm_out_(nf) = A_(na)
                  nf = nf + 1
                  na = na + n_X
               enddo !do nw=0,n_w_max-1
            enddo !do nx=nx_bgn,nx_end-1
            if (flag_fftw.eqv..true.) then
               call dfftw_execute_dft(fpm_back,fpm_out_,fpm_in1_)
            end if ! if (flag_fftw.eqv..true.) then
            if (flag_fftw.eqv..false.) then
               call cp1_c16(n_w_max*nx_tot,fpm_out_,fpm_in1_)
            end if ! if (flag_fftw.eqv..false.) then
            nf = 0;
            do nx=nx_bgn,nx_end-1
               na = nx
               do nw=0,n_w_max-1
c$$$                  nf = nw + nx*n_w_max
c$$$                  na = nx + nw*n_X
                  A_(na) = fpm_in1_(nf)
                  nf = nf + 1
                  na = na + n_X
               enddo !do nw=0,n_w_max-1
            enddo !do nx=nx_bgn,nx_end-1            
         end if !if (nx_tot.gt.0) then
      enddo !do nblock=0,n_block-1      
      if (verbose.gt.0) then
         write(6,'(A)') '[finished dfftw_block_many_0]'
      end if !if (verbose.gt.0) then
      end
