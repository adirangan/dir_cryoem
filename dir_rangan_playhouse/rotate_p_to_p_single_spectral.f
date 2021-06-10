      subroutine rotate_p2p_fx_single(fftw_plan_frwd,fftw_plan_back,n_w
     $     ,fftw_in1_,fftw_out_,S_p_1_,gamma,M_p_1_)
c$$$  Assumes that M_p_1_ is the same size and dimensions as S_p_1_
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_w
      real *8 gamma
      integer *8 fftw_plan_frwd
      integer *8 fftw_plan_back
      complex *16 fftw_in1_(0:n_w-1),fftw_out_(0:n_w-1)
      complex *16 S_p_1_(0:n_w-1),M_p_1_(0:n_w-1)
      complex *16 C
      real *8 pi
      integer ic,nw,nq,q
      real *8 al2_c16_f
      real *8 l2_pre,l2_pos
      pi = 4.0*atan(1.0)

      if (verbose.gt.0) then
         write(6,'(A)') '[entering rotate_p2p_fx_single]'
      end if !if (verbose.gt.0) then

c$$$      We assume that C_ is preallocated as follows: ;
c$$$      allocate(C_(0:n_w_max-1));
c$$$      do nw=0,n_w_max-1
c$$$         q = nw - n_w_max/2
c$$$         C = cmplx( +cos(q*gamma) , -sin(q*gamma) )
c$$$         C_(nw) = C
c$$$      enddo !do nw=0,n_w_max-1
c$$$      This implies that C_(n_w_max/2 + q) = exp(-i*q*gamma)

      ic = 0
      if (n_w.gt.0) then            
         call cp1_c16(n_w,S_p_1_(ic),fftw_in1_(ic))
         if (verbose.gt.0) then
            write(6,'(A,F16.8)') ' fftw_in1_ l2: ' ,
     $           al2_c16_f(n_w,fftw_in1_(ic))
         end if                 !if (verbose.gt.0) then
         call dfftw_execute_(fftw_plan_frwd)
         call af1_c16(n_w,1.0d0*cmplx(1.0d0/dsqrt(1.0d0
     $        *n_w)),1.0d0*cmplx(0.0d0),fftw_out_(ic)
     $        ,fftw_out_(ic))
         if (verbose.gt.0) then
            write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $           al2_c16_f(n_w,fftw_out_(ic))
         end if                 !if (verbose.gt.0) then
         do nq=0,n_w-1
            q = nq
            if (q.gt.n_w/2-1) then 
               q = q - n_w
            end if              !if (q.ge.n_w/2-1) then 
c$$$            C = C_(n_w_max/2 + q)
            C = cmplx( +cos(q*gamma) , -sin(q*gamma) )
            fftw_out_(ic + nq) = fftw_out_(ic + nq) * C
         enddo                  !do nq=0,n_w-1
         if (verbose.gt.0) then
            write(6,'(A,F16.8)') ' fftw_out_ l2: ' ,
     $           al2_c16_f(n_w,fftw_out_(ic))
         end if                 !if (verbose.gt.0) then
         call dfftw_execute_(fftw_plan_back)
         call af1_c16(n_w,1.0d0*cmplx(1.0d0/dsqrt(1.0d0
     $        *n_w)),1.0d0*cmplx(0.0d0),fftw_in1_(ic)
     $        ,fftw_in1_(ic))
         if (verbose.gt.0) then
            write(6,'(A,F16.8)') ' fftw_in1_ l2: ' ,
     $           al2_c16_f(n_w,fftw_in1_(ic))
         end if                 !if (verbose.gt.0) then
         call cp1_c16(n_w,fftw_in1_(ic),M_p_1_(ic))
         if (verbose.gt.2) then
            write(6,*) '% % setting ic ',ic,' to ic ',ic + n_w
         end if                 !if (verbose.gt.2) then
         ic = ic + n_w
      end if                    !if (n_w.gt.0) then

      if (verbose.gt.0) then
         write(6,'(A)') '[finished rotate_p2p_fx_single]'
      end if !if (verbose.gt.0) then

      end
