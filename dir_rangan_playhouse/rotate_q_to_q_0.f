      subroutine rotate_q_to_q_0(n_r,n_w_,n_A,S_q_,gamma,M_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as S_q_
c$$$      Multiplication performed in place
      implicit none
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 gamma
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      real *8 pi
      integer nr,ic,n_w_max,nw,nq,q
      complex *16 C
      complex *16, allocatable :: C_(:)
      pi = 4.0*atan(1.0)

      n_w_max = n_w_(n_r-1)
      allocate(C_(0:n_w_max-1))
      do nw=0,n_w_max-1
         q = nw - n_w_max/2
         C = cmplx( +cos(q*gamma) , -sin(q*gamma) )
         C_(nw) = C
      enddo !do nw=0,n_w_max-1

      ic=0
      do nr=0,n_r-1
         do nq=0,n_w_(nr)-1
            q = nq
            if (q.gt.n_w_(nr)/2-1) then 
               q = q - n_w_(nr)
            end if              !if (q.ge.n_w_(nr)/2-1) then 
            C = C_(n_w_max/2 + q)
            M_q_(ic) = S_q_(ic)*C
            ic = ic + 1
         enddo
      enddo
      deallocate(C_)
      end
