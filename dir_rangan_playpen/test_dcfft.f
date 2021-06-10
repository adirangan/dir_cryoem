      subroutine test_dcfft(n_r,n_w_,n_w_W_,n_A_W,Wsave_,n_A,Wdata_)
      implicit none
      integer n_r,n_w_(0:n_r-1),n_w_W_(0:n_r-1),n_A,n_A_W
      complex *16 Wsave_(0:n_A_W-1)
      complex *16 Wdata_(0:n_A-1)
      complex *16, allocatable :: Z_(:)
      integer ic,nr,ic_W,n_w_max,ng
      real *8 pi,theta
      write(6,*) ' [entering test_dcfft] '
      write(6,*) '% n_r: ',n_r
      write(6,*) '% n_w_: ',(n_w_(nr),nr=0,n_r-1)
      write(6,*) '% n_w_W_: ',(n_w_W_(nr),nr=0,n_r-1)
      write(6,*) '% n_A_W: ',n_A_W
      write(6,*) '% n_A: ',n_A
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
c$$$      write(6,*) n_w_max
c$$$      write(6,*) n_w_W_(n_r-1)
      allocate(Z_(0:n_w_max))
      ic=0
      ic_W=0
      do nr=0,n_r-1
         if (n_w_(nr).gt.0) then
            write(6,*) nr,'filling Z_'
            do ng=0,n_w_(nr)-1
               theta = 2*pi*ng/n_w_(nr)
               Z_(ng) = cmplx( cos(theta) , sin(theta) )
            enddo
            write(6,*) nr,'calling dcfftf using Z_'
            call dcfftf(n_w_(nr),Z_,Wsave_(ic_W))
            if (n_w_(nr).gt.3) then
               write(6,'(A,F8.3,F8.3,F8.3,F8.3)') 'Z_ : ',real(Z_(0)) ,
     $              real(Z_(1))/n_w_(nr) ,real(Z_(2)), real(Z_(3))
            end if
            write(6,*) nr,'filling Wdata_'
            do ng=0,n_w_(nr)-1
               theta = 2*pi*ng/n_w_(nr)
               Wdata_(n_w_max-n_w_(nr) + ng) = cmplx( cos(theta) ,
     $              sin(theta) )
            enddo
            write(6,*) nr,'calling dcfftf using Wdata_'
            call dcfftf(n_w_(nr),Wdata_(n_w_max-n_w_(nr)),Wsave_(ic_W))
            if (n_w_(nr).gt.3) then
               write(6,'(A,F8.3,F8.3,F8.3,F8.3)') 'Wdata_ : '
     $              ,real(Wdata_(n_w_max-n_w_(nr)+0))
     $              ,real(Wdata_(n_w_max -n_w_(nr)+1))/n_w_(nr)
     $              ,real(Wdata_(n_w_max-n_w_(nr)+2)) ,
     $              real(Wdata_(n_w_max-n_w_(nr)+3))
            end if
            ic = ic + n_w_(nr)
         end if
         ic_W = ic_W + n_w_W_(nr)
      enddo
      deallocate(Z_)
      write(6,*) ' [finished test_dcfft] '
      end
