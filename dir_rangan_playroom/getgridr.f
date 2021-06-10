ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create radial grid out to rkmax.
c     
c     INPUT:
c
c     rkmax          outer radius of sphere (r in [0,rkmax])
c     ngridr         number of nodes in r
c     ityper         quadrature scheme in radial direction
c                       1 = Gaussian
c                       0 = uniform grid from 1..ngridr
c                           rad(i) = rkmax*i/ngridr
c
c     OUTPUT:
c
c     xnodesr        nodes in radial direction
c     wtsr           quadrature weights in radial direction
c
c----------------------------------------------------------------------
      implicit none
      integer i,ngridr,ityper
      real *8 u,v,rkmax
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8, allocatable :: b(:)
c
      if (ityper.eq.1) then
         allocate(b(ngridr))
c$$$         call gaussq(0, ngridr, 0, 0, 0, 0, b, xnodesr, wtsr)
         u = rkmax/2.0d0
         v = rkmax/2.0d0
         do i = 1,ngridr
            xnodesr(i) = u*xnodesr(i) + v
            wtsr(i) = xnodesr(i)*xnodesr(i)*wtsr(i)*rkmax/2.0d0
         enddo
      else
         do i = 1,ngridr
            xnodesr(i) = rkmax*(i)/ngridr
            wtsr(i) =  xnodesr(i)*xnodesr(i)*rkmax/ngridr
         enddo
      endif
c
      return
      end
