      subroutine getquadwtsirreg(nout,k,nquad,nquadm,
     1           icols,wts,qwts)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 wts(nquad)
      real *8 qwts(nout)
      integer icols(nout,k*k)
      integer, allocatable :: isnap(:) 
c
c     From icols, determine approximate quadrature weight for each
c     irregular point. The algorithm for this is:
c
c     1) find nearest reg grid point.
c     2) count number N_g of irreg points assoc with that reg grid point.
c     3) assign quadrature weight to be wt for regular grid point divided 
c        by N_g.
c
c     INPUT:
c    
c     nout            :  number of target points
c     k               :  interpolation order
c     nquad           :  number of nodes in theta
c     nquadm          :  number of nodes in phi
c     icols           :  unrolled indices of source points for interp.
c     wts             :  quadrature weights for reg grid points (in theta).
c
c     OUTPUT:
c    
c     qwts            :  quadrature weights for irreg grid points 
c
c
      allocate(isnap(nquad*nquadm)) 
c
c     iss is "snap to grid node" for each irregular point.
c
      iss = (k*k+1)/2
c
      do i = 1,nquad*nquadm
         isnap(i) = 0
      enddo
c
      do i = 1,nout
         jj = icols(i,iss)
         isnap(jj) = isnap(jj)+1
      enddo
c
      do i = 1,nout
         jj = icols(i,iss)
         itheta = (jj-1)/nquadm + 1
         qwts(i) = wts(itheta)/(2*nquadm)
         qwts(i) = qwts(i)/isnap(jj)
      enddo
      return
      end
c
c
c

