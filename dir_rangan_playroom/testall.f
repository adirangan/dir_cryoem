c      Testing suite for K-space discretization routines
c
c      Test 1: call testfourier(...)
c
c      Make Cartesian grid in physical space and spherical
c      grid in Fourier space together with suitable quadratures.
c      Check that Forward transform followed by inverse transform
c      recovers function.
c     
c      Invokes: mkphysgrid, evalgrid, get_ntot, mkfouriergrid
c             getgridr, getspheregrid, getgridth,getgridph,
c             nufft3df90
c
c      Can be called with both regular and sparse grid options and with 
c      regular and Gaussian grids in the radial direction.
c-----------------------------------------------------------------------
c      Test 2: call testprojeval(...)
c
c      Tests projectio n and evaluation of spherical harmonic
c      expansions using tensor product or sparse grid.
c   
c      Invokes: getspheregrid, projshexp, evalshexp.
c
c      Test 3: call greatcircle_test.
c
c      a) makes all great circles for normals corresponding to all 
c      grid points on sphere.  
c      b) Get Cartesian coords of all such target points and evaluate
c         some potential directly.
c      c) Make interpoolation matrix from tensor product grid to all 
c      target points. 
c      Evaluate field at all targets via interpoaltion.
c      Then, evaluate field at all targets via subroutine sheval_greatcircs.
c      Check error

c      Invokes:  mkgreatcircles, mkshinterp, shinterp, sheval_greatcircs
c
c      Test 4: call  testlsq(...)

c      Check that Sinterp_adj * diag * Sinterp = Identity.
c
c      Create a bunch of points in sphere and grab data there.
c      Solve least squares problem for SH coefficients.
c      Tests: lsqsolvespharm
c
c      Invokes: getspheregrid, sheval_spheregrid, quadscale,
c      		sheval_spheregrid_adj, lsqsolvespharm
c
      implicit real *8 (a-h,o-z)
c
      parameter (ngridr=80)
      integer nlats(ngridr)
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8 source(3,10)
      real *8 x0y0z0(3,10),sigma
c
      complex *16 wavek
      complex *16, allocatable :: local(:)
      complex *16 charge(10)
      complex *16 eye
      external multaha
      external fgauss
c
      data eye/(0.0d0,1.0d0)/
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
c
c
c     set sigma and location of Gaussian sources (x0,y0,z0)
c     ns = number of sources
c
      ns = 2
      sigma = 1.0d-1
      x0y0z0(1,1) = 0.1d0
      x0y0z0(2,1) = 0.2d0
      x0y0z0(3,1) = -0.13d0
      x0y0z0(1,2) = -0.1d0
      x0y0z0(2,2) = -0.2d0
      x0y0z0(3,2) = 0.23d0
c
      rmax = 60.0d0
      ityper = 0
      itypep = 1
      ngrid = 60
      ngridth=30
      ngridc=20
c
c     do test1
      iftest1 = 1
      if (iftest1.eq.1) then
         call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)
         do i = 1,ngridr
            nlats(i) = nint(pi*xnodesr(i))
            if (nlats(i).lt.6) nlats(i) = 6
         enddo
c      
         call testfourier(ngrid,ngridr,rmax,ityper,itypep,
     1     nlats,fgauss,ns,x0y0z0,sigma,err,errback)
         call prin2(' l2 err in ffhat is *',err,1)
         call prin2(' l2 err in ff is *',errback,1)
      endif
c
      nterms = 16
      allocate(local((nterms+1)**2))
c
c
c     create source.
c
      wavek=1.1d-1 + eye*0.1d0
      ns = 1
      source(1,1)=1.0d0
      source(2,1)=0.1d0
      source(3,1)=-1.5d0
      charge(1)= 1.0d0 + eye*0.2d0
c
      rk = 1.0d0
      call testprojeval(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms)
c
      call  testgreatcircles(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms,ngridc)
c
      call  testlsq(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms)
c
      stop
      end
c
c

