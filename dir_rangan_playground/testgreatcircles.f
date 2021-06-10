!> Doxygen comment: ;\n
!>     greatcircle_test. ;\n
!>     a) makes all great circles for normals corresponding to all  ;\n
!>     grid points on sphere.   ;\n
!>     b) Get Cartesian coords of all such target points and evaluate ;\n
!>        some potential directly. ;\n
!>     c) Make interpoolation matrix from tensor product grid to all  ;\n
!>     target points.  ;\n
!>     d) Evaluate field at all targets via interpoaltion. ;\n
!>     e) Then, evaluate field at all targets via subroutine  ;\n
!>        shevalspherecircsp. ;\n
!>     Check errors ;\n
!>     Invokes:  mkgreatcircles, mkshinterp, shinterp, shevalspherecircsp ;\n
!>    INPUT: ;\n
!>    ngridth        number of latitude points used on sphere ;\n
!>    itypep         quadrature scheme in phi direction         ;\n
!>                     0 = same on all latitudes (which oversamples poles) ;\n
!>                     1 = adaptive (undersamples toward poles) ;\n
!>    wavek          frequency defining test Helmholtz field ;\n
!>    rk             radius of testing sphere ;\n
!>    source         locations of Helmholtz sources ;\n
!>    charge         strengths of Helmholtz sources ;\n
!>    ns             number of Helmholtz sources ;\n
!>    local          local spherical harmonic expansion of Helmholtz field ;\n
!>    nterms         order of local expansion  ;\n
!>    ngridc         number of points on each great circle ;\n
!> Doxygen comment: ;\n
      subroutine testgreatcircles(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms,ngridc)
c
c      greatcircle_test.
c
c      a) makes all great circles for normals corresponding to all 
c      grid points on sphere.  
c      b) Get Cartesian coords of all such target points and evaluate
c         some potential directly.
c      c) Make interpoolation matrix from tensor product grid to all 
c      target points. 
c      d) Evaluate field at all targets via interpoaltion.
c      e) Then, evaluate field at all targets via subroutine 
c         shevalspherecircsp.
c      Check errors
c
c      Invokes:  mkgreatcircles, mkshinterp, shinterp, shevalspherecircsp
c
c     INPUT:
c     ngridth        number of latitude points used on sphere
c     itypep         quadrature scheme in phi direction        
c                      0 = same on all latitudes (which oversamples poles)
c                      1 = adaptive (undersamples toward poles)
c     wavek          frequency defining test Helmholtz field
c     rk             radius of testing sphere
c     source         locations of Helmholtz sources
c     charge         strengths of Helmholtz sources
c     ns             number of Helmholtz sources
c     local          local spherical harmonic expansion of Helmholtz field
c     nterms         order of local expansion 
c     ngridc         number of points on each great circle
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,ns),center(3),ztrg(3)
      integer, allocatable :: ngridps(:)
      integer, allocatable :: icols(:,:)
      real *8, allocatable :: xx(:),yy(:),zz(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: sthetas(:),phsteps(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      real *8, allocatable :: xg(:,:)
      real *8, allocatable :: yg(:,:)
      real *8, allocatable :: zg(:,:)
      real *8, allocatable :: values(:,:)
      real *8, allocatable :: phis(:)
      real *8, allocatable :: thetas(:)

      complex *16 wavek,eye,charge(ns)
      complex *16 opot,ofld(3)
      complex *16 local((nterms+1)**2)
      complex *16, allocatable :: phiall(:,:)
      complex *16, allocatable :: phiall2(:,:)
      complex *16, allocatable :: phival(:)
      complex *16, allocatable :: phival2(:)
      complex *16, allocatable :: uout(:,:)
c   
      data eye/(0.0d0,1.0d0)/
      pi=4.0d0*datan(1.0d0)
c
c     print output flag (1 turns on)
c
      ifpr = 0
c
c      create grid on sphere of radius rk according to 
c      itypep
c      get Cartesian coordinates for spherical grid
c
      allocate(xnodesth(ngridth))
      allocate(sthetas(ngridth))
      allocate(wtsth(ngridth))
      allocate(ngridps(ngridth))
      allocate(phsteps(ngridth))

      call getspheregrid(ngridth,itypep,xnodesth,
     1             sthetas,wtsth,ngridps,phsteps,numonsphere)
      call prinf(' numonsphere is *',  numonsphere,1)
c
      allocate(xg(ngridc,numonsphere))
      allocate(yg(ngridc,numonsphere))
      allocate(zg(ngridc,numonsphere))
      allocate(phiall(ngridc,numonsphere))
      allocate(phiall2(ngridc,numonsphere))
      allocate(phis(ngridc*numonsphere))
      allocate(thetas(ngridc*numonsphere))
c
      call mkallgreatcircles(rk,ngridth,itypep,
     1           numonsphere,ngridc,xg,yg,zg)
c
c     get phi/theta values for all great circles
c     and define phiall to be exact sol at those points.
c
      iffld=1
      iincr = 0
      do ii = 1,numonsphere
      do ll = 1,ngridc
         ztrg(1) = xg(ll,ii)
         ztrg(2) = yg(ll,ii)
         ztrg(3) = zg(ll,ii)
         call hpotfld3dall(iffld,source,charge,ns,ztrg,wavek,opot,ofld)
         ctheta = ztrg(3)
         phi = datan2(ztrg(2),ztrg(1))
         if (phi.lt.0.0d0) phi = phi + 2.0d0*pi
         phiall(ll,ii) = opot
         iincr = iincr+1
         phis(iincr) = phi
         thetas(iincr) = dacos(ctheta)
         write(76,*) ctheta
         write(77,*) thetas(iincr)
         write(78,*) phis(iincr)
      enddo
      enddo
      call prin2('finished direct evaluation *', dsqrt(err/stot),0)
      nout = ngridc*numonsphere
      call prinf('nout is *', nout,1)
c
c     set order of interpolation
c
      kord = 7
c
c     create interpolation matrix using tensor product grid
c
      ngridph = 2*ngridth
      allocate(uout(ngridc,numonsphere))
      allocate(values(ngridc*numonsphere,kord*kord))
      allocate(icols(ngridc*numonsphere,kord*kord))
      t1 = second()
c
c
      call mkshinterp(thetas,phis,nout,kord,ngridth,ngridph,values,
     1           icols)
      t2 = second()
      call prin2(' time for mkshinterp is *',t2-t1,1)
c
c     get data on tensor product grid
c
      allocate(phival2(ngridth*ngridph))
      itensor = 0
      call evalshexp(local,nterms,ngridth,itensor,
     1          numonsphere,phival2)
ccc      call prin2(' phival2 is *',phival2,2*numonsphere)
c
c     perform checksum
c
      sum = 0.0d0
      do j = 1,kord*kord
         do i = 1,nout
            sum = sum + values(i,j)
         enddo
      enddo
      sum = sum/nout
      call prin2(' sum is is *',sum,1)
      call prin2(' sum should be 0.10000E+01 *',sum,0)
c
c     perform interpolation phival2 -> uout
c
      t1 = second()
      call shinterp(nout,kord,ngridth,ngridph,values,
     1           icols,phival2,uout)
      t2 = second()
      call prin2(' time for shinterp is *',t2-t1,1)
c
      t1 = second()
      ldp = ngridc
      call sheval_greatcircs(local,nterms,
     1           ngridth,itypep,numonsphere,ngridc,phiall2,ldp)
      t2 = second()
      call prin2('evaluate by spherical harmonics, time is *',t2-t1,1)
c
c     compare error for interp (uout) and for eval of SH exp (phiall2)
c
      err = 0.0d0
      stot = 0.0d0
      errsh = 0.0d0
      do ii = 1,numonsphere
      do ll = 1,ngridc
         errsh = errsh + abs(phiall(ll,ii)-phiall2(ll,ii))**2
         err = err + abs(phiall(ll,ii)-uout(ll,ii))**2
         stot = stot + abs(phiall(ll,ii))**2
      enddo
      enddo
      call prin2('err from shinterp is*', dsqrt(err/stot),1)
      call prin2('err from shevalspherecircsp is*', dsqrt(errsh/stot),1)
      return
      end
c
c
