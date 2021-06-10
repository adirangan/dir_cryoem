      subroutine testlsq(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms)

c
c      test shinterp_adj and least squares solver.
c      
c      Check that Sinterp (adjoint) *diag * Sinterp is identity.
c
c      Create a bunch of points on sphere and grab data there.
c      Solve least squares problem for SH coefficients. Check error.
c      Tests: lsqsolvespharm
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,ns),center(3),ztrg(3)
      integer, allocatable :: ngridps(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: sthetas(:),phsteps(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      real *8, allocatable :: phis(:)
      real *8, allocatable :: thetas(:)

      complex *16 wavek,eye,charge(ns)
      complex *16 opot,ofld(3)
      complex *16 local((nterms+1)**2)
      complex *16, allocatable :: phiall(:)
      complex *16, allocatable :: urand(:)
      complex *16, allocatable :: local2(:)
      complex *16, allocatable :: local3(:)
c
ccc      external multaha
c
c     check shinterp_adj
c
      pi=4.0d0*datan(1.0d0)
c
      ngridph = 2*ngridth
      allocate(xnodesth(ngridth))
      allocate(sthetas(ngridth))
      allocate(wtsth(ngridth))
      allocate(ngridps(ngridth))
      allocate(phsteps(ngridth))
      allocate(phiall(ngridth*ngridph))
      allocate(urand(ngridth*ngridph))
      allocate(local2((nterms+1)**2))
      allocate(local3((nterms+1)**2))
      allocate(phis(ngridth*ngridph))
      allocate(thetas(ngridth*ngridph))

      itensor = 0
      call getspheregrid(ngridth,itensor,xnodesth,
     1             sthetas,wtsth,ngridps,phsteps,numonsphere)

      call sheval_spheregrid(local,phiall,nterms,ngridth,ngridph,
     1     xnodesth)
      call prin2(' phiall is *',phiall,10)
      call quadscale(ngridth,ngridph,wtsth,phiall)
      call prin2(' phiall is *',phiall,10)
      call sheval_spheregrid_adj(nterms,ngridth,ngridph,xnodesth,
     1                       phiall,local2)
      np = (nterms+1)**2
      errl = 0.0d0
      stot = 0.0d0
      do i = 1,np
         errl = errl + abs(local(i)-local2(i))**2
         stot = stot + local(i)**2
      enddo
      call prin2(' err in S_adj D S = I *',dsqrt(errl/stot),1)
c
      nout =ngridth*ngridph
      next = 0
      call prinf(' np *',np,1)
      call prinf(' nout *',nout,1)
      do ii = 1,ngridth
      do jj = 1,ngridph
         next = next+1
         thet = pi*rand()
         phi = 2*pi*rand()
         phis(next) = phi
         thetas(next) = thet
         ctheta = dcos(thet)
         stheta = dsin(thet)
         cphi = dcos(phi)
         sphi = dsin(phi)
         ztrg(1) = stheta*cphi
         ztrg(2) = stheta*sphi
         ztrg(3) = ctheta
         call hpotfld3dall(iffld,source,charge,ns,ztrg,wavek,opot,ofld)
         urand(next) = opot
      enddo
      enddo
c
c     random grid of data is now defined by:
c     urand,thetas,phis,
c     Now, solve for coefficients (local3).
c     The parameters ngridth,ngridph,kord,xnodesth,wtshth
c     are used for applying the matrix in the conjugate
c     gradient iteration.
c
c
      t1 = second()
      kord = 5
      ifpr = 0
      if (ifpr.eq.1) call prin2(' urand is *',urand,2*nout)
      if (ifpr.eq.1) call prin2(' thetas is *',thetas,10)
      if (ifpr.eq.1) call prin2(' phis is *',phis,10)
      if (ifpr.eq.1) call prinf(' nout is *',nout,1)
      if (ifpr.eq.1) call prinf(' ngridph is *',ngridph,1)
      if (ifpr.eq.1) call prinf(' ngridth is *',ngridth,1)
      if (ifpr.eq.1) call prinf(' kord is *',kord,1)
      oversamp = 1.0d0
      nquad = oversamp*ngridth
      eps = 1.0d-6
      numit = 1000
      call lsqsolvespharm(urand,thetas,phis,nout,local3,nterms,
     1                    nquad,kord,eps,numit)
      t2 = second()
      call prin2('time for lsq is *',t2-t1,1)
ccc      call prin2(' local is *',local,2*np)
ccc      call prin2(' local3 is *',local3,2*np)
      errl = 0.0d0
      stot = 0.0d0
ccc      call prin2(' errl is *',dsqrt(errl),1)
      do i = 1,np
         errl = errl + abs(local2(i)-local3(i))**2
         stot = stot + local2(i)**2
      enddo
      call prin2(' err in LSQ sol *',dsqrt(errl/stot),1)
c
      return
      end
c
c
c
