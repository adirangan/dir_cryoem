      subroutine testprojeval(ngridth,itypep,wavek,rk,
     1           source,charge,ns,local,nterms)
c
c      Define complex-valed potential U using solution to Helmholtz
c      equation via subroutine hpotfld3dall.
c
c      1)  sample U on sphere.
c      2)  project on spherical harmonic expansions.
c      3)  evaluate back on grid and compare
c          with values from step 1.
c
      implicit real *8 (a-h,o-z)
      real *8 source(3,ns),center(3),ztrg(3)
      integer, allocatable :: ngridps(:)
      real *8, allocatable :: xx(:),yy(:),zz(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: sthetas(:),phsteps(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      complex *16 wavek,eye,charge(ns)
      complex *16 opot,ofld(3)
      complex *16, allocatable :: phival(:)
      complex *16, allocatable :: phival2(:)
      complex *16 local((nterms+1)**2)
c   
      data eye/(0.0d0,1.0d0)/
c
c     print output flag (1 turns on)
c
      ifpr = 0
c
c      create grid on sphere of radius rk
c      and get Cartesian coordinates for spherical grid
c
      allocate(xnodesth(ngridth))
      allocate(sthetas(ngridth))
      allocate(wtsth(ngridth))
      allocate(ngridps(ngridth))
      allocate(phsteps(ngridth))
c
      call getspheregrid(ngridth,itypep,xnodesth,
     1             sthetas,wtsth,ngridps,phsteps,numonsphere)
      if (ifpr.eq.1) call prinf(' numonsphere is *',numonsphere,1)
      if (ifpr.eq.1) call prin2(' xnodesth is *',xnodesth,ngridth)
      if (ifpr.eq.1) call prin2(' sthetas is *',sthetas,ngridth)
      if (ifpr.eq.1) call prin2(' phsteps is *',phsteps,ngridth)
c
      allocate(xx(numonsphere))
      allocate(yy(numonsphere))
      allocate(zz(numonsphere))
      allocate(phival(numonsphere))
      allocate(phival2(numonsphere))
c
      nnn=0
      do kk = 1,ngridth
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
c
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            xx(nnn) = rk*stheta*dcos(phi)
            yy(nnn) = rk*stheta*dsin(phi)
            zz(nnn) = rk*ctheta
         enddo
      enddo
c
c       direct calculation:
c
      iffld=1
      do kk = 1,numonsphere
         ztrg(1) = xx(kk)
         ztrg(2) = yy(kk)
         ztrg(3) = zz(kk)
         call hpotfld3dall(iffld,source,charge,ns,ztrg,wavek,opot,ofld)
         phival(kk) = opot
      enddo
      if (ifpr.eq.1) call prin2(' phival is *',phival,2*numonsphere)
c
c    project on spherical harmonics
c
ccc      allocate(local((nterms+1)**2))
c
      t1 = second()
      call projshexp(phival,numonsphere,ngridth,itypep,
     1       nterms,local)
      t2 = second()
      call prin2('time for proj is *',t2-t1,1)
c
c     evaluate spherical harmonic expansion back at grid points.
c      
      t1 = second()
      call evalshexp(local,nterms,ngridth,itypep,
     1          numonsphere,phival2)
      t2 = second()
      call prin2('time for sheval is *',t2-t1,1)
      if (ifpr.eq.1) call prin2(' phival2 is *',phival2,2*numonsphere)
c
      err = 0.0d0
      stot = 0.0d0
      do kk = 1,numonsphere
         err = err + abs(phival(kk)-phival2(kk))**2
         stot = stot + abs(phival(kk))**2
      enddo
      call prin2('l2 of function *', sqrt(stot),1)
      call prin2('l2 err, forwrd/back SH transform*', sqrt(err/stot),1)
      return
      end
c
