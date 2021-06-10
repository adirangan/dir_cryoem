c
c      TESTING CODE FOR INTERP
c
c      Define simple function of (theta,phi) and collection of 
c      random points on sphere.
c
c      Note: vex, uout, phival are complex BUT
c            values is real
c
      implicit real *8 (a-h,o-z)
c
      integer iscale(0:5000)
      parameter (nmax=2000)
      integer, allocatable ::  icols(:,:)
      real *8, allocatable ::  phis(:)
      real *8, allocatable ::  thetas(:)
      real *8, allocatable ::  values(:,:)
      complex *16, allocatable ::  vex(:)
      complex *16, allocatable ::  uout(:)
      complex *16, allocatable ::  phival(:)
      real *8, allocatable ::  xnodesth(:)
      real *8, allocatable ::  wtsth(:)
      integer, allocatable ::  ngridps(:)
      real *8, allocatable ::  sthetas(:)
      real *8, allocatable ::  phsteps(:)
      real *8, allocatable ::  ynm(:,:)
c
      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      nout = 100
      nterms = 3
      allocate(phis(nout))
      allocate(thetas(nout))
      allocate(vex(nout))
      allocate(ynm(0:nterms,0:nterms))
      do i = 1,nout
         thet = pi*rand()
         phi = 2*pi*rand()
ccc         write(6,*) ' thet =',thet
ccc         write(6,*) ' phi =',phi
         phis(i) = phi
         thetas(i) = thet
         ctheta = dcos(thet)
         call ylgndr(nterms,ctheta,ynm)
         vex(i) = ynm(nterms,nterms)*dcos(nterms*phi)
      enddo
c
      kord = 7
      ngridth = 20
      ngridph = 40
      allocate(icols(nout,kord*kord))
      allocate(values(nout,kord*kord))
      allocate(uout(nout))
      allocate(phival(ngridph*ngridth))
      allocate(xnodesth(ngridth))
      allocate(wtsth(ngridth))
      allocate(ngridps(ngridth))
      allocate(sthetas(ngridth))
      allocate(phsteps(ngridth))
c
      itypep = 0
      call getspheregrid(ngridth,itypep,xnodesth,
     1             sthetas,wtsth,ngridps,phsteps,numonsphere)
c
      call prin2(' xnodesth is *',sthetas,ngridth)
      call prin2(' sthetas is *',sthetas,ngridth)
      call prin2(' wtsth is *',wtsth,ngridth)
      call prinf(' ngridps is *',ngridps,ngridth)
      call prin2(' phsteps is *',phsteps,ngridth)
      nnn = 0
      do jj = 1,ngridth
         ctheta = xnodesth(jj)
         stheta = sthetas(jj)
         call ylgndr(nterms,ctheta,ynm)
         do kk = 1,ngridps(jj)
            nnn=nnn+1
            phi = 2*pi*(kk-1)/ngridps(jj)
            phival(nnn) = ynm(nterms,nterms)*dcos(nterms*phi)
         enddo
      enddo
c
c     create interpolation matrix
c
      t1 = second()
      call mkshinterp(thetas,phis,nout,kord,ngridth,ngridph,values,
     1           icols)
      t2 = second()
      call prin2(' time for mkshinterp is *',t2-t1,1)
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
c
c     perform interpolation
c
      t1 = second()
      call shinterp(nout,kord,ngridth,ngridph,values,
     1           icols,phival,uout)
      t2 = second()
      call prin2(' vex is *',vex,10)
      call prin2(' uout is *',uout,10)
      call prin2(' time for shinterp is *',t2-t1,1)
c
c     compute error
c
      err = 0.0d0
      stot = 0.0d0
      do ii = 1,nout
         err = err + abs(vex(ii)-uout(ii))**2
         stot = stot + abs(vex(ii))**2
      enddo
      call prin2('l2 of function  is*', sqrt(stot),1)
      call prin2('rel l2 err from shinterp is*', sqrt(err/stot),1)
c
      stop
      end
c
c
c
c
c
C
C
      SUBROUTINE PRINM(MPOLE,NTERMS)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,-nterms:NTERMS)
      real *8 MPOLE2(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-nterms:NTERMS)
      INTEGER NTERMS
C
C     print out coefficients of multipole expansion
C
1000  FORMAT(6E12.5)
1001  FORMAT(/)
      DO 100 L = 0,NTERMS
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
100   CONTINUE
        return
C
C
C
C
      ENTRY PRINM0(MPOLE0,NTERMS)
      DO 200 L = 0,NTERMS
         WRITE(6,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
200   CONTINUE
      RETURN
C
C
      ENTRY PRINM2(MPOLE2,NTERMS)
      DO L = 0,NTERMS
         WRITE(6,1000)(MPOLE2(L,M),M=0,L)
         WRITE(13,1000)(MPOLE2(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
      ENDDO
c
c
      RETURN
      end
c
