c  A version of lsqsolvespharm which solves a least squares problem
c  including a complex-valued factor (eg CTF) on each data point
c
c  The measured data f takes the form of 
c  f = A alpha = D T S alpha, where S is spherical harmonic transform,
c                             T is interpolation
c                             and D is CTF weights (diagonal mult by d_i)
c  The Normal equations are 
c
c     S^H T^H D^H D T S alpha = S^H T^H D^H f
c                      
c  which are solved with CG.
c
c  Solves the least-squares problem for a
c  truncated spherical harmonic expansion up to degree nterms, times weights,
c  matching complex data
c  given at arbitrary points (thetas,phis) on the sphere. Useful for cryo-EM to
c  get Fourier data on a k-shell given set of image rings at the same k.
c
c  Let S[c] := \sum_{n=0}^nterms \sum_{m=-n}^n c_{nm} Y_{nm}(z) e^{im\phi}
c
c  define the a function on S^2 that is the spherical harmonic expansion with
c  length-(nterms+1)^2 coefficient vector c. Let (theta,phi) define a point
c  on S^2.
c
c  Then the routine returns
c
c      \hat{c}  =     arg min   \sum_{i=1}^nout | b_i - d_i S[c](z_i,\phi_i)|^2
c                c in C^{(nterms+1)^2}
c
c  ie the least-squares solution to Ac=b where A is the matrix of spherical
c  harmonics evaluated at the data points then multiplied by the CTL weights d_i
c
c  It uses conjugate gradient on the normal equations and a fast algorithm 
c  that requires O(nout + nterms^3) work per iteration, where nout is number of 
c  data points.
c
c     INPUT:
c   
c     phiout:   input data (including CTF contribution), length nout
c     thetas:   input theta values, length nout
c     phis:     input phi values, length nout
c     nout:     number of data points
c     nterms:   order of spherical harmonic approximation (ie P)
c     nquad:    number of theta discretization points on regular grid
c               (used for interpolation)
c     k:        interpolation order of accuracy to irregular points
c     ctfw:     ctf weight for each data point, called d_i above.
c               Note: All points from the same great circle on the same
c               k-space sphere have the same ctf weight but we assume
c               each data point has its own ctf weight on input and pay
c               no attention to which point comes from which reference
c               image
c     multaha   matvec for CG iteration (passed external subroutine)
c
c     OUTPUT:
c
c     localp:   spherical harmonic coefficients in unrolled format
c               ie ordered (0,0), (1,-1),
c               (1,0), (1,1), (2,-2), ... (nterms,nterms).
c               Length is (nterms+1)^2.
c
c Things you can adjust in code & should be made arguments:
c   eps = CG tolerance
c
c Below is a wrapper that omits the multaha argument, which can thus be
c called from MATLAB
c
c Fortran files this code needs:
c   spharmroutsnew.f        (note: not "p" version)
c   ccongrlsq.f
c
c Old doc notes:
c  Matvec:  A * localp ->
c           (1) localp -> vals on grid VALS(nquad*nquadm) vector.
c           (2) vals are INTERP*VALS  (sparse nout * nquad*nquadm matrix)
c
c  MatvecH: a.  INTERP^H,  ( nquad*nquadm x nout sparse matrix)
c           b.  FORMMP from reg grid
c           (2) vals are INTERP*VALS.
c
c Greengard based on Barnett lsqsolvespharm MATLAB code.
c Doc updated ahb 3/1/16

      subroutine lsqctfsolve(phiout,thetas,phis,nout,localp,nterms,
     1           nquad,k,ctfw,multaha,eps,numit,niter)

      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,numit,nitermax,niter
      integer, allocatable :: icols(:,:)
      real *8, allocatable :: xnodes(:),wts(:),twork(:),sthetas(:)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 errs(1000),eps
      real *8, allocatable :: values(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phiout(nout)
      complex *16 ctfw(nout)
      complex *16, allocatable :: y(:)
      complex *16, allocatable :: ugrid(:)
      complex *16, allocatable :: phivals(:)
      complex *16, allocatable :: work(:)
      complex *16 a,z1,z2,z3
      external multaha

c     alex used for print output to matlab shell... (fails with std fortran)
      integer i
      character(len=80) str
c     integer*4, external :: mexPrintf
c     
c
      nquadm = 2*nquad
      allocate(xnodes(nquad))
      allocate(wts(nquad))
      allocate(sthetas(nquad))
      allocate(twork(nquad))
      call getgridth(nquad,twork,xnodes,sthetas,wts)

      np = (nterms+1)*(nterms+1)
      allocate(icols(nout,k*k))
      allocate(values(nout,k*k))
      allocate(y(np))
      allocate(ugrid(nquadm*nquad))
      allocate(phivals(nout))
      allocate(work(5*np+1000))
c
c     solving neq:   
c
      t1 = second()
ccc      call prinf(' calling mkshinterp nout is *',nout,1)
      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
      t2 = second()
c      write(str,*) 'interp build time =',t2-t1,achar(15)
c      i = mexPrintf(str)
c
      do i = 1,nout
         phivals(i) = dconjg(ctfw(i))*phiout(i)    
      enddo

ccc      call prinf(' calling shinterp_adj nout is *',nout,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phivals,ugrid)
ccc      call prinf(' calling shevalspherep_adj nout is *',nout,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
ccc      call prin2(' ugrid is *',ugrid,2*nquad*nquadm)
ccc      call prin2(' y is *',y,2*np)
c
      call ccongr2(ier,np,a,multaha,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,ctfw,y,eps,numit,localp,niter,
     1    errs,work)
c      call prinf(' ier *',ier,1)
c$$$      call prinf(' niter *',niter,1)
c$$$      call prin2(' errs *',errs,niter)

c     alex printing to matlab:
c      write(str,*) 'niter =',niter,achar(5)
c      i = mexPrintf(str)
      return
      end
c
c

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c    This does the A^H A apply. Other matvecs can be used by changing the
c    function name passed into lsqctfsolve:
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine multaha(a,x,y,n,nterms,nquad,nquadm,nout,k,icols,
     1                   xnodes,wts,values,phiout,ugrid,ctfw)
      implicit real *8 (a-h,o-z)
      real *8 icols(nout,k*k)
      real *8 xnodes(nquad)
      real *8 wts(nquad)
      real *8 values(nout,k*k)
      real *8 p3
      complex *16 ugrid(nquadm*nquad)
      complex *16 phiout(nout),z3
      complex *16 a,x(n),y(n),d
      complex *16 ctfw(nout)
c
ccc      call prinf(' calling shevalspherep *',n,1)
      call shevalspherep(x,ugrid,nterms,nquad,nquadm,xnodes)
ccc      call prinf(' calling shinterp *',n,1)
      call shinterp(nout,k,nquad,nquadm,values,icols,ugrid,phiout)
c
      do i = 1,nout
         phiout(i) = dconjg(ctfw(i))*ctfw(i)*phiout(i)    
      enddo
ccc      call prinf(' calling shinterp_adj *',n,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout,ugrid)
ccc      call prinf(' calling shevalspherep_adj *',n,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
c
      return
      end
c
c

ccccccccccccccccccccccccccccc debugging codes ccccccccccccccccccccccccccccccc
c
c     This code sets up spares interp matrix then does phiout=Ax and
c     ahax=A^H Ax on given coeff vecotr x (=cnm)
C     to allow testing against matlab.  Hardwired to CTF=1
C     See spharm.mw for checkmatvec MATLAB driver and comparelsqsolvers.m
c     do run the test.
c
C     barnett 6/22/16
      subroutine checkmatvecs(phiout,thetas,phis,nout,x,ahax,nterms,
     1           nquad,nquadm,k,xnodes,wts,ctfw)

      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,nitermax
      integer, allocatable :: icols(:,:)
      real *8 xnodes(nquad),wts(nquad)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 errs(1000),eps
      real *8, allocatable :: values(:,:)
      complex *16 x((nterms+1)*(nterms+1))
      complex *16 ahax((nterms+1)*(nterms+1))
      complex *16 phiout(nout)
      complex *16 ctfw(nout)
      complex *16, allocatable :: ugrid(:)
      complex *16, allocatable :: work(:)
      complex *16 a,z1,z2,z3
c
      np = (nterms+1)*(nterms+1)
      allocate(icols(nout,k*k))
      allocate(values(nout,k*k))
      allocate(ugrid(nquadm*nquad))
      allocate(work(5*np+1000))
c
c     set up interp matrix... (fills values,icols)
      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
c     do Ax ... (takes x -> ugrid -> phiout)
      call shevalspherep(x,ugrid,nterms,nquad,nquadm,xnodes)
      call shinterp(nout,k,nquad,nquadm,values,icols,ugrid,phiout)
c     do A^H on the Ax from above... (takes phiout -> ugrid -> ahax)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout,ugrid)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,ahax)
c
      return
      end
c
