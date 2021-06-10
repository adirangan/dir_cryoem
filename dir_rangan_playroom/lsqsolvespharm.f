c
c  Fortran version of the MATLAB solver lsqsolvespharm (Alex B.).
c
c  Subsumed by CTF version - see lsqctfsolvespharm.f
c
c  Solves the least-squares problem for a
c  truncated spherical harmonic expansion up to degree nterms,
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
c        \hat{c}  =     arg min     \sum_{i=1}^nout | b_i - S[c](z_i,\phi_i)|^2
c                   c in C^{(nterms+1)^2}
c
c  ie the least-squares solution to Ac=b where A is the matrix of spherical
c  harmonics evaluated at the data points.
c
c  It uses conjugate gradient on the normal equations and a fast algorithm 
c  that requires O(nout + nterms^3) work per iteration, where nout is number of 
c  data points.
c
c     INPUT:
c   
c     phiout:   input data, length nout
c     thetas:   input theta values, length nout
c     phis:     input phi values, length nout
c     nout:     number of data points
c     nterms:   order of spherical harmonic approximation (P)
c     nquad:    number of theta discretization points on regular grid
c               (used for interpolation)
c     nquadm:   number of phi discretization points on regular grid
c               (used for interpolation)
c     k:        interpolation order of accuracy to irregular points
c     xnodes:   grid locations in theta
c     wts:      quadrature weights in theta (currently not used).
c     multaha   matvec function for CG iteration (passed external subroutine)
c
c     OUTPUT:
c
c     localp:   spherical harmonic coefficients in unrolled format
c               ie ordered (0,0), (1,-1),
c               (1,0), (1,1), (2,-2), ... (P,P).  Length is (nterms+1)^2.
c
c
c  Matvec:  A * localp ->
c           (1) localp -> vals on grid VALS(nquad*nquadm) vector.
c           (2) vals are INTERP*VALS  (sparse nout * nquad*nquadm matrix)
c
c  MatvecH: a.  INTERP^H,  ( nquad*nquadm x nout sparse matrix)
c           b.  FORMMP from reg grid
c           (2) vals are INTERP*VALS.
c
c Things you can adjust in code:
c   eps = CG tolerance
c
c Below is a wrapper that omits the multaha argument, which can thus be
c called from MATLAB
c
c Doc updated ahb 1/3/16
      subroutine lsqsolvespharm(phiout,thetas,phis,nout,localp,nterms,
     1           nquad,nquadm,k,xnodes,wts,multaha)
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,nitermax
      integer, allocatable :: icols(:,:)
      real *8 xnodes(nquad),wts(nquad)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 errs(100000),eps
      real *8, allocatable :: values(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phiout(nout)
      complex *16, allocatable :: y(:)
      complex *16, allocatable :: ugrid(:)
      complex *16, allocatable :: phivals(:)
      complex *16, allocatable :: work(:)
      complex *16 a,z1,z2,z3
      external multaha
c
c
      np = (nterms+1)*(nterms+1)
      allocate(icols(nout,k*k))
      allocate(values(nout,k*k))
      allocate(y(np))
      allocate(ugrid(nquadm*nquad))
      allocate(phivals(nout))
      allocate(work(5*np+1000))

c
c     Operator is: A =  T S, where S is spherical harmonic transform
c                            and T is interpolation
c     For CGN: T S alpha = f -> conjg(S') T' T S alpha = conjg(S') T' f
c                      
c     solving neq:   
c
c      call prinf(' calling mkshinterp nout is *',nout,1)
      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
ccc      call prinf(' calling shinterp_adj nout is *',nout,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout,ugrid)
ccc      call prinf(' calling shevalspherep_adj nout is *',nout,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
ccc      call prin2(' ugrid is *',ugrid,2*nquad*nquadm)
ccc      call prin2(' y is *',y,2*np)
c
      eps = 1d-06
      numit = 1000
      call ccongr2(ier,np,a,multaha,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,z3,y,eps,numit,localp,niter,
     1    errs,work)
c      call prinf(' ier *',ier,1)
c      call prinf(' niter *',niter,1)
c      call prin2(' errs *',errs,niter)
      return
      end
c
c
c
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c    This does the A^H A apply. Other matvecs can be used by changing the
c    function name passed into lsqsolvespharm:
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine multaha(a,x,y,n,nterms,nquad,nquadm,nout,k,icols,
     1                   xnodes,wts,values,phiout,ugrid,z3)
c     xnodes are cos(theta_j) grid (used only in shevalspherep_adj in
c       spharmroutsnewp.f)
c     note: wts are not used
      implicit real *8 (a-h,o-z)
      real *8 icols(nout,k*k)
      real *8 xnodes(nquad)
      real *8 wts(nquad)
      real *8 values(nout,k*k)
      real *8 p3
      complex *16 ugrid(nquadm*nquad)
      complex *16 phiout(nout),z3
      complex *16 a,x(n),y(n),d
c
ccc      call prinf(' calling shevalspherep *',n,1)
      call shevalspherep(x,ugrid,nterms,nquad,nquadm,xnodes)
ccc      call prinf(' calling shinterp *',n,1)
      call shinterp(nout,k,nquad,nquadm,values,icols,ugrid,phiout)
c
ccc      call prinf(' calling shinterp_adj *',n,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout,ugrid)
ccc      call prinf(' calling shevalspherep_adj *',n,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
c
      return
      end



ccccccccccccccccccccccccccccccc
c     Wrapper which fixes a particular multaha (ie the one above).
c     Needed since Matlab can't pass in a fortran subroutine.    Barnett 3/1/16

      subroutine lsqsolvespharmwrap(phiout,thetas,phis,nout,localp,
     1    nterms,nquad,nquadm,k,xnodes,wts)
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm
      real *8 xnodes(nquad),wts(nquad)
      real *8 thetas(nout)
      real *8 phis(nout)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phiout(nout)

      external multaha
C     put your choice of multaha func as last argument here...
      call lsqsolvespharm(phiout,thetas,phis,nout,localp,nterms,
     1           nquad,nquadm,k,xnodes,wts,multaha)
      return
      end
