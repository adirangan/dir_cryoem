c
c  Fortran version of the MATLAB solver lsqsolvespharm (Alex B.):
c
c  Solves the least-squares problem for a
c  truncated spherical harmonic expansion up to degree P, matching complex data
c  given at arbitrary points (z,phi) on the sphere. Useful for cryo-EM to
c  get Fourier data on a k-shell given set of image rings at the same k.
c
c  Let S[c] := \sum_{n=0}^P \sum_{m=-n}^n c_{nm} Y_{nm}(z) e^{im\phi}
c
c  define the a function on S^2 that is the spherical harmonic expansion with
c  length-(P+1)^2 coefficient vector c. Let (z,phi) define a point on S^2.
c
c  Then the routine returns
c
c        \hat{c}  =     arg min     \sum_{i=1}^m | b_i - S[c](z_i,\phi_i)|^2
c                   c in C^{(P+1)^2}
c
c  ie the least-squares solution to Ac=b where A is the matrix of spherical
c  harmonics evaluated at the data points.
c
c  It uses conjugate gradient on the normal equations and a fast algorithm 
c  that requires O(m + P^3) work per iteration, where m is number of 
c  data points.
c
c Inputs:
c   u = complex length-m vector of values on points.
c   z,phi = length-m vectors of z and phi locations of data points on S^2.
c   P = maximum degree of spherical harmonic projections.
c   nth, nph = dimensions of regular output grid for SH transform 
c              from which interpolation is carried out.         
c   k = order of interpolation to irregular locations (so that k^2 is number of
c       nonzeros per row in sparse format)
c   eps = CG tolerance
c
c Outputs:
c   cnm = complex coefficients returned as column vector, ordered (0,0), (1,-1),
c         (1,0), (1,1), (2,-2), ... (P,P).  Length is (P+1)^2.
c
c  Matvec:  A * cnm ->
c           (1) cnm -> vals on grid VALS(nth*nph) vector.
c           (2) vals are INTERP*VALS  (sparse m * nth*nph matrix)
c
c  MatvecH: a.  INTERP^H,  ( nth*nph x m sparse matrix)
c           b.  FORMMP from reg grid
c           (2) vals are INTERP*VALS.
c
c
c
      subroutine lsqsolvespharm(phiout,thetas,phis,nout,localp,nterms,
     1           nquad,nquadm,k,xnodes,wts,multaha)
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,nitermax
      integer, allocatable :: icols(:,:)
      real *8 xnodes(nquad),wts(nquad)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 errs(1000),eps
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
c
ccc      call prinf(' calling mkshinterp nout is *',nout,1)
      call mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,icols)
ccc      call prinf(' calling shinterp_adj nout is *',nout,1)
      call shinterp_adj(nout,k,nquad,nquadm,values,icols,phiout,ugrid)
ccc      call prinf(' calling shevalspherep_adj nout is *',nout,1)
      call shevalspherep_adj(nterms,nquad,nquadm,xnodes,ugrid,y)
ccc      call prin2(' ugrid is *',ugrid,2*nquad*nquadm)
ccc      call prin2(' y is *',y,2*np)
c

      eps = 1d-12
      numit = 1000
      call ccongr2(ier,np,a,multaha,nterms,nquad,nquadm,nout,k,icols,
     1    xnodes,wts,values,phivals,ugrid,z3,y,eps,numit,localp,niter,
     1    errs,work)
      call prinf(' ier *',ier,1)
      call prinf(' niter *',niter,1)
      call prin2(' errs *',errs,niter)
      return
      end
c
c
c
c
c
      subroutine multaha(a,x,y,n,nterms,nquad,nquadm,nout,k,icols,
     1                   xnodes,wts,values,phiout,ugrid,z3)
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
c
c
