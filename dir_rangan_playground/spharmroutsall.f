!> Doxygen comment: ;\n
!>  Fortran library for spherical harmonic projection and evaluation ;\n
!>  on spherical grids (tensor product or quasiuniform/sparse grids) ;\n
!>  using packed storage format for expansion coefficients. ;\n
!>   ;\n
!>--------------------------------------------------------------------- ;\n
!>  projshexp               Compute spherical harmonic expansion on unit  ;\n
!>	            	      sphere of function tabulated on sphere. ;\n
!>  packshexp               Packs a complex 2d spherical expansion format  ;\n
!>			      to a 1d list ;\n
!>  unpackshexp             Unpacks a packed 1d list to a 2d spherical  ;\n
!>			      expansion format  ;\n
!>  evalshexp               Evaluates a spherical harmonic expansion on a  ;\n
!>		              unit sphere at discretization points  ;\n
!>  sheval_greatcircs       Evaluates a spherical harmonic expansion at  ;\n
!>   			      the equatorial circles for z-axis rotation to ;\n
!>			      all points on discretized sphere  ;\n
!>  mkallgreatcircles       Create Cartesian coordinates for equatorial  ;\n
!>			      circles on sphere of radius rk with z-axis  ;\n
!>			      rotated to all discretization nodes ;\n
!>  mkonegreatcircle        Create Cartesian coordinates for equatorial  ;\n
!>			      circle on sphere of radius rk with z-axis  ;\n
!>			      rotated to single discretization node ;\n
!>  sheval_greatcirc        Evaluates a spherical harmonic expansion at a  ;\n
!>			      single equatorial circle corresponding to one  ;\n
!>                            z-axis rotation  ;\n
!>  mkshinterp              Create interpolation matrix in sparse format ;\n
!>  getinterpsparse         Create row of interp. matrix in sparse format ;\n
!>  getlagrangewt           Create interpolation weights for getinterpsparse ;\n
!>  shinterp                Apply interpolation matrix in sparse format ;\n
!>  shinterp_adj            Apply adjoint of interp. matrix in sparse format ;\n
!>  sheval_spheregrid_adj   Apply adjoint of evalshexp  ;\n
!>  quadscale               Scales function on sphere by quadrature weights ;\n
!>  sheval_spheregrid       Evaluates a spherical harmonic expansion on a  ;\n
!>		              tensor product grid on unit sphere. Duplicates ;\n
!>                            evalshexp but calling sequence is different ;\n
!>                            and seems best to keep both.  ;\n
!>                            .......Perhaps to be made obsolete .... ;\n
!>  kspacegrid_to_model     takes discrete values on spherical grid and ;\n
!>                          converts to spherical harmonic representation ;\n
!>  model_to_kspacegrid     takes spherical harmonic representation and  ;\n
!>                          converts to values on spherical grid ;\n
!>---------------------------------------------------------------------- ;\n
!>  Alex Barnett, Leslie Greengard, Andras Pataki, Marina Spivak  ;\n
!>  various precursor libraries merged/updated here (12/26/16). ;\n
!> Doxygen comment: ;\n
C   Fortran library for spherical harmonic projection and evaluation
C   on spherical grids (tensor product or quasiuniform/sparse grids)
C   using packed storage format for expansion coefficients.
C
C----------------------------------------------------------------------
C   projshexp               Compute spherical harmonic expansion on unit 
C 	            	      sphere of function tabulated on sphere.
C   packshexp               Packs a complex 2d spherical expansion format 
C			      to a 1d list
C   unpackshexp             Unpacks a packed 1d list to a 2d spherical 
C			      expansion format 
C   evalshexp               Evaluates a spherical harmonic expansion on a 
C		              unit sphere at discretization points 
C   sheval_greatcircs       Evaluates a spherical harmonic expansion at 
C    			      the equatorial circles for z-axis rotation to
C			      all points on discretized sphere 
C   mkallgreatcircles       Create Cartesian coordinates for equatorial 
C			      circles on sphere of radius rk with z-axis 
C			      rotated to all discretization nodes
C   mkonegreatcircle        Create Cartesian coordinates for equatorial 
C			      circle on sphere of radius rk with z-axis 
C			      rotated to single discretization node
C   sheval_greatcirc        Evaluates a spherical harmonic expansion at a 
C			      single equatorial circle corresponding to one 
C                             z-axis rotation 
C   mkshinterp              Create interpolation matrix in sparse format
C   getinterpsparse         Create row of interp. matrix in sparse format
C   getlagrangewt           Create interpolation weights for getinterpsparse
C   shinterp                Apply interpolation matrix in sparse format
C   shinterp_adj            Apply adjoint of interp. matrix in sparse format
C   sheval_spheregrid_adj   Apply adjoint of evalshexp 
C   quadscale               Scales function on sphere by quadrature weights
C   sheval_spheregrid       Evaluates a spherical harmonic expansion on a 
C		              tensor product grid on unit sphere. Duplicates
C                             evalshexp but calling sequence is different
C                             and seems best to keep both. 
C                             .......Perhaps to be made obsolete ....
C   kspacegrid_to_model     takes discrete values on spherical grid and
C                           converts to spherical harmonic representation
C   model_to_kspacegrid     takes spherical harmonic representation and 
C                           converts to values on spherical grid
C----------------------------------------------------------------------
C   Alex Barnett, Leslie Greengard, Andras Pataki, Marina Spivak 
C   various precursor libraries merged/updated here (12/26/16).
C
C
C
C***********************************************************************
      subroutine projshexp(phival,numonsphere,ngridt,itypep,
     1           nterms,shexp)
C***********************************************************************
C     Compute spherical harmonic expansion on unit sphere
C     of function tabulated on sphere.
C---------------------------------------------------------------------
C     INPUT:
C
C     phival()       tabulated function
c     numonsphere    number of points used on sphere, defined by 
c                    ngridt,itypep
c     ngridt         number of latitude points used on sphere
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(ngridt*phi_over)
c                    if itypep=1
c                       nphi = nint(ngridt*phi_over*dsin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
C     nterms         order of spherical harmonic expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     shexp()   = coefficients of s.h. expansion
C
C     NOTE:
C
C     yrecursion.f produces Ynm with a nonstandard scaling:
C     (without the 1/dsqrt(4*pi)). Thus the orthogonality relation
C     is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C     In the first loop below, you see
C
C           marray(jj,m) = sum*phsteps(jj)/(4*pi)
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,ngridt,nquadm,numonsphere
      integer l,m,jj,kk
      integer, allocatable ::  ngridps(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
      real *8, allocatable ::  ynm(:,:)
c
      complex *16 phival(numonsphere)
      complex *16 shexp((nterms+1)*(nterms+1))
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(marray(ngridt,-nterms:nterms))
      allocate(xnodesth(ngridt))
      allocate(sthetas(ngridt))
      allocate(wtsth(ngridt))
      allocate(ngridps(ngridt))
      allocate(phsteps(ngridt))
c
      call getspheregrid(ngridt,itypep,xnodesth,
     1     sthetas,wtsth,ngridps,phsteps,numonsphere)

c            
c     initialize shexp exp to zero
c
      do ix = 1,(nterms+1)*(nterms+1)
         shexp(ix) = 0.0d0
      enddo
c
c     create marray (intermediate array)
c
      do m=-nterms,nterms
         nnn = 0
         do jj=1,ngridt
            emul = cdexp(imag*m*2.0d0*pi/ngridps(jj))
            sum = 0
            ephi = 1.0d0
            do kk = 1,ngridps(jj)
               nnn=nnn+1
               sum = sum + phival(nnn)*dconjg(ephi)
               ephi = ephi*emul
            enddo
            marray(jj,m) = sum*phsteps(jj)/(4*pi)
         enddo
      enddo
c
c     get shexp 
c
      do jj=1,ngridt
         cthetaj = xnodesth(jj)
         call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
            zmul = marray(jj,m)*wtsth(jj)
            do l=abs(m),nterms
               ix = l*(l+1) + m + 1
               shexp(ix) = shexp(ix) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo

      return
      end
C
c**********************************************************************
      subroutine packshexp(nterms,in,out)
C***********************************************************************
c     Packs a complex *16 2d local expansion array in a 1d array.
c     Alex Barnett 3/28/12   Leslie Greengard modified 12/25/16
C---------------------------------------------------------------------
c     INPUT:
c
c     nterms  degree of expansion
c     in(:,:) expansion in 2d matrix format
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     out(:)     coefficients packed in 1D array
C***********************************************************************
      implicit none
      integer i,n,m,nterms,stride
      complex *16 in(0:nterms,-nterms:nterms), out(*)
      
      i = 1
      stride = 1
      do n=0,nterms
         do m=-n,n
            out(i) = in(n,m)
            i = i+stride
         enddo
      enddo
      end

c**********************************************************************
      subroutine unpackshexp(nterms,in,out)
C***********************************************************************
c     Unpacks a complex *16 local expansion stored in 1d array
c     in standard 2d array format.
c     Alex Barnett 2/4/13     Leslie Greengard modified 12/25/16
C---------------------------------------------------------------------
c     INPUT:
c
c     nterms  degree of expansion
c     in(:)      coefficients packed in 1D array
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     out(:,:)     expansion in 2d matrix format
c
c**********************************************************************
      implicit none
      integer i,n,m,nterms,stride
      complex *16 out(0:nterms,-nterms:nterms), in((nterms+1)**2)
      
      i = 1
      stride = 1
      do n=0,nterms
         do m=-n,n
            out(n,m) = in(i)
            i = i+stride
         enddo
      enddo
      end
C
C
C
C***********************************************************************
      subroutine evalshexp(localp,nterms,ngridt,itypep,
     1           numonsphere,phival)
C***********************************************************************
C     This subroutine evaluates a spherical harmonic expansion on a 
C     unit sphere at discretization points defined by ngridt, itypep.
C---------------------------------------------------------------------
C     INPUT:
C
C     localp()       packed coefficients of s.h. expansion
C     nterms         order of spherical harmonic expansion
c     ngridt         number of latitude points used on sphere
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
C     numonsphere    total number of points on sphere
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival()       tabulated function
C***********************************************************************
C     converted to OPENMP, Barnett 6/22/16
C
      implicit real *8 (a-h,o-z)
      integer nterms,ngridt,itypep,numonsphere
      integer jj,m,n,ix,nnn
      integer, allocatable ::  ngridps(:)
      real *8 phi
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
      real *8, allocatable ::  ynm(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phival(numonsphere)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      data imag/(0.0d0,1.0d0)/

c     needed for output that makes it to MATLAB terminal
c      character(len=80) str
c      integer*4, external :: mexPrintf
c     needed for omp
c      integer OMP_GET_THREAD_NUM, omp_get_num_threads
c
c
      pi = 4.0d0*datan(1.0d0)
      allocate(ngridps(ngridt))
      allocate(phsteps(ngridt))
      allocate(sthetas(ngridt))
      allocate(xnodesth(ngridt))
      allocate(wtsth(ngridt))
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,ngridt))
c
      call getspheregrid(ngridt,itypep,xnodesth,
     1     sthetas,wtsth,ngridps,phsteps,numonsphere)
c
c$OMP PARALLEL DO
      do jj=1,ngridt
      do m=-nterms,nterms
         phitemp(m,jj) = 0.0d0
      enddo
      enddo
c$OMP END PARALLEL DO
c
cccc$OMP PARALLEL
c      if (omp_get_thread_num().eq.0) then
c         write(str,*) omp_get_num_threads(), ' threads',achar(10)
c         i=mexPrintf(str)
c      endif
      
c$OMP PARALLEL DO PRIVATE(ctheta,stheta,ynm,m,mabs,ix,n)
      do jj=1,ngridt
         ctheta = xnodesth(jj)
         stheta = sthetas(jj)
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               ix = n*(n+1) + m + 1
               phitemp(m,jj) = phitemp(m,jj) +
     1                localp(ix)*ynm(n,mabs)
            enddo
         enddo
      enddo
c$OMP END PARALLEL DO
c
cccc$OMP PARALLEL DO PRIVATE(kk,ephik,ephi,m)
c
      nnn = 0
      do jj = 1,ngridt
      do kk = 1,ngridps(jj)
         nnn=nnn+1
         phival(nnn) = 0.0d0
         phi = 2.0d0*pi*(kk-1)/ngridps(jj)
         ephik = cdexp(imag*phi)
         ephi = cdexp(-nterms*imag*phi)
         do m = -nterms,nterms
            phival(nnn) = phival(nnn) + phitemp(m,jj)*ephi
            ephi = ephi*ephik
         enddo
      enddo
      enddo
cccc$OMP END PARALLEL DO
      return
      end
C
cc
C***********************************************************************
      subroutine sheval_greatcircs(localp,nterms,
     1           ngridt,itypep,numonsphere,nquadc,phival,ldp)
C***********************************************************************
C     This subroutine evaluates a spherical harmonic expansion at
C     the equatorial circles for z-axis rotation to all points on 
C     discretized sphere with numonsphere points. 
C     nquadc is the number of points desired on each great circle.
C---------------------------------------------------------------------
C     INPUT:
C
C     localp()       coefficients of spherical harmonic exp. (packed)
C     nterms         number of terms in the orig. expansion
c     ngridt         number of latitude points used on sphere
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
C     numonsphere    total number of points on sphere
C
C     nquadc   : number of points on equator for each rotation
C     ldp      : leading dimenson of phival array (must be at least
C                nquadc)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival()  : (l,m) corresponds to lth node on equatorial circle
C                        for normal corresponding to (theta_m,phi_m).
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,ngridt,itypep,numonsphere,nquadc
      integer ii,ll,m,n,ix,nnn
      integer, allocatable ::  ngridps(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
      real *8, allocatable ::  ynm(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
ccc      complex *16 phival(nquadc,numonsphere)
      complex *16 phival(ldp,numonsphere)
      complex *16, allocatable :: phitemp(:,:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(phitemp(-nterms:nterms,nquadc,ngridt))
ccc      allocate(ynm(0:nterms,0:nterms))
      allocate(xnodesth(ngridt))
      allocate(sthetas(ngridt))
      allocate(wtsth(ngridt))
      allocate(ngridps(ngridt))
      allocate(phsteps(ngridt))
c
      call getspheregrid(ngridt,itypep,xnodesth,
     1     sthetas,wtsth,ngridps,phsteps,numonsphere)
c

      do ii=1,ngridt
      do ll=1,nquadc
      do m=-nterms,nterms
         phitemp(m,ll,ii) = 0.0d0
      enddo
      enddo
      enddo
ccc      call prin2(' initalized phitemp *',pi,0)
c
c$OMP PARALLEL PRIVATE(ynm)
c     each thread gets own allocation of ynm
      allocate(ynm(0:nterms,0:nterms))
C$OMP DO PRIVATE(ll,calphai,salphai,phil,cphil,ctheta,m,mabs,ix,n)
      do ii=1,ngridt
         calphai = xnodesth(ii)
         salphai = sthetas(ii)
         do ll=1,nquadc
            phil = 2.0d0*pi*(ll-1)/nquadc
            cphil = dcos(phil)
            ctheta = -salphai*cphil
            call ylgndr(nterms,ctheta,ynm)
            do m=-nterms,nterms
               mabs = abs(m)
               do n=mabs,nterms
                  ix = n*(n+1) + m + 1
                  phitemp(m,ll,ii) = phitemp(m,ll,ii) +
     1                localp(ix)*ynm(n,mabs)
               enddo
            enddo
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL

ccc      call prin2(' computed phitemp *',pi,0)
c

cccC$OMP PARALLEL DO PRIVATE(jj,ll,calphai,salphai,betaj,cbetaj,sbetaj,phil,cphil,&
cccC$OMP& sphil,x,y,phi,ephik,ephi,m)
      nnn = 0
      do ii = 1,ngridt
         calphai = xnodesth(ii)
         salphai = sthetas(ii)
         do jj = 1,ngridps(ii)
            betaj = 2.0d0*pi*(jj-1)/ngridps(ii)
            cbetaj = dcos(betaj)
            sbetaj = dsin(betaj)
            nnn=nnn+1
            do ll = 1,nquadc
               phil = 2.0d0*pi*(ll-1)/nquadc
               cphil = dcos(phil)
               sphil = dsin(phil)
               x = cbetaj*calphai*cphil-sbetaj*sphil
               y = sbetaj*calphai*cphil+cbetaj*sphil
               phi = datan2(y,x)
c         
               phival(ll,nnn) = 0.0d0
               ephik = cdexp(imag*phi)
               ephi = ephik**(-nterms)
               do m = -nterms,nterms
                  phival(ll,nnn)=phival(ll,nnn)+phitemp(m,ll,ii)*ephi
                  ephi = ephi*ephik
               enddo
ccc         if ( ((ii.eq.1).and.(jj.eq.1)).and.(ll.eq.1)) then
ccc         call prin2('111 phival is *', phival(ll,jj,ii),2)
ccc         endif
            enddo
         enddo
      enddo
cccC$OMP END PARALLEL DO

ccc      call prin2(' computed phival *',pi,0)
      return
      end
C
C
C
C***********************************************************************
      subroutine mkallgreatcircles(rk,ngridt,itypep,
     1           numonsphere,ngridc,xg,yg,zg)
C***********************************************************************
c
c     create Cartesian coordinates for equatorial circles on sphere of
c     radius rk with z-axis rotated to sparse grid locations.
c
C---------------------------------------------------------------------
c     INPUT:
c
c     rk             radius of sphere
c     ngridt         number of latitude points used on sphere
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
C     numonsphere    total number of points on sphere
C     ngridc     number of nodes on great circle
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     xg(l,m)  xcoord of lth node in equatorial circle for
c                z-axis rotated to mth sparse sphere point.
c     yg(l,m)  y-coord of lth node in equatorial circle for
c                z-axis rotated to mth sparse sphere point.
c     zg(l,m)  z-coord of lth node in equatorial circle for
c                z-axis rotated to mth sparse sphere point.
C***********************************************************************
      implicit none
      integer i,j,l,ngridc,ngridt,itypep,numonsphere
      integer nnn
      integer, allocatable ::  ngridps(:)
      real *8 pi,rk,calphai,salphai
      real *8 betaj,cbetaj,sbetaj,phil,cphil,sphil
      real *8 xg(ngridc,numonsphere)
      real *8 yg(ngridc,numonsphere)
      real *8 zg(ngridc,numonsphere)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
c
      pi = 4.0d0*datan(1.0d0)
      allocate(ngridps(ngridt))
      allocate(phsteps(ngridt))
      allocate(sthetas(ngridt))
      allocate(xnodesth(ngridt))
      allocate(wtsth(ngridt))
c
      call getspheregrid(ngridt,itypep,xnodesth,
     1     sthetas,wtsth,ngridps,phsteps,numonsphere)
c
      nnn=0
      do i = 1,ngridt
         calphai = xnodesth(i)
         salphai = sthetas(i)
         do j = 1,ngridps(i)
            betaj = 2.0d0*pi*(j-1)/ngridps(i)
            cbetaj = dcos(betaj)
            sbetaj = dsin(betaj)
            nnn=nnn+1
            do l = 1,ngridc
               phil = 2.0d0*pi*(l-1)/ngridc
               cphil = dcos(phil)
               sphil = dsin(phil)
               xg(l,nnn) = rk*(cbetaj*calphai*cphil-sbetaj*sphil)
               yg(l,nnn) = rk*(sbetaj*calphai*cphil+cbetaj*sphil)
               zg(l,nnn) = rk*(-salphai*cphil)
            enddo
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine mkonegreatcircle(rk,cthetause,phiuse,rgamma,
     1           ngridc,xg,yg,zg)
C***********************************************************************
c     Create Cartesian coordinates a single equatorial circle on 
c     sphere of radius rk with z-axis rotated to location 
c     (cthetause, phiuse).
C---------------------------------------------------------------------
c     INPUT:
c
c     rk         radius of sphere
c     ngridt     number of discretization nodes in theta in lab frame
c     xnodesth   discretization nodes in theta in lab frame
c     ngridp     number of discretization nodes in phi in lab frame
c     xnodesph   discretization nodes in phi in lab frame
c     ngridc     number of nodes on great circle
C---------------------------------------------------------------------
c     OUTPUT:
c
c     xg(l)  xcoord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     yg(l)  y-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     zg(l)  z-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c
C***********************************************************************
      implicit none
      integer l,ngridc
      real *8 cthetause, phiuse, rgamma
      real *8 pi,rk,calphai,salphai
      real *8 betaj,cbetaj,sbetaj,phil,cphil,sphil
      real *8 xg(ngridc)
      real *8 yg(ngridc)
      real *8 zg(ngridc)
c
      pi = 4.0d0*datan(1.0d0)
c
      calphai = cthetause
      salphai = dsqrt(1 - calphai**2)
      betaj = phiuse
      cbetaj = dcos(betaj)
      sbetaj = dsin(betaj)
      do l = 1,ngridc
         phil = 2.0d0*pi*(l-1)/ngridc
         phil = phil + rgamma
         cphil = dcos(phil)
         sphil = dsin(phil)
         xg(l) = rk*(cbetaj*calphai*cphil-sbetaj*sphil)
         yg(l) = rk*(sbetaj*calphai*cphil+cbetaj*sphil)
         zg(l) = rk*(-salphai*cphil)
      enddo

      return
      end
C
C***********************************************************************
      subroutine sheval_greatcirc(localp,nterms,cthetause,phiuse,
     1           rgamma,nquadc,phival)
C***********************************************************************
C     This subroutine evaluates a spherical harmonic expansion at
C     a single equatorial circle corresponding to z-axis rotation 
C     to (theta,phi), where ctheta = dcos(theta).
C     nquadc is the number of points desired on each great circle.
C---------------------------------------------------------------------
C     INPUT:
C
C     localp()  : coefficients of spherical harmonic exp. (packed)
C     nterms    : number of terms in the orig. expansion
C     cthetause : cos theta for y-axis rotation.
C     phiuse    : angle for extrinsic z-axis rotation.
C     rgamma    : rotation of reference image
C     nquadc    : number of output points 
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival()  : phival(l) corresponds to lth node on equatorial circle
C                for normal corresponding to (theta,phi) rotated
C                by gamma
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8, allocatable ::  ynm(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phival(nquadc)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m,ix
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquadc))
c
      do ll=1,nquadc
      do m=-nterms,nterms
         phitemp(m,ll) = 0.0d0
      enddo
      enddo
c
      do ll=1,nquadc
         phival(ll) = 0.0d0
      enddo
c
      do ll=1,nquadc
         calphai = cthetause
         salphai = dsqrt(1.0d0 - calphai**2)
         phil = 2.0d0*pi*(ll-1)/nquadc
         phil = phil + rgamma
         cphil = dcos(phil)
         ctheta = -salphai*cphil
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               ix = n*(n+1) + m + 1
               phitemp(m,ll) = phitemp(m,ll) +
     1                localp(ix)*ynm(n,mabs)
            enddo
         enddo
      enddo
c
      do ll = 1,nquadc
         calphai = cthetause
         salphai = dsqrt(1.0d0 - calphai**2)
         betaj = phiuse
         cbetaj = dcos(betaj)
         sbetaj = dsin(betaj)
         phil = 2.0d0*pi*(ll-1)/nquadc
         phil = phil + rgamma
         cphil = dcos(phil)
         sphil = dsin(phil)
         x = cbetaj*calphai*cphil-sbetaj*sphil
         y = (sbetaj*calphai*cphil+cbetaj*sphil)
         phi = datan2(y,x)
c         
         phival(ll) = 0.0d0
         ephik = cdexp(imag*phi)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(ll) = phival(ll) + phitemp(m,ll)*ephi
            ephi = ephi*ephik
         enddo
      enddo
      return
      end
C
C
C
C
C***********************************************************************
      subroutine mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,
     1           icols)
C***********************************************************************
c
c     Create interpolation matrix in simple sparse storage 
c     format. Since we are using a kxk interpolation stencil,
c     we only need to know the column index and interpolation
c     weights.
C---------------------------------------------------------------------
c     INPUT:
c    
c     thetas(),phis() :  coordinates of target points
c     nout            :  number of target points
c     k               :  interpolation order
c     nquad           :  number of nodes in theta
c     nquadm          :  number of nodes in phi
c
C---------------------------------------------------------------------
c     OUTPUT:
c    
c     values(:,:)     :  interpolation weights from reg grid.
c     icols (:,:)     :  unrolled indices of source points for interp.
c
c
c     For each target point, compute the next k*k interpolation
c     weights and increment values,icols.
C***********************************************************************
      implicit real *8 (a-h,o-z)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 values(nout,k*k)
      real *8, allocatable :: wtemp(:)
      integer icols(nout,k*k)
      integer, allocatable :: icoltemp(:) 
c
      pi = 4*datan(1.0d0)
      iwrite = 1
C$OMP PARALLEL PRIVATE(icoltemp,wtemp)
c     note each thread gets own copy of these local temp arrays...
      allocate(icoltemp(k*k)) 
      allocate(wtemp(k*k))
c     multithread the loop over filling sparse matrix rows...
C$OMP DO PRIVATE(j)
      do i = 1,nout
         call getinterpsparse(thetas(i),phis(i),k,nquad,nquadm,
     1                    icoltemp,wtemp)
         do j = 1,k*k
            icols(i,j) = icoltemp(j)
            values(i,j) = wtemp(j)
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine getinterpsparse(theta,phi,k,nquad,nquadm,icol,
     1           wtemp)
C***********************************************************************
c     Create one row of interpolation matrix in simple sparse storage 
c     format. Since we are using a kxk interpolation stencil,
c     we only need to know the column index and interpolation
c     weights
C---------------------------------------------------------------------
c     INPUT:
c    
c     theta,phi  : coordinates of target point
c     k          : interpolation order
c     nquad      : number of equspaced nodes in theta
c     nquadm     : number of equspaced node in phi
C---------------------------------------------------------------------
c     OUTPUT:
c
c     icol()     : source pt on sphere in unrolled vector format.
c                  Assuming phival is given in the 2d fortran format:
c                  phival(nquadm,nquad)
c     wtemp()    : corresponding interpolation weights
c
c     f(x,y) = sum_{j1,j_2} f(j_1h_1,j_2h_2) l_{j_1,j_2}(x,y)
c
c     l_{j_1,j_2}(x,y)=\prod_{k_1\neq j_1} (x_1-k_1h_1)/(j_1h_1-k_1h_1)
c                      * \prod_{k_2\neq j_2} (x_2-k_2h_2)/(j_2h_2-k_2h_2)
c
c     Need to grab correctly wrapped (NP,NT) values for "source pt" 
c     weight is computed without regard for wrapping
C***********************************************************************
      implicit real *8 (a-h,o-z)
      real *8 pi,theta,phi
      real *8 wtemp(k*k)
      integer icol(k*k)
c
      pi = 4*datan(1.0d0)
c
c     get nearest grid coords (np0,nt0)
c
      np0 = 1 + nint(nquadm*phi/(2.0d0*pi))
ccc      nt0 = 1 + nint(nquad*(theta + pi/(2*nquad))/pi)
      nt0 = nint(nquad*(theta + pi/(2*nquad))/pi)
c
c     need to compute offset from base point before thinking about
c     wrapping interpolation points...
c
      ht = pi/nquad
      hp = 2.0d0*pi/nquadm
      dp = phi - (np0-1)*hp
      dth = theta - ((2*nt0-1)*pi)/(2*nquad) 
c
      if (np0 .gt. nquadm) np0 = np0- nquadm
      if (nt0 .gt. nquad) nt0 = nt0-1
c
      if ( (np0.lt.1).or.(np0.gt.nquadm)) then
         call prinf(' np0 is *',np0,1)
      endif
      if ( (nt0.lt.1).or.(nt0.gt.nquad)) then
         call prinf(' nt0 is *',nt0,1)
      endif
c
c
      incr = 0
      do i = -k/2,k/2
      do j = -k/2,k/2
         incr = incr+1
         ip = np0 +i
         it = nt0 +j
         if (it .lt. 1) then
            ip = ip + nquadm/2
            it = 1 - it
         else if (it .gt. nquad) then
            ip = ip + nquadm/2
            it = 2*nquad + 1 - it
         endif
         if (ip .gt. nquadm) then
            ip = ip - nquadm 
         else if (ip .lt. 1) then
            ip = ip + nquadm 
         endif
         if ( (ip.lt.1).or.(ip.gt.nquadm)) then
            call prinf(' ip is *',ip,1)
         endif
         if ( (it.lt.1).or.(it.gt.nquad)) then
            call prinf(' it is *',it,1)
         endif
c
c     now, ip,it are indices of desired contributor to interpolant
c     value ->  lagrange poyl eval.
c
ccc         write(6,*) ' np0, nt0 ',np0, nt0
ccc         write(6,*) ' ip-1 x hp ',(ip-1)*hp
ccc         write(6,*) ' phi,theta ',phi,theta
ccc         if ( (abs(dp).gt.hp/2).or.(abs(dt).gt.ht/2)) then
ccc            write(6,*) ' dp,dth,ht,hp ',dp,dth,ht,hp
ccc         endif
ccc         write(6,*) ' phi,theta ',phi,theta
ccc         write(6,*) ' np0,nt0 ',np0,nt0
ccc         write(6,*) ' dp,dth,ht,hp ',dp,dth,ht,hp
ccc         write(6,*) ' k,i,j,ip,it ',k,i,j,ip,it
         call getlagrangewt(dth,dp,k,i,j,ht,hp,val)
ccc         write(6,*) ' val ',val
         icol(incr) = (it-1)*nquadm + ip
         wtemp(incr) = val               
      enddo
      enddo
ccc      pause
      return
      end
c
c
c
c
C***********************************************************************
      subroutine getlagrangewt(dth,dp,k,i,j,ht,hp,val)
C***********************************************************************
c     Compute interpolation weight at (dth,dp) from i,j node.
C---------------------------------------------------------------------
c     INPUT:
c
c     dth,dp     : offsets from base point (assumed to be origin).
c     k          : interpolation order
c     i,j        : particular grid point whose lagrange polynomial
c                  is evaluated at (dp,dth).
c     ht,hp      : grid spacing in theta and phi, respectively.
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     val        : value of Lagrange interpolant assoc with i,j point.
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer k,i,j
      real *8 dth,dp,ht,hp,val
c
      prodp = 1.0d0
      dploc = dp/hp
      dtloc = dth/ht
      do ip = -k/2,k/2
         if (i.ne.ip) prodp = prodp*(dp-ip*hp)/(i*hp-ip*hp)  
      enddo
c
      prodt = 1.0d0
      do jt = -k/2,k/2
         if (j.ne.jt) prodt = prodt*(dth-jt*ht)/(j*ht-jt*ht)  
      enddo
      val = prodp*prodt
      return
      end        
c
c
c
c
C***********************************************************************
      subroutine shinterp(nout,k,nquad,nquadm,values,
     1           icols,uval,uout)
C***********************************************************************
c     Apply interpolation matrix in sparse formt.
C---------------------------------------------------------------------
c     INPUT:
c
c     nout          : number of targets
c     k             : interpolation order
c     nquad         : number of nodes in theta
c     nquadm        : number of nodes in phi
c     values(:,:)   : sparse format interpolation weights
c     icols(:,:)    : sparse format col indices
c     uval()        : values on regular grid
C---------------------------------------------------------------------
c     OUTPUT:
c
c     uout()        : values at irregular nodes corresp. to 
c                     rows of values/icols
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer icols(nout,k*k)
      real *8 values(nout,k*k)
      complex *16 uout(nout)
      complex *16 uval(nquadm*nquad)
c
c
C$OMP PARALLEL DO PRIVATE(j)
      do i = 1,nout
         uout(i) = 0.0d0
         do j = 1,k*k
            uout(i) = uout(i) + values(i,j)*uval(icols(i,j))
         enddo
      enddo
C$OMP END PARALLEL DO
      return
      end
c
c
c
c
C***********************************************************************
      subroutine shinterp_adj(nout,k,nquad,nquadm,values,
     1           icols,uin,ugrid)
C***********************************************************************
c     Apply adjoint of interpolation matrix in sparse formt.
C---------------------------------------------------------------------
c     INPUT:
c
c     nout          : number of targets
c     k             : interpolation order
c     nquad         : number of nodes in theta
c     nquadm        : number of nodes in phi
c     values(:,:)   : sparse format interpolation weights
c     icols(:,:)    : sparse format col indices
c     uin()         : input vector (dimensioned as values on irregular grid)
c
c     OUTPUT:
c
c     ugrid()       : adjoint of interp matrix applied to uin
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer icols(nout,k*k)
      real *8 values(nout,k*k)
      complex *16 uin(nout)
      complex *16 ugrid(nquadm*nquad)
c
      ngrid = nquadm*nquad
      do i = 1,ngrid
         ugrid(i) = 0.0d0
      enddo
C     Note since ugrid is added to by different threads, need reduction...
C$OMP PARALLEL DO PRIVATE(i,iuse) REDUCTION(+:ugrid)
      do j = 1,nout
      do i = 1,k*k
         iuse = icols(j,i)
         ugrid(iuse) = ugrid(iuse)+values(j,i)*uin(j)
      enddo
      enddo
C$OMP END PARALLEL DO
      return
      end
c
c
C***********************************************************************
      subroutine sheval_spheregrid_adj(nterms,nquad,nquadm,xnodes,
     1           phival,localp)
C***********************************************************************
C     Compute adjoint of spherical harmonic expansion on unit 
C     sphere of function tabulated at nquadm*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C     nterms        : order of spherical harmonic expansion
C     nquad         : number of quadrature nodes in theta direction.
C     nquadm        : number of quadrature nodes in phi direction.
C     xnodes()      : Gauss-Legendre nodes x_j = cos theta_j
C     phival(:,:)   : tabulated function
C                      phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C     NOTE:    We assume phi_i = (i-1)*2.0d0*pi/nquadm, as do the 
C              routines in projection.f. However, we permit
C              different numbers of nodes in theta and phi.
C---------------------------------------------------------------------
C     OUTPUT:
C
C     localp()      : cadjoint of sheval applied to phival
C
C           This is different from projshexp since it omits the 
C           quadrature weights. Thus the lines below are commented out
C
Cccc        marray(jj,m) = sum*2.0d0*pi/nquadm
Cccc        marray(jj,m) = sum/(2*nquadm)
C
ccc         zmul = marray(jj,m)*wts(jj)
c
C     The latter has incorporated the 1/(4*pi) normalization factor
C     into the azimuthal quadrature weight (2.0d0*pi/nquadm).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm
      integer l,m,jj,kk,ix
      real *8 xnodes(nquad)
      real *8, allocatable ::  ynm(:,:)
      complex *16 phival(nquadm,nquad)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(marray(nquad,-nterms:nterms))
c
c     initialize localp exp to zero
c
      do ix = 1,(nterms+1)*(nterms+1)
         localp(ix) = 0.0d0
      enddo
c
c     create marray (intermediate array)
c
      do m=-nterms,nterms
         emul = cdexp(imag*m*2.0d0*pi/nquadm)
c     for simplicity we thread over jj since m loop has emul changing...
C$OMP PARALLEL DO PRIVATE(sum,ephi,kk)
         do jj=1,nquad
            sum = 0
            ephi = 1.0d0
            do kk = 1,nquadm
ccc               ephi = cdexp(imag*m*(kk-1)*2.0d0*pi/nquadm)
               sum = sum + phival(kk,jj)*dconjg(ephi)
               ephi = ephi*emul
            enddo
ccc         marray(jj,m) = sum*2.0d0*pi/nquadm
ccc         marray(jj,m) = sum/(2*nquadm)
            marray(jj,m) = sum
         enddo
C$OMP END PARALLEL DO
c     I don't understand what value emul1 has - it's never set! (ahb):
         emul = emul*emul1
      enddo
c
c     get local exp
c
C$OMP PARALLEL PRIVATE(ynm)
c     note each thread gets own allocation of ynm output from ylgndr
      allocate(ynm(0:nterms,0:nterms))
C$OMP DO PRIVATE(cthetaj,m,zmul,l,ix) REDUCTION(+:localp) 
      do jj=1,nquad
         cthetaj = xnodes(jj)
         call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
ccc            zmul = marray(jj,m)*wts(jj)
            zmul = marray(jj,m)
            do l=abs(m),nterms
               ix = l*(l+1) + m + 1
               localp(ix) = localp(ix) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
C$OMP END DO
C$OMP END PARALLEL
      return
      end
c
c
C***********************************************************************
      subroutine quadscale(nquad,nquadm,wts,phival)
C***********************************************************************
C     scale function phival on sphere by quadrature weights.
C---------------------------------------------------------------------
C     INPUT:
C
C     nquad       : number of quadrature nodes in theta direction.
C     nquadm      : number of quadrature nodes in phi direction.
C     wts()       : quadrature wrights in theta
C     phival(:,:) : tabulated function
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival(:,:) : OVER-WRITTEN by phival*(quadrature weight)
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad,nquadm
      integer jj,kk
      real *8 wts(1)
      complex *16 phival(nquadm,nquad)
C
      pi = 4.0d0*datan(1.0d0)
c
      do jj=1,nquad
         w = wts(jj)/(2*nquadm)
         do kk = 1,nquadm
            phival(kk,jj) = phival(kk,jj)*w
         enddo
      enddo
      return
      end
C
C
C
C
C***********************************************************************
      subroutine sheval_spheregrid(localp,phival,
     1           nterms,nquad,nquadm,xnodes)
C***********************************************************************
C     Evaluates a spherical harmonic expansion on a tensor product
C     nquad x nquadm grid on the unit sphere. 
C---------------------------------------------------------------------
C     INPUT:
C
C     localp      : coefficients of spherical harmonic exp. (packed)
C     nterms      : number of terms in the orig. expansion
C     nquad       : number of quadrature nodes in theta
C     nquadm      : number of quadrature nodes in phi
C     xnodes()    : Legendre nodes in theta (x_j = cos theta_j).
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival(:,:) : function value on tensor product
C                mesh on target sphere. phi is the fast variable, theta slow
C***********************************************************************
C     converted to OPENMP, Barnett 6/22/16
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 xnodes(nquad)
      real *8, allocatable ::  ynm(:,:)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phival(nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m
      data imag/(0.0d0,1.0d0)/

c     needed for output that makes it to MATLAB terminal
c      character(len=80) str
c      integer*4, external :: mexPrintf
c     needed for omp
c      integer OMP_GET_THREAD_NUM, omp_get_num_threads


      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquad))
c
c$OMP PARALLEL DO
      do jj=1,nquad
      do m=-nterms,nterms
         phitemp(m,jj) = 0.0d0
      enddo
      enddo
c$OMP END PARALLEL DO
c
cccc$OMP PARALLEL
c      if (omp_get_thread_num().eq.0) then
c         write(str,*) omp_get_num_threads(), ' threads',achar(10)
c         i=mexPrintf(str)
c      endif
      
c$OMP PARALLEL DO PRIVATE(ctheta,stheta,ynm,m,mabs,ix,n)
      do jj=1,nquad
         ctheta = xnodes(jj)
         stheta = dsqrt(1.0d0 - ctheta**2)
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               ix = n*(n+1) + m + 1
               phitemp(m,jj) = phitemp(m,jj) +
     1                localp(ix)*ynm(n,mabs)
            enddo
         enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(kk,ephik,ephi,m)
      do jj = 1,nquad
      do kk = 1,nquadm
         phival(kk,jj) = 0.0d0
         ephik = cdexp(2.0d0*pi*(kk-1)*imag/nquadm)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(kk,jj) = phival(kk,jj) + phitemp(m,jj)*ephi
            ephi = ephi*ephik
         enddo
      enddo
      enddo
c$OMP END PARALLEL DO
      return
      end
c
c
c
C
C
C***********************************************************************
      subroutine kspacegrid_to_model(ffhat,rkmax,nterms_sph,
     1                  ityper,ngridr,nlats,itypep,numonsphere,modsph,
     2                  nsphstore)
C***********************************************************************
C     Converts from discrete spherical grid to model: spherical harmonic
c     expansions on successive spheres.
C---------------------------------------------------------------------
C     INPUT:
C
C     ffhat      : samples of fucntion in k space on spherical grid
C     rkmax      : maximul radius in k space
c     nterms_sph : array of length ngridr defining the orders of the 
c                  spherical harmonic expansions on successive spheres.
c     ityper       quadrature scheme in radial direction
c                     1 = Gaussian
c                     0 = uniform grid from 1..ngridr
c                         rad(i) = rkmax*i/ngridr
c     ngridr       number of nodes in r
c     nlats()      number of latitude points used on each sphere
c                  permits coarser sampling near origin 
c     itypep       quadrature scheme in phi direction
c                     0 = same on all latitudes (which oversamples poles)
c                     1 = adaptive (undersamples toward poles)
c                  if itypep=0
c                     nphi = nint(nlats(i)*phi_over)
c                  if itypep=1
c                     nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                  phi_over is defined in getgridph and typically
c                  set to 2.
c     numonsphere  number of points used on sphere, defined by 
c                  ngridt,itypep
c     nsphstore    size of spherical harmonic model 
c                  (unrolled on all spheres)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     modsph(:)  : ffhat  expressed in spherical harmonics in successive
C                  spheres
C***********************************************************************

      implicit real *8 (a-h,o-z)
      integer ityper,ngridr,itypep
      integer i,nstart,nsphstart
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      complex *16 ffhat(*)
      complex *16 modsph(0:nsphstore-1)
c
c
      pi = 4.0d0*datan(1.0d0)
      nstart = 1
      nsphstart = 0
      do i = 1,ngridr
         call projshexp(ffhat(nstart),numonsphere(i),nlats(i),itypep,
     1       nterms_sph(i),modsph(nsphstart))
         nstart = nstart + numonsphere(i)
         nsph = (nterms_sph(i)+1)**2
         nsphstart = nsphstart + nsph
      enddo
      return
      end

c
C
C
C***********************************************************************
      subroutine model_to_kspacegrid(ffhat,rkmax,nterms_sph,
     1                  ityper,ngridr,nlats,itypep,numonsphere,modsph,
     2                  nsphstore)
c
C***********************************************************************
C     Evaluates model at discrete points: 
c     spherical harmonic expansions on successive 
c     spheres (nodsph) to discrete spherical grid.
C---------------------------------------------------------------------
C     INPUT:
C
C     modsph(:)  : ffhat  expressed in spherical harmonics in successive
C                  spheres
C
C     rkmax      : maximul radius in k space
c     nterms_sph : array of length ngridr defining the orders of the 
c                  spherical harmonic expansions on successive spheres.
c     ityper       quadrature scheme in radial direction
c                     1 = Gaussian
c                     0 = uniform grid from 1..ngridr
c                         rad(i) = rkmax*i/ngridr
c     ngridr       number of nodes in r
c     nlats()      number of latitude points used on each sphere
c                  permits coarser sampling near origin 
c     itypep       quadrature scheme in phi direction
c                     0 = same on all latitudes (which oversamples poles)
c                     1 = adaptive (undersamples toward poles)
c                  if itypep=0
c                     nphi = nint(nlats(i)*phi_over)
c                  if itypep=1
c                     nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                  phi_over is defined in getgridph and typically
c                  set to 2.
c     numonsphere  number of points used on sphere, defined by 
c                  ngridt,itypep
c     nsphstore    size of spherical harmonic model 
c                  (unrolled on all spheres)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     ffhat      : samples of function in k space on spherical grid
C***********************************************************************

      implicit real *8 (a-h,o-z)
      integer ityper,ngridr,itypep
      integer i,nstart,nsphstart
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      real *8, allocatable :: xnodesr(:),wtsr(:)
      complex *16 ffhat(*)
      complex *16 modsph(0:nsphstore-1)
c
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      call getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
c
      pi = 4.0d0*datan(1.0d0)
      nstart = 1
      nsphstart = 0
      do i = 1,ngridr
         call evalshexp(modsph(nsphstart),nterms_sph(i),nlats(i),
     1          itypep,numonsphere(i),ffhat(nstart))
         nstart = nstart + numonsphere(i)
         nsph = (nterms_sph(i)+1)**2
         nsphstart = nsphstart + nsph
      enddo
      return
      end

