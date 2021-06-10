C     Fortran library for sph harm projection and eval on grids & circles
C
C     Stuff from Leslie/Zydrunas FMMLIB3D, and localexp3d.
C     pulled by Barnett 7/31/15
C     modified by Greengard 8/1/15
c
c     Note there is also a "p" version of this lib with packed cnm storage.
C
C***********************************************************************
C     calling sequences changed from Barnett's spharmrouts.f 
C***********************************************************************
c
C***********************************************************************
      subroutine projloc3d(nterms,ldl,nquad,nquadm,xnodes,wts,
     1           phival,local)
C***********************************************************************
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           ldl    = dimension parameter for local expansion
C           nquad  = number of quadrature nodes in theta direction.
C           nquadm = number of quadrature nodes in phi direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           NOTE:    We assume phi_i = (i-1)*2*pi/nquadm, as do the 
C                    routines in projection.f. However, we permit
C                    different numbers of nodes in theta and phi.
C***********************************************************************
C     OUTPUT:
C
C           local = coefficients of s.h. expansion
C
C     NOTE:
C
C     yrecursion.f produces Ynm with a nonstandard scaling:
C     (without the 1/sqrt(4*pi)). Thus the orthogonality relation
C     is
C             \int_S  Y_nm Y_n'm'*  dA = delta(n) delta(m) * 4*pi. 
C
C     In the first loop below, you see
C
Cccc        marray(jj,m) = sum*2*pi/nquadm
C           marray(jj,m) = sum/(2*nquadm)
C
C     The latter has incorporated the 1/(4*pi) normalization factor
C     into the azimuthal quadrature weight (2*pi/nquadm).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm,ldl
      integer l,m,jj,kk
      real *8 wts(1),xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 phival(nquadm,nquad)
      complex *16 local(0:ldl,-ldl:ldl)
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(marray(nquad,-nterms:nterms))
c
c     initialize local exp to zero
c
      do l = 0,ldl
         do m = -l,l
            local(l,m) = 0.0d0
         enddo
      enddo

c
c     create marray (intermediate array)
c
      do m=-nterms,nterms
         emul = cdexp(imag*m*2*pi/nquadm)
         do jj=1,nquad
            sum = 0
            ephi = 1.0d0
            do kk = 1,nquadm
ccc               ephi = cdexp(imag*m*(kk-1)*2*pi/nquadm)
               sum = sum + phival(kk,jj)*dconjg(ephi)
               ephi = ephi*emul
            enddo
ccc         marray(jj,m) = sum*2*pi/nquadm
            marray(jj,m) = sum/(2*nquadm)
         enddo
      enddo
c
c     get local exp
c
      do jj=1,nquad
         cthetaj = xnodes(jj)
         call ylgndr(nterms,cthetaj,ynm)
         do m=-nterms,nterms
            zmul = marray(jj,m)*wts(jj)
            do l=abs(m),nterms
               local(l,m) = local(l,m) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo

      return
      end



C     from localexp3d/h3dwrappers.f =====================================

c**********************************************************************
      subroutine flattenlocexpz(nterms,in,out,stride)
c     unpacks a complex *16 2d local expansion array (in) into a 1d list (out).
c     Output array is written to using given stride.
c     Alex Barnett 3/28/12
      implicit none
      integer i,n,m,nterms,stride
      complex *16 in(0:nterms,-nterms:nterms), out(*)
      
      i = 1
      do n=0,nterms
         do m=-n,n
            out(i) = in(n,m)
            i = i+stride
         enddo
      enddo
      end

c**********************************************************************
      subroutine stacklocexpz(nterms,in,out,stride)
c     packs a complex *16 1d local expansion list (in) into a 2d array (out).
c     Input array is read using given stride.
c     Alex Barnett 2/4/13
      implicit none
      integer i,n,m,nterms,stride
      complex *16 out(0:nterms,-nterms:nterms), in((nterms+1)**2)
      
      i = 1
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
      subroutine shevalsphere(local,phival,
     1           nterms,lmp,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion on an 
C     nquad x nquadm grid on the unit sphere. 
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local    : coefficients of spherical harmonic exp.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for local
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : function value on tensor product
C                mesh on target sphere. phi is the fast variable, theta slow
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
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
c$OMP PARALLEL PRIVATE(ynm)
c     note each thread gets own allocation of ynm output from ylgndr
      allocate(ynm(0:nterms,0:nterms))
C$OMP DO PRIVATE(ctheta,stheta,m,mabs,ix,n)
      do jj=1,nquad
         ctheta = xnodes(jj)
         stheta = dsqrt(1.0d0 - ctheta**2)
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,jj) = phitemp(m,jj) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
C$OMP END DO
c$OMP END PARALLEL
c
c$OMP PARALLEL DO PRIVATE(kk,ephik,ephi,m)
      do jj = 1,nquad
      do kk = 1,nquadm
         phival(kk,jj) = 0.0d0
         ephik = cdexp(2*pi*(kk-1)*imag/nquadm)
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
C
C
C
C
C***********************************************************************
      subroutine shevalspherecircs(local,phival,
     1           nterms,lmp,nquad,nquadm,xnodes,nquadc)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion at
C     the equatorial circles for z-axis rotation to all standard spherical
C     nodes (on an nquad x nquadm grid). nquadc is the number of points
C     desired on each great circle.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local    : coefficients of spherical harmonic exp.
C     nterms   : number of terms in the orig. expansion
C     lmp      : dimension parameter for local
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C     nquadc   : number of points on equator for each rotation
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : (l,j,i) corresponds to lth node on equatorial circle
C                        for normal corresponding to (theta_i,phi_j).
C     NOTE the ordering of output grids.
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadc,nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquadc,nquad))
c
      do ii=1,nquad
      do ll=1,nquadc
      do m=-nterms,nterms
         phitemp(m,ll,ii) = 0.0d0
      enddo
      enddo
      enddo
ccc      call prin2(' initalized phitemp *',pi,0)
c
      do ii=1,nquad
      do ll=1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         phil = 2*pi*ll/nquadc
         phil = 2*pi*(ll-1)/nquadc
         cphil = dcos(phil)
         ctheta = -salphai*cphil
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,ll,ii) = phitemp(m,ll,ii) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
      enddo
ccc      call prin2(' computed phitemp *',pi,0)
c
      do ii = 1,nquad
      do jj = 1,nquadm
      do ll = 1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         betaj = 2*pi*jj/nquadm
         betaj = 2*pi*(jj-1)/nquadm
         cbetaj = dcos(betaj)
         sbetaj = dsin(betaj)
ccc         phil = 2*pi*ll/nquadc
         phil = 2*pi*(ll-1)/nquadc
         cphil = dcos(phil)
         sphil = dsin(phil)
         x = cbetaj*calphai*cphil-sbetaj*sphil
         y = (sbetaj*calphai*cphil+cbetaj*sphil)
ccc         zg = -salphai*cphil
ccc         if ( ((ii.eq.1).and.(jj.eq.1)).and.(ll.eq.1)) then
ccc         call prin2('111 x is *', x,1)
ccc         call prin2('111 y is *', y,1)
ccc         call prin2('111 zg is *', zg,1)
ccc         endif
         phi = datan2(y,x)
c         
         phival(ll,jj,ii) = 0.0d0
         ephik = cdexp(imag*phi)
         ephi = ephik**(-nterms)
         do m = -nterms,nterms
            phival(ll,jj,ii) = phival(ll,jj,ii) + phitemp(m,ll,ii)*ephi
            ephi = ephi*ephik
         enddo
ccc         if ( ((ii.eq.1).and.(jj.eq.1)).and.(ll.eq.1)) then
ccc         call prin2('111 phival is *', phival(ll,jj,ii),2)
ccc         endif
      enddo
      enddo
      enddo
ccc      call prin2(' computed phival *',pi,0)
      return
      end
C
cccccccccccccccccccccc
c
      subroutine mkgreatcircles(rk,ngridt,xnodesth,ngridp,xnodesph,
     1           ngridc,xg,yg,zg)
c
c     create Cartesian coordinates for equatorial circles on sphere of
c     radius rk with z-axis rotated to location 
c     (xnodesth(i),xnodesph(j)).
c
c     INPUT:
c
c     rk         radius of sphere
c     ngridt     number of discretization nodes in theta in lab frame
c     xnodesth   discretization nodes in theta in lab frame
c     ngridp     number of discretization nodes in phi in lab frame
c     xnodesph   discretization nodes in phi in lab frame
c     ngridc     number of nodes on great circle
c
c     OUTPUT:
c
c     xg(l,j,i)  xcoord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     yg(l,j,i)  y-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     zg(l,j,i)  z-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c
c     NOTE the ordering of output grids.
c
      implicit none
      integer i,j,l,ngridc,ngridt,ngridp
      real *8 pi,rk,calphai,salphai
      real *8 betaj,cbetaj,sbetaj,phil,cphil,sphil
      real *8 xnodesth(ngridt)
      real *8 xnodesph(ngridp)
      real *8 xg(ngridc,ngridp,ngridt)
      real *8 yg(ngridc,ngridp,ngridt)
      real *8 zg(ngridc,ngridp,ngridt)
c
      pi = 4.0d0*datan(1.0d0)
c
      do i = 1,ngridt
         calphai = xnodesth(i)
         salphai = dsqrt(1 - calphai**2)
         do j = 1,ngridp
            betaj = xnodesph(j)
            cbetaj = dcos(betaj)
            sbetaj = dsin(betaj)
            do l = 1,ngridc
ccc               phil = 2*pi*l/ngridc
               phil = 2*pi*(l-1)/ngridc
               cphil = dcos(phil)
               sphil = dsin(phil)
               xg(l,j,i) = rk*(cbetaj*calphai*cphil-sbetaj*sphil)
               yg(l,j,i) = rk*(sbetaj*calphai*cphil+cbetaj*sphil)
               zg(l,j,i) = rk*(-salphai*cphil)
            enddo
         enddo
      enddo
      return
      end
cc
cc
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkonegreatcircle(rk,cthetause,phiuse,rgamma,
     1           ngridc,xg,yg,zg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     create Cartesian coordinates a single equatorial circle on sphere of
c     radius rk with z-axis rotated to location 
c     (cthetause, phiuse).
c
c     INPUT:
c
c     rk         radius of sphere
c     ngridt     number of discretization nodes in theta in lab frame
c     xnodesth   discretization nodes in theta in lab frame
c     ngridp     number of discretization nodes in phi in lab frame
c     xnodesph   discretization nodes in phi in lab frame
c     ngridc     number of nodes on great circle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     OUTPUT:
c
c     xg(l)  xcoord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     yg(l)  y-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c     zg(l)  z-coord of lth node in equatorial circle for
c                z-axis rotated to (xnodesth(i),xnodesph(j)).
c
c     NOTE the ordering of output grids.
c
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
         phil = 2*pi*(l-1)/ngridc
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
C
C***********************************************************************
      subroutine shevalsphereonecirc(local,nterms,lmp,cthetause,phiuse,
     1           rgamma,nquadc,phival)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion at
C     a single equatorial circle corresponding to z-axis rotation 
C     to (theta,phi), where ctheta = cos(theta).
C     nquadc is the number of points desired on each great circle.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     local     : coefficients of spherical harmonic exp.
C     nterms    : number of terms in the orig. expansion
C     lmp       : dimension parameter for local
C     cthetause : cos theta for y-axis rotation.
C     phiuse    : angle for extrinsic z-axis rotation.
C     rgamma    : rotation of reference image
C     nquadc   : number of output points 
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : phival(l) corresponds to lth node on equatorial circle
C                for normal corresponding to (theta,phi) rotated
C                by gamma
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8, allocatable ::  ynm(:,:)
      complex *16 local(0:lmp,-lmp:lmp)
      complex *16 phival(nquadc)
      complex *16, allocatable :: phitemp(:,:)
ccc      complex *16, allocatable :: phivaltemp(:)
ccc      complex *16, allocatable :: wsave(:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
      allocate(phitemp(-nterms:nterms,nquadc))
ccc      allocate(wsave(4*nquadc+16))
ccc      allocate(phivaltemp(nquadc))
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
ccc      call prin2(' initalized phitemp *',pi,0)
c
      do ll=1,nquadc
         calphai = cthetause
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         call prin2(' calphai is *',calphai,1)
         phil = 2*pi*(ll-1)/nquadc
         phil = phil + rgamma
         cphil = dcos(phil)
         ctheta = -salphai*cphil
         call ylgndr(nterms,ctheta,ynm)
         do m=-nterms,nterms
            mabs = abs(m)
            do n=mabs,nterms
               phitemp(m,ll) = phitemp(m,ll) +
     1                local(n,m)*ynm(n,mabs)
            enddo
         enddo
      enddo
ccc      call prin2(' computed phitemp *',phitemp(0,1),2)
ccc      call prin2(' computed phitemp *',phitemp,2*(2*nterms+1)*nquadc)
c
      do ll = 1,nquadc
c         call prinf(' ll= *',ll,1)
         calphai = cthetause
         salphai = dsqrt(1.0d0 - calphai**2)
         betaj = phiuse
         cbetaj = dcos(betaj)
         sbetaj = dsin(betaj)
         phil = 2*pi*(ll-1)/nquadc
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
c         call prin2(' computed phival(ll) = *',phival(ll),2)
         enddo
      enddo
ccc      do ll = 1,nquadc
ccc         phivaltemp(ll) = phival(ll)
ccc      enddo
ccc      
ccc      call dcffti(nquadc,wsave)
ccc      call dcfftf(nquadc,phivaltemp,wsave)
ccc      do ll = 1,nquadc
ccc         phivaltemp(ll) = phivaltemp(ll)/nquadc
ccc      enddo
ccc      call prin2(' phivaltemp *',phivaltemp,2*nquadc)
ccc      call prin2(' rgamma *',rgamma,1)
ccc      call prinf(' nquadc *',nquadc,1)
ccc      do ll = 1,nquadc
ccc         phival(ll) = 0.0d0
ccc         phi = 2*pi*(ll-1)/nquadc + rgamma
ccc         do ii = 1,nquadc/2
ccc            phival(ll) = phival(ll) + 
ccc     1                   cdexp(phi*(ii-1)*imag)*phivaltemp(ii)
ccc         enddo
ccc         do ii = 1,nquadc/2-1
ccc            iuse = -ii
ccc            phival(ll) = phival(ll) + 
ccc     1                  cdexp(phi*iuse*imag)*phivaltemp(nquadc-ii+1)
ccc         enddo
ccc      enddo
ccc      call prin2(' phival *',phival,2*nquadc)
      return
      end
C
C
c-----------------------------------------------------------------------
c
c      sheval3d: computes potential due to SH expansion 
c                f_j =  A * {M_nm}  f_j is value at point j given
c                       by pts (cos theta_j,phi_j)
c
c      shproj3d: adjoint of sheval3d wth diagonal prefactor
c                 a collection of charges  
c                 N_nm = A^H D f_j
c
c**********************************************************************
      subroutine sheval3d(mpole,nterms,ctheta,phi,pot)
c**********************************************************************
c
c     This subroutine evaluates the spherical harmonic expansion
c
c     pot =  sum sum  mpole(n,m) Y_nm(theta,phi) 
c             n   m
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     mpole  :    multipole expansion in 2d matrix format
c     nterms :    order of the multipole expansion
c     ctheta,phi  :    target location
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ier
      integer lpp,ipp,lephi,lused
      real *8 ctheta,phi
      real *8, allocatable :: w(:)
      complex *16 pot
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16, allocatable :: wephi(:)
c
      ier=0
c
c     Carve up workspace for Ynm 
c
      lpp=(nterms+1)**2+5
      ipp=1
      lused = ipp+lpp
      allocate(w(lused))
ccc      iephi = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
ccc      lephi=2*(2*nterms+3)+5 
      lephi=(2*nterms+3)+5 
c
      lused=1+lephi
      allocate(wephi(lused))
c
      call sheval3d0(mpole,nterms,ctheta,phi,pot,w(ipp),wephi)
c
      return
      end
c
c
c
c**********************************************************************
      subroutine sheval3d0(mpole,nterms,ctheta,phi,pot,ynm,ephi)
c**********************************************************************
c
c     See sheval3d for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer nterms
      integer i,n,m
      real *8 ctheta,phi
      real *8 ynm(0:nterms,0:nterms)
      complex *16 pot,ephi1
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
c
      real *8 done,cphi,sphi,stheta
      complex *16 eye
      complex *16 ztmp2
c
      data eye/(0.0d0,1.0d0)/
c
      done=1.0d0
c
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     get the associated Legendre functions
c     and scale by 1/sqrt(2l+1)
c
      call ylgndr(nterms,ctheta,ynm)
c
      pot=mpole(0,0)
c
      do n=1,nterms
         pot=pot+mpole(n,0)*ynm(n,0)
         do m=1,n
            ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
            pot=pot+ynm(n,m)*ztmp2
         enddo
      enddo
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine shproj3d(cthetas,phis,weight,charge,ns,nterms,mpole)
C***********************************************************************
C
C 
C     Compute adjoint of sheval3d: A^H D with D defined by weights.
c-----------------------------------------------------------------------
C     INPUT:
c
C     cthetas(ns)   : source cthetas
C     phis(ns)      : source phis
C     weight(ns)    : diag weights
C     charge(ns)    : function val at ctheta,phi
C     ns              : number of sources
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     mpole           : coeffs of the multipole expansion
c-----------------------------------------------------------------------
      implicit none
      integer nterms,ns,i,l,m, ier
      real *8 cthetas(ns)
      real *8 phis(ns)
      real *8 weight(ns)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,charge(ns)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      ier = 0
      do i = 1, ns
         call shproj3d1(cthetas(i),phis(i),weight(i),charge(i),
     1        nterms,mpole)
      enddo
c
      return
      end
C
c**********************************************************************
      subroutine shproj3d1(ctheta,phi,weight,charge,
     1          nterms,mpole)
c**********************************************************************
c
c     This is the memory management routine. Work is done in the
c     secondary call to l3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     ctheta,phi  : coordinates 
c     weight  : diag weight 
c     charge  : complex charge strength
c     nterms  : order of the expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c     mpole   : coeffs of the expansion
c-----------------------------------------------------------------------
      implicit none
      integer ier,nterms
      integer ipp,lpp,iephi,lephi,lused
      real *8 weight,ctheta,phi
      real *8, allocatable :: w(:)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
      complex *16, allocatable :: wephi(:)
c
c     compute work space components:
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
      lused=ipp+lpp
      allocate(w(lused))
      
c
      iephi=1
      lephi=(2*nterms+1)+7
c
      lused=iephi+lephi
      allocate(wephi(lused))
ccc      call prinf(' in formmp lused is *',lused,1)
c
      call shproj3d0(ctheta,phi,weight,charge,nterms,
     1          mpole,w(ipp),wephi(iephi))
      return
      end
c
c
c
c**********************************************************************
      subroutine shproj3d0(ctheta,phi,weight,charge,nterms,
     1          mpole,pp,ephi)
c**********************************************************************
c
c     See l3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit none
      integer ier,nterms,i,n,m
      real *8 weight,ctheta,phi
      real *8 pp(0:nterms,0:nterms)
      real *8 stheta,cphi,sphi
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 charge
      complex *16 ephi(-nterms:nterms),ephi1
      complex *16  ztmp
c
      ier=0
c
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     get the associated Legendre functions:
c
      call ylgndr(nterms,ctheta,pp)
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     Compute contribution to mpole coefficients.
c
      mpole(0,0)= mpole(0,0) + charge
      do n=1,nterms
         ztmp=pp(n,0)*charge
         mpole(n,0)= mpole(n,0) + ztmp
         do m=1,n
            ztmp=pp(n,m)*charge
            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
         enddo
      enddo
      return
      end
c
c
c
      subroutine mkshinterp(thetas,phis,nout,k,nquad,nquadm,values,
     1           icols)
      real *8 thetas(nout)
      real *8 phis(nout)
      real *8 values(nout,k*k)
      real *8, allocatable :: wtemp(:)
      integer icols(nout,k*k)
      integer, allocatable :: icoltemp(:) 
c
c     Create interpolation matrix in simple sparse storage 
c     format. Since we are using a kxk interpolation stencil,
c     we only need to know the column index and interpolation
c     weights.
c
c     INPUT:
c    
c     thetas,phis     :  coordinates of target points
c     nout            :  number of target points
c     k               :  interpolation order
c     nquad           :  number of nodes in theta
c     nquadm          :  number of nodes in phi
c
c     OUTPUT:
c    
c     values          :  interpolation weights from reg grid.
c     icols           :  unrolled indices of source points for interp.
c
c
c     For each target point, compute the next k*k interpolation
c     weights and increment values,icols.
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
      subroutine getinterpsparse(theta,phi,k,nquad,nquadm,icol,
     1           wtemp)
      implicit real *8 (a-h,o-z)
      real *8 pi,theta,phi
      real *8 wtemp(k*k)
      integer icol(k*k)
c
c     Create one row of interpolation matrix in simple sparse storage 
c     format. Since we are using a kxk interpolation stencil,
c     we only need to know the column index and interpolation
c     weights
c
c     INPUT:
c    
c     theta,phi  : coordinates of target point
c     k          : interpolation order
c     nquad      : number of equspaced nodes in theta
c     nquadm     : number of equspaced node in phi
c
c     OUTPUT:
c
c     icol       : source pt on sphere in unrolled vector format.
c                  Assuming phival is given in the 2d fortran format:
c                  phival(nquadm,nquad)
c     wtemp      : corresponding interpolation weights
c
c     f(x,y) = sum_{j1,j_2} f(j_1h_1,j_2h_2) l_{j_1,j_2}(x,y)
c
c     l_{j_1,j_2}(x,y)=\prod_{k_1\neq j_1} (x_1-k_1h_1)/(j_1h_1-k_1h_1)
c                      * \prod_{k_2\neq j_2} (x_2-k_2h_2)/(j_2h_2-k_2h_2)
c
c     Need to grab correctly wrapped (NP,NT) values for "source pt" 
c     weight is computed without regard for wrapping
c
c----------------------------------------------------------------------    
      pi = 4*datan(1.0d0)
c
c     get nearest grid coords (np0,nt0)
c
      np0 = 1 + nint(nquadm*phi/(2*pi))
ccc      nt0 = 1 + nint(nquad*(theta + pi/(2*nquad))/pi)
      nt0 = nint(nquad*(theta + pi/(2*nquad))/pi)
c
c     need to compute offset from base point before thinking about
c     wrapping interpolation points...
c
      ht = pi/nquad
      hp = 2*pi/nquadm
      dp = phi - (np0-1)*hp
      dth = theta - ((2*nt0-1)*pi)/(2*nquad) 
c
      if (np0 .gt. nquadm) np0 = np0- nquadm
      if (nt0 .gt. nquad) nt0 = nt0-1
c
      if ( (np0.lt.1).or.(np0.gt.nquadm)) then
c         call prinf(' np0 is *',np0,1)
         write(6,*) 'np0 ', np0, theta, phi
      endif
      if ( (nt0.lt.1).or.(nt0.gt.nquad)) then
c         call prinf(' nt0 is *',nt0,1)
         write(6,*) 'np0 ', np0, theta, phi
      endif

      

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
      subroutine getlagrangewt(dth,dp,k,i,j,ht,hp,val)
      implicit real *8 (a-h,o-z)
      integer k,i,j
      real *8 dth,dp,ht,hp,val
c
c     Compute interpolation weight at (dth,dp) from i,j node.
c
c     INPUT:
c     dth,dp are offsets from base point (assumed to be origin).
c     k is interpolation order
c     i,j is the particular grid point whose lagrange polynomial
c     is evaluated at (dp,dth).
c     ht,hp  : grid spacing in theta and phi, respectively.
c
c     OUPUT:
c     val is value of Lagrange interpolant assoc with i,j point.
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
      subroutine shinterp(nout,k,nquad,nquadm,values,
     1           icols,uval,uout)
      implicit real *8 (a-h,o-z)
      integer icols(nout,k*k)
      real *8 values(nout,k*k)
      complex *16 uout(nout)
      complex *16 uval(nquadm*nquad)
c
c     Apply interpolation matrix in sparse formt.
c
c     INPUT:
c     nout     : number of targets
c     k        : interpolation order
c     nquad    : number of nodes in theta
c     nquadm   : number of nodes in phi
c     values   : sparse format interpolation weights
c     icols    : sparse format col indices
c     uval     : values on regular grid
c
c     OUTPUT:
c     uout     : values at irregular nodes corresp. to 
c                rows of values/icols
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shinterp_adj(nout,k,nquad,nquadm,values,
     1           icols,uin,ugrid)
      implicit real *8 (a-h,o-z)
      integer icols(nout,k*k)
      real *8 values(nout,k*k)
      complex *16 uin(nout)
      complex *16 ugrid(nquadm*nquad)
c
c     Apply adjoint of interpolation matrix in sparse formt.
c
c     INPUT:
c     nout     : number of targets
c     k        : interpolation order
c     nquad    : number of nodes in theta
c     nquadm   : number of nodes in phi
c     values   : sparse format interpolation weights
c     icols    : sparse format col indices
c     uin      : input vector (dimensioned as values on irregular grid)
c
c     OUTPUT:
c     ugrid     : adjoint of interp matrix applied to uin
c   
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
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine unrolltoshxp(nterms,sol,sh1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine converts an "unrolled" spherical harmonic
c     expansion to a doubly indexed array where (i,j) refers to 
c     the Y_{i,j} coefficient.
c
c     INPUT:  
c             sol  (unrolled)
c     OUTPUT: 
c             sh1  (matrix format)
c
      implicit real *8 (a-h,o-z)
      complex *16 sol(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
c
      next = 1
      do i = 0,nterms
         do j = -i,i
            sh1(i,j) = sol(next)
	    next = next+1
         enddo
      enddo
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shxptounroll(nterms,sol,sh1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine converts a double indexed spherical harmonic
c     expansion to an "unrolled" one. 
c
c     INPUT: 
c              sh1  (matrix format)
c     OUTPUT:  
c              sol  (unrolled)
c
      implicit real *8 (a-h,o-z)
      complex *16 sol(1)
      complex *16 sh1(0:nterms,-nterms:nterms)
c
      next = 1
      do i = 0,nterms
         do j = -i,i
            sol(next) = sh1(i,j) 
	    next = next+1
         enddo
      enddo
      return
      end
