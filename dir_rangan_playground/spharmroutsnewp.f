!> Doxygen comment: ;\n
!>    Fortran library for sph harm projection and eval on grids & circles ;\n
!>    ;\n
!>    Stuff from Leslie/Zydrunas FMMLIB3D, and localexp3d. ;\n
!>    pulled by Barnett 7/31/15 ;\n
!>    modified by Greengard 8/1/15 ;\n
!>    ;\n
!>    routines that end in 'p' use the packed representation of spherical harmonic coefficients ;\n
!>    stored linearly without gaps in the (n,m) order of ;\n
!>      (0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2), (3,-3), ... ;\n
!>    The total number of coefficients is (nterms+1)*(nterms+1) ;\n
!>    for an expansion up to order nterms (inclusive). ;\n
!>    Using 1 based indexing, the linear index for element (n, m) is ix = n*(n+1)+m+1 ;\n
!>    ;\n
!>    ;\n
!>********************************************************************** ;\n
!>    calling sequences changed from Barnett's spharmrouts.f  ;\n
!>********************************************************************** ;\n
!> Doxygen comment: ;\n
C     Fortran library for sph harm projection and eval on grids & circles
C
C     Stuff from Leslie/Zydrunas FMMLIB3D, and localexp3d.
C     pulled by Barnett 7/31/15
C     modified by Greengard 8/1/15
C
C     routines that end in 'p' use the packed representation of spherical harmonic coefficients
C     stored linearly without gaps in the (n,m) order of
C       (0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2), (3,-3), ...
C     The total number of coefficients is (nterms+1)*(nterms+1)
C     for an expansion up to order nterms (inclusive).
C     Using 1 based indexing, the linear index for element (n, m) is ix = n*(n+1)+m+1
C
C
C***********************************************************************
C     calling sequences changed from Barnett's spharmrouts.f 
C***********************************************************************
c
C***********************************************************************
      subroutine projloc3dp(nterms,nquad,nquadm,xnodes,wts,
     1           phival,localp)
C***********************************************************************
C     Usage:
C
C           compute spherical harmonic expansion on unit sphere
C           of function tabulated at nquad*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           nquad  = number of quadrature nodes in theta direction.
C           nquadm = number of quadrature nodes in phi direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           wts    = Gauss quadrature weights
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           NOTE:    We assume phi_i = (i-1)*2.0d0*pi/nquadm, as do the 
C                    routines in projection.f. However, we permit
C                    different numbers of nodes in theta and phi.
C***********************************************************************
C     OUTPUT:
C
C           localp = coefficients of s.h. expansion (packed)
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
Cccc        marray(jj,m) = sum*2.0d0*pi/nquadm
C           marray(jj,m) = sum/(2*nquadm)
C
C     The latter has incorporated the 1/(4*pi) normalization factor
C     into the azimuthal quadrature weight (2.0d0*pi/nquadm).
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nquad,nquadm
      integer l,m,jj,kk,ix
      real *8 wts(1),xnodes(1)
      real *8, allocatable ::  ynm(:,:)
      complex *16 phival(nquadm,nquad)
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 ephi,imag,emul,sum,zmul,emul1
      complex *16, allocatable :: marray(:,:)
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
      allocate(ynm(0:nterms,0:nterms))
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
         do jj=1,nquad
            sum = 0
            ephi = 1.0d0
            do kk = 1,nquadm
ccc               ephi = cdexp(imag*m*(kk-1)*2.0d0*pi/nquadm)
               sum = sum + phival(kk,jj)*dconjg(ephi)
               ephi = ephi*emul
            enddo
ccc         marray(jj,m) = sum*2.0d0*pi/nquadm
            marray(jj,m) = sum/(2*nquadm)
         enddo
c     I don't understand what value emul1 has - it's never set! (ahb):
         emul = emul*emul1
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
               ix = l*(l+1) + m + 1
               localp(ix) = localp(ix) + zmul*ynm(l,abs(m))
            enddo
         enddo
      enddo
      return
      end



C***********************************************************************
      subroutine shevalspherep(localp,phival,
     1           nterms,nquad,nquadm,xnodes)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion on an 
C     nquad x nquadm grid on the unit sphere. 
C
C---------------------------------------------------------------------
C     INPUT:
C
C     localp   : coefficients of spherical harmonic exp. (packed)
C     nterms   : number of terms in the orig. expansion
C     nquad    : number of quadrature nodes in theta
C     nquadm   : number of quadrature nodes in phi
C     xnodes   : Legendre nodes in theta (x_j = cos theta_j).
C---------------------------------------------------------------------
C     OUTPUT:
C
C     phival   : function value on tensor product
C                mesh on target sphere. phi is the fast variable, theta slow
C***********************************************************************
C     converted to OPENMP, Barnett 6/22/16
C
      implicit real *8 (a-h,o-z)
      integer nterms
      real *8 xnodes(1)
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


cc
C***********************************************************************
      subroutine shevalspherecircsp(localp,phival,
     1           nterms,nquad,nquadm,xnodes,nquadc)
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
C     localp   : coefficients of spherical harmonic exp. (packed)
C     nterms   : number of terms in the orig. expansion
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
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phival(nquadc,nquadm,nquad)
      complex *16, allocatable :: phitemp(:,:,:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m,ix
      data imag/(0.0d0,1.0d0)/
C
      pi = 4.0d0*datan(1.0d0)
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
c$OMP PARALLEL PRIVATE(ynm)
c     each thread gets own allocation of ynm
      allocate(ynm(0:nterms,0:nterms))
C$OMP DO PRIVATE(ll,calphai,salphai,phil,cphil,ctheta,m,mabs,ix,n)
      do ii=1,nquad
      do ll=1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         phil = 2.0d0*pi*ll/nquadc
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

C$OMP PARALLEL DO PRIVATE(jj,ll,calphai,salphai,betaj,cbetaj,sbetaj,phil,cphil,&
C$OMP& sphil,x,y,phi,ephik,ephi,m)
      do ii = 1,nquad
      do jj = 1,nquadm
      do ll = 1,nquadc
         calphai = xnodes(ii)
         salphai = dsqrt(1.0d0 - calphai**2)
ccc         betaj = 2.0d0*pi*jj/nquadm
         betaj = 2.0d0*pi*(jj-1)/nquadm
         cbetaj = dcos(betaj)
         sbetaj = dsin(betaj)
ccc         phil = 2.0d0*pi*ll/nquadc
         phil = 2.0d0*pi*(ll-1)/nquadc
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
C$OMP END PARALLEL DO

ccc      call prin2(' computed phival *',pi,0)
      return
      end


C***********************************************************************
      subroutine shevalsphereonecircp(localp,nterms,cthetause,phiuse,
     1           rgamma,nquadc,phival)
C***********************************************************************
C
C     This subroutine evaluates a spherical harmonic expansion at
C     a single equatorial circle corresponding to z-axis rotation 
C     to (theta,phi), where ctheta = dcos(theta).
C     nquadc is the number of points desired on each great circle.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     localp    : coefficients of spherical harmonic exp. (packed)
C     nterms    : number of terms in the orig. expansion
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
      complex *16 localp((nterms+1)*(nterms+1))
      complex *16 phival(nquadc)
      complex *16, allocatable :: phitemp(:,:)
ccc      complex *16, allocatable :: phivaltemp(:)
ccc      complex *16, allocatable :: wsave(:)
      complex *16 imag
      complex *16 ephi,ephik
C
      integer m,ix
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
ccc         phi = 2.0d0*pi*(ll-1)/nquadc + rgamma
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

C***********************************************************************
      subroutine shevalspherep_adj(nterms,nquad,nquadm,xnodes,
     1           phival,localp)
C***********************************************************************
C     Usage:
C
C           compute adjoint of spherical harmonic expansion on unit 
C           sphere of function tabulated at nquadm*nquad grid points.
C---------------------------------------------------------------------
C     INPUT:
C
C           nterms = order of spherical harmonic expansion
C           nquad  = number of quadrature nodes in theta direction.
C           nquadm = number of quadrature nodes in phi direction.
C           xnodes = Gauss-Legendre nodes x_j = cos theta_j
C           phival = tabulated function
C                    phival(i,j) = phi(sin theta_j cos phi_i,
C                                      sin theta_j sin phi_i,
C                                      cos theta_j).
C
C           NOTE:    We assume phi_i = (i-1)*2.0d0*pi/nquadm, as do the 
C                    routines in projection.f. However, we permit
C                    different numbers of nodes in theta and phi.
C***********************************************************************
C     OUTPUT:
C
C           localp = cadjoint of sheval applied to phival
C           This is different from projloc3dp since it omits the 
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
      real *8 xnodes(1)
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
      subroutine quadscale99(nquad,nquadm,wts,phival)
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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sph_gridconfig199(ngridr, sph_size, sph_nterms, 
     1     sph_sizes,sph_startix)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Configure the grid of spherical harmonics expansions by radius
c
c     Inputs:
c     ngridr: number of radial grid points
c
c     Outputs:
c     sph_size:            total size of the array required to store the spherical harmonics expansions
c     sph_sizes(ngridr):   individual sizes of spherical hardmonics expansion by radius
c     sph_nterms(ngridr):  number of terms of the expansion by radius (nterms)
c     sph_startix(ngridr): starting index for an expansion by radius
c
      implicit none
      integer ngridr
      integer sph_size, sph_sizes(ngridr), sph_startix(ngridr)
      integer sph_nterms(ngridr)
c
      integer i, n
c
      sph_size = 0
      do i = 1,ngridr
c        hard code for all expansions
c         n = 30
c     size of expansions varies with i
c         n = max(5,i)
         n = i+2
         sph_nterms(i) = n
         sph_sizes(i) = (n+1)*(n+1)
         sph_startix(i) = sph_size

         sph_size = sph_size + (n+1)*(n+1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sph_gridconfig99(ngridr, sph_size, sph_nterms, 
     1     sph_sizes,sph_startix, xnodesr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Configure the grid of spherical harmonics expansions by radius
c
c     Inputs:
c     ngridr: number of radial grid points
c
c     Outputs:
c     sph_size:            total size of the array required to store the spherical harmonics expansions
c     sph_sizes(ngridr):   individual sizes of spherical hardmonics expansion by radius
c     sph_nterms(ngridr):  number of terms of the expansion by radius (nterms)
c     sph_startix(ngridr): starting index for an expansion by radius
c
      implicit none
      integer ngridr
      integer sph_size, sph_sizes(ngridr), sph_startix(ngridr)
      integer sph_nterms(ngridr)
      real *8 xnodesr(ngridr)
c
      integer i, n
c
      sph_size = 0
      do i = 1,ngridr
c        hard code for all expansions
c         n = 30
c     size of expansions varies with i
         n = floor(xnodesr(i))+2
c         write(6,*) i, xnodesr(i),n
         sph_nterms(i) = n
         sph_sizes(i) = (n+1)*(n+1)
         sph_startix(i) = sph_size

         sph_size = sph_size + (n+1)*(n+1)
      enddo

      return
      end






