C***********************************************************************
      subroutine rebuild_full_model(cslices,imagesize,nimages,
     1           alphas,icstart,nlats,itypep,ngridc,nlow,ncur,
     2           lmod,isph_start,nterms_sph,oversamp,kord,eps,modelsph)
C***********************************************************************
C     Builds new model (in spherical harmonic basis) 
C     from slices (Fourier transforms of images) with image parameters 
C     given by alphas array.
C---------------------------------------------------------------------
C     INPUT:
C
C     cslices        Fourier transforms of images on polar grids
C     imagesize      dimension of images on quasiunform polar grid
C     nimages        number of images
c     alphas(3,*)    image parameters: 
c                          alphas(1,*) is first Euler angle (alpha)
c                          alphas(2,*) is second Euler angle (beta)
c                          alphas(3,*) is third Euler angle  (gamma)
c
c         alphas(5,*) TO BE IMPLEMENTED
c                          alphas(4,*) is x-offset of center
c                          alphas(5,*) is y-offset of center
c                          alphas(6,*) is nromalization factor
c
c     icstart()      indexing array for points on successive circles
c     nlats()        number of quarature nodes in theta on sphere
c                    defined by index ncur.
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
c     ngridc()       number of output points on successive circles in templates
c     nlow           index of lowest frequency sphere under consideration
c     ncur           index of highest frequency sphere under consideration
c     lmod           length of packed spherical harmonic model up to sphere
c                    ncur
c     isph_start()   array of length ngridr indicating where in the 
c                    modhat_sph vector the coefficients for the 
c                    corresponding sphere begin
c     nterms_sph()   array of length ngridr defining the orders of the
c                    spherical harmonic expansions on successive spheres.
c     oversamp       oversampling parameter for least squares solver
c     kord           interpolation order for least squares solver
c     eps            tolerance for least squares solver
C
C     OUTPUT: 
c
C     modelsph(:)    solution to least squares problem expressed in 
c                    spherical harmonics in successive spheres
c
C***********************************************************************
      implicit none
      integer itypep,nlow,ncur,lmod,imagesize,nimages,kord
      integer numit,is,jj,jstart,nout,noutmax,nquad
      integer nsphmax,nsphtot,nterms,ntmax
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      integer isph_start(ncur),nterms_sph(ncur)
      real *8 oversamp,eps
      real *8 alphas(3,nimages)
      real *8, allocatable :: thetas(:)
      real *8, allocatable :: phis(:)
      complex *16 modelsph(0:lmod-1)
      complex *16 cslices(imagesize,nimages)
      complex *16, allocatable :: localp(:)
      complex *16, allocatable :: phiout(:)
c
      nterms = nterms_sph(ncur)
      nsphmax = (nterms+1)**2
      allocate(localp(0:nsphmax))
c
      ntmax = nint(oversamp*nlats(ncur))
      write(6,*) ' oversamp = ',oversamp
      write(6,*) ' ntmax = ',ntmax
c
      noutmax = nimages*ngridc(ncur)
      allocate(phiout(noutmax))
      allocate(thetas(noutmax))
      allocate(phis(noutmax))
c
c
c
      do is = nlow,ncur
         write(6,*) ' is = ',is
         call get_lsqdata(cslices,imagesize,nimages,icstart(is),
     1           ngridc(is),alphas,phiout,thetas,phis,nout)
         nterms = nterms_sph(is)
c
c        set up parameters for tensor product grid on sphere
c        to be used in least squares solve (for interpolation
c        to arbitrary points)
c
         nquad = nint(oversamp*nlats(is))
         nsphtot = (nterms+1)**2
         numit = 1000
         call lsqsolvespharm(phiout,thetas,phis,nout,localp,nterms,
     1           nquad,kord,eps,numit)
         jstart = isph_start(is)
         do jj = 0,nsphtot-1
            modelsph(jstart+jj) = localp(jj)
         enddo
      enddo
      return
      end
c
c
c
c
c
c
      subroutine get_lsqdata(ximagesk,imagesize,nimages,icstart,
     1           ngridc,alphas,phiout,thetas,phis,nout)
c
      implicit none
      integer imagesize,nimages,icstart,ngridc
      integer nout,ii,iincr,ll
      complex *16 ximagesk(imagesize,nimages)
      complex *16 phiout(ngridc*nimages)
      real *8 alphas(3,nimages)
      real *8 thetas(ngridc*nimages)
      real *8 phis(ngridc*nimages)
      real *8 pi,ctheta,cthetause
      real *8 phi,phiuse,rgammause,rho,rkx,rky,rkz,rone
      real *8, allocatable :: xc(:)
      real *8, allocatable :: yc(:)
      real *8, allocatable :: zc(:)
c
      allocate(xc(ngridc))
      allocate(yc(ngridc))
      allocate(zc(ngridc))
c
      pi=4.0d0*datan(1.0d0)
      iincr = 0
      do ii = 1,nimages
         cthetause = dcos(alphas(1,ii))
         phiuse = alphas(2,ii)
         rgammause = alphas(3,ii)
         rone = 1.0d0
         call mkonegreatcircle(rone,cthetause,phiuse,rgammause,
     1        ngridc,xc,yc,zc)
c
         do ll = 1,ngridc
            rkx = xc(ll)
            rky = yc(ll)
            rkz = zc(ll)
            ctheta = rkz
            phi = datan2(rky,rkx)
            if (phi.lt.0.0d0) phi = phi + 2*pi
            iincr = iincr+1
            phis(iincr) = phi
            rho = dsqrt(rkx**2+ rky**2)
ccc            thetas(iincr) = dacos(ctheta)
            thetas(iincr) = datan2(rho,ctheta)
            phiout(iincr) = ximagesk(icstart+ll-1,ii)
         enddo
      enddo
      nout = ngridc*nimages
      return
      end
