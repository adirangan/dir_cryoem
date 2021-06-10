ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkfouriergrid(rkmax,ngridr,ityper,
     1                    nlats,itypep,ntot,numonsphere,
     2                    kgridx,kgridy,kgridz,wts)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ? add nphis(j,i) = number of azimuthal pts on jth latitude on
c                        ith sphere 
C     ? output theta, phi, r grid 
c
c     utility: given rk -> grid on that sphere, thetas/phis
c              
c
c     PURPOSE: Create spherical grid out to radius rkmax.
c     This routine can create either a tensor product
c     grid, or a discretization which uses different meshes on each
c     successive sphere, properly sampling in theta or phi or both.
c     
c     INPUT:
c
c     rkmax          outer radius of sphere (r in [0,rkmax])
c     ngridr         number of nodes in r
c     ityper         quadrature scheme in radial direction
c                       1 = Gaussian
c                       0 = uniform grid from 1..ngridr
c                           rad(i) = rkmax*i/ngridr
c     nlats()        number of latitude points used on each sphere
c                    permits coarser sampling near origin 
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is defined in getgridph and typically
c                    set to 2.
c     ntot           total number of points in grid
c
c
c     OUTPUT:
c
c     numonsphere()  number of points on ith sphere
c                    sum_i numonsphere(i) should equal ntot (input).
c                    
c     kgridx()       1st Cartesian coordinate of k-space nodes
c     kgridy()       2nd Cartesian coordinate of k-space nodes
c     kgridz()       3rd Cartesian coordinate of k-space nodes
c     wts()          quadrature weights for k-space nodes
c
c----------------------------------------------------------------------
      implicit none
      integer ngridr,ngridp,ngridt,i
      integer ityper,itypep,ntot,nnn
      integer ifprint,jj,kk,ll,ngridtmax
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      real *8 rkmax,u,v,pi,phi,phstep,rr
      real *8 ctheta,stheta
      integer, allocatable ::  ngridps(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
      real *8 kgridx(ntot)
      real *8 kgridy(ntot)
      real *8 kgridz(ntot)
      real *8 wts(ntot)
c
      pi = 4.0d0*datan(1.0d0)
c
c     get max theta array size and allocate local arrays
c
      ngridtmax = nlats(1)
      do i = 2,ngridr
         ngridtmax = max(nlats(i),ngridtmax)
      enddo
c
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
c
      allocate(xnodesth(ngridtmax))
      allocate(sthetas(ngridtmax))
      allocate(wtsth(ngridtmax))
      allocate(ngridps(ngridtmax))
      allocate(phsteps(ngridtmax))
c
      call getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
c      
      ifprint = 0
      if (ifprint.eq.1) call prin2(' xnodesr is *',xnodesr,ngridr)
      if (ifprint.eq.1) call prin2(' wtsr is *',wtsr,ngridr)
      if (ifprint.eq.1) write(6,*) ' ngridtmax is ',ngridtmax
c
      nnn = 0
      do jj = 1,ngridr
        rr = xnodesr(jj)
        ngridt = nlats(jj)
        call getspheregrid(ngridt,itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,numonsphere(jj))
        do kk = 1,ngridt
           ctheta = xnodesth(kk)
           stheta = sthetas(kk)
           phstep = phsteps(kk)
c
           do ll = 1,ngridps(kk)
              phi = (ll-1)*phstep
              nnn = nnn+1
              kgridx(nnn) = rr*stheta*dcos(phi)
              kgridy(nnn) = rr*stheta*dsin(phi)
              kgridz(nnn) = rr*ctheta
              wts(nnn) = wtsr(jj)*wtsth(kk)*phstep
           enddo
        enddo
      enddo
      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create radial grid out to rkmax.
c     
c     INPUT:
c
c     rkmax          outer radius of sphere (r in [0,rkmax])
c     ngridr         number of nodes in r
c     ityper         quadrature scheme in radial direction
c                       1 = Gaussian
c                       0 = uniform grid from 1..ngridr
c                           rad(i) = rkmax*i/ngridr
c
c     OUTPUT:
c
c     xnodesr        nodes in radial direction
c     wtsr           quadrature weights in radial direction
c
c----------------------------------------------------------------------
      implicit none
      integer i,ngridr,ityper
      real *8 u,v,rkmax
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8, allocatable :: b(:)
c
      if (ityper.eq.1) then
         allocate(b(ngridr))
c$$$         call gaussq(0, ngridr, 0, 0, 0, 0, b, xnodesr, wtsr)
         u = rkmax/2.0d0
         v = rkmax/2.0d0
         do i = 1,ngridr
            xnodesr(i) = u*xnodesr(i) + v
            wtsr(i) = xnodesr(i)*xnodesr(i)*wtsr(i)*rkmax/2.0d0
         enddo
      else
         do i = 1,ngridr
            xnodesr(i) = rkmax*(i)/ngridr
            wtsr(i) =  xnodesr(i)*xnodesr(i)*rkmax/ngridr
         enddo
      endif
c
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getgridth(ngridt,t,xnodesth,sthetas,wtsth)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create theta grid on sphere 
c     
c     INPUT:
c
c     ngridt         number of latitudes (equispaced in theta)
c     t              workspace of length ngridt
c
c     OUTPUT:
c
c     xnodesth       nodes in [-1,1] corresponding to cos(theta_j)
c                    where theta_j are equispaced at classical 
c                    Chebyshev nodes
c     sthetas        corresponding sin(theta) values
c     wtsth          quadrature weights in theta direction
c
c----------------------------------------------------------------------
      implicit none
      integer i,ngridt,itype
      real *8 alpha,beta
      real *8 u,v,rkmax
      real *8 t(ngridt)
      real *8 xnodesth(ngridt)
      real *8 sthetas(ngridt)
      real *8 wtsth(ngridt)
c
c    Chebyshev quadrature in theta (incorporates sin(theta) weight)
c    using equispaced (classical) nodes
c
      itype = 1
      alpha = 0.0d0
      beta = 0.0d0
      call chebexps(itype,ngridt,t,alpha,beta,wtsth)
      do i = 1,ngridt
         xnodesth(i) = t(ngridt+1-i)
         sthetas(i) = dsqrt(1.0d0 - xnodesth(i)**2)
      enddo
c
      return
      end
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getgridph(itypep,ngridt,stheta,ngridp,phstep)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create phi grid on sphere at a specific latitude
c              theta.
c     
c     INPUT:
c
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    Typically, phi_over = 2 is sufficient for resolution,
c                    since we are integrating on a circle of radius
c                    [0,2 pi]*gridr(i)*sin(theta) while the
c                    nlats(i) points are integrating on the interval
c                    [0, pi]*gridr(i)
c
c     ngridt         number of grid points used in theta direction
c     stheta         value of sin(theta) on given latitude
c
c     OUTPUT:
c
c     ngridp       number of nodes in phi
c     phstep       quadrature weightin phi direction
c
c----------------------------------------------------------------------
c     IMPORTANT: phi_over is set here and all codes call this routine
c                to find ngridp, phstep... 
c----------------------------------------------------------------------
      implicit none
      integer itypep,ngridt,ngridp
      real *8 stheta,pi,phstep,phi_over
c
      pi = 4.0d0*datan(1.0d0)
      phi_over = 2.0d0
      if (itypep.eq.0) then
         ngridp = nint(ngridt*phi_over)
      else 
         ngridp = max(nint(ngridt*phi_over*stheta),6)
      endif
      phstep = pi*2.0d0/ngridp
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getspheregrid(ngridt,itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,numonsphere)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create phi grid on sphere for each latitude
c     
c     INPUT:
c
c     ngridt         number of grid points used in theta direction
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
c
c     OUTPUT:
c
c     xnodesth()     nodes in theta
c     sthetas()      corresponding sin(theta) values
c     wtsth()        quadrature weights in theta direction
c     ngridps()      number of nodes in phi for each latitude
c     phsteps()      quadrature weightin phi direction
c     numonsphere    total number of points on sphere
c
c----------------------------------------------------------------------
      implicit none
      integer i,ngridt,itypep,numonsphere,kk
      integer ngridps(ngridt)
      real *8 xnodesth(ngridt)
      real *8 sthetas(ngridt)
      real *8 wtsth(ngridt)
      real *8 phsteps(ngridt),pi
      real *8, allocatable :: t(:)
c
      allocate(t(ngridt))
c
      call getgridth(ngridt,t,xnodesth,sthetas,wtsth)
      numonsphere = 0
      do kk = 1,ngridt
           call getgridph(itypep,ngridt,sthetas(kk),
     1                    ngridps(kk),phsteps(kk))
           numonsphere = numonsphere+ngridps(kk)
      enddo
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_ntot(ngridr,nlats,itypep,ntot)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Compute total number of points used in spherical grid.
c     
c     INPUT:
c
c     ngridr         number of grid points used in radial direction
c     nlats()        number of latitude points used on each sphere
c                    permits coarser sampling near origin 
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is set in getgridph (typically = 2). 
c
c     OUTPUT:
c
c     ntot           total number of points in grid
c----------------------------------------------------------------------
      implicit none
      integer ntot,ngridr,ntmax,nlats(ngridr)
      integer itypep,numonsphere,jj
      real *8 pi
      integer, allocatable ::  ngridps(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: wtsth(:)
c
      ntmax = nlats(1)
      do jj = 2,ngridr
         ntmax = max(nlats(jj),ntmax)
      enddo
c
      allocate(ngridps(ntmax))
      allocate(phsteps(ntmax))
      allocate(xnodesth(ntmax))
      allocate(sthetas(ntmax))
      allocate(wtsth(ntmax))

      ntot = 0
      do jj = 1,ngridr
         call getspheregrid(nlats(jj),itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,numonsphere)
         ntot = ntot + numonsphere
      enddo
      return
      end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fun_to_kspacegrid(x,y,z,ff,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat,ffhatnorm,ier)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: compute transform from density in physical space to 
c              spherical grid in k-space
c     
c     INPUT:
c
c     x,y,z          3D arrays of physical space coordinates
c     ff             3D array  of density in physical space
c     h              grid spacing in x,y,z
c     ngrid          number of grid pts in each linear dimension
c     eps            NUFFT tolerance
c     rmax           max k-space radius
c     ngridr         number of grid points used in radial direction
c     ityper         quadrature scheme in radial direction
c                       1 = Gaussian
c                       0 = uniform grid from 1..ngridr
c                           rad(i) = rkmax*i/ngridr
c     nlats()        number of latitude points used on each sphere
c                    permits coarser sampling near origin 
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c     ntot           total number of points in k-space grid
c     numonsphere    number of points on successive spheres
c
c     OUTPUT:
c
c     kgridx         first coordinates of target points in k-space 
c     kgridy         second coordinates of target points in k-space 
c     kgridz         third coordinates of target points in k-space 
c     wts            quad weights in radial direction
c     ffhat          Fourier transform values on spherical grid
c     ffhatnorm      L2 norm of computed Fourier transform 
c     ier            error return code from NUFFT
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer ngrid,ngridr,ityper,itypep,ntot,ier
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      real *8 pi,eps,rmax,ffhatnorm
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      real *8 kgridx(ntot)
      real *8 kgridy(ntot)
      real *8 kgridz(ntot)
      real *8 wts(ntot)
      complex *16 ff(ngrid,ngrid,ngrid)
      complex *16 ffhat(ntot)
c
      pi = 4.0d0*datan(1.0d0)
c
      call mkfouriergrid(rmax,ngridr,ityper,
     1                    nlats,itypep,ntot,numonsphere,
     2                    kgridx,kgridy,kgridz,wts)
      ifpr = 0
      if (ifpr.eq.1) call prinf('ntot is *',ntot,1)
      if (ifpr.eq.1) call prinf('numonsphere is *',numonsphere,ngridr)
      npts = ngrid*ngrid*ngrid
c
c     construct fourier transform of ff
c
      if (ifpr.eq.1) call prinf(' (gr) npts is *',npts,1)
      iflag = 1
      call  finufft3d3_f(npts,x,y,z,ff,iflag,eps,ntot,
     1        kgridx,kgridy,kgridz,ffhat,ier)
c
c     rescale to account for trapezodal rule quadrature weight
c
      ffhatnorm = 0.0d0
      do ii = 1,ntot
         ffhat(ii) = ffhat(ii)*h*h*h
         ffhatnorm = ffhatnorm + wts(ii)*cdabs(ffhat(ii))**2
      enddo
      ffhatnorm = dsqrt(ffhatnorm)/(dsqrt(2*pi)**3)
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_thetas_phis(nlats,itypep,numonsphere,thetas,phis)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     PURPOSE: Create theta/phi values for pts on sphere
c              with nlats latitudes. 
c     
c     INPUT:
c
c     nlats          number of latitude points on sphere
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is defined in getgridph and typically
c                    set to 2.
c
c     OUTPUT:
c
c     numonsphere    number of points on sphere
c     thetas         theta components of pts on sphere                    
c     phis           phi components of pts on sphere                    
c
c----------------------------------------------------------------------
      implicit none
      integer ngridp,ngridt,i
      integer itypep,ntot,nnn
      integer ifprint,jj,kk,ll,ngridtmax
      integer nlats
      integer numonsphere
      real *8 rkmax,u,v,pi,phi,phstep,rr
      real *8 ctheta,stheta
      real *8 phis(numonsphere)
      real *8 thetas(numonsphere)
      integer, allocatable ::  ngridps(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: wtsth(:)
c
      pi = 4.0d0*datan(1.0d0)
c
c     get max theta array size and allocate local arrays
c
      allocate(xnodesth(nlats))
      allocate(sthetas(nlats))
      allocate(wtsth(nlats))
      allocate(ngridps(nlats))
      allocate(phsteps(nlats))
c
      nnn = 0
      ngridt = nlats
      call getspheregrid(ngridt,itypep,xnodesth,
     1         sthetas,wtsth,ngridps,phsteps,numonsphere)
      do kk = 1,ngridt
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
c
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            phis(nnn) = phi
            thetas(nnn) = dacos(ctheta)
         enddo
      enddo
      return
      end
