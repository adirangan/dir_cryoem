      subroutine testfourier(ngrid,ngridr,rmax,ityper,itypep,
     1            nlats,fgauss,ns,x0y0z0,sigma,err,errback)
c
c
c     testing Fourier transform routines:
c
c     create function on physical space grid
c     create spherical grid in k-space (tensor product or sparse)
c     compute Fourier transform on spherical grid using nufft type 3.
c     compare with exact transform
c     convert to spherical harmonic basis and comute norm in that basis
c     convert from spherical harmonic basis to grid
c     compute inverse Fourier transform back to rectangular grid 
c             using nufft type 3.
c     compare with original function
c
      implicit real *8 (a-h,o-z)
      integer nlats(ngridr)
      integer, allocatable :: numonsphere(:)
      integer, allocatable :: nterms_sph(:)
      integer, allocatable :: isph_start(:)

      real *8 x0y0z0(3,ns)
      real *8 fnorm,rmodelnorm,phi_over
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      real *8, allocatable :: xnodesr(:),wtsr(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ff(:,:,:)
      complex *16, allocatable :: ff2(:,:,:)
      complex *16, allocatable :: ffhat(:)
      complex *16, allocatable :: modsph(:)
      complex *16 eye,uex
      external fgauss
c
      eye = dcmplx(0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      call prini(6,13)
c
c     ifpr is printing flag (1 = on)
c
      ifpr = 1
c
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
      allocate(ff(ngrid,ngrid,ngrid))
      allocate(ff2(ngrid,ngrid,ngrid))
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(xnodesth(nlats(ngridr)))
      allocate(wtsth(nlats(ngridr)))
      allocate(xnodesph(2*nlats(ngridr)))
      allocate(wtsph(2*nlats(ngridr)))
c
c     create physical space grid and sample function on it.
c
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,ns,sigma,
     1     fgauss,ff,fnorm)
      ifpr = 1
      if (ifpr.eq.1) write(6,*) '  h is ',h
      if (ifpr.eq.1) write(6,*) '  fnorm is ',fnorm
c
c     Determine size of spherical grid, allocate space and compute
c     Fourier transform at corresponding points in k-space.
c 
      call get_ntot(ngridr,nlats,itypep,ntot)
      if (ifpr.eq.1) write(6,*) '  ntot is ',ntot
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(ffhat(ntot))
c
c     Determine size of spherical harmonic model 
c
      allocate(numonsphere(ngridr))
      allocate(nterms_sph(ngridr))
      allocate(isph_start(ngridr))
      nsphstore = 0
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)
      do i = 1,ngridr
         isph_start(i) = nsphstore
         nterms_sph(i) = nint(xnodesr(i) + 2)
         nsphstore = nsphstore + (nterms_sph(i)+1)**2
      enddo
      call prinf(' nterms_sph is *',nterms_sph,ngridr)
      call prinf(' nsphstore is *',nsphstore,1)
      allocate(modsph(0:nsphstore-1))
c
      eps = 1.0d-6
c
      call fun_to_kspacegrid(x,y,z,ff,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat,ffhatnorm,ier)
      if (ifpr.eq.1) write(6,*) '  ffhatnorm is ',ffhatnorm
c
c     check error in ffhat using exact transform of Gaussian
c
      err = 0 
      stot = 0 
      do ii = 1,ntot
         rkx = kgridx(ii)
         rky = kgridy(ii)
         rkz = kgridz(ii)
         uex = 0.0d0
         do is = 1,ns
            arg = rkx*x0y0z0(1,is)+rky*x0y0z0(2,is)+rkz*x0y0z0(3,is)
            uex = uex+dexp(-(rkx*rkx+rky*rky+rkz*rkz)*sigma*sigma/2)*
     1           cdexp(eye*arg)*(2*pi)*dsqrt(2*pi)*sigma**3
         enddo
         err = err + abs(ffhat(ii)-uex)**2
         stot = stot + abs(uex)**2
      enddo
      err = sqrt(err/stot)
      if (ifpr.eq.1) call prin2(' l2 err in ffhat is *',err,1)
c
      call kspacegrid_to_model(ffhat,rmax,nterms_sph,
     1         ityper,ngridr,nlats,itypep,numonsphere,modsph,nsphstore)
      call kspace_model_norm(modsph,nsphstore,nterms_sph,wtsr,ngridr,
     1     rmodelnorm)
      if (ifpr.eq.1) write(6,*) '  rmodelnorm is ',rmodelnorm
      call model_to_kspacegrid(ffhat,rmax,nterms_sph,
     1         ityper,ngridr,nlats,itypep,numonsphere,modsph,nsphstore)
c
c    now carry out inverse transform and check against original function.
c
      do ii = 1,ntot
         ffhat(ii) = ffhat(ii)*wts(ii)
      enddo
      iflag = -1
      npts = ngrid*ngrid*ngrid
      call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat,iflag,eps,
     1        npts,x,y,z,ff2,ier)
c
      errback = 0 
      stot = 0 
      do i = 1,ngrid
      do j = 1,ngrid
      do k = 1,ngrid
          ff2(i,j,k) = ff2(i,j,k)/(8*pi*pi*pi)
          errback = errback + abs(ff2(i,j,k)-ff(i,j,k))**2
          stot = stot + abs(ff(i,j,k))**2
      enddo
      enddo
      enddo
      errback = sqrt(errback/stot)
      if (ifpr.eq.1) call prin2(' l2 err in ff is *',errback,1)
      return
      end
C
      subroutine createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,ns,sigma,
     1     fgauss,ff,fnorm)
c
      implicit real *8 (a-h,o-z)
      integer ngrid,ns
      real *8 a,h,sigma,pi
      real *8 x0y0z0(3,ns)
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      complex *16 ff(ngrid,ngrid,ngrid)
      external fgauss
c
c
c$$$      write(6,*) ' Calling: ','mkphysgrid'
      call mkphysgrid(a,ngrid,h,x,y,z)
c$$$      write(6,*) 'x: ',(x(i,1,1),i=1,ngrid)
c$$$      write(6,*) 'y: ',(y(1,j,1),j=1,ngrid)
c$$$      write(6,*) 'z: ',(z(1,1,k),k=1,ngrid)
c$$$      do i = 1,ngrid
c$$$         do j = 1,ngrid
c$$$            do k = 1,ngrid
c$$$               write(6,*) 'i,j,k: ',i,j,k,'; x,y,z: ',x(i,j,k),y(i,j,k)
c$$$     $              ,z(i,j,k)
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      write(6,*) ' Calling: ','evalgrid'
      call evalgrid(x,y,z,ngrid,ngrid,ngrid,x0y0z0,ns,sigma,
     1     fgauss,ff)
c 
c$$$      write(6,*) ' Computing: ','fnorm'
      fnorm = 0.0d0
      do i = 1,ngrid
      do j = 1,ngrid
      do k = 1,ngrid
         fnorm = fnorm + h*h*h*cdabs(ff(i,j,k))**2
      enddo
      enddo
      enddo
      fnorm = dsqrt(fnorm)
      return
      end

