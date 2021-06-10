C
c      Testing code for rebuilding model from slices.
c
      implicit real *8 (a-h,o-z)
c
      parameter (ngridr=20)
      parameter (ncur=10)
      parameter (ngrid=60)
      parameter (ntmax=1000)
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      integer isph_start(ngridr)
      integer ngridc(2*ntmax)
      integer ngridc2(2*ntmax)
      integer icstart(ngridr)
      integer ngridps(ntmax)
      real *8 alphas(3,100000)
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8 xnodesth(ntmax)
      real *8 sthetas(ntmax)
      real *8 wtsth(ntmax)
      real *8 phsteps(ntmax)
      real *8 source(3,10)
      real *8 x0y0z0(3,10),sigma
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      complex *16 ff(ngrid,ngrid,ngrid)
ccc      complex *16 cslices(350,1232)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ffhat(:)
      complex *16, allocatable :: modsph(:)
      complex *16, allocatable :: modsph2(:)
      complex *16, allocatable :: cslices(:,:)
      complex *16, allocatable :: templates(:,:)
c
      external fgauss
c
c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
c
      done=1
      pi=4.0d0*atan(done)
      call prini(6,13)
c
c     set sigma and location of Gaussian sources (x0,y0,z0)
c     ns = number of sources
c
      ns = 2
      sigma = 1.0d-1
      x0y0z0(1,1) = 0.1d0
      x0y0z0(2,1) = 0.2d0
      x0y0z0(3,1) = -0.13d0
      x0y0z0(1,2) = -0.1d0
      x0y0z0(2,2) = -0.2d0
      x0y0z0(3,2) = 0.23d0
c
c     create physical space grid and sample function on it.
c
      rmax = 10.0d0
      ityper = 0
      itypep = 1
      ifpr = 1
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,ns,sigma,
     1     fgauss,ff,fnorm)
c
      if (ifpr.eq.1) write(6,*) '  fnorm is ',fnorm
c
c-----------------------------------------------------------------------
c     2) Determine size of spherical grid, allocate space and compute
c     Fourier transform at corresponding points in k-space.
c     Also determine size of spherical harmonic model  (nsphstore)
c-----------------------------------------------------------------------
c 
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)
      nsphstore = 0
      nslice_size = 0
      do i = 1,ngridr
         nlats(i) = nint(pi*xnodesr(i))
         if (nlats(i).lt.6) nlats(i) = 6
         ngridc(i) = nlats(i)*2
         isph_start(i) = nsphstore
         nterms_sph(i) = nint(xnodesr(i) + 2)
         nsphstore = nsphstore + (nterms_sph(i)+1)**2
         nslice_size = nslice_size + ngridc(i)
      enddo
      if (ifpr.eq.1) call prinf('  nlats is *',nlats,ngridr)
      if (ifpr.eq.1) call prinf('  ngridc is *',ngridc,ngridr)
      if (ifpr.eq.1) call prinf('  nslice_size is *',nslice_size,1)
c
      call get_ntot(ngridr,nlats,itypep,ntot)
      if (ifpr.eq.1) write(6,*) '  ntot is ',ntot
      if (ifpr.eq.1) call prinf(' nterms_sph is *',nterms_sph,ngridr)
      if (ifpr.eq.1) call prinf(' nsphstore is *',nsphstore,1)
c
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(ffhat(ntot))
c
c
      allocate(modsph(0:nsphstore-1))
      allocate(modsph2(0:nsphstore-1))
c
c     create kspace grid data from physical space model
c
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat,ffhatnorm,ier)
      if (ifpr.eq.1) write(6,*) '  ffhatnorm is ',ffhatnorm
c
c     convert to model: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(ffhat,rmax,nterms_sph,
     1         ityper,ngridr,nlats,itypep,numonsphere,modsph,nsphstore)
c
c-----------------------------------------------------------------------
c     3) Assign slices to orientations of spherical grid on highest
c        frequency sphere and define gammas.
c-----------------------------------------------------------------------
c
      boxsize = 2.0d0
      nslices = numonsphere(ngridr)
      call getspheregrid(nlats(ngridr),itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,nspherebig)
      if (ifpr.eq.1) write(6,*) '  nslices is ',nslices
      if (ifpr.eq.1) write(6,*) '  nspherebig is ',nspherebig
c
      nnn = 0
      do kk = 1,nlats(ngridr)
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            alphas(1,nnn) = dacos(ctheta)
            alphas(2,nnn) = phi
            alphas(3,nnn) = -rand()
ccc            alphas(3,nnn) = 0.0d0
         enddo
      enddo
      if (ifpr.eq.1) write(6,*) '  nnn is ',nnn
c
      allocate(cslices(nslice_size,nslices))
      allocate(templates(nslice_size,nslices))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     1           icstart)
      if (ifpr.eq.1) write(6,*) '  ntemplatesize is ',ntemplatesize
c
c-----------------------------------------------------------------------
c     4) compute slices directly from physical space function (ff).
c        with orientation vectors defined by alphas.
c-----------------------------------------------------------------------
c
      call mk_simulated_slices(ff,ngrid,boxsize,eps,ngridr,
     1       ityper,ngridc,rmax,alphas,nslices,nslice_size,cslices)
c
c-----------------------------------------------------------------------
c     5) TEST 1: generate templates (should match slices since we took slices
c     from the same grid, so long as gamma (third Euler angle is set to zero!!)
c-----------------------------------------------------------------------
c
      call template_gen(modsph,nsphstore,isph_start,
     1     nterms_sph,ngridr,ngridc,ntemplatesize,icstart,
     2     nlats,itypep,ngridr,nslices,templates)

ccc      write(6,*) ' cslices(2,1) = ',cslices(2,1)
ccc      write(6,*) ' templates(2,1) = ',templates(2,1)
c
      err = 0.0d0
      denom = 0.0d0
      do jj = 1,nslices
      do ii = 1,nslice_size
         denom = denom + abs(cslices(ii,jj))**2
         err = err + abs(cslices(ii,jj) - templates(ii,jj))**2
      enddo
      enddo
      write(6,*) ' err in template_gen vs mk_simulated_slices = ',
     1            dsqrt(err/denom)
c
c-----------------------------------------------------------------------
c     6) TEST 2: rebuild model from slices and Euler angles and compare to 
c     original model.
c-----------------------------------------------------------------------
      oversamp = 2.0d0
      eps = 1.0d-6
      kord = 5
      nlow = 3
      call rebuild_full_model(cslices,nslice_size,nslices,alphas,
     1           icstart,nlats,itypep,ngridc,nlow,ncur,nsphstore,
     2           isph_start,nterms_sph,oversamp,kord,eps,modsph2)
c
      ncheck = 0
      do i = 1,ncur
         ncheck = ncheck + (nterms_sph(i)+1)**2
      enddo
      nlowstart = 0
      do i = 1,nlow-1
         nlowstart = nlowstart + (nterms_sph(i)+1)**2
      enddo
c
      do i = 0,nlowstart-1
         modsph(i) = 0.0d0
      enddo
      do i = nlowstart,ncheck-1
         modsph(i) = modsph(i) - modsph2(i)
      enddo
      call kspace_model_norm(modsph,ncheck,nterms_sph,wtsr,ncur,
     1          rmodelnorm)
      call kspace_model_norm(modsph2,ncheck,nterms_sph,wtsr,ncur,
     1          rmodelnorm2)
      write(6,*) 'rel error in rebuilt model is ',rmodelnorm/rmodelnorm2

      stop
      end
c
c

