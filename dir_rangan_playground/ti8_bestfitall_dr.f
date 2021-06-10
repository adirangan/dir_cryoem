c      Testing code for bestfit
c
      implicit real *8 (a-h,o-z)
c
      parameter (ntmax=100000)
      integer verbose
      integer ngrid
      integer ngridr, ncur
      real *8 source(3,100)
      real *8 x0y0z0(3,100),sigma
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      complex *16, allocatable :: ff(:,:,:)
      complex *16, allocatable :: ff1(:,:,:)
      complex *16, allocatable :: ff2(:,:,:)

      integer, allocatable ::  nlats(:)
      integer, allocatable :: numonsphere(:)
      integer numsphmax
      integer, allocatable :: nterms_sph(:)
      integer, allocatable :: isph_start(:)
      integer ngridc(2*ntmax)
      integer ngridc2(2*ntmax)
      integer, allocatable :: icstart(:)
      integer, allocatable :: icstart2(:)
      integer ngridps(ntmax)
      real *8, allocatable :: alphas(:,:)
      real *8, allocatable :: alphas2(:,:)
      real *8, allocatable :: alphas2_copy(:,:)

      real *8 trange,tstep,h
      real *8, allocatable :: alltrans(:,:)
      integer ndelta

      real *8 smtheta,smphi,smgamma,smdeltax,smdeltay
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      real *8 xnodesth(ntmax)
      real *8 sthetas(ntmax)
      real *8 wtsth(ntmax)
      real *8 phsteps(ntmax)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: thetas(:)
      real *8, allocatable :: phis(:)
      complex *16, allocatable :: ffhat(:)
      complex *16, allocatable :: ffhat1(:)
      complex *16, allocatable :: ffhat2(:)
      complex *16, allocatable :: modsph(:)
      complex *16, allocatable :: modsph1(:)
      complex *16, allocatable :: modsph2(:)
      complex *16, allocatable :: modsphdiff(:)
      integer nslices
      complex *16, allocatable :: cslices(:,:)
      integer ntemplates
      complex *16, allocatable :: templates(:,:)
c
      external fgauss

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     CTF parameters
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer nctf
c     spherical aberation of the lens
      real *8 Cs
      real *8, allocatable :: SphericalAberration(:)
c     voltage
      real *8 Kv
      real *8,allocatable :: Voltage(:)
c     amplitude contrast
      real *8 AmpCnst,w1,w2
      real *8, allocatable :: AmplitudeContrast(:)
c     magnification of microscope
      real *8 xmag
      real *8, allocatable :: Magnification(:)
c     pixel size of scanner in microns
      real *8 dstep
      real *8 bfactor
      real *8, allocatable :: DetectorPixelSize(:)
c     defocus (in Angstroms) and angle of astigmatism
      real *8 df1,df2,angast
      real *8, allocatable :: DefocusU(:)
      real *8, allocatable :: DefocusV(:)
      real *8, allocatable :: DefocusAngle(:)

      complex *16, allocatable :: ctfw(:,:)
      integer, allocatable :: ctf_ind(:)
      integer cn, idx
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer rr
      real *8 sm

      integer ngridtt,ngridpp,ngridcc
      real *8 rcs,csm
      character(len=100) :: outdir
      character(len=100) :: fffile
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     size of the box in which the real-space function lives
      boxsize = 2.0d0
      a = 1.0d0
c     number of grid points for the real-space function
      ngrid = 100
c     how far to go out in frequency, i.e. max of radial component
      rmax = 40.0d0
c     number of radial quadrature nodes
      ngridr = 20
c     number of quadrature nodes in theta - should be >rmax by Nyquist
      ngridtt = 20
c     number of quadrature nodes in phi - should be >2rmax by Nyquist
      ngridpp = 40
c     number of points on radial circles for the projections, >2rmax by Nyquist
      ngridcc = 40
c
      outdir = 'output/'

c     oversampled at poles (0) or not (1)
      itypep = 1
c     Gaussian quadrature (1) or uniform grid (0) in radial direction
      ityper = 0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     allocate some space
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
      allocate(ff(ngrid,ngrid,ngrid))
      allocate(ff1(ngrid,ngrid,ngrid))
      allocate(ff2(ngrid,ngrid,ngrid))

      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(nlats(ngridr))
      allocate(numonsphere(ngridr))
      allocate(nterms_sph(ngridr))
      allocate(isph_start(ngridr))

      allocate(icstart(ngridr))
      allocate(icstart2(ngridr))
c
c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
      done=1
      pi=4.0d0*atan(done)
      call prini(6,13)
c
c     set sigma and location of Gaussian sources (x0,y0,z0)
c     ns = number of sources
c     set sigma and location of Gaussian sources (x0,y0,z0)
      sigma = 0.08d0
c     number of sources, i.e. Gaussians
      ns = 100
      do j = 1,ns
         x0y0z0(1,j) = (rand()-0.5d0)/1.5d0
         x0y0z0(2,j) = (rand()-0.5d0)/1.5d0
         x0y0z0(3,j) = (rand()-0.5d0)/1.5d0
      enddo
c     create physical space grid and sample function on it.
      call createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,ns,sigma,
     1     fgauss,ff,fnorm)
c
c-----------------------------------------------------------------------
c     2) Determine size of spherical grid, allocate space
c     Also determine size of spherical harmonic model  (nsphstore)
c-----------------------------------------------------------------------
c
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)
c     size of the whole array of spherical harmonic coefficients for all spheres
      nsphstore = 0
c     size of a single slice out to ngridr (including all rings at each k)
      nslice_size = 0
      do i = 1,ngridr
c     how many latitutes for each sphere determined by radius xnodesr(i)
         nlats(i) = nint(1.5d0*pi*xnodesr(i))
         if (nlats(i).lt.30) nlats(i) = 30
         if (nlats(i).ge.150) nlats(i) = 150
c     how many points on each ring of a slice
c     (each ring will have a different number of points depending on radius k)
         ngridc(i) = nint(2.0d0*pi*xnodesr(i))
         if(ngridc(i) .lt. 20) ngridc(i) = 20
         sm = ngridc(i)
         ngridc(i) = next235(sm)
c     where spherical harmonic expansion for each sphere starts
         isph_start(i) = nsphstore
         nterms_sph(i) = nint(xnodesr(i) + 2)
         nsphstore = nsphstore + (nterms_sph(i)+1)**2
         nslice_size = nslice_size + ngridc(i)
      enddo

c     where rings on a slice start
      icstart(1) = 1
      do i = 1,ngridr-1
         icstart(i+1) = icstart(i) + ngridc(i)
      enddo

c     get the number of all points in 3D Fourier on all spheres
      call get_ntot(ngridr,nlats,itypep,ntot)

      write(6,*) (nlats(i),i=1,ngridr)
      write(6,*) (ngridc(i),i=1,ngridr)
      write(6,*) 'nslice_size' , nslice_size
c      write(6,*) 'nsphstore', nsphstore
c      write(6,*) (nterms_sph(i),i=1,ngridr)
c      write(6,*) 'ntot', ntot
c      write(6,*) (nterms_sph(rr),rr=1,ngridr)
c      write(6,*) 'nsphstore',nsphstore

c-------------------------------------------------------
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(ffhat(ntot))
      allocate(ffhat1(ntot))
      allocate(ffhat2(ntot))
c
      allocate(modsph(0:nsphstore-1))
      allocate(modsph1(0:nsphstore-1))
      allocate(modsph2(0:nsphstore-1))
      allocate(modsphdiff(0:nsphstore-1))
c-------------------------------------------------------------
c
c      compute Fourier transform at corresponding points in k-space.
      t0=second()
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat,ffhatnorm,ier)
      t1 = second()
      write(6,*) 'time to generate fhat', t1-t0
c     convert to model: spherical harmonic expansions on successive
c     spheres
      t0=second()
      call kspacegrid_to_model(ffhat,rmax,nterms_sph,
     1         ityper,ngridr,nlats,itypep,numonsphere,modsph,nsphstore)
      t1 = second()
      write(6,*) 'time to generate model', t1-t0
c
      t0=second()
      call model_to_kspacegrid(ffhat2,rmax,nterms_sph,
     1         ityper,ngridr,nlats,itypep,numonsphere,modsph,nsphstore)
      t1 = second()
      write(6,*) 'time to go back to ffhat', t1-t0
c
      err = 0.0d0
      denom = 0.0d0
      do kk=1,ntot
         err=err+wts(kk)*abs(ffhat2(kk)-ffhat(kk))**2
         denom=denom+wts(kk)*abs(ffhat(kk))**2
      enddo
      write(6,*) 'ffhat err', sqrt(err/denom)
c     go back to physical space
      call kspacegrid_to_fun(x,y,z,ff2,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat2,ier)

      errback = 0
      stot = 0
      do i = 1,ngrid
      do j = 1,ngrid
      do k = 1,ngrid
           errback = errback + abs(ff2(i,j,k)-ff(i,j,k))**2
          stot = stot + abs(ff(i,j,k))**2
      enddo
      enddo
      enddo
      errback = sqrt(errback/stot)
      write(6,*) ' l2 err in ff is *',errback

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c choose CTF parameters
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nctf = 100

      allocate(Voltage(nctf))
      allocate(DefocusU(nctf))
      allocate(DefocusV(nctf))
      allocate(DefocusAngle(nctf))
      allocate(SphericalAberration(nctf))
      allocate(DetectorPixelSize(nctf))
      allocate(Magnification(nctf))
      allocate(AmplitudeContrast(nctf))

c     the CTF parameters are taken from the rib80s dataset
      Kv = 300.0d0
      Cs = 2.0d0
      dstep = 3.7687499999999998d0
      AmpCnst = 0.1d0
      bfactor = 50
      xmag = 1.0d0

c     maximum possible value of defocus
      max_df = 26462.15039099d0
c     minimum possible value of defocus
      min_df = 8112.450195d0
c     max astigmatism (meaning the difference between df1 and df2)
      max_ast = 457.03125d0
c     min astigmatism
      min_ast = 10.509765d0

      do cn=1,nctf

c     generate a value between the min and max defocus
         df2 = rand()*(max_df-min_df)
         df2 = df2 + min_df
c         write(6,*) df1
c     generate a value between max and min astigmatism
         diff = rand()*(max_ast-min_ast)
         diff = diff + min_ast

c     df1 > df1
         df1 = df2+diff
         df1 = df2

         Voltage(cn) = Kv
         DefocusU(cn) = df1
         DefocusV(cn) = df2
         DefocusAngle(cn) = angast
         SphericalAberration(cn) = Cs
         DetectorPixelSize(cn) = dstep
         Magnification(cn) = xmag
         AmplitudeContrast(cn) = AmpCnst
      enddo

      allocate(ctfw(nslice_size,nctf))
      call compute_ctf(nctf,nslice_size,SphericalAberration,
     1     Voltage,DefocusU,DefocusV,DefocusAngle,AmplitudeContrast,
     2     DetectorPixelSize,ngrid,rmax,ngridr,ityper,ngridc,ctfw)

c
c-----------------------------------------------------------------------
c     3) Assign slices to orientations of spherical grid on highest
c        frequency sphere and define gammas.
c-----------------------------------------------------------------------
c
      nslices = 100
      allocate(cslices(nslice_size,nslices))
      allocate(ctf_ind(nslices))
      allocate(alphas(5,nslices))
      allocate(alphas2(5,nslices))
      allocate(alphas2_copy(5,nslices))

      ncur = ngridr-1

c     choose angles for slices
      call getspheregrid(nlats(ncur),itypep,xnodesth,
     1     sthetas,wtsth,ngridps,phsteps,nspherebig)
      h = 2.0d0/ngrid
      trange = 1.0d0*h
      do nnn=1,nslices
c     values on the grid
         kk = int(rand()*nlats(ncur)-1)+1
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
         ll = int(rand()*ngridps(kk)-1)+1
         phi = (ll-1)*phstep
         jj = int(rand()*ngridc(ncur)-1)+1
         gamma = 2*pi*(jj-1)/ngridc(ncur)
         alphas(1,nnn) = dacos(ctheta)
         alphas(2,nnn) = phi
         alphas(3,nnn) = gamma
         alphas(4,nnn) = trange
         alphas(5,nnn) = -trange

c     values off the grid
         alphas(1,nnn) = dacos(rand()*2.0d0-1.0d0)
         alphas(2,nnn) = rand()*2*pi
         alphas(3,nnn) = rand()*2*pi
         alphas(4,nnn) = rand()*trange
         alphas(5,nnn) = -rand()*trange
      enddo
c
c-----------------------------------------------------------------------
c     4) compute slices directly from physical space function (ff).
c        with orientation vectors defined by alphas.
c-----------------------------------------------------------------------
      eps = 1.0d-6
      t0 = second()
      call mk_simulated_slices(ff,ngrid,boxsize,eps,ngridr,
     1       ityper,ngridc,rmax,alphas,nslices,nslice_size,cslices)
      t1 = second()
      write(6,*) ' time for simulating slices is ',t1-t0

c      do cn = 1,nctf
c         do ir = 1,nslice_size
c            ctfw(ir,cn) = 1.0d0
c         enddo
c      enddo

      do cn = 1,nslices
c     the rounding is to the lower integer, so add 1
         idx = int(rand()*nctf)+1
         ctf_ind(cn) = idx
      do ir = 1,nslice_size
         cslices(ir,cn)=cslices(ir,cn)*ctfw(ir,ctf_ind(cn))
      enddo
      enddo

      call translate_fourier_images(nslices,nslice_size,rmax,
     1     ngridr,ityper,ngridc,alphas,cslices)

c-----------------------------------------------------------------------
c
      ncur =  1

c-------------------------------------------------------------------------
c     starting model
c-------------------------------------------------------------------------

      do nnn=1,nslices
c     true angles
         alphas2(1,nnn) = alphas(1,nnn)
         alphas2(2,nnn) = alphas(2,nnn)
         alphas2(3,nnn) = alphas(3,nnn)
         alphas2(4,nnn) = alphas(4,nnn)
         alphas2(5,nnn) = alphas(5,nnn)
c     random starting angles
c         alphas2(1,nnn) = dacos(rand()*2.0d0-1.0d0)
c         alphas2(2,nnn) = rand()*2*pi
c         alphas2(3,nnn) = rand()*2*pi
      enddo
c     solve for the model
      oversamp = 2.0d0
      eps = 1.0d-2
      kord = 5
      nlow = 1
      call rebuild_full_model(cslices,nslice_size,nslices,ctfw,nctf,
     1     ctf_ind,alphas2,icstart,nlats,xnodesr,ngridc,nlow,ngridr,
     2     nsphstore,isph_start,nterms_sph,oversamp,kord,eps,modsph2)

c     starting model up tp ncur
      ncheck = 0
      do i = 1,ncur
         ncheck = ncheck + (nterms_sph(i)+1)**2
      enddo
      do i = 0,ncheck-1
         modsph1(i) = modsph2(i)
      enddo

c---------------------------------------------------------------------
c     Main loop
c--------------------------------------------------------------------------

      do while(ncur < ngridr)
         write(6,*) 'ncur', ncur

         call getspheregrid(nlats(ncur),itypep,xnodesth,
     1        sthetas,wtsth,ngridps,phsteps,ntemplates)

         ntemplatesize=0
         do i=1,ncur
            ntemplatesize=ntemplatesize+ngridc(i)
         enddo
         write(6,*) 'ntemplates', ntemplates
         write(6,*) 'ntemplatesize', ntemplatesize, nslice_size

         allocate(templates(ntemplatesize,ntemplates))
         allocate(thetas(ntemplates))
         allocate(phis(ntemplates))

         t0 = time()
         call template_gen(modsph1,nsphstore,isph_start,
     1        nterms_sph,ngridr,ngridc,ntemplatesize,icstart,
     2        nlats,itypep,ncur,ntemplates,templates)
         t1 = time()
         write(6,*) ' time for template gen is ',t1-t0
c
         call get_thetas_phis(nlats(ncur),itypep,numsphmax,thetas,phis)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up translations

      tstep = 0.5d0*h
      ndelta = nint(2*trange/tstep+1)**2
      allocate(alltrans(2,ndelta))
      write(6,*) 'ndelta', ndelta
      call get_translations(trange,tstep,ndelta,alltrans)
c$$$      %%%%%%%% ;
c$$$      Maybe reduce number of translations for debugging. ;
c$$$      %%%%%%%% ;
c$$$      ndelta = 1 ; alltrans(1,1) = +0.050d0 ; alltrans(2,1) = -0.030d0
c$$$      ndelta = 1 ; alltrans(1,1) = +0.001d0 ; alltrans(2,1) = -0.002d0
c$$$      ndelta = 1 ; alltrans(1,1) = +0.000d0 ; alltrans(2,1) = -0.000d0

      verbose = 1
c$$$      %%%%%%%%
      t0 = time()
      call cp1_r8(5*nslices,alphas2,alphas2_copy)
      if (verbose.gt.1) then
         write(6,'(A)') ' first and last images: '
         call print_sub_c16(ntemplatesize,cslices(1,1),' first image: ')
         call print_sub_c16(ntemplatesize,cslices(1,nslices)
     $        ,'  last image: ')
         write(6,'(A)') ' first and last templates: '
         call print_sub_c16(ntemplatesize,templates(1,1)
     $        ,' first template: ')
         call print_sub_c16(ntemplatesize,templates(1,ntemplates)
     $        ,'  last template: ')
      end if !if (verbose.gt.1) then
      call adi_getbestfit(
     $     cslices,nslice_size,nslices,ctfw,nctf,
     1     ctf_ind,icstart,ncur,ngridc,alphas2_copy,xnodesr,
     2     wtsr,templates,ntemplatesize,ntemplates,
     3     thetas,phis,ndelta,alltrans)
      t1 = time()
         write(6,*) 'time to match the angles (adi_getbestfit)', t1-t0
         if (verbose.gt.-1) then
            write(6,'(A)') ' alpha_first5_ for first 10 images: '
            do i=1,min(10,nslices)
            call print_all_r8(5,alphas2_copy(1,i)
     $           ,' alphas2_copy(1,i): ')
            enddo !do i=1,min(10,nslices)
         end if !if (verbose.gt.1) then
      do i=1,nslices
         if (verbose.gt.2) then
         if (abs(alphas(2,i)-alphas2_copy(2,i)) > pi) then
            write(6,*) abs(alphas(2,i)-alphas2_copy(2,i)),
     1           alphas(2,i), alphas2_copy(2,i)
         endif
         end if !if (verbose.gt.2) then
      enddo
         smtheta=0.0d0
         smphi=0.0d0
         smgamma=0.0d0
         smdeltax = 0.0d0
         smdeltay = 0.0d0
         do i=1,nslices
            smtheta=smtheta+abs(alphas(1,i)-alphas2_copy(1,i))
            smphi=smphi+abs(alphas(2,i)-alphas2_copy(2,i))
            smgamma=smgamma+abs(alphas(3,i)-alphas2_copy(3,i))
            smdeltax=smdeltax+abs(alphas(4,i)-alphas2_copy(4,i))
            smdeltay=smdeltay+abs(alphas(5,i)-alphas2_copy(5,i))
         enddo
         write(6,*) 'error in theta', smtheta/nslices
         write(6,*) 'error in phi', smphi/nslices
         write(6,*) 'error in gamma', smgamma/nslices
         write(6,*) 'error in deltax', smdeltax/nslices
         write(6,*) 'error in deltay', smdeltay/nslices
c$$$  %%%%%%%%
      t0 = time()
      call cp1_r8(5*nslices,alphas2,alphas2_copy)
      if (verbose.gt.1) then
         write(6,'(A)') ' first and last images: '
         call print_sub_c16(ntemplatesize,cslices(1,1),' first image: ')
         call print_sub_c16(ntemplatesize,cslices(1,nslices)
     $        ,'  last image: ')
         write(6,'(A)') ' first and last templates: '
         call print_sub_c16(ntemplatesize,templates(1,1)
     $        ,' first template: ')
         call print_sub_c16(ntemplatesize,templates(1,ntemplates)
     $        ,'  last template: ')
      end if !if (verbose.gt.1) then
      call ti8_getbestfit(verbose,
     $     cslices,nslice_size,nslices,ctfw,nctf,
     1     ctf_ind,icstart,ncur,ngridc,alphas2_copy,xnodesr,
     2     wtsr,templates,ntemplatesize,ntemplates,
     3     thetas,phis,ndelta,alltrans)
      t1 = time()
         write(6,*) 'time to match the angles (ti8_getbestfit)', t1-t0
         if (verbose.gt.-1) then
            write(6,'(A)') ' alpha_first5_ for first 10 images: '
            do i=1,min(10,nslices)
            call print_all_r8(5,alphas2_copy(1,i)
     $           ,' alphas2_copy(1,i): ')
            enddo !do i=1,min(10,nslices)
         end if !if (verbose.gt.1) then
      do i=1,nslices
         if (verbose.gt.2) then
         if (abs(alphas(2,i)-alphas2_copy(2,i)) > pi) then
            write(6,*) abs(alphas(2,i)-alphas2_copy(2,i)),
     1           alphas(2,i), alphas2_copy(2,i)
         endif
         end if !if (verbose.gt.2) then
      enddo
         smtheta=0.0d0
         smphi=0.0d0
         smgamma=0.0d0
         smdeltax = 0.0d0
         smdeltay = 0.0d0
         do i=1,nslices
            smtheta=smtheta+abs(alphas(1,i)-alphas2_copy(1,i))
            smphi=smphi+abs(alphas(2,i)-alphas2_copy(2,i))
            smgamma=smgamma+abs(alphas(3,i)-alphas2_copy(3,i))
            smdeltax=smdeltax+abs(alphas(4,i)-alphas2_copy(4,i))
            smdeltay=smdeltay+abs(alphas(5,i)-alphas2_copy(5,i))
         enddo
         write(6,*) 'error in theta', smtheta/nslices
         write(6,*) 'error in phi', smphi/nslices
         write(6,*) 'error in gamma', smgamma/nslices
         write(6,*) 'error in deltax', smdeltax/nslices
         write(6,*) 'error in deltay', smdeltay/nslices

c$$$  %%%%%%%%
         call cp1_r8(5*nslices,alphas2_copy,alphas2)

         deallocate(templates)
         deallocate(thetas)
         deallocate(phis)
         deallocate(alltrans)

c      do nnn=1,nslices
cc     true angles
c         alphas2(1,nnn) = alphas(1,nnn)
c         alphas2(2,nnn) = alphas(2,nnn)
c         alphas2(3,nnn) = alphas(3,nnn)
c         alphas2(4,nnn) = alphas(4,nnn)
c         alphas2(5,nnn) = alphas(5,nnn)
c      enddo

c
c-----------------------------------------------------------------------
c     6) rebuild model from slices and Euler angles and compare to
c     original model.
c-----------------------------------------------------------------------
         ncur=ncur+1

         oversamp = 2.0d0
         eps = 1.0d-2
         kord = 5
         nlow = 1
         t0=second()
         call rebuild_full_model(cslices,nslice_size,nslices,ctfw,nctf,
     1        ctf_ind,alphas2,icstart,nlats,xnodesr,ngridc,nlow,ncur,
     2        nsphstore,isph_start,nterms_sph,oversamp,kord,eps,modsph2)
         t1=second()
         write(6,*) 'time to rebuild the model', t1-t0
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
            modsphdiff(i) = 0.0d0
         enddo
         do i = nlowstart,ncheck-1
            modsphdiff(i) = modsph(i) - modsph2(i)
         enddo
         call kspace_model_norm(modsphdiff,ncheck,nterms_sph,wtsr,ncur,
     1        rmodelnorm)
         call kspace_model_norm(modsph2,ncheck,nterms_sph,wtsr,ncur,
     1        rmodelnorm2)
         write(6,*) 'rel error in rebuilt model is ',
     1        rmodelnorm/rmodelnorm2

         do i = 0,ncheck-1
            modsph1(i) = modsph2(i)
         enddo

      enddo
c-------------------------------------------------------------------------------
c     end of the main loop
c-------------------------------------------------------------------------

      call model_to_kspacegrid(ffhat1,rmax,nterms_sph,
     1     ityper,ngridr,nlats,itypep,numonsphere,modsph1,nsphstore)

      err = 0.0d0
      denom = 0.0d0
      nstart = 1
      do rr = 1,ngridr
      do kk=1,numonsphere(rr)
         err=err+abs(ffhat1(nstart+kk-1)-
     1        ffhat(nstart+kk-1))**2
         denom=denom+abs(ffhat(nstart+kk-1))**2
      enddo
         nstart=nstart+numonsphere(rr)
      enddo
      write(6,*) 'ffhat err', sqrt(err/denom)

c    now carry out inverse transform and check against original function.
      call kspacegrid_to_fun(x,y,z,ff1,h,ngrid,eps,rmax,ngridr,
     1                  ityper,nlats,itypep,ntot,numonsphere,
     2                  kgridx,kgridy,kgridz,wts,ffhat1,ier)

      errback = 0
      stot = 0
      do i = 1,ngrid
      do j = 1,ngrid
      do k = 1,ngrid
          errback = errback + abs(ff1(i,j,k)-ff(i,j,k))**2
          stot = stot + abs(ff(i,j,k))**2
      enddo
      enddo
      enddo
      errback = sqrt(errback/stot)
      write(6,*) ' l2 err in ff is *',errback

c      fffile = trim(adjustl(outdir)) // 'ff1'
c      write(6,*) 'fffile', fffile
c      open(17,FILE=fffile)
c      do k=1,ngrid
c      do j=1,ngrid
c      do i=1,ngrid
c         write(17,*) x(i,j,k),y(i,j,k),z(i,j,k),real(ff1(i,j,k))
c      enddo
c      enddo
c      enddo
c      close(17)





      stop
      end
