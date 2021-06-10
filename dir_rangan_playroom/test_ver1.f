C
c      Testing code for rebuilding model from slices.
c
      implicit real *8 (a-h,o-z)
c     
      integer verbose
      data verbose / 1 /
      parameter (ngridr=20)
      parameter (ncur=10)
      parameter (ngrid=60)
      parameter (ntmax=2*ngridr)
      integer nlats(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      integer isph_start(ngridr)
      integer ngridc(ntmax)
      integer ngridc2(ntmax)
      integer icstart(ngridr)
      integer ngridps(ntmax)
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8 xnodesth(ntmax)
      real *8 sthetas(ntmax)
      real *8 wtsth(ntmax)
      real *8 phsteps(ntmax)
      integer n_source,n_source_use,nsource
      parameter (n_source=128)
      real *8 x0y0z0(3,n_source),sigma,source_r
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      complex *16 ff(ngrid,ngrid,ngrid)
      complex *16 ff2(ngrid,ngrid,ngrid)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ffhat(:)
      complex *16, allocatable :: ffhat2(:)
      complex *16, allocatable :: modsph(:)
      complex *16, allocatable :: modsph2(:)
      complex *16, allocatable :: cslices(:,:)
      complex *16, allocatable :: templates(:,:)
      integer d_(0:5)
      character(len=1024) fname,format_string
      real *8 max_x_c
      parameter(max_x_c=1.0d0)
      integer n_M_sample,nM_sample,n_S_sample,nS_sample
      integer *4, allocatable :: I_M_sample_(:)
      integer *4, allocatable :: I_S_sample_(:)
      complex *16, allocatable :: M_sample_(:)
      complex *16, allocatable :: S_sample_(:)
      real *8, allocatable :: alpha5(:,:)
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      real *8, allocatable :: delta_x_sort_(:)
      real *8, allocatable :: delta_y_sort_(:)
      real *8, allocatable :: gamma_z_sort_(:)
      complex *16, allocatable :: C_Z_sort_(:)
      integer *4, allocatable :: I_permute_(:)
      integer *4, allocatable :: I_inverse_(:)
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
c$$$      n_source_use = 32             ! number of gaussian sources
c$$$      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
c$$$      do nsource=0,n_source_use-1
c$$$         source_r = 0.5d0*(nsource)/max(1,n_source_use-1);
c$$$         x0y0z0(1,1+nsource) = source_r*cos(2*pi*2*nsource
c$$$     $        /n_source_use)
c$$$         x0y0z0(2,1+nsource) = source_r*sin(2*pi*2*nsource
c$$$     $        /n_source_use)
c$$$         x0y0z0(3,1+nsource) = 2.0d0*source_r - 0.5 
c$$$      enddo
      n_source_use = 128             ! number of gaussian sources
      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
      do nsource=0,n_source_use-1
         source_r = 2.0d0*pi*(nsource)/max(1,n_source_use)
         x0y0z0(1,1+nsource) = cos(2.0d0*source_r)*(3.0d0+cos(3.0d0
     $        *source_r))/8.0d0
         x0y0z0(2,1+nsource) = sin(2.0d0*source_r)*(3.0d0+cos(3.0d0
     $        *source_r))/8.0d0
         x0y0z0(3,1+nsource) = 2.0d0*sin(3.0d0*source_r)/8.0d0
      enddo

      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './source_xyz_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,x0y0z0,fname)

c
c     create physical space grid and sample function on it.
c
      rmax = 1.0d0*ncur
      ityper = 0
      itypep = 1
      ifpr = 1
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,n_source_use,sigma,
     1     fgauss,ff,fnorm)
c
      if (verbose.gt.1) write(6,*) '  fnorm is ',fnorm

      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
      write(fname,'(A)') './x_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x,fname)
      write(fname,'(A)') './y_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y,fname)
      write(fname,'(A)') './z_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './S_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff,fname)

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
      if (verbose.gt.1) call prinf('  nlats is *',nlats,ngridr)
      if (verbose.gt.1) call prinf('  ngridc is *',ngridc,ngridr)
      if (verbose.gt.1) call prinf('  nslice_size is *',nslice_size,1)
c
      call get_ntot(ngridr,nlats,itypep,ntot)
      if (verbose.gt.1) write(6,*) '  ntot is ',ntot
      if (verbose.gt.1) call prinf(' nterms_sph is *',nterms_sph,ngridr)
      if (verbose.gt.1) call prinf(' nsphstore is *',nsphstore,1)
c
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(ffhat(ntot))
      allocate(ffhat2(ntot))
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
      if (verbose.gt.1) write(6,*) '  ffhatnorm is ',ffhatnorm
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
      if (verbose.gt.0) write(6,*) '  nslices is ',nslices
      if (verbose.gt.1) write(6,*) '  nspherebig is ',nspherebig
c      
      allocate(alpha5(5,nslices))
      nnn = 0
      do kk = 1,nlats(ngridr)
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            alpha5(1,nnn) = dacos(ctheta)
            alpha5(2,nnn) = phi
c$$$            alpha5(3,nnn) = 0.0d0
c$$$            alpha5(4,nnn) = 0.0d0
c$$$            alpha5(5,nnn) = 0.0d0
            alpha5(3,nnn) = 0.5d0 - rand()
            alpha5(4,nnn) = (1.0d0 - 2*rand())/rmax
            alpha5(5,nnn) = (1.0d0 - 2*rand())/rmax
c$$$            alpha5(3,nnn) = 8.0d0*pi/32
c$$$            alpha5(4,nnn) = 1.0d0/rmax
c$$$            alpha5(5,nnn) = 0.0d0/rmax
         enddo
      enddo
      if (verbose.gt.1) write(6,*) '  nnn is ',nnn
c
      allocate(cslices(nslice_size,nslices))
      allocate(templates(nslice_size,nslices))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     1           icstart)
      if (verbose.gt.1) write(6,*) '  ntemplatesize is ',ntemplatesize
c
c-----------------------------------------------------------------------
c     4) compute slices directly from physical space function (ff).
c        with orientation vectors defined by alpha5.
c-----------------------------------------------------------------------
c
      call mk_simulated_slices_alpha5(ff,ngrid,boxsize,eps,ngridr,
     1       ityper,ngridc,rmax,alpha5,nslices,nslice_size,cslices)
c
c-----------------------------------------------------------------------
c     5) TEST 1: generate templates; should match slices since we took slices
c     from the same grid, so long as gamma_z (third Euler angle)
c     and the  delta_x and delta_y (displacements) are set to zero.
c-----------------------------------------------------------------------
c
      call template_gen(modsph,nsphstore,isph_start,
     1     nterms_sph,ngridr,ngridc,ntemplatesize,icstart,
     2     nlats,itypep,ngridr,nslices,templates)

      n_M_sample = 64
      allocate(I_M_sample_(0:n_M_sample-1))
      allocate(M_sample_(0:nslice_size*n_M_sample-1))
      allocate(alpha_tru_(0:5*n_M_sample-1))
      allocate(alpha_est_(0:5*n_M_sample-1))
      do nM_sample=0,n_M_sample-1
         I_M_sample_(nM_sample) = floor(nslices*1.0d0*nM_sample
     $        /(n_M_sample))
         call cp1_c16(nslice_size,cslices(1,1+I_M_sample_(nM_sample))
     $        ,M_sample_(nM_sample*nslice_size))
         alpha_tru_(0 + 5*nM_sample) = alpha5(1,1
     $        +I_M_sample_(nM_sample))
         alpha_tru_(1 + 5*nM_sample) = alpha5(2,1
     $        +I_M_sample_(nM_sample))
         alpha_tru_(2 + 5*nM_sample) = alpha5(3,1
     $        +I_M_sample_(nM_sample))
         alpha_tru_(3 + 5*nM_sample) = alpha5(4,1
     $        +I_M_sample_(nM_sample))
         alpha_tru_(4 + 5*nM_sample) = alpha5(5,1
     $        +I_M_sample_(nM_sample))
         alpha_est_(0 + 5*nM_sample) = 0.0d0
         alpha_est_(1 + 5*nM_sample) = 0.0d0
         alpha_est_(2 + 5*nM_sample) = 0.0d0
         alpha_est_(3 + 5*nM_sample) = 0.0d0
         alpha_est_(4 + 5*nM_sample) = 0.0d0
      enddo

      n_S_sample = 64
      allocate(I_S_sample_(0:n_S_sample-1))
      allocate(S_sample_(0:ntemplatesize*n_S_sample-1))
      do nS_sample=0,n_S_sample-1
         I_S_sample_(nS_sample) = floor(nslices*1.0d0*nS_sample
     $        /(n_S_sample))
         call cp1_c16(ntemplatesize,templates(1,1
     $        +I_S_sample_(nS_sample)),S_sample_(nS_sample
     $        *ntemplatesize))
      enddo

      allocate(delta_x_sort_(0:n_S_sample*n_M_sample-1))
      allocate(delta_y_sort_(0:n_S_sample*n_M_sample-1))
      allocate(gamma_z_sort_(0:n_S_sample*n_M_sample-1))
      allocate(C_Z_sort_(0:n_S_sample*n_M_sample-1))
      allocate(I_permute_(0:n_S_sample*n_M_sample-1))
      allocate(I_inverse_(0:n_S_sample*n_M_sample-1))

      if (verbose.gt.1) then
         write(6,*) 'Analyzing a few figures...'
         write(6,*) 'I_M_sample_: ',(I_M_sample_(nM_sample),nM_sample=0
     $        ,n_M_sample-1)
         write(6,*) 'I_S_sample_: ',(I_S_sample_(nS_sample),nS_sample=0
     $        ,n_S_sample-1)
      end if

c$$$         call Fig_gen_ver0(ncur,nlats,xnodesr,n_M_sample,ntemplatesize
c$$$     $        ,S_sample_,n_M_sample,nslice_size,M_sample_)
c$$$         call test_innerproduct_batch_wrapper0(ncur,nlats,xnodesr
c$$$     $        ,max_x_c,n_S_sample,ntemplatesize,S_sample_,n_M_sample
c$$$     $        ,nslice_size,M_sample_,alpha_tru_,alpha_est_,0.1d0,1.5d0
c$$$     $        ,17,17,max(32,2*ncur),-1)
      call test_innerproduct_batch_wrapper1(verbose,ncur,nlats,xnodesr
     $     ,max_x_c,n_S_sample,ntemplatesize,S_sample_,n_M_sample
     $     ,nslice_size,M_sample_,alpha_est_,0.1d0,1.5d0,17,17,max(32,2
     $     *ncur),-1,delta_x_sort_,delta_y_sort_,gamma_z_sort_
     $     ,C_Z_sort_,I_permute_,I_inverse_)
      if (verbose.gt.1) then
         write(6,'(A)') 'delta_x_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample,'(F8.3,1X))'
         write(6,format_string) (delta_x_sort_(na),na=0,n_S_sample
     $        *n_M_sample-1)
         write(6,'(A)') 'delta_y_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample,'(F8.3,1X))'
         write(6,format_string) (delta_y_sort_(na),na=0,n_S_sample
     $        *n_M_sample-1)
         write(6,'(A)') 'gamma_z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample,'(F8.3,1X))'
         write(6,format_string) (gamma_z_sort_(na),na=0,n_S_sample
     $        *n_M_sample-1)
         write(6,'(A)') 'C_Z_sort^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample*2,'(F8.3,1X))'
         write(6,format_string) (C_Z_sort_(na),na=0,n_S_sample
     $        *n_M_sample-1)
         write(6,'(A)') 'I_permute_^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample,'(I8,1X))'
         write(6,format_string) (I_permute_(na),na=0,n_S_sample
     $        *n_M_sample-1)
         write(6,'(A)') 'I_inverse_^T: '
         write(format_string,'(A,I0,A)') '(',n_S_sample,'(I8,1X))'
         write(6,format_string) (I_inverse_(na),na=0,n_S_sample
     $        *n_M_sample-1)
      end if

      goto 10


      err = 0.0d0
      denom = 0.0d0
      do jj = 1,nslices
      do ii = 1,nslice_size
         denom = denom + abs(cslices(ii,jj))**2
         err = err + abs(cslices(ii,jj) - templates(ii,jj))**2
      enddo
      enddo
      if (verbose.gt.1) then
         write(6,*) ' err in template_gen vs mk_simulated_slices = ',
     $        dsqrt(err/denom)
      end if
c
c-----------------------------------------------------------------------
c     6) TEST 2: rebuild model from slices and Euler angles and compare to 
c     original model.
c-----------------------------------------------------------------------
      oversamp = 2.0d0
      eps = 1.0d-6
      kord = 5
      nlow = 3
      if (verbose.gt.1) then
         write(6,*) 'nslice_size: ',nslice_size
         write(6,*) 'nslices: ',nslices
         write(6,*) 'icstart: ',(icstart(ii),ii=1,ncur)
         write(6,*) 'nlats: ',(nlats(ii),ii=1,ncur)
         write(6,*) 'itypep: ',itypep
         write(6,*) 'xnodesr: ',(xnodesr(ii),ii=1,ncur)
         write(6,*) 'ngridc: ',(ngridc(ii),ii=1,ncur)
         write(6,*) 'nlow: ',nlow
         write(6,*) 'ncur: ',ncur
         write(6,*) 'nsphstore: ',nsphstore
         write(6,*) 'isph_start: ',(isph_start(ii),ii=1,ncur)
         write(6,*) 'nterms_sph: ',(nterms_sph(ii),ii=1,ncur)
         write(6,*) 'oversamp: ',oversamp
         write(6,*) 'kord: ',kord
         write(6,*) 'eps: ',eps
      end if
      call rebuild_model_alpha5(cslices,nslice_size,nslices,alpha5,
     $     icstart,nlats,itypep,xnodesr,ngridc,nlow,ncur,nsphstore,
     $     isph_start,nterms_sph,oversamp,kord,eps,modsph2)

      call model_to_kspacegrid(ffhat2,rmax,nterms_sph, ityper,ngridr
     $     ,nlats,itypep,numonsphere,modsph2,nsphstore)

      do ii = 1,ntot
         ffhat2(ii) = ffhat2(ii)*wts(ii)
      enddo
      iflag = -1
      npts = ngrid*ngrid*ngrid
      call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat2,iflag,eps,
     $     npts,x,y,z,ff2,ier)

      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
      write(fname,'(A)') './x_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x,fname)
      write(fname,'(A)') './y_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y,fname)
      write(fname,'(A)') './z_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './S_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff2,fname)

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
      if (verbose.gt.0) then
         write(6,*) 'rel error in rebuilt model is ',rmodelnorm
     $        /rmodelnorm2
      end if

 10   continue
      stop
      end
