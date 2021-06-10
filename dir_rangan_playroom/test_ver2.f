C
c      Testing code for rebuilding model from slices.
c
      implicit real *8 (a-h,o-z)
c     
      integer verbose
      data verbose / 1 /
      parameter (ngridr=20)
      integer ncur
      parameter (ngrid=60)
      parameter (ntmax=4*ngridr)
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
      integer ngridps_tmp(ntmax)
      real *8 xnodesth_tmp(ntmax)
      real *8 sthetas_tmp(ntmax)
      real *8 wtsth_tmp(ntmax)
      real *8 phsteps_tmp(ntmax)
      integer n_source_max,n_source_use,nsource
      parameter (n_source_max=128)
      real *8 x0y0z0(3,n_source_max),sigma,source_r
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      complex *16 ff_ori(ngrid,ngrid,ngrid)
      complex *16 ff_tru(ngrid,ngrid,ngrid)
      complex *16 ff_est(ngrid,ngrid,ngrid)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ffhat_ori(:)
      complex *16, allocatable :: ffhat_tru(:)
      complex *16, allocatable :: ffhat_est(:)
      complex *16, allocatable :: modsph_ori(:)
      complex *16, allocatable :: modsph_tru(:)
      complex *16, allocatable :: modsph_est(:)
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
      real *8, allocatable :: alpha5_tru_(:,:)
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      real *8, allocatable :: alpha5_est_(:,:)
      real *8, allocatable :: alpha_0_all_(:)
      real *8, allocatable :: alpha_1_all_(:)
      real *8, allocatable :: alpha_0_sample_(:)
      real *8, allocatable :: alpha_1_sample_(:)
      real *8, allocatable :: delta_x_sort_(:)
      real *8, allocatable :: delta_y_sort_(:)
      real *8, allocatable :: gamma_z_sort_(:)
      complex *16, allocatable :: C_Z_sort_(:)
      integer *4, allocatable :: I_permute_(:)
      integer *4, allocatable :: I_inverse_(:)
c
      external fgauss

c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
      done=1
      pi=4.0d0*atan(done)
      call prini(6,13)
      
      n_source_use = min(n_source_max,32) ! number of gaussian sources
      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
      do nsource=0,n_source_use-1
         source_r = 0.5d0*(nsource)/max(1,n_source_use-1);
         x0y0z0(1,1+nsource) = source_r*cos(2*pi*2*nsource
     $        /n_source_use)
         x0y0z0(2,1+nsource) = source_r*sin(2*pi*2*nsource
     $        /n_source_use)
         x0y0z0(3,1+nsource) = 2.0d0*source_r - 0.5 
      enddo
      
c$$$      n_source_use = min(n_source_max,128) ! number of gaussian sources
c$$$      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
c$$$      do nsource=0,n_source_use-1
c$$$         source_r = 2.0d0*pi*(nsource)/max(1,n_source_use)
c$$$         x0y0z0(1,1+nsource) = cos(2.0d0*source_r)*(3.0d0+cos(3.0d0
c$$$     $        *source_r))/8.0d0
c$$$         x0y0z0(2,1+nsource) = sin(2.0d0*source_r)*(3.0d0+cos(3.0d0
c$$$     $        *source_r))/8.0d0
c$$$         x0y0z0(3,1+nsource) = 2.0d0*sin(3.0d0*source_r)/8.0d0
c$$$      enddo

c$$$      n_source_use = min(n_source_max,3) ! number of gaussian sources
c$$$      sigma = 1.0d0/16.0d0      ! isotropic std of each gaussian source
c$$$      do nsource=0,n_source_use-1
c$$$         x0y0z0(1,1+nsource) = 0.5-rand();
c$$$         x0y0z0(2,1+nsource) = 0.5-rand();
c$$$         x0y0z0(3,1+nsource) = 0.5-rand();
c$$$      enddo

      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './dir_mda/source_xyz_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,x0y0z0,fname)
c
c     create physical space grid and sample function on it.
      ityper = 0
      itypep = 1
      ifpr = 1
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,x0y0z0,n_source_use,sigma,
     1     fgauss,ff_ori,fnorm)
      if (verbose.gt.1) write(6,*) '  fnorm is ',fnorm
      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
      write(fname,'(A)') './dir_mda/x_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x,fname)
      write(fname,'(A)') './dir_mda/y_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y,fname)
      write(fname,'(A)') './dir_mda/z_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './dir_mda/S_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff_ori,fname)

c-----------------------------------------------------------------------
c     2) Determine size of spherical grid, allocate space and compute
c     Fourier transform at corresponding points in k-space.
c     Also determine size of spherical harmonic model  (nsphstore)
c-----------------------------------------------------------------------
      rmax = 1.0d0*ngridr
      call getgridr(rmax,ngridr,ityper,xnodesr,wtsr)
      nsphstore = 0
      nimage_size = 0
      do i = 1,ngridr
         nlats(i) = nint(pi*xnodesr(i))
         if (nlats(i).lt.6) nlats(i) = 6
         ngridc(i) = nlats(i)*2
         isph_start(i) = nsphstore
         nterms_sph(i) = nint(xnodesr(i) + 2)
         nsphstore = nsphstore + (nterms_sph(i)+1)**2
         nimage_size = nimage_size + ngridc(i)
      enddo
      if (verbose.gt.1) call prinf('  nlats is *',nlats,ngridr)
      if (verbose.gt.1) call prinf('  ngridc is *',ngridc,ngridr)
      if (verbose.gt.1) call prinf('  nimage_size is *',nimage_size,1)
      ntemplate_size = nimage_size
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
      allocate(ffhat_ori(ntot))
      allocate(ffhat_tru(ntot))
      allocate(ffhat_est(ntot))
c
      allocate(modsph_ori(0:nsphstore-1))
      allocate(modsph_tru(0:nsphstore-1))
      allocate(modsph_est(0:nsphstore-1))
c
c     create kspace grid data from physical space model
c
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff_ori,h,ngrid,eps,rmax,ngridr,
     $     ityper,nlats,itypep,ntot,numonsphere, kgridx,kgridy,kgridz
     $     ,wts,ffhat_ori,ffhatnorm,ier)
      if (verbose.gt.1) then
         write(6,*) 'ffhatnorm: ',ffhatnorm
         write(6,*) 'numonsphere: ',(numonsphere(ii),ii=1,ngridr)
      end if
c
c     convert to model: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(ffhat_ori,rmax,nterms_sph, ityper,ngridr
     $     ,nlats,itypep,numonsphere,modsph_ori,nsphstore)

c-----------------------------------------------------------------------
c     3) Assign slices to orientations of spherical grid on highest
c        frequency sphere and define gammas.
c-----------------------------------------------------------------------
      nimages = numonsphere(ngridr)
      ntemplates = numonsphere(ngridr)
      call getspheregrid(nlats(ngridr),itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,nspherebig)
      if (verbose.gt.0) write(6,*) '  nimages is ',nimages
      if (verbose.gt.1) write(6,*) '  nspherebig is ',nspherebig
c      
      allocate(alpha5_tru_(5,nimages))
      allocate(alpha5_est_(5,nimages))
      nnn = 0
      do kk = 1,nlats(ngridr)
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            alpha5_tru_(1,nnn) = dacos(ctheta)
            alpha5_tru_(2,nnn) = phi
            alpha5_tru_(3,nnn) = (1.0d0 - 2*rand())*pi
            alpha5_tru_(4,nnn) = (1.0d0 - 2*rand())/rmax
            alpha5_tru_(5,nnn) = (1.0d0 - 2*rand())/rmax
c$$$            alpha5_tru_(4,nnn) = 0.0d0/rmax
c$$$            alpha5_tru_(5,nnn) = 0.0d0/rmax
            alpha5_est_(1,nnn) = pi*rand()
            alpha5_est_(2,nnn) = 2*rand()*pi
            alpha5_est_(3,nnn) = (1.0d0 - 2*rand())*pi
            alpha5_est_(4,nnn) = 0.0d0/rmax
            alpha5_est_(5,nnn) = 0.0d0/rmax
c$$$            alpha5_est_(1,nnn) = alpha5_tru_(1,nnn)
c$$$            alpha5_est_(2,nnn) = alpha5_tru_(2,nnn)
c$$$            alpha5_est_(3,nnn) = alpha5_tru_(3,nnn)
c$$$            alpha5_est_(4,nnn) = alpha5_tru_(4,nnn)
c$$$            alpha5_est_(5,nnn) = alpha5_tru_(5,nnn)
         enddo
      enddo
      if (verbose.gt.1) write(6,*) '  nnn is ',nnn
      if (nnn.ne.numonsphere(ngridr)) then
         write(6,'(A,I0,A,I0,A,I0,A)') 'Warning, nnn ',nnn,'.ne.'
     $        ,numonsphere,'(',ngridr,')'
      end if
c
      allocate(cslices(nimage_size,nimages))
      allocate(templates(ntemplate_size,nimages))
      allocate(alpha_0_all_(0:ntemplates-1))
      allocate(alpha_1_all_(0:ntemplates-1))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     1           icstart)
      if (verbose.gt.1) write(6,*) '  ntemplatesize is ',ntemplatesize

c-----------------------------------------------------------------------
c     4) compute slices directly from physical space function (ff_ori).
c        with orientation vectors defined by alpha5_tru_.
c-----------------------------------------------------------------------
c
      boxsize = 2.0d0
      call mk_simulated_slices_alpha5(ff_ori,ngrid,boxsize,eps,ngridr,
     1       ityper,ngridc,rmax,alpha5_tru_,nimages,nimage_size,cslices)

c-----------------------------------------------------------------------
c     5) TEST 1: generate templates; should match slices since we took slices
c     from the same grid, so long as gamma_z (third Euler angle)
c     and the  delta_x and delta_y (displacements) are set to zero.
c-----------------------------------------------------------------------
c
c$$$      n_M_sample_max = nimages
      n_M_sample_max = 512
      allocate(I_M_sample_(0:n_M_sample_max-1))
      allocate(M_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(alpha_tru_(0:5*n_M_sample_max-1))
      allocate(alpha_est_(0:5*n_M_sample_max-1))
c$$$      n_S_sample_max = ntemplates
      n_S_sample_max = 256
      allocate(I_S_sample_(0:n_S_sample_max-1))
      allocate(S_sample_(0:ntemplatesize*n_S_sample_max-1))
      allocate(alpha_0_sample_(0:n_S_sample_max-1))
      allocate(alpha_1_sample_(0:n_S_sample_max-1))
      allocate(delta_x_sort_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_sort_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_sort_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_sort_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_permute_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_inverse_(0:n_S_sample_max*n_M_sample_max-1))
      nlow = 3
      do ncur=nlow,ngridr
         if (verbose.gt.0) write(6,*) 'ncur: ',ncur
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         select images to use
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         n_M_sample = min(nimages,n_M_sample_max)
         if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_M_sample
     $        ,' images'
         do nM_sample=0,n_M_sample-1
            if (n_M_sample.lt.nimages) then
               I_M_sample_(nM_sample) = min(nimages-1,floor(nimages
     $              *1.0d0*nM_sample/(n_M_sample-1)))
            else
               I_M_sample_(nM_sample) = nM_sample
            end if
            nm = I_M_sample_(nM_sample)
            call cp1_c16(nimage_size,cslices(1,1+nm)
     $           ,M_sample_(nM_sample*nimage_size))
         enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         copy alpha5 to alphas
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'copying alpha-->alpha5'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(5,alpha5_tru_(1,1+nm),alpha_tru_(0+5
     $           *nM_sample))
            call cp1_r8(5,alpha5_est_(1,1+nm),alpha_est_(0+5
     $           *nM_sample))
         enddo
         d_(0) = 5
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda/alpha_tru_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_tru_,fname)
         d_(0) = 5
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda/alpha_est_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_est_,fname)
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         rebuild model
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'rebuilding modsph_tru'
         oversamp = 2.0d0
         eps = 1.0d-6
         kord = 5
         if (verbose.gt.1) then
            write(6,*) 'nimage_size: ',nimage_size
            write(6,*) 'n_M_sample: ',n_M_sample
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
         call rebuild_model_alpha5(M_sample_,nimage_size,n_M_sample
     $        ,alpha_tru_,icstart,nlats,itypep,xnodesr,ngridc,nlow,ncur
     $        ,nsphstore,isph_start,nterms_sph,oversamp,kord,eps
     $        ,modsph_tru)
         call model_to_kspacegrid(ffhat_tru,rmax,nterms_sph,ityper
     $        ,ngridr,nlats,itypep,numonsphere,modsph_tru,nsphstore)
         do ii = 1,ntot
            ffhat_tru(ii) = ffhat_tru(ii)*wts(ii)
         enddo
         iflag = -1
         npts = ngrid*ngrid*ngrid
         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_tru,iflag
     $        ,eps,npts,x,y,z,ff_tru,ier)
         d_(0) = ngrid
         d_(1) = ngrid
         d_(2) = ngrid
c$$$         write(fname,'(A)') './dir_mda/x_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,x,fname)
c$$$         write(fname,'(A)') './dir_mda/y_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,y,fname)
c$$$         write(fname,'(A)') './dir_mda/z_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,z,fname)
         write(fname,'(A,I0,A)') './dir_mda/S_tru_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_c16(3,d_,ff_tru,fname)
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'rebuilding modsph_est'
         call rebuild_model_alpha5(M_sample_,nimage_size,n_M_sample
     $        ,alpha_est_,icstart,nlats,itypep,xnodesr,ngridc,nlow,ncur
     $        ,nsphstore,isph_start,nterms_sph,oversamp,kord,eps
     $        ,modsph_est)
         call model_to_kspacegrid(ffhat_est,rmax,nterms_sph,ityper
     $        ,ngridr,nlats,itypep,numonsphere,modsph_est,nsphstore)
         do ii = 1,ntot
            ffhat_est(ii) = ffhat_est(ii)*wts(ii)
         enddo
         iflag = -1
         npts = ngrid*ngrid*ngrid
         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_est,iflag
     $        ,eps,npts,x,y,z,ff_est,ier)
         d_(0) = ngrid
         d_(1) = ngrid
         d_(2) = ngrid
c$$$         write(fname,'(A)') './dir_mda/x_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,x,fname)
c$$$         write(fname,'(A)') './dir_mda/y_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,y,fname)
c$$$         write(fname,'(A)') './dir_mda/z_pos_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$         call MDA_write_r8(3,d_,z,fname)
         write(fname,'(A,I0,A)') './dir_mda/S_est_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_c16(3,d_,ff_est,fname)
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         generate templates
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         numonsphere_cur = numonsphere(ncur)
         if (verbose.gt.0) write(6,'(A,I0,A)')
     $        'generating numonsphere_cur = ',numonsphere_cur
     $        ,' templates'
         call template_gen(modsph_est,nsphstore,isph_start, nterms_sph
     $        ,ngridr,ngridc,ntemplatesize,icstart, nlats,itypep,ncur
     $        ,numonsphere_cur,templates)
         call getspheregrid(nlats(ncur),itypep,xnodesth_tmp, sthetas_tmp
     $        ,wtsth_tmp,ngridps_tmp,phsteps_tmp,nspherebig_tmp)
         nnn = 0
         do kk = 1,nlats(ncur)
            ctheta = xnodesth_tmp(kk)
            stheta = sthetas_tmp(kk)
            phstep = phsteps_tmp(kk)
            do ll = 1,ngridps_tmp(kk)
               phi = (ll-1)*phstep
               alpha_0_all_(nnn) = dacos(ctheta)
               alpha_1_all_(nnn) = phi
               nnn = nnn+1
            enddo
         enddo
         if (verbose.gt.1) write(6,*) '  nnn is ',nnn
         if (nnn.ne.numonsphere_cur) then
            write(6,'(A,I0,A,I0)') 'Warning, nnn ',nnn,'.ne.'
     $           ,numonsphere(ncur)
         end if
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         select templates to use
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         n_S_sample = min(numonsphere_cur,n_S_sample_max)
         if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_S_sample
     $        ,' templates'
         do nS_sample=0,n_S_sample-1
            if (I_S_sample.lt.numonsphere_cur) then
               I_S_sample_(nS_sample) = min(numonsphere_cur-1
     $              ,floor(numonsphere_cur*1.0d0*nS_sample/(n_S_sample
     $              -1)))
            else
               I_S_sample_(nS_sample) = nS_sample
            end if
            ns = I_S_sample_(nS_sample)
            call cp1_c16(ntemplatesize,templates(1,1+ns)
     $           ,S_sample_(nS_sample*ntemplatesize))
            alpha_0_sample_(nS_sample) = alpha_0_all_(ns)
            alpha_1_sample_(nS_sample) = alpha_1_all_(ns)
         enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         calculate innerproducts
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $        ,n_M_sample,')-x-(n_S = ',n_S_sample, ') = (', n_M_sample
     $        *n_S_sample,') innerproducts'
         call test_innerproduct_batch_wrapper1(0,ncur,nlats
     $        ,xnodesr,max_x_c,n_S_sample,ntemplatesize,S_sample_
     $        ,n_M_sample,nimage_size,M_sample_,alpha_est_,0.1d0,0.5d0
     $        ,5,5,max(32,2*ncur),-1,delta_x_sort_,delta_y_sort_
     $        ,gamma_z_sort_,C_Z_sort_,I_permute_,I_inverse_)
c$$$         call test_innerproduct_batch_wrapper1(0,ncur,nlats
c$$$     $        ,xnodesr,max_x_c,n_S_sample,ntemplatesize,S_sample_
c$$$     $        ,n_M_sample,nimage_size,M_sample_,alpha_est_,0.1d0,1.0d0
c$$$     $        ,1,1,max(32,2*ncur),-1,delta_x_sort_,delta_y_sort_
c$$$     $        ,gamma_z_sort_,C_Z_sort_,I_permute_,I_inverse_)
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         update alphas
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'update alphas'
         do nm=0,n_M_sample-1
            ns = I_permute_(n_S_sample-1 + nm*n_S_sample)
            alpha_est_(0 + nm*5) = alpha_0_sample_(ns)
            alpha_est_(1 + nm*5) = alpha_1_sample_(ns)
            alpha_est_(2 + nm*5) = alpha_est_(2 + nm*5) +
     $           gamma_z_sort_(n_S_sample-1 + nm*n_S_sample)
            alpha_est_(3 + nm*5) = alpha_est_(3 + nm*5) +
     $           delta_x_sort_(n_S_sample-1 + nm*n_S_sample)
            alpha_est_(4 + nm*5) = alpha_est_(4 + nm*5) +
     $           delta_y_sort_(n_S_sample-1 + nm*n_S_sample)
         enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         copy alpha to alpha5
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'copying alpha-->alpha5'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(5,alpha_est_(0+5*nM_sample),alpha5_est_(1,1
     $           +nm))
         enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         end loop
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         
      enddo

 10   continue
      stop
      end
