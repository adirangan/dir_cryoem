C
c      Testing code for rebuilding model from slices.
c$$$      Adding ctf and image+template normalization
c
      implicit real *8 (a-h,o-z)
c     
      include 'omp_lib.h'
      integer verbose
      data verbose / 1 /
      parameter (ngridr=20)
      integer ncur
      parameter (ngrid=60)
      parameter (ntmax=4*ngridr)
      integer n_source_max,n_source_use,nsource
      parameter (n_source_max=128)
      real *8 sigma,source_r
      integer, allocatable :: nlats(:)
      integer, allocatable :: numonsphere(:)
      integer, allocatable :: nterms_sph(:)
      integer, allocatable :: isph_start(:)
      integer, allocatable :: ngridc(:)
      integer, allocatable :: ngridc2(:)
      integer, allocatable :: icstart(:)
      integer, allocatable :: ngridps(:)
      integer, allocatable :: ngridps_tmp(:)
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      real *8, allocatable :: xnodesth(:)
      real *8, allocatable :: sthetas(:)
      real *8, allocatable :: wtsth(:)
      real *8, allocatable :: phsteps(:)
      real *8, allocatable :: xnodesr_tmp(:)
      real *8, allocatable :: wtsr_tmp(:)
      real *8, allocatable :: xnodesth_tmp(:)
      real *8, allocatable :: sthetas_tmp(:)
      real *8, allocatable :: wtsth_tmp(:)
      real *8, allocatable :: phsteps_tmp(:)
      real *8, allocatable :: x0y0z0(:,:)
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      complex *16, allocatable :: ff_ori(:,:,:)
      complex *16, allocatable :: ff_tru(:,:,:)
      complex *16, allocatable :: ff_est(:,:,:)
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
      logical displacement_flag
      parameter (displacement_flag=.true.)
c$$$      logical normalization_flag
c$$$      parameter (normalization_flag=.true.)
      integer n_delta_x,n_delta_y
      real *8 N_pixels_in
      real *8 displacement_max
      integer n_M_sample,nM_sample,n_S_sample,nS_sample
      integer *4, allocatable :: I_M_sample_(:)
      integer *4, allocatable :: I_S_sample_(:)
      complex *16, allocatable :: M_sample_(:)
c$$$      complex *16, allocatable :: N_sample_(:)
      complex *16, allocatable :: S_sample_(:)
      integer *4 n_alpha
c$$$      parameter (n_alpha=7)
      parameter (n_alpha=5)
      include 'nalpha_define.f'
c$$$      real *8 l2_norm
      real *8, allocatable :: alpha2d_tru_(:,:)
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      real *8, allocatable :: alpha2d_est_(:,:)
c$$$      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: alpha_polar_a_all_(:)
      real *8, allocatable :: alpha_azimu_b_all_(:)
      real *8, allocatable :: alpha_polar_a_sample_(:)
      real *8, allocatable :: alpha_azimu_b_sample_(:)
      real *8, allocatable :: delta_x_max_(:)
      real *8, allocatable :: delta_y_max_(:)
      real *8, allocatable :: gamma_z_max_(:)
c$$$      complex *16, allocatable :: C_S_max_(:)
      complex *16, allocatable :: C_Z_max_(:)
      real *8, allocatable :: delta_x_sort_SM_(:)
      real *8, allocatable :: delta_y_sort_SM_(:)
      real *8, allocatable :: gamma_z_sort_SM_(:)
c$$$      complex *16, allocatable :: C_S_sort_SM_(:)
      complex *16, allocatable :: C_Z_sort_SM_(:)
      integer *4, allocatable :: I_permute_SM_(:)
      integer *4, allocatable :: I_inverse_SM_(:)
      real *8, allocatable :: delta_x_sort_MS_(:)
      real *8, allocatable :: delta_y_sort_MS_(:)
      real *8, allocatable :: gamma_z_sort_MS_(:)
c$$$      complex *16, allocatable :: C_S_sort_MS_(:)
      complex *16, allocatable :: C_Z_sort_MS_(:)
      integer *4, allocatable :: I_permute_MS_(:)
      integer *4, allocatable :: I_inverse_MS_(:)
      real *8 timing_tic,timing_toc
      external fgauss
c
c$$$      integer *4 n_ctf,nctf
c$$$      integer *4, allocatable :: I_ctf_(:)
c$$$      complex *16, allocatable :: CTF_k_p_(:)
c$$$      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_ver5]'
      end if

c-----------------------------------------------------------------------
c     0) allocate some arrays (dynamic allocation required for openmp)
c-----------------------------------------------------------------------
      allocate(nlats(ngridr))
      allocate(numonsphere(ngridr))
      allocate(nterms_sph(ngridr))
      allocate(isph_start(ngridr))
      allocate(ngridc(ntmax))
      allocate(ngridc2(ntmax))
      allocate(icstart(ngridr))
      allocate(ngridps(ntmax))
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(xnodesth(ntmax))
      allocate(sthetas(ntmax))
      allocate(wtsth(ntmax))
      allocate(phsteps(ntmax))
      allocate(ngridps_tmp(ntmax))
      allocate(xnodesth_tmp(ntmax))
      allocate(sthetas_tmp(ntmax))
      allocate(wtsth_tmp(ntmax))
      allocate(phsteps_tmp(ntmax))
      allocate(x0y0z0(3,n_source_max))
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
      allocate(ff_ori(ngrid,ngrid,ngrid))
      allocate(ff_tru(ngrid,ngrid,ngrid))
      allocate(ff_est(ngrid,ngrid,ngrid))

c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
      done=1
      pi=4.0d0*atan(done)
      call prini(6,13)
      
      n_source_use = min(n_source_max,32) ! number of gaussian sources
      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
      do nsource=0,n_source_use-1
         source_r = 0.5d0*(nsource)/max(1,n_source_use-1)
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
c$$$         x0y0z0(1,1+nsource) = 0.5-rand()
c$$$         x0y0z0(2,1+nsource) = 0.5-rand()
c$$$         x0y0z0(3,1+nsource) = 0.5-rand()
c$$$      enddo

      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './dir_mda/source_xyz_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,x0y0z0,fname)

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
c        frequency sphere and define alphas.
c-----------------------------------------------------------------------
      nimages = numonsphere(ngridr)
      ntemplates = numonsphere(ngridr)
      call getspheregrid(nlats(ngridr),itypep,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,nspherebig)
      if (verbose.gt.0) write(6,*) '  nimages is ',nimages
      if (verbose.gt.1) write(6,*) '  nspherebig is ',nspherebig
c      
      if (displacement_flag) then
         displacement_max = 1.0d0/rmax
      else
         displacement_max = 0.0d0/rmax
      end if
      allocate(alpha2d_tru_(n_alpha,nimages))
      allocate(alpha2d_est_(n_alpha,nimages))
      nnn = 0
      do kk = 1,nlats(ngridr)
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            alpha2d_tru_(1+nalpha_polar_a,nnn) = dacos(ctheta)
            alpha2d_tru_(1+nalpha_azimu_b,nnn) = phi
            alpha2d_tru_(1+nalpha_gamma_z,nnn) = (1.0d0 - 2*rand())*pi
            if (displacement_flag) then
               angle_tmp = 2*pi*rand()
               delta_tmp = rand()*displacement_max
               alpha2d_tru_(1+nalpha_delta_x,nnn) = cos(angle_tmp)
     $              *delta_tmp
               alpha2d_tru_(1+nalpha_delta_y,nnn) = sin(angle_tmp)
     $              *delta_tmp
            else
               alpha2d_tru_(1+nalpha_delta_x,nnn) = 0.0d0/rmax
               alpha2d_tru_(1+nalpha_delta_y,nnn) = 0.0d0/rmax
            end if !displacement_flag
c$$$            if (normalization_flag) then
c$$$c$$$               draw v from (1/tau)exp(-t/tau) via: v = -tau*log(1-rand())
c$$$c$$$               alpha2d_tru_(1+nalpha_l2_norm,nnn) = - 1.0d0 * log(1.0d0 - rand())
c$$$c$$$               draw v from [0.5,1.5]
c$$$               alpha2d_tru_(1+nalpha_l2_norm,nnn) = 0.5d0 + rand()
c$$$            else
c$$$               alpha2d_tru_(1+nalpha_l2_norm,nnn) = 1.0d0
c$$$            end if !normalization_flag
c$$$c$$$            to start with somewhat arbitrary alpha2d_est_
c$$$            alpha2d_tru_(1+nalpha_ctf_ind,nnn) = 0.0d0
            alpha2d_est_(1+nalpha_polar_a,nnn) = pi*rand()
            alpha2d_est_(1+nalpha_azimu_b,nnn) = 2*rand()*pi
            alpha2d_est_(1+nalpha_gamma_z,nnn) = (1.0d0 - 2*rand())*pi
            alpha2d_est_(1+nalpha_delta_x,nnn) = 0.0d0/rmax
            alpha2d_est_(1+nalpha_delta_y,nnn) = 0.0d0/rmax
c$$$            alpha2d_est_(1+nalpha_l2_norm,nnn) = 1.0d0
c$$$            alpha2d_est_(1+nalpha_ctf_ind,nnn) = 0.0d0
c$$$            to start with alpha2d_est_ .eq. alpha2d_tru_
c$$$            call cp1_r8(n_alpha,alpha2d_tru_(1,1+nm),alpha2d_est_(1,1+nm))
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
      allocate(alpha_polar_a_all_(0:ntemplates-1))
      allocate(alpha_azimu_b_all_(0:ntemplates-1))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     $     icstart)
      if (verbose.gt.1) write(6,*) '  ntemplatesize is ',ntemplatesize

c-----------------------------------------------------------------------
c     4) compute slices directly from physical space function (ff_ori).
c        with orientation vectors defined by alpha2d_tru_.
c-----------------------------------------------------------------------
c
      boxsize = 2.0d0
c$$$      call mk_simulated_slices_alpha2d(ff_ori,ngrid,boxsize,eps,ngridr,
c$$$     $     ityper,ngridc,rmax,n_alpha,alpha2d_tru_,nimages,nimage_size
c$$$     $     ,cslices)
      call mk_simulated_slices_alpha5(ff_ori,ngrid,boxsize,eps,ngridr,
     $     ityper,ngridc,rmax,alpha2d_tru_,nimages,nimage_size,cslices)

c-----------------------------------------------------------------------
c     5) generate ctf functions to associate with slices.
c-----------------------------------------------------------------------
c

c$$$c$$$      Number of different ctfs used is n_ctf.
c$$$c$$$      We expect this to be proportional to the number of images,
c$$$c$$$      with batches of roughly 30-300 images each having the same ctf.
c$$$      n_ctf = 3
c$$$      allocate(CTF_k_p_(0:n_ctf*nimage_size-1))
c$$$c$$$      Here we use an unrealistic ctf-function with a high degree
c$$$c$$$      of anisotropy to test our code.
c$$$      do nctf=0,n_ctf-1
c$$$         ctf_p_0 = nctf * pi/max(1,n_ctf)
c$$$         ctf_p_1 = 0.125d0 + 0.0625*nctf
c$$$         ctf_p_2 = 10.0d0
c$$$         ctf_p_3 = 10.0d0
c$$$         ctf_p_4 = 5.0d0
c$$$         call get_ctf_k_p_(ngridr,CTF_k_p_(nctf*nimage_size),ctf_p_0
c$$$     $        ,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4)
c$$$      enddo
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         print out the ctf-functions
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$      write(fname,'(A)') './dir_jpg/Fig_ctf'
c$$$      call Fig_gen_ctf_ver0(ngridr,nlats,xnodesr,n_ctf,nimage_size
c$$$     $     ,CTF_k_p_,fname)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         Randomly determine which images correspond to which ctf-function.
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$      allocate(I_ctf_(0:nimages-1))
c$$$      do nm=0,nimages-1
c$$$         I_ctf_(nm) = mod(nm,n_ctf)
c$$$         alpha2d_tru_(1+nalpha_ctf_ind,nm) = I_ctf_(nm)
c$$$         alpha2d_est_(1+nalpha_ctf_ind,nm) = I_ctf_(nm)
c$$$      enddo

c-----------------------------------------------------------------------
c     6) Run the main loop (i.e., try and recover molecule from images).
c-----------------------------------------------------------------------
c
c$$$      n_M_sample_max = maximum number of images to use at each step.
      n_M_sample_max = 64
      allocate(I_M_sample_(0:n_M_sample_max-1))
      allocate(M_sample_(0:nimage_size*n_M_sample_max-1))
c$$$      allocate(N_sample_(0:nimage_size*n_M_sample_max-1))
c$$$      allocate(C_M_(0:n_M_sample_max-1))
      allocate(alpha_tru_(0:n_alpha*n_M_sample_max-1))
      allocate(alpha_est_(0:n_alpha*n_M_sample_max-1))
c$$$      n_S_sample_max = maximum number of templates to use at each step.
      n_S_sample_max = 48
      allocate(I_S_sample_(0:n_S_sample_max-1))
      allocate(S_sample_(0:ntemplatesize*n_S_sample_max-1))
      allocate(alpha_polar_a_sample_(0:n_S_sample_max-1))
      allocate(alpha_azimu_b_sample_(0:n_S_sample_max-1))
      allocate(delta_x_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_max_(0:n_S_sample_max*n_M_sample_max-1))
c$$$      allocate(C_S_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_x_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
c$$$      allocate(C_S_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_permute_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_inverse_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_x_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
c$$$      allocate(C_S_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_permute_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_inverse_MS_(0:n_S_sample_max*n_M_sample_max-1))
      if (displacement_flag) then
         n_delta_x = 5
         n_delta_y = 5
         N_pixels_in = 0.5d0
      else
         n_delta_x = 1
         n_delta_y = 1
         N_pixels_in = 0.5d0
      end if
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
         enddo ! do nM_sample=0,n_M_sample-1

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         copy alpha2d to alphas
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'copying alpha2d-->alpha'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha2d_tru_(1,1+nm),alpha_tru_(0
     $           +n_alpha*nM_sample))
            call cp1_r8(n_alpha,alpha2d_est_(1,1+nm),alpha_est_(0
     $           +n_alpha*nM_sample))
         enddo
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         write alphas to disk for postmortem analysis.
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         d_(0) = n_alpha
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda/alpha_tru_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_tru_,fname)
         d_(0) = n_alpha
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda/alpha_est_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_est_,fname)
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         apply ctf-convolution to selected images 
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         do nM_sample=0,n_M_sample-1
c$$$            nm = I_M_sample_(nM_sample)
c$$$            nctf = I_ctf_(nm)
c$$$            call xx1_c16(nimage_size,cslices(1,1+nm),CTF_k_p_(nctf
c$$$     $           *nimage_size),N_sample_(nM_sample*nimage_size))
c$$$         enddo ! do nM_sample=0,n_M_sample-1

c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         multiply selected images by true normalization factor 
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         do nM_sample=0,n_M_sample-1
c$$$            nm = I_M_sample_(nM_sample)
c$$$            l2_norm = alpha_tru_(nalpha_l2_norm+n_alpha*nM_sample)
c$$$            call af1_c16(nimage_size,1.0d0*cmplx(l2_norm,0.0d0),1.0d0
c$$$     $           *cmplx(0.0d0,0.0d0),N_sample_(nM_sample*nimage_size)
c$$$     $           ,N_sample_(nM_sample*nimage_size))
c$$$         enddo ! do nM_sample=0,n_M_sample-1

c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         print out a subset of image+ctf pairs
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         write(fname,'(A,I0)') './dir_jpg/Fig_MNC_',ncur
c$$$         call Fig_gen_ctf_ver1(ncur,nlats,xnodesr,n_M_sample,nimage_size
c$$$     $        ,M_sample_,N_sample_,nimage_size,CTF_k_p_,n_alpha
c$$$     $        ,alpha_est_,min(16,n_M_sample),fname)

c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         rebuild model based on current alpha.
c$$$         We rebuild two models.
c$$$         The first, modsph_tru is the model obtained by using the 
c$$$         true alpha-parameters for each image.
c$$$         The second, modsph_est is the model obtained by using the
c$$$         best-guess alpha-parameters for each image (which depends
c$$$         on previous iterations).
c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         First find modsph_tru
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
         timing_tic = second()
c$$$         call getspheregrid(nlats(ncur),itypep,xnodesth_tmp, sthetas_tmp
c$$$     $        ,wtsth_tmp,ngridps_tmp,phsteps_tmp,nspherebig_tmp)
c$$$         call rebuild_model_alpha2d(N_sample_,nimage_size,n_M_sample
c$$$     $        ,n_CTF,nimage_size,CTF_k_p_,I_ctf_,n_alpha,alpha_tru_
c$$$     $        ,icstart ,nlats,itypep,xnodesr,ngridc,nlow,ncur
c$$$     $        ,xnodesth_tmp,wtsth_tmp,nsphstore ,isph_start,nterms_sph
c$$$     $        ,oversamp,kord,modsph_tru)

c$$$         call rebuild_model_alphaX(M_sample_,nimage_size,n_M_sample
c$$$     $        ,n_alpha,alpha_tru_ ,icstart ,nlats,itypep,xnodesr,ngridc
c$$$     $        ,nlow,ncur ,nsphstore ,isph_start,nterms_sph ,oversamp
c$$$     $        ,kord,eps,modsph_tru)
         call rebuild_model_alphaX(M_sample_,nimage_size,n_M_sample ,5
     $        ,alpha_tru_,icstart,nlats,itypep,xnodesr,ngridc,nlow,ncur
     $        ,nsphstore,isph_start,nterms_sph,oversamp,kord,eps
     $        ,modsph_tru)

         timing_toc = second()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'rebuild_model_alpha2d (tru):'
     $           ,' total_time ',timing_toc-timing_tic
         end if

         goto 10

c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         write model to disk for postmortem analysis
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         if (verbose.gt.0) write(6,*) 'modsph_tru to kspacegrid'
c$$$         call model_to_kspacegrid(ffhat_tru,rmax,nterms_sph,ityper
c$$$     $        ,ngridr,nlats,itypep,numonsphere,modsph_tru,nsphstore)
c$$$         do ii = 1,ntot
c$$$            ffhat_tru(ii) = ffhat_tru(ii)*wts(ii)
c$$$         enddo
c$$$         iflag = -1
c$$$         npts = ngrid*ngrid*ngrid
c$$$         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_tru,iflag
c$$$     $        ,eps,npts,x,y,z,ff_tru,ier)
c$$$         d_(0) = ngrid
c$$$         d_(1) = ngrid
c$$$         d_(2) = ngrid
c$$$         write(fname,'(A,I0,A)') './dir_mda/S_tru_',ncur,'_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
c$$$         call MDA_write_c16(3,d_,ff_tru,fname)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         Second find modsph_est
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         if (verbose.gt.0) write(6,*) 'rebuilding modsph_est'
c$$$         timing_tic = second()
c$$$c$$$         call getspheregrid(nlats(ncur),itypep,xnodesth_tmp, sthetas_tmp
c$$$c$$$     $        ,wtsth_tmp,ngridps_tmp,phsteps_tmp,nspherebig_tmp)
c$$$c$$$         call rebuild_model_alpha2d(N_sample_,nimage_size,n_M_sample
c$$$c$$$     $        ,n_CTF,nimage_size,CTF_k_p_,I_ctf_,n_alpha,alpha_est_
c$$$c$$$     $        ,icstart ,nlats,itypep,xnodesr,ngridc,nlow,ncur
c$$$c$$$     $        ,xnodesth_tmp,wtsth_tmp,nsphstore ,isph_start,nterms_sph
c$$$c$$$     $        ,oversamp,kord,modsph_est)
c$$$
c$$$         call rebuild_model_alphaX(N_sample_,nimage_size,n_M_sample
c$$$     $        ,n_alpha,alpha_est_ ,icstart ,nlats,itypep,xnodesr,ngridc
c$$$     $        ,nlow,ncur ,nsphstore ,isph_start,nterms_sph ,oversamp
c$$$     $        ,kord,eps,modsph_est)
c$$$
c$$$         timing_toc = second()
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,A,F8.3)') 'rebuild_model_alpha2d (est):'
c$$$     $           ,' total_time ',timing_toc-timing_tic
c$$$         end if
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         write model to disk for postmortem analysis
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         if (verbose.gt.0) write(6,*) 'modsph_est to kspacegrid'
c$$$         call model_to_kspacegrid(ffhat_est,rmax,nterms_sph,ityper
c$$$     $        ,ngridr,nlats,itypep,numonsphere,modsph_est,nsphstore)
c$$$         do ii = 1,ntot
c$$$            ffhat_est(ii) = ffhat_est(ii)*wts(ii)
c$$$         enddo
c$$$         iflag = -1
c$$$         npts = ngrid*ngrid*ngrid
c$$$         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_est,iflag
c$$$     $        ,eps,npts,x,y,z,ff_est,ier)
c$$$         d_(0) = ngrid
c$$$         d_(1) = ngrid
c$$$         d_(2) = ngrid
c$$$         write(fname,'(A,I0,A)') './dir_mda/S_est_',ncur,'_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
c$$$         call MDA_write_c16(3,d_,ff_est,fname)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         generate templates based on current model (using modsph_est)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         numonsphere_cur = numonsphere(ncur)
c$$$         if (verbose.gt.0) write(6,'(A,I0,A)')
c$$$     $        'generating numonsphere_cur = ',numonsphere_cur
c$$$     $        ,' templates'
c$$$         timing_tic = second()
c$$$         call template_gen(modsph_est,nsphstore,isph_start, nterms_sph
c$$$     $        ,ngridr,ngridc,ntemplatesize,icstart, nlats,itypep,ncur
c$$$     $        ,numonsphere_cur,templates)
c$$$         timing_toc = second()
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,A,F8.3)') 'template_gen:' ,' total_time '
c$$$     $           ,timing_toc-timing_tic
c$$$         end if
c$$$         call getspheregrid(nlats(ncur),itypep,xnodesth_tmp, sthetas_tmp
c$$$     $        ,wtsth_tmp,ngridps_tmp,phsteps_tmp,nspherebig_tmp)
c$$$         nnn = 0
c$$$         do kk = 1,nlats(ncur)
c$$$            ctheta = xnodesth_tmp(kk)
c$$$            stheta = sthetas_tmp(kk)
c$$$            phstep = phsteps_tmp(kk)
c$$$            do ll = 1,ngridps_tmp(kk)
c$$$               phi = (ll-1)*phstep
c$$$               alpha_polar_a_all_(nnn) = dacos(ctheta)
c$$$               alpha_azimu_b_all_(nnn) = phi
c$$$               nnn = nnn+1
c$$$            enddo
c$$$         enddo
c$$$         if (verbose.gt.1) write(6,*) '  nnn is ',nnn
c$$$         if (nnn.ne.numonsphere_cur) then
c$$$            write(6,'(A,I0,A,I0)') 'Warning, nnn ',nnn,'.ne.'
c$$$     $           ,numonsphere(ncur)
c$$$         end if
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         select templates to use
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         n_S_sample = min(numonsphere_cur,n_S_sample_max)
c$$$         if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_S_sample
c$$$     $        ,' templates'
c$$$         do nS_sample=0,n_S_sample-1
c$$$            if (I_S_sample.lt.numonsphere_cur) then
c$$$               I_S_sample_(nS_sample) = min(numonsphere_cur-1
c$$$     $              ,floor(numonsphere_cur*1.0d0*nS_sample/(n_S_sample
c$$$     $              -1)))
c$$$            else
c$$$               I_S_sample_(nS_sample) = nS_sample
c$$$            end if
c$$$            ns = I_S_sample_(nS_sample)
c$$$            call cp1_c16(ntemplatesize,templates(1,1+ns)
c$$$     $           ,S_sample_(nS_sample*ntemplatesize))
c$$$            alpha_polar_a_sample_(nS_sample) = alpha_polar_a_all_(ns)
c$$$            alpha_azimu_b_sample_(nS_sample) = alpha_azimu_b_all_(ns)
c$$$         enddo
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         calculate innerproducts
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
c$$$     $        ,n_M_sample,')-x-(n_S = ',n_S_sample, ') = (', n_M_sample
c$$$     $        *n_S_sample,') innerproducts'
c$$$         timing_tic = omp_get_wtime()
c$$$         call test_innerproduct_batch_wrapper_5(0,ncur,nlats,xnodesr
c$$$     $        ,max_x_c,n_S_sample,ntemplatesize,S_sample_,n_M_sample
c$$$     $        ,nimage_size,N_sample_,nimage_size,CTF_k_p_,n_alpha
c$$$     $        ,alpha_est_,0.1d0,N_pixels_in ,displacement_max,n_delta_x
c$$$     $        ,n_delta_y,max(32,2 *ncur), -1 ,delta_x_max_,delta_y_max_
c$$$     $        ,gamma_z_max_,C_M_,C_S_max_,C_Z_max_)
c$$$         timing_toc = omp_get_wtime()
c$$$         if (verbose.gt.0) then
c$$$            write(6,'(A,A,F8.3)') 'test_innerproduct_batch_wrapper_5:'
c$$$     $           ,' total_time ',timing_toc-timing_tic
c$$$         end if
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         write innerproducts to disk for postmortem analysis
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         d_(0) = n_S_sample
c$$$         d_(1) = n_M_sample
c$$$         write(fname,'(A,I0,A)') './dir_mda/C_Z_max_',ncur,'_.mda'
c$$$         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
c$$$         call MDA_write_c16(2,d_,C_Z_max_,fname)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         sort results
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         call test_innerproduct_batch_wrapper_3a(0,n_S_sample,n_M_sample
c$$$     $        ,delta_x_max_,delta_y_max_,gamma_z_max_,C_S_max_,C_Z_max_
c$$$     $        ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
c$$$     $        ,C_S_sort_SM_ ,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
c$$$     $        ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
c$$$     $        ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         print out a subset of image-template pairs
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         write(fname,'(A,I0)') './dir_jpg/Fig_MTS_',ncur
c$$$         call Fig_gen_ver1a(ncur,nlats,xnodesr,n_S_sample,ntemplatesize
c$$$     $        ,S_sample_,n_M_sample,nimage_size,M_sample_,N_sample_
c$$$     $        ,nimage_size,CTF_k_p_,n_alpha,alpha_est_ ,delta_x_sort_SM_
c$$$     $        ,delta_y_sort_SM_,gamma_z_sort_SM_ ,I_permute_SM_,min(16
c$$$     $        ,n_M_sample),fname)
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         update alphas
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         call test_innerproduct_batch_wrapper_4a(0,n_S_sample,n_M_sample
c$$$     $        ,C_M_,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
c$$$     $        ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
c$$$     $        ,delta_x_sort_MS_,delta_y_sort_MS_,gamma_z_sort_MS_
c$$$     $        ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_,I_inverse_MS_
c$$$     $        ,alpha_polar_a_sample_ ,alpha_azimu_b_sample_,n_alpha
c$$$     $        ,alpha_est_) 
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         copy alpha to alpha2d
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$         if (verbose.gt.0) write(6,*) 'copying alpha-->alpha2d'
c$$$         do nM_sample=0,n_M_sample-1
c$$$            nm = I_M_sample_(nM_sample)
c$$$            call cp1_r8(n_alpha,alpha_est_(0+n_alpha*nM_sample)
c$$$     $           ,alpha2d_est_(1,1+nm))
c$$$         enddo
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$c$$$         end loop
c$$$c$$$         /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         
      enddo

 10   continue
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver5]'
      end if
      stop
      end
