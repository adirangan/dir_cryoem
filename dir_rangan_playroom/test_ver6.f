C
c      Testing code for rebuilding model from slices.
c$$$      Adding ctf and image+template normalization
c
      implicit none
      integer ngridr,ngrid,ntmax,i,ier,iflag,ifpr
      integer itypep,ityper,kk,kord,ll,nnn,ii
      integer nimage_size,nimages,nlow,nbot,nwun,n_loading
      parameter (n_loading=1)
      integer npts,nspherebig,nspherebig_tmp,nsphstore,ntemplate_size
      integer ntemplates,ntemplatesize,ntot,numonsphere_cur
      real *8 a,angle_tmp,boxsize,ctheta,delta_tmp,eps,ffhatnorm
      real *8 fnorm,h,phi,phstep,pi,rmax,stheta,oversamp
c     
      include 'omp_lib.h'
      integer verbose
      data verbose / 1 /
      parameter (ngridr=20)
      integer ncur
      parameter (ngrid=60)
      parameter (ntmax=4*ngridr)
      integer n_source_max,n_source_use,nsource,n_cut
      parameter (n_source_max=128)
      real *8 sigma,source_r,source_w,source_z
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
      real *8, allocatable :: xyz_A(:,:)
      real *8, allocatable :: xyz_B(:,:)
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      complex *16, allocatable :: ff_A_ori(:,:,:)
      complex *16, allocatable :: ff_B_ori(:,:,:)
      complex *16, allocatable :: ff_tru(:,:,:)
      complex *16, allocatable :: ff_est(:,:,:)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ffhat_A_ori(:)
      complex *16, allocatable :: ffhat_B_ori(:)
      complex *16, allocatable :: ffhat_tru(:)
      complex *16, allocatable :: ffhat_est(:)
      complex *16, allocatable :: modsph_A_ori(:)
      complex *16, allocatable :: modsph_B_ori(:)
      complex *16, allocatable :: modsph_tru(:)
      complex *16, allocatable :: modsph_est(:)
      complex *16, allocatable :: cslices_A(:,:)
      complex *16, allocatable :: cslices_B(:,:)
      complex *16, allocatable :: cslices(:,:)
      complex *16, allocatable :: templates(:,:)
      integer d_(0:5)
      character(len=1024) fname,format_string
      real *8 max_x_c
      parameter(max_x_c=1.0d0)
      logical displacement_flag
      parameter (displacement_flag=.false.)
      logical normalization_flag
      parameter (normalization_flag=.false.)
      logical ctf_flag
      parameter (ctf_flag=.false.)
      integer svd_calculation_type
      parameter (svd_calculation_type=-1)
      integer n_delta_x,n_delta_y,n_gamma_z
      logical tesselation_flag
      parameter (tesselation_flag=.false.)
      real *8 distance_req
      integer *4 S_used_all
      real *8 N_pixels_in
      real *8 displacement_max
      integer n_M_sample_max,n_M_sample,nM_sample,nm
      integer n_M_A_sample,n_M_B_sample
      integer n_M_A_p_sample,n_M_A_n_sample
      integer n_M_B_p_sample,n_M_B_n_sample
      integer n_modelsph
      integer n_S_sample_max,n_S_sample,nS_sample,ns
      integer *4, allocatable :: I_M_sample_(:)
      integer *4, allocatable :: I_S_sample_(:)
      complex *16, allocatable :: M_sample_(:)
      complex *16, allocatable :: N_sample_(:)
      complex *16, allocatable :: P_sample_(:) ! ctransf
      complex *16, allocatable :: Q_sample_(:) ! model slices
      complex *16, allocatable :: R_sample_(:) ! residuals = images - slices
      complex *16, allocatable :: L_sample_(:) ! loading vectors
      complex *16, allocatable :: S_sample_(:)
      integer *4 n_alpha
      parameter (n_alpha=7)
      include 'nalpha_define.f'
      real *8 l2_norm
      real *8, allocatable :: alpha2d_tru_(:,:)
      real *8, allocatable :: alpha_tru_(:)
      real *8, allocatable :: alpha_est_(:)
      real *8, allocatable :: alpha2d_est_(:,:)
      complex *16, allocatable :: C_M_(:)
      real *8, allocatable :: S_alpha_polar_a_all_(:)
      real *8, allocatable :: S_alpha_azimu_b_all_(:)
      real *8, allocatable :: S_alpha_polar_a_sample_(:)
      real *8, allocatable :: S_alpha_azimu_b_sample_(:)
      real *8, allocatable :: delta_x_max_(:)
      real *8, allocatable :: delta_y_max_(:)
      real *8, allocatable :: gamma_z_max_(:)
      complex *16, allocatable :: C_S_max_(:)
      complex *16, allocatable :: C_Z_max_(:)
      real *8, allocatable :: delta_x_sort_SM_(:)
      real *8, allocatable :: delta_y_sort_SM_(:)
      real *8, allocatable :: gamma_z_sort_SM_(:)
      complex *16, allocatable :: C_S_sort_SM_(:)
      complex *16, allocatable :: C_Z_sort_SM_(:)
      integer *4, allocatable :: I_permute_SM_(:)
      integer *4, allocatable :: I_inverse_SM_(:)
      real *8, allocatable :: delta_x_sort_MS_(:)
      real *8, allocatable :: delta_y_sort_MS_(:)
      real *8, allocatable :: gamma_z_sort_MS_(:)
      complex *16, allocatable :: C_S_sort_MS_(:)
      complex *16, allocatable :: C_Z_sort_MS_(:)
      integer *4, allocatable :: I_permute_MS_(:)
      integer *4, allocatable :: I_inverse_MS_(:)
      real *8 timing_tic,timing_toc
      external fgauss
c
      integer *4 n_ctf,nctf
      integer *4, allocatable :: I_ctf_(:)
      complex *16, allocatable :: CTF_k_p_(:)
      complex *16, allocatable :: CTI_k_p_(:)
      real *8 ctf_p_0,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4
c
      integer *4, allocatable :: I_model_(:)
      integer *4, allocatable :: I_model_sample_(:)

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_ver6]'
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
      allocate(xyz_A(3,n_source_max))
      allocate(xyz_B(3,n_source_max))
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
      allocate(ff_A_ori(ngrid,ngrid,ngrid))
      allocate(ff_B_ori(ngrid,ngrid,ngrid))
      allocate(ff_tru(ngrid,ngrid,ngrid))
      allocate(ff_est(ngrid,ngrid,ngrid))

c-----------------------------------------------------------------------
c     1) Create function in physical space
c-----------------------------------------------------------------------
      pi=4.0d0*atan(1.0d0)
      call prini(6,13)
      
      n_source_use = min(n_source_max,32) ! number of gaussian sources
      n_cut = floor(0.825d0*n_source_use)
      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
      do nsource=0,n_source_use-1
         source_r = 0.5d0*(nsource)/(n_source_use-1)
         source_w = 2*pi*2*nsource/(n_source_use)
         source_z = 1.0d0*(nsource)/(n_source_use-1) - 0.50d0
         xyz_A(1,1+nsource) = source_r*cos(source_w)
         xyz_A(2,1+nsource) = source_r*sin(source_w)
         xyz_A(3,1+nsource) = source_z
         if (nsource.lt.n_cut) then
            xyz_B(1,1+nsource) = xyz_A(1,1+nsource)
            xyz_B(2,1+nsource) = xyz_A(2,1+nsource)
            xyz_B(3,1+nsource) = xyz_A(3,1+nsource)
         else
            source_r = 0.5d0*(n_cut + 3.5d0*(nsource-n_cut))
     $           /(n_source_use-1)
            source_w = 2*pi*2*(n_cut - 0.15d0*(nsource-n_cut)*(nsource
     $           -n_source_use+1))/(n_source_use)
            source_z = 1.0d0*(n_cut - 0.50d0*(nsource-n_cut)*(nsource
     $           -n_source_use+1))/(n_source_use-1) - 0.50d0
            xyz_B(1,1+nsource) = source_r*cos(source_w)
            xyz_B(2,1+nsource) = source_r*sin(source_w)
            xyz_B(3,1+nsource) = source_z
         end if !if (nsource.lt.n_cut) then
      enddo

      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './dir_mda6/source_xyz_A_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,xyz_A,fname)
      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './dir_mda6/source_xyz_B_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,xyz_B,fname)

c     create physical space grid and sample function_A on it.
      ityper = 0
      itypep = 1
      ifpr = 1
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,xyz_A,n_source_use,sigma,
     1     fgauss,ff_A_ori,fnorm)
      if (verbose.gt.1) write(6,*) '  fnorm is ',fnorm
      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
      write(fname,'(A)') './dir_mda6/x_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x,fname)
      write(fname,'(A)') './dir_mda6/y_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y,fname)
      write(fname,'(A)') './dir_mda6/z_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './dir_mda6/S_A_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff_A_ori,fname)
c     create physical space grid and sample function_B on it.
      ityper = 0
      itypep = 1
      ifpr = 1
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,xyz_B,n_source_use,sigma,
     1     fgauss,ff_B_ori,fnorm)
      if (verbose.gt.1) write(6,*) '  fnorm is ',fnorm
      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
c$$$      write(fname,'(A)') './dir_mda6/x_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,x,fname)
c$$$      write(fname,'(A)') './dir_mda6/y_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,y,fname)
c$$$      write(fname,'(A)') './dir_mda6/z_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './dir_mda6/S_B_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff_B_ori,fname)

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
      if (verbose.gt.1) call prinf(' isph_start *',isph_start,ngridr)
      if (verbose.gt.1) call prinf(' nterms_sph *',nterms_sph,ngridr)
      if (verbose.gt.1) call prinf(' nsphstore is *',nsphstore,1)
c
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(ffhat_A_ori(ntot))
      allocate(ffhat_B_ori(ntot))
      allocate(ffhat_tru(ntot))
      allocate(ffhat_est(ntot))
c
      allocate(modsph_A_ori(0:nsphstore-1))
      allocate(modsph_B_ori(0:nsphstore-1))
      allocate(modsph_tru(0:nsphstore-1))
      allocate(modsph_est(0:nsphstore-1))
c
c     create kspace grid data from physical space model A
c
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff_A_ori,h,ngrid,eps,rmax,ngridr,
     $     ityper,nlats,itypep,ntot,numonsphere, kgridx,kgridy,kgridz
     $     ,wts,ffhat_A_ori,ffhatnorm,ier)
      if (verbose.gt.1) then
         write(6,*) 'ffhatnorm: ',ffhatnorm
         write(6,*) 'numonsphere: ',(numonsphere(ii),ii=1,ngridr)
      end if
c     create kspace grid data from physical space model B
c
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff_B_ori,h,ngrid,eps,rmax,ngridr,
     $     ityper,nlats,itypep,ntot,numonsphere, kgridx,kgridy,kgridz
     $     ,wts,ffhat_B_ori,ffhatnorm,ier)
      if (verbose.gt.1) then
         write(6,*) 'ffhatnorm: ',ffhatnorm
         write(6,*) 'numonsphere: ',(numonsphere(ii),ii=1,ngridr)
      end if
c
c     convert to model_A: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(ffhat_A_ori,rmax,nterms_sph, ityper
     $     ,ngridr,nlats,itypep,numonsphere,modsph_A_ori,nsphstore)
c     convert to model_B: spherical harmonic expansions on successive 
c     spheres
c
      call kspacegrid_to_model(ffhat_B_ori,rmax,nterms_sph, ityper
     $     ,ngridr,nlats,itypep,numonsphere,modsph_B_ori,nsphstore)
c$$$      save both modsph_A_ori and modsph_B_ori to disk.
      d_(0) = ngridr
      d_(1) = 1
      write(fname,'(A)') './dir_mda6/isph_start_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(fname)
      call MDA_write_i4(2,d_,isph_start,fname)
      d_(0) = ngridr
      d_(1) = 1
      write(fname,'(A)') './dir_mda6/nterms_sph_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(fname)
      call MDA_write_i4(2,d_,nterms_sph,fname)
      d_(0) = nsphstore
      d_(1) = 1
      write(fname,'(A)') './dir_mda6/modsph_A_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(2,d_,modsph_A_ori,fname)
      d_(0) = nsphstore
      d_(1) = 1
      write(fname,'(A)') './dir_mda6/modsph_B_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(2,d_,modsph_B_ori,fname)

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
         displacement_max = 2.0d0/rmax
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
            alpha2d_tru_(1+nalpha_gamma_z,nnn) = pi/32.0d0*(mod(nnn,32)
     $           -16)
            if (displacement_flag) then
c$$$               angle_tmp = 2*pi*rand()
c$$$               delta_tmp = rand()*displacement_max
c$$$               angle_tmp = 2*pi*(1.0d0*nnn)/(1.0d0*numonsphere(ngridr))
c$$$               delta_tmp = displacement_max/2*(2*(1.0d0*nnn)/(1.0d0
c$$$     $              *numonsphere(ngridr))-1.0d0)
c$$$               alpha2d_tru_(1+nalpha_delta_x,nnn) = cos(angle_tmp)
c$$$     $              *delta_tmp
c$$$               alpha2d_tru_(1+nalpha_delta_y,nnn) = sin(angle_tmp)
c$$$     $              *delta_tmp
               alpha2d_tru_(1+nalpha_delta_x,nnn) = 0.25d0/rmax*(mod(nnn
     $              ,5)-1)
               alpha2d_tru_(1+nalpha_delta_y,nnn) = 0.50d0/rmax*(mod(nnn
     $              ,3)-1)
            else
               alpha2d_tru_(1+nalpha_delta_x,nnn) = 0.0d0/rmax
               alpha2d_tru_(1+nalpha_delta_y,nnn) = 0.0d0/rmax
            end if !displacement_flag
            if (normalization_flag) then
c$$$               draw v from (1/tau)exp(-t/tau) via: v = -tau*log(1-rand())
c$$$               alpha2d_tru_(1+nalpha_l2_norm,nnn) = - 1.0d0 * log(1.0d0 - rand())
c$$$               draw v from [0.5,1.5]
               alpha2d_tru_(1+nalpha_l2_norm,nnn) = 0.5d0 + rand()
            else
               alpha2d_tru_(1+nalpha_l2_norm,nnn) = 1.0d0
            end if !normalization_flag
c$$$            to start with somewhat arbitrary alpha2d_est_
            alpha2d_tru_(1+nalpha_ctf_ind,nnn) = 0.0d0 + mod(nnn,3)
c$$$            estimated parameters used to initialize least-squares-solve
            alpha2d_est_(1+nalpha_polar_a,nnn) = pi*rand()
            alpha2d_est_(1+nalpha_azimu_b,nnn) = 2*rand()*pi
c$$$            alpha2d_est_(1+nalpha_gamma_z,nnn) = 0.0d0
c$$$            alpha2d_est_(1+nalpha_gamma_z,nnn) = (1.0d0 - 2*rand())*pi
            alpha2d_est_(1+nalpha_gamma_z,nnn) = pi/2.0d0*(2*mod(nnn,2)
     $           -1)
c$$$            alpha2d_est_(1+nalpha_delta_x,nnn) = 0.0d0/rmax
c$$$            alpha2d_est_(1+nalpha_delta_y,nnn) = 0.0d0/rmax
c$$$            alpha2d_est_(1+nalpha_delta_x,nnn) = 0.125d0/rmax*(2*rand()
c$$$     $           -1)
c$$$            alpha2d_est_(1+nalpha_delta_y,nnn) = 0.250d0/rmax*(2*rand()
c$$$     $           -1)

            if (displacement_flag) then
               alpha2d_est_(1+nalpha_delta_x,nnn) = 0.125d0/rmax
     $              *(mod(nnn+2,5)-1)
               alpha2d_est_(1+nalpha_delta_y,nnn) = 0.25d0/rmax*(mod(nnn
     $              +1,3)-1)
            else
               alpha2d_est_(1+nalpha_delta_x,nnn) = 0.0d0
               alpha2d_est_(1+nalpha_delta_y,nnn) = 0.0d0
            end if !displacement_flag
            alpha2d_est_(1+nalpha_l2_norm,nnn) = 1.0d0
            alpha2d_est_(1+nalpha_ctf_ind,nnn) = 0.0d0 + mod(nnn,3)
c$$$            to start with alpha2d_est_ .eq. alpha2d_tru_
c$$$            call cp1_r8(n_alpha,alpha2d_tru_(1,1+nm),alpha2d_est_(1,1
c$$$     $           +nm))
         enddo
      enddo
      if (verbose.gt.1) write(6,*) '  nnn is ',nnn
      if (nnn.ne.numonsphere(ngridr)) then
         write(6,'(A,I0,A,I0,A,I0,A)') 'Warning, nnn ',nnn,'.ne.'
     $        ,numonsphere,'(',ngridr,')'
      end if
c
      allocate(cslices_A(nimage_size,nimages))
      allocate(cslices_B(nimage_size,nimages))
      allocate(cslices(nimage_size,nimages))
      allocate(templates(ntemplate_size,nimages))
      allocate(S_alpha_polar_a_all_(0:ntemplates-1))
      allocate(S_alpha_azimu_b_all_(0:ntemplates-1))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     $     icstart)
      if (verbose.gt.1) write(6,*) '  ntemplatesize is ',ntemplatesize

c-----------------------------------------------------------------------
c     4) compute cslices_A directly from physical space function (ff_A_ori).
c        with orientation vectors defined by alpha2d_tru_.
c-----------------------------------------------------------------------
c
      boxsize = 2.0d0
      call mk_simulated_slices_alpha2d(ff_A_ori,ngrid,boxsize,eps,ngridr
     $     ,ityper,ngridc,rmax,n_alpha,alpha2d_tru_,nimages,nimage_size
     $     ,cslices_A)
c-----------------------------------------------------------------------
c     4) compute cslices_B directly from physical space function (ff_B_ori).
c        with orientation vectors defined by alpha2d_tru_.
c-----------------------------------------------------------------------
c
      boxsize = 2.0d0
      call mk_simulated_slices_alpha2d(ff_B_ori,ngrid,boxsize,eps,ngridr
     $     ,ityper,ngridc,rmax,n_alpha,alpha2d_tru_,nimages,nimage_size
     $     ,cslices_B)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Determine which images correspond to which model.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_model_(0:nimages-1))
      do nm=0,nimages-1
         if (mod(nm,7).le.3) I_model_(nm)=0
         if (mod(nm,7).gt.3) I_model_(nm)=1
      enddo !do nm=0,nimages-1
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use cslices_A and cslices_B to fill cslices for later use
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      do nm=0,nimages-1
         if (I_model_(nm).eq.0) then
            call cp1_c16(nimage_size,cslices_A(1,1+nm),cslices(1,1+nm))
         end if !if (I_model_(nm).eq.0) then
         if (I_model_(nm).eq.1) then
            call cp1_c16(nimage_size,cslices_B(1,1+nm),cslices(1,1+nm))
         end if !if (I_model_(nm).eq.0) then
      enddo !do nm=0,nimages-1
      
c-----------------------------------------------------------------------
c     5) generate ctf functions to associate with slices.
c-----------------------------------------------------------------------
c

c$$$  Number of different ctfs used is n_ctf.
c$$$  We expect this to be proportional to the number of images,
c$$$  with batches of roughly 30-300 images each having the same ctf.
      n_ctf = 3
      allocate(CTF_k_p_(0:n_ctf*nimage_size-1))
      allocate(CTI_k_p_(0:n_ctf*nimage_size-1))
c$$$  Here we use an unrealistic ctf-function with a high degree
c$$$  of anisotropy to test our code.
      do nctf=0,n_ctf-1
         if (ctf_flag) then
c$$$  ctf_p_0 = 128.0d0
c$$$  ctf_p_1 = 0.125d0 + 0.0625*nctf
c$$$  ctf_p_2 = 1.5d0
c$$$  ctf_p_3 = 128.0d0
c$$$  ctf_p_4 = 1.0d0
c$$$  call get_ctf_k_p_(ngridr,CTF_k_p_(nctf*nimage_size),ctf_p_0
c$$$  $           ,ctf_p_1,ctf_p_2,ctf_p_3,ctf_p_4)
c$$$  ctf_p_0 = rmax/4.0d0
            ctf_p_0 = 10.0d0 + (5.0d0*nctf)/n_ctf
            ctf_p_1 = (pi*nctf)/n_ctf
            ctf_p_2 = 2
            call get_ctf_star_k_p_(ngridr,xnodesr,ngridc,nimage_size
     $           ,CTF_k_p_(nctf*nimage_size),ctf_p_0,ctf_p_1,ctf_p_2)
         else
            call get_ctf_ones_k_p_(ngridr,xnodesr,ngridc,nimage_size
     $           ,CTF_k_p_(nctf*nimage_size))
         end if                 ! if (ctf_flag) then
         call get_ctf_ones_k_p_(ngridr,xnodesr,ngridc,nimage_size
     $        ,CTI_k_p_(nctf*nimage_size))
      enddo
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out the ctf-functions
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      write(fname,'(A)') './dir_jpg6/Fig_ctf_c'
      call Fig_gen_ctf_ver0(ngridr,nlats,xnodesr,n_ctf,nimage_size
     $     ,CTF_k_p_,fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Randomly determine which images correspond to which ctf-function.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      allocate(I_ctf_(0:nimages-1))
      do nm=0,nimages-1
         I_ctf_(nm) = mod(nm,n_ctf)
      enddo
      do nm=0,nimages-1
         alpha2d_tru_(1+nalpha_ctf_ind,1+nm) = 1.0d0*I_ctf_(nm)
         alpha2d_est_(1+nalpha_ctf_ind,1+nm) = 1.0d0*I_ctf_(nm)
      enddo

c-----------------------------------------------------------------------
c     6) Run the main loop (i.e., try and recover molecule from images).
c-----------------------------------------------------------------------
c     
c$$$  n_M_sample_max = maximum number of images to use at each step.
      n_M_sample_max = 512*1.0
      allocate(I_M_sample_(0:n_M_sample_max-1))
      allocate(I_model_sample_(0:n_M_sample_max-1))
      allocate(M_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(N_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(P_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(Q_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(R_sample_(0:nimage_size*n_M_sample_max-1))
      allocate(L_sample_(0:n_loading*n_M_sample_max-1))
      allocate(C_M_(0:n_M_sample_max-1))
      allocate(alpha_tru_(0:n_alpha*n_M_sample_max-1))
      allocate(alpha_est_(0:n_alpha*n_M_sample_max-1))
c$$$  n_S_sample_max = maximum number of templates to use at each step.
      n_S_sample_max = 256*2.0
      allocate(I_S_sample_(0:n_S_sample_max-1))
      allocate(S_sample_(0:ntemplatesize*n_S_sample_max-1))
      allocate(S_alpha_polar_a_sample_(0:n_S_sample_max-1))
      allocate(S_alpha_azimu_b_sample_(0:n_S_sample_max-1))
      allocate(delta_x_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_S_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_max_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_x_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_S_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_sort_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_permute_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_inverse_SM_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_x_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(delta_y_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(gamma_z_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_S_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(C_Z_sort_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_permute_MS_(0:n_S_sample_max*n_M_sample_max-1))
      allocate(I_inverse_MS_(0:n_S_sample_max*n_M_sample_max-1))
      if (displacement_flag) then
         n_delta_x = 17
         n_delta_y = 17
         N_pixels_in = 1.5d0
      else
         n_delta_x = 1
         n_delta_y = 1
         N_pixels_in = 0.5d0
      end if !displacement_flag
      nwun = 1
      nlow = 3
      nbot = nlow
      do ncur=nbot,ngridr
         if (verbose.gt.0) write(6,*) 'ncur: ',ncur
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  select images to use 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         n_M_sample = min(nimages,n_M_sample_max)
         n_M_A_sample = 0
         n_M_B_sample = 0
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
            I_model_sample_(nM_sample) = I_model_(nm)
            call cp1_c16(nimage_size,cslices(1,1+nm)
     $           ,M_sample_(nM_sample*nimage_size))
            if (I_model_(nm).eq.0) then
               n_M_A_sample = n_M_A_sample + 1
            end if !if (I_model_(nm).eq.0) then
            if (I_model_(nm).eq.1) then
               n_M_B_sample = n_M_B_sample + 1
            end if !if (I_model_(nm).eq.1) then
         enddo                  ! do nM_sample=0,n_M_sample-1
         if (verbose.gt.0) then
            write(6,'(A,I0,A,I0,A,I0)') 'n_M_sample: ' , n_M_sample ,
     $           '; n_M_A_sample: ' , n_M_A_sample , '; n_M_B_sample: '
     $           , n_M_B_sample
         end if !if (verbose.gt.0) then

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha2d to alphas
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'copying alpha2d-->alpha'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha2d_tru_(1,1+nm),alpha_tru_(0
     $           +n_alpha*nM_sample))
            call cp1_r8(n_alpha,alpha2d_est_(1,1+nm),alpha_est_(0
     $           +n_alpha*nM_sample))
         enddo
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write alphas to disk for postmortem analysis.
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         d_(0) = n_alpha
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda6/alpha_tru_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_tru_,fname)
         d_(0) = n_alpha
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda6/alpha_est_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
         call MDA_write_r8(2,d_,alpha_est_,fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  apply ctf-convolution to selected images 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            nctf = I_ctf_(nm)
            call xx1_c16(nimage_size,cslices(1,1+nm),CTF_k_p_(nctf
     $           *nimage_size),N_sample_(nM_sample*nimage_size))
         enddo                  ! do nM_sample=0,n_M_sample-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  multiply selected images by true normalization factor 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            l2_norm = alpha_tru_(nalpha_l2_norm+n_alpha*nM_sample)
            call af1_c16(nimage_size,1.0d0*cmplx(l2_norm,0.0d0),1.0d0
     $           *cmplx(0.0d0,0.0d0),N_sample_(nM_sample*nimage_size)
     $           ,N_sample_(nM_sample*nimage_size))
         enddo                  ! do nM_sample=0,n_M_sample-1

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image+ctf pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         write(fname,'(A,I0)') './dir_jpg6/Fig_MNC_c_',ncur
         call Fig_gen_ctf_ver1(ncur,nlats,xnodesr,n_M_sample,nimage_size
     $        ,M_sample_,N_sample_,nimage_size,CTF_k_p_,n_alpha
     $        ,alpha_est_,min(16,n_M_sample),fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  rebuild model based on current alpha.
c$$$  We rebuild two models.
c$$$  The first, modsph_tru is the model obtained by using the 
c$$$  true alpha-parameters for each image.
c$$$  The second, modsph_est is the model obtained by using the
c$$$  best-guess alpha-parameters for each image (which depends
c$$$  on previous iterations).
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  First find modsph_tru
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
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
         call rebuild_model_alpha2d(N_sample_,nimage_size,n_M_sample
     $        ,n_CTF,nimage_size,CTF_k_p_,n_alpha,alpha_tru_
     $        ,icstart ,nlats,itypep,xnodesr,ngridc,nlow,ncur
     $        ,nsphstore ,isph_start,nterms_sph
     $        ,oversamp,kord,eps,modsph_tru)

c$$$  call rebuild_model_alphaX(M_sample_,nimage_size,n_M_sample
c$$$  $        ,n_alpha,alpha_tru_ ,icstart ,nlats,itypep,xnodesr,ngridc
c$$$  $        ,nlow,ncur ,nsphstore ,isph_start,nterms_sph ,oversamp
c$$$  $        ,kord,eps,modsph_tru)

         timing_toc = second()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'rebuild_model_alpha2d (tru):'
     $           ,' total_time ',timing_toc-timing_tic
         end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'modsph_tru to kspacegrid'
         call model_to_kspacegrid(ffhat_tru,rmax,nterms_sph,ityper
     $        ,ngridr,nlats,itypep,numonsphere,modsph_tru,nsphstore)
         do ii = 1,ntot
            ffhat_tru(ii) = ffhat_tru(ii)*wts(ii)
         enddo
         iflag = -1
         npts = ngrid*ngrid*ngrid
         if (verbose.gt.0) write(6,*) 'finufft3d3_f'
         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_tru,iflag
     $        ,eps,npts,x,y,z,ff_tru,ier)
         d_(0) = ngrid
         d_(1) = ngrid
         d_(2) = ngrid
         write(fname,'(A,I0,A)') './dir_mda6/S_tru_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_c16(3,d_,ff_tru,fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use modsph_tru and alhpa_tru_ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         n_modelsph = nsphstore
         if (ncur.lt.ngridr) n_modelsph = isph_start(ncur+1)
         call test_resid_alpha2d(N_sample_,nimage_size,n_M_sample ,n_CTF
     $        ,nimage_size,CTF_k_p_,n_alpha,alpha_tru_ ,icstart ,nlats
     $        ,itypep,xnodesr,ngridc,nwun,ncur ,n_modelsph ,isph_start
     $        ,nterms_sph ,oversamp,kord,eps,modsph_tru,P_sample_
     $        ,Q_sample_,R_sample_,n_loading,L_sample_)

         n_M_A_p_sample=0
         n_M_A_n_sample=0
         n_M_B_p_sample=0
         n_M_B_n_sample=0
         do nM_sample=0,n_M_sample-1
            if (I_model_sample_(nM_sample).eq.0 .and.
     $           dreal(L_sample_(nM_sample)).ge.0.0d0) n_M_A_p_sample
     $           = n_M_A_p_sample +1
            if (I_model_sample_(nM_sample).eq.1 .and.
     $           dreal(L_sample_(nM_sample)).ge.0.0d0) n_M_B_p_sample
     $           = n_M_B_p_sample +1
            if (I_model_sample_(nM_sample).eq.0 .and.
     $           dreal(L_sample_(nM_sample)).lt.0.0d0) n_M_A_n_sample
     $           = n_M_A_n_sample +1
            if (I_model_sample_(nM_sample).eq.1 .and.
     $           dreal(L_sample_(nM_sample)).lt.0.0d0) n_M_B_n_sample
     $           = n_M_B_n_sample +1
         enddo !do nM_sample=0,n_M_sample-1
         write(6,'(A,I0,1X,I0)') ' alpha_tru: n_M_A_x_sample: ',
     $        n_M_A_p_sample ,n_M_A_n_sample
         write(6,'(A,I0,1X,I0)') ' alpha_tru: n_M_B_x_sample: ',
     $        n_M_B_p_sample ,n_M_B_n_sample
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for modsph_tru
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         write(fname,'(A,I0)') './dir_jpg6/Fig_NPQR_tru_',ncur
         call Fig_R_ver0(ncur,nlats,xnodesr,n_M_sample ,nimage_size
     $        ,N_sample_,P_sample_,Q_sample_,R_sample_,n_alpha
     $        ,alpha_tru_,min(16 ,n_M_sample),fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Second find modsph_est
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'rebuilding modsph_est'
         timing_tic = second()
         call rebuild_model_alpha2d(N_sample_,nimage_size,n_M_sample
     $        ,n_CTF,nimage_size,CTF_k_p_,n_alpha,alpha_est_
     $        ,icstart ,nlats,itypep,xnodesr,ngridc,nlow,ncur
     $        ,nsphstore ,isph_start,nterms_sph
     $        ,oversamp,kord,eps,modsph_est)

c$$$  call rebuild_model_alphaX(M_sample_,nimage_size,n_M_sample
c$$$  $        ,n_alpha,alpha_est_ ,icstart ,nlats,itypep,xnodesr,ngridc
c$$$  $        ,nlow,ncur ,nsphstore ,isph_start,nterms_sph ,oversamp
c$$$  $        ,kord,eps,modsph_est)

         timing_toc = second()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'rebuild_model_alpha2d (est):'
     $           ,' total_time ',timing_toc-timing_tic
         end if
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write model to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'modsph_est to kspacegrid'
         call model_to_kspacegrid(ffhat_est,rmax,nterms_sph,ityper
     $        ,ngridr,nlats,itypep,numonsphere,modsph_est,nsphstore)
         do ii = 1,ntot
            ffhat_est(ii) = ffhat_est(ii)*wts(ii)
         enddo
         iflag = -1
         npts = ngrid*ngrid*ngrid
         if (verbose.gt.0) write(6,*) 'finufft3d3_f'
         call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_est,iflag
     $        ,eps,npts,x,y,z,ff_est,ier)
         d_(0) = ngrid
         d_(1) = ngrid
         d_(2) = ngrid
         write(fname,'(A,I0,A)') './dir_mda6/S_est_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_c16(3,d_,ff_est,fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  Use modsph_est and alhpa_est_ to calculate residuals 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */         
         n_modelsph = nsphstore
         if (ncur.lt.ngridr) n_modelsph = isph_start(ncur+1)
         call test_resid_alpha2d(N_sample_,nimage_size,n_M_sample ,n_CTF
     $        ,nimage_size,CTF_k_p_,n_alpha,alpha_est_ ,icstart ,nlats
     $        ,itypep,xnodesr,ngridc,nwun,ncur ,n_modelsph ,isph_start
     $        ,nterms_sph ,oversamp,kord,eps,modsph_est,P_sample_
     $        ,Q_sample_,R_sample_,n_loading,L_sample_)

         n_M_A_p_sample=0
         n_M_A_n_sample=0
         n_M_B_p_sample=0
         n_M_B_n_sample=0
         do nM_sample=0,n_M_sample-1
            if (I_model_sample_(nM_sample).eq.0 .and.
     $           dreal(L_sample_(nM_sample)).ge.0.0d0) n_M_A_p_sample
     $           = n_M_A_p_sample +1
            if (I_model_sample_(nM_sample).eq.1 .and.
     $           dreal(L_sample_(nM_sample)).ge.0.0d0) n_M_B_p_sample
     $           = n_M_B_p_sample +1
            if (I_model_sample_(nM_sample).eq.0 .and.
     $           dreal(L_sample_(nM_sample)).lt.0.0d0) n_M_A_n_sample
     $           = n_M_A_n_sample +1
            if (I_model_sample_(nM_sample).eq.1 .and.
     $           dreal(L_sample_(nM_sample)).lt.0.0d0) n_M_B_n_sample
     $           = n_M_B_n_sample +1
         enddo !do nM_sample=0,n_M_sample-1
         write(6,'(A,I0,1X,I0)') ' alpha_est: n_M_A_x_sample: ',
     $        n_M_A_p_sample ,n_M_A_n_sample
         write(6,'(A,I0,1X,I0)') ' alpha_est: n_M_B_x_sample: ',
     $        n_M_B_p_sample ,n_M_B_n_sample
         
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-slice pairs for modsph_est
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         write(fname,'(A,I0)') './dir_jpg6/Fig_NPQR_est_',ncur
         call Fig_R_ver0(ncur,nlats,xnodesr,n_M_sample ,nimage_size
     $        ,N_sample_,P_sample_,Q_sample_,R_sample_,n_alpha
     $        ,alpha_est_,min(16 ,n_M_sample),fname)

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  generate templates based on current model (using modsph_est)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         numonsphere_cur = numonsphere(ncur)
         if (verbose.gt.0) write(6,'(A,I0,A)')
     $        'generating numonsphere_cur = ',numonsphere_cur
     $        ,' templates'
         timing_tic = second()
         call template_gen(modsph_est,nsphstore,isph_start, nterms_sph
     $        ,ngridr,ngridc,ntemplatesize,icstart, nlats,itypep,ncur
     $        ,numonsphere_cur,templates)
         timing_toc = second()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.3)') 'template_gen:' ,' total_time '
     $           ,timing_toc-timing_tic
         end if
         call getspheregrid(nlats(ncur),itypep,xnodesth_tmp, sthetas_tmp
     $        ,wtsth_tmp,ngridps_tmp,phsteps_tmp,nspherebig_tmp)
         nnn = 0
         do kk = 1,nlats(ncur)
            ctheta = xnodesth_tmp(kk)
            stheta = sthetas_tmp(kk)
            phstep = phsteps_tmp(kk)
            do ll = 1,ngridps_tmp(kk)
               phi = (ll-1)*phstep
               S_alpha_polar_a_all_(nnn) = dacos(ctheta)
               S_alpha_azimu_b_all_(nnn) = phi
               nnn = nnn+1
            enddo
         enddo
         if (verbose.gt.1) write(6,*) '  nnn is ',nnn
         if (nnn.ne.numonsphere_cur) then
            write(6,'(A,I0,A,I0)') 'Warning, nnn ',nnn,'.ne.'
     $           ,numonsphere(ncur)
         end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  select templates to use
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         n_S_sample = min(numonsphere_cur,n_S_sample_max)
         if (verbose.gt.0) write(6,'(A,I0,A)') 'selecting ',n_S_sample
     $        ,' templates'
         do nS_sample=0,n_S_sample-1
            if (n_S_sample.lt.numonsphere_cur) then
               I_S_sample_(nS_sample) = min(numonsphere_cur-1
     $              ,floor(numonsphere_cur*1.0d0*nS_sample/(n_S_sample
     $              -1)))
            else
               I_S_sample_(nS_sample) = nS_sample
            end if
            ns = I_S_sample_(nS_sample)
            call cp1_c16(ntemplatesize,templates(1,1+ns)
     $           ,S_sample_(nS_sample*ntemplatesize))
            S_alpha_polar_a_sample_(nS_sample) =
     $           S_alpha_polar_a_all_(ns)
            S_alpha_azimu_b_sample_(nS_sample) =
     $           S_alpha_azimu_b_all_(ns)
         enddo

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  calculate innerproducts
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,'(A,I0,A,I0,A,I0,A)') 'get (n_M = '
     $        ,n_M_sample,')-x-(n_S = ',n_S_sample, ') = (', n_M_sample
     $        *n_S_sample,') innerproducts'
         timing_tic = omp_get_wtime()
         if (tesselation_flag) then
            distance_req = 2.0d0 / max(1,ncur)
         else
            distance_req = 2.0d0
         end if !tesselation_flag
         n_gamma_z = max(32,2*ncur)
         call test_innerproduct_batch_wrapper_2d(0,ncur,nlats,xnodesr
     $        ,max_x_c,n_S_sample,ntemplatesize,S_sample_ ,distance_req
     $        ,S_alpha_polar_a_sample_,S_alpha_azimu_b_sample_
     $        ,n_M_sample ,nimage_size,N_sample_,n_CTF,nimage_size
     $        ,CTF_k_p_,n_alpha ,alpha_est_,0.1d0,N_pixels_in
     $        ,displacement_max,n_delta_x ,n_delta_y,n_gamma_z,
     $        svd_calculation_type ,delta_x_max_,delta_y_max_
     $        ,gamma_z_max_,C_M_,C_S_max_ ,C_Z_max_,S_used_all)
         timing_toc = omp_get_wtime()
         if (verbose.gt.0) then
            write(6,'(A,A,F8.4,A,F8.3,A,I0,A,I0,A,I0,A,F8.4)')
     $           'test_innerproduct_batch_wrapper_2d:', ' distance_req '
     $           , distance_req , ' total_time ' ,timing_toc-timing_tic
     $           , '; calculated ' , S_used_all , ' / ' , n_M_sample
     $           *n_S_sample , ' = ' , floor((S_used_all*100.0d0)
     $           /(1.0d0*n_M_sample*n_S_sample)) , 
     $           '% , ; time for each: ' ,
     $           (timing_toc -timing_tic) / max(1,S_used_all)
         end if

c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  write innerproducts to disk for postmortem analysis
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         d_(0) = n_S_sample
         d_(1) = n_M_sample
         write(fname,'(A,I0,A)') './dir_mda6/C_Z_max_',ncur,'_.mda'
         if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
         call MDA_write_c16(2,d_,C_Z_max_,fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  sort results
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         call test_innerproduct_batch_wrapper_3a(0,n_S_sample,n_M_sample
     $        ,delta_x_max_,delta_y_max_,gamma_z_max_,C_S_max_,C_Z_max_
     $        ,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $        ,C_S_sort_SM_ ,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $        ,delta_x_sort_MS_ ,delta_y_sort_MS_,gamma_z_sort_MS_
     $        ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_ ,I_inverse_MS_)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  print out a subset of image-template pairs
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         write(fname,'(A,I0)') './dir_jpg6/Fig_MTS_c_',ncur
         call Fig_gen_ver1c(ncur,nlats,xnodesr,n_S_sample,ntemplatesize
     $        ,S_sample_,n_M_sample,nimage_size,M_sample_,N_sample_
     $        ,nimage_size,CTF_k_p_,n_alpha,alpha_est_ ,delta_x_sort_SM_
     $        ,delta_y_sort_SM_,gamma_z_sort_SM_ ,I_permute_SM_,min(16
     $        ,n_M_sample),fname)
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  update alphas
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         call test_innerproduct_batch_wrapper_4c(0,n_S_sample,n_M_sample
     $        ,C_M_,delta_x_sort_SM_,delta_y_sort_SM_,gamma_z_sort_SM_
     $        ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_,I_inverse_SM_
     $        ,delta_x_sort_MS_,delta_y_sort_MS_,gamma_z_sort_MS_
     $        ,C_S_sort_MS_,C_Z_sort_MS_,I_permute_MS_,I_inverse_MS_
     $        ,S_alpha_polar_a_sample_ ,S_alpha_azimu_b_sample_,n_alpha
     $        ,alpha_est_) 
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  copy alpha to alpha2d
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         if (verbose.gt.0) write(6,*) 'copying alpha-->alpha2d'
         do nM_sample=0,n_M_sample-1
            nm = I_M_sample_(nM_sample)
            call cp1_r8(n_alpha,alpha_est_(0+n_alpha*nM_sample)
     $           ,alpha2d_est_(1,1+nm))
         enddo
         if (normalization_flag) then
c$$$  do nothing
         else
            do nM_sample=0,n_M_sample-1
               nm = I_M_sample_(nM_sample)
               alpha2d_est_(1+nalpha_l2_norm,1+nm) = 1.0d0
            enddo
         end if                 !normalization_flag
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
c$$$  end loop
c$$$  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
         
      enddo                     ! do ncur = nbot,ngridr

 10   continue
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_ver6]'
      end if
      stop
      end
