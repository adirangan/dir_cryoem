c$$$      tests innerproduct routines (in 3d) for spherical-harmonic-models ;
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
      integer, allocatable :: icstart(:)
      integer, allocatable :: ngridps(:)
      integer, allocatable :: ngridps_tmp(:)
      real *8, allocatable :: xnodesr(:)
      real *8, allocatable :: wtsr(:)
      real *8, allocatable :: xyz_A(:,:)
      real *8, allocatable :: xyz_B(:,:)
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      complex *16, allocatable :: ff_A_ori(:,:,:)
      complex *16, allocatable :: ff_B_ori(:,:,:)
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
      integer d_(0:5)
      character(len=1024) fname,format_string
      real *8 timing_tic,timing_toc
      external fgauss
      integer *4 n_beta,nbeta
      real *8, allocatable :: beta_(:)
      real *8 beta
      integer *4 n_k,nk
      real *8, allocatable :: k_(:)
      integer *4, allocatable :: n_l_(:)
      integer *4 n_l_max,n_m_max
      integer *4 n_A,nm,nl
      real *8, allocatable :: W_(:)
      integer *4, allocatable :: n_W_(:)
      complex *16, allocatable :: X_(:)
      integer *4 max_i4_f

      if (verbose.gt.1) then
         write(6,'(A)') '[entering test_innerproduct_3d_0_dr]'
      end if

c-----------------------------------------------------------------------
c     0) allocate some arrays (dynamic allocation required for openmp)
c-----------------------------------------------------------------------
      allocate(nlats(ngridr))
      allocate(numonsphere(ngridr))
      allocate(nterms_sph(ngridr))
      allocate(isph_start(ngridr))
      allocate(ngridc(ntmax))
      allocate(icstart(ngridr))
      allocate(ngridps(ntmax))
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(xyz_A(3,n_source_max))
      allocate(xyz_B(3,n_source_max))
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
      allocate(ff_A_ori(ngrid,ngrid,ngrid))
      allocate(ff_B_ori(ngrid,ngrid,ngrid))

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
      write(fname,'(A)') './dir_mdaT/source_xyz_A_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,xyz_A,fname)
      d_(0) = 3
      d_(1) = n_source_use
      write(fname,'(A)') './dir_mdaT/source_xyz_B_.mda'
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
      write(fname,'(A)') './dir_mdaT/x_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x,fname)
      write(fname,'(A)') './dir_mdaT/y_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y,fname)
      write(fname,'(A)') './dir_mdaT/z_pre_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './dir_mdaT/S_A_ori_.mda'
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
c$$$      write(fname,'(A)') './dir_mdaT/x_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,x,fname)
c$$$      write(fname,'(A)') './dir_mdaT/y_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,y,fname)
c$$$      write(fname,'(A)') './dir_mdaT/z_pre_.mda'
c$$$      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
c$$$      call MDA_write_r8(3,d_,z,fname)
      write(fname,'(A)') './dir_mdaT/S_B_ori_.mda'
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
      write(fname,'(A)') './dir_mdaT/xnodesr_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,xnodesr,fname)
      d_(0) = ngridr
      d_(1) = 1
      write(fname,'(A)') './dir_mdaT/isph_start_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(fname)
      call MDA_write_i4(2,d_,isph_start,fname)
      d_(0) = ngridr
      d_(1) = 1
      write(fname,'(A)') './dir_mdaT/nterms_sph_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing i4 : ',trim(fname)
      call MDA_write_i4(2,d_,nterms_sph,fname)
      d_(0) = nsphstore
      d_(1) = 1
      write(fname,'(A)') './dir_mdaT/modsph_A_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(2,d_,modsph_A_ori,fname)
      d_(0) = nsphstore
      d_(1) = 1
      write(fname,'(A)') './dir_mdaT/modsph_B_ori_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(2,d_,modsph_B_ori,fname)

      n_beta = 33
      allocate(beta_(0:n_beta-1))
      call linspace(-pi,+pi,n_beta,beta_)

      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_beta,'(F8.4,1X))'
         write(6,format_string) 'beta_: ', (beta_(nbeta),nbeta=0,n_beta
     $        -1)
      end if !if (verbose.gt.0) then
      n_k = ngridr
      allocate(k_(0:n_k-1))
      call cp1_r8(n_k,xnodesr,k_)
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_k,'(F8.4,1X))'
         write(6,format_string) 'k_: ', (k_(nk),nk=0,n_k-1)
      end if !if (verbose.gt.0) then
      allocate(n_l_(0:n_k-1))
      call cp1_i4(n_k,nterms_sph,n_l_)
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_k,'(I0,1X))'
         write(6,format_string) 'n_l_: ', (n_l_(nk),nk=0,n_k-1)
      end if !if (verbose.gt.0) then
      n_l_max = max_i4_f(n_k,n_l_)
      n_m_max = 1 + 2*n_l_max
      n_A = 0
      do nl=0,n_l_max
         nm = 1 + 2*nl
         n_A = n_A + nm*nm
      enddo !do nl=0,n_l_max
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_A: ', n_A
      end if !if (verbose.gt.0) then
      allocate(X_(0:n_m_max*n_m_max*n_beta-1))

      call test_innerproduct_3d_0(verbose,n_beta,beta_,n_k,k_,n_l_
     $     ,modsph_A_ori,modsph_B_ori,X_)

      d_(0) = n_m_max
      d_(1) = n_m_max
      d_(2) = n_beta
      write(fname,'(A)') './dir_mdaT/X_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,X_,fname)

 10   continue
      if (verbose.gt.1) then
         write(6,'(A)') '[finished test_innerproduct_3d_0_dr]'
      end if
      stop
      end

