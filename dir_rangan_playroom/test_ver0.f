      program test_ver0
      implicit none
      integer verbose
      data verbose / 2 /
      integer ngridr
      parameter (ngridr=20)
      integer nlats_max
      parameter (nlats_max=80)
      integer ncur
      integer ngrid
      parameter (ngrid=60)
      integer ntmax
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
      integer n_alpha
      real *8, allocatable :: alpha_(:)
      real *8 xnodesr(ngridr)
      real *8 wtsr(ngridr)
      real *8 xnodesth(ntmax)
      real *8 sthetas(ntmax)
      real *8 wtsth(ntmax)
      real *8 phsteps(ntmax)
      integer nslice_size,nsphstore,ntot
      real *8 rmax

      real *8 source(3,10)
      integer n_source,nsource
      real *8, allocatable :: source_xyz_(:)
      real *8, allocatable :: center_r_(:)
      real *8, allocatable :: center_x_(:)
      real *8, allocatable :: center_y_(:)
      real *8, allocatable :: center_z_(:)
      integer nA,nb,n_A,d_(0:2)
      real *8 r_3d_x_p_half_max
      integer n_3d_x_c
      real *8 h_3d_x_c
      real *8, allocatable :: x_3d_x_c_(:)
      real *8, allocatable :: y_3d_x_c_(:)
      real *8, allocatable :: z_3d_x_c_(:)
      complex *16, allocatable :: S_3d_x_c_(:)
      real *8 S_3d_x_c_norm
      integer n_3d_k_c,n_3d_k_p
      real *8 r_3d_k_p_half_max
      integer r_3d_k_p_quad_flag,b_3d_k_p_quad_flag
      real *8, allocatable :: r_3d_k_p_(:)
      real *8, allocatable :: r_3d_k_p_weight_(:)
      integer, allocatable :: n_b_3d_k_p_(:)
      integer, allocatable :: n_a_3d_k_p_(:)
      integer n_b_3d_k_p_max,n_a_3d_k_p_sum,n_a_2d_k_p_sum
      integer n_3d_k_p_term_sum,n_3d_k_p_point_sum,n_3d_k_p_point_max
      integer, allocatable :: i_3d_k_p_start_(:)
      integer, allocatable :: n_3d_k_p_term_(:)
      real *8 pi
      character(len=1024) fname
      real *8 eps
      real *8, allocatable ::  x_3d_k_p_(:)
      real *8, allocatable ::  y_3d_k_p_(:)
      real *8, allocatable ::  z_3d_k_p_(:)
      real *8, allocatable ::  xyz_3d_k_p_weight_(:)
      complex *16, allocatable ::  S_3d_k_p_(:)
      real *8 S_3d_k_p_norm
      complex *16, allocatable ::  S_3d_k_Y_ver0_(:)
      complex *16, allocatable ::  S_3d_k_Y_ver1_(:)
      complex *16, allocatable ::  S_3d_k_Y_ver2_(:)
      integer, allocatable :: n_3d_k_p_point_(:)
      integer i_error
      real *8, allocatable ::  b_3d_k_p_(:)
      real *8, allocatable ::  b_3d_k_p_sin_(:)
      real *8, allocatable ::  b_3d_k_p_weight_(:)
      integer, allocatable :: n_a_3d_k_p_tmp_(:)
      real *8, allocatable :: a_3d_k_p_weight_(:)
      integer n_a_3d_k_p_tmp_sum
      real *8 b_3d_k_p,b_3d_k_p_sin,a_3d_k_p_weight,a_3d_k_p_tmp
      integer, allocatable :: n_w_2d_k_p_(:)
      integer, allocatable :: i_w_2d_k_p_(:)
      integer n_w_2d_k_p_sum
      complex *16, allocatable :: S_2d_k_p_(:)
      complex *16, allocatable :: M_2d_k_p_(:)
      real *8, allocatable :: alphas_(:)

      real *8 a,h,fnorm
      real *8 x0y0z0(3,10),sigma 
      real *8 x(ngrid,ngrid,ngrid)
      real *8 y(ngrid,ngrid,ngrid)
      real *8 z(ngrid,ngrid,ngrid)
      real *8 x_s(ngrid,ngrid,ngrid)
      real *8 y_s(ngrid,ngrid,ngrid)
      real *8 z_s(ngrid,ngrid,ngrid)
      complex *16 ff(ngrid,ngrid,ngrid)
      complex *16 ff2(ngrid,ngrid,ngrid)
ccc      complex *16 cslices(350,1232)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: kgridx_s(:)
      real *8, allocatable :: kgridy_s(:)
      real *8, allocatable :: kgridz_s(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: ffhat(:)
      complex *16, allocatable :: ffhat_pos(:)
      complex *16, allocatable :: modsph1(:)
      complex *16, allocatable :: modsph2(:)
      complex *16, allocatable :: modsph3(:)
      complex *16, allocatable :: cslices(:,:)
      complex *16, allocatable :: templates(:,:)
      real *8 ffhatnorm
      integer i,j,k,kk,ll,ii,jj,nslices,nnn,nspherebig,ntemplatesize
      real *8 err,denom,boxsize,ctheta,stheta,phstep
      real *8 phi
      external fgauss
      real *8 rmodelnorm2,rmodelnorm3
      real *8 oversamp
      integer kord,nlow,ncheck,nlowstart,npts,stot,iflag,ier
      real *8 errback,lftmp

      pi=4.0d0*atan(1.0d0)
      call prini(6,13)
c$$$      set number, sigma and centers of of Gaussian sources (x0,y0,z0)
      n_source = 32             ! number of gaussian sources
      sigma = 1.0d0/32.0d0      ! isotropic std of each gaussian source
      allocate(center_r_(0:n_source-1))
      allocate(center_x_(0:n_source-1))
      allocate(center_y_(0:n_source-1))
      allocate(center_z_(0:n_source-1))
      do nsource=0,n_source-1
         center_r_(nsource) = 0.5d0*(nsource)/max(1,n_source-1);
         center_x_(nsource) = center_r_(nsource)*cos(2*pi*2*nsource
     $        /n_source)
         center_y_(nsource) = center_r_(nsource)*sin(2*pi*2*nsource
     $        /n_source)
         center_z_(nsource) = 2.0d0*center_r_(nsource) - 0.5 
      enddo
      allocate(source_xyz_(0:3*n_source-1))
      nA=0
      do nsource=0,n_source-1
         source_xyz_(nA) = center_x_(nsource)
         nA = nA+1
         source_xyz_(nA) = center_y_(nsource)
         nA = nA+1
         source_xyz_(nA) = center_z_(nsource)
         nA = nA+1
      enddo

      d_(0) = 3
      d_(1) = n_source
      write(fname,'(A)') './source_xyz_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(2,d_,source_xyz_,fname)

c$$$      create physical space grid and sample function on it.
      n_3d_x_c = 64             ! number of gridpoints across 3d x-space
                                !  (cartesian grid)
      h_3d_x_c = 0d0            ! grid spacing in each direction
      allocate(x_3d_x_c_(0:(n_3d_x_c**3-1))) ! x position of each
                                           ! gridpoint
      allocate(y_3d_x_c_(0:(n_3d_x_c**3-1))) ! y position of each
                                           ! gridpoint
      allocate(z_3d_x_c_(0:(n_3d_x_c**3-1))) ! z position of each
                                           ! gridpoint
      allocate(S_3d_x_c_(0:(n_3d_x_c**3-1))) ! Function evaluated at each
                                           ! gridpoint

      if (verbose.gt.1) write(6,'(A,A)') 'calling : ','createfun_xspace'
      
      r_3d_x_p_half_max = 1.0d0
      call createfun_xspace(r_3d_x_p_half_max,n_3d_x_c,h_3d_x_c
     $     ,x_3d_x_c_,y_3d_x_c_,z_3d_x_c_,source_xyz_,n_source,sigma
     $     ,fgauss,S_3d_x_c_,S_3d_x_c_norm)

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      a = 1.0d0
      call createfun_xspace(a,ngrid,h,x,y,z,source_xyz_,n_source,sigma,
     $     fgauss,ff,fnorm)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      d_(0) = n_3d_x_c
      d_(1) = n_3d_x_c
      d_(2) = n_3d_x_c
      write(fname,'(A)') './x_3d_x_c_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x_3d_x_c_,fname)
      write(fname,'(A)') './y_3d_x_c_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y_3d_x_c_,fname)
      write(fname,'(A)') './z_3d_x_c_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z_3d_x_c_,fname)
      write(fname,'(A)') './S_3d_x_c_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,S_3d_x_c_,fname)

      n_3d_k_c = 20
      n_3d_k_p = n_3d_k_c
      r_3d_k_p_half_max = 0.5d0*n_3d_k_p
      r_3d_k_p_quad_flag = 0
      b_3d_k_p_quad_flag = 1
      allocate(n_b_3d_k_p_(0:n_3d_k_p-1))
      allocate(n_a_3d_k_p_(0:n_3d_k_p-1))
      allocate(r_3d_k_p_(0:n_3d_k_p-1))
      allocate(r_3d_k_p_weight_(0:n_3d_k_p-1))
      allocate(i_3d_k_p_start_(0:n_3d_k_p-1))
      allocate(n_3d_k_p_term_(0:n_3d_k_p-1))
      if (verbose.gt.1) then 
         write(6,*) 'calling : ','getgridr'
         write(6,*) 'r_3d_k_p_half_max: ',r_3d_k_p_half_max
      end if
      call getgridr(r_3d_k_p_half_max,n_3d_k_p,r_3d_k_p_quad_flag
     $     ,r_3d_k_p_,r_3d_k_p_weight_)
      if (verbose.gt.1) write(6,*) 'xnodesr:', (r_3d_k_p_(na),na=0
     $     ,n_3d_k_p-1)
      n_3d_k_p_term_sum = 0
      n_a_3d_k_p_sum = 0
      n_a_2d_k_p_sum = 0
      do na=0,n_3d_k_p-1
         n_b_3d_k_p_(na) = max(6,nint(pi*r_3d_k_p_(na)))
         n_a_3d_k_p_(na) = 2*n_b_3d_k_p_(na)
         i_3d_k_p_start_(na) = n_3d_k_p_term_sum
         n_3d_k_p_term_(na) = nint(r_3d_k_p_(na) + 2)
         n_3d_k_p_term_sum = n_3d_k_p_term_sum + (n_3d_k_p_term_(na)+1)
     $        **2
         n_a_3d_k_p_sum = n_a_3d_k_p_sum + n_a_3d_k_p_(na)
         n_a_2d_k_p_sum = n_a_2d_k_p_sum + n_a_3d_k_p_(na)
      enddo
      n_b_3d_k_p_max = n_b_3d_k_p_(n_3d_k_p-1)
      if (verbose.gt.1) then
         write(6,*) 'n_b_3d_k_p_: ',(n_b_3d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*) 'n_b_3d_k_p_max: ',n_b_3d_k_p_max
         write(6,*) 'n_a_3d_k_p_: ',(n_a_3d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*) 'n_a_3d_k_p_sum: ',n_a_3d_k_p_sum
         write(6,*) 'calling: ','get_ntot'
      end if
      call get_ntot(n_3d_k_p,n_b_3d_k_p_,b_3d_k_p_quad_flag
     $     ,n_3d_k_p_point_sum)
      if (verbose.gt.1) then
         write(6,*) 'n_3d_k_p_point_sum: ',n_3d_k_p_point_sum
         write(6,*) 'n_3d_k_p_term_: ',(n_3d_k_p_term_(na),na=0,n_3d_k_p
     $        -1)
         write(6,*) 'n_3d_k_p_term_sum: ',n_3d_k_p_term_sum
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      rmax = 10.0d0
      call getgridr(rmax,n_3d_k_p,r_3d_k_p_quad_flag,xnodesr,wtsr)
      if (verbose.gt.1) write(6,*) 'xnodesr:', (xnodesr(na),na=1
     $     ,n_3d_k_p)
      nsphstore = 0
      nslice_size = 0
      do na = 1,ngridr
         nlats(na) = nint(pi*xnodesr(na))
         if (nlats(na).lt.6) nlats(na) = 6
         ngridc(na) = nlats(na)*2
         isph_start(na) = nsphstore
         nterms_sph(na) = nint(xnodesr(na) + 2)
         nsphstore = nsphstore + (nterms_sph(na)+1)**2
         nslice_size = nslice_size + ngridc(na)
      enddo
      if (verbose.gt.1) call prinf('  nlats is *',nlats,ngridr)
      if (verbose.gt.1) call prinf('  ngridc is *',ngridc,ngridr)
      if (verbose.gt.1) call prinf('  nslice_size is *',nslice_size,1)
      call get_ntot(ngridr,nlats,b_3d_k_p_quad_flag,ntot)
      if (verbose.gt.1) write(6,*) '  ntot is ',ntot
      if (verbose.gt.1) call prinf(' nterms_sph is *',nterms_sph,ngridr)
      if (verbose.gt.1) call prinf(' nsphstore is *',nsphstore,1)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      call testfourier(ngrid,ngridr,rmax,r_3d_k_p_quad_flag
c$$$     $     ,b_3d_k_p_quad_flag,nlats,fgauss,n_source,source_xyz_,sigma
c$$$     $     ,err,errback)
c$$$      goto 10
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      allocate(x_3d_k_p_(0:n_3d_k_p_point_sum-1))
      allocate(y_3d_k_p_(0:n_3d_k_p_point_sum-1))
      allocate(z_3d_k_p_(0:n_3d_k_p_point_sum-1))
      allocate(xyz_3d_k_p_weight_(0:n_3d_k_p_point_sum-1))
      allocate(S_3d_k_p_(0:n_3d_k_p_point_sum-1))
      eps = 1.0d-6
      allocate(n_3d_k_p_point_(0:n_3d_k_p-1))
      if (verbose.gt.1) then
         write(6,*) 'calling: ','fun_to_kspacegrid'
      end if
      call fun_to_kspacegrid(x_3d_x_c_,y_3d_x_c_,z_3d_x_c_,S_3d_x_c_
     $     ,h_3d_x_c,n_3d_x_c,eps,r_3d_k_p_half_max,n_3d_k_p
     $     ,r_3d_k_p_quad_flag,n_b_3d_k_p_,b_3d_k_p_quad_flag
     $     ,n_3d_k_p_point_sum,n_3d_k_p_point_,x_3d_k_p_,y_3d_k_p_
     $     ,z_3d_k_p_,xyz_3d_k_p_weight_,S_3d_k_p_,S_3d_k_p_norm
     $     ,i_error)
      n_3d_k_p_point_max = n_3d_k_p_point_(n_3d_k_p-1)
      if (verbose.gt.1) then
         write(6,*) 'S_3d_k_p_norm: ',S_3d_k_p_norm
      end if
      if (verbose.gt.1) then
         write(6,*),'S_3d_k_p_: ',(S_3d_k_p_(na),na=0,1)
         write(6,*),'r_3d_k_p_half_max: ',r_3d_k_p_half_max
         write(6,*),'n_3d_k_p_term_: ',(n_3d_k_p_term_(na),na=0,n_3d_k_p
     $        -1)
         write(6,*),'n_3d_k_p: ',n_3d_k_p
         write(6,*),'n_b_3d_k_p_: ',(n_b_3d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*),'n_3d_k_p_point_: ',(n_3d_k_p_point_(na),na=0
     $        ,n_3d_k_p-1)
         write(6,*),'n_3d_k_p_point_max: ',n_3d_k_p_point_max
         write(6,*),'n_3d_k_p_term_sum: ',n_3d_k_p_term_sum
      end if
      allocate(S_3d_k_Y_ver0_(0:n_3d_k_p_term_sum-1))
      allocate(S_3d_k_Y_ver1_(0:n_3d_k_p_term_sum-1))
      allocate(S_3d_k_Y_ver2_(0:n_3d_k_p_term_sum-1))
      if (verbose.gt.1) then
         write(6,*) 'calling: ','kspacegrid_to_model'
      end if
      call kspacegrid_to_model(S_3d_k_p_,r_3d_k_p_half_max
     $     ,n_3d_k_p_term_,r_3d_k_p_quad_flag,n_3d_k_p,n_b_3d_k_p_
     $     ,b_3d_k_p_quad_flag,n_3d_k_p_point_,S_3d_k_Y_ver0_
     $     ,n_3d_k_p_term_sum)
      if (verbose.gt.1) then
         write(6,*),'S_3d_k_Y_ver0_: ',(S_3d_k_Y_ver0_(na),na=0,2)
      end if

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(kgridx_s(ntot))
      allocate(kgridy_s(ntot))
      allocate(kgridz_s(ntot))
      allocate(wts(ntot))
      allocate(ffhat(ntot))
      allocate(ffhat_pos(ntot))
      eps = 1.0d-6
      call fun_to_kspacegrid(x,y,z,ff,h,ngrid,eps,rmax
     $     ,ngridr,r_3d_k_p_quad_flag,nlats,b_3d_k_p_quad_flag,ntot
     $     ,numonsphere,kgridx ,kgridy,kgridz,wts,ffhat,ffhatnorm
     $     ,i_error)
      if (verbose.gt.1) write(6,*) '  ffhatnorm is ',ffhatnorm
      if (verbose.gt.1) then
         lftmp=0;
         do i=1,ntot
            if (kgridx(i).gt.lftmp) lftmp=kgridx(i)
         enddo
         write(6,*),'max(kgridx): ',lftmp
         lftmp=0;
         do i=1,ntot
            if (kgridx(i).lt.lftmp) lftmp=kgridx(i)
         enddo
         write(6,*),'min(kgridx): ',lftmp
c$$$         write(6,*),'kgridx: ',(kgridx(na),na=1,ntot)
c$$$         write(6,*),'kgridy: ',(kgridy(na),na=1,ntot)
c$$$         write(6,*),'kgridz: ',(kgridz(na),na=1,ntot)
         write(6,*),'ffhat: ',(ffhat(na),na=1,2)
         write(6,*),'rmax: ',rmax
         write(6,*),'nterms_sph: ',(nterms_sph(na),na=1,ngridr)
         write(6,*),'ngridr: ',ngridr
         write(6,*),'nlats: ',(nlats(na),na=1,ngridr)
         write(6,*),'numonsphere: ',(numonsphere(na),na=1,ngridr)
         write(6,*),'nsphstore: ',nsphstore
      end if
      allocate(modsph1(0:nsphstore-1))
      allocate(modsph2(0:nsphstore-1))
      allocate(modsph3(0:nsphstore-1))
      call kspacegrid_to_model(ffhat,rmax,nterms_sph,r_3d_k_p_quad_flag
     $     ,ngridr,nlats,b_3d_k_p_quad_flag,numonsphere,modsph1
     $     ,nsphstore)
      if (verbose.gt.1) then
         write(6,*),'modsph1: ',(modsph1(na),na=0,2)
      end if
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      allocate(b_3d_k_p_(0:n_b_3d_k_p_max-1))
      allocate(b_3d_k_p_sin_(0:n_b_3d_k_p_max-1))
      allocate(b_3d_k_p_weight_(0:n_b_3d_k_p_max-1))
      allocate(n_a_3d_k_p_tmp_(0:n_b_3d_k_p_max-1))
      allocate(a_3d_k_p_weight_(0:n_b_3d_k_p_max-1))
      if (verbose.gt.1) write(6,*) 'n_3d_k_p_point_max: '
     $     ,n_3d_k_p_point_max
      call getspheregrid(n_b_3d_k_p_max,b_3d_k_p_quad_flag,b_3d_k_p_,
     $     b_3d_k_p_sin_,b_3d_k_p_weight_,n_a_3d_k_p_tmp_
     $     ,a_3d_k_p_weight_ ,n_a_3d_k_p_tmp_sum)
      if (verbose.gt.1) then
         write(6,*) 'n_b_3d_k_p_max: ',n_b_3d_k_p_max
         write(6,*) 'b_3d_k_p_: ',(b_3d_k_p_(na),na=0,n_b_3d_k_p_max-1)
         write(6,*) 'b_3d_k_p_sin_: ',(b_3d_k_p_sin_(na),na=0
     $        ,n_b_3d_k_p_max-1)
         write(6,*) 'b_3d_k_p_weight_: ',(b_3d_k_p_weight_(na),na=0
     $        ,n_b_3d_k_p_max-1)
         write(6,*) 'n_a_3d_k_p_tmp_: ',(n_a_3d_k_p_tmp_(na),na=0
     $        ,n_b_3d_k_p_max-1)
         write(6,*) 'a_3d_k_p_weight_: ',(a_3d_k_p_weight_(na),na=0
     $        ,n_b_3d_k_p_max-1)
         write(6,*) 'n_a_3d_k_p_tmp_sum: ',n_a_3d_k_p_tmp_sum
      end if
      allocate(alphas_(0:3*n_a_3d_k_p_tmp_sum-1))
      n_A = 0
      do nb = 0,n_b_3d_k_p_max-1
         b_3d_k_p = b_3d_k_p_(nb)
         b_3d_k_p_sin = b_3d_k_p_sin_(nb)
         a_3d_k_p_weight = a_3d_k_p_weight_(nb)
c$$$         write(6,*) '1,2,3,4: ',b_3d_k_p,b_3d_k_p_sin,a_3d_k_p_weight
c$$$     $        ,n_a_3d_k_p_tmp_(nb)
         do na = 0,n_a_3d_k_p_tmp_(nb)-1
            a_3d_k_p_tmp = na*a_3d_k_p_weight
            alphas_(0 + 3*n_A) = dacos(b_3d_k_p)
            alphas_(1 + 3*n_A) = a_3d_k_p_tmp
            alphas_(2 + 3*n_A) = na*2*pi/n_a_3d_k_p_tmp_(nb)
c$$$            alphas_(2 + 3*n_A) = -rand()
c$$$            alphas_(2 + 3*n_A) = 0.0d0
            n_A = n_A+1
         enddo
      enddo
      if (verbose.gt.1) write(6,*) '  n_A is ',n_A
      allocate(n_w_2d_k_p_(0:n_3d_k_p-1))
      allocate(i_w_2d_k_p_(0:n_3d_k_p-1))
      call get_template_size(n_b_3d_k_p_,n_3d_k_p,n_w_2d_k_p_sum
     $     ,n_w_2d_k_p_,i_w_2d_k_p_)
      if (verbose.gt.1) then
         write(6,*) 'n_a_2d_k_p_sum: ',n_a_2d_k_p_sum
         write(6,*) 'n_3d_k_p_point_max: ',n_3d_k_p_point_max
         write(6,*) 'n_b_3d_k_p_: ',(n_b_3d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*) 'n_w_2d_k_p_sum: ',n_w_2d_k_p_sum
         write(6,*) 'n_w_2d_k_p_: ',(n_w_2d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*) 'i_w_2d_k_p_: ',(i_w_2d_k_p_(na),na=0,n_3d_k_p-1)
         write(6,*) 'alphas_: ',(alphas_(na),na=0,2)
         write(6,*) 'alphas_: ',(alphas_(na),na=3,5)
         write(6,*) 'alphas_: ',(alphas_(na),na=6,8)
         write(6,*) 'alphas_: ',(alphas_(na),na=9,11)
      end if
      allocate(S_2d_k_p_(0:n_a_2d_k_p_sum*n_3d_k_p_point_max-1))
      allocate(M_2d_k_p_(0:n_a_2d_k_p_sum*n_3d_k_p_point_max-1))
      if (verbose.gt.1) write(6,*) 'calling: ','mk_simulated slices'
      call mk_simulated_slices(S_3d_x_c_,n_3d_x_c,2.0d0
     $     *r_3d_x_p_half_max,eps,n_3d_k_p,r_3d_k_p_quad_flag
     $     ,n_a_3d_k_p_,r_3d_k_p_half_max,alphas_,n_3d_k_p_point_max
     $     ,n_a_3d_k_p_sum,M_2d_k_p_)

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      boxsize = 2.0d0*r_3d_x_p_half_max
      nslices = numonsphere(ngridr)
      if (verbose.gt.1) write(6,*) '  nslices is ',nslices
      call getspheregrid(nlats(ngridr),b_3d_k_p_quad_flag,xnodesth,
     1           sthetas,wtsth,ngridps,phsteps,nspherebig)
      if (verbose.gt.1) then
         write(6,*) 'nlats(ngridr): ',nlats(ngridr)
         write(6,*) 'xnodesth: ',(xnodesth(na),na=1,nlats(ngridr))
         write(6,*) 'sthetas: ',(sthetas(na),na=1,nlats(ngridr))
         write(6,*) 'wtsth: ',(wtsth(na),na=1,nlats(ngridr))
         write(6,*) 'ngridps: ',(ngridps(na),na=1,nlats(ngridr))
         write(6,*) 'phsteps: ',(phsteps(na),na=1,nlats(ngridr))
         write(6,*) 'nspherebig: ',nspherebig
      end if
      nnn = 0
      do kk = 1,nlats(ngridr)
         ctheta = xnodesth(kk)
         stheta = sthetas(kk)
         phstep = phsteps(kk)
c$$$         write(6,*) '1,2,3: ',ctheta,stheta,phstep,ngridps(kk)
         do ll = 1,ngridps(kk)
            phi = (ll-1)*phstep
            nnn = nnn+1
            alphas(1,nnn) = dacos(ctheta)
            alphas(2,nnn) = phi
            alphas(3,nnn) = (ll-1)*2*pi/ngridps(kk)
c$$$            alphas(3,nnn) = -rand()
c$$$            alphas(3,nnn) = 0.0d0
         enddo
      enddo
      if (verbose.gt.1) write(6,*) '  nnn is ',nnn
      allocate(cslices(nslice_size,nslices))
      allocate(templates(nslice_size,nslices))
      call get_template_size(nlats,ngridr,ntemplatesize,ngridc2,
     1           icstart)
      if (verbose.gt.1) then
         write(6,*) 'nslices: ',nslices
         write(6,*) 'nslice_size: ',nslice_size
         write(6,*) 'nlats: ',(nlats(na),na=1,ngridr)
         write(6,*) 'ntemplatesize is ',ntemplatesize
         write(6,*) 'ngridc2: ',(ngridc2(na),na=1,ngridr)
         write(6,*) 'icstart: ',(icstart(na),na=1,ngridr)
         write(6,*) 'alphas: ',(alphas(na,1),na=1,3)
         write(6,*) 'alphas: ',(alphas(na,2),na=1,3)
         write(6,*) 'alphas: ',(alphas(na,3),na=1,3)
         write(6,*) 'alphas: ',(alphas(na,4),na=1,3)
      end if
      call mk_simulated_slices(ff,ngrid,boxsize,eps,ngridr,
     $     r_3d_k_p_quad_flag,ngridc,rmax,alphas,nslices,nslice_size
     $     ,cslices)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      call template_gen(S_3d_k_Y_ver0_,n_3d_k_p_term_sum,i_3d_k_p_start_
     $     ,n_3d_k_p_term_,n_3d_k_p,n_a_3d_k_p_,n_w_2d_k_p_sum
     $     ,i_w_2d_k_p_,n_b_3d_k_p_,b_3d_k_p_quad_flag,n_3d_k_p
     $     ,n_3d_k_p_point_max,S_2d_k_p_)
      err = 0.0d0
      denom = 0.0d0
      do na=0,n_a_2d_k_p_sum*n_3d_k_p_point_max-1
         denom = denom + abs(M_2d_k_p_(na))**2
         err = err + abs(S_2d_k_p_(na) - M_2d_k_p_(na))**2
      enddo
      write(6,*) 'err in S_2d_k_p_ vs M_2d_k_p_ = ', dsqrt(err/denom)

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call template_gen(modsph1,nsphstore,isph_start,
     1     nterms_sph,ngridr,ngridc,ntemplatesize,icstart,
     2     nlats,b_3d_k_p_quad_flag,ngridr,nslices,templates)
      err = 0.0d0
      denom = 0.0d0
      do jj = 1,nslices
      do ii = 1,nslice_size
         denom = denom + abs(cslices(ii,jj))**2
         err = err + abs(cslices(ii,jj) - templates(ii,jj))**2
      enddo
      enddo
      write(6,*) 'err in template_gen vs mk_simulated_slices = ',
     1            dsqrt(err/denom)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      oversamp = 2.0d0
      eps = 1.0d-6
      kord = 5
      nlow = 3
      ncur = 10
      call rebuild_full_model(M_2d_k_p_,n_a_3d_k_p_sum
     $     ,n_3d_k_p_point_max,alphas_,i_w_2d_k_p_,n_b_3d_k_p_
     $     ,b_3d_k_p_quad_flag,n_a_3d_k_p_,nlow,ncur,n_3d_k_p_term_sum
     $     ,i_3d_k_p_start_,n_3d_k_p_term_,oversamp,kord,eps
     $     ,S_3d_k_Y_ver1_)
      ncheck = 0
      do na = 0,ncur-1
         ncheck = ncheck + (n_3d_k_p_term_(na)+1)**2
      enddo
      nlowstart = 0
      do na = 0,nlow-2
         nlowstart = nlowstart + (n_3d_k_p_term_(na)+1)**2
      enddo
      do na = 0,nlowstart-1
         S_3d_k_Y_ver2_(na) = 0.0d0
      enddo
      do na = nlowstart,ncheck-1
         S_3d_k_Y_ver2_(na) = S_3d_k_Y_ver0_(na) - S_3d_k_Y_ver1_(na)
      enddo
      call kspace_model_norm(S_3d_k_Y_ver2_,ncheck,nterms_sph
     $     ,r_3d_k_p_weight_,ncur,rmodelnorm3)
      call kspace_model_norm(S_3d_k_Y_ver1_,ncheck,nterms_sph
     $     ,r_3d_k_p_weight_,ncur,rmodelnorm2)
      write(6,*) 'rel error in rebuilt model is ',rmodelnorm3
     $     /rmodelnorm2
      

c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      oversamp = 2.0d0
      eps = 1.0d-6
      kord = 5
      nlow = 3
      call rebuild_full_model(cslices,nslice_size,nslices,alphas,
     $     icstart,nlats,b_3d_k_p_quad_flag,ngridc,nlow,ncur,nsphstore,
     $     isph_start,nterms_sph,oversamp,kord,eps,modsph2)
      ncheck = 0
      do ii = 1,ncur
         ncheck = ncheck + (nterms_sph(ii)+1)**2
      enddo
      nlowstart = 0
      do ii = 1,nlow-1
         nlowstart = nlowstart + (nterms_sph(ii)+1)**2
      enddo
      do ii = 0,nlowstart-1
         modsph3(ii) = 0.0d0
      enddo
      do ii = nlowstart,ncheck-1
         modsph3(ii) = modsph1(ii) - modsph2(ii)
      enddo
      call kspace_model_norm(modsph3,ncheck,nterms_sph,wtsr,ncur,
     1          rmodelnorm3)
      call kspace_model_norm(modsph2,ncheck,nterms_sph,wtsr,ncur,
     1          rmodelnorm2)
      write(6,*) 'rel error in rebuilt model is ',rmodelnorm3
     $     /rmodelnorm2
c$$$      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      call model_to_kspacegrid(ffhat_pos,rmax,nterms_sph,
     $     r_3d_k_p_quad_flag,ngridr,nlats,b_3d_k_p_quad_flag
     $     ,numonsphere,modsph2,nsphstore)

      do k = 1,ngrid
      do j = 1,ngrid
      do i = 1,ngrid
          x_s(i,j,k) = i-1-ngrid/2
          y_s(i,j,k) = j-1-ngrid/2
          z_s(i,j,k) = k-1-ngrid/2
      enddo
      enddo
      enddo
      do ii = 1,ntot
         ffhat_pos(ii) = ffhat_pos(ii)*wts(ii)
         kgridx_s(ii) = pi*kgridx(ii)/rmax
         kgridy_s(ii) = pi*kgridy(ii)/rmax
         kgridz_s(ii) = pi*kgridz(ii)/rmax
      enddo
      iflag = -1
      npts = ngrid*ngrid*ngrid
c$$$      call  finufft3d3_f(ntot,kgridx,kgridy,kgridz,ffhat_pos,iflag,eps,
c$$$     1        npts,x,y,z,ff2,ier)
      call  finufft3d1_f(ntot,kgridx_s,kgridy_s,kgridz_s,ffhat_pos,iflag
     $     ,eps,ngrid,ngrid,ngrid,ff2,ier)
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
      if (verbose.gt.1) write(6,*) 'l2 err in ff is: ',errback

      d_(0) = ngrid
      d_(1) = ngrid
      d_(2) = ngrid
      write(fname,'(A)') './x_3d_x_c_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,x_s,fname)
      write(fname,'(A)') './y_3d_x_c_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,y_s,fname)
      write(fname,'(A)') './z_3d_x_c_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing r8 : ',trim(fname)
      call MDA_write_r8(3,d_,z_s,fname)
      write(fname,'(A)') './S_3d_x_c_pos_.mda'
      if (verbose.gt.1) write(6,'(A,A)') 'Writing c16 : ',trim(fname)
      call MDA_write_c16(3,d_,ff2,fname)

      deallocate(center_r_)
      deallocate(center_x_)
      deallocate(center_y_)
      deallocate(center_z_)
      deallocate(source_xyz_)
      deallocate(x_3d_x_c_)
      deallocate(y_3d_x_c_)
      deallocate(z_3d_x_c_)
      deallocate(S_3d_x_c_)
      deallocate(n_b_3d_k_p_)
      deallocate(n_a_3d_k_p_)
      deallocate(r_3d_k_p_)
      deallocate(r_3d_k_p_weight_)
      deallocate(i_3d_k_p_start_)
      deallocate(n_3d_k_p_term_)
      deallocate(x_3d_k_p_)
      deallocate(y_3d_k_p_)
      deallocate(z_3d_k_p_)
      deallocate(xyz_3d_k_p_weight_)
      deallocate(S_3d_k_p_)
      deallocate(n_3d_k_p_point_)
      deallocate(S_3d_k_Y_ver0_)
      deallocate(S_3d_k_Y_ver1_)

 10   continue
      stop
      end

      include 'MDA_write_r8.f'
      include 'MDA_write_c16.f'

