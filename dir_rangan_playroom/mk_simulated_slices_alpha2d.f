c***********************************************************************
      subroutine mk_simulated_slices_alpha2d(density,ngrid,boxsize,eps
     $     ,ngridr,ityper,ngridc,rkmax,n_alpha,alpha2d,nslices
     $     ,nslice_size ,cslices)
c***********************************************************************
c$$$      From density in physical space, generate slices of Fourier 
c$$$      transform on 2D grids at requested parameters :
c$$$      alpha2d -> (cos(beta==theta),alpha==phi,gamma_z,delta_x,delta_y).
c$$$      
c$$$      Note: The parameters gamma and delta correspond to the following
c$$$      perturbation of the slice in real-space (i.e., x-space): 
c$$$      *FIRST* translate by (+delta_x,+delta_y), and
c$$$      *SECOND* rotate by +gamma_z.
c$$$
c$$$      Note: This means that when comparing this slice (as an image) to
c$$$      various templates, one must:
c$$$      *FIRST* rotate by -gamma_z, and
c$$$      *SECOND*  translate by (-delta_x,-delta_y).

c-----------------------------------------------------------------------
c     INPUT:
c
c     density       ngrid x ngrid x ngrid array in REAL space
c     ngrid         dimension of density array
c     boxsize       in physical space [-boxsize/2,boxsize/2]^3.
c     eps           NUFFT tolerance 
c
c     ngridr        the extent of fourier array in the radial direction 
c     ngridc()      number of output points on successive circles
c                   in slices
c     rkmax         max frequency
c     n_alpha       leading dimension of alpha2d
c     alpha2d        n_alpha x nimage array of parameters for each simulated image
c                   parameter list for each image:
c                   1: polar angle (e.g., beta)
c                   2: azimutha angle (e.g., alpha)
c                   3: planar angle (i.e., gamma_z)
c                   4: delta_x
c                   5: delta_y
c                   plus any other parameters
c                   Note that we assume n_alpha is at least 5.
c     nslices       number of slices to generate
c     nslice_size   size of slice
c
c     OUTPUT:
c
c     cslices:      nslice_size x nslices array of central slices.
c-----------------------------------------------------------------------
      implicit none
      integer verbose
      data verbose / 1 /
      integer ngrid,ngridr,nslices,nslice_size,ityper,ier,iflag,ii,ir
      integer istart,npts,nptsk,n_A
      integer ngridc(ngridr)
      complex *16  density(ngrid,ngrid,ngrid)
      real *8  boxsize,cthetause,phiuse,gammause,rkmax,rk
      real *8 delta_x,delta_y,gamma_z
      real *8  a,h,eps,t0,t1,t2
      integer n_alpha
      include 'nalpha_define.f'
      real *8  alpha2d(n_alpha,nslices)
      real *8, allocatable :: rkx(:,:)
      real *8, allocatable :: rky(:,:)
      real *8, allocatable :: rkz(:,:)
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      real *8, allocatable :: xnodesr(:),wtsr(:)
      complex *16 cslices(nslice_size,nslices)
      if (verbose.gt.1) write(6,*)
     $     '[entering mk_simulated_slices_alpha2d]'
c
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(rkx(nslice_size,nslices))
      allocate(rky(nslice_size,nslices))
      allocate(rkz(nslice_size,nslices))
      allocate(x(ngrid,ngrid,ngrid))
      allocate(y(ngrid,ngrid,ngrid))
      allocate(z(ngrid,ngrid,ngrid))
c
c     create physical space grid and sample function on it.
c
      a = boxsize/2.0d0
      call mkphysgrid(a,ngrid,h,x,y,z)
c
c     create k-space output points
c
      call getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
c
      t0 = second()
      istart = 1
      do ir = 1,ngridr
         rk = xnodesr(ir)
         do ii = 1,nslices
            cthetause = dcos(alpha2d(1+nalpha_polar_a,ii))
            phiuse = alpha2d(1+nalpha_azimu_b,ii)
            gammause = 0.0d0
            call mkonegreatcircle(rk,cthetause,phiuse,gammause,
     1         ngridc(ir),rkx(istart,ii),rky(istart,ii),rkz(istart,ii))
         enddo
         istart =  istart + ngridc(ir)
ccc         call prinf(' istart is *',istart,1)
      enddo
      t1 = second()
c
c     compute Fourier transform
c
      npts = ngrid*ngrid*ngrid
      nptsk = nslice_size*nslices
      iflag = 1
      if (verbose.gt.1) then
         write(6,*) ' (si) npts is ',npts
         write(6,*) ' nptsk is ',nptsk
      end if
      call  finufft3d3_f(npts,x,y,z,density,iflag,eps,nptsk,
     1     rkx,rky,rkz,cslices,ier)

      n_A = 0
      do ir = 1,ngridr
         n_A = n_A + ngridc(ir)
      end do
      if (n_A.ne.nslice_size) then
         write(6,*) 'Warning, nslice_size.neq.sum(ngridc)'
     $        ,' in mk_simulated_slices_alpha2d'
      end if

      do ii = 1,nslices
         do ir = 1,nslice_size
            cslices(ir,ii) = cslices(ir,ii)*h*h*h
         enddo
         delta_x = +alpha2d(1+nalpha_delta_x,ii)
         delta_y = +alpha2d(1+nalpha_delta_y,ii)
         call transf_p_to_p(ngridr,xnodesr,ngridc,nslice_size,cslices(1
     $        ,ii),+delta_x,+delta_y,cslices(1,ii))
         gamma_z = +alpha2d(1+nalpha_gamma_z,ii)
         call rotate_p_to_p(ngridr,ngridc,nslice_size,cslices(1,ii)
     $        ,gamma_z,cslices(1,ii))
      enddo
      t2 = second()
      if (verbose.gt.0) then
         write(6,*) ' time for data extraction ',t1-t0
         write(6,*) ' time for nufft ',t2-t1
      end if
c
      if (verbose.gt.1) write(6,*)
     $     '[finished mk_simulated_slices_alpha2d]'
      return
      end
