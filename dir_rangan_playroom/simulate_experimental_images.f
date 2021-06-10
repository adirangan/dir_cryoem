c***********************************************************************
      subroutine mk_simulated_slices(density,ngrid,boxsize,eps,ngridr,
     1       ityper,ngridc,rkmax,alphas,nslices,nslice_size,cslices)
c***********************************************************************
c     From density in physical space, generate slices of Fourier 
c     transform on 2D grids at requested rotations 
c     [alphas -> (alpha,beta,gamma) for now].
c
c     Simulated images are computed in make_experimental_images.
c     which compute the 2D transforms of these slices after 
c     multiplication with CTF.
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
c     alphas        3 x nimage array of rotations for each simulated image
c     nslices       number of slices to generate
c     nslice_size   size of slice
c
c     OUTPUT:
c
c     cslices:      nslice_size x nslices array of central slices.
c-----------------------------------------------------------------------
      implicit none
      integer ngrid,ngridr,nslices,nslice_size,ityper,ier,iflag,ii,ir
      integer istart,npts,nptsk
      integer ngridc(ngridr)
      complex *16  density(ngrid,ngrid,ngrid)
      real *8  boxsize,cthetause,phiuse,gammause,rkmax,rk
      real *8  a,h,eps,t0,t1,t2
      real *8  alphas(3,nslices)
      real *8, allocatable :: rkx(:,:)
      real *8, allocatable :: rky(:,:)
      real *8, allocatable :: rkz(:,:)
      real *8, allocatable :: x(:,:,:)
      real *8, allocatable :: y(:,:,:)
      real *8, allocatable :: z(:,:,:)
      real *8, allocatable :: xnodesr(:),wtsr(:)
      complex *16 cslices(nslice_size,nslices)
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
            cthetause = dcos(alphas(1,ii))
            phiuse = alphas(2,ii)
            gammause = alphas(3,ii)
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
      write(6,*) ' (si) npts is ',npts
      write(6,*) ' nptsk is ',nptsk
      call  finufft3d3_f(npts,x,y,z,density,iflag,eps,nptsk,
     1     rkx,rky,rkz,cslices,ier)

      do ii = 1,nslices
      do ir = 1,nslice_size
         cslices(ir,ii) = cslices(ir,ii)*h*h*h
      enddo
      enddo
      t2 = second()
      write(6,*) ' time for data extraction ',t1-t0
      write(6,*) ' time for nufft ',t2-t1
c
      return
      end
c
ccc      subroutine mk_simulated_images(cslices,nslices,nslice_size,
ccc     1           ngridr,ityper,ngridc,ctf,boxsize,ngrid,rimages)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Generate experimental projections (images) in 2D REAL space from a
c     ground truth 3D density array (density) in 3D REAL space, for purpose
c     of simulations. 

c     input:
c     cslices      nslice_size x nslices array of central slices.
c     nslices      number of slices
c     nslice_size  number of discretization pts in each slice
c     ngridr       number of radial pts in k space
c     ityper       type of radial grid
c     ngridc()     number of output points on successive circles
c                  in slices
c     ctf          nslice_size x nslices array of ctf values for each
c                                        slice in k-space.
c     boxsize      physical space boxsize [-boxsize/2,boxsize/2]^3.
c     ngrid        number of points in linear dimension in phys space
c
c     output:
c
c     rimages      ngrid x ngrid x nslices 
c                  array of 2D projections in REAL space
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
