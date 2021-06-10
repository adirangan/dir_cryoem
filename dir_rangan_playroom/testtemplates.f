c
c      Test template generation on quasiuniform polar grid
c
c     1) define function in k-space
c     2) sample function on concentric spheres
c     3) convert to spherical harmonic expansions on concentric spheres
c     4) call template generator
c     5) compare with direct eval of function at corresponding points.
c
      implicit none
c
      integer itypep,ngridr,ityper
      parameter (ngridr=60)
      integer i,ii,is,nlats(ngridr)
      integer j,l2,ncur,ngridc,ngridpmax,ngridtmax,nn,nnn,ns
      integer nsphstart,nsphstore,nstart,ntot,ntemplatesize
      integer icstart(ngridr),ngridcv(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      integer isph_start(ngridr)
      real *8 x0y0z0(3,10),sigma,done,phi_over,pi,rkmax,rl2,rl2err
      real *8 t2,t1
      real *8, allocatable :: xnodesr(:),wtsr(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      real *8, allocatable :: kgridx(:)
      real *8, allocatable :: kgridy(:)
      real *8, allocatable :: kgridz(:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: xg(:,:)
      real *8, allocatable :: yg(:,:)
      real *8, allocatable :: zg(:,:)
c
      complex *16 zsum,zc
      complex *16, allocatable :: modsph(:)
      complex *16, allocatable :: templates(:,:)
      complex *16, allocatable :: templates_pft(:,:)
      complex *16, allocatable :: phival(:)
      complex *16, allocatable :: phidir(:,:)
      complex *16, allocatable :: phipft(:,:)
      complex *16 eye
      external fgauss
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0d0*datan(done)
      call prini(6,13)
c
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
c
c     set sigma and location of Gaussian sources (x0,y0,z0)
c
      sigma = 1.0d0
      x0y0z0(1,1) = 0.1d0
      x0y0z0(2,1) = 0.2d0
      x0y0z0(3,1) = -0.13d0
      x0y0z0(1,2) = -0.1d0
      x0y0z0(2,2) = -0.2d0
      x0y0z0(3,2) = 0.23d0
      ns = 2
c
      ityper = 0
      itypep = 1
      phi_over = 2.0d0
      rkmax = 60.0d0
      ngridc = 400
c
c     get radial grid parameters
c     set nterms on each shell
c     set nlats on each shell
c     compute total storage needed for s.h. expansions nterms on each shell
c
      call getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
      call prin2(' xnodesr is *',xnodesr,ngridr)
      nsphstore = 0
      do i = 1,ngridr
         isph_start(i) = nsphstore
         nlats(i) = nint(pi*xnodesr(i))
         if (nlats(i).lt.12) nlats(i) = 12
         nterms_sph(i) = nint(xnodesr(i) + 2)
         nsphstore = nsphstore + (nterms_sph(i)+1)**2
      enddo
      allocate(modsph(0:nsphstore-1))
      call prinf(' nterms_sph is *',nterms_sph,ngridr)
      call prinf(' nsphstore is *',nsphstore,1)
c
c     get total number of points on spherical grid and 
c     create it
c
      ngridtmax = nlats(1)
      do i = 2,ngridr
         ngridtmax = max(nlats(i),ngridtmax)
      enddo
      ngridpmax = nint(ngridtmax*phi_over)
c
      allocate(xnodesth(ngridtmax))
      allocate(wtsth(ngridtmax))
      allocate(xnodesph(ngridpmax))
      allocate(wtsph(ngridpmax))
c
      call get_ntot(ngridr,nlats,itypep,ntot)
c
      allocate(kgridx(ntot))
      allocate(kgridy(ntot))
      allocate(kgridz(ntot))
      allocate(wts(ntot))
      allocate(phival(ntot))
c
      call mkfouriergrid(rkmax,ngridr,ityper,
     1                  nlats,itypep,ntot,numonsphere,
     2                    kgridx,kgridy,kgridz,wts)
ccc      call prinf(' after mkfouriergrid, ntot is *',ntot,1)
ccc      call prinf(' after mkf.., nlats is *',nlats,ngridr)
ccc      call prinf(' after mkf.., numonsphere is *',numonsphere,ngridr)
c
c     evaluate function on spherical grid
c     project onto spherical harmonics on successive shells
c
      do i = 1,ntot
         call fgauss(kgridx(i),kgridy(i),kgridz(i),x0y0z0,ns,sigma,
     1        phival(i))
         phival(i) = eye*phival(i)
      enddo
ccc      call prin2(' phival is *',phival,2*ntot)
      nstart = 1
      nsphstart = 0
      do i = 1,ngridr
         call projshexp(phival(nstart),numonsphere(i),nlats(i),itypep,
     1       nterms_sph(i),modsph(nsphstart))
         nstart = nstart + numonsphere(i)
         nsphstart = nsphstart + (nterms_sph(i)+1)**2
      enddo
c
c     set ncur and generate templates to that resolution
c
      ncur = 10
      call get_template_size(nlats,ncur,ntemplatesize,ngridcv,icstart)
      call prinf(' ngridcv is *',ngridcv,ncur)
      call prinf(' icstart is *',icstart,ncur)
      call prinf(' ntemplatesize is *',ntemplatesize,1)
      allocate(templates(ntemplatesize,numonsphere(ncur)))
      allocate(templates_pft(ntemplatesize,numonsphere(ncur)))
      t1 = second()
      call template_gen(modsph,nsphstore,isph_start,
ccc     1     nterms_sph,ngridr,xnodesr,wtsr,ngridcv,ntemplatesize,
     1     nterms_sph,ngridr,ngridcv,ntemplatesize,
     2     icstart,nlats,itypep,ncur,numonsphere(ncur),templates)
      t2 = second()
      call prin2(' time for template gen *',t2-t1,1)
ccc      write(77,*)(templates(i,2),i=1,ntemplatesize)
ccc      call prin2(' phival *',templates,2*ngridc*ncur*numonsphere(ncur))
c
c     compute values of testing function at corresponding grid points
c     directly.
c
      allocate(phidir(ntemplatesize,numonsphere(ncur)))
      allocate(phipft(ntemplatesize,numonsphere(ncur)))
      do nn = 1,ncur
         allocate(xg(ngridcv(nn),numonsphere(ncur)))
         allocate(yg(ngridcv(nn),numonsphere(ncur)))
         allocate(zg(ngridcv(nn),numonsphere(ncur)))
         call mkallgreatcircles(xnodesr(nn),nlats(ncur),itypep,
     1           numonsphere(ncur),ngridcv(nn),xg,yg,zg)
         do j = 1,numonsphere(ncur)
         do i = 1,ngridcv(nn)
            call fgauss(xg(i,j),yg(i,j),zg(i,j),x0y0z0,ns,sigma,
     1           phidir(icstart(nn)+i-1,j))
            phidir(icstart(nn)+i-1,j) = eye*phidir(icstart(nn)+i-1,j)
         enddo
         enddo
         deallocate(xg)
         deallocate(yg)
         deallocate(zg)
      enddo
ccc      write(78,*)(phidir(i,2),i=1,ntemplatesize)
ccc      call prin2(' phdir is *',phidir(1,1,1),
ccc     1     2*ncur*ngridc*numonsphere(ncur))
c
c     check L2 error (without proper quadrature weighting)
c
      nnn =  ntemplatesize*numonsphere(ncur)
      call l2error_c(phidir, templates, nnn, rl2err, rl2)
      call prin2(' l2 is *',rl2,1)
      call prin2(' l2err is *',rl2err,1)
      call prin2(' relerr is *',rl2err/rl2,1)
c
c     Fourier transform templates ("hat hat" representation)
c     and evaluate Fourier series at grid points
c
ccc      do j = 1,numonsphere(ncur)
ccc         call template_trans(ngridcv,ntemplatesize,icstart,
ccc     1           ncur,templates(1,j),templates_pft(1,j))
ccc      enddo
      call template_trans_all(ngridcv,ntemplatesize,icstart,
     1           ncur,numonsphere(ncur),templates,templates_pft)
c
c     evaluate Fourier series at target points
c
      do nn = 1,ncur
         call prinf(' nn is *',nn,1)
         do j = 1,numonsphere(ncur)
         do i = 1,ngridcv(nn)
            zsum = 0.0d0
            do ii = 1,ngridcv(nn)
               zc = templates_pft(icstart(nn)+ii-1,j)
               zsum=zsum+zc*cdexp(eye*2*pi*(i-1)*(ii-1)/ngridcv(nn))
            enddo
            phipft(icstart(nn)+i-1,j) = zsum
         enddo
         enddo
      enddo
      call l2error_c(phidir, phipft, nnn, rl2err, rl2)
      call prin2(' l2 is *',rl2,1)
      call prin2(' l2err is *',rl2err,1)
      call prin2(' relerr is *',rl2err/rl2,1)

      stop
      end
c
c
c
c
c
