c
c      Test inner product on quasiuniform polar grid
c
c     1) define two functions f1, f2 in k-space
c     2) compute inner product at various rotations of f1.
c     3) call innerprod routine
c     4) compare 
c
      implicit none
c
      integer itypep,ngridr,ityper
      parameter (ngridr=40)
      integer i,igamma,ii,is,jj,nlats(ngridr)
      integer j,l2,ncur,ngridc,ngridpmax,ngridtmax,nn,nnn,ns
      integer nsphstart,nsphstore,nstart,ntot
      integer jgamma,imagesize,itemplatesize,nspin
      integer icstart(ngridr),ngridcv(ngridr)
      integer numonsphere(ngridr)
      integer nterms_sph(ngridr)
      integer isph_start(ngridr)
      real *8 x1y1(2,10),sigma1,done,phi_over,pi,rkmax,rl2,rl2err
      real *8 x2y2(2,10),sigma2
      real *8 area,t2,t1,hh,rgamma,theta,xx,yy
      real *8, allocatable :: xnodesr(:),wtsr(:)
      real *8, allocatable :: xnodesth(:),wtsth(:)
      real *8, allocatable :: xnodesph(:),wtsph(:)
      real *8, allocatable :: wts(:)
      complex *16, allocatable :: f1(:)
      complex *16, allocatable :: f2(:)
      complex *16, allocatable :: f1_fourier(:)
      complex *16, allocatable :: f2_fourier(:)
c
      complex *16 zsum,zc
      complex *16 phival
      complex *16 phival1,phival2
      complex *16 eye,zprod,zprods(1000)
      complex *16 zprodd(1000)
      external fgauss2d
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
      sigma1 = 1.0d0
      sigma2 = 1.0d0
      x1y1(1,1) = 5.1d0
      x1y1(2,1) = 0.2d0
      x1y1(1,2) = -0.1d0
      x1y1(2,2) = -0.2d0
c
      x2y2(1,1) = 5.1d0
      x2y2(2,1) = 0.2d0
      x2y2(1,2) = 0.1d0
      x2y2(2,2) = -0.2d0
      ns = 1
c
      ityper = 1
      itypep = 1
      phi_over = 2.0d0
      rkmax = 10.0d0
c
c     get radial grid parameters
c     set nterms on each shell
c     set nlats on each shell
c     compute total storage needed for s.h. expansions nterms on each shell
c
      call getgridr(rkmax,ngridr,ityper,xnodesr,wtsr)
      call prin2(' xnodesr is *',xnodesr,ngridr)
c
      do i = 1,ngridr
         nlats(i) = nint(pi*xnodesr(i))
ccc         if (nlats(i).lt.12) nlats(i) = 12
         if (nlats(i).lt.32) nlats(i) = 32
      enddo
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
      ncur = 40
      area = pi*xnodesr(ncur)**2
      area = pi*rkmax**2
      write(6,*) ' area is ',area
      call get_template_size(nlats,ncur,itemplatesize,ngridcv,icstart)
      allocate(f1(itemplatesize))
      allocate(f2(itemplatesize))
      allocate(f1_fourier(itemplatesize))
      allocate(f2_fourier(itemplatesize))
c
      nspin = ngridcv(ncur)
      do jgamma = 1,nspin
ccc      do jgamma = 10,10
ccc      write(21,*)' x = ['
ccc      do jgamma = 0,0
         rgamma =  (jgamma-1)*2*pi/nspin
         write(6,*) ' rgamma is ',rgamma
         zprod = 0.0d0
         do nn = 1,ncur
            hh = 2*pi/ngridcv(nn)
            do jj = 1,ngridcv(nn)
               theta = (jj-1)*hh
               xx = xnodesr(nn)*dcos(theta-rgamma)
               yy = xnodesr(nn)*dsin(theta-rgamma)
               call fgauss2d(xx,yy,x1y1,ns,sigma1,
     1           phival1)
ccc               phival1 = 1.0d0
               f1(icstart(nn)+jj-1) = phival1 + eye*phival1/2
ccc               write(21,*)xx,yy,abs(phival1)
               xx = xnodesr(nn)*dcos(theta)
               yy = xnodesr(nn)*dsin(theta)
               call fgauss2d(xx,yy,x2y2,ns,sigma2,
     1           phival2)
ccc               phival2 = 1.0d0
               f2(icstart(nn)+jj-1) = phival2 + eye*phival2/3
               zprod = zprod + phival1*dconjg(phival2)*
     1                 hh*wtsr(nn)/xnodesr(nn)
            enddo
         enddo
ccc         write(21,*) ']'
         zprodd(jgamma) = zprod
         write(6,*)' zprod is ',zprod
         write(6,*)' abs(zprod) is ',abs(zprod)
         write(22,*)rgamma,abs(zprod)
      enddo
c
      if (2.ne.3) then
         do nn = 1,ncur
            hh = 2*pi/ngridcv(nn)
            do jj = 1,ngridcv(nn)
               theta = (jj-1)*hh
               xx = xnodesr(nn)*dcos(theta)
               yy = xnodesr(nn)*dsin(theta)
               call fgauss2d(xx,yy,x1y1,ns,sigma1,
     1           phival1)
ccc               phival1 = 1.0d0
               f1(icstart(nn)+jj-1) = phival1 + eye*phival/2
ccc               write(21,*)xx,yy,abs(phival1)
               xx = xnodesr(nn)*dcos(theta)
               yy = xnodesr(nn)*dsin(theta)
               call fgauss2d(xx,yy,x2y2,ns,sigma2,
     1           phival2)
ccc               phival2 = 1.0d0
               f2(icstart(nn)+jj-1) = phival2 + eye*phival/3
               zprod = zprod + phival1*dconjg(phival2)*
     1                 hh*wtsr(nn)/xnodesr(nn)
            enddo
         enddo
      endif 

      call template_trans(ngridcv,itemplatesize,icstart,
     1           ncur,f1,f1_fourier)
      call template_trans(ngridcv,itemplatesize,icstart,
     1           ncur,f2,f2_fourier)
c
      imagesize = itemplatesize
      call innerprod(f1_fourier,imagesize,icstart,
     1           nlats,itypep,ncur,ngridcv,xnodesr,wtsr,
     2           f2_fourier,itemplatesize,igamma,rgamma,zprods)
c
      write(6,*) ' igamma,rgamma are ',igamma,rgamma
      do i = 1,nspin
         write(6,*) i,abs(zprodd(i)),abs(zprods(i)),
     1    abs(zprodd(i)-zprods(i))
      enddo
ccc      write(6,*) ' one pt test'
ccc      write(6,*) abs(zprodd(2)),abs(zprods(1)),
ccc     1    abs(zprodd(2)-zprods(1))
ccc      write(6,*) ' one pt test'
ccc      write(6,*) abs(zprodd(nspin)),abs(zprods(nspin-1)),
ccc     1    abs(zprodd(nspin)-zprods(nspin-1))
ccc      write(6,*) ' one pt test'
ccc      write(6,*) abs(zprodd(1)),abs(zprods(nspin)),
ccc     1    abs(zprodd(1)-zprods(nspin))
      stop
      end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fgauss2d(x,y,x0y0,ns,sigma,ff)
c
c     INPUT:
c
c     x,y     target location
c     x0y0    location of Gaussian sources
c     ns      number of Gaussian sources
c     sigma   number variance of all Gaussian sources
c
c     OUTPUT:
c
c     ff        complex function value at corresponding grid pt.
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real *8 (a-h,o-z)
      real *8 x0y0(2,ns)
      complex *16 ff
c
      ff = 0.0d0
      do i = 1,ns
         xx = x - x0y0(1,i)
         yy = y - x0y0(2,i)
         rr = xx*xx+yy*yy
         ff = ff + dexp(-rr/(2*sigma*sigma))
      enddo
      return
      end
c
