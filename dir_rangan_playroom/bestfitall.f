      subroutine getbestfit(ximagesk,imagesize,nimages,
     1           icstart,nlats,itypep,ncur,ngridc,alphas,xnodesr,
     2           wtsr,templates,itemplatesize,ntemplates,thetas,phis)
C***********************************************************************
C
C     For each image, find best fit using templates from all
c     points on sphere at current k-space resolution (the first two 
c     Euler angles), and a discrete set of rotations for the third 
c     Euler angles.
c
C---------------------------------------------------------------------
C     INPUT:
C
C     ximagesk            images unrolled as point values on successive rings
c                         in k-space from origin outward
C     imagesize           dimension of images (total number of points 
C                         on quasiunform polar grid)
C     nimages             number of images
C     icstart()           indexing array for points on successive circles
c     nlats()             number of quarature nodes in theta on sphere
c                         defined by index ncur.
c     itypep              quadrature scheme in phi direction
c                            0 = same on all latitudes (oversamples poles)
c                            1 = adaptive (undersamples toward poles)
c                         if itypep=0
c                            nphi = nint(nlats(i)*phi_over)
c                         if itypep=1
c                            nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                         phi_over is set in getgridph (typically = 2).
c     ncur                index of highest frequency ring/sphere under 
c                         consideration
c     ngridc()            number of output points/Fourier modes on successive 
c                         circles in templates/images
c     xnodesr             nodes in radial direction
c     wtsr                quadrature weights in radial direction
C     templates           Fourier series on rings in k-space for each image
C     templates           templates unrolled as point values on successive rings
c                         in k-space from origin outward
C     itemplatesize       dimension of templates (series expansions of data
C                         on quasiunform polar grid)
C     ntemplates          number of templates
C     thetas              theta values for template normals
C     phis                phi values for template normals
C
C     OUTPUT:
C
C     alphas              Euler angles for each image
C***********************************************************************

      implicit real *8 (a-h,o-z)
      integer imagesize,nimages
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 alphas(3,nimages)
      real *8 thetas(ntemplates)
      real *8 phis(ntemplates)
      real *8 xnodesr(ncur)
      real *8 wtsr(ncur)
      real *8 xnorm(1)
      real *8, allocatable :: tnorms(:)
      complex *16 ximagesk(imagesize,nimages)
      complex *16 templates(itemplatesize,ntemplates)
      complex *16, allocatable :: templates_fourier(:,:)
      complex *16, allocatable :: ximage_fourier(:)
      complex *16, allocatable :: zprods(:)
c
      allocate(templates_fourier(itemplatesize,ntemplates))
      allocate(ximage_fourier(imagesize))
      allocate(tnorms(ntemplates))
      nzprods = ngridc(ncur)
      allocate(zprods(nzprods))
c
      t0 = second()
      call template_trans_all(ngridc,itemplatesize,icstart,
     1    ncur,ntemplates,templates,templates_fourier)
      t1 = second()
      call templatenorms(templates,itemplatesize,ntemplates,
     1           icstart,nlats,itypep,ncur,ngridc,xnodesr,wtsr,tnorms)
      t2 = second()
      write(6,*)' time for template trans all =',t1-t0
      write(6,*)' time for templatenorms =',t2-t1
c
      do i = 1,nimages
c
         nim1=1
         call templatenorms(ximagesk(1,i),imagesize,nim1,
     1           icstart,nlats,itypep,ncur,ngridc,xnodesr,wtsr,xnorm)
         call template_trans(ngridc,imagesize,icstart,
     1    ncur,ximagesk(1,i),ximage_fourier)
c
         rprodmax = -777.0d0
         do j = 1,ntemplates
ccc            write(6,*)' in bestfit, j is',j
            call innerprod(ximage_fourier,imagesize,icstart,nlats,
     1           itypep,ncur,ngridc,xnodesr,wtsr,templates_fourier(1,j),
     2           itemplatesize,igamma,rgamma,zprods)
ccc            write(6,*)' igamma is',igamma
ccc            write(6,*)' j, rprodmax is',j, abs(zprods(igamma))
ccc     1                      /tnorms(j)/xnorm(1) 
            xtest = dreal(zprods(igamma))/tnorms(j)/xnorm(1)
ccc            xtest = abs(zprods(igamma))/tnorms(j)/xnorm(1)
            if (xtest.gt.rprodmax) then
ccc               write(6,*)' in bestfit, better j is',j
               alphas(1,i) = thetas(j)
               alphas(2,i) = phis(j)
               alphas(3,i) = rgamma
ccc               rprodmax = abs(zprods(igamma))/tnorms(j)/xnorm(1) 
               rprodmax = xtest
            endif
         enddo
      enddo
      return
      end
c
C***********************************************************************
      subroutine innerprod(ximagek_fourier,imagesize,icstart,
     1           nlats,itypep,ncur,ngridc,xnodesr,wtsr,
     2           template_fourier,itemplatesize,igamma,rgamma,zprods)
C***********************************************************************
c
c     Document.....
c     number of modes on ring is Nq = int[(ngridc-1)/2]
c     so -> add up zero mode.
c     go modes 1-Nq and -1,...,-Nq = (ngridc,ngridc-1,...ngridc-nq+1)
c
C---------------------------------------------------------------------
C     INPUT:
C
C     ximagek_fourier     Fourier series on rings in k-space for each image
C     imagesize           dimension of images (series expansions of data
C                         on quasiunform polar grid)
C     icstart()           indexing array for points on successive circles
c     nlats()             number of quarature nodes in theta on sphere
c                         defined by index ncur.
c     itypep              quadrature scheme in phi direction
c                            0 = same on all latitudes (oversamples poles)
c                            1 = adaptive (undersamples toward poles)
c                         if itypep=0
c                            nphi = nint(nlats(i)*phi_over)
c                         if itypep=1
c                            nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                         phi_over is set in getgridph (typically = 2).
c     ncur                index of highest frequency ring/sphere under 
c                         consideration
c     ngridc()            number of output points/Fourier modes on successive 
c                         circles in templates/images
c     xnodesr             nodes in radial direction
c     wtsr                quadrature weights in radial direction
C     template_fourier    Fourier series on rings in k-space for each image
C     itemplatesize       dimension of templates (series expansions of data
C                         on quasiunform polar grid)
C
C     OUTPUT:
C
C     igamma              index of best fit gamma 
C     rgamma              value of best gamma 
C     zprods              inner products with all discrete gammas 
C***********************************************************************


      implicit real *8 (a-h,o-z)
      integer imagesize,nimages
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 pi
      real *8 xnodesr(ncur)
      real *8 wtsr(ncur)
      complex *16 ximagek_fourier(imagesize)
      complex *16 template_fourier(itemplatesize)
      complex *16 zz,zprod
      complex *16 zprods(*)
      complex *16, allocatable :: cn(:)
      complex *16, allocatable :: wsave(:)
      complex *16, allocatable :: wtemp(:)
c
      pi = 4.0d0*datan(1.0d0)
c
      nfftsize = ngridc(ncur)
      allocate(cn(nfftsize))
      allocate(wsave(4*nfftsize+16))
      allocate(wtemp(nfftsize))
c
      do j = 1,nfftsize
         cn(j) = 0.0d0
      enddo
c
      do i = 1,ncur
         cn(1) = cn(1) + ximagek_fourier(icstart(i))*
     1          dconjg(template_fourier(icstart(i)))*
     2          wtsr(i)/xnodesr(i)
         nq = (ngridc(i)-1)/2
         do j = 2,nq+1
            cn(j) = cn(j) + ximagek_fourier(icstart(i)+j-1)*
     1          dconjg(template_fourier(icstart(i)+j-1))*
     2          wtsr(i)/xnodesr(i)
         enddo
         do j = 0,nq-1
            jj = ngridc(i)-j
            cn(nfftsize-j) = cn(nfftsize-j) + 
     1          ximagek_fourier(icstart(i)+jj-1)*
     1          dconjg(template_fourier(icstart(i)+jj-1))*
     2          wtsr(i)/xnodesr(i)
         enddo
      enddo
ccc      write(6,*) 'cn created is ',(cn(j),j=1,nfftsize)
c
c     now evaluate Fourier series with cn coeffs
c
      call dcffti(nfftsize,wsave)
      zz= 0.0d0
      do i = 1,nfftsize
         zz = zz + (2*pi)*cn(i)
      enddo
ccc      write(6,*) 'zz is ',zz
      call dcfftf(nfftsize,cn,wsave)
      do i = 1,nfftsize
         zprods(i) = (2*pi)*cn(i)
      enddo
ccc      write(6,*) 'after transform zprods ',(zprods(j),j=1,nfftsize)
c
c     rgamma is max magnitude inner product
c
      rgamma = 0.0d0
      rmax = 0.0d0
      do i = 1,nfftsize
ccc         if (abs(zprods(i)).gt.rmax) then
         if (real(zprods(i)).gt.rmax) then
ccc            rmax = abs(zprods(i))
            rmax = real(zprods(i))
            igamma = i
            rgamma = 2*pi*(i-1)/nfftsize
         endif
      enddo
      return
      end
C
C***********************************************************************
      subroutine innerprod_dir(ximagek,imagesize,icstart,
     1           nlats,itypep,ncur,ngridc,xnodesr,wtsr,
     2           template,itemplatesize,zprod)
C***********************************************************************
c
c     Document.....
c     number of modes on ring is Nq = int[(ngridc-1)/2]
c     so -> add up zero mode.
c     go modes 1-Nq and -1,...,-Nq = (ngridc,ngridc-1,...ngridc-nq+1)
c
c     Direct calc for fixed rgamma
c
C---------------------------------------------------------------------
C     INPUT:
C
C     ximagek             Fourier series on rings in k-space for each image
C     imagesize           dimension of images (series expansions of data
C                         on quasiunform polar grid)
C     icstart()           indexing array for points on successive circles
c     nlats()             number of quarature nodes in theta on sphere
c                         defined by index ncur.
c     itypep              quadrature scheme in phi direction
c                            0 = same on all latitudes (oversamples poles)
c                            1 = adaptive (undersamples toward poles)
c                         if itypep=0
c                            nphi = nint(nlats(i)*phi_over)
c                         if itypep=1
c                            nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                         phi_over is set in getgridph (typically = 2).
c     ncur                index of highest frequency ring/sphere under 
c                         consideration
c     ngridc()            number of output points/Fourier modes on successive 
c                         circles in templates/images
c     xnodesr             nodes in radial direction
c     wtsr                quadrature weights in radial direction
C     template            Fourier series on rings in k-space for each image
C     itemplatesize       dimension of templates (series expansions of data
C                         on quasiunform polar grid)
C     igamma              0,...,n-1 discrete rotation by gamma 
C
C     OUTPUT:
C
C     zprod               inner product
C***********************************************************************


      implicit real *8 (a-h,o-z)
      integer imagesize,nimages
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 pi
      real *8 xnodesr(ncur)
      real *8 wtsr(ncur)
      complex *16 ximagek(imagesize)
      complex *16 template(itemplatesize)
      complex *16 zz,zprod
      complex *16 phival1,phival2
c
      pi = 4.0d0*datan(1.0d0)
c
      zprod = 0.0d0
      do nn = 1,ncur
         hh = 2*pi/ngridc(nn)
         do jj = 1,ngridc(nn)
            phival1 = ximagek(icstart(nn)+jj-1) 
            phival2 = template(icstart(nn)+jj-1) 
            zprod = zprod + phival1*dconjg(phival2)*
     1                 hh*wtsr(nn)/xnodesr(nn)
            enddo
         enddo
      return
      end
C
C
C
C***********************************************************************
      subroutine templatenorms(templates,itemplatesize,ntemplates,
     1           icstart,nlats,itypep,ncur,ngridc,xnodesr,wtsr,tnorms)
C***********************************************************************
c
c    
c     Document.....
c     number of modes on ring is Nq = int[(ngridc-1)/2]
c     so -> add up zero mode.
c     go modes 1-Nq and -1,...,-Nq = (ngridc,ngridc-1,...ngridc-nq+1)
c
c     Direct calc of norm
c
C---------------------------------------------------------------------
C     INPUT:
C
C     templates           Fourier series on rings in k-space for each image
C     itemplatesize       dimension of templates (series expansions of data
C                         on quasiunform polar grid)
C     ntemplates          number of templates 
C     icstart()           indexing array for points on successive circles
c     nlats()             number of quarature nodes in theta on sphere
c                         defined by index ncur.
c     itypep              quadrature scheme in phi direction
c                            0 = same on all latitudes (oversamples poles)
c                            1 = adaptive (undersamples toward poles)
c                         if itypep=0
c                            nphi = nint(nlats(i)*phi_over)
c                         if itypep=1
c                            nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                         phi_over is set in getgridph (typically = 2).
c     ncur                index of highest frequency ring/sphere under 
c                         consideration
c     ngridc()            number of output points/Fourier modes on successive 
c                         circles in templates/images
c     xnodesr             nodes in radial direction
c     wtsr                quadrature weights in radial direction
C
C     OUTPUT:
C
C     tnorms()            norm out to radius set by ncur
C***********************************************************************


      implicit real *8 (a-h,o-z)
      integer itemplatesize,ntemplates
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 pi,rprod
      real *8 xnodesr(ncur)
      real *8 wtsr(ncur)
      real *8 tnorms(ntemplates)
      complex *16 templates(itemplatesize,ntemplates)
      complex *16 zprod
      complex *16 phival
c
      pi = 4.0d0*datan(1.0d0)
c
      do i = 1,ntemplates
         rprod = 0.0d0
         do nn = 1,ncur
            hh = 2*pi/ngridc(nn)
            do jj = 1,ngridc(nn)
               phival = templates(icstart(nn)+jj-1,i) 
               rprod = rprod + phival*dconjg(phival)*
     1                 hh*wtsr(nn)/xnodesr(nn)
            enddo
         enddo
         tnorms(i) = sqrt(rprod)
      enddo
      return
      end
