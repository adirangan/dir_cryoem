ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_template_size(nlats,ncur,ntemplatesize,ngridc,
     1           icstart)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Generate quasiuniform grid on concentric circles, with
c     the number of points on each circle at radius k corresponding 
c     to the number for the great circle on the corresponding sphere
c     of radius k (encoded n nlats).
c
c     INPUT:
c     nlats()      number of quarature nodes in theta on sphere
c                  defined by index ncur.
c     ncur         index of highest frequency sphere under consideration
c
c     OUT:         
c     ngridc()       number of output points on successive circles 
c                    in templates
c     ntemplatesize  total number of pts in templates 
c     icstart()      indexing array for points on successive circles 
c
      implicit none
      integer ntemplatesize,nlats(ncur),ncur,ngridc(ncur),icstart(ncur)
      integer i
c
      icstart(1) = 1
      do i = 1,ncur-1
         ngridc(i) = nlats(i)*2     
         icstart(i+1) = icstart(i) + ngridc(i)
      enddo
      ngridc(ncur) = nlats(ncur)*2     
      ntemplatesize = icstart(ncur) + ngridc(ncur) - 1
      return
      end
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_gen(model_sph,lmod,sph_start,
     1     nterms_sph,ngridr,xnodesr,wtsr,ngridc,ntemplatesize,
     2     icstart,nlats,itypep,ncur,numonsphere,templates)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Generate all band-limited templates from a model defined by its 
c     spherical harmonic expansion on spheres from radius 
c     xnodesr(1)... xnoodesr(ncur). 
c     A discretization of the outermost sphere with "numonsphere" points 
c     (i.e. normal orientations) is defined
c     by nlats(ncur),itypep, and phiover. 
c     The corresponding central slice intersects
c     the model (discretized on spheres) on a sequence of great circles.
c     For each point on those great circles, we evaluate the spherical
c     harmonic expansion for the corresponding sphere. 
c     This yields a central slice sampled on a polar grid for each 
c     normal orientation.
c
c     INPUT: 
c
c     model_sph()  the model in Fourier space, discretized as a 
c                  spherical harmonic expansion on successive spheres
c                  There are ngridr spheres, of which 1,...,ncur
c                  will be used. This is stored in packed (1D) format
c     lmod         length of modhat_sph
c     sph_start()  array of length ngridr indicating where in the modhat_sph
c                  vector the coefficients for the corresponding 
c                  sphere begin
c     nterms_sph() array of length ngridr defining the orders of the 
c                  spherical harmonic expansions on successive spheres.
c     ngridr       number of points in radial direction
c     xnodesr()    sphere radii 
c     wtsr()       quadrature weights in radial direction
c     ngridc()     number of output points on successive circles in templates
c     ntemplatesize    total number of points on quasiuniformly sampled 
c                      template
c     icstart()    indexing array for points on successive circles 
c     nlats()      number of quarature nodes in theta on sphere
c                  defined by index ncur.
c     itypep       quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
c     ncur          index of largest sphere under consideration
c     numonsphere   number of points on sphere at ncur...
c
c     OUTPUT:
c
c     templates(:,:)  the "numonsphere" templates, in standard ordering 
c
      implicit none
      integer ngridc(ncur),nlats(ncur),ngridr,ncur,itypep,numonsphere
      integer ntemplatesize,lmod,sph_start(ngridr),nterms_sph(ngridr)
      integer icstart(ncur)
      integer ngcmax
      real *8 xnodesr(ngridr),wtsr(ngridr)
      complex *16 model_sph(0:lmod-1)
ccc      complex *16 templates(ngridc,ncur,numonsphere)
      complex *16 templates(ntemplatesize,numonsphere)

      integer nn,jj,ns
      real *8 t1,t2
      complex *16, allocatable ::  values(:,:)
      
      ngcmax = ngridc(ncur)
      allocate(values(ngcmax,numonsphere))
      do nn = 1,ncur
         call sheval_greatcircs(model_sph(sph_start(nn)),
     1        nterms_sph(nn),nlats(ncur),itypep,numonsphere,
     1        ngridc(nn),values,ngcmax)
         do ns = 1,numonsphere
            do jj = 1,ngridc(nn)
ccc              templates(jj,nn,ns) = values(jj,ns)
              templates(icstart(nn)+jj-1,ns) = values(jj,ns)
            enddo
         enddo
      enddo
c
      return
      end

