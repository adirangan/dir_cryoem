ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_gen_old(model_sph,lmod,sph_start,
     1     nterms_sph,ngridr,xnodesr,wtsr,ngridc,nlats,itypep,
     2     ncur,numonsphere,templates)
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
c     model_sph    the model in Fourier space, discretized as a 
c                  spherical harmonic expansion on successive spheres
c                  There are ngridr spheres, of which 1,...,ncur
c                  will be used. This is stored in packed (1D) format
c     lmod         length of modhat_sph
c     sph_size     array of length ngridr specifying the length of the 
c                  expansion on each sphere.  
c     sph_start    array of length ngridr indicating where in the modhat_sph
c                  vector the coefficients for the corresponding 
c                  sphere begin
c     nterms_sph   array of length ngridr defining the orders of the 
c                  spherical harmonic expansions on successive spheres.
c     ngridr       number of points in radial direction
c     xnodesr      sphere radii 
c     wtsr         quadrature weights in radial direction
c     ngridc       number of output points on each circle in templates
c     nlats        number of quarature nodes in theta on sphere
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
c     templates -    the "numonsphere" templates, in standard ordering 
c     thetas    -    the corresponding spherical angles w.r.t. z axis
c     phis      -    the corresponding spherical angles w.r.t. x axis

      implicit none
      integer ngridc,nlats(ncur),ngridr,ncur,itypep,numonsphere
      integer lmod,sph_start(ngridr),nterms_sph(ngridr)
      real *8 xnodesr(ngridr),wtsr(ngridr)
      complex *16 model_sph(0:lmod-1)
      complex *16 templates(ngridc,ncur,numonsphere)

      integer nn,jj,ns
      real *8 t1,t2
      complex *16, allocatable ::  values(:,:)
      
      allocate(values(ngridc,numonsphere))
      do nn = 1,ncur
ccc         call prin2(' inside template_gen modsph *',
ccc     1        model_sph(sph_start(nn)),2*(nterms_sph(nn)+1)**2)
c
         call sheval_greatcircs(model_sph(sph_start(nn)),
     1        nterms_sph(nn),nlats(ncur),itypep,numonsphere,
     1        ngridc,values,ngridc)
ccc         call prin2(' values = *',values,2*ngridc*numonsphere)
         do ns = 1,numonsphere
            do jj = 1,ngridc
              templates(jj,nn,ns) = values(jj,ns)
            enddo
         enddo
      enddo
c
      return
      end

