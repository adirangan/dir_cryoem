!> Doxygen comment: ;\n
!> get_template_size: ;\n
!>    Generate quasiuniform grid on concentric circles, with ;\n
!>    the number of points on each circle at radius k corresponding  ;\n
!>    to the number for the great circle on the corresponding sphere ;\n
!>    of radius k (encoded n nlats). ;\n
!>--------------------------------------------------------------------- ;\n
!>    INPUT ;\n
!>    nlats()        number of quarature nodes in theta on sphere ;\n
!>                   defined by index ncur. ;\n
!>    ncur           index of highest frequency sphere under consideration ;\n
!>    OUTPUT          ;\n
!>    ngridc()       number of output points on successive circles  ;\n
!>                   in templates ;\n
!>    ntemplatesize  total number of pts in templates  ;\n
!>    icstart()      indexing array for points on successive circles  ;\n
!>--------------------------------------------------------------------- ;\n
!> template_gen: ;\n
!>********************************************************************* ;\n
!>    Generate all band-limited templates from a model defined by its  ;\n
!>    spherical harmonic expansion on spheres from radius  ;\n
!>    xnodesr(1)... xnoodesr(ncur).  ;\n
!>    This correspondence is implicit:  ;\n
!>    xnodesr is not used in this subroutine. ;\n
!>    A discretization of the outermost sphere with "numonsphere" points  ;\n
!>    (i.e. normal orientations) is defined ;\n
!>    by nlats(ncur),itypep, and phiover.  ;\n
!>    The corresponding central slice intersects ;\n
!>    the model (discretized on spheres) on a sequence of great circles. ;\n
!>    For each point on those great circles, we evaluate the spherical ;\n
!>    harmonic expansion for the corresponding sphere.  ;\n
!>    This yields a central slice sampled on a polar grid for each  ;\n
!>    normal orientation. ;\n
!>--------------------------------------------------------------------- ;\n
!>    INPUT  ;\n
!>    model_sph()  the model in Fourier space, discretized as a  ;\n
!>                 spherical harmonic expansion on successive spheres ;\n
!>                 There are ngridr spheres, of which 1,...,ncur ;\n
!>                 will be used. This is stored in packed (1D) format ;\n
!>    lmod         length of modhat_sph ;\n
!>    isph_start()  array of length ngridr indicating where in the modhat_sph ;\n
!>                 vector the coefficients for the corresponding  ;\n
!>                 sphere begin ;\n
!>    nterms_sph() array of length ngridr defining the orders of the  ;\n
!>                 spherical harmonic expansions on successive spheres. ;\n
!>    ngridr       number of points in radial direction ;\n
!>    ngridc()     number of output points on successive circles in templates ;\n
!>    ntemplatesize    total number of points on quasiuniformly sampled  ;\n
!>                     template ;\n
!>    icstart()    indexing array for points on successive circles  ;\n
!>    nlats()      number of quarature nodes in theta on sphere ;\n
!>                 defined by index ncur. ;\n
!>    itypep       quadrature scheme in phi direction ;\n
!>                      0 = same on all latitudes (which oversamples poles) ;\n
!>                      1 = adaptive (undersamples toward poles) ;\n
!>                   if itypep=0 ;\n
!>                      nphi = nint(nlats(i)*phi_over) ;\n
!>                   if itypep=1 ;\n
!>                      nphi = nint(nlats(i)*phi_over*dsin(theta)) ;\n
!>                   phi_over is set in getgridph (typically = 2). ;\n
!>    ncur          index of largest sphere under consideration ;\n
!>    numonsphere   number of points on sphere at ncur... ;\n
!>    OUTPUT ;\n
!>    templates(:,:)  the "numonsphere" templates, in standard ordering  ;\n
!> Doxygen comment: ;\n
C**********************************************************************
      subroutine get_template_size(nlats,ncur,ntemplatesize,ngridc,
     1           icstart)
C**********************************************************************
c     Generate quasiuniform grid on concentric circles, with
c     the number of points on each circle at radius k corresponding 
c     to the number for the great circle on the corresponding sphere
c     of radius k (encoded n nlats).
c
C----------------------------------------------------------------------
c     INPUT
c
c     nlats()        number of quarature nodes in theta on sphere
c                    defined by index ncur.
c     ncur           index of highest frequency sphere under consideration
c
c     OUTPUT         
c
c     ngridc()       number of output points on successive circles 
c                    in templates
c     ntemplatesize  total number of pts in templates 
c     icstart()      indexing array for points on successive circles 
C----------------------------------------------------------------------
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
C**********************************************************************
      subroutine template_gen(model_sph,lmod,isph_start,
     1     nterms_sph,ngridr,ngridc,ntemplatesize,icstart,
     2     nlats,itypep,ncur,numonsphere,templates)
C**********************************************************************
c     Generate all band-limited templates from a model defined by its 
c     spherical harmonic expansion on spheres from radius 
c     xnodesr(1)... xnoodesr(ncur). 
c     This correspondence is implicit: 
c     xnodesr is not used in this subroutine.
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
C----------------------------------------------------------------------
c     INPUT 
c
c     model_sph()  the model in Fourier space, discretized as a 
c                  spherical harmonic expansion on successive spheres
c                  There are ngridr spheres, of which 1,...,ncur
c                  will be used. This is stored in packed (1D) format
c     lmod         length of modhat_sph
c     isph_start()  array of length ngridr indicating where in the modhat_sph
c                  vector the coefficients for the corresponding 
c                  sphere begin
c     nterms_sph() array of length ngridr defining the orders of the 
c                  spherical harmonic expansions on successive spheres.
c     ngridr       number of points in radial direction
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
c                       nphi = nint(nlats(i)*phi_over*dsin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
c     ncur          index of largest sphere under consideration
c     numonsphere   number of points on sphere at ncur...
c
c     OUTPUT
c
c     templates(:,:)  the "numonsphere" templates, in standard ordering 
C----------------------------------------------------------------------
      implicit none
      integer ngridc(ncur),nlats(ncur),ngridr,ncur,itypep,numonsphere
      integer ntemplatesize,lmod,isph_start(ngridr),nterms_sph(ngridr)
      integer icstart(ncur)
      integer ngcmax
      complex *16 model_sph(0:lmod-1)
      complex *16 templates(ntemplatesize,numonsphere)
c
      integer nn,jj,ns
      real *8 t1,t2
      complex *16, allocatable ::  values(:,:)
c      
      ngcmax = ngridc(ncur)
cc$OMP PARALLEL PRIVATE(values)
      allocate(values(ngcmax,numonsphere))
cc$OMP DO PRIVATE(nn,ns,jj)
      do nn = 1,ncur
         call sheval_greatcircs(model_sph(isph_start(nn)),
     1        nterms_sph(nn),nlats(ncur),itypep,numonsphere,
     1        ngridc(nn),values,ngcmax)
         do ns = 1,numonsphere
            do jj = 1,ngridc(nn)
              templates(icstart(nn)+jj-1,ns) = values(jj,ns)
            enddo
         enddo
      enddo
cc$OMP END DO
cc$OMP END PARALLEL

      return
      end

