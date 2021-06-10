C***********************************************************************
      subroutine rebuild_model_alpha2d(cslices,imagesize,nimages,n_CTF
     $     ,ld_CTF,CTF_p_,n_alpha,alpha2d,icstart,nlats,itypep
     $     ,xnodesr,ngridc,nlow,ncur,lmod ,isph_start
     $     ,nterms_sph ,oversamp,kord,eps,modelsph)
C***********************************************************************
C     Builds new model (in spherical harmonic basis) 
C     from slices (Fourier transforms of images) with image parameters 
C     given by alpha2d array.
C---------------------------------------------------------------------
C     INPUT:
C
C     cslices        Fourier transforms of images on polar grids
C     imagesize      dimension of images on quasiunform polar grid
C     nimages        number of images
C     n_CTF          number of ctf-functions
C     ld_CTF         leading dimension of CTF-array (usually imagesize)
C     CTF_p_         stack of ctf-functions (size ld_CTF*n_CTF)
c     n_alpha        size of alpha2d
c     alpha2d(3,*)    image parameters (see nalpha_define.f) 
c
c     icstart()      indexing array for points on successive circles
c     nlats()        number of quarature nodes in theta on sphere
c                    defined by index ncur.
c     itypep         quadrature scheme in phi direction
c                       0 = same on all latitudes (which oversamples poles)
c                       1 = adaptive (undersamples toward poles)
c                    if itypep=0
c                       nphi = nint(nlats(i)*phi_over)
c                    if itypep=1
c                       nphi = nint(nlats(i)*phi_over*sin(theta))
c
c                    phi_over is set in getgridph (typically = 2).
c     xnodesr()      radius associated with successive circles in templates
c     ngridc()       number of output points on successive circles in templates
c     nlow           index of lowest frequency sphere under consideration
c     ncur           index of highest frequency sphere under consideration
c     lmod           length of packed spherical harmonic model up to sphere
c                    ncur
c     isph_start()   array of length ngridr indicating where in the 
c                    modhat_sph vector the coefficients for the 
c                    corresponding sphere begin
c     nterms_sph()   array of length ngridr defining the orders of the
c                    spherical harmonic expansions on successive spheres.
c     oversamp       oversampling parameter for least squares solver
c     kord           interpolation order for least squares solver
C
C     OUTPUT: 
c
C     modelsph(:)    solution to least squares problem expressed in 
c                    spherical harmonics in successive spheres
c
C***********************************************************************
      implicit none
      integer verbose
      data verbose / 0 /
      integer itypep,nlow,ncur,lmod,imagesize,nimages,kord
      real *8 eps
      integer n_CTF,ld_CTF
      complex *16 CTF_p_(0:ld_CTF*n_CTF-1)
      integer ii,ir,n_A
      integer numit,niter
      integer is,jj,jstart,nout,noutmax,nquad
      integer nsphmax,nsphtot,nterms,ntmax
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 xnodesr(ncur)
      integer isph_start(ncur),nterms_sph(ncur)
      real *8 oversamp
      integer *4 n_alpha
      real *8 alpha2d(n_alpha,nimages)
      real *8, allocatable :: thetas(:)
      real *8, allocatable :: phis(:)
      complex *16 modelsph(0:lmod-1)
      complex *16 cslices(imagesize,nimages)
      complex *16, allocatable :: localp(:)
      complex *16, allocatable :: phiout(:)
      complex *16, allocatable :: ctfw(:)
      external multaha
      integer nouttmp
      character(len=1024) format_string
      if (verbose.gt.0) then
         write(6,*) '[entering rebuild_model_alpha2d]'
      end if
      if (verbose.gt.1) then
         write(6,*) 'imagesize: ',imagesize
         write(6,*) 'nimages: ',nimages
         write(6,*) 'icstart: ',(icstart(ii),ii=1,ncur)
         write(6,*) 'nlats: ',(nlats(ii),ii=1,ncur)
         write(6,*) 'itypep: ',itypep
         write(6,*) 'xnodesr: ',(xnodesr(ii),ii=1,ncur)
         write(6,*) 'ngridc: ',(ngridc(ii),ii=1,ncur)
         write(6,*) 'nlow: ',nlow
         write(6,*) 'ncur: ',ncur
         write(6,*) 'lmod: ',lmod
         write(6,*) 'isph_start: ',(isph_start(ii),ii=1,ncur)
         write(6,*) 'nterms_sph: ',(nterms_sph(ii),ii=1,ncur)
         write(6,*) 'oversamp: ',oversamp
         write(6,*) 'kord: ',kord
      end if
c
      nterms = nterms_sph(ncur)
      nsphmax = (nterms+1)**2
      if (verbose.gt.1) then
         write(6,*) ' nterms = ',nterms
         write(6,*) ' nsphmax = ',nsphmax
      end if
      allocate(localp(0:nsphmax))
c
      ntmax = nint(oversamp*nlats(ncur))
      if (verbose.gt.1) then
         write(6,*) ' oversamp = ',oversamp
         write(6,*) ' ntmax = ',ntmax
      end if
c
      noutmax = nimages*ngridc(ncur)
      if (verbose.gt.1) then
         write(6,*) ' noutmax = ',noutmax
      end if
      allocate(phiout(noutmax))
      allocate(thetas(noutmax))
      allocate(phis(noutmax))
      allocate(ctfw(noutmax))
c
c
c
      n_A = 0
      do ir = 1,ncur
         n_A = n_A + ngridc(ir)
      end do
      if (n_A.gt.imagesize) then
         write(6,*) 'Warning, imagesize: ',imagesize,' .neq. n_A: ',n_A
     $        ,' in rebuild_model_alpha2d'
      end if
c
      do is = nlow,ncur
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' is = ',is
         end if
         call get_lsqdata_alpha2d(cslices,imagesize,nimages,n_CTF,ld_CTF
     $        ,CTF_p_,icstart(is),xnodesr(is),ngridc(is),n_alpha
     $        ,alpha2d,phiout,thetas,phis,ctfw,nout)

         if (verbose.gt.1) then
c$$$            nouttmp = nimages*ngridc(is)
            nouttmp = 1*ngridc(is)
            write(format_string,'(A,I0,A)') '(A,',nouttmp*2,'F8.3)'
            write(6,format_string) 'phiout: ',(phiout(ir),ir=1,nouttmp)
            write(6,format_string) 'ctfw: ',(ctfw(ir),ir=1,nouttmp)
         end if

c
c        set up parameters for tensor product grid on sphere
c        to be used in least squares solve (for interpolation
c        to arbitrary points)
c
         nterms = nterms_sph(is)
         nquad = nint(oversamp*nlats(is))
         nsphtot = (nterms+1)**2
         numit = 1000
         niter = 0
         call lsqctfsolve(phiout,thetas,phis,nout,localp,nterms, nquad
     $        ,kord,ctfw,multaha,eps,numit,niter)
         jstart = isph_start(is)
         do jj = 0,nsphtot-1
            modelsph(jstart+jj) = localp(jj)
         enddo
      enddo
      if (verbose.gt.0) write(6,*)
     $     '[finished rebuild_model_alpha2d]'
      return
      end

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine get_lsqdata_alpha2d(ximagesk,imagesize,nimages,n_CTF
     $     ,ld_CTF,CTF_p_,icstart,xnodesr,ngridc,n_alpha,alpha2d
     $     ,phiout ,thetas,phis,ctfw,nout)
      implicit none
      integer verbose
      data verbose / 0 /
      integer imagesize,nimages,icstart,ngridc
      integer n_CTF,ld_CTF
      complex *16 CTF_p_(0:ld_CTF*n_CTF-1)
      integer nout,ii,iincr,ll
      complex *16 ximagesk(imagesize,nimages)
      complex *16 phiout(ngridc*nimages)
      complex *16 ctfw(ngridc*nimages)
      real *8 xnodesr
      integer *4 n_alpha
      include 'nalpha_define.f'
      real *8 alpha2d(n_alpha,nimages)
      real *8 thetas(ngridc*nimages)
      real *8 phis(ngridc*nimages)
      real *8 pi,ctheta,cthetause
      real *8 phi,phiuse,rgammause,rho,rkx,rky,rkz,rone
      real *8 delta_x,delta_y,gamma_z,l2_norm
      real *8, allocatable :: xc(:)
      real *8, allocatable :: yc(:)
      real *8, allocatable :: zc(:)
      complex *16, allocatable ::  M_p_1_(:)
      complex *16, allocatable ::  CTf_p_1_(:)
      integer nw,nctf
      if (verbose.gt.0) write(6,*)
     $     '[entering get_lsqdata_alpha2d]'
c
      allocate(xc(ngridc))
      allocate(yc(ngridc))
      allocate(zc(ngridc))
      allocate(M_p_1_(0:ngridc-1))
      allocate(CTF_p_1_(0:ngridc-1))
c
      pi=4.0d0*datan(1.0d0)
      iincr = 0
      do ii = 1,nimages
         nctf = nint(alpha2d(1+nalpha_ctf_ind,ii))
         nctf = max(0,min(n_CTF-1,nctf))
         cthetause = dcos(alpha2d(1+nalpha_polar_a,ii))
         phiuse = alpha2d(1+nalpha_azimu_b,ii)
         rgammause = 0.0d0
         rone = 1.0d0
         call mkonegreatcircle(rone,cthetause,phiuse,rgammause,
     1        ngridc,xc,yc,zc)
c
         if (n_alpha.ge.1+nalpha_l2_norm) then
            l2_norm = alpha2d(1+nalpha_l2_norm,ii)
         else
            l2_norm = 1.0d0
         end if
            
         if (dabs(l2_norm).lt.1.0d-15) then
            l2_norm=1.0d0
         end if
         do nw = 0,ngridc-1
            M_p_1_(nw) = ximagesk(icstart+nw,ii)/l2_norm
            CTF_p_1_(nw) = CTF_p_(icstart+nw-1 + nctf*ld_CTF)
         enddo
c
         gamma_z = -alpha2d(1+nalpha_gamma_z,ii)
         call rotate_p_to_p_1(ngridc,M_p_1_,+gamma_z,M_p_1_)
         call rotate_p_to_p_1(ngridc,CTF_p_1_,+gamma_z,CTF_p_1_)
         delta_x = -alpha2d(1+nalpha_delta_x,ii)
         delta_y = -alpha2d(1+nalpha_delta_y,ii)
         call transf_p_to_p_1(xnodesr,ngridc,M_p_1_,+delta_x,+delta_y
     $        ,M_p_1_)
c
         do ll = 1,ngridc
            rkx = xc(ll)
            rky = yc(ll)
            rkz = zc(ll)
            ctheta = rkz
            phi = datan2(rky,rkx)
            if (phi.lt.0.0d0) phi = phi + 2*pi
            iincr = iincr+1
            phis(iincr) = phi
            rho = dsqrt(rkx**2+ rky**2)
ccc            thetas(iincr) = dacos(ctheta)
            thetas(iincr) = datan2(rho,ctheta)
            phiout(iincr) = M_p_1_(ll-1)
            ctfw(iincr) = CTF_p_1_(ll-1)
         enddo
      enddo
      nout = ngridc*nimages
      if (verbose.gt.0) write(6,*)
     $     '[finished get_lsqdata_alpha2d]'
      return
      end
