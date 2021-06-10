C***********************************************************************
      subroutine test_resid_alpha2d(cslices,imagesize,n_image,n_CTF
     $     ,ld_CTF,CTF_p_,n_alpha,alpha2d,icstart,nlats,itypep ,xnodesr
     $     ,ngridc,nlow,ncur,n_modelsph,isph_start ,nterms_sph ,oversamp
     $     ,kord ,eps,modelsph,ctransf,mslices,residua,n_loading
     $     ,loading)
C***********************************************************************
C     Builds new model (in spherical harmonic basis) 
C     from slices (Fourier transforms of images) with image parameters 
C     given by alpha2d array.
C---------------------------------------------------------------------
C     INPUT:
C
C     cslices        Fourier transforms of images on polar grids
C     imagesize      dimension of images on quasiunform polar grid
C     n_image        number of images
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
c     n_modelsph     length of packed spherical harmonic model up to sphere
c                    ncur
c     isph_start()   array of length ngridr indicating where in the 
c                    modhat_sph vector the coefficients for the 
c                    corresponding sphere begin
c     nterms_sph()   array of length ngridr defining the orders of the
c                    spherical harmonic expansions on successive spheres.
c     oversamp       oversampling parameter for least squares solver
c     kord           interpolation order for least squares solver
C
C     modelsph(:)    solution to least squares problem expressed in 
c                    spherical harmonics in successive spheres
c     n_loading      integer: number of loading vectors per image
c
C     OUTPUT: 
c     ctransf(:)     images transformed
c     mslices(:)     model sliced
c     residua(:)     residuals = ctransf - mslices
c     loading(:)     loading vectors
c
C***********************************************************************
      implicit none
      integer verbose
      data verbose / 0 /
      integer itypep,nlow,ncur,imagesize,n_image,nimage,kord,nctf
      real *8 eps
      integer n_CTF,ld_CTF,n_loading
      complex *16 CTF_p_(0:ld_CTF*n_CTF-1)
      integer ii,ir,n_A,ntmp
      integer numit,niter
      integer is,jj,jstart,noutmax,nquad
      integer nsphmax,nsphtot,nterms,ntmax
      integer icstart(ncur),nlats(ncur),ngridc(ncur)
      real *8 xnodesr(ncur)
      integer isph_start(ncur),nterms_sph(ncur)
      real *8 oversamp
      integer *4 n_alpha
      real *8 alpha2d(n_alpha,n_image)
      include 'nalpha_define.f'
      real *8, allocatable :: thetas(:,:)
      real *8, allocatable :: phis(:,:)
      complex *16 modelsph(0:n_modelsph-1)
      complex *16 cslices(imagesize,n_image)
      complex *16 ctransf(imagesize,n_image)
      complex *16 mslices(imagesize,n_image)
      complex *16 residua(imagesize,n_image)
      complex *16 loading(0:n_loading*n_image-1)
      complex *16, allocatable :: localp(:)
      complex *16, allocatable :: ctfw(:,:)
      integer *4 , allocatable :: nout_(:)
      integer nloading,n_modelsph,nmodelsph
      complex *16, allocatable :: G_(:)
      complex *16, allocatable :: H_(:)
      complex *16 HH,HG,RR
      integer n_iteration,niteration
      external multaha
      integer nouttmp
      character(len=1024) format_string
      if (verbose.gt.0) then
         write(6,*) '[entering test_resid_alpha2d]'
      end if
      if (verbose.gt.1) then
         write(6,*) 'imagesize: ',imagesize
         write(6,*) 'n_image: ',n_image
         write(6,*) 'icstart: ',(icstart(ii),ii=1,ncur)
         write(6,*) 'n_CTF: ',n_CTF
         write(6,*) 'ld_CTF: ',ld_CTF
         write(6,*) 'n_alpha: ',n_alpha
         write(6,*) 'nlats: ',(nlats(ii),ii=1,ncur)
         write(6,*) 'itypep: ',itypep
         write(6,*) 'xnodesr: ',(xnodesr(ii),ii=1,ncur)
         write(6,*) 'ngridc: ',(ngridc(ii),ii=1,ncur)
         write(6,*) 'nlow: ',nlow
         write(6,*) 'ncur: ',ncur
         write(6,*) 'n_modelsph: ',n_modelsph
         write(6,*) 'isph_start: ',(isph_start(ii),ii=1,ncur)
         write(6,*) 'nterms_sph: ',(nterms_sph(ii),ii=1,ncur)
         write(6,*) 'oversamp: ',oversamp
         write(6,*) 'kord: ',kord
         write(6,*) 'n_loading: ',n_loading
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
      noutmax = n_image*ngridc(ncur)
      if (verbose.gt.1) then
         write(6,*) ' noutmax = ',noutmax
      end if
      allocate(thetas(imagesize,n_image))
      allocate(phis(imagesize,n_image))
      allocate(ctfw(imagesize,n_image))
      allocate(nout_(ncur))
c
c
c
      do ir = 1,ncur
         nout_(ir) = 0
      end do
      n_A = 0
      do ir = 1,ncur
         n_A = n_A + ngridc(ir)
      end do
      if (n_A.gt.imagesize) then
         write(6,*) 'Warning, imagesize: ',imagesize,' .neq. n_A: ',n_A
     $        ,' in test_resid_alpha2d'
      end if

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating ctransf'
      end if !if (verbose.gt.0) then
      do nimage=1,n_image
         nctf = nint(alpha2d(1+nalpha_ctf_ind,nimage))
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0)') ' nimage: ' , nimage , ' nctf: '
     $           , nctf, 'imagesize: ' , imagesize
         end if
         do is = nlow,ncur
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I0)') ' is = ' , is , ' icstart(is) ' ,
     $              icstart(is)
            end if
            call get_lsqdata_alpha2d(cslices(1,nimage),imagesize,1,1
     $           ,ld_CTF,CTF_p_(nctf*ld_CTF),icstart(is),xnodesr(is)
     $           ,ngridc(is),n_alpha,alpha2d(1,nimage)
     $           ,ctransf(icstart(is),nimage),thetas(icstart(is),nimage)
     $           ,phis(icstart(is),nimage),ctfw(icstart(is),nimage)
     $           ,nout_(is))
         enddo                  ! do is = nlow,ncur
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,' , imagesize ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'cslices: ' , (cslices(ntmp,nimage)
     $           ,ntmp=1,imagesize)
            write(6,format_string) 'ctransf: ' , (ctransf(ntmp,nimage)
     $           ,ntmp=1,imagesize)
         end if !if (verbose.gt.2) then
      enddo                     ! do nimage=1,n_image

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating mslices'
      end if !if (verbose.gt.0) then
      do nimage=1,n_image
         nctf = nint(alpha2d(1+nalpha_ctf_ind,nimage))
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0)') ' nimage: ' , nimage , ' nctf: '
     $           , nctf, ' imagesize: ' , imagesize
         end if
         do is = nlow,ncur
            if (verbose.gt.1) then
               write(6,'(A,I0,A,I0)') ' is = ' , is , ' icstart(is) ' ,
     $              icstart(is)
            end if
            nterms = nterms_sph(is)
            nquad = nint(oversamp*nlats(is))
c$$$            nsphtot = (nterms+1)**2
c$$$            jstart = isph_start(is)
c$$$            do jj = 0,nsphtot-1
c$$$               localp(jj) = modelsph(jstart+jj)
c$$$            enddo
c$$$            call test_resid_multa(thetas(icstart(is),nimage)
c$$$     $           ,phis(icstart(is),nimage),nout_(is),localp,nterms,nquad
c$$$     $           ,kord ,ctfw(icstart(is),nimage),mslices(icstart(is),nimage))
            call test_resid_multa(thetas(icstart(is),nimage)
     $           ,phis(icstart(is),nimage),nout_(is)
     $           ,modelsph(isph_start(is)),nterms,nquad ,kord
     $           ,ctfw(icstart(is),nimage),mslices(icstart(is),nimage))
         enddo                  ! do is = nlow,ncur
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,' , imagesize ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'mslices: ' , (mslices(ntmp,nimage)
     $           ,ntmp=1,imagesize)
         end if !if (verbose.gt.2) then
      enddo                     ! do nimage=1,n_image      

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating residua'
      end if !if (verbose.gt.0) then
      do nimage=1,n_image
         do jj=1,imagesize
            residua(jj,nimage) = cmplx(0.0d0,0.0d0)
         enddo !do jj=1,imagesize
         do is=3,ncur
            do jj=0,ngridc(is)-1
               residua(icstart(is)+jj,nimage) = ctransf(icstart(is)+jj
     $              ,nimage)-mslices(icstart(is)+jj,nimage)
            enddo !do jj=0,ngridc(is)-1
         enddo !         do is=3,ncur
      enddo                     ! do nimage=1,n_image      

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating loading'
      end if !if (verbose.gt.0) then
      do nloading=0,n_loading*n_image-1
         loading(nloading) = cmplx(rand()-0.5d0,0.0d0)
      enddo !do nloading=0,n_loading*n_image-1
      allocate(G_(0:n_loading*n_modelsph-1))
      allocate(H_(0:n_modelsph-1))
      n_iteration = 3
      do niteration=0,n_iteration-1
         write(6,*) ' niteration ' , niteration
         if (niteration.eq.0) then
            do nloading=0,n_loading-1
               do nmodelsph=0,n_modelsph-1
                  G_(nmodelsph+nloading*n_modelsph) = cmplx(rand()-0.5
     $                 ,0.0d0)
               enddo            !do nmodelsph=0,n_modelsph-1
               call normalize_c16(n_modelsph,G_(nloading*n_modelsph))
            enddo               !do nloading=0,n_loading*n_modelsph-1
         end if                 !if (niteration.eq.0) then
         if (niteration.gt.0) then
            call cl1_c16(n_modelsph*n_loading,G_)
            do nimage=1,n_image
               do is = nlow,ncur
                  nterms = nterms_sph(is)
                  nquad = nint(oversamp*nlats(is))
                  call test_resid_multah(thetas(icstart(is),nimage)
     $                 ,phis(icstart(is),nimage),nout_(is)
     $                 ,residua(icstart(is),nimage),nterms,nquad
     $                 ,kord ,ctfw(icstart(is),nimage)
     $                 ,H_(isph_start(is)))
               enddo            ! do is = nlow,ncur
               do nloading=0,n_loading-1
                  call af1_c16(n_modelsph,loading(nloading + n_loading
     $                 *(nimage-1)),cmplx(0.0d0,0.0d0),H_,G_(0 +
     $                 nloading*n_modelsph))
               enddo            !do nloading=0,n_loading-1
            enddo               ! do nimage=1,n_image
            call gramschmidt_c16(n_modelsph,n_loading,G_)
         end if                 !if (niteration.gt.0) then
         do nimage=1,n_image
            do is = nlow,ncur
               nterms = nterms_sph(is)
               nquad = nint(oversamp*nlats(is))
               call test_resid_multah(thetas(icstart(is),nimage)
     $              ,phis(icstart(is),nimage),nout_(is)
     $              ,residua(icstart(is),nimage),nterms,nquad
     $              ,kord ,ctfw(icstart(is),nimage)
     $              ,H_(isph_start(is)))
            enddo               ! do is = nlow,ncur
            call dot_c16(n_modelsph,H_,H_,HH)
c$$$            call dot_c16(imagesize,residua(1,nimage),residua(1,nimage)
c$$$     $           ,RR)
            do nloading=0,n_loading-1
               call dot_c16(n_modelsph,H_,G_(0+nloading*n_modelsph),HG)
               loading(nloading + n_loading*(nimage-1)) = dreal(HG)
     $              /max(1.0d-15,dreal(HH))
c$$$               loading(nloading + n_loading*(nimage-1)) = dreal(HG)
c$$$     $              /max(1.0d-15,dreal(RR))
            enddo               !do nloading=0,n_loading-1
         enddo                  ! do nimage=1,n_image
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,I5,' , n_loading ,
     $           '(2F10.4,1X))'
            do nimage=1,n_image
               write(6,format_string) 'loading: image ' , nimage-1 ,
     $              (loading(nloading + n_loading*(nimage -1)),nloading
     $              =0,n_loading-1)
            enddo               !do nimage=1,n_image
         end if                 !if (verbose.gt.2) then

         if (verbose.gt.3) then
            write(format_string,'(A,I0,A)') '(' , n_loading ,
     $           '(2F10.4,1X))'
            do nimage=1,n_image
               write(6,format_string) (loading(nloading + n_loading
     $              *(nimage -1)),nloading=0,n_loading-1)
            enddo               !do nimage=1,n_image
         end if                 !if (verbose.gt.3) then

      enddo                     !do niteration=0,n_iteration-1
      deallocate(H_)
      deallocate(G_)

      if (verbose.gt.0) write(6,*)
     $     '[finished test_resid_alpha2d]'
      return
      end

