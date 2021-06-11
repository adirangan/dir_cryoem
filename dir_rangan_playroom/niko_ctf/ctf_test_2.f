      implicit none
c
      real *8 pi
      complex *16 eye

      integer  ngrid
      integer  ngridr, ngridt, ngridp, ngridc
c
      real *8 a, h
c
      real *8, allocatable :: x1(:)
      real *8, allocatable :: xnodesr(:), wtsr(:)
      real *8, allocatable :: gauss_xnodesr(:),gauss_wtsr(:)
      real *8, allocatable :: xnodesc(:), wtsc(:)

      real *8 rmax
      integer nexpim

c     CTF parameters
c     spherical aberation of the lens
      real *8 Cs
      real *8, allocatable :: SphericalAberration(:)
c     voltage
      real *8 Kv
      real *8,allocatable :: Voltage(:)
c     amplitude contrast
      real *8 AmpCnst
      real *8, allocatable :: AmplitudeContrast(:)
c     magnification of microscope
      real *8 xmag
      real *8, allocatable :: Magnification(:)
c     pixel size of scanner in microns
      real *8 dstep
      real *8, allocatable :: DetectorPixelSize(:)
c     defocus (in Angstroms) and angle of astigmatism
      real *8 df1,df2,angast
      real *8, allocatable :: DefocusU(:)
      real *8, allocatable :: DefocusV(:)
      real *8, allocatable :: DefocusAngle(:)

c     electron wavelength in Angstoms
      real *8 lambda
c     weights for the amplitude and phase contrasts in CTF
      real *8 w1,w2
c     pixel size of the scanner in physical space (not magnified) in Angstroms
      real *8 stepr
c     box size in physical space
      real *8 D
      real *8 angspt,rad
c      
      real *8 thetatr
      real *8 ctfv,CTF,ctfv1,ctf1,ctfval

      complex *16, allocatable :: ctf1d(:)
      complex *16, allocatable :: ctfw(:,:)
      complex *16, allocatable :: ctfw1(:,:)
      real *8, allocatable :: envelop(:)

      real *8 sm
      integer i,j,k,l,m
      real *8 ix,iy
      integer ll,rr,cn
      real *8 t1,t2
      character(len=100) :: outdir
      character(len=100) :: fffile

cccccccccccccccccccccccccccccccccccccccccccccccccc


      data eye/(0.0d0,1.0d0)/
      pi=4.0*atan(1.0d0)
c     get the output directory from the command line
      call get_command_argument(1,outdir)
      outdir = trim(adjustl(outdir)) // '/'


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize the random seed
      call srand(time())
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     size of the box in which the real-space function lives
      a = 1.0d0
c     number of grid points for the real-space function
      ngrid = 256
c     how far to go out in frequency, i.e. max of radial component
c      rmax = 300
      rmax = 300
c     number of radial quadrature nodes
c      ngridr = 1000
      ngridr = 100
c     number of quadrature nodes in theta - should be >rmax by Nyquist
c      ngridt = 450
      ngridt = 450
c     number of quadrature nodes in phi - should be >2rmax by Nyquist
c      ngridp = 900
      ngridp = 900
c     number of points on the equatorial circles
c      ngridc = 1000
      ngridc = 100
c     number of experimental images
      nexpim = 10
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      allocate(x1(ngrid))
      allocate(xnodesr(ngridr))
      allocate(wtsr(ngridr))
      allocate(xnodesc(ngridc))
      allocate(wtsc(ngridc))

      allocate(Voltage(nexpim))
      allocate(DefocusU(nexpim))
      allocate(DefocusV(nexpim))
      allocate(DefocusAngle(nexpim))
      allocate(SphericalAberration(nexpim))
      allocate(DetectorPixelSize(nexpim))
      allocate(Magnification(nexpim))
      allocate(AmplitudeContrast(nexpim))
      allocate(ctf1d(ngridr))
      allocate(ctfw(ngridc,ngridr))
      allocate(ctfw1(ngridc,ngridr))
      allocate(envelop(ngridr))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up some CTF parameters

c     spherical aberation of the lens in mm
      Cs = 2.70d0
c     into Angstroms
      Cs=Cs*(10.0d0**7.0d0)
c     voltage in Kvolts
      Kv = 300.0d0
c     into Volts
      Kv=Kv*1000.0 
c     Amplitude Contrast
      AmpCnst = 0.07
c     magnification of the microscope
      xmag = 46627.0d0
c     pixel size on scanner in microns
      dstep = 5.0d0 
c     into Angstroms
      dstep = dstep*(10.0d0**4.0d0)

c     electron wavelength in Angstroms
c      lambda=12.3d0/sqrt(Kv+Kv**2.0d0/(10.0d0**6.0d0))
      lambda = 12.2643247/sqrt(Kv+Kv**2*0.978466d-6) 

      write(6,*) 'lambda ', lambda
c     weights for the amplitude and phase contrasts in CTF
      w1=sqrt(1.0d0-AmpCnst**2)
      w2=AmpCnst
c     defocus values (probably in Angstroms)
c      df1 =13682.0d0 
c      df2 = 13600.0d0
c      df2 = 13600.0d0
      df1 = 20000.0d0 
      df2 = 1000.0d0
c     angle of astigmatism
      angast = 45.0d0
c     into radians
      angast=angast*pi/180.0d0

c     pixel size of the scanner in physical space (not magnified) in Angstroms
      stepr = dstep/xmag

c     box size in Angstroms
      D = stepr*ngrid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Niko's CTF
c     angle of deflection?
      thetatr = lambda/D
      do l=1,20
         m = l+20
         write(6,*) l,m
         ix = l
         iy = m
         call niko_ctf(Cs,lambda,w1,w2,df1,df2,angast,
     1        thetatr,ix,iy,ctfv)
         angspt=atan2(iy,ix)
         rad = m**2+l**2
         rad = sqrt(rad)
         call niko_ctf1(Cs,lambda,w1,w2,df1,df2,angast,
     1        D,angspt,rad,ctfv1)
         call ctf_asym(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ctfval)
         write(6,*) ctfv,ctfv1,ctfval
         enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     r and theta grids in 2D
      do i = 1,ngridr
         xnodesr(i) = rmax*(i)/ngridr
      enddo
      do i = 1,ngridc
         xnodesc(i) = 2.0d0*pi*i/ngridc
      enddo

      call get_ctf_asym_full(Cs,lambda,AmpCnst,df1,df2,angast,D,
     1     ngridr,xnodesr,ngridc,xnodesc,ctfw)

      call our_envelop_fxn(lambda,D,ngridr,xnodesr,envelop)

c     write it out for plotting
      open(20,FILE='ctfw');
      do ll = 1,ngridc
         write(20,*) (real(ctfw(ll,rr)),rr=1,ngridr)
      enddo
      close(20)

      do rr = 1,ngridr
      do ll = 1,ngridc
         ctfw1(ll,rr) = ctfw(ll,rr)*envelop(rr)
      enddo
      enddo

c     write it out for plotting
      open(20,FILE='ctfw1');
      do ll = 1,ngridc
         write(20,*) (real(ctfw1(ll,rr)),rr=1,ngridr)
      enddo
      close(20)


      stop
      end


