!> Doxygen comment: ;\n
!>     Utilities creating Cartesian and spherical grids for ;\n
!>     physical and Fourier space, respectively. ;\n
!> Doxygen comment: ;\n
C     Utilities creating Cartesian and spherical grids for
C     physical and Fourier space, respectively.
C
C
      subroutine writevol(iout,ngrid,x1,x,y,z,ff)
c
c     dump abs(ff) to output file
c
      implicit real *8 (a-h,o-z)
      dimension x1(ngrid)
      dimension x(ngrid,ngrid,ngrid)
      dimension y(ngrid,ngrid,ngrid)
      dimension z(ngrid,ngrid,ngrid)
      complex *16 ff(ngrid,ngrid,ngrid)
c
      do i = 1,ngrid
      do j = 1,ngrid
      do k = 1,ngrid
         write(iout,*) abs(ff(i,j,k))
      enddo
      enddo
      enddo
      return
      end
c
c
      subroutine mkphysgrid(a,ngrid,h,xg,yg,zg)
c
c     create uniform ngrid x ngrid x ngrid grid on [-a,a]^3.
c
c     INPUT:
c
c     a       box dimension
c     ngrid   number of grid points in each dimension
c   
c     OUTPUT:
c
c     h         grid spacing
c     xg        x component of all grid points
c     yg        y component of all grid points
c     zg        z component of all grid points
c     
      implicit none
      integer ngrid
      real *8 a, h
      real *8 xg(ngrid,ngrid,ngrid)
      real *8 yg(ngrid,ngrid,ngrid)
      real *8 zg(ngrid,ngrid,ngrid)

      integer i, j, k
      real *8 x, y, z
c     
c     create physical space grid 
c     grab data in physical space on given grid.
c
      h = 2*a/ngrid
      do i = 1,ngrid
         x = -a + i*h
         do j = 1,ngrid
            y = -a + j*h
            do k = 1,ngrid
               z = -a + k*h
               xg(i,j,k) = x
               yg(i,j,k) = y
               zg(i,j,k) = z
c$$$               write(6,*) 'i,j,k: ',i,j,k,'; x,y,z: ',x,y,z
            enddo
         enddo
      enddo
c
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fgauss(x,y,z,x0y0z0,ns,sigma,ff)
c
c     INPUT:
c
c     x,y,z   target location
c     x0y0z0  location of Gaussian sources
c     ns      number of Gaussian sources
c     sigma   number variance of all Gaussian sources
c   
c     OUTPUT:
c
c     ff        complex function value at corresponding grid pt.
ccccccccccccccccccccccccccccccccccccccccccccccccc     
c
      implicit real *8 (a-h,o-z)
      real *8 x0y0z0(3,ns)
      complex *16 ff
c
c$$$      write(6,*) 'Entering fgauss'
      ff = 0.0d0
      do i = 1,ns
c$$$         write(6,*) 'i: ',i
         xx = x - x0y0z0(1,i)
         yy = y - x0y0z0(2,i)
         zz = z - x0y0z0(3,i)
         rr = xx*xx+yy*yy+zz*zz
         ff = ff + dexp(-rr/(2*sigma*sigma))
      enddo
c$$$      write(6,*) 'finishing fgauss'
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fgauss_many_sigma(x,y,z,x0y0z0,ns,sigmas,ff)
c
c      INPUT:
c
c     x,y,z   target location
c     x0y0z0  location of Gaussian sources
c     ns      number of Gaussian sources
c     sigmas   different variance for each Gaussian sources
c   
c     OUTPUT:
c
c     ff        complex function value at corresponding grid pt.
ccccccccccccccccccccccccccccccccccccccccccc     
c
      implicit none
      integer ns
      real *8 x,y,z
      real *8 x0y0z0(3,ns)
      real *8 sigmas(ns)
      complex *16 ff
c
      integer i
      real *8 sigma,rr,xx,yy,zz
c     
      ff = 0.0d0
      do i = 1,ns
         sigma = sigmas(i)
         xx = x - x0y0z0(1,i)
         yy = y - x0y0z0(2,i)
         zz = z - x0y0z0(3,i)
         rr = xx*xx+yy*yy+zz*zz
         ff = ff + dexp(-rr/(2*sigma*sigma))
      enddo
      return
      end
c
      subroutine fgaussft_many_sigma(kx,ky,kz,x0y0z0,ns,sigmas,ffhat)
c
c     Compute the Fourier transform of Gaussians (as defined in fgauss)
c
c     INPUT:
c
c     kx,ky,kz target location (in Fourier space)
c     x0y0z0   location of Gaussian sources
c     ns       number of Gaussian sources
c     sigmas    number variance for every Gaussian
c   
c     OUTPUT:
c
c     ffhat    complex function value at corresponding grid pt.
c     
c
      implicit none
      integer ns
      real *8 kx, ky, kz, sigma
      real *8 x0y0z0(3,ns),sigmas(ns)
      complex *16 ffhat
c
      integer i
      real *8 arg1, arg2
      complex *16 eye
      real *8 pi
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      ffhat = 0.0d0
      do i = 1,ns
         sigma = sigmas(i)
         arg1 = kx*x0y0z0(1,i) + ky*x0y0z0(2,i) + kz*x0y0z0(3,i)
         arg2 = -(kx*kx + ky*ky + kz*kz) * sigma*sigma/2
         ffhat = ffhat + cdexp(eye*arg1 + arg2)*(2.0d0*pi)*
     1        dsqrt(2.0d0*pi)*sigma**3
      enddo
      return
      end
c
      subroutine evalgrid_many_sigma(x,y,z,n1,n2,n3,x0y0z0,ns,sigmas,
     1     fneval,ff)
c
      implicit none
      integer ns,n1,n2,n3
      real *8 x(n1,n2,n3), y(n1,n2,n3), z(n1,n2,n3)
      real *8 x0y0z0(3,ns),sigmas(ns)
      complex *16 ff(n1,n2,n3)
      external fneval
c
      integer i1, i2, i3
c
      do i3 = 1,n3
         do i2 = 1,n2
            do i1 = 1,n1
               call fneval(x(i1,i2,i3), y(i1,i2,i3), z(i1,i2,i3),
     1              x0y0z0, ns, sigmas, ff(i1,i2,i3))
            enddo
         enddo
      enddo
c
      return
      end
c
c     
      subroutine fgaussft(kx,ky,kz,x0y0z0,ns,sigma,ffhat)
c
c     Compute the Fourier transform of Gaussians (as defined in fgauss)
c
c     INPUT:
c
c     kx,ky,kz target location (in Fourier space)
c     x0y0z0   location of Gaussian sources
c     ns       number of Gaussian sources
c     sigma    number variance of all Gaussian sources
c   
c     OUTPUT:
c
c     ffhat    complex function value at corresponding grid pt.
c     
      implicit none
      integer ns
      real *8 kx, ky, kz, sigma
      real *8 x0y0z0(3,ns)
      complex *16 ffhat
c
      integer i
      real *8 arg1, arg2
      complex *16 eye
      real *8 pi
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
      ffhat = 0.0d0
      do i = 1,ns
         arg1 = kx*x0y0z0(1,i) + ky*x0y0z0(2,i) + kz*x0y0z0(3,i)
         arg2 = -(kx*kx + ky*ky + kz*kz) * sigma*sigma/2
         ffhat = ffhat + cdexp(eye*arg1 + arg2)
      enddo
      ffhat = ffhat * (2.0d0*pi)*dsqrt(2.0d0*pi)*sigma**3
      return
      end
c 
c
      subroutine evalgrid(x,y,z,n1,n2,n3,x0y0z0,ns,sigma,fneval,ff)
c
      implicit none
      integer ns,n1,n2,n3
      real *8 x(n1,n2,n3), y(n1,n2,n3), z(n1,n2,n3), sigma
      real *8 x0y0z0(3,ns)
      complex *16 ff(n1,n2,n3)
      complex *16 ftmp
      external fneval
c 
      integer i1, i2, i3
c$$$      write(6,*) 'entering evalgrid'
c$$$      write(6,*) 'ns: ',ns
c$$$      do i1=1,3
c$$$         write(6,*) (x0y0z0(i1,i2),i2=1,ns)
c$$$      enddo
c$$$      write(6,*) 'sigma: ',sigma
      call fneval(0d0,0d0,0d0,x0y0z0,32,1.0d0/32.0d0,ftmp)
c
      do i3 = 1,n3
         do i2 = 1,n2
            do i1 = 1,n1
c$$$               write(6,*) 'i1,i2,i3: ',i1,i2,i3
c$$$               write(6,*) 'x: ',x(i1,i2,i3)
c$$$               write(6,*) 'y: ',x(i1,i2,i3)
c$$$               write(6,*) 'z: ',x(i1,i2,i3)
c$$$               write(6,*) 'S: ',ff(i1,i2,i3)
c$$$               write(6,*) 'Calling fneval'
               call fneval(x(i1,i2,i3),y(i1,i2,i3),z(i1,i2,i3),
     1               x0y0z0,ns,sigma,ff(i1,i2,i3))
            enddo
         enddo
      enddo

c$$$      write(6,*) 'finished evalgrid'
c
      return
      end
c
c
c
      subroutine l2error_c(a1, a2, n, l2err, l2)
c
c     Compute L2 norm of a1-a2
c
c     INPUT:
c     a1      First array of complex numbers
c     a2      Second array of complex numbers
c     n       Number of elements
c
c     OUTPUT:
c     l2err   L2 norm of (a1 - a2)
c     l2      L2 norm of a1
c
      implicit none
      complex *16 a1(n), a2(n)
      integer n
      real *8 l2err, l2
c
      integer i
c
      l2 = 0.0d0
      l2err = 0.0d0
      do i=1,n
         l2 = l2 + abs(a1(i))**2
         l2err = l2err + abs(a1(i) - a2(i))**2
      enddo
c
      l2 = dsqrt(l2)
      l2err = dsqrt(l2err)
c 
      return
      end
c
c

C
C***********************************************************************
      subroutine sph_gridconfig1(ngridr, sph_size, sph_nterms, 
     1     sph_sizes,sph_startix)
C***********************************************************************
c     Set spherical harmonic expansion degrees according to radius
c     and determine size of various arrays
C---------------------------------------------------------------------
C     INPUT:
c
c     ngridr 		: number of radial grid points
c
C---------------------------------------------------------------------
C     OUTPUT:
C
c     sph_size:            size of the array required to store all 
c                          spherical harmonics expansions
c     sph_sizes(ngridr):   size of spherical harmonic expansion for 
c                          radius r(ngridr)
c     sph_nterms(ngridr):  spherical harmonic degree for radius r(ngridr)
c     sph_startix(ngridr): starting index for expansion at radius r(ngridr)
C***********************************************************************
      implicit none
      integer ngridr
      integer sph_size, sph_sizes(ngridr), sph_startix(ngridr)
      integer sph_nterms(ngridr)
c
      integer i, n
c
      sph_size = 0
      do i = 1,ngridr
c        hard code for all expansions
c         n = 30
c     size of expansions varies with i
c         n = max(5,i)
         n = i+2
         sph_nterms(i) = n
         sph_sizes(i) = (n+1)*(n+1)
         sph_startix(i) = sph_size

         sph_size = sph_size + (n+1)*(n+1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sph_gridconfig(ngridr, sph_size, sph_nterms, 
     1     sph_sizes,sph_startix, xnodesr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Configure the grid of spherical harmonics expansions by radius
c
c     Inputs:
c     ngridr: number of radial grid points
c
c     Outputs:
c     sph_size:            total size of the array required to store the spherical harmonics expansions
c     sph_sizes(ngridr):   individual sizes of spherical hardmonics expansion by radius
c     sph_nterms(ngridr):  number of terms of the expansion by radius (nterms)
c     sph_startix(ngridr): starting index for an expansion by radius
c
      implicit none
      integer ngridr
      integer sph_size, sph_sizes(ngridr), sph_startix(ngridr)
      integer sph_nterms(ngridr)
      real *8 xnodesr(ngridr)
c
      integer i, n
c
      sph_size = 0
      do i = 1,ngridr
c        hard code for all expansions
c         n = 30
c     size of expansions varies with i
         n = floor(xnodesr(i))+2
c         write(6,*) i, xnodesr(i),n
         sph_nterms(i) = n
         sph_sizes(i) = (n+1)*(n+1)
         sph_startix(i) = sph_size

         sph_size = sph_size + (n+1)*(n+1)
      enddo

      return
      end
C
c
c
c
C
C
      SUBROUTINE PRINM(MPOLE,NTERMS)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,-nterms:NTERMS)
      real *8 MPOLE2(0:NTERMS,0:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-nterms:NTERMS)
      INTEGER NTERMS
C
C     print out coefficients of multipole expansion
C
1000  FORMAT(6E12.5)
1001  FORMAT(/)
      DO 100 L = 0,NTERMS
         WRITE(6,1000)(MPOLE(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
100   CONTINUE
        return
C
C
C
C
      ENTRY PRINM0(MPOLE0,NTERMS)
      DO 200 L = 0,NTERMS
         WRITE(6,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(13,1000)(MPOLE0(L,M),M=-L,L)
         WRITE(6,1001)
         WRITE(13,1001)
200   CONTINUE
      RETURN
C
C
      ENTRY PRINM2(MPOLE2,NTERMS)
      DO L = 0,NTERMS
         WRITE(6,1000)(MPOLE2(L,M),M=0,L)
         WRITE(13,1000)(MPOLE2(L,M),M=0,L)
         WRITE(6,1001)
         WRITE(13,1001)
      ENDDO
c
c
      RETURN
      end
c
c
c
c**********************************************************************
      subroutine hpotfld3dall(iffld,sources,charge,ns,
     1                   target,wavek,pot,fld)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a collection of charges at 
c     SOURCE(3,ns). The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = dexp(i*k*r)/r
c		fld = grad(pot)
c
c     INPUT:
c
c     iffld  (integer *4)      : flag for computing gradient
c	                 	   ffld = 0 -> don't compute 
c		                   ffld = 1 -> do compute 
c     sources(3,ns)  (real *8) : location of the sources
c     charge(ns) (complex *16) : charge strength
c     target(3)  (real *8)     : location of the target
c     wavek  (complex *16)     : helmholtz parameter
c
c     OUTPUT:
c
c     pot   (real *8)        : calculated potential
c     fld   (real *8)        : calculated gradient
c
      real *8 sources(3,ns),target(3)
      complex *16 wavek,pot,fld(3),potloc,fldloc(3)
      complex *16 h0,h1,cd,eye,z,ewavek
      complex *16 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      do i = 1,ns
         call hpotfld3d(iffld,sources(1,i),charge(i),target,wavek,
     1        potloc,fldloc)
         pot = pot + potloc
         fld(1) = fld(1) + fldloc(1)
         fld(2) = fld(2) + fldloc(2)
         fld(3) = fld(3) + fldloc(3)
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine hpotfld3d(iffld,source,charge,target,wavek,pot,fld)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT and field FLD
c     at the target point TARGET, due to a charge at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = dexp(i*k*r)/r
c		fld = grad(pot)
c
c     INPUT:
c
c     iffld  (integer *4)    : flag for computing gradient
c	                 	ffld = 0 -> don't compute 
c		                ffld = 1 -> do compute 
c     source(3)  (real *8)   : location of the source 
c     charge  (complex *16)  : charge strength
c     target(3)  (real *8)   : location of the target
c     wavek  (complex *16)   : helmholtz parameter
c
c     OUTPUT:
c
c     pot   (real *8)        : calculated potential
c     fld   (real *8)        : calculated gradient
c
      real *8 source(3),target(3)
      complex *16 wavek,pot,fld(3)
      complex *16 h0,h1,cd,eye,z,ewavek
      complex *16 charge
c
      data eye/(0.0d0,1.0d0)/
c
c ... Caculate offsets and distance
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      zdiff=target(3)-source(3)
      dd=xdiff*xdiff+ydiff*ydiff+zdiff*zdiff
      d=dsqrt(dd)
c
c ... Calculate the potential and field in the regular case:
c
      z=d*wavek
      call h3d01(z,h0,h1)
c
c ... Get potential and field as per required
c
c     Field is - grad(pot).
c
      ewavek=eye*wavek
      pot=h0*ewavek*charge
      if (iffld.eq.1) then
ccc         dinv=-1.0d0/d
         dinv=1.0d0/d
         cd=h1*dinv*ewavek*wavek*charge
         fld(1)=cd*xdiff
         fld(2)=cd*ydiff
         fld(3)=cd*zdiff
      endif
c
      return
      end
c
c
c**********************************************************************
      subroutine h3d01(z,h0,h1)
c**********************************************************************
c
c     Compute spherical Hankel functions of order 0 and 1.
c
c     h0 = dexp(i*z)/(i*z),
c     h1 = - h0' = -h0*(i-1/z) = h0*(1/z-i)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c	z    argument of Hankel functions
c            if abs(z)<1.0d-15, returns zero.
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c	h0 =  h0(z)    (spherical Hankel function of order 0).
c	h1 =  -h0'(z)  (spherical Hankel function of order 1).
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 z,zinv,eye,cd,h0,h1
      data eye/(0.0d0,1.0d0)/, thresh/1.0d-15/, done/1.0d0/
c
      if (abs(z).lt.thresh) then
         h0=0.0d0
         h1=0.0d0
         return
      endif
c
c     Otherwise, use formula
c
      cd = eye*z
      h0=exp(cd)/cd
      h1=h0*(done/z - eye)
c
      return
      end
c


