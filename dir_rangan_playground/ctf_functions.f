ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_ctf_asym_full(Cs,lambda,AmpCnst,df1,df2,
     1     angast,D,ngridr,xnodesr,ngridc,xnodesc,ctfw)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     produces 2D radially asymmetric CTF function

c     input:
c     Cs - spherical aberation of the lens (Angstroms)
c     lambda - electron wavelength (Angstroms) 
c     AmpCnst - amplitude contrast 
c     df1,df2,angast - defocus (Angstroms) and angle of astigmatism
c     D - box size D (Angstroms)
c     ngridr - number of points in the radial direction where to evaluate CTF
c     ngridc - number of discretization points in angle where to evaluate CTF
c     xnodesr - quadrature points in radial direction (in mathematical units,
c     to be rescaled to physical units based on the box size)
c     xnodesc - quadrature points in angle (in mathematical units, to be rescaled
c     to physical units base on the box size)

c     output:
c     ctfw - full 2d radially asymmetric CTF evaluated at quadature points 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      implicit none
c     
      real *8 Cs
      real *8 lambda
      real *8 AmpCnst
      real *8 df1,df2,angast
      real *8 D
      integer ngridr,ngridc
      real *8 xnodesr(ngridr),xnodesc(ngridc)
      complex *16 ctfw(ngridc,ngridr)

      integer ll,rr
      real *8 rad,angspt
      real *8 ctfval
      real *8 pi

      pi = 4.0d0*datan(1.0d0)
      
      do rr = 1,ngridr
      do ll = 1,ngridc
         rad = xnodesr(rr)
         angspt = xnodesc(ll)
      call ctf_asym(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ctfval)
      ctfw(ll,rr) = ctfval
      enddo
      enddo

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_ctf_asym1d(Cs,lambda,AmpCnst,df1,df2,
     1     angast,D,ngridr,xnodesr,angspt,ctfw)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     produces 1D CTF along one radial line (same for every angle)

c     input:
c     Cs - spherical aberation of the lens (Angstroms)
c     lambda - electron wavelength (Angstroms) 
c     AmpCnst - amplitude contrast 
c     df - defocus (Angstroms) and angle of astigmatism
c     D - box size D (Angstroms)
c     ngridr - number of points in the radial direction where to evaluate CTF
c     ngridc - number of discretization points in angle where to evaluate CTF
c     xnodesr - quadrature points in radial direction (in mathematical units,
c     to be rescaled to physical units based on the box size)
c     xnodesc - quadrature points in angle (in mathematical units, to be rescaled
c     to physical units base on the box size)

c     output:
c     ctf1d - 1D CTF along one radial line (same for every angle)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     
      real *8 Cs
      real *8 lambda
      real *8 AmpCnst
      real *8 df1,df2,angast
      real *8 D
      integer ngridr
      real *8 xnodesr(ngridr),angspt
      complex *16 ctfw(ngridr)

      integer ll,rr
      real *8 rad
      real *8 ctfval
      real *8 pi

      pi = 4.0d0*datan(1.0d0)
      
      do rr = 1,ngridr
         rad = xnodesr(rr)
      call ctf_asym(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ctfval)
      ctfw(rr) = ctfval
      enddo

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine our_envelop_fxn(lambda,D,ngridr,xnodesr,envelop)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     produces 2D radially asymmetric CTF function

c     input:
c     lambda - electron wavelength (Angstroms) 
c     D - box size D (Angstroms)
c     ngridr - number of points in the radial direction where to evaluate CTF
c     xnodesr - quadrature points in radial direction (in mathematical units,
c     to be rescaled to physical units based on the box size)

c     output:
c     envelop - 1d envelop function
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      implicit none
c     
      real *8 lambda
      real *8 D
      integer ngridr
      real *8 xnodesr(ngridr)
      real *8 envelop(ngridr)

      integer rr
      real *8 k, theta0,theta(ngridr)
      real *8 pi

      pi = 4.0d0*datan(1.0d0)
      
      k = 2.0d0*pi/lambda
      theta0 = 2.0d-3
      do rr = 1,ngridr
         theta(rr) = xnodesr(rr)/(D*k)
      enddo
c     compute the envelope function
      do rr = 1,ngridr
         envelop(rr) = dexp(-(theta(rr)**2/theta0**2))
      enddo

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_ctf_asym_rad(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ngridc,xnodesc,ctfw)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluates asymmetric CTF at single radial distance for all angles

c     input:
c     Cs - spherical aberation of the lens (Angstroms)
c     lambda - electron wavelength (Angstroms) 
c     AmpCnst - amplitude contrast 
c     df1,df2,angast - defocus (Angstroms) and angle of astigmatism
c     D - box size D (Angstroms)
c     rad - radial distance (in mathematical units)
c     ngridc - number of discretization points in angle where to evaluate CTF
c     xnodesc - quadrature points in angle (in mathematical units, to be rescaled
c     to physical units base on the box size)


c     output:
c     ctfw - asymmetric CTF evaluated at a single radial distance, but many angles 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      implicit none
      integer rr
      real *8 Cs
      real *8 lambda
      real *8 AmpCnst
      real *8 df1,df2,angast
      real *8 D
      real *8 angspt,rad
      integer ngridc
      real *8 xnodesc(ngridc)
      complex *16 ctfw(ngridc)

      integer ll
      real *8 ctfval
      real *8 k, theta0,theta,Bth
      real *8 pi

      pi = 4.0d0*datan(1.0d0)
      
      k = 2.0d0*pi/lambda
      theta0 = 2.0d-3
      theta = rad/(D*k)
c     compute the envelope function
      Bth = dexp(-(theta**2/theta0**2))
      
      do ll = 1,ngridc
         angspt = xnodesc(ll)
      call ctf_asym(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ctfval)
c      ctfw(ll) = Bth*ctfval
      ctfw(ll) = ctfval
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ctf_asym(Cs,lambda,AmpCnst,df1,df2,angast,
     1     D,angspt,rad,ctfval)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluates asymmetric CTF at single point at radial distance rad and angle angspt

c     input:
c     Cs - spherical aberation of the lens (Angstroms)
c     lambda - electron wavelength (Angstroms) 
c     AmpCnst - amplitude contrast 
c     df1,df2,angast - defocus (Angstroms) and angle of astigmatism
c     D - box size D (Angstroms)
c     angspt - angle with x-axis 
c     rad - radial distance (in mathematical units)
c     ngridc - number of discretization points in angle where to evaluate CTF
c     xnodesc - quadrature points in angle (in mathematical units, to be rescaled
c     to physical units base on the box size)

c     output:
c     ctfval - asymmetric CTF evaluated at radial distance rad and angle angspt 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      implicit none
      real *8 Cs
      real *8 lambda
      real *8 AmpCnst,w1,w2
      real *8 df1,df2,angast
      real *8 D
      real *8 angspt,rad

      real *8 ctfval

      real *8 pi,angle,C1,C2,angdif,ccos,df,chi

      pi = 4.0d0*datan(1.0d0)

      w1=dsqrt(1.0d0-AmpCnst**2)
      w2=AmpCnst
c     scale to physical space
      angle = rad/D
      C1 = pi*angle*angle*lambda
      C2 = -C1*Cs*angle*angle*lambda*lambda/2.0d0
      angdif = angspt-angast
      ccos = dcos(2.0d0*angdif)
      df = 0.5d0*(df1+df2+ccos*(df1-df2))
      chi = C1*df+C2
      ctfval=-w1*dsin(chi)-w2*dcos(chi)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine niko_ctf(Cs,lambda,w1,w2,df1,df2,angast,
     1     thetatr,l,m,ctfv)

      implicit none
c     spherical aberation of the lens
      real *8 Cs
c     electron wavelength
      real *8 lambda
      real *8 w1,w2
c     defocus
      real *8 df1,df2,angast
      real *8 thetatr
      real *8 l,m
c     box size D
      real *8 D

      real *8 ctfv

      real *8 pi,rad,angle,angspt,c1,c2,angdif,ccos,df,chi

      pi=4.0d0*datan(1.0d0)
      rad = l**2+m**2
      rad = dsqrt(rad)
      angle = rad*thetatr
      angspt=datan2(m,l)
      c1=2.0d0*pi*angle*angle/(2.0d0*lambda)
      c2=-c1*Cs*angle*angle/2.0d0
      angdif=angspt-angast
      ccos=dcos(2.0d0*angdif)
      df = 0.5d0*(df1+df2+ccos*(df1-df2))
      chi=c1*df+c2

      ctfv=-w1*dsin(chi)-w2*dcos(chi)
      

      return
      end


      subroutine niko_ctf1(Cs,lambda,w1,w2,df1,df2,angast,
     1     D,angspt,rad,ctfv1)

      implicit none
c     spherical aberation of the lens
      real *8 Cs
c     electron wavelength
      real *8 lambda
      real *8 w1,w2
c     defocus
      real *8 df1,df2,angast
c     box size
      real *8 D,angspt
c     distance in mathematical space
      real *8 rad

      real *8 ctfv1

      real *8 pi,angle,c1,c2,angdif,ccos,df,chi 

      pi=4.0d0*datan(1.0d0)

      angle=rad/D
      c1 = pi*angle*angle*lambda
      c2=-c1*Cs*angle*angle*lambda*lambda/2.0d0
      angdif=angspt-angast
      ccos=dcos(2.0d0*angdif)
      df = 0.5d0*(df1+df2+ccos*(df1-df2))
      chi=c1*df+c2
      ctfv1=-w1*dsin(chi)-w2*dcos(chi)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C OUR OLD CTF FUNCTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_ctf_old(ngridr,xnodesr,nrefim,defocus,
     1     ctf1d, addCTFflag)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer addCTFflag
      integer ngridr,nrefim
      real *8 xnodesr(ngridr)
      real *8 defocus(nrefim)

      real *8 lambda,D,k,kD,Cs,w1,w2,theta0
      real *8 theta(ngridr), Bth(ngridr)
      real *8 xi(ngridr,nrefim)
      complex *16 ctf1d(ngridr,nrefim)
      integer rr,cn
      real *8 pi,df
      

      if ( addCTFflag .eq. 0) then
c     set CTF to one

         do cn = 1,nrefim
         do rr = 1,ngridr
               ctf1d(rr,cn)=1.0d0
         enddo
         enddo
         
      else
c     compute CTF
      pi=4.0d0*datan(1.0d0)
c     lambda is the wavelength of electron in A, it is determined by the Voltage
      lambda = 0.0251;
c     the scaling ratio between theta and the mathematical k-space  
c     unit box in mathematical space corresponds to D=100 ( for ex) A in physical space 
      D = 100;
c     wave number
      k = 2.0d0*pi/lambda;
c     scaling factor used in the calculation of theta
      kD = k*D;
c     Cs is the aberation of the lense, Cs = 2mm = 2e7 A
      Cs = 2e7;
c     w
      w2 = 0.07; 
      w1 = dsqrt(1-w2**2);
c     scale theta so it corresponds to the mathematic scale
      theta0 = 2.0d-3
      do rr = 1,ngridr
         theta(rr) = xnodesr(rr)/kD
      enddo
c     compute the envelope function
      do rr = 1,ngridr
         Bth(rr) = dexp(-(theta(rr)**2/theta0**2))
      enddo

      do cn = 1,nrefim
         df = defocus(cn)
      do rr = 1,ngridr
         xi(rr,cn) = -k*df*(theta(rr)**2)/2 + k*Cs*(theta(rr)**4)/4;
      enddo
      enddo
c     the CTF
      do cn = 1,nrefim
         df = defocus(cn)
      do rr = 1,ngridr
         ctf1d(rr,cn)=Bth(rr)*(-w1*dsin(xi(rr,cn))-w2*dcos(xi(rr,cn)))
      enddo
      enddo

      endif

      return
      end




      
