ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine kspace_model_norm(model_sph,lmod,
     1     nterms_sph,wtsr,ncur,rmodelnorm)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute norm of model in k-space in L2 using spherical 
c     harmonic basis on successive spheres.
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
c     ncur         number of points in radial direction
c
c     OUTPUT:
c
      implicit none
      integer lmod,ncur
      integer nterms_sph(ncur)
      integer i,j,nsphstart,nsph
      real *8 wtsr(ncur)
      real *8 rmodelnorm,rmod,pi
      complex *16 model_sph(0:lmod-1)
c      
      pi = 4.0d0*datan(1.0d0)
      nsphstart = 0
      rmodelnorm = 0
      do i = 1,ncur
         nsph = (nterms_sph(i)+1)**2
         rmod = 0.0d0
         do j = nsphstart,nsphstart+nsph-1
            rmod = rmod + 4*pi*abs(model_sph(j))**2
         enddo
         rmod = rmod*wtsr(i)
         rmodelnorm = rmodelnorm + rmod
         nsphstart = nsphstart + nsph
      enddo
      rmodelnorm = dsqrt(rmodelnorm)/(dsqrt(2*pi)**3)
c
      return
      end

