ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_trans(ngridcv,ntemplatesize,icstart,
     1           ncur,template,template_fourier)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Fourier transform templates on rings in polar grid.
c
c     INPUT: 
c
c     ngridcv         number of points on successive rings
c     ntemplatesize   total number of points on quasiuniform polar grid
c     icstart         indexing array for points on successive rings
c     ncur            number of rings
c     template        template
c
c     OUTPUT:
c
c     template_fourier  FFT of template on each ring in standard FFT
c                        format.
c
      implicit none
      integer ik,j
      integer ngridcv(ncur),ncur,ntemplatesize,icstart(ncur),numonsphere
      real *8 rscale
      complex *16 template(ntemplatesize)
      complex *16 template_fourier(ntemplatesize)
      complex *16, allocatable :: wsave(:)
      complex *16, allocatable :: wtemp(:)

c
      do ik = 1,ncur
         allocate(wsave(4*ngridcv(ik)+16))
         allocate(wtemp(ngridcv(ik)))
         call dcffti(ngridcv(ik),wsave)
         do j = 1,ngridcv(ik)
            wtemp(j) = template(icstart(ik)+j-1)
         enddo
         call dcfftf(ngridcv(ik),wtemp,wsave)
ccc         rscale = 1.0d0/dsqrt(ngridcv(ik)+0.0d0)
         rscale = 1.0d0/(ngridcv(ik)+0.0d0)
         do j = 1,ngridcv(ik)
            template_fourier(icstart(ik)+j-1) =
     1              wtemp(j)*rscale
         enddo
         deallocate(wsave)
         deallocate(wtemp)
      enddo
      return
      end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine template_trans_all(ngridcv,ntemplatesize,icstart,
     1           ncur,numonsphere,templates,templates_fourier)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fourier transform templates on rings in polar grid.
c
c     INPUT: 
c
c     ngridcv         number of points on successive rings
c     ntemplatesize   total number of points on quasiuniform polar grid
c     icstart         indexing array for points on successive rings
c     ncur            number of rings
c     numonsphere     number of templates
c     templates       templates
c
c     OUTPUT:
c
c     templates_fourier  FFT of templates on each ring in standard FFT
c                        format.
c
      implicit none
      integer ncur,ik,j,itemp,numonsphere
      integer ngridcv(ncur),ntemplatesize,icstart(ncur)
      real *8 rscale
      complex *16 templates(ntemplatesize,numonsphere)
      complex *16 templates_fourier(ntemplatesize,numonsphere)
      complex *16, allocatable :: wsave(:)
      complex *16, allocatable :: wtemp(:)

c
      do ik = 1,ncur
         allocate(wsave(4*ngridcv(ik)+16))
         allocate(wtemp(ngridcv(ik)))
         call dcffti(ngridcv(ik),wsave)
         do itemp = 1,numonsphere
            do j = 1,ngridcv(ik)
               wtemp(j) = templates(icstart(ik)+j-1,itemp)
            enddo
            call dcfftf(ngridcv(ik),wtemp,wsave)
            rscale = 1.0d0/(ngridcv(ik)+0.0d0)
            do j = 1,ngridcv(ik)
               templates_fourier(icstart(ik)+j-1,itemp) =
     1                 wtemp(j)*rscale
            enddo
         enddo
         deallocate(wsave)
         deallocate(wtemp)
      enddo
      return
      end
