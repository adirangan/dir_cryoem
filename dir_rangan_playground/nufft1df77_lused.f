      subroutine nufft1df77_1d1_lused(eps,ms,lused)
c$$$      returns lused as computed by nufft1d1 within nufft1df77.f ;
      implicit none
      real *8 eps !real *8: accuracy. ;
      integer ms !integer: number of fourier modes given. ;
      real *8 rat !temporary: real *8: oversampling ratio. ;
      integer nspread !temporary: integer: number of points in spreader. ;
      integer nf1 !temporary: integer: number of points in oversampled regular mesh. ;
      integer iw1,iwsav !temporary: integer. ;
      integer lused !integer: workspace requested (number of elements). ;
      integer next235 !integer: function output. ;
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif !if (eps.le.1d-11) then
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif !if (2*nspread.gt.nf1) then
      iw1 = 2*nf1
      iwsav = iw1 + nspread + 1
      lused = iwsav + 4*nf1 + 15
      end !subroutine

      subroutine nufft1df77_1d2_lused(eps,ms,lused)
c$$$      returns lused as computed by nufft1d2 within nufft1df77.f ;
      implicit none
      real *8 eps !real *8: accuracy. ;
      integer ms !integer: number of fourier modes given. ;
      real *8 rat !temporary: real *8: oversampling ratio. ;
      integer nspread !temporary: integer: number of points in spreader. ;
      integer nf1 !temporary: integer: number of points in oversampled regular mesh. ;
      integer iw1,iwsav !temporary: integer. ;
      integer lused !integer: workspace requested (number of elements). ;
      integer next235 !integer: function output. ;
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif !if (eps.le.1d-11) then
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif !if (2*nspread.gt.nf1) then
      iw1 = 2*nf1
      iwsav = iw1 + nspread + 1
      lused = iwsav + 4*nf1 + 15
      end !subroutine

      subroutine nufft1df77_1d3_lused(nj,xj,nk,sk,eps,ms,lused)
c$$$      returns lused as computed by nufft1d3 within nufft1df77.f ;
      implicit none
      integer nj,nk !integer: number of sources and targets. ;
      real *8 xj(nj) !real *8 array (length nj): location of sources. ;
      real *8 sk(nk) !real *8 array (length nk): location of targets. ;
      real *8 eps !real *8: accuracy. ;
      integer ms !integer: number of fourier modes given. ;
      real *8 rat !temporary: real *8: oversampling ratio. ;
      integer j,k1 !temporary: integer. ;
      real *8 t1,t2,xb,xm,sb,sm !temporary: real *8. ;
      integer nspread !temporary: integer: number of points in spreader. ;
      integer nf1 !temporary: integer: number of points in oversampled regular mesh. ;
      integer iw1,iwsave !temporary: integer. ;
      integer lused !integer: workspace requested (number of elements). ;
      integer next235 !integer: function output. ;
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      t1 = xj(1)
      t2 = xj(1)
      do j = 2, nj
         if (xj(j).gt.t2) then
             t2=xj(j)
         else if (xj(j).lt.t1) then
             t1=xj(j)
         endif
      enddo
      xb = (t1+t2) / 2d0
      xm = max(t2-xb,-t1+xb)  ! max(abs(t2-xb),abs(t1-xb))
c
      t1 = sk(1)
      t2 = sk(1)
      do k1 = 2, nk
         if (sk(k1).gt.t2) then
             t2=sk(k1)
         else if (sk(k1).lt.t1) then
             t1=sk(k1)
         endif
      enddo
      sb = (t1+t2) / 2d0
      sm = max(t2-sb,-t1+sb)
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif !if (eps.le.1d-11) then
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      t1 = 2d0/pi * xm*sm
      nf1 = next235(rat*max(rat*t1+2*nspread,2*nspread/(rat-1)))
      iw1 = 2*nf1
      iwsave = iw1 + nspread + 1
      lused = iwsave + 16 + 4*nf1
      end !subroutine

      
      
      
