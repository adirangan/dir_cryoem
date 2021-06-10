c     gfortran -o WritePNMfile_dr.out WritePNMfile_dr.f ; ./WritePNMfile_dr.out ; 
      program WritePNMfile_dr
      real *8 pi,sigma_x,sigma_y,lmax,lmin
      integer *4 i,imax,j,jmax
      real *8, allocatable :: x_(:,:)
      character *1, allocatable :: y_(:,:)
      character(len=20) fmtf,fmtc
      pi = 4*atan(1.0)
      imax = 128
      jmax = 256
      write(fmtf,'(A,I0,A)') '(',jmax,'(F6.5,1X))'
      write(fmtc,'(A,I0,A)') '(',jmax,'(A,1X))'
      sigma_x = imax/8;
      sigma_y = jmax/6;
      allocate(x_(0:imax-1,0:jmax-1))
      do i=0,imax-1
         do j=0,jmax-1
            x_(i,j) = 1.0/(2*pi*sigma_x*sigma_y)*exp(-((i-imax/2.0)**2)
     $           /2.0/sigma_x**2)*exp(-((j-jmax/2.0)**2)/2.0/sigma_y**2)
     $           + 1.0/(2*pi*sigma_x*sigma_y)*exp(-((i-3.0*imax/4.0)**2)
     $           /2.0/sigma_x**2)*exp(-((j-9.0*jmax/10.0)**2)/2.0
     $           /sigma_y**2)
         enddo
      enddo
      if (imax .lt. 16 .and. jmax .lt. 16) then
         write(6,fmtf) ((x_(i,j),j=0,jmax-1),i=0,imax-1)
      end if
      lmax = 2.0/(2*pi*sigma_x*sigma_y)
      lmin = 0.0000
      allocate(y_(0:imax-1,0:jmax-1))
      do i=0,imax-1
         do j=0,jmax-1
            y_(i,j) = char(int(dmax1(0.0,dmin1(255.0,255.0*(x_(i,j)
     $           -lmin)/(lmax-lmin))),4))
         enddo
      enddo
      if (imax .lt. 16 .and. jmax .lt. 16) then
         write(6,'(A)') ''
         write(6,fmtc) ((y_(i,j),j=0,jmax-1),i=0,imax-1)
      end if
      open(7,file='WritePNMfile.pnm',status='replace',form='formatted')
      write(7,'(A)') 'P5'
      write(7,'(A,1X,F6.5,1X,F6.5)') '#',lmin,lmax
      write(7,'(I3,1X,I3)') jmax,imax
      write(7,'(I3)') 255
      close(7,status='keep')
      open(7,file='WritePNMfile.pnm',status='old',access='append',form
     $     ='unformatted')
      write(7) ((y_(i,j),j=0,jmax-1),i=0,imax-1)
      close(7,status='keep')
      stop
      end
