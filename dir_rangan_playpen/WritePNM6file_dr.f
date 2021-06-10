c     gfortran -o WritePNM6file_dr.out WritePNM6file_dr.f ; ./WritePNM6file_dr.out ; 
      program WritePNM6file_dr
      real *8 pi,sigma_x,sigma_y,lmax,lmin,rcolor,gcolor,bcolor
      integer *4 i,imax,j,jmax
      real *8, allocatable :: x_(:,:)
      character *1, allocatable :: y_(:)
      character(len=20) fmtf,fmtc
      pi = 4*atan(1.0)
      imax = 128
      jmax = 256
      write(fmtf,'(A,I0,A)') '(',jmax,'(F6.5,1X))'
      write(fmtc,'(A,A,A)') '(','3','(I0,1X))'
      sigma_x = imax/8;
      sigma_y = jmax/6;
c$$$      allocate(x_(0:imax-1,0:jmax-1))
      allocate(x_(imax,jmax))
      do i=0,imax-1
         do j=0,jmax-1
c$$$            x_(i,j) = 1.0/(2*pi*sigma_x*sigma_y)*exp(-((i-imax/2.0)**2)
c$$$     $           /2.0/sigma_x**2)*exp(-((j-jmax/2.0)**2)/2.0/sigma_y**2)
c$$$     $           + 1.0/(2*pi*sigma_x*sigma_y)*exp(-((i-3.0*imax/4.0)**2)
c$$$     $           /2.0/sigma_x**2)*exp(-((j-9.0*jmax/10.0)**2)/2.0
c$$$     $           /sigma_y**2)
            x_(1+i,1+j) = j
         enddo
      enddo
      if (imax .lt. 16 .and. jmax .lt. 16) then
         write(6,fmtf) ((x_(i,j),j=1,jmax),i=1,imax)
      end if
c$$$      lmax = 1.0/(2*pi*sigma_x*sigma_y)
c$$$      lmin = 0.1/(2*pi*sigma_x*sigma_y)
      lmax = 1.0*jmax;
      lmin = 0.0;
c$$$      allocate(y_(0:(3*imax*jmax-1)))
      allocate(y_(3*imax*jmax))
      do i=0,imax-1
         do j=0,jmax-1
            call colorscale(x_(1+i,1+j),lmin,lmax,rcolor,gcolor,bcolor)
            y_(1 + 2 + j*3 + i*3*jmax) = char(int(dmax1(0.0,dmin1(255.0
     $           ,255.0*rcolor))))
            y_(1 + 0 + j*3 + i*3*jmax) = char(int(dmax1(0.0,dmin1(255.0
     $           ,255.0*gcolor))))
            y_(1 + 1 + j*3 + i*3*jmax) = char(int(dmax1(0.0,dmin1(255.0
     $           ,255.0*bcolor))))
         enddo
      enddo
      if (imax .lt. 16 .and. jmax .lt. 16) then
         write(6,fmtc) (y_(i),i=1,(3*imax*jmax))
      end if
      open(7,file='WritePNM6file.pnm',status='replace',form='formatted')
      write(7,'(A)') 'P6'
      write(7,'(A,1X,F6.5,1X,F6.5)') '#',lmin,lmax
      write(7,'(I3,1X,I3)') jmax,imax
      write(7,'(I3)') 255
      close(7,status='keep')
      open(7,file='WritePNM6file.pnm',status='old',access='append',form
     $     ='unformatted')
      write(7) (y_(i),i=1,(3*imax*jmax))
      close(7,status='keep')
      stop
      end

      include 'hsv2rgb.f'
      include 'colorscale.f'

