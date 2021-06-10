c     gfortran -o WriteFigfile_dr.out WriteFigfile_dr.f ; ./WriteFigfile_dr.out ; 
      program WriteFigfile_dr
      real *8 pi,lmax,lmin,xval_tmp,yval_tmp
      real *8 rval_pre,tval_pre,rval_pos,tval_pos,tval_mid
      integer *4 diameter,arclength
      integer *4 nr,nrmax,nt,nxmax,nx,narc
      integer *4, allocatable :: ntmax_(:)
      real *8, allocatable :: rvals_(:)
      real *8, allocatable :: tvals_(:)
      real *8, allocatable :: x_(:)
      integer *4, allocatable :: y_(:)
      integer *4, allocatable :: xyval_(:)
      character(len=20) fmtc
      character(len=1024) fig2dev_system_call
      integer system,system_error
      pi = 4*atan(1.0)
      nrmax = 32
      allocate(rvals_(0:nrmax-1))
      allocate(ntmax_(0:nrmax-1))
      nxmax=0
      do nr=0,nrmax-1
         rvals_(nr) = 1.0*nr
         ntmax_(nr) = 2*nr+3
         nxmax = nxmax + ntmax_(nr)
      enddo
      allocate(tvals_(0:nxmax-1))
      allocate(x_(0:nxmax-1));
      nx=0
      do nr=0,nrmax-1
         do nt=0,ntmax_(nr)-1
            tvals_(nx) = 2.0*pi*(1.0*nt)/(1.0*ntmax_(nr))
            nx = nx + 1;
         enddo
      enddo
      nx=0
      do nr=0,nrmax-1
         do nt=0,ntmax_(nr)-1
            xval_tmp = rvals_(nr)*cos(tvals_(nx))
            yval_tmp = rvals_(nr)*sin(tvals_(nx))
c$$$            x_(nx) = exp(-rvals_(nr)/(nrmax/3.0))*(cos(2*pi*xval_tmp
c$$$     $           /5.0))
c$$$            x_(nx) = 0.5*xval_tmp/rvals_(nrmax-1)
            x_(nx) = 0.5*yval_tmp/rvals_(nrmax-1)
            nx = nx + 1;
         enddo
      enddo
      lmax = +0.5
      lmin = -0.5
      allocate(y_(0:nxmax-1));
      nx=0
      do nr=0,nrmax-1
         do nt=0,ntmax_(nr)-1
            y_(nx) = 32 + int(dmax1(0.0,dmin1(511.0,511.0*(x_(nx)-lmin)
     $           /(lmax-lmin))))
            nx = nx + 1;
         enddo
      enddo
      open(7,file='WriteFigfile.fig',status='replace',form='formatted')
      call Figheader(7);
      diameter = 1000
      arclength = 8
      allocate(xyval_(0:4*arclength+1))
      write(fmtc,'(A,I0,A)') '(A,1X,',(4*arclength+2),'(I0,1X))'
      nx=0
      do nr=0,nrmax-1
         if (nr.gt.0) then
            rval_pre = 0.5*(rvals_(nr) + rvals_(nr-1))
         else
            rval_pre = rvals_(0)
         end if
         if (nr.lt.nrmax-1) then
            rval_pos = 0.5*(rvals_(nr+1) + rvals_(nr))
         else 
            rval_pos = rvals_(nrmax-1)
         end if
         do nt=0,ntmax_(nr)-1
            if (nt.gt.0) then
               tval_pre = 0.5*(tvals_(nx) + tvals_(nx-1))
            else
               tval_pre = 0.5*(tvals_(nx) + tvals_(nx+ntmax_(nr)-1) - 2
     $              *pi)
            end if
            if (nt.lt.ntmax_(nr)-1) then
               tval_pos = 0.5*(tvals_(nx+1) + tvals_(nx))
            else
               tval_pos = 0.5*(tvals_(nx-ntmax_(nr)+1) + 2*pi +
     $              tvals_(nx))
            end if
            do narc=0,arclength-1
               tval_mid = (narc*tval_pos + (arclength-1-narc)*tval_pre)
     $              /(arclength-1)
               xyval_(0 + 2*narc) = +int(diameter*rval_pre
     $              *cos(tval_mid))
               xyval_(1 + 2*narc) = -int(diameter*rval_pre
     $              *sin(tval_mid))
            enddo
            do narc=0,arclength-1
               tval_mid = (narc*tval_pre + (arclength-1-narc)*tval_pos)
     $              /(arclength-1)
               xyval_(0 + 2*narc + 2*arclength) = +int(diameter*rval_pos
     $              *cos(tval_mid))
               xyval_(1 + 2*narc + 2*arclength) = -int(diameter*rval_pos
     $              *sin(tval_mid))
            enddo
            xyval_(4*arclength+0) = +int(diameter*rval_pre
     $           *cos(tval_pre))
            xyval_(4*arclength+1) = -int(diameter*rval_pre
     $           *sin(tval_pre))
            write(7,'(A,I0,A,I0)') '2 3 0 0 0 ',y_(nx)
     $           ,' 50 -1 20 0.000 0 0 -1 0 0 ',(2*arclength+1)
            write(7,fmtc) char(9),(xyval_(j),j=0,4*arclength+1)
            nx = nx + 1;
         enddo
      enddo
      close(7,status='keep')
      write(fig2dev_system_call,'(A)')
     $     'fig2dev -Ljpeg -q 10 WriteFigfile.fig WriteFigfile.jpg;'
      system_error = system(fig2dev_system_call)
      stop
      end

      include 'hsv2rgb.f'
      include 'colorscale.f'
      include 'Figheader.f'
