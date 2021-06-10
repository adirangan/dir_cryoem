      subroutine readmrc_size(filename,n_x,n_y,n_M)
c$$$      Note: this function assumes that the "mode" of the mrc_file
c$$$            (read in as the fourth integer in the header)
c$$$            is equal to 2. (i.e., mode.eq.2 --> real *4 data)
c$$$      INPUTS:
c$$$      character filename(*) : string containing name of mrc_file
c$$$      OUTPUTS:
c$$$      integer *4 n_x : number of pixels in x-direction of each image
c$$$      integer *4 n_y : number of pixels in y-direction of each image
c$$$      integer *4 n_M : number of images
      implicit none
      integer verbose
      data verbose / 1 /
      character(len=1024) filename
      integer n_x,n_y,n_M,n_M_all
      integer unitnumber
      integer np
      integer n_h0,nh0
      parameter (n_h0=10)
      integer *4 h0_(0:n_h0-1)
      integer n_h1,nh1
      parameter (n_h1=12)
      real *4 h1_(0:n_h1-1)
      integer n_h2,nh2
      parameter (n_h2=30)
      integer *4 h2_(0:n_h2-1)
      integer n_h3,nh3
      parameter (n_h3=8)
      character h3_(0:n_h3-1)
      integer n_h4,nh4
      parameter (n_h4=2)
      integer *4 h4_(0:n_h4-1)
      integer n_tmp
      integer n_h5,nh5
      parameter (n_h5=80)
      character h5_(0:n_h5-1)
      character(len=64) format_string
      if (verbose.gt.0) then
         write(6,'(A,A)') 'opening: ',trim(filename)
      end if
      unitnumber = 7
      open(unitnumber,file=filename,access='stream',form
     $     ='unformatted',status='old')
      np=1
      do nh0=0,n_h0-1
         read(unitnumber,pos=np) h0_(nh0)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h0,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h0_: ',(h0_(nh0),nh0=0,n_h0-1)
      end if
      n_x = h0_(0)
      n_y = h0_(1)
      n_M_all = h0_(2)
      n_M = n_M_all
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0,A,I0,A,I0)') 'n_x: ',n_x,'; n_y: ',n_y
     $        ,'; n_M: ',n_M,'; Mode: ',h0_(3)
      end if
      if (h0_(3).ne.2) then
         write(6,'(A,I0,A)') 'Warning. Mode ',h0_(3),' in readmrc'
      end if
      do nh1=0,n_h1-1
         read(unitnumber,pos=np) h1_(nh1)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h1,'(F8.2,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h1_: ',(h1_(nh1),nh1=0,n_h1-1)
      end if
      do nh2=0,n_h2-1
         read(unitnumber,pos=np) h2_(nh2)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h2,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h2_: ',(h2_(nh2),nh2=0,n_h2-1)
      end if
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'Extended header length: ',h2_(1)
      end if;
      do nh3=0,n_h3-1
         read(unitnumber,pos=np) h3_(nh3)
         np = np+1;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h3,'(A,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h3_: ',(h3_(nh3),nh3=0,n_h3-1)
      end if
      do nh4=0,n_h4-1
         read(unitnumber,pos=np) h4_(nh4)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h4,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h4_: ',(h4_(nh4),nh4=0,n_h4-1)
      end if
      do n_tmp=0,max(9,h4_(1)-1)
         do nh5=0,n_h5-1
            read(unitnumber,pos=np) h5_(nh5)
            np = np+1;
         enddo
         write(format_string,'(A,I0,A)') '(A,',n_h5,'(A,2X))'
         if (verbose.gt.1) then
            write(6,format_string) 'h5_: ',(h5_(nh5),nh5=0,n_h5-1)
         end if
      enddo
      np=np+h2_(1)
      close(unitnumber,status='keep')
      end;

      subroutine readmrc_data(filename,n_M,M_)
c$$$      Note: this function assumes that the "mode" of the mrc_file
c$$$            (read in as the fourth integer in the header)
c$$$            is equal to 2. (i.e., mode.eq.2 --> real *4 data)
c$$$      INPUTS:
c$$$      character filename(*) : string containing name of mrc_file
c$$$      integer n_M : number of images to read
c$$$      OUTPUTS:
c$$$      complex *16 M_ : array storing the images ; 
c$$$                       The nx,ny entry of image nm is stored in:
c$$$                       M_(nx + ny*n_x + nm*n_x*n_y)
c$$$                       where n_x and n_y correspond to image-size.
      implicit none
      integer verbose
      data verbose / 1 /
      character(len=1024) filename
      integer n_M
      complex *16 M_(0:0)
      integer n_x,n_y,n_M_all
      integer unitnumber
      integer np
      integer n_h0,nh0
      parameter (n_h0=10)
      integer *4 h0_(0:n_h0-1)
      integer n_h1,nh1
      parameter (n_h1=12)
      real *4 h1_(0:n_h1-1)
      integer n_h2,nh2
      parameter (n_h2=30)
      integer *4 h2_(0:n_h2-1)
      integer n_h3,nh3
      parameter (n_h3=8)
      character h3_(0:n_h3-1)
      integer n_h4,nh4
      parameter (n_h4=2)
      integer *4 h4_(0:n_h4-1)
      integer n_tmp
      integer n_h5,nh5
      parameter (n_h5=80)
      character h5_(0:n_h5-1)
      real *4 f_tmp
      integer nc
      character(len=64) format_string
      if (verbose.gt.0) then
         write(6,'(A,A)') 'opening: ',trim(filename)
      end if
      unitnumber = 7
      open(unitnumber,file=filename,access='stream',form
     $     ='unformatted',status='old')
      np=1
      do nh0=0,n_h0-1
         read(unitnumber,pos=np) h0_(nh0)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h0,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h0_: ',(h0_(nh0),nh0=0,n_h0-1)
      end if
      n_x = h0_(0)
      n_y = h0_(1)
      n_M_all = h0_(2)
      if (verbose.gt.1) then
         write(6,'(A,I0,A,I0,A,I0,A,I0)') 'n_x: ',n_x,'; n_y: ',n_y
     $        ,'; n_M: ',n_M_all,'; Mode: ',h0_(3)
      end if
      if (h0_(3).ne.2) then
         write(6,'(A,I0,A)') 'Warning. Mode ',h0_(3),' in readmrc'
      end if
      do nh1=0,n_h1-1
         read(unitnumber,pos=np) h1_(nh1)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h1,'(F8.2,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h1_: ',(h1_(nh1),nh1=0,n_h1-1)
      end if
      do nh2=0,n_h2-1
         read(unitnumber,pos=np) h2_(nh2)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h2,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h2_: ',(h2_(nh2),nh2=0,n_h2-1)
      end if
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'Extended header length: ',h2_(1)
      end if;
      do nh3=0,n_h3-1
         read(unitnumber,pos=np) h3_(nh3)
         np = np+1;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h3,'(A,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h3_: ',(h3_(nh3),nh3=0,n_h3-1)
      end if
      do nh4=0,n_h4-1
         read(unitnumber,pos=np) h4_(nh4)
         np = np+4;
      enddo
      write(format_string,'(A,I0,A)') '(A,',n_h4,'(I0,2X))'
      if (verbose.gt.1) then
         write(6,format_string) 'h4_: ',(h4_(nh4),nh4=0,n_h4-1)
      end if
      do n_tmp=0,max(9,h4_(1)-1)
         do nh5=0,n_h5-1
            read(unitnumber,pos=np) h5_(nh5)
            np = np+1;
         enddo
         write(format_string,'(A,I0,A)') '(A,',n_h5,'(A,2X))'
         if (verbose.gt.1) then
            write(6,format_string) 'h5_: ',(h5_(nh5),nh5=0,n_h5-1)
         end if
      enddo
      np=np+h2_(1)
      nc = 0
      do n_tmp=0,n_x*n_y*n_M
         read(unitnumber,pos=np) f_tmp
         np = np+4
         M_(nc) = cmplx(f_tmp,0.0d0)
         nc = nc+1
      enddo
      close(unitnumber,status='keep')
      end;
