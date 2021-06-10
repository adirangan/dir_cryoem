      subroutine Figheader(unitnumber)
      integer *4 unitnumber,nx
      real *8 rcolor,gcolor,bcolor
      character(len=20) fmtcccccc,fmtr,fmtg,fmtb,fmtc_,fmt_c
      integer *4 tmpcc,tmp_c,tmpc_
      write(unitnumber,'(A)') '#FIG 3.2'
      write(unitnumber,'(A)') 'Landscape'
      write(unitnumber,'(A)') 'Center'
      write(unitnumber,'(A)') 'Inches'
      write(unitnumber,'(A)') 'Letter'
      write(unitnumber,'(A)') '100.00'
      write(unitnumber,'(A)') 'Single'
      write(unitnumber,'(A)') '-2'
      write(unitnumber,'(A)') '1200 2'
      do nx=0,511
         call colorscale(1.0d0*nx,0.0d0,511.0d0,rcolor,gcolor,bcolor)
c$$$         if (nx.eq.0 .or. nx.eq.511) then
c$$$            write(6,'(I0,2X,F6.3,F6.3,F6.3)') nx,rcolor,gcolor,bcolor
c$$$         end if
         tmpcc = int(dmax1(0.0,dmin1(255.0,255.0*rcolor)))
         tmpc_ = tmpcc/16
         tmp_c = mod(tmpcc,16)
         if (tmpc_.lt.10) then
            write(fmtc_,'(I0)') tmpc_
         else if (tmpc_.ge.10) then
            write(fmtc_,'(A)') char(97+tmpc_-10)
         end if
         if (tmp_c.lt.10) then
            write(fmt_c,'(I0)') tmp_c
         else if (tmp_c.ge.10) then
            write(fmt_c,'(A)') char(97+tmp_c-10)
         end if
         fmtr = trim(fmtc_)//trim(fmt_c)
         tmpcc = int(dmax1(0.0,dmin1(255.0,255.0*gcolor)))
         tmpc_ = tmpcc/16
         tmp_c = mod(tmpcc,16)
         if (tmpc_.lt.10) then
            write(fmtc_,'(I0)') tmpc_
         else if (tmpc_.ge.10) then
            write(fmtc_,'(A)') char(97+tmpc_-10)
         end if
         if (tmp_c.lt.10) then
            write(fmt_c,'(I0)') tmp_c
         else if (tmp_c.ge.10) then
            write(fmt_c,'(A)') char(97+tmp_c-10)
         end if
         fmtg = trim(fmtc_)//trim(fmt_c)
         tmpcc = int(dmax1(0.0,dmin1(255.0,255.0*bcolor)))
         tmpc_ = tmpcc/16
         tmp_c = mod(tmpcc,16)
         if (tmpc_.lt.10) then
            write(fmtc_,'(I0)') tmpc_
         else if (tmpc_.ge.10) then
            write(fmtc_,'(A)') char(97+tmpc_-10)
         end if
         if (tmp_c.lt.10) then
            write(fmt_c,'(I0)') tmp_c
         else if (tmp_c.ge.10) then
            write(fmt_c,'(A)') char(97+tmp_c-10)
         end if
         fmtb = trim(fmtc_)//trim(fmt_c)
         fmtcccccc = trim(fmtr)//trim(fmtg)//trim(fmtb)
         write(unitnumber,'(A,I0,1X,A,A)') '0 ',(nx+32),'#',fmtcccccc
      enddo
      end
