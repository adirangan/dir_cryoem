!> Doxygen comment: ;\n
!> Simple colorscale to map interval into RGB vector. ;\n
      subroutine colorscale(val,valmin,valmax,rcolor,gcolor,bcolor)
      real *8 val,valmin,valmax,rcolor,gcolor,bcolor,valmin2,valmax2,v,h
      valmin2=valmin
      valmax2=valmax
      if (valmax2.le.valmin2) then
         valmax2 = valmin2+1.0;
      end if
      v=val
      if (v.le.valmin2) then
         v=valmin2
      end if
      if (v.ge.valmax2) then
         v=valmax2
      end if
c$$$      h = 300.0*(valmax2-v)/(valmax2-valmin2) - 180.0
      h = 300.0*(valmax2-v)/(valmax2-valmin2) - 60.0
      if (h.gt.360.0) then
         h=h-360.0
      end if
      if (h.lt.0.0) then
         h=h+360.0
      end if
      call hsv2rgb(h,1.0d0,1.0d0,rcolor,gcolor,bcolor)
      end
