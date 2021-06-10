!> Doxygen comment: ;\n
!> colormap_beach to map interval into RGB vector. ;\n
      subroutine colormap_beach(val,valmin,valmax,rcolor,gcolor,bcolor)
      implicit none
      real *8 val,valmin,valmax,rcolor,gcolor,bcolor,valmin2,valmax2,v,h
      real *8 gamma1,gamma2,v1,v2,v3
      gamma1 = 0.50d0
      gamma2 = 0.25d0
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
      v = (v-valmin2)/(valmax2-valmin2)
      if (v.lt.0.5) then
         v1 = v*2.0d0
         v1 = v1**gamma1
         rcolor = v1**6
         gcolor = 0.5d0 + 0.5d0*v1
         bcolor = 1
      end if !if (v.lt.0.5) then
      if (v.eq.0.5) then
         rcolor = 1.0d0
         gcolor = 1.0d0
         bcolor = 1.0d0
      end if !if (v.eq.0.5) then
      if (v.gt.0.5) then
         v2 = 1.0 - 2*(v-0.5d0)
         v2 = v2**gamma2
         rcolor = 1.0d0
         gcolor = v2**6
         bcolor = 0.5d0 + 0.5d0*v2
      end if !if (v.gt.0.5) then
      v3 = 2*(v-0.5)
      bcolor = bcolor*dsqrt(abs(v3))
      end
