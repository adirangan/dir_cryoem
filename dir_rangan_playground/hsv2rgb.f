!> Doxygen comment: ;\n
!> map from hsv coordinates to rgb coordinates. ;\n
      subroutine hsv2rgb(h,s,v,r,g,b)
      real *8 h,s,v,r,g,b,f,p,q,t
      integer *4 i
      if (s.eq.0.0) then
         r=v
         g=v
         b=v
      else
         i=int(h/60.0)
         f=(h/60.0)-i
         p=v*(1.0-s)
         q=v*(1.0-s*f)
         t=v*(1-s*(1-f))
         if (i.eq.0) then
            r=v
            g=t
            b=p
         else if (i.eq.1) then
            r=q
            g=v
            b=p
         else if (i.eq.2) then
            r=p
            g=v
            b=t
         else if (i.eq.3) then
            r=p
            g=q
            b=v
         else if (i.eq.4) then
            r=t
            g=p
            b=v;
         else if (i.eq.5) then
            r=v
            g=p
            b=q
         end if
      end if
      end
