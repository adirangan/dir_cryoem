!> Doxygen comment: ;\n
!> Calculates angle between complex *16 vectors x_ and y_  ;\n
      real *8 function costheta_c16_f(n_a,x_,y_)
      integer n_a,na
      complex *16 x_(0:n_a-1)
      complex *16 y_(0:n_a-1)
      real *8 xy
      real *8 xx
      real *8 yy
      if (n_a.le.0) then
         xy=0.0d0
         xx=1.0d0
         yy=1.0d0
      else
         xy=dreal(dconjg(x_(0))*y_(0))
         xx=dreal(dconjg(x_(0))*x_(0))
         yy=dreal(dconjg(y_(0))*y_(0))
         do na=1,n_a-1
            xy = xy + dreal(dconjg(x_(na))*y_(na))
            xx = xx + dreal(dconjg(x_(na))*x_(na))
            yy = yy + dreal(dconjg(y_(na))*y_(na))
         enddo
         if (xx.le.0.0d0) xx=1.0d0
         if (yy.le.0.0d0) yy=1.0d0
      end if
      costheta_c16_f = xy/(dsqrt(xx)*dsqrt(yy))
      return
      end
