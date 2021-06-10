      real *8 function al2_c16_f(n_x,x_)
      integer n_x,nx
      complex *16 x_(0:n_x-1)
      real *8 m
      if (n_x.le.0) then
         m=0.0d0
      else
         m=dreal(dconjg(x_(0))*x_(0))
         do nx=1,n_x-1
            m = m + dreal(dconjg(x_(nx))*x_(nx))
         enddo
         m = dsqrt(m)
      end if
      al2_c16_f = m
      return
      end
