!> Doxygen comment: ;\n
!> Calculates frobenius norm of the difference between complex *16 vectors x_ and y_  ;\n
      real *8 function fnod_c16_f(n_a,x_,y_)
      integer n_a,na
      complex *16 x_(0:n_a-1)
      complex *16 y_(0:n_a-1)
      complex *16 tmp_d
      real *8 tmp_f
      if (n_a.le.0) then
         tmp_f=0.0d0
      else
         tmp_d = x_(0)-y_(0)
         tmp_f=dreal(dconjg(tmp_d)*tmp_d)
         do na=1,n_a-1
            tmp_d = x_(na)-y_(na)
            tmp_f = tmp_f + dreal(dconjg(tmp_d)*tmp_d)
         enddo
      end if
      fnod_c16_f = dsqrt(tmp_f)
      return
      end
