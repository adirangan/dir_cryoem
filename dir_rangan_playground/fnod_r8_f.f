!> Doxygen comment: ;\n
!> Calculates frobenius norm of the difference between real *8 vectors x_ and y_  ;\n
      real *8 function fnod_r8_f(n_a,x_,y_)
      integer n_a,na
      real *8 x_(0:n_a-1)
      real *8 y_(0:n_a-1)
      real *8 tmp_d
      real *8 tmp_f
      if (n_a.le.0) then
         tmp_f=0.0d0
      else
         tmp_d = x_(0)-y_(0)
         tmp_f=tmp_d*tmp_d
         do na=1,n_a-1
            tmp_d = x_(na)-y_(na)
            tmp_f = tmp_f + tmp_d*tmp_d
         enddo
      end if
      fnod_r8_f = dsqrt(tmp_f)
      return
      end
