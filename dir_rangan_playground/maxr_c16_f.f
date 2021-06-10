!> Doxygen comment: ;\n
!> maximum of complex *16 array, with maximum taken over real entries. ;\n
      complex *16 function maxr_c16_f(n_x,x_)
      integer n_x,nx
      complex *16 x_(0:n_x-1),m
      if (n_x.le.0) then
         m=dcmplx(0.0d0,0.0d0)
      else
         m=x_(0)
         do nx=0,n_x-1
            if (dreal(x_(nx)).gt.dreal(m)) then
               m = x_(nx)
            end if
         enddo
      end if
      maxr_c16_f = m
      return
      end
