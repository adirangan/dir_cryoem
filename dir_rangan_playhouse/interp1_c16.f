      subroutine interp1_c16(n_x,min_x,max_x,F_,x_loc,output)
c$$$      assumes regular grid for x
c$$$      periodic boundary conditions
      integer *4 n_x
      real *8 min_x,max_x
      complex *16 F_(0:n_x-1)
      real *8 x_loc
      complex *16 output
      real *8 nx_d
      integer *4 nx_pre,nx_pos
      real *8 dx_pre,dx_pos
      real *8 d_pre,d_pos
      real *8 A,B,alpha,beta
      call periodize_r8(x_loc,min_x,max_x,x_loc)
      nx_d = n_x*(x_loc-min_x)/(max_x-min_x)
      nx_pre = floor(nx_d)
      dx_pre = dabs(nx_pre - nx_d)
      nx_pos = ceiling(nx_d)
      dx_pos = dabs(nx_pos - nx_d)
      if (nx_pos.ge.n_x) then
         nx_pos = nx_pos - n_x
      end if
      d_pre = (dx_pre)
      d_pos = (dx_pos)
      A = d_pre
      B = d_pos
      if (A+B.le.0) then
         output = F_(nx_pre)
      else
         alpha = A/(A+B)
         beta = B/(A+B)
         output = beta*F_(nx_pre) + alpha*F_(nx_pos)
      end if
      end
