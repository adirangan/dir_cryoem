      subroutine interp2_c16(n_x,min_x,max_x,n_y,min_y,max_y,F_,x_loc
     $     ,y_loc,output)
c$$$      assumes regular grid for x and y
c$$$      periodic boundary conditions
      integer *4 n_x,n_y
      real *8 min_x,max_x,min_y,max_y
      complex *16 F_(0:n_x*n_y-1)
      real *8 x_loc,y_loc
      complex *16 output
      real *8 nx_d,ny_d
      integer *4 nx_pre,ny_pre,nx_pos,ny_pos
      real *8 dx_pre,dy_pre,dx_pos,dy_pos
      real *8 d_prepre,d_prepos,d_pospre,d_pospos
      real *8 A,B,C,D,alpha,beta
      logical E
      call periodize_r8(x_loc,min_x,max_x,x_loc)
      nx_d = n_x*(x_loc-min_x)/(max_x-min_x)
      nx_pre = floor(nx_d)
      dx_pre = dabs(nx_pre - nx_d)
      nx_pos = ceiling(nx_d)
      dx_pos = dabs(nx_pos - nx_d)
      if (nx_pos.ge.n_x) then
         nx_pos = nx_pos - n_x
      end if
      call periodize_r8(y_loc,min_y,max_y,y_loc)
      ny_d = n_y*(y_loc-min_y)/(max_y-min_y)
      ny_pre = floor(ny_d)
      dy_pre = dabs(ny_pre - ny_d)
      ny_pos = ceiling(ny_d)
      dy_pos = dabs(ny_pos - ny_d)
      if (ny_pos.ge.n_y) then
         ny_pos = ny_pos - n_y
      end if
      d_prepre = dsqrt((dx_pre)**2 + (dy_pre)**2)
      d_prepos = dsqrt((dx_pre)**2 + (dy_pos)**2)
      d_pospre = dsqrt((dx_pos)**2 + (dy_pre)**2)
      d_pospos = dsqrt((dx_pos)**2 + (dy_pos)**2)
      E=.true.
c$$$  case 1
      A = d_prepre
      B = d_prepos
      C = d_pospre
      D = d_pospos
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pre + ny_pre*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pre + ny_pre*n_x) + alpha*F_(nx_pre +
     $           ny_pos*n_x)
         end if
      end if
c$$$  case 2
      A = d_prepre
      B = d_pospre
      C = d_prepos
      D = d_pospos
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pre + ny_pre*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pre + ny_pre*n_x) + alpha*F_(nx_pos +
     $           ny_pre*n_x)
         end if
      end if
c$$$  case 3
      A = d_prepre
      B = d_pospos
      C = d_prepos
      D = d_pospre
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pre + ny_pre*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pre + ny_pre*n_x) + alpha*F_(nx_pos +
     $           ny_pos*n_x)
         end if
      end if
c$$$  case 4
      A = d_prepos
      B = d_pospre
      C = d_prepre
      D = d_pospos
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pre + ny_pos*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pre + ny_pos*n_x) + alpha*F_(nx_pos +
     $           ny_pre*n_x)
         end if
      end if
c$$$  case 5
      A = d_prepos
      B = d_pospos
      C = d_prepre
      D = d_pospre
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pre + ny_pos*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pre + ny_pos*n_x) + alpha*F_(nx_pos +
     $           ny_pos*n_x)
         end if
      end if
c$$$  case 6
      A = d_pospre
      B = d_pospos
      C = d_prepre
      D = d_prepos
      if (E .and. A.le.C .and. A.le.D .and. B.le.C .and. B.le.D) then
         E=.false.
         if (A+B.le.0) then
            output = F_(nx_pos + ny_pre*n_x)
         else
            alpha = A/(A+B)
            beta = B/(A+B)
            output = beta*F_(nx_pos + ny_pre*n_x) + alpha*F_(nx_pos +
     $           ny_pos*n_x)
         end if
      end if      
      end
