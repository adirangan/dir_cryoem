      subroutine interp2(n_x,min_x,max_x,n_y,min_y,max_y,F_,x_loc,y_loc
     $     ,output)
c$$$      assumes regular grid for x and y
      integer *4 n_x,n_y
      real *8 min_x,max_x,min_y,max_y
      real *8 F_(0:n_x*n_y-1),x_loc,y_loc
      real *8 output
      integer *4 nx_pre,ny_pre,nx_pos,ny_pos
      real *8 x_pre,y_pre,x_pos,y_pos
      real *8 d_prepre,d_prepos,d_pospre,d_pospos
      real *8 A,B,C,D,alpha,beta
      logical E
      if (x_loc.ge.max_x) then
         nx_pre = n_x-1
         nx_pos = n_x-1
      else if (x_loc.le.min_x) then
         nx_pre = 0
         nx_pos = 0
      else 
         nx_pre = floor(n_x*(x_loc-min_x)/(max_x-min_x))
         nx_pre = max(0,min(n_x-1,nx_pre))
         nx_pos = ceiling(n_x*(x_loc-min_x)/(max_x-min_x))
         nx_pos = max(0,min(n_x-1,nx_pos))
      end if
      x_pre = min_x + nx_pre*(max_x-min_x)/n_x
      x_pos = min_x + nx_pos*(max_x-min_x)/n_x
      if (y_loc.ge.max_y) then
         ny_pre = n_y-1
         ny_pos = n_y-1
      else if (y_loc.le.min_y) then
         ny_pre = 0
         ny_pos = 0
      else 
         ny_pre = floor(n_y*(y_loc-min_y)/(max_y-min_y))
         ny_pre = max(0,min(n_y-1,ny_pre))
         ny_pos = ceiling(n_y*(y_loc-min_y)/(max_y-min_y))
         ny_pos = max(0,min(n_y-1,ny_pos))
      end if
      y_pre = min_y + ny_pre*(max_y-min_y)/n_y
      y_pos = min_y + ny_pos*(max_y-min_y)/n_y
      d_prepre = dsqrt((x_loc-x_pre)**2 + (y_loc-y_pre)**2)
      d_prepos = dsqrt((x_loc-x_pre)**2 + (y_loc-y_pos)**2)
      d_pospre = dsqrt((x_loc-x_pos)**2 + (y_loc-y_pre)**2)
      d_pospos = dsqrt((x_loc-x_pos)**2 + (y_loc-y_pos)**2)
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
