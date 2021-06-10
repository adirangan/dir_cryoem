      subroutine Fig_c16_polar(unitnumber,n_r,grid_x_p_,n_t_
     $     ,real_or_imag,F_ ,Fmin,Fmax,x_offset,y_offset,x_side,y_side)
      implicit none
      real *8 pi
      integer *4 unitnumber,n_r,n_t_(0:n_r-1)
      real *8 grid_x_p_(0:n_r-1)
      logical real_or_imag
      complex *16 F_(0:0)
      real *8 Fmax,Fmin
      real *8 x_offset,y_offset,x_side,y_side
      real *8 F_avg,F_std,Fmin_use,Fmax_use
      real *8 xval_tmp,yval_tmp
      real *8 rval_pre,tval_pre,rval_pos,tval_pos,tval_mid
      integer *4 diameter,arclength
      integer *4 nr,nt,n_A,na,narc,nd
      real *8, allocatable :: tvals_(:)
      integer *4, allocatable :: y_(:)
      integer *4, allocatable :: xyval_(:)
      character(len=20) format_string
      pi = 4*atan(1.0)
      n_A=0
      do nr=0,n_r-1
         n_A = n_A + n_t_(nr)
      enddo
      Fmin_use = Fmin
      Fmax_use = Fmax
      if (Fmin_use.ge.Fmax_use) then
         call stdlim_c16(n_A,real_or_imag,F_,F_avg,F_std)
         Fmin_use = F_avg - 2.0*F_std
         Fmax_use = F_avg + 2.0*F_std
      end if
      allocate(tvals_(0:n_A-1))
      na=0
      do nr=0,n_r-1
         do nt=0,n_t_(nr)-1
            tvals_(na) = 2.0*pi*(1.0*nt)/(1.0*n_t_(nr))
            na = na + 1;
         enddo
      enddo
      allocate(y_(0:n_A-1));
      na=0
      do nr=0,n_r-1
         do nt=0,n_t_(nr)-1
            if (real_or_imag.eqv..true.) then
               y_(na) = 32 + int(dmax1(0.0,dmin1(511.0,511.0
     $              *(real(F_(na))-Fmin_use)/(Fmax_use-Fmin_use))))
            else
               y_(na) = 32 + int(dmax1(0.0,dmin1(511.0,511.0
     $              *(aimag(F_(na))-Fmin_use)/(Fmax_use-Fmin_use))))
            end if
            na = na + 1;
         enddo
      enddo
      diameter = 1000
      arclength = 8
      allocate(xyval_(0:4*arclength+1))
      write(format_string,'(A,I0,A)') '(A,1X,',(4*arclength+2)
     $     ,'(I0,1X))'
      na=0
      do nr=0,n_r-1
         if (nr.gt.0) then
            rval_pre = 0.5*(grid_x_p_(nr) + grid_x_p_(nr-1))
         else
            rval_pre = grid_x_p_(0)
         end if
         if (nr.lt.n_r-1) then
            rval_pos = 0.5*(grid_x_p_(nr+1) + grid_x_p_(nr))
         else 
            rval_pos = grid_x_p_(n_r-1)
         end if
         do nt=0,n_t_(nr)-1
            if (nt.gt.0) then
               tval_pre = 0.5*(tvals_(na) + tvals_(na-1))
            else
               tval_pre = 0.5*(tvals_(na) + tvals_(na+n_t_(nr)-1) - 2
     $              *pi)
            end if
            if (nt.lt.n_t_(nr)-1) then
               tval_pos = 0.5*(tvals_(na+1) + tvals_(na))
            else
               tval_pos = 0.5*(tvals_(na-n_t_(nr)+1) + 2*pi +
     $              tvals_(na))
            end if
            do narc=0,arclength-1
               tval_mid = (narc*tval_pos + (arclength-1-narc)*tval_pre)
     $              /(arclength-1)
               xyval_(0 + 2*narc) = +int(diameter*(rval_pre
     $              *cos(tval_mid)*x_side + x_offset))
               xyval_(1 + 2*narc) = -int(diameter*(rval_pre
     $              *sin(tval_mid)*y_side + y_offset))
            enddo
            do narc=0,arclength-1
               tval_mid = (narc*tval_pre + (arclength-1-narc)*tval_pos)
     $              /(arclength-1)
               xyval_(0 + 2*narc + 2*arclength) = +int(diameter
     $              *(rval_pos*cos(tval_mid)*x_side + x_offset))
               xyval_(1 + 2*narc + 2*arclength) = -int(diameter
     $              *(rval_pos*sin(tval_mid)*y_side + y_offset))
            enddo
            xyval_(4*arclength+0) = +int(diameter*(rval_pre
     $           *cos(tval_pre)*x_side + x_offset))
            xyval_(4*arclength+1) = -int(diameter*(rval_pre
     $           *sin(tval_pre)*y_side + y_offset))
            write(unitnumber,'(A,I0,A,I0)') '2 3 0 0 0 ',y_(na)
     $           ,' 50 -1 20 0.000 0 0 -1 0 0 ',(2*arclength+1)
            write(unitnumber,format_string) char(9),(xyval_(nd),nd=0,4
     $           *arclength+1)
            na = na + 1;
         enddo
      enddo
      deallocate(xyval_)
      deallocate(y_)
      deallocate(tvals_)
      end

