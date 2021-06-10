      subroutine Fig_c16_bess2(unitnumber,n_r,grid_x_p_,n_t_
     $     ,real_or_imag,F_ ,Fmin,Fmax,x_offset,y_offset,x_side,y_side)
c$$$      expanding q-axis for visibility
c$$$      rows no longer align equivalent q-values as r varies
      implicit none
      integer *4 unitnumber,n_r,n_t_(0:n_r-1)
      real *8 grid_x_p_(0:n_r-1)
      logical real_or_imag
      complex *16 F_(0:0)
      real *8 Fmax,Fmin
      real *8 x_offset,y_offset,x_side,y_side
      real *8 F_avg,F_std,Fmin_use,Fmax_use
      real *8 r_pre,r_pos,q_pre,q_pos
      integer *4 diameter
      integer *4 nr,nt,n_A,na,ntmax,nd
      integer *4 xyval_(0:9)
      integer *4, allocatable :: y_(:)
      character(len=20) format_string
      n_A=0
c$$$      ntmax is not actually used
      ntmax=0
      do nr=0,n_r-1
         n_A = n_A + n_t_(nr)
         if (n_t_(nr).gt.ntmax) then
            ntmax = n_t_(nr)
         end if
      enddo
      Fmin_use = Fmin
      Fmax_use = Fmax
      if (Fmin_use.ge.Fmax_use) then
         call stdlim_c16(n_A,real_or_imag,F_,F_avg,F_std)
         Fmin_use = F_avg - 2.0*F_std
         Fmax_use = F_avg + 2.0*F_std
      end if
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
      write(format_string,'(A,I0,A)') '(A,1X,',10,'(I0,1X))'
      na=0
      do nr=0,n_r-1
         if (nr.gt.0) then
            r_pre = 2.0*0.5*(grid_x_p_(nr) + grid_x_p_(nr-1))
         else
            r_pre = 2.0*grid_x_p_(0)
         end if
         if (nr.lt.n_r-1) then
            r_pos = 2.0*0.5*(grid_x_p_(nr+1) + grid_x_p_(nr))
         else 
            r_pos = 2.0*grid_x_p_(n_r-1)
         end if
         do nt=0,n_t_(nr)-1
            if (nt.ge.n_t_(nr)/2) then
               q_pre = 0.5 + 0.5*(1.0d0*(nt+0-1*n_t_(nr)))/(0.5d0
     $              *n_t_(nr))
               q_pos = 0.5 + 0.5*(1.0d0*(nt+1-1*n_t_(nr)))/(0.5d0
     $              *n_t_(nr))
            else
               q_pre = 0.5 + 0.5*(1.0d0*(nt+0-0*n_t_(nr)))/(0.5d0
     $              *n_t_(nr))
               q_pos = 0.5 + 0.5*(1.0d0*(nt+1-0*n_t_(nr)))/(0.5d0
     $              *n_t_(nr))
            end if
            xyval_(0 + 2*0) = +int(diameter*(r_pre*x_side + x_offset))
            xyval_(1 + 2*0) = -int(diameter*(q_pre*y_side + y_offset))
            xyval_(0 + 2*1) = +int(diameter*(r_pos*x_side + x_offset))
            xyval_(1 + 2*1) = -int(diameter*(q_pre*y_side + y_offset))
            xyval_(0 + 2*2) = +int(diameter*(r_pos*x_side + x_offset))
            xyval_(1 + 2*2) = -int(diameter*(q_pos*y_side + y_offset))
            xyval_(0 + 2*3) = +int(diameter*(r_pre*x_side + x_offset))
            xyval_(1 + 2*3) = -int(diameter*(q_pos*y_side + y_offset))
            xyval_(0 + 2*4) = +int(diameter*(r_pre*x_side + x_offset))
            xyval_(1 + 2*4) = -int(diameter*(q_pre*y_side + y_offset))
            write(unitnumber,'(A,I0,A,I0)') '2 3 0 0 0 ',y_(na)
     $           ,' 50 -1 20 0.000 0 0 -1 0 0 ',5
            write(unitnumber,format_string) char(9),(xyval_(nd),nd=0,9)
            na = na + 1;
         enddo
      enddo
      deallocate(y_)
      end

