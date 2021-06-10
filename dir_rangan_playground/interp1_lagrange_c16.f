!> Doxygen comment: ;\n
!> assumes regular grid for x
!> periodic boundary conditions
!> uses n_q==6 nodes.
      subroutine interp1_lagrange_c16(n_x,min_x,max_x,F_,n_x_loc,x_loc_
     $     ,output_)
c$$$      assumes regular grid for x ;
c$$$      periodic boundary conditions ;
c$$$      uses n_q==6 nodes. ;
c$$$      Note that, if n_x_loc is 0, ;
c$$$      then x_loc_ is read as a single real *8, ;
c$$$      and is treated as a shift to the regular periodic grid ;
c$$$      defined by n_x, min_x and max_x. ;
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 n_x,n_x_loc
      real *8 min_x,max_x
      complex *16 F_(0:n_x-1)
      real *8 x_loc_(0:0)
      complex *16 output_(0:0)
      integer n_x_loc_use
      integer nt_pre,nt_pos,nx_loc
      integer n_q,nq
      parameter(n_q=6)
      integer *4 nt_(0:n_q-1)
      complex *16 G_(0:n_q-1)
      real *8 q_(0:n_q-1)
      real *8 x_(0:n_q-1)
      real *8 d_(0:n_q-1)
      real *8 e_(0:n_q-1)
      real *8 x_loc,nt_r8,x_tmp
      x_(0) = 0.0d0
      x_(1) = 1.0d0
      x_(2) = 2.0d0
      x_(3) = 3.0d0
      x_(4) = 4.0d0
      x_(5) = 5.0d0
      d_(0) = ((x_(0)-x_(1))*(x_(0)-x_(2)) *(x_(0)-x_(3))*(x_(0)-x_(4))
     $     *(x_(0)-x_(5)))
      d_(1) = ((x_(1)-x_(0))*(x_(1)-x_(2)) *(x_(1)-x_(3))*(x_(1)-x_(4))
     $     *(x_(1)-x_(5)))
      d_(2) = ((x_(2)-x_(0))*(x_(2)-x_(1)) *(x_(2)-x_(3))*(x_(2)-x_(4))
     $     *(x_(2)-x_(5)))
      d_(3) = ((x_(3)-x_(0))*(x_(3)-x_(1)) *(x_(3)-x_(2))*(x_(3)-x_(4))
     $     *(x_(3)-x_(5)))
      d_(4) = ((x_(4)-x_(0))*(x_(4)-x_(1)) *(x_(4)-x_(2))*(x_(4)-x_(3))
     $     *(x_(4)-x_(5)))
      d_(5) = ((x_(5)-x_(0))*(x_(5)-x_(1)) *(x_(5)-x_(2))*(x_(5)-x_(3))
     $     *(x_(5)-x_(4)))
      e_(0) = 1.0d0/d_(0)
      e_(1) = 1.0d0/d_(1)
      e_(2) = 1.0d0/d_(2)
      e_(3) = 1.0d0/d_(3)
      e_(4) = 1.0d0/d_(4)
      e_(5) = 1.0d0/d_(5)
      n_x_loc_use = n_x_loc
      if (n_x_loc.le.0) then
         n_x_loc_use = n_x
      end if !if (n_x_loc.le.0) then
         
      do nx_loc=0,n_x_loc_use-1
         if (n_x_loc.le.0) then
            x_loc = min_x + (1.0d0*nx_loc*(max_x-min_x))/(1.0d0*n_x) -
     $           x_loc_(0)
         end if !if (n_x_loc.le.0) then
         if (n_x_loc.gt.0) then
         call periodize_r8(x_loc_(nx_loc),min_x,max_x,x_loc)
         end if !if (n_x_loc.gt.0) then
         nt_r8 = (x_loc-min_x)/(max_x-min_x)*n_x
         nt_pre = floor(nt_r8)
         nt_(0) = nt_pre - ceiling(n_q*0.5d0) + 1 + 0
         nt_(1) = nt_pre - ceiling(n_q*0.5d0) + 1 + 1
         nt_(2) = nt_pre - ceiling(n_q*0.5d0) + 1 + 2
         nt_(3) = nt_pre - ceiling(n_q*0.5d0) + 1 + 3
         nt_(4) = nt_pre - ceiling(n_q*0.5d0) + 1 + 4
         nt_(5) = nt_pre - ceiling(n_q*0.5d0) + 1 + 5
         x_tmp = nt_r8 - nt_(0)
         call periodize_i(nt_(0),0,n_x,nt_(0))
         call periodize_i(nt_(1),0,n_x,nt_(1))
         call periodize_i(nt_(2),0,n_x,nt_(2))
         call periodize_i(nt_(3),0,n_x,nt_(3))
         call periodize_i(nt_(4),0,n_x,nt_(4))
         call periodize_i(nt_(5),0,n_x,nt_(5))
         G_(0) = F_(nt_(0))
         G_(1) = F_(nt_(1))
         G_(2) = F_(nt_(2))
         G_(3) = F_(nt_(3))
         G_(4) = F_(nt_(4))
         G_(5) = F_(nt_(5))
         q_(0) = ((x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp-x_(3))*(x_tmp
     $        -x_(4))*(x_tmp-x_(5)))*e_(0)
         q_(1) = ((x_tmp-x_(0))*(x_tmp-x_(2))*(x_tmp-x_(3))*(x_tmp
     $        -x_(4))*(x_tmp-x_(5)))*e_(1)
         q_(2) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(3))*(x_tmp
     $        -x_(4))*(x_tmp-x_(5)))*e_(2)
         q_(3) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $        -x_(4))*(x_tmp-x_(5)))*e_(3)
         q_(4) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $        -x_(3))*(x_tmp-x_(5)))*e_(4)
         q_(5) = ((x_tmp-x_(0))*(x_tmp-x_(1))*(x_tmp-x_(2))*(x_tmp
     $        -x_(3))*(x_tmp-x_(4)))*e_(5)
         output_(nx_loc) = G_(0)*q_(0) + G_(1)*q_(1) + G_(2)*q_(2) +
     $        G_(3)*q_(3) + G_(4)*q_(4) + G_(5)*q_(5)
      enddo !do nx_loc=0,n_x_loc_use-1

      end !subroutine
