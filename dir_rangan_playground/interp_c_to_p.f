!> Doxygen comment: ;\n
!> Interpolate from cartesian to polar coordinates. ;\n
      subroutine interp_c_to_p(n_x,max_x_c,n_y,max_y_c,S_c_,n_r
     $     ,grid_p_,n_w_,n_A,S_p_)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_x,n_y,n_r,n_A
      real *8 max_x_c,max_y_c
      complex *16 S_c_(0:n_x*n_y-1),S_p_(0:n_A-1)
      real *8 grid_p_(0:n_r-1)
      integer n_w_(0:n_r-1)
      complex *16, allocatable :: T_c_(:)
      real *8 pi,max_r_c
      integer ic,nr,nw
      real *8 X_c,Y_c,R_c,W_c
      complex *16 C_c
      allocate(T_c_(0:n_x*n_y-1));
      call recenter_c16(n_x,n_y,S_c_,T_c_)
      pi = 4.0d0*datan(1.0d0)
      max_r_c = dsqrt(max_x_c**2 + max_y_c**2)
      ic=0
      do nr=0,n_r-1
         R_c = grid_p_(nr)
         if (verbose.gt.0) then
            write(6,*) 'R_c(',nr,') = ',R_c
         end if
         do nw=0,n_w_(nr)-1
            W_c = 0.0 + nw*(2.0d0*pi)/(n_w_(nr))
            X_c = R_c*dcos(W_c)
            Y_c = R_c*dsin(W_c)
            call interp2_c16(n_x,-max_x_c/2.0,max_x_c/2.0,n_y,-max_y_c
     $           /2.0,max_y_c/2.0,T_c_,X_c,Y_c ,C_c)
            S_p_(ic) = C_c
            ic = ic + 1
         enddo
      enddo
      deallocate(T_c_)
      end
