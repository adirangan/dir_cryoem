      subroutine rotate_c_to_c(n_x,max_x_c,n_y,max_y_c,S_c_,gamma
     $     ,M_c_)
c$$$      Assumes that M_c_ is the same size and dimensions as S_c_
c$$$      Rotates around (0,0)
c$$$      Does not periodize boundary
      implicit none
      integer n_x,n_y
      real *8 max_x_c,max_y_c,gamma
      complex *16 S_c_(0:n_x*n_y-1),M_c_(0:n_x*n_y-1)
      complex *16, allocatable :: T_c_(:)
      complex *16, allocatable :: N_c_(:)
      real *8 pi
      integer nx,ny
      real *8 X_c,Y_c,R_c,W_c
      complex *16 C_c
      allocate(T_c_(0:n_x*n_y-1))
      allocate(N_c_(0:n_x*n_y-1))
      call recenter_c16(n_x,n_y,S_c_,T_c_)
      pi = 4.0*atan(1.0)
      do ny=0,n_y-1
         do nx=0,n_x-1
            X_c = 0.0d0 + nx*max_x_c/n_x - max_x_c/2.0
            Y_c = 0.0d0 + ny*max_y_c/n_y - max_y_c/2.0
            R_c = dsqrt(X_c**2 + Y_c**2)
            W_c = atan2(Y_c,X_c) - gamma
            X_c = R_c*cos(W_c)
            Y_c = R_c*sin(W_c)
            X_c = max(0.0d0,min(max_x_c,X_c + max_x_c/2.0))
            Y_c = max(0.0d0,min(max_y_c,Y_c + max_y_c/2.0))
c$$$            call periodize_r8(X_c,0.0d0,max_x_c-d_x_c,X_c)
c$$$            call periodize_r8(Y_c,0.0d0,max_y_c-d_y_c,Y_c)
            call interp2_c16(n_x,0.0d0,max_x_c,n_y,0.0d0,max_y_c,T_c_
     $           ,X_c,Y_c,C_c)
            N_c_(nx+ny*n_x) = C_c
         enddo
      enddo
      call recenter_c16(n_x,n_y,N_c_,M_c_)
      deallocate(N_c_)
      deallocate(T_c_)      
      end
