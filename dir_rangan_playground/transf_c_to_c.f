!> Doxygen comment: ;\n
!> Applies translation in k_c (i.e., fourier-space cartesian) coordinates. ;\n
!> This amounts to a multiplication by a plane-wave. ;\n
!> Assumes that M_c_ is the same size and dimensions as S_c_. ;\n
!> Note that this recenters the image before applying the plane-wave. ;\n
!> i.e., we assume that the 'origin' is in the center of the image. ;\n
      subroutine transf_c_to_c(n_x,max_x_c,n_y,max_y_c,S_c_,delta_x
     $     ,delta_y,M_c_)
c$$$      Assumes that M_c_ is the same size and dimensions as S_c_
c$$$      Multiplication *should* be performed in place 
c$$$      (i.e., without temporary workspace)
      implicit none
      integer n_x,n_y
      real *8 max_x_c,max_y_c,delta_x,delta_y
      complex *16 S_c_(0:n_x*n_y-1),M_c_(0:n_x*n_y-1)
      complex *16, allocatable :: T_c_(:)
      complex *16, allocatable :: N_c_(:)
      real *8 pi
      integer nx,ny
      real *8 X_c,Y_c,L_c
      complex *16 C_c
      pi = 4.0d0*datan(1.0d0)
      allocate(T_c_(0:n_x*n_y-1))
      allocate(N_c_(0:n_x*n_y-1))
      call recenter_c16(n_x,n_y,S_c_,T_c_)
      do ny=0,n_y-1
         do nx=0,n_x-1
            X_c = 0.0d0 + nx*max_x_c/n_x - max_x_c/2.0
            Y_c = 0.0d0 + ny*max_y_c/n_y - max_y_c/2.0
            L_c = (X_c * delta_x) + (Y_c * delta_y)
            C_c = dcmplx(+dcos(2.0d0*pi*L_c),-dsin(2.0d0*pi*L_c))
            N_c_(nx+ny*n_x) = C_c*T_c_(nx+ny*n_x)
         enddo
      enddo
      call decenter_c16(n_x,n_y,N_c_,M_c_)
      deallocate(N_c_)
      deallocate(T_c_)
      end
