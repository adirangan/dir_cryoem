!> Doxygen comment: ;\n
!> Applies translation in x_c (i.e., real-space cartesian) coordinates. ;\n
!> This amounts to a shift by some (fractional) number of pixels. ;\n
!> Uses linear interpolation. ;\n
!> Assumes that M_c_ is the same size and dimensions as S_c_. ;\n
      subroutine transl_c_to_c(n_x,max_x_c,n_y,max_y_c,S_c_,delta_x
     $     ,delta_y,M_c_)
c$$$      Assumes that M_c_ is the same size and dimensions as S_c_
      implicit none
      integer n_x,n_y
      real *8 max_x_c,max_y_c,delta_x,delta_y
      complex *16 S_c_(0:n_x*n_y-1),M_c_(0:n_x*n_y-1)
      complex *16, allocatable :: T_c_(:)
      complex *16, allocatable :: N_c_(:)
      integer nx,ny
      real *8 X_c,Y_c
      complex *16 C_c
      allocate(T_c_(0:n_x*n_y-1))
      allocate(N_c_(0:n_x*n_y-1))
      call recenter_c16(n_x,n_y,S_c_,T_c_)
      do ny=0,n_y-1
         do nx=0,n_x-1
            X_c = 0.0d0 + nx*max_x_c/n_x - delta_x
            Y_c = 0.0d0 + ny*max_y_c/n_y - delta_y
            call interp2_c16(n_x,0.0d0,max_x_c,n_y,0.0d0,max_y_c,T_c_
     $           ,X_c,Y_c,C_c)
            N_c_(nx+ny*n_x) = C_c
         enddo
      enddo
      call decenter_c16(n_x,n_y,N_c_,M_c_)
      deallocate(N_c_)
      deallocate(T_c_)
      end
