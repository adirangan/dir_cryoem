      subroutine test_innerproduct_batch_excerpt_4(n_delta_x,delta_x_
     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r,grid_p_,n_w_,n_A
     $     ,C_M,C_S_,gamma_z_est ,S_p_,M_p_,Z_p_)
      implicit none
      integer n_delta_x,n_delta_y,n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      real *8 gamma_z_est
      complex *16 C_M,C_S_(0:0),Z_p,C_S_use,C_Z_use
      complex *16 S_p_(0:0),M_p_(0:0),Z_p_(0:0)
      complex *16, allocatable :: T_p_(:)
      integer nr,nC,ndx,ndy,ngz
      real *8 delta_x,delta_y,gamma_z
      real *8 pi
      complex *16, allocatable :: C_S_use_(:)
      pi = 4.0*atan(1.0)
      allocate(C_S_use_(0:n_gamma_z-1))
      call get_C_S_use_(gamma_z_est,n_gamma_z,C_S_,C_S_use_)
      allocate(T_p_(0:n_A-1))
      do ngz=0,n_gamma_z-1
         gamma_z = gamma_z_(ngz)
         C_S_use = C_S_use_(ngz)
         if (zabs(C_M*C_S_use).le.1.0d-15) then
            C_Z_use = 1.0d0
         else
            C_Z_use = C_M*C_S_use
         end if
         do ndy=0,n_delta_y-1
            delta_y = delta_y_(ndy)
            do ndx=0,n_delta_x-1
               delta_x = delta_x_(ndx)
               call cp1_c16(n_A,S_p_,T_p_)
               call transf_p_to_p(n_r,grid_p_,n_w_,n_A,T_p_,+delta_x,
     $              +delta_y,T_p_)
               call rotate_p2p_fz(n_r,n_w_,n_A,T_p_,+gamma_z,T_p_)
               call innerproduct_p(n_r,grid_p_,n_w_,n_A,T_p_,M_p_,Z_p)
               nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               Z_p_(nC) = Z_p/(n_r**4)/C_Z_use
            enddo
         enddo
      enddo
      deallocate(T_p_)
      end
