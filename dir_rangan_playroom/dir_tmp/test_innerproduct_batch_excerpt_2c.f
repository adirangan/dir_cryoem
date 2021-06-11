      subroutine test_innerproduct_batch_excerpt_2c(n_delta_x,delta_x_
     $     ,n_delta_y,delta_y_,n_gamma_z,gamma_z_,fftw_plan_frwd_
     $     ,fftw_in1_,fftw_out_,fftw_plan_back_last ,fftw_in1_last_
     $     ,fftw_out_last_,n_r,grid_p_,n_w_,n_A,C_M,C_S_,C_S_use_
     $     ,gamma_z_est ,S_p_,M_q_,Z_q_)
      implicit none
      include '/usr/include/fftw3.f'
      integer n_delta_x,n_delta_y,n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      real *8 gamma_z_est
      integer *8 fftw_plan_frwd_(0:n_r-1),fftw_plan_back_last
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 fftw_in1_last_(0:0),fftw_out_last_(0:0)
      complex *16 C_M,C_S_(0:0),Z_q,C_S_use_(0:0),C_Z_use
      complex *16 S_p_(0:0),M_q_(0:0),Z_q_(0:0)
      complex *16, allocatable :: T_p_(:)
      complex *16, allocatable :: T_q_(:)
      complex *16, allocatable :: Z_tmp_(:)
      integer n_w_max,nr,nC,ndx,ndy,ngz,nw,l
      real *8 delta_x,delta_y,gamma_z
      real *8 pi
      pi = 4.0*atan(1.0)
      call get_C_S_use_(gamma_z_est,n_gamma_z,C_S_,C_S_use_)
      n_w_max = n_w_(n_r-1)
      allocate(T_p_(0:n_A-1))
      allocate(T_q_(0:n_A-1))
      allocate(Z_tmp_(0:n_w_max-1))
      do ndy=0,n_delta_y-1
         delta_y = delta_y_(ndy)
         do ndx=0,n_delta_x-1
            delta_x = delta_x_(ndx)
            call cp1_c16(n_A,S_p_,T_p_)
            call transf_p_to_p(n_r,grid_p_,n_w_
     $           ,n_A,T_p_,+delta_x,+delta_y,T_p_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $           ,n_A,fftw_in1_,fftw_out_,T_p_,T_q_)
            call innerproduct_q__k_only(n_r,grid_p_,n_w_
     $           ,n_A,T_q_,M_q_,Z_tmp_)
            call cp1_c16(n_w_max,Z_tmp_,fftw_out_last_)
            call dfftw_execute_(fftw_plan_back_last)
            call cp1_c16(n_w_max,fftw_in1_last_,Z_tmp_)
            do ngz=0,n_gamma_z-1
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmp_,+gamma_z
     $              ,Z_q)
               nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               if (zabs(C_M*C_S_use_(ngz)).le.1.0d-15) then
                  C_Z_use = 1.0d0
               else
                  C_Z_use = C_M*C_S_use_(ngz)
               end if
               Z_q_(nC) = Z_q/(n_r**4)/C_Z_use
            enddo
         enddo
      enddo
      deallocate(Z_tmp_)
      deallocate(T_q_)
      deallocate(T_p_)
      end
