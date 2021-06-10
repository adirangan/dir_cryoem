      subroutine get_C_S_use_(gamma_z_est,n_gamma_z,C_S_,C_S_use_)
      implicit none
      real *8 gamma_z_est
      integer n_gamma_z
      complex *16 C_S_(0:0),C_S_use_(0:0)
      real *8 ngz_est_v,dz_pre,dz_pos,alpha,beta
      integer ngz_est_pre,ngz_est_pos,ngz_pre,ngz_pos
      integer ngz
      real *8 pi
      pi = 4.0*atan(1.0)
      include 'get_C_S_use_excerpt.f'
      do ngz=0,n_gamma_z-1
         ngz_pre = ngz + ngz_est_pre
         call periodize_i(ngz_pre,0,n_gamma_z,ngz_pre)
         ngz_pos = ngz + ngz_est_pos
         call periodize_i(ngz_pos,0,n_gamma_z,ngz_pos)
         C_S_use_(ngz) = beta*C_S_(ngz_pre) + alpha*C_S_(ngz_pos)
      enddo !do ngz=0,n_gamma_z-1
      end
