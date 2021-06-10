!> Doxygen comment: ;\n
!> Used multiple times when calling get_CTF_R_S_use. ;\n
!> Sets up linear interpolation for gamma_z. ;\n
      subroutine get_CTF_R_S_periodize_0(n_gamma_z,gamma_z_ ,CTF_R_S_
     $     ,gamma_z_est,ngz,CTF_R_S_use)
      implicit none
      integer n_gamma_z,ngz
      real *8 gamma_z_(0:n_gamma_z-1),gamma_z_est
      complex *16 CTF_R_S_(0:0),CTF_R_S_use
      real *8 ngz_est_v,dz_pre,dz_pos,alpha,beta
      integer ngz_est_pre,ngz_est_pos,ngz_pre,ngz_pos
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      include 'get_CTF_R_S_use_excerpt.f'
      ngz_pre = ngz + ngz_est_pre
      call periodize_i(ngz_pre,0,n_gamma_z,ngz_pre)
      ngz_pos = ngz + ngz_est_pos
      call periodize_i(ngz_pos,0,n_gamma_z,ngz_pos)
      CTF_R_S_use = beta*CTF_R_S_(ngz_pre) + alpha*CTF_R_S_(ngz_pos)
      end
