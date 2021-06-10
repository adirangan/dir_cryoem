      subroutine get_angle_to_vp_(alpha_polar_a,alpha_azimu_b,vp_)
      implicit none
      real *8 alpha_polar_a
      real *8 alpha_azimu_b
      real *8 vp_(0:2)
      vp_(0) = cos(alpha_azimu_b)*sin(alpha_polar_a)
      vp_(1) = sin(alpha_azimu_b)*sin(alpha_polar_a)
      vp_(2) = cos(alpha_polar_a)
      end
