!> Doxygen comment: ;\n
!> maps polar and azumithal angles to a vector on sphere. ;\n
      subroutine get_angle_to_vp_(alpha_polar_a,alpha_azimu_b,vp_)
      implicit none
      real *8 alpha_polar_a
      real *8 alpha_azimu_b
      real *8 vp_(0:2)
      vp_(0) = dcos(alpha_azimu_b)*dsin(alpha_polar_a)
      vp_(1) = dsin(alpha_azimu_b)*dsin(alpha_polar_a)
      vp_(2) = dcos(alpha_polar_a)
      end
