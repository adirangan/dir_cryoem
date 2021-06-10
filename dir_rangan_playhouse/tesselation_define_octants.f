      subroutine tesselation_define_octants(v_n00_,v_0n0_,v_00n_,v_p00_
     $     ,v_0p0_,v_00p_)
      implicit none
      real *8 v_n00_(0:2)
      real *8 v_0n0_(0:2)
      real *8 v_00n_(0:2)
      real *8 v_p00_(0:2)
      real *8 v_0p0_(0:2)
      real *8 v_00p_(0:2)
      v_n00_(0) = -1
      v_n00_(1) = 0
      v_n00_(2) = 0
      v_0n0_(0) = 0
      v_0n0_(1) = -1
      v_0n0_(2) = 0
      v_00n_(0) = 0
      v_00n_(1) = 0
      v_00n_(2) = -1
      v_p00_(0) = +1
      v_p00_(1) = 0
      v_p00_(2) = 0
      v_0p0_(0) = 0
      v_0p0_(1) = +1
      v_0p0_(2) = 0
      v_00p_(0) = 0
      v_00p_(1) = 0
      v_00p_(2) = +1
      end
