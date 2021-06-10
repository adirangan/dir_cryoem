      subroutine rotate_delta(gamma,delta_0in_,delta_out_)
c$$$      This function simply applies the rotation matrix R_{+\gamma} to delta_0in_. ;
c$$$      That is to say, delta_out_ = R_{+\gamma} * delta_0in_. ;
c$$$      One aspect of this transformation that might be slightly confusing ;
c$$$      is that multiplying by R_{+\gamma} is actually the same as setting ;
c$$$      the angle of delta_out_ (say delta_out_p_w) to equal (delta_0in_p_w + gamma), ;
c$$$      which is the same thing as applying the transformation R_{-\gamma} to the ;
c$$$      distribution defined by delta_0in_ ;
c$$$      (i.e., the distribution \delta(x-\delta_0in_), which is peaked at delta_0in_). ;
c$$$      In summary, the important thing to keep straight is that, if: ;
c$$$      \delta(x-\tau_1) = R_{-\gamma} \delta(x-\tau_0) ;
c$$$      then ;
c$$$      \tau_1 = R_{+\gamma} \tau_0. ;
      implicit none
      real *8 gamma
      real *8 delta_0in_(0:1)
      real *8 delta_out_(0:1)
      real *8 c,s,d0,d1
      c = dcos(gamma)
      s = dsin(gamma)
      d0 = +c*delta_0in_(0) - s*delta_0in_(1)
      d1 = +s*delta_0in_(0) + c*delta_0in_(1)
      delta_out_(0) = d0
      delta_out_(1) = d1
      end !subroutine
      
     
