      subroutine get_interchange_delta(delta_x_pre,delta_y_pre
     $     ,gamma_z,delta_x_pos,delta_y_pos)
c$$$      This function simply applies the rotation ;
c$$$      R_{-gamma_z} to the vector delta_pre, ;
c$$$      where: ;
c$$$      R_{-gamma_z} = [ +cos(gamma_z) , +sin(gamma_z) ] ;
c$$$                     [ -sin(gamma_z) , +cos(gamma_z) ] ;
c$$$      and
c$$$      delta_pre = [ delta_x_pre ; delta_y_pre ]. ;
c$$$      ;
c$$$      if we define T_{delta} to be translation by delta, ;
c$$$      we can see that the above function is used when converting: ;
c$$$      R_{gamma2}*T_{delta2}*R_{gamma1}*T_{delta1} ;
c$$$      into ;
c$$$      R_{gamma2}*R_{gamma1}T_{R_{-gamma1}*delta2}*T_{delta1} ;
c$$$      which equals ;
c$$$      R_{gamma2 + gamma1}T_{R_{-gamma1}*delta2 + delta1}. ;
c$$$      ;
c$$$      A conversion of this kind is useful when updating image-parameters.

      implicit none
      real *8 delta_x_pre,delta_y_pre,gamma_z
      real *8 delta_x_pos,delta_y_pos
      real *8 delta_x_tmp,delta_y_tmp
      real *8 dc,ds
      dc = dcos(gamma_z)
      ds = dsin(gamma_z)
      delta_x_tmp = delta_x_pre
      delta_y_tmp = delta_y_pre
      delta_x_pos = +dc*delta_x_tmp +ds*delta_y_tmp
      delta_y_pos = -ds*delta_x_tmp +dc*delta_y_tmp
      end
