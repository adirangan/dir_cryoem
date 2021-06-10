      subroutine get_interchange_delta(delta_x_sort,delta_y_sort
     $     ,gamma_z_est,delta_x_tmp,delta_y_tmp)
      implicit none
      real *8 delta_x_sort,delta_y_sort,gamma_z_est
      real *8 delta_x_tmp,delta_y_tmp
      real *8 dc,ds
      dc = dcos(gamma_z_est)
      ds = dsin(gamma_z_est)
      delta_x_tmp = +dc*delta_x_sort +ds*delta_y_sort
      delta_y_tmp = -ds*delta_x_sort +dc*delta_y_sort
      end
