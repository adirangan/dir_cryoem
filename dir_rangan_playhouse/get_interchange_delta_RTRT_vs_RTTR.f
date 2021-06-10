      subroutine get_interchange_delta_RTRT_vs_RTTR(flag_RTRT_vs_RTTR
     $     ,delta_x_est,delta_y_est,gamma_z_est,delta_x_upd,delta_y_upd
     $     ,gamma_z_upd)
      implicit none
      logical flag_RTRT_vs_RTTR
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 delta_x_upd,delta_y_upd,gamma_z_upd
      real *8 delta_x_tmp,delta_y_tmp,gamma_z_tmp
      if (flag_RTRT_vs_RTTR.eqv..true.) then
         call get_interchange_delta(delta_x_est,delta_y_est,gamma_z_upd
     $        ,delta_x_tmp ,delta_y_tmp)
         gamma_z_upd = gamma_z_upd + gamma_z_est
         delta_x_upd = delta_x_upd + delta_x_tmp
         delta_y_upd = delta_y_upd + delta_y_tmp
      end if                    !if (flag_RTRT_vs_RTTR.eqv..true.) then
      if (flag_RTRT_vs_RTTR.eqv..false.) then
         delta_x_tmp = delta_x_est + delta_x_upd
         delta_y_tmp = delta_y_est + delta_y_upd
         call get_interchange_delta(delta_x_tmp
     $        ,delta_y_tmp,gamma_z_upd,delta_x_upd
     $        ,delta_y_upd)
         gamma_z_upd = gamma_z_upd + gamma_z_est
      end if                    !if (flag_RTRT_vs_RTTR.eqv..false.) then
      end
