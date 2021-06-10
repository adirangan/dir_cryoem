      ngz_est_v = n_gamma_z*(gamma_z_est)/(2*pi)
      ngz_est_pre = floor(ngz_est_v)
      dz_pre = dabs(ngz_est_pre - ngz_est_v)
      ngz_est_pos = ceiling(ngz_est_v)
      dz_pos = dabs(ngz_est_pos - ngz_est_v)
      if (ngz_est_pos.eq.n_gamma_z) then
         ngz_est_pos = 0
      end if
      if (dz_pre+dz_pos.le.0.0d0) then
         alpha = 0.0d0
         beta = 1.0d0
      else
         alpha = dz_pre/(dz_pre+dz_pos)
         beta =  dz_pos/(dz_pre+dz_pos)
      end if !dz
