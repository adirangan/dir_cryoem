c$$$      reading in estimated translations and rotations
      if (verbose.gt.1) then
         write(6,'(A)') ' extracting from alpha_est__ '
      end if !if (verbose.gt.1) then
      call cl1_r8(nM_0_per,polar_a_est_)
      call cl1_r8(nM_0_per,azimu_b_est_)
      call cl1_r8(nM_0_per,gamma_z_est_)
      call cl1_r8(nM_0_per,delta_x_est_)
      call cl1_r8(nM_0_per,delta_y_est_)
      call cl1_r8(nM_0_per,l2_norm_est_)
      call cl1_r8(nM_0_per,ctf_ind_est_)
      call cl1_r8(nM_0_per,S_index_est_)
      call cl1_r8(nM_0_per,M_index_est_)
      do nm=0,nM_0_per-1
         polar_a_est_(nm) = alpha_est__(nalpha_polar_a + (nM_0_sum + nm)
     $        *n_alpha)
         azimu_b_est_(nm) = alpha_est__(nalpha_azimu_b + (nM_0_sum + nm)
     $        *n_alpha)
         gamma_z_est_(nm) = alpha_est__(nalpha_gamma_z + (nM_0_sum + nm)
     $        *n_alpha)
         delta_x_est_(nm) = alpha_est__(nalpha_delta_x + (nM_0_sum + nm)
     $        *n_alpha)
         delta_y_est_(nm) = alpha_est__(nalpha_delta_y + (nM_0_sum + nm)
     $        *n_alpha)
         l2_norm_est_(nm) = alpha_est__(nalpha_l2_norm + (nM_0_sum + nm)
     $        *n_alpha)
         ctf_ind_est_(nm) = alpha_est__(nalpha_ctf_ind + (nM_0_sum + nm)
     $        *n_alpha)
         S_index_est_(nm) = alpha_est__(nalpha_S_index + (nM_0_sum + nm)
     $        *n_alpha)
         M_index_est_(nm) = alpha_est__(nalpha_M_index + (nM_0_sum + nm)
     $        *n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,9F8.3)') (nM_0_sum + nm)
     $           ,'<-- nm: pa,ab,gz,dx,dy,l2,ci,ns,nm -->'
     $           ,polar_a_est_(nm) ,azimu_b_est_(nm),gamma_z_est_(nm)
     $           ,delta_x_est_(nm) ,delta_y_est_(nm),l2_norm_est_(nm)
     $           ,ctf_ind_est_(nm) ,S_index_est_(nm),M_index_est_(nm)
         end if !if (verbose.gt.1) then
      enddo !do nm=0,nM_0_per-1
      if (verbose.gt.1) then
         write(6,'(A)') ' finished extracting from alpha_est__ '
      end if !if (verbose.gt.1) then
