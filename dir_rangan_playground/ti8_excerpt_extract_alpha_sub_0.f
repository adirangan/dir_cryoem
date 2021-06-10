!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Extracts image-parameters from stack of image-parameters. ;\n
!> Extracts only across images nM_0_sum  --> nM_0_sum + nM_0_per ;\n
!> I.e., extracts parameters for the current 0-level block of images. ;\n
c$$$      reading in estimated translations and rotations
      if (verbose.gt.1) then
         write(6,'(A)') ' extracting from alpha_est__ '
         write(6,'(A)') ' (extracting across nM_0_per) '
      end if !if (verbose.gt.1) then
      call cl1_r8(nM_0_per_max,polar_a_est_)
      call cl1_r8(nM_0_per_max,azimu_b_est_)
      call cl1_r8(nM_0_per_max,gamma_z_est_)
      call cl1_r8(nM_0_per_max,delta_x_est_)
      call cl1_r8(nM_0_per_max,delta_y_est_)
      call cl1_r8(nM_0_per_max,l2_norm_est_)
      call cl1_r8(nM_0_per_max,ctf_ind_est_)
      call cl1_r8(nM_0_per_max,S_index_est_)
      call cl1_r8(nM_0_per_max,M_index_est_)
      do nm0=0,nM_0_per-1
         nm1 = nM_0_sum + nm0
         polar_a_est_(nm0) = alpha_est__(nalpha_polar_a + nm1*n_alpha)
         azimu_b_est_(nm0) = alpha_est__(nalpha_azimu_b + nm1*n_alpha)
         gamma_z_est_(nm0) = alpha_est__(nalpha_gamma_z + nm1*n_alpha)
         delta_x_est_(nm0) = alpha_est__(nalpha_delta_x + nm1*n_alpha)
         delta_y_est_(nm0) = alpha_est__(nalpha_delta_y + nm1*n_alpha)
         l2_norm_est_(nm0) = alpha_est__(nalpha_l2_norm + nm1*n_alpha)
         ctf_ind_est_(nm0) = alpha_est__(nalpha_ctf_ind + nm1*n_alpha)
         S_index_est_(nm0) = alpha_est__(nalpha_S_index + nm1*n_alpha)
         M_index_est_(nm0) = alpha_est__(nalpha_M_index + nm1*n_alpha)
         if (verbose.gt.2) then
            write(6,'(I2,A,9F8.3,I0)') (nm0)
     $           ,'<-- nm0: pa,ab,gz,dx,dy,l2,ci,ns,nm,nm1 -->'
     $           ,polar_a_est_(nm0) ,azimu_b_est_(nm0),gamma_z_est_(nm0)
     $           ,delta_x_est_(nm0) ,delta_y_est_(nm0),l2_norm_est_(nm0)
     $           ,ctf_ind_est_(nm0) ,S_index_est_(nm0),M_index_est_(nm0)
     $           ,nm1
         end if                 !if (verbose.gt.1) then
      enddo                     !do nM0=0,nM_0_per-1
      if (verbose.gt.1) then
         write(6,'(A)') ' finished extracting from alpha_est__ '
         write(6,'(A)') ' (extracting across nM_0_per) '
      end if                    !if (verbose.gt.1) then
