      if (verbose.gt.0) then

         timing_tot = timing_total_fftw_plan
         gnump_tot = gnump_total_fftw_plan
         write(timing_string,'(A)') 'fftw_plan'
         include 'ti5_timing_display_0.f'

         timing_tot = timing_total_CTF_R_S
         gnump_tot = gnump_total_CTF_R_S
         write(timing_string,'(A)') 'CTF_R_S'
         include 'ti5_timing_display_0.f'

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_O_S_q
         gnump_tot = gnump_total_O_S_q
         write(timing_string,'(A)') 'S_q'
         include 'ti5_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.))then
         timing_tot = timing_total_T_S_q
         gnump_tot = gnump_total_T_S_q
         write(timing_string,'(A)') 'T_S_q'
         include 'ti5_timing_display_0.f'
         end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.)) then
         timing_tot = timing_total_Z_S_q
         gnump_tot = gnump_total_Z_S_q
         write(timing_string,'(A)') 'Z_S_q'
         include 'ti5_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if (flag_RTRT_vs_RTTR.eqv..true.) then
         timing_tot = timing_total_O_T_R_CTF_M_q
         gnump_tot = gnump_total_O_T_R_CTF_M_q
         write(timing_string,'(A)') 'O_T_R_CTF_M_q'
         include 'ti5_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_T_T_R_CTF_M_q
         gnump_tot = gnump_total_T_T_R_CTF_M_q
         write(timing_string,'(A)') 'T_T_R_CTF_M_q'
         include 'ti5_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_Z_T_R_CTF_M_q
         gnump_tot = gnump_total_Z_T_R_CTF_M_q
         write(timing_string,'(A)') 'Z_T_R_CTF_M_q'
         include 'ti5_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         timing_tot = timing_total_zgemm
         gnump_tot = gnump_total_zgemm
         write(timing_string,'(A)') 'zgemm'
         include 'ti5_timing_display_0.f'

         timing_tot = timing_total_fpm_fill
         gnump_tot = gnump_total_fpm_fill
         write(timing_string,'(A)') 'fpm_fill'
         include 'ti5_timing_display_0.f'

         timing_tot = timing_total_fpm_fftw
         gnump_tot = gnump_total_fpm_fftw
         write(timing_string,'(A)') 'fpm_fftw'
         include 'ti5_timing_display_0.f'

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_transpose
         gnump_tot = gnump_total_transpose
         write(timing_string,'(A)') 'transpose'
         include 'ti5_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then

         timing_tot = timing_total_Zstore
         gnump_tot = gnump_total_Zstore
         write(timing_string,'(A)') 'Zstore'
         include 'ti5_timing_display_0.f'

         if (flag_time_Zstore) then
         timing_tot = timing_total_Zstore_a
         gnump_tot = gnump_total_Zstore_a
         write(timing_string,'(A)') 'Zstore_a'
         include 'ti5_timing_display_0.f'
         timing_tot = timing_total_Zstore_b
         gnump_tot = gnump_total_Zstore_b
         write(timing_string,'(A)') 'Zstore_b'
         include 'ti5_timing_display_0.f'
         timing_tot = timing_total_Zstore_c
         gnump_tot = gnump_total_Zstore_c
         write(timing_string,'(A)') 'Zstore_c'
         include 'ti5_timing_display_0.f'
         timing_tot = timing_total_Zstore_x
         gnump_tot = gnump_total_Zstore_x
         write(timing_string,'(A)') 'Zstore_x'
         include 'ti5_timing_display_0.f'
         timing_tot = timing_total_Zstore_d
         gnump_tot = gnump_total_Zstore_d
         write(timing_string,'(A)') 'Zstore_d'
         include 'ti5_timing_display_0.f'
         timing_tot = timing_total_Zstore_e
         gnump_tot = gnump_total_Zstore_e
         write(timing_string,'(A)') 'Zstore_e'
         include 'ti5_timing_display_0.f'
         end if !if (flag_time_Zstore) then

      end if !if (verbose.gt.1) then

 10   continue
      
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished test_innerproduct_5]'
      end if !if (verbose.gt.0) then

      end
