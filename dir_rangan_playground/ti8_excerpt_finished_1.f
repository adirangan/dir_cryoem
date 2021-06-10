!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Prints out several parameters at the end of test_innerproduct_8. ;\n
      if (verbose.gt.0) then

         timing_tot = timing_total_fftw_plan
         gnump_tot = gnump_total_fftw_plan
         write(timing_string,'(A)') 'fftw_plan'
         include 'ti8_excerpt_timing_display_0.f'

         timing_tot = timing_total_CTF_R_S
         gnump_tot = gnump_total_CTF_R_S
         write(timing_string,'(A)') 'CTF_R_S'
         include 'ti8_excerpt_timing_display_0.f'

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_O_S_q
         gnump_tot = gnump_total_O_S_q
         write(timing_string,'(A)') 'S_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.))then
         timing_tot = timing_total_T_S_q
         gnump_tot = gnump_total_T_S_q
         write(timing_string,'(A)') 'T_S_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.)) then
         timing_tot = timing_total_Z_S_q
         gnump_tot = gnump_total_Z_S_q
         write(timing_string,'(A)') 'Z_S_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if (flag_RTRT_vs_RTTR.eqv..true.) then
         timing_tot = timing_total_O_T_R_CTF_M_q
         gnump_tot = gnump_total_O_T_R_CTF_M_q
         write(timing_string,'(A)') 'O_T_R_CTF_M_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_T_T_R_CTF_M_q
         gnump_tot = gnump_total_T_T_R_CTF_M_q
         write(timing_string,'(A)') 'T_T_R_CTF_M_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         timing_tot = timing_total_Z_T_R_CTF_M_q
         gnump_tot = gnump_total_Z_T_R_CTF_M_q
         write(timing_string,'(A)') 'Z_T_R_CTF_M_q'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         timing_tot = timing_total_zgemm
         gnump_tot = gnump_total_zgemm
         write(timing_string,'(A)') 'zgemm'
         include 'ti8_excerpt_timing_display_0.f'

         if (flag_fill.eqv..true.) then
         timing_tot = timing_total_fpm_fill
         gnump_tot = gnump_total_fpm_fill
         write(timing_string,'(A)') 'fpm_fill'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_fill.eqv..true.) then

         timing_tot = timing_total_fpm_fftw
         gnump_tot = gnump_total_fpm_fftw
         write(timing_string,'(A)') 'fpm_fftw'
         include 'ti8_excerpt_timing_display_0.f'

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         timing_tot = timing_total_transpose
         gnump_tot = gnump_total_transpose
         write(timing_string,'(A)') 'transpose'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then

         timing_tot = timing_total_Zstore
         gnump_tot = gnump_total_Zstore
         write(timing_string,'(A)') 'Zstore'
         include 'ti8_excerpt_timing_display_0.f'

         if (flag_time_Zstore) then
         timing_tot = timing_total_Zstore_a
         gnump_tot = gnump_total_Zstore_a
         write(timing_string,'(A)') 'Zstore_a'
         include 'ti8_excerpt_timing_display_0.f'
         timing_tot = timing_total_Zstore_b
         gnump_tot = gnump_total_Zstore_b
         write(timing_string,'(A)') 'Zstore_b'
         include 'ti8_excerpt_timing_display_0.f'
         timing_tot = timing_total_Zstore_c
         gnump_tot = gnump_total_Zstore_c
         write(timing_string,'(A)') 'Zstore_c'
         include 'ti8_excerpt_timing_display_0.f'
         timing_tot = timing_total_Zstore_x
         gnump_tot = gnump_total_Zstore_x
         write(timing_string,'(A)') 'Zstore_x'
         include 'ti8_excerpt_timing_display_0.f'
         timing_tot = timing_total_Zstore_d
         gnump_tot = gnump_total_Zstore_d
         write(timing_string,'(A)') 'Zstore_d'
         include 'ti8_excerpt_timing_display_0.f'
         timing_tot = timing_total_Zstore_e
         gnump_tot = gnump_total_Zstore_e
         write(timing_string,'(A)') 'Zstore_e'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_time_Zstore) then

         if (flag_MS_vs_SM.eqv..true.) then
         timing_tot = timing_total_merge_alpha_MS
         gnump_tot = gnump_total_merge_alpha_MS
         write(timing_string,'(A)') 'merge_alpha_MS'
         include 'ti8_excerpt_timing_display_0.f'
         end if !if (flag_MS_vs_SM.eqv..true.) then

      end if !if (verbose.gt.1) then

 10   continue
      
      if (verbose.gt.0) then
         if (verbose.gt.1) then
         if (flag_tesselation.eqv..true.) then
            call print_all_i4(n_M_9_sub_use,n_S_use_sum_
     $           ,' n_S_use_sum_: ')
            call print_all_i4(n_M,n_SM_use_local_
     $           ,' n_SM_use_local_: ')
         end if !if (flag_tesselation.eqv..true.) then
         end if !if (verbose.gt.1) then
         if (flag_tesselation.eqv..true.) then
            write(6,'(A,I0,A,I0,A)') ' n_SM_use_local_ in range: [' ,
     $           min_i4_f(n_M,n_SM_use_local_) , ',' , max_i4_f(n_M
     $           ,n_SM_use_local_) , ']'
         end if !if (flag_tesselation.eqv..true.) then
         write(6,'(A,I0,A,I0,A,F6.2,A)')
     $        ' [finished test_innerproduct_8]: calculating ' ,
     $        n_S_use_sum , ' out of ' , n_S*n_M ,
     $        ' image-template pairs: ' , 100.0d0*n_S_use_sum / (1.0d0
     $        *n_S*n_M) , '%'
      end if !if (verbose.gt.0) then

      end
