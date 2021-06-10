      if (verbose.gt.0) then

         write(6,'(A,F8.4)') ' timing_total_fftw_plan: ' ,
     $        timing_total_fftw_plan
         write(6,'(A,F8.4)') ' gnump_total_fftw_plan: ' ,
     $        gnump_total_fftw_plan/max(1.0d0,timing_total_fftw_plan
     $        *1.0d9)

         write(6,'(A,F8.4)') ' timing_total_CTF_R_S: ' ,
     $        timing_total_CTF_R_S
         write(6,'(A,F8.4)') ' gnump_total_CTF_R_S: ' ,
     $        gnump_total_CTF_R_S/max(1.0d0,timing_total_CTF_R_S*1.0d9)

         if (flag_RTRT_vs_RTTR.eqv..false.) then
         write(6,'(A,F8.4)') ' timing_total_S_q: ' ,
     $        timing_total_S_q
         write(6,'(A,F8.4)') ' gnump_total_S_q: ' , gnump_total_S_q
     $        /max(1.0d0,timing_total_S_q*1.0d9)
         end if !if (flag_RTRT_vs_RTTR.eqv..false.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.))then
         write(6,'(A,F8.4)') ' timing_total_T_S_q: ' ,
     $        timing_total_T_S_q
         write(6,'(A,F8.4)') ' gnump_total_T_S_q: ' , gnump_total_T_S_q
     $        /max(1.0d0,timing_total_T_S_q*1.0d9)
         end if ! if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..true.)) then
         write(6,'(A,F8.4)') ' timing_total_Z_S_q: ' ,
     $        timing_total_Z_S_q
         write(6,'(A,F8.4)') ' gnump_total_Z_S_q: ' , gnump_total_Z_S_q
     $        /max(1.0d0,timing_total_Z_S_q*1.0d9)
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..true.)) then

         if (flag_RTRT_vs_RTTR.eqv..true.) then
         write(6,'(A,F8.4)') ' timing_total_T_R_CTF_M_q: ' ,
     $        timing_total_T_R_CTF_M_q
         write(6,'(A,F8.4)') ' gnump_total_T_R_CTF_M_q: ' ,
     $        gnump_total_T_R_CTF_M_q/max(1.0d0,timing_total_T_R_CTF_M_q
     $        *1.0d9)
         end if !if (flag_RTRT_vs_RTTR.eqv..true.) then

         if ((svd_calculation_type.eq.2) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         write(6,'(A,F8.4)') ' timing_total_T_T_R_CTF_M_q: ' ,
     $        timing_total_T_T_R_CTF_M_q
         write(6,'(A,F8.4)') ' gnump_total_T_T_R_CTF_M_q: ' ,
     $        gnump_total_T_T_R_CTF_M_q/max(1.0d0
     $        ,timing_total_T_T_R_CTF_M_q *1.0d9)
         end if !if ((svd_calculation_type.eq.2) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         if ((svd_calculation_type.eq.1) .and.
     $        (flag_RTRT_vs_RTTR.eqv..false.)) then
         write(6,'(A,F8.4)') ' timing_total_Z_T_R_CTF_M_q: ' ,
     $        timing_total_Z_T_R_CTF_M_q
         write(6,'(A,F8.4)') ' gnump_total_Z_T_R_CTF_M_q: ' ,
     $        gnump_total_Z_T_R_CTF_M_q/max(1.0d0
     $        ,timing_total_Z_T_R_CTF_M_q *1.0d9)
         end if !if ((svd_calculation_type.eq.1) .and. (flag_RTRT_vs_RTTR.eqv..false.)) then

         write(6,'(A,F8.4)') ' timing_total_zgemm: ' ,
     $        timing_total_zgemm
         write(6,'(A,F8.4)') ' gnump_total_zgemm: ' , gnump_total_zgemm
     $        /max(1.0d0,timing_total_zgemm*1.0d9)

         write(6,'(A,F8.4)') ' timing_total_fpm_fill: ' ,
     $        timing_total_fpm_fill
         write(6,'(A,F8.4)') ' gnump_total_fpm_fill: ' ,
     $        gnump_total_fpm_fill /max(1.0d0,timing_total_fpm_fill
     $        *1.0d9)

         write(6,'(A,F8.4)') ' timing_total_fpm_fftw: ' ,
     $        timing_total_fpm_fftw
         write(6,'(A,F8.4)') ' gnump_total_fpm_fftw: ' ,
     $        gnump_total_fpm_fftw /max(1.0d0,timing_total_fpm_fftw
     $        *1.0d9)

         write(6,'(A,F8.4)') ' timing_total_transpose: ' , 
     $        timing_total_transpose
         write(6,'(A,F8.4)') ' gnump_total_transpose: ' ,
     $        gnump_total_transpose /max(1.0d0,timing_total_transpose
     $        *1.0d9)

         write(6,'(A,F8.4)') ' timing_total_Zstore: ' , 
     $        timing_total_Zstore
         write(6,'(A,F8.4)') ' gnump_total_Zstore: ' ,
     $        gnump_total_Zstore /max(1.0d0,timing_total_Zstore
     $        *1.0d9)


      end if !if (verbose.gt.1) then

 10   end
