!> Doxygen comment: ;\n
!> checks memory allocation for test_residual_5.f ;\n
      flag_memory_checkset = .true.
c$$$      %%%%%%%%
      n = ld_M*n_M
      call cxs_r8(n,polar_a__,'polar_a__',flag_memory_checkset)
      n = ld_M*n_M
      call cxs_r8(n,azimu_b__,'azimu_b__',flag_memory_checkset)
      n = ld_M*n_M
      call cxs_c16(n,weight_CTF_k_c__,'weight_CTF_k_c__'
     $     ,flag_memory_checkset)
      n = n_k_p_r_cur
      call cxs_i4(n,n_w_M_,'n_w_M_',flag_memory_checkset)
      n = n_A
      call cxs_c16(n,fftw_0in_,'fftw_0in_',flag_memory_checkset)
      n = n_A
      call cxs_c16(n,fftw_out_,'fftw_out_',flag_memory_checkset)
      n = ld_M
      call cxs_c16(n,M_transform_,'M_transform_',flag_memory_checkset)
      n = ld_M
      call cxs_c16(n,Y_slice_,'Y_slice_',flag_memory_checkset)
      n = n_Y_lm_csum_cur*n_M
      call cxs_c16(n,H__,'H__',flag_memory_checkset)
      n = n_M
      call cxs_c16(n,HH_,'HH_',flag_memory_checkset)
      n = n_residual_loading*n_Y_lm_csum_cur
      call cxs_c16(n,G_,'G_',flag_memory_checkset)
c$$$      %%%%%%%%
      if ((flag_memory_checkset.eqv..true.).and.(verbose.gt.1)) then
         write(6,'(A)') '[checkset passed]'
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then
