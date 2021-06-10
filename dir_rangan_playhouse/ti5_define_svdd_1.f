      if (svd_calculation_type.eq.1) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Depending on the value of flag_RTRT_vs_RTTR,'
         write(6,'(A)') ' either set up displacement-operator Z_S_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to S,'
         write(6,'(A)') ' or set up displacement-operator Z_M_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to M.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then

      timing_tic = omp_get_wtime()
      if (flag_RTRT_vs_RTTR.eqv..true.) then
         call cl1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
         call innerproduct_q_k_svdd_2(flag_RTRT_vs_RTTR,n_svd_d,svd_d_
     $        ,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $        ,Z_S_svdd_)
c$$$      %%%%%%%%%%%%%%%%
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 2
         MDA_d_(0) = n_svd_l
         MDA_d_(1) = n_delta_v
         write(MDA_string,'(A)') './dir_mda/Z_S_svdd_.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,Z_S_svdd_,MDA_string)
      end if                    !if (flag_MDA.eqv..true.) then
c$$$      %%%%%%%%%%%%%%%%
      end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
      if (flag_RTRT_vs_RTTR.eqv..false.) then
         call cl1_c16(n_svd_l*n_delta_v,Z_M_svdd_)
         call innerproduct_q_k_svdd_2(flag_RTRT_vs_RTTR,n_svd_d,svd_d_
     $        ,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_
     $        ,Z_M_svdd_)
c$$$      %%%%%%%%%%%%%%%%
      flag_MDA = .false.
      if (flag_MDA.eqv..true.) then
         MDA_n_d = 2
         MDA_d_(0) = n_svd_l
         MDA_d_(1) = n_delta_v
         write(MDA_string,'(A)') './dir_mda/Z_M_svdd_.mda'
         call MDA_write_c16(MDA_n_d,MDA_d_,Z_M_svdd_,MDA_string)
      end if                    !if (flag_MDA.eqv..true.) then
c$$$      %%%%%%%%%%%%%%%%
      end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q_k_svdd_1: total time '
     $        ,timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_v*1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' innerproduct_q_k_svdd_1 total Gnumps: '
     $        ,timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      end if ! if (svd_calculation_type.eq.1) then
