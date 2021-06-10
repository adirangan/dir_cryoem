!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calls innerproduct_q_k_svdd_?. ;\n
      if (svd_calculation_type.eq.1) then
      if (verbose.gt.1) then
         write(6,'(A)') ' Depending on the value of flag_RTRT_vs_RTTR,'
         write(6,'(A)') ' either set up displacement-operator Z_S_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to S,'
         write(6,'(A)') ' or set up displacement-operator Z_M_svdd_'
         write(6,'(A)') ' associated with svd-expansion applied to M.'
         write(6,'(A)') ' '
      end if !if (verbose.gt.1) then
c$$$      svd_d_max = 0.0d0;
c$$$      do ndv=0,n_delta_v-1
c$$$         delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
c$$$         if (svd_d_max.lt.delta) then
c$$$            svd_d_max = delta
c$$$         end if !if (svd_d_max.lt.delta) then
c$$$      enddo !do ndv=0,n_delta_v-1
c$$$      if (verbose.gt.0) then
c$$$         write(6,'(A,F8.6)') ' setting svd_d_max: ' , svd_d_max
c$$$      end if !if (verbose.gt.0) then
      timing_tic = omp_get_wtime()
      if (flag_RTRT_vs_RTTR.eqv..true.) then
         call cl1_c16(n_svd_l*n_delta_v,Z_S_svdd_)
c$$$         call innerproduct_q_k_svdd_polyval_0on_0(flag_RTRT_vs_RTTR,svd_d_max
c$$$     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_
c$$$     $        ,delta_y_,Z_S_svdd_)
         call innerproduct_q_k_svdd_polyval_off_0(flag_RTRT_vs_RTTR
     $        ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $        ,svd_polyval_U_d_,n_delta_v,delta_x_ ,delta_y_,Z_S_svdd_)
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
c$$$         call innerproduct_q_k_svdd_polyval_0on_0(flag_RTRT_vs_RTTR,svd_d_max
c$$$     $        ,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_
c$$$     $        ,delta_y_,Z_M_svdd_)
         call innerproduct_q_k_svdd_polyval_off_0(flag_RTRT_vs_RTTR
     $        ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $        ,svd_polyval_U_d_,n_delta_v,delta_x_ ,delta_y_,Z_M_svdd_)
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
      if (verbose.gt.1) then
         write(6,'(A,F8.5)')
     $        ' finished ti8_define_svdd_3: total time '
     $        ,timing_tot
         timing_tmp = (n_svd_l*1.0d0)*(n_delta_v*1.0d0)/timing_tot/1e9
         write(6,'(A,F8.4)') ' ti8_define_svdd_3 total Gnumps: '
     $        ,timing_tmp
         write(6,'(A)') ' '
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      end if ! if (svd_calculation_type.eq.1) then


