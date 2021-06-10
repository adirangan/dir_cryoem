      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
      end if ! if (verbose.gt.1) then
      allocate(fftw_plan_frwd__(0:n_r*n_omp_sub__in-1))
      allocate(fftw_plan_back__(0:n_r*n_omp_sub__in-1))
      allocate(fftw_in1__(0:n_A*n_omp_sub__in-1))
      allocate(fftw_out__(0:n_A*n_omp_sub__in-1))
      timing_tic = omp_get_wtime()
      do nomp_sub=0,n_omp_sub__in-1
         n_r_tmp = nomp_sub*n_r
         n_A_tmp = nomp_sub*n_A
         na = 0
         do nr=0,n_r-1
            call dfftw_plan_dft_1d_(fftw_plan_frwd__(n_r_tmp+nr)
     $           ,n_w_(nr),fftw_in1__(n_A_tmp+na),fftw_out__(n_A_tmp
     $           +na),FFTW_FORWARD,FFTW_ESTIMATE) 
            call dfftw_plan_dft_1d_(fftw_plan_back__(n_r_tmp+nr)
     $           ,n_w_(nr),fftw_out__(n_A_tmp+na),fftw_in1__(n_A_tmp
     $           +na),FFTW_BACKWARD,FFTW_ESTIMATE) 
            na = na + n_w_(nr)
         enddo ! do nr=0,n_r-1
      enddo ! do nomp_sub=0,n_omp_sub__in-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_A)*(1.0d0*n_omp_sub__in)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_omp:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose_timing.gt.0) then
