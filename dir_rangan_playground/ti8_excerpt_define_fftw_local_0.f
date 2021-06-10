!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Generating fftw_plans for local use. ;\n
      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for local use'
         write(6,'(A)') ' '
      end if
      allocate(fftw_plan_frwd_(0:n_r-1))
      allocate(fftw_plan_back_(0:n_r-1))
      allocate(fftw_0in_(0:n_w_sum-1))
      allocate(fftw_out_(0:n_w_sum-1))
      timing_tic = omp_get_wtime()
      na = 0
      do nr=0,n_r-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nr),n_w_(nr)
     $        ,fftw_0in_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nr),n_w_(nr)
     $        ,fftw_out_(na),fftw_0in_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nr)
      enddo
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_w_sum)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_1d:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose_timing.gt.0) then
      p_fftw_plan_frwd_last = loc(fftw_plan_frwd_(n_r-1))
      p_fftw_plan_back_last = loc(fftw_plan_back_(n_r-1))
      p_fftw_0in_last_ = loc(fftw_0in_(n_w_sum-n_w_max))
      p_fftw_out_last_ = loc(fftw_out_(n_w_sum-n_w_max))
