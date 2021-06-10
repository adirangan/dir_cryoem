      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_many_plan for omp sub-blocks'
         write(6,'(A)') ' '
      end if
      fpm_rank = 1
      fpm_howmany = min(fpm_howmany_max,n_transf*nS_1_per_max
     $     *nM_1_per_max)
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*fpm_rank*1.0d0*n_omp_sub__in*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_n__ requires ' , d_mem*1.0d-6
     $        ,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
      allocate(fpm_n__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_n__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub__in-1
      allocate(fpm_inembed__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_inembed__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub__in-1
      allocate(fpm_onembed__(0:fpm_rank*n_omp_sub__in-1))
      do nomp_sub=0,n_omp_sub__in-1
         fpm_onembed__(nomp_sub) = n_w_max
      enddo !do nomp_sub=0,n_omp_sub__in-1
      fpm_istride = 1
      fpm_idist = n_w_max
      fpm_ostride = 1
      fpm_odist = n_w_max
      n_fpm = fpm_howmany*fpm_n__(0)
      d_mem = 0.0d0
      d_mem = d_mem + 1.0d0*n_fpm*1.0d0*n_omp_sub__in*16.0d0
      if (verbose.gt.1 .or. verbose_mem.gt.0) then
         write (6,'(A,2(F16.8,A))') ' fpm_in1__ requires ' , d_mem*1.0d
     $        -6,' MB, ' , d_mem*1.0d-9 , ' GB'
      end if !if (verbose.gt.0) then
      allocate(fpm_in1__(0:n_fpm*n_omp_sub__in-1))
      allocate(fpm_out__(0:n_fpm*n_omp_sub__in-1))     
      allocate(fpm_back_(0:n_omp_sub__in-1))
      timing_tic = omp_get_wtime()
c$$$      We do not actually use the forward plan
      do nomp_sub=0,n_omp_sub__in-1      
         call dfftw_plan_many_dft(fpm_back_(nomp_sub),fpm_rank
     $        ,fpm_n__(nomp_sub),fpm_howmany,fpm_in1__(n_fpm*nomp_sub)
     $        ,fpm_inembed__(nomp_sub) ,fpm_istride,fpm_idist
     $        ,fpm_out__(n_fpm*nomp_sub),fpm_onembed__(nomp_sub)
     $        ,fpm_ostride ,fpm_odist ,FFTW_BACKWARD,FFTW_MEASURE)
      enddo !do nomp_sub=0,n_omp_sub__in-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc - timing_tic
      gnump_tot = (1.0d0*n_fpm)*(1.0d0*n_omp_sub__in)
      timing_total_fftw_plan = timing_total_fftw_plan + timing_tot
      gnump_total_fftw_plan = gnump_total_fftw_plan + gnump_tot
      if (verbose_timing.gt.0) then
         write(6,'(A,A,F8.3)') ' fftw_plan_many:'
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose_timing.gt.0) then
