!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Uses omp to call blocks of ti8_merge_alpha_MS_0. ;\n
      if (flag_MS_vs_SM.eqv..true.) then
          if (verbose.gt.1) then
            write(6,'(A)') ' Merge alpha_MS_omp__ from threads. '
            write(6,'(A)') ' Given thread nM_9_sub, the array in'
            write(6,'(A)')
     $           ' alpha_MS_omp__(n_alpha*n_MS_max*n_S*nM_9_sub)'
            write(6,'(A)') ' must be merged with the output array'
            write(6,'(A)') ' alpha_MS__(0).'
            write(6,'(A)') ' We also must merge the values in'
            write(6,'(A)') ' n_MS_omp__(n_S*nM_9_sub)'
            write(6,'(A)') ' with those in'
            write(6,'(A)') ' n_MS_(0).'
         end if !if (verbose.gt.1) then
         timing_tic = omp_get_wtime()
         if (n_S_9_sub_use.gt.1) then
c$OMP PARALLEL PRIVATE(nS_9_per,nS_9_sum)
c$OMP DO 
         do nS_9_sub=0,n_S_9_sub_use-1
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_merge_alpha_MS_0'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_merge_alpha_MS_0(
     $           verbose-1
     $           ,n_S
     $           ,nS_9_per
     $           ,nS_0_sum + nS_9_sum
     $           ,n_MS_max
     $           ,n_omp_sub_0in
     $           ,n_MS_omp__
     $           ,n_MS_
     $           ,alpha_MS_omp__
     $           ,alpha_MS__
     $           )
         enddo !do nS_9_sub=0,n_S_9_sub_use-1
c$OMP END DO
c$OMP END PARALLEL
         else !if (n_S_9_sub_use.gt.1) then
            nS_9_sub = 0
            nS_9_per = n_S_9_per_(nS_9_sub)
            nS_9_sum = n_S_9_sum_(nS_9_sub)
            if (verbose.gt.1) then
               write(6,'(A,A,A,I0)') 'calling'
     $              ,' ti8_merge_alpha_MS_0'
     $              ,' with thread: ',omp_get_thread_num()
            end if
            call ti8_merge_alpha_MS_0(
     $           verbose-1
     $           ,n_S
     $           ,nS_9_per
     $           ,nS_0_sum + nS_9_sum
     $           ,n_MS_max
     $           ,n_omp_sub_0in
     $           ,n_MS_omp__
     $           ,n_MS_
     $           ,alpha_MS_omp__
     $           ,alpha_MS__
     $           )
         end if !if (n_S_9_sub_use.gt.1) then
         timing_toc = omp_get_wtime()
         timing_tot = timing_toc - timing_tic
         gnump_tot = (n_MS_max*1.0d0)*(nS_9_per*1.0d0)
         timing_total_merge_alpha_MS = timing_total_merge_alpha_MS +
     $        timing_tot
         gnump_total_merge_alpha_MS = gnump_total_merge_alpha_MS +
     $        gnump_tot
         timing_tmp = gnump_tot/timing_tot/1e9
         if (verbose_timing.gt.0) then
            write(6,'(A,A,F8.5)')
     $           ' finished calculating merge_alpha_MS__ for each S. ' ,
     $           ' total time ' , timing_tot
            write(6,'(A,F8.4)') ' merge_alpha_MS__ total Gnumps: ' ,
     $           timing_tmp
            write(6,'(A)') ' '
         end if ! if (verbose_timing.gt.0) then
      end if !if (flag_MS_vs_SM.eqv..true.) then

