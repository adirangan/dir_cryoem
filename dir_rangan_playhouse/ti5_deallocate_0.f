      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nr=0,n_r-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nr))
         call dfftw_destroy_plan(fftw_plan_back_(nr))
      enddo !do nr=0,n_r-1
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_in1_)
      deallocate(fftw_out_)

      if (verbose.gt.1) then
         write(6,'(A)') 'destroying fftw_plans for omp sub-blocks.'
      end if ! if (verbose.gt.1) then
      do nomp_sub=0,n_omp_sub__in-1
         n_r_tmp = nomp_sub*n_r
         do nr=0,n_r-1
            call dfftw_destroy_plan(fftw_plan_frwd__(n_r_tmp+nr))
            call dfftw_destroy_plan(fftw_plan_back__(n_r_tmp+nr))
         enddo ! do nr=0,n_r-1
      enddo ! do nomp_sub=0,n_omp_sub__in-1
      deallocate(fftw_plan_frwd__)
      deallocate(fftw_plan_back__)
      deallocate(fftw_in1__)
      deallocate(fftw_out__)

      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_many_plans for omp.'
      end if
      do nomp_sub=0,n_omp_sub__in-1
         call dfftw_destroy_plan_(fpm_back_(nomp_sub))
c$$$         call dfftw_destroy_plan_(fpm_frwd_(nomp_sub))
      enddo !do nomp_sub=0,n_omp_sub__in-1
      deallocate(fpm_in1__)
      deallocate(fpm_out__)
      deallocate(fpm_n__)
      deallocate(fpm_inembed__)
      deallocate(fpm_onembed__)

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_5]'
      end if
