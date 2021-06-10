      subroutine test_innerproduct_3d_0(verbose,n_beta,beta_,n_k,k_,n_l_
     $     ,a_,b_,X_)
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      integer *4 n_beta,n_k,n_l_(0:n_k-1)
      real *8 beta_(0:n_beta-1),k_(0:n_k-1)
      complex *16 a_(0:0),b_(0:0),X_(0:0)
      integer *4 n_l_max,n_m_max
      integer *4 max_i4_f
c$$$      fftw for omp
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1__(:)
      complex *16, allocatable :: fftw_out__(:)
      integer *4 n_beta_sub,nbeta_sub
      integer *4 l_base,f_base,i_base,n_base,nl,nm,na
      integer *4, allocatable :: n_beta_sub_(:)
      integer *4, allocatable :: I_beta_sub_(:)
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
c$$$      display
      character(len=64) format_string

      if (verbose.gt.0) then
         write (6,'(A)') ' [entering test_innerproduct_3d_0]'
      end if

      n_l_max = max_i4_f(n_k,n_l_)
      n_m_max = 1 + 2*n_l_max

      if (verbose.gt.0) then
         write (6,'(A,I0)') ' n_l_max: ' , n_l_max
         write (6,'(A,I0)') ' n_m_max: ' , n_m_max
      end if

      allocate(n_beta_sub_(0:n_beta-1))
      allocate(I_beta_sub_(0:n_beta-1))
      n_beta_sub = min(8,n_beta)
      if (verbose.gt.0) then
         write (6,'(A,I0,A)') 'Setting up n_beta_sub = ',n_beta_sub
     $        ,' sub-blocks for omp parallelization'
      end if
      do nbeta_sub=0,n_beta_sub-1
         if (nbeta_sub.eq.0) then
            n_beta_sub_(0) = n_beta/n_beta_sub
            I_beta_sub_(0) = 0
         else if (nbeta_sub.gt.0 .and. nbeta_sub.lt.n_beta_sub-1) then
            n_beta_sub_(nbeta_sub) = n_beta /n_beta_sub
            I_beta_sub_(nbeta_sub) = I_beta_sub_(nbeta_sub-1) +
     $           n_beta_sub_(nbeta_sub-1)
         else if (nbeta_sub.eq.n_beta_sub-1) then
            I_beta_sub_(n_beta_sub-1) = I_beta_sub_(n_beta_sub-2) +
     $           n_beta_sub_(n_beta_sub-2)
            n_beta_sub_(n_beta_sub-1) = n_beta - I_beta_sub_(n_beta_sub
     $           -1)
         end if
      enddo
      if (verbose.gt.0) then
         write(format_string,'(A,I0,A)') '(A,',n_beta_sub
     $        ,'(I0,1X))'
         write(6,format_string) 'n_beta_sub_: '
     $        ,(n_beta_sub_(nbeta_sub),nbeta_sub=0
     $        ,n_beta_sub-1)
         write(6,format_string) 'I_beta_sub_: '
     $        ,(I_beta_sub_(nbeta_sub),nbeta_sub=0
     $        ,n_beta_sub-1)
      end if

      if (verbose.gt.1) then
         write(6,'(A)') 'Generating fftw_plans for omp sub-blocks.'
      end if
      allocate(fftw_plan_frwd_(0:n_beta_sub-1))
      allocate(fftw_plan_back_(0:n_beta_sub-1))
      allocate(fftw_in1__(0:n_m_max*n_m_max*n_beta_sub-1))
      allocate(fftw_out__(0:n_m_max*n_m_max*n_beta_sub-1))
      do nbeta_sub=0,n_beta_sub-1
         l_base = nbeta_sub
         f_base = nbeta_sub*n_m_max*n_m_max
         call dfftw_plan_dft_2d_(fftw_plan_frwd_(l_base),n_m_max
     $        ,n_m_max,fftw_in1__(f_base),fftw_out__(f_base)
     $        ,FFTW_FORWARD,FFTW_ESTIMATE) 
         call dfftw_plan_dft_2d_(fftw_plan_back_(l_base),n_m_max
     $        ,n_m_max,fftw_out__(f_base),fftw_in1__(f_base)
     $        ,FFTW_BACKWARD,FFTW_ESTIMATE) 
      enddo
      
      timing_tic = omp_get_wtime()
c$OMP PARALLEL PRIVATE(n_base,i_base,l_base,f_base)
c$OMP DO 
      do nbeta_sub=0,n_beta_sub-1
         n_base = n_beta_sub_(nbeta_sub)
         i_base = I_beta_sub_(nbeta_sub)
         l_base = nbeta_sub
         f_base = nbeta_sub*n_m_max*n_m_max
         if (verbose.gt.1) then
            write (6,'(A,I0)') ' n_base: ' , n_base
            write (6,'(A,I0)') ' i_base: ' , i_base
            write (6,'(A,I0)') ' l_base: ' , l_base
            write (6,'(A,I0)') ' f_base: ' , f_base
         end if
         if (verbose.gt.1) then
            write(6,'(A,A,A,I0)') 'calling'
     $           ,' test_innerproduct_3d_stage_0'
     $           ,' with thread: ',omp_get_thread_num()
         end if
         call test_innerproduct_3d_stage_0(verbose,n_base
     $        ,beta_(i_base),n_k,k_,n_l_,a_,b_,X_(n_m_max*n_m_max
     $        *i_base),fftw_plan_frwd_(l_base),fftw_plan_back_(l_base)
     $        ,fftw_in1__(f_base) ,fftw_out__(f_base))
      enddo !do nbeta_sub=0,n_M_sub-1
c$OMP END DO
c$OMP END PARALLEL
      timing_toc = omp_get_wtime()

      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for omp.'
      end if
      do nbeta_sub=0,n_beta_sub-1
         l_base = nbeta_sub
         call dfftw_destroy_plan_(fftw_plan_back_(l_base))
         call dfftw_destroy_plan_(fftw_plan_frwd_(l_base))
      enddo
      deallocate(fftw_out__)
      deallocate(fftw_in1__)
      deallocate(fftw_plan_back_)
      deallocate(fftw_plan_frwd_)
      deallocate(n_beta_sub_)
      deallocate(I_beta_sub_)

      if (verbose.gt.0) then
         write(6,'(A,A,A,F8.5)') ' [finished',
     $        ' test_innerproduct_batch_stage_1]', 
     $        ' total time ',
     $        timing_toc-timing_tic
      end if

      end
