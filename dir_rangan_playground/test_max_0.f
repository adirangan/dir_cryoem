!> Doxygen comment: ;\n
!> This program tests simple maximum timing. ;\n
      subroutine test_max_0(verbose,n_x,n_iteration)
      implicit none
      include 'omp_lib.h'
      integer *4 verbose
      integer *4 n_x,nx,nx_max
      integer *4 n_iteration,niteration
      real *8, allocatable :: x_(:)
      complex *16 x,x_max
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_max_0]'
      end if !if (verbose.gt.0) then
      allocate(x_(0:n_x-1))
      do nx=0,n_x-1
         x_(nx) = dcmplx( dsin(1.0d0*nx) , 0.0d0 )
      enddo !do nx=0,n_x-1
c$$$      %%%%%%%%%%%%%%%%
      timing_tic = omp_get_wtime()
      do niteration=0,n_iteration-1
      nx_max = 0
      x_max = x_(nx_max)
      do nx=0,n_x-1
         x = x_(nx)
         if (real(x).gt.real(x_max)) then
            nx_max = nx
            x_max = x
         end if !if (x.gt.x_max) then
      enddo !do nx=0,n_x-1
      enddo ! do niteration=0,n_iteration-1
      timing_toc = omp_get_wtime()
      if (verbose.gt.1) then
         write(6,'(A,I0,A,F8.3)') ' max brute force: ' , nx_max
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%%%%%%%%%
      timing_tic = omp_get_wtime()
      do niteration=0,n_iteration-1
      nx_max = 0
      x_max = x_(nx_max)
      nx_max = maxloc(x_,1)
      enddo ! do niteration=0,n_iteration-1
      timing_toc = omp_get_wtime()
      if (verbose.gt.1) then
         write(6,'(A,I0,A,F8.3)') ' max location: ' , nx_max
     $        ,' total_time ',timing_toc-timing_tic
      end if !if (verbose.gt.0) then
c$$$      %%%%%%%%%%%%%%%%

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_max_0]'
      end if !if (verbose.gt.0) then
      end
