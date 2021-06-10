      program test_zgemm
c$$$      gfortran -o test_zgemm.out test_zgemm.f -fopenmp -lopenblas ; ./test_zgemm.out
      implicit none
      include 'omp_lib.h'
      complex *16, allocatable :: A_(:)
      complex *16, allocatable :: B_(:)
      complex *16, allocatable :: C1_(:)
      complex *16, allocatable :: C2_(:)
      complex *16 Z
      integer n_L,n_M,n_N,nl,nm,nn,na,nb,nc
      real *8 timing_tic,timing_toc,timing_tot
      real *8 timing_tmp

      write(6,'(A)') '[entering test_zgemm]'

      n_L = 130
      n_M = 160
      n_N = 120
      allocate(A_(0:n_L*n_M-1));
      allocate(B_(0:n_M*n_N-1));
      allocate(C1_(0:n_L*n_N-1));
      allocate(C2_(0:n_L*n_N-1));

      write(6,'(A)') 'initializing...'
c$$$      n_L = 3; n_M = 6; n_N = 2;
c$$$      A_ = zeros(n_M,n_L); for (na=0:n_M*n_L-1); A_(1+na) = mod(na,3)-1 + i*(mod(na,5)-2); end;
c$$$      B_ = zeros(n_M,n_N); for (nb=0:n_M*n_N-1); B_(1+nb) = mod(nb,7)-3 + i*(mod(nb,11)-5); end;
      
      timing_tic = omp_get_wtime()
      do na=0,n_L*n_M-1
         A_(na) = cmplx(0.0d0 + mod(na,3) - 1 , 0.0d0  + mod(na,5) - 2)
      enddo
      do nb=0,n_M*n_N-1
         B_(nb) = cmplx(0.0d0 + mod(nb,7) - 3 , 0.0d0  + mod(nb,11) - 5)
      enddo
      do nc=0,n_L*n_N-1
         C1_(nc) = (0.0d0,0.0d0)
         C2_(nc) = (0.0d0,0.0d0)
      enddo
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') 'total time: ' , timing_tot
c$$$      write(6,*) (A_(na),na=0,n_L*n_M-1)
c$$$      write(6,*) (B_(nb),nb=0,n_M*n_N-1)

c$$$      write(6,'(A)') 'zgemm multiplying...'
c$$$      timing_tic = omp_get_wtime()
c$$$      call zgemm('N','N',n_L,n_N,n_M,1.0d0,A_,n_L,B_,n_M,0.0d0,C1_,n_L)
c$$$      timing_toc = omp_get_wtime()
c$$$      timing_tot = timing_toc-timing_tic
c$$$      write(6,'(A,F8.4)') 'total time: ' , timing_tot
c$$$      timing_tmp = (n_L*1.0d0)*(n_M*1.0d0)*(n_N*1.0d0)/timing_tot/1e9
c$$$      write(6,'(A,F8.4)') 'total Gflops: ' , timing_tmp

      write(6,'(A)') 'zgemm multiplying...'
      timing_tic = omp_get_wtime()
      call zgemm('T','N',n_L,n_N,n_M,1.0d0*cmplx(1.0d0,0.0d0),A_,n_M,B_
     $     ,n_M,0.0d0*cmplx(1.0d0,0.0d0),C1_,n_L)
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') 'total time: ' , timing_tot
      timing_tmp = (n_L*1.0d0)*(n_M*1.0d0)*(n_N*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') 'total Gflops: ' , timing_tmp
c$$$      write(6,*) (C1_(nc),nc=0,n_L*n_N-1)

      write(6,'(A)') 'slow multiplying...'
      timing_tic = omp_get_wtime()
      nc = 0
      do nn=0,n_N-1
         do nl=0,n_L-1
            na = 0 + nl*n_M
            nb = 0 + nn*n_M
            Z = (0.0d0,0.0d0)
            do nm=0,n_M-1
               Z = Z+A_(na)*B_(nb)
               na = na+1
               nb = nb+1
            enddo !do nm=0,n_M-1
            C2_(nc) = Z
            nc = nc+1
         enddo !do nn=0,n_N-1
      enddo !do nl=0,n_L-1
      timing_toc = omp_get_wtime()
      timing_tot = timing_toc-timing_tic
      write(6,'(A,F8.4)') 'total time: ' , timing_tot
      timing_tmp = (n_L*1.0d0)*(n_M*1.0d0)*(n_N*1.0d0)/timing_tot/1e9
      write(6,'(A,F8.4)') 'total Gflops: ' , timing_tmp
c$$$      write(6,*) (C2_(nc),nc=0,n_L*n_N-1)

      Z = (0.0d0,0.0d0)
      do nc=0,n_L*n_N-1
         Z = Z + zabs(C1_(nc)-C2_(nc))
      enddo
      write(6,'(A,F16.8)') 'error: ' , zabs(Z)
      
      write(6,'(A)') '[finished test_zgemm]'
      stop
      end;
