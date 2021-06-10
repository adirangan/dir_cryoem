      subroutine interp1_nufft1df77_c16(n_x,min_x,max_x,F_
     $     ,fftw_plan_frwd,fftw_plan_back,fftw_0in_,fftw_out_
     $     ,nufft1df77_fw_,nufft1df77_lw,nufft1df77_lused,n_x_loc,x_loc_
     $     ,output_)
c$$$      assumes regular grid for x
c$$$      periodic boundary conditions
      implicit none
      include '/usr/include/fftw3.f'
      include 'omp_lib.h'
      integer verbose
      data verbose / 0 /
      integer *4 n_x,n_x_loc
      real *8 min_x,max_x
      complex *16 F_(0:n_x-1)
      real *8 x_loc_(0:0)
      complex *16 output_(0:0)
      integer *4 n_x_2,nx,nx_loc
      real *8, allocatable :: y_loc_(:)
      complex *16, allocatable :: G_(:)
      complex *16, allocatable :: H_(:)
      complex *16, allocatable :: J_(:)
      integer *8 fftw_plan_frwd
      integer *8 fftw_plan_back
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      real *8 nufft1df77_fw_(0:0)
      integer *4 nufft1df77_lw,nufft1df77_lused
      real *8 pi
      integer iflag,ier
      real *8 eps
      logical flag_memory_checkset
      pi = 4.0d0*datan(1.0d0)
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering interp1_nufft1df77_c16] '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' n_x: ' , n_x
         write(6,'(A,F8.4)') ' min_x: ' , min_x
         write(6,'(A,F8.4)') ' max_x: ' , max_x
         write(6,'(A,I0)') ' n_x_loc: ' , n_x_loc
         call print_mid_r8(n_x_loc,x_loc_,' x_loc_: ')
         call print_sub_c16(n_x,F_,' F_: ')
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A)') ' allocating memory'
      end if !if (verbose.gt.0) then
      allocate(y_loc_(0:1+n_x_loc-1))
      call cs1_r8(n_x_loc,y_loc_)
      allocate(G_(0:1+n_x-1))
      call cs1_c16(n_x,G_)
      allocate(H_(0:1+n_x-1))
      call cs1_c16(n_x,H_)
      allocate(J_(0:1+n_x_loc-1))
      call cs1_c16(n_x_loc,J_)
      call periodize_r8_(n_x_loc,x_loc_,min_x,max_x,y_loc_) !wrap y_loc_ into [x_min,x_max). ;
      if (verbose.gt.0) then
         write(6,'(A)') ' y_loc_ in [x_min,x_max) : '
         call print_mid_r8(n_x_loc,y_loc_,' y_loc_: ')
      end if !if (verbose.gt.0) then
      call af1_r8(n_x_loc,1.0d0/(max_x-min_x),-min_x/(max_x-min_x)
     $     ,y_loc_,y_loc_) !map y_loc_ to [0,1) interval. 
      if (verbose.gt.0) then
         write(6,'(A)') ' y_loc_ in [0,1) : '
         call print_mid_r8(n_x_loc,y_loc_,' y_loc_: ')
      end if !if (verbose.gt.0) then
      do nx_loc=0,n_x_loc-1
         if (y_loc_(nx_loc).ge.0.5d0) then
            y_loc_(nx_loc) = y_loc_(nx_loc)-1.0d0 !wrap y_loc_ into [-0.5,0.5) interval. ;
         end if !if (y_loc_(nx_loc).ge.0.5d0) then
      enddo !do nx_loc=0,n_x_loc-1
      if (verbose.gt.0) then
         write(6,'(A)') ' y_loc_ in [-0.5d0,+0.5d0) : '
         call print_mid_r8(n_x_loc,y_loc_,' y_loc_: ')
      end if !if (verbose.gt.0) then
      call af1_r8(n_x_loc,2.0d0*pi,0.0d0
     $     ,y_loc_,y_loc_) !map y_loc_ to [-pi,+pi) interval. ;
      if (verbose.gt.0) then
         write(6,'(A)') ' y_loc_ in [-pi,+pi) : '
         call print_mid_r8(n_x_loc,y_loc_,' y_loc_: ')
      end if !if (verbose.gt.0) then
      call cp1_c16(n_x,F_,fftw_0in_)
      call dfftw_execute(fftw_plan_frwd)
      call af1_c16(n_x,1.0d0*dcmplx(1.0d0/dsqrt(1.0d0*n_x)),1.0d0
     $     *dcmplx(0.0d0),fftw_out_,G_) ! frequency range 2.0d0*pi*[0,n_x-1]
      if (verbose.gt.0) then
         write(6,'(A)') ' G_ frequencies in 2.0d0*pi*[0,n_x-1] '
         call print_mid_r8(n_x,G_,' G_: ')
      end if !if (verbose.gt.0) then
      n_x_2 = ceiling(n_x/2.0d0)
      do nx=0,n_x_2-1
         H_(n_x - n_x_2 + nx) = G_(nx) !wrap G_ to run from frequency range pi*[-n_x/2,+n_x/2-1]. ;
      enddo !do nx=0,n_x_2-1
      do nx=n_x_2,n_x-1
         H_(nx-n_x_2) = G_(nx) !wrap G_ to run from frequency range pi*[-n_x/2,+n_x/2-1]. ;
      enddo !do nx=0,n_x_2-1
      if (verbose.gt.0) then
         write(6,'(A)') ' H_ frequencies in pi*[-n_x/2,+n_x/2-1] '
         call print_mid_r8(n_x,G_,' G_: ')
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
      call print_mid_r8(n_x_loc,y_loc_,' y_loc_: ')
      call print_sub_c16(n_x,H_,' H_: ')
      end if !if (verbose.gt.0) then
      iflag = +1
      eps = 1.0d-12
      call nufft1d2(n_x_loc,y_loc_,J_,iflag,eps,n_x,H_,nufft1df77_fw_
     $     ,nufft1df77_lw,nufft1df77_lused,ier)
      if (ier.ne.0) then
         write(6,'(A)') 'Warning! ier.ne.0 in interp1_nufft1df77_c16.f'
      end if !if (ier.ne.0) then
      call af1_c16(n_x_loc,1.0d0*dcmplx(1.0d0/dsqrt(1.0d0*n_x)),1.0d0
     $     *dcmplx(0.0d0),J_,output_) 
      if (verbose.gt.0) then
         call print_sub_c16(n_x_loc,J_,' J_: ')
         call print_sub_c16(n_x_loc,output_,' output_: ')
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A)') ' checking memory allocation'
      end if !if (verbose.gt.0) then
      flag_memory_checkset=.true.
      call cxs_r8(n_x_loc,y_loc_,'y_loc_',flag_memory_checkset)
      call cxs_c16(n_x,G_,'G_',flag_memory_checkset)
      call cxs_c16(n_x,H_,'H_',flag_memory_checkset)
      call cxs_c16(n_x_loc,J_,'J_',flag_memory_checkset)
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error
      end if !if (flag_memory_checkset.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished interp1_nufft1df77_c16] '
      end if !if (verbose.gt.0) then
      end !subroutine
