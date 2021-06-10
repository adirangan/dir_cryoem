      subroutine tesselation_get_nl_nm_ll(n_point,L_,tradius_min,nl_max
     $     ,nm_sum,ll_sum)
      implicit none
      include 'omp_lib.h'
      integer verbose
      data verbose / 0 /
      integer *4 n_point,npoint
      real *8 L_(0:3*n_point-1)
      integer *4, allocatable :: LT_(:) ! LT_(j) = j = index of point j within L_
      real *8 tradius_min
      real *8, allocatable :: v_n00_(:)
      real *8, allocatable :: v_0n0_(:)
      real *8, allocatable :: v_00n_(:)
      real *8, allocatable :: v_p00_(:)
      real *8, allocatable :: v_0p0_(:)
      real *8, allocatable :: v_00p_(:)
      integer *4 nl_nnn,nl_nnp,nl_npn,nl_npp
      integer *4 nl_pnn,nl_pnp,nl_ppn,nl_ppp
      integer *4 nl_(0:7),nl_max
      integer *4 nm_nnn,nm_nnp,nm_npn,nm_npp
      integer *4 nm_pnn,nm_pnp,nm_ppn,nm_ppp
      integer *4 nm_(0:7),nm_sum
      integer *4 ll_nnn,ll_nnp,ll_npn,ll_npp
      integer *4 ll_pnn,ll_pnp,ll_ppn,ll_ppp
      integer *4 ll_(0:7),ll_sum
      integer *4 nl,nm,ll
      integer *4 max_i4_f,sum_i4_f
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write (6,'(A)') ' [entering tesselation_get_nl_nm_ll]'
      end if !if (verbose.gt.0) then

      allocate(LT_(0:n_point-1)) ! LT_(j) = j = index of L_(j)
      do npoint=0,n_point-1
         LT_(npoint) = npoint
      enddo !do npoint=0,n_point-1

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A)') 'defining octants'
      end if ! if (verbose.gt.0) then
      allocate(v_n00_(0:2))
      allocate(v_0n0_(0:2))
      allocate(v_00n_(0:2))
      allocate(v_p00_(0:2))
      allocate(v_0p0_(0:2))
      allocate(v_00p_(0:2))
      call tesselation_define_octants(v_n00_,v_0n0_,v_00n_,v_p00_
     $     ,v_0p0_,v_00p_)
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'defining octants:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      timing_tic = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A)') 'calling tesselation_size'
      end if ! if (verbose.gt.0) then
      call tesselation_size(tradius_min,0,v_n00_,v_0n0_,v_00n_,.false.
     $     ,n_point,L_,LT_,nl_nnn,nm_nnn,ll_nnn)
      call tesselation_size(tradius_min,0,v_n00_,v_0n0_,v_00p_,.true.
     $     ,n_point,L_,LT_,nl_nnp,nm_nnp,ll_nnp)
      call tesselation_size(tradius_min,0,v_n00_,v_0p0_,v_00n_,.true.
     $     ,n_point,L_,LT_,nl_npn,nm_npn,ll_npn)
      call tesselation_size(tradius_min,0,v_n00_,v_0p0_,v_00p_,.false.
     $     ,n_point,L_,LT_,nl_npp,nm_npp,ll_npp)
      call tesselation_size(tradius_min,0,v_p00_,v_0n0_,v_00n_,.true.
     $     ,n_point,L_,LT_,nl_pnn,nm_pnn,ll_pnn)
      call tesselation_size(tradius_min,0,v_p00_,v_0n0_,v_00p_,.false.
     $     ,n_point,L_,LT_,nl_pnp,nm_pnp,ll_pnp)
      call tesselation_size(tradius_min,0,v_p00_,v_0p0_,v_00n_,.false.
     $     ,n_point,L_,LT_,nl_ppn,nm_ppn,ll_ppn)
      call tesselation_size(tradius_min,0,v_p00_,v_0p0_,v_00p_,.true.
     $     ,n_point,L_,LT_,nl_ppp,nm_ppp,ll_ppp)
      nl_(0) = nl_nnn
      nl_(1) = nl_nnp
      nl_(2) = nl_npn
      nl_(3) = nl_npp
      nl_(4) = nl_pnn
      nl_(5) = nl_pnp
      nl_(6) = nl_ppn
      nl_(7) = nl_ppp
      nl_max = max_i4_f(8,nl_)
      nm_(0) = nm_nnn
      nm_(1) = nm_nnp
      nm_(2) = nm_npn
      nm_(3) = nm_npp
      nm_(4) = nm_pnn
      nm_(5) = nm_pnp
      nm_(6) = nm_ppn
      nm_(7) = nm_ppp
      nm_sum = sum_i4_f(8,nm_)
      ll_(0) = ll_nnn
      ll_(1) = ll_nnp
      ll_(2) = ll_npn
      ll_(3) = ll_npp
      ll_(4) = ll_pnn
      ll_(5) = ll_pnp
      ll_(6) = ll_ppn
      ll_(7) = ll_ppp
      ll_sum = sum_i4_f(8,ll_)
      if (verbose.gt.0) then
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' nl_: '
     $        ,(nl_(nl),nl=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' nl_max: ',nl_max
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' nm_: '
     $        ,(nm_(nm),nm=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' nm_sum: ',nm_sum
         write(6,'(A,I0,A,8(I0,1X))') 'n_point: ',n_point,' ll_: '
     $        ,(ll_(ll),ll=0,7)
         write(6,'(A,I0,A,I0)') 'n_point: ',n_point,' ll_sum: ',ll_sum
      end if                    !verbose
      timing_toc = omp_get_wtime()
      if (verbose.gt.0) then
         write(6,'(A,A,F8.4)') 'calling tesselation_size:'
     $        ,' total_time ',timing_toc-timing_tic
      end if

      if (verbose.gt.0) then
         write (6,'(A)') ' [finished tesselation_get_nl_nm_ll]'
      end if !if (verbose.gt.0) then

      deallocate(v_n00_)
      deallocate(v_0n0_)
      deallocate(v_00n_)
      deallocate(v_p00_)
      deallocate(v_0p0_)
      deallocate(v_00p_)

      deallocate(LT_) ! LT_(j) = j = index of L_(j)

      end
