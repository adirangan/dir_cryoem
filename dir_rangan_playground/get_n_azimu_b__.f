!> Doxygen comment: ;\n
!> Defines n_azimu_b__, n_azimu_b__sum_ and n_azimu_b__sum_sum. ;\n
!> Doxygen comment: ;\n
      subroutine get_n_azimu_b__(
     $     n_k_p_max
     $     ,n_polar_a_upperbound
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,n_azimu_b__
     $     ,n_azimu_b__sum_
     $     ,n_azimu_b__sum_sum
     $     )
      implicit none
      integer n_k_p_max !integer: maximum number of k-values (in k-space polar coordinates). sometimes named ngridr. ;
      integer n_polar_a_upperbound !integer: upper bound on the number of k-values (in k-space polar coordinates). sometimes named ntmax. ;
      integer *4 n_polar_a_(0:0) !integer *4 array (size at least n_k_p_max): number of polar_a values on sphere of radius grid_k_p_r_(nk). also called nlats(nk). ;
      integer *4 quadrature_type_azimu_b !integer: quadrature type in the azimuthal direction. sometimes named itypep. 0=same number of points for all n_polar_a_ (oversampling poles). 1=adaptive (samples with roughly uniform density). ;
      integer *4 n_azimu_b__(0:0) !integer *4 arrary (size at least n_polar_a_upperbound*n_k_p_max): n_azimu_b__(np + nk*n_polar_a_upperbound) = number of azimu_b values on sphere of radius grid_k_p_r_(nk) and polar_a value polar_a_(np). ;
      integer *4 n_azimu_b__sum_(0:0) !integer *4 array (size at least n_k_p_max): n_azimu_b__sum_(nk) = sum of n_azimu_b__ over np. sometimes called numonsphere(nk). ;
      integer *4 n_azimu_b__sum_sum !integer: sum of azimu_b__sum_. total number of points on all spheres (out to radius grid_k_p_r_(n_k_p_max-1)). sometimes called ntot. ;
      real *8, allocatable :: grid_cos_polar_a_tmp_(:) !temporary: real *8 array (size at least n_k_p_max): values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called xnodesth(nk). ;
      real *8, allocatable :: grid_sin_polar_a_tmp_(:) !temporary: real *8 array (size at least n_k_p_max): values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called sthetas(nk). ;
      real *8, allocatable :: weight_cos_polar_a_tmp_(:) !temporary real *8 array (size at least n_k_p_max): weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called wtsth_(np). ;
      integer *4, allocatable :: n_azimu_b_tmp_(:) !temporary integer *4 array (size at least n_k_p_max): number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_r_(nk). also called ngridps(np). ;
      real *8, allocatable :: azimu_b_step_tmp_(:) !temporary real *8 array (size at least n_k_p_max): grid-spacing for azimu_b. ;
      integer *4 n_azimu_b__sum_tmp !temporary integer *4: total number of points on sphere at a particular radius grid_k_p_r_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
      integer max_i4_f
      integer n_polar_a_max,nk
      logical flag_memcheck
c$$$      %%%%%%%%
      n_polar_a_max = max_i4_f(n_k_p_max,n_polar_a_)
      allocate(grid_cos_polar_a_tmp_(0:1+n_polar_a_max-1))
      call cs1_r8(n_polar_a_max,grid_cos_polar_a_tmp_)
      allocate(grid_sin_polar_a_tmp_(0:1+n_polar_a_max-1))
      call cs1_r8(n_polar_a_max,grid_sin_polar_a_tmp_)
      allocate(weight_cos_polar_a_tmp_(0:1+n_polar_a_max-1))
      call cs1_r8(n_polar_a_max,weight_cos_polar_a_tmp_)
      allocate(n_azimu_b_tmp_(0:1+n_polar_a_max-1))
      call cs1_i4(n_polar_a_max,n_azimu_b_tmp_)
      allocate(azimu_b_step_tmp_(0:1+n_polar_a_max-1))
      call cs1_r8(n_polar_a_max,azimu_b_step_tmp_)
c$$$      %%%%%%%%
      call cl1_i4(n_polar_a_max*n_k_p_max,n_azimu_b__)
      call cl1_i4(n_k_p_max,n_azimu_b__sum_)
c$$$      %%%%%%%%
      n_azimu_b__sum_sum = 0
      do nk=0,n_k_p_max-1
         call getspheregrid(
     $     n_polar_a_(nk) !integer *4: number of polar_a values on sphere of radius grid_k_p_r_(nk). ;
     $     ,quadrature_type_azimu_b !integer: quadrature type in the azimuthal direction. sometimes named itypep. 0=same number of points for all n_polar_a_ (oversampling poles). 1=adaptive (samples with roughly uniform density). ;
     $     ,grid_cos_polar_a_tmp_ !temporary: real *8 array (size at least n_k_p_max): values for dcos(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called xnodesth(nk). ;
     $     ,grid_sin_polar_a_tmp_ !temporary: real *8 array (size at least n_k_p_max): values for dsin(polar_a) associated with n_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called sthetas(nk). ;
     $     ,weight_cos_polar_a_tmp_ !temporary real *8 array (size at least n_k_p_max): weight associated with grid_cos_polar_a_(nk) in k-space polar coordinates for a particular grid_k_p_r_(nk). sometimes called wtsth_(np). ;
     $     ,n_azimu_b_tmp_ !temporary integer *4 array (size at least n_k_p_max): number of nodes in azimu_b for each polar_a in k-space polar coordinates for a particular grid_k_p_r_(nk). also called ngridps(np). ;
     $     ,azimu_b_step_tmp_ !temporary real *8 array (size at least n_k_p_max): grid-spacing for azimu_b. ;
     $     ,n_azimu_b__sum_tmp !temporary integer *4: total number of points on sphere at a particular radius grid_k_p_r_(nk) in k-space polar coordinates. sometimes called nspherebig_tmp. ;
     $     )
         call cp1_i4(n_polar_a_(nk),n_azimu_b_tmp_
     $        ,n_azimu_b__(0+ nk*n_polar_a_upperbound))
         n_azimu_b__sum_(nk) = n_azimu_b__sum_tmp
         n_azimu_b__sum_sum = n_azimu_b__sum_sum + n_azimu_b__sum_tmp
      enddo !do nk=0,n_k_p_max-1
      flag_memcheck = .true.
      call cxs_r8(n_polar_a_max,grid_cos_polar_a_tmp_
     $     ,'grid_cos_polar_a_tmp_',flag_memcheck)
      call cxs_r8(n_polar_a_max,grid_sin_polar_a_tmp_
     $     ,'grid_sin_polar_a_tmp_',flag_memcheck)
      call cxs_r8(n_polar_a_max,weight_cos_polar_a_tmp_
     $     ,'weight_cos_polar_a_tmp_',flag_memcheck)
      call cxs_i4(n_polar_a_max,n_azimu_b_tmp_,'n_azimu_b_tmp_'
     $     ,flag_memcheck)
      call cxs_r8(n_polar_a_max,azimu_b_step_tmp_,'azimu_b_step_tmp_'
     $     ,flag_memcheck)
      if (flag_memcheck.eqv..false.) then
         write(6,'(A)') 'Warning, memcheck failed in get_n_azimu_b__'
         stop !due to error. ;
      end if !if (flag_memcheck.eqv..false.) then
      end; !subroutine
