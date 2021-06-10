!> Doxygen comment: ;\n
!>    Builds new model Y_ (i.e., in spherical harmonic basis)  ;\n
!>    from images M_k_p_ (in k-space polar coordinates), ;\n
!>    with image parameters given by alpha2d__ array. ; ;\n
!>--------------------------------------------------------------------- ;\n
!>    INPUT: ;\n
!>     ;\n
!>    n_M            integer: number of images ; ;\n
!>    I_M_sample_    integer: array indexing the relevant nM entries of M_k_p__. ; ;\n
!>    ld_M           integer: leading dimension of images in k-space polar coordinates. ; ;\n
!>    M_k_p__        complex *16: stack of images in k-space polar coordinates (size ld_M*n_M). ; ;\n
!>    n_CTF          integer: number of ctf-functions. ; ;\n
!>    ld_CTF         integer: leading dimension of CTF-array (usually ld_M). ; ;\n
!>    CTF_k_p_       complex *16: stack of ctf-functions in k-space polar coordinates (size ld_CTF*n_CTF). ; ;\n
!>    n_alpha        integer: size of alpha2d__. ; ;\n
!>    alpha2d__(n_alpha,*) real *8: 2d array of image parameters (see excerpt_define_nalpha.f). ; ;\n
!>    n_w_csum_()     integer: indexing array for points on successive circles within each image. ;  ;\n
!>                   n_w_csum_(nk) represents the starting index of points (in angle-w) on ring at radius grid_k_p_r_(nk). ; ;\n
!>                   Note that n_w_csum_ is base 1 (instead of base 0) for compatibility with get_template_size. ; ;\n
!>    n_polar_a_()   integer:  number of quadrature nodes in polar_a on sphere defined by index n_k_p_r_cur. ; ;\n
!>    quadrature_type_azimu_b integer: flag determining quadrature scheme in azimu_b direction. ; ;\n
!>                      0 = same on all latitudes (which oversamples poles). ; ;\n
!>                      1 = adaptive (appropriately samples toward poles). ; ;\n
!>                   if quadrature_type_azimu_b=0, then n_azimu_b = nint(n_polar_a_(nk)*phi_over). ; ;\n
!>                   if quadrature_type_azimu_b=1, then n_azimu_b = nint(n_polar_a_(nk)*phi_over*dsin(polar_a)). ; ;\n
!>                   Note: phi_over is set in getgridph (typically = 2). ; ;\n
!>    grid_k_p_r_()    real *8: radius associated with successive circles in templates in k-space polar coordinates. ; ;\n
!>    n_w_()         integer: number of output points on successive circles in templates in k-space polar coordinates. ; ;\n
!>    n_k_p_r_low        integer: index of lowest frequency sphere under consideration. Note that this runs from 1 to n_k_p_r_max. ; ;\n
!>    n_k_p_r_cur        integer: index of highest frequency sphere under consideration. Note that this runs from 1 to n_k_p_r_max. ; ;\n
!>    n_Y_lm_csum_()  integer: array of length n_k_p_r_max indicating where the Y_ for that shell begins. ; ;\n
!>    n_Y_l_()       integer: array of length n_k_p_r_max defining the orders of the various Y_ on successive shells. ; ;\n
!>    lsq_oversample integer: lsq_oversampling parameter for least-squares solver. ; ;\n
!>    lsq_interpolation_order integer: interpolation order for least-squares solver. ; ;\n
!>    lsq_eps        real *8: tolerance epsilon passed into least-squares solver ;\n
!>     ;\n
!>    OUTPUT:  ;\n
!>     ;\n
!>    Y_(:)         complex *16: solution to least-squares problem expressed in spherical harmonics in successive shells. ; ;\n
!>     ;\n
!> Doxygen comment: ;\n
C***********************************************************************
      subroutine rebuild_model_2(
     $     n_M
     $     ,I_M_sample_
     $     ,ld_M
     $     ,M_k_p__
     $     ,n_CTF
     $     ,ld_CTF
     $     ,CTF_k_p_
     $     ,alpha2d__
     $     ,n_w_csum_
     $     ,n_polar_a_
     $     ,quadrature_type_azimu_b
     $     ,grid_k_p_r_
     $     ,n_w_
     $     ,n_k_p_r_low
     $     ,n_k_p_r_cur
     $     ,n_Y_lm_csum_
     $     ,n_Y_l_
     $     ,lsq_oversample
     $     ,lsq_interpolation_order
     $     ,lsq_eps
     $     ,Y_
     $     )
C***********************************************************************
C     Builds new model Y_ (i.e., in spherical harmonic basis) 
C     from images M_k_p_ (in k-space polar coordinates),
C     with image parameters given by alpha2d__ array. ;
C---------------------------------------------------------------------
C     INPUT:
C
C     n_M            integer: number of images ;
C     I_M_sample_    integer: array indexing the relevant nM entries of M_k_p__. ;
C     ld_M           integer: leading dimension of images in k-space polar coordinates. ;
C     M_k_p__        complex *16: stack of images in k-space polar coordinates (size ld_M*n_M). ;
C     n_CTF          integer: number of ctf-functions. ;
C     ld_CTF         integer: leading dimension of CTF-array (usually ld_M). ;
C     CTF_k_p_       complex *16: stack of ctf-functions in k-space polar coordinates (size ld_CTF*n_CTF). ;
c     n_alpha        integer: size of alpha2d__. ;
c     alpha2d__(n_alpha,*) real *8: 2d array of image parameters (see excerpt_define_nalpha.f). ;
c     n_w_csum_()     integer: indexing array for points on successive circles within each image. ; 
c                    n_w_csum_(nk) represents the starting index of points (in angle-w) on ring at radius grid_k_p_r_(nk). ;
c                    Note that n_w_csum_ is base 1 (instead of base 0) for compatibility with get_template_size. ;
c     n_polar_a_()   integer:  number of quadrature nodes in polar_a on sphere defined by index n_k_p_r_cur. ;
c     quadrature_type_azimu_b integer: flag determining quadrature scheme in azimu_b direction. ;
c                       0 = same on all latitudes (which oversamples poles). ;
c                       1 = adaptive (appropriately samples toward poles). ;
c                    if quadrature_type_azimu_b=0, then n_azimu_b = nint(n_polar_a_(nk)*phi_over). ;
c                    if quadrature_type_azimu_b=1, then n_azimu_b = nint(n_polar_a_(nk)*phi_over*dsin(polar_a)). ;
c                    Note: phi_over is set in getgridph (typically = 2). ;
c     grid_k_p_r_()    real *8: radius associated with successive circles in templates in k-space polar coordinates. ;
c     n_w_()         integer: number of output points on successive circles in templates in k-space polar coordinates. ;
c     n_k_p_r_low        integer: index of lowest frequency sphere under consideration. Note that this runs from 1 to n_k_p_r_max. ;
c     n_k_p_r_cur        integer: index of highest frequency sphere under consideration. Note that this runs from 1 to n_k_p_r_max. ;
c     n_Y_lm_csum_()  integer: array of length n_k_p_r_max indicating where the Y_ for that shell begins. ;
c     n_Y_l_()       integer: array of length n_k_p_r_max defining the orders of the various Y_ on successive shells. ;
c     lsq_oversample integer: lsq_oversampling parameter for least-squares solver. ;
c     lsq_interpolation_order integer: interpolation order for least-squares solver. ;
c     lsq_eps        real *8: tolerance epsilon passed into least-squares solver
C
C     OUTPUT: 
c
C     Y_(:)         complex *16: solution to least-squares problem expressed in spherical harmonics in successive shells. ;
c
C***********************************************************************
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF
      integer n_w_csum_(0:0),n_polar_a_(0:0)
      integer quadrature_type_azimu_b
      include '/usr/include/fftw3.f'
      include 'excerpt_define_nalpha.f'
      integer n_w_(0:0)
      integer n_k_p_r_low,n_k_p_r_cur
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_0in_(:)
      complex *16, allocatable :: fftw_out_(:)
      integer n_Y_lm_csum_(0:0),n_Y_l_(0:0)
      integer lsq_interpolation_order
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      real *8 grid_k_p_r_(0:n_k_p_r_cur-1)
      real *8 lsq_oversample,lsq_eps
      integer *4 lsq_n_quad
      complex *16 M_k_p__(0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      complex *16 Y_(0:0)
      integer nk,n_A,na,n_w,nw
      integer lsq_n_iteration,lsq_niteration
      integer n_w_M,n_w_M_max
      integer n_Y_lm_max,n_Y_lm,n_Y_l,n_Y_lm_sum
      integer n_polar_a_oversample
      real *8, allocatable :: polar_a_(:)
      real *8, allocatable :: azimu_b_(:)
      real *8, allocatable :: alpha2d_adi__(:,:)
      real *8, allocatable :: alpha2d_marina__(:,:)
      logical flag_adi_vs_marina
      parameter (flag_adi_vs_marina=.false.)
      complex *16, allocatable :: Y_cur_(:)
      complex *16, allocatable :: M_k_c_(:)
      complex *16, allocatable :: weight_CTF_k_c_(:)
      external multaha
      real *8 fnod_r8_f !function output
      logical flag_memory_checkset !temporary flag for checking memory setting. ;
      character(len=1024) format_string
      character(len=1024) tmp_dir
      character(len=1024) tmp_str
      integer *4 tmp_tab

c$$$         %%%%%%%%
      if (verbose.gt.0) then
         write(6,'(A)') '[entering rebuild_model_2]'
      end if !if (verbose.gt.0) then
      n_Y_lm_sum = n_Y_lm_csum_(n_k_p_r_cur-1) + (n_Y_l_(n_k_p_r_cur-1)
     $     +1)**2
      if (verbose.gt.1) then
         write(6,*) 'ld_M: ',ld_M
         write(6,*) 'n_M: ',n_M
         write(6,*) 'n_w_csum_: ',(n_w_csum_(nk),nk=0,n_k_p_r_cur-1)
         write(6,*) 'n_polar_a_: ',(n_polar_a_(nk),nk=0,n_k_p_r_cur-1)
         write(6,*) 'quadrature_type_azimu_b: ',quadrature_type_azimu_b
         write(6,*) 'grid_k_p_r_: ',(grid_k_p_r_(nk),nk=0,n_k_p_r_cur-1)
         write(6,*) 'n_w_: ',(n_w_(nk),nk=0,n_k_p_r_cur-1)
         write(6,*) 'n_k_p_r_low: ',n_k_p_r_low
         write(6,*) 'n_k_p_r_cur: ',n_k_p_r_cur
         write(6,*) 'n_Y_lm_csum_: ',(n_Y_lm_csum_(nk),nk=0,n_k_p_r_cur
     $        -1)
         write(6,*) 'n_Y_l_: ',(n_Y_l_(nk),nk=0,n_k_p_r_cur-1)
         write(6,*) 'lsq_oversample: ',lsq_oversample
         write(6,*) 'lsq_interpolation_order: ',lsq_interpolation_order
      end if
c$$$         %%%%%%%%

c$$$         %%%%%%%%
      n_A = 0
      do nk = 0,n_k_p_r_cur-1
         n_A = n_A + n_w_(nk)
      end do !do nk = 0,n_k_p_r_cur-1
      if (n_A.gt.ld_M) then
         write(6,*) 'Warning, ld_M: ',ld_M
     $        ,' .neq. n_A: ',n_A,' in rebuild_model_2'
      end if
      if (verbose.gt.1) then
         write(6,'(2(A,I0))') ' n_k_p_r_cur ' , n_k_p_r_cur , ' n_A ' ,
     $        n_A
      end if !if (verbose.gt.1) then
c$$$         %%%%%%%%
      allocate(fftw_plan_frwd_(0:n_k_p_r_cur-1))
      allocate(fftw_plan_back_(0:n_k_p_r_cur-1))
      allocate(fftw_0in_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nk=0,n_k_p_r_cur-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nk),n_w_(nk)
     $        ,fftw_0in_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nk),n_w_(nk)
     $        ,fftw_out_(na),fftw_0in_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nk)
      enddo !do nk=0,n_k_p_r_cur-1
c$$$         %%%%%%%%
      n_Y_l = n_Y_l_(n_k_p_r_cur-1)
      n_Y_lm_max = (n_Y_l+1)**2
      if (verbose.gt.1) then
         write(6,*) ' n_Y_l = ',n_Y_l
         write(6,*) ' n_Y_lm_max = ',n_Y_lm_max
      end if
      allocate(Y_cur_(0:1+n_Y_lm_max))
      call cs1_c16(n_Y_lm_max,Y_cur_)
c$$$         %%%%%%%%
      n_polar_a_oversample = nint(lsq_oversample*n_polar_a_(n_k_p_r_cur
     $     -1))
      if (verbose.gt.1) then
         write(6,*) ' lsq_oversample = ',lsq_oversample
         write(6,*) ' n_polar_a_oversample = ',n_polar_a_oversample
      end if
c$$$         %%%%%%%%
      n_w_M_max = n_M*n_w_(n_k_p_r_cur-1)
      if (verbose.gt.1) then
         write(6,*) ' n_w_M_max = ',n_w_M_max
      end if
      allocate(M_k_c_(0:1+n_w_M_max-1))
      call cs1_c16(n_w_M_max,M_k_c_)
      allocate(polar_a_(0:1+n_w_M_max-1))
      call cs1_r8(n_w_M_max,polar_a_)
      allocate(azimu_b_(0:1+n_w_M_max-1))
      call cs1_r8(n_w_M_max,azimu_b_)
      allocate(weight_CTF_k_c_(0:1+n_w_M_max-1))
      call cs1_c16(n_w_M_max,weight_CTF_k_c_)
      allocate(alpha2d_adi__(0:n_alpha-1,0:1+n_M-1))
      call cs1_r8(n_alpha*n_M,alpha2d_adi__)
      allocate(alpha2d_marina__(0:n_alpha-1,0:1+n_M-1))
      call cs1_r8(n_alpha*n_M,alpha2d_marina__)
      if (verbose.gt.2) then
      call print_all_r8(n_alpha,alpha2d_marina__(0,0)
     $     ,' alpha2d_marina__: ')
      call print_all_r8(n_alpha,alpha2d_marina__(0,1)
     $     ,' alpha2d_marina__: ')
      call print_all_r8(n_alpha,alpha2d_marina__(0,2)
     $     ,' alpha2d_marina__: ')
      end if !if (verbose.gt.2) then
      call convert_alpha2d_adi_to_marina_0(
     $     n_M
     $     ,alpha2d__
     $     ,alpha2d_marina__
     $     )
      call convert_alpha2d_marina_to_adi_0(
     $     n_M
     $     ,alpha2d_marina__
     $     ,alpha2d_adi__
     $     )
      if (verbose.gt.1) then
      write(6,'(A,F24.16)') ' fnorm(alpha2d__ - alpha2d_adi__): ' ,
     $     fnod_r8_f(n_alpha*n_M,alpha2d__,alpha2d_adi__)
      end if !if (verbose.gt.1) then
c$$$         %%%%%%%%

c$$$         %%%%%%%%
      na = 0
      do nk = 0,n_k_p_r_cur-1
c$$$         %%%%%%%%
         if (verbose.gt.1) then
            write(6,'(2(A,I0))') ' nk ' , nk , ' n_k_p_r_cur ' ,
     $           n_k_p_r_cur
         end if !if (verbose.gt.1) then
         if (nk.ge.0) then
            if (flag_adi_vs_marina.eqv..true.) then
            call get_lsqdata_2(
     $           n_M
     $           ,I_M_sample_
     $           ,ld_M
     $           ,M_k_p__
     $           ,n_CTF
     $           ,ld_CTF
     $           ,CTF_k_p_
     $           ,grid_k_p_r_(nk)
     $           ,fftw_plan_frwd_(nk)
     $           ,fftw_plan_back_(nk)
     $           ,n_w_csum_(nk)
     $           ,n_w_(nk)
     $           ,fftw_0in_(na)
     $           ,fftw_out_(na)
     $           ,alpha2d_adi__
     $           ,M_k_c_
     $           ,polar_a_
     $           ,azimu_b_
     $           ,weight_CTF_k_c_
     $           ,n_w_M
     $           )
            end if !if (flag_adi_vs_marina.eqv..true.) then
            if (flag_adi_vs_marina.eqv..false.) then
            call get_lsqdata_alamarina(
     $           n_M
     $           ,I_M_sample_
     $           ,ld_M
     $           ,M_k_p__
     $           ,n_CTF
     $           ,ld_CTF
     $           ,CTF_k_p_
     $           ,grid_k_p_r_(nk)
     $           ,fftw_plan_frwd_(nk)
     $           ,fftw_plan_back_(nk)
     $           ,n_w_csum_(nk)
     $           ,n_w_(nk)
     $           ,fftw_0in_(na)
     $           ,fftw_out_(na)
     $           ,alpha2d_marina__
     $           ,M_k_c_
     $           ,polar_a_
     $           ,azimu_b_
     $           ,weight_CTF_k_c_
     $           ,n_w_M
     $           )
            end if !if (flag_adi_vs_marina.eqv..false.) then
c$$$         %%%%%%%%
            n_Y_l = n_Y_l_(nk)
            n_Y_lm = (n_Y_l+1)**2
            lsq_n_quad = nint(lsq_oversample*n_polar_a_(nk))
            lsq_n_iteration = 1000
            lsq_niteration = 0
            if (verbose.gt.1) then
               write(6,'(7(A,I0))')
     $              ' nk: ' , nk 
     $              , ' n_w_M: ' , n_w_M
     $              , ' n_Y_l: ' , n_Y_l
     $              , ' n_Y_lm: ' , n_Y_lm
     $              , ' lsq_n_quad: ' , lsq_n_quad
     $              , ' lsq_interpolation_order: ' 
     $              , lsq_interpolation_order
     $              , ' lsq_n_iteration: ' , lsq_n_iteration
               call print_sub_c16(n_w_M,M_k_c_,' M_k_c_: ')
               call print_sub_c16(n_w_M,weight_CTF_k_c_
     $              ,' C_k_c_: ')
               call print_sub_r8(n_w_M,polar_a_,' polar_a_: ')
               call print_sub_r8(n_w_M,azimu_b_,' azimu_b_: ')
            end if !if (verbose.gt.1) then
c$$$c$$$            %%%%%%%%
c$$$         if (nk.eq.8) then
c$$$            include 'tmp_lsq_write.f'
c$$$         end if !if (nk.eq.8) then
c$$$c$$$            %%%%%%%%
            call lsqctfsolve(
     $           M_k_c_
     $           ,polar_a_
     $           ,azimu_b_
     $           ,n_w_M
     $           ,Y_cur_
     $           ,n_Y_l
     $           ,lsq_n_quad
     $           ,lsq_interpolation_order
     $           ,weight_CTF_k_c_
     $           ,multaha
     $           ,lsq_eps
     $           ,lsq_n_iteration
     $           ,lsq_niteration
     $           )
            call cp1_c16(n_Y_lm,Y_cur_,Y_(n_Y_lm_csum_(nk)))
         end if !if (nk.ge.0) then
         if (verbose.gt.1) then
            call print_sub_c16(n_Y_lm,Y_cur_,' Y_cur_: ')
         end if !if (verbose.gt.1) then
         na = na + n_w_(nk)
c$$$         %%%%%%%%
      enddo !do nk = n_k_p_r_low-1,n_k_p_r_cur-1
c$$$         %%%%%%%%
      if (verbose.gt.2) then
         call print_sub_c16(n_Y_lm_sum,Y_,' Y_: ')         
         call print_mid_c16(n_Y_lm_sum,Y_,' Y_: ')         
      end if !if (verbose.gt.-2) then

c$$$         %%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if !if (verbose.gt.1) then
      do nk=0,n_k_p_r_cur-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nk))
         call dfftw_destroy_plan(fftw_plan_back_(nk))
      enddo !do nk=0,n_k_p_r_cur-1
      if (verbose.gt.1) then
         write(6,'(A)') ' deallocating. '
      end if
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_0in_)
      deallocate(fftw_out_)
c$$$         %%%%%%%%

c$$$         %%%%%%%%
      if (verbose.gt.1) then
         write(6,'(A)') ' checking memory allocation. '
      end if !if (verbose.gt.1) then
      flag_memory_checkset = .true.
      call cxs_c16(n_Y_lm_max,Y_cur_,'Y_cur_',flag_memory_checkset)
      call cxs_c16(n_w_M_max,M_k_c_,'M_k_c_',flag_memory_checkset)
      call cxs_r8(n_w_M_max,polar_a_,'polar_a_',flag_memory_checkset)
      call cxs_r8(n_w_M_max,azimu_b_,'azimu_b_',flag_memory_checkset)
      call cxs_c16(n_w_M_max,weight_CTF_k_c_,'weight_CTF_k_c_'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_M,alpha2d_adi__,'alpha2d_adi__'
     $     ,flag_memory_checkset)
      call cxs_r8(n_alpha*n_M,alpha2d_marina__,'alpha2d_marina__'
     $     ,flag_memory_checkset)
      if (flag_memory_checkset.eqv..true.) then
         if (verbose.gt.1) then 
            write(6,'(A)') '[checkset passed]'
         end if !if (verbose.gt.1) then 
      end if !if (flag_memory_checkset.eqv..true.) then
      if (flag_memory_checkset.eqv..false.) then
         write(6,'(A)') '[checkset failed] <-- WARNING'
         stop !exit program due to error.
      end if !if (flag_memory_checkset.eqv..false.) then
c$$$         %%%%%%%%

      if (verbose.gt.0) then
         write(6,'(A)') '[finished rebuild_model_2]'
      end if !if (verbose.gt.0) then
      end !subroutine

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
