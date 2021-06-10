!> Doxygen comment: ;\n
!>    Calculates residuals and estimates cross-correlation. ;\n
!>    Uses model Y_ (in spherical harmonic basis)  ;\n
!>    as well as images M_k_p_ (in k-space polar coordinates), ;\n
!>    with image parameters given by alpha2d__ array. ; ;\n
!> -------------------------------------------------------------------- ;\n
!>    INPUT: ;\n
!>  ;\n
!>    rseed          integer *4: random seed ; ;\n
!>    n_M            integer: number of images ; ;\n
!>    I_M_sample_    integer: indexing array used to reference n_M images in M_k_p__ ;\n
!>    ld_M           integer: leading dimension of images in k-space polar coordinates. ; ;\n
!>    M_k_p__        complex *16: stack of images in k-space polar coordinates. ; ;\n
!>    n_CTF          integer: number of ctf-functions. ; ;\n
!>    ld_CTF         integer: leading dimension of CTF-array (usually ld_M). ; ;\n
!>    CTF_k_p_       complex *16: stack of ctf-functions in k-space polar coordinates (size ld_CTF*n_CTF). ; ;\n
!>    n_alpha        integer: size of alpha2d__. ; ;\n
!>    alpha2d__(3,*) real *8: 2d array of image parameters (see nalpha_define.f). ; ;\n
!>    n_w_csum_()     integer: indexing array for points on successive circles within each image. ;  ;\n
!>                   n_w_csum_(nk) represents the number of points (in angle-w) on ring at radius grid_k_p_(nk). ; ;\n
!>    n_polar_a_()   integer:  number of quadrature nodes in polar_a on sphere defined by index n_k_cur. ; ;\n
!>    quadrature_type_azimu_b integer: flag determining quadrature scheme in azimu_b direction. ; ;\n
!>                      0 = same on all latitudes (which oversamples poles). ; ;\n
!>                      1 = adaptive (appropriately samples toward poles). ; ;\n
!>                   if quadrature_type_azimu_b=0, then n_azimu_b = nint(n_polar_a_(nk)*phi_over). ; ;\n
!>                   if quadrature_type_azimu_b=1, then n_azimu_b = nint(n_polar_a_(nk)*phi_over*dsin(polar_a)). ; ;\n
!>                   Note: phi_over is set in getgridph (typically = 2). ; ;\n
!>    grid_k_p_()    real *8: radius associated with successive circles in templates in k-space polar coordinates. ; ;\n
!>    n_w_()         integer: number of output points on successive circles in templates in k-space polar coordinates. ; ;\n
!>    n_k_low        integer: index of lowest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ; ;\n
!>    n_k_cur        integer: index of highest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ; ;\n
!>    n_Y_lm_csum_()  integer: array of length n_k_p_max indicating where the Y_ for that shell begins. ; ;\n
!>    n_Y_l_()       integer: array of length n_k_p_max defining the orders of the various Y_ on successive shells. ; ;\n
!>    lsq_oversample integer: lsq_oversampling parameter for least-squares solver. ; ;\n
!>    lsq_interpolation_order integer: interpolation order for least-squares solver. ; ;\n
!>    Y_(:)          complex *16: solution to least squares problem expressed in spherical harmonics in successive shells. ; ;\n
!>    n_residual_loading    integer: number of M_loading_ vectors per image. ; ;\n
!>    n_residual_iteration  integer: number of iterations to perform. ; ;\n
!>  ;\n
!>    OUTPUT:  ;\n
!>  ;\n
!>    M_residual__(:)  complex *16: M_residual__ = M_transform_ - Y_slice_. ; ;\n
!>    M_loading_(:)    complex *16: M_loading_ vectors. ; ;\n
!>  ;\n

C***********************************************************************
      subroutine test_residual_4(rseed,n_M,I_M_sample_,ld_M,M_k_p__
     $     ,n_CTF,ld_CTF,CTF_k_p_,alpha2d__,n_w_csum_,n_polar_a_
     $     ,quadrature_type_azimu_b ,grid_k_p_ ,n_w_,n_k_low,n_k_cur
     $     ,n_Y_lm_csum_cur,n_Y_lm_csum_ ,n_Y_l_ ,lsq_oversample
     $     ,lsq_interpolation_order ,Y_ ,M_residual__,n_residual_loading
     $     ,n_residual_iteration,M_loading_)
C***********************************************************************
C     Calculates residuals and estimates cross-correlation.
C     Uses model Y_ (in spherical harmonic basis) 
C     as well as images M_k_p_ (in k-space polar coordinates),
C     with image parameters given by alpha2d__ array. ;
C---------------------------------------------------------------------
C     INPUT:
C
C     rseed          integer *4: random seed ;
C     n_M            integer: number of images ;
C     I_M_sample_    integer: indexing array used to reference n_M images in M_k_p__
C     ld_M           integer: leading dimension of images in k-space polar coordinates. ;
C     M_k_p__        complex *16: stack of images in k-space polar coordinates. ;
C     n_CTF          integer: number of ctf-functions. ;
C     ld_CTF         integer: leading dimension of CTF-array (usually ld_M). ;
C     CTF_k_p_       complex *16: stack of ctf-functions in k-space polar coordinates (size ld_CTF*n_CTF). ;
c     n_alpha        integer: size of alpha2d__. ;
c     alpha2d__(3,*) real *8: 2d array of image parameters (see nalpha_define.f). ;
c     n_w_csum_()     integer: indexing array for points on successive circles within each image. ; 
c                    n_w_csum_(nk) represents the number of points (in angle-w) on ring at radius grid_k_p_(nk). ;
c     n_polar_a_()   integer:  number of quadrature nodes in polar_a on sphere defined by index n_k_cur. ;
c     quadrature_type_azimu_b integer: flag determining quadrature scheme in azimu_b direction. ;
c                       0 = same on all latitudes (which oversamples poles). ;
c                       1 = adaptive (appropriately samples toward poles). ;
c                    if quadrature_type_azimu_b=0, then n_azimu_b = nint(n_polar_a_(nk)*phi_over). ;
c                    if quadrature_type_azimu_b=1, then n_azimu_b = nint(n_polar_a_(nk)*phi_over*dsin(polar_a)). ;
c                    Note: phi_over is set in getgridph (typically = 2). ;
c     grid_k_p_()    real *8: radius associated with successive circles in templates in k-space polar coordinates. ;
c     n_w_()         integer: number of output points on successive circles in templates in k-space polar coordinates. ;
c     n_k_low        integer: index of lowest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ;
c     n_k_cur        integer: index of highest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ;
c     n_Y_lm_csum_()  integer: array of length n_k_p_max indicating where the Y_ for that shell begins. ;
c     n_Y_l_()       integer: array of length n_k_p_max defining the orders of the various Y_ on successive shells. ;
c     lsq_oversample integer: lsq_oversampling parameter for least-squares solver. ;
c     lsq_interpolation_order integer: interpolation order for least-squares solver. ;
C     Y_(:)          complex *16: solution to least squares problem expressed in spherical harmonics in successive shells. ;
c     n_residual_loading    integer: number of M_loading_ vectors per image. ;
c     n_residual_iteration  integer: number of iterations to perform. ;
cC
C     OUTPUT: 
c
c     M_residual__(:)  complex *16: M_residual__ = M_transform_ - Y_slice_. ;
c     M_loading_(:)    complex *16: M_loading_ vectors. ;
c
C***********************************************************************
      implicit none
      integer verbose
      data verbose / 4 /
      integer *4 rseed
      real *8 adi_rand_f
      include 'adi_rand_f_define.f'
      integer n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF
      complex *16 M_k_p__(0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      include '/usr/include/fftw3.f'
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_0in_(:)
      complex *16, allocatable :: fftw_out_(:)
      include 'nalpha_define.f'
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      integer n_w_csum_(0:n_k_cur-1)
      integer n_polar_a_(0:n_k_cur-1)
      integer n_w_(0:n_k_cur-1)
      integer quadrature_type_azimu_b
      real *8 grid_k_p_(0:n_k_cur-1)
      integer n_k_low,n_k_cur
      integer n_Y_lm_csum_cur
      integer n_Y_lm_csum_(0:n_k_cur-1)
      integer n_Y_l_(0:n_k_cur-1)
      complex *16 Y_(0:n_Y_lm_csum_cur-1)
      real *8 lsq_oversample
      integer lsq_interpolation_order
      complex *16  M_residual__(0:0)
      integer n_residual_loading,n_residual_iteration
      complex *16 M_loading_(0:n_residual_loading*n_M-1)
      integer nM,nk,nA,n_A,nw,ntmp,nctf
      integer n_w_M_max,n_polar_a_oversample
      integer n_Y_lm_max,n_Y_l,n_polar_a_oversample_max
      real *8, allocatable :: polar_a__(:,:)
      real *8, allocatable :: azimu_b__(:,:)
      complex *16, allocatable :: M_transform_(:)
      complex *16, allocatable :: Y_slice_(:)
      complex *16, allocatable :: weight_CTF_k_c__(:,:)
      integer *4 , allocatable :: n_w_M_(:)
      integer nM_loading,nY_lm_sum
c$$$      Below the array G stores a collection of vectors G_{1}, G_{2}, etc, ;
c$$$      which are constructed according to: ;
c$$$      G^{j} = \sum_{nM} Et_{nM}*R_{nM}*v_{nM}^{j}, ;
c$$$      where: ;
c$$$      nM sums over all the images. ;
c$$$      Et_{nM} refers to the adjoint of the evaluation operator En_{nM}, ;
c$$$      where En_{nM} serves to evaluate a spherical harmonic expansion ;
c$$$      on the points (in k-space cartesian coordinates) associated with '
c$$$      image nM (which depend on the image parameters alpha_{nM}). ;
c$$$      R_{nM} refers to the residual associated with image nM. ;
c$$$      v_{nM}^{j} refers to the nM-element of the j^{th} 'loading' vector. ;
c$$$      The vector H is a temporary variable storing Et_{nM}*R_{nM}. ;
      complex *16, allocatable :: G_(:) !temporary: stores array of \sum Et*R*v^{j} (one for each loading, each in terms of spherical harmonic coefficients). ;
c$$$      complex *16, allocatable :: H_(:) !temporary: stores single Et*R (in terms of spherical harmonic coefficients). ;
      complex *16, allocatable :: H__(:) !temporary: stores all Et*R (in terms of spherical harmonic coefficients) across all images. ;
      complex *16 HH,HG,RR
      integer n_iteration,niteration
      external multaha
      real *8 al2_c16_f
      complex *16 dot_c16_f
      character(len=1024) format_string
      if (verbose.gt.0) then
         write(6,*) '[entering test_residual_4]'
      end if
      if (verbose.gt.1) then
         write(6,*) 'ld_M: ',ld_M
         write(6,*) 'n_M: ',n_M
         write(6,*) 'n_w_csum_: ',(n_w_csum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_CTF: ',n_CTF
         write(6,*) 'ld_CTF: ',ld_CTF
         write(6,*) 'n_alpha: ',n_alpha
         write(6,*) 'n_polar_a_: ',(n_polar_a_(nk),nk=0,n_k_cur-1)
         write(6,*) 'quadrature_type_azimu_b: ',quadrature_type_azimu_b
         write(6,*) 'grid_k_p_: ',(grid_k_p_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_w_: ',(n_w_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_k_low: ',n_k_low
         write(6,*) 'n_k_cur: ',n_k_cur
         write(6,*) 'n_Y_lm_csum_cur: ',n_Y_lm_csum_cur
         write(6,*) 'n_Y_lm_csum_: ',(n_Y_lm_csum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_Y_l_: ',(n_Y_l_(nk),nk=0,n_k_cur-1)
         write(6,*) 'lsq_oversample: ',lsq_oversample
         write(6,*) 'lsq_interpolation_order: ',lsq_interpolation_order
         write(6,*) 'n_residual_loading: ',n_residual_loading
      end if
c
      n_Y_l = n_Y_l_(n_k_cur-1)
      n_Y_lm_max = (n_Y_l+1)**2
      if (verbose.gt.1) then
         write(6,*) ' n_Y_l = ',n_Y_l
         write(6,*) ' n_Y_lm_max = ',n_Y_lm_max
      end if
c
      n_polar_a_oversample_max = nint(lsq_oversample*n_polar_a_(n_k_cur
     $     -1))
      if (verbose.gt.1) then
         write(6,*) ' lsq_oversample = ',lsq_oversample
         write(6,*) ' n_polar_a_oversample_max = '
     $        ,n_polar_a_oversample_max
      end if
c
      n_w_M_max = n_M*n_w_(n_k_cur-1)
      if (verbose.gt.1) then
         write(6,*) ' n_w_M_max = ',n_w_M_max
      end if
      allocate(polar_a__(0:ld_M-1,0:n_M-1))
      allocate(azimu_b__(0:ld_M-1,0:n_M-1))
      allocate(weight_CTF_k_c__(0:ld_M-1,0:n_M-1))
      allocate(n_w_M_(0:n_k_cur-1))
c
c
c
      do nk = 0,n_k_cur-1
         n_w_M_(nk) = 0
      end do
      n_A = 0
      do nk = 0,n_k_cur-1
         n_A = n_A + n_w_(nk)
      end do
      if (n_A.gt.ld_M) then
         write(6,*) 'Warning, ld_M: ',ld_M,' .neq. n_A: ',n_A
     $        ,' in test_residual_4'
      end if

      allocate(fftw_plan_frwd_(0:n_k_cur-1))
      allocate(fftw_plan_back_(0:n_k_cur-1))
      allocate(fftw_0in_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nk=0,n_k_cur-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nk),n_w_(nk)
     $        ,fftw_0in_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nk),n_w_(nk)
     $        ,fftw_out_(na),fftw_0in_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nk)
      enddo !do nk=0,n_k_cur-1

      allocate(M_transform_(0:ld_M-1))
      allocate(Y_slice_(0:ld_M-1))
      if (verbose.gt.0) then
         write(6,'(A)') ' M_transform_ - Y_slice_ = M_residual__'
      end if !if (verbose.gt.0) then
      do nM=0,n_M-1
         call cl1_c16(ld_M,M_transform_)
         nctf = nint(alpha2d__(nalpha_ctf_ind,nM))
         if ((verbose.gt.1).and.(nM.eq.0)) then
            write(6,'(A,I0,A,I0,A,I0)') ' nM: ' , nM , ' nctf: ' , nctf,
     $           ' ld_M: ' , ld_M
         end if !if ((verbose.gt.1).and.(nM.eq.0)) then
         na = 0
         do nk = n_k_low-1,n_k_cur-1
            if ((verbose.gt.2).and.(nM.eq.0)) then
               write(6,'(A,I0,A,I0)') ' nk = ' , nk , ' n_w_csum_(nk) '
     $              ,n_w_csum_(nk)
            end if !if verbose
            call get_lsqdata_2(1,I_M_sample_(nM),ld_M,M_k_p__,n_CTF
     $           ,ld_CTF,CTF_k_p_,grid_k_p_(nk),fftw_plan_frwd_(nk)
     $           ,fftw_plan_back_(nk) ,n_w_csum_(nk),n_w_(nk)
     $           ,fftw_0in_(na) ,fftw_out_(na),alpha2d__(0,nM)
     $           ,M_transform_(n_w_csum_(nk) -1),polar_a__(n_w_csum_(nk)
     $           -1 ,nM) ,azimu_b__(n_w_csum_(nk)-1,nM)
     $           ,weight_CTF_k_c__(n_w_csum_(nk)-1,nM) ,n_w_M_(nk))
            na = na + n_w_(nk)
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
         if ((verbose.gt.2).and.(nM.eq.0)) then
            write(format_string,'(A,I0,A)') '(A,' , ld_M ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'M_k_p__: ' , (M_k_p__(ntmp
     $           +I_M_sample_(nM)*ld_M),ntmp=0,ld_M-1)
            write(6,format_string) 'M_transform_: ' ,
     $           (M_transform_(ntmp),ntmp=0,ld_M-1)
         end if !if verbose
         call cl1_c16(ld_M,Y_slice_)
         do nk = n_k_low-1,n_k_cur-1
            if ((verbose.gt.3).and.(nM.eq.0)) then
               write(6,'(A,I0,A,I0)') ' nk = ' , nk , ' n_w_csum_(nk) '
     $              ,n_w_csum_(nk)
            end if !if verbose
            n_Y_l = n_Y_l_(nk)
            n_polar_a_oversample = nint(lsq_oversample*n_polar_a_(nk))
            call test_residual_multa(polar_a__(n_w_csum_(nk)-1,nM)
     $           ,azimu_b__(n_w_csum_(nk)-1,nM),n_w_M_(nk)
     $           ,Y_(n_Y_lm_csum_(nk)),n_Y_l,n_polar_a_oversample
     $           ,lsq_interpolation_order
     $           ,weight_CTF_k_c__(n_w_csum_(nk) -1,nM)
     $           ,Y_slice_(n_w_csum_(nk)-1))
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
         if ((verbose.gt.2).and.(nM.eq.0)) then
            write(format_string,'(A,I0,A)') '(A,' , ld_M ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'Y_slice_: ' , (Y_slice_(ntmp)
     $           ,ntmp=0,ld_M-1)
         end if !if verbose
         call cl1_c16(ld_M,M_residual__(nM*ld_M))
         do nk=max(2,n_k_low-1),n_k_cur-1
            do nw=0,n_w_(nk)-1
               M_residual__(n_w_csum_(nk)-1+nw + nM*ld_M) =
     $              M_transform_(n_w_csum_(nk)-1+nw)
     $              -Y_slice_(n_w_csum_(nk)-1+nw)
            enddo !do nw=0,n_w_(nk)-1
         enddo !         do nk=max(2,n_k_low-1),n_k_cur-1
         if ((verbose.gt.2).and.(nM.eq.0)) then
            call print_sub_c16(n_A,M_transform_
     $           ,' M_transform_')
            call print_sub_c16(n_A,Y_slice_
     $           ,' Y_slice_')
            call print_sub_c16(n_A,M_residual__(nM*ld_M)
     $           ,' M_residual__')
         end if !if verbose
      enddo                     ! do nM=0,n_M-1      

      if ((verbose.gt.0).and.(nM.eq.0)) then
         write(6,'(A)') ' initializing H__ = Et*R across images. '
      end if !if verbose
      if ((verbose.gt.0).and.(nM.eq.0)) then
         write(6,'(A,I0,A)') ' H__ requires ' , n_Y_lm_csum_cur*n_M ,
     $        ' elements'
         write(6,'(A,F6.3,A)') ' H__ requires ' , n_Y_lm_csum_cur*n_M*16
     $        *1.0d-9 ,' GB'
      end if !if verbose
c$$$      allocate(H_(0:n_Y_lm_csum_cur-1))
      allocate(H__(0:n_Y_lm_csum_cur*n_M-1))
      do nM=0,n_M-1
         do nk = n_k_low-1,n_k_cur-1
            n_Y_l = n_Y_l_(nk)
            n_polar_a_oversample = nint(lsq_oversample
     $           *n_polar_a_(nk))
            call test_residual_multah(polar_a__(n_w_csum_(nk)-1
     $           ,nM),azimu_b__(n_w_csum_(nk)-1,nM),n_w_M_(nk)
     $           ,M_residual__(n_w_csum_(nk)-1+nM*ld_M),n_Y_l
     $           ,n_polar_a_oversample ,lsq_interpolation_order
     $           ,weight_CTF_k_c__(n_w_csum_(nk)-1,nM)
     $           ,H__(n_Y_lm_csum_(nk) + nM*n_Y_lm_csum_cur))
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
      enddo                     ! do nM=0,n_M-1

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating M_loading_'
      end if !if verbose
      do nM_loading=0,n_residual_loading*n_M-1
         M_loading_(nM_loading) = dcmplx(adi_rand_f(rseed)-0.5d0,0.0d0)
      enddo !do nM_loading=0,n_residual_loading*n_M-1
      allocate(G_(0:n_residual_loading*n_Y_lm_csum_cur-1))

      n_iteration = n_residual_iteration
      do niteration=0,n_iteration-1
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' niteration ' , niteration
         end if !if verbose
         if (niteration.eq.0) then
            do nM_loading=0,n_residual_loading-1
               do nY_lm_sum=0,n_Y_lm_csum_cur-1
                  G_(nY_lm_sum+nM_loading*n_Y_lm_csum_cur) =
     $                 dcmplx(adi_rand_f(rseed)-0.5,0.0d0)
               enddo            !do nY_lm_sum=0,n_Y_lm_csum_cur-1
               call normalize_c16(n_Y_lm_csum_cur,G_(nM_loading
     $              *n_Y_lm_csum_cur))
            enddo               !do nM_loading=0,n_residual_loading*n_Y_lm_csum_cur-1
         end if                 !if (niteration.eq.0) then
         if (niteration.gt.0) then
            call cl1_c16(n_Y_lm_csum_cur*n_residual_loading,G_)
            do nM=0,n_M-1
               do nM_loading=0,n_residual_loading-1
                  call ac1_c16(n_Y_lm_csum_cur,M_loading_(nM_loading +
     $                 n_residual_loading*(nM)),dcmplx(0.0d0,0.0d0)
     $                 ,H__(0 + nM *n_Y_lm_csum_cur),G_(0 +nM_loading
     $                 *n_Y_lm_csum_cur))
               enddo            !do nM_loading=0,n_residual_loading-1
            enddo               ! do nM=0,n_M-1
            call gramschmidt_c16(n_Y_lm_csum_cur,n_residual_loading,G_)
            if ((verbose.gt.2) .and. (n_residual_loading.ge.2)) then
               write(6,*) ' iteration: ' , niteration
               nM_loading = 0
               write(6,*) ' loading 0: ' , al2_c16_f(n_Y_lm_csum_cur
     $              ,G_(0+nM_loading*n_Y_lm_csum_cur))
               nM_loading = 1
               write(6,*) ' loading 1: ' , al2_c16_f(n_Y_lm_csum_cur
     $              ,G_(0+nM_loading*n_Y_lm_csum_cur))
               write(6,*) ' <0,0>: ' , dot_c16_f(n_Y_lm_csum_cur,G_(0 +
     $              0*n_Y_lm_csum_cur),G_(0 + 0 *n_Y_lm_csum_cur))
               write(6,*) ' <0,1>: ' , dot_c16_f(n_Y_lm_csum_cur,G_(0 +
     $              0*n_Y_lm_csum_cur),G_(0 + 1 *n_Y_lm_csum_cur))
               write(6,*) ' <1,0>: ' , dot_c16_f(n_Y_lm_csum_cur,G_(0 +
     $              1*n_Y_lm_csum_cur),G_(0 + 0 *n_Y_lm_csum_cur))
               write(6,*) ' <1,1>: ' , dot_c16_f(n_Y_lm_csum_cur,G_(0 +
     $              1*n_Y_lm_csum_cur),G_(0 + 1 *n_Y_lm_csum_cur))
            end if !if verbose
         end if                 !if (niteration.gt.0) then
         do nM=0,n_M-1
            call dot_c16(n_Y_lm_csum_cur,H__(0 + nM*n_Y_lm_csum_cur)
     $           ,H__(0+ nM*n_Y_lm_csum_cur),HH)
            do nM_loading=0,n_residual_loading-1
               call dot_c16(n_Y_lm_csum_cur,H__(0 + nM*n_Y_lm_csum_cur)
     $              ,G_(0+nM_loading*n_Y_lm_csum_cur),HG)
               M_loading_(nM_loading + n_residual_loading*(nM)) =
     $              dreal(HG)/max(1.0d-15,dreal(HH))
            enddo               !do nM_loading=0,n_residual_loading-1
         enddo                  ! do nM=0,n_M-1
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,I5,' ,
     $           n_residual_loading ,'(2F10.4,1X))'
            do nM=0,n_M-1
               write(6,format_string) 'M_loading_: image ' , nM ,
     $              (M_loading_(nM_loading + n_residual_loading*(nM))
     $              ,nM_loading =0,n_residual_loading-1)
            enddo               !do nM=0,n_M-1
         end if                 !if verbosen

         if (verbose.gt.3) then
            write(format_string,'(A,I0,A)') '(' , n_residual_loading ,
     $           '(2F10.4,1X))'
            do nM=0,n_M-1
               write(6,format_string) (M_loading_(nM_loading +
     $              n_residual_loading*(nM)),nM_loading=0
     $              ,n_residual_loading-1)
            enddo               !do nM=0,n_M-1
         end if                 !if (verbose.gt.3) then

      enddo                     !do niteration=0,n_iteration-1
c$$$      deallocate(H_)
      deallocate(H__)
      deallocate(G_)
      deallocate(polar_a__)
      deallocate(azimu_b__)
      deallocate(weight_CTF_k_c__)
      deallocate(n_w_M_)

      if (verbose.gt.-1) then
         write(6,'(A)') ' Destroying fftw_plans for local use.'
      end if
      do nk=0,n_k_cur-1
         call dfftw_destroy_plan(fftw_plan_frwd_(nk))
         call dfftw_destroy_plan(fftw_plan_back_(nk))
      enddo !do nk=0,n_k_cur-1
      if (verbose.gt.-1) then
         write(6,'(A)') ' deallocating. '
      end if
      deallocate(fftw_plan_frwd_)
      deallocate(fftw_plan_back_)
      deallocate(fftw_0in_)
      deallocate(fftw_out_)

      if (verbose.gt.0) write(6,*)
     $     '[finished test_residual_4]'
      return
      end

