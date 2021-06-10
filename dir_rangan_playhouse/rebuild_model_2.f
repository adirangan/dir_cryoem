C***********************************************************************
      subroutine rebuild_model_2(n_M,I_M_sample_,ld_M,M_k_p__,n_CTF
     $     ,ld_CTF,CTF_k_p_,alpha2d__,n_w_sum_,n_polar_a_
     $     ,quadrature_type_azimu_b ,grid_k_p_,n_w_,n_k_low,n_k_cur
     $     ,n_Y_lm_sum_ ,n_Y_l_ ,lsq_oversample,lsq_interpolation_order
     $     ,lsq_eps ,Y_)
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
c     alpha2d__(3,*) real *8: 2d array of image parameters (see nalpha_define.f). ;
c     n_w_sum_()     integer: indexing array for points on successive circles within each image. ; 
c                    n_w_sum_(nk) represents the number of points (in angle-w) on ring at radisu grid_k_p_(nk). ;
c     n_polar_a_()   integer:  number of quadrature nodes in polar_a on sphere defined by index n_k_cur. ;
c     quadrature_type_azimu_b integer: flag determining quadrature scheme in azimu_b direction. ;
c                       0 = same on all latitudes (which oversamples poles). ;
c                       1 = adaptive (appropriately samples toward poles). ;
c                    if quadrature_type_azimu_b=0, then n_azimu_b = nint(n_polar_a_(nk)*phi_over). ;
c                    if quadrature_type_azimu_b=1, then n_azimu_b = nint(n_polar_a_(nk)*phi_over*sin(polar_a)). ;
c                    Note: phi_over is set in getgridph (typically = 2). ;
c     grid_k_p_()    real *8: radius associated with successive circles in templates in k-space polar coordinates. ;
c     n_w_()         integer: number of output points on successive circles in templates in k-space polar coordinates. ;
c     n_k_low        integer: index of lowest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ;
c     n_k_cur        integer: index of highest frequency sphere under consideration. Note that this runs from 1 to n_k_p_max. ;
c     n_Y_lm_sum_()  integer: array of length n_k_p_max indicating where the Y_ for that shell begins. ;
c     n_Y_l_()       integer: array of length n_k_p_max defining the orders of the various Y_ on successive shells. ;
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
      integer n_w_sum_(0:0),n_polar_a_(0:0)
      integer quadrature_type_azimu_b
      include '/usr/include/fftw3.f'
      include 'nalpha_define.f'
      integer n_w_(0:0)
      integer n_k_low,n_k_cur
      integer *8, allocatable :: fftw_plan_frwd_(:)
      integer *8, allocatable :: fftw_plan_back_(:)
      complex *16, allocatable :: fftw_in1_(:)
      complex *16, allocatable :: fftw_out_(:)
      integer n_Y_lm_sum_(0:0),n_Y_l_(0:0)
      integer lsq_interpolation_order
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      real *8 grid_k_p_(0:n_k_cur-1)
      real *8 lsq_oversample,lsq_eps
      complex *16 M_k_p__(0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      complex *16 Y_(0:0)
      integer nk,n_A,na,n_w,nw
      integer n_iteration,niteration
      integer n_w_M,n_w_M_max
      integer n_Y_lm_max,n_Y_lm,n_Y_l,n_polar_a_oversample
      real *8, allocatable :: polar_a_(:)
      real *8, allocatable :: azimu_b_(:)
      complex *16, allocatable :: Y_cur_(:)
      complex *16, allocatable :: M_k_c_(:)
      complex *16, allocatable :: weight_CTF_k_c_(:)
      external multaha
      character(len=1024) format_string
      if (verbose.gt.0) then
         write(6,*) '[entering rebuild_model_2]'
      end if
      if (verbose.gt.1) then
         write(6,*) 'ld_M: ',ld_M
         write(6,*) 'n_M: ',n_M
         write(6,*) 'n_w_sum_: ',(n_w_sum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_polar_a_: ',(n_polar_a_(nk),nk=0,n_k_cur-1)
         write(6,*) 'quadrature_type_azimu_b: ',quadrature_type_azimu_b
         write(6,*) 'grid_k_p_: ',(grid_k_p_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_w_: ',(n_w_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_k_low: ',n_k_low
         write(6,*) 'n_k_cur: ',n_k_cur
         write(6,*) 'n_Y_lm_sum_: ',(n_Y_lm_sum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_Y_l_: ',(n_Y_l_(nk),nk=0,n_k_cur-1)
         write(6,*) 'lsq_oversample: ',lsq_oversample
         write(6,*) 'lsq_interpolation_order: ',lsq_interpolation_order
      end if

      n_A = 0
      do nk = 0,n_k_cur-1
         n_A = n_A + n_w_(nk)
      end do !do nk = 0,n_k_cur-1
      if (n_A.gt.ld_M) then
         write(6,*) 'Warning, ld_M: ',ld_M
     $        ,' .neq. n_A: ',n_A,' in rebuild_model_2'
      end if
      if (verbose.gt.1) then
         write(6,'(2(A,I0))') ' n_k_cur ' , n_k_cur , ' n_A ' , n_A
      end if !if (verbose.gt.1) then
c
      allocate(fftw_plan_frwd_(0:n_k_cur-1))
      allocate(fftw_plan_back_(0:n_k_cur-1))
      allocate(fftw_in1_(0:n_A-1))
      allocate(fftw_out_(0:n_A-1))
      na = 0
      do nk=0,n_k_cur-1
         call dfftw_plan_dft_1d_(fftw_plan_frwd_(nk),n_w_(nk)
     $        ,fftw_in1_(na),fftw_out_(na),FFTW_FORWARD,FFTW_MEASURE) 
         call dfftw_plan_dft_1d_(fftw_plan_back_(nk),n_w_(nk)
     $        ,fftw_out_(na),fftw_in1_(na),FFTW_BACKWARD,FFTW_MEASURE) 
         na = na + n_w_(nk)
      enddo !do nk=0,n_k_cur-1
c
      n_Y_l = n_Y_l_(n_k_cur-1)
      n_Y_lm_max = (n_Y_l+1)**2
      if (verbose.gt.1) then
         write(6,*) ' n_Y_l = ',n_Y_l
         write(6,*) ' n_Y_lm_max = ',n_Y_lm_max
      end if
      allocate(Y_cur_(0:1+n_Y_lm_max))
      call cs1_c16(n_Y_lm_max,Y_cur_)
c
      n_polar_a_oversample = nint(lsq_oversample*n_polar_a_(n_k_cur-1))
      if (verbose.gt.1) then
         write(6,*) ' lsq_oversample = ',lsq_oversample
         write(6,*) ' n_polar_a_oversample = ',n_polar_a_oversample
      end if
c
      n_w_M_max = n_M*n_w_(n_k_cur-1)
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
c
      na = 0
      do nk = 0,n_k_cur-1
         if (verbose.gt.1) then
            write(6,'(2(A,I0))') ' nk ' , nk , ' n_k_cur ' , n_k_cur
         end if
         if (nk.ge.0) then
            call get_lsqdata_2(n_M,I_M_sample_,ld_M,M_k_p__,n_CTF,ld_CTF
     $           ,CTF_k_p_,grid_k_p_(nk),fftw_plan_frwd_(nk)
     $           ,fftw_plan_back_(nk),n_w_sum_(nk),n_w_(nk)
     $           ,fftw_in1_(na),fftw_out_(na),alpha2d__,M_k_c_ ,polar_a_
     $           ,azimu_b_ ,weight_CTF_k_c_ ,n_w_M)
            if (verbose.gt.1) then
               n_w = 1*n_w_(nk)
               write(format_string,'(A,I0,A)') '(A,',n_w*2,'F8.3)'
               write(6,format_string) 'M_k_c_: ',(M_k_c_(nw),nw=0,n_w
     $              -1)
               write(6,format_string) 'weight_CTF_k_c_: '
     $              ,(weight_CTF_k_c_(nw),nw=0,n_w-1)
            end if              !verbose
c
c     set up parameters for tensor product grid on sphere
c     to be used in least-squares solve (for interpolation
c     to arbitrary points)
c     
            n_Y_l = n_Y_l_(nk)
            n_Y_lm = (n_Y_l+1)**2
            n_iteration = 1000
            niteration = 0
            call lsqctfsolve(M_k_c_,polar_a_,azimu_b_,n_w_M,Y_cur_,n_Y_l
     $           ,nint(lsq_oversample*n_polar_a_(nk))
     $           ,lsq_interpolation_order,weight_CTF_k_c_,multaha
     $           ,lsq_eps,n_iteration ,niteration)
            call cp1_c16(n_Y_lm,Y_cur_,Y_(n_Y_lm_sum_(nk)))
         end if                 !if (nk.ge.0) then
         na = na + n_w_(nk)
      enddo                     !do nk = n_k_low-1,n_k_cur-1

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
      deallocate(fftw_in1_)
      deallocate(fftw_out_)

      if (verbose.gt.0) write(6,*)
     $     '[finished rebuild_model_2]'
      return
      end

c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
