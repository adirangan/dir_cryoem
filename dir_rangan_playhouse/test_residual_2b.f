C***********************************************************************
      subroutine test_residual_2b(rseed,n_M,I_M_sample_,ld_M,M_k_p__
     $     ,n_CTF,ld_CTF,CTF_k_p_,n_alpha,alpha2d__,n_w_sum_,n_polar_a_
     $     ,quadrature_type_azimu_b ,grid_k_p_ ,n_w_,n_k_low,n_k_cur
     $     ,n_Y_lm_sum_cur,n_Y_lm_sum_ ,n_Y_l_ ,lsq_oversample
     $     ,lsq_interpolation_order ,Y_ ,M_transform__ ,Y_slice__
     $     ,M_residual__ ,n_M_loading ,n_M_iteration,M_loading_)
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
C     Y_(:)          complex *16: solution to least squares problem expressed in spherical harmonics in successive shells. ;
c     n_M_loading    integer: number of M_loading_ vectors per image. ;
c     n_M_iteration  integer: number of iterations to perform. ;
cC
C     OUTPUT: 
c
c     M_transform__(:) complex *16: images transformed. ;
c     Y_slice__(:)     complex *16: model sliced. ;
c     M_residual__(:)  complex *16: M_residual__ = M_transform__ - Y_slice__. ;
c     M_loading_(:)    complex *16: M_loading_ vectors. ;
c
C***********************************************************************
      implicit none
      integer verbose
      data verbose / 0 /
      integer *4 rseed
      real *8 adi_rand_f
      include 'adi_rand_f_define.f'
      integer n_M,I_M_sample_(0:n_M-1),ld_M,n_CTF,ld_CTF
      complex *16 M_k_p__(0:ld_M-1,0:0)
      complex *16 CTF_k_p_(0:ld_CTF*n_CTF-1)
      integer *4 n_alpha
      include 'nalpha_define.f'
      real *8 alpha2d__(0:n_alpha-1,0:n_M-1)
      integer n_w_sum_(0:n_k_cur-1)
      integer n_polar_a_(0:n_k_cur-1)
      integer n_w_(0:n_k_cur-1)
      integer quadrature_type_azimu_b
      real *8 grid_k_p_(0:n_k_cur-1)
      integer n_k_low,n_k_cur
      integer n_Y_lm_sum_cur
      integer n_Y_lm_sum_(0:n_k_cur-1)
      integer n_Y_l_(0:n_k_cur-1)
      complex *16 Y_(0:n_Y_lm_sum_cur-1)
      real *8 lsq_oversample
      integer lsq_interpolation_order
      complex *16 M_transform__(0:ld_M-1,0:n_M-1)
      complex *16     Y_slice__(0:ld_M-1,0:n_M-1)
      complex *16  M_residual__(0:ld_M-1,0:n_M-1)
      integer n_M_loading,n_M_iteration
      complex *16 M_loading_(0:n_M_loading*n_M-1)
      integer nM,nk,nA,n_A,nw,ntmp,nctf
      integer n_w_M_max,n_polar_a_oversample
      integer n_Y_lm_max,n_Y_l,n_polar_a_oversample_max
      real *8, allocatable :: polar_a__(:,:)
      real *8, allocatable :: azimu_b__(:,:)
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
         write(6,*) '[entering test_residual_2b]'
      end if
      if (verbose.gt.1) then
         write(6,*) 'ld_M: ',ld_M
         write(6,*) 'n_M: ',n_M
         write(6,*) 'n_w_sum_: ',(n_w_sum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_CTF: ',n_CTF
         write(6,*) 'ld_CTF: ',ld_CTF
         write(6,*) 'n_alpha: ',n_alpha
         write(6,*) 'n_polar_a_: ',(n_polar_a_(nk),nk=0,n_k_cur-1)
         write(6,*) 'quadrature_type_azimu_b: ',quadrature_type_azimu_b
         write(6,*) 'grid_k_p_: ',(grid_k_p_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_w_: ',(n_w_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_k_low: ',n_k_low
         write(6,*) 'n_k_cur: ',n_k_cur
         write(6,*) 'n_Y_lm_sum_cur: ',n_Y_lm_sum_cur
         write(6,*) 'n_Y_lm_sum_: ',(n_Y_lm_sum_(nk),nk=0,n_k_cur-1)
         write(6,*) 'n_Y_l_: ',(n_Y_l_(nk),nk=0,n_k_cur-1)
         write(6,*) 'lsq_oversample: ',lsq_oversample
         write(6,*) 'lsq_interpolation_order: ',lsq_interpolation_order
         write(6,*) 'n_M_loading: ',n_M_loading
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
     $        ,' in test_residual_2b'
      end if

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating M_transform__'
      end if !if (verbose.gt.0) then
      do nM=0,n_M-1
         nctf = nint(alpha2d__(nalpha_ctf_ind,nM))
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0)') ' nM: ' , nM , ' nctf: ' , nctf,
     $           'ld_M: ' , ld_M
         end if
         do nk = n_k_low-1,n_k_cur-1
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I0)') ' nk = ' , nk , ' n_w_sum_(nk) '
     $              ,n_w_sum_(nk)
            end if
            call get_lsqdata_0(1,ld_M,M_k_p__(0,I_M_sample_(nM)),1
     $           ,ld_CTF,CTF_k_p_(nctf*ld_CTF),grid_k_p_(nk)
     $           ,n_w_sum_(nk),n_w_(nk),n_alpha,alpha2d__(0,nM)
     $           ,M_transform__(n_w_sum_(nk)-1,nM)
     $           ,polar_a__(n_w_sum_(nk)-1,nM),azimu_b__(n_w_sum_(nk)-1
     $           ,nM),weight_CTF_k_c__(n_w_sum_(nk)-1,nM) ,n_w_M_(nk))
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,' , ld_M ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'M_k_p__: ' , (M_k_p__(ntmp
     $           ,I_M_sample_(nM)),ntmp=0,ld_M-1)
            write(6,format_string) 'M_transform__: ' ,
     $           (M_transform__(ntmp,nM),ntmp=0,ld_M-1)
         end if !if (verbose.gt.2) then
      enddo                     ! do nM=0,n_M-1

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating Y_slice__'
      end if !if (verbose.gt.0) then
      do nM=0,n_M-1
         nctf = nint(alpha2d__(nalpha_ctf_ind,nM))
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0)') ' nM: ' , nM , ' nctf: '
     $           , nctf, ' ld_M: ' , ld_M
         end if
         do nk = n_k_low-1,n_k_cur-1
            if (verbose.gt.1) then
               write(6,'(A,I0,A,I0)') ' nk = ' , nk , ' n_w_sum_(nk) '
     $              ,n_w_sum_(nk)
            end if
            n_Y_l = n_Y_l_(nk)
            n_polar_a_oversample = nint(lsq_oversample*n_polar_a_(nk))
            call test_residual_multa(polar_a__(n_w_sum_(nk)-1,nM)
     $           ,azimu_b__(n_w_sum_(nk)-1,nM),n_w_M_(nk)
     $           ,Y_(n_Y_lm_sum_(nk)),n_Y_l,n_polar_a_oversample
     $           ,lsq_interpolation_order ,weight_CTF_k_c__(n_w_sum_(nk)
     $           -1,nM) ,Y_slice__(n_w_sum_(nk)-1,nM))
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,' , ld_M ,
     $           '(2F8.4,1X))'
            write(6,format_string) 'Y_slice__: ' , (Y_slice__(ntmp,nM)
     $           ,ntmp=0,ld_M-1)
         end if !if (verbose.gt.2) then
      enddo                     ! do nM=0,n_M-1      

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating M_residual__'
      end if !if (verbose.gt.0) then
      do nM=0,n_M-1
         do nA=0,ld_M-1
            M_residual__(nA,nM) = cmplx(0.0d0,0.0d0)
         enddo !do nA=0,ld_M-1
         do nk=max(2,n_k_low-1),n_k_cur-1
            do nw=0,n_w_(nk)-1
               M_residual__(n_w_sum_(nk)-1+nw,nM) =
     $              M_transform__(n_w_sum_(nk)-1+nw,nM)
     $              -Y_slice__(n_w_sum_(nk)-1+nw,nM)
            enddo !do nw=0,n_w_(nk)-1
         enddo !         do nk=max(2,n_k_low-1),n_k_cur-1
      enddo                     ! do nM=0,n_M-1      

      if (verbose.gt.0) then
         write(6,'(A)') ' initializing H__ = Et*R across images. '
      end if !if (verbose.gt.0) then
      if (verbose.gt.0) then
         write(6,'(A,I0,A)') ' H__ requires ' , n_Y_lm_sum_cur*n_M ,
     $        ' elements'
         write(6,'(A,F6.3,A)') ' H__ requires ' , n_Y_lm_sum_cur*n_M*16
     $        *1.0d-9 ,' GB'
      end if !if (verbose.gt.0) then
c$$$      allocate(H_(0:n_Y_lm_sum_cur-1))
      allocate(H__(0:n_Y_lm_sum_cur*n_M-1))
      do nM=0,n_M-1
         do nk = n_k_low-1,n_k_cur-1
            n_Y_l = n_Y_l_(nk)
            n_polar_a_oversample = nint(lsq_oversample
     $           *n_polar_a_(nk))
            call test_residual_multah(polar_a__(n_w_sum_(nk)-1
     $           ,nM),azimu_b__(n_w_sum_(nk)-1,nM),n_w_M_(nk)
     $           ,M_residual__(n_w_sum_(nk)-1,nM),n_Y_l
     $           ,n_polar_a_oversample ,lsq_interpolation_order
     $           ,weight_CTF_k_c__(n_w_sum_(nk)-1,nM)
     $           ,H__(n_Y_lm_sum_(nk) + nM*n_Y_lm_sum_cur))
         enddo                  ! do nk = n_k_low-1,n_k_cur-1
      enddo                     ! do nM=0,n_M-1

      if (verbose.gt.0) then
         write(6,'(A)') ' calculating M_loading_'
      end if !if (verbose.gt.0) then
      do nM_loading=0,n_M_loading*n_M-1
         M_loading_(nM_loading) = cmplx(adi_rand_f(rseed)-0.5d0,0.0d0)
      enddo !do nM_loading=0,n_M_loading*n_M-1
      allocate(G_(0:n_M_loading*n_Y_lm_sum_cur-1))

      n_iteration = n_M_iteration
      do niteration=0,n_iteration-1
         if (verbose.gt.1) then
            write(6,'(A,I0)') ' niteration ' , niteration
         end if !if (verbose.gt.1) then
         if (niteration.eq.0) then
            do nM_loading=0,n_M_loading-1
               do nY_lm_sum=0,n_Y_lm_sum_cur-1
                  G_(nY_lm_sum+nM_loading*n_Y_lm_sum_cur) =
     $                 cmplx(adi_rand_f(rseed)-0.5,0.0d0)
               enddo            !do nY_lm_sum=0,n_Y_lm_sum_cur-1
               call normalize_c16(n_Y_lm_sum_cur,G_(nM_loading
     $              *n_Y_lm_sum_cur))
            enddo               !do nM_loading=0,n_M_loading*n_Y_lm_sum_cur-1
         end if                 !if (niteration.eq.0) then
         if (niteration.gt.0) then
            call cl1_c16(n_Y_lm_sum_cur*n_M_loading,G_)
            do nM=0,n_M-1
               do nM_loading=0,n_M_loading-1
                  call ac1_c16(n_Y_lm_sum_cur,M_loading_(nM_loading +
     $                 n_M_loading*(nM)),cmplx(0.0d0,0.0d0),H__(0 + nM
     $                 *n_Y_lm_sum_cur),G_(0 +nM_loading
     $                 *n_Y_lm_sum_cur))
               enddo            !do nM_loading=0,n_M_loading-1
            enddo               ! do nM=0,n_M-1
            call gramschmidt_c16(n_Y_lm_sum_cur,n_M_loading,G_)
            if (verbose.gt.2 .and. n_M_loading.ge.2) then
               write(6,*) ' iteration: ' , niteration
               nM_loading = 0
               write(6,*) ' loading 0: ' , al2_c16_f(n_Y_lm_sum_cur,G_(0
     $              +nM_loading*n_Y_lm_sum_cur))
               nM_loading = 1
               write(6,*) ' loading 1: ' , al2_c16_f(n_Y_lm_sum_cur,G_(0
     $              +nM_loading*n_Y_lm_sum_cur))
               write(6,*) ' <0,0>: ' , dot_c16_f(n_Y_lm_sum_cur,G_(0 + 0
     $              *n_Y_lm_sum_cur),G_(0 + 0 *n_Y_lm_sum_cur))
               write(6,*) ' <0,1>: ' , dot_c16_f(n_Y_lm_sum_cur,G_(0 + 0
     $              *n_Y_lm_sum_cur),G_(0 + 1 *n_Y_lm_sum_cur))
               write(6,*) ' <1,0>: ' , dot_c16_f(n_Y_lm_sum_cur,G_(0 + 1
     $              *n_Y_lm_sum_cur),G_(0 + 0 *n_Y_lm_sum_cur))
               write(6,*) ' <1,1>: ' , dot_c16_f(n_Y_lm_sum_cur,G_(0 + 1
     $              *n_Y_lm_sum_cur),G_(0 + 1 *n_Y_lm_sum_cur))
            end if !if (verbose.gt.2) then
         end if                 !if (niteration.gt.0) then
         do nM=0,n_M-1
            call dot_c16(n_Y_lm_sum_cur,H__(0 + nM*n_Y_lm_sum_cur),H__(0
     $           + nM*n_Y_lm_sum_cur),HH)
            do nM_loading=0,n_M_loading-1
               call dot_c16(n_Y_lm_sum_cur,H__(0 + nM*n_Y_lm_sum_cur)
     $              ,G_(0+nM_loading*n_Y_lm_sum_cur),HG)
               M_loading_(nM_loading + n_M_loading*(nM)) = dreal(HG)
     $              /max(1.0d-15,dreal(HH))
            enddo               !do nM_loading=0,n_M_loading-1
         enddo                  ! do nM=0,n_M-1
         if (verbose.gt.2) then
            write(format_string,'(A,I0,A)') '(A,I5,' , n_M_loading ,
     $           '(2F10.4,1X))'
            do nM=0,n_M-1
               write(6,format_string) 'M_loading_: image ' , nM ,
     $              (M_loading_(nM_loading + n_M_loading*(nM))
     $              ,nM_loading =0,n_M_loading-1)
            enddo               !do nM=0,n_M-1
         end if                 !if (verbose.gt.2) then

         if (verbose.gt.3) then
            write(format_string,'(A,I0,A)') '(' , n_M_loading ,
     $           '(2F10.4,1X))'
            do nM=0,n_M-1
               write(6,format_string) (M_loading_(nM_loading +
     $              n_M_loading*(nM)),nM_loading=0,n_M_loading-1)
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

      if (verbose.gt.0) write(6,*)
     $     '[finished test_residual_2b]'
      return
      end

