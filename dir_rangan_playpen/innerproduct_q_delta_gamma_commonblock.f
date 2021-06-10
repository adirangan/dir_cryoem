      subroutine innerproduct_q_delta_gamma_commonblock(type_flag
     $     ,eps_target,n_delta_x,delta_x_,n_delta_y,delta_y_,n_gamma_z
     $     ,gamma_z_,n_r,grid_p_,n_w_,Wlast_,n_A,S_q_,M_q_,C_Z,Z_q_)
      implicit none
      include '/usr/include/fftw3.f'
c$$$      Uses svd-expansion defined via n_svd_r, .. , svd_V_r_
c$$$      common block to hold svd-expansions
      include 'gen_Jsvd_blockdecl.txt'
      integer verbose
      data verbose / 1 /
      logical warning_flag
      data warning_flag / .true. /
      pointer (p_n_svd_r,n_svd_r)
      pointer (p_n_svd_d,n_svd_d)
      pointer (p_n_svd_l,n_svd_l)
      pointer (p_svd_l_,svd_l_)
      pointer (p_svd_r_,svd_r_)
      pointer (p_svd_d_,svd_d_)
      pointer (p_svd_s_,svd_s_)
      pointer (p_svd_U_d_,svd_U_d_)
      pointer (p_svd_V_r_,svd_V_r_)
      integer n_svd_r,n_svd_d,n_svd_l,svd_l_(0:0)
      real *8 svd_r_(0:0),svd_d_(0:0)
      real *8 svd_s_(0:0)
      real *8 svd_U_d_(0:0)
      real *8 svd_V_r_(0:0)      
      integer type_flag,n_delta_x,n_delta_y,n_gamma_z
      real *8 eps_target
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 Wlast_(0:4*n_w_(n_r-1)+16-1)
      complex *16 S_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16 C_Z,Z_q_(0:n_delta_x*n_delta_y*n_gamma_z-1)
      integer n_fftw
      parameter (n_fftw=1024)
      double complex fftw_in1_(0:n_fftw-1)
      double complex fftw_in2_(0:n_fftw-1)
      double complex fftw_out_(0:n_fftw-1)
      integer*8 fftw_plan_back
      integer*8 fftw_plan_frwd      
      real *8 pi
      real *8 R_max,K_max,delta,delta_max,N_pixels
      complex *16, allocatable :: Z_svds_(:)
      complex *16, allocatable :: Z_svdd_(:)
      complex *16, allocatable :: Z_tmpC_(:)
      integer ndx,ndy,ngz,nw,n_w_max,l
      real *8 gamma_z
      integer nC1,nC2,nC3,nC4,nC5
      complex *16 Z_q
c$$$      parameters for timing
      real *8 timing_tic,timing_toc
      if (verbose.gt.0) then
         write(6,'(A,I0)') '[entering innerproduct_q_delta_gamma_commonb
     $lock], type_flag: ',type_flag
      end if
      pi = 4.0*atan(1.0)
      R_max = 2*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
            delta = dsqrt(delta_x_(ndx)**2 + delta_y_(ndy)**2)
            if (delta.gt.delta_max) then
               delta_max = delta
            end if
         enddo
      enddo
      N_pixels = delta_max/dsqrt(2.0d0)*2*K_max
      if (verbose.gt.0) then
         write(6,'(A,F8.3,A,F8.3)') ' R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') ' K_max: ',K_max,'; N_pixels: '
     $        ,N_pixels
      end if
      if (K_max.gt.48 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, K_max ',K_max
     $        ,' too large in innerproduct_q_delta_gamma_commonblock'
      end if
      if (N_pixels.gt.3 .and. warning_flag) then
         write(6,'(A,F8.3,A)') 'Warning, N_pixels ',N_pixels
     $        ,' too large in innerproduct_q_delta_gamma_commonblock'
      end if
      if (eps_target.lt.0.1d0 .and. warning_flag) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in innerproduct_q_delta_gamma_commonblock'
      end if
      if (.false.) then
      include 'gen_Jsvd_blockpick.txt'
      end if
      n_w_max = n_w_(n_r-1)
      call dfftw_plan_dft_1d_(fftw_plan_frwd,n_w_max,fftw_in1_,fftw_out_
     $     ,FFTW_FORWARD ,FFTW_ESTIMATE)
      call dfftw_plan_dft_1d_(fftw_plan_back,n_w_max,fftw_out_,fftw_in2_
     $     ,FFTW_BACKWARD,FFTW_ESTIMATE)
      allocate(Z_svds_(0:n_svd_l*n_w_max-1))
      allocate(Z_svdd_(0:n_svd_l*n_delta_x*n_delta_y-1))
      allocate(Z_tmpC_(0:n_w_max*n_delta_x*n_delta_y-1))

c$$$      multiply first, fft second
      if (type_flag.eq.0) then
      timing_tic = second()
      call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_
     $     ,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_ ,M_q_,Z_svds_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svds, total time ',timing_toc
     $        -timing_tic
      end if
      timing_tic = second()
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svdd, total time ',timing_toc
     $        -timing_tic
      end if
      timing_tic = second()
      nC1 = 0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
c$$$            nC4 = ndx*n_w_max + ndy*n_w_max*n_delta_x
            nC4 = nC1
            do nw=0,n_w_max-1
c$$$               nC1 = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
               Z_tmpC_(nC1) = cmplx( 0.0 , 0.0 )
               nC2 = ndx*n_svd_l + ndy*n_svd_l*n_delta_x
               nC3 = nw*n_svd_l
               do l=0,n_svd_l-1
                  Z_tmpC_(nC1) = Z_tmpC_(nC1) + Z_svdd_(nC2)
     $                 *Z_svds_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo
               nC1 = nC1 + 1
            enddo
            call cp1_c16(n_w_max,Z_tmpC_(nC4),fftw_out_)
            call dfftw_execute_(fftw_plan_back)
            call cp1_c16(n_w_max,fftw_in2_,Z_tmpC_(nC4))
            do ngz=0,n_gamma_z-1
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmpC_(nC4),+gamma_z
     $              ,Z_q)
               nC5 = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               Z_q_(nC5) = Z_q/(n_r**4)/C_Z
c$$$               Z_q_(nC5) = Z_q
            enddo
         enddo
      enddo
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculation of Z_q_, total time ',timing_toc
     $        -timing_tic
      end if

c$$$      fft first, multiply second
      else if (type_flag.eq.1) then
      timing_tic = second()
      call innerproduct_q__k_svds(n_svd_r,svd_r_,n_svd_l,svd_l_,svd_s_
     $     ,svd_V_r_,n_r,grid_p_,n_w_,n_A,S_q_ ,M_q_,Z_svds_)
      do l=0,n_svd_l-1
         call cps_c16(n_w_max,Z_svds_(l),n_svd_l,fftw_out_,1)
         call dfftw_execute_(fftw_plan_back)
         call cps_c16(n_w_max,fftw_in2_,1,Z_svds_(l),n_svd_l)
      enddo
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svds, total time ',timing_toc
     $        -timing_tic
      end if
      timing_tic = second()
      call innerproduct_q__k_svdd(n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_
     $     ,n_delta_x,delta_x_,n_delta_y,delta_y_,Z_svdd_)
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished innerproduct_q__k_svdd, total time ',timing_toc
     $        -timing_tic
      end if
      timing_tic = second()
      nC1 = 0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
c$$$            nC4 = ndx*n_w_max + ndy*n_w_max*n_delta_x
            nC4 = nC1
            do nw=0,n_w_max-1
c$$$               nC1 = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
               Z_tmpC_(nC1) = cmplx( 0.0 , 0.0 )
               nC2 = ndx*n_svd_l + ndy*n_svd_l*n_delta_x
               nC3 = nw*n_svd_l
               do l=0,n_svd_l-1
                  Z_tmpC_(nC1) = Z_tmpC_(nC1) + Z_svdd_(nC2)
     $                 *Z_svds_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo
               nC1 = nC1 + 1
            enddo
            do ngz=0,n_gamma_z-1
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmpC_(nC4),+gamma_z
     $              ,Z_q)
               nC5 = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               Z_q_(nC5) = Z_q/(n_r**4)/C_Z
c$$$               Z_q_(nC5) = Z_q
            enddo
         enddo
      enddo
      timing_toc = second()
      if (verbose.gt.0) then
         write(6,'(A,F8.5)')
     $        ' finished calculation of Z_q_, total time ',timing_toc
     $        -timing_tic
      end if
      end if

      call dfftw_destroy_plan_(fftw_plan_back)
      call dfftw_destroy_plan_(fftw_plan_frwd)
      deallocate(Z_svds_)
      deallocate(Z_svdd_)
      deallocate(Z_tmpC_)
      if (verbose.gt.0) then
         write(6,*) '[finished innerproduct_q_delta_gamma_commonblock]'
      end if
      end
