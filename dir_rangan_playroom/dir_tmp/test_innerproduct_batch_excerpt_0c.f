      subroutine test_innerproduct_batch_excerpt_0c(n_delta_x,n_delta_y
     $     ,n_gamma_z,gamma_z_,n_w_max,n_svd_l,Z_svdd_,Z_svds_,Z_tmpC_
     $     ,fftw_plan_back_last,fftw_in1_last_,fftw_out_last_,n_r,C_M
     $     ,C_S_,C_S_use_,gamma_z_est,Z_q_)
      implicit none
      include '/usr/include/fftw3.f'
      integer n_delta_x,n_delta_y,n_gamma_z,n_w_max,n_svd_l,n_r
      real *8 gamma_z_(0:n_gamma_z-1),gamma_z_est
      complex *16 Z_svdd_(0:0),Z_svds_(0:0),Z_tmpC_(0:0)
      integer *8 fftw_plan_back_last
      complex *16 fftw_in1_last_(0:0),fftw_out_last_(0:0)
      complex *16 C_M,C_S_(0:0),Z_q,C_S_use_(0:0),C_S_use,C_Z_use
      complex *16 Z_q_(0:0)
      integer nC1,nC2,nC3,nC4,nC5,ndx,ndy,ngz,nw,l
      real *8 ngz_est_v,dz_pre,dz_pos,alpha,beta
      integer ngz_est_pre,ngz_est_pos,ngz_pre,ngz_pos
      real *8 gamma_z
      real *8 pi
      pi = 4.0*atan(1.0)
c$$$      call get_C_S_use_(gamma_z_est,n_gamma_z,C_S_,C_S_use_)
      ngz_est_v = n_gamma_z*(gamma_z_est)/(2*pi)
      ngz_est_pre = floor(ngz_est_v)
      dz_pre = dabs(ngz_est_pre - ngz_est_v)
      ngz_est_pos = ceiling(ngz_est_v)
      dz_pos = dabs(ngz_est_pos - ngz_est_v)
      if (ngz_est_pos.eq.n_gamma_z) then
         ngz_est_pos = 0
      end if
      if (dz_pre+dz_pos.le.0.0d0) then
         alpha = 0.0d0
         beta = 1.0d0
      else
         alpha = dz_pre/(dz_pre+dz_pos)
         beta =  dz_pos/(dz_pre+dz_pos)
      end if
      nC1 = 0
      do ndy=0,n_delta_y-1
         do ndx=0,n_delta_x-1
c$$$        nC4 = ndx*n_w_max + ndy*n_w_max*n_delta_x
            nC4 = nC1
            do nw=0,n_w_max-1
c$$$           nC1 = nw + ndx*n_w_max + ndy*n_w_max*n_delta_x
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
            call cp1_c16(n_w_max,Z_tmpC_(nC4),fftw_out_last_)
            call dfftw_execute_(fftw_plan_back_last)
            call cp1_c16(n_w_max,fftw_in1_last_,Z_tmpC_(nC4))
            do ngz=0,n_gamma_z-1
               ngz_pre = ngz + ngz_est_pre
               call periodize_i(ngz_pre,0,n_gamma_z,ngz_pre)
               ngz_pos = ngz + ngz_est_pos
               call periodize_i(ngz_pos,0,n_gamma_z,ngz_pos)
               C_S_use = beta*C_S_(ngz_pre) + alpha*C_S_(ngz_pos)
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmpC_(nC4),
     $              +gamma_z,Z_q)
               nC5 = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
               if (zabs(C_M*C_S_use).le.1.0d-15) then
                  C_Z_use = 1.0d0
               else
                  C_Z_use = C_M*C_S_use
               end if
               Z_q_(nC5) = Z_q/(n_r**4)/C_Z_use
            enddo
         enddo
      enddo
      end
