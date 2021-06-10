      subroutine test_innerproduct_batch_excerpt_2(n_delta_x,n_delta_y
     $     ,n_gamma_z,gamma_z_,n_w_max,n_svd_l,Z_svdd_,Z_svdr_,Z_tmpC_
     $     ,n_r,C_M,C_S_,gamma_z_est,Z_q_)
      implicit none
      integer n_delta_x,n_delta_y,n_gamma_z,n_w_max,n_svd_l,n_r
      real *8 gamma_z_(0:n_gamma_z-1),gamma_z_est
      complex *16 Z_svdd_(0:0),Z_svdr_(0:0),Z_tmpC_(0:0)
      complex *16 C_M,C_S_(0:0),Z_q,C_S_use,C_Z_use
      complex *16 Z_q_(0:0)
      integer nC1,nC2,nC3,nC4,nC5,ndx,ndy,ngz,nw,l
      real *8 gamma_z
      real *8 pi
      complex *16, allocatable :: C_S_use_(:)
      pi = 4.0*atan(1.0)
      allocate(C_S_use_(0:n_gamma_z-1))
      call get_C_S_use_(gamma_z_est,n_gamma_z,C_S_,C_S_use_)
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
     $                 *Z_svdr_(nC3)
                  nC2 = nC2 + 1
                  nC3 = nC3 + 1
               enddo
               nC1 = nC1 + 1
            enddo
            do ngz=0,n_gamma_z-1
               C_S_use = C_S_use_(ngz)
               gamma_z = gamma_z_(ngz)
               call interp1_c16(n_w_max,0.0d0,2*pi,Z_tmpC_(nC4),+gamma_z
     $              ,Z_q)
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
      deallocate(C_S_use_)
      end
