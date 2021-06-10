      subroutine get_svd_1(eps_target,n_svd_r_out,n_svd_d_out
     $     ,n_svd_l_out,svd_r_out_,svd_d_out_,svd_l_out_,svd_U_d_out_
     $     ,svd_s_out_,svd_V_r_out_,svd_unitnumber_out,svd_fname_out
     $     ,grid_p_,n_r,n_delta_v,delta_x_,delta_y_ ,flag_warning,R_max
     $     ,K_max,delta_max ,n_pixels)
      implicit none
      integer verbose
      data verbose / 1 /
      real *8 eps_target
      integer *4 n_svd_r_out,n_svd_d_out,n_svd_l_out
      real *8 svd_r_out_(0:0)
      real *8 svd_d_out_(0:0)
      integer *4 svd_l_out_(0:0)
      real *8 svd_U_d_out_(0:0)
      real *8 svd_s_out_(0:0)
      real *8 svd_V_r_out_(0:0)
      integer *4 svd_unitnumber_out
      character(len=1024) svd_fname_out
      include './dir_gen_Jsvd_1/gen_Jsvd_svddecl.txt'
      logical flag_warning
      real *8 grid_p_(0:0)
      integer n_r
      integer n_delta_v,ndv
      real *8 delta_x_(0:0),delta_y_(0:0)
      real *8 R_max,K_max,delta,delta_max,n_pixels
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_svd_1]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      R_max = 2*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndv=0,n_delta_v-1
         delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
         if (delta.gt.delta_max) then
            delta_max = delta
         end if
      enddo !do ndv=0,n_delta_v-1
      n_pixels = delta_max/dsqrt(2.0d0)*2*K_max
      if (verbose.gt.1) then
         write(6,'(A,F8.3,A,F8.3)') 'R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') 'K_max: ',K_max,'; n_pixels: '
     $        ,n_pixels
      end if
      if (K_max.gt.96 .and. flag_warning) then
         write(6,'(A,F8.3,A)') 'Warning, K_max ',K_max
     $        ,' too large in get_svd_1'
      end if
      if (n_pixels.gt.3 .and. flag_warning) then
         write(6,'(A,F8.3,A)') 'Warning, n_pixels ',n_pixels
     $        ,' too large in get_svd_1'
      end if
      if (eps_target.lt.0.1d0 .and. flag_warning) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in get_svd_1'
      end if
      if (.false.) then ! do nothing and continue to next line ;
      include './dir_gen_Jsvd_1/gen_Jsvd_svdpick.txt'
      end if
      include './dir_gen_Jsvd_1/gen_Jsvd_svdload.txt'
      if (verbose.gt.1) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
         write(6,'(A,I0)') 'svd_unitnumber: ',svd_unitnumber
      end if
      n_svd_r_out = n_svd_r
      n_svd_d_out = n_svd_d
      n_svd_l_out = n_svd_l
      call cp1_r8(n_svd_r,svd_r_,svd_r_out_)
      call cp1_r8(n_svd_d,svd_d_,svd_d_out_)
      call cp1_i4(n_svd_l,svd_l_,svd_l_out_)
      call cp1_r8(n_svd_d*n_svd_l,svd_U_d_,svd_U_d_out_)
      call cp1_r8(n_svd_l,svd_s_,svd_s_out_)
      call cp1_r8(n_svd_r*n_svd_l,svd_V_r_,svd_V_r_out_)
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_svd_1]'
      end if !if (verbose.gt.0) then
      end
