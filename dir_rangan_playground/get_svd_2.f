!> Doxygen comment: ;\n
!> load svd_U_d_ and svd_V_r_ from disc. ;\n
!> takes in certain parameters (e.g., eps_target, n_pixel). ;\n
      subroutine get_svd_2(eps_target,n_svd_r_out,n_svd_d_out
     $     ,n_svd_l_out,svd_r_out_,svd_d_out_,svd_l_out_,svd_U_d_out_
     $     ,svd_s_out_,svd_V_r_out_,svd_unitnumber_out,svd_fname_out
     $     ,grid_p_,n_r,n_delta_v,delta_x_,delta_y_ ,flag_warning,R_max
     $     ,K_max,delta_max ,n_pixel)
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 eps_target
      integer *4 n_svd_r_out,n_svd_d_out,n_svd_l_out
      real *8 svd_r_out_(0:0)
      real *8 svd_d_out_(0:0)
      integer *4 svd_l_out_(0:0)
      real *8 svd_U_d_out_(0:0)
      real *8 svd_s_out_(0:0)
      real *8 svd_V_r_out_(0:0)
      logical, allocatable :: flag_s_(:)
      integer *4 svd_unitnumber_out
      character(len=1024) svd_fname_out
      include './dir_gen_Jsvd_6/gen_Jsvd_svddecl.txt'
      logical flag_warning
      real *8 grid_p_(0:0)
      integer n_r
      integer n_delta_v,ndv
      real *8 delta_x_(0:0),delta_y_(0:0)
      real *8 R_max,K_max,delta,delta_max,n_pixel
      integer nl
      integer *4 sum_l2_f
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering get_svd_2]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
      R_max = 2.0d0*pi*grid_p_(n_r-1)
      K_max = grid_p_(n_r-1)
      delta_max = 0.0d0
      do ndv=0,n_delta_v-1
         delta = dsqrt(delta_x_(ndv)**2 + delta_y_(ndv)**2)
         if (delta.gt.delta_max) then
            delta_max = delta
         end if
      enddo !do ndv=0,n_delta_v-1
      n_pixel = delta_max/dsqrt(2.0d0)*2.0d0*K_max
      if (verbose.gt.0) then
         write(6,'(A,F8.3,A,F8.3)') 'R_max: ',R_max,'; delta_max: '
     $        ,delta_max
         write(6,'(A,F8.3,A,F8.3)') 'K_max: ',K_max,'; n_pixel: '
     $        ,n_pixel
      end if
      if (n_pixel.gt.5 .and. flag_warning) then
         write(6,'(A,F8.3,A)') 'Warning, n_pixel ',n_pixel
     $        ,' too large in get_svd_2'
      end if
      if (eps_target.lt.1.0d-6 .and. flag_warning) then
         write(6,'(A,F8.5,A)') 'Warning, eps_target ',eps_target
     $        ,' too small in get_svd_2'
      end if
      if (.false.) then ! do nothing and continue to next line ;
      include './dir_gen_Jsvd_6/gen_Jsvd_svdpick.txt'
      end if
      if (verbose.gt.0) then
         write(6,'(A,A)') 'svd_fname: ' , svd_fname
      end if !if (verbose.gt.0) then
      include './dir_gen_Jsvd_6/gen_Jsvd_svdload.txt'
      if (verbose.gt.0) then
         write(6,'(A,I0)') 'n_svd_r: ',n_svd_r
         write(6,'(A,I0)') 'n_svd_d: ',n_svd_d
         write(6,'(A,I0)') 'n_svd_l: ',n_svd_l
         write(6,'(A,I0)') 'svd_unitnumber: ',svd_unitnumber
      end if
      n_svd_r_out = n_svd_r
      n_svd_d_out = n_svd_d
      allocate(flag_s_(0:n_svd_l-1))
      n_svd_l_out = 0
      do nl=0,n_svd_l-1
         flag_s_(nl) = .false.
         if (svd_s_(nl).lt.eps_target) then
            if (verbose.gt.1) then
               write(6,'(A,I0,A,F8.4,A)') ' nl: ' , nl , ' svd_s: ' ,
     $              svd_s_(nl) , ' skipping '
            end if !if (verbose.gt.1) then
         end if !if (svd_s_(nl).lt.eps_target) then
         if (svd_s_(nl).ge.eps_target) then
            flag_s_(nl) = .true.
            if (verbose.gt.1) then
               write(6,'(A,I0,A,F8.4,A,I0)') ' nl: ' , nl , ' svd_s: ' ,
     $              svd_s_(nl) , ' retaining; n_svd_l_out: ' ,
     $              n_svd_l_out
            end if !if (verbose.gt.1) then
            svd_l_out_(n_svd_l_out) = svd_l_(nl)
            call cp1_r8(n_svd_d,svd_U_d_(nl*n_svd_d)
     $           ,svd_U_d_out_(n_svd_l_out*n_svd_d))
            call cp1_r8(1,svd_s_(nl),svd_s_out_(n_svd_l_out))
            call cp1_r8(n_svd_r,svd_V_r_(nl*n_svd_r)
     $           ,svd_V_r_out_(n_svd_l_out*n_svd_r))
            n_svd_l_out = n_svd_l_out + 1
         end if !if (svd_s_(nl).ge.eps_target) then
      enddo !do nl=0,n_svd_l-1
      call cp1_r8(n_svd_r,svd_r_,svd_r_out_)
      call cp1_r8(n_svd_d,svd_d_,svd_d_out_)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' sum(flag_s_) ' , sum_l2_f(n_svd_l,flag_s_)
         if (verbose.gt.1) then
            call print_all_l2(n_svd_l,flag_s_,' flag_s_: ')
         end if                 !if (verbose.gt.1) then
      end if !if (verbose.gt.0) then
      deallocate(flag_s_)
      if (verbose.gt.0) then
         write(6,'(A)') '[finished get_svd_2]'
      end if !if (verbose.gt.0) then
      end
