!> Doxygen comment: ;\n
!> This function updates a stack of image-parameters by matching templates to images. ;\n
!> (This uses a parameter 'f', stored as alpha_update_f, as described in leslie's frequency marching paper). ;\n
      subroutine test_alpha_update_SM_3(verbose,rseed,n_SM_max,n_SM
     $     ,alpha_SM_,flag_RTRT_vs_RTTR,alpha_est_,alpha_update_f
     $     ,alpha_upd_)
      implicit none
      include 'excerpt_define_nalpha.f'
      integer *4 verbose
      integer *4 rseed
      integer *4 n_SM_max,n_SM
      real *8 alpha_SM_(0:n_alpha*n_SM_max-1)
      logical flag_RTRT_vs_RTTR
      real *8 alpha_est_(0:n_alpha-1)
      real *8 alpha_update_f
      real *8 alpha_upd_(0:n_alpha-1)
      real *8 adi_rand_f
      integer *4 ns_use
      real *8 delta_x_upd,delta_y_upd,gamma_z_upd
      real *8 delta_x_est,delta_y_est,gamma_z_est

      if (verbose.gt.1) then
         write(6,'(A)') ' [entering test_alpha_update_SM_3].'
      end if !if (verbose.gt.1) then

      if (verbose.gt.2) then
         write(6,'(A,I0)') ' verbose: ' , verbose
         write(6,'(A,I0)') ' rseed: ' , rseed
         write(6,'(A,I0)') ' n_SM_max: ' , n_SM_max
         write(6,'(A,I0)') ' n_SM: ' , n_SM
         call print_sub_r8(n_alpha*n_SM_max,alpha_SM_,12,' alpha_SM_: ')
         write(6,'(A,L1)') ' flag_RTRT_vs_RTTR: ' , flag_RTRT_vs_RTTR
         call print_all_r8(n_alpha,alpha_est_,13,' alpha_est_: ')
         call print_all_r8(n_alpha,alpha_upd_,13,' alpha_upd_: ')
      end if !if (verbose.gt.2) then
      
      if (n_SM.lt.n_SM_max) then
         call alpha_SM_sort_0(n_SM_max,n_SM,alpha_SM_)
      end if !if (n_SM.lt.n_SM_max) then

c$$$      In this particular context, alpha_update_f is used to 
c$$$      choose an entry from 0,n_SM-1
      ns_use = max(0,min(n_SM-1,floor((1.0d0-adi_rand_f(rseed)
     $     *alpha_update_f)*n_SM)))
      delta_x_upd = alpha_SM_(nalpha_delta_x + n_alpha*ns_use)
      delta_y_upd = alpha_SM_(nalpha_delta_y + n_alpha*ns_use)
      gamma_z_upd = alpha_SM_(nalpha_gamma_z + n_alpha*ns_use)
      delta_x_est = alpha_est_(nalpha_delta_x)
      delta_y_est = alpha_est_(nalpha_delta_y)
      gamma_z_est = alpha_est_(nalpha_gamma_z)
      call get_interchange_delta_RTRT_vs_RTTR(flag_RTRT_vs_RTTR
     $     ,delta_x_est,delta_y_est,gamma_z_est,delta_x_upd,delta_y_upd
     $     ,gamma_z_upd)
      call cp1_r8(n_alpha,alpha_SM_(n_alpha*ns_use),alpha_upd_)
      alpha_upd_(nalpha_delta_x) = delta_x_upd
      alpha_upd_(nalpha_delta_y) = delta_y_upd
      alpha_upd_(nalpha_gamma_z) = gamma_z_upd

      if (verbose.gt.1) then
         write(6,'(A)') ' [finished test_alpha_update_SM_3].'
      end if !if (verbose.gt.1) then
      end
      
