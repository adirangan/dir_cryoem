      subroutine test_innerproduct_fast_ZZ_2(verbose,rseed,n_delta_x
     $     ,delta_x_,n_delta_y,delta_y_,n_gamma_z,gamma_z_,n_r ,grid_p_
     $     ,n_w_,n_A,n_S,S_alpha_polar_a_,S_alpha_azimu_b_,n_M_sub
     $     ,polar_a_est_,azimu_b_est_,gamma_z_est_ ,delta_x_est_
     $     ,delta_y_est_,l2_norm_est_,ctf_ind_est_ ,S_index_est_
     $     ,polar_a_upd_,azimu_b_upd_,gamma_z_upd_ ,delta_x_upd_
     $     ,delta_y_upd_,l2_norm_upd_,ctf_ind_upd_ ,S_index_upd_
     $     ,alpha_update_f,displacement_max,C_M_,CTF_R_S_ ,n_M_tot
     $     ,S_T_T_R_CTF_M_q__,ZZ_)
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer *4 rseed !random seed. ;
      real *8 adi_rand_f
      integer n_delta_x,n_delta_y,n_gamma_z
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 S_alpha_polar_a_(0:0),S_alpha_azimu_b_(0:0)
      real *8 delta_x_(0:n_delta_x-1),delta_y_(0:n_delta_y-1)
      real *8 gamma_z_(0:n_gamma_z-1),grid_p_(0:n_r-1)
      real *8 polar_a_est_(0:n_M_sub-1),azimu_b_est_(0:n_M_sub-1)
     $     ,gamma_z_est_(0:n_M_sub-1)
      real *8 delta_x_est_(0:n_M_sub-1),delta_y_est_(0:n_M_sub-1)
     $     ,l2_norm_est_(0:n_M_sub-1)
      real *8 ctf_ind_est_(0:n_M_sub-1)
      real *8 S_index_est_(0:n_M_sub-1)
      real *8 polar_a_upd_(0:n_M_sub-1),azimu_b_upd_(0:n_M_sub-1)
     $     ,gamma_z_upd_(0:n_M_sub-1)
      real *8 delta_x_upd_(0:n_M_sub-1),delta_y_upd_(0:n_M_sub-1)
     $     ,l2_norm_upd_(0:n_M_sub-1)
      real *8 ctf_ind_upd_(0:n_M_sub-1)
      real *8 S_index_upd_(0:n_M_sub-1)
      real *8 displacement_max
      real *8 alpha_update_f
      integer n_S,n_M_sub,n_M_tot
      complex *16 C_M_(0:0),C_M
      complex *16 CTF_R_S_(0:0)
      complex *16 S_T_T_R_CTF_M_q__(0:0)
      complex *16 ZZ_(0:0)
      complex *16 ZZ
      complex *16, allocatable :: CTF_R_S_use_(:)
      complex *16, allocatable :: ZZ_sub_(:)
      complex *16 CTF_R_S_use,C_Z_use
      integer n_w_max,nr,nw,nC
      integer ndx,ndy,ngz
      integer ns,nm,nctf,ns_use,nr_use
      integer nx2,nx3
      real *8 delta_x,delta_y,gamma_z,gamma_z_est
      integer ndx_optimal,ndy_optimal,ngz_optimal
      real *8 delta_x_optimal,delta_y_optimal,gamma_z_optimal
      real *8 delta_x_tmp,delta_y_tmp,gamma_z_tmp
      complex *16 C_Z_optimal,CTF_R_S_optimal
      real *8, allocatable ::  delta_x_SM_(:)
      real *8, allocatable ::  delta_y_SM_(:)
      real *8, allocatable ::  gamma_z_SM_(:)
      complex *16, allocatable ::  C_S_SM_(:)
      complex *16, allocatable ::  C_Z_SM_(:)
      real *8, allocatable ::  delta_x_sort_SM_(:)
      real *8, allocatable ::  delta_y_sort_SM_(:)
      real *8, allocatable ::  gamma_z_sort_SM_(:)
      complex *16, allocatable ::  C_S_sort_SM_(:)
      complex *16, allocatable ::  C_Z_sort_SM_(:)
      integer, allocatable :: I_permute_SM_(:)
      integer, allocatable :: I_inverse_SM_(:)
      real *8 pi
      character(len=1024) format_string
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_fast_ZZ_2]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)

c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      allocate(C_Z_SM_(0:n_S-1))
      allocate(C_S_SM_(0:n_S-1))
      allocate(delta_x_SM_(0:n_S-1))
      allocate(delta_y_SM_(0:n_S-1))
      allocate(gamma_z_SM_(0:n_S-1))
      allocate(C_Z_sort_SM_(0:n_S-1))
      allocate(C_S_sort_SM_(0:n_S-1))
      allocate(delta_x_sort_SM_(0:n_S-1))
      allocate(delta_y_sort_SM_(0:n_S-1))
      allocate(gamma_z_sort_SM_(0:n_S-1))
      allocate(I_permute_SM_(0:n_S-1))
      allocate(I_inverse_SM_(0:n_S-1))
      allocate(CTF_R_S_use_(0:n_gamma_z-1))
      call cl1_c16(n_gamma_z,CTF_R_S_use_)
      allocate(ZZ_sub_(0:n_w_max-1))
      call cl1_c16(n_w_max,ZZ_sub_)
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do nm=0,n_M_sub-1
         C_M = C_M_(nm)
         nctf = nint(ctf_ind_est_(nm))
         gamma_z_est = gamma_z_est_(nm)
         do ns=0,n_S-1
c$$$            call cp1_c16(n_gamma_Z,CTF_R_S_(0 + ns *n_gamma_z + nctf
c$$$     $           *n_gamma_z*n_S),CTF_R_S_use_)
            call get_CTF_R_S_use_(gamma_z_est,n_gamma_z,CTF_R_S_(0 + ns
     $           *n_gamma_z + nctf*n_gamma_z*n_S),CTF_R_S_use_)
            call cl1_c16(n_delta_x*n_delta_y*n_gamma_z,ZZ_)
            do ndy=0,n_delta_y-1
               delta_y = delta_y_(ndy)
               do ndx=0,n_delta_x-1
                  delta_x = delta_x_(ndx)
                  call cl1_c16(n_w_max,ZZ_sub_)
                  do nw=0,n_w_max-1
                     nx2 = ndx + n_delta_x*(ndy + n_delta_y*(ns + n_S
     $                    *nm))
                     nx3 = nw*n_delta_x*n_delta_y*n_S*n_M_tot
                     ZZ_sub_(nw) = S_T_T_R_CTF_M_q__(nx2 + nx3)
                  enddo !do nw
                  do ngz=0,n_gamma_z-1
                     CTF_R_S_use = CTF_R_S_use_(ngz)
                     gamma_z = gamma_z_(ngz)
                     call interp1_c16(n_w_max,0.0d0,2*pi,ZZ_sub_,
     $                    +gamma_z,ZZ)
                     nC = ndx + ndy*n_delta_x + ngz*n_delta_x*n_delta_y
                     if (zabs(C_M*CTF_R_S_use).le.1.0d-15) then
                        C_Z_use = 1.0d0
                     else
                        C_Z_use = C_M*CTF_R_S_use
                     end if
                     ZZ_(nC) = ZZ/(n_r**4)/C_Z_use
                  enddo !do ngz
               enddo !do ndx
            enddo ! do ndy
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call test_innerproduct_fast_ZZ_scan_0(delta_x_est_(nm)
     $           ,delta_y_est_(nm),gamma_z_est_(nm),displacement_max
     $           ,n_delta_x ,delta_x_ ,n_delta_y,delta_y_,n_gamma_z
     $           ,gamma_z_,ZZ_ ,ndx_optimal ,ndy_optimal ,ngz_optimal
     $           ,C_Z_optimal)
            delta_x_optimal = delta_x_(ndx_optimal)
            delta_y_optimal = delta_y_(ndy_optimal)
            gamma_z_optimal = gamma_z_(ngz_optimal)
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call test_innerproduct_batch_excerpt_6(n_gamma_z,gamma_z_
     $           ,CTF_R_S_(0+ns*n_gamma_z + nctf*n_gamma_z*n_S),
     $           +gamma_z_est_(nm),ngz_optimal,CTF_R_S_optimal)
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I0,A,I2,1X,I2,1X,I2,A,2F8.3,A,2F8.3)')
     $              ' Match between M_q_(' , nm , ') and S_q_(',ns
     $              ,') at ndx_optimal,ndy_optimal,ngz_optimal: '
     $              ,ndx_optimal ,ndy_optimal ,ngz_optimal
     $              ,'; CTF_R_S_optimal: ',CTF_R_S_optimal ,
     $              '; C_Z_optimal: ' ,C_Z_optimal
            end if              !verbose
            delta_x_SM_(ns) = delta_x_optimal
            delta_y_SM_(ns) = delta_y_optimal
            gamma_z_SM_(ns) = gamma_z_optimal
            C_S_SM_(ns) = CTF_R_S_optimal
            C_Z_SM_(ns) = C_Z_optimal
         enddo !do ns=0,n_S-1
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (0+verbose.gt.1) then
            write(6,'(A)') ' Sort innerproducts.'
         end if !if (verbose.gt.1) then
         if (0+verbose.gt.1) then
            write(6,'(A)') 'delta_x_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_x_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'delta_y_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_y_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'gamma_z_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (gamma_z_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_S_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_S_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_Z_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_Z_SM_(ns),ns=0,n_S-1)
         end if !if (verbose.gt.1) then
         call test_innerproduct_batch_sort_a(n_S,1,delta_x_SM_
     $        ,delta_y_SM_,gamma_z_SM_,C_S_SM_,C_Z_SM_
     $        ,delta_x_sort_SM_ ,delta_y_sort_SM_,gamma_z_sort_SM_
     $        ,C_S_sort_SM_,C_Z_sort_SM_,I_permute_SM_ ,I_inverse_SM_)
         if (0+verbose.gt.1) then
            write(6,'(A)') 'delta_x_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_x_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'delta_y_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (delta_y_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'gamma_z_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(F8.3,1X))'
            write(6,format_string) (gamma_z_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_S_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_S_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'C_Z_sort_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S*2,'(F8.3,1X))'
            write(6,format_string) (C_Z_sort_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'I_permute_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
            write(6,format_string) (I_permute_SM_(ns),ns=0,n_S-1)
            write(6,'(A)') 'I_inverse_SM_^T: '
            write(format_string,'(A,I0,A)') '(',n_S,'(I8,1X))'
            write(6,format_string) (I_inverse_SM_(ns),ns=0,n_S-1)
         end if !if (verbose.gt.1) then
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (verbose.gt.1) then
            write(6,'(A)') ' Update image parameters.'
         end if !if (verbose.gt.1) then
         nr_use = max(0,min(n_S-1,floor((1.0d0-adi_rand_f(rseed)
     $        *alpha_update_f)*n_S)))
         ns_use = I_permute_SM_(nr_use)
         call get_interchange_delta(delta_x_est_(nm) ,delta_y_est_(nm)
     $        ,gamma_z_SM_(ns_use),delta_x_tmp ,delta_y_tmp)
         polar_a_upd_(nm) = S_alpha_polar_a_(ns_use)
         azimu_b_upd_(nm) = S_alpha_azimu_b_(ns_use)
         gamma_z_upd_(nm) = gamma_z_SM_(ns_use) + gamma_z_est_(nm)
         delta_x_upd_(nm) = delta_x_SM_(ns_use) + delta_x_tmp
         delta_y_upd_(nm) = delta_y_SM_(ns_use) + delta_y_tmp
         S_index_upd_(nm) = ns_use
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,2F8.3,A,2F8.3,A,2F8.3)') 'nm ' ,nm
     $           ,'; ns_use ',ns_use,'; C_S ' ,C_S_SM_(ns_use),'; C_M '
     $           ,C_M_(nm) ,'; C_Z ', C_Z_SM_(ns_use)
         end if !if (verbose.gt.1) then
         if (zabs(C_S_SM_(ns_use)).le.1.0d-15) then
            l2_norm_upd_(nm) = 1.0d0
         else 
            l2_norm_upd_(nm) = zabs(C_M_(nm)) *zabs(C_Z_SM_(ns_use)) /
     $           zabs(C_S_SM_(ns_use))
         end if ! if (zabs(C_S_SM_(ns_use)).le.1.0d-15) then
c$$$            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      enddo !do nm=0,n_M_sub-1

      deallocate(C_Z_sort_SM_)
      deallocate(C_S_sort_SM_)
      deallocate(delta_x_sort_SM_)
      deallocate(delta_y_sort_SM_)
      deallocate(gamma_z_sort_SM_)
      deallocate(C_Z_SM_)
      deallocate(C_S_SM_)
      deallocate(delta_x_SM_)
      deallocate(delta_y_SM_)
      deallocate(gamma_z_SM_)
      deallocate(ZZ_sub_)
      deallocate(CTF_R_S_use_)

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_fast_ZZ_2]'
      end if !if (verbose.gt.0) then

      end
