!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calculate bessel coefficients of T(S). ;\n
!> T = translation by delta_x,delta_y. ;\n
!> S = template in k-space polar coord. ;\n
!> T_S_q__(nr + n_r*(ndv + n_delta_v*(ns + n_S_tot*nw))) ;\n
!> is equal to the bessel-coefficients ;\n
!> dconjg(S_q_(ic)) ;\n
!> for ic = nw + n_w_csum_(nr), ;\n
!> where S_q_ = T(S), with ;\n
!> T <-- delta_x_(ndv) and delta_y_(ndv) ;\n
!> and S <-- S_p__(ns*ld_S). ;\n
      subroutine ti8_build_T_S_q_4(
     $     verbose
     $     ,n_delta_v
     $     ,delta_x_
     $     ,delta_y_
     $     ,fftw_plan_frwd_
     $     ,fftw_0in_
     $     ,fftw_out_ 
     $     ,n_r
     $     ,grid_p_
     $     ,weight_p_
     $     ,n_w_
     $     ,n_w_sum
     $     ,S_p_
     $     ,S_q_
     $     ,n_S_sub
     $     ,I_S_sample_
     $     ,ld_S
     $     ,S_p__ 
     $     ,n_S_tot 
     $     ,T_S_q__
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of T(S).
c$$$      T = translation by delta_x,delta_y.
c$$$      S = template in k-space polar coord.
c$$$      T_S_q__(nr + n_r*(ndv + n_delta_v*(ns + n_S_tot*nw)))
c$$$      is equal to the bessel-coefficients
c$$$      dconjg(S_q_(ic))
c$$$      for ic = nw + n_w_csum_(nr),
c$$$      where S_q_ = T(S), with
c$$$      T <-- delta_x_(ndv) and delta_y_(ndv)
c$$$      and S <-- S_p__(ns*ld_S).
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_delta_v
      integer n_r,n_w_(0:n_r-1),n_w_sum,n_S_sub,I_S_sample_(0:0)
     $     ,n_S_tot,ld_S
      real *8 delta_x_(0:n_delta_v-1)
      real *8 delta_y_(0:n_delta_v-1)
      real *8 grid_p_(0:n_r-1)
      real *8 weight_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 T_S_q__(0:0)
      complex *16 Z_tmp
      integer n_w_max,nr,nw,ns
      integer ndv,nx,ld_X
      logical flag_dconjg
      real *8 delta_x,delta_y
      real *8 pi
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_build_T_S_q_4]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)

      do ns=0,n_S_sub-1
         if (verbose.gt.2) then
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,S_p__(I_S_sample_(ns)*ld_S) ,S_p__(I_S_sample_(ns)
     $           *ld_S),Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_p_(',ns,')|: ' , Z_tmp
            call cp1_c16(n_w_sum,S_p__(I_S_sample_(ns)*ld_S),S_p_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_w_sum
     $           ,fftw_0in_,fftw_out_,S_p_,S_q_)
            call print_sub_c16(n_w_sum,S_q_,'S_q_: ')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,S_q_,S_q_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_q_(',ns,')|: ' , Z_tmp
         end if                 !if (verbose.gt.2) then
         do ndv=0,n_delta_v-1
            delta_x = delta_x_(ndv)
            delta_y = delta_y_(ndv)
            call cp1_c16(n_w_sum,S_p__(I_S_sample_(ns)*ld_S),S_p_)
            call transf_p_to_p(n_r,grid_p_,n_w_
     $           ,n_w_sum,S_p_,+delta_x,+delta_y,S_p_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_
     $           ,n_w_sum,fftw_0in_,fftw_out_,S_p_,S_q_)
            flag_dconjg = .true.
            nx = 0 + n_r*(ndv + n_delta_v*ns)
            ld_X = n_r*n_delta_v*n_S_tot
            call innerproduct_q_stretch_quad_0(n_r,grid_p_,weight_p_
     $           ,n_w_,n_w_sum,S_q_,flag_dconjg,ld_X,T_S_q__(nx))
            if (verbose.gt.2) then
               write(6,'(A,I0,A,I0)') ' ns: ' , ns , ' ndv: ' , ndv
               call print_sub_c16(n_w_sum,S_q_,'S_q_: ')
               call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_
     $              ,n_w_sum,S_q_,S_q_,Z_tmp)
               Z_tmp = zsqrt(Z_tmp)
               write(6,'(A,F8.4,1X,F8.4)') '|S_q_|: ' , Z_tmp
            end if              ! if (verbose.gt.2) then
         enddo                  !do ndv=0,n_delta_v-1
      enddo                     ! do ns=0,n_S_sub-1

      if (verbose.gt.0) then
         write(6,'(A)') '[finished ti8_build_T_S_q_4]'
      end if !if (verbose.gt.0) then
      end
