!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calculate bessel coefficients of T(R(CTF.*M)). ;\n
!> T = translation by -delta_x_est,-delta_y_est.  ;\n
!> R = rotation by -gamma_est. ;\n
!> CTF = dconjg(CTF) in k-space polar coord. ;\n
!> M = template in k-space polar coord. ;\n
!> O_T_R_CTF_M_q__(nr + n_r*(nm + n_M_tot*nw))  ;\n
!> is equal to the bessel-coefficients ;\n
!> M_q_(ic))  ;\n
!> for ic = nw + n_w_csum_(nr), ;\n
!> where M_q_ = T(R(CTF.*M)), with  ;\n
!> T <-- -delta_est ;\n
!> R <-- -gamma_est ;\n
!> and CTF <-- dconjg(CTF_(nctf)) ;\n
!> and M <-- M_p__(nm*ld_M). ;\n
      subroutine ti8_build_O_T_R_CTF_M_q_3(
     $     verbose
     $     ,delta_x_est_
     $     ,delta_y_est_
     $     ,gamma_z_est_
     $     ,ctf_ind_est_ 
     $     ,fftw_plan_frwd_
     $     ,fftw_plan_back_ 
     $     ,fftw_0in_
     $     ,fftw_out_
     $     ,n_r 
     $     ,grid_p_
     $     ,weight_p_
     $     ,n_w_
     $     ,n_w_sum 
     $     ,M_p_ 
     $     ,M_q_
     $     ,n_M_sub
     $     ,I_M_sample_
     $     ,ld_M 
     $     ,M_p__ 
     $     ,CTF_p_
     $     ,n_CTF
     $     ,ld_CTF
     $     ,CTF_p__ 
     $     ,C_M_ 
     $     ,n_M_tot
     $     ,O_T_R_CTF_M_q__
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of T(R(CTF.*M)).
c$$$      T = translation by -delta_x_est,-delta_y_est. 
c$$$      R = rotation by -gamma_est.
c$$$      CTF = dconjg(CTF) in k-space polar coord.
c$$$      M = template in k-space polar coord.
c$$$      O_T_R_CTF_M_q__(nr + n_r*(nm + n_M_tot*nw)) 
c$$$      is equal to the bessel-coefficients
c$$$      M_q_(ic)) 
c$$$      for ic = nw + n_w_csum_(nr),
c$$$      where M_q_ = T(R(CTF.*M)), with 
c$$$      T <-- -delta_est
c$$$      R <-- -gamma_est
c$$$      and CTF <-- dconjg(CTF_(nctf))
c$$$      and M <-- M_p__(nm*ld_M).
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      integer n_r,n_w_(0:n_r-1),n_w_sum,n_M_sub,I_M_sample_(0:0)
     $     ,n_M_tot,ld_M,n_CTF,ld_CTF
      real *8 delta_x_est_(0:0)
      real *8 delta_y_est_(0:0)
      real *8 gamma_z_est_(0:0)
      real *8 ctf_ind_est_(0:0)
      real *8 grid_p_(0:n_r-1)
      real *8 weight_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      integer *8 fftw_plan_back_(0:n_r-1)
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      complex *16 M_p_(0:0),M_q_(0:0),M_p__(0:0)
      complex *16 CTF_p_(0:0),CTF_q_(0:0),CTF_p__(0:0)
      complex *16 C_M_(0:0)
      complex *16 O_T_R_CTF_M_q__(0:0)
      complex *16 C_M
      complex *16 Z_tmp
      integer n_w_max,nr,nw,nm,nctf
      integer ndx,ndy,nx,ld_X
      logical flag_dconjg
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 pi

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[entering ti8_build_O_T_R_CTF_M_q_3]'
      end if !if (verbose.gt.0) then

      if (verbose.gt.1) then
         write(6,'(A,I0)') ' verbose: ' , verbose
         write(6,'(A,I0)') ' n_M_sub: ' , n_M_sub         
         call print_sub_r8(n_M_sub,delta_x_est_,' delta_x_est_: ')
         call print_sub_r8(n_M_sub,delta_y_est_,' delta_y_est_: ')
         call print_sub_r8(n_M_sub,gamma_z_est_,' gamma_z_est_: ')
         call print_sub_r8(n_M_sub,ctf_ind_est_,' ctf_ind_est_: ')
         write(6,'(A,I0)') ' n_r: ' , n_r
         call print_sub_r8(n_r,grid_p_,' grid_p_: ')
         call print_sub_r8(n_r,weight_p_,' weight_p_: ')
         write(6,'(A,I0)') ' n_w_sum: ' , n_w_sum
         call print_sub_c16(n_w_sum,M_p_,' M_p_: ')
         call print_sub_c16(n_w_sum,M_q_,' M_q_: ')
         write(6,'(A,I0)') ' n_M_sub: ' , n_M_sub
         call print_all_i4(n_M_sub,I_M_sample_,' I_M_sample_: ')
         write(6,'(A,I0)') ' ld_M: ' , ld_M
         call print_sub_c16(n_w_sum,M_p__(I_M_sample_(0)*ld_M)
     $        ,' M_p__(0): ')
         call print_sub_c16(n_w_sum,M_p__(I_M_sample_(n_M_sub-1)*ld_M)
     $        ,' M_p__(end): ')
         call print_sub_c16(n_w_sum,CTF_p_,' CTF_p_: ')
         write(6,'(A,I0)') ' n_CTF: ' , n_CTF
         write(6,'(A,I0)') ' ld_CTF: ' , ld_CTF
         call print_sub_c16(n_w_sum*n_CTF,CTF_p__,' CTF_p__: ')
         call print_sub_c16(n_M_sub,C_M_,' C_M_: ')
         write(6,'(A,I0)') ' n_M_tot: ' , n_M_tot
      end if !if (verbose.gt.1) then

      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)

      do nm=0,n_M_sub-1
         delta_x_est = +delta_x_est_(nm)
         delta_y_est = +delta_y_est_(nm)
         gamma_z_est = +gamma_z_est_(nm)
         nctf = nint(ctf_ind_est_(nm))
         if (0+verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,I0,A,I0)') 'Processing image ',nm,'/'
     $           ,n_M_sub,' associated with ctf ',nctf,'/',n_CTF
            write(6,'(A)') ' First we calculate l2_norm C_M.'
         end if
         call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $        ,M_p__(I_M_sample_(nm)*ld_M) ,M_p__(I_M_sample_(nm) *ld_M)
     $        ,C_M)
         C_M = zsqrt(C_M)/dsqrt(2.0d0)
         C_M_(nm) = C_M
         if (verbose.gt.2) then
            write(6,'(A,I0,A,2F16.3)') ' nm: ',nm,'; C_M: ',C_M
         end if
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply ctf-star to M_p_ '
         end if
         call cp1_c16(n_w_sum,M_p__(I_M_sample_(nm)*ld_M),M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call print_sub_c16(n_w_sum,M_p_,'original M_p_: ')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,M_p_,M_p_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,F8.4,1X,F8.4)') '|original M_p_|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         call xc1_c16(n_w_sum,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call print_sub_c16(n_w_sum,M_p_,'  CTF *  M_p_: ')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,M_p_,M_p_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,F8.4,1X,F8.4)') '|  CTF *  M_p_|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply R_{-gamma_est} to M_p_ '
         end if
         call rotate_p_to_p_fftw(n_r,fftw_plan_frwd_,fftw_plan_back_
     $        ,n_w_,n_w_sum,fftw_0in_,fftw_out_,M_p_,-gamma_z_est,M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call print_sub_c16(n_w_sum,M_p_,'  R(CTF*M_p_): ')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,M_p_,M_p_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,F8.4,1X,F8.4)') '|  R(CTF*M_p_)|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply T_{-delta_est} to M_p_ '
         end if
         call transf_p_to_p(n_r,grid_p_,n_w_,n_w_sum,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call print_sub_c16(n_w_sum,M_p_,'T(R(CTF*M_p_)):')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,M_p_,M_p_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,F8.4,1X,F8.4)') '|T(R(CTF*M_p_))|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we convert M_p_ to M_q_ '
         end if
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_w_sum
     $        ,fftw_0in_,fftw_out_,M_p_,M_q_)
         flag_dconjg = .false.
         nx = 0 + n_r*nm
         ld_X = n_r*n_M_tot
         call innerproduct_q_stretch_quad_0(n_r,grid_p_,weight_p_ ,n_w_
     $        ,n_w_sum,M_q_,flag_dconjg,ld_X,O_T_R_CTF_M_q__(nx))
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call print_sub_c16(n_w_sum,M_q_,'M_q_: ')
            call innerproduct_p_quad(n_r,grid_p_,weight_p_,n_w_,n_w_sum
     $           ,M_q_,M_q_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)
            write(6,'(A,F8.4,1X,F8.4)') '|M_q_|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
      enddo ! do nm=0,n_M_sub-1

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished ti8_build_O_T_R_CTF_M_q_3]'
      end if !if (verbose.gt.0) then
      end
