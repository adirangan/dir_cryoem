      subroutine test_innerproduct_5_Z_T_R_CTF_M_q_1(verbose,n_svd_r
     $     ,svd_r_,n_svd_l,svd_l_,svd_s_,svd_V_r_,delta_x_est_
     $     ,delta_y_est_,gamma_z_est_,ctf_ind_est_ ,fftw_plan_frwd_
     $     ,fftw_in1_,fftw_out_,n_r,grid_p_,n_w_,n_A ,M_p_ ,M_q_,n_M_sub
     $     ,I_M_sample_,ld_M ,M_p__,CTF_p_ ,n_CTF,ld_CTF,CTF_p__
     $     ,C_M_ ,n_M_tot,Z_T_R_CTF_M_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of Z(T(R(CTF.*M))).
c$$$      Z = svd of translation operator. 
c$$$      T = translation by -delta_x_est,-delta_y_est. 
c$$$      R = rotation by -gamma_est.
c$$$      CTF = conjg(CTF) in k-space polar coord.
c$$$      M = template in k-space polar coord.
c$$$      Z_T_R_CTF_M_q__(nr + n_r*(nl + n_svd_l*(nm + n_M_tot*nw))) 
c$$$      is equal to the bessel-coefficients
c$$$      M_q_(ic)) 
c$$$      for ic = nw + n_w_sum_(nr),
c$$$      where M_q_ = Z(T(R(CTF.*M))), with 
c$$$      Z representing the left-hand factor 
c$$$      of the translation-operator
c$$$      associated with delta,
c$$$      T <-- -delta_est
c$$$      R <-- -gamma_est
c$$$      and CTF <-- conjg(CTF_(nctf))
c$$$      and M <-- M_p__(nm*ld_M).
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      logical warning_flag
      data warning_flag / .true. /
      integer n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      integer n_r,n_w_(0:n_r-1),n_A,n_M_sub,I_M_sample_(0:0)
     $     ,n_M_tot,ld_M,n_CTF,ld_CTF
      real *8 delta_x_est_(0:0)
      real *8 delta_y_est_(0:0)
      real *8 gamma_z_est_(0:0)
      real *8 ctf_ind_est_(0:0)
      real *8 grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 Z_q
      complex *16 M_p_(0:0),M_q_(0:0),M_p__(0:0)
      complex *16 CTF_p_(0:0),CTF_q_(0:0),CTF_p__(0:0)
      complex *16 C_M_(0:0)
      complex *16 Z_T_R_CTF_M_q__(0:0)
      complex *16 C_M
      complex *16 Z_tmp,C_q
      integer n_w_max,ic_store,ic,ict,icr,nr,nw,nwt,nwr,nw_fix,nw_C
      integer nm,nctf
      integer nx,ld_X
      real *8 R_q,R_pos,R_pre,dr,dw,dA,dAn,dsqrt_dAn
      real *8 svd_r_max,n_svd_r_v
      integer n_svd_r_pre,n_svd_r_pos
      real *8 d_svd_r_pre,d_svd_r_pos
      real *8 alpha_r,beta_r
      real *8 D_V_r,D_s
      integer nl,I_l
      real *8 delta_x_est,delta_y_est,gamma_z_est
      real *8 pi
      integer, allocatable :: n_X_(:)
      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[entering test_innerproduct_5_Z_T_R_CTF_M_q_1]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      svd_r_max = svd_r_(n_svd_r-1)
      R_q = 2*pi*grid_p_(n_r-1)
      if (R_q.gt.svd_r_max .and. warning_flag) then
         write(6,'(A,F6.3,A,F6.3,A)') 'Warning, 2*pi*r ',R_q,'>'
     $        ,svd_r_max,' in innerproduct_q_k_svdr'
      end if
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      if (verbose.gt.2) then
         allocate(n_X_(0:n_r*n_svd_l*n_w_max-1))
         call cl1_i4(n_r*n_svd_l*n_w_max,n_X_)
      end if !if (verbose.gt.2) then

      ld_X = n_r*n_svd_l*n_M_tot
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
         call innerproduct_p(n_r,grid_p_,n_w_,n_A
     $        ,M_p__(I_M_sample_(nm)*ld_M)
     $        ,M_p__(I_M_sample_(nm)*ld_M),C_M)
         C_M = zsqrt(C_M)/(n_r*n_r)
         C_M_(nm) = C_M
         if (verbose.gt.2) then
            write(6,'(A,I0,A,2F16.3)') ' nm: ',nm,'; C_M: ',C_M
         end if
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply ctf-star to M_p_ '
         end if
         call cp1_c16(n_A,M_p__(I_M_sample_(nm)*ld_M),M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call write_sub_c16(n_A,M_q_,14,'original M_p_: ')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_p_,M_p_
     $           ,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,F8.4,1X,F8.4)') '|original M_p_|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         call xc1_c16(n_A,M_p_,CTF_p__(nctf*ld_CTF),M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call write_sub_c16(n_A,M_q_,14,'  CTF *  M_p_: ')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_p_,M_p_
     $           ,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,F8.4,1X,F8.4)') '|  CTF *  M_p_|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply R_{-gamma_est} to M_p_ '
         end if
         call rotate_p2p_fz(n_r,n_w_,n_A,M_p_,-gamma_z_est,M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call write_sub_c16(n_A,M_q_,14,'  R(CTF*M_p_): ')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_p_,M_p_
     $           ,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,F8.4,1X,F8.4)') '|  R(CTF*M_p_)|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we apply T_{-delta_est} to M_p_ '
         end if
         call transf_p_to_p(n_r,grid_p_,n_w_,n_A,M_p_,-delta_x_est,
     $        -delta_y_est,M_p_)
         if (verbose.gt.2) then
            write(6,'(A,I0,A,I0,A,I0)') ' nm: ' , nm 
            call write_sub_c16(n_A,M_q_,14,'T(R(CTF*M_p_)):')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,M_p_,M_p_
     $           ,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,F8.4,1X,F8.4)') '|T(R(CTF*M_p_))|: ' , Z_tmp
         end if                 ! if (verbose.gt.2) then
         if (verbose.gt.2) then
            write(6,'(A)') ' Now we convert M_p_ to M_q_ '
         end if
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_,n_A,fftw_in1_
     $        ,fftw_out_,M_p_,M_q_)
         if (verbose.gt.2) then
            call cl1_i4(n_r*n_svd_l*n_w_max,n_X_)
         end if !if (verbose.gt.2) then
         ic = 0
         do nr=0,n_r-1
            R_q = 2*pi*grid_p_(nr)
            n_svd_r_v = (n_svd_r-1)*R_q/svd_r_max
            n_svd_r_pre = max(0,min(n_svd_r-1,floor(n_svd_r_v)))
            n_svd_r_pos = max(0,min(n_svd_r-1,ceiling(n_svd_r_v)))
            d_svd_r_pre = abs(n_svd_r_v - n_svd_r_pre)
            d_svd_r_pos = abs(n_svd_r_pos - n_svd_r_v)
            if (d_svd_r_pre+d_svd_r_pos.le.0.0d0) then
               alpha_r = 0.0d0
               beta_r = 1.0d0
            else
               alpha_r = d_svd_r_pre / (d_svd_r_pre + d_svd_r_pos)
               beta_r = d_svd_r_pos / (d_svd_r_pre + d_svd_r_pos)
            end if !if (d_svd_r_pre+d_svd_r_pos.le.0.0d0) then
            if (verbose.gt.1) then
               write(6,'(A,I0,A,F6.3,A,I0,1X,I0)') ' % nr ',nr,'; R_q '
     $              ,R_q,'; n_svd_r_pre n_svd_r_pos: ',n_svd_r_pre
     $              ,n_svd_r_pos
               write(6,'(A,F6.3,1X,F6.3)')
     $              ' % % d_svd_r_pre d_svd_r_pos: ',d_svd_r_pre
     $              ,d_svd_r_pos
               write(6,'(A,F6.3,1X,F6.3)') ' % % alpha_r beta_r: '
     $              ,alpha_r,beta_r
            end if !if (verbose.gt.1) then
            if (nr.gt.0) then
               R_pre = 0.5*(grid_p_(nr-1) + grid_p_(nr))
            else
               R_pre = grid_p_(0)
            end if !if (nr.gt.0) then
            if (nr.lt.n_r-1) then
               R_pos = 0.5*(grid_p_(nr+1) + grid_p_(nr))
            else 
               R_pos = grid_p_(n_r-1)
            end if !if (nr.lt.n_r-1) then
            dr = R_pos - R_pre
c$$$  We set the zero-mode to zero
            if (grid_p_(nr).le.0.0d0) then
               dr = 0.0d0
            end if !if (grid_p_(nr).le.0.0d0) then
            if (verbose.gt.1) then
               write(6,'(A,I0,A,I0,A,F6.3,1X,F6.3,1X,F6.3)') ' % nr ',nr
     $              ,'; n_w_(nr) ',n_w_(nr),'; R_pre,R_pos,dr: ',R_pre
     $              ,R_pos,dr
            end if !if (verbose.gt.1) then
            dw = 2*pi/(1.0d0*max(1,n_w_(nr)))
            dA = (R_pre*dr + (dr**2)/2)*dw
c$$$  We assume that the fourier basis is orthonormal (not merely orthogonal)
            dAn = dA
            dsqrt_dAn = dsqrt(dAn)
            if (verbose.gt.1) then
               write(6,'(A,I0,A,F6.3,1X,F6.3,1X,F6.3)') ' % nr ',nr
     $              ,'; dr dw dA: ',dr,dw,dA
            end if !if (verbose.gt.1) then
            ic_store = ic
            do nl=0,n_svd_l-1
               D_V_r = beta_r*svd_V_r_(n_svd_r_pre + nl*n_svd_r) +
     $              alpha_r*svd_V_r_(n_svd_r_pos + nl*n_svd_r)
               D_s = svd_s_(nl)
               I_l = svd_l_(nl)
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(2(A,I3),1X,2(A,F16.3))') ' % % nl:' , nl ,
     $                 '; I_l:' , I_l , '; D_V_r:' , D_V_r , '; D_s:' ,
     $                 D_s
               end if !if (verbose.gt.2 .and. nr.lt.5) then
               nx = 0 + n_r*(nl + n_svd_l*nm);
               ic = ic_store
               do nw=0,n_w_(nr)-1
                  call periodize_i(nw+I_l,0,n_w_(nr),nwt)
                  call periodize_i(nw-I_l,0,n_w_(nr),nwr)
                  ict = ic-nw+nwt
                  icr = ic-nw+nwr
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(10(A,I0))') ' % % % nl:' , nl , '; I_l:'
     $                    , I_l , '; ic_store:' , ic_store,
     $                    '; n_w_(nr):' , n_w_(nr) , '; nw:', nw ,
     $                    '; nwt:' ,nwt , '; nwr:' , nwr ,'; ic:' , ic ,
     $                    '; ict:' , ict , '; icr:', icr
                  end if !if (verbose.gt.3 .and. nr.lt.5) then
                  if (nw.gt.n_w_(nr)/2) then
                     nw_fix = nw - n_w_(nr) + n_w_max
                     if (verbose.gt.4 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (full loop)'
                     end if
c$$$                     no conjugate, since we are acting on M and not S
                     C_q = (D_s*D_V_r*M_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_T_R_CTF_M_q__(nx + nr + ld_X*nw_C) =
     $                    Z_T_R_CTF_M_q__(nx + nr +ld_X*nw_C) +C_q
     $                    *dsqrt_dAn
                  else if (nw.eq.n_w_(nr)/2) then
                     nw_fix = nw
                     if (verbose.gt.4 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (half orig)'
                     end if
c$$$                     no conjugate, since we are acting on M and not S
                     C_q = dsqrt(0.5d0)*(D_s*D_V_r*M_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_T_R_CTF_M_q__(nx + nr + ld_X*nw_C) =
     $                    Z_T_R_CTF_M_q__(nx + nr +ld_X*nw_C) +C_q
     $                    *dsqrt_dAn
                     nw_fix = nw - n_w_(nr) + n_w_max
                     if (verbose.gt.4 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (half loop)'
                     end if
c$$$                     no conjugate, since we are acting on M and not S
                     C_q = dsqrt(0.5d0)*(D_s*D_V_r*M_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_T_R_CTF_M_q__(nx + nr + ld_X*nw_C) =
     $                    Z_T_R_CTF_M_q__(nx + nr +ld_X*nw_C) +C_q
     $                    *dsqrt_dAn
                  else
                     nw_fix = nw
                     if (verbose.gt.4 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (full orig)'
                     end if
c$$$                     no conjugate, since we are acting on M and not S
                     C_q = (D_s*D_V_r*M_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_T_R_CTF_M_q__(nx + nr + ld_X*nw_C) =
     $                    Z_T_R_CTF_M_q__(nx + nr +ld_X*nw_C) +C_q
     $                    *dsqrt_dAn
                  end if !if nw
                  ic = ic + 1
               enddo !do nw=0,n_w_(nr)-1
            enddo !do nl=0,n_svd_l-1
         enddo !do nr=0,n_r-1
         if (verbose.gt.2) then
            write(6,'(A,I0)') ' nm: ' , nm
            call write_all_i4__(n_r,n_svd_l*n_w_max,n_X_,6,'n_X_: ')
         end if !if (verbose.gt.2) then
      enddo ! do nm=0,n_M_sub-1
      if (verbose.gt.2) then
         deallocate(n_X_)
      end if !if (verbose.gt.2) then

      if (verbose.gt.0) then
         write(6,'(A)')
     $        '[finished test_innerproduct_5_Z_T_R_CTF_M_q_1]'
      end if !if (verbose.gt.0) then
      end
