      subroutine test_innerproduct_5_Z_S_q_3(verbose,n_svd_r,svd_r_
     $     ,n_svd_l,svd_l_,svd_s_,svd_V_r_,fftw_plan_frwd_,fftw_in1_
     $     ,fftw_out_ ,n_r,grid_p_,n_w_,n_A,S_p_,S_q_,n_S_sub
     $     ,I_S_sample_,ld_S,S_p__ ,n_S_tot,Z_S_q__)
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of Z(S).
c$$$      Z = svd of translation operator. 
c$$$      S = template in k-space polar coord.
c$$$      Z_S_q__(nr + n_r*(nl + n_svd_l*(ns + n_S_tot*nw))) 
c$$$      is equal to the bessel-coefficients
c$$$      conjg(S_q_(ic)) 
c$$$      for ic = nw + n_w_sum_(nr),
c$$$      where S_q_ = Z(S), with 
c$$$      Z representing the left-hand factor 
c$$$      of the translation-operator
c$$$      associated with delta,
c$$$      and S <-- S_p__(ns*ld_S).
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      logical warning_flag
      data warning_flag / .true. /
      integer n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      integer n_r,n_w_(0:n_r-1),n_A,n_S_sub,I_S_sample_(0:0)
     $     ,n_S_tot,ld_S
      real *8 grid_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      complex *16 fftw_in1_(0:0),fftw_out_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 Z_S_q__(0:0)
      complex *16 Z_tmp,C_q
      integer n_w_max,ic_store,ic,ict,nr,nw,nwt,nw_fix,nw_C,ns
      integer nx,ld_X
      real *8 R_q,R_pos,R_pre,dr,dw,dA,dAn,dsqrt_dAn
      real *8 svd_r_max,n_svd_r_v
      integer n_svd_r_pre,n_svd_r_pos
      real *8 d_svd_r_pre,d_svd_r_pos
      real *8 alpha_r,beta_r
      real *8 D_V_r,D_s
      integer nl,I_l
      real *8 pi
      integer, allocatable :: n_X_(:)
      if (verbose.gt.0) then
         write(6,'(A)') '[entering test_innerproduct_5_Z_S_q_3]'
      end if !if (verbose.gt.0) then
      pi = 4.0*atan(1.0)
      svd_r_max = svd_r_(n_svd_r-1)
      R_q = 2*pi*grid_p_(n_r-1)
      if (R_q.gt.svd_r_max .and. warning_flag) then
         write(6,'(A,F6.3,A,F6.3,A)') 'Warning, 2*pi*r ',R_q,'>'
     $        ,svd_r_max,' in innerproduct_q_k_svdr'
      end if
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      if (verbose.gt.2) then
         allocate(n_X_(0:n_r*n_svd_l*n_w_max-1))
         call cl1_i4(n_r*n_svd_l*n_w_max,n_X_)
      end if !if (verbose.gt.2) then

      ld_X = n_r*n_svd_l*n_S_tot
      do ns=0,n_S_sub-1
         if (verbose.gt.2) then
            call innerproduct_p(n_r,grid_p_,n_w_,n_A
     $           ,S_p__(I_S_sample_(ns)*ld_S)
     $           ,S_p__(I_S_sample_(ns)*ld_S),Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_p_(',ns,')|: ' , Z_tmp
            call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_)
            call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $           ,fftw_in1_,fftw_out_,S_p_,S_q_)
            call write_sub_c16(n_A,S_q_,5,'S_q_: ')
            call innerproduct_p(n_r,grid_p_,n_w_,n_A,S_q_,S_q_,Z_tmp)
            Z_tmp = zsqrt(Z_tmp)/(n_r*n_r)
            write(6,'(A,I0,A,F8.4,1X,F8.4)') '|S_q_(',ns,')|: ' , Z_tmp
         end if                 !if (verbose.gt.2) then
         call cp1_c16(n_A,S_p__(I_S_sample_(ns)*ld_S),S_p_) 
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_A
     $        ,fftw_in1_,fftw_out_,S_p_,S_q_)         
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
                  write(6,'(A,I3,1X,I3,1X,F16.3,1X,F16.3)')
     $                 ' % % l I_l D_V_r D_s: ',nl,I_l,D_V_r,D_s
               end if !if (verbose.gt.2 .and. nr.lt.5) then
               nx = 0 + n_r*(nl + n_svd_l*ns);
               ic = ic_store
               do nw=0,n_w_(nr)-1
                  call periodize_i(nw+I_l,0,n_w_(nr),nwt)
                  ict = ic-nw+nwt
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(6(A,I0))') ' % % % nw:' , nw , '; I_l:' ,
     $                    I_l , '; n_w_(nr):' , n_w_(nr) , '; nwt:' ,
     $                    nwt , '; ic:' , ic , '; ict:' , ict
                  end if !if (verbose.gt.3 .and. nr.lt.5) then
                  if (nw.gt.n_w_(nr)/2) then
                     nw_fix = nw - n_w_(nr) + n_w_max
                     if (verbose.gt.3 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (full loop)'
                     end if
                     C_q = conjg(D_s*D_V_r*S_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_S_q__(nx + nr + ld_X*nw_C) = Z_S_q__(nx + nr +
     $                    ld_X*nw_C) +C_q*dsqrt_dAn
                  else if (nw.eq.n_w_(nr)/2) then
                     nw_fix = nw
                     if (verbose.gt.3 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (half orig)'
                     end if
                     C_q = dsqrt(0.5d0)*conjg(D_s*D_V_r*S_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_S_q__(nx + nr + ld_X*nw_C) = Z_S_q__(nx + nr +
     $                    ld_X*nw_C) +C_q*dsqrt_dAn
                     nw_fix = nw - n_w_(nr) + n_w_max
                     if (verbose.gt.3 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (half loop)'
                     end if
                     C_q = dsqrt(0.5d0)*conjg(D_s*D_V_r*S_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_S_q__(nx + nr + ld_X*nw_C) = Z_S_q__(nx + nr +
     $                    ld_X*nw_C) +C_q*dsqrt_dAn
                  else
                     nw_fix = nw
                     if (verbose.gt.3 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (full orig)'
                     end if
                     C_q = conjg(D_s*D_V_r*S_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_S_q__(nx + nr + ld_X*nw_C) = Z_S_q__(nx + nr +
     $                    ld_X*nw_C) +C_q*dsqrt_dAn
                  end if !if nw
                  ic = ic + 1
               enddo !do nw=0,n_w_(nr)-1
            enddo !do nl=0,n_svd_l-1
         enddo !do nr=0,n_r-1
         if (verbose.gt.2) then
            write(6,'(A,I0)') ' ns: ' , ns
            call write_all_i4__(n_r,n_svd_l*n_w_max,n_X_,6,'n_X_: ')
         end if !if (verbose.gt.2) then
      enddo !do ns=0,n_S_sub-1
      if (verbose.gt.2) then
         deallocate(n_X_)
      end if !if (verbose.gt.2) then

      if (verbose.gt.0) then
         write(6,'(A)') '[finished test_innerproduct_5_Z_S_q_3]'
      end if !if (verbose.gt.0) then
      end
