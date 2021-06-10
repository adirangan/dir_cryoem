!> Doxygen comment: ;\n
!> test_innerproduct_8 excerpt: ;\n
!> Calculate bessel coefficients of Z(S). ; ;\n
!> Z = svd of translation operator.  ; ;\n
!> S = template in k-space polar coord. ; ;\n
!> Z_S_q__(nr + n_r*(nl + n_svd_l*(ns + n_S_tot*nw)))  ; ;\n
!> is equal to the bessel-coefficients ; ;\n
!> dconjg(S_q_(ic))  ; ;\n
!> for ic = nw + n_w_csum_(nr), ; ;\n
!> where S_q_ = Z(S), with  ; ;\n
!> Z representing the left-hand factor  ; ;\n
!> of the translation-operator ; ;\n
!> associated with delta, ; ;\n
!> and S <-- S_p__(ns*ld_S). ; ;\n
!> Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ; ;\n
      subroutine ti8_build_Z_S_q_quad_5(
     $     verbose
     $     ,svd_r_max
     $     ,n_svd_r 
     $     ,svd_r_
     $     ,n_svd_l
     $     ,svd_l_
     $     ,svd_s_
     $     ,svd_V_r_
     $     ,svd_polyval_V_r_
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
     $     ,Z_S_q__
     $     )
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
c$$$      Calculate bessel coefficients of Z(S). ;
c$$$      Z = svd of translation operator.  ;
c$$$      S = template in k-space polar coord. ;
c$$$      Z_S_q__(nr + n_r*(nl + n_svd_l*(ns + n_S_tot*nw)))  ;
c$$$      is equal to the bessel-coefficients ;
c$$$      dconjg(S_q_(ic))  ;
c$$$      for ic = nw + n_w_csum_(nr), ;
c$$$      where S_q_ = Z(S), with  ;
c$$$      Z representing the left-hand factor  ;
c$$$      of the translation-operator ;
c$$$      associated with delta, ;
c$$$      and S <-- S_p__(ns*ld_S). ;
c$$$      Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ;
c$$$      %%%%%%%%%%%%%%%%%%%%%%%
      implicit none
      include '/usr/include/fftw3.f'
      integer verbose
      logical warning_flag
      data warning_flag / .true. /
      integer n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      real *8 svd_polyval_V_r_(0:n_svd_l*n_r-1)
      integer n_r,n_w_(0:n_r-1),n_w_sum,n_S_sub,I_S_sample_(0:0)
     $     ,n_S_tot,ld_S
      real *8 grid_p_(0:n_r-1)
      real *8 weight_p_(0:n_r-1)
      integer *8 fftw_plan_frwd_(0:n_r-1)
      complex *16 fftw_0in_(0:0),fftw_out_(0:0)
      complex *16 Z_q
      complex *16 S_p_(0:0),S_q_(0:0),S_p__(0:0)
      complex *16 Z_S_q__(0:0)
      complex *16 Z_tmp,C_q
      integer n_w_max,ic_store,ic,ict,icr,nr,nw,nwt,nwr,nw_fix,nw_C,ns
      integer n_w_t ! positive threshold for overflow. ;
      integer nwc ! centered nw. ;
      integer nwd ! displaced nw. ;
      logical flag_ic0_overflow ! notes whether or not M_q_ coefficient should be set to 0. ;
      logical flag_ict_overflow ! notes whether or not M_q_ coefficient should be set to 0. ;
      logical flag_icr_overflow ! notes whether or not M_q_ coefficient should be set to 0. ;
      integer nx,ld_X
      real *8 dw,dA,dAn,dsqrt_dAn
      real *8 svd_r_max,svd_r_m,svd_r_c,svd_r(0:0)
      real *8 D_V_r,D_s
      integer nl,I_l
      real *8 pi
      integer, allocatable :: n_X_(:)
c$$$      real *8, allocatable :: V_r_(:)
      if (verbose.gt.0) then
         write(6,'(A)') '[entering ti8_build_Z_S_q_quad_5]'
      end if !if (verbose.gt.0) then
      pi = 4.0d0*datan(1.0d0)
c$$$      allocate(V_r_(0:n_svd_l*n_r-1))
      svd_r_m = svd_r_max / 2.0d0
      svd_r_c = svd_r_m
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.1) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      if (verbose.gt.2) then
         allocate(n_X_(0:n_r*n_svd_l*n_w_max-1))
         call cl1_i4(n_r*n_svd_l*n_w_max,n_X_)
      end if !if (verbose.gt.2) then
      do nr=0,n_r-1
         if (grid_p_(nr).gt.svd_r_max .and. warning_flag) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)')
     $           'Warning, grid_p_(nr) ',grid_p_(nr),'>',svd_r_max
     $           ,'; ratio = ',grid_p_(nr)/svd_r_max
     $           ,' in ti8_build_Z_S_q_quad_5.f'
         end if
         svd_r(0) = (grid_p_(nr) - svd_r_m)/svd_r_c
c$$$         do nl=0,n_svd_l-1
c$$$            call polyval_r8_reverse_0(n_svd_r,svd_V_r_(0+nl*n_svd_r),1
c$$$     $           ,svd_r(0),V_r_(nl+nr*n_svd_l))
c$$$         enddo !do nl=0,n_svd_l-1         
      enddo !do nr=0,n_r-1
      ld_X = n_r*n_svd_l*n_S_tot
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
         call cp1_c16(n_w_sum,S_p__(I_S_sample_(ns)*ld_S),S_p_) 
         call interp_p_to_q_fftw(n_r,fftw_plan_frwd_,n_w_ ,n_w_sum
     $        ,fftw_0in_,fftw_out_,S_p_,S_q_)         
         if (verbose.gt.2) then
            call cl1_i4(n_r*n_svd_l*n_w_max,n_X_)
         end if !if (verbose.gt.2) then
         ic = 0
         do nr=0,n_r-1
            dw = 2.0d0*pi/(1.0d0*max(1,n_w_(nr)))
            dA = weight_p_(nr)*dw
c$$$  We assume that the fourier basis is orthonormal (not merely orthogonal)
            dAn = dA
            dsqrt_dAn = dsqrt(dAn)
            if (verbose.gt.1) then
               write(6,'(A,I0,A,F6.3,1X,F6.3)') ' % nr ',nr
     $              ,'; dw dA: ',dw,dA
            end if !if (verbose.gt.1) then
            n_w_t = floor(1.0d0*n_w_(nr)/2.0d0)
            ic_store = ic
            do nl=0,n_svd_l-1
               D_V_r = svd_polyval_V_r_(nl+nr*n_svd_l)
               D_s = svd_s_(nl)
               I_l = svd_l_(nl)
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,1X,I3,1X,F16.3,1X,F16.3)')
     $                 ' % % l I_l D_V_r D_s: ',nl,I_l,D_V_r,D_s
               end if !if (verbose.gt.2 .and. nr.lt.5) then
               nx = 0 + n_r*(nl + n_svd_l*ns);
               ic = ic_store
               do nw=0,n_w_(nr)-1
                  nwc = nw
                  if (nwc.ge.n_w_t) then
                     nwc = nwc - n_w_(nr)
                  end if !if (nwc.ge.n_w_t) then
                  if (abs(nwc).lt.n_w_t) then
                     flag_ic0_overflow = .false.
                  else
                     flag_ic0_overflow = .true.
                  end if !if (abs(nwc).lt.n_w_t) then
                  nwc = nw
                  if (nwc.ge.n_w_t) then
                     nwc = nwc - n_w_(nr)
                  end if !if (nwc.ge.n_w_t) then
                  flag_ict_overflow = .false.
                  nwd = nwc + I_l
                  if (abs(nwd).lt.n_w_t) then
                     call periodize_i(nwd,0,n_w_(nr),nwt)
                  else
                     nwt = 0
                     flag_ict_overflow = .true.
                  end if !if (abs(nwd).lt.n_w_t) then
                  flag_icr_overflow = .false.
                  nwd = nwc - I_l
                  if (abs(nwd).lt.n_w_t) then
                     call periodize_i(nwd,0,n_w_(nr),nwr)
                  else
                     nwr = 0
                     flag_icr_overflow = .true.
                  end if !if (abs(nwd).lt.n_w_t) then
                  ict = ic-nw+nwt
                  icr = ic-nw+nwr
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(10(A,I0))') ' % % % nl:' , nl , '; I_l:'
     $                    , I_l , '; ic_store:' , ic_store,
     $                    '; n_w_(nr):' , n_w_(nr) , '; nw:', nw ,
     $                    '; nwt:' ,nwt , '; nwr:' , nwr ,'; ic:' , ic ,
     $                    '; ict:' , ict , '; icr:', icr 
                     write(6,'(2(A,L1))') ' % % % flag_ict: ' ,
     $                    flag_ict_overflow , '; flag_icr: ' ,
     $                    flag_icr_overflow
                  end if !if (verbose.gt.3 .and. nr.lt.5) then
                  if ((flag_ic0_overflow.eqv..false.) .and.
     $                 (flag_ict_overflow.eqv..false.)) then
                  if (nw.gt.n_w_(nr)/2) then
                     nw_fix = nw - n_w_(nr) + n_w_max
                     if (verbose.gt.3 .and. nr.lt.5) then
                        write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                       ,'; nw_fix ',nw_fix,'; (full loop)'
                     end if
                     C_q = dconjg(D_s*D_V_r*S_q_(ict))
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
     $                       ,'; nw_fix ',nw_fix,'; (first orig)'
                     end if
                     C_q = dsqrt(0.0d0*0.5d0)*dconjg(D_s*D_V_r
     $                    *S_q_(ict))
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
     $                       ,'; nw_fix ',nw_fix,'; (then loop)'
                     end if
                     C_q = dsqrt(0.0d0*0.5d0)*dconjg(D_s*D_V_r
     $                    *S_q_(ict))
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
                     C_q = dconjg(D_s*D_V_r*S_q_(ict))
                     nw_C = nw_fix
                     if (verbose.gt.2) then
                        n_X_(nr + n_r*(nl + n_svd_l*nw_fix)) = n_X_(nr +
     $                       n_r*(nl + n_svd_l*nw_fix)) + 1
                     end if     !if (verbose.gt.2) then
                     Z_S_q__(nx + nr + ld_X*nw_C) = Z_S_q__(nx + nr +
     $                    ld_X*nw_C) +C_q*dsqrt_dAn
                  end if !if nw
                  end if !if ((flag_ic0_overflow.eqv..false.) .and. (flag_ict_overflow.eqv..false.)) then
                  ic = ic + 1
               enddo !do nw=0,n_w_(nr)-1
            enddo !do nl=0,n_svd_l-1
         enddo !do nr=0,n_r-1
         if (verbose.gt.2) then
            write(6,'(A,I0)') ' ns: ' , ns
            call print_all_i4__(n_r,n_svd_l*n_w_max,n_X_,'n_X_: ')
         end if !if (verbose.gt.2) then
      enddo !do ns=0,n_S_sub-1
c$$$      deallocate(V_r_);
      if (verbose.gt.2) then
         deallocate(n_X_)
      end if !if (verbose.gt.2) then
      if (verbose.gt.0) then
         write(6,'(A)') '[finished ti8_build_Z_S_q_quad_5]'
      end if !if (verbose.gt.0) then
      end
