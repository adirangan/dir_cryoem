      subroutine innerproduct_q_flip_stretch(n_r,grid_p_,n_w_,n_A,T_q_
     $     ,flag_conjg,ld_X,X_q_)
c$$$      Assumes T_q_ corresponds to bessel-coefficients associated 
c$$$      with quasi-uniform polar-grid.
c$$$      Thus, the (nr,nw) value of T_q_ is stored at:
c$$$      T_q_(ic) = T_q_(nw + n_w_sum_(nr)).
c$$$      Depending on flag_conjg, we write T_q_(ic) or conjg(T_q_(ic))
c$$$      to the (nr,nw) element of X_q_, stored at:
c$$$      X_q_(nr + nw*ld_X).
c$$$      We assume that X_q_ is large enough to hold all the values of T_q_.
c$$$      This means X_q_ has to be at least size n_r * n_w_max, 
c$$$      with ld_X .ge. n_r.
c$$$      We store each column of X_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 T_q_(0:n_A-1)
      logical flag_conjg
      integer ld_X
      complex *16 X_q_(0:0)
      complex *16 C_q
      real *8 pi
      integer n_w_max,ic_store,ic,ict,icr,nr,nw,nwt,nwr,nw_fix,nw_C
      integer n_w_t ! positive threshold for overflow. ;
      integer nwc ! centered nw. ;
      integer nwd ! displaced nw. ;
      logical flag_ict_overflow ! notes whether or not M_q_ coefficient should be set to 0. ;
      logical flag_icr_overflow ! notes whether or not M_q_ coefficient should be set to 0. ;
      real *8 R_pos,R_pre,dr,dw,dA,dAn,dsqrt_dAn
      integer nq,nk
      integer, allocatable :: n_X_(:)
      if (verbose.gt.0) then
         write(6,'(A,I0)')
     $        ' % [entering innerproduct_q_flip_stretch] n_r ',n_r
      end if
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      if (verbose.gt.2) then
         allocate(n_X_(0:n_r*n_w_max-1))
         call cl1_i4(n_r*n_w_max,n_X_)
      end if !if (verbose.gt.2) then
      C_q = cmplx( 0.0 , 0.0 )
      ic = 0
      do nr=0,n_r-1
         if (nr.gt.0) then
            R_pre = 0.5*(grid_p_(nr-1) + grid_p_(nr))
         else
            R_pre = grid_p_(0)
         end if
         if (nr.lt.n_r-1) then
            R_pos = 0.5*(grid_p_(nr+1) + grid_p_(nr))
         else
            R_pos = grid_p_(n_r-1)
         end if
         dr = R_pos - R_pre
c$$$         We set the zero-mode to zero
         if (grid_p_(nr).le.0.0d0) then
            dr = 0.0d0
         end if
         if (verbose.gt.1) then
            write(6,'(A,I0,A,I0,A,F6.3,1X,F6.3,1X,F6.3)') ' % nr ',nr
     $           ,'; n_w_(nr) ',n_w_(nr),'; R_pre,R_pos,dr: ',R_pre
     $           ,R_pos,dr
         end if
         dw = 2*pi/(1.0d0*max(1,n_w_(nr)))
         dA = (R_pre*dr + (dr**2)/2)*dw
c$$$         We assume that the fourier basis is orthonormal (not merely orthogonal)
         dAn = dA
         dsqrt_dAn = dsqrt(dAn)
         n_w_t = floor(1.0d0*n_w_(nr)/2.0d0)
         do nw=0,n_w_(nr)-1
            nwc = nw
            if (nwc.ge.n_w_t) then
               nwc = nwc - n_w_(nr)
            end if              !if (nwc.ge.n_w_t) then
            flag_ict_overflow = .false.
            nwd = nwc + I_l
            if (abs(nwd).lt.n_w_t) then
               call periodize_i(nwd,0,n_w_(nr),nwt)
            else
               nwt = 0
               flag_ict_overflow = .true.
            end if              !if (abs(nwd).lt.n_w_t) then
            flag_icr_overflow = .false.
            nwd = nwc - I_l
            if (abs(nwd).lt.n_w_t) then
               call periodize_i(nwd,0,n_w_(nr),nwr)
            else
               nwr = 0
               flag_icr_overflow = .true.
            end if              !if (abs(nwd).lt.n_w_t) then
            ict = ic-nw+nwt
            icr = ic-nw+nwr
            if (verbose.gt.3 .and. nr.lt.5) then
               write(6,'(10(A,I0))') ' % % % nl:' , nl , '; I_l:'
     $              , I_l , '; ic_store:' , ic_store,
     $              '; n_w_(nr):' , n_w_(nr) , '; nw:', nw ,
     $              '; nwt:' ,nwt , '; nwr:' , nwr ,'; ic:' , ic ,
     $              '; ict:' , ict , '; icr:', icr 
               write(6,'(2(A,L1))') ' % % % flag_ict: ' ,
     $              flag_ict_overflow , '; flag_icr: ' ,
     $              flag_icr_overflow
            end if              !if (verbose.gt.3 .and. nr.lt.5) then
            if (flag_ict_overflow.eqv..false.) then
            if (nw.gt.n_w_(nr)/2) then
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full loop)'
               end if
               C_q = T_q_(ic)
               if (flag_conjg) C_q = conjg(T_q_(ict))
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
            else if (nw.eq.n_w_(nr)/2) then
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (half orig)'
               end if
               C_q = dsqrt(0.5d0)*T_q_(ic)
               if (flag_conjg) C_q = conjg(dsqrt(0.5d0)*T_q_(ict))
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (half loop)'
               end if
               C_q = dsqrt(0.5d0)*T_q_(ic)
               if (flag_conjg) C_q = conjg(dsqrt(0.5d0)*T_q_(ict))
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
            else
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full orig)'
               end if
               C_q = T_q_(ic)
               if (flag_conjg) C_q = conjg(T_q_(ict))
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
            end if !if (nw.gt.n_w_(nr)/2) then
            end if !if (flag_ict_overflow.eqv..false.) then
            ic = ic + 1
         enddo
      enddo
      if (verbose.gt.2) then
         call write_all_i4__(n_r,n_w_max,n_X_,6,'n_X_: ')
         deallocate(n_X_)
      end if !if (verbose.gt.2) then
      if (verbose.gt.0) then
         write(6,'(A)') ' % [finished innerproduct_q_flip_stretch]'
      end if
      end
