      subroutine innerproduct_q_k_svdr(n_svd_r,svd_r_,n_svd_l,svd_l_
     $     ,svd_s_,svd_V_r_,n_r,grid_p_,n_w_,n_A,T_q_,M_q_,C_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as T_q_.
c$$$      Assumes quasi-uniform polar-grid defined via n_r , .. , n_A.
c$$$      Uses svd-expansion defined via n_svd_r , .. , svd_V_r_.
c$$$      Assumes that C_q_ is large enough to hold all n_w_max = n_w_(n_r-1) 
c$$$      modes for each of the n_svd_l terms in the svd-expansion
c$$$      (assuming of course that n_w_(n_r-1) is the largest value within n_w_).
c$$$      modes in C_q_ are stored in the order: (mode 0 , mode 1 , ... , mode -1).
c$$$      The mode-k for term-l is stored in C_q_(l + k*n_svd_l).
      implicit none
      integer verbose
      data verbose / 0 /
      logical warning_flag
      data warning_flag / .true. /
      integer n_svd_r,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_r_(0:n_svd_r-1),svd_s_(0:n_svd_l-1)
      real *8 svd_V_r_(0:n_svd_r*n_svd_l-1)
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 T_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16 C_q_(0:n_svd_l*n_w_(n_r-1)-1),C_q
      real *8 pi
      integer n_w_max,ic,ic_store,ict,nr,nw,nwt,nw_fix,nw_C
      real *8 R_q,R_pos,R_pre,dr,dw,dA,dAn
      real *8 svd_r_max,n_svd_r_v
      integer n_svd_r_pre,n_svd_r_pos
      real *8 d_svd_r_pre,d_svd_r_pos
      real *8 alpha_r,beta_r
      real *8 D_V_r,D_s
      integer l,I_l
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % [entering innerproduct_q_k_svdr] n_r '
     $        ,n_r
      end if
      pi = 4.0*atan(1.0)
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
      do l=0,n_svd_l-1
         do nw=0,n_w_max-1
            C_q_(l+nw*n_svd_l) = cmplx( 0.0 , 0.0 )
         enddo
      enddo
      C_q = cmplx( 0.0 , 0.0 )
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
         end if
         if (verbose.gt.1) then
            write(6,'(A,I0,A,F6.3,A,I0,1X,I0)') ' % nr ',nr,'; R_q ',R_q
     $           ,'; n_svd_r_pre n_svd_r_pos: ',n_svd_r_pre,n_svd_r_pos
            write(6,'(A,F6.3,1X,F6.3)') ' % % d_svd_r_pre d_svd_r_pos: '
     $           ,d_svd_r_pre,d_svd_r_pos
            write(6,'(A,F6.3,1X,F6.3)') ' % % alpha_r beta_r: ',alpha_r
     $           ,beta_r
         end if
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
         if (verbose.gt.1) then
            write(6,'(A,I0,A,F6.3,1X,F6.3,1X,F6.3)') ' % nr ',nr
     $           ,'; dr dw dA: ',dr,dw,dA
         end if
         ic_store = ic
         do l=0,n_svd_l-1
            D_V_r = beta_r*svd_V_r_(n_svd_r_pre + l*n_svd_r) +
     $           alpha_r*svd_V_r_(n_svd_r_pos + l*n_svd_r)
            D_s = svd_s_(l)
            I_l = svd_l_(l)
            if (verbose.gt.2 .and. nr.lt.5) then
               write(6,'(A,I3,1X,I3,1X,F16.3,1X,F16.3)')
     $              ' % % l I_l D_V_r D_s: ',l,I_l,D_V_r,D_s
            end if
            ic = ic_store
            do nw=0,n_w_(nr)-1
               call periodize_i(nw+I_l,0,n_w_(nr),nwt)
               ict = ic-nw+nwt
               if (verbose.gt.3 .and. nr.lt.5) then
                  write(6,'(A,I0,1X,I0,1X,I0)') ' % % % nw I_l nwt: ',nw
     $                 ,I_l,nwt
               end if
               if (nw.gt.n_w_(nr)/2) then
                  nw_fix = nw - n_w_(nr) + n_w_max
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                    ,'; nw_fix ',nw_fix,'; (full loop)'
                  end if
                  C_q = conjg(D_s*D_V_r*T_q_(ict))*M_q_(ic)
                  nw_C = l + nw_fix*n_svd_l
                  C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
               else if (nw.eq.n_w_(nr)/2) then
                  nw_fix = nw
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                    ,'; nw_fix ',nw_fix,'; (half orig)'
                  end if
                  C_q = 0.5*conjg(D_s*D_V_r*T_q_(ict))*M_q_(ic)
                  nw_C = l + nw_fix*n_svd_l
                  C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
                  nw_fix = nw - n_w_(nr) + n_w_max
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                    ,'; nw_fix ',nw_fix,'; (half loop)'
                  end if
                  C_q = 0.5*conjg(D_s*D_V_r*T_q_(ict))*M_q_(ic)
                  nw_C = l + nw_fix*n_svd_l
                  C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
               else
                  nw_fix = nw
                  if (verbose.gt.3 .and. nr.lt.5) then
                     write(6,'(A,I3,A,I3,A)') ' % % % nw ',nw
     $                    ,'; nw_fix ',nw_fix,'; (full orig)'
                  end if
                  C_q = conjg(D_s*D_V_r*T_q_(ict))*M_q_(ic)
                  nw_C = l + nw_fix*n_svd_l
                  C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
               end if
               ic = ic + 1
            enddo
         enddo         
      enddo
      if (verbose.gt.0) then
         write(6,'(A)') ' % [finished innerproduct_q_k_svdr]'
      end if
      end
