      subroutine innerproduct_q_k_only(n_r,grid_p_,n_w_,n_A,T_q_,M_q_
     $     ,C_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as T_q_.
c$$$      Assumes quasi-uniform polar-grid.
c$$$      Assumes that C_q_ is large enough to hold all n_w_(n_r-1) modes
c$$$      (assuming of course that n_w_(n_r-1) is the largest value within n_w_).
c$$$      Stores C_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1)
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 T_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16 C_q_(0:n_w_(n_r-1)-1)
      complex *16 C_q
      real *8 pi
      integer n_w_max,ic,nr,nw,nw_fix,nw_C
      real *8 R_pos,R_pre,dr,dw,dA,dAn
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % [entering innerproduct_q_k_only] n_r '
     $        ,n_r
      end if
      pi = 4.0*atan(1.0)
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      do nw=0,n_w_max-1
         C_q_(nw) = cmplx( 0.0 , 0.0 )
      enddo
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
         do nw=0,n_w_(nr)-1
            if (nw.gt.n_w_(nr)/2) then
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full loop)'
               end if
               C_q = conjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            else if (nw.eq.n_w_(nr)/2) then
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (half orig)'
               end if
               C_q = 0.5*conjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (half loop)'
               end if
               C_q = 0.5*conjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            else
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full orig)'
               end if
               C_q = conjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            end if
            ic = ic + 1
         enddo
      enddo
      if (verbose.gt.0) then
         write(6,'(A)') ' % [finished innerproduct_q_k_only]'
      end if
      end
