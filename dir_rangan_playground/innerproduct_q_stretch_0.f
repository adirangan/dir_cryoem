!> Doxygen comment: ;\n
!> Reorganizes an image in k_q_ (i.e., fourier-space bessel) coordinates ;\n
!> to allow for a subsequent multiplication to perform an innerproduct. ;\n
!> Assumes T_q_ corresponds to bessel-coefficients associated ;\n
!> with quasi-uniform polar-grid. ;\n
!> Thus, the (nr,nw) value of T_q_ is stored at: ;\n
!> T_q_(ic) = T_q_(nw + n_w_csum_(nr)). ;\n
!> Depending on flag_dconjg, we write T_q_(ic) or dconjg(T_q_(ic)) ;\n
!> to the (nr,nw) element of X_q_, stored at: ;\n
!> X_q_(nr + nw*ld_X). ;\n
!> We assume that X_q_ is large enough to hold all the values of T_q_. ;\n
!> This means X_q_ has to be at least size n_r * n_w_max, ;\n
!> with ld_X .ge. n_r. ;\n
!> We store each column of X_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;\n
!> Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ;\n
      subroutine innerproduct_q_stretch_0(n_r,grid_p_,n_w_,n_A,T_q_
     $     ,flag_dconjg,ld_X,X_q_)
c$$$      Assumes T_q_ corresponds to bessel-coefficients associated ;
c$$$      with quasi-uniform polar-grid. ;
c$$$      Thus, the (nr,nw) value of T_q_ is stored at: ;
c$$$      T_q_(ic) = T_q_(nw + n_w_csum_(nr)). ;
c$$$      Depending on flag_dconjg, we write T_q_(ic) or dconjg(T_q_(ic)) ;
c$$$      to the (nr,nw) element of X_q_, stored at: ;
c$$$      X_q_(nr + nw*ld_X). ;
c$$$      We assume that X_q_ is large enough to hold all the values of T_q_. ;
c$$$      This means X_q_ has to be at least size n_r * n_w_max, ;
c$$$      with ld_X .ge. n_r. ;
c$$$      We store each column of X_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;
c$$$      Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ;
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      complex *16 T_q_(0:n_A-1)
      logical flag_dconjg
      integer ld_X
      complex *16 X_q_(0:0)
      complex *16 C_q
      real *8 pi
      integer n_w_max,ic,nr,nw,nw_fix,nw_C
      real *8 R_pos,R_pre,dr,dw,dA,dAn,dsqrt_dAn
      integer nq,nk
      integer, allocatable :: n_X_(:)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % [entering innerproduct_q_stretch_0] n_r '
     $        ,n_r
      end if
      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      if (verbose.gt.2) then
         allocate(n_X_(0:n_r*n_w_max-1))
         call cl1_i4(n_r*n_w_max,n_X_)
      end if !if (verbose.gt.2) then
      C_q = dcmplx(0.0d0,0.0d0)
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
         dw = 2.0d0*pi/(1.0d0*max(1,n_w_(nr)))
         dA = (R_pre*dr + (dr**2)/2)*dw
c$$$         We assume that the fourier basis is orthonormal (not merely orthogonal)
         dAn = dA
         dsqrt_dAn = dsqrt(dAn)
         do nw=0,n_w_(nr)-1
            if (nw.gt.n_w_(nr)/2) then
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full loop)'
               end if
               C_q = T_q_(ic)
               if (flag_dconjg) C_q = dconjg(C_q)
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
     $                 ,nw_fix,'; (first orig)'
               end if
               C_q = dsqrt(0.0d0*0.5d0)*T_q_(ic)
               if (flag_dconjg) C_q = dconjg(C_q)
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (then loop)'
               end if
               C_q = dsqrt(0.0d0*0.5d0)*T_q_(ic)
               if (flag_dconjg) C_q = dconjg(C_q)
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
               if (flag_dconjg) C_q = dconjg(C_q)
               nw_C = nw_fix
               if (verbose.gt.2) then
                  n_X_(nr + nw_C*n_r) = n_X_(nr + nw_C*n_r) + 1
               end if !if (verbose.gt.2) then
               X_q_(nr + ld_X*nw_C) = X_q_(nr + ld_X*nw_C) + C_q
     $              *dsqrt_dAn
            end if
            ic = ic + 1
         enddo
      enddo
      if (verbose.gt.2) then
         call print_all_i4__(n_r,n_w_max,n_X_,'n_X_: ')
         deallocate(n_X_)
      end if !if (verbose.gt.2) then
      if (verbose.gt.0) then
         write(6,'(A)') ' % [finished innerproduct_q_stretch_0]'
      end if
      end
