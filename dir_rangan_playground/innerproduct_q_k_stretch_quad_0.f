!> Doxygen comment: ;\n
!> Calculates in innerproduct in k_q_ (i.e., fourier-space bessel) coordinates. ;\n
!> Assumes that M_q_ is the same size and dimensions as T_q_. ;\n
!> Assumes quasi-uniform polar-grid. ;\n
!> Assumes that C_q_ is large enough to hold all n_w_(n_r-1) modes ;\n
!> (assuming of course that n_w_(n_r-1) is the largest value within n_w_). ;\n
!> Stores C_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;\n
!> Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ;\n
!> Upgraded to use radial quadrature weight. ;\n
      subroutine innerproduct_q_k_stretch_quad_0(n_r,grid_p_,weight_p_
     $     ,n_w_,n_A,T_q_,M_q_,C_q_)
c$$$      Assumes that M_q_ is the same size and dimensions as T_q_. ;
c$$$      Assumes quasi-uniform polar-grid. ;
c$$$      Assumes that C_q_ is large enough to hold all n_w_(n_r-1) modes ;
c$$$      (assuming of course that n_w_(n_r-1) is the largest value within n_w_). ;
c$$$      Stores C_q_ in fft ordering (mode 0 , mode 1 , ... , mode -1) ;
c$$$      Upgraded to ignore frequencies of magnitude n_w_(nr)/2 or larger. ;
c$$$      Upgraded to use radial quadrature weight. ;\n
      implicit none
      integer verbose
      data verbose / 0 /
      integer n_r,n_w_(0:n_r-1),n_A
      real *8 grid_p_(0:n_r-1)
      real *8 weight_p_(0:n_r-1)
      complex *16 T_q_(0:n_A-1),M_q_(0:n_A-1)
      complex *16 C_q_(0:n_w_(n_r-1)-1)
      complex *16 C_q
      real *8 pi
      integer n_w_max,ic,nr,nw,nw_fix,nw_C
      real *8 dA,dAn,dw
      if (verbose.gt.0) then
         write(6,'(A,I0)')
     $        ' % [entering innerproduct_q_k_stretch_quad_0] n_r ',n_r
      end if
      pi = 4.0d0*datan(1.0d0)
      n_w_max = n_w_(n_r-1)
      if (verbose.gt.0) then
         write(6,'(A,I0)') ' % n_w_max ',n_w_max
      end if
      do nw=0,n_w_max-1
         C_q_(nw) = dcmplx(0.0d0,0.0d0)
      enddo
      C_q = dcmplx(0.0d0,0.0d0)
      ic = 0
      do nr=0,n_r-1
         dw = 2.0d0*pi/(1.0d0*max(1,n_w_(nr)))
         dA = weight_p_(nr)
c$$$         We assume that the fourier basis is orthonormal (not merely orthogonal)
         dAn = dA*dw
         do nw=0,n_w_(nr)-1
            if (nw.gt.n_w_(nr)/2) then
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full loop)'
               end if
               C_q = dconjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            else if (nw.eq.n_w_(nr)/2) then
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (first orig)'
               end if
               C_q = 0.0d0*0.5d0*dconjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
               nw_fix = nw - n_w_(nr) + n_w_max
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (then loop)'
               end if
               C_q = 0.0d0*0.5d0*dconjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            else
               nw_fix = nw
               if (verbose.gt.2 .and. nr.lt.5) then
                  write(6,'(A,I3,A,I3,A)') ' % % nw ',nw ,'; nw_fix '
     $                 ,nw_fix,'; (full orig)'
               end if
               C_q = dconjg(T_q_(ic))*M_q_(ic)
               nw_C = nw_fix
               C_q_(nw_C) = C_q_(nw_C) + C_q*dAn
            end if
            ic = ic + 1
         enddo
      enddo
      if (verbose.gt.0) then
         write(6,'(A)') ' % [finished innerproduct_q_k_stretch_quad_0]'
      end if
      end
