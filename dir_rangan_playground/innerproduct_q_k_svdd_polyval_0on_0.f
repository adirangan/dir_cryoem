!> Doxygen comment: ;\n
!> Sets up an array C_q_ to calculate the 'left side' ;\n
!> (i.e., translation-specific component) to the svd-expansion ;\n
!> of the translation operator. ;\n
!> Creates a matrix for translations.  \n
!> Uses svd-expansion defined via n_svd_d, .. , svd_U_d_.\n
!> Term-l of Translation-d is stored in \n
!> C_q_(l + ndv*n_svd_l).\n
!> The logical flag_S_vs_M determines the sign of the complex exponential.\n
!> flag_S_vs_M .eqv. .true. --> transformation applied to S, use +.\n
!> flag_S_vs_M .eqv. .false. --> transformation applied to M, use -.\n
!> Performs polynomial evaluation internally. ;\n
      subroutine innerproduct_q_k_svdd_polyval_0on_0(flag_S_vs_M
     $     ,svd_d_max,n_svd_d,svd_d_,n_svd_l,svd_l_,svd_U_d_,n_delta_v
     $     ,delta_x_,delta_y_,C_q_)
c$$$      Creates a matrix for translations.  
c$$$      Uses svd-expansion defined via n_svd_d, .. , svd_U_d_.
c$$$      Term-l of Translation-d is stored in 
c$$$      C_q_(l + ndv*n_svd_l).
c$$$      The logical flag_S_vs_M determines the sign of the complex exponential.
c$$$      flag_S_vs_M .eqv. .true. --> transformation applied to S, use +.
c$$$      flag_S_vs_M .eqv. .false. --> transformation applied to M, use -.
      implicit none
      integer verbose
      data verbose / 0 /
      logical warning_flag
      data warning_flag / .true. /
      logical flag_S_vs_M
      real *8 svd_d_max
      integer n_svd_d,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_d_(0:n_svd_d-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      integer n_delta_v
      real *8 delta_x_(0:n_delta_v-1),delta_x
      real *8 delta_y_(0:n_delta_v-1),delta_y
      complex *16 C_q_(0:n_svd_l*n_delta_v-1)
      real *8 pi
      real *8 svd_d_m,svd_d_c,svd_d(0:0)
      real *8 delta,omega,theta
      integer nl,I_l,ndv,nC,nl_pre
      integer *4 max_i4_f
      integer l_max
c$$$      l_max must equal or exceed order of bessel-function expansion
      complex *16 C_q
      complex *16, allocatable :: C_q_pre_(:)
      real *8, allocatable :: U_d_(:)
      if (verbose.gt.0) then
         write (6,'(A,I0,1X,I0,1X,I0)')
     $        ' % [entering innerproduct_q_k_svdd_polyval_0on_0] '
     $        ,n_svd_d ,n_svd_l ,n_delta_v
      end if
      l_max = max_i4_f(n_svd_l,svd_l_)
      allocate(C_q_pre_(-l_max:+l_max))
      allocate(U_d_(0:n_svd_l-1))
      pi = 4.0d0*datan(1.0d0)
      nC = 0
      svd_d_m = svd_d_max / 2.0d0
      svd_d_c = svd_d_m
      do ndv=0,n_delta_v-1
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         if (verbose.gt.1) then
            write (6,'(A,I0,1X,F6.3,1X,F6.3)')
     $           ' % ndv dx dy : ',ndv,delta_x,delta_y
         end if                 !if (verbose.gt.1) then
         delta = dsqrt(delta_x**2 + delta_y**2)
         if (delta.gt.svd_d_max .and. warning_flag) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)') 'Warning, delta '
     $           ,delta,'>',svd_d_max,'; ratio = ',delta/svd_d_max,
     $           ' in innerproduct_q_k_svdd_polyval_0on_0'
         end if
         omega = datan2(delta_y,delta_x)
         if (verbose.gt.1) then
            write (6,'(A,F6.3,1X,F6.3)') ' % delta omega : ',delta
     $           ,omega
         end if
c$$$            here we initialize C_q_pre using omega and l_max
c$$$            The goal is for 
c$$$            C_q_pre_(nl_pre) 
c$$$            to equal:
c$$$            S-type --> dcmplx( +dcos(theta) , +dsin(theta) )
c$$$            M-type --> dcmplx( +dcos(theta) , -dsin(theta) )
c$$$            with
c$$$            theta = nl_pre*(pi/2 - omega).
         C_q_pre_(0) = dcmplx(1.0d0,0.0d0)
         if (flag_S_vs_M.eqv..true.) theta = (pi/2.0d0 - omega)
         if (flag_S_vs_M.eqv..false.) theta = (pi/2.0d0 - omega + pi)
         if (flag_S_vs_M.eqv..true.) C_q = dcmplx( +dcos(theta) ,
     $        +dsin(theta) )
         if (flag_S_vs_M.eqv..false.) C_q = dcmplx( +dcos(theta) ,
     $        -dsin(theta) )
         do nl_pre=1,l_max
            C_q_pre_(+nl_pre) = C_q_pre_(+nl_pre-1)*C_q
         end do
         if (flag_S_vs_M.eqv..true.) C_q = dcmplx( +dcos(theta) ,
     $        -dsin(theta) )
         if (flag_S_vs_M.eqv..false.) C_q = dcmplx( +dcos(theta) ,
     $        +dsin(theta) )
         do nl_pre=1,l_max
            C_q_pre_(-nl_pre) = C_q_pre_(-nl_pre+1)*C_q
         end do
         svd_d(0) = (delta - svd_d_m)/svd_d_c
         do nl=0,n_svd_l-1
            call polyval_r8_reverse_0(n_svd_d,svd_U_d_(0+nl*n_svd_d),1
     $           ,svd_d(0),U_d_(nl))
            if (verbose.gt.1) then
               write (6,'(A,F6.3,A,F6.3,A,I0,A,F8.5)') ' % delta ' ,
     $              delta , ' svd_d ' ,svd_d(0) , ' U_d_(' , nl , ') ',
     $              U_d_(nl)
            end if !if (verbose.gt.1) then
         enddo !do nl=0,n_svd_l-1
         do nl=0,n_svd_l-1
            I_l = svd_l_(nl)
            if (abs(I_l).gt.l_max) then
               write(6,'(A,I0,A)') 'Warning, I_l: ',I_l
     $              ,'>l_max in innerproduct_q_k_svdd_polyval_0on_0'
            end if              !if (abs(I_l).gt.l_max) then
            if (verbose.gt.2) then
               write(6,'(A,I3,1X,I3)') ' % % nl I_l: ',nl,I_l
            end if              !if (verbose.gt.2) then
c$$$               theta = I_l*(pi/2 - omega)
c$$$               If transformation were applied to S we would use:
c$$$               C_q = dcmplx( +dcos(theta) , +dsin(theta) )
c$$$               If transformation were applied to M we would use:
c$$$               C_q = dcmplx( +dcos(theta) , -dsin(theta) )
            C_q = C_q_pre_(I_l)
c$$$               nC = nl + ndv*n_svd_l
            C_q_(nC) = U_d_(nl)*C_q
            nC = nC+1
         enddo                  !do nl=0,n_svd_l-1
      enddo                     !do ndv=0,n_delta_v-1
      deallocate(U_d_)
      deallocate(C_q_pre_)
      if (verbose.gt.0) then
         write (6,'(A)')
     $        ' % [finished innerproduct_q_k_svdd_polyval_0on_0]'
      end if
      end
