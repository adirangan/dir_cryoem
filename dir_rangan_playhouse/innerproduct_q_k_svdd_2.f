      subroutine innerproduct_q_k_svdd_2(flag_S_vs_M,n_svd_d,svd_d_
     $     ,n_svd_l,svd_l_,svd_U_d_,n_delta_v,delta_x_,delta_y_,C_q_)
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
      integer n_svd_d,n_svd_l,svd_l_(0:n_svd_l-1)
      real *8 svd_d_(0:n_svd_d-1)
      real *8 svd_U_d_(0:n_svd_d*n_svd_l-1)
      integer n_delta_v
      real *8 delta_x_(0:n_delta_v-1),delta_x
      real *8 delta_y_(0:n_delta_v-1),delta_y
      complex *16 C_q_(0:n_svd_l*n_delta_v-1)
      real *8 pi
      real *8 svd_d_max,n_svd_d_v
      integer n_svd_d_pre,n_svd_d_pos
      real *8 d_svd_d_pre,d_svd_d_pos
      real *8 alpha_d,beta_d
      real *8 D_U_d
      real *16 delta,omega,theta
      integer l,I_l,ndv,nC,nl_pre
      integer l_max
c$$$      l_max must equal or exceed order of bessel-function expansion
      parameter (l_max = 15)
      complex *16 C_q,C_q_pre_(-l_max:l_max)
      if (verbose.gt.0) then
         write (6,'(A,I0,1X,I0,1X,I0)')
     $        ' % [entering innerproduct_q_k_svdd_2] ',n_svd_d,n_svd_l
     $        ,n_delta_v
      end if
      pi = 4.0*atan(1.0)
      nC = 0
      do ndv=0,n_delta_v-1
         delta_x = delta_x_(ndv)
         delta_y = delta_y_(ndv)
         if (verbose.gt.1) then
            write (6,'(A,I0,1X,F6.3,1X,F6.3)')
     $           ' % ndv dx dy : ',ndv,delta_x,delta_y
         end if                 !if (verbose.gt.1) then
         delta = dsqrt(delta_x**2 + delta_y**2)
         svd_d_max = svd_d_(n_svd_d-1)
         if (delta.gt.svd_d_max*1.0625 .and. warning_flag) then
            write(6,'(A,F6.3,A,F6.3,A,F6.3,A)') 'Warning, delta '
     $           ,delta,'>',svd_d_max,'; ratio = ',delta/svd_d_max,
     $           ' in transf_svd_q_k_svdd_2'
         end if
         omega = atan2(delta_y,delta_x)
         if (verbose.gt.1) then
            write (6,'(A,F6.3,1X,F6.3)') ' % delta omega : ',delta
     $           ,omega
         end if
c$$$            here we initialize C_q_pre using omega and l_max
c$$$            The goal is for 
c$$$            C_q_pre_(nl_pre) 
c$$$            to equal:
c$$$            S-type --> cmplx( +cos(theta) , +sin(theta) )
c$$$            M-type --> cmplx( +cos(theta) , -sin(theta) )
c$$$            with
c$$$            theta = nl_pre*(pi/2 - omega).
         C_q_pre_(0) = cmplx( 1.0 , 0.0 )
         if (flag_S_vs_M.eqv..true.) theta = (pi/2 - omega)
         if (flag_S_vs_M.eqv..false.) theta = (pi/2 - omega + pi)
         if (flag_S_vs_M.eqv..true.) C_q = cmplx( +cos(theta) ,
     $        +sin(theta) )
         if (flag_S_vs_M.eqv..false.) C_q = cmplx( +cos(theta) ,
     $        -sin(theta) )
         do nl_pre=1,l_max
            C_q_pre_(+nl_pre) = C_q_pre_(+nl_pre-1)*C_q
         end do
         if (flag_S_vs_M.eqv..true.) C_q = cmplx( +cos(theta) ,
     $        -sin(theta) )
         if (flag_S_vs_M.eqv..false.) C_q = cmplx( +cos(theta) ,
     $        +sin(theta) )
         do nl_pre=1,l_max
            C_q_pre_(-nl_pre) = C_q_pre_(-nl_pre+1)*C_q
         end do
         n_svd_d_v = (n_svd_d-1)*delta/svd_d_max
         n_svd_d_pre = max(0,min(n_svd_d-1,floor(n_svd_d_v)))
         n_svd_d_pos = max(0,min(n_svd_d-1,ceiling(n_svd_d_v)))
         d_svd_d_pre = abs(n_svd_d_v - n_svd_d_pre)
         d_svd_d_pos = abs(n_svd_d_pos - n_svd_d_v)
         if (d_svd_d_pre+d_svd_d_pos.le.0.0d0) then
            alpha_d = 0.0d0
            beta_d = 1.0d0
         else
            alpha_d = d_svd_d_pre / (d_svd_d_pre + d_svd_d_pos)
            beta_d = d_svd_d_pos / (d_svd_d_pre + d_svd_d_pos)
         end if
         if (verbose.gt.1) then
            write (6,'(A,I0,1X,I0)') ' % n_svd_d_pre n_svd_d_pos : '
     $           ,n_svd_d_pre,n_svd_d_pos
            write (6,'(A,F6.3,1X,F6.3)')
     $           ' % d_svd_d_pre d_svd_d_pos : ' ,d_svd_d_pre
     $           ,d_svd_d_pos
            write (6,'(A,F6.3,1X,F6.3)') ' % alpha_d beta_d : '
     $           ,alpha_d,beta_d
         end if
         do l=0,n_svd_l-1
            I_l = svd_l_(l)
            if (abs(I_l).gt.l_max) then
               write(6,'(A,I0,A)') 'Warning, I_l: ',I_l
     $              ,'>l_max in innerproduct_q_k_svdd_2'
            end if              !if (abs(I_l).gt.l_max) then
            if (verbose.gt.2) then
               write(6,'(A,I3,1X,I3)') ' % % l I_l: ',l,I_l
            end if              !if (verbose.gt.2) then
c$$$               theta = I_l*(pi/2 - omega)
c$$$               If transformation were applied to S we would use:
c$$$               C_q = cmplx( +cos(theta) , +sin(theta) )
c$$$               If transformation were applied to M we would use:
c$$$               C_q = cmplx( +cos(theta) , -sin(theta) )
            C_q = C_q_pre_(I_l)
            D_U_d = beta_d*svd_U_d_(n_svd_d_pre + l*n_svd_d) +
     $           alpha_d*svd_U_d_(n_svd_d_pos + l*n_svd_d)
c$$$               nC = l + ndv*n_svd_l
            C_q_(nC) = D_U_d*C_q
            nC = nC+1
         enddo                  !do l=0,n_svd_l-1
      enddo                     !do ndv=0,n_delta_v-1
      if (verbose.gt.0) then
         write (6,'(A)') ' % [finished innerproduct_q_k_svdd_2]'
      end if
      end
