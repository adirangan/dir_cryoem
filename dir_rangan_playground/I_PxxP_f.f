      real *8 function I_PxRTRTxP_f(
     $     K
     $     ,phi_M
     $     ,d_M
     $     ,omega_M
     $     ,phi_S
     $     ,d_S
     $     ,omega_S
     $     ,flag_RTRT_vs_RTTR
     $     ,gamma_z_est
     $     ,delta_x_est
     $     ,delta_y_est
     $     ,gamma_z_upd 
     $     ,delta_x_upd
     $     ,delta_y_upd
     $     )
c$$$      Calculates the integral of a product of: ;
c$$$      one rotated and translated version of a plane-wave multiplied by a planar-function
c$$$      with another rotated and translated verions of a plane-wave multiplied by a planar-function,
c$$$      where the domain of integration is the disc of radius K. ;
c$$$      The integral is either [flag_RTRT_vs_RTTR.eqv..true.]: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} conj(R_{+gamma_z_upd} T_{+delta_upd_} S) * (T_{-delta_est} R_{-gamma_est} conj(C) M) * dpsi * kdk ;
c$$$      or [flag_RTRT_vs_RTTR.eqv..false.]: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} conj(T_{+gamma_z_upd} R_{+delta_upd_} S) * (T_{-delta_est} R_{-gamma_est} conj(C) M) * dpsi * kdk ;
c$$$      with: ;
c$$$      C = 2*k*cos(psi-phi_M) ;
c$$$      S = 2*k*cos(psi-phi_S) * exp(+i*2*pi*k*d_S*cos(psi-omega_S)) ;
c$$$      M = exp(+i*2*pi*k*d_M*cos(psi-omega_M)) ;
c$$$      ;
c$$$      Note: ;
c$$$      % TS = exp(-i*2*pi*dot(k_,delta_upd_)) * 2*tmp_k*cos(tmp_w-fix_phi_S) * exp(+i*2*pi*dot(k_,delta_S_)) ;
c$$$      %    = 2*tmp_k*cos(tmp_w-fix_phi_S) * exp(+i*2*pi*dot(k_,delta_S_-delta_upd_));
c$$$      % RTS = 2*dot(R_{-gamma_upd}k_,n_fix_phi_S) * exp(+i*2*pi*dot(R_{-gamma_upd}k_,delta_S_-delta_upd_)) ;
c$$$      %     = 2*dot(k_,R_{+gamma_upd}n_fix_phi_S) * exp(+i*2*pi*dot(k_,R_{+gamma_upd}(delta_S_-delta_upd_))) ;
c$$$      %     = 2*tmp_k*cos(tmp_w-(fix_phi_S+gamma_upd)) * exp(+i*2*pi*dot(k_,R_{+gamma_upd}*(delta_S_-delta_upd_))) ;
c$$$      % CM = 2*tmp_k*cos(tmp_w-fix_phi_M) * exp(+i*2*pi*dot(k_,delta_M_)) ;
c$$$      % RCM = 2*dot(R_{+gamma_est}k_,n_phi_M) * exp(+i*2*pi*dot(R_{+gamma_est}k_,delta_M_)) ;
c$$$      %     = 2*dot(k_,R_{-gamma_est}n_phi_M) * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_)) ;
c$$$      %     = 2*tmp_k*cos(tmp_w - (fix_phi_M-gamma_est)) * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_)) ;
c$$$      % TRCM = exp(+i*2*pi*dot(k_,delta_est_)) * 2*tmp_k*cos(tmp_w - (fix_phi_M-gamma_est)) * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_)) ;
c$$$      %      = 2*tmp_k*cos(tmp_w - (fix_phi_M-gamma_est)) * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_ + delta_est_)) ;
c$$$      % <RTS,TRCM> = 2*tmp_k*cos(tmp_w-(fix_phi_S+gamma_upd)) * exp(+i*2*pi*dot(k_,-R_{+gamma_upd}*(delta_S_-delta_upd_))) ;
c$$$      %              * 2*tmp_k*cos(tmp_w - (fix_phi_M-gamma_est)) * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_ + delta_est_)) ;
c$$$      %            = 2*tmp_k*cos(tmp_w-(fix_phi_S+gamma_upd)) * 2*tmp_k*cos(tmp_w - (fix_phi_M-gamma_est)) 
c$$$      %              * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_ + delta_est_ - R_{+gamma_upd}*(delta_S_-delta_upd_))) ;
c$$$      %            = tmp_k^2 * (2*cos(2*tmp_w-(fix_phi_S+gamma_upd)-(fix_phi_M-gamma_est)) + 2*cos((fix_phi_M-gamma_est) - (fix_phi_S+gamma_upd))) 
c$$$      %              * exp(+i*2*pi*dot(k_,R_{-gamma_est}delta_M_ + delta_est_ - R_{+gamma_upd}*(delta_S_-delta_upd_))) ;
c$$$      %            = I_PxxP( K , fix_phi_M-gamma_est , R_{-gamma_est}delta_M_+delta_est_ , fix_phi_S+gamma_upd , R_{+gamma_upd}(delta_S_-delta_upd_) ) ; %<-- use temporary variables for delta and omega. ;
c$$$      % RS = 2*dot(R_{-gamma_upd}k_,n_fix_phi_S) * exp(+i*2*pi*dot(R_{-gamma_upd}k_,delta_S_)) ;
c$$$      %    = 2*dot(k_,R_{+gamma_upd}n_fix_phi_S) * exp(+i*2*pi*dot(k_,R_{+gamma_upd}delta_S_)) ;
c$$$      % TRS = exp(-i*2*pi*dot(k_,delta_upd_)) * 2*dot(k_,R_{+gamma_upd}n_fix_phi_S) * exp(+i*2*pi*dot(k_,R_{+gamma_upd}delta_S_)) ;
c$$$      %     = 2*tmp_k*cos(tmp_w-(fix_phi_S+gamma_upd)) * exp(+i*2*pi*dot(k_,R_{+gamma_upd}delta_S_-delta_upd_)) ;
c$$$      % <TRS,TRCM> = I_PxxP( K , fix_phi_M-gamma_est , R_{-gamma_est}delta_M_+delta_est_ , fix_phi_S+gamma_upd , R_{+gamma_upd}delta_S_-delta_upd_ ) ; %<-- use temporary variables for delta and omega. ;
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 K ! real *8: radius of disc ;
      real *8 phi_M ! real *8: orientation of planar-function ;
      real *8 d_M ! real *8: frequency of plane-wave M ;
      real *8 omega_M ! real *8: orientation of plane-wave M ;
      real *8 phi_S ! real *8: orientation of planar-function ;
      real *8 d_S ! real *8: frequency of plane-wave S ;
      real *8 omega_S ! real *8: orientation of plane-wave S ;
      logical flag_RTRT_vs_RTTR ! logical: <RTS,RTCM> (.true.) or <TRS,RTCM> (.false.) ;
      real *8 gamma_z_est,delta_x_est,delta_y_est ! real *8: estimated (original) image parameters. ;
      real *8 gamma_z_upd,delta_x_upd,delta_y_upd ! real *8: updates (increments) to image parameters. ;
      real *8 fix_phi_S,fix_phi_M
      real *8 delta_est_x_c_(0:1)
      real *8 delta_upd_x_c_(0:1)
      real *8 fix_delta_S_x_c_(0:1)
      real *8 fix_delta_M_x_c_(0:1)
      real *8 tmp_delta_S_x_c_(0:1)
      real *8 tmp_delta_M_x_c_(0:1)
      real *8 tmp_delta_S_x_p_r,tmp_delta_S_x_p_w
      real *8 tmp_delta_M_x_p_r,tmp_delta_M_x_p_w
      real *8 tmp_phi_S,tmp_phi_M
      real *8 I_SX,I_MX
      real *8 I_PxxP_f ! real *8 function evaluation ;
      if (verbose.gt.0) then
         write(6,'(A)') ' [entering I_PxRTRTxP_f]'
      end if !if (verbose.gt.0) then
      fix_phi_S = phi_S
      fix_phi_M = phi_M
      delta_est_x_c_(0) = delta_x_est
      delta_est_x_c_(1) = delta_y_est
      delta_upd_x_c_(0) = delta_x_upd
      delta_upd_x_c_(1) = delta_y_upd
      fix_delta_S_x_c_(0) = d_S*dcos(omega_S)
      fix_delta_S_x_c_(1) = d_S*dsin(omega_S)
      fix_delta_M_x_c_(0) = d_M*dcos(omega_M)
      fix_delta_M_x_c_(1) = d_M*dsin(omega_M)
      tmp_delta_S_x_c_(0) = fix_delta_S_x_c_(0) - delta_upd_x_c_(0)
      tmp_delta_S_x_c_(1) = fix_delta_S_x_c_(1) - delta_upd_x_c_(1)
      call rotate_delta(+gamma_z_upd,tmp_delta_S_x_c_
     $     ,tmp_delta_S_x_c_)
      tmp_delta_S_x_p_r = dsqrt(tmp_delta_S_x_c_(0)**2 +
     $     tmp_delta_S_x_c_(1)**2)
      tmp_delta_S_x_p_w = datan2(tmp_delta_S_x_c_(1)
     $     ,tmp_delta_S_x_c_(0))
      tmp_delta_M_x_c_(0) = fix_delta_M_x_c_(0)
      tmp_delta_M_x_c_(1) = fix_delta_M_x_c_(1)
      call rotate_delta(-gamma_z_est,tmp_delta_M_x_c_
     $     ,tmp_delta_M_x_c_)
      tmp_delta_M_x_c_(0) = tmp_delta_M_x_c_(0) + delta_est_x_c_(0)
      tmp_delta_M_x_c_(1) = tmp_delta_M_x_c_(1) + delta_est_x_c_(1)
      tmp_delta_M_x_p_r = dsqrt(tmp_delta_M_x_c_(0)**2 +
     $     tmp_delta_M_x_c_(1)**2)
      tmp_delta_M_x_p_w = datan2(tmp_delta_M_x_c_(1)
     $     ,tmp_delta_M_x_c_(0))
      tmp_phi_M = fix_phi_M - gamma_z_est
      tmp_phi_S = fix_phi_S + gamma_z_upd
      if (verbose.gt.0) then
         write(6,'(6(A,F16.8))')
     $     ' K: ' , K
     $     ,' tmp_phi_M: ' , tmp_phi_M
     $     ,' tmp_delta_M_x_p_r: ' , tmp_delta_M_x_p_r
     $     ,' tmp_delta_M_x_p_w: ' , tmp_delta_M_x_p_w
     $     ,' tmp_phi_S: ' , tmp_phi_S
     $     ,' tmp_delta_S_x_p_r: ' , tmp_delta_S_x_p_r
     $     ,' tmp_delta_S_x_p_w: ' , tmp_delta_S_x_p_w
      end if !if (verbose.gt.0) then
      I_SX = I_PxxP_f(
     $     K
     $     ,tmp_phi_M
     $     ,tmp_delta_M_x_p_r
     $     ,tmp_delta_M_x_p_w
     $     ,tmp_phi_S
     $     ,tmp_delta_S_x_p_r 
     $     ,tmp_delta_S_x_p_w
     $     )
      tmp_delta_S_x_c_(0) = fix_delta_S_x_c_(0)
      tmp_delta_S_x_c_(1) = fix_delta_S_x_c_(1)
      call rotate_delta(+gamma_z_upd,tmp_delta_S_x_c_
     $     ,tmp_delta_S_x_c_)
      tmp_delta_S_x_c_(0) = tmp_delta_S_x_c_(0) - delta_upd_x_c_(0)
      tmp_delta_S_x_c_(1) = tmp_delta_S_x_c_(1) - delta_upd_x_c_(1)
      tmp_delta_S_x_p_r = dsqrt(tmp_delta_S_x_c_(0)**2 +
     $     tmp_delta_S_x_c_(1)**2)
      tmp_delta_S_x_p_w = datan2(tmp_delta_S_x_c_(1)
     $     ,tmp_delta_S_x_c_(0))
      if (verbose.gt.0) then
         write(6,'(6(A,F16.8))')
     $     ' K: ' , K
     $     ,' tmp_phi_M: ' , tmp_phi_M
     $     ,' tmp_delta_M_x_p_r: ' , tmp_delta_M_x_p_r
     $     ,' tmp_delta_M_x_p_w: ' , tmp_delta_M_x_p_w
     $     ,' tmp_phi_S: ' , tmp_phi_S
     $     ,' tmp_delta_S_x_p_r: ' , tmp_delta_S_x_p_r
     $     ,' tmp_delta_S_x_p_w: ' , tmp_delta_S_x_p_w
      end if !if (verbose.gt.0) then
      I_MX = I_PxxP_f(
     $     K
     $     ,tmp_phi_M
     $     ,tmp_delta_M_x_p_r
     $     ,tmp_delta_M_x_p_w
     $     ,tmp_phi_S
     $     ,tmp_delta_S_x_p_r 
     $     ,tmp_delta_S_x_p_w
     $     )
      I_PxRTRTxP_f = dcmplx(0.0d0,0.0d0)
      if (flag_RTRT_vs_RTTR.eqv..true.) then
         I_PxRTRTxP_f = I_SX
      end if !if (flag_RTRT_vs_RTTR.eqv..true.) then
      if (flag_RTRT_vs_RTTR.eqv..false.) then
         I_PxRTRTxP_f = I_MX
      end if !if (flag_RTRT_vs_RTTR.eqv..false.) then
      if (verbose.gt.0) then
         write(6,'(A)') ' [finished I_PxRTRTxP_f]'
      end if !if (verbose.gt.0) then
      end !function

      real *8 function l2_xRPx_f(K,phi_M,phi_S,gamma_z)
c$$$      calculates the l2-norm of a product of a planar-function ;
c$$$      and a rotated plane-wave multiplied by a planar-function, ;
c$$$      where the domain of integration is the disc of radius K: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} conj(C*(R*S)) * (C*(R*S)) * dpsi * kdk ;
c$$$      with: ;
c$$$      C = 2*k*cos(psi-phi_M) ;
c$$$      S = 2*k*cos(psi-phi_S) * exp(+i*2*pi*k*d_S*cos(psi-omega_S)) ;
c$$$      R = rotation by gamma_z ;
c$$$      R*S = 2*k*cos(psi-gamma_z-phi_S) * exp(+i*2*pi*k*d_S*cos(psi-gamma_z-omega_S)) ;
      implicit none
      real *8 K ! real *8: radius of disc ;
      real *8 phi_M ! real *8: orientation of planar-function ;
      real *8 phi_S ! real *8: orientation of planar-function ;
      real *8 gamma_z ! real *8: rotation. ;
      real *8 I ! real *8: integral ;
      real *8 phi_n ! real *8: phase. ;
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      phi_n = phi_S + gamma_z - phi_M
      I = 2.0d0*pi * K**6 / 3.0d0 * (1.0d0 + 2.0d0*dcos(phi_n)
     $     *dcos(phi_n))
      l2_xRPx_f = dsqrt(I)
      end !function

      real *8 function I_PxxP_f(K,phi_M,d_M,omega_M,phi_S,d_S,omega_S)
c$$$      calculates the integral of a product of two plane-waves each multiplied by a planar-function, ;
c$$$      where the domain of integration is the disc of radius K: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} conj(C*S) * M * dpsi * kdk ;
c$$$      with: ;
c$$$      C = 2*k*cos(psi-phi_M) ;
c$$$      S = 2*k*cos(psi-phi_S) * exp(+i*2*pi*k*d_S*cos(psi-omega_S)) ;
c$$$      M = exp(+i*2*pi*k*d_M*cos(psi-omega_M)) ;
      implicit none
      real *8 K ! real *8: radius of disc ;
      real *8 phi_M ! real *8: orientation of planar-function ;
      real *8 d_M ! real *8: frequency of plane-wave M ;
      real *8 omega_M ! real *8: orientation of plane-wave M ;
      real *8 phi_S ! real *8: orientation of planar-function ;
      real *8 d_S ! real *8: frequency of plane-wave S ;
      real *8 omega_S ! real *8: orientation of plane-wave S ;
      real *8 I ! real *8: integral ;
      real *8 d_S_(0:1) ! real *8 array (length 2): delta_S_ ;
      real *8 d_M_(0:1) ! real *8 array (length 2): delta_M_ ;
      real *8 d_T_(0:1) ! real *8 array (length 2): delta_M_ - delta_S_ ;
      real *8 d_T ! real *8: frequency of plane-wave T ;
      real *8 omega_T ! real *8: orientation of plane-wave T ;
      real *8 phi_p,phi_n ! real *8: phase. ;
      real *8 t,tK ! real *8: exponents ;
      real *8 j0,j2,j4 ! real *8: bessel function evaluations
      real *8 mode_2,mode_0 ! real *8: components of integral
      real *8 pi
      pi = 4.0d0*datan(1.0d0)
      d_S_(0) = d_S*dcos(omega_S)
      d_S_(1) = d_S*dsin(omega_S)
      d_M_(0) = d_M*dcos(omega_M)
      d_M_(1) = d_M*dsin(omega_M)
      d_T_(0) = d_M_(0) - d_S_(0)
      d_T_(1) = d_M_(1) - d_S_(1)
      d_T = dsqrt(d_T_(0)**2 + d_T_(1)**2)
      omega_T = datan2(d_T_(1),d_T_(0))
      phi_p = phi_S + phi_M
      phi_n = phi_S - phi_M
      t = 2*pi*d_T
      tK = t*K
      call jbess(0.0d0,1,tK,j0)
      call jbess(2.0d0,1,tK,j2)
      call jbess(4.0d0,1,tK,j4)
      mode_2 = -4*pi*cos(phi_p - 2*omega_T) * K**4 * (j2+j4)/6.0d0
      mode_0 = 4*pi*cos(phi_n) * K**4 * (j0/4.0d0+j2/6.0d0-j4/12.0d0)
      I_PxxP_f = mode_0 + mode_2
      end !function

      real *8 function I_PxP_f(K,phi,d_M,omega_M,d_S,omega_S)
c$$$      calculates the integral of a product of two plane-waves multiplied by a planar-function, ;
c$$$      where the domain of integration is the disc of radius K: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} conj(C*S) * M * dpsi * kdk ;
c$$$      with: ;
c$$$      C = -2*i*k*cos(psi-phi) ;
c$$$      S = exp(+i*2*pi*k*d_S*cos(psi-omega_S)) ;
c$$$      M = exp(+i*2*pi*k*d_M*cos(psi-omega_M)) ;
      implicit none
      real *8 K ! real *8: radius of disc ;
      real *8 phi ! real *8: orientation of planar-function ;
      real *8 d_M ! real *8: frequency of plane-wave M ;
      real *8 omega_M ! real *8: orientation of plane-wave M ;
      real *8 d_S ! real *8: frequency of plane-wave S ;
      real *8 omega_S ! real *8: orientation of plane-wave S ;
      real *8 I ! real *8: integral ;
      real *8 d_S_(0:1) ! real *8 array (length 2): delta_S_ ;
      real *8 d_M_(0:1) ! real *8 array (length 2): delta_M_ ;
      real *8 d_T_(0:1) ! real *8 array (length 2): delta_M_ - delta_S_ ;
      real *8 d_T ! real *8: frequency of plane-wave T ;
      real *8 omega_T ! real *8: orientation of plane-wave T ;
      real *8 I_xP_f ! real *8 function evaluation ;
      d_S_(0) = d_S*dcos(omega_S)
      d_S_(1) = d_S*dsin(omega_S)
      d_M_(0) = d_M*dcos(omega_M)
      d_M_(1) = d_M*dsin(omega_M)
      d_T_(0) = d_M_(0) - d_S_(0)
      d_T_(1) = d_M_(1) - d_S_(1)
      d_T = dsqrt(d_T_(0)**2 + d_T_(1)**2)
      omega_T = datan2(d_T_(1),d_T_(0))
      I_PxP_f = -I_xP_f(K,phi,d_T,omega_T)
      end !function

      real *8 function I_xP_f(K,phi,d,omega)
c$$$      calculates the integral of a plane-wave multiplied by a planar-function, ;
c$$$      where the domain of integration is the disc of radius K: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} -2*i*k*cos(psi-phi) * exp(+i*2*pi*k*d*cos(psi-omega)) * dpsi * kdk ;
c$$$      which equals: ;
c$$$      (pi*K^2) * K*cos(omega-phi) * ( jbess(1,2*pi*K*d) + jbess(3,2*pi*K*d) ) ;
      implicit none
      real *8 K ! real *8: radius of disc ;
      real *8 phi ! real *8: orientation of planar-function ;
      real *8 d ! real *8: frequency of plane-wave ;
      real *8 omega ! real *8: orientation of plane-wave ;
      real *8 I ! real *8: integral ;
      real *8 pi
      real *8 tK,j1,j3 !temporary: real *8: bessel function evaluations ;
      pi = 4.0d0*datan(1.0d0)
      tK = 2*pi*K*d
      call jbess(1.0d0,1,tK,j1)
      call jbess(3.0d0,1,tK,j3)
      I = pi*K**2 * K*cos(phi-omega) * ( j1 + j3 )
      I_xP_f = I
      end !function

      real *8 function I_P_f(K,d)
c$$$      calculates the integral of a plane-wave, ;
c$$$      where the domain of integration is the disc of radius K: ;
c$$$      \int_{0}^{K} \int_{0}^{2*pi} exp(+i*2*pi*k*d*cos(psi-omega)) * dpsi * kdk ;
c$$$      which equals: ;
c$$$      (pi*K^2) * ( jbess(0,2*pi*K*d) + jbess(2,2*pi*K*d) ) ;
      implicit none
      real *8 K ! real *8: radius of disc ;
      real *8 d ! real *8: frequency of plane-wave ;
      real *8 I ! real *8: integral ;
      real *8 pi
      real *8 tK,j0,j2 !temporary: real *8: bessel function evaluations ;
      pi = 4.0d0*datan(1.0d0)
      tK = 2*pi*K*d
      call jbess(0.0d0,1,tK,j0)
      call jbess(2.0d0,1,tK,j2)
      I = pi*K**2 * ( j0 + j2 )
      I_P_f = I
      end !function
