      subroutine get_ctf_k_c(max_k_c,X_k_c,Y_k_c,C_k_c,param_1,param_2
     $     ,param_3,param_4,param_5)
c$$$      real *8 param_1 sigma_E
c$$$      real *8 param_2 omega_E
c$$$      real *8 param_3 tau_G
c$$$      real *8 param_4 omega_G
c$$$      real *8 param_5 eccentricity
      real *8 max_k_c,X_k_c,Y_k_c
      real *8 param_1,param_2,param_3,param_4,param_5
      complex *16 C_k_c
      real *8 pi,R_k_c,W_k_c
      real *8 omega_E,sigma_E,cw,sw,omega_G,tau_G,eccentricity
      real *8 rw00,rw10,rw01,rw11
      real *8 rm00,rm10,rm01,rm11
      real *8 dw00,dw10,dw01,dw11
      real *8 dm00,dm10,dm01,dm11
      real *8 cc00,cc10,cc01,cc11
      real *8 xt,yt,et,rt,gt
      pi = 4*atan(1.0)
      R_k_c = dsqrt(X_k_c**2 + Y_k_c**2)
      W_k_c = atan2(Y_k_c,X_k_c)
      omega_E = param_1
      cw = cos(omega_E)
      sw = sin(omega_E)
      rw00 = cw
      rw01 = +sw
      rw10 = -sw
      rw11 = cw
      rm00 = cw
      rm01 = -sw
      rm10 = +sw
      rm11 = cw
      sigma_E = param_2
      eccentricity = param_5
      dw00 = 1 / (sigma_E)
      dw01 = 0 / (sigma_E)
      dw10 = 0 / (sigma_E)
      dw11 = sqrt(eccentricity) / (sigma_E)
      dm00 = 1 * (sigma_E)
      dm01 = 0 * (sigma_E)
      dm10 = 0 * (sigma_E)
      dm11 = (sigma_E) / sqrt(eccentricity)
      cc00 = dw00*rm00 + dw01*rm10
      cc01 = dw00*rm01 + dw01*rm11
      cc10 = dw10*rm00 + dw11*rm10
      cc11 = dw10*rm01 + dw11*rm11
      xt = c00*X_k_c + c01*Y_k_c
      yt = c10*X_k_c + c11*Y_k_c
      et = xt*xt + yt*yt
      rt = sqrt(xt**2 + yt**2)
      tau_G = param_3
      omega_G = param_4
      gt = exp(-rt/tau_G)*cos(2*pi*omega_G*rt)
      C_k_c = cmplx(exp(-et)*gt,0.0)
      end
