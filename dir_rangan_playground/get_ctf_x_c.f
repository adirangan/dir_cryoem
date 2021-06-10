!> Doxygen comment: ;\n
!> sets up a particular class of (very anisotropic) ctf functions. ;\n
      subroutine get_ctf_x_c(max_x_c,X_x_c,Y_x_c,C_x_c,param_0,param_1
     $     ,param_2,param_3,param_4)
c$$$      Inputs:
c$$$      real *8 max_x_c not used
c$$$      real *8 X_x_c x-value (x-space cartesian-grid)
c$$$      real *8 Y_x_c y-value (x-space cartesian-grid)
c$$$      real *8 param_0 sigma_E, e.g., 0.25
c$$$      real *8 param_1 omega_E, e.g., pi/4
c$$$      real *8 param_2 eccentricity, e.g., 10
c$$$      real *8 param_3 tau_G, e.g., 10
c$$$      real *8 param_4 omega_G, e.g., 5
c$$$      Outputs:
c$$$      complex *16 C_x_c value of ctf on x-space cartesian-grid
c$$$
c$$$      See following matlab code for ctf-function:
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$N = 128;
c$$$[X,Y] = meshgrid(linspace(-1,1,N));
c$$$R = dsqrt(X.^2+Y.^2);W = datan2(Y,X);
c$$$sigma_E = 1/4; eccentricity = 10;
c$$$omega_E = pi/4; cw = dcos(omega_E); sw = dsin(omega_E); 
c$$$rw = [cw +sw ; -sw cw] ; rm = [cw -sw ; +sw cw]; 
c$$$dw = (1/sigma_E).*[ 1 0 ; 0 dsqrt(eccentricity/1) ]; 
c$$$dm = (sigma_E/1).*[ 1 0 ; 0 dsqrt(1/eccentricity) ]; 
c$$$C = dw * rm ;
c$$$Xt = C(1).*X + C(3).*Y; Yt = C(2).*X + C(4).*Y; 
c$$$E1 = Xt.*Xt + Yt.*Yt;
c$$$C = rw * dw * dw * rm;
c$$$E2 = X.*C(1).*X + X.*C(3).*Y + Y.*C(2).*X + Y.*C(4).*Y; 
c$$$disp(sprintf(' %% norm E1-E2 %f',norm(E1-E2)));
c$$$tau_G = 10;
c$$$Rt = dsqrt(Xt.^2 + Yt.^2)/tau_G;
c$$$omega_G = 5;
c$$$G1 = dexp(-Rt).*dsin(min(pi,2.0d0*pi*omega_G*Rt.^2));
c$$$F = dexp(-E1).*G1;
c$$$%F = dexp(-E1);
c$$$K = recenter2(fft2(recenter2(F)));
c$$$subplot(2,2,1); imagesc(F); axis square; colorbar; 
c$$$subplot(2,2,2); imagesc(real(K)); axis square; colorbar;
c$$$subplot(2,2,3); imagesc(imag(K)); axis square; colorbar;
c$$$subplot(2,2,4); imagesc(abs(K)); axis square; colorbar;
c$$$%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c$$$      real *8 param_0 sigma_E
c$$$      real *8 param_1 omega_E
c$$$      real *8 param_2 eccentricity
c$$$      real *8 param_3 tau_G
c$$$      real *8 param_4 omega_G
      implicit none
      integer verbose
      data verbose / 0 /
      real *8 max_x_c,X_x_c,Y_x_c
      real *8 param_0,param_1,param_2,param_3,param_4
      complex *16 C_x_c
      real *8 pi,R_x_c,W_x_c
      real *8 omega_E,sigma_E,cw,sw,omega_G,tau_G,eccentricity
      real *8 rw00,rw10,rw01,rw11
      real *8 rm00,rm10,rm01,rm11
      real *8 dw00,dw10,dw01,dw11
      real *8 dm00,dm10,dm01,dm11
      real *8 cc00,cc10,cc01,cc11
      real *8 xt,yt,et,rt,gt
      if (verbose.gt.0) then
         write(6,'(A,F8.5)') 'param_0: ',param_0
         write(6,'(A,F8.5)') 'param_1: ',param_1
         write(6,'(A,F8.5)') 'param_2: ',param_2
         write(6,'(A,F8.5)') 'param_3: ',param_3
         write(6,'(A,F8.5)') 'param_4: ',param_4
      end if
      pi = 4*datan(1.0d0)
      R_x_c = dsqrt(X_x_c**2 + Y_x_c**2)
      W_x_c = datan2(Y_x_c,X_x_c)
      omega_E = param_0
      sigma_E = param_1
      eccentricity = param_2
      tau_G = param_3
      omega_G = param_4
      cw = dcos(omega_E)
      sw = dsin(omega_E)
      rw00 = cw
      rw01 = +sw
      rw10 = -sw
      rw11 = cw
      rm00 = cw
      rm01 = -sw
      rm10 = +sw
      rm11 = cw
      dw00 = 1 / (sigma_E)
      dw01 = 0 / (sigma_E)
      dw10 = 0 / (sigma_E)
      dw11 = dsqrt(eccentricity) / (sigma_E)
      dm00 = 1 * (sigma_E)
      dm01 = 0 * (sigma_E)
      dm10 = 0 * (sigma_E)
      dm11 = (sigma_E) / dsqrt(eccentricity)
      cc00 = dw00*rm00 + dw01*rm10
      cc01 = dw00*rm01 + dw01*rm11
      cc10 = dw10*rm00 + dw11*rm10
      cc11 = dw10*rm01 + dw11*rm11
      xt = cc00*X_x_c + cc01*Y_x_c
      yt = cc10*X_x_c + cc11*Y_x_c
      et = xt*xt + yt*yt
      rt = dsqrt(xt**2 + yt**2)/tau_G
      gt = dexp(-rt)*dsin(min(pi,2.0d0*pi*omega_G*rt**2))
      C_x_c = dcmplx(dexp(-et)*gt,0.0d0)
      if (verbose.gt.0) then
         write(6,*) 'dm ',dm00,dm01,dm10,dm11
         write(6,*) 'dw ',dw00,dw01,dw10,dw11
         write(6,*) 'rm ',rm00,rm01,rm10,rm11
         write(6,*) 'rw ',rw00,rw01,rw10,rw11
         write(6,*) 'cc ',cc00,cc01,cc10,cc11
         write(6,*) 'X_x_c: ',X_x_c,'Y_x_c: ',Y_x_c
         write(6,*) 'xt: ',xt, 'yt: ',yt
         write(6,*) 'C_x_c: ',C_x_c
      end if
      end
