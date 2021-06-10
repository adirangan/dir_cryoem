      subroutine get_F2_x_c(max_x_c,X_x_c,Y_x_c,F_x_c,param_1,param_2)
      real *8 max_x_c,X_x_c,Y_x_c,param_1,param_2
      complex *16 F_x_c
      real *8 pi,R_x_c,W_x_c
      real *8 sx,mx,sy,my,tmp_a,tmp_b,tmp_c,tmp_d
      real *8 tmp_crop
      pi = 4*atan(1.0)
      R_x_c = dsqrt(X_x_c**2 + Y_x_c**2)
      W_x_c = atan2(Y_x_c,X_x_c)
      sx=0.1
      mx=0.2-1*0.5 + param_1*1
      call periodize_r8(mx,-0.5d0,+0.5d0,mx)
      sy=0.2
      my=0.2-1*0.5 + param_1*3
      call periodize_r8(my,-0.5d0,+0.5d0,my)
      tmp_a = 1/(2*pi*max_x_c*sx*max_x_c*sy)*exp(-(X_x_c - max_x_c *mx)
     $     **2/2/(max_x_c*sx)**2)*exp(-(Y_x_c - max_x_c*my)**2 /2
     $     /(max_x_c*sy)**2)
      sx=0.1
      mx=0.6-1*0.5 + param_1*5
      call periodize_r8(mx,-0.5d0,+0.5d0,mx)
      sy=0.1
      my=0.3-1*0.5 + param_1*7
      call periodize_r8(my,-0.5d0,+0.5d0,my)
      tmp_b = 1/(2*pi*max_x_c*sx*max_x_c*sy)*exp(-(X_x_c - max_x_c
     $     *mx)**2/2/(max_x_c*sx)**2)*exp(-(Y_x_c - max_x_c*my)**2
     $     /2 /(max_x_c*sy)**2)
      sx=0.1
      mx=0.7-1*0.5 + param_1*11
      call periodize_r8(mx,-0.5d0,+0.5d0,mx)
      sy=0.2
      my=0.8-1*0.5 + param_1*13
      call periodize_r8(my,-0.5d0,+0.5d0,my)
      tmp_c = 1/(2*pi*max_x_c*sx*max_x_c*sy)*exp(-(X_x_c - max_x_c
     $     *mx)**2/2/(max_x_c*sx)**2)*exp(-(Y_x_c - max_x_c*my)**2
     $     /2 /(max_x_c*sy)**2)
      sx=0.1
      mx=0.2-1*0.5 + param_1*17
      call periodize_r8(mx,-0.5d0,+0.5d0,mx)
      sy=0.2
      my=0.8-1*0.5 + param_1*19
      call periodize_r8(my,-0.5d0,+0.5d0,my)
      tmp_d = 1/(2*pi*max_x_c*sx*max_x_c*sy)*exp(-(X_x_c - max_x_c
     $     *mx)**2/2/(max_x_c*sx)**2)*exp(-(Y_x_c - max_x_c*my)**2
     $     /2 /(max_x_c*sy)**2)
      tmp_crop = 0.5*(1+erf(-16.0*(R_x_c-0.25*(1 + 0.5*cos(param_2
     $     *W_x_c)))))
      F_x_c = cmplx(tmp_crop*(tmp_a+tmp_b+tmp_c+tmp_d),0.0)
      end
