%%%%%%%%;
% Paraxial approximation to shrodinger-equation. ;
% Beam headed in z-direction. ;
% Using polar (i.e., radial) description of x-y-plane for now. ;
%%%%%%%%;

flag_verbose=1; nf=0;
flag_disp = 1; nf=0; flag_replot = 1;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_besselbeam/dir_jpg',string_root);

%%%%%%%%;
% Simple test of J0 integral. ;
%%%%%%%%;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(w_));
I_quad = integral(g,0,2*pi);
I_form = besselj(0,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero.');
%%%%%%%%;

%%%%%%%%;
% Here we test a version of the integral in Gradshteyn and Ryzhik: 2.33 (just to make sure there are no typos). ;
%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(b^2/max(1e-12,4*a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Here we test another integral. ;
%%%%%%%%;
eta = 2.5; k = 0.1; x_min = -10*4; x_max = +10*4;
g0 = @(x_) exp(-i*(eta.^2*x_.^2 + 2*pi*k*x_));
%g0 = @(x_) exp(-i*(-eta.^2*x_.^2 + 2*pi*k*x_));
I_quad = integral(g0,x_min,x_max);
g1 = @(x_) exp(+i*(pi*k/eta)^2) * exp(-i*(eta*x_ + pi*k/eta).^2);
I_qua1 = integral(g1,x_min,x_max);
fnorm_disp(flag_verbose,'I_qua1',I_qua1,'I_quad',I_quad,' %<-- should be zero (change of variables).');
g2 = @(y_) exp(+i*(pi*k/eta)^2) * exp(-i*y_.^2);
y_min = x_min*eta - (pi*k/eta); y_max = x_max*eta - (pi*k/eta);
I_qua2 = integral(g2,y_min,y_max) / eta;
fnorm_disp(flag_verbose,'I_qua2',I_qua2,'I_quad',I_quad,' %<-- should be zero (change of variables).');
I_form = exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k/eta)^2); %<-- principal-value of oscillatory integral. ;
%I_form = i*exp(-i*pi/4)*sqrt(pi)/eta * exp(-i*(pi*k/eta)^2); %<-- principal-value of oscillatory integral. ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%x_ = linspace(x_min,x_max,1024*8); plot(x_,g0(x_),'o-');
%%%%%%%%;

%%%%%%%%;
% Consequently, the (full-space) integral of: ;
% exp(-i*(eta.^2*x_0_.^2 + 2*pi*k_0*x_0_)) .* exp(-i*(eta.^2*x_1_.^2 + 2*pi*k_1*x_1_)) dx_0 dx_1 ;
% --> exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k_0/eta)^2) .* exp(-i*pi/4)*sqrt(pi)/eta * exp(+i*(pi*k_1/eta)^2) ;
% --> exp(-i*pi/2) * pi/eta^2 * exp( +i * (pi/eta)^2 * (k_0^2 + k_1^2) ) ;
%%%%%%%%;
eta = 2.50; k = 0.1; x_min = 0; x_max = 19;
h0 = @(x_) x_ .* exp(-i*eta^2*x_.^2) .* besselj(0,2*pi*k*x_) * (2*pi);
I_quad = integral(h0,x_min,x_max);
%x_ = linspace(x_min,x_max,1024*8); plot(x_,h0(x_),'o-');
h1 = @(x_) +1/(2*-i*eta^2) * exp(-i*eta^2*x_.^2) .* besselj(0,2*pi*k*x_) * (2*pi);
h2 = @(x_) -1/(2*-i*eta^2) * exp(-i*eta^2*x_.^2) .* besselj(1,2*pi*k*x_) * (2*pi*k)*(2*pi);
I_qua2 = h1(x_max) - h1(x_min) - integral(h2,x_min,x_max);
%x_ = linspace(x_min,x_max,1024*8); plot(x_,h1(x_),'o-');
fnorm_disp(flag_verbose,'I_qua2',I_qua2,'I_quad',I_quad,' %<-- should be zero (integration by parts).');
I_form = exp(-i*pi/2) * pi/eta^2 * exp( +i * (pi/eta)^2 * (k^2) ) ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% This is basically the same as test_basis_3.m ;
% (using the same scale in fourier-space). ;
%%%%%%%%;
f_res = 16;
k_int = 48*f_res;
k_eq_d_double = 0.5;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64*f_res;
n_x_M_u = n_x_u_pack;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1);
x_u_0_ = transpose(x_u_0_(1:n_x_M_u));
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1);
x_u_1_ = transpose(x_u_1_(1:n_x_M_u));
[x_u_0_xx__,x_u_1_xx__] = ndgrid(x_u_0_,x_u_1_);
x_u_r_xx__ = sqrt(x_u_0_xx__.^2 + x_u_1_xx__.^2);
x_u_w_xx__ = atan2(x_u_1_xx__,x_u_0_xx__);
dx = diameter_x_c/n_x_M_u;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi);
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,k_p_r_max ...
,k_eq_d ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_p_r_)*4*pi - volume: %0.16f',sum(weight_3d_k_p_r_)*4*pi - 4/3*pi*k_p_r_max^3)); end;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
l_max_max = max(l_max_);
n_w_max = 2*(l_max_max+1);
template_k_eq_d = -1; %<-- default value of n_w_max. ;
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_wk_ ...
] = ...
get_weight_2d_1( ...
 max(0,flag_verbose-1) ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_k_p_r_) - area %0.16f',sum(weight_2d_k_p_r_) - pi*k_p_r_max^2)); end;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
k_p_r_wk_ = reshape(transpose(repmat(k_p_r_,[1,n_w_max])),[n_w_sum,1]);

%%%%%%%%;
% simple plane-wave incident on lens. ;
%%%%%%%%;
eta = 2.50; %<-- kz/(2*f) = z-frequency / (2*focal_length). ;
%%%%%%%%;
str_type = 'ringlens';
%%%%;
if strcmp(str_type,'quadlens');
M_k_p_r_ = exp(-i*pi/2) * pi/max(1e-12,eta^2) * exp( + i * (pi/max(1e-12,eta))^2 * k_p_r_.^2 );
M_k_p_wk_ = reshape(transpose(repmat(M_k_p_r_,[1,n_w_max])),[n_w_sum,1]);
end;%if strcmp(str_type,'quadlens');
%%%%;
if strcmp(str_type,'ringlens');
ring_inn = 0.35; ring_ext = 0.36; %<-- interior and exterior radii of annulus. ;
M_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
x_min = ring_inn; x_max = ring_ext;
h0 = @(x_) x_ .* exp(-i*eta^2*x_.^2) .* besselj(0,2*pi*k_p_r*x_) * (2*pi);
I_quad = integral(h0,x_min,x_max);
M_k_p_r_(1+nk_p_r) = I_quad;
end;%for nk_p_r=0:n_k_p_r-1;
M_k_p_wk_ = reshape(transpose(repmat(M_k_p_r_,[1,n_w_max])),[n_w_sum,1]);
end;%if strcmp(str_type,'quadlens');
%%%%;

%%%%%%%%;
% initial radial waveform. ;
%%%%%%%%;
flag_check=0;
if flag_check;
n_x_p_r = n_k_p_r; x_p_r_ = linspace(0,x_p_r_max,n_x_p_r);
M_x_p_r_ = zeros(n_x_p_r,1);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
tmp_I = sum(2*pi*besselj(0,2*pi*k_p_r_*x_p_r).*M_k_p_r_.*k_p_r_.*weight_2d_k_p_r_);
M_x_p_r_(1+nx_p_r) = tmp_I;
end;%for nx_p_r=0:n_x_p_r-1;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); cla; hold on;
plot(x_p_r_,abs(M_x_p_r_),'ko-');
plot(x_p_r_,real(M_x_p_r_),'r-');
plot(x_p_r_,imag(M_x_p_r_),'g-');
title('abs(M_x_p_r_)','Interpreter','none');
subplot(1,2,2); cla; hold on;
plot(k_p_r_,abs(M_k_p_r_),'ko-');
plot(k_p_r_,real(M_k_p_r_),'r-');
plot(k_p_r_,imag(M_k_p_r_),'g-');
title('abs(M_k_p_r_)','Interpreter','none');
disp('returning'); return;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% evolve radial waveform. ;
%%%%%%%%;
flag_check=1;
if flag_check;
%%%%;
n_x_p_r = 1024; x_p_r_ = linspace(0,x_p_r_max,n_x_p_r);
x_p_r_from_k_p_r_xk__ = zeros(n_x_p_r,n_k_p_r);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
x_p_r_from_k_p_r_xk__(1+nx_p_r,:) = 2*pi*besselj(0,2*pi*k_p_r_*x_p_r).*k_p_r_.*weight_2d_k_p_r_;
end;%for nx_p_r=0:n_x_p_r-1;
%%%%;
lambda = 1.0; %<-- here measure z in terms of lambda. ;
dzl = 1.0/8192.0;
zl_max = 128*dzl;
n_zl = zl_max/dzl;
P_k_p_r_ = exp(-i*pi*lambda*dzl*k_p_r_.^2); %<-- radial fourier propagator ;
M_k_p_r_kt__ = zeros(n_k_p_r,n_zl);
M_k_p_r_kt__(:,1+0) = M_k_p_r_;
for nzl=1:n_zl-1;
M_k_p_r_kt__(:,1+nzl+0) = M_k_p_r_kt__(:,1+nzl-1).*P_k_p_r_;
end;%for nzl=1:n_zl-1;
%%%%;
M_x_p_r_xt__ = x_p_r_from_k_p_r_xk__*M_k_p_r_kt__;
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot_{1} = subplot(1,2,1);
tmp_xt__ = abs(M_x_p_r_xt__);
tmp_xt__ = [flipud(tmp_xt__) ; tmp_xt__];
imagesc(tmp_xt__); set(gca,'ydir','normal'); axisnotick; title('abs(M_x_p_r_xt__)','Interpreter','none');
subplot_{2} = subplot(1,2,2);
tmp_xt__ = angle(M_x_p_r_xt__);
tmp_xt__ = [flipud(tmp_xt__) ; tmp_xt__];
imagesc(tmp_xt__,[-pi,+pi]); set(gca,'ydir','normal'); axisnotick; title('angle(M_x_p_r_xt__)','Interpreter','none');
colormap(subplot_{2},colormap('hsv'));
colormap(subplot_{1},colormap_81s());
%%%%;
disp('returning'); return;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
lambda = 1.0; %<-- here measure z in terms of lambda. ;
dzl = 1.0/64.0;
zl_max = 2.0;
n_zl = zl_max/dzl;
tmp_M_k_p_wk_ = M_k_p_wk_;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 18;
gamma_use = 2.0;
colormap(colormap_pm().^gamma_use);
clim_ = [];
%%%%;
for nzl=0:n_zl-1;
tmp_M_k_p_wk_ = tmp_M_k_p_wk_.*exp(-i*pi*lambda*dzl*k_p_r_wk_.^2); %<-- paraxial fourier propagator ;
tmp_N_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
if isempty(clim_); clim_ = [-1,+1]*prctile(abs(real(tmp_N_x_c_xx__)),100,'all'); end;
subplot(1,1,1); cla;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(tmp_N_x_c_xx__),clim_,colormap_pm.^gamma_use); axis image;
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$','Interpreter','latex');
drawnow();
end;%for nzl=0:n_zl-1;
%%%%;



