%%%%%%%%;
% Paraxial approximation to shrodinger-equation. ;
% Beam headed in z-direction. ;
% Using polar (i.e., radial) description of x-y-plane for now. ;
%%%%%%%%;

flag_verbose=1; nf=0;
flag_disp = 1; nf=0; flag_replot = 1;
str_thisfunction = 'test_bb_4';

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
% definition of erf. ;
%%%%%%%%;
sigma = 2.5;
u = 3;
I_form = 0.5*erf(u/(sqrt(2)*sigma));
I_quad = integral(@(x) (1/sqrt(2*pi))/sigma *exp(-x.^2/(2*sigma^2)),0,u);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero.');

flag_check=0;
if flag_check;
%%%%%%%%;
% test fourier-transform on chebyshev grid. ;
%%%%%%%%;
tolerance_nufft = 1e-12;
lx_max = 1;
lx_lim_ = lx_max*[-1,+1];
%%%%%%%%;
%n_lx = 256+1;
for n_lx = 1+2.^[5:10];
%%%%%%%%;
[node_lx_,weight_lx_] = chebpts(n_lx,lx_lim_);
node_lx_ = reshape(node_lx_,[n_lx,1]);
weight_lx_ = reshape(weight_lx_,[n_lx,1]);
k_val_max = round(0.5*n_lx);
k_val_ = transpose(-k_val_max:1:+k_val_max); n_k = numel(k_val_);
nk0 = efind(k_val_==0);
FF__ = zeros(n_k,n_k);
for nk=0:n_k-1;
if (flag_verbose>0); if (mod(nk,128)==0); disp(sprintf(' %% nk %.4d/%.4d',nk,n_k)); end; end;
k_val = k_val_(1+nk);
F_ = finufft1d3(node_lx_,exp(i*2*pi*k_val*node_lx_).*weight_lx_,-1,tolerance_nufft,2*pi*k_val_);
FF__(:,1+nk) = F_;
end;%for nk=0:n_k-1;
%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figmed; fig81s;
subplot(1,2,1);
imagesc(log(abs(FF__)),[-12,0]); axis image; axisnotick; colorbar;
title('FF__','Interpreter','none');
subplot(1,2,2); hold on;
plot(k_val_,log(abs(FF__(:,1+nk0))),'ko-');
nk1 = efind(k_val_==-round(n_lx/3)); plot(k_val_(1+nk1)*[1;1],[-16;+1],'r-','LineWidth',3);
nk1 = efind(k_val_==+round(n_lx/3)); plot(k_val_(1+nk1)*[1;1],[-16;+1],'r-','LineWidth',3);
end;%if flag_disp>1;
%%%%%%%%;
end;%for n_lx = 1+2.^[5,10];
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% repeat test, estimating typical k-threshold for spectral-accuracy. ;
% Note that the threshold readily decreases below 4, even for small/moderate values of n_lx. ;
% Note that the threshold approaches 3 or so for large values of n_lx. ;
%%%%%%%%;
tolerance_nufft = 1e-12;
lx_max = 1;
lx_lim_ = lx_max*[-1,+1];
%%%%%%%%;
n_lx_ = transpose(1+round(2.^[4:0.25:16])); n_n_lx = numel(n_lx_);
dnk_cut_ = zeros(n_n_lx,1);
k_cut_ = zeros(n_n_lx,1);
for nn_lx=0:n_n_lx-1;
n_lx = n_lx_(1+nn_lx);
[node_lx_,weight_lx_] = chebpts(n_lx,lx_lim_);
node_lx_ = reshape(node_lx_,[n_lx,1]);
weight_lx_ = reshape(weight_lx_,[n_lx,1]);
k_val_max = round(0.5*n_lx);
k_val_ = transpose(-k_val_max:1:+k_val_max); n_k = numel(k_val_);
nk0 = efind(k_val_==0);
F_ = zeros(n_k,1); k_val = k_val_(1+nk0);
F_ = finufft1d3(node_lx_,exp(i*2*pi*k_val*node_lx_).*weight_lx_,-1,tolerance_nufft,2*pi*k_val_);
nk1 = min(efind((k_val_>0) & (log(abs(F_))>-12)));
dnk_cut_(1+nn_lx) = nk1-nk0;
k_cut_(1+nn_lx) = k_val_(1+nk1) - k_val_(1+nk0);
end;%for nn_lx=0:n_n_lx-1;
%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figsml;
plot(1:n_n_lx,n_lx_./max(1,k_cut_),'o');
xlabel('n_n_lx','Interpreter','none');
ylabel('k_cut_./n_lx_','Interpreter','none');
end;%if flag_disp>1;
%%%%%%%%;
end;%if flag_check;

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
% Simple test of Jq integral (frwd fft). ;
%%%%%%%%;
for nq=0:8-1;
q = nq; omega = pi/3;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(omega - w_)) .* exp(-i*q*w_) ;
I_quad = integral(g,0,2*pi);
I_form = (-i).^q .* exp(-i*q*omega) .* besselj(q,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- q %d: error should be zero.',q));
end;%for nq=0:8-1;
%%%%%%%%;

%%%%%%%%;
% Simple test of Jq integral (back fft). ;
%%%%%%%%;
for nq=0:8-1;
q = nq; omega = pi/3;
kx = 1.75;
g = @(w_) exp(-i*2*pi*kx*cos(omega - w_)) .* exp(+i*q*w_) ;
I_quad = integral(g,0,2*pi);
I_form = (-i).^q .* exp(+i*q*omega) .* besselj(q,2*pi*kx)*(2*pi);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- q %d: error should be zero.',q));
end;%for nq=0:8-1;
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

%{
%%%%%%%%;
% Now testing integral from Talman 1978. ;
% Q_{q,p}(t) = \int_{0}^{\infty} r^{-p - i*t} J_{q}(r) dr ;
% Q_{q,p}(t) = 2^{-p-i*t} \frac{\Gamma(0.5(n-p-i*t+1))}{\Gamma(0.5(n+p+i*t+1))} ;
% with -1/2 < p < q + 1 ;
%%%%%%%%;
%%%%%%%%;
for q_val = -3:+3;
p_min = -1/2; p_max = q_val + 0.5;
p_val_ = p_min:0.5:p_max;
for p_val = p_val_;
t_val = randn();
g = @(r_) r_.^(-p_val-i*t_val) .* besselj(q_val,r_) ;
I_quad = integral(g,0,+Inf);
I_form = exp((-p_val-i*t_val)*log(2) + gamma_godfrey_1(0.5*(q_val-p_val-i*t_val+1)) - gamma_godfrey_1(0.5*(q_val+p_val+i*t_val+1)));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- p_val %0.2f q_val %.2d: error should be zero.',p_val,q_val));
end;%for p_val = p_val_;
end;%for q_val = -3:+3;
%%%%%%%%;
 %}

%{
%%%%%%%%;
% Now testing integral from 6.561.14 in TISPISGIMR. ;
% Q_{mu_r,nu}(mu_i) = \int_{0}^{\infty} r^{mu_r + i*mu_i} J_{nu}(a*r) dr ;
% Q_{mu_r,nu}(mu_i) = 2.^{mu}.*a.^{1+mu} .* \frac{\Gamma(0.5(1+nu+mu))}{\Gamma(0.5(1+nu-mu))} ;
% with -\Re(nu)-1 < mu_r < 1/2 ;
%%%%%%%%;
a = 2*pi;
for q_val = -3:+3;
if abs(q_val)==0; mu_r = -0.75; end;
if abs(q_val)> 0; mu_r = -1.00; end;
nu = q_val;
for mu_i = -2:+2;
mu = mu_r + i*mu_i;
%nu_sgn = +1; if (nu<0); nu_sgn = (-1).^nu; end;
g = @(r_) r_.^(mu) .* besselj(nu,a*r_) ;
I_quad = integral(g,0,+Inf,'RelTol',0,'AbsTol',1e-2);
I_form = (-1).^(nu .* (nu<0)) .* exp( mu*log(2) - (1+mu)*log(a) + gamma_godfrey_1(0.5*(1 + abs(nu) + mu)) - gamma_godfrey_1(0.5*(1 + abs(nu) - mu)) ) ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- mu_r %+0.2f mu_i %+0.2f q_val %+.2d: error should be zero.',mu_r,mu_i,q_val));
end;%for mu_i = -2:+2;
end;%for q_val = -3:+3;
%%%%%%%%;
Q = @(q_val,mu_val) (-1).^(q_val .* (q_val<0)) .* exp( mu_val*log(2) - (1+mu_val)*log(2*pi) + gamma_godfrey_1(0.5*(1 + abs(q_val) + mu_val)) - gamma_godfrey_1(0.5*(1 + abs(q_val) - mu_val)) ) ;
 %}

%{
%%%%%%%%;
% Now testing integral from 6.561.14 in TISPISGIMR, but this time with approximation to bessel-function from 8.402 and 3.326.2 in TISPISGIMR. ;
% Q_{mu_r,nu}(mu_i) = \int_{0}^{\infty} r^{mu_r + i*mu_i} J_{nu}(a*r) dr ;
% Q_{mu_r,nu}(mu_i) = 2.^{mu}.a.^{1+mu} .frac{\Gamma(0.5(1+nu+mu))}{\Gamma(0.5(1+nu-mu))} ;
% with -Re(nu)-1 < mu_r < 1/2 ;
% J_{nu}(z) \approx (z/2)^{nu} * (1/Gamma(1+nu) - z^{2}/(4Gamma(2+nu)) + z^{4}/(32Gamma(3+nu))) * exp(-\beta z^{4}) ;
% \int_{0}^{\infty} x^{mu} exp(-\beta x^{nu}) dx = Gamma(\delta) / (nu*\beta^{\delta}), with \delta = (mu+1)/nu . ;
%%%%%%%%;
b = 1.0; nu = 4.0; a = 2*pi;
for mu_r=[-0.5,0,1,2];
mu_i = 2*pi*rand(); mu = mu_r + i*mu_i;
g = @(r_) (a*r_).^(mu) .* exp(-b*(a*r_).^(nu)) ;
I_quad = integral(g,0,+Inf,'RelTol',0,'AbsTol',1e-12);
d = (mu+1)/nu;
I_form = exp(gamma_godfrey_1(d)) / max(1e-12,nu.*b.^(d)) / a ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- mu_r %+0.2f mu_i %+0.2f: error should be zero.',mu_r,mu_i));
end;%for mu_r=[-0.5,0,1,2];
%%%%%%%%;
J_approx = @(q_val,x_p_r_) (-1).^(q_val .* (q_val<0)) .* (x_p_r_/2).^(abs(q_val)) .* (1/gamma(1+abs(q_val)) - x_p_r_.^2/(4*gamma(2+abs(q_val))) + x_p_r_.^4/(16*2*gamma(3+abs(q_val)))) .* exp(-b.*x_p_r_.^6) ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
hold on;
for q_val=-6:+6;
plot(log(x_p_r_),besselj(q_val,x_p_r_),'k.-');
plot(log(x_p_r_),J_approx(q_val,x_p_r_),'ro-');
end;%for q_val=-6:+6;
%%%%%%%%;
a = 2*pi;
for q_val = -3:+3;
if abs(q_val)==0; mu_r = -0.75+0.5; end;
if abs(q_val)> 0; mu_r = -1.00+0.5; end;
nu = q_val;
for mu_i = -2:+2;
mu = mu_r + i*mu_i;
g_appr = @(r_) r_.^(mu) .* J_approx(q_val,a*r_) ;
I_quad = integral(g_appr,0,+Inf,'RelTol',0,'AbsTol',1e-12);
R = @(mu,nu) exp(gamma_godfrey_1((mu+1)./nu)) ./ max(1e-12,nu.*b.^((mu+1)./nu)) ./ a ;
I_form = ...
 (-1).^(q_val .* (q_val<0)) ...
.* a.^(-mu) .* (1/2).^(abs(q_val)) ...
.* ( ...
     +(1/( 1*gamma(1+abs(q_val)))) .* R(mu+abs(q_val)+0,6) ...
     -(1/( 4*gamma(2+abs(q_val)))) .* R(mu+abs(q_val)+2,6) ...
     +(1/(32*gamma(3+abs(q_val)))) .* R(mu+abs(q_val)+4,6) ...
     ) ...
;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,sprintf(' %%<-- mu_r %+0.2f mu_i %+0.2f q_val %+.2d: error should be zero.',mu_r,mu_i,q_val));
end;%for mu_i = -2:+2;
end;%for q_val = -3:+3;
%%%%%%%%;
%}

%%%%%%%%;
% Set up spatial-grid for testing. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
%%%%%%%%;
n_x_M_u = 128;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u+1)); x_c_0_ = x_c_0_(1:n_x_M_u);
x_c_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u+1)); x_c_1_ = x_c_1_(1:n_x_M_u);
%x_c_0_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u));
%x_c_1_ = transpose(linspace(-half_diameter_x_c,+half_diameter_x_c,n_x_M_u));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
dx = mean(diff(x_c_0_));
n_xx_M_u = n_x_M_u^2;

%%%%%%%%;
% Now set up some polar grids in real-space and fourier-space. ;
%%%%%%%%;
x_res = 2;
x_int = 48*x_res;
x_eq_d_double = 0.5*(x_res/2);
n_t_int = 1;
%%%%%%%%;
k_res = 2; %<-- be careful with the aliasing associated with n_x_M_u. ;
k_int = 48*k_res;
k_eq_d_double = 0.5*(k_res/2);
n_w_int = n_t_int;
%%%%%%%%;
x_p_r_max = half_diameter_x_c;
k_p_r_max = k_int/(2*pi);
x_eq_d = x_eq_d_double/(2*pi)/16 * 1.05; %<-- to check dimensions. ;
k_eq_d = k_eq_d_double/(2*pi)/ 8 * 1.00; %<-- to check dimensions. ;
%%%%%%%%;

%%%%%%%%;
[ ...
 n_x_p_r ...
,x_p_r_ ...
,weight_3d_x_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,x_p_r_max ...
,x_eq_d ...
);
fnorm_disp(flag_verbose,'sum(weight_3d_x_p_r_)*4*pi',sum(weight_3d_x_p_r_)*4*pi,'volume',4/3*pi*x_p_r_max^3,' %<-- should be small.');
%%%%%%%%;
n_t_max = n_t_int*2*(x_int+1) + 2;  %<-- to check dimensions. ;
n_t_0in_ = n_t_max*ones(n_x_p_r,1);
%%%%%%%%;
[ ...
 n_t_ ...
,weight_2d_x_p_r_ ...
,weight_2d_x_p_tx_ ...
,x_p_r_tx_ ...
,x_p_t_tx_ ...
,x_c_0_tx_ ...
,x_c_1_tx_ ...
] = ...
get_weight_2d_2( ...
 max(0,flag_verbose-1) ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,-1 ...
,n_t_0in_ ...
,weight_3d_x_p_r_ ...
);
fnorm_disp(flag_verbose,'sum(weight_2d_x_p_r_)',sum(weight_2d_x_p_r_),'area',pi*x_p_r_max^2,' %<-- should be small.');
n_t_sum = sum(n_t_); n_t_csum_ = cumsum([0;n_t_]);
%%%%%%%%;

%%%%%%%%;
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
fnorm_disp(flag_verbose,'sum(weight_3d_k_p_r_)*4*pi',sum(weight_3d_k_p_r_)*4*pi,'volume',4/3*pi*k_p_r_max^3,' %<-- should be small.');
%%%%%%%%;
%n_w_max = 12; 
n_w_max = n_w_int*2*(k_int+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 max(0,flag_verbose-1) ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
fnorm_disp(flag_verbose,'sum(weight_2d_k_p_r_)',sum(weight_2d_k_p_r_),'area',pi*k_p_r_max^2,' %<-- should be small.');
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%;
% Now test gaussian. ;
%%%%%%%%;
sigma_x_c = 0.0625;
sigma_k_p = 1/sigma_x_c;
delta_ = 0.25*0.75*[+0.1,-0.2];
M_x_c_form_ = 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0__-delta_(1+0)).^2 + (x_c_1__-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
M_x_c_form_l2 = sum(M_x_c_form_.^2,'all')*dx^2;
disp(sprintf(' %% sum(M_x_c_form_*dx^2,''all'') = %0.16f',sum(M_x_c_form_*dx^2,'all')));
disp(sprintf(' %% M_x_c_form_l2 = %0.16f',M_x_c_form_l2));
M_x_p_form_ = 1/(sqrt(2*pi)*sigma_x_c)^2 * exp( -( (x_c_0_tx_-delta_(1+0)).^2 + (x_c_1_tx_-delta_(1+1)).^2 ) / (2*sigma_x_c^2) );
M_x_p_form_l2 = sum(abs(M_x_p_form_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_p_form_l2 = %0.16f',M_x_p_form_l2));
M_k_p_quad_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_form_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2;
M_k_p_quad_l2 = sum(abs(M_k_p_quad_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_quad_l2 = %0.16f',M_k_p_quad_l2));
M_k_p_form_ = exp( -( (2*pi*k_c_0_wk_).^2 + (2*pi*k_c_1_wk_).^2 ) / (2/sigma_x_c^2) ) .* exp( - 2*pi*i*( k_c_0_wk_*delta_(1+0) + k_c_1_wk_*delta_(1+1) ) );
M_k_p_form_l2 = sum(abs(M_k_p_form_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_form_l2 = %0.16f',M_k_p_form_l2));
fnorm_disp(flag_verbose,'M_k_p_form_',M_k_p_form_,'M_k_p_quad_',M_k_p_quad_,' %<-- should be small if no aliasing');
M_k_p_reco_ = interp_x_p_to_k_p_xxnufft(n_k_p_r,k_p_r_,k_p_r_max,n_w_,weight_2d_k_p_wk_,n_x_p_r,x_p_r_,x_p_r_max,n_t_,weight_2d_x_p_tx_,M_x_p_form_)*(2*pi)^2;
M_k_p_reco_l2 = sum(abs(M_k_p_reco_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_p_reco_l2 = %0.16f',M_k_p_reco_l2));
M_x_c_reco_ = interp_k_p_to_x_c_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_quad_.*weight_2d_k_p_wk_*(2*pi)^2)*sqrt(n_x_M_u^2) * n_w_sum;
M_x_c_reco_l2 = sum(abs(M_x_c_reco_).^2,'all')*dx^2;
disp(sprintf(' %% M_x_c_reco_l2 = %0.16f',M_x_c_reco_l2));
fnorm_disp(flag_verbose,'M_x_c_form_',M_x_c_form_,'M_x_c_reco_',M_x_c_reco_,' %<-- should be small if no aliasing');
M_x_p_reco_ = interp_k_p_to_x_p_xxnufft(n_x_p_r,x_p_r_,x_p_r_max,n_t_,weight_2d_x_p_tx_,n_k_p_r,k_p_r_,k_p_r_max,n_w_,weight_2d_k_p_wk_,M_k_p_form_)*(2*pi)^2;
M_x_p_reco_l2 = sum(abs(M_x_p_reco_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_p_reco_l2 = %0.16f',M_x_p_reco_l2));

%%%%%%%%;
% Now test fourier-bessel representation. ;
%%%%%%%%;
M_x_q_quad_ = interp_p_to_q(n_x_p_r,n_t_,n_t_sum,M_x_p_form_);
M_x_q_quad_l2 = sum(abs(M_x_q_quad_).^2 .* weight_2d_x_p_tx_) * (2*pi)^2;
disp(sprintf(' %% M_x_q_quad_l2 = %0.16f',M_x_q_quad_l2));
M_k_q_quad_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_form_);
M_k_q_quad_l2 = sum(abs(M_k_q_quad_).^2 .* weight_2d_k_p_wk_) * (2*pi)^2;
disp(sprintf(' %% M_k_q_quad_l2 = %0.16f',M_k_q_quad_l2));
%%%%%%%%;
%imagesc_q(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_q_quad_));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_q_quad_));
%%%%%%%%;
M_x_q_quad__ = reshape(M_x_q_quad_,[n_t_max,n_x_p_r]);
M_k_q_quad__ = reshape(M_k_q_quad_,[n_w_max,n_k_p_r]);
%%%%%%%%;
for q_val=-6:+6;
nt = q_val; if (nt<0); nt=nt+n_t_max; end;
nw = q_val; if (nw<0); nw=nw+n_w_max; end;
M_x_q_x_ = reshape(M_x_q_quad__(1+nt,:),[n_x_p_r,1]);
M_k_q_k_ = reshape(M_k_q_quad__(1+nw,:),[n_k_p_r,1]);
N_x_q_x_ = zeros(n_x_p_r,1);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
N_x_q_x_(1+nx_p_r) = sum( (+i).^q_val .* besselj(q_val,2*pi*x_p_r*k_p_r_) .* M_k_q_k_ .* weight_2d_k_p_r_ ); %<-- note that weight_2d_k_p_r_ accounts for radial weighting. ;
end;%for nx_p_r=0:n_x_p_r-1;
fnorm_disp(flag_verbose,'M_x_q_x_',M_x_q_x_,'N_x_q_x_',N_x_q_x_,sprintf(' %%<-- q_val %+d: should be small.',q_val));
N_k_q_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
N_k_q_k_(1+nk_p_r) = sum( (-i).^q_val .* besselj(q_val,2*pi*k_p_r*x_p_r_) .* M_x_q_x_ .* weight_2d_x_p_r_ ); %<-- note that weight_2d_x_p_r_ accounts for radial weighting. ;
end;%for nk_p_r=0:n_k_p_r-1;
fnorm_disp(flag_verbose,'M_k_q_k_',M_k_q_k_,'N_k_q_k_',N_k_q_k_,sprintf(' %%<-- q_val %+d: should be small.',q_val));
end;%for q_val=-6:+6;
%%%%%%%%;

%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;figbig;fig80s;
p_row = 2; p_col = 4; np=0;
M_x_p_lim_ = reshape(prctile(real(M_x_p_form_),[  0,100],'all'),[1,2]);
M_k_p_lim_ = reshape(prctile(real(M_k_p_form_),[  0,100],'all'),[1,2]);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_form_),M_x_p_lim_);axis image;axisnotick; title('M_x_c_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_p_form_),M_x_p_lim_);axis image;axisnotick; title('M_x_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_x_p_r,x_p_r_,n_t_,n_t_sum,real(M_x_p_reco_),M_x_p_lim_);axis image;axisnotick; title('M_x_p_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_quad_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_quad_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_form_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_form_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_reco_),M_k_p_lim_);axis image;axisnotick; title('M_k_p_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_),M_x_p_lim_);axis image;axisnotick; title('M_x_c_reco_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_M_u,x_c_0_,n_x_M_u,x_c_1_,real(M_x_c_reco_-M_x_c_form_),M_x_p_lim_);axis image;axisnotick; title('difference','Interpreter','none');
disp(sprintf(' %% returning')); return;
end;%if flag_disp>1;
%%%%%%%%;

%%%%%%%%;
% Now try exponential change of variables. ;
%%%%%%%%;
alpha = 1.0; tmp_b = 1.0;
q_val = 0;
Q = @(q_val,mu_val) (-1).^(q_val .* (q_val<0)) .* exp( mu_val*log(2) - (1+mu_val)*log(2*pi) + gamma_godfrey_1(0.5*(1 + abs(q_val) + mu_val)) - gamma_godfrey_1(0.5*(1 + abs(q_val) - mu_val)) ) ;
J_approx = @(q_val,x_p_r_) (-1).^(q_val .* (q_val<0)) .* (x_p_r_/2).^(abs(q_val)) .* (1/gamma(1+abs(q_val)) - x_p_r_.^2/(4*gamma(2+abs(q_val))) + x_p_r_.^4/(16*2*gamma(3+abs(q_val)))) .* exp(-tmp_b.*x_p_r_.^6) ;
R = @(mu,nu) exp(gamma_godfrey_1((mu+1)./nu)) ./ max(1e-12,nu.*tmp_b.^((mu+1)./nu)) ./ (2*pi) ;
S = @(q_val,mu_val) ...
 (-1).^(q_val .* (q_val<0)) ...
.* (2*pi).^(-mu_val) .* (1/2).^(abs(q_val)) ...
.* ( ...
     +(1/( 1*exp(gamma_godfrey_1(1+abs(q_val)))) .* R(mu_val+abs(q_val)+0,6)) ...
     -(1/( 4*exp(gamma_godfrey_1(2+abs(q_val)))) .* R(mu_val+abs(q_val)+2,6)) ...
     +(1/(32*exp(gamma_godfrey_1(3+abs(q_val)))) .* R(mu_val+abs(q_val)+4,6)) ...
     ) ...
;
%%%%%%%%;
nt = q_val; if (nt<0); nt=nt+n_t_max; end;
nw = q_val; if (nw<0); nw=nw+n_w_max; end;
M_x_q_x_ = reshape(M_x_q_quad__(1+nt,:),[n_x_p_r,1]);
M_k_q_k_ = reshape(M_k_q_quad__(1+nw,:),[n_k_p_r,1]);
%%%%%%%%;

%%%%%%%%;
% calculate M-integral. ;
%%%%%%%%;
lk_lim_ = [log(min(k_p_r_)),log(max(k_p_r_))]/alpha; lk_dia = diff(lk_lim_);
lx_lim_ = [log(min(x_p_r_)),log(max(x_p_r_))]/alpha; lx_dia = diff(lx_lim_);
nyquist_factor = 2.0;
n_gamma_z_pre = max(1024+1,n_x_p_r);
gamma_z_max = max(1,n_gamma_z_pre)/max(1e-12,diff(lx_lim_))*nyquist_factor; %<-- first ensure that gamma_z_max is sufficient to resolve lx_. ;
n_gamma_z_pos = round(gamma_z_max * max(1e-12,diff(lk_lim_)) / nyquist_factor) ; %<-- now ensure that n_gamma_z is sufficient to resolve lk_. ;
n_gamma_z = max(n_gamma_z_pre,n_gamma_z_pos) ;
if (flag_verbose>0); disp(sprintf(' %% n_gamma_z_pre %d gamma_z_max %0.2f n_gamma_z_pos %d ',n_gamma_z_pre,gamma_z_max,n_gamma_z_pos)); end;
%%%%;
if (gamma_z_max >= 1.05*n_gamma_z/max(1e-12,lx_dia)*nyquist_factor);
disp(sprintf(' %% Warning, lx_dia %0.2f vs n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor %0.2f in %s',lx_dia,n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor,str_thisfunction));
end;%if (gamma_z_max >= n_gamma_z/max(1e-12,lx_dia)*nyquist_factor);
%%%%;
if (gamma_z_max >= 1.05*n_gamma_z/max(1e-12,lk_dia)*nyquist_factor);
disp(sprintf(' %% Warning, lk_dia %0.2f vs n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor %0.2f in %s',lk_dia,n_gamma_z/max(1e-12,gamma_z_max)*nyquist_factor,str_thisfunction));
end;%if (gamma_z_max >= n_gamma_z/max(1e-12,lk_dia)*nyquist_factor);
%%%%%%%%;
[gamma_z_,weight_1d_gamma_z_] = chebpts(n_gamma_z,gamma_z_max*[-1,+1]);
gamma_z_ = reshape(gamma_z_,[n_gamma_z,1]);
weight_1d_gamma_z_ = reshape(weight_1d_gamma_z_,[n_gamma_z,1]);
%%%%%%%%;
if abs(q_val)==0; n_val = -0.75; n_va2 = -0.75 + 1*0.50; end;
if abs(q_val)> 0; n_val = -1.00; n_va2 = -1.00 + 1*0.50; end;
%%%%%%%%;
fB_nuff_z_ = finufft1d3(log(x_p_r_)/alpha,M_x_q_x_.*x_p_r_.^(n_val).*weight_2d_x_p_r_,-1,1e-12,gamma_z_);
fB_nuf2_z_ = finufft1d3(log(x_p_r_)/alpha,M_x_q_x_.*x_p_r_.^(n_va2).*weight_2d_x_p_r_,-1,1e-12,gamma_z_);
%%%%%%%%;
flag_check=0;
if flag_check;
fB_quad_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_.*x_p_r_.^(n_val).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_quad_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuff_z_',fB_nuff_z_,'fB_quad_z_',fB_quad_z_,'%<-- should be zero');
fB_qua2_z_ = zeros(n_gamma_z,1);
for ngamma_z=0:n_gamma_z-1;
gamma_z = gamma_z_(1+ngamma_z);
fB = sum(M_x_q_x_.*x_p_r_.^(n_va2).*x_p_r_.^(-i*gamma_z/alpha).*weight_2d_x_p_r_);
fB_qua2_z_(1+ngamma_z) = fB;
end;%for ngamma_z=0:n_gamma_z-1;
fnorm_disp(flag_verbose,'fB_nuf2_z_',fB_nuf2_z_,'fB_qua2_z_',fB_qua2_z_,'%<-- should be zero');
end;%flag_check;
%%%%%%%%;
% Calculate J-integral. ;
%%%%%%%%;
fA_form_z_ = (-i).^(q_val) / alpha .* Q(q_val,-(1+n_val+i*gamma_z_/alpha));
fA_for2_z_ = (-i).^(q_val) / alpha .* S(q_val,-(1+n_va2+i*gamma_z_/alpha));
%%%%%%%%;
fC_comb_z_ = fA_form_z_.*flipud(fB_nuff_z_);
bC_comb_y_ = finufft1d3(gamma_z_,fC_comb_z_.*weight_1d_gamma_z_,+1,1e-12,log(k_p_r_)/alpha);
O_k_q_k_ = k_p_r_.^n_val .* bC_comb_y_  / (2*pi) ;
fC_com2_z_ = fA_for2_z_.*flipud(fB_nuf2_z_);
bC_com2_y_ = finufft1d3(gamma_z_,fC_com2_z_.*weight_1d_gamma_z_,+1,1e-12,log(k_p_r_)/alpha);
P_k_q_k_ = k_p_r_.^n_va2 .* bC_com2_y_  / (2*pi) ;
%%%%%%%%;
tmp_O_vs_M = fnorm_disp(0*flag_verbose,'O_k_q_k_',O_k_q_k_,'M_k_q_k_',M_k_q_k_,'%<-- should be small.');
disp(sprintf(' %% n_gamma_z %.6d gamma_z_max %6.2f tmp_O_vs_M %+0.6f',n_gamma_z,gamma_z_max,tmp_O_vs_M));
tmp_P_vs_M = fnorm_disp(0*flag_verbose,'P_k_q_k_',P_k_q_k_,'M_k_q_k_',M_k_q_k_,'%<-- should be small.');
disp(sprintf(' %% n_gamma_z %.6d gamma_z_max %6.2f tmp_P_vs_M %+0.6f',n_gamma_z,gamma_z_max,tmp_P_vs_M));
%%%%%%%%;
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(2,2,1); plot(log(k_p_r_),real(O_k_q_k_),'ko-',log(k_p_r_),real(M_k_q_k_),'rx-'); title('real(O_k_q_k_)','Interpreter','none');
subplot(2,2,2); plot(log(k_p_r_),imag(O_k_q_k_),'ko-',log(k_p_r_),imag(M_k_q_k_),'rx-'); title('imag(O_k_q_k_)','Interpreter','none');
subplot(2,2,3); plot(log(k_p_r_),real(P_k_q_k_),'ko-',log(k_p_r_),real(M_k_q_k_),'rx-'); title('real(P_k_q_k_)','Interpreter','none');
subplot(2,2,4); plot(log(k_p_r_),imag(P_k_q_k_),'ko-',log(k_p_r_),imag(M_k_q_k_),'rx-'); title('imag(P_k_q_k_)','Interpreter','none');
end;%if flag_disp>1;
%%%%%%%%;
%end;%for gamma_z_max = 50*[1,2,3]; end;%for n_gamma_z = 1024*[1,2,4,8];
%%%%%%%%;
 
%%%%%%%%;
% Compare with fht_2. ;
%%%%%%%%;
flag_sign = -1;
[ ...
 ~ ...
,L_k_q_k_ ...
] = ...
fht_2( ...
 [] ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_ ...
,q_val ...
,flag_sign ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_ ...
);
fnorm_disp(flag_verbose,'L_k_q_k_',L_k_q_k_,'M_k_q_k_',M_k_q_k_,'%<-- should be small.');
%%%%%%%%;
flag_sign = +1;
[ ...
 ~ ...
,L_x_q_x_ ...
] = ...
fht_2( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_ ...
,q_val ...
,flag_sign ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_ ...
);
fnorm_disp(flag_verbose,'L_x_q_x_',L_x_q_x_,'M_x_q_x_',M_x_q_x_,'%<-- should be small.');
if flag_disp>1;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(2,1,1); plot(x_p_r_,real(L_x_q_x_),'ko-',x_p_r_,real(M_x_q_x_),'rx-'); title('real');
subplot(2,1,2); plot(x_p_r_,imag(L_x_q_x_),'ko-',x_p_r_,imag(M_x_q_x_),'rx-'); title('imag');
end;%if flag_disp>1;
%%%%%%%%;

%%%%%%%%;
N_x_q_x_ = zeros(n_x_p_r,1);
for nx_p_r=0:n_x_p_r-1;
x_p_r = x_p_r_(1+nx_p_r);
N_x_q_x_(1+nx_p_r) = sum( (+i).^q_val .* besselj(q_val,2*pi*x_p_r*k_p_r_) .* M_k_q_k_ .* weight_2d_k_p_r_ ); %<-- note that weight_2d_k_p_r_ accounts for radial weighting. ;
end;%for nx_p_r=0:n_x_p_r-1;
fnorm_disp(flag_verbose,'M_x_q_x_',M_x_q_x_,'N_x_q_x_',N_x_q_x_,sprintf(' %%<-- q_val %+d: should be small.',q_val));
N_k_q_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
N_k_q_k_(1+nk_p_r) = sum( (-i).^q_val .* besselj(q_val,2*pi*k_p_r*x_p_r_) .* M_x_q_x_ .* weight_2d_x_p_r_ ); %<-- note that weight_2d_x_p_r_ accounts for radial weighting. ;
end;%for nk_p_r=0:n_k_p_r-1;
fnorm_disp(flag_verbose,'M_k_q_k_',M_k_q_k_,'N_k_q_k_',N_k_q_k_,sprintf(' %%<-- q_val %+d: should be small.',q_val));
%%%%%%%%;

%%%%%%%%;
% test rest of q_val for fht_2. ;
%%%%%%%%;
for q_val = -6:+6;
if (flag_verbose>0); disp(sprintf(' %% q_val %+d',q_val)); end;
nt = q_val; if (nt<0); nt=nt+n_t_max; end;
nw = q_val; if (nw<0); nw=nw+n_w_max; end;
M_x_q_x_ = reshape(M_x_q_quad__(1+nt,:),[n_x_p_r,1]);
M_k_q_k_ = reshape(M_k_q_quad__(1+nw,:),[n_k_p_r,1]);
%%%%%%%%;
flag_sign = -1;
[ ...
 ~ ...
,L_k_q_k_ ...
] = ...
fht_2( ...
 [] ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_ ...
,q_val ...
,flag_sign ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
);
fnorm_disp(flag_verbose,'L_k_q_k_',L_k_q_k_,'M_k_q_k_',M_k_q_k_,'%<-- should be small.');
%%%%%%%%;
flag_sign = +1;
[ ...
 ~ ...
,L_x_q_x_ ...
] = ...
fht_2( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_ ...
,q_val ...
,flag_sign ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
);
fnorm_disp(flag_verbose,'L_x_q_x_',L_x_q_x_,'M_x_q_x_',M_x_q_x_,'%<-- should be small.');
%%%%%%%%;
end;%for q_val = -6:+6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp('returning'); return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


flag_verbose=1;
flag_disp=1;nf=0;
%%%%%%%%;
% Redefine grids. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
%%%%%%%%;
x_res = 2;
x_int = 48*x_res;
x_eq_d_double = 0.5*(x_res/2);
%%%%%%%%;
k_res = 2*8*2; %<-- be careful with the aliasing associated with n_x_M_u. ;
k_int = 48*k_res;
k_eq_d_double = 0.5*(k_res/2);
%%%%%%%%;
x_p_r_max = half_diameter_x_c;
k_p_r_max = k_int/(2*pi);
x_eq_d = x_eq_d_double/(2*pi)/32 * 1.05; %<-- to check dimensions. ;
k_eq_d = k_eq_d_double/(2*pi)/ 4 * 1.00; %<-- to check dimensions. ;
%%%%%%%%;
[ ...
 n_x_p_r ...
,x_p_r_ ...
,weight_3d_x_p_r_ ...
,weight_2d_x_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,x_p_r_max ...
,x_eq_d ...
);
fnorm_disp(flag_verbose,'sum(weight_3d_x_p_r_)*4*pi',sum(weight_3d_x_p_r_)*4*pi,'volume',4/3*pi*x_p_r_max^3,' %<-- should be small.');
fnorm_disp(flag_verbose,'sum(weight_2d_x_p_r_)',sum(weight_2d_x_p_r_),'area',pi*x_p_r_max^2,' %<-- should be small.');
%%%%%%%%;
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,k_p_r_max ...
,k_eq_d ...
);
fnorm_disp(flag_verbose,'sum(weight_3d_k_p_r_)*4*pi',sum(weight_3d_k_p_r_)*4*pi,'volume',4/3*pi*k_p_r_max^3,' %<-- should be small.');
fnorm_disp(flag_verbose,'sum(weight_2d_k_p_r_)',sum(weight_2d_k_p_r_),'area',pi*k_p_r_max^2,' %<-- should be small.');
%%%%%%%%;

%%%%%%%%;
% set units. ;
%%%%%%%%;
unit = struct('type','unit');
unit.m = 1.0;
unit.cm=1e-2*unit.m;
unit.mm=1e-3*unit.m;
unit.um=1e-6*unit.m;
unit.nm=1e-9*unit.m;
unit.rad=1.0;
unit.mrad=1e-3*unit.rad;
unit.urad=1e-6*unit.rad;
unit.deg=(2*pi*unit.rad)/360;
unit.W = 1.0;
unit.mW = 1e-3*unit.W;
%%%%%%%%;
lightpipe = struct('type','lightpipe');
lightpipe.wavelength = 405 * unit.nm ; % wavelength ;
lightpipe.lam = lightpipe.wavelength ; % lambda = wavelength. ;
lightpipe.size       = 20  * unit.mm ; % simulation window size ;
lightpipe.N          = 5096 ; % grid resolution ;
lightpipe.Kmax       = k_p_r_max/max(1e-12,lightpipe.size) ; % bandlimit (max spatial-frequency) for x-y-plane in 1/unit.m. i.e., exp(\pm i*2*pi*k_p_r_max*x_p_r_) are the highest-frequency components;
lightpipe.z          = 1  * unit.m ; % long distance for Fraunhofer approx ;
lightpipe.a          = 2.5  * unit.mm ; % slit width ;
%%%%%%%%;
% circular aperture. ;
% Note (6.561.5 in TISPISMIGR): ;
% \int_{x=0}^{x=1} x^{q+1} besselj(q,a*x) dx = besselj(q+1,a)/a ;
% \int_{x=0}^{x=1}dx (\pm i).^(q0) .* x^{q0+1} .* besselj(q0,2*pi*k*x) = besselj(q0+1,(2*pi*k))/(2*pi*k) ;
% \int_{x=0}^{x=x_cut}dx x^{1} .* besselj(0,2*pi*k*x) = ...
% \int_{y=0}^{y=1}dy x_cut*x_cut*y^{1} .* besselj(0,2*pi*k*x_cut*y) = x_cut^2 * besselj(1,(2*pi*k*x_cut))/(2*pi*k*x_cut) ;
%%%%%%%%;

x_cut = ( (lightpipe.a/2) / (lightpipe.size/2) ) * x_p_r_max ;
M_x_q_x_form_ = 1.0*(x_p_r_ <= x_cut) ;
M_k_q_k_form_ = 2*pi*x_cut * besselj(1,(2*pi*k_p_r_*x_cut)) ./ max(1e-12,2*pi*k_p_r_) ;
x_step = x_cut/8;
M_x_q_x_orig_ = 1 - 0.5*(1 + erf((x_p_r_ - x_cut)/x_step));
lightpipe.focal_length = 10^8*lightpipe.wavelength / 4 ;
%M_x_q_x_orig_ = M_x_q_x_orig_.*exp(-i*2*pi/lightpipe.wavelength .* (x_p_r_*lightpipe.size/2).^2 / (2*lightpipe.focal_length));
tmp_t = tic();
flag_sign = -1; q_val = 0;
[ ...
 ~ ...
,M_k_q_k_fht2_ ...
] = ...
fht_2( ...
 [] ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_orig_ ...
,q_val ...
,flag_sign ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% fht_2: time %0.6fs',tmp_t)); end;
tmp_t = tic();
[ ...
 ~ ...
,M_k_q_k_sht1_ ...
] = ...
sht_1( ...
 [] ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
,M_x_q_x_orig_ ...
,q_val ...
,flag_sign ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% sht_1: time %0.6fs',tmp_t)); end;
tmp_t = tic();
flag_sign = +1; q_val = 0;
[ ...
 ~ ...
,M_x_q_x_fht2_ ...
] = ...
fht_2( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_sht1_ ...
,q_val ...
,flag_sign ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% fht_2: time %0.6fs',tmp_t)); end;
tmp_t = tic();
flag_sign = +1; q_val = 0;
[ ...
 ~ ...
,M_x_q_x_sht1_ ...
] = ...
sht_1( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_k_sht1_ ...
,q_val ...
,flag_sign ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% sht_1: time %0.6fs',tmp_t)); end;
fnorm_disp(flag_verbose,'M_k_q_k_sht1_',M_k_q_k_sht1_,'M_k_q_k_fht2_',M_k_q_k_fht2_,'%<-- should be small.');
fnorm_disp(flag_verbose,'M_x_q_x_orig_',M_x_q_x_orig_,'M_x_q_x_sht1_',M_x_q_x_sht1_,'%<-- should be small.');
fnorm_disp(flag_verbose,'M_x_q_x_orig_',M_x_q_x_orig_,'M_x_q_x_fht2_',M_x_q_x_fht2_,'%<-- should be small.');
%%%%;
if flag_disp>0;
M_x_lim_ = prctile(real(M_x_q_x_form_),[ 1,99]); M_x_lim_ = mean(M_x_lim_) + 1.25*0.5*diff(M_x_lim_)*[-1,+1];
M_k_lim_ = prctile(real(M_k_q_k_form_),[ 1,99]); M_k_lim_ = mean(M_k_lim_) + 1.25*0.5*diff(M_k_lim_)*[-1,+1];
subplot(1,2,1);
plot(x_p_r_,real(M_x_q_x_orig_),'kx',x_p_r_,real(M_x_q_x_sht1_),'go',x_p_r_,real(M_x_q_x_fht2_),'ro');
xlim([0,x_p_r_max]); ylim(M_x_lim_);
subplot(1,2,2);
plot(k_p_r_,real(M_k_q_k_form_),'kx',k_p_r_,real(M_k_q_k_sht1_),'go',k_p_r_,real(M_k_q_k_fht2_),'ro');
xlim([0,k_p_r_max]); ylim(M_k_lim_);
end;%if flag_disp>0;
%%%%%%%%;

flag_check=1;
if flag_check;
%%%%%%%%;
% evolve radial waveform. ;
%%%%%%%%;
lambda = lightpipe.wavelength;
z_max = 10^8*lambda ; % here we think of z in terms of number of wavelengths. ;
n_z = 1+0.5*1024; z_ = transpose(linspace(0,z_max,n_z)); dz = mean(diff(z_));
P_k_p_r_ = exp(-i*pi*lambda*dz*(k_p_r_/lightpipe.size).^2); %<-- radial fourier propagator ;
M_x_q_xz__ = zeros(n_x_p_r,n_z);
M_x_q_xz__(:,1+0) = M_x_q_x_orig_;
M_k_q_kz__ = zeros(n_k_p_r,n_z);
M_k_q_kz__(:,1+0) = M_k_q_k_fht2_;
for nz=1:n_z-1;
if (flag_verbose>0); if mod(nz,32)==0; disp(sprintf(' %% nz %.4d/%.4d',nz,n_z)); end; end;
M_k_q_kz__(:,1+nz+0) = M_k_q_kz__(:,1+nz-1).*P_k_p_r_;
tmp_t = tic();
flag_sign = +1; q_val = 0;
[ ...
 ~ ...
,M_x_q_xz__(:,1+nz+0) ...
] = ...
fht_2( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_2d_k_p_r_ ...
,M_k_q_kz__(:,1+nz+0) ...
,q_val ...
,flag_sign ...
,n_x_p_r ...
,x_p_r_ ...
,x_p_r_max ...
,weight_2d_x_p_r_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% fht_2: time %0.6fs',tmp_t)); end;
end;%for nz=1:n_z-1;
%%%%%%%%;
if flag_disp>0;
figure(1+nf);nf=nf+1;clf;figmed;
x_p_r_2x_ = [-flipud(x_p_r_) ; x_p_r_];
n_y = 2*n_x_p_r; y_p_r_2y_ = transpose(linspace(min(x_p_r_2x_),max(x_p_r_2x_),n_y));
subplot_{1} = subplot(2,1,1);
tmp_xz__ = M_x_q_xz__;
tmp_2xz__ = [+flipud(tmp_xz__) ; +tmp_xz__];
tmp_2yz__ = log(abs(interp1(x_p_r_2x_,tmp_2xz__,y_p_r_2y_)));
imagesc(tmp_2yz__,[-8,+2]);
set(gca,'ydir','normal'); axisnotick; title('abs(M_x_q_xz__)','Interpreter','none');
subplot_{2} = subplot(2,1,2);
tmp_xz__ = M_x_q_xz__;
tmp_2xz__ = [+flipud(tmp_xz__) ; +tmp_xz__];
tmp_2yz__ = angle(interp1(x_p_r_2x_,tmp_2xz__,y_p_r_2y_));
imagesc(tmp_2yz__,[-pi,+pi]);
set(gca,'ydir','normal'); axisnotick; title('angle(M_x_q_xz__)','Interpreter','none');
colormap(subplot_{2},colormap('hsv'));
colormap(subplot_{1},colormap_81s());
end;%if flag_disp>0;
%%%%%%%%;
if flag_disp>2;
figure(1+nf);nf=nf+1;clf;figmed;
x_p_r_2x_ = [-flipud(x_p_r_) ; x_p_r_];
subplot_{1} = subplot(2,1,1);
tmp_xz__ = log(abs(M_x_q_xz__));
tmp_2xz__ = [+flipud(tmp_xz__) ; +tmp_xz__];
imagesc_c(n_z,z_,2*n_x_p_r,x_p_r_2x_,transpose(tmp_2xz__),[-8,+2],colormap_81s);
xlim([0,z_max]); ylim(x_p_r_max*[-1,+1]); axisnotick;
set(gca,'ydir','normal'); axisnotick; title('abs(M_x_q_xz__)','Interpreter','none');
subplot_{2} = subplot(2,1,2);
tmp_xz__ = angle(M_x_q_xz__);
tmp_2xz__ = [+flipud(tmp_xz__) ; +tmp_xz__];
imagesc_c(n_z,z_,2*n_x_p_r,x_p_r_2x_,transpose(tmp_2xz__),[-pi,+pi],colormap('hsv'));
xlim([0,z_max]); ylim(x_p_r_max*[-1,+1]); axisnotick;
set(gca,'ydir','normal'); axisnotick; title('angle(M_x_q_xz__)','Interpreter','none');
%colormap(subplot_{2},colormap('hsv'));
%colormap(subplot_{1},colormap_81s());
end;%if flag_disp>2;
%%%%;
disp('returning'); return;
end;%if flag_check;
%%%%%%%%;



