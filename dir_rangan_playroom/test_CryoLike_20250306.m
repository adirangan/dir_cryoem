%%%%%%%%;
% Updated to include full likelihood calculation in likelihood_fourier <-- likelihood.py ;
% References include: I_PxxP_g.tex, ;
% Gradshteyn and Ryzhik: 5.5 ;
% Gradshteyn and Ryzhik: 6.5-6.7 ;
%%%%%%%%;
flag_verbose=1;
flag_disp=0; nf=0;
rng(0);

%%%%%%%%;
% test inner-product. ;
%%%%%%%%;
n_pixel = 7;
y_ = randn(n_pixel,1) + i*randn(n_pixel,1);
x_ = randn(n_pixel,1) + i*randn(n_pixel,1);
sigma = 0.3;
mu_r = -1.2; mu_i = 3.4;
mu = mu_r + i*mu_i;
nu = 0.7;
e_ = randn(n_pixel,1) + i*randn(n_pixel,1);
z_ = nu*x_ + mu*e_;
I_all_0 = sum(abs(y_ - z_).^2,'all');
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
I_all_1 = I_yy + nu^2*I_xx - 2*nu*real(I_xy) + mu_r^2*I_11 + mu_i^2*I_11 - 2*mu_r*(real(I_1y) - nu*real(I_1x)) - 2*mu_i*(imag(I_1y) - nu*imag(I_1x));
fnorm_disp(flag_verbose,'I_all_0',I_all_0,'I_all_1',I_all_1);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(b^2/max(1e-12,4*a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + 2*b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(b^2/max(1e-12,a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ; 
%%%%%%%%;
a = 3.0;
b = -0.5;
sigma = 1.75;
g = @(x_) exp(-(a*x_.^2 + 2*b*x_)/max(1e-12,2*sigma^2));
I_quad = integral(g,-10,+10);
I_form = sqrt(2*pi*sigma^2/max(1e-12,a)) * exp(b^2/max(1e-12,a*2*sigma^2));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
a = 3.0;
b = 2.25*i;
h = @(x_) exp(-a*x_.^2 + i*b*x_);
I_quad = integral(h,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(-b^2/max(1e-12,4*a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
n = -5.5; a = 3.25;
egamma = @(sigma) sigma.^n .* exp(-a./sigma.^2);
x_max = 64;
I_quad = integral(egamma,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:1:x_max]);
I_form = 0.5 * a.^(+n/2+1/2) * gamma(-n/2 - 1/2);
%x_ = linspace(0,x_max,1024); plot(x_,egamma(x_),'r-','LineWidth',3);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
n_pixel = 7; %<-- should be more than 3+1. ;
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_2B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = -real(I_1y)^2 - imag(I_1y)^2 + I_yy*I_11;
tmp_factor = (I_2B^2/I_A - I_C)/(2*I_11);
tmp_integrand = @(sigma_) (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * (sigma_).^(3-n_pixel) .* exp( tmp_factor ./ sigma_.^2 );
x_max = 64;
I_quad = integral(tmp_integrand,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:0.25:x_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (-tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(-tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral. ;
%%%%%%%%;
n_pixel = 9; %<-- should be more than 3+1. ;
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_2B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = -real(I_1y)^2 - imag(I_1y)^2 + I_yy*I_11;
tmp_factor = -(I_2B^2/I_A - I_C)/(2*I_11);
tmp_integrand = @(sigma_) (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * (sigma_).^(3-n_pixel) .* exp( -tmp_factor ./ sigma_.^2 );
x_max = 64;
I_quad = integral(tmp_integrand,1e-9,x_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:0.25:x_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test Gradshteyn and Ryzhik: 5.55. ;
%%%%%%%%;
p = 1+0.2; q = 1-0.2;
x_min = 3; x_max = 5;
tmp_f_ = @(x_) besselj(p,x_).*besselj(q,x_)./x_ ;
tmp_F_ = @(x_) x_.*( besselj(p,x_).*besselj(q+1,x_) - besselj(p+1,x_).*besselj(q,x_) )/(p^2-q^2) + besselj(p,x_).*besselj(q,x_)/(p+q) ;
I_quad = integral(tmp_f_,x_min,x_max,'AbsTol',1e-9,'RelTol',1e-9);
I_form = tmp_F_(x_max)-tmp_F_(x_min);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad);
%%%%%%%%;

%%%%%%%%;
% test integral using a very small image (denoted by y) and ctf*template (denoted by x). ;
%%%%%%%%;
n_pixel = 7;
y_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1); %<-- image. ;
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
y3 = y_(1+3);
y4 = y_(1+4);
y5 = y_(1+5);
y6 = y_(1+6);
x_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1); %<-- ctf*template. ;
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
x3 = x_(1+3);
x4 = x_(1+4);
x5 = x_(1+5);
x6 = x_(1+6);
sigma = 0.3; %<-- temperature. ;
mu_r = -0.2; mu_i = 0.45;
mu = mu_r + i*mu_i; %<-- offest. ;
nu = 0.37; %<-- intensity. ;
%e_ = ones(n_pixel,1); %<-- mask. ;
e_ = randn(n_pixel,1) + i*randn(n_pixel,1); %<-- mask. ;
e0 = e_(1+0);
e1 = e_(1+1);
e2 = e_(1+2);
e3 = e_(1+3);
e4 = e_(1+4);
e5 = e_(1+5);
e6 = e_(1+6);
ssnll_mu_nu_function_0 = @(mu_r_,mu_i_,nu_) ...
  + 0.5*conj(y0 - nu_*x0 - (mu_r_ + i*mu_i_)*e0).*(y0 - nu_*x0 - (mu_r_ + i*mu_i_)*e0) ...
  + 0.5*conj(y1 - nu_*x1 - (mu_r_ + i*mu_i_)*e1).*(y1 - nu_*x1 - (mu_r_ + i*mu_i_)*e1) ...
  + 0.5*conj(y2 - nu_*x2 - (mu_r_ + i*mu_i_)*e2).*(y2 - nu_*x2 - (mu_r_ + i*mu_i_)*e2) ...
  + 0.5*conj(y3 - nu_*x3 - (mu_r_ + i*mu_i_)*e3).*(y3 - nu_*x3 - (mu_r_ + i*mu_i_)*e3) ...
  + 0.5*conj(y4 - nu_*x4 - (mu_r_ + i*mu_i_)*e4).*(y4 - nu_*x4 - (mu_r_ + i*mu_i_)*e4) ...
  + 0.5*conj(y5 - nu_*x5 - (mu_r_ + i*mu_i_)*e5).*(y5 - nu_*x5 - (mu_r_ + i*mu_i_)*e5) ...
  + 0.5*conj(y6 - nu_*x6 - (mu_r_ + i*mu_i_)*e6).*(y6 - nu_*x6 - (mu_r_ + i*mu_i_)*e6) ...
  ;
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_2B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = -real(I_1y)^2 - imag(I_1y)^2 + I_yy*I_11;
ssnll_nu_function_0 = @(nu_) ...
  0.5*(I_A*nu_.^2 -2*I_2B*nu_ + I_C)./I_11 ...
  ;
ybar_ = I_1y*e_/max(1e-12,I_11); %<-- image average. ;
xbar_ = I_1x*e_/max(1e-12,I_11); %<-- ctf*template average. ;
ycen_ = y_ - ybar_; %<-- value-centered image. ;
xcen_ = x_ - xbar_; %<-- value-centered ctf*template. ;
I_ycenycen = sum(conj(ycen_).*ycen_,'all');
I_xcenycen = sum(conj(xcen_).*ycen_,'all');
I_ycenxcen = sum(conj(ycen_).*xcen_,'all');
I_xcenxcen = sum(conj(xcen_).*xcen_,'all');
ssnll_nu_function_1 = @(nu_) ...
  0.5*(I_ycenycen + nu_.^2*I_xcenxcen - nu_*I_xcenycen - nu_*I_ycenxcen) ...
  ;
ynrm_ = ycen_/max(1e-12,sqrt(I_ycenycen)); %<-- normalized image. ;
xnrm_ = xcen_/max(1e-12,sqrt(I_xcenxcen)); %<-- normalized ctf*template. ;
I_ynrmxnrm = sum(conj(ynrm_).*xnrm_,'all');
I_xnrmynrm = sum(conj(xnrm_).*ynrm_,'all');
I_costheta = 0.5*(I_xnrmynrm + I_ynrmxnrm); %<-- similarity measure. ;
I_sinthetasquared = 1 - I_costheta.^2; %<-- dis-similarity measure. ;
ssnll_2 = 0.5*I_ycenycen*I_sinthetasquared ;
%%%%%%%%;
% set up quadrature weights for mu_ := [mu*cos(psi);mu*sin(psi)]. ;
%%%%%%%%;
mu_int = 8;
mu_eq_d_double = 1/16;
n_psi_int = 8;
mu_p_r_max = mu_int/(2*pi); mu_eq_d = mu_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_mu_p_r ...
,mu_p_r_ ...
,weight_3d_mu_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,mu_p_r_max ...
,mu_eq_d ...
,str_L ...
);
%%%%;
n_psi_max = n_psi_int*2*(mu_int+1); n_psi_0in = n_psi_max; n_psi_0in_ = n_psi_max*ones(n_mu_p_r,1);
[ ...
 n_psi_ ...
,weight_2d_mu_p_r_ ...
,weight_2d_psimu_ ...
,mu_p_r_psimu_ ...
,mu_p_psi_psimu_ ...
,mu_c_0_psimu_ ...
,mu_c_1_psimu_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_mu_p_r ...
,mu_p_r_ ...
,mu_p_r_max ...
,-1 ...
,n_psi_0in_ ...
);
n_psi_max = max(n_psi_);
n_psi_sum = sum(n_psi_);
n_psi_csum_ = cumsum([0;n_psi_]);
mu_psimu_ = mu_c_0_psimu_ + i*mu_c_1_psimu_ ;
%%%%%%%%;
likelihood_mu_nu_sigma_function_0 = @(mu_r_,mu_i_,nu_,sigma) exp(-ssnll_mu_nu_function_0(mu_r_,mu_i_,nu_)./sigma.^2) .* (2*pi*sigma.^2).^(-n_pixel/2) ;
likelihood_mu_nu_sigma_psimu_ = likelihood_mu_nu_sigma_function_0(mu_c_0_psimu_,mu_c_1_psimu_,nu,sigma);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
likelihood_lim_ = [0,max(abs(likelihood_mu_nu_sigma_psimu_))];
imagesc_p(n_mu_p_r,mu_p_r_,n_psi_,n_psi_sum,likelihood_mu_nu_sigma_psimu_,likelihood_lim_,colormap_81s);
axis image; axisnotick;
title('likelihood_mu_nu_sigma_psimu_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
likelihood_nu_sigma_function_0 = @(nu_,sigma_) exp(-ssnll_nu_function_0(nu_)./sigma_.^2) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11) ;
likelihood_nu_sigma_function_1 = @(nu_,sigma_) exp(-ssnll_nu_function_1(nu_)./sigma_.^2) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11) ;
likelihood_nu_sigma_auto_0 = integral2( ...
 @(mu_r,mu_i) likelihood_mu_nu_sigma_function_0(mu_r,mu_i,nu,sigma) ...
,-20,+20,-20,+20 ...
,'AbsTol',1e-9 ...
,'RelTol',1e-9 ...
); %<-- built-in adaptive quadrature is not very accurate. ;
likelihood_nu_sigma_quad_0 = sum(likelihood_mu_nu_sigma_psimu_.*weight_2d_psimu_,'all')*(2*pi)^2; %<-- homebrew quadrature is much more accurate. ;
likelihood_nu_sigma_form_0 = likelihood_nu_sigma_function_0(nu,sigma);
likelihood_nu_sigma_form_1 = likelihood_nu_sigma_function_1(nu,sigma);
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_auto_0',likelihood_nu_sigma_auto_0,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_quad_0',likelihood_nu_sigma_quad_0,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_form_1',likelihood_nu_sigma_form_1,' %<-- should be zero');
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
likelihood_sigma_auto_1 = integral(@(nu_) likelihood_nu_sigma_function_1(nu_,sigma),-10,+10);
likelihood_sigma_quad_1 = sum(bsxfun(@times,likelihood_nu_sigma_function_1(nu_value_q_,sigma),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form_1 = sqrt(2*pi*sigma^2*I_11/I_A)*exp((I_2B^2/(2*sigma^2*I_11*I_A)))*exp(-I_C/(2*sigma^2*I_11)) .* (2*pi*sigma.^2).^(-n_pixel/2) .* (2*pi*sigma.^2./I_11);
likelihood_sigma_form_2 = sqrt(2*pi*sigma^2/I_xcenxcen)*exp(-ssnll_2/(sigma^2)) .* (2*pi*sigma.^2).^(-n_pixel/2) .* (2*pi*sigma.^2./I_11);
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_auto_1',likelihood_sigma_auto_1,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_sigma_form_2',likelihood_sigma_form_2,'likelihood_sigma_auto_1',likelihood_sigma_auto_1,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_quad_1',likelihood_sigma_quad_1,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_sigma_form_2',likelihood_sigma_form_2,'likelihood_sigma_quad_1',likelihood_sigma_quad_1,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_form_2',likelihood_sigma_form_2,' %<-- should be zero');
%%%%%%%%;
sigma_lim_ = 16*[0,+1];
n_sigma_q = 256;
[sigma_value_q_,sigma_weight_q_] = legpts(n_sigma_q,sigma_lim_);
likelihood_sigma_function_2 = @(sigma_) sqrt(2*pi*sigma_.^2./I_xcenxcen).*exp(-ssnll_2./(sigma_.^2)) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11);
likelihood_auto_2 = integral(@(sigma_) likelihood_sigma_function_2(sigma_),0,10);
likelihood_quad_2 = sum(bsxfun(@times,likelihood_sigma_function_2(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form_2 = (2*pi)^(3/2-n_pixel/2) * (I_xcenxcen).^(-1/2) * 0.5 * (ssnll_2).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2) / I_11;
fnorm_disp(flag_verbose,'likelihood_form_2',likelihood_form_2,'likelihood_auto_2',likelihood_auto_2,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_form_2',likelihood_form_2,'likelihood_quad_2',likelihood_quad_2,' %<-- should be small');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now set up a similar test, but with a fixed image and template. ;
% Here the image will be denoted by "M", the template by "S" and the CTF by "CTF". ;
% The CTF*S will be denoted by "T". ;
% The mask will be denoted by "E". ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First set up polar-grid quadrature. ;
%%%%%%%%;
k_int = 8;
k_eq_d_double = 0.250;
n_w_int = 4;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in = n_w_max; n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% Set up template S, CTF, image M and mask E. ;
%%%%%%%%;
rng(2);
K = k_p_r_max;
phi_S = 2*pi*rand();
delta_S_ = 1e-1*2*pi*rand(2,1); delta_S = fnorm(delta_S_); omega_S = atan2(delta_S_(1+1),delta_S_(1+0));
phi_M = 2*pi*rand();
delta_M_ = 1e-1*2*pi*rand(2,1); delta_M = fnorm(delta_M_); omega_M = atan2(delta_M_(1+1),delta_M_(1+0));
S_k_p_wk_ = exp(2*pi*i*k_p_r_wk_.*delta_S.*cos(k_p_w_wk_ - omega_S));
CTF_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - phi_S);
T_k_p_wk_ = CTF_k_p_wk_ .* S_k_p_wk_ ;
M_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - phi_M).*exp(2*pi*i*k_p_r_wk_.*delta_M.*cos(k_p_w_wk_ - omega_M));
delta_E_ = 1e-1*2*pi*rand(2,1); delta_E = fnorm(delta_E_); omega_E = atan2(delta_E_(1+1),delta_E_(1+0));
E_k_p_wk_ = exp(2*pi*i*k_p_r_wk_.*delta_E.*cos(k_p_w_wk_ - omega_E));
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig; p_row=2;p_col=3;np=0;
Tlim_ = max(abs(T_k_p_wk_),[],'all')*[-1,+1];
Mlim_ = max(abs(M_k_p_wk_),[],'all')*[-1,+1];
Elim_ = max(abs(E_k_p_wk_),[],'all')*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_),Tlim_); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_wk_),Tlim_); axis image; axisnotick; title('imag(T_k_p_wk_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_wk_),Mlim_); axis image; axisnotick; title('real(M_k_p_wk_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_wk_),Mlim_); axis image; axisnotick; title('imag(M_k_p_wk_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(E_k_p_wk_),Elim_); axis image; axisnotick; title('real(E_k_p_wk_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1; imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(E_k_p_wk_),Elim_); axis image; axisnotick; title('imag(E_k_p_wk_)','Interpreter','none');
end;%if flag_disp;
%%%%;
I_T_T_quad = sum(conj(T_k_p_wk_).*(T_k_p_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_T_T_form ...
] = ...
I_xPPx_0( ...
 K ...
,phi_S ...
,delta_S_ ...
,phi_S ...
,delta_S_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_T_T_quad',I_T_T_quad,'I_T_T_form',I_T_T_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
I_M_M_quad = sum(conj(M_k_p_wk_).*(M_k_p_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_M_M_form ...
] = ...
I_xPPx_0( ...
 K ...
,phi_M ...
,delta_M_ ...
,phi_M ...
,delta_M_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_M_M_quad',I_M_M_quad,'I_M_M_form',I_M_M_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
I_E_E_quad = sum(conj(E_k_p_wk_).*(E_k_p_wk_).*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_E_E_form ...
] = ...
I_PP_0( ...
 K ...
,delta_E_ ...
,delta_E_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_E_E_quad',I_E_E_quad,'I_E_E_form',I_E_E_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
I_T_M_quad = sum(conj(T_k_p_wk_).*M_k_p_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_T_M_form ...
] = ...
I_xPPx_0( ...
 K ...
,phi_S ...
,delta_S_ ...
,phi_M ...
,delta_M_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_T_M_quad',I_T_M_quad,'I_T_M_form',I_T_M_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
I_T_E_quad = sum(conj(T_k_p_wk_).*E_k_p_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_T_E_form ...
] = ...
I_xPP_0( ...
 K ...
,phi_S ...
,delta_S_ ...
,delta_E_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_T_E_quad',I_T_E_quad,'I_T_E_form',I_T_E_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
I_M_E_quad = sum(conj(M_k_p_wk_).*E_k_p_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
[ ...
 I_M_E_form ...
] = ...
I_xPP_0( ...
 K ...
,phi_M ...
,delta_M_ ...
,delta_E_ ...
);
%%%%;
fnorm_disp(flag_verbose,'I_M_E_quad',I_M_E_quad,'I_M_E_form',I_M_E_form,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%;
%%%%%%%%;
% Quick test of standard l2-distance. ;
%%%%%%%%;
nu = 1.3; mu_r = 0.4; mu_i = -0.7; mu = mu_r + i*mu_i; lu = 1.0; lu_ = 1.0;
tmp_wk_ = nu*T_k_p_wk_ + mu*E_k_p_wk_ - M_k_p_wk_ ;
ssnll_quad_0 = 0.5*sum(conj(tmp_wk_).*tmp_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
ssnll_mu_nu_function_0 = ...
@(mu_,nu_) ...
 + 0.5*conj(nu_).*nu_.*I_T_T_form ...
 + 0.5*conj(mu_).*mu_.*I_E_E_form ...
 + 0.5*conj(lu_).*lu_.*I_M_M_form ...
 + 1.0*real(conj(nu_).*mu_.*I_T_E_form) ...
 - 1.0*real(conj(nu_).*lu_.*I_T_M_form) ...
 - 1.0*real(conj(lu_).*mu_.*I_M_E_form) ...
;
ssnll_form_0 = ssnll_mu_nu_function_0(mu,nu);
fnorm_disp(flag_verbose,'ssnll_quad_0',ssnll_quad_0,'ssnll_form_0',ssnll_form_0,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%%%%%;

%%%%%%%%;
% Now set up and test quadrature for integration over complex mu, real nu and positive sigma. ;
%%%%%%%%;

%%%%%%%%;
% First define the inner-products that involve formulas. ;
%%%%%%%%;
n_pixel = I_E_E_form; %<-- the number of relevant pixels is given by the l2-norm of the mask. ;
I_M_T_form = conj(I_T_M_form);
I_E_T_form = conj(I_T_E_form);
I_E_M_form = conj(I_M_E_form);
I_Tcen_Tcen_form = I_T_T_form - I_T_E_form*I_E_T_form/max(1e-12,I_E_E_form);
I_Mcen_Mcen_form = I_M_M_form - I_M_E_form*I_E_M_form/max(1e-12,I_E_E_form);
I_Tcen_Mcen_form = I_T_M_form - I_T_E_form*I_E_M_form/max(1e-12,I_E_E_form);
I_Mcen_Tcen_form = I_M_T_form - I_M_E_form*I_E_T_form/max(1e-12,I_E_E_form);
I_Tnrm_Tnrm_form = I_Tcen_Tcen_form/max(1e-12,sqrt(I_Tcen_Tcen_form*I_Tcen_Tcen_form));
I_Tnrm_Mnrm_form = I_Tcen_Mcen_form/max(1e-12,sqrt(I_Tcen_Tcen_form*I_Mcen_Mcen_form));
I_Mnrm_Tnrm_form = I_Mcen_Tcen_form/max(1e-12,sqrt(I_Mcen_Mcen_form*I_Tcen_Tcen_form));
I_Mnrm_Mnrm_form = I_Mcen_Mcen_form/max(1e-12,sqrt(I_Mcen_Mcen_form*I_Mcen_Mcen_form));
I_costheta = 0.5*(I_Mnrm_Tnrm_form + I_Tnrm_Mnrm_form); %<-- similarity measure. ;
I_sinthetasquared = 1 - I_costheta.^2; %<-- dis-similarity measure. ;
%%%%%%%%;
% Now define quadrature-grid for mu. ;
%%%%%%%%;
mu_int = 24;
mu_eq_d_double = 1/16;
n_psi_int = 8;
mu_p_r_max = mu_int/(2*pi); mu_eq_d = mu_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_mu_p_r ...
,mu_p_r_ ...
,weight_3d_mu_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*flag_verbose ...
,mu_p_r_max ...
,mu_eq_d ...
,str_L ...
);
%%%%;
n_psi_max = n_psi_int*2*(mu_int+1); n_psi_0in = n_psi_max; n_psi_0in_ = n_psi_max*ones(n_mu_p_r,1);
[ ...
 n_psi_ ...
,weight_2d_mu_p_r_ ...
,weight_2d_psimu_ ...
,mu_p_r_psimu_ ...
,mu_p_psi_psimu_ ...
,mu_c_0_psimu_ ...
,mu_c_1_psimu_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_mu_p_r ...
,mu_p_r_ ...
,mu_p_r_max ...
,-1 ...
,n_psi_0in_ ...
);
n_psi_max = max(n_psi_);
n_psi_sum = sum(n_psi_);
n_psi_csum_ = cumsum([0;n_psi_]);
mu_psimu_ = mu_c_0_psimu_ + i*mu_c_1_psimu_ ;
%%%%%%%%;
% Now define formula for integration over mu. ;
%%%%%%%%;
ssnll_nu_sigma_function = @(nu_,sigma_) ...
  0.5*(I_Mcen_Mcen_form + nu_.^2*I_Tcen_Tcen_form - nu_*I_Tcen_Mcen_form - nu_*I_Mcen_Tcen_form) ...
  ;
likelihood_nu_sigma_function = @(nu_,sigma_) exp(-ssnll_nu_sigma_function(nu_,sigma_)./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_E_E_form) ;
%%%%%%%%;
% Now test integration over mu. ;
%%%%%%%%;
sigma_value = 1.3; nu_value = 0.8;
likelihood_nu_sigma_psimu_ = (2*pi*sigma_value^2).^(-n_pixel/2) * exp(-ssnll_mu_nu_function_0(mu_psimu_,nu_value)/max(1e-12,sigma_value^2)) ;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
l_lim_ = [0,max(abs(likelihood_nu_sigma_psimu_))];
imagesc_p(n_mu_p_r,mu_p_r_,n_psi_,n_psi_sum,likelihood_nu_sigma_psimu_,l_lim_,colormap_81s);
axis image; axisnotick;
title('likelihood_nu_sigma_psimu_','Interpreter','none');
end;%if flag_disp;
likelihood_nu_sigma_quad = sum(bsxfun(@times,likelihood_nu_sigma_psimu_,reshape(weight_2d_psimu_,[n_psi_sum,1])),'all')*(2*pi)^2 ;
likelihood_nu_sigma_form = likelihood_nu_sigma_function(nu_value,sigma_value);
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form',likelihood_nu_sigma_form,'likelihood_nu_sigma_quad',likelihood_nu_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
% Now define quadrature-grid for nu. ;
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
%%%%%%%%;
% Now define formula for integration over nu. ;
%%%%%%%%;
ssnll_sigma_function = ...
  0.5*I_Mcen_Mcen_form*I_sinthetasquared ...
  ;
likelihood_sigma_function = @(sigma_) sqrt(2*pi*sigma_.^2./max(1e-12,I_Tcen_Tcen_form)) .* exp(-ssnll_sigma_function./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* (2*pi*sigma_.^2./max(1e-12,I_E_E_form));
%%%%%%%%;
% Now test integration over nu. ;
%%%%%%%%;
sigma_value = 1.5;
likelihood_sigma_quad = sum(bsxfun(@times,likelihood_nu_sigma_function(nu_value_q_,sigma_value),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form = likelihood_sigma_function(sigma_value);
fnorm_disp(flag_verbose,'likelihood_sigma_form',likelihood_sigma_form,'likelihood_sigma_quad',likelihood_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
% Now define quadrature-grid for sigma. ;
%%%%%%%%;
sigma_lim_ = 1024*[0,+1];
n_sigma_q = 8192;
[sigma_value_q_,sigma_weight_q_] = legpts(n_sigma_q,sigma_lim_);
%%%%%%%%;
% Now define formula for integration over sigma. ;
%%%%%%%%;
likelihood_function = (2*pi)^(3/2-n_pixel/2) * (max(1e-12,I_Tcen_Tcen_form)).^(-1/2) * 0.5 * (max(1e-12,ssnll_sigma_function)).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2) / max(1e-12,I_E_E_form) ;
%%%%%%%%;
% Now test integration over sigma. ;
%%%%%%%%;
likelihood_quad = sum(bsxfun(@times,likelihood_sigma_function(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form = likelihood_function;
fnorm_disp(flag_verbose,'likelihood_form',likelihood_form,'likelihood_quad',likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');
%%%%%%%%;
% Now test formula for log-likelihood. ;
%%%%%%%%;
negative_log_likelihood_function = ...
-(3/2 - n_pixel/2) * log(2*pi) ...
+(1/2) * log(max(1e-12,I_Tcen_Tcen_form)) ...
+(1.0) * log(2) ...
-(2 - n_pixel/2) * log(max(1e-12,ssnll_sigma_function)) ...
- gammaln(n_pixel/2 - 2) ...
+ log(max(1e-12,I_E_E_form)) ...
;
negative_log_likelihood_quad = -log(likelihood_quad);
negative_log_likelihood_form = negative_log_likelihood_function;
fnorm_disp(flag_verbose,'negative_log_likelihood_form',negative_log_likelihood_form,'negative_log_likelihood_quad',negative_log_likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');


