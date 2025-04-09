%%%%%%%%;
% This script is intended to help set up an analytical test of the likelihood-calculations in CryoLike. ;
% The general notation is: ;
% n_pixel: (integer) size of each image and/or template, although later this will simply be a (double) serving as a parameter in the likelihood. ;
% y_ or T_: (complex) array of size (n_pixel,1) storing the appropriately oriented CTF*template-values. ;
% x_ or M_: (complex) array of size (n_pixel,1) storing the image-values. ;
% e_ or E_: (complex) array of size (n_pixel,1) storing the mask-values. ;
% mu: (complex) offset used as prefactor for the mask. ;
% mu_r: (real) real part of mu. ;
% mu_i: (real) imaginary part of mu (i.e., mu = mu_r + i*mu_i). ;
% nu: (real) prefactor used for the image. ;
% z_: (complex) array of size (n_pixel,1) storing the shifted and scaled image: z_ := nu*x_ + mu*e_. ;
% sigma: (real) estimated temperature or noise-level used in likelihood calculations. Assumed to be positive. ;
%%%%%%%%;
% Note that, for this script we will assume that we are comparing the CTF*template to the shifted-and-scaled image z_ := nu*x_ + mu*e_. ;
% This is different from earlier versions of the script, and also different from the current appendix to the paper, ;
% which shift and scale the ctf*template, rather than the image. ;
% This impacts the final scaling
% (Note, we should probably fix the appendix afterwards to be consistent with this script). ;
%%%%%%%%;
% We will also refer to:
% likelihood_mu_nu_sigma: (function) the likelihood associated with the image and ctf*template as a function of mu, nu and sigma. ;
% ssnll_mu_nu: (function) the dominant term in the numerator of the log of likelihood_mu_nu_sigma (which itself does not depend on sigma). ;
% likelihood_nu_sigma: (function) the likelihood associated with the image and ctf*template as a function of nu and sigma (after marginalizing over mu). ;
% ssnll_nu: (function) the dominant term in the numerator of the log of likelihood_nu_sigma (which itself does not depend on sigma). ;
% likelihood_sigma: (function) the likelihood associated with the image and ctf*template as a function of sigma (after marginalizing over mu and nu). ;
% ssnll: (function, but also just a constant) the dominant term in the numerator of the log of likelihood_sigma (which itself does not depend on sigma). ;
% likelihood: (function, but also just a constant) the likelihood associated with the image and ctf*template (after marginalizing over mu, nu and sigma). ;
%%%%%%%%;
% In performing these likelihood calculations we will refer to: ;
% ybar_ or Tbar_: (complex) array storing the projection of the ctf*template onto the mask (e.g., in standard usage the 'value-averaged' ctf*template). ;
% ycen_ or Tcen_: (complex) array storing the portion of the ctf*template perpendicular to the mask (e.g., the 'value-centered' ctf*template). ;
% ynrm_ or Tnrm_: (complex) array storing a normalized ctf*template perpendicular to the mask (e.g., the 'normalized' ctf*template). ;
% xbar_ or Mbar_: (complex) array storing the projection of the image onto the mask (e.g., in standard usage the 'value-averaged' image). ;
% xcen_ or Mcen_: (complex) array storing the portion of the image perpendicular to the mask (e.g., the 'value-centered' image). ;
% xnrm_ or Mnrm_: (complex) array storing a normalized image perpendicular to the mask (e.g., the 'normalized' image). ;
%%%%%%%%;
% Finally, we will refer to: ;
% costheta: (real) This is the real-part of the pearson-correlation between the ctf*template and the image. This is a measure of similarity. ;
% sinthetasquared: (real) This is 1 - costheta^2, which is a measure of dis-similarity. ;
%%%%%%%%;
% Ultimately, we should set up a similar test in CryoLike. ;
% This will require specifying the following: ;
% The image should be a plane-wave times a planar-function. ;
% The template should be a plane-wave. ;
% The CTF should be a planar-function. ;
% The mask should be a plane-wave. ;
% Note that, right now, the mask is hard-coded (as a product of sincs). ;
% This will have to be altered to allow for this test. ;
% Note also that a full test of this kind will require sorting out the effects of the translation and in-plane rotation, ;
% in terms of how they affect the relative orientation of the CTF, template and image (see test_cross_correlation_2.py). ;
%%%%%%%%;

%%%%%%%%;
% References include: 
% I_PxxP_g.tex, I_xPPx_0.m, I_xPP_0.m, I_PP_0.m, test_cross_correlation_2.py ;
%%%%;
% Gradshteyn and Ryzhik: 2.33 ;
% Gradshteyn and Ryzhik: 3.326.2 ;
% Gradshteyn and Ryzhik: 5.5 ;
% Gradshteyn and Ryzhik: 6.5-6.7 ;
% Note: Gradshteyn and Ryzhik can be accessed via: ;
% http://fisica.ciens.ucv.ve/~svincenz/TISPISGIMR.pdf ;
%%%%%%%%;
flag_verbose=1;
flag_disp=0; nf=0;
rng(0);

%%%%%%%%;
% Here we just test to make sure that I have the correct formula for the redistribution of an l2-norm into a sum of inner-products. ;
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
I_all_0 = sum(abs(y_ - z_).^2,'all'); %<-- this is the standard l2-norm of (y_-z_). ;
I_yy = sum(conj(y_).*y_,'all'); %<-- this is the inner-product of y_ and y_. ;
I_xy = sum(conj(x_).*y_,'all'); %<-- this is the inner-product of x_ and y_, and so forth. ;
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
I_all_1 = I_yy + nu^2*I_xx - 2*nu*real(I_xy) + mu_r^2*I_11 + mu_i^2*I_11 - 2*mu_r*(real(I_1y) - nu*real(I_1x)) - 2*mu_i*(imag(I_1y) - nu*imag(I_1x)); %<-- this is also the l2-norm of (y_-z_), written out in terms of inner-products. ;
fnorm_disp(flag_verbose,'I_all_0',I_all_0,'I_all_1',I_all_1,' %<-- should be zero'); %<-- now compare the two. ;
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
% Here we test another version of the integral in Gradshteyn and Ryzhik: 2.33 (just to make sure there are no typos). ;
%%%%%%%%;
a = 3.0;
b = -0.5;
g = @(x_) exp(-a*x_.^2 + 2*b*x_);
I_quad = integral(g,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(b^2/max(1e-12,a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Here we test another version of the integral in Gradshteyn and Ryzhik: 2.33 (just to make sure there are no typos). ;
%%%%%%%%;
a = 3.0;
b = -0.5;
sigma = 1.75;
g = @(x_) exp(-(a*x_.^2 + 2*b*x_)/max(1e-12,2*sigma^2));
I_quad = integral(g,-10,+10);
I_form = sqrt(2*pi*sigma^2/max(1e-12,a)) * exp(b^2/max(1e-12,a*2*sigma^2));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Here we test another version of the integral in Gradshteyn and Ryzhik: 2.33 (just to make sure there are no typos). ;
%%%%%%%%;
a = 3.0;
b = 2.25*i;
h = @(x_) exp(-a*x_.^2 + i*b*x_);
I_quad = integral(h,-10,+10);
I_form = sqrt(pi/max(1e-12,a)) * exp(-b^2/max(1e-12,4*a));
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved.');
%%%%%%%%;

%%%%%%%%;
% Here we test a version of the integral in Gradshteyn and Ryzhik: 3.326.2 (just to make sure there are no typos). ;
%%%%%%%%;
n = -5.5; a = 3.25;
egamma = @(sigma) sigma.^n .* exp(-a./sigma.^2);
sigma_max = 64;
I_quad = integral(egamma,1e-9,sigma_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:1:sigma_max]);
I_form = 0.5 * a.^(+n/2+1/2) * gamma(-n/2 - 1/2);
%sigma_ = linspace(0,sigma_max,1024); subplot(1,2,1);plot(sigma_,egamma(sigma_),'r-','LineWidth',3); subplot(1,2,2);semilogy(sigma_,egamma(sigma_),'r-','LineWidth',3); %<-- you can plot the integrand to make sure the numerical integral is resolved. ;
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- can be difficult to resolve due to algebraic decay of integrand.');
%%%%%%%%;
if flag_disp;
tmp_n = -2;
sigma_max = 10; n_sigma = 256; sigma_ = transpose(linspace(0,sigma_max,n_sigma)); dsigma = mean(diff(sigma_));
x_max = sigma_max; n_x = 1*n_sigma; x_0_ = transpose(linspace(-x_max,+x_max,n_x)); x_1_ = transpose(linspace(-x_max,+x_max,n_x));
[x_0_xx__,x_1_xx__] = ndgrid(x_0_,x_1_);
[x_0_xxs___,x_1_xxs___,s_xxs___] = ndgrid(x_0_,x_1_,sigma_);
g_xxs = @(x_0_,x_1_,sigma_) 1/(2*pi) ./max(1e-12,sigma_.^2) .* exp(-(x_0_.^2 + x_1_.^2)./max(1e-12,(2*sigma_.^2))) ;
g_xxs___ = g_xxs(x_0_xxs___,x_1_xxs___,s_xxs___);
g_xx = @(x_0_,x_1_) 1/(2*pi) .* 0.5 * (0.5*(x_0_.^2 + x_1_.^2)).^(+tmp_n/2+1/2) .* gamma(-tmp_n/2-1/2) ;
g_xx__ = g_xx(x_0_xx__,x_1_xx__);
figure(1+nf);nf=nf+1;clf;figbig;
isosurface_f_x_u_1(struct('vval_',5e-3),g_xxs___(:));
zlabel('sigma');
fname_fig_pre = '../dir_CryoBIFE_MD/test_CryoLike_GR.3.326.2';
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
print('-djpeg',fname_fig_jpg);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Here we test another version of the integral in Gradshteyn and Ryzhik: 3.326.2 (just to make sure there are no typos). ;
% This is closer to the structure we will see later when calculating the likelihood after marginalizing over mu and nu. ;
%%%%%%%%;
n_pixel = 7; %<-- should be more than 3+1. ;
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_2B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = -real(I_1y)^2 - imag(I_1y)^2 + I_yy*I_11;
tmp_factor = (I_2B^2 - I_A*I_C)/max(1e-12,2*I_11*I_A);
tmp_integrand = @(sigma_) (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * (sigma_).^(3-n_pixel) .* exp( tmp_factor ./ sigma_.^2 );
sigma_max = 64;
I_quad = integral(tmp_integrand,1e-9,sigma_max,'AbsTol',1e-9,'RelTol',1e-9,'WayPoints',[1e-9:0.25:sigma_max]);
I_form = (2*pi)^(3/2-n_pixel/2) * (I_11*I_A).^(-1/2) * 0.5 * (-tmp_factor).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- can be difficult to resolve due to algebraic decay of integrand.');
lI_form = (3/2-n_pixel/2)*log(2*pi) -0.5*log(I_11*I_A) - log(2) - (n_pixel/2-2)*log(-tmp_factor) + gammaln(n_pixel/2 - 2);
I_form = exp(lI_form);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- can be difficult to resolve due to algebraic decay of integrand.');
%%%%%%%%;

%%%%%%%%;
% Here we test a version of the integral in Gradshteyn and Ryzhik: 5.55 (just to make sure there are no typos). ;
%%%%%%%%;
p = 1+0.2; q = 1-0.2;
x_min = 3; x_max = 5;
tmp_f_ = @(x_) besselj(p,x_).*besselj(q,x_)./x_ ;
tmp_F_ = @(x_) x_.*( besselj(p,x_).*besselj(q+1,x_) - besselj(p+1,x_).*besselj(q,x_) )/(p^2-q^2) + besselj(p,x_).*besselj(q,x_)/(p+q) ;
I_quad = integral(tmp_f_,x_min,x_max,'AbsTol',1e-9,'RelTol',1e-9);
I_form = tmp_F_(x_max)-tmp_F_(x_min);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %<-- should be zero if integral is resolved');
%%%%%%%%;

%%%%%%%%;
% Now we test out the various likelihood-formulae using a very small image (denoted by x_) and ctf*template (denoted by y_). ;
% Recall that: ;
% likelihood_mu_nu_sigma refers to the likelihood associated with the image and ctf*template as a function of mu, nu and sigma. ;
% likelihood_nu_sigma refers to the likelihood associated with the image and ctf*template as a function of nu and sigma (after marginalizing over mu). ;
% likelihood_sigma refers to the likelihood associated with the image and ctf*template as a function of sigma (after marginalizing over mu and nu). ;
% likelihood refers to the likelihood associated with the image and ctf*template (after marginalizing over mu, nu and sigma). ;
% The 'ssnll' terms are used as sub-steps towards calculating the likelihood.
% Also recall that: ;
% ybar_ stores the projection of the ctf*template onto the mask. ;
% ycen_ stores the portion of the ctf*template perpendicular to the mask. ;
% ynrm_ stores a normalized ctf*template perpendicular to the mask. ;
% xbar_ stores the projection of the image onto the mask. ;
% xcen_ stores the portion of the image perpendicular to the mask. ;
% xnrm_ stores a normalized image perpendicular to the mask. ;
%%%%%%%%;
n_pixel = 7;
%%%%%%%%;
% First we define the ctf*template, and extract the individual pixel-values (for the inline-function ssnll_mu_nu_function). ;
%%%%%%%%;
y_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1); %<-- ctf*template. ;
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
y3 = y_(1+3);
y4 = y_(1+4);
y5 = y_(1+5);
y6 = y_(1+6);
%%%%%%%%;
% Now we define the image. ;
%%%%%%%%;
x_ = 0.25*randn(n_pixel,1) + i*0.25*randn(n_pixel,1); %<-- image. ;
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
x3 = x_(1+3);
x4 = x_(1+4);
x5 = x_(1+5);
x6 = x_(1+6);
%%%%%%%%;
% Now we fix some of the parameters for integration. ;
%%%%%%%%;
sigma = 0.3; %<-- temperature. ;
mu_r = -0.2; mu_i = 0.45;
mu = mu_r + i*mu_i; %<-- offset. ;
nu = 0.37; %<-- intensity. ;
%%%%%%%%;
% Now we define the mask. ;
%%%%%%%%;
%e_ = ones(n_pixel,1); %<-- typical mask, but we will not use this. ;
e_ = randn(n_pixel,1) + i*randn(n_pixel,1); %<-- more complicated mask used for testing. ;
e0 = e_(1+0);
e1 = e_(1+1);
e2 = e_(1+2);
e3 = e_(1+3);
e4 = e_(1+4);
e5 = e_(1+5);
e6 = e_(1+6);
%%%%%%%%;
% The reason we extracted the individual pixel values was so that we could define the following inline-function easily. ;
%%%%%%%%;
ssnll_mu_nu_function_0 = @(mu_r_,mu_i_,nu_) ... %<-- a standard definition for the ssnll_mu_nu_function. ;
  + 0.5*conj(y0 - nu_*x0 - (mu_r_ + i*mu_i_)*e0).*(y0 - nu_*x0 - (mu_r_ + i*mu_i_)*e0) ...
  + 0.5*conj(y1 - nu_*x1 - (mu_r_ + i*mu_i_)*e1).*(y1 - nu_*x1 - (mu_r_ + i*mu_i_)*e1) ...
  + 0.5*conj(y2 - nu_*x2 - (mu_r_ + i*mu_i_)*e2).*(y2 - nu_*x2 - (mu_r_ + i*mu_i_)*e2) ...
  + 0.5*conj(y3 - nu_*x3 - (mu_r_ + i*mu_i_)*e3).*(y3 - nu_*x3 - (mu_r_ + i*mu_i_)*e3) ...
  + 0.5*conj(y4 - nu_*x4 - (mu_r_ + i*mu_i_)*e4).*(y4 - nu_*x4 - (mu_r_ + i*mu_i_)*e4) ...
  + 0.5*conj(y5 - nu_*x5 - (mu_r_ + i*mu_i_)*e5).*(y5 - nu_*x5 - (mu_r_ + i*mu_i_)*e5) ...
  + 0.5*conj(y6 - nu_*x6 - (mu_r_ + i*mu_i_)*e6).*(y6 - nu_*x6 - (mu_r_ + i*mu_i_)*e6) ...
  ;
%%%%%%%%;
% Here we calculate the inner-products between the ctf*template, image and mask. ;
%%%%%%%%;
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
%%%%%%%%;
% We can use these inner-products to define ssnll_nu_function... ;
%%%%%%%%;
I_A = -real(I_1x)^2 - imag(I_1x)^2 + I_xx*I_11;
I_2B = -real(I_1x)*real(I_1y) - imag(I_1x)*imag(I_1y) + real(I_xy)*I_11;
I_C = -real(I_1y)^2 - imag(I_1y)^2 + I_yy*I_11;
ssnll_nu_function_0 = @(nu_) ... %<-- one version of the ssnll_nu_function. ;
  0.5*(I_A*nu_.^2 -2*I_2B*nu_ + I_C)./I_11 ...
  ;
%%%%%%%%;
% ... but this definition is not elegant, and I prefer the one below, ;
% involving the centered ctf*template and image. ;
%%%%%%%%;
ybar_ = I_1y*e_/max(1e-12,I_11); %<-- ctf*template average. ;
xbar_ = I_1x*e_/max(1e-12,I_11); %<-- image average. ;
ycen_ = y_ - ybar_; %<-- value-centered ctf*template. ;
xcen_ = x_ - xbar_; %<-- value-centered image. ;
I_ycenycen = sum(conj(ycen_).*ycen_,'all');
I_xcenycen = sum(conj(xcen_).*ycen_,'all');
I_ycenxcen = sum(conj(ycen_).*xcen_,'all');
I_xcenxcen = sum(conj(xcen_).*xcen_,'all');
ssnll_nu_function_1 = @(nu_) ... %<-- another (more expressive) representation of the ssnll_nu_function. ;
  0.5*(I_ycenycen + nu_.^2*I_xcenxcen - nu_*I_xcenycen - nu_*I_ycenxcen) ...
  ;
%%%%%%%%;
% Now we normalize the centered ctf*template and image... ;
%%%%%%%%;
ynrm_ = ycen_/max(1e-12,sqrt(I_ycenycen)); %<-- normalized ctf*template. ;
xnrm_ = xcen_/max(1e-12,sqrt(I_xcenxcen)); %<-- normalized image. ;
I_ynrmxnrm = sum(conj(ynrm_).*xnrm_,'all');
I_xnrmynrm = sum(conj(xnrm_).*ynrm_,'all');
%%%%%%%%;
% ... and calculate measures of similarity and dis-similarity. ;
%%%%%%%%;
I_costheta = 0.5*(I_xnrmynrm + I_ynrmxnrm); %<-- similarity measure. ;
I_sinthetasquared = 1 - I_costheta.^2; %<-- dis-similarity measure. ;
%%%%%%%%;
% This dis-similarity measure can be used to represent the ssnll_function. ;
%%%%%%%%;
ssnll_2 = 0.5*I_ycenycen*I_sinthetasquared ; %<-- a simple representation of the ssnll_function. ;
%%%%%%%%;
% Now we set up a polar-quadrature grid and weights for mu_ := [mu*cos(psi);mu*sin(psi)]. ;
% The range and precision of this quadrature-grid is chosen so that the likelihood_mu_nu_sigma is accurately integrated. ;
%%%%%%%%;
mu_int = 8; %<-- determines the maximum radius in mu-space. ;
mu_eq_d_double = 1/16; %<-- determines the local-resolution in mu-space. ;
n_psi_int = 8; %<-- determines the angular-resolution in mu-space. ;
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
% Here we extract the integration-weights associated with the polar-quadrature-grid above. ;
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
mu_psimu_ = mu_c_0_psimu_ + i*mu_c_1_psimu_ ; %<-- Now we list the mu-values on the quadrature-grid. ;
%%%%%%%%;
% Now we define the likelihood_mu_nu_sigma_functions in the abstract, as well as on the quadrature-grid (the latter for subsequent numerical integration). ;
%%%%%%%%;
likelihood_mu_nu_sigma_function_0 = @(mu_r_,mu_i_,nu_,sigma) exp(-ssnll_mu_nu_function_0(mu_r_,mu_i_,nu_)./sigma.^2) .* (2*pi*sigma.^2).^(-n_pixel/2) ; %<-- here is the standard version of the likelihood_mu_nu_sigma_function. ;
likelihood_mu_nu_sigma_psimu_ = likelihood_mu_nu_sigma_function_0(mu_c_0_psimu_,mu_c_1_psimu_,nu,sigma);
%%%%%%%%;
% If we want, we can check to make sure that the integrand likelihood_mu_nu_sigma is captured by the quadrature-grid. ;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
likelihood_lim_ = [0,max(abs(likelihood_mu_nu_sigma_psimu_))];
imagesc_p(n_mu_p_r,mu_p_r_,n_psi_,n_psi_sum,likelihood_mu_nu_sigma_psimu_,likelihood_lim_,colormap_81s); %<-- The integrand integrand likelihood_mu_nu_sigma should look well-resolved to the eye. ;
axis image; axisnotick;
title('likelihood_mu_nu_sigma_psimu_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
% Now we use numerical-quadrature to calculate the likelihood after marginalizing over mu (for a particular value of nu and sigma). ;
%%%%%%%%;
likelihood_nu_sigma_function_0 = @(nu_,sigma_) exp(-ssnll_nu_function_0(nu_)./sigma_.^2) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11) ; %<-- Here is one version of the likelihood_nu_sigma_function. ;
likelihood_nu_sigma_function_1 = @(nu_,sigma_) exp(-ssnll_nu_function_1(nu_)./sigma_.^2) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11) ; %<-- Here is a slightly different version of the likelihood_nu_sigma_function which uses the more expressive ssnll_nu_function. ;
likelihood_nu_sigma_auto_0 = integral2( ...
 @(mu_r,mu_i) likelihood_mu_nu_sigma_function_0(mu_r,mu_i,nu,sigma) ...
,-20,+20,-20,+20 ...
,'AbsTol',1e-9 ...
,'RelTol',1e-9 ...
); %<-- built-in adaptive 2d-quadrature is not very accurate. ;
likelihood_nu_sigma_quad_0 = sum(likelihood_mu_nu_sigma_psimu_.*weight_2d_psimu_,'all')*(2*pi)^2; %<-- our hand-picked polar quadrature is much more accurate. ;
likelihood_nu_sigma_form_0 = likelihood_nu_sigma_function_0(nu,sigma);
likelihood_nu_sigma_form_1 = likelihood_nu_sigma_function_1(nu,sigma);
%%%%%%%%;
% Finally, we compare the results of numerical-quadrature to the analytical formulae. ;
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_auto_0',likelihood_nu_sigma_auto_0,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_quad_0',likelihood_nu_sigma_quad_0,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_0',likelihood_nu_sigma_form_0,'likelihood_nu_sigma_form_1',likelihood_nu_sigma_form_1,' %<-- should be zero');
%%%%%%%%;
% Now we set up a simple 1d-quadrature to integrate over nu (for a particular value of sigma). ;
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
%%%%%%%%;
% We use numerical-quadrature to calculate the likelihood after marginalizing over nu. ;
%%%%%%%%;
likelihood_sigma_auto_1 = integral(@(nu_) likelihood_nu_sigma_function_1(nu_,sigma),-10,+10); %<-- Here we have a one-dimensional integral, which matlab should be able to handle. ;
likelihood_sigma_quad_1 = sum(bsxfun(@times,likelihood_nu_sigma_function_1(nu_value_q_,sigma),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form_1 = sqrt(2*pi*sigma^2*I_11/I_A)*exp((I_2B^2/(2*sigma^2*I_11*I_A)))*exp(-I_C/(2*sigma^2*I_11)) .* (2*pi*sigma.^2).^(-n_pixel/2) .* (2*pi*sigma.^2./I_11);
likelihood_sigma_form_2 = sqrt(2*pi*sigma^2/I_xcenxcen)*exp(-ssnll_2/(sigma^2)) .* (2*pi*sigma.^2).^(-n_pixel/2) .* (2*pi*sigma.^2./I_11);
%%%%%%%%;
% We compare the results of numerical-quadrature to the analytical formulae. ;
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_auto_1',likelihood_sigma_auto_1,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_sigma_form_2',likelihood_sigma_form_2,'likelihood_sigma_auto_1',likelihood_sigma_auto_1,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_quad_1',likelihood_sigma_quad_1,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_sigma_form_2',likelihood_sigma_form_2,'likelihood_sigma_quad_1',likelihood_sigma_quad_1,' %<-- should be small');
fnorm_disp(flag_verbose,'likelihood_sigma_form_1',likelihood_sigma_form_1,'likelihood_sigma_form_2',likelihood_sigma_form_2,' %<-- should be zero');
%%%%%%%%;
% Now we set up a simple 1d-quadrature to integrate over sigma. ;
%%%%%%%%;
sigma_lim_ = 16*[0,+1];
n_sigma_q = 256;
[sigma_value_q_,sigma_weight_q_] = legpts(n_sigma_q,sigma_lim_);
%%%%%%%%;
% We use numerical-quadrature to calculate the likelihood after marginalizing over sigma. ;
%%%%%%%%;
likelihood_sigma_function_2 = @(sigma_) sqrt(2*pi*sigma_.^2./I_xcenxcen).*exp(-ssnll_2./(sigma_.^2)) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_11); %<-- Here is an expressive version of likelihood_sigma. ;
likelihood_auto_2 = integral(@(sigma_) likelihood_sigma_function_2(sigma_),0,10); %<-- Here we have a one-dimensional integral, which matlab should be able to handle. ;
likelihood_quad_2 = sum(bsxfun(@times,likelihood_sigma_function_2(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form_2 = (2*pi)^(3/2-n_pixel/2) * (I_xcenxcen).^(-1/2) * 0.5 * (ssnll_2).^(2 - n_pixel/2) * gamma(n_pixel/2 - 2) / I_11;
%%%%%%%%;
% We compare the results of numerical-quadrature to the analytical formulae. ;
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_form_2',likelihood_form_2,'likelihood_auto_2',likelihood_auto_2,' %<-- can be large');
fnorm_disp(flag_verbose,'likelihood_form_2',likelihood_form_2,'likelihood_quad_2',likelihood_quad_2,' %<-- should be small');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Here we test the formulae when y_, x_, e_ and mu are strictly real. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_pixel = 7;
%%%%%%%%;
y_ = 0.25*randn(n_pixel,1); %<-- ctf*template. ;
y0 = y_(1+0);
y1 = y_(1+1);
y2 = y_(1+2);
y3 = y_(1+3);
y4 = y_(1+4);
y5 = y_(1+5);
y6 = y_(1+6);
%%%%%%%%;
x_ = 0.25*randn(n_pixel,1); %<-- image. ;
x0 = x_(1+0);
x1 = x_(1+1);
x2 = x_(1+2);
x3 = x_(1+3);
x4 = x_(1+4);
x5 = x_(1+5);
x6 = x_(1+6);
%%%%%%%%;
sigma = 0.3; %<-- temperature. ;
mu_r = -0.2;
mu = mu_r; %<-- offset. ;
nu = 0.37; %<-- intensity. ;
%%%%%%%%;
e_ = randn(n_pixel,1); %<-- more complicated mask used for testing. ;
e0 = e_(1+0);
e1 = e_(1+1);
e2 = e_(1+2);
e3 = e_(1+3);
e4 = e_(1+4);
e5 = e_(1+5);
e6 = e_(1+6);
%%%%%%%%;
ssnll_mu_nu_function_0 = @(mu_r_,nu_) ... %<-- a standard definition for the ssnll_mu_nu_function. ;
  + 0.5*conj(y0 - nu_*x0 - (mu_r_)*e0).*(y0 - nu_*x0 - (mu_r_)*e0) ...
  + 0.5*conj(y1 - nu_*x1 - (mu_r_)*e1).*(y1 - nu_*x1 - (mu_r_)*e1) ...
  + 0.5*conj(y2 - nu_*x2 - (mu_r_)*e2).*(y2 - nu_*x2 - (mu_r_)*e2) ...
  + 0.5*conj(y3 - nu_*x3 - (mu_r_)*e3).*(y3 - nu_*x3 - (mu_r_)*e3) ...
  + 0.5*conj(y4 - nu_*x4 - (mu_r_)*e4).*(y4 - nu_*x4 - (mu_r_)*e4) ...
  + 0.5*conj(y5 - nu_*x5 - (mu_r_)*e5).*(y5 - nu_*x5 - (mu_r_)*e5) ...
  + 0.5*conj(y6 - nu_*x6 - (mu_r_)*e6).*(y6 - nu_*x6 - (mu_r_)*e6) ...
  ;
%%%%%%%%;
I_yy = sum(conj(y_).*y_,'all');
I_xy = sum(conj(x_).*y_,'all');
I_1y = sum(conj(e_).*y_,'all');
I_yx = sum(conj(y_).*x_,'all');
I_xx = sum(conj(x_).*x_,'all');
I_1x = sum(conj(e_).*x_,'all');
I_y1 = sum(conj(y_).*e_,'all');
I_x1 = sum(conj(x_).*e_,'all');
I_11 = sum(conj(e_).*e_,'all');
%%%%%%%%;
ybar_ = I_1y*e_/max(1e-12,I_11); %<-- ctf*template average. ;
xbar_ = I_1x*e_/max(1e-12,I_11); %<-- image average. ;
ycen_ = y_ - ybar_; %<-- value-centered ctf*template. ;
xcen_ = x_ - xbar_; %<-- value-centered image. ;
I_ycenycen = sum(conj(ycen_).*ycen_,'all');
I_xcenycen = sum(conj(xcen_).*ycen_,'all');
I_ycenxcen = sum(conj(ycen_).*xcen_,'all');
I_xcenxcen = sum(conj(xcen_).*xcen_,'all');
ssnll_nu_function_1 = @(nu_) ... %<-- another (more expressive) representation of the ssnll_nu_function. ;
  0.5*(I_ycenycen + nu_.^2*I_xcenxcen - nu_*I_xcenycen - nu_*I_ycenxcen) ...
  ;
%%%%%%%%;
ynrm_ = ycen_/max(1e-12,sqrt(I_ycenycen)); %<-- normalized ctf*template. ;
xnrm_ = xcen_/max(1e-12,sqrt(I_xcenxcen)); %<-- normalized image. ;
I_ynrmxnrm = sum(conj(ynrm_).*xnrm_,'all');
I_xnrmynrm = sum(conj(xnrm_).*ynrm_,'all');
%%%%%%%%;
I_costheta = 0.5*(I_xnrmynrm + I_ynrmxnrm); %<-- similarity measure. ;
I_sinthetasquared = 1 - I_costheta.^2; %<-- dis-similarity measure. ;
%%%%%%%%;
ssnll_2 = 0.5*I_ycenycen*I_sinthetasquared ; %<-- a simple representation of the ssnll_function. ;
%%%%%%%%;
% Now we set up a polar-quadrature grid and weights for mu_ := mu_r_. ;
%%%%%%%%;
mu_lim_ = 5.0*[-1,+1];
n_mu_q = 256*8;
[mu_value_q_,mu_weight_q_] = legpts(n_mu_q,mu_lim_);
%%%%%%%%;
likelihood_mu_nu_sigma_function_0 = @(mu_r_,nu_,sigma) exp(-ssnll_mu_nu_function_0(mu_r_,nu_)./sigma.^2) .* (2*pi*sigma.^2).^(-n_pixel/2) ; %<-- here is the standard version of the likelihood_mu_nu_sigma_function. ;
likelihood_mu_nu_sigma_q_ = likelihood_mu_nu_sigma_function_0(mu_value_q_,nu,sigma);
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
plot(mu_value_q_,likelihood_mu_nu_sigma_q_,'o');
xlabel('mu');ylabel('likelihood');
title('likelihood_mu_nu_sigma_q_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
likelihood_nu_sigma_function_1 = @(nu_,sigma_) exp(-ssnll_nu_function_1(nu_)./sigma_.^2) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* sqrt(2*pi*sigma_.^2./I_11) ; %<-- Here is a slightly different version of the likelihood_nu_sigma_function which uses the more expressive ssnll_nu_function. ;
likelihood_nu_sigma_quad_0 = sum(likelihood_mu_nu_sigma_q_.*reshape(mu_weight_q_,[n_mu_q,1]),'all');
likelihood_nu_sigma_form_1 = likelihood_nu_sigma_function_1(nu,sigma);
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form_1',likelihood_nu_sigma_form_1,'likelihood_nu_sigma_quad_0',likelihood_nu_sigma_quad_0,' %<-- should be small');
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
%%%%%%%%;
likelihood_sigma_quad_1 = sum(bsxfun(@times,likelihood_nu_sigma_function_1(nu_value_q_,sigma),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form_2 = sqrt(2*pi*sigma^2/I_xcenxcen)*exp(-ssnll_2/(sigma^2)) .* (2*pi*sigma.^2).^(-n_pixel/2) .* sqrt(2*pi*sigma.^2./I_11);
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_sigma_form_2',likelihood_sigma_form_2,'likelihood_sigma_quad_1',likelihood_sigma_quad_1,' %<-- should be small');
%%%%%%%%;
sigma_lim_ = 16*[0,+1];
n_sigma_q = 256;
[sigma_value_q_,sigma_weight_q_] = legpts(n_sigma_q,sigma_lim_);
%%%%%%%%;
likelihood_sigma_function_2 = @(sigma_) sqrt(2*pi*sigma_.^2./I_xcenxcen).*exp(-ssnll_2./(sigma_.^2)) .* (2*pi*sigma_.^2).^(-n_pixel/2) .* sqrt(2*pi*sigma_.^2./I_11); %<-- Here is an expressive version of likelihood_sigma. ;
likelihood_quad_2 = sum(bsxfun(@times,likelihood_sigma_function_2(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form_2 = (2*pi)^(2/2-n_pixel/2) * (I_xcenxcen).^(-1/2) * 0.5 * (ssnll_2).^(3/2 - n_pixel/2) * gamma(n_pixel/2 - 3/2) / sqrt(I_11);
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_form_2',likelihood_form_2,'likelihood_quad_2',likelihood_quad_2,' %<-- should be small');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now we set up a test with a fixed image and template. ;
% Here the image will be denoted by "M", the template by "S" and the CTF by "CTF". ;
% The CTF*S will be denoted by "T". ;
% The mask will be denoted by "E". ;
% We first step through the formulae for integration over mu complex. ;
%%%%%%%%;
% Note that, in this example, we define the image to be a plane-wave multiplied by a planar-function. ;
% We also define the template to be a plane-wave, and the CTF to be a planar-function. ;
% This means that T will be a plane-wave multiplied by a planar-function. ;
% The mask is also a plane-wave. ;
% These choices mean that we will be able to use the formulae in: ;
% I_xPPx_0, I_xPP_0 and I_PP_0 ;
% to calculate the inner-products between T, M and E. ;
% These analytically-calculated values can then be used in the formula for the likelihood (ssnll_function below). ;
%%%%%%%%;
% Ultimately, we should set up a similar test in CryoLike. ;
% Note that a full test of this kind will require sorting out the effects of the translation and in-plane rotation, ;
% in terms of how they affect the relative orientation of the CTF, template and image (see test_cross_correlation_2.py). ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
% First we set up polar-grid for quadrature over the image-domain (i.e., typically frequency-space). ;
% The range and precision of this quadrature-grid is chosen so that the plane-waves below can be accurately integrated. ;
%%%%%%%%;
k_int = 8; %<-- this determines the maximum radius in frequency-space. ;
k_eq_d_double = 0.250; %<-- this determines the local-resolution in frequency-space. ;
n_w_int = 4; %<-- this determines the angular-resolution in frequency-space. ;
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
% Here we extract the integration-weights associated with the polar-quadrature-grid above. ;
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
% Now we define a particular set of template S, CTF, image M and mask E. ;
% For this definition we make sure to pick a relatively low spatial-frequency for the plane-waves involved. ;
% This is because we want these plane-waves to be resolved on the quadrature-grid above (involve frequency-space). ;
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
% If we want, we can check to see that the objects T_, M_ and E_ are well-resolved in frequency-space. ;
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
% Now we use numerical-quadrature, as well as the formulae in I_xPPx_0, I_xPP_0 and I_PP_0 to (analytically) calculate the inner-products between T_, M_ and E_. ;
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
% Now we quickly check to see if our formula for ssnll_mu_nu is correct. ;
% Note that, in this example, mu is complex. ;
%%%%%%%%;
nu = 1.3; mu_r = 0.4; mu_i = -0.7; mu = mu_r + i*mu_i; lu = 1.0; lu_ = 1.0;
tmp_wk_ = nu*M_k_p_wk_ + mu*E_k_p_wk_ - T_k_p_wk_ ;
ssnll_mu_nu_quad_0 = 0.5*sum(conj(tmp_wk_).*tmp_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
ssnll_mu_nu_function_0 = ...
@(mu_,nu_) ...
 + 0.5*conj(nu_).*nu_.*I_M_M_form ...
 + 0.5*conj(mu_).*mu_.*I_E_E_form ...
 + 0.5*conj(lu_).*lu_.*I_T_T_form ...
 + 1.0*real(conj(nu_).*mu_.*I_M_E_form) ...
 - 1.0*real(conj(nu_).*lu_.*I_T_M_form) ...
 - 1.0*real(conj(lu_).*mu_.*I_T_E_form) ...
;
ssnll_mu_nu_form_0 = ssnll_mu_nu_function_0(mu,nu);
fnorm_disp(flag_verbose,'ssnll_mu_nu_quad_0',ssnll_mu_nu_quad_0,'ssnll_mu_nu_form_0',ssnll_mu_nu_form_0,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now we set up and test quadrature for integration over complex mu, real nu and positive sigma. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
% First we define the inner-products that involve formulas. ;
% Note that here quantities such as I_Mcen_Mcen_form and I_sinthetasquared can be computed analytically from the prescription of M_, S_, CTF_ and E_. ;
%%%%%%%%;
n_pixel = I_E_E_form; %<-- It is not unreasonable to expect that the number of relevant pixels is given by the l2-norm of the mask. ;
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
% Once again, we set up a polar-quadrature grid and weights for mu_ := [mu*cos(psi);mu*sin(psi)]. ;
% The range and precision of this quadrature-grid is again chosen so that the likelihood_mu_nu_sigma is accurately integrated. ;
%%%%%%%%;
mu_int = 24; %<-- Note the increased maximum-radius in mu-space. ;
mu_eq_d_double = 1/16; %<-- determines the local-resolution in mu-space. ;
n_psi_int = 8; %<-- determines the angular-resolution in mu-space. ;
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
% Again we extract the integration-weights associated with the polar-quadrature-grid above. ;
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
% Now we define formula for integration over mu. ;
%%%%%%%%;
ssnll_nu_function = @(nu_,sigma_) ... %<-- This is the expressive version of ssnll_nu. ;
  0.5*(I_Tcen_Tcen_form + nu_.^2*I_Mcen_Mcen_form - nu_*I_Tcen_Mcen_form - nu_*I_Mcen_Tcen_form) ...
  ;
likelihood_nu_sigma_function = @(nu_,sigma_) exp(-ssnll_nu_function(nu_,sigma_)./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* (2*pi*sigma_.^2./I_E_E_form) ; %<-- This is the expressive version of likelihood_nu_sigma. ;
%%%%%%%%;
% Now we test integration over mu for a particular value of nu and sigma. ;
%%%%%%%%;
sigma_value = 1.3; nu_value = 0.8;
likelihood_mu_nu_sigma_psimu_ = (2*pi*sigma_value^2).^(-n_pixel/2) * exp(-ssnll_mu_nu_function_0(mu_psimu_,nu_value)/max(1e-12,sigma_value^2)) ;
%%%%;
% if we like, we can make sure that the integrand likelihood_mu_nu_sigma is well resolved by the quadrature-grid in mu-space. ;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
l_lim_ = [0,max(abs(likelihood_mu_nu_sigma_psimu_))];
imagesc_p(n_mu_p_r,mu_p_r_,n_psi_,n_psi_sum,likelihood_mu_nu_sigma_psimu_,l_lim_,colormap_81s);
axis image; axisnotick;
title('likelihood_mu_nu_sigma_psimu_','Interpreter','none');
end;%if flag_disp;
likelihood_nu_sigma_quad = sum(bsxfun(@times,likelihood_mu_nu_sigma_psimu_,reshape(weight_2d_psimu_,[n_psi_sum,1])),'all')*(2*pi)^2 ;
likelihood_nu_sigma_form = likelihood_nu_sigma_function(nu_value,sigma_value);
%%%%;
% Now we compare the numerical-quadrature to our formula. ;
%%%%;
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form',likelihood_nu_sigma_form,'likelihood_nu_sigma_quad',likelihood_nu_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
% Now we test integration over nu. ;
% We first define a quadrature-grid for nu... ;
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
%%%%%%%%;
% ... then define formula for integration over nu. ;
%%%%%%%%;
ssnll_function = ...
  0.5*I_Tcen_Tcen_form*I_sinthetasquared ...
  ;
likelihood_sigma_function = @(sigma_) sqrt(2*pi*sigma_.^2./max(1e-12,I_Mcen_Mcen_form)) .* exp(-ssnll_function./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* (2*pi*sigma_.^2./max(1e-12,I_E_E_form)); %<-- This is the expressive version of likelihood_sigma. ;
%%%%%%%%;
% Now we test integration over nu for a particular value of sigma... ;
%%%%%%%%;
sigma_value = 1.5;
likelihood_sigma_quad = sum(bsxfun(@times,likelihood_nu_sigma_function(nu_value_q_,sigma_value),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form = likelihood_sigma_function(sigma_value);
%%%%%%%%;
% ... and compare the numerical-quadrature to the formula. ;
%%%%%%%%;
fnorm_disp(flag_verbose,'likelihood_sigma_form',likelihood_sigma_form,'likelihood_sigma_quad',likelihood_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
% Now we test integration over sigma. ;
% We first define a quadrature-grid for sigma... ;
%%%%%%%%;
n_sigma_q_sub = 8; [sigma_value_q_sub_,sigma_weight_q_sub_] = legpts(n_sigma_q_sub,[0,1]);
n_sigma_max = 1024*8; sigma_lim_ = n_sigma_max*[0,+1];
n_sigma_q = n_sigma_q_sub*n_sigma_max;
sigma_value_q_ = reshape(bsxfun(@plus,sigma_value_q_sub_,[0:n_sigma_max-1]),[n_sigma_q,1]);
sigma_weight_q_ = repmat(sigma_weight_q_sub_,[1,n_sigma_max]);
%%%%%%%%;
% ... then define formula for integration over sigma. ;
%%%%%%%%;
likelihood_function = ...
1 ...
* (2*pi)^(3/2-n_pixel/2) ...
* (max(1e-12,I_Mcen_Mcen_form)).^(-1/2) ...
* 0.5 ...
* (max(1e-12,ssnll_function)).^(2 - n_pixel/2) ...
* gamma(n_pixel/2 - 2) ...
/ max(1e-12,I_E_E_form) ...
;
%%%%%%%%;
% Now we test integration over sigma. ;
%%%%%%%%;
likelihood_quad = sum(bsxfun(@times,likelihood_sigma_function(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form = likelihood_function;
fnorm_disp(flag_verbose,'likelihood_form',likelihood_form,'likelihood_quad',likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');
%%%%%%%%;
% And we also test the formula for the log-likelihood. ;
%%%%%%%%;
negative_log_likelihood_function = ...
-(3/2 - n_pixel/2) * log(2*pi) ...
+(1/2) * log(max(1e-12,I_Mcen_Mcen_form)) ...
+(1.0) * log(2) ...
-(2 - n_pixel/2) * log(max(1e-12,ssnll_function)) ...
- gammaln(n_pixel/2 - 2) ...
+ log(max(1e-12,I_E_E_form)) ...
;
negative_log_likelihood_quad = -log(likelihood_quad);
negative_log_likelihood_form = negative_log_likelihood_function;
fnorm_disp(flag_verbose,'negative_log_likelihood_form',negative_log_likelihood_form,'negative_log_likelihood_quad',negative_log_likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now we set up a similar test, but with a fixed image and template and real mask (with real mu). ;
% We re-use the settings established above for the image-template test. ;
% Note that we can keep the same mask as the plane-wave mask used above, ;
% since can still integrate over nu real (indeed, E_x_c_xx__ should be real). ;
%%%%%%%%;
% However, we have to redefine the likelihood to account for real mu. ;
%%%%%%%%;
nu = 1.3; mu_r = 0.4; mu_i = 0; mu = mu_r + i*mu_i; lu = 1.0; lu_ = 1.0;
tmp_wk_ = nu*M_k_p_wk_ + mu*E_k_p_wk_ - T_k_p_wk_ ;
ssnll_mu_nu_quad_0 = 0.5*sum(conj(tmp_wk_).*tmp_wk_.*weight_2d_wk_,'all')*(2*pi)^2;
ssnll_mu_nu_function_0 = ...
@(mu_,nu_) ...
 + 0.5*conj(nu_).*nu_.*I_M_M_form ...
 + 0.5*conj(mu_).*mu_.*I_E_E_form ...
 + 0.5*conj(lu_).*lu_.*I_T_T_form ...
 + 1.0*real(conj(nu_).*mu_.*I_M_E_form) ...
 - 1.0*real(conj(nu_).*lu_.*I_T_M_form) ...
 - 1.0*real(conj(lu_).*mu_.*I_T_E_form) ...
;
ssnll_mu_nu_form_0 = ssnll_mu_nu_function_0(mu,nu);
fnorm_disp(flag_verbose,'ssnll_mu_nu_quad_0',ssnll_mu_nu_quad_0,'ssnll_mu_nu_form_0',ssnll_mu_nu_form_0,' %<-- should be small (e.g., <1e-6 if quadrature resolved)');
%%%%%%%%;
mu_int = 24;
mu_lim_ = mu_int*[-1,+1];
n_mu_q = 256*8;
[mu_value_q_,mu_weight_q_] = legpts(n_mu_q,mu_lim_);
%%%%%%%%;
ssnll_nu_function = @(nu_,sigma_) ... %<-- This is the expressive version of ssnll_nu. ;
  0.5*(I_Tcen_Tcen_form + nu_.^2*I_Mcen_Mcen_form - nu_*I_Tcen_Mcen_form - nu_*I_Mcen_Tcen_form) ...
  ;
likelihood_nu_sigma_function = @(nu_,sigma_) exp(-ssnll_nu_function(nu_,sigma_)./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* sqrt(2*pi*sigma_.^2./I_E_E_form) ; %<-- This is the expressive version of likelihood_nu_sigma. ;
%%%%%%%%;
sigma_value = 1.3; nu_value = 0.8;
likelihood_mu_nu_sigma_q_ = (2*pi*sigma_value^2).^(-n_pixel/2) * exp(-ssnll_mu_nu_function_0(mu_value_q_,nu_value)/max(1e-12,sigma_value^2)) ;
%%%%;
likelihood_nu_sigma_quad = sum(bsxfun(@times,likelihood_mu_nu_sigma_q_,reshape(mu_weight_q_,[n_mu_q,1])),'all')*(2*pi)^2 ;
likelihood_nu_sigma_form = likelihood_nu_sigma_function(nu_value,sigma_value);
fnorm_disp(flag_verbose,'likelihood_nu_sigma_form',likelihood_nu_sigma_form,'likelihood_nu_sigma_quad',likelihood_nu_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
nu_lim_ = 16*[-1,+1];
n_nu_q = 256;
[nu_value_q_,nu_weight_q_] = legpts(n_nu_q,nu_lim_);
%%%%%%%%;
ssnll_function = ...
  0.5*I_Tcen_Tcen_form*I_sinthetasquared ...
  ;
likelihood_sigma_function = @(sigma_) sqrt(2*pi*sigma_.^2./max(1e-12,I_Mcen_Mcen_form)) .* exp(-ssnll_function./max(1e-12,sigma_.^2)) .* (2*pi*max(1e-12,sigma_.^2)).^(-n_pixel/2) .* sqrt(2*pi*sigma_.^2./max(1e-12,I_E_E_form)); %<-- This is the expressive version of likelihood_sigma. ;
%%%%%%%%;
sigma_value = 1.5;
likelihood_sigma_quad = sum(bsxfun(@times,likelihood_nu_sigma_function(nu_value_q_,sigma_value),reshape(nu_weight_q_,[n_nu_q,1])),'all');
likelihood_sigma_form = likelihood_sigma_function(sigma_value);
fnorm_disp(flag_verbose,'likelihood_sigma_form',likelihood_sigma_form,'likelihood_sigma_quad',likelihood_sigma_quad,' %<-- should be small as long as quadrature includes support');
%%%%%%%%;
n_sigma_q_sub = 8; [sigma_value_q_sub_,sigma_weight_q_sub_] = legpts(n_sigma_q_sub,[0,1]);
n_sigma_max = 1024*8; sigma_lim_ = n_sigma_max*[0,+1];
n_sigma_q = n_sigma_q_sub*n_sigma_max;
sigma_value_q_ = reshape(bsxfun(@plus,sigma_value_q_sub_,[0:n_sigma_max-1]),[n_sigma_q,1]);
sigma_weight_q_ = repmat(sigma_weight_q_sub_,[1,n_sigma_max]);
%%%%%%%%;
likelihood_function = ...
1 ...
* (2*pi)^(2/2-n_pixel/2) ...
* (max(1e-12,I_Mcen_Mcen_form)).^(-1/2) ...
* 0.5 ...
* (max(1e-12,ssnll_function)).^(3/2 - n_pixel/2) ...
* gamma(n_pixel/2 - 3/2) ...
/ max(1e-12,sqrt(I_E_E_form)) ...
;
%%%%%%%%;
likelihood_quad = sum(bsxfun(@times,likelihood_sigma_function(sigma_value_q_),reshape(sigma_weight_q_,[n_sigma_q,1])),'all');
likelihood_form = likelihood_function;
fnorm_disp(flag_verbose,'likelihood_form',likelihood_form,'likelihood_quad',likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');
%%%%%%%%;
negative_log_likelihood_function = ...
-(2/2 - n_pixel/2) * log(2*pi) ...
+(1/2) * log(max(1e-12,I_Mcen_Mcen_form)) ...
+(1.0) * log(2) ...
-(3/2 - n_pixel/2) * log(max(1e-12,ssnll_function)) ...
- gammaln(n_pixel/2 - 3/2) ...
+ log(max(1e-12,sqrt(I_E_E_form))) ...
;
negative_log_likelihood_quad = -log(likelihood_quad);
negative_log_likelihood_form = negative_log_likelihood_function;
fnorm_disp(flag_verbose,'negative_log_likelihood_form',negative_log_likelihood_form,'negative_log_likelihood_quad',negative_log_likelihood_quad,' %<-- can be difficult to resolve due to algebraic decay');
%%%%%%%%;
