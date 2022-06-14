function ...
[ ...
 X0_ ...
,n_op ...
,n_mult ...
,alpha_ ...
,gamma_ ...
] = ...
register_spharm_to_spharm_single_beta_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_ ...
,b_ ...
,beta ...
,n_alpha ...
,alpha_ ...
,n_gamma ...
,gamma_ ...
,a___ ...
,b___ ...
,ab__ ...
,d__ ...
);
% registration between molecule_A and molecule_B using a single beta (fast only);
% ;
% verbose = integer verbosity_level ;
% n_k_p_r = integer maximum k ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk) = k_p_r_value for shell nk ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% beta = beta angle ;
% n_alpha = integer number of alpha angles (optional);
% alpha_ = real array of alpha angles (optional);
% n_gamma = integer number of gamma angles (optional);
% gamma_ = real array of gamma angles (optional);
% ;
% If arrays alpha_ and gamma_ are not provided we use the standard arrays: ;
% alpha_ = linspace(0,2*pi,n_m_max+1); alpha_ = alpha_(1:end-1);
% gamma_ = linspace(0,2*pi,n_m_max+1); gamma_ = gamma_(1:end-1);
% and use the standard fft to calculate X0_.; 
% ;
% However, if n_alpha, alpha_ and n_gamma, gamma_ are provided, we use these arrays instead, ;
% and use the nufft to calculate X0_. ;
% ;
% X_ = complex array of size (n_alpha,n_gamma,1);
% The default values of n_alpha and n_gamma are n_m_max. ;
% X_(nalpha,ngamma,nbeta) corresponds to the innerproduct between molecule_A and molecule_B, ;
% where one of the molecules has been rotated by euler-angles alpha,beta,gamma. ;
% Note that in this case we specifically mean that either: ;
% (i) rotate_spharm_to_spharm_2 has been applied to b_ ;
% with euler_ = -[gamma , beta , alpha], or ;
% (ii) rotate_spharm_to_spharm_2 has been applied to a_ ;
% with euler_ = +[alpha , beta , gamma] ;
% Note that alpha_ and gamma_ are arrays from 0 to 2*pi, ;
% whereas beta_ is an array from -pi to pi. ;
%%%%%%%%;
% a___ = a_ rearranged so that mn varies quickly, then k, then l. ;
% b___ = b_ rearranged so that mp varies quickly, then k, then l. ;
% ab__ = \sum_{k} conj(a___)*b___ rearranged so that mn varies quickly, then mp, then l. ;
% d__ = W_ rearranged so that mn varies quickly, then mp, then l. ;

if (nargin<1);
verbose=1; nf=0;
if (verbose); disp(sprintf(' %% testing register_spharm_to_spharm_single_beta_3')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
%%%%%%%%;
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = half_diameter_x_c;
n_x_c = 64;
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
xxx_c_weight = (2*x_p_r_max/n_x_c)^3;
%%%%%%%%;
% Now set up and test k-quadrature on sphere. ;
%%%%%%%%;
verbose=0; k_p_r_max = 48.0/(2*pi); k_eq_d = 1.0/(2*pi); TorL = 'L';
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; n_plot = p_row*p_col;
for nplot=0:n_plot-1;
nk_p_r = max(0,min(n_k_p_r-1,round(n_k_p_r*nplot/n_plot)));
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
subplot(p_row,p_col,1+nplot);
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
plot3(k_c_0_all_(1+tmp_index_),k_c_1_all_(1+tmp_index_),k_c_2_all_(1+tmp_index_),'.');
axis vis3d; axisnotick;
end;%for nplot=0:n_plot-1;
end;%if flag_plot;
%%%%%%%%;
n_source = 4;
rng(0);
delta_a_c__ = zeros(3,n_source);
delta_b_c__ = zeros(3,n_source);
for nsource=0:n_source-1;
delta_a_c_ = 0.25*(2*rand(3,1)-1);
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = 0.25*(2*rand(3,1)-1);
delta_b_c__(:,1+nsource) = delta_b_c_;
end;%for nsource=0:n_source-1;
a_k_p_form_ = zeros(n_k_all,1);
b_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1) + k_c_2_all_*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
I_a_quad = sum(a_k_p_form_.*weight_3d_k_all_);
I_b_quad = sum(b_k_p_form_.*weight_3d_k_all_);
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h3d_(tmp_kd)*k_p_r_max^3;
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h3d_(tmp_kd)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
disp(sprintf(' %% I_a_form vs I_a_quad %0.16f',fnorm(I_a_form-I_a_quad)/fnorm(I_a_form)));
disp(sprintf(' %% I_b_form vs I_b_quad %0.16f',fnorm(I_b_form-I_b_quad)/fnorm(I_b_form)));
%%%%%%%%;
% now for nufft3d3: ;
% as per: https://github.com/zgimbutas/nufft3d/blob/master/nufft3d3.m ;
%%%%%%%%;
%  [FK,IER] = NUFFT3D3(NJ,XJ,YJ,ZJ,CJ,IFLAG,NK,SK,TK,UK);
%
%                  nj
%     fk(k)    =  SUM cj(j) exp(+/-i (s(k),t(k),u(k))*(xj(j),yj(j),zj(j)))
%                 j=1
%
%     If (iflag .ge.0) the + sign is used in the exponential.
%     If (iflag .lt.0) the - sign is used in the exponential.
%%%%%%%%;
% Thus, if we let: ;
% eta = pi/k_p_r_max, ;
% and the sum over j run over all quadrature nodes, ;
% xj, yj, zj = 2*pi*k_c_?_all_*eta, ;
% s, t, u = x_c_?___/eta, ;
% cj = a_k_p_form_*weight_3d_k_all_, ;
% then the sum approximates: ;
% f = \integral_{sphere of radius 2*pi*k_p_r_max} a_k_p_form_ .* exp(+i * < 2*pi*k_ , x_ > ) ;
% or, if a_k_p_form_ = exp(+i * <2*pi*k_ , delta_ > ), then ;
% f = \integral_{sphere of radius 2*pi*k_p_r_max} exp(+i * < 2*pi*k_ , x_ + delta_ > ) ;
% f = h3d_(tmp_kd)*k_p_r_max^3, ;
% where tmp_kd = 2*pi*k_p_r_max*fnorm(x_ + delta_). ;
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_form_.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_quad_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
b_x_c_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,b_k_p_form_.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: b_x_c_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
a_x_c_form_ = zeros(n_x_c,n_x_c,n_x_c);
b_x_c_form_ = zeros(n_x_c,n_x_c,n_x_c);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd___ = 2*pi * k_p_r_max * sqrt( (x_c_0___ + delta_a_c_(1+0)).^2 + (x_c_1___ + delta_a_c_(1+1)).^2 + (x_c_2___ + delta_a_c_(1+2)).^2 ) ;
a_x_c_form_ = a_x_c_form_ + h3d_(tmp_kd___)*k_p_r_max^3;
tmp_kd___ = 2*pi * k_p_r_max * sqrt( (x_c_0___ + delta_b_c_(1+0)).^2 + (x_c_1___ + delta_b_c_(1+1)).^2 + (x_c_2___ + delta_b_c_(1+2)).^2 ) ;
b_x_c_form_ = b_x_c_form_ + h3d_(tmp_kd___)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
disp(sprintf(' %% a_x_c_form_(:) vs a_x_c_quad_(:): %0.16f',fnorm(a_x_c_form_(:) - a_x_c_quad_(:))/fnorm(a_x_c_form_(:))));
disp(sprintf(' %% b_x_c_form_(:) vs b_x_c_quad_(:): %0.16f',fnorm(b_x_c_form_(:) - b_x_c_quad_(:))/fnorm(b_x_c_form_(:))));
%%%%%%%%;
% Now if we let: ;
% eta = pi/x_p_r_max, ;
% and the sum over j run over all x_c_, ;
% xj, yj, zj = x_c_?___*eta, ;
% s, t, u = 2*pi*k_c_?_all_/eta, ;
% cj = a_x_c_quad_ * xxx_c_weight, ;
% then the sum approximates: ;
% f = \integral_{cube of side-length 2} a_x_c_quad_ .* exp(-1 * < x_ , 2*pi*k_ > ) ;
% f = a_k_p_. ;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,a_x_c_quad_(:).*xxx_c_weight,-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_form_ vs a_k_p_quad_: %0.16f',fnorm(a_k_p_form_(:)-a_k_p_quad_)/fnorm(a_k_p_form_(:))));
eta = pi/x_p_r_max; tmp_t = tic;
b_k_p_quad_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,b_x_c_quad_(:).*xxx_c_weight,-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: b_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: b_k_p_form_ vs b_k_p_quad_: %0.16f',fnorm(b_k_p_form_(:)-b_k_p_quad_)/fnorm(b_k_p_form_(:))));
%%%%%%%%;
% Now convert a_k_p_form_ into a_k_Y_quad_. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
flag_Y_quad = input(' test convert_k_p_to_spharm_1 and convert_spharm_to_k_p_1 (default no): ');
if (isempty(flag_Y_quad)); flag_Y_quad = 0; end;
if flag_Y_quad;
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[b_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_quad_ time %0.2fs',tmp_t));
end;%if flag_Y_quad;
%%%%%%%%;
a_k_Y_form_ = zeros(n_lm_sum,1);
b_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c_,l_max_);
b_k_Y_form_ = b_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_b_c_,l_max_);
end;%for nsource=0:n_source-1;
if flag_Y_quad;
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_form_ - a_k_Y_quad_)/fnorm(a_k_Y_form_)));
disp(sprintf(' %% b_k_Y_form_ vs b_k_Y_quad_: %0.16f',fnorm(b_k_Y_form_ - b_k_Y_quad_)/fnorm(b_k_Y_form_)));
end;%if flag_Y_quad;
%%%%%%%%;
if flag_Y_quad;
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_form_-a_k_p_reco_)/fnorm(a_k_p_form_))); %<-- this should be 2-3 digits. ;
tmp_t = tic;
[b_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_quad_ --> b_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: b_k_p_reco error: %0.16f',fnorm(b_k_p_form_-b_k_p_reco_)/fnorm(b_k_p_form_))); %<-- this should be 2-3 digits. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;ns=0;
clf;
subplot(2,3,1+ns); ns=ns+1; plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(2,3,1+ns); ns=ns+1; plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(2,3,1+ns); ns=ns+1; plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(2,3,1+ns); ns=ns+1; plot(Y_l_val_,log10(abs(b_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(b))');
subplot(2,3,1+ns); ns=ns+1; plot(Y_m_val_,log10(abs(b_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(b))');
subplot(2,3,1+ns); ns=ns+1; plot(Y_k_val_,log10(abs(b_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(b))');
figbig;
end;%if flag_plot;
end;%if flag_Y_quad;
%%%%%%%%;
% Now register a against b. ;
%%%%%%%%;
beta = +pi/3;
[ ...
 X_quad_ ...
,n_op ...
,n_mult ...
,alpha_ ...
,gamma_ ...
] = ...
register_spharm_to_spharm_single_beta_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,b_k_Y_form_ ...
,beta ...
);
n_m_max = 1+2*l_max_max;
%%%%%%%%;
% test a few values. ;
%%%%%%%%;
n_test=3;
for ntest=0:n_test-1;
rng(ntest);
nalpha = max(0,min(n_m_max-1,floor(n_m_max*rand())));
alpha = alpha_(1+nalpha);
ngamma = max(0,min(n_m_max-1,floor(n_m_max*rand())));
gamma = gamma_(1+ngamma);
X_quad = X_quad_(1+nalpha,1+ngamma);
%%%%;
euler_ = -[gamma , beta , alpha];
[ ...
 b_k_Y_rota_ ...
] = ...
rotate_spharm_to_spharm_2( ...
 verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,b_k_Y_form_ ...
,euler_ ...
);
X_form = 0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm-1];
X_form = X_form + sum(conj(a_k_Y_form_(1+tmp_index_)).*b_k_Y_rota_(1+tmp_index_)) * weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
disp(sprintf(' %% ntest %d/%d: X_form - X_quad: %0.16f',ntest,n_test,fnorm(X_form - X_quad)/fnorm(X_form)));
%%%%;
euler_ = +[alpha , beta , gamma];
[ ...
 a_k_Y_rota_ ...
] = ...
rotate_spharm_to_spharm_2( ...
 verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,euler_ ...
);
X_form = 0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm-1];
X_form = X_form + sum(conj(a_k_Y_rota_(1+tmp_index_)).*b_k_Y_form_(1+tmp_index_)) * weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
disp(sprintf(' %% ntest %d/%d: X_form - X_quad: %0.16f',ntest,n_test,fnorm(X_form - X_quad)/fnorm(X_form)));
%%%%;
end;%for ntest=0:n_test-1;
%%%%%%%%;
% Now display the correlations X_quad_. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;figbeach();
Xlim_ = [min(abs(X_quad_)),max(abs(X_quad_))];
subplot(1,1,1); imagesc(real(X_quad_));
xlabel('alpha'); set(gca,'XTick',1:n_m_max,'XTickLabel',alpha_);
ylabel('gamma'); set(gca,'YTick',1:n_m_max,'YTickLabel',gamma_);
axis image;
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',6);
end;%if flag_plot;
%%%%%%%%;
% test again with xxnufft2d2.; 
%%%%%%%%;
n_alpha = 16;
alpha_ = 2*pi*rand(n_alpha,1);
n_gamma = 12;
gamma_ = 2*pi*rand(n_gamma,1);
beta = +pi/3;
[ ...
 X_quad_ ...
,n_op ...
,n_mult ...
,alpha_ ...
,gamma_ ...
] = ...
register_spharm_to_spharm_single_beta_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,b_k_Y_form_ ...
,beta ...
,n_alpha ...
,alpha_ ...
,n_gamma ...
,gamma_ ...
);
%%%%%%%%;
% test a few values. ;
%%%%%%%%;
n_test=3;
for ntest=0:n_test-1;
rng(ntest);
nalpha = max(0,min(n_alpha-1,floor(n_alpha*rand())));
alpha = alpha_(1+nalpha);
ngamma = max(0,min(n_gamma-1,floor(n_gamma*rand())));
gamma = gamma_(1+ngamma);
X_quad = X_quad_(1+nalpha,1+ngamma);
%%%%;
euler_ = -[gamma , beta , alpha];
[ ...
 b_k_Y_rota_ ...
] = ...
rotate_spharm_to_spharm_2( ...
 verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,b_k_Y_form_ ...
,euler_ ...
);
X_form = 0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm-1];
X_form = X_form + sum(conj(a_k_Y_form_(1+tmp_index_)).*b_k_Y_rota_(1+tmp_index_)) * weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
disp(sprintf(' %% ntest %d/%d: X_form - X_quad: %0.16f',ntest,n_test,fnorm(X_form - X_quad)/fnorm(X_form)));
%%%%;
euler_ = +[alpha , beta , gamma];
[ ...
 a_k_Y_rota_ ...
] = ...
rotate_spharm_to_spharm_2( ...
 verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,euler_ ...
);
X_form = 0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm-1];
X_form = X_form + sum(conj(a_k_Y_rota_(1+tmp_index_)).*b_k_Y_form_(1+tmp_index_)) * weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
disp(sprintf(' %% ntest %d/%d: X_form - X_quad: %0.16f',ntest,n_test,fnorm(X_form - X_quad)/fnorm(X_form)));
%%%%;
end;%for ntest=0:n_test-1;
%%%%%%%%
flag_plot=1;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;figbeach();
Xlim_ = [min(abs(X_quad_)),max(abs(X_quad_))];
subplot(1,1,1); imagesc(real(X_quad_));
xlabel('gamma'); set(gca,'XTick',1:n_gamma,'XTickLabel',gamma_);
ylabel('alpha'); set(gca,'YTick',1:n_alpha,'YTickLabel',alpha_);
axis image;
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',6);
end;%if flag_plot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp('returning'); return;
end;%if (nargin<1);

na=0;
if (nargin<1+na); verbose = []; end; na=na+1;
if (nargin<1+na); n_k_p_r = []; end; na=na+1;
if (nargin<1+na); k_p_r_ = []; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_ = []; end; na=na+1;
if (nargin<1+na); l_max_ = []; end; na=na+1;
if (nargin<1+na); a_ = []; end; na=na+1;
if (nargin<1+na); b_ = []; end; na=na+1;
if (nargin<1+na); beta = []; end; na=na+1;
if (nargin<1+na); n_alpha = []; end; na=na+1;
if (nargin<1+na); alpha_ = []; end; na=na+1;
if (nargin<1+na); n_gamma = []; end; na=na+1;
if (nargin<1+na); gamma_ = []; end; na=na+1;
if (nargin<1+na); a___ = []; end; na=na+1;
if (nargin<1+na); b___ = []; end; na=na+1;
if (nargin<1+na); ab__ = []; end; na=na+1;
if (nargin<1+na); d__ = []; end; na=na+1;

if (isempty(n_alpha)); n_alpha = 0; end;
if (isempty(n_gamma)); n_gamma = 0; end;

n_lm_ = (l_max_+1).^2; n_lm_csum_ = cumsum([0;n_lm_(:)]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);

if (n_alpha==0 | n_gamma==0);
fft_flag = 0; 
n_alpha = n_m_max;
n_gamma = n_m_max;
alpha_ = linspace(0,2*pi,n_m_max+1); alpha_ = alpha_(1:end-1);
gamma_ = linspace(0,2*pi,n_m_max+1); gamma_ = gamma_(1:end-1);
 else;
fft_flag = 1;
end;%else;

X0_ = zeros(n_alpha,n_gamma,1);
% Note that rotating a molecule by [+alpha,+beta,+gamma] ;
% corresponds to rotating the coordinate-frame by [-gamma,-beta,-alpha] ;
if (isempty(d__));
W_ = wignerd_b(l_max_max,-beta);
C_ = zeros(n_m_max,n_m_max);
n_W_ = zeros(1,1+l_max_max); for (l_val=0:l_max_max); n_W_(1+l_val) = numel(W_{1+l_val}); end;
d__ = zeros(n_m_max,n_m_max,1+l_max_max);
for l_val=0:l_max_max;
d__(1+l_max_max + [-l_val:+l_val],1+l_max_max + [-l_val:+l_val],1+l_val) = W_{1+l_val};
end;%for l_val=0:l_max_max;
end%if (isempty(d__));
if ((isempty(a___) | isempty(b___)) & isempty(ab__));
a___ = zeros(n_m_max,n_k_p_r,1+l_max_max); b___ = zeros(n_m_max,n_k_p_r,1+l_max_max);
ixk=0;
for nk=0:n_k_p_r-1;
ixl=0;
for l_val=0:l_max_(1+nk);
%ixl = l_val*(l_val+1);
%a___(1+l_max_max + [-l_val:+l_val],1+nk,1+l_val) = a_(1+ixk+ixl+-[-l_val:+l_val]);
%b___(1+l_max_max + [-l_val:+l_val],1+nk,1+l_val) = b_(1+ixk+ixl+-[-l_val:+l_val]);
a___(1+l_max_max + [-l_val:+l_val],1+nk,1+l_val) = a_(1+ixk+ixl+l_val+-[-l_val:+l_val]);
b___(1+l_max_max + [-l_val:+l_val],1+nk,1+l_val) = b_(1+ixk+ixl+l_val+-[-l_val:+l_val]);
ixl = ixl+(2*l_val+1);
end;%for l_val=0:l_max_(1+nk);
assert(ixl==n_lm_(1+nk));
ixk = ixk+n_lm_(1+nk);
end;%for nk=0:n_k_p_r-1;
assert(ixk==sum(n_lm_));
end;%if ((isempty(a___) | isempty(b___)) & isempty(ab__));
if (isempty(ab__));
t_0in = tic;
ab__ = zeros(n_m_max,n_m_max,1+l_max_max);
for l_val=0:l_max_max;
ab__(:,:,1+l_val) = conj(squeeze(a___(:,:,1+l_val)))*diag(weight_3d_k_p_r_)*transpose(b___(:,:,1+l_val));
end;%for l_val=0:l_max_max;
t_out = toc(t_0in);
if (verbose>1); disp(sprintf(' %% sum over k: t %0.6f',t_out)); end;
end;%if (isempty(ab__));
t_0in = tic;
C_ = zeros(n_m_max,n_m_max);
C_ = sum(ab__.*d__,3);
t_out = toc(t_0in);
if (verbose>1); disp(sprintf(' %% sum over l: t %0.6f',t_out)); end;
n_op = (n_m_max*n_m_max)*n_k_p_r*(1+l_max_max) + (n_m_max*n_m_max)*(1+l_max_max);
n_mult = (n_m_max*n_m_max)*n_k_p_r*(1+l_max_max) + (n_m_max*n_m_max)*(1+l_max_max);

if (fft_flag==0);
tmp_C = recenter2(squeeze(C_(:,:)));
X0_(:,:,1) = fft2(tmp_C);
end;%if (fft_flag==0);
if (fft_flag==1);
[gamma__,alpha__] = meshgrid(gamma_,alpha_);
tmp_C = squeeze(C_(:,:));
if (verbose>1); disp(sprintf(' %% calling xxnufft2d2')); end;
tmp = xxnufft2d2(n_alpha*n_gamma,alpha__(:),gamma__(:),-1,1e-12,n_m_max,n_m_max,tmp_C);
X0_(:,:,1) = reshape(tmp,n_alpha,n_gamma);
end;%if (fft_flag==1);

