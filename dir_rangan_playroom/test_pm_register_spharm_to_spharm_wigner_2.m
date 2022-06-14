clear;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
flag_center = 1;
flag_invert = 0;
tolerance_master = 1e-2;
nf=0;

dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

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
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[b_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
a_k_Y_form_ = zeros(n_lm_sum,1);
b_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c_,l_max_);
b_k_Y_form_ = b_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_b_c_,l_max_);
end;%for nsource=0:n_source-1;
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_form_ - a_k_Y_quad_)/fnorm(a_k_Y_form_)));
disp(sprintf(' %% b_k_Y_form_ vs b_k_Y_quad_: %0.16f',fnorm(b_k_Y_form_ - b_k_Y_quad_)/fnorm(b_k_Y_form_)));
%%%%%%%%;
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
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
b_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = b_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
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

%%%%%%%%;
% Now register a against b, testing fft. ;
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
%%%%%%%%
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;figbeach();
Xlim_ = [min(abs(X_quad_)),max(abs(X_quad_))];
subplot(1,1,1); imagesc(real(X_quad_));
xlabel('gamma'); set(gca,'XTick',1:n_m_max,'XTickLabel',gamma_);
ylabel('alpha'); set(gca,'YTick',1:n_m_max,'YTickLabel',alpha_);
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
flag_plot=0;
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

%%%%%%%%;
% Now calculate principal-modes for a_k_Y_form_. ;
%%%%%%%%;
a_k_Y_form__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_form__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_form_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
b_k_Y_form__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_k_Y_form__(1:n_lm_(1+nk_p_r),1+nk_p_r) = b_k_Y_form_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
CTF_ones_k_p_r_ = ones(n_k_p_r,1);
CTF_ones_k_p_r__ = ones(n_k_p_r,n_k_p_r);
[ ...
 X_3d_ones_d0__ ...
,X_3d_ones_d0_weight_r_ ...
] = ...
principled_marching_cost_matrix_3( ...
 n_k_p_r ...
,weight_3d_k_p_r_ ...
,l_max_max ...
,a_k_Y_form__ ...
,CTF_ones_k_p_r__ ...
); 
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX_3d_ones_d0__,SX_3d_ones_d0__,VX_3d_ones_d0__] = svds(X_3d_ones_d0__,n_UX_rank); SX_3d_ones_d0_ = diag(SX_3d_ones_d0__);
X__ = X_3d_ones_d0__ ;
X_weight_r_ = X_3d_ones_d0_weight_r_ ;
UX__ = UX_3d_ones_d0__ ;
SX__ = diag(SX_3d_ones_d0_) ;
SX_ = SX_3d_ones_d0_ ;
VX__ = VX_3d_ones_d0__ ;
%%%%%%%%;
% Now transform into principal-volumes. ;
%%%%%%%%;
a_CTF_ones_UX_Y_form__ = zeros(n_lm_max,n_UX_rank);
b_CTF_ones_UX_Y_form__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
assert(tmp_n_lm==n_lm_(1+nk_p_r));
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_ones_UX_Y_form__(1:tmp_n_lm,1+nUX_rank) = a_CTF_ones_UX_Y_form__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_form_(1+tmp_index_)*CTF_ones_k_p_r_(1+nk_p_r);
b_CTF_ones_UX_Y_form__(1:tmp_n_lm,1+nUX_rank) = b_CTF_ones_UX_Y_form__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*b_k_Y_form_(1+tmp_index_)*CTF_ones_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
% Now register principal-volumes against one another. ;
%%%%%%%%;
beta = +pi/3;
n_m_max = 1+2*l_max_max;
%%%%%%%%;
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
%%%%%%%%;
pm_n_k_p_r = n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2;
pm_n_lm_sum = sum(pm_n_lm_);
pm_n_lm_csum_ = cumsum([0;pm_n_lm_]);
[ ...
 pm_X_quad_ ...
,pm_n_op ...
,pm_n_mult ...
] = ...
register_spharm_to_spharm_single_beta_3( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_3d_k_p_r_ ...
,pm_l_max_ ...
,reshape(a_CTF_ones_UX_Y_form__(:,1:pm_n_k_p_r),[pm_n_lm_sum,1]) ...
,reshape(b_CTF_ones_UX_Y_form__(:,1:pm_n_k_p_r),[pm_n_lm_sum,1]) ...
,beta ...
);
%%%%%%%%;
disp(sprintf(' %% X_quad_ vs pm_X_quad_: %0.16f',fnorm(X_quad_ - pm_X_quad_)/fnorm(X_quad_)));
%%%%%%%%;

%%%%%%%%;
% Now calculate accuracy as a function of tmp_pm_n_UX_rank. ;
%%%%%%%%;
pm_E_ = zeros(n_UX_rank,1);
pm_n_op_ = zeros(n_UX_rank,1);
for tmp_pm_n_UX_rank=0:n_UX_rank-1;
if (mod(tmp_pm_n_UX_rank,8)==0); disp(sprintf(' %% tmp_pm_n_UX_rank %d/%d',tmp_pm_n_UX_rank,n_UX_rank)); end;
pm_n_k_p_r = 1+tmp_pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2;
pm_n_lm_sum = sum(pm_n_lm_);
pm_n_lm_csum_ = cumsum([0;pm_n_lm_]);
[ ...
 pm_X_quad_ ...
,pm_n_op ...
,pm_n_mult ...
] = ...
register_spharm_to_spharm_single_beta_3( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_3d_k_p_r_ ...
,pm_l_max_ ...
,reshape(a_CTF_ones_UX_Y_form__(:,1:pm_n_k_p_r),[pm_n_lm_sum,1]) ...
,reshape(b_CTF_ones_UX_Y_form__(:,1:pm_n_k_p_r),[pm_n_lm_sum,1]) ...
,beta ...
);
pm_E = fnorm(X_quad_ - pm_X_quad_)/fnorm(X_quad_);
pm_E_(1+tmp_pm_n_UX_rank) = pm_E;
pm_n_op_(1+tmp_pm_n_UX_rank) = pm_n_op;
end;%for tmp_pm_n_UX_rank=0:n_UX_rank-1;


