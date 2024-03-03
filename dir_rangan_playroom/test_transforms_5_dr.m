%%%%%%%%;
% Sets up a simple volume (in spherical-harmonic-coordinates), ;
% then tests pm_template_2.m. ;
%%%%%%%%;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_verbose = 1;
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
rng(1+nsource);
delta_a_c_ = 2*randn(3,1)/k_p_r_max;
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = 2*randn(3,1)/k_p_r_max;
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
disp(sprintf(' %% I_a_form vs I_a_quad %0.16f %%<-- should be <1e-6',fnorm(I_a_form-I_a_quad)/fnorm(I_a_form)));
disp(sprintf(' %% I_b_form vs I_b_quad %0.16f %%<-- should be <1e-6',fnorm(I_b_form-I_b_quad)/fnorm(I_b_form)));

%%%%%%%%;
% Now set up and test polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = 2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
tmp_S_delta_x_c_ = [cos(pi/4);sin(pi/4)]/k_p_r_max;
tmp_T_delta_x_c_ = [cos(pi/3);sin(pi/3)]/k_p_r_max;
S_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_c_1_wk_*tmp_S_delta_x_c_(1+1)));
T_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_T_delta_x_c_(1+0) + k_c_1_wk_*tmp_T_delta_x_c_(1+1)));
S_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
T_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
if flag_plot;
figure(1+nf);nf=nf+1;figbig; np=0;
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_)); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_c_xx_)); axis image; axisnotick; title('real(S_x_c_xx_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_x_c_xx_)); axis image; axisnotick; title('real(T_x_c_xx_)','Interpreter','none');
end;%if flag_plot;
I_quad = sum(conj(S_k_p_wk_).*T_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
tmp_kd = 2*pi*k_p_r_max*fnorm(tmp_S_delta_x_c_ - tmp_T_delta_x_c_);
I_form = h2d_(tmp_kd)/(2*pi)^2 * (pi*k_p_r_max^2);
disp(sprintf(' %% I_form vs I_quad %0.16f %%<-- should be <1e-6',fnorm(I_form-I_quad)/fnorm(I_form)));
%%%%%%%%;
S_k_p_wk_ = zeros(n_w_sum,1);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
S_k_p_wk_ = S_k_p_wk_ + exp(+2*pi*i*(k_c_0_wk_*delta_a_c__(1+0,1+nsource) + k_c_1_wk_*delta_a_c__(1+1,1+nsource)));
T_k_p_wk_ = T_k_p_wk_ + exp(+2*pi*i*(k_c_0_wk_*delta_b_c__(1+0,1+nsource) + k_c_1_wk_*delta_b_c__(1+1,1+nsource)));
end;%for nsource=0:n_source-1;
S_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
T_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
if flag_plot;
figure(1+nf);nf=nf+1;figbig; np=0;
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_)); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_c_xx_)); axis image; axisnotick; title('real(S_x_c_xx_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_x_c_xx_)); axis image; axisnotick; title('real(T_x_c_xx_)','Interpreter','none');
end;%if flag_plot;
I_quad = sum(conj(S_k_p_wk_).*T_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
I_form = 0;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c__(1:2,1+nsource0) - delta_b_c__(1:2,1+nsource1));
I_form = I_form + h2d_(tmp_kd)/(2*pi)^2 * (pi*k_p_r_max^2);
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
disp(sprintf(' %% I_form vs I_quad %0.16f %%<-- should be <1e-2',fnorm(I_form-I_quad)/fnorm(I_form)));
%%%%%%%%;

%%%%%%%%;
% Now set up spherical-harmonics. ;
%%%%%%%%;
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
a_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c__(:,1+nsource),l_max_);
end;%for nsource=0:n_source-1;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_k_p_to_spharm_4( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_p_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f %%<-- should be <1e-2',fnorm(a_k_Y_form_-a_k_Y_quad_)/fnorm(a_k_Y_form_)));
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_quad_ ...
] = ...
convert_spharm_to_k_p_4( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_p_form_ vs a_k_p_quad_: %0.16f %%<-- should be <1e-2',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_reco_ ...
] = ...
convert_spharm_to_k_p_4( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_p_form_ vs a_k_p_reco_: %0.16f %%<-- should be <1e-2',fnorm(a_k_p_form_-a_k_p_reco_)/fnorm(a_k_p_form_)));
disp(sprintf(' %% a_k_p_quad_ vs a_k_p_reco_: %0.16f %%<-- should be <1e-2',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_)));
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_Y_reco_ ...
] = ...
convert_k_p_to_spharm_4( ...
 verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_reco_: %0.16f %%<-- should be <1e-2',fnorm(a_k_Y_form_-a_k_Y_reco_)/fnorm(a_k_Y_form_)));
disp(sprintf(' %% a_k_Y_quad_ vs a_k_Y_reco_: %0.16f %%<-- should be <1e-2',fnorm(a_k_Y_quad_-a_k_Y_reco_)/fnorm(a_k_Y_quad_)));
%%%%%%%%;
a_k_Y_form_lmk__ = zeros(n_lm_max,n_k_p_r);
a_k_Y_form_lmk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
tmp_a_k_Y_form_lm_ = a_k_Y_form_(1+tmp_index_);
a_k_Y_form_lmk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = tmp_a_k_Y_form_lm_;
tmp_a_k_Y_form_lm__ = zeros(1+l_max_max,n_m_max);
l_max = l_max_(1+nk_p_r);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(na==l_val*(l_val+1)+m_val);
tmp_a_k_Y_form_lm__(1+l_val,1+l_max_max+m_val) = tmp_a_k_Y_form_lm_(1+na);
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
a_k_Y_form_lmk___(:,:,1+nk_p_r) = tmp_a_k_Y_form_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
flag_check=1;
if flag_check;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max=l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_lmk__(1+l_val*(l_val+1)+m_val,1+nk_p_r));
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_lmk___(1+l_val,1+l_max_max+m_val,1+nk_p_r));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
end;%if flag_check;
%%%%%%%%;
a_k_Y_form_lkm___ = permute(a_k_Y_form_lmk___,[1,3,2]);
%%%%%%%%;

%%%%%%%%;
% define rotations. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;

%%%%%%%%;
% Now generate templates from the volumes. ;
%%%%%%%%;
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_form_lmk__,[n_lm_max,n_k_p_r]) ...
,1.0/k_p_r_max ...
,-1 ...
,n_w_max ...
);
n_S = n_viewing_S;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
%%%%%%%%;
% Now step through and reconstitute the templates. ;
%%%%%%%%;
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
if flag_plot;
clf;figmed;
subplot(1,2,1);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_)); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(1,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
end;%if flag_plot;
T_k_p_wkS__(:,1+nS) = T_k_p_wk_;
end;%for nS=0:n_S-1;
disp(sprintf(' %% T_k_p_wkS__ vs S_k_p_wkS__: %0.16f %%<-- should be <1e-2',fnorm(T_k_p_wkS__-S_k_p_wkS__)/fnorm(T_k_p_wkS__)));
if flag_plot;
clf;figsml; plot(0:n_S-1,sum(abs(T_k_p_wkS__ - S_k_p_wkS__).^2,1)./sum(abs(T_k_p_wkS__).^2,1),'o');
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now test out direct spherical-harmonic contruction of S_k_q_wk_. ;
% Note that below we introduce a (unnecessary) tmp_gamma_z to check consistency. ;
%%%%%%%%;
%nS = 185;
nS = 238;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = pi/12;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
disp(sprintf(' %% S_k_p_wk_ vs T_k_p_wk_: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_wk_-T_k_p_wk_)/fnorm(S_k_p_wk_)));
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
disp(sprintf(' %% S_k_q_wk_ vs T_k_q_wk_: %0.16f %%<-- should be <1e-2',fnorm(S_k_q_wk_-T_k_q_wk_)/fnorm(S_k_q_wk_)));
%%%%;
tmp_euler_ = [-tmp_azimu_b,-tmp_polar_a,-tmp_gamma_z];
W_beta__ = wignerd_b(l_max_max,+tmp_euler_(1+1));
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
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
,tmp_euler_ ...
);
%%%%;
Q_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm_(1+nk_p_r)-1];
a_k_Y_sub_ = a_k_Y_rota_(1+tmp_index_);
l_max = l_max_(1+nk_p_r);
for nw=0:n_w_max-1;
nq = periodize(nw,-n_w_max/2,+n_w_max/2); m0_val = nq;
for l_val=0:l_max;
if abs(m0_val)<=l_val;
Q_k_q_wk_(1+nw+nk_p_r*n_w_max) = Q_k_q_wk_(1+nw+nk_p_r*n_w_max) + ...
 (pi^2) ...
*a_k_Y_sub_(1+(l_val^2 + l_val + m0_val)) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
;
end;%if abs(m0_val)<=l_val;
end;%for l_val=0:l_max;
end;%for nw=0:n_w_max-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm_(1+nk_p_r)-1];
a_k_Y_sub_ = a_k_Y_form_(1+tmp_index_);
l_max = l_max_(1+nk_p_r);
for nw=0:n_w_max-1;
nq = periodize(nw,-n_w_max/2,+n_w_max/2); m0_val = nq;
for l_val=0:l_max;
for m1_val=-l_val:+l_val;
if abs(m0_val)<=l_val;
R_k_q_wk_(1+nw+nk_p_r*n_w_max) = R_k_q_wk_(1+nw+nk_p_r*n_w_max) + ...
 (pi^2) ...
*exp(-i*m0_val*(+tmp_euler_(1+2))) ...
*W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*exp(-i*m1_val*(+tmp_euler_(1+0))) ...
*a_k_Y_sub_(1+(l_val^2 + l_val + m1_val)) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
;
end;%if abs(m0_val)<=l_val;
end;%for m1_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nw=0:n_w_max-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
disp(sprintf(' %% R_k_q_wk_ vs Q_k_q_wk_: %0.16f %%<-- should be <1e-6',fnorm(R_k_q_wk_-Q_k_q_wk_)/fnorm(R_k_q_wk_)));
disp(sprintf(' %% S_k_q_wk_ vs Q_k_q_wk_: %0.16f %%<-- should be <1e-2',fnorm(S_k_q_wk_-Q_k_q_wk_)/fnorm(S_k_q_wk_)));
disp(sprintf(' %% S_k_q_wk_ vs R_k_q_wk_: %0.16f %%<-- should be <1e-2',fnorm(S_k_q_wk_-R_k_q_wk_)/fnorm(S_k_q_wk_)));
disp(sprintf(' %% T_k_q_wk_ vs Q_k_q_wk_: %0.16f %%<-- should be <1e-2',fnorm(T_k_q_wk_-Q_k_q_wk_)/fnorm(T_k_q_wk_)));
disp(sprintf(' %% T_k_q_wk_ vs R_k_q_wk_: %0.16f %%<-- should be <1e-2',fnorm(T_k_q_wk_-R_k_q_wk_)/fnorm(T_k_q_wk_)));
%%%%%%%%;

%%%%%%%%;
% Now we try and set up a template-operator ;
% for a collection of azimu_b associated with a single polar_a. ;
%%%%%%%%;
polar_a_use = viewing_polar_a_S_(round(n_S/4));
tmp_index_ = efind(abs(viewing_polar_a_S_-polar_a_use)<1e-6);
n_azimu_b_use = numel(tmp_index_);
azimu_b_use_ = viewing_azimu_b_S_(1+tmp_index_);
S_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
S_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
nS = tmp_index_(1+nazimu_b_use);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = -pi/9;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
S_k_p_sub_wkb__(:,1+nazimu_b_use) = S_k_p_wk_;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_sub_wkb__(:,1+nazimu_b_use) = T_k_p_wk_;
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
S_k_q_sub_wkb__(:,1+nazimu_b_use) = S_k_q_wk_;
T_k_q_sub_wkb__(:,1+nazimu_b_use) = T_k_q_wk_;
clear S_k_p_wk_ T_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
disp(sprintf(' %% S_k_p_sub_wkb__ vs T_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__-T_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% S_k_q_sub_wkb__ vs T_k_q_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_q_sub_wkb__-T_k_q_sub_wkb__)/fnorm(S_k_q_sub_wkb__)));
%%%%;
flag_plot=0;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
%title(sprintf('real(S_k_p_sub_wkb__(:,1+%d))',nazimu_b_use),'Interpreter','none');
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
end;%if flag_plot;
%%%%;
tmp_t = tic();
W_beta__ = wignerd_b(l_max_max,-polar_a_use);
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_betazeta_mlm___: %0.2fs',tmp_t)); end;
%%%%;
flag_check=1;
if flag_check;
tmp_w_ = crandn(n_m_max);
tmp_azimu_b_use_ = 2*pi*rand(n_azimu_b_use,1);
tmp_f__ = exp(-i*(-reshape(tmp_azimu_b_use_,[n_azimu_b_use,1]))*reshape(m_max_,[1,n_m_max]));
tmp_fw_0_ = tmp_f__*tmp_w_;
tmp_fw_1_ = xxnufft1d2(n_azimu_b_use,tmp_azimu_b_use_,+1,1e-6,n_m_max,tmp_w_);
disp(sprintf(' %% tmp_fw_0_ vs tmp_fw_1_: %0.16f %%<-- should be <1e-6',fnorm(tmp_fw_0_-tmp_fw_1_)/fnorm(tmp_fw_0_)));
clear tmp_azimu_b_use_ tmp_w_ tmp_f__ tmp_fw_0_ tmp_fw_1_ ;
end;%if flag_check;
%%%%;
tmp_t = tic();
W_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag(exp(+i*m0_val_*tmp_gamma_z))*W_betazeta_ml__*a_k_Y_form_lk__ for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
W_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 (pi^2) ...
 *diag(exp(-i*m_max_*(-tmp_gamma_z))) ...
 *reshape(W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_form_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_caza_mkm___: %0.2fs',tmp_t)); end;
tmp_t = tic();
W_caza_mmk___ = permute(W_caza_mkm___,[3,1,2]);
W_caza_bmk___ = reshape(xxnufft1d2(n_azimu_b_use,azimu_b_use_,+1,1e-6,n_m_max,reshape(W_caza_mmk___,[n_m_max,n_m_max*n_k_p_r])),[n_azimu_b_use,n_m_max,n_k_p_r]);
W_caza_mkb___ = permute(W_caza_bmk___,[2,3,1]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% W_caza_mkb___: %0.2fs',tmp_t)); end;
%%%%%%%%;
R_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
R_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
R_k_p_wk_ = zeros(n_w_sum,1);
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq = m_val; if (nq<0); nq=nq+n_w_max; end;
R_k_q_wk_(1+nq+nk_p_r*n_w_max) = W_caza_mkb___(1+l_max_max+m_val,1+nk_p_r,1+nazimu_b_use);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
R_k_q_sub_wkb__(:,1+nazimu_b_use) = R_k_q_wk_;
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
R_k_p_sub_wkb__(:,1+nazimu_b_use) = R_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
disp(sprintf(' %% S_k_p_sub_wkb__ vs T_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__ - T_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% S_k_p_sub_wkb__ vs R_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__ - R_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% T_k_p_sub_wkb__ vs R_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(T_k_p_sub_wkb__ - R_k_p_sub_wkb__)/fnorm(T_k_p_sub_wkb__)));
%%%%;
flag_plot=0;
if flag_plot;
%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
%%;
end;%if flag_plot;
%%%%%%%%;
[ ...
 U_k_p_sub_wkb__ ...
,W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_form_lkm___ ...
,n_w_max ...
,polar_a_use ...
,n_azimu_b_use ...
,azimu_b_use_ ...
,tmp_gamma_z ...
,W_betazeta_mlm___ ...
);
disp(sprintf(' %% S_k_p_sub_wkb__ vs U_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_sub_wkb__ - U_k_p_sub_wkb__)/fnorm(S_k_p_sub_wkb__)));
disp(sprintf(' %% T_k_p_sub_wkb__ vs U_k_p_sub_wkb__: %0.16f %%<-- should be <1e-2',fnorm(T_k_p_sub_wkb__ - U_k_p_sub_wkb__)/fnorm(T_k_p_sub_wkb__)));
disp(sprintf(' %% R_k_p_sub_wkb__ vs U_k_p_sub_wkb__: %0.16f %%<-- should be <1e-9',fnorm(R_k_p_sub_wkb__ - U_k_p_sub_wkb__)/fnorm(R_k_p_sub_wkb__)));
%%%%%%%%;

%%%%%%%%;
% We check the following sum: ;
% \sum_{nq} zeta_lm__(1+l0,1+nq) * zeta_lm__(1+l1,1+nq) * W_beta__{1+l0}(1+nq,1+m0) * W_beta__{1+l1}(1+nq,1+m1) ;
% Verdict: not too sparse. ;
%%%%%%%%;
flag_check=0;
if flag_check;
sum0_lmlm____ = zeros(1+l_max_max,n_m_max,1+l_max_max,n_m_max);
for l0_val=0:l_max_max;
disp(sprintf(' %% slow calculation: l0_val %d/%d',l0_val,l_max_max));
for l1_val=0:l_max_max;
lm_val = min(l0_val,l1_val);
for m0_val=-l0_val:+l0_val;
for m1_val=-l1_val:+l1_val;
for q_val=-lm_val:+lm_val;
sum0_lmlm____(1+l0_val,1+l_max_max+m0_val,1+l1_val,1+l_max_max+m1_val) = ...
 sum0_lmlm____(1+l0_val,1+l_max_max+m0_val,1+l1_val,1+l_max_max+m1_val) ...
+zeta_lm__(1+l0_val,1+l_max_max+q_val) ...
*zeta_lm__(1+l1_val,1+l_max_max+q_val) ...
*W_beta__{1+l0_val}(1+l0_val+q_val,1+l0_val+m0_val) ...
*W_beta__{1+l1_val}(1+l1_val+q_val,1+l1_val+m1_val) ...
;
end;%for q_val=-lm_val:+lm_val;
end;%for m1_val=-l1_val:+l1_val;
end;%for m0_val=-l0_val:+l0_val;
end;%for l1_val=0:l_max_max;
end;%for l0_val=0:l_max_max;
%%%%;
%figure(1);clf;figsml;imagesc(reshape(log10(abs(sum0_lmlm____)),[(1+l_max_max)*n_m_max,(1+l_max_max)*n_m_max]),[-5,0]);colorbar; axis image;
%figure(2);clf;figsml;imagesc(reshape(log10(abs(permute(sum0_lmlm____,[1,3,2,4]))),[(1+l_max_max)^2,n_m_max^2]),[-5,0]);colorbar; axis image;
%%%%;
disp(sprintf(' %% nnz sum0_lmlm____: %d / %d --> %0.2f',nnz(sum0_lmlm____),numel(sum0_lmlm____),nnz(sum0_lmlm____)/numel(sum0_lmlm____)));
%%%%;
sum1_lmlm____ = zeros(1+l_max_max,n_m_max,1+l_max_max,n_m_max);
tmp_mlm__ = reshape(W_betazeta_mlm___,[n_m_max,(1+l_max_max)*n_m_max]);
sum1_lmlm____ = reshape(transpose(tmp_mlm__)*tmp_mlm__,[1+l_max_max,n_m_max,1+l_max_max,n_m_max]);
disp(sprintf(' %% nnz sum1_lmlm____: %d / %d --> %0.2f',nnz(sum1_lmlm____),numel(sum1_lmlm____),nnz(sum1_lmlm____)/numel(sum1_lmlm____)));
disp(sprintf(' %% sum0_lmlm____ vs sum1_lmlm____: %0.16f %%<-- should be <1e-2',fnorm(sum0_lmlm____ - sum1_lmlm____)/fnorm(sum0_lmlm____)));
%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% integral of plane-wave times bessel times plane-wave. ;
% (i.e., isotropic CTF). ;
%%%%%%%%;
% C = J_{0}(\alpha k) ;
% S = exp(i 2\pi \vk \dot \vd_{S}) = exp(i 2\pi k \delta_{S} cos(\psi - \omega_{S})) ;
% M = exp(i 2\pi \vk \dot \vd_{M}) = exp(i 2\pi k \delta_{M} cos(\psi - \omega_{M})) ;
% \vd_{T} = \vd_{M} - \vd_{S} = \delta_{T} [ \cos(\omega_{T}) ; \sin(\omega_{T}) ] ;
% a = \alpha K ;
% b = 2\pi K \delta_{T} ;
% \int_{0}^{K}\int_{0}^{2\pi} [C\cdot S]^{\dagger} M d\psi kdk ... ;
%                            = 2\pi K^2 \frac{ bJ_{-1}(b)J_{0}(a) - aJ_{-1}(a)J_{0}(b) }{ a^2 - b^2 } ;
%%%%%%%%;
k_p_r_max = 48/(2*pi);
delta_max = 0.1;
tmp_delta_S = rand()*delta_max;
tmp_omega_S = 2*pi*rand();
tmp_delta_S_ = tmp_delta_S * [cos(tmp_omega_S) ; sin(tmp_omega_S)];
tmp_delta_M = rand()*delta_max;
tmp_omega_M = 2*pi*rand();
tmp_delta_M_ = tmp_delta_M * [cos(tmp_omega_M) ; sin(tmp_omega_M)];
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_alpha = rand();
tmp_a = tmp_alpha*k_p_r_max;
tmp_b = 2*pi*k_p_r_max*tmp_delta_T;
tmp_I = integral2( ...
@(k_p_r,psi) k_p_r ...
 .* conj(besselj(0,tmp_alpha*k_p_r).*exp(2*pi*i*k_p_r.*tmp_delta_S.*cos(psi-tmp_omega_S))) ...
 .* exp(2*pi*i*k_p_r.*tmp_delta_M.*cos(psi-tmp_omega_M)) ...
,0,k_p_r_max,0,2*pi);
tmp_c = max(tmp_a,tmp_b); tmp_d = min(tmp_a,tmp_b);
tmp_J = 2*pi*k_p_r_max^2 * (tmp_d*besselj(-1,tmp_d)*besselj(0,tmp_c) - tmp_c*besselj(-1,tmp_c)*besselj(0,tmp_d))/max(1e-12,tmp_c^2-tmp_d^2);
if (flag_verbose>0); disp(sprintf(' %% tmp_I vs tmp_J: %0.16f %%<-- should be <1e-6 ',fnorm(tmp_I-tmp_J)/fnorm(tmp_I))); end;
%%%%%%%%;

%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
S_k_q_wkS__(:,1+nS) = S_k_q_wk_;
end;%for nS=0:n_S-1;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
%%%%%%%%;

%%%%%%%%;
% Now we generate images for the likelihood-calculation. ;
%%%%%%%%;
n_CTF = 3;
%CTF_alpha_C_ = transpose(linspace(0.5,1.0,n_CTF)); %<-- changes sign. ;
CTF_alpha_C_ = transpose(linspace(0.05,2.4048/k_p_r_max,n_CTF)); %<-- stays positive, since first root of besselj(0,.) is ~2.4048 ;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
CTF_k_p_r_kC__(:,1+nCTF) = besselj(0,CTF_alpha_C_(1+nCTF)*k_p_r_);
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
pm_n_UX_rank = n_k_p_r; %<-- keep the same for now. ;
UX_knC___ = zeros(n_k_p_r,pm_n_UX_rank,n_CTF);
X_weight_rC__ = zeros(n_k_p_r,n_CTF);
UX_kn__ = eye(n_k_p_r);
X_weight_r_ = sqrt(weight_2d_k_p_r_);
for nCTF=0:n_CTF-1;
UX_knC___(:,:,1+nCTF) = UX_kn__;
X_weight_rC__(:,1+nCTF) = X_weight_r_;
end;%for nCTF=0:n_CTF-1;
clear nCTF UX_kn__ X_weight_r_ ;
%%%%%%%%;
tmp_t = tic();
n_M = 1024;
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
index_nCTF_from_nM_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
index_nw_from_nM_ = zeros(n_M,1);
euler_azimu_b_true_M_ = zeros(n_M,1);
euler_polar_a_true_M_ = zeros(n_M,1);
euler_gamma_z_true_M_ = zeros(n_M,1);
rng(0);
for nM=0:n_M-1;
nCTF = mod(nM,n_CTF);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
index_nCTF_from_nM_(1+nM) = nCTF;
nS = max(0,min(n_S-1,floor(n_S*rand())));
index_nS_from_nM_(1+nM) = nS;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
azimu_b = viewing_azimu_b_S_(1+nS);
polar_a = viewing_polar_a_S_(1+nS);
index_gamma_z = periodize(nM,0,n_gamma_z);
index_nw_from_nM_(1+nM) = periodize(-index_gamma_z,0,n_w_max);
gamma_z = gamma_z_(1+index_gamma_z);
euler_azimu_b_true_M_(1+nM) = azimu_b;
euler_polar_a_true_M_(1+nM) = polar_a;
euler_gamma_z_true_M_(1+nM) = gamma_z;
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-gamma_z); %<-- note here we rotate S_k_p_wk_ by -gamma_z to form M_k_p_wk_. ;
%M_k_p_wk_ = reshape(bsxfun(@times,reshape(T_k_p_wk_,[n_w_max,n_k_p_r]),reshape(CTF_k_p_r_k_,[1,n_k_p_r])),[n_w_sum,1]); %<-- do not multiply the images by CTF. ;
M_k_p_wk_ = T_k_p_wk_;
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_wkM__(:,1+nM) = M_k_q_wk_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_p_wkM__: %0.2fs',tmp_t)); end;
clear nCTF CTF_k_p_r_k_ nS azimu_b polar_a index_gamma_z gamma_z S_k_p_wk_ T_k_p_wk_ M_k_p_wk_ M_k_q_wk_ ;
%%%%%%%%;

%%%%%%%%;
delta_r_max = 0; svd_eps = 1e-2; n_delta_v_requested = 0;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
%%%%%%%%;

%%%%%%%%;
% Now calculate correlation X_wSM___ and innerproduct Z_wSM___
%%%%%%%%;
X_wSM___ = zeros(n_w_max,n_S,n_M);
gamma_z_wSM___ = zeros(n_w_max,n_S,n_M);
Z_wSM___ = zeros(n_w_max,n_S,n_M);
UX_M_l2_M_ = zeros(n_M,1);
CTF_UX_S_l2_SC__ = zeros(n_S,n_CTF);
for nCTF=0:n_CTF-1;
if (verbose); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
UX_kn__ = UX_knC___(:,:,1+nCTF);
X_weight_r_ = X_weight_rC__(:,1+nCTF);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
tmp_t = tic();
svd_VUXM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q_wkM__(:,1+index_M_sub_),pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% tmpmh_VUXM_lwnM____3: %0.2fs',tmp_t)); end;
tmp_t = tic();
UX_M_sub_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_sub_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% ampmh_UX_M_sub_l2_dm__1: %0.2fs',tmp_t)); end;
tmp_index_d_ = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); %<-- should be zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_M_sub_l2_dM__(1+tmp_index_d_(1+0),:),[n_M_sub,1]);
[UX_M_sub_k_q_wnM___,UX_M_sub_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M_sub,svd_VUXM_sub_lwnM____);
%%%%;
tmp_t = tic();
CTF_UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% CTF_UX_S_k_q_wnS__: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
parameter_ampmh = struct('type','parameter');
[ ...
 parameter_ampmh ...
,X_sub_wSM___ ...
,~ ...
,~ ...
,gamma_z_sub_wSM___ ...
,~ ...
] = ...
ampmh_X_wSM___8( ...
 parameter_ampmh ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_S_ ...
,n_M_sub ...
,svd_VUXM_sub_lwnM____ ...
,UX_M_sub_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% ampmh_X_wSM___8: %0.2fs',tmp_t)); end;
UX_M_l2_M_(1+index_M_sub_) = UX_M_sub_l2_M_;
CTF_UX_S_l2_SC__(:,1+nCTF) = CTF_UX_S_l2_S_;
Z_sub_wSM___ = bsxfun(@times,bsxfun(@times,X_sub_wSM___,reshape(sqrt(CTF_UX_S_l2_S_),[1,n_S,1])),reshape(sqrt(UX_M_sub_l2_M_),[1,1,n_M_sub]));
X_wSM___(:,:,1+index_M_sub_) = X_sub_wSM___;
gamma_z_wSM___(:,:,1+index_M_sub_) = gamma_z_sub_wSM___;
Z_wSM___(:,:,1+index_M_sub_) = Z_sub_wSM___;
%%%%;
clear UX_kn__ X_weight_r_ CTF_k_p_r_k_ index_M_sub_ svd_VUXM_sub_lwnM____ UX_M_sub_l2_dM__ UX_M_sub_l2_M_ ;
clear UX_M_sub_k_q_wnM___ UX_M_sub_k_p_wnM___ ;
clear CTF_UX_S_k_q_wnS__ CTF_UX_S_l2_S_ ;
clear X_sub_wSM___ gamma_z_sub_wSM___ Z_sub_wSM___ ;
end;%for nCTF=0:n_CTF-1;
clear nCTF;
%%%%%%%%;

%%%%%%%%;
% Now test Z_wSM___ for a few images and templates. ;
%%%%%%%%;
n_inpla_z = n_w_max; inpla_z_ = transpose(linspace(0,2*pi,n_inpla_z+1)); inpla_z_ = inpla_z_(1:n_inpla_z);
flag_plot = 1;
if flag_plot;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 5; p_col = 5; np=0;
end;%if flag_plot;
for nM = [58,128,250,377,800];
nS_from_nM = index_nS_from_nM_(1+nM);
nCTF_from_nM = index_nCTF_from_nM_(1+nM);
CTF_alpha = CTF_alpha_C_(1+nCTF_from_nM);
nw_from_nM = index_nw_from_nM_(1+nM);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF_from_nM);
azimu_b = euler_azimu_b_true_M_(1+nM);
polar_a = euler_polar_a_true_M_(1+nM);
gamma_z = euler_gamma_z_true_M_(1+nM);
tmp_R_M__ = Rz(-gamma_z)*Ry(-polar_a)*Rz(-azimu_b);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
N_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS_from_nM),-gamma_z);
disp(sprintf(' %% %% M_k_p_wk_ vs N_k_p_wk_: %0.16f %%<-- should be <1e-2',fnorm(M_k_p_wk_-N_k_p_wk_)/fnorm(M_k_p_wk_)));
O_k_p_wk_ = zeros(n_w_sum,1);
for nsource_M=0:n_source-1;
tmp_delta_M_ = tmp_R_M__*delta_a_c__(:,1+nsource_M);
O_k_p_wk_ = O_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_M_(1+0) + k_c_1_wk_*tmp_delta_M_(1+1)));
end;%for nsource_M=0:n_source-1;
disp(sprintf(' %% %% M_k_p_wk_ vs O_k_p_wk_: %0.16f %%<-- should be <1e-2',fnorm(M_k_p_wk_-O_k_p_wk_)/fnorm(M_k_p_wk_)));
%%%%%%%%;
for nS = [nS_from_nM,65,82,129,900];
disp(sprintf(' %% nS %.4d nM %.4d',nS,nM));
%%%%;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R_S__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource_S=0:n_source-1;
tmp_delta_S_ = tmp_R_S__*delta_a_c__(:,1+nsource_S);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_S_(1+0) + k_c_1_wk_*tmp_delta_S_(1+1)));
end;%for nsource_S=0:n_source-1;
disp(sprintf(' %% %% S_k_p_wk_ vs T_k_p_wk_: %0.16f %%<-- should be <1e-2',fnorm(S_k_p_wk_-T_k_p_wk_)/fnorm(S_k_p_wk_)));
%%%%;
CTF_S_k_p_wk_ = bsxfun(@times,reshape(S_k_p_wk_,[n_w_max,n_k_p_r]),reshape(CTF_k_p_r_k_,[1,n_k_p_r]));
CTF_UX_S_l2_0_w_ = zeros(n_w_max,1);
CTF_UX_S_l2_1_w_ = zeros(n_w_max,1);
UX_M_l2_0_w_ = zeros(n_w_max,1);
UX_M_l2_1_w_ = zeros(n_w_max,1);
I0_w_ = zeros(n_w_max,1);
I1_w_ = zeros(n_w_max,1);
I2_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
inpla_z = (2*pi*nw)/n_w_max;
R_CTF_S_k_p_wk_ = rotate_p_to_p(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_,+inpla_z);
CTF_UX_S_l2_0_w_(1+nw) = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,R_CTF_S_k_p_wk_,R_CTF_S_k_p_wk_))/(2*pi);
CTF_UX_S_l2_1_w_(1+nw) = CTF_UX_S_l2_SC__(1+nS,1+nCTF_from_nM);
UX_M_l2_0_w_(1+nw) = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_wk_,M_k_p_wk_))/(2*pi);
UX_M_l2_1_w_(1+nw) = UX_M_l2_M_(1+nM);
I0_w_(1+nw) = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,R_CTF_S_k_p_wk_,M_k_p_wk_))/(2*pi);
assert( abs( Z_wSM___(1+nw,1+nS,1+nM) - X_wSM___(1+nw,1+nS,1+nM)*sqrt(UX_M_l2_M_(1+nM))*sqrt(CTF_UX_S_l2_SC__(1+nS,1+nCTF_from_nM)) )<1e-6 );
I1_w_(1+nw) = Z_wSM___(1+nw,1+nS,1+nM);
%%;
tmp_I = 0;
for nsource_S=0:n_source-1;
for nsource_M=0:n_source-1;
tmp_delta_S_ = Rz(+inpla_z)*tmp_R_S__*delta_a_c__(:,1+nsource_S);
tmp_delta_M_ = tmp_R_M__*delta_a_c__(:,1+nsource_M);
tmp_delta_T_ = tmp_delta_M_ - tmp_delta_S_;
tmp_delta_T_ = tmp_delta_T_(1:2);
tmp_delta_T = fnorm(tmp_delta_T_);
tmp_omega_T = atan2(tmp_delta_T_(1+1),tmp_delta_T_(1+0));
tmp_alpha = CTF_alpha_C_(1+nCTF_from_nM);
tmp_a = tmp_alpha*k_p_r_max;
tmp_b = 2*pi*k_p_r_max*tmp_delta_T;
tmp_c = max(tmp_a,tmp_b); tmp_d = min(tmp_a,tmp_b);
tmp_J = 2*pi*k_p_r_max^2 * (tmp_d*besselj(-1,tmp_d)*besselj(0,tmp_c) - tmp_c*besselj(-1,tmp_c)*besselj(0,tmp_d))/max(1e-12,tmp_c^2-tmp_d^2);
tmp_I = tmp_I + tmp_J;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
%%;
I2_w_(1+nw) = tmp_I;
end;%for nw=0:n_w_max-1;
%%%%;
if flag_plot;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
plot(inpla_z_,I0_w_,'k.');
plot(inpla_z_,I1_w_,'mx-','LineWidth',2);
plot(inpla_z_,I2_w_,'co','LineWidth',2);
plot(inpla_z_(1+nw_from_nM),I0_w_(1+nw_from_nM),'kh','MarkerSize',24);
title(sprintf('nM %d nS %d[%d] (K\alpha %0.2f)',nM,nS,nS_from_nM,CTF_alpha*k_p_r_max),'Interpreter','none');
hold off;
xlim([0,2*pi]); xlabel('inpla_z','Interpreter','none');
ylabel('C','Interpreter','none');
end;%if flag_plot;
%%%%;
disp(sprintf(' %% %% CTF_UX_S_l2_0_w_ vs CTF_UX_S_l2_1_w_: %0.16f %%<-- should be <1e-6',fnorm(CTF_UX_S_l2_0_w_-CTF_UX_S_l2_1_w_)/fnorm(CTF_UX_S_l2_0_w_)));
disp(sprintf(' %% %% UX_M_l2_0_w_ vs UX_M_l2_1_w_: %0.16f %%<-- should be <1e-6',fnorm(UX_M_l2_0_w_-UX_M_l2_1_w_)/fnorm(UX_M_l2_0_w_)));
disp(sprintf(' %% %% I0_w_ vs I1_w_: %0.16f %%<-- should be <1e-6',fnorm(I0_w_-I1_w_)/fnorm(I0_w_)));
disp(sprintf(' %% %% I0_w_ vs I2_w_: %0.16f %%<-- should be <1e-2',fnorm(I0_w_-I2_w_)/fnorm(I0_w_)));
disp(sprintf(' %% %% I1_w_ vs I2_w_: %0.16f %%<-- should be <1e-2',fnorm(I1_w_-I2_w_)/fnorm(I1_w_)));
%%%%%%%%;
clear S_k_p_wk_ tmp_azimu_b tmp_polar_a tmp_gamma_z tmp_R_S__ T_k_p_wk_ tmp_delta_S_ T_k_p_wk_ CTF_S_k_p_wk_ CTF_UX_S_l2_0_w_ CTF_UX_S_l2_1_w_ UX_M_l2_0_w_ UX_M_l2_1_w_ I0_w_ I1_w_ I2_w_ ;
end;%for nS = [nS_from_nM,65,82,129,900];
clear nS_from_nM nCTF_from_nM CTF_alpha nw_from_nM CTF_k_p_r_k_ azimu_b polar_a gamma_z tmp_R_M__ M_k_p_wk_ N_k_p_wk_ O_k_p_wk_ tmp_delta_M_ O_k_p_wk_ ;
end;%for nM = [58,128,250,377,800];
%%%%%%%%;

%%%%%%%%;
% Now calculate the likelihood directly from the images and the templates. ;
%%%%%%%%;
ssnll_0_M_ = zeros(n_M,1);
for nM=0:n_M-1;
nCTF_from_nM = index_nCTF_from_nM_(1+nM);
nS_from_nM = index_nS_from_nM_(1+nM);
nw_from_nM = index_nw_from_nM_(1+nM);
ssnll_0_M_(1+nM) = 0.5 * (CTF_UX_S_l2_SC__(1+nS_from_nM,1+nCTF_from_nM) + UX_M_l2_M_(1+nM) - 2*Z_wSM___(1+nw_from_nM,1+nS_from_nM,1+nM));
clear nCTF_from_nM nS_from_nM nw_from_nM ;
end;%for nM=0:n_M-1;
ssnll_0 = sum(ssnll_0_M_);

%%%%%%%%;
% Redefine zeta_lm__ (in case it was overwritten). ;
%%%%%%%%;
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%%%%%;

%%%%%%%%;
% Express likelihood using spherical-harmonics and fourier-bessel coordinates. ;
% Test on only a few images. ;
%%%%%%%%;
ssnll_1_M_ = zeros(n_M,1);
%%%%%%%%;
n_M_sub = 32;
for nM=0:n_M_sub-1;%for nM=0:n_M-1;
if (verbose); if (mod(nM,32)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
%%%%%%%%;
nCTF_from_nM = index_nCTF_from_nM_(1+nM);
nS_from_nM = index_nS_from_nM_(1+nM);
nw_from_nM = index_nw_from_nM_(1+nM);
tmp_azimu_b = euler_azimu_b_true_M_(1+nM);
assert(tmp_azimu_b==viewing_azimu_b_S_(1+nS_from_nM));
tmp_polar_a = euler_polar_a_true_M_(1+nM);
assert(tmp_polar_a==viewing_polar_a_S_(1+nS_from_nM));
tmp_gamma_z = euler_gamma_z_true_M_(1+nM);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF_from_nM);
%%%%;
W_beta__ = wignerd_b(l_max_max,-tmp_polar_a);
W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
W_cazba_mk__ = zeros(n_m_max,n_k_p_r); %<-- exp(+i*m0_val_*tmp_gamma_z)*W_betazeta_mlm___*exp(+i*m1_val_*tmp_azimu_b)*a_k_Y_form_mlk___ for each m0_val. ;
for nk_p_r=0:n_k_p_r-1;
for m0_val=-l_max_max:+l_max_max;
W_cazba = ...
 (pi^2) ...
*exp(+i*m0_val*(+tmp_gamma_z)) ...
*sum( ...
      bsxfun(@times ...
	     ,bsxfun(@times ...
		     ,reshape(W_betazeta_mlm___(1+l_max_max+m0_val,:,:),[1+l_max_max,n_m_max]) ...
		     ,reshape(exp(+i*m_max_*(+tmp_azimu_b)),[1,n_m_max]) ...
		     ) ...
	     ,reshape(a_k_Y_form_lmk___(:,:,1+nk_p_r),[1+l_max_max,n_m_max]) ...
	     ) ...
      ,'all' ...
      ) ...
  ;
W_cazba_mk__(1+l_max_max+m0_val,1+nk_p_r) = W_cazba;
end;%for m0_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
clear W_beta__ clear W_betazeta_mlm___;
%%%%;
R_k_p_wk_ = zeros(n_w_sum,1);
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq = m_val; if (nq<0); nq=nq+n_w_max; end;
R_k_q_wk_(1+nq+nk_p_r*n_w_max) = W_cazba_mk__(1+l_max_max+m_val,1+nk_p_r);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS_from_nM);
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z); %<-- note here we rotate S_k_p_wk_ by -tmp_gamma_z to form M_k_p_wk_. ;
if (verbose); disp(sprintf(' %% R_k_p_wk_ vs T_k_p_wk_: %0.16f %%<-- should be <1e-2',fnorm(R_k_p_wk_-T_k_p_wk_)/fnorm(R_k_p_wk_))); end;
CTF_R_k_p_wk_ = reshape(bsxfun(@times,reshape(CTF_k_p_r_k_,[1,n_k_p_r]),reshape(R_k_p_wk_,[n_w_max,n_k_p_r])),[n_w_sum,1]);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = M_k_q_wkM__(:,1+nM);
CTF_UX_R_l2_0 = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,CTF_R_k_p_wk_,CTF_R_k_p_wk_))/(2*pi);
CTF_UX_R_l2_1 = CTF_UX_S_l2_SC__(1+nS_from_nM,1+nCTF_from_nM);
if (verbose); disp(sprintf(' %% CTF_UX_R_l2_0 vs CTF_UX_R_l2_1: %0.16f %%<-- should be <1e-2',fnorm(CTF_UX_R_l2_0-CTF_UX_R_l2_1)/fnorm(CTF_UX_R_l2_0))); end;
UX_M_l2_0 = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_wk_,M_k_p_wk_))/(2*pi);
UX_M_l2_1 = UX_M_l2_M_(1+nM);
if (verbose); disp(sprintf(' %% UX_M_l2_0 vs UX_M_l2_1: %0.16f %%<-- should be <1e-2',fnorm(UX_M_l2_0-UX_M_l2_1)/fnorm(UX_M_l2_0))); end;
Z_0 = real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,CTF_R_k_p_wk_,M_k_p_wk_))/(2*pi);
Z_1 = Z_wSM___(1+nw_from_nM,1+nS_from_nM,1+nM);
if (verbose); disp(sprintf(' %% Z_0 vs Z_1: %0.16f %%<-- should be <1e-2',fnorm(Z_0-Z_1)/fnorm(Z_0))); end;
ssnll_1_M_(1+nM) = 0.5 * real(innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,CTF_R_k_p_wk_-M_k_p_wk_,CTF_R_k_p_wk_-M_k_p_wk_))/(2*pi);
ssnll_1_M_(1+nM) = 0.5 * (CTF_UX_R_l2_0 + UX_M_l2_0 - 2*Z_0);
%%%%%%%%;
clear R_k_p_wk_ R_k_q_wk_ S_k_p_wk_ T_k_p_wk_ CTF_R_k_p_wk_ M_k_p_wk_ M_k_q_wk_ CTF_UX_R_l2_0 CTF_UX_R_l2_1 UX_M_l2_0 UX_M_l2_1 Z_0 Z_1 ;
clear nCTF_from_nM nS_from_nM nw_from_nM ;
clear W_cazba_mk__
end;%for nM=0:n_M-1;
%%%%%%%%;
disp(sprintf(' %% ssnll_0_M_ vs ssnll_1_M_: %0.16f %%<-- should be <1e-2',fnorm(ssnll_0_M_(1:n_M_sub) - ssnll_1_M_(1:n_M_sub))/fnorm(ssnll_0_M_(1:n_M_sub))));

%%%%%%%%;
% Express likelihood using spherical-harmonics and fourier-bessel coordinates. ;
% Process all images. ;
%%%%%%%%;
M_k_q_qkM___ = zeros(n_m_max,n_k_p_r,n_M);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq=m_val; if (nq<0); nq=nq+n_w_max; end;
M_k_q_qkM___(1+l_max_max+m_val,1+nk_p_r,:) = M_k_q_wkM__(1+nq+nk_p_r*n_w_max,:);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
u_polar_a_ = unique(euler_polar_a_true_M_); n_u_polar_a = numel(u_polar_a_);
index_nM_assoc_ua_aM__ = cell(n_u_polar_a);
n_index_nM_assoc_ua_a_ = zeros(n_u_polar_a,1);
for nu_polar_a=0:n_u_polar_a-1;
u_polar_a = u_polar_a_(1+nu_polar_a);
tmp_index_ = efind(euler_polar_a_true_M_==u_polar_a);
n_index_nM_assoc_ua_a_(1+nu_polar_a) = numel(tmp_index_);
index_nM_assoc_ua_aM__{1+nu_polar_a} = tmp_index_;
clear u_polar_a tmp_index_ ;
end;%for nu_polar_a=0:n_u_polar_a-1;
%%%%%%%%;
W_betazeta_amlm____ = cell(n_u_polar_a,1);
for nu_polar_a=0:n_u_polar_a-1;
u_polar_a = u_polar_a_(1+nu_polar_a);
W_beta__ = wignerd_b(l_max_max,-u_polar_a);
V_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
tmp_index_ = l_max_max+[-l_val:+l_val];
V_betazeta_mlm___(1+tmp_index_,1+l_val,1+tmp_index_) = bsxfun(@times,W_beta__{1+l_val},reshape(zeta_lm__(1+l_val,1+tmp_index_),[1+2*l_val,1]));
end;%for l_val=0:l_max_max;
flag_check=0;
if flag_check;
W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
disp(sprintf(' %% V_betazeta_mlm___ vs W_betazeta_mlm___: %0.16f %%<-- should be <1e-6',fnorm(V_betazeta_mlm___-W_betazeta_mlm___)/fnorm(V_betazeta_mlm___)));
end;%if flag_check;
W_betazeta_amlm____{1+nu_polar_a} = V_betazeta_mlm___;
clear u_polar_a W_beta__ W_betazeta_mlm___ V_betazeta_mlm___ ;
end;%for nu_polar_a=0:n_u_polar_a-1;
%%%%%%%%;
%ssnll_ddPdFdF_k_ = zeros(n_k_p_r,1);
ssnll_2 = 0;
ssnll_2_M_ = zeros(n_M,1);
for nu_polar_a=0:n_u_polar_a-1;
u_polar_a = u_polar_a_(1+nu_polar_a);
index_nM_assoc_ua_M_ = index_nM_assoc_ua_aM__{1+nu_polar_a};
n_index_nM_assoc_ua = numel(index_nM_assoc_ua_M_);
assert(numel(index_nM_assoc_ua_M_)==n_index_nM_assoc_ua);
if (verbose); disp(sprintf(' %% nu_polar_a %d/%d n_index_nM_assoc_ua %d',nu_polar_a,n_u_polar_a,n_index_nM_assoc_ua)); end;
W_betazeta_mlm___ = W_betazeta_amlm____{1+nu_polar_a};
for nl=0:n_index_nM_assoc_ua-1;
nM = index_nM_assoc_ua_M_(1+nl);
nCTF_from_nM = index_nCTF_from_nM_(1+nM);
tmp_polar_a = euler_polar_a_true_M_(1+nM);
assert(tmp_polar_a==u_polar_a);
tmp_azimu_b = euler_azimu_b_true_M_(1+nM);
tmp_gamma_z = euler_gamma_z_true_M_(1+nM);
tmp_M_qk__ = M_k_q_qkM___(:,:,1+nM);
tmp_S_qk__ = zeros(n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF_from_nM);
a_k_Y_form_lm__ = a_k_Y_form_lmk___(:,:,1+nk_p_r);
ba_lm_ = reshape(bsxfun(@times,reshape(exp(+i*m_max_*tmp_azimu_b),[1,n_m_max]),a_k_Y_form_lm__),[(1+l_max_max)*n_m_max,1]);
tmp_S_qk__(:,1+nk_p_r) = (pi^2)*CTF_k_p_r*bsxfun(@times,reshape(exp(+i*m_max_*tmp_gamma_z),[n_m_max,1]),reshape( W_betazeta_mlm___,[n_m_max,(1+l_max_max)*n_m_max]) * ba_lm_ );
clear CTF_k_p_r  a_k_Y_form_lm__ ba_lm_;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_T_qk__ = tmp_M_qk__ - tmp_S_qk__;
tmp_ssnll = 0.5*sum(abs(tmp_T_qk__).^2*weight_2d_k_p_r_)/max(1,n_w_max); %<-- average over q. ;
ssnll_2_M_(1+nM) = tmp_ssnll;
ssnll_2 = ssnll_2 + tmp_ssnll;
clear nM nCTF_from_nM tmp_azimu_b tmp_gamma_z tmp_polar_a tmp_T_qk__ tmp_M_qk__ tmp_S_qk__ tmp_ssnll ;
end;%for nl=0:n_index_nM_assoc_ua-1;
clear W_betazeta_mlm___ ;
end;%for nu_polar_a=0:n_u_polar_a-1;
%%%%%%%%;
disp(sprintf(' %% ssnll_0_M_ vs ssnll_2_M_: %0.16f %%<-- should be <1e-2',fnorm(ssnll_0_M_(1:n_M) - ssnll_2_M_(1:n_M))/fnorm(ssnll_0_M_(1:n_M))));
%%%%%%%%;





