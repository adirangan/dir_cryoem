%%%%%%%%;
% Sets up a simple volume (in spherical-harmonic-coordinates), ;
% then tests the various template-generators. ;
% These calculations consider anisotropic CTF. ;
%%%%%%%%;

str_thisfunction = 'test_transforms_7_dr';

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
k_int = 16;
k_eq_d_double = 1.0;
template_k_eq_d_double = 0.5;
n_w_int = 1.0;
n_x_c = max(64,2*k_int);
tolerance_master = 1e-2;
flag_disp=1; nf=0;

dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_jpg = sprintf('%s/dir_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
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
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
weight_xxx_c = (2*x_p_r_max/n_x_c)^3;

%%%%%%%%;
% Now set up and test k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
flag_uniform_over_n_k_p_r = 1;
flag_uniform_over_polar_a = 0;
[ ...
 n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ;
%%%%%%%%;
% For eventual reconstruction we generate quadrature on shell. ;
%%%%%%%%;
qref_k_eq_d = k_eq_d/max(1e-12,k_p_r_max);
[ ...
 qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,qref_k_eq_d ...
,str_L ...
,flag_uniform_over_polar_a ...
);
tmp_index_ = [n_qk_csum_(end-1):n_qk_csum_(end)-1];
fnorm_disp(flag_verbose,'k_p_azimu_b_qk_(1+tmp_index_)',k_p_azimu_b_qk_(1+tmp_index_),'qref_azimu_b_shell_',qref_azimu_b_shell_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'k_p_polar_a_qk_(1+tmp_index_)',k_p_polar_a_qk_(1+tmp_index_),'qref_polar_a_shell_',qref_polar_a_shell_,' %%<-- should be zero');
%%%%%%%%;
% Visualization. ;
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; n_plot = p_row*p_col;
for nplot=0:n_plot-1;
nk_p_r = max(0,min(n_k_p_r-1,round(n_k_p_r*nplot/n_plot)));
tmp_index_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
subplot(p_row,p_col,1+nplot);
plot3(k_c_0_qk_(1+tmp_index_),k_c_1_qk_(1+tmp_index_),k_c_2_qk_(1+tmp_index_),'.');
axis equal; axis vis3d; axisnotick3d;
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
end;%for nplot=0:n_plot-1;
end;%if (flag_disp>1);
%%%%%%%%;
delta_a_c_3s__ = [  ...
 +1.5 , -0.5 ...
;-0.5 , -1.5 ...
;+0.3 , +2.0 ...
] / 2 / k_p_r_max ;
delta_b_c_3s__ = [  ...
 -0.5 , +0.8 ...
;-1.0 , +0.2 ...
;+1.2 , -0.7 ...
] / 2 / k_p_r_max ;
n_source = size(delta_a_c_3s__,2);
a_k_p_form_ = zeros(n_qk,1);
b_k_p_form_ = zeros(n_qk,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_a_c_(1+0) + k_c_1_qk_*delta_a_c_(1+1) + k_c_2_qk_*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_b_c_(1+0) + k_c_1_qk_*delta_b_c_(1+1) + k_c_2_qk_*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
I_a_quad = sum(a_k_p_form_.*weight_3d_k_p_qk_);
I_b_quad = sum(b_k_p_form_.*weight_3d_k_p_qk_);
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h3d_(tmp_kd)*k_p_r_max^3;
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h3d_(tmp_kd)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'I_a_form',I_a_form,'I_a_quad',I_a_quad,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'I_b_form',I_b_form,'I_b_quad',I_b_quad,' %%<-- should be <1e-6');
%%%%%%%%;
a_k_p_l2_quad = sum(conj(a_k_p_form_).*a_k_p_form_.*weight_3d_k_p_qk_);
a_k_p_l2_form = 0;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
delta_a_c_0_ = delta_a_c_3s__(:,1+nsource0);
delta_a_c_1_ = delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_0_ - delta_a_c_1_);
tmp_h3d = 4*pi/3; if abs(tmp_kd)>1e-12; tmp_h3d = h3d_(tmp_kd); end;
a_k_p_l2_form = a_k_p_l2_form + tmp_h3d*k_p_r_max^3;
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
fnorm_disp(flag_verbose,'a_k_p_l2_form',a_k_p_l2_form,'a_k_p_l2_quad',a_k_p_l2_quad,' %%<-- should be <1e-6');
%%%%%%%%;

%%%%%%%%;
% Now set up and test polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
tmp_S_delta_x_c_ = [cos(pi/4);sin(pi/4)]/max(1e-12,k_p_r_max);
tmp_S_phi = +1*pi/5;
tmp_T_delta_x_c_ = [cos(-pi/3);sin(-pi/3)]/max(1e-12,k_p_r_max);
tmp_T_phi = -3*pi/7;
tmp_S_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_c_1_wk_*tmp_S_delta_x_c_(1+1)));
tmp_plane_S_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - tmp_S_phi);
tmp_T_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_T_delta_x_c_(1+0) + k_c_1_wk_*tmp_T_delta_x_c_(1+1)));
tmp_plane_T_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - tmp_T_phi);
tmp_S_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
tmp_T_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_T_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
if (flag_disp>1);
figure(1+nf);nf=nf+1;figbig; np=0;
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_S_k_p_wk_)); axis image; axisnotick; title('real(tmp_S_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_S_x_c_xx_)); axis image; axisnotick; title('real(tmp_S_x_c_xx_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_T_k_p_wk_)); axis image; axisnotick; title('real(tmp_T_k_p_wk_)','Interpreter','none');
subplot(2,2,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_T_x_c_xx_)); axis image; axisnotick; title('real(tmp_T_x_c_xx_)','Interpreter','none');
end;%if (flag_disp>1);
I_quad = sum(conj(tmp_plane_S_k_p_wk_.*tmp_S_k_p_wk_).*(tmp_plane_T_k_p_wk_.*tmp_T_k_p_wk_).*weight_2d_wk_)*(2*pi)^2;
I_form = I_xPPx_0(k_p_r_max,tmp_S_phi,tmp_S_delta_x_c_,tmp_T_phi,tmp_T_delta_x_c_);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %%<-- should be <1e-6');
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
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c_3s__(:,1+nsource),l_max_);
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
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_quad_',a_k_Y_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_quad_ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_quad_',a_k_p_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_reco_ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_reco_',a_k_p_reco_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'a_k_p_quad_',a_k_p_quad_,'a_k_p_reco_',a_k_p_reco_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_Y_reco_ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_reco_',a_k_Y_reco_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'a_k_Y_quad_',a_k_Y_quad_,'a_k_Y_reco_',a_k_Y_reco_,' %%<-- should be <1e-2');
%%%%%%%%;
a_k_Y_form_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_form_);
a_k_Y_form_yk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
tmp_a_k_Y_form_lm_ = a_k_Y_form_(1+tmp_index_);
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
a_k_Y_form_yk___(:,:,1+nk_p_r) = tmp_a_k_Y_form_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
flag_check=1;
if flag_check;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max=l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_yk__(1+l_val*(l_val+1)+m_val,1+nk_p_r));
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_yk___(1+l_val,1+l_max_max+m_val,1+nk_p_r));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
end;%if flag_check;
%%%%%%%%;
a_k_Y_form_lkm___ = permute(a_k_Y_form_yk___,[1,3,2]);
%%%%%%%%;
a_k_Y_l2_quad = sum(conj(a_k_Y_form_yk__).*a_k_Y_form_yk__*reshape(weight_3d_k_p_r_,[n_k_p_r,1]));

%%%%%%%%;
% define rotations. ;
%%%%%%%%;
R2 = @(gamma_z) ...
[ +cos(gamma_z) -sin(gamma_z) ; ...
  +sin(gamma_z) +cos(gamma_z) ; ...
] ;
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
 0*flag_verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_form_yk__,[n_lm_max,n_k_p_r]) ...
,template_k_eq_d_double/max(1e-12,k_p_r_max) ...
,-1 ...
,n_w_max ...
);
n_S = n_viewing_S;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
S_k_p_l2_quad_S_ = reshape(sum(bsxfun(@times,conj(S_k_p_wkS__).*S_k_p_wkS__,reshape(weight_2d_wk_,[n_w_sum,1])),1),[n_S,1])*(2*pi)^2;
%%%%%%%%;
% Now step through and reconstitute the templates. ;
%%%%%%%%;
T_k_p_l2_quad_S_ = zeros(n_S,1);
T_k_p_l2_form_S_ = zeros(n_S,1);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_l2_quad_S_(1+nS) = sum(conj(T_k_p_wk_).*T_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
tmp_delta_a_c_0_ = tmp_R__*delta_a_c_3s__(:,1+nsource0);
tmp_delta_a_c_1_ = tmp_R__*delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(tmp_delta_a_c_0_(1:2) - tmp_delta_a_c_1_(1:2));
tmp_h2d = (2*pi)^2; if abs(tmp_kd)>1e-12; tmp_h2d = h2d_(tmp_kd); end;
T_k_p_l2_form_S_(1+nS) = T_k_p_l2_form_S_(1+nS) + tmp_h2d/(2*pi)^2 * (pi*k_p_r_max^2);
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
if (flag_disp>1);
if mod(nS,128)==0;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_)); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(1,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
sgtitle(sprintf('nS %d/%d',nS,n_S));
end;%if mod(nS,128)==0;
end;%if (flag_disp>1);
T_k_p_wkS__(:,1+nS) = T_k_p_wk_;
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'S_k_p_wkS__',S_k_p_wkS__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_l2_form_S_',T_k_p_l2_form_S_,'T_k_p_l2_quad_S_',T_k_p_l2_quad_S_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_l2_form_S_',T_k_p_l2_form_S_,'S_k_p_l2_quad_S_',S_k_p_l2_quad_S_,' %%<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%;
% Now test out direct spherical-harmonic contruction of S_k_q_wk_. ;
% Note that below we introduce a (unnecessary) tmp_gamma_z to check consistency. ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*2/3)));
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
%tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = pi/12;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'T_k_p_wk_',T_k_p_wk_,' %%<-- should be <1e-2');
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'T_k_q_wk_',T_k_q_wk_,' %%<-- should be <1e-2');
%%%%;
S_k_p_l2_qua2 = sum(conj(S_k_p_wk_).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
S_k_p_l2_quad = S_k_p_l2_quad_S_(1+nS);
fnorm_disp(flag_verbose,'S_k_p_l2_quad',S_k_p_l2_quad,'S_k_p_l2_qua2',S_k_p_l2_qua2,' %%<-- should be <1e-9');
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
 0*flag_verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,tmp_euler_ ...
);
%%%%;
a_k_Y_rota_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_rota_);
a_k_Y_l2_rota = sum(conj(a_k_Y_rota_yk__).*a_k_Y_rota_yk__*reshape(weight_3d_k_p_r_,[n_k_p_r,1]));
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
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
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
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
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
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
Q_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,Q_k_q_wk_);
R_k_p_l2_qua2 = sum(conj(R_k_p_wk_).*R_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
Q_k_p_l2_qua2 = sum(conj(Q_k_p_wk_).*Q_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
%%%%;
fnorm_disp(flag_verbose,'R_k_q_wk_',R_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'R_k_q_wk_',R_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_q_wk_',T_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_q_wk_',T_k_q_wk_,'R_k_q_wk_',R_k_q_wk_,' %%<-- should be <1e-2');
%%%%%%%%;
fnorm_disp(flag_verbose,'R_k_p_wk_',R_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'R_k_p_wk_',R_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_wk_',T_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_wk_',T_k_p_wk_,'R_k_p_wk_',R_k_p_wk_,' %%<-- should be <1e-2');
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 4; np=0;
Slim_ = max(abs(S_k_p_wk_))*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(Q_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('real(Q_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(Q_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('imag(Q_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('real(R_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(R_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('imag(R_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('real(S_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('imag(S_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('real(T_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_wk_),Slim_,colormap_beach()); axisnotick; axis image; title(sprintf('imag(T_k_p_wk_)'),'Interpreter','none');
sgtitle(sprintf(' nS %d/%d',nS,n_S),'Interpreter','none');
end;%if (flag_disp>0);
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 4; np=0;
Slim_ = max(abs(S_k_q_wk_))*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(Q_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('real(Q_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(Q_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('imag(Q_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('real(R_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(R_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('imag(R_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('real(S_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('imag(S_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('real(T_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_q_wk_),Slim_,colormap_beach()); axisnotick; axis tight; title(sprintf('imag(T_k_q_wk_)'),'Interpreter','none');
sgtitle(sprintf(' nS %d/%d',nS,n_S),'Interpreter','none');
end;%if (flag_disp>0);
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
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_sub_wkb__(:,1+nazimu_b_use) = T_k_p_wk_;
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
S_k_q_sub_wkb__(:,1+nazimu_b_use) = S_k_q_wk_;
T_k_q_sub_wkb__(:,1+nazimu_b_use) = T_k_q_wk_;
clear S_k_p_wk_ T_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_q_sub_wkb__',S_k_q_sub_wkb__,'T_k_q_sub_wkb__',T_k_q_sub_wkb__,' %%<-- should be <1e-2');
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use))); axis image; axisnotick;
%title(sprintf('real(S_k_p_sub_wkb__(:,1+%d))',nazimu_b_use),'Interpreter','none');
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
end;%if (flag_disp>1);
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
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_betazeta_mlm___: %0.2fs',tmp_t)); end;
%%%%;
flag_check=1;
if flag_check;
tmp_w_ = randn(n_m_max) + i*randn(n_m_max);
tmp_azimu_b_use_ = 2*pi*rand(n_azimu_b_use,1);
tmp_f__ = exp(-i*(-reshape(tmp_azimu_b_use_,[n_azimu_b_use,1]))*reshape(m_max_,[1,n_m_max]));
tmp_fw_0_ = tmp_f__*tmp_w_;
tmp_fw_1_ = xxnufft1d2(n_azimu_b_use,tmp_azimu_b_use_,+1,1e-6,n_m_max,tmp_w_);
fnorm_disp(flag_verbose,'tmp_fw_0_',tmp_fw_0_,'tmp_fw_1_',tmp_fw_1_,' %%<-- should be <1e-6');
clear tmp_azimu_b_use_ tmp_w_ tmp_f__ tmp_fw_0_ tmp_fw_1_ ;
end;%if flag_check;
%%%%;
tmp_t = tic();
W_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag(exp(+i*m0_val_*tmp_gamma_z))*W_betazeta_ml__*a_k_Y_form_lk__ for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
W_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
 *diag(exp(-i*m_max_*(-tmp_gamma_z))) ...
 *reshape(W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_form_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_caza_mkm___: %0.2fs',tmp_t)); end;
tmp_t = tic();
W_caza_mmk___ = permute(W_caza_mkm___,[3,1,2]);
W_caza_bmk___ = reshape(xxnufft1d2(n_azimu_b_use,azimu_b_use_,+1,1e-6,n_m_max,reshape(W_caza_mmk___,[n_m_max,n_m_max*n_k_p_r])),[n_azimu_b_use,n_m_max,n_k_p_r]);
W_caza_mkb___ = permute(W_caza_bmk___,[2,3,1]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_caza_mkb___: %0.2fs',tmp_t)); end;
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
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,' %%<-- should be <1e-2');
%%%%;
if (flag_disp>1);
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
end;%if (flag_disp>1);
%%%%%%%%;
[ ...
 U_k_p_sub_wkb__ ...
,W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 0*flag_verbose ...
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
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%;
% We check the following sum: ;
% \sum_{nq} zeta_lm__(1+l0,1+nq) * zeta_lm__(1+l1,1+nq) * W_beta__{1+l0}(1+nq,1+m0) * W_beta__{1+l1}(1+nq,1+m1) ;
% Verdict: not too sparse. ;
%%%%%%%%;
flag_check=1;
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
fnorm_disp(flag_verbose,'sum0_lmlm____',sum0_lmlm____,'sum1_lmlm____',sum1_lmlm____,' %%<-- should be <1e-2');
%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
S_k_q_wkS__(:,1+nS) = S_k_q_wk_;
end;%for nS=0:n_S-1;
tmp_ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__),[n_w_sum,n_S]); %<-- this accomplishes the same thing. ;
fnorm_disp(flag_verbose,'S_k_q_wkS__',S_k_q_wkS__,'tmp_',tmp_,' %%<-- should be zero');
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
%%%%%%%%;

flag_delta = 1;
%%%%%%%%;
if ~flag_delta; delta_r_max = 0.0/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested =  1; end;
if  flag_delta; delta_r_max = 0.5/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested = 128; end;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
%%%%%%%%;
n_delta_v = FTK.n_delta_v;

flag_CTF = 1;
%%%%%%%%;
% Now we generate images for the likelihood-calculation. ;
% This test accounts for anisotropic CTF. ;
% We also shift each image by an 'on-grid' displacement. ;
%%%%%%%%;
n_CTF = 3;
CTF_phi_C_ = pi*[+2/3;-1/5;+6/7];
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if ~flag_CTF; CTF_k_p_wkC__(:,1+nCTF) = ones(n_w_sum,1); end;
if  flag_CTF; CTF_k_p_wkC__(:,1+nCTF) = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTF_phi_C_(1+nCTF)); end;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
tmp_t = tic();
n_M = 2*n_S; %<-- pick something bigger than n_S;
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
index_nCTF_from_nM_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
euler_azimu_b_true_M_ = zeros(n_M,1);
euler_polar_a_true_M_ = zeros(n_M,1);
index_nw_from_nM_ = zeros(n_M,1);
euler_gamma_z_true_M_ = zeros(n_M,1);
index_nd_from_nM_ = zeros(n_M,1);
image_delta_x_true_M_ = zeros(n_M,1);
image_delta_y_true_M_ = zeros(n_M,1);
rng(0);
for nM=0:n_M-1;
nCTF = mod(nM,n_CTF);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
index_nCTF_from_nM_(1+nM) = nCTF;
nS = max(0,min(n_S-1,mod(nM,n_S)));
index_nS_from_nM_(1+nM) = nS;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
azimu_b = viewing_azimu_b_S_(1+nS);
polar_a = viewing_polar_a_S_(1+nS);
index_gamma_z = periodize(nM,0,n_gamma_z);
gamma_z = gamma_z_(1+index_gamma_z);
index_nw_from_nM_(1+nM) = index_gamma_z;
euler_azimu_b_true_M_(1+nM) = azimu_b;
euler_polar_a_true_M_(1+nM) = polar_a;
euler_gamma_z_true_M_(1+nM) = gamma_z;
index_nd_from_nM_(1+nM) = periodize(nM,0,n_delta_v);
delta_x = FTK.delta_x_(1+index_nd_from_nM_(1+nM)); delta_y = FTK.delta_y_(1+index_nd_from_nM_(1+nM));
image_delta_x_true_M_(1+nM) = delta_x;
image_delta_y_true_M_(1+nM) = delta_y;
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z); %<-- note here we rotate S_k_p_wk_ by +gamma_z to form M_k_p_wk_. ;
T_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_wk_,-delta_x,-delta_y); %<-- note here we translate T_k_p_wk_ by -[delta_x,delta_y]. ;
M_k_p_wk_ = CTF_k_p_wk_.*T_k_p_wk_;
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_wkM__(:,1+nM) = M_k_q_wk_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_p_wkM__: %0.2fs',tmp_t)); end;
clear nCTF CTF_k_p_wk_ nS azimu_b polar_a index_gamma_z gamma_z S_k_p_wk_ T_k_p_wk_ M_k_p_wk_ M_k_q_wk_ ;
%%%%%%%%;

%%%%%%%%;
% quickly test inner-product: <R(+gamma_z)*M_k_p_wk_,S_k_p_wk_> ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*2/3)));
nM = max(0,min(n_M-1,round(n_M*1/5)));
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM));
M_k_p_wk_ = CTF_k_p_wk_.*M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
tmp_1_w_ = ifft(sum(bsxfun(@times,reshape(conj(M_k_q_wk_).*S_k_q_wk_,[n_w_max,n_k_p_r]),reshape(weight_2d_k_p_r_,[1,n_k_p_r])),2));
tmp_0_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
tmp_0_w_(1+nw) = sum(conj(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_,+gamma_z)).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'tmp_0_w_',tmp_0_w_,'tmp_1_w_',tmp_1_w_,' %%<-- should be zero');
%%%%%%%%;

%%%%%%%%;
% Now construct the template-norms: ;
% <R(+gamma_z)*CTF_k_p_wk_.*S_k_p_wk_,R(+gamma_z)*CTF_k_p_wk_.*S_k_p_wk_> ;
% Note that this does not involve collapsing onto principal-modes. ;
%%%%%%%%;
CTF_S_l2_wSC_quad___ = zeros(n_w_max,n_S,n_CTF);
SS_k_p_wkS__ = conj(S_k_p_wkS__).*S_k_p_wkS__;
SS_k_q_wkS__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,SS_k_p_wkS__),[n_w_sum,n_S]);
CC_k_p_wkC__ = conj(CTF_k_p_wkC__).*CTF_k_p_wkC__;
CC_k_q_wkC__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,CC_k_p_wkC__),[n_w_sum,n_CTF]);
for nCTF=0:n_CTF-1;
CTF_S_l2_wSC_quad___(:,:,1+nCTF) = ifft(squeeze(sum(bsxfun(@times,reshape(bsxfun(@times,conj(CC_k_q_wkC__(:,1+nCTF)),SS_k_q_wkS__),[n_w_max,n_k_p_r,n_S]),reshape(weight_2d_k_p_r_,[1,n_k_p_r,1])),2)));
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% Now check one template. ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*1/5)));
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R_S__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
nCTF = max(0,min(n_CTF-1,round(n_CTF*1/3)));
tmp_CTF_phi = CTF_phi_C_(1+nCTF);
CTF_S_l2_w_quad_ = CTF_S_l2_wSC_quad___(:,1+nS,1+nCTF);
CTF_S_l2_w_qua2_ = zeros(n_w_max,1);
CTF_S_l2_w_form_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
CS_k_p_wk_ = S_k_p_wkS__(:,1+nS).*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),+gamma_z);
CTF_S_l2_w_qua2_(1+nw) = sum(conj(CS_k_p_wk_).*CS_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
tmp_S_delta_0_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource0); tmp_S_delta_0_2_ = tmp_S_delta_0_3_(1:2);
for nsource1=0:n_source-1;
tmp_S_delta_1_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource1); tmp_S_delta_1_2_ = tmp_S_delta_1_3_(1:2);
CTF_S_l2_w_form_(1+nw) = CTF_S_l2_w_form_(1+nw) + I_xPPx_0(k_p_r_max,tmp_CTF_phi+gamma_z,tmp_S_delta_0_2_,tmp_CTF_phi+gamma_z,tmp_S_delta_1_2_);
end;%for nsource0=0:n_source-1;
end;%for nsource1=0:n_source-1;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'CTF_S_l2_w_qua2_',CTF_S_l2_w_qua2_,'CTF_S_l2_w_quad_',CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'CTF_S_l2_w_qua2_',CTF_S_l2_w_qua2_,'CTF_S_l2_w_form_',CTF_S_l2_w_form_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'CTF_S_l2_w_quad_',CTF_S_l2_w_quad_,'CTF_S_l2_w_form_',CTF_S_l2_w_form_,' %%<-- should be <1e-6');
%%%%%%%%;

flag_pm = 1;
%%%%%%%%%%%%%%%%;
% Now calculate innerproduct Z_dwSM____. ;
% Here we account for anisotropic CTF. ;
%%%%%%%%%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
Z_dwSM_ampm____ = zeros(n_delta_v,n_w_max,n_S,n_M);
UX_M_l2_M_ = zeros(n_M,1);
UX_M_l2_dM__ = zeros(n_delta_v,n_M);
UX_knC___ = zeros(n_k_p_r,n_UX_rank,n_CTF);
X_weight_rC__ = zeros(n_k_p_r,n_CTF);
%%%%%%%%;
for nCTF=0:n_CTF-1;
if (flag_verbose); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
%%%%;
% Prepare principal-modes. ;
%%%%;
UX_kn__ = eye(n_k_p_r,n_UX_rank);
X_weight_r_ = sqrt(weight_2d_k_p_r_);
if  flag_pm;
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M_sub ...
,M_k_p_wkM__(:,1+index_M_sub_) ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
end;%if  flag_pm;
UX_knC___(:,:,1+nCTF) = UX_kn__;
X_weight_rC__(:,1+nCTF) = X_weight_r_;
%%%%;
% Prepare quasi-images. ;
%%%%;
tmp_t = tic();
CTF_M_sub_k_p_wkM__ = bsxfun(@times,CTF_k_p_wk_,M_k_p_wkM__(:,1+index_M_sub_));
CTF_M_sub_k_q_wkM__ = zeros(n_w_sum,n_M_sub);
for nM_sub=0:n_M_sub-1;
CTF_M_sub_k_q_wkM__(:,1+nM_sub) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_M_sub_k_p_wkM__(:,1+nM_sub));
end;%for nM_sub=0:n_M_sub-1;
M_sub_k_q_wkM__ = M_k_q_wkM__(:,1+index_M_sub_);
svd_VUXCTFM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,CTF_M_sub_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
svd_VUXM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_sub_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmpmh_VUXM_lwnM____3: %0.2fs',tmp_t)); end;
%%%%;
% Now calculate norms of the translated images. ;
%%%%;
tmp_t = tic();
UX_M_sub_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_sub_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_UX_M_sub_l2_dm__1: %0.2fs',tmp_t)); end;
tmp_index_d0 = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); assert(numel(tmp_index_d0)==1); %<-- should be zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_M_sub_l2_dM__(1+tmp_index_d0,:),[n_M_sub,1]);
%%%%;
% Prepare UX_S_k_q_wnS__. ;
%%%%;
tmp_t = tic();
UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% UX_S_k_q_wnS__: %0.2fs',tmp_t)); end;
%%%%;
% Calculate Z_sub_dwSM____. ;
%%%%;
tmp_t = tic();
UX_S_k_q_nSw___ = permute(reshape(UX_S_k_q_wnS__,[n_w_max,n_UX_rank,n_S]),[2,3,1]);
tmp_t_sub = tic();
svd_VUXCTFM_sub_nMwl____ = permute(svd_VUXCTFM_sub_lwnM____,[3,4,2,1]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose); disp(sprintf(' %% svd_VUXCTFM_sub_nMwl____: %0.2fs',tmp_t_sub)); end;
tmp_t_sub = tic();
svd_SVUXCTFM_sub_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_SVUXCTFM_sub_SMwl____(:,:,1+nw,1+nl) = ctranspose(UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXCTFM_sub_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_SMwl____: %0.6f',tmp_t_sub)); end;
tmp_t_sub = tic();
svd_SVUXCTFM_sub_lwSM____ = permute(ifft(permute(svd_SVUXCTFM_sub_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_sub_lwSM____: %0.6f',tmp_t_sub)); end;
tmp_t_sub = tic(); nop=0;
svd_USESVUXCTFM_sub_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXCTFM_sub_lwSM____,[FTK.n_svd_l,n_w_max*n_S*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]));
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_USESVUXCTFM_sub_dwSM____: %0.6f',tmp_t_sub)); end;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% Z_sub_dwSM____: %0.2fs',tmp_t)); end;
%%%%;
% Store results. ;
%%%%;
UX_M_l2_M_(1+index_M_sub_) = UX_M_sub_l2_M_;
UX_M_l2_dM__(:,1+index_M_sub_) = UX_M_sub_l2_dM__;
Z_dwSM_ampm____(:,:,:,1+index_M_sub_) = svd_USESVUXCTFM_sub_dwSM____ ;
%%%%;
clear UX_kn__ X_weight_r_ CTF_k_p_wk_ index_M_sub_ UX_M_sub_l2_dM__ UX_M_sub_l2_M_ ;
clear UX_S_k_q_wnS__ ;
clear svd_VUXCTFM_sub_lwnM____ svd_VUXCTFM_sub_nMwl____ ;
clear svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_lwSM____ ;
clear svd_USESVUXCTFM_sub_dwSM____ ;
%%%%%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
clear nCTF;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now estimate landscape of innerproducts across delta_x_ and delta_y_. ;
% Limited to a single image-template pair and a fixed nw. ;
%%%%%%%%;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
nS = max(0,min(n_S-1,round(n_S*2/5)));
azimu_b_S = viewing_azimu_b_S_(1+nS);
polar_a_S = viewing_polar_a_S_(1+nS);
gamma_z_S = 0.0;
R_S__ = Rz(+gamma_z_S)*Ry(-polar_a_S)*Rz(-azimu_b_S);
nM = max(0,min(n_M-1,round(n_M*4/5)));
index_nS_from_nM = index_nS_from_nM_(1+nM);
azimu_b_M = euler_azimu_b_true_M_(1+nM);
polar_a_M = euler_polar_a_true_M_(1+nM);
gamma_z_M = euler_gamma_z_true_M_(1+nM);
index_nd_from_nM = index_nd_from_nM_(1+nM);
delta_x_M = image_delta_x_true_M_(1+nM);
delta_y_M = image_delta_y_true_M_(1+nM);
R_M__ = Rz(+gamma_z_M)*Ry(-polar_a_M)*Rz(-azimu_b_M);
nw = max(0,min(n_w_max-1,round(n_w_max*3/5)));
Z_d_ampm_ = Z_dwSM_ampm____(:,1+nw,1+nS,1+nM);
Z_d_quad_ = zeros(n_delta_v,1);
Z_d_form_ = zeros(n_delta_v,1);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_phi = CTF_phi_C_(1+nCTF);
%%%%;
gamma_z = (2*pi*nw)/n_w_max;
RS_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
CTF_RS_k_p_wk_ = RS_k_p_wk_.*CTF_k_p_wk_;
CTF_RS_k_p_l2 = sum(conj(CTF_RS_k_p_wk_).*CTF_RS_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
T0M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+0.0,+0.0);
T0M_k_p_l2 = sum(conj(T0M_k_p_wk_).*T0M_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
%%;
for ndelta_v=0:n_delta_v-1;
%%;
delta_x = FTK.delta_x_(1+ndelta_v); delta_y = FTK.delta_y_(1+ndelta_v);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
CTF_RS_k_p_TM_k_p = sum(conj(CTF_RS_k_p_wk_).*TM_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
Z_d_quad_(1+ndelta_v) = CTF_RS_k_p_TM_k_p;
%%;
%%%%;
Z_d_form = 0;
for nsource_S=0:n_source-1;
delta_S_3_ = R_S__*delta_a_c_3s__(:,1+nsource_S);
delta_S_2_ = R2(+gamma_z) * delta_S_3_(1:2);
for nsource_M=0:n_source-1;
delta_M_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M);
delta_M_2_ = delta_M_3_(1:2) - [delta_x;delta_y] + [delta_x_M;delta_y_M];
tmp_I = I_xPPx_0(k_p_r_max,CTF_phi,delta_S_2_,CTF_phi,delta_M_2_);
Z_d_form = Z_d_form + tmp_I;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
Z_d_form_(1+ndelta_v) = Z_d_form;
%%;
end;%for ndelta_v=0:n_delta_v-1;
%%;
fnorm_disp(flag_verbose,'Z_d_quad_',Z_d_quad_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_quad_',Z_d_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
Zlim_ = prctile(real(Z_d_form_),[ 5,95],'all');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_form_),Zlim_); axisnotick; axis image; title('real(Z_d_form_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_quad_),Zlim_); axisnotick; axis image; title('real(Z_d_quad_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_ampm_),Zlim_); axisnotick; axis image; title('real(Z_d_ampm_)','Interpreter','none');
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now try and pick out optimal poses from displacements. ;
%%%%%%%%;
euler_azimu_b_opti_M_ = zeros(n_M,1);
euler_polar_a_opti_M_ = zeros(n_M,1);
euler_gamma_z_opti_M_ = zeros(n_M,1);
image_delta_x_opti_M_ = zeros(n_M,1);
image_delta_y_opti_M_ = zeros(n_M,1);
for nM=0:n_M-1;
if (flag_verbose>0); if mod(nM,32)==0; disp(sprintf(' %% nM %.3d/%.3d',nM,n_M)); end; end;
nCTF = index_nCTF_from_nM_(1+nM);
tmp_UX_M_l2_d_ = UX_M_l2_dM__(:,1+nM);
tmp_CTF_S_l2_wS__ = flipud(CTF_S_l2_wSC_quad___(:,:,1+nCTF));
tmp_Z_dwS_ampm___ = Z_dwSM_ampm____(:,:,:,1+nM);
tmp_X_dwS_ampm___ = bsxfun(@rdivide,bsxfun(@rdivide,tmp_Z_dwS_ampm___,max(1e-12,reshape(sqrt(tmp_CTF_S_l2_wS__),[1,n_w_max,n_S]))),max(1e-12,reshape(sqrt(tmp_UX_M_l2_d_),[n_delta_v,1,1])));
[~,tmp_ij] = max(tmp_X_dwS_ampm___,[],'all','linear'); tmp_index = tmp_ij-1;
tmp_index_tmp = tmp_index;
tmp_nd = mod(tmp_index_tmp,n_delta_v); tmp_index_tmp = (tmp_index_tmp - tmp_nd)/max(1,n_delta_v);
tmp_nw = mod(tmp_index_tmp,  n_w_max); tmp_index_tmp = (tmp_index_tmp - tmp_nw)/max(1,  n_w_max);
tmp_nS = mod(tmp_index_tmp,      n_S); tmp_index_tmp = (tmp_index_tmp - tmp_nS)/max(1,      n_S);
assert(fnorm(tmp_index_tmp)==0);
euler_azimu_b_opti_M_(1+nM) = viewing_azimu_b_S_(1+tmp_nS);
euler_polar_a_opti_M_(1+nM) = viewing_polar_a_S_(1+tmp_nS);
euler_gamma_z_opti_M_(1+nM) = (2*pi*tmp_nw)/max(1,n_w_max);
image_delta_x_opti_M_(1+nM) = FTK.delta_x_(1+tmp_nd);
image_delta_y_opti_M_(1+nM) = FTK.delta_y_(1+tmp_nd);
end;%for nM=0:n_M-1;
%%%%%%%%;
flag_check=0;
%%%%%%%%;
if flag_check;
euler_azimu_b_pose_M_ = zeros(n_M,1);
euler_polar_a_pose_M_ = zeros(n_M,1);
euler_gamma_z_pose_M_ = zeros(n_M,1);
image_delta_x_pose_M_ = zeros(n_M,1);
image_delta_y_pose_M_ = zeros(n_M,1);
index_true_pose_prct_M_ = zeros(n_M,1);
for nM=0:n_M-1;
if (flag_verbose>0); if mod(nM,32)==0; disp(sprintf(' %% nM %.3d/%.3d',nM,n_M)); end; end;
nCTF = index_nCTF_from_nM_(1+nM);
tmp_UX_M_l2_d_ = UX_M_l2_dM__(:,1+nM);
tmp_CTF_S_l2_wS__ = flipud(CTF_S_l2_wSC_quad___(:,:,1+nCTF));
tmp_Z_dwS_ampm___ = Z_dwSM_ampm____(:,:,:,1+nM);
tmp_X_dwS_ampm___ = bsxfun(@rdivide,bsxfun(@rdivide,tmp_Z_dwS_ampm___,max(1e-12,reshape(sqrt(tmp_CTF_S_l2_wS__),[1,n_w_max,n_S]))),max(1e-12,reshape(sqrt(tmp_UX_M_l2_d_),[n_delta_v,1,1])));
[tmp_X_dwS_sort_,tmp_ij_] = sort(tmp_X_dwS_ampm___(:),'ascend'); [~,tmp_ji_] = sort(tmp_ij_,'ascend'); tmp_index_ = tmp_ij_-1;
tmp_index_tmp_ = tmp_index_;
tmp_nd_ = mod(tmp_index_tmp_,n_delta_v); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nd_)/max(1,n_delta_v);
tmp_nw_ = mod(tmp_index_tmp_,  n_w_max); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nw_)/max(1,  n_w_max);
tmp_nS_ = mod(tmp_index_tmp_,      n_S); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nS_)/max(1,      n_S);
assert(fnorm(tmp_index_tmp_)==0);
opt_nd = tmp_nd_(end); opt_delta_x = FTK.delta_x_(1+opt_nd); opt_delta_y = FTK.delta_y_(1+opt_nd);
opt_nw = tmp_nw_(end); opt_gamma_z = gamma_z_(1+opt_nw);
opt_nS = tmp_nS_(end); opt_polar_a = viewing_polar_a_S_(1+opt_nS); opt_azimu_b = viewing_azimu_b_S_(1+opt_nS);
euler_azimu_b_pose_M_(1+nM) = opt_azimu_b;
euler_polar_a_pose_M_(1+nM) = opt_polar_a;
euler_gamma_z_pose_M_(1+nM) = opt_gamma_z;
image_delta_x_pose_M_(1+nM) = opt_delta_x;
image_delta_y_pose_M_(1+nM) = opt_delta_y;
tru_nd = index_nd_from_nM_(1+nM); tru_delta_x = image_delta_x_true_M_(1+nM); tru_delta_y = image_delta_y_true_M_(1+nM);
tru_nw = index_nw_from_nM_(1+nM); tru_gamma_z = gamma_z_(1+tru_nw);
tru_nS = index_nS_from_nM_(1+nM); tru_polar_a = euler_polar_a_true_M_(1+nM); tru_azimu_b = euler_azimu_b_true_M_(1+nM);
index_true_pose = tru_nd + (tru_nw + (tru_nS)*n_w_max)*n_delta_v ;
index_true_pose_prct = tmp_ji_(1+index_true_pose)/max(1,n_delta_v*n_w_max*n_S);
index_true_pose_prct_M_(1+nM) = index_true_pose_prct;
if (flag_verbose>1);
disp(sprintf(' %% nM %.3d/%.3d: ',nM,n_M));
disp(sprintf(' %% opt_nS %.3d tru_nS %.3d opt_nw %.3d tru_nw %.3d opt_nd %.3d tru_nd %.3d',opt_nS,tru_nS,opt_nw,tru_nw,opt_nd,tru_nd));
disp(sprintf(' %% opt_azimu_b %+0.2f opt_polar_a %+0.2f opt_gamma_z %+0.2f opt_delta_x %+0.2f opt_delta_y %+0.2f',opt_azimu_b,opt_polar_a,opt_gamma_z,opt_delta_x,opt_delta_y));
disp(sprintf(' %% tru_azimu_b %+0.2f tru_polar_a %+0.2f tru_gamma_z %+0.2f tru_delta_x %+0.2f tru_delta_y %+0.2f',tru_azimu_b,tru_polar_a,tru_gamma_z,tru_delta_x,tru_delta_y));
end;%if (flag_verbose>0);
end;%for nM=0:n_M-1;
%%%%%%%%;
fnorm_disp(flag_verbose,'euler_azimu_b_opti_M_',euler_azimu_b_opti_M_,'euler_azimu_b_pose_M_',euler_azimu_b_pose_M_);
fnorm_disp(flag_verbose,'euler_polar_a_opti_M_',euler_polar_a_opti_M_,'euler_polar_a_pose_M_',euler_polar_a_pose_M_);
fnorm_disp(flag_verbose,'euler_gamma_z_opti_M_',euler_gamma_z_opti_M_,'euler_gamma_z_pose_M_',euler_gamma_z_pose_M_);
fnorm_disp(flag_verbose,'image_delta_x_opti_M_',image_delta_x_opti_M_,'image_delta_x_pose_M_',image_delta_x_pose_M_);
fnorm_disp(flag_verbose,'image_delta_y_opti_M_',image_delta_y_opti_M_,'image_delta_y_pose_M_',image_delta_y_pose_M_);
%%%%%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% Now perform reconstruction. ;
%%%%%%%%;
parameter_kappa = struct('type','parameter');
parameter_kappa.kernel_basic_l_max_use = l_max_max;
parameter_kappa.kernel_basic_qref_k_eq_d_double = k_eq_d/max(1e-12,k_p_r_max);
parameter_kappa.flag_recalc_qref_from_data = 1;
parameter_kappa.flag_recalc_dtau_qref_from_data = 0;
parameter_kappa.flag_recalc_dtau_dtau_qref_from_data = 0;
KAPPA = [];
weight_imagecount_M_ = ones(n_M,1);
TM_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
delta_x = image_delta_x_true_M_(1+nM); delta_y = image_delta_y_true_M_(1+nM);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
TM_k_p_wkM__(:,1+nM) = TM_k_p_wk_;
end;%for nM=0:n_M-1;
tmp_t = tic();
[ ...
 parameter_kappa ...
,KAPPA ...
,a_from_M_C2M0_k_p_qk__ ...
,a_from_M_C1M1_k_p_qk__ ...
,a_from_M_C0M2_k_p_qk__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,+euler_polar_a_true_M_ ...
,+euler_azimu_b_true_M_ ...
,+euler_gamma_z_true_M_ ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,TM_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
a_k_p_from_M_ = reshape(a_from_M_C1M1_k_p_qk__./max(1e-12,a_from_M_C2M0_k_p_qk__),[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_M_',a_k_p_from_M_);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter_kappa ...
,~ ...
,a_from_S_C2M0_k_p_qk__ ...
,a_from_S_C1M1_k_p_qk__ ...
,a_from_S_C0M2_k_p_qk__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,[] ...
,n_w_max ...
,n_S ...
,ones(n_S,1) ...
,+viewing_polar_a_S_ ...
,+viewing_azimu_b_S_ ...
,+zeros(n_S,1) ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,[] ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
a_k_p_from_S_ = reshape(a_from_S_C1M1_k_p_qk__./max(1e-12,a_from_S_C2M0_k_p_qk__),[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_S_',a_k_p_from_S_);
%%%%%%%%;
qbp_eps = 1e-6;
[ ...
 a_k_Y_from_M_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_true_M_ ...
,euler_azimu_b_true_M_ ...
,euler_gamma_z_true_M_ ...
,+image_delta_x_true_M_ ...
,+image_delta_y_true_M_ ...
);
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_from_M_',a_k_Y_from_M_);
%%%%%%%%;
qbp_eps = 1e-6;
[ ...
 a_k_Y_from_S_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,zeros(n_S,1) ...
);
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_from_S_',a_k_Y_from_S_);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

return;





