function ...
[ ...
 parameter ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
test_slice_vs_volume_integral_helper_eig_helper_figure_5( ...
 parameter ...
,n_x_u_pack ...
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
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,k_p_r_max ...
,X_weight_r_ ...
,pm_a_k_Y_use_yk__ ...
,from_pm_UX_kn__...
,n_viewing_S_use ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,euler_polar_a_tavg_ ...
,euler_azimu_b_tavg_ ...
);

str_thisfunction = 'test_slice_vs_volume_integral_helper_eig_helper_figure_5';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if (~isfield(parameter,'dir_ssnll')); parameter.dir_ssnll = ''; end; %<-- parameter_bookmark. ;
dir_ssnll = parameter.dir_ssnll;
if (~isfield(parameter,'str_fname_nopath_prefix')); parameter.str_fname_nopath_prefix = ''; end; %<-- parameter_bookmark. ;
str_fname_nopath_prefix = parameter.str_fname_nopath_prefix;
if (~isfield(parameter,'str_fname_nopath_sub_prefix')); parameter.str_fname_nopath_sub_prefix = ''; end; %<-- parameter_bookmark. ;
str_fname_nopath_sub_prefix = parameter.str_fname_nopath_sub_prefix;
if (~isfield(parameter,'tmp_fname_sub_mat')); parameter.tmp_fname_sub_mat = ''; end; %<-- parameter_bookmark. ;
tmp_fname_sub_mat = parameter.tmp_fname_sub_mat;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_implicit_dtau')); parameter.flag_implicit_dtau = 1; end; %<-- parameter_bookmark. ;
flag_implicit_dtau = parameter.flag_implicit_dtau;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
flag_replot = parameter.flag_replot; nf=0;
if (~isfield(parameter,'index_lambda')); parameter.index_lambda = 0; end; %<-- parameter_bookmark. ;
index_lambda = parameter.index_lambda;
if (~isfield(parameter,'lambda_mid')); parameter.lambda_mid = 0; end; %<-- parameter_bookmark. ;
lambda_mid = parameter.lambda_mid;
if (~isfield(parameter,'lambda_dif')); parameter.lambda_dif = 0; end; %<-- parameter_bookmark. ;
lambda_dif = parameter.lambda_dif;
if (~isfield(parameter,'lambda_lsq')); parameter.lambda_lsq = 0; end; %<-- parameter_bookmark. ;
lambda_lsq = parameter.lambda_lsq;
if (~isfield(parameter,'RR_bar')); parameter.RR_bar = 0; end; %<-- parameter_bookmark. ;
RR_bar = parameter.RR_bar;
if (~isfield(parameter,'val_zoom_use')); parameter.val_zoom_use = 0; end; %<-- parameter_bookmark. ;
val_zoom_use = parameter.val_zoom_use;
if (~isfield(parameter,'prct_use')); parameter.prct_use = 0; end; %<-- parameter_bookmark. ;
prct_use = parameter.prct_use;
if (~isfield(parameter,'vlim_g_max')); parameter.vlim_g_max = 0; end; %<-- parameter_bookmark. ;
vlim_g_max = parameter.vlim_g_max;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

str_dir_mat = sprintf('%s_mat',dir_ssnll);
if ~exist(str_dir_mat,'dir'); disp(sprintf(' %% mkdir %s',str_dir_mat)); mkdir(str_dir_mat); end;
str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if ~exist(str_dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg)); mkdir(str_dir_jpg); end;
str_dir_jpg_stripped = sprintf('%s_jpg_stripped',dir_ssnll);
if ~exist(str_dir_jpg_stripped,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg_stripped)); mkdir(str_dir_jpg_stripped); end;
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
if ~exist(str_dir_sub_mat,'dir'); disp(sprintf(' %% mkdir %s',str_dir_sub_mat)); mkdir(str_dir_sub_mat); end;
str_dir_sub_jpg = sprintf('%s/dir_%s',str_dir_jpg,str_fname_nopath_prefix);
if ~exist(str_dir_sub_jpg,'dir'); disp(sprintf(' %% mkdir %s',str_dir_sub_jpg)); mkdir(str_dir_sub_jpg); end;

x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;

k_c_0_all_ = k_p_r_all_.*cos(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_1_all_ = k_p_r_all_.*sin(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_2_all_ = k_p_r_all_.*cos(k_p_polar_a_all_);

l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = load(tmp_fname_sub_mat);
%%%%;
pm_n_k_p_r = tmp_.pm_n_k_p_r;
pm_k_p_r_ = tmp_.pm_k_p_r_;
pm_k_p_r_max = tmp_.pm_k_p_r_max;
pm_l_max_ = tmp_.pm_l_max_;
pm_weight_3d_k_p_r_ = tmp_.pm_weight_3d_k_p_r_;
n_M_use = tmp_.n_M_use;
weight_imagecount_M_use_ = tmp_.weight_imagecount_M_use_;
euler_polar_a_M_use_ = tmp_.euler_polar_a_M_use_;
euler_azimu_b_M_use_ = tmp_.euler_azimu_b_M_use_;
euler_gamma_z_M_use_ = tmp_.euler_gamma_z_M_use_;
n_M_imp = tmp_.n_M_imp;
weight_imagecount_M_imp_ = tmp_.weight_imagecount_M_imp_;
scaling_volumetric = tmp_.scaling_volumetric;
pm_weight_3d_riesz_k_p_r_ = tmp_.pm_weight_3d_riesz_k_p_r_;
pm_weight_3d_riesz_weight_imagecount_ykabc_ = tmp_.pm_weight_3d_riesz_weight_imagecount_ykabc_;
pm_v_tilde_eig_ykabc_ = tmp_.pm_v_tilde_eig_ykabc_;
pm_v_eig_ykabc_ = tmp_.pm_v_eig_ykabc_;
%%%%;

pm_l_max_max = max(pm_l_max_);
pm_n_lm_ = (1+pm_l_max_).^2;
pm_n_lm_sum = sum(pm_n_lm_);
pm_n_lm_csum_ = cumsum([0;pm_n_lm_]);

%%%%%%%%;
% reconstruct tmp_a_x_u_reco_frompm_. ;
%%%%%%%%;
a_k_Y_reco_frompm_yk_ = zeros(n_lm_sum,1);
for pm_nk_p_r=0:pm_n_k_p_r-1;
for nk_p_r=0:n_k_p_r-1;
if (flag_verbose>1); disp(sprintf(' %% adding pm_nk_p_r %d/%d nk_p_r %d/%d',pm_nk_p_r,pm_n_k_p_r,nk_p_r,n_k_p_r)); end;
tmp_l_max = l_max_(1+nk_p_r);
pm_tmp_l_max = pm_l_max_(1+pm_nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
pm_tmp_n_lm = (pm_tmp_l_max+1).^2;
pm_tmp_index_ = pm_n_lm_csum_(1+pm_nk_p_r) + (0:pm_tmp_n_lm-1);
a_k_Y_reco_frompm_yk_(1+tmp_index_) = a_k_Y_reco_frompm_yk_(1+tmp_index_) + from_pm_UX_kn__(1+nk_p_r,1+pm_nk_p_r)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_Y_use_yk__(1:tmp_n_lm,1+pm_nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%;
tmp_a_k_p_reco_frompm_ = zeros(n_k_all,1);
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
tmp_a_k_p_reco_frompm_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
0*flag_verbose ...
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
,a_k_Y_reco_frompm_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_a_k_p_reco_frompm_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
tmp_a_x_u_reco_frompm_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_a_k_p_reco_frompm_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_a_x_u_reco_frompm_ time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
% visualize pm_v_eig_dvol_yk_;
%%%%%%%%;
[pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_polar_a_M_use_,pm_v_tilde_eig_azimu_b_M_use_,pm_v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(pm_n_k_p_r,pm_l_max_,n_M_imp,pm_v_tilde_eig_ykabc_);
[pm_v_eig_dvol_yk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(pm_n_k_p_r,pm_l_max_,n_M_imp,pm_v_eig_ykabc_);
pm_v_k_Y_use_yk_ = pm_v_eig_dvol_yk_;
pm_v_k_Y_use_yk__ = local_yk__from_yk_(pm_n_k_p_r,pm_l_max_,pm_v_k_Y_use_yk_);
v_k_Y_reco_yk_ = zeros(n_lm_sum,1);
pm_n_UX_rank = pm_n_k_p_r;
for nUX_rank=0:pm_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
v_k_Y_reco_yk_(1+tmp_index_) = v_k_Y_reco_yk_(1+tmp_index_) + from_pm_UX_kn__(1+nk_p_r,1+nUX_rank)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_v_k_Y_use_yk__(1:tmp_n_lm,1+nUX_rank);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
%%%%;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
v_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
0*flag_verbose ...
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
,v_k_Y_reco_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% v_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
v_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,v_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: v_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
tmp_a_x_u_reco_ = tmp_a_x_u_reco_frompm_;

%%%%%%%%;
% estimate magnitude of perturbation. ;
%%%%%%%%;
if ~exist('pm_a_k_Y_use_yk_','var'); pm_a_k_Y_use_yk_ = local_yk_from_yk__(pm_n_k_p_r,pm_l_max_,pm_a_k_Y_use_yk__); end;
[~,~,tmp_v_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_v_eig_dvol_yk_);
[~,~,tmp_a_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_a_k_Y_use_yk_);
if (flag_verbose>0); disp(sprintf(' %% tmp_v_std/tmp_a_std %0.6f',tmp_v_std/tmp_a_std)); end;
[pm_a_k_Y_use_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_a_k_Y_use_yk_,pm_a_k_Y_use_yk_));
[pm_v_tilde_eig_dvol_lr] = sqrt(local_weightless_f_bar_dot_g_(pm_n_k_p_r,pm_l_max_,0,pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_dvol_yk_));
[pm_v_eig_dvol_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_v_eig_dvol_yk_,pm_v_eig_dvol_yk_));
if (flag_verbose>0); disp(sprintf(' %% pm_v_tilde_eig_dvol_lr %0.6f pm_v_eig_dvol_lr %0.6f',pm_v_tilde_eig_dvol_lr,pm_v_eig_dvol_lr)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_v_eig_dvol_lr/pm_a_k_Y_use_lr %0.6f',pm_v_eig_dvol_lr/pm_a_k_Y_use_lr)); end;
if (flag_verbose>0); disp(sprintf(' %% pm_a_k_Y_use_lr/pm_v_eig_dvol_lr %0.6f',pm_a_k_Y_use_lr/pm_v_eig_dvol_lr)); end;
% So we can imagine multiplying pm_v_eig_dvol_yk_ by (tmp_a_std/tmp_v_std) to produce something with equivalent norm to pm_a_k_Y_use_yk_. ;
% This then multiplies the value of H by (tmp_a_std/tmp_v_std)^2 --> lambda_mid*(tmp_a_std/tmp_v_std)^2 (in image-driven units). ;
% Note that sum(weight_imagecount_M_use_) --> 1.0, ;
% so our likelihood takes the form of:
% ssnll(a+delta*dvol) - ssnll(a) = N_image * (0.5 * H_scaled * delta^2) ; %<-- here a and dvol have comparable l2-norm. ;
% nll(a+delta*dvol) - nll(a) = sigma^-2 * N_image * (0.5 * H_scaled * delta^2) ; %<-- delta could be thought of as dimensionless. ;
%%%%%%%%;
tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_/max(1e-12,fnorm(tmp_a_x_u_reco_));
v_x_u_reco_fnrm_ = v_x_u_reco_/max(1e-12,fnorm(v_x_u_reco_));
tmp_H_scale = min(lambda_dif,lambda_lsq)*(tmp_a_std/tmp_v_std)^2;
tmp_N_image = 1024;
tmp_inv_sigma_squared = inv(RR_bar/n_x_u_pack^2);
tmp_log20 = 3.0;
tmp_dvol_mag = real(sqrt(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / tmp_N_image ));
tmp_dvol_mag = real( sqrt(RR_bar)/n_x_u_pack * sqrt(2.0 * tmp_log20 / (tmp_H_scale) / tmp_N_image ));
if (flag_verbose>0); disp(sprintf(' %% min(lambda_dif,lambda_lsq) %+0.6f (tmp_a_std/tmp_v_std)=(%+0.6f/%+0.6f) %+0.6f tmp_H_scale %+0.6f tmp_dvol_mag %+0.6f',min(lambda_dif,lambda_lsq),tmp_a_std,tmp_v_std,(tmp_a_std/tmp_v_std),tmp_H_scale,tmp_dvol_mag)); end;
if tmp_dvol_mag>=0.125;
tmp_N_image = ceil(real(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / 0.125^2 ));
tmp_N_image = max(tmp_N_image,1000*ceil(tmp_N_image/1000));
tmp_dvol_mag = real(sqrt(tmp_log20 / (tmp_inv_sigma_squared * 0.5 * tmp_H_scale) / tmp_N_image ));
if (flag_verbose>0); disp(sprintf(' %% adjusting: tmp_dvol_mag %+0.6f tmp_N_image %d',tmp_dvol_mag,tmp_N_image)); end;
end;%if tmp_dvol_mag>=0.125;
%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGF',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
vmax = max(abs(v_x_u_reco_),[],'all');
amax = max(abs(tmp_a_x_u_reco_),[],'all');
dvol_ = 3e-1*[-1:+1]; p_col = numel(dvol_);
prct_ = [98.5,99.0]; p_row = numel(prct_);
np=0;
for prow=0:p_row-1;
prct = prct_(1+prow);
pcol=0;
for dvol=0:numel(dvol_)-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
dvol = dvol_(1+dvol);
tmp_tmp_a_x_u_reco_ = tmp_a_x_u_reco_/max(1e-12,amax) + dvol*v_x_u_reco_/max(1e-12,vmax);
isosurface_f_x_u_1(struct('percent_threshold_',[prct]),tmp_tmp_a_x_u_reco_);
title(sprintf('prct %5.2f dvol %0.2f',prct,dvol));
pcol=pcol+1;
end;%for dvol=0:numel(dvol_)-1;
end;%for prow=0:p_row-1;
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGG',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
prct_ = [98.5,99.0]; p_row = numel(prct_);
p_col = 2; np=0;
for prow=0:p_row-1;
prct = prct_(1+prow);
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
title(sprintf('prct [%5.2f,%5.2f]:',[100-prct,prct]));
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_2(struct('percent_threshold_',prct),tmp_a_x_u_reco_,v_x_u_reco_);
title(sprintf('prct [%5.2f]:',[prct]));
end;%for prow=0:p_row-1;
set(gcf,'Position',1+[0,0,1024*1.5,1024*1.5]);
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGH',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
figbig;
subplot(1,1,1);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/2*sqrt(0.5),'flag_2d_vs_3d',1,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); ylabel('polar_a','Interpreter','none');
axisnotick;
title('dtau_euler_','Interpreter','none');
str_sgtitle = sprintf('%s index_lambda %d lambda %+0.2f[%+0.2f] = exp(%+0.2f)[exp(%+0.2f)]',str_fname_nopath_prefix,index_lambda,lambda_dif,lambda_lsq,log(abs(lambda_dif)),log(abs(lambda_lsq)));
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGI',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
figbig;
for np=0:1;
subplot(1,2,1+np);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',0,'flag_normalize',1) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d; axisnotick3d;
if np==0; view(0,+90); title('dtau_euler_ view(0,+90)','Interpreter','none'); end;
if np==1; view(0,-90); title('dtau_euler_ view(0,-90)','Interpreter','none'); end;
end;%for np=0:1;
%%%%;
sgtitle(str_sgtitle,'Interpreter','none');
%sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%;
end;%if flag_disp;
 
flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
dvol_ = tmp_dvol_mag*[-1:0.5:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
dvol = dvol_(1+ndvol);
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 16;
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_tmp_a_x_u_reco_fnrm_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),v_x_u_reco_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dvol %+.04f',prct,dvol),'Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
set(gca,'FontSize',fontsize_use);
str_sgtitle = sprintf('%s index_lambda %d lambda_dif %0.4f tmp_dvol_mag %0.6f',str_fname_nopath_prefix,index_lambda,lambda_dif,tmp_dvol_mag);
sgtitle(str_sgtitle,'Interpreter','none');
str_subplot = sprintf('p%.4d_dvol%d',100*prct,ndvol);
fname_fig_pre = sprintf('%s/%s_FIGK_%s',str_dir_sub_jpg,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGK_%s_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%for ndvol=0:n_dvol-1;
%%%%%%%%;
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
[ ...
 tmp_h_raw_ab_ ...
 tmp_h_w3d_ab_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,euler_polar_a_tavg_(:) ...
,euler_azimu_b_tavg_(:) ...
,[] ...
,[] ...
,-1 ...
);
factor_imagecount_M_use_ = tmp_h_w3d_ab_;
tmp_f = sum(factor_imagecount_M_use_.*viewing_weight_S_use_);
factor_imagecount_M_use_ = factor_imagecount_M_use_./max(1e-12,tmp_f);
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
subplot(1,2,[1]);
%%%%;
flim_ = 4.0*[0,2.0/(4*pi)];
flag_2d_vs_3d=0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_S_use_ ... 
,viewing_azimu_b_S_use_ ... 
,factor_imagecount_M_use_ ...
,flim_ ... 
,colormap_beach ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
title('$\mu(\tau)$','Interpreter','latex');
axisnotick3d; axis equal; axis vis3d;
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,2,[2]);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
val_zoom = val_zoom_use;
isosurface_f_x_u_1( ...
 struct('vval_',[vval]) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_a_x_u_reco_fnrm_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGL',str_dir_sub_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGL_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;  
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 4.0; compass_r_base_use = factor_amplify*1.0/(2*pi)/2*sqrt(0.5);
dpolar_a = min(diff(unique(tmp_.euler_polar_a_M_use_)));
dazimu_b = min(diff(unique(tmp_.euler_azimu_b_M_use_)));
dsurf = min(dpolar_a,dazimu_b)*2*pi*0.5/8;
compass_r_base_use = dsurf;
%%%%;
for np=0;%for np=0:1;
subplot(1,2,[1]);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',compass_r_base_use,'flag_2d_vs_3d',0,'flag_normalize',0) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlabel(''); ylabel(''); zlabel(''); axis vis3d; axisnotick3d; axis([-1,+1,-1,+1,-1,+1]);
if np==0; view(0,+90); title('$\Delta\tau$: azimuthal $+\pi/2$','Interpreter','latex'); end;
if np==1; view(0,-90); title('$\Delta\tau$: azimuthal $-\pi/2$','Interpreter','latex'); end;
set(gca,'FontSize',fontsize_use);
end;%for np=0:1;
%%%%;
subplot(1,2,[2]);
prct = prct_use;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGM',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGM_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;  
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_.tmp_dtau_euler_polar_a_M_use_,tmp_.tmp_dtau_euler_azimu_b_M_use_,tmp_.tmp_dtau_euler_gamma_z_M_use_].^2,tmp_.weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_.tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_.tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_.tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
%%%%;
%{
deleteme = struct();
deleteme.n_M_use = tmp_.n_M_use;
deleteme.euler_polar_a_M_use_ = tmp_.euler_polar_a_M_use_;
deleteme.euler_azimu_b_M_use_ = tmp_.euler_azimu_b_M_use_;
deleteme.dtau_euler_nrm_polar_a_M_ = tmp_dtau_euler_nrm_polar_a_M_;
deleteme.dtau_euler_nrm_azimu_b_M_ = tmp_dtau_euler_nrm_azimu_b_M_;
deleteme.dtau_euler_nrm_gamma_z_M_ = tmp_dtau_euler_nrm_gamma_z_M_;
save('/data/rangan/dir_cryoem/dir_rangan_playroom/deleteme.mat','deleteme');
%%%%;
load('/data/rangan/dir_cryoem/dir_rangan_playroom/deleteme.mat','deleteme');
fontsize_use = 16;
hold on;
factor_amplify = 4.0; compass_r_base_use = factor_amplify*1.0/(2*pi)/2*sqrt(0.5);
dpolar_a = min(diff(unique(deleteme.euler_polar_a_M_use_)));
dazimu_b = min(diff(unique(deleteme.euler_azimu_b_M_use_)));
dsurf = min(dpolar_a,dazimu_b)*2*pi*0.5/8;
compass_r_base_use = dsurf;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',compass_r_base_use,'flag_2d_vs_3d',1,'flag_normalize',0) ...
,deleteme.n_M_use ...
,deleteme.euler_polar_a_M_use_ ...
,deleteme.euler_azimu_b_M_use_ ...
,deleteme.dtau_euler_nrm_polar_a_M_ ...
,deleteme.dtau_euler_nrm_azimu_b_M_ ...
,deleteme.dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimuthal','Interpreter','latex');
set(gca,'XTick',0:pi/4:2*pi,'XTickLabel',[]);
ylim([0,1*pi]); ylabel('polar','Interpreter','latex');
set(gca,'YTick',0:pi/4:1*pi,'YTickLabel',[]);
grid on;
title('$\Delta\tau$: equatorial','Interpreter','latex');
set(gca,'FontSize',fontsize_use);
return;
%}
%%%%;
factor_amplify = 4.0; compass_r_base_use = factor_amplify*1.0/(2*pi)/2*sqrt(0.5);
dpolar_a = min(diff(unique(tmp_.euler_polar_a_M_use_)));
dazimu_b = min(diff(unique(tmp_.euler_azimu_b_M_use_)));
dsurf = min(dpolar_a,dazimu_b)*2*pi*0.5/8;
compass_r_base_use = dsurf;
subplot(1,3,[1,2]);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',compass_r_base_use,'flag_2d_vs_3d',1,'flag_normalize',0) ...
,tmp_.n_M_use ...
,tmp_.euler_polar_a_M_use_ ...
,tmp_.euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimuthal','Interpreter','latex');
set(gca,'XTick',0:pi/4:2*pi,'XTickLabel',[]);
ylim([0,1*pi]); ylabel('polar','Interpreter','latex');
set(gca,'YTick',0:pi/4:1*pi,'YTickLabel',[]);
grid on;
title('$\Delta\tau$: equatorial','Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,[3]);
prct = prct_use;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.75,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGN',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGN_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); end;
if ndvol==1; subplot(1,2,[2]); end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),tmp_tmp_a_x_u_reco_fnrm_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),v_x_u_reco_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$',dvol),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGO',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGO_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
interp_k = 1;
%%%%;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
tmp_tmp_a_x_u_reco_mini_ = min(tmp_a_x_u_reco_fnrm_ + min(dvol_)*v_x_u_reco_fnrm_,tmp_a_x_u_reco_fnrm_ + max(dvol_)*v_x_u_reco_fnrm_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); c_use__ = flipud(colormap_pm); end;
if ndvol==1; subplot(1,2,[2]); c_use__ =        colormap_pm ; end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
hold on;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'c_use__',0.95*[1,1,1],'flag_projection',0,'flag_collapse',0,'flag_boxgrid',1,'v_alpha',1.000) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom*1.00),interp3(reshape(tmp_tmp_a_x_u_reco_mini_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'c_use__',     c_use__,'flag_projection',0,'flag_collapse',0,'flag_boxgrid',1,'v_alpha',0.500) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom*0.99),interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
hold off;
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$',dvol),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGP',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGP_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
interp_k = 1;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = [0,vlim_g_max];
val_zoom = val_zoom_use;
dvol = max(dvol_);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
figure(1+nf);nf=nf+1;clf;figsml;
[~,d_v_] = ...
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_,'flag_plot',0) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
close(gcf);
vlim_g_ = [0,prctile(d_v_,90)];
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
val_zoom = val_zoom_use;
for ndvol=0:n_dvol-1;
%if ndvol==0; subplot(2,5,[ 2, 3, 7, 8]); end;
%if ndvol==1; subplot(2,5,[ 4, 5, 9,10]); end;
if ndvol==0; subplot(1,2,[1]); end;
if ndvol==1; subplot(1,2,[2]); end;
dvol = dvol_(1+ndvol);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$, $J$%: $%d$',dvol,tmp_N_image),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGQ',str_dir_sub_jpg,str_fname_nopath_sub_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGQ_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
interp_k = 1;
dvol_ = tmp_dvol_mag*[-1:0.5:+1]; n_dvol = numel(dvol_);
prct_ = [prct_use]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = [0,vlim_g_max];
val_zoom = val_zoom_use;
dvol = max(dvol_);
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
figure(1+nf);nf=nf+1;clf;figsml;
[~,d_v_] = ...
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_,'flag_plot',0) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
close(gcf);
vlim_g_ = [0,prctile(d_v_,90)];
%%%%%%%%;
for ndvol=0:n_dvol-1;
dvol = dvol_(1+ndvol);
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 16;
interp_k = 2;
tmp_tmp_a_x_u_reco_fnrm_ = tmp_a_x_u_reco_fnrm_ + dvol*v_x_u_reco_fnrm_;
isosurface_f_x_u_3( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom) ...
			 ,interp3(reshape(tmp_a_x_u_reco_fnrm_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline') ...
			 ) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dvol %+.04f',prct,dvol),'Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*2,2*(768-128)]);
set(gca,'FontSize',fontsize_use);
str_sgtitle = sprintf('%s index_lambda %d lambda_dif %0.4f tmp_dvol_mag %0.6f tmp_N_image %d',str_fname_nopath_prefix,index_lambda,lambda_dif,tmp_dvol_mag,tmp_N_image);
sgtitle(str_sgtitle,'Interpreter','none');
str_subplot = sprintf('p%.4d_dvol%d',100*prct,ndvol);
fname_fig_pre = sprintf('%s/%s_FIGK2_%s',str_dir_sub_jpg,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGK2_%s_stripped',str_dir_jpg_stripped,str_fname_nopath_sub_prefix,str_subplot);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%for ndvol=0:n_dvol-1;
%%%%%%%%;
end;%if flag_disp;

%%%%;
clear tmp_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

