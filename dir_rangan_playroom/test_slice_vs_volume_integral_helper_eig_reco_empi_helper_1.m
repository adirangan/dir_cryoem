%%%%%%%%;
% First collect/collate data. ;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
str_dir_mat = sprintf('%s_mat',dir_ssnll);
if ~exist(str_dir_mat,'dir'); disp(sprintf(' %% mkdir %s',str_dir_mat)); mkdir(str_dir_mat); end;
str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if ~exist(str_dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg)); mkdir(str_dir_jpg); end;
str_dir_jpg_stripped = sprintf('%s_jpg_stripped',dir_ssnll);
if ~exist(str_dir_jpg_stripped,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg_stripped)); mkdir(str_dir_jpg_stripped); end;
%%%%%%%%;

%%%%%%%%;
flag_implicit_dtau = 1;
viewing_weight_M_use_ = ones(n_M_use,1);
factor_imagecount_M_use_ = 1/max(1e-12,sum(viewing_weight_M_use_));
weight_imagecount_M_use_ = viewing_weight_M_use_ .* factor_imagecount_M_use_ ;
nlsigma_dist = 0; n_lsigma_dist = 1; lsigma_dist = NaN;
str_infix = sprintf('p_reco_empi');
if (flag_verbose>0); disp(sprintf(' %% %s',str_infix)); end;
lanczos_n_iteration_max = 32;
ddssnll_mid_q2d_i_ = zeros(1,lanczos_n_iteration_max);
ddssnll_dif_q2d_i_ = zeros(1,lanczos_n_iteration_max);
ddssnll_lsq_q2d_i_ = zeros(1,lanczos_n_iteration_max);
nfound = 0;
str_fname_nopath_prefix = sprintf('eig_i1_from_image_%s_%s',str_tolerance_pm,str_infix);
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
for index_lambda=0:lanczos_n_iteration_max-1;
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,lanczos_n_iteration_max-1,index_lambda);
tmp_fname_sub_mat = sprintf('%s/%s.mat',str_dir_sub_mat,str_fname_nopath_sub_prefix);
if ~exist(tmp_fname_sub_mat,'file');
disp(sprintf(' %% Warning, %s not found',tmp_fname_sub_mat));
end;%if ~exist(tmp_fname_sub_mat,'file');
if  exist(tmp_fname_sub_mat,'file');
tmp_ = load(tmp_fname_sub_mat,'tmp_ddssnll_mid_q2d','tmp_ddssnll_dif_q2d','tmp_ddssnll_lsq_q2d');
disp(sprintf(' %% %% %s: mid %+16.4f dif %+16.4f lsq %+16.4f',str_fname_nopath_sub_prefix,tmp_.tmp_ddssnll_mid_q2d,tmp_.tmp_ddssnll_dif_q2d,tmp_.tmp_ddssnll_lsq_q2d));
ddssnll_mid_q2d_i_(1,1+index_lambda) = tmp_.tmp_ddssnll_mid_q2d;
ddssnll_dif_q2d_i_(1,1+index_lambda) = tmp_.tmp_ddssnll_dif_q2d;
ddssnll_lsq_q2d_i_(1,1+index_lambda) = tmp_.tmp_ddssnll_lsq_q2d;
nfound  = nfound+1;
clear tmp_;
end;%if  exist(tmp_fname_sub_mat,'file');
end;%for index_lambda=0:lanczos_n_iteration_max-1;
disp(sprintf(' %% '));
if (flag_verbose>0); disp(sprintf(' %% found %d/%d',nfound,lanczos_n_iteration_max)); end;
%%%%%%%%;

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
a_k_Y_reco_frompm_yk_(1+tmp_index_) = a_k_Y_reco_frompm_yk_(1+tmp_index_) + UX_kn__(1+nk_p_r,1+pm_nk_p_r)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_Y_reco_empi_yk__(1:tmp_n_lm,1+pm_nk_p_r);
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

lambda_cut = 3.5;
%%%%%%%%;
% Run through each single example and try to visualize the alignment perturbation. ;
npick_ = { ...
  ,[15] ...
  ,[16] ...
  ,[17] ...
};
n_pick = numel(npick_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
npick = npick_{1+npick};
index_lambda = npick(1+0);

%%%%%%%%;
lambda_mid = real(ddssnll_mid_q2d_i_(1,1+index_lambda));
lambda_dif = real(ddssnll_dif_q2d_i_(1,1+index_lambda));
lambda_lsq = real(ddssnll_lsq_q2d_i_(1,1+index_lambda));
lambda_mid_min = min(real(ddssnll_mid_q2d_i_(1,:)));
lambda_dif_min = min(real(ddssnll_dif_q2d_i_(1,:)));
lambda_lsq_min = min(real(ddssnll_lsq_q2d_i_(1,:)));
flag_ismin = (lambda_dif==lambda_dif_min) | (lambda_lsq==lambda_lsq_min);
disp(sprintf(' %% index_lambda %d lambda_mid +%0.6f lambda_dif +%0.6f lambda_lsq +%0.6f',index_lambda,lambda_mid,lambda_dif,lambda_lsq));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%if (flag_ismin | (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nlsigma_dist = 0; n_lsigma_dist = 1; lsigma_dist = NaN;
str_infix = sprintf('p_reco_empi');
if (flag_verbose>0); disp(sprintf(' %% %s',str_infix)); end;
str_dir_sub_mat = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
str_dir_sub_jpg = sprintf('%s/dir_%s',str_dir_jpg,str_fname_nopath_prefix);
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,lanczos_n_iteration_max-1,index_lambda);
tmp_fname_sub_mat = sprintf('%s/%s.mat',str_dir_sub_mat,str_fname_nopath_sub_prefix);
if ~exist(tmp_fname_sub_mat,'file');
disp(sprintf(' %% Warning, %s not found',tmp_fname_sub_mat));
end;%if ~exist(tmp_fname_sub_mat,'file');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = load(tmp_fname_sub_mat);
%%%%;
from_pm_UX_kn__ = UX_kn__;
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
[~,~,tmp_v_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_v_eig_dvol_yk_);
[~,~,tmp_a_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_a_k_Y_reco_empi_yk_);
if (flag_verbose>0); disp(sprintf(' %% tmp_v_std/tmp_a_std %0.6f',tmp_v_std/tmp_a_std)); end;
[pm_a_k_Y_reco_empi_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_a_k_Y_reco_empi_yk_,pm_a_k_Y_reco_empi_yk_));
[pm_v_tilde_eig_dvol_lr] = sqrt(local_weightless_f_bar_dot_g_(pm_n_k_p_r,pm_l_max_,0,pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_dvol_yk_));
[pm_v_eig_dvol_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_v_eig_dvol_yk_,pm_v_eig_dvol_yk_));
if (flag_verbose>0); disp(sprintf(' %% pm_v_eig_dvol_lr/pm_a_k_Y_reco_empi_lr %0.6f',pm_v_eig_dvol_lr/pm_a_k_Y_reco_empi_lr)); end;
% So we can imagine multiplying pm_v_eig_dvol_yk_ by (tmp_a_std/tmp_v_std) to produce something with equivalent norm to pm_a_k_Y_reco_empi_yk_. ;
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
if (flag_verbose>0); disp(sprintf(' %% tmp_dvol_mag %+0.6f',tmp_dvol_mag)); end;
%%%%%%%%;
 
flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
dvol_ = tmp_dvol_mag*[-1:0.5:+1]; n_dvol = numel(dvol_);
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = 2.0;
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

flag_disp=0;
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
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
val_zoom = 2.0;
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

flag_disp=0;
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
factor_amplify = 4.0;
%%%%;
for np=0;%for np=0:1;
subplot(1,2,[1]);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',factor_amplify*1.0/(2*pi)/4,'flag_2d_vs_3d',0,'flag_normalize',1) ...
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
prct = 98.5;
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

flag_disp=0;
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
factor_amplify = 4.0;
%%%%;
subplot(1,3,[1,2]);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',factor_amplify*1.0/(2*pi)/2*sqrt(0.5),'flag_2d_vs_3d',1,'flag_normalize',1) ...
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
prct = 98.5;
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

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
dvol_ = tmp_dvol_mag*[-1:2.0:+1]; n_dvol = numel(dvol_);
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = 2.0;
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
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(tmp_a_x_u_reco_fnrm_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = 2.0;
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

%%%%;
clear tmp_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if  exist(tmp_fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;%if (lambda_dif<=lambda_cut & lambda_lsq<=lambda_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npick=0:n_pick-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
