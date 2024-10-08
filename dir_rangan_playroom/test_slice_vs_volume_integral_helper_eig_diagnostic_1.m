str_type = 'equa_band';
str_thissubfunction = 'test_slice_vs_volume_integral_helper_eig_diagnostic_1';

scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
%%%%%%%%;
% Set up weights. ;
%%%%%%%%;
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
weight_3d_riesz_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
weight_3d_riesz_yk_(1+tmp_index_) = weight_3d_riesz_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
n_ykabc = n_lm_sum + n_M_use*3;
weight_3d_riesz_ykabc_ = cat(1,weight_3d_riesz_yk_/scaling_volumetric,ones(3*n_M_use,1));
numerator_root_weight_3d_riesz_ykabc_ = reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]);
denomator_root_weight_3d_riesz_ykabc_ = 1./max(1e-12,reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]));
%%%%%%%%;

if strcmp(str_type,'reco_empi');
%%%%%%%%;
% visualize the empirical fit. ;
%%%%%%%%;
a_k_Y_use_yk_ = a_k_Y_reco_empi_yk_;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_p_use_ ...
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
,a_k_Y_use_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_use_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_use_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_use_ time %0.2fs',tmp_t));
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[05:10:95];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),a_x_u_use_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle('a_k_Y_reco_empi_yk_','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
end;%if strcmp(str_type,'reco_empi');

%%%%%%%%;
if strcmp(str_type,'reco_empi');
str_M = sprintf('M%.4d',n_M_use); str_fname_pre = sprintf('eig_from_reco_empi_%s',str_M);
str_title = sprintf('reco_empi_%s',str_M);
fname_mat = sprintf('%s_mat/%s.mat',dir_ssnll,str_fname_pre);
end;%if strcmp(str_type,'reco_empi');
%%%%%%%%;
if strcmp(str_type,'equa_band');
str_fname_pre = sprintf('eig_from_synth_equa_band');
str_title = sprintf('equa_band');
fname_mat = sprintf('%s_mat/%s.mat',dir_ssnll,str_fname_pre);
end;%if strcmp(str_type,'equa_band');
%%%%%%%%;
if strcmp(str_type,'polar_cap');
str_fname_pre = sprintf('eig_from_synth_polar_cap');
str_title = sprintf('polar_cap');
fname_mat = sprintf('%s_mat/%s.mat',dir_ssnll,str_fname_pre);
end;%if strcmp(str_type,'polar_cap');
%%%%%%%%;
if (~exist(fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found',fname_mat));
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found',fname_mat));
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
tmp_ = load(fname_mat);
tmp_.n_iteration = numel(tmp_.alph_tilde_i_); tmp_.T_tilde__ = real(spdiags([circshift(tmp_.beta_tilde_i_,-1),tmp_.alph_tilde_i_,tmp_.beta_tilde_i_],[-1,0,+1],tmp_.n_iteration,tmp_.n_iteration));
tmp_.lambda_xi__ = -Inf*ones(tmp_.n_iteration,tmp_.n_iteration);
for niteration=0:tmp_.n_iteration-1;
tmp_.T_tilde_sub__ = tmp_.T_tilde__(1:1+niteration,1:1+niteration);
tmp_.lambda_sub_ = eigs(tmp_.T_tilde_sub__,[],1+niteration);
tmp_.lambda_xi__(1:1+niteration,1+niteration) = sort(tmp_.lambda_sub_,'ascend');
end;%for niteration=0:tmp_.n_iteration-1;
tmp_.S_x_ = sort(eigs(tmp_.T_tilde__,[],tmp_.n_iteration),'ascend');
n_iteration = tmp_.n_iteration;
S_x_ = tmp_.S_x_;
S_x_min = min(tmp_.S_x_);
S_x_max = max(tmp_.S_x_);
lambda_xi__ = tmp_.lambda_xi__;
T_tilde__ = tmp_.T_tilde__;
%%%%;
if isfield(tmp_,'n_M_use'); tmp_n_M_use = tmp_.n_M_use; end;
if isfield(tmp_,'n_synth_M'); tmp_n_M_use = tmp_.n_synth_M; end;
vv_n4__ = zeros(tmp_.n_iteration,4);
ee_ns4___ = zeros(tmp_.n_iteration,tmp_.n_iteration,4);
for niteration=0:tmp_.n_iteration-1;
v_tilde_ykabc_ = tmp_.v_tilde_ykabci__(:,1+niteration);
[v_tilde_dvol_yk_,v_tilde_polar_a_M_use_,v_tilde_azimu_b_M_use_,v_tilde_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_weightless_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_ykabc_,v_tilde_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_vv)); end;
vv_n4__(1+niteration,:) = [tmp_vv_dvol;tmp_vv_a;tmp_vv_b;tmp_vv_c];
tmp_.T_tilde_sub__ = tmp_.T_tilde__(1:1+niteration,1:1+niteration);
[tmp_.TV_tilde_sub__,tmp_.lambda_sub__] = eigs(tmp_.T_tilde_sub__,[],1+niteration);
tmp_.lambda_sub_ = diag(tmp_.lambda_sub__);
[lambda_srt_,ij_srt_] = sort(tmp_.lambda_sub_,'ascend');
for index_lambda=0:1+niteration-1;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_tilde_eig_ = tmp_.TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = tmp_.v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ee)); end;
ee_ns4___(1+niteration,1+index_lambda,:) = [tmp_ee_dvol;tmp_ee_a;tmp_ee_b;tmp_ee_c];
end;%for index_lambda=0:1+niteration-1;
end;%for niteration=0:tmp_.n_iteration-1;
%%%%%%%%;

%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/%s_FIGB',dir_ssnll,str_fname_pre);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);fig80s();
p_row = 1; p_col = 5; np=0;
fontsize_use = 12;
ilim_ = [-0.125,+1.125];
%%%%;
for pcol=0:4-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(squeeze(ee_ns4___(:,:,1+pcol)),ilim_);
xlabel('index'); ylabel('iteration'); %axisnotick;
title(sprintf('ee_ns4___(:,:,1+%d)',pcol),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
end;%for pcol=0:4-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(vv_n4__,ilim_);
xlabel('1-4'); ylabel('iteration'); %axisnotick;
title('vv_n4__','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/%s_FIGA',dir_ssnll,str_fname_pre);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figmed;fig81s;
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
%%%%;
hold on;
plot(repmat([0;n_iteration],[1,n_iteration]),repmat(reshape(S_x_,[1,n_iteration]),[2,1]),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_sml);
ni_xi__ = repmat([1:n_iteration],[n_iteration,1]);
tmp_index_ = efind(isfinite(lambda_xi__));
plot(ni_xi__(1+tmp_index_),lambda_xi__(1+tmp_index_),'r.','MarkerSize',markersize_use);
hold off;
xlabel('iteration'); ylabel('sigma');
xlim([0,1+n_iteration]);
ylim([S_x_min-0.25,S_x_max+0.25]);
title(str_title,'Interpreter','none');
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
niteration = n_iteration-1;
index_lambda = 0; %<-- 0 is lowest, niteration is highest. ;
tmp_.T_tilde_sub__ = tmp_.T_tilde__(1:1+niteration,1:1+niteration);
%%%%%%%%;
% use T_tilde_sub__ to estimate minimum eigenvector. ;
%%%%%%%%;
[tmp_.TV_tilde_sub__,tmp_.lambda_sub__] = eigs(tmp_.T_tilde_sub__,[],1+niteration);
tmp_.lambda_sub_ = diag(tmp_.lambda_sub__);
[lambda_srt_,ij_srt_] = sort(tmp_.lambda_sub_,'ascend');
%[~,index_lambda] = min(abs(lambda_srt_)); index_lambda = index_lambda - 1; %<-- pick min(abs). ;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_tilde_eig_ = tmp_.TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = tmp_.v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;
%%%%%%%%;

[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;
v_eig_ykabc_ = bsxfun(@times,denomator_root_weight_3d_riesz_ykabc_,v_tilde_eig_ykabc_);
[v_eig_dvol_yk_,v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_n_M_use,v_eig_ykabc_,v_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;

%%%%%%%%;
% visualize v_eig_dvol_yk_;
%%%%%%%%;
v_k_Y_use_yk_ = v_eig_dvol_yk_;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 v_k_p_use_ ...
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
,v_k_Y_use_yk_ ...
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
v_x_u_use_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,v_k_p_use_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: v_x_u_use_ time %0.2fs',tmp_t));
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[05:10:95];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),v_x_u_use_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle(sprintf('v_eig_dvol_yk_: %s niteration %d lambda %0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_title,niteration,lambda_use,log(lambda_use),index_lambda,niteration,str_ee),'Interpreter','none');
end;%if flag_disp;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
v_x_u_use___ = reshape(v_x_u_use_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
subplot(1,3,1);
imagesc(squeeze(real(v_x_u_use___(:,:,end/2))));
axis image; axisnotick; colorbar;
subplot(1,3,2);
imagesc(squeeze(real(v_x_u_use___(:,end/2,:))));
axis image; axisnotick; colorbar;
subplot(1,3,3);
imagesc(squeeze(real(v_x_u_use___(end/2,:,:))));
axis image; axisnotick; colorbar;
sgtitle(sprintf('%s niteration %d lambda %0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_title,niteration,lambda_use,log(lambda_use),index_lambda,niteration,str_ee),'Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if strcmp(str_type,'reco_empi') & strcmp(platform,'access1'); %<-- Now evaluate ddssnll. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
parameter_ddssnll = struct('type','ddssnll');
parameter_ddssnll.flag_verbose = flag_verbose;
parameter_ddssnll.flag_check = 1;
parameter_ddssnll.flag_disp = 1;
parameter_ddssnll.flag_kernel_qpro_d0 = 1;
parameter_ddssnll.flag_kernel_qpro_d1 = 1;
parameter_ddssnll.kernel_qpro_polar_a_pole_north=KAPPA_pole_north_double;
parameter_ddssnll.kernel_qpro_polar_a_pole_south=KAPPA_pole_south_double;
parameter_ddssnll.kernel_qpro_qref_k_eq_d_double=KAPPA_qref_k_eq_d_double;
parameter_ddssnll.kernel_qpro_l_max_use = l_max;
S_k_p_q2d_wkS__ = S_k_p_wkS__;
%%%%;
euler_polar_a_use_M_ = euler_polar_a_empi_(1:n_M_use);
euler_azimu_b_use_M_ = euler_azimu_b_empi_(1:n_M_use);
euler_gamma_z_use_M_ = euler_gamma_z_empi_(1:n_M_use);
image_delta_x_use_M_ = image_delta_x_empi_(1:n_M_use);
image_delta_y_use_M_ = image_delta_y_empi_(1:n_M_use);
M_use_k_p_wkM__ = M_k_p_wkM__(:,1:n_M_use);
for nM_use=0:n_M_use-1;
M_use_k_p_wkM__(:,1+nM_use) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_use_k_p_wkM__(:,1+nM_use),+image_delta_x_use_M_(1+nM_use),+image_delta_y_use_M_(1+nM_use));
end;%for nM_use=0:n_M-1;
%%%%;

if ~exist('KAPPA','var'); KAPPA=[]; end;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlma____','var'); d0W_betazeta_mlma____=[]; end;
if ~exist('d1W_betazeta_mlma____','var'); d1W_betazeta_mlma____=[]; end;
if ~exist('d2W_betazeta_mlma____','var'); d2W_betazeta_mlma____=[]; end;

[ ...
 ~ ...
,Hvt_ykabc_ ...
,Hv_q3d_k_Y_quad_yk_ ...
,Hv_q3d_k_Y_quad_yk__ ...
,Hv_q3d_k_p_quad_ ...
,Ht_q2d_M3__ ...
,a_restore_C2M0_k_Y_lmk_ ...
,a_restore_C2M0_k_p_quad_ ...
] = ...
ddssnll_1( ...
 parameter_ddssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_use_yk_ ...
,[] ...
,v_eig_dvol_yk_ ...
,[] ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M_use ...
,M_use_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_use_M_ ...
,euler_azimu_b_use_M_ ...
,euler_gamma_z_use_M_ ...
,v_eig_polar_a_M_use_ ...
,v_eig_azimu_b_M_use_ ...
,v_eig_gamma_z_M_use_ ...
,KAPPA ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
);

%%%%%%%%;
% check difference quotient. ;
%%%%%%%%;
ddssnll_mid_q2d = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M_use,Hvt_ykabc_,v_eig_ykabc_);
%%%%;
dvt = 1e-3;
index_neta_from_nM_ =[];
eta_k_p_wke__=[];
parameter_ssnll = parameter_ddssnll;
%%%%;
[ ...
 ~ ...
,~ ...
,ssnll_mid_q2d ...
] = ...
ssnll_from_a_k_Y_12( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_use_yk_ + 0.0*dvt*v_eig_dvol_yk_...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M_use ...
,M_use_k_p_wkM__ ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_use_M_ + 0.0*dvt*real(v_eig_polar_a_M_use_) ...
,euler_azimu_b_use_M_ + 0.0*dvt*real(v_eig_azimu_b_M_use_) ...
,euler_gamma_z_use_M_ + 0.0*dvt*real(v_eig_gamma_z_M_use_) ...
,[] ...
,[] ...
,[] ...
);
%%%%;
[ ...
 ~ ...
,~ ...
,ssnll_pos_q2d ...
] = ...
ssnll_from_a_k_Y_12( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_use_yk_ + 1.0*dvt*v_eig_dvol_yk_...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M_use ...
,M_use_k_p_wkM__ ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_use_M_ + 1.0*dvt*real(v_eig_polar_a_M_use_ ) ...
,euler_azimu_b_use_M_ + 1.0*dvt*real(v_eig_azimu_b_M_use_) ...
,euler_gamma_z_use_M_ + 1.0*dvt*real(v_eig_gamma_z_M_use_) ...
,[] ...
,[] ...
,[] ...
);
%%%%;
[ ...
 ~ ...
,~ ...
,ssnll_neg_q2d ...
] = ...
ssnll_from_a_k_Y_12( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_use_yk_ - 1.0*dvt*v_eig_dvol_yk_...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M_use ...
,M_use_k_p_wkM__ ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_use_M_ - 1.0*dvt*real(v_eig_polar_a_M_use_ ) ...
,euler_azimu_b_use_M_ - 1.0*dvt*real(v_eig_azimu_b_M_use_) ...
,euler_gamma_z_use_M_ - 1.0*dvt*real(v_eig_gamma_z_M_use_) ...
,[] ...
,[] ...
,[] ...
);
%%%%;
ddssnll_dif_q2d = (ssnll_pos_q2d - 2*ssnll_mid_q2d + ssnll_neg_q2d)/max(1e-12,dvt^2);
if (flag_verbose>0); disp(sprintf(' %% ddssnll_dif_q2d vs ddssnll_mid_q2d: %0.16f',fnorm(ddssnll_dif_q2d -ddssnll_mid_q2d)/max(1e-12,fnorm(ddssnll_dif_q2d)))); end;
%%%%%%%%;

%%%%%%%%;
% display the rayleigh-quotients. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% lambda_use = %0.6f',lambda_use)); end;
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_n_M_use,v_eig_ykabc_,v_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;
[tmp_hv,tmp_hv_dvol,tmp_hv_a,tmp_hv_b,tmp_hv_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_n_M_use,Hvt_ykabc_,v_eig_ykabc_);
str_hv = sprintf('tmp_hv %0.2f,tmp_hv_dvol %0.2f,tmp_hv_a %0.2f,tmp_hv_b %0.2f,tmp_hv_c %0.2f',tmp_hv,tmp_hv_dvol,tmp_hv_a,tmp_hv_b,tmp_hv_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_hv)); end;
if (flag_verbose>0); disp(sprintf(' %% corr(Hvt_ykabc_,v_eig_ykabc_) = %0.6f',corr(Hvt_ykabc_,v_eig_ykabc_))); end;
if (flag_verbose>0); disp(sprintf(' %% corr(Hv_q3d_k_Y_quad_yk_,v_eig_dvol_yk_) = %0.6f',corr(Hv_q3d_k_Y_quad_yk_,v_eig_dvol_yk_))); end;
v_eig_abc_M3__ = [v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_];
if (flag_verbose>0); disp(sprintf(' %% corr(Ht_q2d_M3__,v_eig_abc_M3__) = %0.6f',corr(Ht_q2d_M3__(:),v_eig_abc_M3__(:)))); end;
%%%%%%%%;

%%%%%%%%;
% visualize a_restore_C2M0_k_p_quad_;
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_restore_C2M0_x_u_use_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_restore_C2M0_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_restore_C2M0_x_u_use_ time %0.2fs',tmp_t));
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[05:10:95];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),a_restore_C2M0_x_u_use_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle(sprintf('v_eig_dvol_yk_: %s niteration %d lambda %0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_title,niteration,lambda_use,log(lambda_use),index_lambda,niteration,str_ee),'Interpreter','none');
end;%if flag_disp;
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
a_restore_C2M0_x_u_use___ = reshape(a_restore_C2M0_x_u_use_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
subplot(1,3,1);
imagesc(squeeze(real(a_restore_C2M0_x_u_use___(:,:,end/2))));
axis image; axisnotick; colorbar;
subplot(1,3,2);
imagesc(squeeze(real(a_restore_C2M0_x_u_use___(:,end/2,:))));
axis image; axisnotick; colorbar;
subplot(1,3,3);
imagesc(squeeze(real(a_restore_C2M0_x_u_use___(end/2,:,:))));
axis image; axisnotick; colorbar;
sgtitle(sprintf('%s a_restore_C2M0_x_u_use_',str_title));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if strcmp(str_type,'reco_empi') & strcmp(platform,'access1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
