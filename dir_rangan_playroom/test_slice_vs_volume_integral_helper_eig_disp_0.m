n_trial = n_T + 3;
fname_mat_t_ = cell(n_trial,1);
str_title_t_ = cell(n_trial,1);
for ntrial=0:n_trial-1;
if ntrial==n_T+0;
str_M = sprintf('M%.4d',n_M_use); str_eig = sprintf('eig_from_reco_empi_%s',str_M);
fname_mat = sprintf('%s_mat/%s.mat',dir_ssnll,str_eig);
str_title = sprintf('reco_empi_%s',str_M);
end;%if;
if ntrial==n_T+1;
fname_mat = sprintf('%s_mat/eig_from_synth_polar_cap.mat',dir_ssnll);
str_title = sprintf('synth_polar_cap');
end;%if;
if ntrial==n_T+2;
fname_mat = sprintf('%s_mat/eig_from_synth_equa_band.mat',dir_ssnll);
str_title = sprintf('synth_equa_band');
end;%if;
if ntrial<=n_T-1;
nT = ntrial;
str_T = sprintf('T%s',num2str(floor(100*nT*T_MAX/max(n_T-1)),'%.3d'));
fname_mat = sprintf('%s_mat/eig_from_synth_%s.mat',dir_ssnll,str_T);
str_title = sprintf('synth_%s',str_T);
end;%if;
fname_mat_t_{1+ntrial} = fname_mat;
str_title_t_{1+ntrial} = str_title;
clear fname_mat str_title;
end;%for ntrial=0:n_trial-1;

n_iteration_t_ = zeros(n_trial,1);
S_x_t__ = cell(n_trial,1);
lambda_txi___ = cell(n_trial,1);
S_x_min = +Inf; S_x_max = -Inf;
vv_tns4____ = cell(n_trial,1);
for ntrial=0:n_trial-1;
fname_mat = fname_mat_t_{1+ntrial};
if (~exist(fname_mat,'file'));
disp(sprintf(' %% Warning, %s not found',fname_mat));
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found',fname_mat));
tmp_ = load(fname_mat);
tmp_.n_iteration = numel(tmp_.alph_i_); tmp_.T__ = real(spdiags([circshift(tmp_.beta_i_,-1),tmp_.alph_i_,tmp_.beta_i_],[-1,0,+1],tmp_.n_iteration,tmp_.n_iteration));
tmp_.lambda_xi__ = -Inf*ones(tmp_.n_iteration,tmp_.n_iteration);
for niteration=0:tmp_.n_iteration-1;
tmp_.T_sub__ = tmp_.T__(1:1+niteration,1:1+niteration);
tmp_.lambda_sub_ = eigs(tmp_.T_sub__,[],1+niteration);
tmp_.lambda_xi__(1:1+niteration,1+niteration) = sort(tmp_.lambda_sub_,'ascend');
end;%for niteration=0:tmp_.n_iteration-1;
tmp_.S_x_ = sort(eigs(tmp_.T__,[],tmp_.n_iteration),'ascend');
n_iteration_t_(1+ntrial) = tmp_.n_iteration;
S_x_t__{1+ntrial} = tmp_.S_x_;
S_x_min = min(S_x_min,min(tmp_.S_x_));
S_x_max = max(S_x_max,max(tmp_.S_x_));
lambda_txi___{1+ntrial} = tmp_.lambda_xi__;
%%%%;
vv_ns4___ = zeros(tmp_.n_iteration,tmp_.n_iteration,4);
for niteration=0:tmp_.n_iteration-1;
tmp_.T_sub__ = tmp_.T__(1:1+niteration,1:1+niteration);
[tmp_.TV_sub__,tmp_.lambda_sub__] = eigs(tmp_.T_sub__,[],1+niteration);
tmp_.lambda_sub_ = diag(tmp_.lambda_sub__);
[lambda_srt_,ij_srt_] = sort(tmp_.lambda_sub_,'ascend');
for index_lambda=0:1+niteration-1;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_min_ = tmp_.TV_sub__(:,ij_use);
v_min_ykabc_ = tmp_.v_ykabci__(:,1:1+niteration)*TV_min_;
if isfield(tmp_,'n_M_use'); tmp_n_M_use = tmp_.n_M_use; end;
if isfield(tmp_,'n_synth_M'); tmp_n_M_use = tmp_.n_synth_M; end;
[v_min_dvol_yk_,v_min_polar_a_synth_M_,v_min_azimu_b_synth_M_,v_min_gamma_z_synth_M_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_n_M_use,v_min_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_n_M_use,v_min_ykabc_,v_min_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_vv)); end;
vv_ns4___(1+niteration,1+index_lambda,:) = [tmp_vv_dvol;tmp_vv_a;tmp_vv_b;tmp_vv_c];
end;%for index_lambda=0:1+niteration-1;
end;%for niteration=0:tmp_.n_iteration-1;
vv_tns4____{1+ntrial} = vv_ns4___;
%%%%;
clear vv_ns4___ tmp_;
end;%if ( exist(fname_mat,'file'));
end;%for ntrial=0:n_trial-1;

if flag_disp;
fname_fig_pre = sprintf('%s_jpg/eig_from_x_FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;fig81s;
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
p_row = 3; p_col = ceil(n_trial/p_row); np=0;
%%%%;
for ntrial=0:n_trial-1;
str_title = str_title_t_(1+ntrial);
n_iteration = n_iteration_t_(1+ntrial);
S_x_ = S_x_t__{1+ntrial};
lambda_xi__ = lambda_txi___{1+ntrial};
subplot(p_row,p_col,1+np);np=np+1;
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
end;%for ntrial=0:n_trial-1;
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
end;%if flag_disp;

tmp_logmin = exp(-10);
lS_x_min = log(max(tmp_logmin,S_x_min));
lS_x_max = log(max(tmp_logmin,S_x_max));
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/eig_from_x_FIGB',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;fig81s;
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
p_row = 3; p_col = ceil(n_trial/p_row); np=0;
%%%%;
for ntrial=0:n_trial-1;
str_title = str_title_t_(1+ntrial);
n_iteration = n_iteration_t_(1+ntrial);
S_x_ = S_x_t__{1+ntrial};
lambda_xi__ = lambda_txi___{1+ntrial};
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(repmat([0;n_iteration],[1,n_iteration]),repmat(reshape(log(max(tmp_logmin,S_x_)),[1,n_iteration]),[2,1]),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_sml);
ni_xi__ = repmat([1:n_iteration],[n_iteration,1]);
tmp_index_ = efind(isfinite(lambda_xi__));
plot(ni_xi__(1+tmp_index_),log(max(tmp_logmin,lambda_xi__(1+tmp_index_))),'r.','MarkerSize',markersize_use);
hold off;
xlabel('iteration'); ylabel('log(sigma)');
xlim([0,1+n_iteration]);
ylim([lS_x_min-0.25,lS_x_max+0.25]);
title(str_title,'Interpreter','none');
end;%for ntrial=0:n_trial-1;
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
end;%if flag_disp;

%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nT=0; %<-- 0 is the lowest, n_T-1 is the nighest. ;
index_lambda = 3; %<-- 0 is the lowest, niteration-1 is the highest. ;
str_T = sprintf('T%s',num2str(floor(100*nT*T_MAX/max(n_T-1)),'%.3d'));
fname_mat = sprintf('%s_mat/eig_from_synth_%s.mat',dir_ssnll,str_T);
%%%%%%%%%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct tmp_.weight_3d_riesz')); end;
%%%%%%%%;
tmp_.weight_3d_riesz_k_p_r_ = tmp_.weight_3d_k_p_r_;
tmp_.weight_3d_riesz_k_all_ = tmp_.weight_3d_k_all_;
for nk_p_r=0:tmp_.n_k_p_r-1;
k_p_r = tmp_.k_p_r_(1+nk_p_r);
tmp_.weight_3d_k_p_r = tmp_.weight_3d_k_p_r_(1+nk_p_r);
tmp_.weight_2d_k_p_r = tmp_.weight_2d_k_p_r_(1+nk_p_r);
tmp_.weight_3d_riesz_k_p_r_(1+nk_p_r) = tmp_.weight_3d_k_p_r_(1+nk_p_r) * tmp_.weight_2d_k_p_r / max(1e-16,tmp_.weight_3d_k_p_r);
tmp_index_ = tmp_.n_k_all_csum_(1+nk_p_r):tmp_.n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(tmp_.weight_3d_k_all_(1+tmp_index_))/(tmp_.weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(tmp_.weight_3d_k_all_(1+tmp_index_))/(4*pi*tmp_.weight_3d_k_p_r))); end;
tmp_.weight_3d_riesz_k_all_(1+tmp_index_) = tmp_.weight_3d_k_all_(1+tmp_index_) * tmp_.weight_2d_k_p_r / max(1e-16,tmp_.weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(tmp_.weight_3d_riesz_k_all_(1+tmp_index_))/(tmp_.weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(tmp_.weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*tmp_.weight_2d_k_p_r))); end;
end;%for nk_p_r=0:tmp_.n_k_p_r-1;
%%%%%%%%;
tmp_.n_iteration = numel(tmp_.alph_i_);
tmp_.T__ = real(spdiags([circshift(tmp_.beta_i_,-1),tmp_.alph_i_,tmp_.beta_i_],[-1,0,+1],tmp_.n_iteration,tmp_.n_iteration));
tmp_.lambda_xi__ = -Inf*ones(tmp_.n_iteration,tmp_.n_iteration);
%%%%%%%%;
% pick niteration to define T_sub__. ;
%%%%%%%%;
niteration = tmp_.n_iteration-1; %<-- yes final niteration. ;
%niteration = tmp_.n_iteration-3; %<-- not final niteration. ;
tmp_.T_sub__ = tmp_.T__(1:1+niteration,1:1+niteration);
%%%%%%%%;
% use T_sub__ to estimate minimum eigenvector. ;
%%%%%%%%;
[tmp_.TV_sub__,tmp_.lambda_sub__] = eigs(tmp_.T_sub__,[],1+niteration);
tmp_.lambda_sub_ = diag(tmp_.lambda_sub__);
[lambda_srt_,ij_srt_] = sort(tmp_.lambda_sub_,'ascend');
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_min_ = tmp_.TV_sub__(:,ij_use);
v_min_ykabc_ = tmp_.v_ykabci__(:,1:1+niteration)*TV_min_;
[v_min_dvol_yk_,v_min_polar_a_synth_M_,v_min_azimu_b_synth_M_,v_min_gamma_z_synth_M_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_.n_synth_M,v_min_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_.n_synth_M,v_min_ykabc_,v_min_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_vv)); end;
%%%%%%%%;
% visualize volumetric perturbation as well as viewing-angle perturbation. ;
%%%%%%%%;
tmp_t = tic();
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 v_min_dvol_k_p_quad_ ...
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
,v_min_dvol_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% convert_spharm_to_k_p_4: time %0.6fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
v_min_dvol_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,v_min_dvol_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
subplot(1,3,1);
isosurface_f_x_u_1([],v_min_dvol_x_u_reco_); axisnotick3d; title('v_min_dvol_','Interpreter','none');
%%%%;
subplot(1,3,2);
v_min_nrm_polar_a_synth_M_ = real(v_min_polar_a_synth_M_);
v_min_nrm_azimu_b_synth_M_ = real(v_min_azimu_b_synth_M_);
v_min_nrm_gamma_z_synth_M_ = real(v_min_gamma_z_synth_M_);
v_min_nrm_fnorm = fnorm([v_min_nrm_polar_a_synth_M_,v_min_nrm_azimu_b_synth_M_,v_min_nrm_gamma_z_synth_M_]);
v_min_nrm_polar_a_synth_M_ = v_min_nrm_polar_a_synth_M_/max(1e-12,v_min_nrm_fnorm);
v_min_nrm_azimu_b_synth_M_ = v_min_nrm_azimu_b_synth_M_/max(1e-12,v_min_nrm_fnorm);
v_min_nrm_gamma_z_synth_M_ = v_min_nrm_gamma_z_synth_M_/max(1e-12,v_min_nrm_fnorm);
factor_amplify = 2.0;
plot_sphere_grid_0;
hold on;
sphere_post__0( ...
 struct('type','sphere_post','post_r_base',1.0/(2*pi)/8) ...
,tmp_.n_synth_M ...
,tmp_.euler_polar_a_synth_M_ ...
,tmp_.euler_azimu_b_synth_M_ ...
,factor_amplify*v_min_nrm_polar_a_synth_M_ ...
,factor_amplify*v_min_nrm_azimu_b_synth_M_ ...
,factor_amplify*v_min_nrm_gamma_z_synth_M_ ...
);
hold off;
axis equal; axis vis3d; axisnotick3d;
title('v_min_dtau_','Interpreter','none');
%%%%;
subplot(1,3,3);
plot_sphere_grid_0;
hold on;
sphere_compass__0( ...
 struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4) ...
,tmp_.n_synth_M ...
,tmp_.euler_polar_a_synth_M_ ...
,tmp_.euler_azimu_b_synth_M_ ...
,v_min_nrm_polar_a_synth_M_ ...
,v_min_nrm_azimu_b_synth_M_ ...
,v_min_nrm_gamma_z_synth_M_ ...
);
hold off;
axis equal; axis vis3d; axisnotick3d;
title('v_min_dtau_','Interpreter','none');
%%%%;
sgtitle(sprintf('%s niteration %d lambda %0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_T,niteration,lambda_use,log(lambda_use),index_lambda,niteration,str_vv),'Interpreter','none');
%%%%%%%%%%%%%%%%;
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%}


%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nT=0; %<-- 0 is the lowest, n_T-1 is the nighest. ;
index_lambda = 0; %<-- 0 is the lowest, niteration-1 is the highest. ;
str_T = sprintf('T%s',num2str(floor(100*nT*T_MAX/max(n_T-1)),'%.3d'));
fname_mat = sprintf('%s_mat/eig_from_synth_%s.mat',dir_ssnll,str_T);
%%%%%%%%%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct tmp_.weight_3d_riesz')); end;
%%%%%%%%;
tmp_.weight_3d_riesz_k_p_r_ = tmp_.weight_3d_k_p_r_;
tmp_.weight_3d_riesz_k_all_ = tmp_.weight_3d_k_all_;
for nk_p_r=0:tmp_.n_k_p_r-1;
k_p_r = tmp_.k_p_r_(1+nk_p_r);
tmp_.weight_3d_k_p_r = tmp_.weight_3d_k_p_r_(1+nk_p_r);
tmp_.weight_2d_k_p_r = tmp_.weight_2d_k_p_r_(1+nk_p_r);
tmp_.weight_3d_riesz_k_p_r_(1+nk_p_r) = tmp_.weight_3d_k_p_r_(1+nk_p_r) * tmp_.weight_2d_k_p_r / max(1e-16,tmp_.weight_3d_k_p_r);
tmp_index_ = tmp_.n_k_all_csum_(1+nk_p_r):tmp_.n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(tmp_.weight_3d_k_all_(1+tmp_index_))/(tmp_.weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(tmp_.weight_3d_k_all_(1+tmp_index_))/(4*pi*tmp_.weight_3d_k_p_r))); end;
tmp_.weight_3d_riesz_k_all_(1+tmp_index_) = tmp_.weight_3d_k_all_(1+tmp_index_) * tmp_.weight_2d_k_p_r / max(1e-16,tmp_.weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(tmp_.weight_3d_riesz_k_all_(1+tmp_index_))/(tmp_.weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(tmp_.weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*tmp_.weight_2d_k_p_r))); end;
end;%for nk_p_r=0:tmp_.n_k_p_r-1;
%%%%%%%%;
tmp_.n_iteration = numel(tmp_.alph_i_);
tmp_.T__ = real(spdiags([circshift(tmp_.beta_i_,-1),tmp_.alph_i_,tmp_.beta_i_],[-1,0,+1],tmp_.n_iteration,tmp_.n_iteration));
tmp_.lambda_xi__ = -Inf*ones(tmp_.n_iteration,tmp_.n_iteration);
%%%%%%%%;
% pick niteration to define T_sub__. ;
%%%%%%%%;
niteration = tmp_.n_iteration-1; %<-- yes final niteration. ;
%niteration = tmp_.n_iteration-3; %<-- not final niteration. ;
tmp_.T_sub__ = tmp_.T__(1:1+niteration,1:1+niteration);
%%%%%%%%;
% use T_sub__ to estimate minimum eigenvector. ;
%%%%%%%%;
[tmp_.TV_sub__,tmp_.lambda_sub__] = eigs(tmp_.T_sub__,[],1+niteration);
tmp_.lambda_sub_ = diag(tmp_.lambda_sub__);
[lambda_srt_,ij_srt_] = sort(tmp_.lambda_sub_,'ascend');
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
 TV_min_ = tmp_.TV_sub__(:,ij_use);
v_min_ykabc_ = tmp_.v_ykabci__(:,1:1+niteration)*TV_min_;
[v_min_dvol_yk_,v_min_polar_a_synth_M_,v_min_azimu_b_synth_M_,v_min_gamma_z_synth_M_] = local_yk_a_b_c_from_ykabc_(tmp_.n_k_p_r,tmp_.l_max_,tmp_.n_synth_M,v_min_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_f_bar_dot_g_(tmp_.n_k_p_r,tmp_.weight_3d_riesz_k_p_r_,tmp_.l_max_,tmp_.n_synth_M,v_min_ykabc_,v_min_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_vv)); end;
%%%%%%%%;
% visualize viewing-angle perturbation. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
hold on;
sphere_compass__0( ...
 struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',1) ...
,tmp_.n_synth_M ...
,tmp_.euler_polar_a_synth_M_ ...
,tmp_.euler_azimu_b_synth_M_ ...
,v_min_nrm_polar_a_synth_M_ ...
,v_min_nrm_azimu_b_synth_M_ ...
,v_min_nrm_gamma_z_synth_M_ ...
);
hold off;
axis([0,2*pi,0,1*pi]); axisnotick;
title('v_min_dtau_','Interpreter','none');
%%%%;
sgtitle(sprintf('%s niteration %d lambda %0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_T,niteration,lambda_use,log(lambda_use),index_lambda,niteration,str_vv),'Interpreter','none');
%%%%%%%%%%%%%%%%;
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%}
