% clear; test_pm_trpv1_2;

%%%%%%%%;
% Load one of the runs which exhibits divergence. ;
%%%%%%%%;
n_iteration = 32;
str_strategy = '';
delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [1,2,4,6,8,10,12,14,16];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
n_strategy = 4;
str_strategy_ = cell(n_strategy,1);
for flag_alternate_MS_vs_SM = [0:1];
for flag_local_exclusion = [0:1];
str_strategy = ''; 
if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
if (flag_local_exclusion); str_strategy = sprintf('%se1',str_strategy); end;
nstrategy = flag_local_exclusion + 2*flag_alternate_MS_vs_SM;
str_strategy_{1+nstrategy} = str_strategy;
end;%for flag_local_exclusion = [0:1];
end;%for flag_alternate_MS_vs_SM = [0:1];
%%%%%%%%;
XA_align_a_CTF_avg_UX_Y_ = zeros(n_iteration,1);
XA_align_a_k_Y_ = zeros(n_iteration,1);
XB_align_a_CTF_avg_UX_Y_ = zeros(n_iteration,1);
XB_align_a_k_Y_ = zeros(n_iteration,1);
XC_align_a_CTF_avg_UX_Y_ = zeros(2*n_iteration,1);
XC_align_a_k_Y_ = zeros(2*n_iteration,1);
%%%%%%%%;
nfound=0; n_all=0;
ndat_rseed=1;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
ndelta_r_max_factor=0;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
ndat_n_UX_rank=7;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
flag_alternate_MS_vs_SM=0;%for flag_alternate_MS_vs_SM = [0:1];
flag_local_exclusion=0;%for flag_local_exclusion = [0:1];
str_strategy = ''; 
if (flag_alternate_MS_vs_SM); str_strategy = sprintf('%sa1',str_strategy); end;
if (flag_local_exclusion); str_strategy = sprintf('%se1',str_strategy); end;
nstrategy = flag_local_exclusion + 2*flag_alternate_MS_vs_SM;
%%%%;
XA_fname_pre = sprintf('%s_mat/X_2d_Memp_d1_%st%.4dn%.2dr%d',dir_pm,str_strategy,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_XA_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',XA_fname_pre);
if ( exist(fname_XA_align_a_CTF_avg_UX_Y_mat,'file'));
tmp_ = load(fname_XA_align_a_CTF_avg_UX_Y_mat);
XA_align_a_CTF_avg_UX_Y_(:) = tmp_.X_best_;
XC_align_a_CTF_avg_UX_Y_(0*n_iteration + [1:n_iteration]) = tmp_.X_best_;
end;%if ( exist(fname_XA_align_a_CTF_avg_UX_Y_mat,'file'));
fname_XA_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',XA_fname_pre);
if ( exist(fname_XA_align_a_k_Y_mat,'file'));
tmp_ = load(fname_XA_align_a_k_Y_mat);
XA_align_a_k_Y_(:) = tmp_.X_best_;
XC_align_a_k_Y_(0*n_iteration + [1:n_iteration]) = tmp_.X_best_;
end;%if ( exist(fname_XA_align_a_k_Y_mat,'file'));
XB_fname_pre = sprintf('%s_mat/X_2d_xcor_d0_%st%.4dn%.2dr%d',dir_pm,str_strategy,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_XB_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',XB_fname_pre);
if ( exist(fname_XB_align_a_CTF_avg_UX_Y_mat,'file'));
tmp_ = load(fname_XB_align_a_CTF_avg_UX_Y_mat);
XB_align_a_CTF_avg_UX_Y_(:) = tmp_.X_best_;
XC_align_a_CTF_avg_UX_Y_(1*n_iteration + [1:n_iteration]) = tmp_.X_best_;
end;%if ( exist(fname_XB_align_a_CTF_avg_UX_Y_mat,'file'));
fname_XB_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',XB_fname_pre);
if ( exist(fname_XB_align_a_k_Y_mat,'file'));
tmp_ = load(fname_XB_align_a_k_Y_mat);
XB_align_a_k_Y_(:) = tmp_.X_best_;
XC_align_a_k_Y_(1*n_iteration + [1:n_iteration]) = tmp_.X_best_;
end;%if ( exist(fname_XB_align_a_k_Y_mat,'file'));
%%%%;
%end;%for flag_local_exclusion = [0:1];
%end;%for flag_alternate_MS_vs_SM = [0:1];
%end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
%end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;

%%%%%%%%;
% Now examine the reconstructed molecule. ;
%%%%%%%%;
tmp_ = load(fname_XA_align_a_k_Y_mat);
[~,niteration] = max(tmp_.X_best_); niteration = niteration-1;
%%%%%%%%;
fname_XA_mat = sprintf('%s.mat',XA_fname_pre);
tmp_ = load(fname_XA_mat);
tmp_n_order = 5; tmp_n_M = n_M;
tmp_euler_polar_a_ = +tmp_.euler_polar_a__(:,1+niteration);
tmp_euler_azimu_b_ = +tmp_.euler_azimu_b__(:,1+niteration);
tmp_euler_gamma_z_ = +tmp_.euler_gamma_z__(:,1+niteration);
tmp_image_delta_x_ = +tmp_.image_delta_x_acc__(:,1+niteration) + tmp_.image_delta_x_upd__(:,1+niteration);
tmp_image_delta_y_ = +tmp_.image_delta_y_acc__(:,1+niteration) + tmp_.image_delta_y_upd__(:,1+niteration);
%%%%%%%%;
tmp_t = tic;
d_k_Y_reco_ = ...
cg_lsq_4( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> d_k_Y_reco_ time %0.2fs',tmp_t));
[ ...
 tmp_X_best_reco ...
,tmp_X_best_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,d_k_Y_reco_ ...
);
disp(sprintf(' %% tmp_X_best_reco %0.3f flag_flip %d',tmp_X_best_reco,tmp_X_best_flag_flip));
%%%%%%%%;
tmp_t = tic;
[d_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,d_k_Y_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% d_k_Y_reco_ --> d_k_p_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; 
d_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,d_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: d_x_u_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
fname_XA_a_k_Y_mat = sprintf('%s_a_k_Y_.mat',XA_fname_pre);
tmp_ = load(fname_XA_a_k_Y_mat);
e_k_Y_reco_ = tmp_.a_k_Y_reco_;
tmp_t = tic;
[e_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,e_k_Y_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% e_k_Y_reco_ --> e_k_p_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; n_X_u = x_u_res^3;
e_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,e_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: e_x_u_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
fname_fig = sprintf('%s_visualize_divergence_FIGA',XA_fname_pre);
if ( 1 | ~exist(sprintf('%s.jpg',fname_fig),'file') );
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1);clf;figmed;
tmp_p_ = [90,95,99];
tmp_p_ = 98.5;
subplot(1,3,1);
isosurface_f_x_u_0(reshape(real(c_x_u_reco_),x_u_res,x_u_res,x_u_res),tmp_p_);
title(sprintf('c_x_u_: corr %0.4f',X_best_reco),'Interpreter','none');
subplot(1,3,2);
isosurface_f_x_u_0(reshape(real(d_x_u_reco_),x_u_res,x_u_res,x_u_res),tmp_p_);
title(sprintf('d_x_u_: (ni%d) corr %0.4f',niteration,XA_align_a_k_Y_(1+niteration)),'Interpreter','none');
subplot(1,3,3);
isosurface_f_x_u_0(reshape(real(e_x_u_reco_),x_u_res,x_u_res,x_u_res),tmp_p_);
title(sprintf('e_x_u_: corr %0.4f',XA_align_a_k_Y_(end)),'Interpreter','none');
sgtitle(sprintf('%s',XA_fname_pre),'Interpreter','none');
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if ( 1 | ~exist(sprintf('%s.jpg',fname_fig),'file') );
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
