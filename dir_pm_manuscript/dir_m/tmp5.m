
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% repeat using clusters. ;
% this time with idealized principal-modes. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat = sprintf('%s_mat/pm_fig_X_2d_cluster_d0_FIGI__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_t = tic();
delta_sigma_base = 0;
a_k_Y_base_yk_ = a_k_Y_quad_;
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
tmp_CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
tmp_CTF_k_p_r_xavg_kk__ = tmp_CTF_k_p_r_xavg_k_*transpose(tmp_CTF_k_p_r_xavg_k_);
[ ...
 X_2d_xavg_dx_kk__ ...
,X_2d_xavg_dx_weight_r_ ...
] = ...
principled_marching_cost_matrix_6( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_ ...
,[] ...
,[] ...
,a_k_Y_base_yk_ ...
,tmp_CTF_k_p_r_xavg_kk__ ...
,delta_sigma_base ...
);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
X_kkc___ = X_2d_xavg_dx_kkc___;
X_weight_rc__ = X_2d_xavg_dx_weight_rc__;
clear X_2d_xavg_dx_kkc__ X_2d_xavg_dx_weight_rc__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_kkc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
delta_r_max = 0*delta_sigma; svd_eps = tolerance_master; n_delta_v_requested = 32;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
%%%%%%%%;
tolerance_pm_ = 0.1.^[0.5:0.5:18]; n_tolerance_pm = numel(tolerance_pm_);
X_SMp__ = zeros(n_S,n_M,n_tolerance_pm);
euler_gamma_z_SMp__ = zeros(n_S,n_M,n_tolerance_pm);
X_t_ = zeros(n_tolerance_pm,1);
X_t0_ = zeros(n_tolerance_pm,1);
X_t1_ = zeros(n_tolerance_pm,1);
X_t2_ = zeros(n_tolerance_pm,1);
X_t3_ = zeros(n_tolerance_pm,1);
pm_n_UX_rank_cp__ = zeros(n_cluster,n_tolerance_pm);
%%%%%%%%;
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
X_kk__ = X_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
if isempty(pm_n_UX_rank); pm_n_UX_rank = 1; end;
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
pm_n_UX_rank_cp__(:,1+ntolerance_pm) = pm_n_UX_rank_c_;
%%%%%%%%;
pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
M_k_q_wkM__ = M_k_q__;
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M);
%%%%%%%%;
tmp_t = tic();
tmp_M_index_ = 0:n_M-1; tmp_n_M = n_M;
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
tmp_M_index_sub_ = intersect(tmp_M_index_,index_nM_from_ncluster_);
tmp_n_M_sub = numel(tmp_M_index_sub_);
if (tmp_n_M_sub> 0);
svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+tmp_M_index_sub_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M_sub,M_k_q_wkM__(:,1+tmp_M_index_sub_),pm_n_UX_rank,UX_kn__,X_weight_r_);
end;%if (tmp_n_M_sub> 0);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
a_k_Y_reco_yk_ = a_k_Y_quad_;
a_k_Y_reco_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_reco_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_reco_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
S_k_p_wkS__ = S_k_p__;
S_k_q_wkS__ = S_k_q__;
%%%%%%%%;
%%%%%%%%;
% storing all the principal-templates takes lots of memory. ;
%%%%%%%%;
CTF_UX_S_k_q_wnSc___ = [];
CTF_UX_S_l2_Sc__ = [];
%%%%%%%%;
% Use given volume to align principal-images. ;
% Groups principal-images by cluster. ;
% Calculates principal-templates associated with each cluster. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.tolerance_cluster = parameter.tolerance_cluster;
tmp_t = tic();
[ ...
 tmp_parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_cluster_wrap_SM__10( ...
 tmp_parameter ...
,FTK ...
,n_w_max ...
,n_k_p_r ...
,n_S ...
,S_k_q_wkS__ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,index_nCTF_from_nM_ ...
,n_M ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,CTF_k_p_r_xavg_kc__ ...
,CTF_UX_S_k_q_wnSc___ ...
,CTF_UX_S_l2_Sc__ ...
,index_ncluster_from_nM_ ...
,index_nM_from_ncluster__ ...
,n_index_nM_from_ncluster_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__10',tmp_t);
X_SMp___(:,:,1+ntolerance_pm) = X_SM__;
euler_gamma_z_SMp___(:,:,1+ntolerance_pm) = gamma_z_SM__;
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10'));
X_t_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_VUXM_nMwl____'));
X_t0_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_SVUXM_SMwl____'));
X_t1_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_SVUXM_lwSM____'));
X_t2_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_USESVUXM_dwSM____'));
X_t3_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
end;%for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','delta_sigma_base','n_UX_rank' ...
     ,'FTK' ...
     ,'tolerance_pm_','n_tolerance_pm' ...
     ,'pm_n_UX_rank_cp__' ...
     ,'X_SMp___','X_t_','X_t0_','X_t1_','X_t2_','X_t3_' ...
     ,'euler_gamma_z_SMp___' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% Calculate error in euler_gamma_z. ;
% mean absolute error (in terms of angle). ;
% as well as percentile (over SM). ;
%%%%%%%%;
prctile_ = [90,95,99]; n_prctile = numel(prctile_);
error_euler_gamma_z_avg_p_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_prc_pp__ = zeros(n_tolerance_pm,n_prctile);
for ntolerance_pm=0:n_tolerance_pm-1;
tmp_ = abs(periodize(euler_gamma_z_SMp___(:,:,end) - euler_gamma_z_SMp___(:,:,1+ntolerance_pm),-pi,pi));
error_euler_gamma_z_avg_p_(1+ntolerance_pm) = mean(tmp_,'all');
error_euler_gamma_z_prc_pp__(1+ntolerance_pm,:) = prctile(tmp_,prctile_,'all');
clear tmp_;
end%;for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
% fraction of image-template-pairs within discretization error in gamma_z. ;
%%%%%%%%;
error_euler_gamma_z_f1_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_f2_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_f4_ = zeros(n_tolerance_pm,1);
dgamma = 2*pi/n_w_max;
for ntolerance_pm=0:n_tolerance_pm-1;
tmp_ = abs(periodize(euler_gamma_z_SMp___(:,:,end) - euler_gamma_z_SMp___(:,:,1+ntolerance_pm),-pi,pi));
error_euler_gamma_z_f1_(1+ntolerance_pm) = mean( tmp_> 1*dgamma + 1e-12 , 'all' );
error_euler_gamma_z_f2_(1+ntolerance_pm) = mean( tmp_> 2*dgamma + 1e-12 , 'all' );
error_euler_gamma_z_f4_(1+ntolerance_pm) = mean( tmp_> 4*dgamma + 1e-12 , 'all' );
clear tmp_;
end;%for ntolerance_pm=0:n_tolerance_pm-1;

figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,512,768]);
subplot(2,1,1);
hold on;
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_avg_p_,'k.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_prc_pp__,'.-')
hold off;
subplot(2,1,2);
hold on;
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f1_,'r.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f2_,'r.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f4_,'r.-');
hold off;


pm_n_UX_rank_avg_p_ = mean(pm_n_UX_rank_cp__(1+index_ncluster_from_nM_,:),1);
tmp_ = reshape(X_SMp___,[n_S*n_M,n_tolerance_pm]);
corr_X_pp__ = corr(tmp_,tmp_);
corr_X_pe_ = corr_X_pp__(:,end);
%%%%%%%%;
prctile_ = [5,15,50,85,95]; n_prctile = numel(prctile_);
cap_pn__ = zeros(n_prctile,n_tolerance_pm);
X_SM_end__ = X_SMp___(:,:,end);
for nprctile=0:n_prctile-1;
cut_X_SM_end__ = zeros(size(X_SM_end__));
prc = prctile_(1+nprctile);
cut_X_SM_end__ = X_SM_end__ > repmat(prctile(X_SM_end__,prc),[n_S,1]);
for ntolerance_pm=0:n_tolerance_pm-1;
X_SM_sub__ = X_SMp___(:,:,1+ntolerance_pm);
cut_X_SM_sub__ = zeros(size(X_SM_sub__));
cut_X_SM_sub__ = X_SM_sub__ > repmat(prctile(X_SM_sub__,prc),[n_S,1]);
cap_pn__(1+nprctile,1+ntolerance_pm) = sum(cut_X_SM_end__.*cut_X_SM_sub__,'all');
end;%for ntolerance_pm=0:n_tolerance_pm-1;
end;%for nprctile=0:n_prctile-1;
cap_nrm_pn__ = cap_pn__./repmat(cap_pn__(:,end),[1,n_tolerance_pm]);
%%%%%%%%;

fname_mat = sprintf('%s_mat/pm_fig_X_2d_cluster_d0_timing_FIGI__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_t_mult_ = zeros(n_tolerance_pm,1);
tmp_t_ifft_ = zeros(n_tolerance_pm,1);
tmp_o_mult_ = zeros(n_tolerance_pm,1);
tmp_o_ifft_ = zeros(n_tolerance_pm,1);
for ntolerance_pm=0:n_tolerance_pm-1;%for ntolerance_pm=[0,5,n_tolerance_pm-1];
n_n = round(pm_n_UX_rank_avg_p_(1+ntolerance_pm));
tmp_n_w_max = n_w_max;
tmp_n_S = 1024; tmp_n_M = 1; 
%tmp_n_S = 1; tmp_n_M = 1024; 
tmp_n_S = min(tmp_n_S,n_S); tmp_n_M = min(tmp_n_M,n_M);
n_M_batch = ceil(n_M/tmp_n_M); n_S_batch = ceil(n_S/tmp_n_S);
tmp_t_mult = 0; tmp_t_ifft = 0;
tmp_o_mult = 0; tmp_o_ifft = 0;
for nS_batch=0:n_S_batch-1;for nM_batch=0:n_M_batch-1;
tmp_CTF_UX_S_k_q_Snw___ = randn(tmp_n_S,n_n,tmp_n_w_max) + i;
tmp_svd_VUXM_nMwl____ = randn(n_n,tmp_n_M,tmp_n_w_max) + i;
tmp_svd_SVUXM_SMwl____ = zeros(tmp_n_S,tmp_n_M,tmp_n_w_max);
tmp_t = tic();
for nw=0:tmp_n_w_max-1;
%tmp_svd_SVUXM_SMwl____(:,:,1+nw) = tmp_svd_SVUXM_SMwl____(:,:,1+nw) + tmp_CTF_UX_S_k_q_Snw___(:,:,1+nw)*tmp_svd_VUXM_nMwl____(:,:,1+nw);
tmp_svd_SVUXM_SMwl____(:,:,1+nw) = tmp_CTF_UX_S_k_q_Snw___(:,:,1+nw)*tmp_svd_VUXM_nMwl____(:,:,1+nw);
end;%for nw=0:tmp_n_w_max-1;
%tmp_svd_SVUXM_SMwl____ = multiprod(tmp_CTF_UX_S_k_q_Snw___,tmp_svd_VUXM_nMwl____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% tmp_svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
tmp_t_mult = tmp_t_mult + tmp_t;
tmp_o = tmp_n_w_max*tmp_n_S*tmp_n_M*n_n;
tmp_o_mult = tmp_o_mult + tmp_o;
tmp_t = tic();
tmp_svd_SVUXM_wSM___ = ifft(permute(tmp_svd_SVUXM_SMwl____,[3,1,2]),[],1)*tmp_n_w_max;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_wSM___: %0.6f',tmp_t)); end;
tmp_t_ifft = tmp_t_ifft + tmp_t;
tmp_o = tmp_n_w_max*log(tmp_n_w_max)*tmp_n_S*tmp_n_M;
tmp_o_ifft = tmp_o_ifft + tmp_o;
end;end;%for nS_batch=0:n_S_batch-1;for nM_batch=0:n_M_batch-1;
tmp_t_mult_(1+ntolerance_pm) = tmp_t_mult;
tmp_t_ifft_(1+ntolerance_pm) = tmp_t_ifft;
tmp_o_mult_(1+ntolerance_pm) = tmp_o_mult;
tmp_o_ifft_(1+ntolerance_pm) = tmp_o_ifft;
disp(sprintf(' %% ntolerance_pm %d n_n %d tmp_t_mult %0.5fs tmp_s_mult %0.5fs tmp_t_ifft %0.5fs tmp_s_ifft %0.5fs',ntolerance_pm,n_n,tmp_t_mult,tmp_o_mult/tmp_t_mult/1e9,tmp_t_ifft,tmp_o_ifft/tmp_t_ifft/1e9));
end;%for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
save(fname_mat ...
     ,'n_tolerance_pm','tolerance_pm_','n_w_max','n_S','n_M','tmp_n_w_max','tmp_n_S','tmp_n_M' ...
     ,'tmp_t_mult_','tmp_o_mult_','tmp_t_ifft_','tmp_o_ifft_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_fig_X_2d_cluster_d0_FIGI__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,512,768]);
p_row=3;p_col=1;ns=0;
tmp_ij=1:1:n_tolerance_pm; 
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
plot(pm_n_UX_rank_avg_p_(tmp_ij),-log10(tolerance_pm_(tmp_ij)),'ko-','MarkerFaceColor',0.85*[1,1,1]);
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylabel('$-\log_{10}(\mbox{tolerance})$','Interpreter','latex'); set(gca,'Ytick',0:2:20);
legend({'tolerance'},'Location','NorthWest');
grid on;
title('user-specified tolerance','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(pm_n_UX_rank_avg_p_(tmp_ij),corr_X_pe_(tmp_ij),'kd-','MarkerFaceColor',0.65*[1,1,1]);
plot(pm_n_UX_rank_avg_p_(tmp_ij),cap_nrm_pn__(end,tmp_ij),'ko-','MarkerFaceColor',[0.85,0.25,0.15]);
hold off;
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylabel('value','Interpreter','latex');ylim([0.5,1.0]);set(gca,'Ytick',0.5:0.1:1.0);
legend({'correlation','top 5%'},'Location','SouthEast');
grid on;
title('accuracy','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(pm_n_UX_rank_avg_p_(tmp_ij),max(tmp_t_mult_)./tmp_t_mult_(tmp_ij),'ks-','MarkerFaceColor',[0.15,0.35,0.85]);
plot(pm_n_UX_rank_avg_p_(tmp_ij),max(tmp_o_mult_)./tmp_o_mult_(tmp_ij),'k^-','MarkerFaceColor',[0.15,0.95,0.15]);
hold off;
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylim([0.0,8.0]);set(gca,'Ytick',0:1:8); ylabel('factor');
legend({'timing','operations'},'Location','NorthEast');
grid on;
title('efficiency','Interpreter','latex');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

disp('returning'); return; 
