
load('ampmh_X_wSM___5_debug.mat');
if (verbose); disp(sprintf(' %% [entering ampmh_X_wSM___5_experiment]')); end;

tmp_index_M_ = 0:n_M-1; %<-- consider all images. ;

%%%%%%%%;
% Find CTF-averaged templates.; 
%%%%%%%%;
tmp_t = tic();
VSCTF_avg_ = mean(VSCTF_Mc__(1+tmp_index_M_,:),1);
VSCTF_std_ = std(VSCTF_Mc__(1+tmp_index_M_,:),1,1);
disp(sprintf(' %% VSCTF std/avg:')); disp(num2str(VSCTF_std_./max(1e-12,abs(VSCTF_avg_))));
CTF_UX_S_k_p_wnS__ = zeros(pm_n_w_sum,n_S);
for nCTF_rank=0:n_CTF_rank-1;
CTF_UX_S_k_p_wnS__ = CTF_UX_S_k_p_wnS__ + UCTF_UX_S_k_p_wSc___(:,:,1+nCTF_rank) * VSCTF_avg_(1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
%%%%%%%%;
CTF_UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
CTF_UX_S_l2_(1+nS) = ...
innerproduct_p_quad( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_2d_k_p_r_/(2*pi) ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
);
end;%for nS=0:n_S-1;
%%%%%%%%;
CTF_UX_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
CTF_UX_S_k_q_wnS__(:,1+nS) = ...
interp_p_to_q( ...
 pm_n_k_p_r ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wnS__(:,1+nS) ...
); 
end;%for nS=0:n_S-1; 
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
SS_k_q_ = svd(CTF_UX_S_k_q_wnS__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(CTF_UX_S_k_q_wnS__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(CTF_UX_S_k_q_wnS__,n_S_rank);
if (verbose>1); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% combining UCTF_UX_S_k_p_wSc___ to form templates: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% visualise. ;
%%%%%%%%;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(svd_VUXM_lwnM____)),[1,3,4,2]),[FTK.n_svd_l*pm_n_UX_rank*n_M,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(svd_VUXM_lwnM____).^2,[1,3,4,2]),[FTK.n_svd_l*pm_n_UX_rank*n_M,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(CTF_UX_S_k_q_wnS__)),[n_w_max,pm_n_UX_rank,n_S]),[2,3,1]),[pm_n_UX_rank*n_S,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(CTF_UX_S_k_q_wnS__).^2,[n_w_max,pm_n_UX_rank,n_S]),[2,3,1]),[pm_n_UX_rank*n_S,n_w_max]),1))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_mat = sprintf('ampmh_X_wSM___5_experiment_2.mat');
if (~exist(fname_mat,'file')); disp(sprintf(' %% %s not found, creating',fname_mat));

tmp_index_M_ = 0:2:128; %<-- take a subset of images. ;
tmp_n_M = numel(tmp_index_M_);

%%%%%%%%;
X_SM_full__ = zeros(n_S,tmp_n_M);
delta_x_SM_full__ = zeros(n_S,tmp_n_M);
delta_y_SM_full__ = zeros(n_S,tmp_n_M);
gamma_z_SM_full__ = zeros(n_S,tmp_n_M);
I_value_SM_full__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.svd_eps_use = 0;
parameter.n_svd_l_use = 0;
parameter.n_delta_v_use = 0;
parameter.pm_n_UX_rank_use = 0;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 24;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_SM_full__ ...
,delta_x_SM_full__ ...
,delta_y_SM_full__ ...
,gamma_z_SM_full__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_full__: %0.3fs',tmp_t)); end;
t_full = tmp_t;

%%%%%%%%;
% test svd_eps_use: ;
% verdict: svd_eps_use = 0.01 manages to retain about 90% of the top 1% of matches. ;
%%%%%%%%;
svd_eps_use_ = 10.^[-3:0.25:0]; n_svd_eps_use = numel(svd_eps_use_);
t_rsvd_ = zeros(n_svd_eps_use,1);
c_XX_rsvd_ = zeros(n_svd_eps_use,1);
c_99_rsvd_ = zeros(n_svd_eps_use,1);
c_95_rsvd_ = zeros(n_svd_eps_use,1);
c_90_rsvd_ = zeros(n_svd_eps_use,1);
c_80_rsvd_ = zeros(n_svd_eps_use,1);
for nsvd_eps_use=0:n_svd_eps_use-1;
svd_eps_use = svd_eps_use_(1+nsvd_eps_use);
%%%%%%%%;
X_SM_rsvd__ = zeros(n_S,tmp_n_M);
delta_x_SM_rsvd__ = zeros(n_S,tmp_n_M);
delta_y_SM_rsvd__ = zeros(n_S,tmp_n_M);
gamma_z_SM_rsvd__ = zeros(n_S,tmp_n_M);
I_value_SM_rsvd__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.svd_eps_use = svd_eps_use;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 24;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_SM_rsvd__ ...
,delta_x_SM_rsvd__ ...
,delta_y_SM_rsvd__ ...
,gamma_z_SM_rsvd__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_rsvd__: %0.3fs',tmp_t)); end;
t_rsvd_(1+nsvd_eps_use) = tmp_t;
c_XX_rsvd_(1+nsvd_eps_use) = corr(X_SM_rsvd__(:),X_SM_full__(:));
c_99_rsvd_(1+nsvd_eps_use) = corr(X_SM_rsvd__(:)>prctile(X_SM_rsvd__(:),99),X_SM_full__(:)>prctile(X_SM_full__(:),99));
c_95_rsvd_(1+nsvd_eps_use) = corr(X_SM_rsvd__(:)>prctile(X_SM_rsvd__(:),95),X_SM_full__(:)>prctile(X_SM_full__(:),95));
c_90_rsvd_(1+nsvd_eps_use) = corr(X_SM_rsvd__(:)>prctile(X_SM_rsvd__(:),90),X_SM_full__(:)>prctile(X_SM_full__(:),90));
c_80_rsvd_(1+nsvd_eps_use) = corr(X_SM_rsvd__(:)>prctile(X_SM_rsvd__(:),80),X_SM_full__(:)>prctile(X_SM_full__(:),80));
end;%for nsvd_eps_use=0:n_svd_eps_use-1;

%%%%%%%%;
% test n_delta_v_use: ;
% verdict: n_delta_v_use = 27 manages to retain about 90% of the top 1% of matches. ;
% but n_delta_v_use = 20 or so still manages to retain about 90% of the top 5% of matches. ;
%%%%%%%%;
n_delta_v_use_ = max(1,round(FTK.n_delta_v*linspace(1,0,13))); n_n_delta_v_use = numel(n_delta_v_use_);
t_rdvu_ = zeros(n_n_delta_v_use,1);
c_XX_rdvu_ = zeros(n_n_delta_v_use,1);
c_99_rdvu_ = zeros(n_n_delta_v_use,1);
c_95_rdvu_ = zeros(n_n_delta_v_use,1);
c_90_rdvu_ = zeros(n_n_delta_v_use,1);
c_80_rdvu_ = zeros(n_n_delta_v_use,1);
for nn_delta_v_use=0:n_n_delta_v_use-1;
n_delta_v_use = n_delta_v_use_(1+nn_delta_v_use);
%%%%%%%%;
X_SM_rdvu__ = zeros(n_S,tmp_n_M);
delta_x_SM_rdvu__ = zeros(n_S,tmp_n_M);
delta_y_SM_rdvu__ = zeros(n_S,tmp_n_M);
gamma_z_SM_rdvu__ = zeros(n_S,tmp_n_M);
I_value_SM_rdvu__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.n_delta_v_use = n_delta_v_use;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 24;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_SM_rdvu__ ...
,delta_x_SM_rdvu__ ...
,delta_y_SM_rdvu__ ...
,gamma_z_SM_rdvu__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_rdvu__: %0.3fs',tmp_t)); end;
t_rdvu_(1+nn_delta_v_use) = tmp_t;
c_XX_rdvu_(1+nn_delta_v_use) = corr(X_SM_rdvu__(:),X_SM_full__(:));
c_99_rdvu_(1+nn_delta_v_use) = corr(X_SM_rdvu__(:)>prctile(X_SM_rdvu__(:),99),X_SM_full__(:)>prctile(X_SM_full__(:),99));
c_95_rdvu_(1+nn_delta_v_use) = corr(X_SM_rdvu__(:)>prctile(X_SM_rdvu__(:),95),X_SM_full__(:)>prctile(X_SM_full__(:),95));
c_90_rdvu_(1+nn_delta_v_use) = corr(X_SM_rdvu__(:)>prctile(X_SM_rdvu__(:),90),X_SM_full__(:)>prctile(X_SM_full__(:),90));
c_80_rdvu_(1+nn_delta_v_use) = corr(X_SM_rdvu__(:)>prctile(X_SM_rdvu__(:),80),X_SM_full__(:)>prctile(X_SM_full__(:),80));
end;%for nn_delta_v_use=0:n_n_delta_v_use-1;

%%%%%%%%;
% test pm_n_UX_rank_use: ;
% verdict: decreasing pm_n_UX_rank_use even slightly causes a dramatic reduction in the quality of the output, ;
% while only reducing the computation time by a little bit! ;
% do not recommend!. ;
%%%%%%%%;
pm_n_UX_rank_use_ = pm_n_UX_rank:-1:1; n_pm_n_UX_rank_use = numel(pm_n_UX_rank_use_);
t_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
c_XX_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
c_99_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
c_95_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
c_90_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
c_80_rpmn_ = zeros(n_pm_n_UX_rank_use,1);
for npm_n_UX_rank_use=0:n_pm_n_UX_rank_use-1;
pm_n_UX_rank_use = pm_n_UX_rank_use_(1+npm_n_UX_rank_use);
%%%%%%%%;
X_SM_rpmn__ = zeros(n_S,tmp_n_M);
delta_x_SM_rpmn__ = zeros(n_S,tmp_n_M);
delta_y_SM_rpmn__ = zeros(n_S,tmp_n_M);
gamma_z_SM_rpmn__ = zeros(n_S,tmp_n_M);
I_value_SM_rpmn__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.pm_n_UX_rank_use = pm_n_UX_rank_use;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 24;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_SM_rpmn__ ...
,delta_x_SM_rpmn__ ...
,delta_y_SM_rpmn__ ...
,gamma_z_SM_rpmn__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_rpmn__: %0.3fs',tmp_t)); end;
t_rpmn_(1+npm_n_UX_rank_use) = tmp_t;
c_XX_rpmn_(1+npm_n_UX_rank_use) = corr(X_SM_rpmn__(:),X_SM_full__(:));
c_99_rpmn_(1+npm_n_UX_rank_use) = corr(X_SM_rpmn__(:)>prctile(X_SM_rpmn__(:),99),X_SM_full__(:)>prctile(X_SM_full__(:),99));
c_95_rpmn_(1+npm_n_UX_rank_use) = corr(X_SM_rpmn__(:)>prctile(X_SM_rpmn__(:),95),X_SM_full__(:)>prctile(X_SM_full__(:),95));
c_90_rpmn_(1+npm_n_UX_rank_use) = corr(X_SM_rpmn__(:)>prctile(X_SM_rpmn__(:),90),X_SM_full__(:)>prctile(X_SM_full__(:),90));
c_80_rpmn_(1+npm_n_UX_rank_use) = corr(X_SM_rpmn__(:)>prctile(X_SM_rpmn__(:),80),X_SM_full__(:)>prctile(X_SM_full__(:),80));
end;%for npm_n_UX_rank_use=0:n_pm_n_UX_rank_use-1;

%%%%%%%%;
% test n_w_max_use: ;
% verdict: decreasing n_w_max_use down to 50 or 60 retains around 90% of the top 1% of matches. ;
% for this example this is limited by the information in the templates. ;
%%%%%%%%;
n_w_max_use_ = max(1,round(n_w_max*linspace(1,0,13))); n_n_w_max_use = numel(n_w_max_use_);
t_rnwm_ = zeros(n_n_w_max_use,1);
c_XX_rnwm_ = zeros(n_n_w_max_use,1);
c_99_rnwm_ = zeros(n_n_w_max_use,1);
c_95_rnwm_ = zeros(n_n_w_max_use,1);
c_90_rnwm_ = zeros(n_n_w_max_use,1);
c_80_rnwm_ = zeros(n_n_w_max_use,1);
for nn_w_max_use=0:n_n_w_max_use-1;
n_w_max_use = n_w_max_use_(1+nn_w_max_use);
%%%%%%%%;
X_SM_rnwm__ = zeros(n_S,tmp_n_M);
delta_x_SM_rnwm__ = zeros(n_S,tmp_n_M);
delta_y_SM_rnwm__ = zeros(n_S,tmp_n_M);
gamma_z_SM_rnwm__ = zeros(n_S,tmp_n_M);
I_value_SM_rnwm__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.n_w_max_use = n_w_max_use;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 24;
n_S_per_Sbatch = 24;
tmp_t = tic();
[ ...
 X_SM_rnwm__ ...
,delta_x_SM_rnwm__ ...
,delta_y_SM_rnwm__ ...
,gamma_z_SM_rnwm__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_rnwm__: %0.3fs',tmp_t)); end;
t_rnwm_(1+nn_w_max_use) = tmp_t;
c_XX_rnwm_(1+nn_w_max_use) = corr(X_SM_rnwm__(:),X_SM_full__(:));
c_99_rnwm_(1+nn_w_max_use) = corr(X_SM_rnwm__(:)>prctile(X_SM_rnwm__(:),99),X_SM_full__(:)>prctile(X_SM_full__(:),99));
c_95_rnwm_(1+nn_w_max_use) = corr(X_SM_rnwm__(:)>prctile(X_SM_rnwm__(:),95),X_SM_full__(:)>prctile(X_SM_full__(:),95));
c_90_rnwm_(1+nn_w_max_use) = corr(X_SM_rnwm__(:)>prctile(X_SM_rnwm__(:),90),X_SM_full__(:)>prctile(X_SM_full__(:),90));
c_80_rnwm_(1+nn_w_max_use) = corr(X_SM_rnwm__(:)>prctile(X_SM_rnwm__(:),80),X_SM_full__(:)>prctile(X_SM_full__(:),80));
end;%for nn_w_max_use=0:n_n_w_max_use-1;

profile clear;profile on;
%%%%%%%%;
% Finally, test a multiple reduction. ;
% verdict: preserves about 90% of the top 1% of matches. ;
%%%%%%%%;
t_rall = 0;
c_XX_rall = 0;
c_99_rall = 0;
c_95_rall = 0;
c_90_rall = 0;
c_80_rall = 0;
%%%%%%%%;
X_SM_rall__ = zeros(n_S,tmp_n_M);
delta_x_SM_rall__ = zeros(n_S,tmp_n_M);
delta_y_SM_rall__ = zeros(n_S,tmp_n_M);
gamma_z_SM_rall__ = zeros(n_S,tmp_n_M);
I_value_SM_rall__ = zeros(n_S,tmp_n_M);
parameter = struct('type','parameter');
parameter.svd_eps_use = 0.01;
parameter.n_delta_v_use = FTK.n_delta_v/2;
parameter.pm_n_UX_rank_use = pm_n_UX_rank;
parameter.n_w_max_use = n_w_max/2;
parameter.flag_optimize_over_gamma_z = 1;
n_M_per_Mbatch = 16;
n_S_per_Sbatch = 4;
tmp_t = tic();
[ ...
 X_SM_rall__ ...
,delta_x_SM_rall__ ...
,delta_y_SM_rall__ ...
,gamma_z_SM_rall__ ...
,~ ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_per_Sbatch ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM_rall__: %0.3fs',tmp_t)); end;
t_rall = tmp_t;
c_XX_rall = corr(X_SM_rall__(:),X_SM_full__(:));
c_99_rall = corr(X_SM_rall__(:)>prctile(X_SM_rall__(:),99),X_SM_full__(:)>prctile(X_SM_full__(:),99));
c_95_rall = corr(X_SM_rall__(:)>prctile(X_SM_rall__(:),95),X_SM_full__(:)>prctile(X_SM_full__(:),95));
c_90_rall = corr(X_SM_rall__(:)>prctile(X_SM_rall__(:),90),X_SM_full__(:)>prctile(X_SM_full__(:),90));
c_80_rall = corr(X_SM_rall__(:)>prctile(X_SM_rall__(:),80),X_SM_full__(:)>prctile(X_SM_full__(:),80));
disp(sprintf(' %% t_rall/t_full %0.3f/%0.3f; c_99_rall %0.2f',t_rall,t_full,c_99_rall));
profile viewer;
profile off;

save( fname_mat ...
      ,'t_full' ...
      ,'svd_eps_use_','n_svd_eps_use' ...
      ,'t_rsvd_' ...
      ,'c_XX_rsvd_' ...
      ,'c_99_rsvd_' ...
      ,'c_95_rsvd_' ...
      ,'c_90_rsvd_' ...
      ,'c_80_rsvd_' ...
      ,'n_delta_v_use_','n_n_delta_v_use' ...
      ,'t_rdvu_' ...
      ,'c_XX_rdvu_' ...
      ,'c_99_rdvu_' ...
      ,'c_95_rdvu_' ...
      ,'c_90_rdvu_' ...
      ,'c_80_rdvu_' ...
      ,'pm_n_UX_rank_use_','n_pm_n_UX_rank_use' ...
      ,'t_rpmn_' ...
      ,'c_XX_rpmn_' ...
      ,'c_99_rpmn_' ...
      ,'c_95_rpmn_' ...
      ,'c_90_rpmn_' ...
      ,'c_80_rpmn_' ...
      ,'n_w_max_use_','n_n_w_max_use' ...
      ,'t_rnwm_' ...
      ,'c_XX_rnwm_' ...
      ,'c_99_rnwm_' ...
      ,'c_95_rnwm_' ...
      ,'c_90_rnwm_' ...
      ,'c_80_rnwm_' ...
      ,'t_rall' ...
      ,'c_XX_rall' ...
      ,'c_99_rall' ...
      ,'c_95_rall' ...
      ,'c_90_rall' ...
      ,'c_80_rall' ...
      );

end;%if (~exist(fname_mat,'file')); disp(sprintf(' %% %s not found, creating',fname_mat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

load(fname_mat);

figure(1);clf;figbig;
%%%%%%%%;
subplot(2,2,1);
hold on;
plot(t_rsvd_,c_XX_rsvd_,'ko-');
plot(t_rsvd_,c_99_rsvd_,'mo-');
plot(t_rsvd_,c_95_rsvd_,'ro-');
plot(t_rsvd_,c_90_rsvd_,'go-');
plot(t_rsvd_,c_80_rsvd_,'co-');
plot(t_full*[1,1],[0,1],'k-');
xlim([0,max(t_full,max(t_rsvd_))]);
ylim([0,1]); grid on;
xlabel('time');
ylabel('recovery');
title('rsvd');
hold off;
%%%%%%%%;
subplot(2,2,2);
hold on;
plot(t_rdvu_,c_XX_rdvu_,'ko-');
plot(t_rdvu_,c_99_rdvu_,'mo-');
plot(t_rdvu_,c_95_rdvu_,'ro-');
plot(t_rdvu_,c_90_rdvu_,'go-');
plot(t_rdvu_,c_80_rdvu_,'co-');
plot(t_full*[1,1],[0,1],'k-');
xlim([0,max(t_full,max(t_rdvu_))]);
ylim([0,1]); grid on;
xlabel('time');
ylabel('recovery');
title('rdvu');
hold off;
%%%%%%%%;
subplot(2,2,3);
hold on;
plot(t_rpmn_,c_XX_rpmn_,'ko-');
plot(t_rpmn_,c_99_rpmn_,'mo-');
plot(t_rpmn_,c_95_rpmn_,'ro-');
plot(t_rpmn_,c_90_rpmn_,'go-');
plot(t_rpmn_,c_80_rpmn_,'co-');
plot(t_full*[1,1],[0,1],'k-');
xlim([0,max(t_full,max(t_rpmn_))]);
ylim([0,1]); grid on;
xlabel('time');
ylabel('recovery');
title('rpmn');
hold off;
%%%%%%%%;
subplot(2,2,4);
hold on;
plot(t_rnwm_,c_XX_rnwm_,'ko-');
plot(t_rnwm_,c_99_rnwm_,'mo-');
plot(t_rnwm_,c_95_rnwm_,'ro-');
plot(t_rnwm_,c_90_rnwm_,'go-');
plot(t_rnwm_,c_80_rnwm_,'co-');
plot(t_full*[1,1],[0,1],'k-');
xlim([0,max(t_full,max(t_rnwm_))]);
ylim([0,1]); grid on;
xlabel('time');
ylabel('recovery');
title('rnwm');
hold off;
disp(sprintf(' %% t_rall/t_full %0.3f/%0.3f; c_XX_rall %0.2f , c_99_rall %0.2f , c_95_rall %0.2f , c_90_rall %0.2f',t_rall,t_full,c_XX_rall,c_99_rall,c_95_rall,c_90_rall));

if (verbose); disp(sprintf(' %% [finished ampmh_X_wSM___5_experiment]')); end;

return;
