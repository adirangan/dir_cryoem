
a_k_Y_reco_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_reco_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_reco_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
S_k_p_wkS__ = S_k_p__;
S_k_q_wkS__ = S_k_q__;

fname_mat = sprintf('%s_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now calculate principal-modes across all image-clusters. ;
%%%%%%%%;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+index_nCTF_from_nM_(1:n_M)),n_k_p_r);
SCTF_c_ = diag(SCTF_c__);
n_CTF_rank = min(efind(SCTF_c_/max(SCTF_c_)<1e-2));
n_CTF_rank = 1;
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+index_nCTF_from_nM_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% VSCTF_Mc__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_ = M_k_p__(:,1+nM);
M_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[X_2d_Semp_d1__,X_2d_Semp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,T_k_p__);
[UX_2d_Semp_d1__,SX_2d_Semp_d1__,VX_2d_Semp_d1__] = svds(X_2d_Semp_d1__,n_UX_rank); SX_2d_Semp_d1_ = diag(SX_2d_Semp_d1__);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_2d_Semp_d1__: %0.3fs',tmp_t)); end;
%%%%%%%%;
UX__ = UX_2d_Semp_d1__;
X_weight_r_ = X_2d_Semp_d1_weight_r_;
%%%%%%%%;
pm_n_UX_rank = 16;
pm_n_k_p_r = pm_n_UX_rank;
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now form principal-images. ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,a_UCTF_UX_Y_ync__ ... 
] = ...
a_UCTF_UX_Y_wrap_ync__0( ...
 parameter ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___,[n_w_max*pm_n_k_p_r,n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_...
,euler_gamma_z_true_...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_ync__: %0.3fs',tmp_t)); end;
%%%%%%%%;
a_UCTF_UX_Y_ync__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_ync__);
%%%%%%%%;
a_UCTF_UX_Y_ync___ = reshape(a_UCTF_UX_Y_ync__,[n_lm_max,pm_n_UX_rank,n_CTF_rank]);
%%%%%%%%;
verbose=2;
X_SMn___ = zeros(n_S,n_M,pm_n_UX_rank);
X_t_ = zeros(pm_n_UX_rank,1);
for pm_nUX_rank=0:pm_n_UX_rank-1;
%%%%%%%%;
% Use current principal-model to align principal-images. ;
% Groups principal-images by micrograph (i.e., inefficient if there are only a few images per micrograph). ;
% Calculates principal-templates associated with each micrograph. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
parameter.flag_compress_S = 0;
tmp_a_UCTF_UX_Y_ync__ = reshape(a_UCTF_UX_Y_ync___(:,1:1+pm_nUX_rank,:),[n_lm_max*(1+pm_nUX_rank),n_CTF_rank]);
tmp_t = tic();
[ ...
 parameter ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_wrap_wrap_SM__8( ...
 parameter ...
,FTK ...
,n_w_max ...
,l_max_max ...
,1+pm_nUX_rank ...
,n_CTF_rank ...
,tmp_a_UCTF_UX_Y_ync__ ...
,n_M ...
,index_nCTF_from_nM_ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____(:,:,1:1+pm_nUX_rank,:) ...
,UX_M_l2_dM__ ...
,[] ...
,zeros(n_M,1) ...
,zeros(n_M,1) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
X_SMn___(:,:,1+pm_nUX_rank) = X_SM__;
X_t_(1+pm_nUX_rank) = tmp_t;
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','n_UX_rank' ...
     ,'FTK','n_CTF_rank','SCTF_c_','VSCTF_Mc__','X_2d_Semp_d1__','X_2d_Semp_d1_weight_r_' ...
     ,'UX_2d_Semp_d1__','SX_2d_Semp_d1_' ...
     ,'pm_n_UX_rank' ...
     ,'a_UCTF_UX_Y_ync__' ...
     ,'X_SMn___','X_t_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_fig_X_2d_Semp_d1_FIGH__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));

figure(1+nf);nf=nf+1;figsml;

prctile_ = [5,15,50,85,95]; n_prctile = numel(prctile_);
cap_pn__ = zeros(n_prctile,pm_n_UX_rank);
X_SM_end__ = X_SMn___(:,:,end);
for nprctile=0:n_prctile-1;
cut_X_SM_end__ = zeros(size(X_SM_end__));
prc = prctile_(1+nprctile);
cut_X_SM_end__ = X_SM_end__ > repmat(prctile(X_SM_end__,prc),[n_S,1]);
for pm_nUX_rank=0:pm_n_UX_rank-1;
X_SM_sub__ = X_SMn___(:,:,1+pm_nUX_rank);
cut_X_SM_sub__ = zeros(size(X_SM_sub__));
cut_X_SM_sub__ = X_SM_sub__ > repmat(prctile(X_SM_sub__,prc),[n_S,1]);
cap_pn__(1+nprctile,1+pm_nUX_rank) = sum(cut_X_SM_end__.*cut_X_SM_sub__,'all');
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
end;%for nprctile=0:n_prctile-1;
cap_nrm_pn__ = cap_pn__./repmat(cap_pn__(:,end),[1,pm_n_UX_rank]);

subplot(1,3,1);
imagesc(X_SMn___(:,:,1+15)); colorbar;
subplot(1,3,2);
imagesc(X_SMn___(:,:,1+14)); colorbar;
subplot(1,3,3);
imagesc(X_SMn___(:,:,1+15) - X_SMn___(:,:,1+14)); colorbar;


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
