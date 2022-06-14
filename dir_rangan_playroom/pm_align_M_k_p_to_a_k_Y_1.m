function ...
[ ...
 parameter ...
,a_k_Y_reco_yki__ ...
,corr_a_k_Y_ ...
,euler_polar_a__ ...
,euler_azimu_b__ ...
,euler_gamma_z__ ...
,image_delta_x_acc__ ...
,image_delta_y_acc__ ...
,image_delta_x_upd__ ...
,image_delta_y_upd__ ...
,flag_image_delta_upd__ ...
,image_I_value__ ...
,image_X_value__ ...
,image_S_index__ ...
] = ...
pm_align_M_k_p_to_a_k_Y_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
);

verbose=1;
if (verbose); disp(sprintf(' %% [entering pm_align_M_k_p_to_a_k_Y_1]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); n_cluster=[]; end; na=na+1;
if (nargin<1+na); index_ncluster_from_nCTF_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_true_yk_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_pm')); parameter.tolerance_pm = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_upb')); parameter.delta_r_upb = 2*delta_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_delta_v_requested')); parameter.n_delta_v_requested = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_MS_vs_SM')); parameter.flag_MS_vs_SM = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'template_viewing_k_eq_d')); parameter.template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max); end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
tolerance_pm = parameter.tolerance_pm;
rseed = parameter.rseed;
n_iteration = parameter.n_iteration;
delta_r_max = parameter.delta_r_max;
delta_r_upb = parameter.delta_r_upb;
n_delta_v_requested = parameter.n_delta_v_requested;
flag_MS_vs_SM = parameter.flag_MS_vs_SM;
template_viewing_k_eq_d = parameter.template_viewing_k_eq_d;
svd_eps = tolerance_master;
%%%%%%%%;

%%%%%%%%;
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
l_max_max = max(l_max_);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%%%%%%%%;
n_cluster = 1+max(index_ncluster_from_nCTF_);
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;

%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;

%%%%%%%%;
% First calculate average CTFs for each cluster. ;
%%%%%%%%;
CTF_k_p_r_xavg_kc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_k_p_r_xavg_kc__(:,1+ncluster) = CTF_k_p_r_xavg_k_;
end;%for ncluster=0:n_cluster-1;

%%%%%%%%;
% Now calculate idealized principal-modes for each nCTF. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
a_k_Y_true_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_true_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_true_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
CTF_k_p_r_xavg_kk__ = CTF_k_p_r_xavg_k_*transpose(CTF_k_p_r_xavg_k_);
tmp_delta_sigma = delta_r_max;
if (tmp_delta_sigma> 0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_true_yk_,CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
end;%if (tmp_delta_sigma> 0);
if (tmp_delta_sigma==0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max,a_k_Y_true_yk__,CTF_k_p_r_xavg_kk__);
end;%if (tmp_delta_sigma==0);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_2d_xavg_dx_kkc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
tmp_X__ = X_2d_xavg_dx_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(tmp_X__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'X__',tmp_t);
%%%%%%%%;

pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
if (verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (verbose);

FTK = [];
if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

flag_check=0;
if flag_check;
%%%%%%%%;
% Form all principal-volumes a_CTF_UX_Y_true_ync__. ;
%%%%%%%%;
tmp_t = tic();
a_CTF_UX_Y_true_ync___ = zeros(pm_n_lm_max_max,pm_n_UX_rank_max,n_cluster);
for ncluster=0:n_cluster-1;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_2d_xavg_dx_weight_r_ = X_2d_xavg_dx_weight_rc__(:,1+ncluster);
a_CTF_UX_Y_true_yn__ = zeros(pm_n_lm_max,pm_n_UX_rank);
for pm_nUX_rank=0:pm_n_UX_rank-1;
a_CTF_UX_Y_true_y_ = zeros(pm_n_lm_max,1);
for nk_p_r=0:n_k_p_r-1;
a_CTF_UX_Y_true_y_ = a_CTF_UX_Y_true_y_ + a_k_Y_true_yk__(:,1+nk_p_r) * UX_kn__(1+nk_p_r,1+pm_nUX_rank) * CTF_k_p_r_xavg_k_(1+nk_p_r) * X_2d_xavg_dx_weight_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
a_CTF_UX_Y_true_yn__(:,1+pm_nUX_rank) = a_CTF_UX_Y_true_y_;
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
%a_CTF_UX_Y_true_yn__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_CTF_UX_Y_true_yn__(:));
a_CTF_UX_Y_true_ync___(1:pm_n_lm_max,1:pm_n_UX_rank,1+ncluster) = reshape(a_CTF_UX_Y_true_yn__,[pm_n_lm_max,pm_n_UX_rank]);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_CTF_UX_Y_true_ync___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now generate all principal-templates (storing requires lots of memory). ;
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
[ ...
 CTF_UX_S_k_p_wnSc___ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 tmp_verbose ...
,l_max_max ...
,pm_n_UX_rank_max*n_cluster ...
,reshape(a_CTF_UX_Y_true_ync___,[pm_n_lm_max_max,pm_n_UX_rank_max*n_cluster]) ...
,template_viewing_k_eq_d ...
,-1 ...
,n_w_max ...
);
CTF_UX_S_k_p_wnSc___ = permute(reshape(CTF_UX_S_k_p_wnSc___,[n_w_max*pm_n_UX_rank_max,n_cluster,n_S]),[1,3,2]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_template_2 (n_S %d): %0.3fs',n_S,tmp_t)); end;
parameter = parameter_timing_update(parameter,'pm_template_2',tmp_t);
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
[ ...
 S_k_p_wkS__ ...
,~ ...
,~ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
get_template_1( ...
 tmp_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,+1 ...
,-1 ...
,n_w_max*ones(n_k_p_r,1) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1 (n_S %d): %0.3fs',n_S,tmp_t)); end;
%%%%%%%%;
S_k_p_wSk___ = permute(reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
tmp_t = tic();
tmp_CTF_UX_S_k_p_wnSc___ = zeros(n_w_max*pm_n_UX_rank_max,n_S,n_cluster);
for ncluster=0:n_cluster-1;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_2d_xavg_dx_weight_r_ = X_2d_xavg_dx_weight_rc__(:,1+ncluster);
tmp_CTF_UX_S_k_p_wnSc___(1:pm_n_w_sum,1:n_S,1+ncluster) = reshape(permute(reshape(reshape(S_k_p_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_2d_xavg_dx_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% CTF_UX_S_k_p_wnSc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
disp(sprintf(' %% CTF_UX_S_k_p_wnSc___ vs tmp_CTF_UX_S_k_p_wnSc___: %0.16f',fnorm(CTF_UX_S_k_p_wnSc___ - tmp_CTF_UX_S_k_p_wnSc___)/fnorm(CTF_UX_S_k_p_wnSc___)));
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 tmp_verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_true_yk__ ...
,1/k_p_r_max ...
,-1 ...
,n_w_max ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_template_2 (n_S %d): %0.3fs',n_S,tmp_t)); end;
%%%%%%%%;
S_k_p_wSk___ = permute(reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
tmp_t = tic();
tmp_CTF_UX_S_k_p_wnSc___ = zeros(n_w_max*pm_n_UX_rank_max,n_S,n_cluster);
for ncluster=0:n_cluster-1;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_2d_xavg_dx_weight_r_ = X_2d_xavg_dx_weight_rc__(:,1+ncluster);
tmp_CTF_UX_S_k_p_wnSc___(1:pm_n_w_sum,1:n_S,1+ncluster) = reshape(permute(reshape(reshape(S_k_p_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_2d_xavg_dx_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% CTF_UX_S_k_p_wnSc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
disp(sprintf(' %% CTF_UX_S_k_p_wnSc___ vs tmp_CTF_UX_S_k_p_wnSc___: %0.16f',fnorm(CTF_UX_S_k_p_wnSc___ - tmp_CTF_UX_S_k_p_wnSc___)/fnorm(CTF_UX_S_k_p_wnSc___)));
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
] = ...
pm_template_2( ...
 tmp_verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_true_yk__ ...
,1/k_p_r_max ...
,-1 ...
,n_w_max ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_template_2 (n_S %d): %0.3fs',n_S,tmp_t)); end;
%%%%%%%%;
S_k_q_wk__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
% storing all the principal-templates takes lots of memory. ;
%%%%%%%%;
CTF_UX_S_k_q_wnSc___ = [];
CTF_UX_S_l2_Sc__ = [];
flag_check=0;
if flag_check;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
tmp_t = tic();
CTF_UX_S_k_q_wnSc___ = zeros(n_w_max*pm_n_UX_rank_max,n_S,n_cluster);
for ncluster=0:n_cluster-1;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_2d_xavg_dx_weight_r_ = X_2d_xavg_dx_weight_rc__(:,1+ncluster);
CTF_UX_S_k_q_wnSc___(1:pm_n_w_sum,1:n_S,1+ncluster) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_2d_xavg_dx_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% CTF_UX_S_k_q_wnSc___: %0.3fs',tmp_t)); end;
end;%if flag_check;
%%%%%%%%;

M_k_q_wkM__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];

%%%%%%%%;
euler_polar_a__ = zeros(n_M,n_iteration);
euler_azimu_b__ = zeros(n_M,n_iteration);
euler_gamma_z__ = zeros(n_M,n_iteration);
image_delta_x_acc__ = zeros(n_M,n_iteration);
image_delta_y_acc__ = zeros(n_M,n_iteration);
image_delta_x_upd__ = zeros(n_M,n_iteration);
image_delta_y_upd__ = zeros(n_M,n_iteration);
flag_image_delta_upd__ = zeros(n_M,n_iteration);
image_I_value__ = zeros(n_M,n_iteration);
image_X_value__ = zeros(n_M,n_iteration);
image_S_index__ = zeros(n_M,n_iteration);
a_k_Y_reco_yki__ = zeros(n_lm_sum,n_iteration);
corr_a_k_Y_ = zeros(n_iteration,1);
%%%%%%%%;
rng(rseed);
euler_polar_a_ = 1*pi*rand(n_M,1);
euler_azimu_b_ = 2*pi*rand(n_M,1);
euler_gamma_z_ = 2*pi*rand(n_M,1);
image_delta_x_acc_ = zeros(n_M,1); %<-- accumulated displacement (i.e., current image center). ;
image_delta_y_acc_ = zeros(n_M,1); %<-- accumulated displacement (i.e., current image center). ;
image_delta_x_upd_ = zeros(n_M,1); %<-- update to displacement (i.e., current image shift). ;
image_delta_y_upd_ = zeros(n_M,1); %<-- update to displacement (i.e., current image shift). ;
image_delta_x_bit_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
image_delta_y_bit_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
flag_image_delta_upd_ = ones(n_M,1); %<-- flag identifying principal-images that need to be recalculated. ;
image_I_value_ = ones(n_M,1);
M_k_q_wkM__ = zeros(n_w_sum,n_M);
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M);
UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
% Construct M_k_q_wkM__ while taking into account the translations. ;
%%%%%%%%;
tmp_M_index_ = efind(flag_image_delta_upd_); tmp_n_M = numel(tmp_M_index_);
if (verbose>0); disp(sprintf(' %% updating M_k_q_wkM__ for tmp_n_M %d/%d images',tmp_n_M,n_M)); end;
tmp_t = tic();
M_k_q_wkM__(:,1+tmp_M_index_) = zeros(n_w_sum,tmp_n_M);
for tmp_nM=0:tmp_n_M-1;
nM = tmp_M_index_(1+tmp_nM);
M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p_wkM__(:,1+nM) ...
,+image_delta_x_acc_(1+nM) ...
,+image_delta_y_acc_(1+nM) ...
);
M_k_q_wkM__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for tmp_nM=0:tmp_n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q_wkM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'M_k_q_wkM__',tmp_t);
%%%%%%%%;
% Now form svd_VUXM_lwnM____ using these translated images. ;
%%%%%%%%;
tmp_t = tic();
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_2d_xavg_dx_weight_r_ = X_2d_xavg_dx_weight_rc__(:,1+ncluster);
tmp_M_index_sub_ = intersect(tmp_M_index_,index_nM_from_ncluster_);
tmp_n_M_sub = numel(tmp_M_index_sub_);
if (tmp_n_M_sub> 0);
svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+tmp_M_index_sub_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M_sub,M_k_q_wkM__(:,1+tmp_M_index_sub_),pm_n_UX_rank,UX_kn__,X_2d_xavg_dx_weight_r_);
end;%if (tmp_n_M_sub> 0);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images (using the displacement-updates). ;
% If we had not included the accumulated-displacements +image_delta_x_acc_ and +image_delta_y_acc_ above, ;
% we would add them to the displacement-updates below (also with a positive-sign). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____,+image_delta_x_upd_,+image_delta_y_upd_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
flag_image_delta_upd_ = zeros(n_M,1);
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
qbp_eps = tolerance_master;
a_k_Y_reco_yk_ ...
= ...
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
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_k_Y_reco_yk_: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'qbp_6',tmp_t);
%%%%%%%%;
% Now correlate the two models. ;
%%%%%%%%;
[~,corr_a_k_Y] = ...
register_spharm_to_spharm_3( ...
 0*verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_true_yk_ ...
,a_k_Y_reco_yk_ ...
);
if (verbose); disp(sprintf(' %% niteration %d/%d, correlation %0.4f',niteration,n_iteration,corr_a_k_Y)); end;
%%%%%%%%;
% Now store image-parameters. ;
%%%%%%%%;
if (nargout>1);
a_k_Y_reco_yki__(:,1+niteration) = a_k_Y_reco_yk_;
corr_a_k_Y_(1+niteration)  = corr_a_k_Y;
euler_polar_a__(:,1+niteration) = euler_polar_a_;
euler_azimu_b__(:,1+niteration) = euler_azimu_b_;
euler_gamma_z__(:,1+niteration) = euler_gamma_z_;
image_delta_x_acc__(:,1+niteration) = image_delta_x_acc_;
image_delta_y_acc__(:,1+niteration) = image_delta_y_acc_;
image_delta_x_upd__(:,1+niteration) = image_delta_x_upd_;
image_delta_y_upd__(:,1+niteration) = image_delta_y_upd_;
image_I_value__(:,1+niteration) = image_I_value_;
end;%if (nargout>1);
%%%%%%%%;
% Now if we are not at the final iteration, align the principal images to the principal volume. ;
%%%%%%%%;
if (niteration<n_iteration-1);
%%%%%%%%;
% Use given volume to align principal-images. ;
% Groups principal-images by cluster. ;
% Calculates principal-templates associated with each cluster. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_cluster_wrap_SM__10( ...
 parameter ...
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
,X_2d_xavg_dx_weight_rc__ ...
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
%%%%%%%%;
% Use current correlations to update current euler-angles. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_bit_ ...
,image_delta_y_bit_ ...
,image_I_value_ ...
,image_X_value_ ...
,image_S_index_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter ...
,n_w_max ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,n_M ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% MS_vs_SM: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_MS_vs_SM_2',tmp_t);
%%%%%%%%;
% update translations. ;
%%%%%%%%;
if (~isfield(parameter,'delta_r_upd_threshold')); parameter.delta_r_upd_threshold = 0.0*delta_r_max; end; %<-- parameter_bookmark. ;
delta_r_upd_threshold = parameter.delta_r_upd_threshold;
image_delta_x_upd_ = image_delta_x_upd_ + image_delta_x_bit_;
image_delta_y_upd_ = image_delta_y_upd_ + image_delta_y_bit_;
image_delta_r_upd_prenorm_ = sqrt(image_delta_x_upd_.^2 + image_delta_y_upd_.^2);
image_delta_x_tot_ = image_delta_x_acc_ + image_delta_x_upd_;
image_delta_y_tot_ = image_delta_y_acc_ + image_delta_y_upd_;
image_delta_r_tot_ = sqrt(image_delta_x_tot_.^2 + image_delta_y_tot_.^2);
image_delta_x_nrm_ = image_delta_x_tot_;
image_delta_y_nrm_ = image_delta_y_tot_;
tmp_index_ = efind(image_delta_r_tot_> delta_r_upb);
if (numel(tmp_index_)> 0);
if (verbose>1); disp(sprintf(' %% normalizing %d/%d image displacements',numel(tmp_index_),n_M)); end;
image_delta_x_nrm_(1+tmp_index_) = image_delta_x_tot_(1+tmp_index_)*delta_r_upb./image_delta_r_tot_(1+tmp_index_);
image_delta_y_nrm_(1+tmp_index_) = image_delta_y_tot_(1+tmp_index_)*delta_r_upb./image_delta_r_tot_(1+tmp_index_);
image_delta_x_upd_(1+tmp_index_) = image_delta_x_nrm_(1+tmp_index_) - image_delta_x_acc_(1+tmp_index_);
image_delta_y_upd_(1+tmp_index_) = image_delta_y_nrm_(1+tmp_index_) - image_delta_y_acc_(1+tmp_index_);
end;%if (numel(tmp_index_)> 0);
flag_image_delta_upd_ = zeros(n_M,1);
image_delta_r_upd_posnorm_ = sqrt(image_delta_x_upd_.^2 + image_delta_y_upd_.^2);
tmp_index_ = efind( (image_delta_r_upd_prenorm_>=delta_r_upd_threshold) | (image_delta_r_upd_posnorm_>=delta_r_upd_threshold) );
if (numel(tmp_index_)> 0);
if (verbose>1); disp(sprintf(' %% accumulating %d/%d image displacements',numel(tmp_index_),n_M)); end;
flag_image_delta_upd_(1+tmp_index_) = 1;
image_delta_x_acc_(1+tmp_index_) = image_delta_x_acc_(1+tmp_index_) + image_delta_x_upd_(1+tmp_index_);
image_delta_y_acc_(1+tmp_index_) = image_delta_y_acc_(1+tmp_index_) + image_delta_y_upd_(1+tmp_index_);
image_delta_x_upd_(1+tmp_index_) = 0;
image_delta_y_upd_(1+tmp_index_) = 0;
end;%if (numel(tmp_index_)> 0);
%%%%%%%%;
flag_image_delta_upd__(:,1+niteration) = flag_image_delta_upd_; %<-- these actually refer to end of iteration. ;
image_X_value__(:,1+niteration) = image_X_value_; %<-- these actually refer to end of iteration. ;
image_S_index__(:,1+niteration) = image_S_index_; %<-- these actually refer to end of iteration. ;
end;%if niteration<n_iteration-1;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished pm_align_M_k_p_to_a_k_Y_1]')); end;
