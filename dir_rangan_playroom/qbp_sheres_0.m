function ...
[ ...
 parameter ...
,a_k_Y_zero_yk_ ...
,a_k_Y_zero_yk__ ...
,a_k_Y_0qbp_yk_ ...
,a_k_Y_0qbp_yk__ ...
,a_k_Y_sher_yk_ ...
,a_k_Y_sher_yk__ ...
,X_base_wSM___ ...
,delta_x_base_wSM___ ...
,delta_y_base_wSM___ ...
,gamma_z_base_wSM___ ...
,I_value_base_wSM___ ...
,frac_qk__ ...
,quad_from_data_T_CTF_normalized_qk__ ...
,expR2_sum_M_ ...
,numerator_qkM___ ...
,denomator_qkM___ ...
,Ylm_w_yq__ ...
] = ...
qbp_sheres_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,n_M ...
,M_k_p_wkM__ ...
,image_delta_x_prev_M_ ...
,image_delta_y_prev_M_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,index_ncluster_from_nCTF_ ...
,a_k_Y_base_yk_ ...
,X_base_wSM___ ...
,delta_x_base_wSM___ ...
,delta_y_base_wSM___ ...
,gamma_z_base_wSM___ ...
,I_value_base_wSM___ ...
,euler_polar_a_base_M_ ...
,euler_azimu_b_base_M_ ...
,euler_gamma_z_base_M_ ...
,image_delta_x_base_M_ ...
,image_delta_y_base_M_ ...
,image_I_value_base_M_ ...
,FTK ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (nargin<1);
%%%%%%%%;
test_pm_trpv1c_9b;
parameter = struct('type','parameter');
parameter.flag_verbose = 1;
parameter.sigma_sheres = 0.001;
CTF_k_p_r_kC__ = CTF_k_p_r__;
a_k_Y_base_yk_ = a_k_Y_quad_; %<-- start with oracle for now. ;
tolerance_pm = tolerance_master;
FTK = [];
delta_r_max = 0.1; %<-- default. ;
M_k_p_wkM__ = M_k_p__;
n_M = 256;
tmp_index_nM_ = 0:n_M-1;
M_k_p_wkM__ = M_k_p_wkM__(:,1+tmp_index_nM_);
CTF_k_p_r_kC__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+tmp_index_nM_));
index_nCTF_from_nM_ = tmp_index_nM_;
n_CTF = n_M;
%%%%%%%%;
[ ...
 parameter ...
,a_k_Y_zero_yk_ ...
,a_k_Y_zero_yk__ ...
,a_k_Y_0qbp_yk_ ...
,a_k_Y_0qbp_yk__ ...
,a_k_Y_sher_yk_ ...
,a_k_Y_sher_yk__ ...
,X_base_wSM___ ...
,delta_x_base_wSM___ ...
,delta_y_base_wSM___ ...
,gamma_z_base_wSM___ ...
,I_value_base_wSM___ ...
,frac_qk__ ...
,quad_from_data_T_CTF_normalized_qk__ ...
,expR2_sum_M_ ...
,numerator_qkM___ ...
,denomator_qkM___ ...
,Ylm_w_yq__ ...
] = ...
qbp_sheres_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,n_M ...
,M_k_p_wkM__ ...
,[] ...
,[] ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,[] ...
,a_k_Y_base_yk_ ...
);
%%%%%%%%;
disp(sprintf('returning')); return;
end;%if (nargin<1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_thisfunction = 'qbp_sheres_0';

a_k_Y_zero_yk_ = [];
a_k_Y_zero_yk__ = [];
a_k_Y_0qbp_yk_ = [];
a_k_Y_0qbp_yk__ = [];
a_k_Y_sher_yk_ = [];
a_k_Y_sher_yk__ = [];
X_base_wSM___ = [];
delta_x_base_wSM___ = [];
delta_y_base_wSM___ = [];
gamma_z_base_wSM___ = [];
I_value_base_wSM___ = [];
frac_qk__ = [];
quad_from_data_T_CTF_normalized_qk__ = [];
expR2_sum_M_ = [];
numerator_qkM___ = [];
denomator_qkM___ = [];
Ylm_w_yq__ = [];

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_prev_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_prev_M_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); index_ncluster_from_nCTF_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_base_yk_=[]; end; na=na+1;
if (nargin<1+na); X_base_wSM___=[]; end; na=na+1;
if (nargin<1+na); delta_x_base_wSM___=[]; end; na=na+1;
if (nargin<1+na); delta_y_base_wSM___=[]; end; na=na+1;
if (nargin<1+na); gamma_z_base_wSM___=[]; end; na=na+1;
if (nargin<1+na); I_value_base_wSM___=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_base_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_base_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_base_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_base_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_base_M_=[]; end; na=na+1;
if (nargin<1+na); image_I_value_base_M_=[]; end; na=na+1;
if (nargin<1+na); FTK=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_rank_vs_tolerance')); parameter.flag_rank_vs_tolerance = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'qbp_sheres_stop')); parameter.qbp_sheres_stop = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_pm')); parameter.tolerance_pm = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rank_pm')); parameter.rank_pm = 10; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'k_p_r_max')); parameter.k_p_r_max = k_p_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'template_viewing_k_eq_d')); parameter.template_viewing_k_eq_d = 1.0/max(1e-12,parameter.k_p_r_max); end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_xcor_vs_Memp')); parameter.flag_xcor_vs_Memp = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_M_per_Mbatch')); parameter.n_M_per_Mbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_optimize_over_gamma_z')); parameter.flag_optimize_over_gamma_z = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_compute_I_value')); parameter.flag_compute_I_value = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'sigma_sheres')); parameter.sigma_sheres = 0.1; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_rank_vs_tolerance = parameter.flag_rank_vs_tolerance;
qbp_sheres_stop = parameter.qbp_sheres_stop; %<-- 0 = no stop ; 1 = stop after innerproduct calculation. ; 2 = stop after standard qbp. ;
tolerance_pm = parameter.tolerance_pm;
rank_pm = parameter.rank_pm;
k_p_r_max = parameter.k_p_r_max;
template_viewing_k_eq_d = parameter.template_viewing_k_eq_d;
delta_r_max = parameter.delta_r_max;
flag_xcor_vs_Memp = parameter.flag_xcor_vs_Memp;
n_M_per_Mbatch = parameter.n_M_per_Mbatch;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
flag_optimize_over_gamma_z = parameter.flag_optimize_over_gamma_z;
flag_compute_I_value = parameter.flag_compute_I_value;
sigma_sheres = parameter.sigma_sheres;

if (flag_verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
flag_disp = flag_verbose>1; nf=0;

l_max_max = n_w_max/2 - 1;
if (l_max_max~=max(l_max_));
disp(sprintf(' %% Warning, l_max_max %d max(l_max_) %d in %s',l_max_max,max(l_max_),str_thisfunction));
end;%if (l_max_max~=max(l_max_));
assert(l_max_max==max(l_max_));
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

if isempty(FTK);
svd_eps = tolerance_master;
n_delta_v_requested = 0;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

%%%%%%%%;
% center images. ;
%%%%%%%%;
if (~isempty(image_delta_x_prev_M_) & ~isempty(image_delta_y_prev_M_));
tmp_t = tic();
for nM=0:n_M-1;
M_k_p_ = M_k_p_wkM__(:,1+nM);
image_delta_x_prev = image_delta_x_prev_M_(1+nM);
image_delta_y_prev = image_delta_y_prev_M_(1+nM);
if ((image_delta_x_prev~=0) | (image_delta_y_prev~=0));
M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
,+image_delta_x_prev ...
,+image_delta_y_prev ...
);
M_k_p_wkM__(:,1+nM) = M_k_p_;
end;%if ((image_delta_x_prev~=0) | (image_delta_y_prev~=0));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% transf_p_to_p: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'transf_p_to_p',tmp_t);
end;%if (~isempty(image_delta_x_prev_M_) & ~isempty(image_delta_y_prev_M_));
%%%%%%%%;
% construct M_k_q_wkM__. ;
%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_ = M_k_p_wkM__(:,1+nM);
M_k_q_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
M_k_q_wkM__(:,1+nM) = M_k_q_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% M_k_q_wkM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'M_k_q_wkM__',tmp_t);

%%%%%%%%;
% Find CTF-clusters. ;
%%%%%%%%;
if isempty(index_ncluster_from_nCTF_);
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
end;%if isempty(index_ncluster_from_nCTF_);
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
if (flag_verbose); disp(sprintf(' %% found %d CTF-clusters',n_cluster)); end;

%%%%%%%%;
% Now determine the principal-modes for each cluster. ;
%%%%%%%%;
if flag_xcor_vs_Memp==0;
%%%%;
% Empirical principal modes (using translated images). ;
%%%%;
X_2d_Memp_d1_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_Memp_d1_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
[ ...
 X_2d_Memp_d1_kk__ ...
,X_2d_Memp_d1_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_index_nM_from_ncluster ...
,M_k_p_wkM__(:,1+index_nM_from_ncluster_) ...
);
X_2d_Memp_d1_kkc___(:,:,1+ncluster) = X_2d_Memp_d1_kk__;
X_2d_Memp_d1_weight_rc__(:,1+ncluster) = X_2d_Memp_d1_weight_r_;
clear X_2d_Memp_d1_kk__ X_2d_Memp_d1_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
X_kkc___ = X_2d_Memp_d1_kkc___;
X_weight_rc__ = X_2d_Memp_d1_weight_rc__;
clear X_2d_Memp_d1_kkc__ X_2d_Memp_d1_weight_rc__;
end;%if flag_xcor_vs_Memp==0;
%%%%%%%%;
if flag_xcor_vs_Memp==1;
%%%%;
% Volumetric principal modes (assuming centered volume). ;
%%%%;
delta_sigma_base = 0.0; %<-- no translation as default. ;
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
end;%if flag_xcor_vs_Memp==1;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
X_kk__ = X_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
if flag_rank_vs_tolerance==0;
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
end;%if flag_rank_vs_tolerance==0;
if flag_rank_vs_tolerance==1;
pm_n_UX_rank = max(0,min(n_UX_rank,rank_pm));
end;%if flag_rank_vs_tolerance==1;
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (flag_verbose); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'UX_knc___',tmp_t);
%%%%%%%%;

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
% then calculate average CTFs for each cluster. ;
%%%%%%%%;
CTF_k_p_r_xavg_kc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_k_p_r_xavg_kc__(:,1+ncluster) = CTF_k_p_r_xavg_k_;
end;%for ncluster=0:n_cluster-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
c_80s__ = colormap_80s; n_c_80s = size(c_80s__,1);
markersize_use = 8;
hold on;
for ncluster=0:n_cluster-1;
nc_80s = max(0,min(n_c_80s-1,floor(n_c_80s*ncluster/n_cluster)));
plot(k_p_r_,CTF_k_p_r_xavg_kc__(:,1+ncluster),'ko-','MarkerFaceColor',c_80s__(1+nc_80s,:),'MarkerSize',markersize_use);
end;%for ncluster=0:n_cluster-1;
xlabel('k'); ylabel('CTF'); axis tight; 
end;%if flag_disp;

%%%%%%%%;
% Now form svd_VUXM_lwnM____ using the translated images. ;
%%%%%%%%;
tmp_t = tic();
pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_ncluster_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M,M_k_q_wkM__(:,1+index_nM_from_ncluster_),pm_n_UX_rank,UX_kn__,X_weight_r_);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images (using the displacement-updates). ;
% If we had not included the accumulated-displacements +image_delta_x_base_M_ and +image_delta_y_base_M_ above, ;
% we would add them to the displacement-updates below (also with a positive-sign). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;

%%%%%%%%;
% Construct templates using the volume. ;
%%%%%%%%;
a_k_Y_base_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_base_yk__(1:n_lm,1+nk_p_r) = a_k_Y_base_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
tmp_flag_verbose=0;
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
,template_viewing_weight_all_ ...
] = ...
pm_template_2( ...
 tmp_flag_verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_base_yk__ ...
,template_viewing_k_eq_d ...
,-1 ...
,n_w_max ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% pm_template_2 (n_S %d): %0.3fs',n_S,tmp_t)); end;
parameter = parameter_timing_update(parameter,'pm_template_2',tmp_t);
%%%%%%%%;
gamma_z_all_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_val_ = transpose(linspace(0,2*pi,n_w_(1+nk_p_r)+1));
gamma_z_all_(1+tmp_index_) = tmp_val_(1:n_w_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;
cc_ = cos(gamma_z_all_);
sc_ = sin(gamma_z_all_);
%%%%%%%%;
template_k_p_r_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
template_k_p_r_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
n_viewing_all = n_S;
template_k_c_0__ = zeros(n_w_sum,n_viewing_all);
template_k_c_1__ = zeros(n_w_sum,n_viewing_all);
template_k_c_2__ = zeros(n_w_sum,n_viewing_all);
for nviewing_all=0:n_viewing_all-1;
template_viewing_polar_a = template_viewing_polar_a_all_(1+nviewing_all); ca = cos(template_viewing_polar_a); sa = sin(template_viewing_polar_a);
template_viewing_azimu_b = template_viewing_azimu_b_all_(1+nviewing_all); cb = cos(template_viewing_azimu_b); sb = sin(template_viewing_azimu_b);
template_k_c_0__(:,1+nviewing_all) = (+cb*ca*cc_ - sb*sc_).*template_k_p_r_;
template_k_c_1__(:,1+nviewing_all) = (+sb*ca*cc_ + cb*sc_).*template_k_p_r_;
template_k_c_2__(:,1+nviewing_all) = (-sa*cc_            ).*template_k_p_r_;
end;%for nviewing_all=0:n_viewing_all-1;
%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;

%%%%%%%%;
% Now estimate innerproducts. ;
%%%%%%%%;
tmp_flag_isempty = ...
  0 ...
| isempty(X_base_wSM___) ...
| isempty(delta_x_base_wSM___) ...
| isempty(delta_y_base_wSM___) ...
| isempty(gamma_z_base_wSM___) ...
| isempty(I_value_base_wSM___) ...
;
if tmp_flag_isempty;
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
tmp_t = tic();
[ ...
 parameter ...
,X_base_wSM___ ...
,delta_x_base_wSM___ ...
,delta_y_base_wSM___ ...
,gamma_z_base_wSM___ ...
,I_value_base_wSM___ ...
] = ...
ampmh_X_cluster_wrap_wSM___11( ...
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
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_base_wSM___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__11',tmp_t);
%%%%%%%%;
end;%if tmp_flag_isempty;

if (qbp_sheres_stop>=1);
if (flag_verbose); disp(sprintf(' %% %s: qbp_sheres_stop %d, returning',str_thisfunction,qbp_sheres_stop)); end;
return;
end;%if (qbp_sheres_stop>=1);

%%%%%%%%;
if numel(X_base_wSM___)==n_S*n_M;
X_base_SM__ = X_base_wSM___;
delta_x_base_SM__ = delta_x_base_wSM___;
delta_y_base_SM__ = delta_y_base_wSM___;
gamma_z_base_SM__ = gamma_z_base_wSM___;
I_value_base_SM__ = I_value_base_wSM___;
end;%if numel(X_base_wSM___)==n_S*n_M;
if numel(X_base_wSM___)==n_w_max*n_S*n_M;
X_base_SM__ = zeros(n_S,n_M);
delta_x_base_SM__ = zeros(n_S,n_M);
delta_y_base_SM__ = zeros(n_S,n_M);
gamma_z_base_SM__ = zeros(n_S,n_M);
I_value_base_SM__ = zeros(n_S,n_M);
for nM=0:n_M-1;
for nS=0:n_S-1;
[X_base,tmp_ij_nw] = max(X_base_wSM___(:,1+nS,1+nM));
X_base_SM__(1+nS,1+nM) = X_base;
delta_x_base_SM__(1+nS,1+nM) = delta_x_base_wSM___(tmp_ij_nw,1+nS,1+nM);
delta_y_base_SM__(1+nS,1+nM) = delta_y_base_wSM___(tmp_ij_nw,1+nS,1+nM);
gamma_z_base_SM__(1+nS,1+nM) = gamma_z_base_wSM___(tmp_ij_nw,1+nS,1+nM);
I_value_base_SM__(1+nS,1+nM) = I_value_base_wSM___(tmp_ij_nw,1+nS,1+nM);
end;%for nS=0:n_S-1;
end;%for nM=0:n_M-1;
end;%if numel(X_base_wSM___)==n_w_max*n_S*n_M;
R2_base_SM__ = 2 - 2*X_base_SM__;
R2_min_M_ = reshape(min(R2_base_SM__,[],1),[n_M,1]);
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:n_M,R2_min_M_,'ko'); xlim([1,n_M]); ylim([0,4]); xlabel('image index'); ylabel('R2'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%;

%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
tmp_flag_isempty = ...
  0 ...
| isempty(euler_polar_a_base_M_) ...
| isempty(euler_azimu_b_base_M_) ...
| isempty(euler_gamma_z_base_M_) ...
| isempty(image_delta_x_base_M_) ...
| isempty(image_delta_y_base_M_) ...
| isempty(image_I_value_base_M_) ...
;
if tmp_flag_isempty;
parameter.flag_MS_vs_SM = 0; %<-- ensure maximum-likelihood. ;
parameter.f_rand = 0.00; %<-- ensure maximum-likelihood. ;
tmp_t = tic();
[ ...
 parameter ...
,euler_polar_a_base_M_ ...
,euler_azimu_b_base_M_ ...
,euler_gamma_z_base_M_ ...
,image_delta_x_base_M_ ...
,image_delta_y_base_M_ ...
,image_I_value_base_M_ ...
,image_X_value_base_M_ ...
,image_S_index_base_M_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter ...
,n_w_max ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,n_M ...
,X_base_SM__ ...
,delta_x_base_SM__ ...
,delta_y_base_SM__ ...
,gamma_z_base_SM__ ...
,I_value_base_SM__ ...
);
tmp_str = 'SM'; if (parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_MS_vs_SM_2',tmp_t);
end;%if tmp_flag_isempty;

%%%%%%%%;
% Set up quadrature on the sphere for qbp. ;
%%%%%%%%;
quad_k_eq_d_ = sqrt(4*pi./n_lm_);
flag_unique_n = 0;
if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1); flag_unique_n = 1; end;
if ~flag_unique_n; if (flag_verbose); disp(sprintf(' %% Note, flag_unique_n %d',flag_unique_n)); end; end;
nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
%%%%;
quad_k_eq_d = quad_k_eq_d_(1+nk_p_r);
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
%%%%;
[ ...
 quad_n_all ...
,quad_azimu_b_all_ ...
,quad_polar_a_all_ ...
,quad_weight_all_ ...
,quad_k_c_0_all_ ...
,quad_k_c_1_all_ ...
,quad_k_c_2_all_ ...
,~ ...
,~ ...
,~ ...
] = ...
sample_shell_5( ...
 1.0 ...
,quad_k_eq_d ...
,'L' ...
) ;
quad_k_c_qd__ = [ quad_k_c_0_all_ , quad_k_c_1_all_ , quad_k_c_2_all_ ];
%%%%;
Ylm__ = get_Ylm__(1+l_max,0:l_max,quad_n_all,quad_azimu_b_all_,quad_polar_a_all_);
Ylm_yq__ = zeros(n_lm,quad_n_all);
nml=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Ylm_yq__(1+nml,:) = Ylm__{1+l_val}(1+l_val+m_val,:);
nml=nml+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
Ylm_w_yq__ = Ylm_yq__ * sparse(1:quad_n_all,1:quad_n_all,quad_weight_all_,quad_n_all,quad_n_all);
%%%%;
[ ...
 data_k_p_polar_a_wS__ ...
,data_k_p_azimu_b_wS__ ...
,data_k_c_0_wS__ ...
,data_k_c_1_wS__ ...
,data_k_c_2_wS__ ...
] = ...
cg_rhs_1( ...
 n_S ...
,n_w_max ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,+zeros(n_S,1) ...
);
data_k_c_wSd__ = [ data_k_c_0_wS__(:) , data_k_c_1_wS__(:) , data_k_c_2_wS__(:) ];
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wSd__,'K',1)-1;
quad_from_data_qwS__ = sparse(1+index_quad_from_data_,1:n_w*n_S,1,quad_n_all,n_w*n_S);
n_quad_from_data_q_ = quad_from_data_qwS__*ones(n_w*n_S,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wSq__ = bsxfun(@rdivide,transpose(quad_from_data_qwS__),max(1,transpose(n_quad_from_data_q_)));
%%%%%%%%;

%%%%%%%%;
% Compute maximum-likelihood qbp. ;
% Here we use the same translations and principal-modes as for sheres. ;
% Thus, this should serve as a zero-temperature-limit. ;
%%%%%%%%;
quad_from_data_T_CTF_numerator_qk__ = zeros(quad_n_all,n_k_p_r);
quad_from_data_T_CTF_denomator_qk__ = zeros(quad_n_all,n_k_p_r);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster);
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
euler_polar_a_sub_M_ = euler_polar_a_base_M_(1+index_nM_from_ncluster_);
euler_azimu_b_sub_M_ = euler_azimu_b_base_M_(1+index_nM_from_ncluster_);
euler_gamma_z_sub_M_ = euler_gamma_z_base_M_(1+index_nM_from_ncluster_);
image_delta_x_sub_M_ = image_delta_x_base_M_(1+index_nM_from_ncluster_);
image_delta_y_sub_M_ = image_delta_y_base_M_(1+index_nM_from_ncluster_);
image_I_value_sub_M_ = image_I_value_base_M_(1+index_nM_from_ncluster_);
[ ...
 data_sub_k_p_polar_a_wM__ ...
,data_sub_k_p_azimu_b_wM__ ...
,data_sub_k_c_0_wM__ ...
,data_sub_k_c_1_wM__ ...
,data_sub_k_c_2_wM__ ...
] = ...
cg_rhs_1( ...
 tmp_n_M ...
,n_w ...
,euler_polar_a_sub_M_ ...
,euler_azimu_b_sub_M_ ...
,+euler_gamma_z_sub_M_ ...
);
data_sub_k_c_wMd__ = [ data_sub_k_c_0_wM__(:) , data_sub_k_c_1_wM__(:) , data_sub_k_c_2_wM__(:) ];
%%%%;
index_quad_from_data_sub_ = knnsearch(quad_k_c_qd__,data_sub_k_c_wMd__,'K',1)-1;
quad_from_data_sub_qwM__ = sparse(1+index_quad_from_data_sub_,1:n_w*tmp_n_M,1,quad_n_all,n_w*tmp_n_M);
n_quad_from_data_sub_q_ = quad_from_data_sub_qwM__*ones(n_w*tmp_n_M,1); %<-- number of data-points per quadrature-point. ;
data_sub_from_quad_wMq__ = bsxfun(@rdivide,transpose(quad_from_data_sub_qwM__),max(1,transpose(n_quad_from_data_sub_q_)));
%%%%%%%%;
T_k_p_sub_wkM__ = zeros(n_w_sum,tmp_n_M);
for nM_sub=0:tmp_n_M-1;
nM = index_nM_from_ncluster_(1+nM_sub);
image_delta_index = knnsearch([FTK.delta_x_,FTK.delta_y_],[+image_delta_x_sub_M_(1+nM_sub),+image_delta_y_sub_M_(1+nM_sub)],'K',1)-1;
tmp_svd_VUXM_lwn___ = svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM);
UX_T_k_q_ = reshape(FTK.svd_U_d_expiw_s__(1+image_delta_index,:)*reshape(tmp_svd_VUXM_lwn___,[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
UX_T_k_p_ = reshape(ifft(UX_T_k_q_,[],1)*sqrt(n_w_max),[n_w_max*pm_n_UX_rank,1]);
XU_UX_T_k_p_ = reshape(reshape(UX_T_k_p_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
T_k_p_sub_wkM__(:,1+nM_sub) = image_I_value_sub_M_(1+nM_sub) * XU_UX_T_k_p_;
end;%for nM_sub=0:tmp_n_M-1;
CTF_sub_wMn__ = reshape(permute(reshape(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),[n_w,n_k_p_r,tmp_n_M]),[1,3,2]),[n_w*tmp_n_M,n_k_p_r]);
CTF2_sub_qk__ = quad_from_data_sub_qwM__*abs(CTF_sub_wMn__).^2;
quad_from_data_sub_T_CTF_sub_numerator_qk__ = zeros(quad_n_all,n_k_p_r);
quad_from_data_sub_T_CTF_sub_denomator_qk__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
T_sub_wM__ = T_k_p_sub_wkM__(1+index_nw_,:);
CTF_sub_wM_ = CTF_sub_wMn__(:,1+nk_p_r);
CTF2_sub_q_ = CTF2_sub_qk__(:,1+nk_p_r);
T_CTF_sub_wM__ = T_sub_wM__.*reshape(CTF_sub_wM_,[n_w,tmp_n_M]);
quad_from_data_sub_T_CTF_sub_numerator_q_ = quad_from_data_sub_qwM__ * T_CTF_sub_wM__(:);
quad_from_data_sub_T_CTF_sub_numerator_qk__(:,1+nk_p_r) = quad_from_data_sub_T_CTF_sub_numerator_q_;
quad_from_data_sub_T_CTF_sub_denomator_q_ = CTF2_sub_q_;
quad_from_data_sub_T_CTF_sub_denomator_qk__(:,1+nk_p_r) = quad_from_data_sub_T_CTF_sub_denomator_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
quad_from_data_T_CTF_numerator_qk__ = quad_from_data_T_CTF_numerator_qk__ + quad_from_data_sub_T_CTF_sub_numerator_qk__;
quad_from_data_T_CTF_denomator_qk__ = quad_from_data_T_CTF_denomator_qk__ + quad_from_data_sub_T_CTF_sub_denomator_qk__;
end;%for ncluster=0:n_cluster-1;
quad_from_data_T_CTF_normalized_qk__ = bsxfun(@rdivide,quad_from_data_T_CTF_numerator_qk__,max(1e-12,quad_from_data_T_CTF_denomator_qk__));
a_k_Y_zero_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
quad_from_data_T_CTF_normalized_q_ = quad_from_data_T_CTF_normalized_qk__(:,1+nk_p_r);
a_k_Y_zero_yk__(:,1+nk_p_r) = conj(Ylm_w_yq__)*quad_from_data_T_CTF_normalized_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
a_k_Y_zero_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_zero_yk_(1+tmp_index_) = a_k_Y_zero_yk__(1:n_lm,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;

%%%%%%%%;
% also construct the standard qbp. ;
%%%%%%%%;
qbp_eps = tolerance_master;
[ ...
 a_k_Y_0qbp_yk_ ...
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
,euler_polar_a_base_M_ ...
,euler_azimu_b_base_M_ ...
,euler_gamma_z_base_M_ ...
,image_delta_x_base_M_ ...
,image_delta_y_base_M_ ...
,image_I_value_base_M_ ...
);
%%%%;
a_k_Y_0qbp_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_0qbp_yk__(1:n_lm,1+nk_p_r) = a_k_Y_0qbp_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;

if sigma_sheres<=0;
frac_qk__ = [];
numerator_qkM___ = [];
denomator_qkM___ = [];
expR2_sum_M_ = [];
a_k_Y_sher_yk__ = a_k_Y_zero_yk__;
a_k_Y_sher_yk_ = a_k_Y_zero_yk_;
end;%if sigma_sheres<=0;

if (qbp_sheres_stop>=2);
if (flag_verbose); disp(sprintf(' %% %s: qbp_sheres_stop %d, returning',str_thisfunction,qbp_sheres_stop)); end;
return;
end;%if (qbp_sheres_stop>=2);

if sigma_sheres> 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now put together sheres-style qbp. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
numerator_qkM___ = zeros(quad_n_all,n_k_p_r,n_M);
denomator_qkM___ = zeros(quad_n_all,n_k_p_r,n_M);
expR2_sum_M_ = zeros(n_M,1);
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
npair = 0; n_pair = n_M*n_S;
tmp_t_sher_pre = tic();
for ncluster=0:n_cluster-1;
%%%%%%%%%%%%%%%%;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: tmp_n_M %d',ncluster,n_cluster,tmp_n_M)); end;
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_k_p_r)); end;
pm_n_k_p_r = pm_n_UX_rank;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
CTF_UX_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S);
CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
CTF_UX_S_k_q_wnS___ = reshape(CTF_UX_S_k_q_wnS__,[n_w_max,pm_n_UX_rank,n_S]);
%%%%%%%%;
R2_sub_min_M_ = R2_min_M_(1+index_nM_from_ncluster_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now taking small temperature sigma_sheres %0.6f. ;',sigma_sheres)); end;
%%%%%%%%;
frac_qk__ = zeros(quad_n_all,n_k_p_r);
numerator_sub_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
denomator_sub_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
expR2_sub_sum_M_ = zeros(tmp_n_M,1);
%%%%%%%%;
tmp_n_Mbatch = ceil(tmp_n_M/n_M_per_Mbatch);
if (flag_verbose>1); disp(sprintf(' %% tmp_n_Mbatch %d',tmp_n_Mbatch)); end;
for nMbatch=0:tmp_n_Mbatch-1;
index_tmp_nM_in_Mbatch_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
index_tmp_nM_in_Mbatch_ = index_tmp_nM_in_Mbatch_(find(index_tmp_nM_in_Mbatch_<tmp_n_M)); tmp_n_M_sub = numel(index_tmp_nM_in_Mbatch_);
if (flag_verbose>1); disp(sprintf(' %% nMbatch %d/%d index_tmp_nM_in_Mbatch_ %d-->%d',nMbatch,tmp_n_Mbatch,index_tmp_nM_in_Mbatch_(1+0),index_tmp_nM_in_Mbatch_(1+tmp_n_M_sub-1))); end;
%if (flag_verbose>0 & mod(nMbatch,32)==0); disp(sprintf(' %% nMbatch %d/%d index_tmp_nM_in_Mbatch_ %d-->%d',nMbatch,tmp_n_Mbatch,index_tmp_nM_in_Mbatch_(1+0),index_tmp_nM_in_Mbatch_(1+tmp_n_M_sub-1))); end;
if (tmp_n_M_sub>0);
%%%%%%%%;
svd_VUXM_sub_lwnM____ = svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_(1+index_tmp_nM_in_Mbatch_));
UX_M_sub_l2_dM__ = UX_M_l2_dM__(:,1+index_nM_from_ncluster_(1+index_tmp_nM_in_Mbatch_));
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nSbatch=0:n_Sbatch-1;
index_nS_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_nS_in_Sbatch_ = index_nS_in_Sbatch_(find(index_nS_in_Sbatch_<n_S)); n_S_sub = numel(index_nS_in_Sbatch_);
if (flag_verbose>1); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
%if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
CTF_UX_S_k_q_sub_wnS___ = reshape(CTF_UX_S_k_q_wnS___(:,:,1+index_nS_in_Sbatch_),[pm_n_w_max,pm_n_UX_rank,n_S_sub]);
CTF_UX_S_sub_l2_S_ = CTF_UX_S_l2_S_(1+index_nS_in_Sbatch_);
%%%%%%%%;
tmp_X_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M_sub);
tmp_delta_x_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M_sub);
tmp_delta_y_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M_sub);
tmp_gamma_z_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M_sub);
tmp_I_value_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M_sub);
tmp_X_sub_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S_sub,tmp_n_M_sub);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,tmp_X_sub_wSM___ ...
,tmp_delta_x_sub_wSM___ ...
,tmp_delta_y_sub_wSM___ ...
,tmp_gamma_z_sub_wSM___ ...
,tmp_I_value_sub_wSM___ ...
,tmp_X_sub_dwSM____ ...
] = ...
ampmh_X_single_cluster_wSM___11( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S_sub ...
,CTF_UX_S_k_q_sub_wnS___ ...
,CTF_UX_S_sub_l2_S_ ...
,tmp_n_M_sub ...
,svd_VUXM_sub_lwnM____ ...
,UX_M_sub_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_single_cluster_wSM___11 sub',tmp_t);
%%%%%%%%;
tmp_index_ = reshape(repmat(transpose(0:n_w_max-1),[1,n_S_sub]) + n_w_max*repmat(reshape(index_nS_in_Sbatch_,[1,n_S_sub]),[n_w_max,1]),[n_w_max*n_S_sub,1]);
quad_from_data_sub_qwS__ = quad_from_data_qwS__(:,1+tmp_index_);
%%%%%%%%;
tmp_t = tic();
tmp_R2_wSdM____ = bsxfun(@minus,permute(2 - 2*tmp_X_sub_dwSM____,[2,3,1,4]),reshape(R2_sub_min_M_(1+index_tmp_nM_in_Mbatch_),[1,1,1,tmp_n_M_sub]));
tmp_expR2_p_wSdM____ = exp(-tmp_R2_wSdM____/(2*sigma_sheres^2));
tmp_p_p_wSdMk_____ = bsxfun(@times,reshape(tmp_expR2_p_wSdM____,[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M_sub,1]),reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]))/sigma_sheres^2;
tmp_p_q_wSdMk_____ = reshape(interp_p_to_q_block_0(n_S_sub*FTK.n_delta_v*tmp_n_M_sub*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*tmp_n_M_sub*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*tmp_n_M_sub*n_k_p_r,tmp_p_p_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r]);
expR2_sub_sum_M_(1+index_tmp_nM_in_Mbatch_) = expR2_sub_sum_M_(1+index_tmp_nM_in_Mbatch_) + reshape(sum(reshape(tmp_expR2_p_wSdM____,[n_w_max*n_S_sub*FTK.n_delta_v,tmp_n_M_sub]),1),[tmp_n_M_sub,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_sub_lwnM____,[FTK.n_svd_l,n_w_max*pm_n_UX_rank*tmp_n_M_sub]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,tmp_n_M_sub])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*tmp_n_M_sub,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r]),[n_w_max,1,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r]),reshape(conj(tmp_p_q_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*tmp_n_M_sub*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M_sub*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(bsxfun(@times,reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdMk_____,1),[n_w_max,1,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M_sub*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M_sub,n_k_p_r]);
numerator_sub_qkM___(:,:,1+index_tmp_nM_in_Mbatch_) = numerator_sub_qkM___(:,:,1+index_tmp_nM_in_Mbatch_) + permute(sum(tmp_numerator_qdMk____,2),[1,4,3,2]);
denomator_sub_qkM___(:,:,1+index_tmp_nM_in_Mbatch_) = denomator_sub_qkM___(:,:,1+index_tmp_nM_in_Mbatch_) + permute(sum(tmp_denomator_qdMk____,2),[1,4,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% numerator_sub_qkM___ denomator_sub_qkM___: %0.6fs',tmp_t)); end
%tmp_f = (n_M*n_S)/(tmp_n_M_sub*n_S_sub);
%if (flag_verbose>2); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
npair = npair + tmp_n_M_sub*n_S_sub;
tmp_t_sher_pos = toc(tmp_t_sher_pre);
tmp_f_sher = (n_pair)/max(1,npair);
if (flag_verbose>0); disp(sprintf(' %% calculation at %0.6fs so far <-- full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t_sher_pos,tmp_t_sher_pos*tmp_f_sher,tmp_t_sher_pos*tmp_f_sher/60,tmp_t_sher_pos*tmp_f_sher/3600)); end;
clear tmp_X_sub_wSM___;
clear tmp_delta_x_sub_wSM___;
clear tmp_delta_y_sub_wSM___;
clear tmp_gamma_z_sub_wSM___;
clear tmp_I_value_sub_wSM___;
clear tmp_X_sub_dwSM____;
clear tmp_R2_wSdM____;
clear tmp_expR2_p_wSdM____;
clear tmp_p_p_wSdMk_____;
clear tmp_p_q_wSdMk_____;
clear tmp_numerator_qdMk____;
clear tmp_denomator_qdMk____;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
%%%%%%%%;
end;%if (tmp_n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;
expR2_sum_M_(1+index_nM_from_ncluster_) = expR2_sub_sum_M_;
numerator_qkM___(:,:,1+index_nM_from_ncluster_) = numerator_sub_qkM___;
denomator_qkM___(:,:,1+index_nM_from_ncluster_) = denomator_sub_qkM___;
%%%%%%%%%%%%%%%%;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% now use numerator and denomator to construct data on the sphere. ;
%%%%%%%%;
tolerance_denom = 1e-12;
frac_qk__ = bsxfun(@rdivide,sum(bsxfun(@rdivide,numerator_qkM___,reshape(expR2_sum_M_,[1,1,n_M])),3),max(sum(bsxfun(@rdivide,denomator_qkM___,reshape(expR2_sum_M_,[1,1,n_M])),3),tolerance_denom));
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:n_M,expR2_sum_M_,'ko'); xlim([1,n_M]); ylim([0,4]); xlabel('image index'); ylabel('expR2_sum_M_','Interpreter','none'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
a_k_Y_sher_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
frac_q_ = frac_qk__(:,1+nk_p_r);
a_k_Y_sher_yk__(:,1+nk_p_r) = conj(Ylm_w_yq__)*frac_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
a_k_Y_sher_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_sher_yk_(1+tmp_index_) = a_k_Y_sher_yk__(1:n_lm,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;figsml;
plot(real(quad_from_data_T_CTF_normalized_qk__(:)),real(frac_qk__(:)),'.');
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Finished with sheres-style qbp. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if sigma_sheres> 0;

if (flag_verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

   



