function ...
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
,M_alig_k_p_wkMS___ ...
,S_k_p_wkS__ ...
] = ...
image_align_to_template_0( ...
 parameter ....
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,index_nCTF_from_nM_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,n_M ...
,M_k_p_wkM__ ...
,n_S ...
,S_k_p_wkS__ ...
,image_delta_x_acc_M_0in_ ...
,image_delta_y_acc_M_0in_ ...
,image_delta_x_upd_M_0in_ ...
,image_delta_y_upd_M_0in_ ...
);

str_thisfunction = 'image_align_to_template_0';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_acc_M_0in_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_acc_M_0in_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_upd_M_0in_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_upd_M_0in_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_rank_vs_tolerance')); parameter.flag_rank_vs_tolerance = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_clump_vs_cluster')); parameter.flag_clump_vs_cluster = parameter.flag_rank_vs_tolerance; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_cluster')); parameter.tolerance_cluster = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_pm')); parameter.tolerance_pm = parameter.tolerance_master; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rank_pm')); parameter.rank_pm = 10; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rank_CTF'));
SCTF_ = svd(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1:n_M)));
rank_CTF = max(find(SCTF_/max(SCTF_)>parameter.tolerance_master));
parameter.rank_CTF = rank_CTF; %<-- parameter_bookmark. ;
end;%if (~isfield(parameter,'rank_CTF')); 
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_delta_v_requested')); parameter.n_delta_v_requested = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_upb')); parameter.delta_r_upb = 2*parameter.delta_r_max; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
flag_verbose = parameter.flag_verbose;
flag_rank_vs_tolerance = parameter.flag_rank_vs_tolerance;
parameter.flag_clump_vs_cluster = parameter.flag_rank_vs_tolerance; flag_clump_vs_cluster = parameter.flag_clump_vs_cluster; %<-- force;
tolerance_cluster = parameter.tolerance_cluster;
tolerance_pm = parameter.tolerance_pm;
rank_pm = parameter.rank_pm;
rank_CTF = parameter.rank_CTF;
delta_r_max = parameter.delta_r_max;
n_delta_v_requested = parameter.n_delta_v_requested;
delta_r_upb = parameter.delta_r_upb;
svd_eps = tolerance_master;
%%%%%%%%;

%%%%%%%%;
% index bounds. ;
%%%%%%%%;
n_w_max = n_w_max + mod(n_w_max,2); %<-- round up to nearest even number. ;
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% Initialize some simple templates. ;
% nS==0: one blob. ;
% nS==1: two blob. ;
%%%%%%%%;
if isempty(n_S);
weight_2d_k_all_ = reshape(ones(n_w_max,1)*reshape(weight_2d_k_p_r_,[1,n_k_p_r])/n_w_max/(4*pi^2),[n_w_sum,1]);
n_x = 128; x_ = transpose(linspace(-1,+1,n_x)); diameter_x_c = 2.0;
[x_0__,x_1__] = ndgrid(x_,x_); x_r__ = sqrt(x_0__.^2 + x_1__.^2);
tmp_g = @(x_0,x_1,mu_0,mu_1,tmp_sigma) 1./(2*pi*tmp_sigma.^2) .* exp(-( (x_0 - mu_0).^2 + (x_1 - mu_1).^2 )./(2*tmp_sigma.^2));
nS=0;
S_x_c_xxS___ = zeros(n_x,n_x,n_S);
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,1.50*tmp_sigma,1.50*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,1.75*tmp_sigma,1.75*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,2.00*tmp_sigma,2.00*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,2.25*tmp_sigma,2.25*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,2.50*tmp_sigma,2.50*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
tmp_sigma = diameter_x_c*sqrt(2)/16;
S_x_c__ = tmp_g(x_0__,x_1__,0,0,tmp_sigma) + tmp_g(x_0__,x_1__,2.75*tmp_sigma,2.75*tmp_sigma,tmp_sigma);
S_x_c_xxS___(:,:,1+nS) = S_x_c__; nS=nS+1;
%%%%;
n_S = nS;
clear S_x_c__;
S_k_p_wkS__ = zeros(n_w_sum,n_S);
T_x_c_xxS___ = zeros(n_x,n_x,n_S);
for nS=0:n_S-1;
S_k_p_wkS__(:,1+nS) = interp_x_c_to_k_p_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,S_x_c_xxS___(:,:,1+nS),n_k_p_r,k_p_r_,n_w_);
T_x_c_xxS___(:,:,1+nS) = real(interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wkS__(:,1+nS).*weight_2d_k_all_));
end;%for nS=0:n_S-1;
end;%if isempty(n_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_clump_vs_cluster==1);
disp(sprintf(' %% Warning, flag_clump_vs_cluster==1 not yet implemented in %s',str_thisfunction));
end;%if (flag_clump_vs_cluster==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_clump_vs_cluster==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Find clusters. ;
%%%%%%%%;
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
% Now determine the principal-modes for each cluster. ;
%%%%%%%%;
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
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
X_kk__ = X_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (flag_verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'X__',tmp_t);
%%%%%%%%;

%%%%%%%%;
% Now draw from ampmut_5. ;
%%%%%%%%;
FTK = [];
pm_n_UX_rank = rank_pm;
M_k_q_wkM__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];
VSCTF_Mc__ = [];
if (~isempty(image_delta_x_acc_M_0in_)); image_delta_x_acc_M_ = image_delta_x_acc_M_0in_; else image_delta_x_acc_M_= zeros(n_M,1); end;
if (~isempty(image_delta_y_acc_M_0in_)); image_delta_y_acc_M_ = image_delta_y_acc_M_0in_; else image_delta_y_acc_M_= zeros(n_M,1); end;
image_delta_x_upd_M_ = [];
image_delta_y_upd_M_ = [];
flag_image_delta_upd_M_ = [];
image_I_value_M_ = [];
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
%%%%%%%%;
pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
if (flag_verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (flag_verbose);
%%%%%%%%;
if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);
%%%%%%%%;
if isempty(image_delta_x_acc_M_); image_delta_x_acc_M_ = zeros(n_M,1); end; %<-- accumulated displacement (i.e., current image center). ;
if isempty(image_delta_y_acc_M_); image_delta_y_acc_M_ = zeros(n_M,1); end; %<-- accumulated displacement (i.e., current image center). ;
if isempty(image_delta_x_upd_M_); image_delta_x_upd_M_ = zeros(n_M,1); end; %<-- update to displacement (i.e., current image shift). ;
if isempty(image_delta_y_upd_M_); image_delta_y_upd_M_ = zeros(n_M,1); end; %<-- update to displacement (i.e., current image shift). ;
image_delta_x_bit_M_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
image_delta_y_bit_M_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
if isempty(flag_image_delta_upd_M_); flag_image_delta_upd_M_ = ones(n_M,1); end; %<-- flag identifying principal-images that need to be recalculated. ;
if isempty(image_I_value_M_); image_I_value_M_ = ones(n_M,1); end;
if isempty(M_k_q_wkM__); M_k_q_wkM__ = zeros(n_w_sum,n_M); end;
if isempty(svd_VUXM_lwnM____); svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M); end;
if isempty(UX_M_l2_dM__); UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M); end;
%%%%%%%%;
% Construct M_k_q_wkM__ while taking into account the translations. ;
%%%%%%%%;
tmp_M_index_ = efind(flag_image_delta_upd_M_); tmp_n_M = numel(tmp_M_index_);
if (flag_verbose>0); disp(sprintf(' %% updating M_k_q_wkM__ for tmp_n_M %d/%d images',tmp_n_M,n_M)); end;
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
,+image_delta_x_acc_M_(1+nM) ...
,+image_delta_y_acc_M_(1+nM) ...
);
M_k_q_wkM__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for tmp_nM=0:tmp_n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% M_k_q_wkM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'M_k_q_wkM__',tmp_t);
%%%%%%%%;
% Now form svd_VUXM_lwnM____ using these translated images. ;
%%%%%%%%;
tmp_t = tic();
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
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images (using the displacement-updates). ;
% If we had not included the accumulated-displacements +image_delta_x_acc_M_ and +image_delta_y_acc_M_ above, ;
% we would add them to the displacement-updates below (also with a positive-sign). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____,+image_delta_x_upd_M_,+image_delta_y_upd_M_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
flag_image_delta_upd_M_ = zeros(n_M,1);
%%%%%%%%;
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
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__10',tmp_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (flag_clump_vs_cluster==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% construct aligned images. ;
%%%%%%%%;
tmp_n_S = n_S;
M_alig_k_p_wkMS___ = zeros(n_w_sum,n_M,tmp_n_S);
for nS=0:tmp_n_S-1;
for nM=0:n_M-1;
delta_x = delta_x_SM__(1+nS,1+nM);
delta_y = delta_y_SM__(1+nS,1+nM);
gamma_z = gamma_z_SM__(1+nS,1+nM);
M_k_p_ = M_k_p_wkM__(:,1+nM);
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+delta_x,+delta_y);
M_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,M_k_p_,-gamma_z);
M_alig_k_p_wkMS___(:,1+nM,1+nS) = M_k_p_;
end;%for nM=0:n_M-1;
end;%for nS=0:tmp_n_S-1;

if (flag_verbose);
%%%%%%%%;
% display output. ;
%%%%%%%%;
n_x = 128; x_ = transpose(linspace(-1,+1,n_x)); diameter_x_c = 2.0;
weight_2d_k_all_ = reshape(ones(n_w_max,1)*reshape(weight_2d_k_p_r_,[1,n_k_p_r])/n_w_max/(4*pi^2),[n_w_sum,1]);
figure(1);clf;figmed;
p_row = 2; n_S_sub = min(6,tmp_n_S); p_col = 3*n_S_sub;
for nS=0:n_S_sub-1;
M_alig_k_p_ = mean(M_alig_k_p_wkMS___(:,:,1+nS),2);
M_alig_x_c_ = real(interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_alig_k_p_.*weight_2d_k_all_));
S_k_p_ = S_k_p_wkS__(:,1+nS);
S_x_c_ = real(interp_k_p_to_x_c_xxnufft(n_x,diameter_x_c,n_x,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_));
%%%%;
subplot(p_row,p_col,1+0*p_col+3*nS+0);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_));
axis image; axisnotick; title(sprintf('Real(S_k_p_) %d',nS),'Interpreter','none');
subplot(p_row,p_col,1+0*p_col+3*nS+1);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_));
axis image; axisnotick; title(sprintf('Imag(S_k_p_) %d',nS),'Interpreter','none');
subplot(p_row,p_col,1+0*p_col+3*nS+2);
imagesc_c(n_x,x_,n_x,x_,S_x_c_);
axis image; axisnotick; title(sprintf('Real(S_x_c_) %d',nS),'Interpreter','none');
%%%%;
subplot(p_row,p_col,1+1*p_col+3*nS+0);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_alig_k_p_));
axis image; axisnotick; title(sprintf('Real(M_alig_k_p_) %d',nS),'Interpreter','none');
subplot(p_row,p_col,1+1*p_col+3*nS+1);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_alig_k_p_));
axis image; axisnotick; title(sprintf('Imag(M_alig_k_p_) %d',nS),'Interpreter','none');
subplot(p_row,p_col,1+1*p_col+3*nS+2);
imagesc_c(n_x,x_,n_x,x_,M_alig_x_c_);
axis image; axisnotick; title(sprintf('Real(M_alig_x_c_) %d',nS),'Interpreter','none');
%%%%;
end;%for nS=0:n_S_sub-1;
%%%%%%%%;
end;%if (flag_verbose);




