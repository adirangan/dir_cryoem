function ...
[ ...
 parameter ...
,flag_image_delta_upd_M_ ...
,FTK ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,M_k_q_wkM__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,a_k_Y_reco_yki__ ...
,euler_polar_a_Mi__ ...
,euler_azimu_b_Mi__ ...
,euler_gamma_z_Mi__ ...
,image_delta_x_acc_Mi__ ...
,image_delta_y_acc_Mi__ ...
,image_delta_x_upd_Mi__ ...
,image_delta_y_upd_Mi__ ...
,flag_image_delta_upd_Mi__ ...
,image_I_value_Mi__ ...
,image_X_value_Mi__ ...
,image_S_index_Mi__ ...
] = ...
ampmut_5( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,FTK ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,n_M ...
,M_k_p_wkM__ ...
,M_k_q_wkM__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,image_delta_x_acc_M_ ...
,image_delta_y_acc_M_ ...
,image_delta_x_upd_M_ ...
,image_delta_y_upd_M_ ...
,flag_image_delta_upd_M_ ...
,image_I_value_M_ ...
);
%%%%%%%%;
% simple alternating minimization using principal modes. ;
% displacements drawn from disc (e.g., of radius delta_r_max). ;
% updates displacements (i.e., updates translations) each iteration. ;
% Assumes that images have been batched into groups corresponding to specific CTFs. ;
% ;
% Input: ;
% n_k_p_r: integer number of shells. ;
% k_p_r_: real array of size n_k_p_r. k_p_r for each shell. ;
% k_p_r_max: real maximum radius (of sphere) to be integrated. ;
% weight_3d_k_p_r_: real array of size n_k_p_r. radial quadrature-weights for volume integration. ;
% n_w_max: integer number of inplane_gamma_z values recorded on each image-ring. ;
%          Should be even. ;
%          Defined such that n_w_max = 2*(l_max_max+1);
% n_w_: integer array of size n_k_p_r. ;
%       fixed to be n_w_max*ones(n_k_p_r,1). ;
% FTK: structure produced by ampmh_FTK_1.m ;
% n_cluster: integer number of image-clusters (sorted by CTF). ;
% index_ncluster_from_nCTF_: integer array of size n_CTF. ;
%                            index_ncluster_from_nCTF_(1+nCTF) is the ncluster used for ctf CTF_k_p_r_kC__(:,1+nCTF). ;
% n_CTF: integer number of CTFs. ;
% index_nCTF_from_nM_: integer array of size n_M. nCTF for each images. ;
%                      index_nCTF_from_nM_(1+nM) is the nCTF used for image M_k_p_wkM__(:,1+nM). ;
% CTF_k_p_r_kC__: real array of size (n_k_p_r,n_CTF). CTF_k_p_r_ for each nCTF. ;
% pm_n_UX_rank_c_: integer array of size (n_cluster). pm rank used for each ncluster. ;
% UX_knc___: real cell array of size (n_k_p_r,pm_n_UX_rank_max,n_cluster). principal-modes for each image-cluster. ;
% X_weight_rc__: real array of size (n_k_p_r,n_cluster). radial weights associated with each image-cluster. ;
% n_M: integer number of images. ;
% M_k_p_wkM__: complex array of size (n_w_sum,n_M). stack of images in k_p_ format. ;
% M_k_q_wkM__: complex array of size (n_w_sum,n_M). stack of images in k_q_ format. ;
% svd_VUXM_lwnM____: complex array of size(FTK.n_svd_l,n_w_max,pm_n_UX_rank,n_M). Principal-images combined with FTK-modes. ;
%                    Note that, for ampmut_5, we assume nCTF-specific principal-modes from UX_Ckn___. ;
% UX_M_l2_dM__: real array of size(FTK.n_delta_v,n_M). ;
%               norms of principal-images (translated in accordance with the FTK.delta_x_ and FTK.delta_y_). ; 
%               Note that, for ampmut_5, we assume nCTF-specific principal-modes from UX_Ckn___. ;
% l_max_max: integer order used for the volumes a_k_Y_reco_yk_. ;
%            fixed to be n_w_max/2 - 1. ;
% euler_polar_a_M_: real array of size n_M. initial polar_a used for each image (randomized if empty). ;
% euler_azimu_b_M_: real array of size n_M. initial azimu_b used for each image (randomized if empty). ;
% euler_gamma_z_M_: real array of size n_M. initial gamma_z used for each image (randomized if empty). ;
% image_delta_x_acc_M_: real array of size n_M. accumulated delta_x used for each image (zero if empty). ;
% image_delta_y_acc_M_: real array of size n_M. accumulated delta_y used for each image (zero if empty). ;
% image_delta_x_upd_M_: real array of size n_M. update to delta_x used for each image (zero if empty). ;
% image_delta_y_upd_M_: real array of size n_M. update to delta_y used for each image (zero if empty). ;
% flag_image_delta_upd_M_: integer array of size n_M. flag identifying principal-images that need to be recalculated. (one if empty);
% image_I_value_M_: real array of size n_M. initial I_value used for each image (one if empty). ;
%%%%%%%%;
% Parameters: ;
% rseed: integer random seed. default 0. ;
% n_iteration: integer number of iterations for alternating minimization. default 16. ;
% delta_r_max: real number maximum displacement used for FTK. ;
%              default 0.1. ;
% delta_r_upb: real number maximum displacement used to bound accumulated translation. default 2*delta_r_max ;
% delta_r_upd_threshold: real number used to determine when to combine displacement-update with accumulated-displacement. ;
%                        default 0*delta_r_max (i.e., always combine). ;
%%%%;
% flag_I_value_vs_1: integer 0 or 1. ;
%                    if set to 1 will use I_value estimated from alignment. ;
%                    if set to 0 will ignore I_value. ;
%                    default 0. ;
% flag_qbp_vs_lsq: integer 0 or 1. ;
%                  if set to 1 will use quadrature-back-projection. ;
%                  if set to 0 will use least-squares ;
%                  default 1 (hardcoded). ;
% n_order: integer polynomial interpolation order used for conjugate-gradient interpolation operator. default 5. ;
% order_limit_MS: integer greater or equal to 0. This is the order to which the MS-phase model is limited during the first half of the iterations. ;
%                 default -1. (i.e., not used). ;
%%%%;
% svd_eps: real number FTK svd error. ;
%          default is tolerance_master. ;
% n_delta_v_requested: integer number of displacements requested (actual number of displacements will be stored in FTK.n_delta_v). ;
%                      default is 2*FTK.n_svd_l. ;
% flag_optimize_over_gamma: integer equal to 1. ; optimizes over gamma when aligning images (hardcoded). ;
% n_M_per_Mbatch: integer number of images to consider in each batch during alignment. default 24. ;
% n_S_per_Sbatch: integer number of templates to consider in each batch during alignment. default 24. ;
% flag_compute_I_value: integer set to 0. Will not compute I-values during alignment (hardcoded). ;
%%%%;
% flag_MS_vs_SM: integer 0 or 1. ;
%                if set to 1 will enforce uniform distribution of viewing angles. ;
%                if set to 0 will match templates to images. ;
%                default 1. ;
% f_rand: real number. fraction of best templates to select from when matching templates to images. default 0.05. ;
%%%%;
% flag_X_local_vs_global: integer 0 or 1. ;
%                         if set to 1 will apply local search to align images. ;
%                         if set to 0 will apply global search to align images. ;
%                         default 0 (hardcoded). ;
%%%%%%%%;
% ;
% Output: ;
% Many of the inputs are listed as outputs as well, so that subsequent iterations of this function can be called without needless recalculation. ;
% In addition, the following are returned: ;
% a_k_Y_reco_yki__: complex array of size (n_lm_sum,n_iteration). stack of output volumes in k_Y_ format. ;
% euler_polar_a_Mi__: real array of size (n_M,n_iteration). list of polar_a used for each image at each iteration. ;
% euler_azimu_b_Mi__: real array of size (n_M,n_iteration). list of azimu_b used for each image at each iteration. ;
% euler_gamma_z_Mi__: real array of size (n_M,n_iteration). list of gamma_z used for each image at each iteration. ;
% image_delta_x_acc_Mi__: real array of size (n_M,n_iteration). list of accumulated delta_x used for each image at each iteration. ;
% image_delta_y_acc_Mi__: real array of size (n_M,n_iteration). list of accumulated delta_y used for each image at each iteration. ;
% image_delta_x_upd_Mi__: real array of size (n_M,n_iteration). list of update to delta_x used for each image at each iteration. ;
% image_delta_y_upd_Mi__: real array of size (n_M,n_iteration). list of update to delta_y used for each image at each iteration. ;
% flag_image_delta_upd_Mi__: integer array of size (n_M,n_iteration). list of flag identifying principal-images that need to be recalculated. ;
% image_I_value_Mi__: real array of size (n_M,n_iteration). list of I_value (i.e., image intensity) used for each image at each iteration. ;
% image_X_value_Mi__: real array of size (n_M,n_iteration). list of X_S (i.e., template-image correlation) used for each image at each iteration. ;
% image_S_index_Mi__: real array of size (n_M,n_iteration). list of S_index (i.e., template index) used for each image at each iteration. ;
%%%%%%%%;

if nargin<1;
disp(sprintf(' %% testing ampmut_5'));
disp(sprintf(' %% Warning, not yet implemented'));
disp('returning'); return;
end;%if nargin<1;

verbose=2;
if (verbose); disp(sprintf(' %% [entering ampmut_5]')); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); FTK=[]; end; na=na+1;
if (nargin<1+na); n_cluster=[]; end; na=na+1;
if (nargin<1+na); index_ncluster_from_nCTF_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); pm_n_UX_rank_c_=[]; end; na=na+1;
if (nargin<1+na); UX_knc___=[]; end; na=na+1;
if (nargin<1+na); X_weight_rc__=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); M_k_q_wkM__=[]; end; na=na+1;
if (nargin<1+na); svd_VUXM_lwnM____=[]; end; na=na+1;
if (nargin<1+na); UX_M_l2_dM__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_acc_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_acc_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_upd_M_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_upd_M_=[]; end; na=na+1;
if (nargin<1+na); flag_image_delta_upd_M_=[]; end; na=na+1;
if (nargin<1+na); image_I_value_M_=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_delta_v_requested')); parameter.n_delta_v_requested = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_MS_vs_SM')); parameter.flag_MS_vs_SM = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'order_limit_MS')); parameter.order_limit_MS = -1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_upb')); parameter.delta_r_upb = 2*delta_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'template_viewing_k_eq_d')); parameter.template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max); end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_alternate_MS_vs_SM')); parameter.flag_alternate_MS_vs_SM = 1; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
rseed = parameter.rseed;
n_iteration = parameter.n_iteration;
n_delta_v_requested = parameter.n_delta_v_requested;
flag_MS_vs_SM = parameter.flag_MS_vs_SM;
order_limit_MS = parameter.order_limit_MS;
delta_r_max = parameter.delta_r_max;
delta_r_upb = parameter.delta_r_upb;
svd_eps = tolerance_master;
template_viewing_k_eq_d = parameter.template_viewing_k_eq_d;
flag_alternate_MS_vs_SM = parameter.flag_alternate_MS_vs_SM;

%%%%%%%%;
% index bounds. ;
%%%%%%%%;
n_w_max = n_w_max + mod(n_w_max,2); %<-- round up to nearest even number. ;
l_max_max = n_w_max/2 - 1; assert(l_max_max==max(l_max_));
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
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

pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
if (verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (verbose);

if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

if (nargout>1);
euler_polar_a_Mi__ = zeros(n_M,n_iteration);
euler_azimu_b_Mi__ = zeros(n_M,n_iteration);
euler_gamma_z_Mi__ = zeros(n_M,n_iteration);
image_delta_x_acc_Mi__ = zeros(n_M,n_iteration);
image_delta_y_acc_Mi__ = zeros(n_M,n_iteration);
image_delta_x_upd_Mi__ = zeros(n_M,n_iteration);
image_delta_y_upd_Mi__ = zeros(n_M,n_iteration);
flag_image_delta_upd_Mi__ = zeros(n_M,n_iteration);
image_I_value_Mi__ = zeros(n_M,n_iteration);
image_X_value_Mi__ = zeros(n_M,n_iteration);
image_S_index_Mi__ = zeros(n_M,n_iteration);
a_k_Y_reco_yki__ = zeros(n_lm_sum,n_iteration);
end;%if (nargout>1);
%%%%%%%%;
% initialize current euler-angles randomly. ;
%%%%%%%%;
if isempty(euler_polar_a_M_); euler_polar_a_M_ = 1*pi*rand(n_M,1); end;
if isempty(euler_azimu_b_M_); euler_azimu_b_M_ = 2*pi*rand(n_M,1); end;
if isempty(euler_gamma_z_M_); euler_gamma_z_M_ = 2*pi*rand(n_M,1); end;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
if flag_alternate_MS_vs_SM==0;
parameter.flag_MS_vs_SM = (niteration< floor(n_iteration/2)); %<-- 1,1,1,...,0,0,0,... ;
end;%if flag_alternate_MS_vs_SM==0;
if flag_alternate_MS_vs_SM==1;
parameter.flag_MS_vs_SM = mod(1+niteration,2); %<-- 1,0,1,0,1,0,... ;
end;%if flag_alternate_MS_vs_SM==1;
flag_MS_vs_SM = parameter.flag_MS_vs_SM;
%%%%%%%%;
% Construct M_k_q_wkM__ while taking into account the translations. ;
%%%%%%%%;
tmp_M_index_ = efind(flag_image_delta_upd_M_); tmp_n_M = numel(tmp_M_index_);
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
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images (using the displacement-updates). ;
% If we had not included the accumulated-displacements +image_delta_x_acc_M_ and +image_delta_y_acc_M_ above, ;
% we would add them to the displacement-updates below (also with a positive-sign). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____,+image_delta_x_upd_M_,+image_delta_y_upd_M_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
flag_image_delta_upd_M_ = zeros(n_M,1);
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
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_k_Y_reco_yk_: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'qbp_6',tmp_t);
%%%%%%%%;
% Now normalize a_k_Y_reco_yk_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
a_k_Y_reco_yk_ = spharm__normalize_1(n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_reco_yk_);
%%%%%%%%;
% Now store image-parameters. ;
%%%%%%%%;
if (nargout>1);
a_k_Y_reco_yki__(:,1+niteration) = a_k_Y_reco_yk_;
euler_polar_a_Mi__(:,1+niteration) = euler_polar_a_M_;
euler_azimu_b_Mi__(:,1+niteration) = euler_azimu_b_M_;
euler_gamma_z_Mi__(:,1+niteration) = euler_gamma_z_M_;
image_delta_x_acc_Mi__(:,1+niteration) = image_delta_x_acc_M_;
image_delta_y_acc_Mi__(:,1+niteration) = image_delta_y_acc_M_;
image_delta_x_upd_Mi__(:,1+niteration) = image_delta_x_upd_M_;
image_delta_y_upd_Mi__(:,1+niteration) = image_delta_y_upd_M_;
image_I_value_Mi__(:,1+niteration) = image_I_value_M_;
end;%if (nargout>1);
%%%%%%%%;
% Construct templates using the volume. ;
%%%%%%%%;
a_k_Y_reco_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_reco_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_reco_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
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
,a_k_Y_reco_yk__ ...
,template_viewing_k_eq_d ...
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
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__10',tmp_t);
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,image_delta_x_bit_M_ ...
,image_delta_y_bit_M_ ...
,image_I_value_M_ ...
,image_X_value_M_ ...
,image_S_index_M_ ...
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
tmp_str = 'SM'; if (parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_MS_vs_SM_2',tmp_t);
%%%%%%%%;
% update translations. ;
%%%%%%%%;
if (~isfield(parameter,'delta_r_upd_threshold')); parameter.delta_r_upd_threshold = 0.0*delta_r_max; end; %<-- parameter_bookmark. ;
delta_r_upd_threshold = parameter.delta_r_upd_threshold;
image_delta_x_upd_M_ = image_delta_x_upd_M_ + image_delta_x_bit_M_;
image_delta_y_upd_M_ = image_delta_y_upd_M_ + image_delta_y_bit_M_;
image_delta_r_upd_prenorm_M_ = sqrt(image_delta_x_upd_M_.^2 + image_delta_y_upd_M_.^2);
image_delta_x_tot_M_ = image_delta_x_acc_M_ + image_delta_x_upd_M_;
image_delta_y_tot_M_ = image_delta_y_acc_M_ + image_delta_y_upd_M_;
image_delta_r_tot_M_ = sqrt(image_delta_x_tot_M_.^2 + image_delta_y_tot_M_.^2);
image_delta_x_nrm_M_ = image_delta_x_tot_M_;
image_delta_y_nrm_M_ = image_delta_y_tot_M_;
tmp_index_ = efind(image_delta_r_tot_M_> delta_r_upb);
if (numel(tmp_index_)> 0);
if (verbose>1); disp(sprintf(' %% normalizing %d/%d image displacements',numel(tmp_index_),n_M)); end;
image_delta_x_nrm_M_(1+tmp_index_) = image_delta_x_tot_M_(1+tmp_index_)*delta_r_upb./image_delta_r_tot_M_(1+tmp_index_);
image_delta_y_nrm_M_(1+tmp_index_) = image_delta_y_tot_M_(1+tmp_index_)*delta_r_upb./image_delta_r_tot_M_(1+tmp_index_);
image_delta_x_upd_M_(1+tmp_index_) = image_delta_x_nrm_M_(1+tmp_index_) - image_delta_x_acc_M_(1+tmp_index_);
image_delta_y_upd_M_(1+tmp_index_) = image_delta_y_nrm_M_(1+tmp_index_) - image_delta_y_acc_M_(1+tmp_index_);
end;%if (numel(tmp_index_)> 0);
flag_image_delta_upd_M_ = zeros(n_M,1);
image_delta_r_upd_posnorm_M_ = sqrt(image_delta_x_upd_M_.^2 + image_delta_y_upd_M_.^2);
tmp_index_ = efind( (image_delta_r_upd_prenorm_M_>=delta_r_upd_threshold) | (image_delta_r_upd_posnorm_M_>=delta_r_upd_threshold) );
if (numel(tmp_index_)> 0);
if (verbose>1); disp(sprintf(' %% accumulating %d/%d image displacements',numel(tmp_index_),n_M)); end;
flag_image_delta_upd_M_(1+tmp_index_) = 1;
image_delta_x_acc_M_(1+tmp_index_) = image_delta_x_acc_M_(1+tmp_index_) + image_delta_x_upd_M_(1+tmp_index_);
image_delta_y_acc_M_(1+tmp_index_) = image_delta_y_acc_M_(1+tmp_index_) + image_delta_y_upd_M_(1+tmp_index_);
image_delta_x_upd_M_(1+tmp_index_) = 0;
image_delta_y_upd_M_(1+tmp_index_) = 0;
end;%if (numel(tmp_index_)> 0);
%%%%%%%%;
if (nargout>1);
if niteration<n_iteration-1;
flag_image_delta_upd_Mi__(:,1+niteration) = flag_image_delta_upd_M_; %<-- these actually refer to end of iteration. ;
image_X_value_Mi__(:,1+niteration) = image_X_value_M_; %<-- these actually refer to end of iteration. ;
image_S_index_Mi__(:,1+niteration) = image_S_index_M_; %<-- these actually refer to end of iteration. ;
end;%if niteration<n_iteration-1;
end;%if (nargout>1);
%%%%%%%%;
% Now return to beginning of loop. ;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished ampmut_5]')); end;
