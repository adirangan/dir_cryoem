function ...
[ ...
 parameter ...
,a_CTF_avg_UX_Y_reco_hyki___ ...
,corr_a_CTF_avg_UX_Y_hi__ ...
,euler_polar_a_hMi___ ...
,euler_azimu_b_hMi___ ...
,euler_gamma_z_hMi___ ...
,image_delta_x_acc_hMi___ ...
,image_delta_y_acc_hMi___ ...
,image_delta_x_upd_hMi___ ...
,image_delta_y_upd_hMi___ ...
,flag_image_delta_upd_hMi___ ...
,image_I_value_hMi___ ...
,image_X_value_hMi___ ...
,image_S_index_hMi___ ...
] = ...
pm_align_M_k_p_to_a_het_k_Y_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
,n_CTF ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_cmb_k_Y_true_ ...
,n_volume ...
,a_het_k_Y_true_hyk__ ...
);

verbose=2;
if (verbose); disp(sprintf(' %% [entering pm_align_M_k_p_to_a_het_k_Y_0]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_upb')); parameter.delta_r_upb = 2*delta_r_max; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_delta_v_requested')); parameter.n_delta_v_requested = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_MS_vs_SM')); parameter.flag_MS_vs_SM = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'template_viewing_k_eq_d')); parameter.template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max); end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
rseed = parameter.rseed;
n_iteration = parameter.n_iteration;
delta_r_max = parameter.delta_r_max;
delta_r_upb = parameter.delta_r_upb;
n_delta_v_requested = parameter.n_delta_v_requested;
flag_MS_vs_SM = parameter.flag_MS_vs_SM;
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

%%%%%%%%;
CTF_k_p_r__ = zeros(n_k_p_r,n_CTF);
for nctf=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_k_p__(1+tmp_index_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_CTF-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xcor__ = CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1))) * transpose(CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1)))) / n_M;
%%%%%%%%;
SCTF_ = svd(CTF_k_p_r__(:,1+CTF_index_(1:n_M)));
n_CTF_rank = max(find(SCTF_/max(SCTF_)>tolerance_master));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_index_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;

%%%%%%%%;
% First calculate idealized principal-modes for combined molecule. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
a_cmb_k_Y_true__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_cmb_k_Y_true__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_cmb_k_Y_true_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = tic();
[ ...
 X__ ...
,X_weight_r_ ...
] = ...
principled_marching_cost_matrix_3( ...
 n_k_p_r ...
,weight_2d_k_p_r_ ...
,l_max_max ...
,a_cmb_k_Y_true__ ...
,CTF_k_p_r_xcor__ ...
);
[UX__,SX__,VX__] = svds(X__,n_UX_rank); SX_ = diag(SX__);
pm_n_UX_rank = max(find(SX_/max(SX_)> tolerance_master));
if (verbose); disp(sprintf(' %% pm_n_UX_rank %d/%d',pm_n_UX_rank,numel(SX_))); end;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'X__',tmp_t);

FTK = [];
pm_n_UX_rank = pm_n_UX_rank;
M_k_q__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];
VSCTF_Mc__ = [];

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

if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

%%%%%%%%;
% Form VSCTF_Mc__. ;
%%%%%%%%;
if isempty(VSCTF_Mc__);
if (n_CTF_rank<=0 | isempty(CTF_index_));
n_CTF_rank = 1;
USCTF_kc__ = ones(n_k_p_r,1); SCTF_C__ = 1; VCTF_Mc__ = ones(n_M,1);
else;
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_index_(1:n_M)),n_CTF_rank);
end;%if (n_CTF_rank<=0 | isempty(CTF_index_));
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
end;%if isempty(VSCTF_Mc__);

%%%%%%%%;
euler_polar_a_hMi___ = zeros(n_volume,n_M,n_iteration);
euler_azimu_b_hMi___ = zeros(n_volume,n_M,n_iteration);
euler_gamma_z_hMi___ = zeros(n_volume,n_M,n_iteration);
image_delta_x_acc_hMi___ = zeros(n_volume,n_M,n_iteration);
image_delta_y_acc_hMi___ = zeros(n_volume,n_M,n_iteration);
image_delta_x_upd_hMi___ = zeros(n_volume,n_M,n_iteration);
image_delta_y_upd_hMi___ = zeros(n_volume,n_M,n_iteration);
flag_image_delta_upd_hMi___ = zeros(n_volume,n_M,n_iteration);
image_I_value_hMi___ = zeros(n_volume,n_M,n_iteration);
image_X_value_hMi___ = zeros(n_volume,n_M,n_iteration);
image_S_index_hMi___ = zeros(n_volume,n_M,n_iteration);
a_CTF_avg_UX_Y_reco_hyki___ = zeros(n_volume,pm_n_lm_sum,n_iteration);
corr_a_CTF_avg_UX_Y_hi__ = zeros(n_volume,n_iteration,1);
%%%%%%%%;

for nvolume=0:n_volume-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (verbose); disp(sprintf(' %% nvolume %d/%d',nvolume,n_volume)); end;
a_k_Y_true_ = reshape(a_het_k_Y_true_hyk__{1+nvolume},[n_lm_sum,1]);
a_k_Y_true__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_true__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_true_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%;
% Form a_UCTF_UX_Y_true_ync__. ;
%%%%%%%%;
a_UCTF_UX_Y_true_ync__ = zeros(pm_n_lm_sum,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
a_UCTF_UX_Y_true_yn__ = zeros(pm_n_lm_max,pm_n_UX_rank);
for pm_nUX_rank=0:pm_n_UX_rank-1;
a_UCTF_UX_Y_true_y_ = zeros(pm_n_lm_max,1);
for nk_p_r=0:n_k_p_r-1;
a_UCTF_UX_Y_true_y_ = a_UCTF_UX_Y_true_y_ + a_k_Y_true__(:,1+nk_p_r) * UX__(1+nk_p_r,1+pm_nUX_rank) * UCTF_kc__(1+nk_p_r,1+nCTF_rank);
end;%for nk_p_r=0:n_k_p_r-1;
a_UCTF_UX_Y_true_yn__(:,1+pm_nUX_rank) = a_UCTF_UX_Y_true_y_;
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
a_UCTF_UX_Y_true_ync__(:,1+nCTF_rank) = a_UCTF_UX_Y_true_yn__(:);
end;%for nCTF_rank=0:n_CTF_rank-1;

%%%%%%%%;
% Now normalize a_UCTF_UX_Y_true_ync__. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
a_UCTF_UX_Y_true_ync__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_true_ync__);
%%%%%%%%;
% Use a_UCTF_UX_Y_true_ync__ as well VSCTF_Mc__ to approximate the image-averaged a_CTF_avg_UX_Y_. ;
% This is not actually used in the calculation, but can be useful for postprocessing. ;
%%%%%%%%;
a_CTF_avg_UX_Y_true_ = spharm_normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,mean(a_UCTF_UX_Y_true_ync__*transpose(VSCTF_Mc__),2));

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
a_CTF_avg_UX_Y_reco__ = zeros(pm_n_lm_sum,n_iteration);
corr_a_CTF_avg_UX_Y_ = zeros(n_iteration,1);
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
M_k_q__ = zeros(n_w_sum,n_M);
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank,n_M);
UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
% Construct M_k_q__ while taking into account the translations. ;
%%%%%%%%;
tmp_M_index_ = efind(flag_image_delta_upd_); tmp_n_M = numel(tmp_M_index_);
if (verbose>0); disp(sprintf(' %% updating M_k_q__ for tmp_n_M %d/%d images',tmp_n_M,n_M)); end;
tmp_t = tic();
M_k_q__(:,1+tmp_M_index_) = zeros(n_w_sum,tmp_n_M);
for tmp_nM=0:tmp_n_M-1;
nM = tmp_M_index_(1+tmp_nM);
M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,n_w_sum ...
,M_k_p__(:,1+nM) ...
,+image_delta_x_acc_(1+nM) ...
,+image_delta_y_acc_(1+nM) ...
);
M_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for tmp_nM=0:tmp_n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'M_k_q__',tmp_t);
%%%%%%%%;
% Now form svd_VUXM_lwnM____ using these translated images. ;
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M,M_k_q__(:,1+tmp_M_index_),pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images (using the displacement-updates). ;
% If we had not included the accumulated-displacements +image_delta_x_acc_ and +image_delta_y_acc_ above, ;
% we would add them to the displacement-updates below (also with a positive-sign). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M,svd_VUXM_lwnM____,+image_delta_x_upd_,+image_delta_y_upd_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
flag_image_delta_upd_ = zeros(n_M,1);
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,a_UCTF_UX_Y_reco_ync__ ... 
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
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_I_value_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_reco_ync__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
%%%%%%%%;
% Now normalize a_CTF_avg_UX_Y_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
a_UCTF_UX_Y_reco_ync__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_reco_ync__);
%%%%%%%%;
% Use a_UCTF_UX_Y_reco_ync__ as well VSCTF_Mc__ to approximate the image-averaged a_CTF_avg_UX_Y_. ;
% This is not actually used in the calculation, but can be useful for postprocessing. ;
%%%%%%%%;
a_CTF_avg_UX_Y_reco_ = spharm_normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,mean(a_UCTF_UX_Y_reco_ync__*transpose(VSCTF_Mc__),2));
%%%%%%%%;
% Now correlate the two models. ;
%%%%%%%%;
corr_a_CTF_avg_UX_Y = corr(a_CTF_avg_UX_Y_true_,a_CTF_avg_UX_Y_reco_);
if (verbose); disp(sprintf(' %% niteration %d/%d, correlation %0.4f',niteration,n_iteration,corr_a_CTF_avg_UX_Y)); end;
%%%%%%%%;
% Now store image-parameters. ;
%%%%%%%%;
if (nargout>1);
a_CTF_avg_UX_Y_reco__(:,1+niteration) = a_CTF_avg_UX_Y_reco_;
corr_a_CTF_avg_UX_Y_(1+niteration)  = corr_a_CTF_avg_UX_Y;
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
% Use given principal-model to align principal-images. ;
% Groups principal-images by micrograph (i.e., inefficient if there are only a few images per micrograph). ;
% Calculates principal-templates associated with each micrograph. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
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
,pm_n_UX_rank ...
,n_CTF_rank ...
,a_UCTF_UX_Y_true_ync__ ...
,n_M ...
,CTF_index_ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,[] ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_wrap_wrap_SM__8',tmp_t);
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
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

a_CTF_avg_UX_Y_reco_hyki___(1+nvolume,:,:) = a_CTF_avg_UX_Y_reco__;
corr_a_CTF_avg_UX_Y_hi__(1+nvolume,:) = corr_a_CTF_avg_UX_Y_;
euler_polar_a_hMi___(1+nvolume,:,:) = euler_polar_a__;
euler_azimu_b_hMi___(1+nvolume,:,:) = euler_azimu_b__;
euler_gamma_z_hMi___(1+nvolume,:,:) = euler_gamma_z__;
image_delta_x_acc_hMi___(1+nvolume,:,:) = image_delta_x_acc__;
image_delta_y_acc_hMi___(1+nvolume,:,:) = image_delta_y_acc__;
image_delta_x_upd_hMi___(1+nvolume,:,:) = image_delta_x_upd__;
image_delta_y_upd_hMi___(1+nvolume,:,:) = image_delta_y_upd__;
flag_image_delta_upd_hMi___(1+nvolume,:,:) = flag_image_delta_upd__;
image_I_value_hMi___(1+nvolume,:,:) = image_I_value__;
image_X_value_hMi___(1+nvolume,:,:) = image_X_value__;
image_S_index_hMi___(1+nvolume,:,:) = image_S_index__;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nvolume=0:n_volume-1;

if (verbose); disp(sprintf(' %% [finished pm_align_M_k_p_to_a_het_k_Y_0]')); end;
