function ...
[X_best_...
,a_UX_Y_0lsq__...
,euler_polar_a__...
,euler_azimu_b__...
,euler_gamma_z__...
,image_delta_x__...
,image_delta_y__...
,image_X_value__...
,image_S_index__...
] = ...
ampm_3(...
 rseed...
,n_iteration...
,n_iteration_register...
,viewing_k_eq_d...
,n_order...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_k_p_r_...
,weight_2d_k_p_r_...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,FTK...
,n_w_...
,pm_n_UX_rank...
,UX__...
,X_weight_r_...
,n_M...
,M_k_p__...
,M_k_q__...
,CTF_index_...
,CTF_k_p__...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
,l_max_...
,a_UX_Y_true_...
,euler_polar_a_...
,euler_azimu_b_...
,euler_gamma_z_...
,image_delta_x_...
,image_delta_y_...
,flag_MS_vs_SM...
,f_rand...
);
%%%%%%%%;
% simple alternating minimization using principal modes. ;
% displacements drawn from disc of radius delta_r_max. ;
% ;
% Input: ;
% rseed: integer random seed. ;
% n_iteration: integer number of iterations for alternating minimization. ;
% n_iteration_register: integer number of iterations to wait before registering result against a_k_Y_true_. ;
% viewing_k_eq_d: real equatorial distance for building viewing angles (i.e., each viewing angle is associated with a principal-template). ;
% n_order: integer polynomial interpolation order used for conjugate-gradient interpolation operator. ;
% n_k_p_r: integer number of shells. ;
% k_p_r_: real array of size n_k_p_r. k_p_r for each shell. ;
% k_p_r_max: real maximum radius (of sphere) to be integrated. ;
% weight_k_p_r_: real array of size n_k_p_r. radial quadrature weighting for sphere. ;
% weight_2d_k_p_r_: real array of size n_k_p_r. radial quadrature weighting for disk. ;
% delta_r_max: real number maximum displacement. ;
% svd_eps: real number FTK svd error. ;
% n_delta_v_requested: integer number of displacements requested (actual number of displacements will be stored in FTK.n_delta_v). ;
% FTK: structure produced by ampmh_FTK_1.m ;
% n_w_: integer array of size n_k_p_r. n_w_max = max(n_w_) is the number of inplane_gamma_z values recorded at that ring. ;
% pm_n_UX_rank: integer number of principal-modes used in UX__. ;
% UX__: real array of size (n_w_max,pm_n_UX_rank). princpal-mode right-singular-vectors. ;
% X_weight_r_: real array of size n_w_max. variance scaling (weight) for different principle-image rings. ;
% n_M: integer number of images. ;
% M_k_p__: complex array of size (n_w_sum,n_M). stack of images in k_p_ format. ;
% M_k_q__: complex array of size (n_w_sum,n_M). stack of images in k_q_ format. ;
%          will calculate if empty. ;
% (unused): CTF_index_: integer array of size n_M. CTF_index_(1+nM) is the CTF_index used for image M_k_p__(:,1+nM). ;
%             This can be empty or set to 1, in which case the same CTF_k_p_ will be used for each image. ;
% (unused): CTF_k_p__: complex array of size(n_w_sum,n_CTF). stack of ctf-functions in k_p_ format. ;
%            If CTF_index_ is empty or set to 1, then we assume this contains only a single CTF_k_p_, ;
%            which will then be used for all images. ;
% svd_VUXM_lwnM____: complex array of size (FTK.n_svd_l,n_w_max,pm_n_UX_rank,n_M): contains combination of FTK.V_r_ and princpal-weights UX and images M. ;
% UX_M_l2_dM__: complex array of size (FTK.n_delta_v,n_M): contains image-norms for each displacement in FTK. ;
% l_max_: integer array of size n_k_p_r. l_max_(1+nk_p_r) is the order used for a_UX_Y_0lsq_ for each principal mode. ;
% a_UX_Y_true_: complex array of size pm_n_lm_sum = sum((1+l_max_).^2). The 'ground-truth' to which a_UX_Y_0lsq_ is compared. ;
% euler_polar_a_: real array of size n_M. initial polar_a used for each image (randomized if empty). ;
% euler_azimu_b_: real array of size n_M. initial azimu_b used for each image (randomized if empty). ;
% euler_gamma_z_: real array of size n_M. initial gamma_z used for each image (randomized if empty). ;
% image_delta_x_: real array of size n_M. initial delta_x used for each image (randomized if empty). ;
% image_delta_y_: real array of size n_M. initial delta_y used for each image (randomized if empty). ;
% flag_MS_vs_SM: integer 0 or 1. ;
%                if set to 1 will enforce uniform distribution of viewing angles. ;
%                if set to 0 will match templates to images. ;
%                default 1. ;
% f_rand: real number. fraction of best templates to select from when matching templates to images. ;
%         default 0.15. ;
% ;
% Output: ;
% X_best_: real array of size n_iteration. correlation with a_k_Y_true_ at every n_iteration_register iterations. ;
% a_UX_Y_0lsq__: complex array of size (pm_n_lm_sum,n_iteration). stack of output functions in k_Y_ format. ;
% euler_polar_a__: real array of size (n_M,n_iteration). list of polar_a used for each image at each iteration. ;
% euler_azimu_b__: real array of size (n_M,n_iteration). list of azimu_b used for each image at each iteration. ;
% euler_gamma_z__: real array of size (n_M,n_iteration). list of gamma_z used for each image at each iteration. ;
% image_delta_x__: real array of size (n_M,n_iteration). list of delta_x used for each image at each iteration. ;
% image_delta_y__: real array of size (n_M,n_iteration). list of delta_y used for each image at each iteration. ;
% image_X_value__: real array of size (n_M,n_iteration). list of X_S (i.e., template-image correlation) used for each image at each iteration. ;
% image_S_index__: real array of size (n_M,n_iteration). list of S_index (i.e., template index) used for each image at each iteration. ;
%%%%%%%%;

verbose=1;
if isempty(flag_MS_vs_SM); flag_MS_vs_SM=1; end;
if isempty(f_rand); f_rand = 0.15; end;

rng(rseed);

l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

pm_n_order = n_order;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ;
pm_CTF_index_ = 1; pm_CTF_k_p__ = ones(pm_n_w_sum,1);

if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
end;%if isempty(FTK);

if isempty(M_k_q__);
tmp_t = tic();
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;
end;%if isempty(M_k_q__);

if isempty(svd_VUXM_lwnM____);
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
end;%if isempty(svd_VUXM_lwnM____);

if isempty(UX_M_l2_dM__);
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
end;%if isempty(UX_M_l2_dM__);

X_best_ = zeros(n_iteration,1);
if (nargout>1);
euler_polar_a__ = zeros(n_M,n_iteration);
euler_azimu_b__ = zeros(n_M,n_iteration);
euler_gamma_z__ = zeros(n_M,n_iteration);
image_delta_x__ = zeros(n_M,n_iteration);
image_delta_y__ = zeros(n_M,n_iteration);
image_X_value__ = zeros(n_M,n_iteration);
image_S_index__ = zeros(n_M,n_iteration);
a_UX_Y_0lsq__ = zeros(pm_n_lm_sum,n_iteration);
end;%if (nargout>1);
%%%%%%%%;
% initialize current euler-angles randomly. ;
%%%%%%%%;
if isempty(euler_polar_a_); euler_polar_a_ = 1*pi*rand(n_M,1); end;
if isempty(euler_azimu_b_); euler_azimu_b_ = 2*pi*rand(n_M,1); end;
if isempty(euler_gamma_z_); euler_gamma_z_ = 2*pi*rand(n_M,1); end;
if isempty(image_delta_x_); image_delta_x_ = zeros(n_M,1); end;
if isempty(image_delta_y_); image_delta_y_ = zeros(n_M,1); end;
for niteration=0:n_iteration-1;
%%%%%%%%;
% Find out which on-grid displacements correspond most closely to current displacements. ;
% Then use current displacements to form principal-images. ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M,svd_VUXM_lwnM____,image_delta_x_,image_delta_y_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
a_UX_Y_0lsq_ = ...
cg_lsq_3( ...
 pm_n_order ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___,[n_w_max*pm_n_k_p_r,n_M]) ...
,pm_CTF_index_ ...
,pm_CTF_k_p__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% cg_lsq_3 for a_UX_Y_0lsq_: %0.3fs',tmp_t)); end;
if (nargout>1);
a_UX_Y_0lsq__(:,1+niteration) = a_UX_Y_0lsq_;
euler_polar_a__(:,1+niteration) = euler_polar_a_;
euler_azimu_b__(:,1+niteration) = euler_azimu_b_;
euler_gamma_z__(:,1+niteration) = euler_gamma_z_;
image_delta_x__(:,1+niteration) = image_delta_x_;
image_delta_y__(:,1+niteration) = image_delta_y_;
end;%if (nargout>1);
%%%%%%%%;
if (mod(niteration,n_iteration_register)==n_iteration_register-1);
%%%%%%%%;
% Compare current model to a_UX_Y_true_. ;
%%%%%%%%;
tmp_t = tic();
[X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(pm_n_k_p_r,pm_k_p_r_,pm_k_p_r_max,pm_weight_k_p_r_,0,pm_l_max_,a_UX_Y_true_,a_UX_Y_0lsq_);
[X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(pm_n_k_p_r,pm_k_p_r_,pm_k_p_r_max,pm_weight_k_p_r_,0,pm_l_max_,a_UX_Y_true_,flipY(pm_n_k_p_r,pm_l_max_,a_UX_Y_0lsq_));
X_best = max(real(X_best_orig),real(X_best_flip));
X_best_(1+niteration) = X_best;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% register_spharm_to_spharm_wigner_0: %0.3fs',tmp_t)); end;
disp(sprintf(' %% a_UX_Y_true_ vs a_UX_Y_lsq0_: correlation %+0.6f',X_best));
end;%if (mod(niteration,n_iteration_register)==n_iteration_register-1);
%%%%%%%%;
% use current model to generate current principal-templates. ;
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
[tmp_S_k_p__,~,~,~,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,pm_n_k_p_r,pm_k_p_r_,pm_k_p_r_max,pm_weight_k_p_r_,pm_l_max_,a_UX_Y_0lsq_,viewing_k_eq_d,-1,pm_n_w_);
if (tmp_verbose>0); disp(sprintf(' %% viewing_k_eq_d %0.3f, n_viewing_all %d',viewing_k_eq_d,n_viewing_all)); end;
n_S = n_viewing_all;
UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
UX_S_l2_(1+nS) = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
tmp_S_k_q__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
tmp_S_k_q__(:,1+nS) = interp_p_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS)); 
end;%for nS=0:n_S-1; 
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_template_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
tmp_t = tic();
SS_k_q_ = svd(tmp_S_k_q__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(tmp_S_k_q__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(tmp_S_k_q__,n_S_rank);
if (verbose); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
% Batches images into batches of size n_M_Mbatch. ;
% Batches templates into batches of size n_S_Sbatch. ;
% Only stores the optimal translation for each image. ;
%%%%%%%%;
tmp_t = tic();
n_M_Mbatch = 24;
n_S_Sbatch = 24;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___1(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_Sbatch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M...
,n_M_Mbatch...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
tmp_t = tic();
[...
 euler_polar_a_...
,euler_azimu_b_...
,euler_gamma_z_...
,image_delta_x_...
,image_delta_y_...
,image_X_value_...
,image_S_index_...
] = ...
ampmh_MS_0(...
 n_w_...
,n_S...
,viewing_polar_a_all_...
,viewing_azimu_b_all_...
,n_M...
,X_wSM___...
,FTK.delta_x_(1+delta_j_wSM___)...
,FTK.delta_y_(1+delta_j_wSM___)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% MS: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
%%%%%%%%;
if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
tmp_t = tic();
[...
 euler_polar_a_...
,euler_azimu_b_...
,euler_gamma_z_...
,image_delta_x_...
,image_delta_y_...
,image_X_value_...
,image_S_index_...
] = ...
ampmh_SM_0(...
 f_rand...
,n_w_...
,n_S...
,viewing_polar_a_all_...
,viewing_azimu_b_all_...
,n_M...
,X_wSM___...
,FTK.delta_x_(1+delta_j_wSM___)...
,FTK.delta_y_(1+delta_j_wSM___)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% SM: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
if (nargout>1);
if niteration<n_iteration-1;
image_X_value__(:,1+niteration+1) = image_X_value_;
image_S_index__(:,1+niteration+1) = image_S_index_;
end;%if niteration<n_iteration-1;
end;%if (nargout>1);

%%%%%%%%;
% Now return to beginning of loop. ;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;

