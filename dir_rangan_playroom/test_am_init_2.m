% attempts volume-alignment after random sampling. ;
% accounts for some level of translation. ;
verbose=1;


date_diff_threshold = 1.5;
n_S = n_viewing_all;
n_M = n_image_sub;

delta_sigma_use = delta_sigma;
svd_eps = 1e-2;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.00,0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;

str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
ndat_rseed=0;%for ndat_rseed=0:n_dat_rseed-1;
am_init_rseed = dat_rseed_(1+ndat_rseed);

n_w_uni_ = n_w_max*ones(n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(n_w_uni_),:);
init_a_CTF_UX_Y_quad__ = zeros(n_lm_max,dat_n_UX_rank);
for dat_nUX_rank=0:dat_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = (l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
init_a_CTF_UX_Y_quad__(1:n_lm,1+dat_nUX_rank) = ...
init_a_CTF_UX_Y_quad__(1:n_lm,1+dat_nUX_rank) ...
+ UX__(1+nk_p_r,1+dat_nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r);
%<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
dat_f_rand = 0.05; dat_flag_plot=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for am_init_nUX_rank=0:dat_n_UX_rank-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_am_init = sprintf('am_init_nUX%.3drng%.3d',am_init_nUX_rank,am_init_rseed);
infix_am_init = sprintf('%s_%s',dat_infix,str_am_init);

fname_mat_A = sprintf('%s_mat/%s_A.mat',dir_trunk,infix_am_init);
fname_mat_B = sprintf('%s_mat/%s_B.mat',dir_trunk,infix_am_init);
fname_mat_C = sprintf('%s_mat/%s_C.mat',dir_trunk,infix_am_init);
fname_tmp_ABC = sprintf('%s_mat/%s_ABC.tmp',dir_trunk,infix_am_init);

flag_exist = ( exist(fname_mat_A,'file') &  exist(fname_mat_B,'file') &  exist(fname_mat_C,'file'));
if ( flag_exist);
disp(sprintf(' %% %s found, not creating',fname_mat_A));
disp(sprintf(' %% %s found, not creating',fname_mat_B));
disp(sprintf(' %% %s found, not creating',fname_mat_C));
flag_skip=0;
end;%if ( flag_exist);
if (~flag_exist);
flag_skip=0;
if ( exist(fname_tmp_ABC,'file'));
tmp_date_diff = datenum(clock) - datenum(dir(fname_tmp_ABC).date);
if (tmp_date_diff< date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = recent, skipping',fname_tmp_ABC,tmp_date_diff));
flag_skip=1;
end;%if (tmp_date_diff< date_diff_threshold);
if (tmp_date_diff>=date_diff_threshold);
disp(sprintf(' %% %s found, tmp_date_diff = %0.2f = stale, deleting',fname_tmp_ABC,tmp_date_diff));
delete(fname_tmp_ABC);
flag_skip=0;
end;%if (tmp_date_diff>=date_diff_threshold);
end;%if ( exist(fname_tmp_ABC,'file'));
end;%if (~flag_exist);

if (~flag_skip);
save(fname_tmp_ABC,'fname_tmp_ABC');

delta_r_max = delta_r_max_use;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d, n_delta_v_requested %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l,n_delta_v_requested));

n_w_uni_max = max(n_w_uni_);
n_w_uni_sum = sum(n_w_uni_);
n_w_uni_csum_ = cumsum([0;n_w_uni_]);

pm_n_UX_rank = dat_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 

init_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
init_n_order = dat_n_order;
init_pm_n_UX_rank = 1+am_init_nUX_rank;
init_pm_n_k_p_r = init_pm_n_UX_rank;
init_pm_k_p_r_ = ones(init_pm_n_k_p_r,1);
init_pm_k_p_r_max = 1;
init_pm_n_w_ = n_w_max*ones(init_pm_n_k_p_r,1);
init_pm_n_w_max = n_w_max;
init_pm_n_w_sum = sum(init_pm_n_w_);
init_pm_n_w_csum_ = cumsum([0;init_pm_n_w_]);
init_pm_l_max_ = l_max_max*ones(init_pm_n_k_p_r,1);
init_pm_n_lm_ = (1+init_pm_l_max_).^2; init_pm_n_lm_sum = sum(init_pm_n_lm_);
init_pm_weight_3d_k_p_r_ = ones(init_pm_n_k_p_r,1);
init_pm_weight_2d_k_p_r_ = ones(init_pm_n_k_p_r,1);
flag_MS_vs_SM = 1;
init_a_CTF_UX_Y_true_ = reshape(init_a_CTF_UX_Y_quad__(:,1:init_pm_n_UX_rank),[n_lm_max*init_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;

[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_idx_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;

%%%%%%%%;
% Prepare for order-limit. ;
%%%%%%%%;
init_Y_l_val_ = zeros(init_pm_n_lm_sum,1);
init_Y_m_val_ = zeros(init_pm_n_lm_sum,1);
na=0;
for nk_p_r=0:init_pm_n_k_p_r-1;
l_max = init_pm_l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
init_Y_l_val_(1+na) = l_val;
init_Y_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:init_pm_n_k_p_r-1;

%%%%%%%%;
% Construct M_k_q__ without taking into account the translations. ;
%%%%%%%%;
disp(sprintf(' %% not-centering the images. '));
tmp_t = tic();
dat_M_k_q__ = zeros(n_w_uni_sum,n_M);
for nM=0:n_M-1;
dat_M_k_p_ = ...
transf_p_to_p( ...
 n_k_p_r ...
,k_p_r_ ...
,n_w_uni_ ...
,n_w_uni_sum ...
,dat_M_k_p__(:,1+nM) ...
,+0.0*delta_read_x_(1+nM) ...
,+0.0*delta_read_y_(1+nM) ...
);
dat_M_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_uni_ ...
,n_w_uni_sum ...
,dat_M_k_p_ ...
);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now form svd_VUXM_lwnM____ using these translated images. ;
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_uni_,n_M,dat_M_k_q__,init_pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the non-translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_uni_,n_M,init_pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now, form principal-images (using the zero-displacement). ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_uni_,init_pm_n_UX_rank,n_M,svd_VUXM_lwnM____,zeros(n_M,1),zeros(n_M,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;

init_Y_lm_cut_threshold_ = [2,4,6,8,12,24];%init_Y_lm_cut_threshold_ = 1:1:1+2*11; %init_Y_lm_cut_threshold_ = 1:init_pm_l_max_max; 
n_init_Y_lm_cut_threshold = numel(init_Y_lm_cut_threshold_);
n_iteration = 16;
init_a_CTF_UX_Y_0lsq_yncit____ = zeros(init_pm_n_lm_sum,n_CTF_rank,n_iteration,n_init_Y_lm_cut_threshold);
init_a_CTF_UX_Y_0lsq_Yit___ = zeros(init_pm_n_lm_sum,n_iteration,n_init_Y_lm_cut_threshold);
init_X_best_it__ = zeros(n_iteration,n_init_Y_lm_cut_threshold);
init_euler_polar_a_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_euler_azimu_b_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_euler_gamma_z_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_delta_x_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_delta_y_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_I_value_MS_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_euler_polar_a_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_euler_azimu_b_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_euler_gamma_z_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_delta_x_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_delta_y_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);
init_image_I_value_SM_Mit___ = zeros(n_M,n_iteration,n_init_Y_lm_cut_threshold);

if ( exist(fname_mat_A,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_A)); end;
if (~exist(fname_mat_A,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_A));
%%%%%%%%;
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
for niteration=0:n_iteration-1;
if (verbose>0); disp(sprintf(' %% ninit_Y_lm_cut_threshold %d/%d; niteration %d/%d',ninit_Y_lm_cut_threshold,n_init_Y_lm_cut_threshold,niteration,n_iteration)); end;
%init_Y_lm_cut_ = init_Y_l_val_+abs(init_Y_m_val_)<=init_Y_lm_cut_threshold;
init_Y_lm_cut_ = init_Y_l_val_+0*abs(init_Y_m_val_)<=init_Y_lm_cut_threshold;
flag_true_vs_rand=0;
if (verbose>1); disp(sprintf(' %% flag_true_vs_rand %d',flag_true_vs_rand)); end;
if flag_true_vs_rand==0;
%%%%%%%%;
% initialize current euler-angles using 'true' angles. ;
%%%%%%%%;
euler_polar_a_ = +euler_angle_marina_(1,1:n_M);
euler_azimu_b_ = +euler_angle_marina_(2,1:n_M);
euler_gamma_z_ = -euler_angle_marina_(3,1:n_M);
image_delta_x_ = +1.0*delta_read_x_(1:n_M);
image_delta_y_ = +1.0*delta_read_y_(1:n_M);
image_I_value_ = ones(n_M,1);
end;%if flag_true_vs_rand==1;
if flag_true_vs_rand==0;
%%%%%%%%;
% initialize current euler-angles using random angles. ;
% do not center. ;
%%%%%%%%;
rng(1024*am_init_rseed + 729*niteration);
euler_polar_a_ = 1*pi*rand(n_M,1);
euler_azimu_b_ = 2*pi*rand(n_M,1);
euler_gamma_z_ = 2*pi*rand(n_M,1);
image_delta_x_ = +0.0*delta_read_x_(1:n_M); %image_delta_x_ = zeros(n_M,1);
image_delta_y_ = +0.0*delta_read_y_(1:n_M); %image_delta_y_ = zeros(n_M,1);
image_I_value_ = ones(n_M,1);
end;%if flag_true_vs_rand==0;
%%%%%%%%;
% use current euler-angles and displacements to solve for current model (i.e., a random function). ;
%%%%%%%%;
tmp_t = tic();
init_a_UCTF_UX_Y_0lsq_ync__ = ...
cg_lsq_pm_0( ...
 init_n_order ...
,init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_l_max_ ...
,init_pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___(:,1:init_pm_n_k_p_r,:),[n_w_uni_max*init_pm_n_k_p_r,n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
,image_I_value_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_pm_0 for init_a_UCTF_UX_Y_0lsq_ync__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now normalize init_a_CTF_UX_Y_0lsq_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
init_a_UCTF_UX_Y_0lsq_ync__ = spharm__normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,init_a_UCTF_UX_Y_0lsq_ync__);
%%%%%%%%;
% use init_a_UCTF_UX_Y_0lsq_ync__ as well VSCTF_Mc__ to approximate the image-averaged init_a_CTF_UX_Y_0lsq_. ;
%%%%%%%%;
init_a_CTF_UX_Y_0lsq_ = spharm_normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,mean(init_a_UCTF_UX_Y_0lsq_ync__*transpose(VSCTF_Mc__),2));
%init_a_CTF_UX_Y_0lsq_Yit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = init_a_CTF_UX_Y_0lsq_; %<-- do not store, will recalculate later. ;
%%%%%%%%;
% Compare current model to init_a_CTF_UX_Y_true_. ;
%%%%%%%%;
tmp_t = tic();
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,init_a_CTF_UX_Y_0lsq_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,flipY(init_pm_n_k_p_r,init_pm_l_max_,init_a_CTF_UX_Y_0lsq_));
tmp_X_best = max(real(tmp_X_best_orig),real(tmp_X_best_flip));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% register_spharm_to_spharm_wigner_0: %0.3fs',tmp_t)); end;
if (verbose>-1); disp(sprintf(' %% init_a_CTF_UX_Y_true_ vs init_a_CTF_UX_Y_0lsq_: correlation %+0.6f',tmp_X_best)); end;
%%%%%%%%;
% bandlimit init_a_UCTF_UX_Y_0lsq_ync__. ;
%%%%%%%%;
blim_a_UCTF_UX_Y_0lsq_ync__ = init_a_UCTF_UX_Y_0lsq_ync__.*repmat(init_Y_lm_cut_,[1,n_CTF_rank]);
%%%%%%%%;
% Groups images by micrograph. ;
% Calculates principal-templates associated with each micrograph. ;
% Takes svd of principal-templates. ;
% Batches images into batches of size n_M_Mbatch. ;
% Batches templates into batches of size n_S_Sbatch. ;
% Only stores the optimal translation for each image. ;
%%%%%%%%;
tmp_t = tic();
n_M_Mbatch = 24;
n_S_Sbatch = 24;
[ ...
 n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,X_wSM___ ...
,delta_j_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___5( ...
 FTK ...
,init_viewing_k_eq_d ...
,n_w_max ...
,l_max_max ...
,init_pm_n_UX_rank ...
,blim_a_UCTF_UX_Y_0lsq_ync__ ...
,CTF_idx_ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,n_S_Sbatch ...
,n_M ...
,n_M_Mbatch ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
flag_MS_vs_SM=1;
if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
tmp_t = tic();
[ ...
 euler_polar_a_MS_ ...
,euler_azimu_b_MS_ ...
,euler_gamma_z_MS_ ...
,image_delta_x_upd_MS_ ...
,image_delta_y_upd_MS_ ...
,image_I_value_MS_ ...
,image_X_value_MS_ ...
,image_S_index_MS_ ...
] = ...
ampmh_MS_1( ...
 n_w_max ...
,n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,n_M ...
,X_wSM___ ...
,FTK.delta_x_(1+delta_j_wSM___) ...
,FTK.delta_y_(1+delta_j_wSM___) ...
,I_value_wSM___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% MS: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
image_delta_x_MS_ = image_delta_x_ + image_delta_x_upd_MS_;
image_delta_y_MS_ = image_delta_y_ + image_delta_y_upd_MS_;
flag_MS_vs_SM=0;
f_rand = dat_f_rand;
if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
tmp_t = tic();
[ ...
 euler_polar_a_SM_ ...
,euler_azimu_b_SM_ ...
,euler_gamma_z_SM_ ...
,image_delta_x_upd_SM_ ...
,image_delta_y_upd_SM_ ...
,image_I_value_SM_ ...
,image_X_value_SM_ ...
,image_S_index_SM_ ...
] = ...
ampmh_SM_1( ...
 f_rand ...
,n_w_max ...
,n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,n_M ...
,X_wSM___ ...
,FTK.delta_x_(1+delta_j_wSM___) ...
,FTK.delta_y_(1+delta_j_wSM___) ...
,I_value_wSM___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% SM: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
image_delta_x_SM_ = image_delta_x_ + image_delta_x_upd_SM_;
image_delta_y_SM_ = image_delta_y_ + image_delta_y_upd_SM_;
%%%%%%%%;
% Now, with these updated MS euler-angles, recalculate the (now nonrandom) current model. ;
% Note that we use displacments, but still set the intensity to 1. ;
%%%%%%%%;
tmp_t = tic();
init_a_UCTF_UX_Y_0lsq_ync__ = ...
cg_lsq_pm_0( ...
 init_n_order ...
,init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_l_max_ ...
,init_pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___(:,1:init_pm_n_k_p_r,:),[n_w_uni_max*init_pm_n_k_p_r,n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_MS_ ...
,euler_azimu_b_MS_ ...
,euler_gamma_z_MS_ ...
,image_delta_x_MS_ ...
,image_delta_y_MS_ ...
,image_I_value_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_pm_0 for init_a_UCTF_UX_Y_0lsq_ync__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now normalize init_a_CTF_UX_Y_0lsq_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
init_a_UCTF_UX_Y_0lsq_ync__ = spharm__normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,init_a_UCTF_UX_Y_0lsq_ync__);
init_a_CTF_UX_Y_0lsq_yncit____(:,:,1+niteration,1+ninit_Y_lm_cut_threshold) = init_a_UCTF_UX_Y_0lsq_ync__; %<-- store. ;
%%%%%%%%;
% use init_a_UCTF_UX_Y_0lsq_ync__ as well VSCTF_Mc__ to approximate the image-averaged init_a_CTF_UX_Y_0lsq_. ;
%%%%%%%%;
init_a_CTF_UX_Y_0lsq_ = spharm_normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,mean(init_a_UCTF_UX_Y_0lsq_ync__*transpose(VSCTF_Mc__),2));
init_a_CTF_UX_Y_0lsq_Yit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = init_a_CTF_UX_Y_0lsq_; %<-- store. ;
%%%%%%%%;
% Compare current model to init_a_CTF_UX_Y_true_. ;
%%%%%%%%;
tmp_t = tic();
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,init_a_CTF_UX_Y_0lsq_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,flipY(init_pm_n_k_p_r,init_pm_l_max_,init_a_CTF_UX_Y_0lsq_));
tmp_X_best = max(real(tmp_X_best_orig),real(tmp_X_best_flip));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% register_spharm_to_spharm_wigner_0: %0.3fs',tmp_t)); end;
if (verbose>-1); disp(sprintf(' %% init_a_CTF_UX_Y_true_ vs init_a_CTF_UX_Y_0lsq_: correlation %+0.6f',tmp_X_best)); end;
init_X_best_it__(1+niteration,1+ninit_Y_lm_cut_threshold) = tmp_X_best;
%%%%%%%%;
init_euler_polar_a_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_polar_a_MS_ ;
init_euler_azimu_b_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_azimu_b_MS_ ;
init_euler_gamma_z_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_gamma_z_MS_ ;
init_image_delta_x_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_delta_x_MS_ ;
init_image_delta_y_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_delta_y_MS_ ;
init_image_I_value_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_I_value_MS_ ;
init_euler_polar_a_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_polar_a_SM_ ;
init_euler_azimu_b_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_azimu_b_SM_ ;
init_euler_gamma_z_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = euler_gamma_z_SM_ ;
init_image_delta_x_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_delta_x_SM_ ;
init_image_delta_y_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_delta_y_SM_ ;
init_image_I_value_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = image_I_value_SM_ ;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
save(fname_mat_A ...
     ,'n_iteration' ...
     ,'n_init_Y_lm_cut_threshold' ...
     ,'init_Y_lm_cut_threshold_' ...
     ,'init_a_CTF_UX_Y_0lsq_yncit____' ...
     ,'init_a_CTF_UX_Y_0lsq_Yit___' ...
     ,'init_X_best_it__' ...
     ,'init_euler_polar_a_MS_Mit___' ...
     ,'init_euler_azimu_b_MS_Mit___' ...
     ,'init_euler_gamma_z_MS_Mit___' ...
     ,'init_image_delta_x_MS_Mit___' ...
     ,'init_image_delta_y_MS_Mit___' ...
     ,'init_image_I_value_MS_Mit___' ...
     ,'init_euler_polar_a_SM_Mit___' ...
     ,'init_euler_azimu_b_SM_Mit___' ...
     ,'init_euler_gamma_z_SM_Mit___' ...
     ,'init_image_delta_x_SM_Mit___' ...
     ,'init_image_delta_y_SM_Mit___' ...
     ,'init_image_I_value_SM_Mit___' ...
     );
end;%if (~exist(fname_mat_A,'file')); 
load(fname_mat_A);

if ( exist(fname_mat_B,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_B)); end;
if (~exist(fname_mat_B,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_B));
%%%%%%%%;
% Now align the init_a_CTF_UX_Y_0lsq_Yit___(:,:,1+ninit_Y_lm_cut_threshold);
%%%%%%%%;
init_a_CTF_UX_Y_0lsq_ynct___ = zeros(init_pm_n_lm_sum,n_CTF_rank,n_init_Y_lm_cut_threshold);
init_a_CTF_UX_Y_0lsq_flip_flag__ = zeros(n_iteration,n_init_Y_lm_cut_threshold);
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
if (verbose>0); disp(sprintf(' %% ninit_Y_lm_cut_threshold %d/%d;',ninit_Y_lm_cut_threshold,n_init_Y_lm_cut_threshold)); end;
[ ...
 tmp_a_ync_ ...
, ~ ...
, init_a_CTF_UX_Y_0lsq_flip_flag__(:,1+ninit_Y_lm_cut_threshold) ...
] = ...
am3d_0( ...
 init_pm_n_k_p_r*n_CTF_rank ...
,repmat(init_pm_k_p_r_,[n_CTF_rank,1]) ...
,init_pm_k_p_r_max ...
,repmat(init_pm_weight_3d_k_p_r_,[n_CTF_rank,1]) ...
,n_iteration ...
,repmat(init_pm_l_max_,[n_CTF_rank,1]) ...
,reshape(init_a_CTF_UX_Y_0lsq_yncit____(:,:,:,1+ninit_Y_lm_cut_threshold),[init_pm_n_lm_sum*n_CTF_rank,n_iteration]) ...
);
init_a_CTF_UX_Y_0lsq_ynct___(:,:,1+ninit_Y_lm_cut_threshold) = reshape(tmp_a_ync_,[init_pm_n_lm_sum,n_CTF_rank]);
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
%%%%%%%%;
save(fname_mat_B ...
     ,'init_a_CTF_UX_Y_0lsq_ynct___' ...
     ,'init_a_CTF_UX_Y_0lsq_flip_flag__' ...
     );
end;%if (~exist(fname_mat_B,'file')); 
load(fname_mat_B);

if ( exist(fname_mat_C,'file')); disp(sprintf(' %% %s found, not creating',fname_mat_C)); end;
if (~exist(fname_mat_C,'file')); disp(sprintf(' %% %s not found, creating',fname_mat_C));
alig_a_CTF_UX_Y_0lsq_ynct___ = zeros(init_pm_n_lm_sum,n_CTF_rank,n_init_Y_lm_cut_threshold);
alig_a_CTF_UX_Y_0lsq_Yt__ = zeros(init_pm_n_lm_sum,n_init_Y_lm_cut_threshold);
alig_X_best_t_ = zeros(n_init_Y_lm_cut_threshold,1);
alig_euler_polar_a_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_euler_azimu_b_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_euler_gamma_z_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_delta_x_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_delta_y_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_I_value_MS_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_euler_polar_a_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_euler_azimu_b_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_euler_gamma_z_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_delta_x_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_delta_y_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
alig_image_I_value_SM_Mt__ = zeros(n_M,n_init_Y_lm_cut_threshold);
%%%%%%%%;
% Check to see if aligned function is any better than the (order-limited) random initializations. ;
%%%%%%%%;
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
if (verbose>0); disp(sprintf(' %% ninit_Y_lm_cut_threshold %d/%d;',ninit_Y_lm_cut_threshold,n_init_Y_lm_cut_threshold)); end;
%%%%%%%%;
% First perform alignment of the images to the aligned function. ;
%%%%%%%%;
tmp_t = tic();
n_M_Mbatch = 24;
n_S_Sbatch = 24;
[ ...
 n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,X_wSM___ ...
,delta_j_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___5( ...
 FTK ...
,init_viewing_k_eq_d ...
,n_w_max ...
,l_max_max ...
,init_pm_n_UX_rank ...
,init_a_CTF_UX_Y_0lsq_ynct___(:,:,1+ninit_Y_lm_cut_threshold) ...
,CTF_idx_ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,n_S_Sbatch ...
,n_M ...
,n_M_Mbatch ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
flag_MS_vs_SM=1;
if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
tmp_t = tic();
[ ...
 euler_polar_a_MS_ ...
,euler_azimu_b_MS_ ...
,euler_gamma_z_MS_ ...
,image_delta_x_upd_MS_ ...
,image_delta_y_upd_MS_ ...
,image_I_value_MS_ ...
,image_X_value_MS_ ...
,image_S_index_MS_ ...
] = ...
ampmh_MS_1( ...
 n_w_max ...
,n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,n_M ...
,X_wSM___ ...
,FTK.delta_x_(1+delta_j_wSM___) ...
,FTK.delta_y_(1+delta_j_wSM___) ...
,I_value_wSM___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% MS: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==1; %<-- assign images to templates, ensuring a uniform distribution of viewing angles. ;
image_delta_x_MS_ = image_delta_x_ + image_delta_x_upd_MS_;
image_delta_y_MS_ = image_delta_y_ + image_delta_y_upd_MS_;
flag_MS_vs_SM=0;
f_rand = dat_f_rand;
if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
tmp_t = tic();
[ ...
 euler_polar_a_SM_ ...
,euler_azimu_b_SM_ ...
,euler_gamma_z_SM_ ...
,image_delta_x_upd_SM_ ...
,image_delta_y_upd_SM_ ...
,image_I_value_SM_ ...
,image_X_value_SM_ ...
,image_S_index_SM_ ...
] = ...
ampmh_SM_1( ...
 f_rand ...
,n_w_max ...
,n_S ...
,tmp_viewing_polar_a_all_ ...
,tmp_viewing_azimu_b_all_ ...
,n_M ...
,X_wSM___ ...
,FTK.delta_x_(1+delta_j_wSM___) ...
,FTK.delta_y_(1+delta_j_wSM___) ...
,I_value_wSM___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% SM: update euler_polar_a_ euler_azimu_b_ euler_gamma_z_ : %0.3fs',tmp_t)); end;
end;%if flag_MS_vs_SM==0; %<-- assign templates to images (uses leslie f_rand). ;
image_delta_x_SM_ = image_delta_x_ + image_delta_x_upd_SM_;
image_delta_y_SM_ = image_delta_y_ + image_delta_y_upd_SM_;
%%%%%%%%;
% Now, with these updated MS euler-angles, recalculate the (now nonrandom) current model. ;
% Note that we still set the displacements to zero and the intensity to 1. ;
%%%%%%%%;
tmp_t = tic();
alig_a_UCTF_UX_Y_0lsq_ync__ = ...
cg_lsq_pm_0( ...
 init_n_order ...
,init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_l_max_ ...
,init_pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___(:,1:init_pm_n_k_p_r,:),[n_w_uni_max*init_pm_n_k_p_r,n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_MS_ ...
,euler_azimu_b_MS_ ...
,euler_gamma_z_MS_ ...
,image_delta_x_MS_ ...
,image_delta_y_MS_ ...
,image_I_value_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% cg_lsq_pm_0 for alig_a_UCTF_UX_Y_0lsq_ync__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now normalize alig_a_CTF_UX_Y_0lsq_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
alig_a_UCTF_UX_Y_0lsq_ync__ = spharm__normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,alig_a_UCTF_UX_Y_0lsq_ync__);
alig_a_CTF_UX_Y_0lsq_ynct___(:,:,1+ninit_Y_lm_cut_threshold) = alig_a_UCTF_UX_Y_0lsq_ync__; %<-- store. ;
%%%%%%%%;
% use alig_a_UCTF_UX_Y_0lsq_ync__ as well VSCTF_Mc__ to approximate the image-averaged alig_a_CTF_UX_Y_0lsq_. ;
%%%%%%%%;
alig_a_CTF_UX_Y_0lsq_ = spharm_normalize_1(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_weight_3d_k_p_r_,init_pm_l_max_,mean(alig_a_UCTF_UX_Y_0lsq_ync__*transpose(VSCTF_Mc__),2));
alig_a_CTF_UX_Y_0lsq_Yt__(:,1+ninit_Y_lm_cut_threshold) = alig_a_CTF_UX_Y_0lsq_; %<-- store. ;
%%%%%%%%;
% Compare current model to init_a_CTF_UX_Y_true_. ;
%%%%%%%%;
tmp_t = tic();
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,alig_a_CTF_UX_Y_0lsq_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,flipY(init_pm_n_k_p_r,init_pm_l_max_,alig_a_CTF_UX_Y_0lsq_));
tmp_X_best = max(real(tmp_X_best_orig),real(tmp_X_best_flip));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% register_spharm_to_spharm_wigner_0: %0.3fs',tmp_t)); end;
if (verbose>-1); disp(sprintf(' %% init_a_CTF_UX_Y_true_ vs alig_a_CTF_UX_Y_0lsq_: correlation %+0.6f',tmp_X_best)); end;
alig_X_best_t_(1+ninit_Y_lm_cut_threshold) = tmp_X_best;
%%%%%%%%;
alig_euler_polar_a_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_polar_a_MS_ ;
alig_euler_azimu_b_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_azimu_b_MS_ ;
alig_euler_gamma_z_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_gamma_z_MS_ ;
alig_image_delta_x_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_delta_x_MS_ ;
alig_image_delta_y_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_delta_y_MS_ ;
alig_image_I_value_MS_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_I_value_MS_ ;
alig_euler_polar_a_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_polar_a_SM_ ;
alig_euler_azimu_b_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_azimu_b_SM_ ;
alig_euler_gamma_z_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = euler_gamma_z_SM_ ;
alig_image_delta_x_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_delta_x_SM_ ;
alig_image_delta_y_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_delta_y_SM_ ;
alig_image_I_value_SM_Mt__(:,1+ninit_Y_lm_cut_threshold) = image_I_value_SM_ ;
%%%%%%%%;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
%%%%%%%%;
save(fname_mat_C ...
     ,'alig_a_CTF_UX_Y_0lsq_ynct___' ...
     ,'alig_a_CTF_UX_Y_0lsq_Yt__' ...
     ,'alig_X_best_t_' ...
     ,'alig_euler_polar_a_MS_Mt__' ...
     ,'alig_euler_azimu_b_MS_Mt__' ...
     ,'alig_euler_gamma_z_MS_Mt__' ...
     ,'alig_image_delta_x_MS_Mt__' ...
     ,'alig_image_delta_y_MS_Mt__' ...
     ,'alig_image_I_value_MS_Mt__' ...
     ,'alig_euler_polar_a_SM_Mt__' ...
     ,'alig_euler_azimu_b_SM_Mt__' ...
     ,'alig_euler_gamma_z_SM_Mt__' ...
     ,'alig_image_delta_x_SM_Mt__' ...
     ,'alig_image_delta_y_SM_Mt__' ...
     ,'alig_image_I_value_SM_Mt__' ...
     );
end;%if (~exist(fname_mat_C,'file')); 
load(fname_mat_C);

k_c_d = @(e0,e1) fnorm([sin(e0(1))*cos(e0(2));sin(e0(1))*sin(e0(2));cos(e0(1))]-[sin(e1(1))*cos(e1(2));sin(e1(1))*sin(e1(2));cos(e1(1))]);
dlim_ = [0,2];
k_c_w = @(e0,e1) acos(1 - k_c_d(e0,e1).^2/2);
wlim_=[0,pi];
init_euler_k_c_d_true__ = zeros(n_M,n_M);
for nM0=0:n_M-1;
for nM1=1+nM0:n_M-1;
tmp_k_c_d = k_c_d(euler_angle_marina_(:,1+nM0),euler_angle_marina_(:,1+nM1));
init_euler_k_c_d_true__(1+nM0,1+nM1) = tmp_k_c_d;
init_euler_k_c_d_true__(1+nM1,1+nM0) = tmp_k_c_d;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
%%%%%%%%;
init_euler_k_c_d_rand_ = zeros(n_M,n_M);
n_h = 24;
init_h2_rand__ = zeros(n_h,n_h);
for niteration=0:n_iteration-1;
for nM0=0:n_M-1; tmp_e0_ = [1*pi*rand;2*pi*rand;2*pi*rand];
for nM1=1+nM0:n_M-1; tmp_e1_ = [1*pi*rand;2*pi*rand;2*pi*rand];
tmp_k_c_d = k_c_d(tmp_e0_,tmp_e1_);
init_euler_k_c_d_rand__(1+nM0,1+nM1) = tmp_k_c_d;
init_euler_k_c_d_rand__(1+nM1,1+nM0) = tmp_k_c_d;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
tmp_ = ones(n_M,n_M)-eye(n_M,n_M);
init_index_offdiag_ = efind(tmp_(:));
init_h2_rand__ = init_h2_rand__ + hist2d_0(init_euler_k_c_d_true__(1+init_index_offdiag_),init_euler_k_c_d_rand__(1+init_index_offdiag_),n_h,n_h,dlim_,dlim_)/n_iteration;
end;%for niteration=0:n_iteration-1;
init_h2_rand_norm__ = init_h2_rand__./repmat(sum(init_h2_rand__,1),[n_h,1]);

flag_plot=1;
if flag_plot;
fname_fig_A = sprintf('%s_jpg/%s_FIGA',dir_trunk,infix_am_init);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig_A),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_A));
%%%%%%%%;
% plot viewing-angle-spread. ;
%%%%%%%%;
figure(4); clf; figbeach(); colormap(colormap_pm(64));
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
subplot(2,3,1+ninit_Y_lm_cut_threshold);
tmp_h2_0lsq__ = zeros(n_h,n_h);
for niteration=0:n_iteration-1;
tmp_euler_k_c_d_0lsq_ = zeros(n_M,n_M);
for nM0=0:n_M-1;
tmp_e0_ = [ ...
 init_euler_polar_a_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_azimu_b_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_gamma_z_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
	    ];
for nM1=1+nM0:n_M-1;
tmp_e1_ = [ ...
 init_euler_polar_a_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_azimu_b_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_gamma_z_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
	    ];
tmp_k_c_d = k_c_d(tmp_e0_,tmp_e1_);
tmp_euler_k_c_d_0lsq__(1+nM0,1+nM1) = tmp_k_c_d;
tmp_euler_k_c_d_0lsq__(1+nM1,1+nM0) = tmp_k_c_d;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
tmp_h2_0lsq__ = tmp_h2_0lsq__ + hist2d_0(init_euler_k_c_d_true__(1+init_index_offdiag_),tmp_euler_k_c_d_0lsq__(1+init_index_offdiag_),n_h,n_h,dlim_,dlim_)/n_iteration;
end;%for niteration=0:n_iteration-1;
tmp_h2_0lsq_norm__ = tmp_h2_0lsq__./repmat(sum(tmp_h2_0lsq__,1),[n_h,1]);
hold on;
imagesc(tmp_h2_0lsq_norm__ - init_h2_rand_norm__,0.1*[-1,+1]);
%imagesc(tmp_h2_0lsq__);
%plot(dlim_,dlim_,'k-');
plot([0,n_h+1],[0,n_h+1],'k-');
hold off;
%axis([0,2,0,2]); axisnotick;
axis([0,n_h+1,0,n_h+1]); axisnotick;
xlabel('true'); ylabel('0lsq');
title(sprintf('MS cut %0.2f',init_Y_lm_cut_threshold));
drawnow;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig_A));
print('-djpeg',sprintf('%s.jpg',fname_fig_A));
print('-depsc',sprintf('%s.eps',fname_fig_A));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig_A),'file'));
if ( exist(sprintf('%s.jpg',fname_fig_A),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_A));
end;%if ( exist(sprintf('%s.jpg',fname_fig_A),'file'));
end;%if flag_plot;

flag_plot=1;
if flag_plot;
fname_fig_B = sprintf('%s_jpg/%s_FIGB',dir_trunk,infix_am_init);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig_B),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_B));
%%%%%%%%;
% plot viewing-angle-spread for those viewing-angles that are close. ;
% accumulatd over iterations. ;
%%%%%%%%;
figure(5); clf; figbeach(); colormap(colormap_pm(64));
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
subplot(2,3,1+ninit_Y_lm_cut_threshold);
tmp_h2_0lsq__ = zeros(n_h,n_h);
for niteration=0:n_iteration-1;
tmp_euler_k_c_d_0lsq_ = zeros(n_M,n_M);
for nM0=0:n_M-1;
tmp_e0_ = [ ...
 init_euler_polar_a_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_azimu_b_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_gamma_z_MS_Mit___(1+nM0,1+niteration,1+ninit_Y_lm_cut_threshold) ...
	    ];
for nM1=1+nM0:n_M-1;
tmp_e1_ = [ ...
 init_euler_polar_a_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_azimu_b_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
 init_euler_gamma_z_MS_Mit___(1+nM1,1+niteration,1+ninit_Y_lm_cut_threshold) ...
	    ];
tmp_k_c_d = k_c_d(tmp_e0_,tmp_e1_);
tmp_euler_k_c_d_0lsq__(1+nM0,1+nM1) = tmp_k_c_d;
tmp_euler_k_c_d_0lsq__(1+nM1,1+nM0) = tmp_k_c_d;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
tmp_h2_0lsq__ = tmp_h2_0lsq__ + hist2d_0(init_euler_k_c_d_true__(1+init_index_offdiag_),tmp_euler_k_c_d_0lsq__(1+init_index_offdiag_),n_h,n_h,dlim_,dlim_)/n_iteration;
end;%for niteration=0:n_iteration-1;
tmp_h2_0lsq_norm__ = tmp_h2_0lsq__./repmat(sum(tmp_h2_0lsq__,1),[n_h,1]);
hold on;
stairs(linspace(0,2,n_h),tmp_h2_0lsq_norm__);
%stairs(linspace(0,2,n_h),tmp_h2_0lsq_norm__ - init_h2_rand_norm__);
%imagesc(tmp_h2_0lsq_norm__ - init_h2_rand_norm__,0.1*[-1,+1]);
%imagesc(tmp_h2_0lsq__);
%plot(dlim_,dlim_,'k-');
%plot([0,n_h+1],[0,n_h+1],'k-');
hold off;
xlim(dlim_);
ylim(0.2*[0,+1]);
%ylim(0.1*[-1,+1]);
%axis([0,2,0,2]); axisnotick;
%axis([0,n_h+1,0,n_h+1]); axisnotick;
xlabel('true distance'); ylabel('histogram (#)');
title(sprintf('avg MS cut %0.2f',init_Y_lm_cut_threshold));
drawnow;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig_B));
print('-djpeg',sprintf('%s.jpg',fname_fig_B));
print('-depsc',sprintf('%s.eps',fname_fig_B));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig_B),'file'));
if ( exist(sprintf('%s.jpg',fname_fig_B),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_B));
end;%if ( exist(sprintf('%s.jpg',fname_fig_B),'file'));
end;%if flag_plot;

flag_plot=1;
if flag_plot;
fname_fig_C = sprintf('%s_jpg/%s_FIGC',dir_trunk,infix_am_init);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig_C),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_C));
%%%%%%%%;
% plot viewing-angle-spread for those viewing-angles that are close. ;
% usign aligned function. ;
%%%%%%%%;
figure(5); clf; figbeach(); colormap(colormap_pm(64));
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
subplot(2,3,1+ninit_Y_lm_cut_threshold);
tmp_h2_0lsq__ = zeros(n_h,n_h);
tmp_euler_k_c_d_0lsq_ = zeros(n_M,n_M);
for nM0=0:n_M-1;
tmp_e0_ = [ ...
 alig_euler_polar_a_MS_Mt__(1+nM0,1+ninit_Y_lm_cut_threshold) ...
 alig_euler_azimu_b_MS_Mt__(1+nM0,1+ninit_Y_lm_cut_threshold) ...
 alig_euler_gamma_z_MS_Mt__(1+nM0,1+ninit_Y_lm_cut_threshold) ...
	    ];
for nM1=1+nM0:n_M-1;
tmp_e1_ = [ ...
 alig_euler_polar_a_MS_Mt__(1+nM1,1+ninit_Y_lm_cut_threshold) ...
 alig_euler_azimu_b_MS_Mt__(1+nM1,1+ninit_Y_lm_cut_threshold) ...
 alig_euler_gamma_z_MS_Mt__(1+nM1,1+ninit_Y_lm_cut_threshold) ...
	    ];
tmp_k_c_d = k_c_d(tmp_e0_,tmp_e1_);
tmp_euler_k_c_d_0lsq__(1+nM0,1+nM1) = tmp_k_c_d;
tmp_euler_k_c_d_0lsq__(1+nM1,1+nM0) = tmp_k_c_d;
end;%for nM1=0:n_M-1;
end;%for nM0=0:n_M-1;
tmp_h2_0lsq__ = tmp_h2_0lsq__ + hist2d_0(init_euler_k_c_d_true__(1+init_index_offdiag_),tmp_euler_k_c_d_0lsq__(1+init_index_offdiag_),n_h,n_h,dlim_,dlim_);
tmp_h2_0lsq_norm__ = tmp_h2_0lsq__./repmat(sum(tmp_h2_0lsq__,1),[n_h,1]);
hold on;
stairs(linspace(0,2,n_h),tmp_h2_0lsq_norm__);
%stairs(linspace(0,2,n_h),tmp_h2_0lsq_norm__ - init_h2_rand_norm__);
%imagesc(tmp_h2_0lsq_norm__ - init_h2_rand_norm__,0.1*[-1,+1]);
%imagesc(tmp_h2_0lsq__);
%plot(dlim_,dlim_,'k-');
%plot([0,n_h+1],[0,n_h+1],'k-');
hold off;
xlim(dlim_);
ylim(0.2*[0,+1]);
%ylim(0.1*[-1,+1]);
%axis([0,2,0,2]); axisnotick;
%axis([0,n_h+1,0,n_h+1]); axisnotick;
xlabel('true distance'); ylabel('histogram (#)');
title(sprintf('alig MS cut %0.2f',init_Y_lm_cut_threshold));
drawnow;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig_C));
print('-djpeg',sprintf('%s.jpg',fname_fig_C));
print('-depsc',sprintf('%s.eps',fname_fig_C));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig_C),'file'));
if ( exist(sprintf('%s.jpg',fname_fig_C),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_C));
end;%if ( exist(sprintf('%s.jpg',fname_fig_C),'file'));
end;%if flag_plot;

true_viewing_k_c_0_ = transpose(sin(euler_angle_marina_(1+0,:)).*cos(euler_angle_marina_(1+1,:)));
true_viewing_k_c_1_ = transpose(sin(euler_angle_marina_(1+0,:)).*sin(euler_angle_marina_(1+1,:)));
true_viewing_k_c_2_ = transpose(cos(euler_angle_marina_(1+0,:)));
true_viewing_k_c__ = [true_viewing_k_c_0_ , true_viewing_k_c_1_ , true_viewing_k_c_2_ ];
true_viewing_k_c__ = true_viewing_k_c__(1:n_M,:);
n_knn = 8;%n_knn = ceil(sqrt(n_M));
true_viewing_index_knn__ = knnsearch(true_viewing_k_c__,true_viewing_k_c__,'K',n_knn) - 1;

flag_plot=1;
if (flag_plot);
for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
fname_fig_D = sprintf('%s_jpg/%s_FIGD_ncut%d',dir_trunk,infix_am_init,ninit_Y_lm_cut_threshold);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig_D),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_D));
%%%%%%%%;
% Now step through each iteration and measure the angular location of the 32=sqrt(n_M) nearest-neighbors to each image. ;
%%%%%%%%;
tmp_markersize = 8;
c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_marina_(1,1:n_M),+euler_angle_marina_(2,1:n_M),0.35);
n_h = 64; lh_lim_ = [0,10];
c__ = colormap_beach(); n_c = size(c__,1); 
figure(3); np=0;
figbeach();
figbig;
subplot(3,6,1+np); np=np+1;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),tmp_markersize,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
axisnotick; axis image;
%%%%%%%%;
for niteration=0:n_iteration-1;
tmp_euler_azimu_b_ = init_euler_azimu_b_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold);
tmp_euler_polar_a_ = init_euler_polar_a_MS_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp_k_c__ = [tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ];
tmp_d_avg_ = zeros(n_M,1);
tmp_w_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_index_ = true_viewing_index_knn__(1+nM,:);
tmp_d_avg_(1+nM) = mean(sqrt(sum((repmat(tmp_k_c__(1+nM,:),[n_knn,1]) - tmp_k_c__(1+tmp_index_,:)).^2,2)));
tmp_w_avg_(1+nM) = mean( acos(1 - sum((repmat(tmp_k_c__(1+nM,:),[n_knn,1]) - tmp_k_c__(1+tmp_index_,:)).^2,2)/2) );
end;%for nM=0:n_M-1;
%tmp_nc_ = max(0,min(n_c-1,floor(n_c*real(tmp_d_avg_)/2)));
tmp_nc_ = max(0,min(n_c-1,floor(n_c*real(tmp_w_avg_)/pi)));
%%%%%%%%;
figure(3);
subplot(3,6,1+np); np = np+1; cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),tmp_markersize,c__(1+tmp_nc_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('ni%0.2d MS',niteration));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
tmp_euler_azimu_b_ = alig_euler_azimu_b_MS_Mt__(:,1+ninit_Y_lm_cut_threshold);
tmp_euler_polar_a_ = alig_euler_polar_a_MS_Mt__(:,1+ninit_Y_lm_cut_threshold);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp_k_c__ = [tmp_k_c_0_ , tmp_k_c_1_ , tmp_k_c_2_ ];
tmp_d_avg_ = zeros(n_M,1);
tmp_w_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_index_ = true_viewing_index_knn__(1+nM,:);
tmp_d_avg_(1+nM) = mean(sqrt(sum((repmat(tmp_k_c__(1+nM,:),[n_knn,1]) - tmp_k_c__(1+tmp_index_,:)).^2,2)));
tmp_w_avg_(1+nM) = mean( acos(1 - sum((repmat(tmp_k_c__(1+nM,:),[n_knn,1]) - tmp_k_c__(1+tmp_index_,:)).^2,2)/2) );
end;%for nM=0:n_M-1;
%tmp_nc_ = max(0,min(n_c-1,floor(n_c*real(tmp_d_avg_)/2)));
tmp_nc_ = max(0,min(n_c-1,floor(n_c*real(tmp_w_avg_)/pi)));
%%%%%%%%;
figure(3);
subplot(3,6,1+np); np = np+1; cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),tmp_markersize,c__(1+tmp_nc_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('alig MS'));
%%%%%%%%;
figbig;
disp(sprintf(' %% writing %s',fname_fig_D));
print('-djpeg',sprintf('%s.jpg',fname_fig_D));
print('-depsc',sprintf('%s.eps',fname_fig_D));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig_D),'file'));
if ( exist(sprintf('%s.jpg',fname_fig_D),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_D));
end;%if ( exist(sprintf('%s.jpg',fname_fig_D),'file'));
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
end;%if flag_plot;

clear tmp_*;

delete(fname_tmp_ABC);
end;%if (~flag_skip);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for am_init_nUX_rank=0:dat_n_UX_rank-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
