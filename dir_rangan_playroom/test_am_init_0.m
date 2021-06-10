verbose=1;

n_S = n_viewing_all;
n_M = n_image_sub;

%%%%%%%%;
% First create 'noiseless' images. ;
% Saved as T_k_p__. ;
%%%%%%%%;
%{
X_TM_ = zeros(n_M,1);
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
T_k_p__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
if (mod(nM,100)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
tmp_euler_polar_a = +euler_angle_marina_(1,1+nM);
tmp_euler_azimu_b = +euler_angle_marina_(2,1+nM);
tmp_euler_gamma_z = -euler_angle_marina_(3,1+nM);
tmp_image_delta_x = +1.0*delta_read_x_(1+nM);
tmp_image_delta_y = +1.0*delta_read_y_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_idx_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_);
X_TM_(1+nM) = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
T_k_p__(:,1+nM) = T_k_p_;
end;%for nM=0:n_M-1;
%%%%%%%%;
% Plot example of one image. ;
%%%%%%%%;
flag_plot=0;
nM = 688; %<-- low correlation for rib80s. ;
M_k_p_ = M_k_p__(:,1+nM);
T_k_p_ = T_k_p__(:,1+nM);
if flag_plot;
figure(1); clf;
figbeach();
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(M(k))');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(M(k))');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(T(k))');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(T(k))');
figbig;
drawnow;
end;%if flag_plot;
 %}

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
ndelta_r_max_factor=0;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
ndat_rseed=0;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
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

tmp_t = tic();
dat_M_k_q__ = zeros(n_w_uni_sum,dat_n_M);
for nM=0:dat_n_M-1;
dat_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_uni_,n_w_uni_sum,dat_M_k_p__(:,1+nM));
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% dat_M_k_q__: %0.3fs',tmp_t)); end;

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

nUX_rank = 2; %<-- increase up to n_UX_rank. ;
init_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
init_n_order = dat_n_order;
init_pm_n_UX_rank = 1+nUX_rank;
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
f_rand = 0;
init_a_CTF_UX_Y_true_ = reshape(init_a_CTF_UX_Y_quad__(:,1:init_pm_n_UX_rank),[n_lm_max*init_pm_n_UX_rank,1]); %<-- use as ground truth for this set of principled-image-rings. ;

[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_idx_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;

init_Y_lm_cut_threshold_ = [4];%init_Y_lm_cut_threshold_ = [2,4,6,8,12,24];%init_Y_lm_cut_threshold_ = 1:1:1+2*11; %init_Y_lm_cut_threshold_ = 1:init_pm_l_max_max; 
n_init_Y_lm_cut_threshold = numel(init_Y_lm_cut_threshold_);
n_iteration = 16;
init_a_CTF_UX_Y_0lsq_Yit___ = zeros(init_pm_n_lm_sum,n_iteration,n_init_Y_lm_cut_threshold);
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

for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
init_Y_lm_cut_threshold = init_Y_lm_cut_threshold_(1+ninit_Y_lm_cut_threshold);
for niteration=0:n_iteration-1;

if (verbose>0); disp(sprintf(' %% ninit_Y_lm_cut_threshold %d/%d; niteration %d/%d',ninit_Y_lm_cut_threshold,n_init_Y_lm_cut_threshold,niteration,n_iteration)); end;

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
%init_Y_lm_cut_ = init_Y_l_val_+abs(init_Y_m_val_)<=init_Y_lm_cut_threshold;
init_Y_lm_cut_ = init_Y_l_val_+0*abs(init_Y_m_val_)<=init_Y_lm_cut_threshold;

flag_true_vs_rand=0;
disp(sprintf(' %% flag_true_vs_rand %d',flag_true_vs_rand));

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
%%%%%%%%;
euler_polar_a_ = 1*pi*rand(n_M,1);
euler_azimu_b_ = 2*pi*rand(n_M,1);
euler_gamma_z_ = 2*pi*rand(n_M,1);
disp(sprintf(' %% pre-centering the images. '));
image_delta_x_ = +1.0*delta_read_x_(1:n_M); %image_delta_x_ = zeros(n_M,1);
image_delta_y_ = +1.0*delta_read_y_(1:n_M); %image_delta_y_ = zeros(n_M,1);
image_I_value_ = ones(n_M,1);
end;%if flag_true_vs_rand==0;
%%%%%%%%;
% Construct M_k_q__ while taking into account the translations. ;
%%%%%%%%;
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
,+image_delta_x_(1+nM) ...
,+image_delta_y_(1+nM) ...
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
% Now calculate norms of the translated images. ;
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
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
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
,[] ...
,[] ...
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
init_a_CTF_UX_Y_0lsq_Yit___(:,1+niteration,1+ninit_Y_lm_cut_threshold) = init_a_CTF_UX_Y_0lsq_;
%%%%%%%%;
% Compare current model to init_a_CTF_UX_Y_true_. ;
%%%%%%%%;
tmp_t = tic();
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,init_a_CTF_UX_Y_0lsq_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(init_pm_n_k_p_r,init_pm_k_p_r_,init_pm_k_p_r_max,init_pm_weight_3d_k_p_r_,0,init_pm_l_max_,init_a_CTF_UX_Y_true_,flipY(init_pm_n_k_p_r,init_pm_l_max_,init_a_CTF_UX_Y_0lsq_));
tmp_X_best = max(real(tmp_X_best_orig),real(tmp_X_best_flip));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% register_spharm_to_spharm_wigner_0: %0.3fs',tmp_t)); end;
if (verbose>-1); disp(sprintf(' %% init_a_CTF_UX_Y_true_ vs a_CTF_UX_Y_lsq0_: correlation %+0.6f',tmp_X_best)); end;

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

if flag_plot;
%%%%%%%%;
% Now plot displacements. ;
%%%%%%%%;
c2d__ = colormap_gaussian_2d(delta_read_x_,delta_read_y_,delta_sigma,0.35);
n_h = 64; lh_lim_ = [0,10];
figure(2);
figbeach();
figbig;
subplot(2,3,1);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_M),delta_read_y_(1:n_M),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
subplot(2,3,4); cla;
imagesc(log2(1+hist2d_0(delta_read_x_(1:n_M),delta_read_y_(1:n_M),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
figure(2);
subplot(2,3,2); cla;
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(image_delta_x_SM_(:),image_delta_y_SM_(:),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square;
%xlabel('dx'); ylabel('dy');
title(sprintf('SM'));
subplot(2,3,5); cla;
imagesc(log2(1+hist2d_0(image_delta_x_SM_(1:n_M),image_delta_y_SM_(1:n_M),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
subplot(2,3,3); cla;
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_M) - image_delta_x_SM_(:),delta_read_y_(1:n_M) - image_delta_y_SM_(:),16,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square;
%xlabel('dx'); ylabel('dy');
title(sprintf('true - SM'));
subplot(2,3,6); cla;
imagesc(log2(1+hist2d_0(delta_read_x_(1:n_M) - image_delta_x_SM_(1:n_M),delta_read_y_(1:n_M) - image_delta_y_SM_(1:n_M),n_h,n_h,[-1,+1],[-1,+1])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
end;%if flag_plot;

if (flag_plot);
%%%%%%%%;
% Now plot viewing angles. ;
%%%%%%%%;
c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_marina_(1,1:n_M),+euler_angle_marina_(2,1:n_M),0.35);
n_h = 64; lh_lim_ = [0,10];
figure(3);
figbeach();
figbig;
subplot(2,3,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
subplot(2,3,4); cla;
imagesc(log2(1+hist2d_0(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
%%%%%%%%;
N_wavelength = 0.0;
[ ...
 X_best_orig ...
,polar_a_best_orig ...
,azimu_b_best_orig ...
,gamma_z_best_orig ...
,delta_best_orig_ ...
,a_k_Y_best_orig_ ...
,b_k_Y_best_orig_ ...
] = ...
register_spharm_to_spharm_wigner_0( ...
 init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_k_p_r_max ...
,init_pm_weight_3d_k_p_r_ ...
,N_wavelength ...
,init_pm_l_max_ ...
,init_a_CTF_UX_Y_true_ ...
,init_a_CTF_UX_Y_0lsq_ ...
);
[ ...
 X_best_flip ...
,polar_a_best_flip ...
,azimu_b_best_flip ...
,gamma_z_best_flip ...
,delta_best_flip_ ...
,a_k_Y_best_flip_ ...
,b_k_Y_best_flip_ ...
] = ...
register_spharm_to_spharm_wigner_0( ...
 init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_k_p_r_max ...
,init_pm_weight_3d_k_p_r_ ...
,N_wavelength ...
,init_pm_l_max_ ...
,init_a_CTF_UX_Y_true_ ...
,flipY(init_pm_n_k_p_r,init_pm_l_max_,init_a_CTF_UX_Y_0lsq_) ...
);
%%%%%%%%;
X_best = X_best_orig;
polar_a_best = polar_a_best_orig;
azimu_b_best = azimu_b_best_orig;
gamma_z_best = gamma_z_best_orig;
delta_best_ = delta_best_orig_;
a_k_Y_best_ = a_k_Y_best_orig_;
b_k_Y_best_ = b_k_Y_best_orig_;
if (X_best_flip>X_best_orig);
X_best = X_best_flip;
polar_a_best = polar_a_best_flip;
azimu_b_best = azimu_b_best_flip;
gamma_z_best = gamma_z_best_flip;
delta_best_ = delta_best_flip_;
a_k_Y_best_ = a_k_Y_best_flip_;
b_k_Y_best_ = b_k_Y_best_flip_;
end;%if (X_best_flip>X_best_orig);
euler_best_ = [+azimu_b_best,+polar_a_best,+gamma_z_best];
R_best_ = euler_to_R(euler_best_);
%%%%%%%%;
tmp_euler_azimu_b_ = euler_azimu_b_MS_(:);
tmp_euler_polar_a_ = euler_polar_a_MS_(:);
tmp_k_c_0_ = sin(tmp_euler_polar_a_).*cos(tmp_euler_azimu_b_);
tmp_k_c_1_ = sin(tmp_euler_polar_a_).*sin(tmp_euler_azimu_b_);
tmp_k_c_2_ = cos(tmp_euler_polar_a_);
tmp__ = inv(R_best_) * transpose([tmp_k_c_0_(:) , tmp_k_c_1_(:) , tmp_k_c_2_(:)]);
tmp_k_c_0_(:) = transpose(tmp__(1+0,:));
tmp_k_c_1_(:) = transpose(tmp__(1+1,:));
tmp_k_c_2_(:) = transpose(tmp__(1+2,:));
tmp_k_c_01_ = sqrt(tmp_k_c_0_.^2 + tmp_k_c_1_.^2);
tmp_euler_azimu_b_ = periodize(atan2(tmp_k_c_1_,tmp_k_c_0_),0,2*pi);
tmp_euler_polar_a_ = atan2(tmp_k_c_01_,tmp_k_c_2_);
tmp_euler_azimu_b_dif_ = pi/1 + periodize((0*pi+transpose(euler_angle_marina_(2,1:n_M))) - (0*pi+tmp_euler_azimu_b_),-pi,+pi);
tmp_euler_polar_a_dif_ = pi/2 + (1*pi-transpose(euler_angle_marina_(1,1:n_M))) - (1*pi-tmp_euler_polar_a_);
%%%%%%%%;
figure(3);
subplot(2,3,2); cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,64,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('MS'));
subplot(2,3,5); cla;
imagesc(log2(1+hist2d_0(0*pi+tmp_euler_azimu_b_,1*pi-tmp_euler_polar_a_,2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
subplot(2,3,3); cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter( ...
 tmp_euler_azimu_b_dif_ ...
,tmp_euler_polar_a_dif_ ...
,64 ...
,c2d__(1:n_M,:) ...
,'filled' ...
);
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('true - MS'));
subplot(2,3,6); cla;
imagesc(log2(1+hist2d_0( ...
 tmp_euler_azimu_b_dif_ ...
,tmp_euler_polar_a_dif_ ...
,2*n_h,n_h,[0,2*pi],[0,1*pi])),lh_lim_);
set(gca,'Ydir','normal');
axisnotick; axis image;
drawnow();
clear tmp_*;
end;%if flag_plot;

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

end;%for niteration=0:n_iteration-1;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;

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
end;%if flag_plot;

flag_plot=1;
if flag_plot;
%%%%%%%%;
% plot viewing-angle-spread for those viewing-angles that are close. ;
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
title(sprintf('MS cut %0.2f',init_Y_lm_cut_threshold));
drawnow;
end;%for ninit_Y_lm_cut_threshold=0:n_init_Y_lm_cut_threshold-1;
figbig;
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

%%%%%%%%;
% Now step through each iteration and measure the angular location of the 32=sqrt(n_M) nearest-neighbors to each image. ;
%%%%%%%%;
ninit_Y_lm_cut_threshold=0; %<-- single threshold. ;
tmp_markersize = 8;
c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_marina_(1,1:n_M),+euler_angle_marina_(2,1:n_M),0.35);
n_h = 64; lh_lim_ = [0,10];
c__ = colormap_beach(); n_c = size(c__,1); 
figure(3);
figbeach();
figbig;
subplot(3,6,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),tmp_markersize,c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
axisnotick; axis image;
%%%%%%%%;
for niteration=0:n_iteration-1;
tmp_euler_azimu_b_ = init_euler_azimu_b_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold);
tmp_euler_polar_a_ = init_euler_polar_a_SM_Mit___(:,1+niteration,1+ninit_Y_lm_cut_threshold);
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
subplot(3,6,1+niteration+1); cla;
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),tmp_markersize,c__(1+tmp_nc_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
title(sprintf('SM'));
end;%for niteration=0:n_iteration-1;

clear tmp_*;
end;%if flag_plot;

disp('returning'); return;

%%%%%%%%;
% plot the spherical-harmonic expansion. ;
%%%%%%%%;
[ ...
 tmp_n_all ...
,tmp_n_sub_ ...
,tmp_k_p_r_all_ ...
,tmp_azimu_b_all_ ...
,tmp_polar_a_all_ ...
,tmp_weight_all_ ...
,tmp_a_true_all_ ...
,tmp_k_c_0_all_ ...
,tmp_k_c_1_all_ ...
,tmp_k_c_2_all_ ...
] = ...
convert_spharm_to_k_p_0( ...
 verbose ...
,init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_l_max_ ...
,init_a_CTF_UX_Y_true_ ...
,init_viewing_k_eq_d/8 ...
);
alim_ = 0*mean(abs(tmp_a_true_all_)) + 2.5*std(abs(tmp_a_true_all_),1)*[-1,+1];
subplot(2,3,1);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,real(tmp_a_true_all_),alim_);
xlim([0,2*pi]);ylim([0,1*pi]); axisnotick; title('real true');
subplot(2,3,2);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,imag(tmp_a_true_all_),alim_);
xlim([0,2*pi]);ylim([0,1*pi]); axisnotick; title('imag true');
subplot(2,3,3);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,abs(tmp_a_true_all_),[0,2*max(alim_)],[],0);
xlim([-1,+1]);ylim([-1,+1]);zlim([-1,+1]); axisnotick; title('abs true');
[ ...
 tmp_n_all ...
,tmp_n_sub_ ...
,tmp_k_p_r_all_ ...
,tmp_azimu_b_all_ ...
,tmp_polar_a_all_ ...
,tmp_weight_all_ ...
,tmp_a_0lsq_all_ ...
,tmp_k_c_0_all_ ...
,tmp_k_c_1_all_ ...
,tmp_k_c_2_all_ ...
] = ...
convert_spharm_to_k_p_0( ...
 verbose ...
,init_pm_n_k_p_r ...
,init_pm_k_p_r_ ...
,init_pm_l_max_ ...
,0*init_a_CTF_UX_Y_0lsq_ + 1*crandn(size(init_a_CTF_UX_Y_0lsq_)).*init_Y_lm_cut_ ...
,init_viewing_k_eq_d/8 ...
);
alim_ = 0*mean(abs(tmp_a_0lsq_all_)) + 2.5*std(abs(tmp_a_0lsq_all_),1)*[-1,+1];
subplot(2,3,4);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,real(tmp_a_0lsq_all_),alim_);
xlim([0,2*pi]);ylim([0,1*pi]); axisnotick; title('real 0lsq');
subplot(2,3,5);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,imag(tmp_a_0lsq_all_),alim_);
xlim([0,2*pi]);ylim([0,1*pi]); axisnotick; title('imag 0lsq');
subplot(2,3,6);
imagesc_polar_a_azimu_b_0(tmp_polar_a_all_,tmp_azimu_b_all_,abs(tmp_a_0lsq_all_),[0,2*max(alim_)],[],0);
xlim([-1,+1]);ylim([-1,+1]);zlim([-1,+1]); axisnotick; title('abs 0lsq');
figbig;
clear tmp_*;

%%%%%%%%;
% Now repeat, except this time limiting the order+degree of the random function. ;
%%%%%%%%;
