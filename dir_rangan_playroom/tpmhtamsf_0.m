function [...
 est_prctile__...
,est_vdist__...
,tmp_n_viewing_all...
,tmp_viewing_azimu_b_all_...
,tmp_viewing_polar_a_all_...
,true_viewing_polar_a_...
,true_viewing_azimu_b_...
,true_viewing_gamma_z_...
,true_viewing_delta_x_...
,true_viewing_delta_y_...
,est_M_S_index_...
] = ...
tpmhtamsf_0(...
 est_rseed...
,est_n_M...
,est_snr...
,est_n_UX_rank...
,dir_trunk...
,n_x_u...
,diameter_x_c...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,delta_r_max...
,svd_eps...
,n_delta_v_requested...
,n_w_max...
,CTF_uni_avg_k_p_...
,l_max_...
,c_k_Y_reco_...
,UX__...
,X_weight_r_...
);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% test_principled_marching_helper_trpv1_alternating_minimization_sample_frand_0. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

rng(0);

n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
verbose=1;

%%%%%%%%;
% First use the estimated volume to generate templates. ;
%%%%%%%%;
tmp_t = tic;
[S_k_p__,tmp_n_w_,tmp_weight_2d_k_p_r_,tmp_weight_2d_k_all_,tmp_n_viewing_all,tmp_viewing_azimu_b_all_,tmp_viewing_polar_a_all_,tmp_n_viewing_polar_a,tmp_viewing_polar_a_,tmp_n_viewing_azimu_b_,tmp_template_k_c_0__,tmp_template_k_c_1__,tmp_template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,l_max_,c_k_Y_reco_,viewing_k_eq_d,0,n_w_);
n_S = tmp_n_viewing_all;
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_reco_ --> S_k_p__ time %0.2fs',tmp_t));
%%%%%%%%;
% Use templates S_k_p___ to estimate typical signal strength (per unit area in k_p). ;
%%%%%%%%;
tmp_t = tic();
template_area_element_ = tmp_weight_2d_k_all_*(2*pi)^2; %<-- now sum(template_area_element_)==pi*k_p_r_max^2;
S_k_p_l1_avg = mean(transpose(template_area_element_)*abs(S_k_p__.*CTF_uni_avg_k_p_))/(pi*k_p_r_max^2);
if (est_snr<=0); est_sigma = 0; end;
if (est_snr> 0); est_sigma = S_k_p_l1_avg/est_snr; end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% est_sigma %0.6f: %0.3fs',est_sigma,tmp_t)); end;
%%%%%%%%
% Set up translation distribution. ;
%%%%%%%%;
% Note that a 2d isotropic gaussian with std delta_sigma has the property ;
% that a fraction (1-exp(-R^2/(2*delta_sigma^2))) is contained within radius R. ;
% For example, try out: ;
% tmp_delta_sigma = 1.45; tmp_R = 3.3; tmp_A_ = randn(1024*64,2)*tmp_delta_sigma; tmp_r_ = sqrt(sum(tmp_A_.^2,2)); tmp_0in = numel(find(tmp_r_<tmp_R))/numel(tmp_r_); tmp_1in = 1-exp(-tmp_R^2/(2*tmp_delta_sigma^2)); disp(sprintf(' %% I_monte: %0.6f I_form: %0.6f',tmp_0in,tmp_1in));
% Thus, to ensure that 95-percent of the true-displacements lie within delta_r_max, ;
% we require that exp(-delta_r_max^2/(2*delta_sigma^2)) = 0.05 ;
% or that delta_sigma = sqrt(delta_r_max^2/log(20^2)). ;
% Roughly speaking, delta_sigma should approximately equal delta_r_max/2.5. ;
% Conversely: if we know delta_sigma, we can conclude that delta_r_max = delta_sigma*sqrt(log(20^2)). ; 
%%%%%%%%;
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
true_viewing_delta_x_ = delta_r_s*randn(est_n_M,1);
true_viewing_delta_y_ = delta_r_s*randn(est_n_M,1);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f',delta_r_p,delta_r_max,delta_r_s,delta_r_N));
%%%%%%%%;
% Generate synthetic images from templates S_k_p__. ;
%%%%%%%%;
tmp_t = tic();
true_viewing_polar_a_ = zeros(est_n_M,1);
true_viewing_azimu_b_ = zeros(est_n_M,1);
true_viewing_gamma_z_ = 2*pi*rand(est_n_M,1);
est_M_k_p__ = zeros(n_w_sum,est_n_M);
est_UX_M_k_p___ = zeros(n_w_max,est_n_M,est_n_UX_rank);
est_UX_M_k_q___ = zeros(n_w_max,est_n_M,est_n_UX_rank);
est_S_permutation_ = randperm(n_S)-1;
est_M_S_index_ = zeros(est_n_M,1);
nS=0;
for nM=0:est_n_M-1;
if (mod(nM,100)==0); disp(sprintf(' %% nM %d/%d',nM,est_n_M)); end;
est_M_S_index = est_S_permutation_(1+nS);
est_M_S_index_(1+nM) = est_M_S_index;
tmp_M_k_p_ = S_k_p__(:,1+est_M_S_index).*CTF_uni_avg_k_p_; %<-- extract random template. ;
true_viewing_polar_a_(1+nM) = tmp_viewing_polar_a_all_(1+est_M_S_index);
true_viewing_azimu_b_(1+nM) = tmp_viewing_azimu_b_all_(1+est_M_S_index);
tmp_M_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_,+true_viewing_gamma_z_(1+nM));
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-true_viewing_delta_x_(1+nM),-true_viewing_delta_y_(1+nM));
tmp_M_k_p_ = tmp_M_k_p_ + est_sigma*randn(n_w_sum,1)./sqrt( template_area_element_ * n_w_sum / (pi*k_p_r_max^2) );
est_M_k_p__(:,1+nM) = tmp_M_k_p_;
tmp_M_k_p__ = reshape(tmp_M_k_p_,[n_w_max,n_k_p_r]);
for nUX_rank=0:est_n_UX_rank-1;
tmp_UX_M_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
est_UX_M_k_p___(:,1+nM,1+nUX_rank) = tmp_UX_M_k_p_;
est_UX_M_k_q___(:,1+nM,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_M_k_p_);
end;%for nUX_rank=0:est_n_UX_rank-1;
nS=nS+1; if (nS>=n_S); nS=0; end;
end;%for nM=0:est_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% est_M_k_p__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% set est_M_k_q__. ;
%%%%%%%%;
tmp_t = tic();
est_M_k_q__ = zeros(n_w_sum,est_n_M);
for nM=0:est_n_M-1;
est_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,est_M_k_p__(:,1+nM));
end;%for nM=0:est_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% est_M_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% prepare precomputation for ampm. ;
%%%%%%%%;
pm_n_UX_rank = est_n_UX_rank;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,est_n_M,est_M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,est_n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
pm_n_UX_rank = est_n_UX_rank;
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
%%%%%%%%;
% Now generate principal-templates. ;
%%%%%%%%;
tmp_t = tic();
UX_S_CTF_k_p_wSn___ = zeros(n_w_max,n_S,pm_n_UX_rank);
UX_S_CTF_k_q_wSn___ = zeros(n_w_max,n_S,pm_n_UX_rank);
for nS=0:n_S-1;
if (mod(nS,100)==0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
tmp_S_CTF_k_p_ = S_k_p__(:,1+nS).*CTF_uni_avg_k_p_; %<-- use average templates here under assumption that templates are used alone. ;
tmp_S_CTF_k_p__ = reshape(tmp_S_CTF_k_p_,[n_w_max,n_k_p_r]);
for nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_S_CTF_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_UX_S_CTF_k_p_ = tmp_UX_S_CTF_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_S_CTF_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
UX_S_CTF_k_p_wSn___(:,1+nS,1+nUX_rank) = tmp_UX_S_CTF_k_p_;
UX_S_CTF_k_q_wSn___(:,1+nS,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_S_CTF_k_p_);
end;%for nUX_rank=0:pm_n_UX_rank-1;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
%%%%%%%%;
% Now use precomputation to calculate X_wSM___. ; 
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
% Batches images into batches of size n_M_Mbatch. ;
% Batches templates into batches of size n_S_Sbatch. ;
% Only stores the optimal translation for each image. ;
%%%%%%%%;
est_prctile__ = zeros(est_n_M,pm_n_UX_rank);
est_vdist__ = zeros(est_n_M,pm_n_UX_rank);
for tmp_pm_n_UX_rank=1:pm_n_UX_rank;
disp(sprintf(' %% tmp_pm_n_UX_rank %d/%d',tmp_pm_n_UX_rank,pm_n_UX_rank));
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
tmp_t = tic();
UX_S_CTF_k_p__ = reshape(permute(UX_S_CTF_k_p_wSn___(:,:,1:tmp_pm_n_UX_rank),[1,3,2]),[n_w_max*tmp_pm_n_UX_rank,n_S]);
UX_S_CTF_k_q__ = reshape(permute(UX_S_CTF_k_q_wSn___(:,:,1:tmp_pm_n_UX_rank),[1,3,2]),[n_w_max*tmp_pm_n_UX_rank,n_S]);
UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
UX_S_l2_(1+nS) = innerproduct_p_quad(tmp_pm_n_k_p_r,tmp_pm_k_p_r_,tmp_pm_weight_2d_k_p_r_/(2*pi),tmp_pm_n_w_,tmp_pm_n_w_sum,UX_S_CTF_k_p__(:,1+nS),UX_S_CTF_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
SS_k_q_ = svd(UX_S_CTF_k_q__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(UX_S_CTF_k_q__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(UX_S_CTF_k_q__,n_S_rank);
if (verbose); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
tmp_t = tic();
tmp_UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,est_n_M,tmp_pm_n_UX_rank,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% tmp_UX_M_l2_dM__: %0.3fs',tmp_t)); end;
tmp_t = tic();
n_M_Mbatch = 24;
n_S_Sbatch = 24;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___1(...
 FTK...
,n_w_...
,tmp_pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_Sbatch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,est_n_M...
,n_M_Mbatch...
,svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,tmp_UX_M_l2_dM__...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
est_prctile_ = zeros(est_n_M,1);
for nM=0:est_n_M-1;
tmp_X_S_ = max(real(X_wSM___(:,:,1+nM)),[],1);
tmp_X_S = tmp_X_S_(1+est_M_S_index_(1+nM));
est_prctile_(1+nM) = sum(tmp_X_S_>=tmp_X_S)/n_S;
end;%for nM=0:est_n_M-1;
est_prctile__(:,tmp_pm_n_UX_rank) = est_prctile_;
%%%%%%%%;
est_vdist_ = zeros(est_n_M,1);
for nM=0:est_n_M-1;
tmp_X_S_ = max(real(X_wSM___(:,:,1+nM)),[],1);
[~,tmp_index] = max(tmp_X_S_); tmp_index = tmp_index-1;
est_viewing_azimu_b = tmp_viewing_azimu_b_all_(1+tmp_index);
est_viewing_polar_a = tmp_viewing_polar_a_all_(1+tmp_index);
tru_viewing_azimu_b = true_viewing_azimu_b_(1+nM);
tru_viewing_polar_a = true_viewing_polar_a_(1+nM);
est_vdist_(1+nM) = sphere_distance_0(est_viewing_polar_a,est_viewing_azimu_b,tru_viewing_polar_a,tru_viewing_azimu_b);
end;%for nM=0:est_n_M-1;
est_vdist__(:,tmp_pm_n_UX_rank) = est_vdist_;
%%%%%%%%;
end;%for tmp_pm_n_UX_rank=1:pm_n_UX_rank;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;





