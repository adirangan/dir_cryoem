%clear;
%test_pm_trpv1_1;

delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [1,2,4,8,16];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
ndat_rseed = 0;%for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
ndelta_r_max_factor = 2;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
ndat_n_UX_rank=n_dat_n_UX_rank-1;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
%%%%%%%%;
ut_fname_pre = sprintf('%s_mat/ut%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
ut_fname_mat = sprintf('%s.mat',ut_fname_pre);
tmp_ut_ = load(ut_fname_mat);

euler_polar_a_true_ = +euler_polar_a_star_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_star_(1+(0:n_M-1));
euler_gamma_z_true_ = -euler_gamma_z_star_(1+(0:n_M-1)); %<-- note sign change. ;
image_delta_x_true_ = +image_delta_x_star_plus_M_abs_x_c_0_avg_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_star_plus_M_abs_x_c_1_avg_(1+(0:n_M-1));

parameter = tmp_ut_.parameter;
delta_r_max = parameter.delta_r_max;
pm_n_UX_rank = dat_n_UX_rank;
FTK = [];
VSCTF_Mc__ = [];
flag_image_delta_upd_ = [];
image_I_value_ = [];
M_k_q__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];
%%%%%%%%;
a_CTF_avg_UX_Y_true__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_avg_UX_Y_true__(1:tmp_n_lm,1+nUX_rank) = a_CTF_avg_UX_Y_true__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UX_Y_true_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
a_CTF_avg_UX_Y_true_ = reshape(a_CTF_avg_UX_Y_true__(:,1:pm_n_UX_rank),[n_lm_max*pm_n_UX_rank,1]);
%%%%%%%%;

verbose=2;
if (verbose); disp(sprintf(' %% [entering test_ampmleut_3]: alternating-minimization with principal-modes, local-exclusion and updating-translations')); end;

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

%%%%;
a_CTF_avg_UX_Y_ = tmp_ut_.a_CTF_avg_UX_Y__(:,n_iteration);
euler_polar_a_ = tmp_ut_.euler_polar_a__(:,n_iteration);
euler_azimu_b_ = tmp_ut_.euler_azimu_b__(:,n_iteration);
euler_gamma_z_ = tmp_ut_.euler_gamma_z__(:,n_iteration);
image_delta_x_acc_ = tmp_ut_.image_delta_x_acc__(:,n_iteration);
image_delta_y_acc_ = tmp_ut_.image_delta_y_acc__(:,n_iteration);
image_delta_x_upd_ = tmp_ut_.image_delta_x_upd__(:,n_iteration);
image_delta_y_upd_ = tmp_ut_.image_delta_y_upd__(:,n_iteration);
image_delta_x_ = image_delta_x_acc_ + image_delta_x_upd_ ;
image_delta_y_ = image_delta_y_acc_ + image_delta_y_upd_ ;
%%%%%%%%;

%%%%%%%%;
flag_plot=0;
if flag_plot;
ampmut_align_to_reference_0( ...
 parameter ...
,l_max_max ...
,pm_n_UX_rank ...
,a_CTF_avg_UX_Y_true_ ...
,n_M ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
,1 ...
,a_CTF_avg_UX_Y_ ...
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_delta_x_ ...
,image_delta_y_ ...
);
end;%if flag_plot;
%%%%%%%%;

%{
%%%%%%%%;
% Typical definitions. ;
%%%%%%%%;
l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
%}

n_w_max = n_w_max + mod(n_w_max,2); %<-- round up to nearest even number. ;
l_max_max = n_w_max/2 - 1;
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

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

if (order_limit_MS>-1);
Y_lm_cut_ = ones(pm_n_lm_sum,1);
Y_l_val_ = zeros(pm_n_lm_sum,1);
Y_m_val_ = zeros(pm_n_lm_sum,1);
na=0;
for pm_nk_p_r=0:pm_n_k_p_r-1;
l_max = pm_l_max_(1+pm_nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Y_l_val_(1+na) = l_val;
Y_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%Y_lm_cut_ = Y_l_val_+abs(Y_m_val_)<=Y_lm_cut_threshold;
Y_lm_cut_ = Y_l_val_+0*abs(Y_m_val_)<=order_limit_MS;
end;%if (order_limit_MS>-1);

if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

if isempty(VSCTF_Mc__);
if (n_CTF_rank<=0 | isempty(CTF_index_));
n_CTF_rank = 1;
USCTF_kc__ = ones(n_k_p_r,1); SCTF_C__ = 1; VCTF_Mc__ = ones(n_M,1);
else;
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_index_(1:n_M)),n_CTF_rank);
end;%if (n_CTF_rank<=0 | isempty(CTF_index_));
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
end;%if isempty(VSCTF_Mc__);

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
a_CTF_avg_UX_Y__ = zeros(pm_n_lm_sum,n_iteration);
%%%%%%%%;
% initialize current euler-angles randomly. ;
%%%%%%%%;
if isempty(euler_polar_a_); euler_polar_a_ = 1*pi*rand(n_M,1); end;
if isempty(euler_azimu_b_); euler_azimu_b_ = 2*pi*rand(n_M,1); end;
if isempty(euler_gamma_z_); euler_gamma_z_ = 2*pi*rand(n_M,1); end;
if isempty(image_delta_x_acc_); image_delta_x_acc_ = zeros(n_M,1); end; %<-- accumulated displacement (i.e., current image center). ;
if isempty(image_delta_y_acc_); image_delta_y_acc_ = zeros(n_M,1); end; %<-- accumulated displacement (i.e., current image center). ;
if isempty(image_delta_x_upd_); image_delta_x_upd_ = zeros(n_M,1); end; %<-- update to displacement (i.e., current image shift). ;
if isempty(image_delta_y_upd_); image_delta_y_upd_ = zeros(n_M,1); end; %<-- update to displacement (i.e., current image shift). ;
image_delta_x_bit_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
image_delta_y_bit_ = zeros(n_M,1); %<-- increment to displacement update (calculated each iteration). ;
if isempty(flag_image_delta_upd_); flag_image_delta_upd_ = ones(n_M,1); end; %<-- flag identifying principal-images that need to be recalculated. ;
if isempty(image_I_value_); image_I_value_ = ones(n_M,1); end;
if isempty(M_k_q__); M_k_q__ = zeros(n_w_sum,n_M); end;
if isempty(svd_VUXM_lwnM____); svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank,n_M); end;
if isempty(UX_M_l2_dM__); UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M); end;

ampmleut_template_tree = [];
if ( isempty(ampmleut_template_tree));
if (~isfield(parameter,'template_viewing_k_eq_d_min')); parameter.template_viewing_k_eq_d_min = template_viewing_k_eq_d; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'ampmleut_template_tree_n_level')); parameter.ampmleut_template_tree_n_level = 3; end; %<-- parameter_bookmark. ;
ampmleut_template_tree = get_template_tree_0(parameter.template_viewing_k_eq_d_min,parameter.ampmleut_template_tree_n_level);
ampmleut_n_level = parameter.ampmleut_template_tree_n_level;
end;%if ( isempty(ampmleut_template_tree));

%%%%%%%%;
% Visualize second-to-finest-level. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(1);clf;hold on; figsml;
tmp_nlevel = ampmleut_n_level-2;
tmp_n_S = ampmleut_template_tree.n_S_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_use = 25;
hold on;
plot3( ...
 ampmleut_template_tree.viewing_pole_k_c_0_all__{ampmleut_n_level}(:) ...
,ampmleut_template_tree.viewing_pole_k_c_1_all__{ampmleut_n_level}(:) ...
,ampmleut_template_tree.viewing_pole_k_c_2_all__{ampmleut_n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
for tmp_nS=0:tmp_n_S-1;
str_symbol = str_symbol_{1+mod(tmp_nS,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nS/tmp_n_S)));
tmp_index_ = ampmleut_template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS};
plot3( ...
 ampmleut_template_tree.viewing_pole_k_c_0_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_1_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_2_all__{ampmleut_n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_use ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
end;%for tmp_nS=0:tmp_n_S-1;
hold off;
axis equal; axis vis3d;
end;%if flag_plot;
%%%%%%%%;
% Visualize expanded neighborhoods at second-to-finest-level. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(2);clf;hold on; figsml;
tmp_nlevel = ampmleut_n_level-2;
tmp_n_S = ampmleut_template_tree.n_S_(1+tmp_nlevel);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_big = 25; markersize_sml = 15;
tmp_nS0 = floor(tmp_n_S/2);
str_symbol = str_symbol_{1+mod(tmp_nS0,5)};
nc = max(0,min(n_c_beach-1,floor(n_c_beach*tmp_nS0/tmp_n_S)));
hold on;
plot3( ...
 ampmleut_template_tree.viewing_pole_k_c_0_all__{ampmleut_n_level}(:) ...
,ampmleut_template_tree.viewing_pole_k_c_1_all__{ampmleut_n_level}(:) ...
,ampmleut_template_tree.viewing_pole_k_c_2_all__{ampmleut_n_level}(:) ...
,'.' ...
,'Color' , 'k' ...
);
tmp_index_ = setdiff(ampmleut_template_tree.index_neighbor_fine_from_coarse_exp___{1+tmp_nlevel}{1+tmp_nS0},ampmleut_template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS0});
plot3( ...
 ampmleut_template_tree.viewing_pole_k_c_0_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_1_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_2_all__{ampmleut_n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_sml ...
,'MarkerFaceColor', 0.5*c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
tmp_index_ = ampmleut_template_tree.index_neighbor_fine_from_coarse___{1+tmp_nlevel}{1+tmp_nS0};
plot3( ...
 ampmleut_template_tree.viewing_pole_k_c_0_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_1_all__{ampmleut_n_level}(1+tmp_index_) ...
,ampmleut_template_tree.viewing_pole_k_c_2_all__{ampmleut_n_level}(1+tmp_index_) ...
,str_symbol ...
,'Color' , c_beach__(1+nc,:) ...
,'MarkerSize', markersize_big ...
,'MarkerFaceColor', c_beach__(1+nc,:) ...
,'MarkerEdgeColor', 'k' ...
);
hold off;
axis equal; axis vis3d;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
niteration = 0;%for niteration=0:n_iteration-1;
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

%%%%%%%%;
% Now use ampmleut_template_tree to sort images into groups. ;
%%%%%%%%;
pole_k_c_0_ = cos(euler_azimu_b_).*sin(euler_polar_a_);
pole_k_c_1_ = sin(euler_azimu_b_).*sin(euler_polar_a_);
pole_k_c_2_ = cos(euler_polar_a_);
pole_k_Md__ = [ pole_k_c_0_ , pole_k_c_1_ , pole_k_c_2_ ];
tmp_nlevel = ampmleut_n_level-2;
tree_k_Td__ = ampmleut_template_tree.viewing_pole_k_c_Sd___{1+tmp_nlevel};
index_tree_from_pole_ = knnsearch(tree_k_Td__,pole_k_Md__);
index_tree_from_pole_ = index_tree_from_pole_ - 1;
n_neighborhood = ampmleut_template_tree.n_S_(1+tmp_nlevel);
index_neighborhood_ori_MT__ = cell(n_neighborhood,1);
index_neighborhood_exp_MT__ = cell(n_neighborhood,1);
for nneighborhood=0:n_neighborhood-1;
tmp_nT = nneighborhood;
tmp_avg_d_ = ampmleut_template_tree.index_neighbor_fine_from_coarse_avg___{1+tmp_nlevel}(1+tmp_nT,:);
tmp_rad = ampmleut_template_tree.index_neighbor_fine_from_coarse_rad__{1+tmp_nlevel}(1+tmp_nT);
index_neighborhood_ori_MT_ = efind(index_tree_from_pole_==tmp_nT);
index_neighborhood_ori_MT__{1+nneighborhood} = index_neighborhood_ori_MT_;
pole_k_cnt_Md__ = bsxfun(@minus,pole_k_Md__,tmp_avg_d_);
pole_rad_M_ = sqrt(sum(pole_k_cnt_Md__.^2,2));
index_neighborhood_exp_MT_ = efind(pole_rad_M_<=2*tmp_rad);
index_neighborhood_exp_MT__{1+nneighborhood} = index_neighborhood_exp_MT_;
end;%for nneighborhood=0:n_neighborhood-1;

%%%%%%%%;
% Visualize. ;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(1);clf;figsml;hold on;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_use = 15;
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_use);
end;%for nM=0:n_M-1;
hold off;
axis equal;axis vis3d;
end;%if flag_plot;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(2);clf;figsml;hold on;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_big = 25; markersize_med = 15; markersize_sml = 5;
tmp_nT = floor(1*n_neighborhood/2);
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_sml);
end;%for nM=0:n_M-1;
tmp_index_ = setdiff(index_neighborhood_exp_MT__{1+tmp_nT},index_neighborhood_ori_MT__{1+tmp_nT});
for nl=0:numel(tmp_index_)-1;
nM = tmp_index_(1+nl);
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_med);
end;%for nl=0:numel(tmp_index_)-1;
tmp_index_ = index_neighborhood_ori_MT__{1+tmp_nT};
for nl=0:numel(tmp_index_)-1;
nM = tmp_index_(1+nl);
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_big);
end;%for nl=0:numel(tmp_index_)-1;
hold off;
axis equal;axis vis3d;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% use current euler-angles and displacements to solve for one current model for each complementary neighborhood. ;
%%%%%%%%;
a_UCTF_UX_Y_0qbp_yncx___ = zeros(pm_n_lm_sum,n_CTF_rank,n_neighborhood);
parameter.flag_qbp_vs_lsq = 1;
tmp_t = tic();
for nneighborhood=0:n_neighborhood-1;
if (verbose>2); disp(sprintf(' %% nneighborhood %d/%d',nneighborhood,n_neighborhood)); end;
tmp_nT = nneighborhood;
tmp_index_ = setdiff(0:n_M-1,index_neighborhood_exp_MT__{1+tmp_nT}); %<-- use complement of expanded neighborhood. ;
tmp_n_M = numel(tmp_index_);
[ ...
 parameter ...
,a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) ... 
] = ...
a_UCTF_UX_Y_wrap_ync__0( ...
 parameter ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,tmp_n_M ...
,reshape(UX_M_k_p_wnM___(:,:,1+tmp_index_),[n_w_max*pm_n_k_p_r,tmp_n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__(1+tmp_index_,:) ...
,euler_polar_a_(1+tmp_index_) ...
,euler_azimu_b_(1+tmp_index_) ...
,euler_gamma_z_(1+tmp_index_) ...
,image_I_value_(1+tmp_index_) ...
);
end;%for nneighborhood=0:n_neighborhood-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_0qbp_yncx___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);

%{
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
parameter.flag_qbp_vs_lsq = 1;
tmp_t = tic();
[ ...
 parameter ...
,a_UCTF_UX_Y_0qbp_ync__ ... 
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
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_0qbp_ync__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
%%%%%%%%;
parameter.flag_qbp_vs_lsq = 0;
tmp_t = tic();
[ ...
 parameter ...
,a_UCTF_UX_Y_0lsq_ync__ ... 
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
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_0lsq_ync__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
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
,euler_polar_a_ ...
,euler_azimu_b_ ...
,euler_gamma_z_ ...
,image_I_value_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_ync__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_wrap_ync__0',tmp_t);
%%%%%%%%;
% Now normalize a_CTF_avg_UX_Y_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
a_UCTF_UX_Y_ync__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_ync__);
%}

%%%%%%%%;
% Now normalize a_CTF_avg_UX_Y_. ;
% This step is necessary to prevent the intensity from diverging over successive iterations. ;
%%%%%%%%;
for nneighborhood=0:n_neighborhood-1;
a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood));
end;%for nneighborhood=0:n_neighborhood-1;
%%%%%%%%;
% Use a_UCTF_UX_Y_ync__ as well VSCTF_Mc__ to approximate the image-averaged a_CTF_avg_UX_Y_. ;
% This is not actually used in the calculation, but can be useful for postprocessing. ;
%%%%%%%%;
a_CTF_avg_UX_Y_ = zeros(pm_n_lm_sum,1);
for nneighborhood=0:n_neighborhood-1;
tmp_nT = nneighborhood;
tmp_index_ = setdiff(0:n_M-1,index_neighborhood_exp_MT__{1+tmp_nT}); %<-- use complement of expanded neighborhood. ;
a_CTF_avg_UX_Y_ = a_CTF_avg_UX_Y_ + spharm_normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,mean(a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood)*transpose(VSCTF_Mc__(1+tmp_index_,:)),2));
end;%for nneighborhood=0:n_neighborhood-1;
a_CTF_avg_UX_Y_ = spharm_normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_CTF_avg_UX_Y_);

%%%%%%%%;
% Now store image-parameters. ;
% Be sure to update the translations using the previous step. ;
%%%%%%%%;
a_CTF_avg_UX_Y__(:,1+niteration) = a_CTF_avg_UX_Y_;
euler_polar_a__(:,1+niteration) = euler_polar_a_;
euler_azimu_b__(:,1+niteration) = euler_azimu_b_;
euler_gamma_z__(:,1+niteration) = euler_gamma_z_;
image_delta_x_acc__(:,1+niteration) = image_delta_x_acc_;
image_delta_y_acc__(:,1+niteration) = image_delta_y_acc_;
image_delta_x_upd__(:,1+niteration) = image_delta_x_upd_;
image_delta_y_upd__(:,1+niteration) = image_delta_y_upd_;
image_I_value__(:,1+niteration) = image_I_value_;
%%%%%%%%;
% if flag_MS_vs_SM & order_limit_MS>-1 & niteration<n_iteration/2, then bandlimit a_UCTF_UX_Y_ync__. ;
%%%%%%%%;
if ( (flag_MS_vs_SM==1) & (order_limit_MS>-1) & (niteration<n_iteration/2) );
for nneighborhood=0:n_neighborhood-1;
a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) = a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood).*repmat(Y_lm_cut_,[1,n_CTF_rank]);
end;%for nneighborhood=0:n_neighborhood-1;
end;%if ( (flag_MS_vs_SM==1) & (order_limit_MS>-1) & (niteration<n_iteration/2) );
%%%%%%%%;
% Use each of the individual (excluded) principal-models to align the associated principal-images. ;
% Combines micrographs (i.e., considering n_CTF_rank_use = 1). ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;

tmp_t = tic();
parameter.n_CTF_rank_use = 1;
for nneighborhood=0:n_neighborhood-1;
if (verbose); disp(sprintf(' %% nneighborhood %d/%d',nneighborhood,n_neighborhood)); end;
tmp_nT = nneighborhood;
tmp_index_ = index_neighborhood_ori_MT__{1+tmp_nT}; %<-- use original neighborhood. ;
tmp_n_M = numel(tmp_index_);
[ ...
 parameter ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,X_SM__(:,1+tmp_index_) ...
,delta_x_SM__(:,1+tmp_index_) ...
,delta_y_SM__(:,1+tmp_index_) ...
,gamma_z_SM__(:,1+tmp_index_) ...
,I_value_SM__(:,1+tmp_index_) ...
] = ...
ampmh_X_wrap_wrap_SM__8( ...
 parameter ...
,FTK ...
,n_w_max ...
,l_max_max ...
,pm_n_UX_rank ...
,n_CTF_rank ...
,a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood) ...
,tmp_n_M ...
,CTF_index_(1+tmp_index_) ...
,VSCTF_Mc__(1+tmp_index_,:) ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_) ...
,UX_M_l2_dM__(:,1+tmp_index_) ...
,[] ...
,euler_polar_a_(1+tmp_index_) ...
,euler_azimu_b_(1+tmp_index_) ...
);
end;%for nneighborhood=0:n_neighborhood-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_wrap_wrap_SM__8',tmp_t);
parameter.n_CTF_rank_use = 0; %<-- revert back to full n_CTF_rank. ;

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
if niteration<n_iteration-1;
flag_image_delta_upd__(:,1+niteration) = flag_image_delta_upd_; %<-- these actually refer to end of iteration. ;
image_X_value__(:,1+niteration) = image_X_value_; %<-- these actually refer to end of iteration. ;
image_S_index__(:,1+niteration) = image_S_index_; %<-- these actually refer to end of iteration. ;
end;%if niteration<n_iteration-1;

%%%%%%%%;
% Now return to beginning of loop. ;
%%%%%%%%;
%end;%for niteration=0:n_iteration-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
