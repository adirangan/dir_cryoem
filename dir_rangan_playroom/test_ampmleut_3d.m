% alternating-minimization principled-marching local-exclusion. ;
% defining neighborhoods to include both poles. ;
%clear; test_pm_trpv1_1;

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
if (verbose); disp(sprintf(' %% [entering test_ampmleut_3c]: \n %% alternating-minimization with principal-modes, local-exclusion and updating-translations \n %% defining neighborhoods to include both poles. ;')); end;

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
% Now use ampmleut_template_tree to sort images into groups based on nearest pnode. ;
%%%%%%%%;
pole_k_c_0_ = cos(euler_azimu_b_).*sin(euler_polar_a_);
pole_k_c_1_ = sin(euler_azimu_b_).*sin(euler_polar_a_);
pole_k_c_2_ = cos(euler_polar_a_);
pole_k_Md__ = [ pole_k_c_0_ , pole_k_c_1_ , pole_k_c_2_ ];
tmp_nlevel = ampmleut_n_level-2;
tmp_n_P = ampmleut_template_tree.n_P_(1+tmp_nlevel);
tmp_index_nS_from_nP_ = ampmleut_template_tree.index_nS_from_nP__{1+tmp_nlevel};
tree_k_Pd__ = ampmleut_template_tree.viewing_pole_k_c_Sd___{1+tmp_nlevel}(1+tmp_index_nS_from_nP_,:);
index_tree_from_pole_ = knnsearch([+tree_k_Pd__;-tree_k_Pd__],pole_k_Md__);
index_tree_from_pole_ = index_tree_from_pole_ - 1;
index_tree_from_pole_ = periodize(index_tree_from_pole_,0,tmp_n_P);
n_neighborhood = tmp_n_P;
index_neighborhood_ori_MP__ = cell(n_neighborhood,1);
index_neighborhood_exp_MP__ = cell(n_neighborhood,1);
for nneighborhood=0:n_neighborhood-1;
tmp_nP = nneighborhood;
tmp_nS = tmp_index_nS_from_nP_(1+tmp_nP);
tmp_avg_d_ = ampmleut_template_tree.index_neighbor_fine_from_coarse_avg___{1+tmp_nlevel}(1+tmp_nS,:);
tmp_rad = ampmleut_template_tree.index_neighbor_fine_from_coarse_rad__{1+tmp_nlevel}(1+tmp_nS);
index_neighborhood_ori_MP_ = efind(index_tree_from_pole_==tmp_nP);
index_neighborhood_ori_MP__{1+nneighborhood} = index_neighborhood_ori_MP_;
index_neighborhood_exp_MP_ = [];
pole_k_cnt_Md__ = bsxfun(@minus,pole_k_Md__,+tmp_avg_d_);
pole_rad_M_ = sqrt(sum(pole_k_cnt_Md__.^2,2));
index_neighborhood_exp_MP_ = union(index_neighborhood_exp_MP_,efind(pole_rad_M_<=2*tmp_rad));
pole_k_cnt_Md__ = bsxfun(@minus,pole_k_Md__,-tmp_avg_d_);
pole_rad_M_ = sqrt(sum(pole_k_cnt_Md__.^2,2));
index_neighborhood_exp_MP_ = union(index_neighborhood_exp_MP_,efind(pole_rad_M_<=2*tmp_rad));
index_neighborhood_exp_MP__{1+nneighborhood} = index_neighborhood_exp_MP_;
end;%for nneighborhood=0:n_neighborhood-1;
%%%%%%%%;
% Now define neighborhoods. ;
%%%%%%%%;
index_neighborhood_exc_MP__ = cell(n_neighborhood,1);
for nneighborhood=0:n_neighborhood-1;
index_neighborhood_exc_MP__{1+nneighborhood} = setdiff(0:n_M-1,index_neighborhood_exp_MP__{1+nneighborhood});
end;%for nneighborhood=0:n_neighborhood-1;

%%%%%%%%;
% Visualize. ;
%%%%%%%%;
flag_plot=0;
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
flag_plot=0;
if flag_plot;
figure(2);clf;figbig;
for np=0:6-1;
subplot(2,3,1+np);
tmp_nP = max(0,min(tmp_n_P-1,floor(tmp_n_P*np/(6-1))));
hold on;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
str_symbol_ = {'o','^','s','p','h'}; markersize_big = 25; markersize_med = 15; markersize_sml = 5;
for nM=0:n_M-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_sml);
end;%for nM=0:n_M-1;
tmp_index_ = setdiff(index_neighborhood_exp_MP__{1+tmp_nP},index_neighborhood_ori_MP__{1+tmp_nP});
for nl=0:numel(tmp_index_)-1;
nM = tmp_index_(1+nl);
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_med);
end;%for nl=0:numel(tmp_index_)-1;
tmp_index_ = index_neighborhood_ori_MP__{1+tmp_nP};
for nl=0:numel(tmp_index_)-1;
nM = tmp_index_(1+nl);
nc = max(0,min(n_c_beach-1,floor(n_c_beach*index_tree_from_pole_(1+nM)/n_neighborhood)));
str_symbol = str_symbol_{1+mod(index_tree_from_pole_(1+nM),5)};
plot3(pole_k_c_0_(1+nM),pole_k_c_1_(1+nM),pole_k_c_2_(1+nM),str_symbol,'MarkerFaceColor',c_beach__(1+nc,:),'MarkerEdgeColor','k','MarkerSize',markersize_big);
end;%for nl=0:numel(tmp_index_)-1;
hold off;
axis equal;axis vis3d;
end;%for np=0:6-1;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% use current euler-angles and displacements to solve for one current model for each complementary neighborhood. ;
%%%%%%%%;
c_UCTF_UX_Y_0qbp_yncx___ = zeros(pm_n_lm_sum,n_CTF_rank,n_neighborhood);
tmp_t = tic();
[ ...
 parameter ...
,c_UCTF_UX_Y_0qbp_yncx___ ... 
] = ...
a_UCTF_UX_Y_0qbp_yncx___1( ...
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
,n_neighborhood ...
,index_neighborhood_exc_MP__ ... 
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_0qbp_yncx___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'a_UCTF_UX_Y_0qbp_yncx___1',tmp_t);
%%%%%%%%;
%disp(sprintf(' %% a_UCTF_UX_Y_0qbp_yncx___ vs b_UCTF_UX_Y_0qbp_yncx___: %0.16f',fnorm(a_UCTF_UX_Y_0qbp_yncx___(:,2,:) - b_UCTF_UX_Y_0qbp_yncx___(:,2,:))));
%disp(sprintf(' %% a_UCTF_UX_Y_0qbp_yncx___ vs c_UCTF_UX_Y_0qbp_yncx___: %0.16f',fnorm(a_UCTF_UX_Y_0qbp_yncx___(:,2,:) - c_UCTF_UX_Y_0qbp_yncx___(:,2,:))));
a_UCTF_UX_Y_0qbp_yncx___ = c_UCTF_UX_Y_0qbp_yncx___;

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
tmp_nP = nneighborhood;
tmp_index_ = index_neighborhood_ori_MP__{1+tmp_nP}; %<-- use original neighborhood. ;
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

%%%%%%%%;
% Generate all the templates (for all local neighborhoods) simultaneously, ;
% since generating templates one neighborhood at a time is too slow. ;
% Warning, this can be quite memory intensive if n_neighborhood*n_S is large. ;
% Consider blocking by neighborhood. ;
%%%%%%%%;
tmp_t = tic();
tmp_verbose=1;
[ ...
 UCTF_UX_S_k_q_wncSx____ ...
,~ ...
,n_S ...
] = ...
pm_template_2( ...
 tmp_verbose ...
,l_max_max ...
,pm_n_UX_rank*n_CTF_rank*n_neighborhood ...
,reshape(a_UCTF_UX_Y_0qbp_yncx___(:,:,1:n_neighborhood),[n_lm_max,pm_n_UX_rank*n_CTF_rank*n_neighborhood]) ...
,template_viewing_k_eq_d ...
,-1 ...
,n_w_max ...
);
UCTF_UX_S_k_q_wncSx____ = permute(reshape(UCTF_UX_S_k_q_wncSx____,[n_w_max*pm_n_UX_rank,n_CTF_rank,n_neighborhood,n_S]),[1,2,4,3]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% pm_template_2: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'pm_template_2',tmp_t);
flag_plot=1;
if (flag_plot);
plot(svd(reshape(UCTF_UX_S_k_q_wncSx____,[pm_n_w_sum,n_CTF_rank*n_S*n_neighborhood])),'.');
end;%if (flag_plot);
%%%%%%%%;
% Now convert S_k_p_ to S_k_q_. ;
%%%%%%%%;
tmp_t = tic();
UCTF_UX_S_k_q_wncSx____ = reshape(fft(reshape(UCTF_UX_S_k_q_wncSx____,[pm_n_w_max,pm_n_UX_rank*n_CTF_rank*n_S*n_neighborhood]),[],1)/sqrt(pm_n_w_max),[pm_n_w_sum,n_CTF_rank,n_S,n_neighborhood]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UCTF_UX_S_k_q_wncxS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now compare against a few direct calculations. ;
%%%%%%%%;
flag_check=0;
if flag_check;
n_test=9;
for ntest=0:n_test-1;
nl = max(0,min(n_neighborhood-1,floor(n_neighborhood*ntest/(n_test-1))));
nneighborhood = nl;
disp(sprintf(' %% ntest %d/%d: nl %d nneighborhood %d',ntest,n_test,nl,nneighborhood));
[ ...
 tmp_UCTF_UX_S_k_q_wncS___ ...
] = ...
get_template_1( ...
 tmp_verbose ...
,pm_n_UX_rank*n_CTF_rank ...
,ones(pm_n_UX_rank*n_CTF_rank,1) ...
,1 ...
,ones(pm_n_UX_rank*n_CTF_rank,1) ...
,l_max_max*ones(pm_n_UX_rank*n_CTF_rank,1) ...
,reshape(a_UCTF_UX_Y_0qbp_yncx___(:,:,1+nneighborhood),[pm_n_lm_sum*n_CTF_rank,1]) ...
,template_viewing_k_eq_d ...
,-1 ...
,pm_n_w_max*ones(pm_n_UX_rank*n_CTF_rank,1) ...
);
tmp_UCTF_UX_S_k_q_wncS___ = reshape(tmp_UCTF_UX_S_k_q_wncS___,[pm_n_w_sum,n_CTF_rank,n_S]);
tmp_UCTF_UX_S_k_q_wncS___ = reshape(fft(reshape(tmp_UCTF_UX_S_k_q_wncS___,[pm_n_w_max,pm_n_UX_rank*n_CTF_rank*n_S]),[],1)/sqrt(pm_n_w_max),[pm_n_w_sum,n_CTF_rank,n_S]);
disp(sprintf(' %% tmp_UCTF_UX_S_k_q_wncS___ vs UCTF_UX_S_k_q_wncSx____(:,:,:,1+nneighborhood): %0.16f',fnorm(tmp_UCTF_UX_S_k_q_wncS___ - UCTF_UX_S_k_q_wncSx____(:,:,:,1+nneighborhood))/fnorm(tmp_UCTF_UX_S_k_q_wncS___)));
end;%for ntest=0:n_test-1;
end;%if flag_check;
%%%%%%%%;

tmp_t = tic();
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
] = ...
ampmh_X_single_neighborhood_wrap_SM__9( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_CTF_rank ...
,n_neighborhood ...
,n_S ...
,UCTF_UX_S_k_q_wncSx____ ...
,n_M ...
,CTF_index_ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,index_neighborhood_ori_MP__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% ampmh_X_single_neighborhood_wrap_SM__9: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_single_neighborhood_wrap_SM__9',tmp_t);

%%%%%%%%;
% Finally, compare with a single call to ampmh_X_wrap_wrap_SM__8. ;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.template_viewing_k_eq_d = parameter.template_viewing_k_eq_d;
tmp_t = tic();
ampmh_X_wrap_wrap_SM__8( ...
 tmp_parameter ...
,FTK ...
,n_w_max ...
,l_max_max ...
,pm_n_UX_rank ...
,n_CTF_rank ...
,a_UCTF_UX_Y_0qbp_yncx___(:,:,1+0) ...
,n_M ...
,CTF_index_ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,[] ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% ampmh_X_wrapwrap_SM__8: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_wrap_wrap_SM__8',tmp_t);

%%%%%%%%;
% Verdict: ~160s vs 100s or so. ;
% Not too bad. ;
%%%%%%%%;

%{
for nneighborhood=0:n_neighborhood-1;
index_neighborhood_ori_M_ = index_neighborhood_ori_MP__{1+nneighborhood};
tmp_n_M = numel(index_neighborhood_ori_M_);
if (tmp_n_M>0);
if (verbose); disp(sprintf(' %% nneighborhood %d/%d: tmp_n_M %d',nneighborhood,n_neighborhood,tmp_n_M)); end;
%%%%%%%%;
UCTF_UX_S_k_q_wnSc___ = permute(UCTF_UX_S_k_q_wncSx____(:,:,:,1+nneighborhood),[1,3,2]);
CTF_UX_S_l2_SM__ = transpose(CTF_UX_S_l2_cSx___(1+index_nu_CTF_index_from_nM_(1+index_neighborhood_ori_M_),:,1+nneighborhood));
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_SM__(:,1+index_neighborhood_ori_M_) ...
,delta_x_SM__(:,1+index_neighborhood_ori_M_) ...
,delta_y_SM__(:,1+index_neighborhood_ori_M_) ...
,gamma_z_SM__(:,1+index_neighborhood_ori_M_) ...
] = ...
ampmh_X_single_neighborhood_SM__9( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_CTF_rank ...
,reshape(UCTF_UX_S_k_q_wnSc___,[pm_n_w_max,pm_n_UX_rank,n_S,n_CTF_rank]) ...
,tmp_n_M ...
,CTF_UX_S_l2_SM__ ...
,VSCTF_Mc__(1+index_neighborhood_ori_M_,:) ...
,svd_VUXM_lwnM____(:,:,:,1+index_neighborhood_ori_M_) ...
,UX_M_l2_dM__(:,1+index_neighborhood_ori_M_) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% ampmh_X_single_neighborhood_SM__9: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_single_neighborhood_SM__9',tmp_t);
%%%%%%%%;
flag_check=1;
if (flag_check);
n_test=9;
for ntest=0:n_test-1;
nl = max(0,min(tmp_n_M-1,floor(tmp_n_M*ntest/(n_test-1))));
tmp_nM = index_neighborhood_ori_M_(1+nl);
disp(sprintf(' %% ntest %d/%d: nl %d tmp_nM %d',ntest,n_test,nl,tmp_nM));
CTF_UX_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S);
for nCTF_rank=0:n_CTF_rank-1;
CTF_UX_S_k_q_wnS__ = CTF_UX_S_k_q_wnS__ + UCTF_UX_S_k_q_wnSc___(:,:,1+nCTF_rank)*VSCTF_Mc__(1+tmp_nM,1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_optimize_over_gamma_z = 1;
[ ...
 tmp_parameter ...
,tmp_X_Sm_ ...
,tmp_delta_x_Sm_ ...
,tmp_delta_y_Sm_ ...
,tmp_gamma_z_Sm_ ...
] = ...
ampmh_X_wSM___8( ...
 tmp_parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_SM__(:,1+nl) ...
,1 ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_nM) ...
,UX_M_l2_dM__(:,1+tmp_nM) ...
);
disp(sprintf(' %% tmp_X_Sm_ vs X_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_X_Sm_ - X_SM__(:,1+tmp_nM))/fnorm(tmp_X_Sm_)));
disp(sprintf(' %% tmp_delta_x_Sm_ vs delta_x_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_delta_x_Sm_ - delta_x_SM__(:,1+tmp_nM))/fnorm(tmp_delta_x_Sm_)));
disp(sprintf(' %% tmp_delta_y_Sm_ vs delta_y_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_delta_y_Sm_ - delta_y_SM__(:,1+tmp_nM))/fnorm(tmp_delta_y_Sm_)));
disp(sprintf(' %% tmp_gamma_z_Sm_ vs gamma_z_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_gamma_z_Sm_ - gamma_z_SM__(:,1+tmp_nM))/fnorm(tmp_gamma_z_Sm_)));
end;%for ntest=0:n_test-1;
end;% if (flag_check);
%%%%%%%%;
clear UCTF_UX_S_k_q_wnSc___ CTF_UX_S_l2_SM__;
%%%%%%%%;
end;%if (tmp_n_M>0);
end;%for nneighborhood=0:n_neighborhood-1;
%}

