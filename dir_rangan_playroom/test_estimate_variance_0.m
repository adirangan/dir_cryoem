%%%%%%%%;
% trying to estimate the variance used in scheres bayesian approach. ;
%%%%%%%%;

% clear; test_pm_trpv1c_9b;
CTF_k_p_r_kC__ = CTF_k_p_r__;
a_k_Y_base_yk_ = a_k_Y_quad_; %<-- start with oracle for now. ;
delta_sigma_base = []; %<-- set to 0 for now. ;
tolerance_pm = tolerance_master;
FTK = [];
delta_r_max = 0.1; %<-- default. ;
svd_eps = tolerance_master;
n_delta_v_requested = 0;
flag_image_delta_upd_M_ = ones(n_M,1);
M_k_p_wkM__ = M_k_p__;
image_delta_x_acc_M_ = zeros(n_M,1);
image_delta_y_acc_M_ = zeros(n_M,1);
image_delta_x_upd_M_ = zeros(n_M,1);
image_delta_y_upd_M_ = zeros(n_M,1);
template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% We need to construct an X_wSM__, 
% which records the innerproduct for each template-image-pair, ;
% rather than just for the best pairing. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% First set up CTF-clusters, as in ampmut_wrap_5.m. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

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
if  isempty(a_k_Y_base_yk_);
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
end;%if  isempty(a_k_Y_base_yk_);
%%%%%%%%;
if ~isempty(a_k_Y_base_yk_);
if isempty(delta_sigma_base); delta_sigma_base = 0.0; end; %<-- no translation as default. ;
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
end;%if ~isempty(a_k_Y_base_yk_);
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
if ( isfield(parameter,'fname_pre'));
fname_fig = sprintf('%s_UX_FIGA',parameter.fname_pre);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;figbeach();
ns=0;
%%%%;
subplot(2,3,1+ns);ns=ns+1;
[tmp_UCTF__,tmp_SCTF__,tmp_VCTF__] = svds(CTF_k_p_r_kC__,min(n_CTF,min(n_k_p_r,2)));
tmp_VSCTF__ = tmp_VCTF__*tmp_SCTF__;
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
markersize_use = 16;
str_symbol_ = {'o','^','s','p','h'};
hold on;
for nCTF=0:n_CTF-1;
ncluster = index_ncluster_from_nCTF_(1+nCTF);
nc_beach = max(0,min(n_c_beach,floor(n_c_beach*ncluster/n_cluster)));
str_symbol = str_symbol_{1+mod(ncluster,5)};
plot(tmp_VSCTF__(1+nCTF,1),tmp_VSCTF__(1+nCTF,min(2,size(tmp_VSCTF__,2))),str_symbol,'MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc_beach,:),'MarkerEdgeColor',0.85*[1,1,1]);
end;%for nCTF=0:n_CTF-1;
hold off;
axis equal;
xlabel('pc0');ylabel('pc1');
title('CTF-space colored by cluster');
%%%%;
for nl=0:min(5,n_cluster)-1;
ncluster = nl;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
subplot(2,3,1+ns);ns=ns+1;
tmp_X_ = diag(X_kkc___(:,:,1+ncluster));
tmp_X_ = tmp_X_/fnorm(tmp_X_);
tmp_X_U__ = [ tmp_X_ , UX_knc___(:,1:pm_n_UX_rank,1+ncluster) ];
tmp_lim_ = 1.5*std(tmp_X_U__,1,'all')*[-1,+1];
imagesc(tmp_X_U__,tmp_lim_);
xlabel('mode'); set(gca,'XTick',1:size(tmp_X_U__,2),'XTickLabel',{'X',1:pm_n_UX_rank});
ylabel('k');
title(sprintf('ncluster %d/%d',ncluster,n_cluster));
end;%for nl=0:min(5,n_cluster)-1;
%%%%;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
%%%%;
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%;
end;%if ( isfield(parameter,'fname_pre'));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now find X_wSM___ as in ampmut_5. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

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
if (flag_verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (flag_verbose);

if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

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

a_k_Y_reco_yk_ = a_k_Y_base_yk_; %<-- start with oracle for now. ;

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
,a_k_Y_reco_yk__ ...
,template_viewing_k_eq_d ...
,-1 ...
,n_w_max ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% pm_template_2 (n_S %d): %0.3fs',n_S,tmp_t)); end;
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
,X_wSM___ ...
,delta_x_wSM___ ...
,delta_y_wSM___ ...
,gamma_z_wSM___ ...
,I_value_wSM___ ...
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
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__10',tmp_t);
%%%%%%%%;
[X_SM__,tmp_ij__] = max(reshape(X_wSM___,[n_w_max,n_S*n_M]),[],1);
X_SM__ = reshape(X_SM__,[n_S,n_M]);
tmp_ij__ = reshape(tmp_ij__,[n_S,n_M]);
delta_x_SM__ = delta_x_wSM___(tmp_ij__);
delta_y_SM__ = delta_y_wSM___(tmp_ij__);
gamma_z_SM__ = gamma_z_wSM___(tmp_ij__);
I_value_SM__ = I_value_wSM___(tmp_ij__);
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
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_MS_vs_SM_2',tmp_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now estimate the actual variance, iterating as in scheres 2012. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Note that we assume the variance is fixed across all degrees of freedom within an image. ;
% I.e., fixed across $k$. ;
% But we still imagine $\sigma$ varying across images. ;
%%%%%%%%;
tmp_t = tic();
scheres_sigma_M_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_R2_wS__ = 2 - 2*X_wSM___(:,:,1+nM); %<-- assume images and templates are normalized to have unit norm. ;
tmp_sigma_pre = 1.0;
n_iteration = 32; niteration=0; tmp_tolerance = 1e-6;
flag_continue=1;
while flag_continue;
tmp_p_wS__ = exp(-tmp_R2_wS__/(2*max(1e-12,tmp_sigma_pre)^2));
tmp_p_sum = 2*pi*sum(tmp_p_wS__*template_viewing_weight_all_)/n_w_max;
tmp_p_wS__ = tmp_p_wS__/tmp_p_sum;
if (flag_verbose>3); disp(sprintf(' %% 2*pi*mean(tmp_p_wS__*template_viewing_weight_all_): %0.16f; <-- should be 1',2*pi*mean(tmp_p_wS__*template_viewing_weight_all_))); end;
tmp_sigma_pos = sqrt(0.5*2*pi*sum((tmp_p_wS__.*tmp_R2_wS__)*template_viewing_weight_all_)/n_w_max);
tmp_error = fnorm(tmp_sigma_pos-tmp_sigma_pre)/max(1e-12,fnorm(tmp_sigma_pre));
if (flag_verbose>2); disp(sprintf(' %% niteration %d/%d, tmp_sigma_pre %0.6f tmp_error %0.6f',niteration,n_iteration,tmp_sigma_pre,tmp_error)); end;
flag_continue = (tmp_error> tmp_tolerance) & (niteration<n_iteration);
niteration = niteration+1;
tmp_sigma_pre = tmp_sigma_pos;
end;%while flag_continue;
if (flag_verbose>1); disp(sprintf(' %% nM %d/%d: niteration %d/%d, tmp_sigma_pre %0.6f tmp_error %0.6f',nM,n_M,niteration,n_iteration,tmp_sigma_pre,tmp_error)); end;
scheres_sigma_M_(1+nM) = tmp_sigma_pre;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% scheres_sigma_M_: %0.3fs',tmp_t)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now reconstruct one of the inner products. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
nM = floor(n_M/10);
nS = floor(n_S*2/5);
gamma_z_ = (2*pi)*transpose([0:n_w_max-1])/n_w_max;
%%%%;
M_k_q_ = M_k_q__(:,1+nM);
M_k_p_ = M_k_p__(:,1+nM);
S_k_p_ = S_k_p__(:,1+nS);
tmp_X_0_ = X_wSM___(:,1+nS,1+nM);
%%%%;
tmp_X_1_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_wSM___(1+nw,1+nS,1+nM);
tmp_delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
tmp_delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_delta_x,-tmp_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_);
tmp_X_1 = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
tmp_X_1_(1+nw) = tmp_X_1;
end;%for nw=0:n_w_max-1;
%%%%;
ncluster = index_ncluster_from_nM_(1+nM);
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
%pm_n_UX_rank = n_k_p_r-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
tmp_CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
tmp_CTF_k_p_r_xavg_wk_ = reshape(repmat(transpose(tmp_CTF_k_p_r_xavg_k_),[n_w_max,1]),[n_w_max*n_k_p_r,1]);
tmp_X_2_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_wSM___(1+nw,1+nS,1+nM);
tmp_delta_x = delta_x_wSM___(1+nw,1+nS,1+nM);
tmp_delta_y = delta_y_wSM___(1+nw,1+nS,1+nM);
%T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_gamma_z);
T_k_p_ = S_k_p_;
%T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,tmp_CTF_k_p_r_xavg_wk_,n_w_sum,n_w_sum)*T_k_p_;
N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+tmp_delta_x,+tmp_delta_y);
N_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,N_k_p_,-tmp_gamma_z);
UX_T_k_p_ = reshape(reshape(T_k_p_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
UX_N_k_p_ = reshape(reshape(N_k_p_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
tmp_TT = dot(UX_T_k_p_,UX_T_k_p_);
tmp_MM = dot(UX_N_k_p_,UX_N_k_p_);
tmp_TM = dot(UX_T_k_p_,UX_N_k_p_);
tmp_X_2 = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
tmp_X_2_(1+nw) = tmp_X_2;
end;%for nw=0:n_w_max-1;
%%%%;
flag_disp=1; if flag_disp; figure(1);clf;figsml; plot(gamma_z_,tmp_X_0_,'r^',gamma_z_,tmp_X_1_,'kx',gamma_z_,tmp_X_2_,'go'); end;
%%%%%%%%;
% Verdict: This looks pretty good; tmp_X_0_ is close to tmp_X_2_. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% tmp_X_0_ vs tmp_X_2_: %0.16f',fnorm(tmp_X_0_-tmp_X_2_)/fnorm(tmp_X_0_))); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now try and put one ring of that particular image back onto the sphere. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+tmp_delta_x,+tmp_delta_y);
quad_k_eq_d_ = sqrt(4*pi./n_lm_);
flag_unique_n = 0;
if (numel(unique(l_max_))==1 & numel(unique(n_lm_))==1 & numel(unique(n_w_))==1); flag_unique_n = 1; end;
if ~flag_unique_n; disp(sprintf(' %% Note, flag_unique_n %d',flag_unique_n)); end;
nk_p_r = floor(n_k_p_r/2);
k_p_r = k_p_r_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
N_w_ = N_k_p_(1+index_nw_,:);
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
 data_k_p_polar_a_ ...
,data_k_p_azimu_b_ ...
,data_k_c_0_ ...
,data_k_c_1_ ...
,data_k_c_2_ ...
] = ...
cg_rhs_1( ...
 1 ...
,n_w ...
,template_viewing_polar_a_all_(1+nS) ...
,template_viewing_azimu_b_all_(1+nS) ...
,+0 ...
);
data_k_c_wd__ = [ data_k_c_0_(:) , data_k_c_1_(:) , data_k_c_2_(:) ];
if (flag_verbose); disp(sprintf(' %% data_k_c_0_ vs template_k_c_0__(1+index_nw_,1+nS)/k_p_r: %0.16f',fnorm(data_k_c_0_(:)-template_k_c_0__(1+index_nw_,1+nS)/k_p_r))); end;
if (flag_verbose); disp(sprintf(' %% data_k_c_1_ vs template_k_c_1__(1+index_nw_,1+nS)/k_p_r: %0.16f',fnorm(data_k_c_1_(:)-template_k_c_1__(1+index_nw_,1+nS)/k_p_r))); end;
if (flag_verbose); disp(sprintf(' %% data_k_c_1_ vs template_k_c_1__(1+index_nw_,1+nS)/k_p_r: %0.16f',fnorm(data_k_c_1_(:)-template_k_c_1__(1+index_nw_,1+nS)/k_p_r))); end;
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wd__,'K',1); index_quad_from_data_ = index_quad_from_data_ - 1;
quad_from_data_qw__ = sparse(1+index_quad_from_data_,1:n_w,1,quad_n_all,n_w);
n_quad_from_data_q_ = quad_from_data_qw__*ones(n_w,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wq__ = bsxfun(@rdivide,transpose(quad_from_data_qw__),max(1,transpose(n_quad_from_data_q_)));
CTF_w_ = reshape(CTF_k_p_wkC__(1+index_nw_,1+index_nCTF_from_nM_(1+nM)),[n_w,1]);
CTF2_q_ = quad_from_data_qw__*abs(CTF_w_).^2;

flag_disp=0;
if flag_disp;
figure(1);clf;figsml;
c_use__ = colormap('hsv'); n_c_use = size(c_use__,1);
markersize_sml = 3;
markersize_med = 6;
markersize_big = 9;
hold on;
plot3( quad_k_c_0_all_ , quad_k_c_1_all_ , quad_k_c_2_all_ , 'ko' , 'MarkerFaceColor',0.85*[1,1,1],'MarkerSize',markersize_sml);
for nw=0:size(data_k_c_wd__,1)-1;
nc_use = max(0,min(n_c_use-1,floor(n_c_use*nw/size(data_k_c_wd__,1))));
plot3( data_k_c_0_(1+nw) , data_k_c_1_(1+nw) , data_k_c_2_(1+nw) , 'ko' , 'MarkerFaceColor',c_use__(1+nc_use,:),'MarkerSize',markersize_med);
tmp_index = index_quad_from_data_(1+nw);
plot3( quad_k_c_0_all_(1+tmp_index) , quad_k_c_1_all_(1+tmp_index) , quad_k_c_2_all_(1+tmp_index) , 'ko' , 'MarkerFaceColor',c_use__(1+nc_use,:),'MarkerSize',markersize_big);
end;%for nw=0:size(data_k_c_wd__,1)-1;
hold off;
title(sprintf('nM %d/%d nS %d/%d nk_p_r %d/%d',nM,n_M,nS,n_S,nk_p_r,n_k_p_r),'Interpreter','none');
axis equal;
axis vis3d;
end;%if flag_disp;

N_CTF_w_ = N_w_.*reshape(CTF_w_,[n_w,1]);
quad_from_data_N_CTF_normalized_q_ = (quad_from_data_qw__ * N_CTF_w_)./max(1e-12,CTF2_q_);
a_k_Y_(1+index_Y_) = conj(Ylm_w_yq__)*quad_from_data_N_CTF_normalized_q_;

%%%%%%%%;
% test FTK for displacements. ;
%%%%%%%%;
N_k_p_0_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+tmp_delta_x,+tmp_delta_y);
UX_N_k_p_0_ = reshape(reshape(N_k_p_0_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
image_delta_index = knnsearch([FTK.delta_x_,FTK.delta_y_],[tmp_delta_x,tmp_delta_y]) - 1;
if (flag_verbose); disp(sprintf(' %% (+tmp_delta_x,+tmp_delta_y) (%+0.6f,%+0.6f)',+tmp_delta_x,+tmp_delta_y)); end;
if (flag_verbose); disp(sprintf(' %% (+FTK.delta_x_(1+image_delta_index),+FTK.delta_y_(1+image_delta_index)) (%+0.6f,%+0.6f)',+FTK.delta_x_(1+image_delta_index),+FTK.delta_y_(1+image_delta_index))); end;
svd_VUXM_lwn___ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,1,M_k_q_,pm_n_UX_rank,UX_kn__,X_weight_r_);
svd_VUXM_lwn___ = svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM);
UX_N_k_q_1_ = reshape(FTK.svd_U_d_expiw_s__(1+image_delta_index,:)*reshape(svd_VUXM_lwn___,[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
UX_N_k_p_1_ = reshape(ifft(UX_N_k_q_1_,[],1)*sqrt(n_w_max),[n_w_max*pm_n_UX_rank,1]);
if (flag_verbose); disp(sprintf(' %% UX_N_k_p_0_ vs UX_N_k_p_1_: %0.16f',fnorm(UX_N_k_p_0_ - UX_N_k_p_1_)/fnorm(UX_N_k_p_0_))); end;
XU_UX_N_k_p_0_ = reshape(reshape(UX_N_k_p_0_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
UX_XU_UX_N_k_p_0_ = reshape(reshape(XU_UX_N_k_p_0_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
if (flag_verbose); disp(sprintf(' %% UX_N_k_p_0_ vs UX_XU_UX_N_k_p_0_: %0.16f',fnorm(UX_N_k_p_0_ - UX_XU_UX_N_k_p_0_)/fnorm(UX_N_k_p_0_))); end;
XU_UX_N_k_p_1_ = reshape(reshape(UX_N_k_p_1_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
if (flag_verbose); disp(sprintf(' %% XU_UX_N_k_p_0_ vs XU_UX_N_k_p_1_: %0.16f',fnorm(XU_UX_N_k_p_0_ - XU_UX_N_k_p_1_)/fnorm(XU_UX_N_k_p_0_))); end;
flag_disp=0;
if flag_disp;
figure(1);clf;figmed;
subplot(1,3,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(      N_k_p_0_),[],colormap_beach()); axis image; axisnotick; title('      N_k_p_0_','Interpreter','none');
subplot(1,3,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(XU_UX_N_k_p_0_),[],colormap_beach()); axis image; axisnotick; title('XU_UX_N_k_p_0_','Interpreter','none');
subplot(1,3,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(XU_UX_N_k_p_1_),[],colormap_beach()); axis image; axisnotick; title('XU_UX_N_k_p_1_','Interpreter','none');
end;%if flag_disp;

%%%%%%%%;
% Now, assuming only a single translation is used, ;
% as well as a single ring and a single image and a single template, ;
% practice placing the (rotated) image-values back into the appropriate locations. ;
%%%%%%%%;
tmp_N_w_ = randn(n_w_max,1) + i*randn(n_w_max,1); %<-- random image. ;
tmp_N_q_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_N_w_);
tmp_p_p_ = randn(n_w_max,1); %<-- random weights. ;
tmp_p_q_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_p_p_);
tmp_0_ = zeros(quad_n_all,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_p_p = tmp_p_p_(1+nw);
tmp_0_ = tmp_0_ + quad_from_data_qw__*rotate_p_to_p_fftw(1,n_w_max,n_w_max,tmp_N_w_,-tmp_gamma_z)*tmp_p_p;
end;%for nw=0:n_w_max-1;
tmp_1_ = quad_from_data_qw__*ifft(tmp_N_q_.*conj(tmp_p_q_))*n_w_max;
if (flag_verbose); disp(sprintf(' %% tmp_0_ vs tmp_1_: %0.16f',fnorm(tmp_0_-tmp_1_)/fnorm(tmp_0_))); end;

%%%%%%%%;
% Now try no translation, 1 image, 1 template, ;
% but all rings. ;
%%%%%%%%;
tmp_p_p_wk_ = randn(n_w_sum,1); %<-- random weights. ;
tmp_p_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_p_p_wk_);
tmp_N_k_p_ = M_k_p__(:,1+nM);
tmp_N_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_);
tmp_0__ = zeros(quad_n_all,n_k_p_r);
tmp_1__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
tmp_N_w_ = tmp_N_k_p_(1+index_nw_);
tmp_N_q_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_N_w_);
tmp_p_p_w_ = tmp_p_p_wk_(1+index_nw_);
tmp_p_q_w_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_p_p_w_);
tmp_0_ = zeros(quad_n_all,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_p_p = tmp_p_p_w_(1+nw);
tmp_0_ = tmp_0_ + quad_from_data_qw__*rotate_p_to_p_fftw(1,n_w_max,n_w_max,tmp_N_w_,-tmp_gamma_z)*tmp_p_p;
end;%for nw=0:n_w_max-1;
tmp_0__(:,1+nk_p_r) = tmp_0_;
tmp_1_ = quad_from_data_qw__*ifft(tmp_N_q_.*conj(tmp_p_q_w_))*n_w_max;
tmp_1__(:,1+nk_p_r) = tmp_1_;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_2__ = quad_from_data_qw__*ifft(reshape(tmp_N_k_q_.*conj(tmp_p_q_wk_),[n_w_max,n_k_p_r]),[],1)*n_w_max;
if (flag_verbose); disp(sprintf(' %% tmp_0__ vs tmp_1__: %0.16f',fnorm(tmp_0__-tmp_1__)/fnorm(tmp_0__))); end;
if (flag_verbose); disp(sprintf(' %% tmp_0__ vs tmp_2__: %0.16f',fnorm(tmp_0__-tmp_2__)/fnorm(tmp_0__))); end;

%%%%%%%%;
% Now try 1 translation, 1 image, 1 template, all rings. ;
%%%%%%%%;
index_delta_v = floor(FTK.n_delta_v/2);
tmp_delta_x = FTK.delta_x_(1+index_delta_v);
tmp_delta_y = FTK.delta_y_(1+index_delta_v);
%%%%;
tmp_p_p_w_ = randn(n_w_max,1); %<-- random weights across a single ring. ;
tmp_p_q_w_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_p_p_w_);
tmp_M_k_p_ = M_k_p__(:,1+nM);
tmp_N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,+tmp_delta_x,+tmp_delta_y);
tmp_N_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_);
tmp_UX_N_k_p_0_ = reshape(reshape(tmp_N_k_p_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
tmp_UX_N_k_q_1_ = reshape(FTK.svd_U_d_expiw_s__(1+index_delta_v,:)*reshape(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
tmp_UX_N_k_p_1_ = reshape(ifft(tmp_UX_N_k_q_1_,[],1)*sqrt(n_w_max),[n_w_max*pm_n_UX_rank,1]);
if (flag_verbose); disp(sprintf(' %% tmp_UX_N_k_p_0_ vs tmp_UX_N_k_p_1_: %0.16f',fnorm(tmp_UX_N_k_p_0_-tmp_UX_N_k_p_1_)/fnorm(tmp_UX_N_k_p_0_))); end;
tmp_XU_UX_N_k_p_0_ = reshape(reshape(tmp_UX_N_k_p_0_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
tmp_XU_UX_N_k_p_1_ = reshape(reshape(tmp_UX_N_k_p_1_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
if (flag_verbose); disp(sprintf(' %% tmp_XU_UX_N_k_p_0_ vs tmp_XU_UX_N_k_p_1_: %0.16f',fnorm(tmp_XU_UX_N_k_p_0_-tmp_XU_UX_N_k_p_1_)/fnorm(tmp_XU_UX_N_k_p_0_))); end;
tmp_0__ = zeros(quad_n_all,n_k_p_r);
tmp_1__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
tmp_N_w_ = tmp_XU_UX_N_k_p_1_(1+index_nw_);
tmp_N_q_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_N_w_);
tmp_0_ = zeros(quad_n_all,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_p_p = tmp_p_p_w_(1+nw);
tmp_0_ = tmp_0_ + quad_from_data_qw__*rotate_p_to_p_fftw(1,n_w_max,n_w_max,tmp_N_w_,-tmp_gamma_z)*tmp_p_p;
end;%for nw=0:n_w_max-1;
tmp_0__(:,1+nk_p_r) = tmp_0_;
tmp_1_ = quad_from_data_qw__*ifft(tmp_N_q_.*conj(tmp_p_q_w_))*n_w_max;
tmp_1__(:,1+nk_p_r) = tmp_1_;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_2__ = quad_from_data_qw__*ifft(reshape(tmp_UX_N_k_q_1_.*conj(repmat(tmp_p_q_w_,[1,pm_n_UX_rank])),[n_w_max,pm_n_UX_rank]),[],1)*n_w_max;
tmp_3__ = tmp_2__*transpose(UX_kn__)*diag(X_weight_r_.^(-1));
if (flag_verbose); disp(sprintf(' %% tmp_0__ vs tmp_1__: %0.16f',fnorm(tmp_0__-tmp_1__)/fnorm(tmp_0__))); end;
if (flag_verbose); disp(sprintf(' %% tmp_0__ vs tmp_3__: %0.16f',fnorm(tmp_0__-tmp_3__)/fnorm(tmp_0__))); end;

%%%%%%%%;
% Now try all translations, 1 image, 1 template, all rings. ;
% using k-dependent probabilities. ;
%%%%%%%%;
tmp_M_k_p_ = M_k_p__(:,1+nM);
tmp_p_p_wkd___ = randn(n_w_max,n_k_p_r,FTK.n_delta_v); %<-- random weights across a single ring. ;
tmp_p_q_wkd___ = reshape(interp_p_to_q(n_k_p_r*FTK.n_delta_v,n_w_max*ones(n_k_p_r*FTK.n_delta_v,1),n_w_max*n_k_p_r*FTK.n_delta_v,reshape(tmp_p_p_wkd___,[n_w_max,n_k_p_r*FTK.n_delta_v])),[n_w_max,n_k_p_r,FTK.n_delta_v]);
tmp_0___ = zeros(quad_n_all,n_k_p_r,FTK.n_delta_v);
tmp_1___ = zeros(quad_n_all,n_k_p_r,FTK.n_delta_v);
tmp_3___ = zeros(quad_n_all,n_k_p_r,FTK.n_delta_v);
for ndelta_v=0:FTK.n_delta_v-1;
tmp_delta_x = FTK.delta_x_(1+ndelta_v);
tmp_delta_y = FTK.delta_y_(1+ndelta_v);
%%%%;
tmp_p_p_wk__ = tmp_p_p_wkd___(:,:,1+ndelta_v);
tmp_p_q_wk__ = reshape(interp_p_to_q(n_k_p_r,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,reshape(tmp_p_p_wk__,[n_w_max,n_k_p_r])),[n_w_max,n_k_p_r]);
tmp_N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,+tmp_delta_x,+tmp_delta_y);
tmp_N_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_N_k_p_);
tmp_UX_N_k_p_0_ = reshape(reshape(tmp_N_k_p_,[n_w_max,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max*pm_n_UX_rank,1]);
tmp_UX_N_k_q_1_ = reshape(FTK.svd_U_d_expiw_s__(1+ndelta_v,:)*reshape(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
tmp_UX_N_k_p_1_ = reshape(ifft(tmp_UX_N_k_q_1_,[],1)*sqrt(n_w_max),[n_w_max*pm_n_UX_rank,1]);
if (flag_verbose); disp(sprintf(' %% tmp_UX_N_k_p_0_ vs tmp_UX_N_k_p_1_: %0.16f',fnorm(tmp_UX_N_k_p_0_-tmp_UX_N_k_p_1_)/fnorm(tmp_UX_N_k_p_0_))); end;
tmp_XU_UX_N_k_p_0_ = reshape(reshape(tmp_UX_N_k_p_0_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
tmp_XU_UX_N_k_p_1_ = reshape(reshape(tmp_UX_N_k_p_1_,[n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max*n_k_p_r,1]);
if (flag_verbose); disp(sprintf(' %% tmp_XU_UX_N_k_p_0_ vs tmp_XU_UX_N_k_p_1_: %0.16f',fnorm(tmp_XU_UX_N_k_p_0_-tmp_XU_UX_N_k_p_1_)/fnorm(tmp_XU_UX_N_k_p_0_))); end;
tmp_0__ = zeros(quad_n_all,n_k_p_r);
tmp_1__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
tmp_N_w_ = tmp_XU_UX_N_k_p_1_(1+index_nw_);
tmp_N_q_ = interp_p_to_q(1,n_w_max,n_w_max,tmp_N_w_);
tmp_p_p_w_ = tmp_p_p_wk__(:,1+nk_p_r);
tmp_p_q_w_ = tmp_p_q_wk__(:,1+nk_p_r);
tmp_0_ = zeros(quad_n_all,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_p_p = tmp_p_p_w_(1+nw);
tmp_0_ = tmp_0_ + quad_from_data_qw__*rotate_p_to_p_fftw(1,n_w_max,n_w_max,tmp_N_w_,-tmp_gamma_z)*tmp_p_p;
end;%for nw=0:n_w_max-1;
tmp_0__(:,1+nk_p_r) = tmp_0_;
tmp_1_ = quad_from_data_qw__*ifft(tmp_N_q_.*conj(tmp_p_q_w_))*n_w_max;
tmp_1__(:,1+nk_p_r) = tmp_1_;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_0___(:,:,1+ndelta_v) = tmp_0__;
tmp_1___(:,:,1+ndelta_v) = tmp_1__;
tmp_3__ = quad_from_data_qw__*ifft(reshape((tmp_UX_N_k_q_1_*transpose(UX_kn__)*diag(X_weight_r_.^(-1))).*conj(tmp_p_q_wk__),[n_w_max,n_k_p_r]),[],1)*n_w_max;
tmp_3___(:,:,1+ndelta_v) = tmp_3__;
end;%for ndelta_v=0:FTK.n_delta_v-1;
svd_VUXM_lwn___ = squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM));
svd_VUXM_wnl___ = permute(svd_VUXM_lwn___,[2,3,1]);
svd_U_d_expiw_s_VUXM_dwn___ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_lwn___,[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max;
svd_U_d_expiw_s_VUXM_dwk___ = reshape(reshape(svd_U_d_expiw_s_VUXM_dwn___,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]);
svd_U_d_expiw_s_VUXM_wkd___ = permute(svd_U_d_expiw_s_VUXM_dwk___,[2,3,1]);
tmp_4_qkd___ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,svd_U_d_expiw_s_VUXM_wkd___,reshape(conj(tmp_p_q_wkd___),[n_w_max,n_k_p_r,FTK.n_delta_v])),[n_w_max,n_k_p_r*FTK.n_delta_v]),[],1)*n_w_max,[quad_n_all,n_k_p_r,FTK.n_delta_v]);
tmp_4___ = tmp_4_qkd___;
if (flag_verbose); disp(sprintf(' %% tmp_0___ vs tmp_1___: %0.16f',fnorm(tmp_0___-tmp_1___)/fnorm(tmp_0___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_0___ vs tmp_3___: %0.16f',fnorm(tmp_0___-tmp_3___)/fnorm(tmp_0___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_0___ vs tmp_4___: %0.16f',fnorm(tmp_0___-tmp_4___)/fnorm(tmp_0___))); end;

%%%%;
tmp_t = tic();
svd_VUXM_lwn___ = squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM));
svd_U_d_expiw_s_VUXM_dwn___ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_lwn___,[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max;
svd_U_d_expiw_s_VUXM_dwk___ = reshape(reshape(svd_U_d_expiw_s_VUXM_dwn___,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]);
svd_U_d_expiw_s_VUXM_wkd___ = permute(svd_U_d_expiw_s_VUXM_dwk___,[2,3,1]);
tmp_4_qkd___ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,svd_U_d_expiw_s_VUXM_wkd___,reshape(conj(tmp_p_q_wkd___),[n_w_max,n_k_p_r,FTK.n_delta_v])),[n_w_max,n_k_p_r*FTK.n_delta_v]),[],1)*n_w_max,[quad_n_all,n_k_p_r,FTK.n_delta_v]);
tmp_4___ = tmp_4_qkd___;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_4___: %0.6fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% tmp_0___ vs tmp_4___: %0.16f',fnorm(tmp_0___-tmp_4___)/fnorm(tmp_0___))); end;
%%%%;
tmp_p_q_wdk___ = permute(tmp_p_q_wkd___,[1,3,2]);
tmp_t = tic();
svd_VUXM_lwn___ = squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM));
svd_U_d_expiw_s_VUXM_dwn___ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_lwn___,[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max;
svd_U_d_expiw_s_VUXM_dwk___ = reshape(reshape(svd_U_d_expiw_s_VUXM_dwn___,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]);
svd_U_d_expiw_s_VUXM_wdk___ = permute(svd_U_d_expiw_s_VUXM_dwk___,[2,1,3]);
tmp_5_qdk___ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,svd_U_d_expiw_s_VUXM_wdk___,reshape(conj(tmp_p_q_wdk___),[n_w_max,FTK.n_delta_v,n_k_p_r])),[n_w_max,FTK.n_delta_v*n_k_p_r]),[],1)*n_w_max,[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_5___ = tmp_5_qdk___;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_5___: %0.6fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% tmp_4___ vs tmp_5___: %0.16f',fnorm(tmp_4___-permute(tmp_5___,[1,3,2]))/fnorm(tmp_4___))); end;
%%%%;
tmp_p_q_wdk___ = permute(tmp_p_q_wkd___,[1,3,2]);
tmp_t = tic();
tmp_6___ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,permute(reshape(reshape(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]),[2,1,3]),reshape(conj(tmp_p_q_wdk___),[n_w_max,FTK.n_delta_v,n_k_p_r])),[n_w_max,FTK.n_delta_v*n_k_p_r]),[],1)*n_w_max,[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_6___: %0.6fs',tmp_t)); end;
if (flag_verbose); disp(sprintf(' %% tmp_4___ vs tmp_6___: %0.16f',fnorm(tmp_4___-permute(tmp_6___,[1,3,2]))/fnorm(tmp_4___))); end;
%%%%;
% Verdict: tmp_6___ seems the fastest. ;
%%%%;

 
%%%%%%%%;
% Now treat the denominator. ;
%%%%%%%%;
tmp_p_p_wdk___ = permute(tmp_p_p_wkd___,[1,3,2]);
tmp_8___ = zeros(quad_n_all,n_k_p_r,FTK.n_delta_v);
for ndelta_v=0:FTK.n_delta_v-1;
tmp_p_p_wk__ = tmp_p_p_wkd___(:,:,1+ndelta_v);
tmp_8__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_p_p_w_ = tmp_p_p_wk__(:,1+nk_p_r);
tmp_8_ = zeros(quad_n_all,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_p_p = tmp_p_p_w_(1+nw);
tmp_8_ = tmp_8_ + quad_from_data_qw__*ones(n_w_max,1)*tmp_CTF_k_p_r_xavg_k_(1+nk_p_r)*tmp_p_p;
end;%for nw=0:n_w_max-1;
tmp_8__(:,1+nk_p_r) = tmp_8_;
end;%for nk_p_r=0:n_k_p_r-1;
tmp_8___(:,:,1+ndelta_v) = tmp_8__;
end;%for ndelta_v=0:FTK.n_delta_v-1;
tmp_7___ = reshape(sum(quad_from_data_qw__,2)*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,n_k_p_r,1]),reshape(sum(tmp_p_p_wkd___,1),[1,n_k_p_r,FTK.n_delta_v])),[1,n_k_p_r*FTK.n_delta_v]),[quad_n_all,n_k_p_r,FTK.n_delta_v]);
tmp_9___ = reshape(sum(quad_from_data_qw__,2)*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,n_k_p_r]),reshape(sum(tmp_p_p_wdk___,1),[1,FTK.n_delta_v,n_k_p_r])),[1,FTK.n_delta_v*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_k_p_r]);
if (flag_verbose); disp(sprintf(' %% tmp_8___ vs tmp_7___: %0.16f',fnorm(tmp_8___-tmp_7___)/fnorm(tmp_8___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_8___ vs tmp_9___: %0.16f',fnorm(tmp_8___-permute(tmp_9___,[1,3,2]))/fnorm(tmp_8___))); end;


%%%%%%%%;
% Going forward, need to accumulate the following weights: ;
% w_p = exp(-R2/(2*sigma^2)) ;
% w_n = exp(-R2/(2*sigma^2)) * CTF^1 / sigma^2 ;
% w_d = exp(-R2/(2*sigma^2)) * CTF^2 / sigma^2 ;
% for each image individually (across all templates),
% as these will eventually be normalized per image before being normalized (once more) per quadrature-point. ;
%%%%%%%%;
% trying all translations, some images and templates, all rings. ;
%%%%%%%%;
ncluster=floor(n_cluster/2);
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
tmp_CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
tmp_CTF_k_p_r_xavg_wk_ = reshape(repmat(transpose(tmp_CTF_k_p_r_xavg_k_),[n_w_max,1]),[n_w_max*n_k_p_r,1]);
n_M_sub = min(24,numel(index_nM_from_ncluster_)); index_nM_from_nM_sub_ = index_nM_from_ncluster_(1:n_M_sub);
n_S_sub = min(25,n_S); index_nS_from_nS_sub_ = randperm(n_S,n_S_sub)-1;
tmp_CTF_k_p_wk_ = mean(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+index_nM_from_nM_sub_)),2);
if (flag_verbose); disp(sprintf(' %% tmp_CTF_k_p_r_xavg_wk_ vs tmp_CTF_k_p_wk_: %0.16f',fnorm(tmp_CTF_k_p_r_xavg_wk_ - tmp_CTF_k_p_wk_)/fnorm(tmp_CTF_k_p_wk_))); end;
%%%%;
[ ...
 data_k_p_polar_a_wS__ ...
,data_k_p_azimu_b_wS__ ...
,data_k_c_0_wS__ ...
,data_k_c_1_wS__ ...
,data_k_c_2_wS__ ...
] = ...
cg_rhs_1( ...
 n_S_sub ...
,n_w ...
,template_viewing_polar_a_all_(1+index_nS_from_nS_sub_) ...
,template_viewing_azimu_b_all_(1+index_nS_from_nS_sub_) ...
,+zeros(n_S_sub,1) ...
);
data_k_c_wSd__ = [ data_k_c_0_wS__(:) , data_k_c_1_wS__(:) , data_k_c_2_wS__(:) ];
%%%%;
index_quad_from_data_ = knnsearch(quad_k_c_qd__,data_k_c_wSd__,'K',1)-1;
quad_from_data_qwS__ = sparse(1+index_quad_from_data_,1:n_w*n_S_sub,1,quad_n_all,n_w*n_S_sub);
n_quad_from_data_q_ = quad_from_data_qwS__*ones(n_w*n_S_sub,1); %<-- number of data-points per quadrature-point. ;
data_from_quad_wSq__ = bsxfun(@rdivide,transpose(quad_from_data_qwS__),max(1,transpose(n_quad_from_data_q_)));
nk_p_r = n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
if (flag_verbose); disp(sprintf(' %% data_k_c_0_wS__vs template_k_c_0__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r: %0.16f',fnorm(data_k_c_0_wS__-template_k_c_0__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r))); end;
if (flag_verbose); disp(sprintf(' %% data_k_c_1_wS__vs template_k_c_1__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r: %0.16f',fnorm(data_k_c_1_wS__-template_k_c_1__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r))); end;
if (flag_verbose); disp(sprintf(' %% data_k_c_1_wS__vs template_k_c_1__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r: %0.16f',fnorm(data_k_c_1_wS__-template_k_c_1__(1+index_nw_,1+index_nS_from_nS_sub_)/k_p_r))); end;
%%%%;

tmp_sigma = 1.5;
tmp_R2_dwSM____ = rand(FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub);

%%%%%%%%;
% calculation one image+template at a time. ;
%%%%%%%%;
tmp_numerator_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
tmp_denomator_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
expR2_sum_M_ = zeros(n_M_sub,1);
tmp_t = tic();
for nS_sub=0:n_S_sub-1;
quad_from_data_qw__ = quad_from_data_qwS__(:,1 + nS_sub*n_w_max + [0:n_w_max-1]);
for nM_sub=0:n_M_sub-1;
nM = index_nM_from_nM_sub_(1+nM_sub);
tmp_R2_wd__ = permute(tmp_R2_dwSM____(:,:,1+nS_sub,1+nM_sub),[2,1,3,4]);
tmp_expR2_p_wd__ = exp(-tmp_R2_wd__/(2*tmp_sigma^2));
tmp_p_p_wdk___ = bsxfun(@times,reshape(tmp_expR2_p_wd__,[n_w_max,FTK.n_delta_v,1]),reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wdk___ = reshape(interp_p_to_q(FTK.n_delta_v*n_k_p_r,n_w_max*ones(FTK.n_delta_v*n_k_p_r,1),n_w_max*FTK.n_delta_v*n_k_p_r,tmp_p_p_wdk___),[n_w_max,FTK.n_delta_v,n_k_p_r]);
expR2_sum_M_(1+nM_sub) = expR2_sum_M_(1+nM_sub) + sum(tmp_expR2_p_wd__,'all');
tmp_numerator_qdk___ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,permute(reshape(reshape(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]),[2,1,3]),reshape(conj(tmp_p_q_wdk___),[n_w_max,FTK.n_delta_v,n_k_p_r])),[n_w_max,FTK.n_delta_v*n_k_p_r]),[],1)*n_w_max,[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_denomator_qdk___ = reshape(sum(quad_from_data_qw__,2)*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,n_k_p_r]),reshape(sum(tmp_p_p_wdk___,1),[1,FTK.n_delta_v,n_k_p_r])),[1,FTK.n_delta_v*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_numerator_qkM___(:,:,1+nM_sub) = tmp_numerator_qkM___(:,:,1+nM_sub) + permute(sum(tmp_numerator_qdk___,2),[1,3,2]);
tmp_denomator_qkM___(:,:,1+nM_sub) = tmp_denomator_qkM___(:,:,1+nM_sub) + permute(sum(tmp_denomator_qdk___,2),[1,3,2]);
end;%for nM_sub=0:n_M_sub-1;
end;%for nS_sub=0:n_S_sub-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ tmp_denomator_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(n_M_sub*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;

%%%%%%%%;
% try and make this more efficient by combining images. ;
%%%%%%%%;
tmp_numerator_0_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
tmp_denomator_0_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
tmp_numerator_0_qMk___ = zeros(quad_n_all,n_M_sub,n_k_p_r);
tmp_denomator_0_qMk___ = zeros(quad_n_all,n_M_sub,n_k_p_r);
expR2_sum_0_M_ = zeros(n_M_sub,1);
tmp_t = tic();
for nS_sub=0:n_S_sub-1;
quad_from_data_qw__ = quad_from_data_qwS__(:,1 + nS_sub*n_w_max + [0:n_w_max-1]);
tmp_R2_wdM___ = permute(tmp_R2_dwSM____(:,:,1+nS_sub,:),[2,1,4,3]);
tmp_expR2_p_wdM___ = exp(-tmp_R2_wdM___/(2*tmp_sigma^2));
tmp_p_p_wdMk____ = bsxfun(@times,reshape(tmp_expR2_p_wdM___,[n_w_max,FTK.n_delta_v,n_M_sub,1]),reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wdMk____ = reshape(interp_p_to_q(FTK.n_delta_v*n_M_sub*n_k_p_r,n_w_max*ones(FTK.n_delta_v*n_M_sub*n_k_p_r,1),n_w_max*FTK.n_delta_v*n_M_sub*n_k_p_r,tmp_p_p_wdMk____),[n_w_max,FTK.n_delta_v,n_M_sub,n_k_p_r]);
expR2_sum_0_M_ = expR2_sum_0_M_ + reshape(sum(reshape(tmp_expR2_p_wdM___,[n_w_max*FTK.n_delta_v,n_M_sub]),1),[n_M_sub,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_qw__*ifft(reshape(bsxfun(@times,reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_nM_sub_),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M_sub]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M_sub])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*n_M_sub,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,n_M_sub,n_k_p_r]),reshape(conj(tmp_p_q_wdMk____),[n_w_max,FTK.n_delta_v,n_M_sub,n_k_p_r])),[n_w_max,FTK.n_delta_v*n_M_sub*n_k_p_r]),[],1)*n_w_max,[quad_n_all,FTK.n_delta_v,n_M_sub,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(sum(quad_from_data_qw__,2)*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,n_k_p_r]),reshape(sum(tmp_p_p_wdMk____,1),[1,FTK.n_delta_v,n_M_sub,n_k_p_r])),[1,FTK.n_delta_v*n_M_sub*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_M_sub,n_k_p_r]);
tmp_numerator_0_qMk___ = tmp_numerator_0_qMk___ + squeeze(sum(tmp_numerator_qdMk____,2));
tmp_denomator_0_qMk___ = tmp_denomator_0_qMk___ + squeeze(sum(tmp_denomator_qdMk____,2));
end;%for nS_sub=0:n_S_sub-1;
tmp_numerator_0_qkM___ = permute(tmp_numerator_0_qMk___,[1,3,2]);
tmp_denomator_0_qkM___ = permute(tmp_denomator_0_qMk___,[1,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_qdMk____ tmp_denomator_qdMk____: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(n_M_sub*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
if (flag_verbose); disp(sprintf(' %% expR2_sum_M_ vs expR2_sum_0_M_: %0.16f',fnorm(expR2_sum_M_ - expR2_sum_0_M_)/fnorm(expR2_sum_M_))); end;
if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ vs tmp_numerator_0_qkM___: %0.16f',fnorm(tmp_numerator_qkM___ - tmp_numerator_0_qkM___)/fnorm(tmp_numerator_qkM___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_denomator_qkM___ vs tmp_denomator_0_qkM___: %0.16f',fnorm(tmp_denomator_qkM___ - tmp_denomator_0_qkM___)/fnorm(tmp_denomator_qkM___))); end;

%%%%%%%%;
% try and make this more efficient by combining templates. ;
%%%%%%%%;
tmp_numerator_1_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
tmp_denomator_1_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
expR2_sum_1_M_ = zeros(n_M_sub,1);
tmp_t = tic();
for nM_sub=0:n_M_sub-1;
nM = index_nM_from_nM_sub_(1+nM_sub);
tmp_R2_wSd___ = permute(tmp_R2_dwSM____(:,:,:,1+nM_sub),[2,3,1]);
tmp_expR2_p_wSd___ = exp(-tmp_R2_wSd___/(2*tmp_sigma^2));
tmp_p_p_wSdk____ = bsxfun(@times,reshape(tmp_expR2_p_wSd___,[n_w_max,n_S_sub,FTK.n_delta_v,1]),reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wSdk____ = reshape(interp_p_to_q(n_S_sub*FTK.n_delta_v*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*n_k_p_r,tmp_p_p_wSdk____),[n_w_max,n_S_sub,FTK.n_delta_v,n_k_p_r]);
expR2_sum_1_M_(1+nM_sub) = expR2_sum_1_M_(1+nM_sub) + sum(tmp_expR2_p_wSd___,'all');
tmp_numerator_qdk___ = reshape(quad_from_data_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(permute(reshape(reshape(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+nM)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank])*n_w_max,[FTK.n_delta_v*n_w_max,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[FTK.n_delta_v,n_w_max,n_k_p_r]),[2,1,3]),[n_w_max,1,FTK.n_delta_v,n_k_p_r]),reshape(conj(tmp_p_q_wSdk____),[n_w_max,n_S_sub,FTK.n_delta_v,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_denomator_qdk___ = reshape(quad_from_data_qwS__*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdk____,1),[n_w_max,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_k_p_r]);
tmp_numerator_1_qkM___(:,:,1+nM_sub) = tmp_numerator_1_qkM___(:,:,1+nM_sub) + permute(sum(tmp_numerator_qdk___,2),[1,3,2]);
tmp_denomator_1_qkM___(:,:,1+nM_sub) = tmp_denomator_1_qkM___(:,:,1+nM_sub) + permute(sum(tmp_denomator_qdk___,2),[1,3,2]);
end;%for nM_sub=0:n_M_sub-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_1_qkM___ tmp_denomator_1_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(n_M_sub*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
if (flag_verbose); disp(sprintf(' %% expR2_sum_M_ vs expR2_sum_1_M_: %0.16f',fnorm(expR2_sum_M_ - expR2_sum_1_M_)/fnorm(expR2_sum_M_))); end;
if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ vs tmp_numerator_1_qkM___: %0.16f',fnorm(tmp_numerator_qkM___ - tmp_numerator_1_qkM___)/fnorm(tmp_numerator_qkM___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_denomator_qkM___ vs tmp_denomator_1_qkM___: %0.16f',fnorm(tmp_denomator_qkM___ - tmp_denomator_1_qkM___)/fnorm(tmp_denomator_qkM___))); end;

%%%%%%%%;
% try and make this more efficient by combining templates and images. ;
%%%%%%%%;
tmp_numerator_2_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
tmp_denomator_2_qkM___ = zeros(quad_n_all,n_k_p_r,n_M_sub);
expR2_sum_2_M_ = zeros(n_M_sub,1);
tmp_t = tic();
tmp_R2_wSdM____ = permute(tmp_R2_dwSM____,[2,3,1,4]);
tmp_expR2_p_wSdM____ = exp(-tmp_R2_wSdM____/(2*tmp_sigma^2));
tmp_p_p_wSdMk_____ = bsxfun(@times,reshape(tmp_expR2_p_wSdM____,[n_w_max,n_S_sub,FTK.n_delta_v,n_M_sub,1]),reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wSdMk_____ = reshape(interp_p_to_q_block_0(n_S_sub*FTK.n_delta_v*n_M_sub*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*n_M_sub*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*n_M_sub*n_k_p_r,tmp_p_p_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,n_M_sub,n_k_p_r]);
expR2_sum_2_M_ = expR2_sum_2_M_ + reshape(sum(reshape(tmp_expR2_p_wSdM____,[n_w_max*n_S_sub*FTK.n_delta_v,n_M_sub]),1),[n_M_sub,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_nM_sub_)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M_sub]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M_sub])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*n_M_sub,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,n_M_sub,n_k_p_r]),[n_w_max,1,FTK.n_delta_v,n_M_sub,n_k_p_r]),reshape(conj(tmp_p_q_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,n_M_sub,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*n_M_sub*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*n_M_sub*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_M_sub,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(quad_from_data_qwS__*reshape(bsxfun(@times,reshape(tmp_CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdMk_____,1),[n_w_max,1,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*n_M_sub*n_k_p_r]),[quad_n_all,FTK.n_delta_v,n_M_sub,n_k_p_r]);
tmp_numerator_2_qkM___ = tmp_numerator_2_qkM___ + permute(sum(tmp_numerator_qdMk____,2),[1,4,3,2]);
tmp_denomator_2_qkM___ = tmp_denomator_2_qkM___ + permute(sum(tmp_denomator_qdMk____,2),[1,4,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_2_qkM___ tmp_denomator_2_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(n_M_sub*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
if (flag_verbose); disp(sprintf(' %% expR2_sum_M_ vs expR2_sum_2_M_: %0.16f',fnorm(expR2_sum_M_ - expR2_sum_2_M_)/fnorm(expR2_sum_M_))); end;
if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ vs tmp_numerator_2_qkM___: %0.16f',fnorm(tmp_numerator_qkM___ - tmp_numerator_2_qkM___)/fnorm(tmp_numerator_qkM___))); end;
if (flag_verbose); disp(sprintf(' %% tmp_denomator_qkM___ vs tmp_denomator_2_qkM___: %0.16f',fnorm(tmp_denomator_qkM___ - tmp_denomator_2_qkM___)/fnorm(tmp_denomator_qkM___))); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now try and extract actual X_dwSM____. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_sigma = 0.1;
%%%%%%%%;
[ ...
 data_k_p_polar_a_wS__ ...
,data_k_p_azimu_b_wS__ ...
,data_k_c_0_wS__ ...
,data_k_c_1_wS__ ...
,data_k_c_2_wS__ ...
] = ...
cg_rhs_1( ...
 n_S ...
,n_w ...
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
n_S_per_Sbatch = 24;
%%%%%%%%;
ncluster = floor(n_cluster/2);
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: tmp_n_M %d',ncluster,n_cluster,tmp_n_M)); end;
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
%pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
%%%%%%%%;
X_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
delta_x_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
delta_y_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
gamma_z_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
I_value_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
tmp_numerator_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
tmp_denomator_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
tmp_expR2_sum_M_ = zeros(tmp_n_M,1);
%%%%%%%%;
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>0); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
%if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_index_ = reshape(repmat(transpose(0:n_w_max-1),[1,n_S_sub]) + n_w_max*repmat(reshape(index_S_in_Sbatch_,[1,n_S_sub]),[n_w_max,1]),[n_w_max*n_S_sub,1]);
quad_from_data_sub_qwS__ = quad_from_data_qwS__(:,1+tmp_index_);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_sub_wSM___(:,1+index_S_in_Sbatch_,:) ...
,delta_x_sub_wSM___(:,1+index_S_in_Sbatch_,:) ...
,delta_y_sub_wSM___(:,1+index_S_in_Sbatch_,:) ...
,gamma_z_sub_wSM___(:,1+index_S_in_Sbatch_,:) ...
,I_value_sub_wSM___(:,1+index_S_in_Sbatch_,:) ...
,tmp_X_sub_dwSM____ ...
] = ...
ampmh_X_single_cluster_wSM___11( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S_sub ...
,reshape(CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1+index_S_in_Sbatch_),[pm_n_w_max,pm_n_UX_rank,n_S_sub]) ...
,CTF_UX_S_l2_S_(1+index_S_in_Sbatch_) ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
tmp_R2_wSdM____ = permute(2 - 2*tmp_X_sub_dwSM____,[2,3,1,4]);
tmp_expR2_p_wSdM____ = exp(-tmp_R2_wSdM____/(2*tmp_sigma^2));
tmp_p_p_wSdMk_____ = bsxfun(@times,reshape(tmp_expR2_p_wSdM____,[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,1]),reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wSdMk_____ = reshape(interp_p_to_q_block_0(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,tmp_p_p_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_expR2_sum_M_ = tmp_expR2_sum_M_ + reshape(sum(reshape(tmp_expR2_p_wSdM____,[n_w_max*n_S_sub*FTK.n_delta_v,tmp_n_M]),1),[tmp_n_M,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_ncluster_)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*tmp_n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,tmp_n_M])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*tmp_n_M,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,tmp_n_M,n_k_p_r]),[n_w_max,1,FTK.n_delta_v,tmp_n_M,n_k_p_r]),reshape(conj(tmp_p_q_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(bsxfun(@times,reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdMk_____,1),[n_w_max,1,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_numerator_qkM___ = tmp_numerator_qkM___ + permute(sum(tmp_numerator_qdMk____,2),[1,4,3,2]);
tmp_denomator_qkM___ = tmp_denomator_qkM___ + permute(sum(tmp_denomator_qdMk____,2),[1,4,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ tmp_denomator_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(tmp_n_M*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
clear tmp_X_sub_dwSM____;
clear tmp_R2_wSdM____;
clear tmp_expR2_p_wSdM____;
clear tmp_p_p_wSdMk_____;
clear tmp_p_q_wSdMk_____;
clear tmp_numerator_qdMk____;
clear tmp_denomator_qdMk____;
%%%%%%%%;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% now use numerator and denomator to construct data on the sphere. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
size(tmp_numerator_qkM___);
size(tmp_denomator_qkM___);
size(tmp_expR2_sum_M_);
tmp_frac_qk__ = bsxfun(@rdivide,sum(bsxfun(@rdivide,tmp_numerator_qkM___,reshape(tmp_expR2_sum_M_,[1,1,tmp_n_M])),3),max(sum(bsxfun(@rdivide,tmp_denomator_qkM___,reshape(tmp_expR2_sum_M_,[1,1,tmp_n_M])),3),tolerance_denom));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Compare with standard maximum-likelihood estimate. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nS=0:n_S-1; for nM_sub=0:tmp_n_M-1;
[~,tmp_ij] = max(X_sub_wSM___(:,1+nS,1+nM_sub)); tmp_nw = tmp_ij-1;
X_sub_SM__(1+nS,1+nM_sub) = X_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
delta_x_sub_SM__(1+nS,1+nM_sub) = delta_x_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
delta_y_sub_SM__(1+nS,1+nM_sub) = delta_y_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
gamma_z_sub_SM__(1+nS,1+nM_sub) = gamma_z_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
I_value_sub_SM__(1+nS,1+nM_sub) = I_value_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
end;end;%for nS=0:n_S-1; for nM_sub=0:tmp_n_M-1;
%%%%%%%%;
tmp_t = tic();
parameter.flag_MS_vs_SM = 0;
[ ...
 parameter ...
,euler_polar_a_sub_M_ ...
,euler_azimu_b_sub_M_ ...
,euler_gamma_z_sub_M_ ...
,image_delta_x_sub_M_ ...
,image_delta_y_sub_M_ ...
,image_I_value_sub_M_ ...
,image_X_value_sub_M_ ...
,image_S_index_sub_M_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter ...
,n_w_max ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,tmp_n_M ...
,X_sub_SM__ ...
,delta_x_sub_SM__ ...
,delta_y_sub_SM__ ...
,gamma_z_sub_SM__ ...
,I_value_sub_SM__ ...
);
tmp_str = 'SM'; if (parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
%%%%%%%%;
T_k_p_sub_wkM__ = M_k_p_wkM__(:,1+index_nM_from_ncluster_);
for nM_sub=0:tmp_n_M-1;
T_k_p_sub_wkM__(:,1+nM_sub) = image_I_value_sub_M_(1+nM_sub) * transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_sub_wkM__(:,1+nM_sub),+image_delta_x_sub_M_(1+nM_sub),+image_delta_y_sub_M_(1+nM_sub));
end;%for nM_sub=0:tmp_n_M-1;
%%%%%%%%;
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
CTF_sub_wMn__ = reshape(permute(reshape(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),[n_w,n_k_p_r,tmp_n_M]),[1,3,2]),[n_w*tmp_n_M,n_k_p_r]);
CTF2_sub_qk__ = quad_from_data_sub_qwM__*abs(CTF_sub_wMn__).^2;
quad_from_data_sub_T_CTF_sub_normalized_qk__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
T_sub_wM__ = T_k_p_sub_wkM__(1+index_nw_,:);
CTF_sub_wM_ = CTF_sub_wMn__(:,1+nk_p_r);
CTF2_sub_q_ = CTF2_sub_qk__(:,1+nk_p_r);
T_CTF_sub_wM__ = T_sub_wM__.*reshape(CTF_sub_wM_,[n_w,tmp_n_M]);
quad_from_data_sub_T_CTF_sub_normalized_q_ = (quad_from_data_sub_qwM__ * T_CTF_sub_wM__(:))./max(1e-12,CTF2_sub_q_);
quad_from_data_sub_T_CTF_sub_normalized_qk__(:,1+nk_p_r) = quad_from_data_sub_T_CTF_sub_normalized_q_;
end;%for nk_p_r=0:n_k_p_r-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Verdict: Hard to tell. ;
% Test again using only, say, a few images and 4 templates. ;
% Need to estimate the X_SM__ first ;
% in order to set an image-dependent lower-bound for R2_. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_verbose = 1; flag_disp = flag_verbose; nf=0;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% testing sheres-style qbp. ;')); end;
ncluster = floor(n_cluster/2);
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
%index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster}(1); %<-- pick a single image. ;
%tmp_n_M = numel(index_nM_from_ncluster_); %<-- pick a single image. ;
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: tmp_n_M %d',ncluster,n_cluster,tmp_n_M)); end;
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_k_p_r)); end;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
%pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
n_S_sub = min(4,n_S);
index_nS_from_nS_sub_ = 0 + transpose([0:n_S_sub-1]);
S_sub_k_q_wSk___ = permute(reshape(S_k_q_wkS__(:,1+index_nS_from_nS_sub_),[n_w_max,n_k_p_r,n_S_sub]),[1,3,2]);
CTF_UX_S_sub_k_q_wnS__ = zeros(pm_n_w_sum,n_S_sub);
CTF_UX_S_sub_k_q_wnS__(1:pm_n_w_sum,1:n_S_sub) = reshape(permute(reshape(reshape(S_sub_k_q_wSk___,[n_w_max*n_S_sub,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S_sub,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S_sub]);
CTF_UX_S_sub_l2_S_ = sum(abs(CTF_UX_S_sub_k_q_wnS__).^2,1)/n_w_max;
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(k_p_r_,CTF_k_p_r_xavg_k_,'ko-'); xlabel('k');ylabel('CTF'); axis tight; title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%;
X_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
delta_x_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
delta_y_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
gamma_z_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
I_value_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
%%%%%%%%;
tmp_index_ = reshape(repmat(transpose(0:n_w_max-1),[1,n_S_sub]) + n_w_max*repmat(reshape(index_nS_from_nS_sub_,[1,n_S_sub]),[n_w_max,1]),[n_w_max*n_S_sub,1]);
quad_from_data_sub_qwS__ = quad_from_data_qwS__(:,1+tmp_index_);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_sub_wSM___ ...
,delta_x_sub_wSM___ ...
,delta_y_sub_wSM___ ...
,gamma_z_sub_wSM___ ...
,I_value_sub_wSM___ ...
] = ...
ampmh_X_single_cluster_wSM___11( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S_sub ...
,reshape(CTF_UX_S_sub_k_q_wnS__,[pm_n_w_max,pm_n_UX_rank,n_S_sub]) ...
,CTF_UX_S_sub_l2_S_ ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
R2_sub_wSM___ = 2-2*X_sub_wSM___;
R2_sub_min_M_ = reshape(min(R2_sub_wSM___,[],[1,2]),[tmp_n_M,1]);
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:tmp_n_M,R2_sub_min_M_,'ko'); xlim([1,tmp_n_M]); ylim([0,4]); xlabel('image index'); ylabel('R2'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_sigma = 0.001;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now taking small temperature tmp_sigma %0.6f. ;',tmp_sigma)); end;
%%%%%%%%;
[ ...
 data_sub_k_p_polar_a_wS__ ...
,data_sub_k_p_azimu_b_wS__ ...
,data_sub_k_c_0_wS__ ...
,data_sub_k_c_1_wS__ ...
,data_sub_k_c_2_wS__ ...
] = ...
cg_rhs_1( ...
 n_S_sub ...
,n_w ...
,template_viewing_polar_a_all_(1+index_nS_from_nS_sub_) ...
,template_viewing_azimu_b_all_(1+index_nS_from_nS_sub_) ...
,+zeros(n_S_sub,1) ...
);
data_sub_k_c_wSd__ = [ data_sub_k_c_0_wS__(:) , data_sub_k_c_1_wS__(:) , data_sub_k_c_2_wS__(:) ];
%%%%;
index_quad_from_data_sub_ = knnsearch(quad_k_c_qd__,data_sub_k_c_wSd__,'K',1)-1;
quad_from_data_sub_qwS__ = sparse(1+index_quad_from_data_sub_,1:n_w*n_S_sub,1,quad_n_all,n_w*n_S_sub);
n_quad_from_data_sub_q_ = quad_from_data_sub_qwS__*ones(n_w*n_S_sub,1); %<-- number of data-points per quadrature-point. ;
data_sub_from_quad_wSq__ = bsxfun(@rdivide,transpose(quad_from_data_sub_qwS__),max(1,transpose(n_quad_from_data_sub_q_)));
%%%%%%%%;
X_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
delta_x_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
delta_y_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
gamma_z_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
I_value_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
tmp_numerator_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
tmp_denomator_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
tmp_expR2_sum_M_ = zeros(tmp_n_M,1);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_sub_wSM___ ...
,delta_x_sub_wSM___ ...
,delta_y_sub_wSM___ ...
,gamma_z_sub_wSM___ ...
,I_value_sub_wSM___ ...
,tmp_X_sub_dwSM____ ...
] = ...
ampmh_X_single_cluster_wSM___11( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S_sub ...
,reshape(CTF_UX_S_sub_k_q_wnS__,[pm_n_w_max,pm_n_UX_rank,n_S_sub]) ...
,CTF_UX_S_sub_l2_S_ ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
tmp_R2_wSdM____ = bsxfun(@minus,permute(2 - 2*tmp_X_sub_dwSM____,[2,3,1,4]),reshape(R2_sub_min_M_,[1,1,1,tmp_n_M]));
tmp_expR2_p_wSdM____ = exp(-tmp_R2_wSdM____/(2*tmp_sigma^2));
tmp_p_p_wSdMk_____ = bsxfun(@times,reshape(tmp_expR2_p_wSdM____,[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,1]),reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wSdMk_____ = reshape(interp_p_to_q_block_0(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,tmp_p_p_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_expR2_sum_M_ = tmp_expR2_sum_M_ + reshape(sum(reshape(tmp_expR2_p_wSdM____,[n_w_max*n_S_sub*FTK.n_delta_v,tmp_n_M]),1),[tmp_n_M,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_ncluster_)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*tmp_n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,tmp_n_M])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*tmp_n_M,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,tmp_n_M,n_k_p_r]),[n_w_max,1,FTK.n_delta_v,tmp_n_M,n_k_p_r]),reshape(conj(tmp_p_q_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(bsxfun(@times,reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdMk_____,1),[n_w_max,1,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_numerator_qkM___ = tmp_numerator_qkM___ + permute(sum(tmp_numerator_qdMk____,2),[1,4,3,2]);
tmp_denomator_qkM___ = tmp_denomator_qkM___ + permute(sum(tmp_denomator_qdMk____,2),[1,4,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmp_numerator_qkM___ tmp_denomator_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(tmp_n_M*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
clear tmp_X_sub_dwSM____;
clear tmp_R2_wSdM____;
clear tmp_expR2_p_wSdM____;
clear tmp_p_p_wSdMk_____;
clear tmp_p_q_wSdMk_____;
clear tmp_numerator_qdMk____;
clear tmp_denomator_qdMk____;
%%%%%%%%;
% now use numerator and denomator to construct data on the sphere. ;
%%%%%%%%;
size(tmp_numerator_qkM___);
size(tmp_denomator_qkM___);
size(tmp_expR2_sum_M_);
tolerance_denom = 1e-12;
tmp_frac_qk__ = bsxfun(@rdivide,sum(bsxfun(@rdivide,tmp_numerator_qkM___,reshape(tmp_expR2_sum_M_,[1,1,tmp_n_M])),3),max(sum(bsxfun(@rdivide,tmp_denomator_qkM___,reshape(tmp_expR2_sum_M_,[1,1,tmp_n_M])),3),tolerance_denom));
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:tmp_n_M,tmp_expR2_sum_M_,'ko'); xlim([1,tmp_n_M]); ylim([0,4]); xlabel('image index'); ylabel('tmp_expR2_sum_M_','Interpreter','none'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now calculating maximum-likelihood qbp. ;')); end;
X_sub_SM__ = zeros(n_S_sub,tmp_n_M);
delta_x_sub_SM__ = zeros(n_S_sub,tmp_n_M);
delta_y_sub_SM__ = zeros(n_S_sub,tmp_n_M);
gamma_z_sub_SM__ = zeros(n_S_sub,tmp_n_M);
I_value_sub_SM__ = zeros(n_S_sub,tmp_n_M);
for nS_sub=0:n_S_sub-1; for nM_sub=0:tmp_n_M-1;
[~,tmp_ij] = max(X_sub_wSM___(:,1+nS_sub,1+nM_sub)); tmp_nw = tmp_ij-1;
X_sub_SM__(1+nS_sub,1+nM_sub) = X_sub_wSM___(1+tmp_nw,1+nS_sub,1+nM_sub);
delta_x_sub_SM__(1+nS_sub,1+nM_sub) = delta_x_sub_wSM___(1+tmp_nw,1+nS_sub,1+nM_sub);
delta_y_sub_SM__(1+nS_sub,1+nM_sub) = delta_y_sub_wSM___(1+tmp_nw,1+nS_sub,1+nM_sub);
gamma_z_sub_SM__(1+nS_sub,1+nM_sub) = gamma_z_sub_wSM___(1+tmp_nw,1+nS_sub,1+nM_sub);
I_value_sub_SM__(1+nS_sub,1+nM_sub) = I_value_sub_wSM___(1+tmp_nw,1+nS_sub,1+nM_sub);
end;end;%for nS_sub=0:n_S_sub-1; for nM_sub=0:tmp_n_M-1;
%%%%%%%%;
tmp_t = tic();
parameter.flag_MS_vs_SM = 0;
[ ...
 parameter ...
,euler_polar_a_sub_M_ ...
,euler_azimu_b_sub_M_ ...
,euler_gamma_z_sub_M_ ...
,image_delta_x_sub_M_ ...
,image_delta_y_sub_M_ ...
,image_I_value_sub_M_ ...
,image_X_value_sub_M_ ...
,image_S_index_sub_M_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter ...
,n_w_max ...
,n_S_sub ...
,template_viewing_polar_a_all_(1+index_nS_from_nS_sub_) ...
,template_viewing_azimu_b_all_(1+index_nS_from_nS_sub_) ...
,tmp_n_M ...
,X_sub_SM__ ...
,delta_x_sub_SM__ ...
,delta_y_sub_SM__ ...
,gamma_z_sub_SM__ ...
,I_value_sub_SM__ ...
);
tmp_str = 'SM'; if (parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
%%%%%%%%;
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
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% We use FTK and UX for qbp. ;')); end;
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
quad_from_data_sub_T_CTF_sub_normalized_qk__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
T_sub_wM__ = T_k_p_sub_wkM__(1+index_nw_,:);
CTF_sub_wM_ = CTF_sub_wMn__(:,1+nk_p_r);
CTF2_sub_q_ = CTF2_sub_qk__(:,1+nk_p_r);
T_CTF_sub_wM__ = T_sub_wM__.*reshape(CTF_sub_wM_,[n_w,tmp_n_M]);
quad_from_data_sub_T_CTF_sub_normalized_q_ = (quad_from_data_sub_qwM__ * T_CTF_sub_wM__(:))./max(1e-12,CTF2_sub_q_);
quad_from_data_sub_T_CTF_sub_normalized_qk__(:,1+nk_p_r) = quad_from_data_sub_T_CTF_sub_normalized_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(real(tmp_frac_qk__(:)),real(quad_from_data_sub_T_CTF_sub_normalized_qk__(:)),'ko'); xlabel('tmp_frac_qk__','Interpreter','none'); ylabel('quad_from_data_sub_T_CTF_sub_normalized_qk__','Interpreter','none'); axis equal; axis tight; title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;

%{
[~,tmp_ij] = max(tmp_expR2_p_wSdM____,[],'all','linear'); [ij_nw,ij_nd] = ind2sub([n_w_max,n_S_sub,FTK.n_delta_v],tmp_ij),;
%%%%;
subplot(1,1,1);
hold on;
plot3(data_sub_k_c_wMd__(:,1+0),data_sub_k_c_wMd__(:,1+1),data_sub_k_c_wMd__(:,1+2),'ro');
plot3(data_sub_k_c_wSd__(:,1+0),data_sub_k_c_wSd__(:,1+1),data_sub_k_c_wSd__(:,1+2),'gx');
hold off;
xlim([-1,+1]);
ylim([-1,+1]);
zlim([-1,+1]);
axis equal; axis vis3d;
%%%%;
plot(circshift(data_sub_k_c_wMd__,-(ij_nw-1)),data_sub_k_c_wSd__,'.')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Verdict: looks good so far. ;
% Test again using only one cluster of images and all templates. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_verbose = 1; flag_disp = flag_verbose; nf=0;
n_S_per_Sbatch = 24;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% First setting up templates on sphere. ;')); end;
%%%%%%%%;
[ ...
 data_k_p_polar_a_wS__ ...
,data_k_p_azimu_b_wS__ ...
,data_k_c_0_wS__ ...
,data_k_c_1_wS__ ...
,data_k_c_2_wS__ ...
] = ...
cg_rhs_1( ...
 n_S ...
,n_w ...
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
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% testing sheres-style qbp. ;')); end;
ncluster = floor(n_cluster/2);
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
tmp_n_M = numel(index_nM_from_ncluster_); assert(tmp_n_M==n_index_nM_from_ncluster_(1+ncluster));
%index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster}(1); %<-- pick a single image. ;
%tmp_n_M = numel(index_nM_from_ncluster_); %<-- pick a single image. ;
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: tmp_n_M %d',ncluster,n_cluster,tmp_n_M)); end;
%%%%%%%%;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
if (flag_verbose); disp(sprintf(' %% ncluster %d/%d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_k_p_r)); end;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
%pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
%%%%%%%%;
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
CTF_UX_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S);
CTF_UX_S_k_q_wnS__(1:pm_n_w_sum,1:n_S) = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_.*CTF_k_p_r_xavg_k_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
CTF_UX_S_l2_S_ = sum(abs(CTF_UX_S_k_q_wnS__).^2,1)/n_w_max;
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(k_p_r_,CTF_k_p_r_xavg_k_,'ko-'); xlabel('k');ylabel('CTF'); axis tight; title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%;
X_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
delta_x_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
delta_y_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
gamma_z_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
I_value_sub_wSM___ = zeros(n_w_max,n_S,tmp_n_M);
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,X_sub_wSM___ ...
,delta_x_sub_wSM___ ...
,delta_y_sub_wSM___ ...
,gamma_z_sub_wSM___ ...
,I_value_sub_wSM___ ...
] = ...
ampmh_X_single_cluster_wSM___11( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,reshape(CTF_UX_S_k_q_wnS__,[pm_n_w_max,pm_n_UX_rank,n_S]) ...
,CTF_UX_S_l2_S_ ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
R2_sub_wSM___ = 2-2*X_sub_wSM___;
R2_sub_min_M_ = reshape(min(R2_sub_wSM___,[],[1,2]),[tmp_n_M,1]);
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:tmp_n_M,R2_sub_min_M_,'ko'); xlim([1,tmp_n_M]); ylim([0,4]); xlabel('image index'); ylabel('R2'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_sigma = 0.001;
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now taking small temperature tmp_sigma %0.6f. ;',tmp_sigma)); end;
%%%%%%%%;
numerator_sub_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
denomator_sub_qkM___ = zeros(quad_n_all,n_k_p_r,tmp_n_M);
expR2_sub_sum_M_ = zeros(tmp_n_M,1);
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>0); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nSbatch=0:n_Sbatch-1;
index_nS_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_nS_in_Sbatch_ = index_nS_in_Sbatch_(find(index_nS_in_Sbatch_<n_S)); n_S_sub = numel(index_nS_in_Sbatch_);
if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
%if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
%%%%%%%%;
tmp_X_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
tmp_delta_x_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
tmp_delta_y_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
tmp_gamma_z_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
tmp_I_value_sub_wSM___ = zeros(n_w_max,n_S_sub,tmp_n_M);
%%%%%%%%;
CTF_UX_S_k_q_wnS___ = reshape(CTF_UX_S_k_q_wnS__,[n_w_max,pm_n_UX_rank,n_S]);
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
,reshape(CTF_UX_S_k_q_wnS___(:,:,1+index_nS_in_Sbatch_),[pm_n_w_max,pm_n_UX_rank,n_S_sub]) ...
,CTF_UX_S_l2_S_(1+index_nS_in_Sbatch_) ...
,tmp_n_M ...
,svd_VUXM_lwnM____(:,1:n_w_max,1:pm_n_UX_rank,1+index_nM_from_ncluster_) ...
,UX_M_l2_dM__(:,1+index_nM_from_ncluster_) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_X_single_cluster_wSM___11 sub: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_index_ = reshape(repmat(transpose(0:n_w_max-1),[1,n_S_sub]) + n_w_max*repmat(reshape(index_nS_in_Sbatch_,[1,n_S_sub]),[n_w_max,1]),[n_w_max*n_S_sub,1]);
quad_from_data_sub_qwS__ = quad_from_data_qwS__(:,1+tmp_index_);
%%%%%%%%;
tmp_t = tic();
tmp_R2_wSdM____ = bsxfun(@minus,permute(2 - 2*tmp_X_sub_dwSM____,[2,3,1,4]),reshape(R2_sub_min_M_,[1,1,1,tmp_n_M]));
tmp_expR2_p_wSdM____ = exp(-tmp_R2_wSdM____/(2*tmp_sigma^2));
tmp_p_p_wSdMk_____ = bsxfun(@times,reshape(tmp_expR2_p_wSdM____,[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,1]),reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]))/tmp_sigma^2;
tmp_p_q_wSdMk_____ = reshape(interp_p_to_q_block_0(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,n_w_max*ones(n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,1),n_w_max*n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r,tmp_p_p_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
expR2_sub_sum_M_ = expR2_sub_sum_M_ + reshape(sum(reshape(tmp_expR2_p_wSdM____,[n_w_max*n_S_sub*FTK.n_delta_v,tmp_n_M]),1),[tmp_n_M,1]);
tmp_numerator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(ifft(reshape(bsxfun(@times,reshape(reshape(reshape(permute(reshape(FTK.svd_U_d_expiw_s__*reshape(squeeze(svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+index_nM_from_ncluster_)),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*tmp_n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,tmp_n_M])*n_w_max,[2,1,4,3]),[n_w_max*FTK.n_delta_v*tmp_n_M,pm_n_UX_rank])*transpose(UX_kn__)*diag(X_weight_r_.^(-1)),[n_w_max,FTK.n_delta_v,tmp_n_M,n_k_p_r]),[n_w_max,1,FTK.n_delta_v,tmp_n_M,n_k_p_r]),reshape(conj(tmp_p_q_wSdMk_____),[n_w_max,n_S_sub,FTK.n_delta_v,tmp_n_M,n_k_p_r])),[n_w_max,n_S_sub*FTK.n_delta_v*tmp_n_M*n_k_p_r]),[],1)*n_w_max,[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
tmp_denomator_qdMk____ = reshape(quad_from_data_sub_qwS__*reshape(bsxfun(@times,reshape(CTF_k_p_r_xavg_k_,[1,1,1,1,n_k_p_r]),repmat(sum(tmp_p_p_wSdMk_____,1),[n_w_max,1,1,1,1])),[n_w_max*n_S_sub,FTK.n_delta_v*tmp_n_M*n_k_p_r]),[quad_n_all,FTK.n_delta_v,tmp_n_M,n_k_p_r]);
numerator_sub_qkM___ = numerator_sub_qkM___ + permute(sum(tmp_numerator_qdMk____,2),[1,4,3,2]);
denomator_sub_qkM___ = denomator_sub_qkM___ + permute(sum(tmp_denomator_qdMk____,2),[1,4,3,2]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% numerator_sub_qkM___ denomator_sub_qkM___: %0.6fs',tmp_t)); end
tmp_f = (n_M*n_S)/(tmp_n_M*n_S_sub);
if (flag_verbose); disp(sprintf(' %% full calculation estimated at %0.6fs <-- %0.6fm <-- %0.6fh ',tmp_t*tmp_f,tmp_t*tmp_f/60,tmp_t*tmp_f/3600)); end;
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
% now use numerator and denomator to construct data on the sphere. ;
%%%%%%%%;
size(numerator_sub_qkM___);
size(denomator_sub_qkM___);
size(expR2_sub_sum_M_);
tolerance_denom = 1e-12;
tmp_frac_qk__ = bsxfun(@rdivide,sum(bsxfun(@rdivide,numerator_sub_qkM___,reshape(expR2_sub_sum_M_,[1,1,tmp_n_M])),3),max(sum(bsxfun(@rdivide,denomator_sub_qkM___,reshape(expR2_sub_sum_M_,[1,1,tmp_n_M])),3),tolerance_denom));
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(1:tmp_n_M,expR2_sub_sum_M_,'ko'); xlim([1,tmp_n_M]); ylim([0,4]); xlabel('image index'); ylabel('expR2_sub_sum_M_','Interpreter','none'); title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% Now calculating maximum-likelihood qbp. ;')); end;
X_sub_SM__ = zeros(n_S,tmp_n_M);
delta_x_sub_SM__ = zeros(n_S,tmp_n_M);
delta_y_sub_SM__ = zeros(n_S,tmp_n_M);
gamma_z_sub_SM__ = zeros(n_S,tmp_n_M);
I_value_sub_SM__ = zeros(n_S,tmp_n_M);
for nS=0:n_S-1; for nM_sub=0:tmp_n_M-1;
[~,tmp_ij] = max(X_sub_wSM___(:,1+nS,1+nM_sub)); tmp_nw = tmp_ij-1;
X_sub_SM__(1+nS,1+nM_sub) = X_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
delta_x_sub_SM__(1+nS,1+nM_sub) = delta_x_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
delta_y_sub_SM__(1+nS,1+nM_sub) = delta_y_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
gamma_z_sub_SM__(1+nS,1+nM_sub) = gamma_z_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
I_value_sub_SM__(1+nS,1+nM_sub) = I_value_sub_wSM___(1+tmp_nw,1+nS,1+nM_sub);
end;end;%for nS=0:n_S-1; for nM_sub=0:tmp_n_M-1;
%%%%%%%%;
tmp_t = tic();
parameter.flag_MS_vs_SM = 0; %<-- ensure maximum-likelihood. ;
parameter.f_rand = 0.00; %<-- ensure maximum-likelihood. ;
[ ...
 parameter ...
,euler_polar_a_sub_M_ ...
,euler_azimu_b_sub_M_ ...
,euler_gamma_z_sub_M_ ...
,image_delta_x_sub_M_ ...
,image_delta_y_sub_M_ ...
,image_I_value_sub_M_ ...
,image_X_value_sub_M_ ...
,image_S_index_sub_M_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter ...
,n_w_max ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,tmp_n_M ...
,X_sub_SM__ ...
,delta_x_sub_SM__ ...
,delta_y_sub_SM__ ...
,gamma_z_sub_SM__ ...
,I_value_sub_SM__ ...
);
tmp_str = 'SM'; if (parameter.flag_MS_vs_SM); tmp_str = 'MS'; end;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% %s: update euler_polar_a_M_ euler_azimu_b_M_ euler_gamma_z_M_ : %0.3fs',tmp_str,tmp_t)); end;
%%%%%%%%;
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
if (flag_verbose); disp(sprintf(' %% ;')); end;
if (flag_verbose); disp(sprintf(' %% We use FTK and UX for qbp. ;')); end;
if (flag_verbose); disp(sprintf(' %% Note: normalization only works because we have a single cluster. ;')); end;
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
quad_from_data_sub_T_CTF_sub_normalized_qk__ = zeros(quad_n_all,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
T_sub_wM__ = T_k_p_sub_wkM__(1+index_nw_,:);
CTF_sub_wM_ = CTF_sub_wMn__(:,1+nk_p_r);
CTF2_sub_q_ = CTF2_sub_qk__(:,1+nk_p_r);
T_CTF_sub_wM__ = T_sub_wM__.*reshape(CTF_sub_wM_,[n_w,tmp_n_M]);
quad_from_data_sub_T_CTF_sub_normalized_q_ = (quad_from_data_sub_qwM__ * T_CTF_sub_wM__(:))./max(1e-12,CTF2_sub_q_);
quad_from_data_sub_T_CTF_sub_normalized_qk__(:,1+nk_p_r) = quad_from_data_sub_T_CTF_sub_normalized_q_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
if flag_disp; figure(1+nf);nf=nf+1;clf;figsml; plot(real(tmp_frac_qk__(:)),real(quad_from_data_sub_T_CTF_sub_normalized_qk__(:)),'ko'); xlabel('tmp_frac_qk__','Interpreter','none'); ylabel('quad_from_data_sub_T_CTF_sub_normalized_qk__','Interpreter','none'); axis equal; axis tight; title(sprintf('ncluster %d/%d',ncluster,n_cluster)); end;

