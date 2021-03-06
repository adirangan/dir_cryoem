function ...
[ ...
 parameter ...
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
ampmut_wrap_5( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,index_nCTF_from_nM_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,n_M ...
,M_k_p_wkM__ ...
,euler_polar_a_M_0in_ ...
,euler_azimu_b_M_0in_ ...
,euler_gamma_z_M_0in_ ...
,image_delta_x_acc_M_0in_ ...
,image_delta_y_acc_M_0in_ ...
,a_k_Y_base_yk_ ...
,delta_sigma_base ...
);

str_thisfunction = 'ampmut_wrap_5';
verbose=2;
if (verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); n_w_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_0in_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_0in_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_0in_=[]; end; na=na+1;
if (nargin<1+na); image_delta_x_acc_M_0in_=[]; end; na=na+1;
if (nargin<1+na); image_delta_y_acc_M_0in_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_base_yk_=[]; end; na=na+1;
if (nargin<1+na); delta_sigma_base=[]; end; na=na+1;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'fname_pre')); if (verbose); disp(sprintf(' %% Warning, fname_pre not set in %s',str_thisfunction)); end; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
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
if (~isfield(parameter,'rseed')); parameter.rseed = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'delta_r_max')); parameter.delta_r_max = 0.1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = 32; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_alternate_MS_vs_SM')); parameter.flag_alternate_MS_vs_SM = 1; end; %<-- parameter_bookmark. ;
%%%%%%%%;
tolerance_master = parameter.tolerance_master;
flag_rank_vs_tolerance = parameter.flag_rank_vs_tolerance;
parameter.flag_clump_vs_cluster = parameter.flag_rank_vs_tolerance; flag_clump_vs_cluster = parameter.flag_clump_vs_cluster; %<-- force;
tolerance_cluster = parameter.tolerance_cluster;
tolerance_pm = parameter.tolerance_pm;
rank_pm = parameter.rank_pm;
rank_CTF = parameter.rank_CTF;
rseed = parameter.rseed;
delta_r_max = parameter.delta_r_max;
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
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
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

%%%%%%%%;
% Now run ampmut_5. ;
%%%%%%%%;
a_k_Y_reco_yki__ = []; a_CTF_avg_UX_Y_yni__= [];
FTK = [];
M_k_q_wkM__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];
rng(rseed);
if (~isempty(euler_polar_a_M_0in_)); euler_polar_a_M_ = euler_polar_a_M_0in_; else euler_polar_a_M_= 1*pi*rand(n_M,1); end;
if (~isempty(euler_azimu_b_M_0in_)); euler_azimu_b_M_ = euler_azimu_b_M_0in_; else euler_azimu_b_M_= 2*pi*rand(n_M,1); end;
if (~isempty(euler_gamma_z_M_0in_)); euler_gamma_z_M_ = euler_gamma_z_M_0in_; else euler_gamma_z_M_= 2*pi*rand(n_M,1); end;
if (~isempty(image_delta_x_acc_M_0in_)); image_delta_x_acc_M_ = image_delta_x_acc_M_0in_; else image_delta_x_acc_M_= zeros(n_M,1); end;
if (~isempty(image_delta_y_acc_M_0in_)); image_delta_y_acc_M_ = image_delta_y_acc_M_0in_; else image_delta_y_acc_M_= zeros(n_M,1); end;
image_delta_x_upd_M_ = [];
image_delta_y_upd_M_ = [];
flag_image_delta_upd_M_ = [];
image_I_value_M_ = [];
%%%%;
tmp_t = tic;
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
tmp_t = toc(tmp_t); disp(sprintf(' %% ampmut_5: time %0.2fs',tmp_t));
parameter = parameter_timing_update(parameter,'ampmut_5',tmp_t);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (flag_clump_vs_cluster==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_clump_vs_cluster==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
%%%%%%%%;
if ~isempty(a_k_Y_base_yk_);
%%%%%%%%;
CTF_k_p_r_xcor_kk__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))) * transpose(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1)))) / n_M;
[ ...
 X_2d_xcor_d0__ ...
,X_2d_xcor_d0_weight_r_ ...
] = ...
principled_marching_cost_matrix_6( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_ ...
,1 ...
,[] ...
,a_k_Y_base_yk_ ...
,CTF_k_p_r_xcor_kk__ ...
);
[UX_2d_xcor_d0__,SX_2d_xcor_d0__,VX_2d_xcor_d0__] = svds(X_2d_xcor_d0__,n_UX_rank); SX_2d_xcor_d0_ = diag(SX_2d_xcor_d0__);
if (verbose); disp(sprintf(' %% cumsum(SX_2d_xcor_d0_): ')); disp(sprintf(' %% %0.4f',cumsum(SX_2d_xcor_d0_/sum(SX_2d_xcor_d0_)))); end;
UX__ = UX_2d_xcor_d0__; X_weight_r_ = X_2d_xcor_d0_weight_r_;
%%%%%%%%;
end;%if ~isempty(a_k_Y_base_yk_);
%%%%%%%%;

%%%%%%%%;
if  isempty(a_k_Y_base_yk_);
%%%%%%%%;
% Calculate empirical principal-modes. ;
%%%%%%%%;
[ ...
 X_2d_Memp_d1_kk__ ...
,X_2d_Memp_d1_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
%%%%%%%%;
[UX_2d_Memp_d1__,SX_2d_Memp_d1__,VX_2d_Memp_d1__] = svds(X_2d_Memp_d1_kk__,n_UX_rank); SX_2d_Memp_d1_ = diag(SX_2d_Memp_d1__);
if (verbose); disp(sprintf(' %% cumsum(SX_2d_Memp_d1_): ')); disp(sprintf(' %% %0.4f',cumsum(SX_2d_Memp_d1_/sum(SX_2d_Memp_d1_)))); end;
UX__ = UX_2d_Memp_d1__; X_weight_r_ = X_2d_Memp_d1_weight_r_;
%%%%%%%%;
end;%if  isempty(a_k_Y_base_yk_);
%%%%%%%%;

%%%%%%%%;
% Now run ampmut_4. ;
%%%%%%%%;
a_k_Y_reco_yki__ = []; a_CTF_avg_UX_Y_yni__= [];
FTK = [];
pm_n_UX_rank = rank_pm;
M_k_q_wkM__ = [];
svd_VUXM_lwnM____ = [];
UX_M_l2_dM__ = [];
VSCTF_Mc__ = [];
rng(rseed);
if (~isempty(euler_polar_a_M_0in_)); euler_polar_a_M_ = euler_polar_a_M_0in_; else euler_polar_a_M_= 1*pi*rand(n_M,1); end;
if (~isempty(euler_azimu_b_M_0in_)); euler_azimu_b_M_ = euler_azimu_b_M_0in_; else euler_azimu_b_M_= 2*pi*rand(n_M,1); end;
if (~isempty(euler_gamma_z_M_0in_)); euler_gamma_z_M_ = euler_gamma_z_M_0in_; else euler_gamma_z_M_= 2*pi*rand(n_M,1); end;
if (~isempty(image_delta_x_acc_M_0in_)); image_delta_x_acc_M_ = image_delta_x_acc_M_0in_; else image_delta_x_acc_M_= zeros(n_M,1); end;
if (~isempty(image_delta_y_acc_M_0in_)); image_delta_y_acc_M_ = image_delta_y_acc_M_0in_; else image_delta_y_acc_M_= zeros(n_M,1); end;
image_delta_x_upd_M_ = [];
image_delta_y_upd_M_ = [];
flag_image_delta_upd_M_ = [];
image_I_value_M_ = [];
%%%%;
tmp_t = tic;
[ ...
 parameter ...
,flag_image_delta_upd_M_ ...
,FTK ...
,M_k_q_wkM__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,VSCTF_Mc__ ...
,a_CTF_avg_UX_Y_yni__ ...
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
ampmut_4( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_max ...
,FTK ...
,pm_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p_wkM__ ...
,M_k_q_wkM__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,rank_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,VSCTF_Mc__ ...
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
tmp_t = toc(tmp_t); disp(sprintf(' %% ampmut_4: time %0.2fs',tmp_t));
parameter = parameter_timing_update(parameter,'ampmut_4',tmp_t);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (flag_clump_vs_cluster==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
if ( isfield(parameter,'fname_pre'));
fname_mat = sprintf('%s.mat',parameter.fname_pre);
disp(sprintf(' %% writing %s',fname_mat));
save(fname_mat ...
,'parameter' ...
,'a_k_Y_reco_yki__' ...
,'a_CTF_avg_UX_Y_yni__' ...
,'euler_polar_a_Mi__' ...
,'euler_azimu_b_Mi__' ...
,'euler_gamma_z_Mi__' ...
,'image_delta_x_acc_Mi__' ...
,'image_delta_y_acc_Mi__' ...
,'image_delta_x_upd_Mi__' ...
,'image_delta_y_upd_Mi__' ...
,'flag_image_delta_upd_Mi__' ...
,'image_I_value_Mi__' ...
,'image_X_value_Mi__' ...
,'image_S_index_Mi__' ...
);
end;%if ( isfield(parameter,'fname_pre'));
%%%%%%%%;

if (verbose); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
