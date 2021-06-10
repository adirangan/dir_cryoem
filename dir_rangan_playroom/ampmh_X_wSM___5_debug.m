function ...
[ ...
 n_S ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,X_wSM___ ...
,delta_j_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___5_debug( ...
 FTK ...
,viewing_k_eq_d ...
,n_w_max ...
,l_max_max ...
,pm_n_UX_rank ...
,a_UCTF_UX_Y_0lsq_ync__ ...
,CTF_index_ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,n_S_per_Sbatch ...
,n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
%%%%%%%%;
% calls ampmh_X_wSM___4_debug.m for each micrograph. ;
%%%%%%%%;
verbose = 2;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_wSM___5_debug]')); end;

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

%%%%%%%%;
[n_S] = sample_shell_5(pm_k_p_r_max,viewing_k_eq_d,'L') ; %<-- obtain viewing angles on outer shell. ;
tmp_t = tic();
tmp_verbose=0;
UCTF_UX_S_k_p_wSc___ = zeros(pm_n_w_sum,n_S,n_CTF_rank);
for nCTF_rank=0:n_CTF_rank-1;
[ ...
 UCTF_UX_S_k_p_wSc___(:,:,1+nCTF_rank) ...
,~ ...
,~ ...
,~ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_template_1( ...
 tmp_verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_k_p_r_ ...
,pm_l_max_ ...
,a_UCTF_UX_Y_0lsq_ync__(:,1+nCTF_rank) ...
,viewing_k_eq_d ...
,-1 ...
,pm_n_w_ ...
);
assert(n_S==n_viewing_all);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% group images by micrograph (i.e., sort by CTF_index_). ;
%%%%%%%%;
u_CTF_index_ = unique(CTF_index_(1:n_M)); n_u_CTF_index = numel(u_CTF_index_);
index_M_CTF_index__ = cell(n_u_CTF_index,1);
n_u_CTF_index_ = zeros(n_u_CTF_index,1);
for nu_CTF_index=0:n_u_CTF_index-1;
u_CTF_index = u_CTF_index_(1+nu_CTF_index);
index_M_CTF_index__{1+nu_CTF_index} = efind(CTF_index_(1:n_M)==u_CTF_index);
n_u_CTF_index_(1+nu_CTF_index) = numel(index_M_CTF_index__{1+nu_CTF_index});
end;%for nu_CTF_index=0:n_u_CTF_index-1;
if (verbose); disp(sprintf(' %% n_u_CTF_index %d',n_u_CTF_index)); end;
%%%%%%%%;
if (verbose);
for nu_CTF_index=0:n_u_CTF_index-1;
disp(sprintf(' %% nu_CTF_index %.3d/%.3d (%.8d) <-- n %.3d',nu_CTF_index,n_u_CTF_index,u_CTF_index_(1+nu_CTF_index),n_u_CTF_index_(1+nu_CTF_index)));
end;%for nu_CTF_index=0:n_u_CTF_index-1;
disp(sprintf(' %% n_M %.4d sum(n_u_CTF_index_) = %.4d',n_M,sum(n_u_CTF_index_)));
end;%if (verbose);

X_wSM___ = zeros(n_w_max,n_S,n_M);
delta_j_wSM___ = zeros(n_w_max,n_S,n_M);
I_value_wSM___ = zeros(n_w_max,n_S,n_M);
%%%%%%%%;
% step through each micrograph, ;
% calculate the templates associated with that particular CTF-function, ;
% and then calculate innerproducts between those templates and all the images. ;
%%%%%%%%;
for nu_CTF_index=0:n_u_CTF_index-1;
tmp_index_M_ = index_M_CTF_index__{1+nu_CTF_index};
tmp_n_M = n_u_CTF_index_(1+nu_CTF_index);
if (verbose); disp(sprintf(' %% nu_CTF_index %d/%d --> tmp_n_M %d [%d,..,%d] ',nu_CTF_index,n_u_CTF_index,tmp_n_M,tmp_index_M_(0+1),tmp_index_M_(tmp_n_M-1+1))); end;
%%%%%%%%;
% Find templates. ;
%%%%%%%%;
tmp_t = tic();
VSCTF_avg_ = mean(VSCTF_Mc__(1+tmp_index_M_,:),1);
VSCTF_std_ = std(VSCTF_Mc__(1+tmp_index_M_,:),1,1);
assert(max(VSCTF_std_./max(1e-12,abs(VSCTF_avg_)))<1e-6); %<-- consider lowering threshold to 1e-3. ;
CTF_UX_S_k_p_wS__ = zeros(pm_n_w_sum,n_S);
for nCTF_rank=0:n_CTF_rank-1;
CTF_UX_S_k_p_wS__ = CTF_UX_S_k_p_wS__ + UCTF_UX_S_k_p_wSc___(:,:,1+nCTF_rank) * VSCTF_avg_(1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
%%%%%%%%;
CTF_UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
CTF_UX_S_l2_(1+nS) = ...
innerproduct_p_quad( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_2d_k_p_r_/(2*pi) ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wS__(:,1+nS) ...
,CTF_UX_S_k_p_wS__(:,1+nS) ...
);
end;%for nS=0:n_S-1;
%%%%%%%%;
CTF_UX_S_k_q_wS__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
CTF_UX_S_k_q_wS__(:,1+nS) = ...
interp_p_to_q( ...
 pm_n_k_p_r ...
,pm_n_w_ ...
,pm_n_w_sum ...
,CTF_UX_S_k_p_wS__(:,1+nS) ...
); 
end;%for nS=0:n_S-1; 
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
SS_k_q_ = svd(CTF_UX_S_k_q_wS__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(CTF_UX_S_k_q_wS__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(CTF_UX_S_k_q_wS__,n_S_rank);
if (verbose>1); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% combining UCTF_UX_S_k_p_wSc___ to form templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
save('ampmh_X_wSM___5_debug.mat');
error('exiting early');
disp('returning');return;
%%%%%%%%;
tmp_t = tic();
[ ...
 X_wSM___(:,:,1+tmp_index_M_) ...
,delta_j_wSM___(:,:,1+tmp_index_M_) ...
,I_value_wSM___(:,:,1+tmp_index_M_) ...
] = ...
ampmh_X_wSM___4_debug( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_S_rank ...
,n_S_per_Sbatch ...
,US_k_q__ ...
,SS_k_q__ ...
,VS_k_q__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,n_M_per_Mbatch ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
end;%for nu_CTF_index=0:n_u_CTF_index-1;

if (verbose); disp(sprintf(' %% [finished ampmh_X_wSM___5_debug]')); end;

return;

%%%%%%%%;
% profile ampmh_X_wSM___6, which uses raw templates. ;
%%%%%%%%;
profile on;
parameter = struct('type','parameter');
parameter.svd_eps_use = 0.1;
parameter.n_svd_l_use = 0;
parameter.n_delta_v_use = 1;
parameter.pm_n_UX_rank_use = floor(pm_n_UX_rank/4);
parameter.flag_optimize_over_gamma_z = 1;
tmp_t = tic();
[ ...
 X_wSM___(:,:,1+tmp_index_M_) ...
,delta_x_wSM___(:,:,1+tmp_index_M_) ...
,delta_y_wSM___(:,:,1+tmp_index_M_) ...
,gamma_z_wSM___(:,:,1+tmp_index_M_) ...
] = ...
ampmh_X_wSM___6( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,16 ...
,CTF_UX_S_k_q_wS__ ...
,CTF_UX_S_l2_ ...
,tmp_n_M ...
,8 ...
,svd_VUXM_lwnM____(:,:,:,1+tmp_index_M_) ...
,UX_M_l2_dM__(:,1+tmp_index_M_) ...
,parameter ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
profile viewer;
profile off;
