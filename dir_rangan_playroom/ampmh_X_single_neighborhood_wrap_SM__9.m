function ...
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
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
%%%%%%%%;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_single_neighborhood_wrap_SM__9]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;

pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%{
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
%}
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);

%%%%%%%%;
% group images by micrograph (i.e., sort by CTF_index_). ;
%%%%%%%%;
u_CTF_index_ = unique(CTF_index_(1:n_M)); n_u_CTF_index = numel(u_CTF_index_);
index_nM_from_nu_CTF_index__ = cell(n_u_CTF_index,1);
n_u_CTF_index_ = zeros(n_u_CTF_index,1);
for nu_CTF_index=0:n_u_CTF_index-1;
u_CTF_index = u_CTF_index_(1+nu_CTF_index);
index_nM_from_nu_CTF_index__{1+nu_CTF_index} = efind(CTF_index_(1:n_M)==u_CTF_index);
n_u_CTF_index_(1+nu_CTF_index) = numel(index_nM_from_nu_CTF_index__{1+nu_CTF_index});
end;%for nu_CTF_index=0:n_u_CTF_index-1;
if (verbose); disp(sprintf(' %% n_u_CTF_index %d',n_u_CTF_index)); end;
%%%%%%%%;
index_nu_CTF_index_from_nM_ = zeros(n_M,1);
for nu_CTF_index=0:n_u_CTF_index-1;
index_nM_from_nu_CTF_index_ = index_nM_from_nu_CTF_index__{1+nu_CTF_index};
index_nu_CTF_index_from_nM_(1+index_nM_from_nu_CTF_index_) = nu_CTF_index;
end;%for nu_CTF_index=0:n_u_CTF_index-1;
%%%%%%%%;

%%%%%%%%;
% Calculate l2-norm of each template. ;
%%%%%%%%;
tmp_t = tic();
tmp_UCTF_UX_S_k_q_cwnSx__ = reshape(permute(UCTF_UX_S_k_q_wncSx____,[2,1,3,4]),[n_CTF_rank,pm_n_w_sum*n_S*n_neighborhood]);
CTF_UX_S_l2_cSx__ = zeros(n_u_CTF_index,n_S*n_neighborhood);
for nu_CTF_index=0:n_u_CTF_index-1;
tmp_index_nM_from_nu_CTF_index_ = index_nM_from_nu_CTF_index__{1+nu_CTF_index};
VSCTF_avg_ = mean(VSCTF_Mc__(1+tmp_index_nM_from_nu_CTF_index_,:),1);
VSCTF_std_ = std(VSCTF_Mc__(1+tmp_index_nM_from_nu_CTF_index_,:),1,1);
if (n_CTF_rank>1); assert(max(VSCTF_std_./max(1e-12,abs(VSCTF_avg_)))<1e-6); end; %<-- consider lowering threshold to 1e-3. ;
tmp_CTF_UX_S_k_q_wnSx__ = zeros(pm_n_w_sum,n_S*n_neighborhood);
tmp_CTF_UX_S_k_q_wnSx__ = reshape(VSCTF_avg_*tmp_UCTF_UX_S_k_q_cwnSx__,[pm_n_w_sum,n_S*n_neighborhood]);
tmp_CTF_UX_S_l2_Sx_ = zeros(n_S*n_neighborhood,1);
tmp_CTF_UX_S_l2_Sx_ = sum(abs(tmp_CTF_UX_S_k_q_wnSx__).^2,1)/pm_n_w_max;
CTF_UX_S_l2_cSx__(1+nu_CTF_index,:) = tmp_CTF_UX_S_l2_Sx_;
end;%for nu_CTF_index=0:n_u_CTF_index-1;
clear tmp_UCTF_UX_S_k_q_cwnSx__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% CTF_UX_S_l2_cSx__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'CTF_UX_S_l2_cSx__',tmp_t);
%%%%%%%%;
% Now compare against a few direct calculations. ;
%%%%%%%%;
flag_check=0;
if (flag_check);
n_test=9;
for nu_CTF_index=0:n_u_CTF_index-1;
tmp_index_nM_from_nu_CTF_index_ = index_nM_from_nu_CTF_index__{1+nu_CTF_index};
VSCTF_avg_ = mean(VSCTF_Mc__(1+tmp_index_nM_from_nu_CTF_index_,:),1);
VSCTF_std_ = std(VSCTF_Mc__(1+tmp_index_nM_from_nu_CTF_index_,:),1,1);
if (n_CTF_rank>1); assert(max(VSCTF_std_./max(1e-12,abs(VSCTF_avg_)))<1e-6); end; %<-- consider lowering threshold to 1e-3. ;
disp(sprintf(' %% nu_CTF_index %d/%d tmp_n_M %d/%d',nu_CTF_index,n_u_CTF_index,numel(tmp_index_nM_from_nu_CTF_index_),n_M));
for ntest=0:n_test-1;
nl = max(0,min(n_S*n_neighborhood-1,floor(n_S*n_neighborhood*ntest/(n_test-1))));
[nS,nneighborhood] = ind2sub([n_S,n_neighborhood],1+nl); nS=nS-1; nneighborhood=nneighborhood-1;
disp(sprintf(' %% ntest %d/%d: nl %d nS %d nneighborhood %d',ntest,n_test,nl,nS,nneighborhood));
tmp_CTF_UX_S_k_q_wn_ = zeros(pm_n_w_sum,1);
for nCTF_rank=0:n_CTF_rank-1;
tmp_CTF_UX_S_k_q_wn_ = tmp_CTF_UX_S_k_q_wn_ + UCTF_UX_S_k_q_wncSx____(:,1+nCTF_rank,1+nS,1+nneighborhood) * VSCTF_avg_(1+nCTF_rank);
end;%for nCTF_rank=0:n_CTF_rank-1;
tmp_CTF_UX_S_l2 = ...
innerproduct_p_quad( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_weight_2d_k_p_r_/(2*pi) ...
,pm_n_w_ ...
,pm_n_w_sum ...
,tmp_CTF_UX_S_k_q_wn_ ...
,tmp_CTF_UX_S_k_q_wn_ ...
);
disp(sprintf(' %% tmp_CTF_UX_S_l2 vs CTF_UX_S_l2_cSx__(1+nu_CTF_index,1+nl): %0.16f',fnorm(tmp_CTF_UX_S_l2 - CTF_UX_S_l2_cSx__(1+nu_CTF_index,1+nl))/fnorm(tmp_CTF_UX_S_l2)));
end;%for ntest=0:n_test-1;
end;%for nu_CTF_index=0:n_u_CTF_index-1;
end;%if (flag_check);
%%%%%%%%;
CTF_UX_S_l2_cSx___ = reshape(CTF_UX_S_l2_cSx__,[n_u_CTF_index,n_S,n_neighborhood]);

X_SM__ = zeros(n_S,n_M);
delta_x_SM__ = zeros(n_S,n_M);
delta_y_SM__ = zeros(n_S,n_M);
gamma_z_SM__ = zeros(n_S,n_M);
I_value_SM__ = zeros(n_S,n_M);

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
,I_value_SM__(:,1+index_neighborhood_ori_M_) ...
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
flag_check=0;
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
,tmp_I_value_Sm_ ...
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
disp(sprintf(' %% tmp_I_value_Sm_ vs I_value_SM__(:,1+tmp_nM): %0.16f',fnorm(tmp_I_value_Sm_ - I_value_SM__(:,1+tmp_nM))/fnorm(tmp_I_value_Sm_)));
end;%for ntest=0:n_test-1;
end;% if (flag_check);
%%%%%%%%;
clear UCTF_UX_S_k_q_wnSc___ CTF_UX_S_l2_SM__;
%%%%%%%%;
end;%if (tmp_n_M>0);
end;%for nneighborhood=0:n_neighborhood-1;

if ( (nargout>5) & (isempty(I_value_SM__)) ); I_value_SM__ = ones(n_S,n_M); end;

if (verbose>0); disp(sprintf(' %% [finished ampmh_X_single_neighborhood_wrap_SM__9]')); end;
