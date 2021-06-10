function ...
[ ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_use_ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM_use__ ...
,parameter ...
] = ...
ampmh_X_wSM_reduce_0( ...
 FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,parameter ...
);

verbose = 0;
if (verbose); disp(sprintf(' %% [entering ampmh_X_wSM_reduce_0]')); end;
CTF_UX_S_k_q_wnS___ = reshape(CTF_UX_S_k_q_wnS__(:,1:n_S),[n_w_max,pm_n_UX_rank,n_S]); %<-- used later. ;

if nargin<10;
parameter = [];
end;%if nargin<13;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);

%%%%%%%%;
if (~isfield(parameter,'n_w_max_use')); parameter.n_w_max_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'svd_eps_use')); parameter.svd_eps_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_svd_l_use')); parameter.n_svd_l_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_delta_v_use')); parameter.n_delta_v_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'pm_n_UX_rank_use')); parameter.pm_n_UX_rank_use = 0; end; %<-- parameter_bookmark. ;
%%%%%%%%;
n_w_max_use = max(0,min(n_w_max,parameter.n_w_max_use));
svd_eps_use = max(0,parameter.svd_eps_use);
n_svd_l_use = max(0,min(FTK.n_svd_l,parameter.n_svd_l_use));
n_delta_v_use = max(0,min(FTK.n_delta_v,parameter.n_delta_v_use));
pm_n_UX_rank_use = max(0,min(pm_n_UX_rank,parameter.pm_n_UX_rank_use));
%%%%%%%%;
flag_n_w_max_subselect = 0;
if (n_w_max_use> 0 & n_w_max_use< n_w_max);
if (verbose); disp(sprintf(' %% reducing n_w_max %d --> %d',n_w_max,n_w_max_use)); end;
flag_n_w_max_subselect = 1;
index_nw_ = 0:n_w_max-1;
index_nw_ = index_nw_(find(abs(periodize(index_nw_,-floor(n_w_max/2),+floor(n_w_max/2)))<=floor(n_w_max_use/2)));
n_w_max_use = numel(index_nw_);
svd_VUXM_lwnM____ = svd_VUXM_lwnM____(:,1+index_nw_,:,:); %<-- replace input with reduced version. ;
CTF_UX_S_k_q_wnS___ = CTF_UX_S_k_q_wnS___(1+index_nw_,:,:);
parameter.n_w_max_use = n_w_max_use; %<-- replace input with reduced version. ;
n_w_max = n_w_max_use; %<-- replace input with reduced version. ;
index_nw_ = 0:n_w_max-1;
end;%if (n_w_max_use> 0 & n_w_max_use< n_w_max);
if (n_w_max_use<=0 | n_w_max_use>=n_w_max);
n_w_max_use = n_w_max;
index_nw_ = 0:n_w_max-1;
end;%if (n_w_max_use<=0 | n_w_max_use>=n_w_max);
%%%%%%%%;
flag_svd_s_subselect = 0;
FTK_n_svd_l = FTK.n_svd_l;
FTK_index_svd_s_sort_from_all_ = 0:FTK_n_svd_l-1;
if (svd_eps_use>0 | n_svd_l_use>0);
if (verbose); disp(sprintf(' %% reducing svd_eps_use to %0.3f, n_svd_l_use to %d',svd_eps_use,n_svd_l_use)); end;
flag_svd_s_subselect = 1;
if (n_svd_l_use<=0); n_svd_l_use = FTK.n_svd_l; end;
[svd_s_sort_,index_svd_s_sort_from_all_] = sort(FTK.svd_s_,'descend'); index_svd_s_sort_from_all_ = index_svd_s_sort_from_all_ - 1;
FTK_n_svd_l = min(n_svd_l_use,max(find(svd_s_sort_>=svd_eps_use)));
FTK_index_svd_s_sort_from_all_ = index_svd_s_sort_from_all_(1:FTK_n_svd_l);
if (verbose); disp(sprintf(' %% FTK_n_svd_l = %d',FTK_n_svd_l)); for nl=0:FTK_n_svd_l-1; disp(sprintf(' %% svd_s_ %0.4f',FTK.svd_s_(1+FTK_index_svd_s_sort_from_all_(1+nl)))); end; end;
end;%if (svd_eps_use>0 | n_svd_l_use>0);
%%%%;
flag_n_delta_v_subselect = 0;
FTK_n_delta_v = FTK.n_delta_v;
FTK_delta_x_ = FTK.delta_x_;
FTK_delta_y_ = FTK.delta_y_;
FTK_svd_U_d_expiw_s__ = FTK.svd_U_d_expiw_s__;
if (n_delta_v_use>0);
if (verbose); disp(sprintf(' %% reducing n_delta_v %d --> %d',n_delta_v,n_delta_v_use)); end;
flag_n_delta_v_subselect = 1;
[ ...
 FTK_n_delta_v ...
,FTK_delta_x_ ...
,FTK_delta_y_ ...
,~ ...
,~ ...
,FTK_svd_U_d_expiw_s__ ...
] = ...
ampmh_FTK_subselect_1( ...
 FTK ...
,n_delta_v_use ...
);
end;%if (n_delta_v_use>0);
%%%%;
if (flag_svd_s_subselect | flag_n_delta_v_subselect);
FTK_use = FTK;
FTK_use.n_svd_l = FTK_n_svd_l;
FTK_use.n_delta_v = FTK_n_delta_v;
FTK_use.delta_x_ = FTK_delta_x_;
FTK_use.delta_y_ = FTK_delta_y_;
FTK_use.svd_s_ = FTK.svd_s_(1+FTK_index_svd_s_sort_from_all_); %<-- replace input with reduced version. ;
FTK_use.svd_U_d_expiw_s__ = FTK_svd_U_d_expiw_s__(:,1+FTK_index_svd_s_sort_from_all_); %<-- replace input with reduced version. ;
svd_VUXM_lwnM____ = svd_VUXM_lwnM____(1+FTK_index_svd_s_sort_from_all_,:,:,:); %<-- replace input with reduced version. ;
FTK = FTK_use; %<-- replace input with reduced version. ;
svd_eps_use = 0; n_svd_l_use = 0; n_delta_v_use = 0; %<-- replace input with reduced version. ;
parameter.svd_eps_use = svd_eps_use; %<-- replace input with reduced version. ;
parameter.n_svd_l_use = n_svd_l_use; %<-- replace input with reduced version. ;
parameter.n_delta_v_use = n_delta_v_use; %<-- replace input with reduced version. ;
end;%if (flag_svd_s_subselect | flag_n_delta_v_subselect);
%%%%%%%%;
flag_pm_n_UX_rank_subselect = 0;
if (pm_n_UX_rank_use> 0 & pm_n_UX_rank_use< pm_n_UX_rank);
if (verbose); disp(sprintf(' %% reducing pm_n_UX_rank %d --> %d',pm_n_UX_rank,pm_n_UX_rank_use)); end;
flag_pm_n_UX_rank_subselect = 1;
pm_n_UX_rank_use = min(pm_n_UX_rank,pm_n_UX_rank_use);
end;%if (pm_n_UX_rank_use> 0 & pm_n_UX_rank_use< pm_n_UX_rank);
if (pm_n_UX_rank_use<=0 | pm_n_UX_rank_use>=pm_n_UX_rank);
pm_n_UX_rank_use = pm_n_UX_rank;
end;%if (pm_n_UX_rank_use<=0 | pm_n_UX_rank_use>=pm_n_UX_rank);
%%%%;
if (flag_pm_n_UX_rank_subselect);
svd_VUXM_lwnM____ = svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank_use,:); %<-- replace input with reduced version. ;
CTF_UX_S_k_q_wnS___ = CTF_UX_S_k_q_wnS___(:,1:pm_n_UX_rank_use,:);
pm_n_UX_rank = pm_n_UX_rank_use;
parameter.pm_n_UX_rank_use = pm_n_UX_rank;
end;%if (flag_pm_n_UX_rank_subselect);

if (parameter.n_w_max_use>0); assert(parameter.n_w_max_use==n_w_max); end;
if (parameter.n_svd_l_use>0); assert(parameter.n_svd_l_use==FTK.n_svd_l); end;
if (parameter.n_delta_v_use>0); assert(parameter.n_delta_v_use==FTK.n_delta_v); end;
if (parameter.pm_n_UX_rank_use>0); assert(parameter.pm_n_UX_rank_use==pm_n_UX_rank); end;
assert(size(svd_VUXM_lwnM____,1)==FTK.n_svd_l);
assert(size(svd_VUXM_lwnM____,2)==n_w_max);
assert(size(svd_VUXM_lwnM____,3)==pm_n_UX_rank);
assert(size(svd_VUXM_lwnM____,4)<=n_M);
assert(size(CTF_UX_S_k_q_wnS___,1)==n_w_max);
assert(size(CTF_UX_S_k_q_wnS___,2)==pm_n_UX_rank);
assert(size(CTF_UX_S_k_q_wnS___,3)<=n_S);

%%%%%%%%;
UX_M_l2_dM_use__ = UX_M_l2_dM__;
if (flag_pm_n_UX_rank_subselect | flag_n_delta_v_subselect | flag_svd_s_subselect | flag_n_w_max_subselect);
UX_M_l2_dM_use__ = zeros(FTK.n_delta_v,n_M);
tmp_UX_M_k_q_dwnM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_VUXM_lwnM____(:,:,:,:),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M]);
UX_M_l2_dM_use__(:,:) = sum(reshape(permute(abs(tmp_UX_M_k_q_dwnM____).^2,[2,3,1,4]),[n_w_max*pm_n_UX_rank,FTK.n_delta_v,n_M]),1) * n_w_max;
end;%if (flag_pm_n_UX_rank_subselect | flag_n_delta_v_subselect | flag_svd_s_subselect | flag_n_w_max_subselect);
assert(size(UX_M_l2_dM_use__,1)==FTK.n_delta_v);
assert(size(UX_M_l2_dM_use__,2)<=n_M);
%%%%%%%%;
CTF_UX_S_l2_use_ = CTF_UX_S_l2_;
if (flag_pm_n_UX_rank_subselect | flag_n_w_max_subselect);
for nS=0:n_S-1;
CTF_UX_S_l2_use_(1+nS) = ...
innerproduct_p_quad( ...
 pm_n_UX_rank_use ...
,ones(pm_n_UX_rank,1) ...
,ones(pm_n_UX_rank,1)/(2*pi) ...
,n_w_max*ones(pm_n_UX_rank,1) ...
,n_w_max*pm_n_UX_rank ...
,reshape(CTF_UX_S_k_q_wnS___(:,:,1+nS),[n_w_max*pm_n_UX_rank,1]) ...
,reshape(CTF_UX_S_k_q_wnS___(:,:,1+nS),[n_w_max*pm_n_UX_rank,1]) ...
);
end;%for nS=0:n_S-1;
CTF_UX_S_k_q_wnS__ = reshape(CTF_UX_S_k_q_wnS___,[n_w_max*pm_n_UX_rank,n_S]);
end;%if (flag_pm_n_UX_rank_subselect | flag_n_w_max_subselect);
%%%%%%%%;
assert(size(CTF_UX_S_k_q_wnS__,1)==n_w_max*pm_n_UX_rank);
assert(size(CTF_UX_S_k_q_wnS__,2)<=n_S);

if (verbose); disp(sprintf(' %% [finished ampmh_X_wSM_reduce_0]')); end;
