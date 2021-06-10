function ...
[ ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
] = ...
cg_lsq_pm_reduce_1( ...
 parameter ...
,pm_n_UX_rank ...
,pm_n_w_ ...
,n_M ...
,UX_M_k_p_wnM__ ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
);

verbose = 0;
if (verbose); disp(sprintf(' %% [entering cg_lsq_pm_reduce_1]')); end;

pm_n_w_max = max(pm_n_w_);
flag_unique_pm_n = (numel(unique(pm_n_w_))==1);
if (~flag_unique_pm_n); error(sprintf(' Error: cg_lsq_pm_reduce_1 expects all pm_n_w_ to be the same.')); end;
UX_M_k_p_wnM___ = reshape(UX_M_k_p_wnM__,[pm_n_w_max,pm_n_UX_rank,n_M]); %<-- used later. ;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'pm_n_w_max_use')); parameter.pm_n_w_max_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'pm_n_UX_rank_use')); parameter.pm_n_UX_rank_use = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_CTF_rank_use')); parameter.n_CTF_rank_use = 0; end; %<-- parameter_bookmark. ;
%%%%%%%%;
pm_n_w_max_use = max(0,min(pm_n_w_max,parameter.pm_n_w_max_use));
pm_n_UX_rank_use = max(0,min(pm_n_UX_rank,parameter.pm_n_UX_rank_use));
n_CTF_rank_use = max(0,min(n_CTF_rank,parameter.n_CTF_rank_use));
%%%%%%%%;
flag_pm_n_w_max_subselect = 0;
if (pm_n_w_max_use> 0 & pm_n_w_max_use< pm_n_w_max);
if (verbose); disp(sprintf(' %% reducing pm_n_w_max %d --> %d',pm_n_w_max,pm_n_w_max_use)); end;
flag_pm_n_w_max_subselect = 1;
index_pm_nw_ = 0:pm_n_w_max-1;
index_pm_nw_ = index_pm_nw_(find(abs(periodize(index_pm_nw_,-floor(pm_n_w_max/2),+floor(pm_n_w_max/2)))<=floor(pm_n_w_max_use/2)));
pm_n_w_max_use = numel(index_pm_nw_);
UX_M_k_p_wnM___ = UX_M_k_p_wnM___(1+index_pm_nw_,:,:); %<-- replace input with reduced version. ;
parameter.pm_n_w_max_use = pm_n_w_max_use; %<-- replace input with reduced version. ;
pm_n_w_max = pm_n_w_max_use; %<-- replace input with reduced version. ;
index_pm_nw_ = 0:pm_n_w_max-1;
end;%if (pm_n_w_max_use> 0 & pm_n_w_max_use< pm_n_w_max);
if (pm_n_w_max_use<=0 | pm_n_w_max_use>=pm_n_w_max);
pm_n_w_max_use = pm_n_w_max;
index_pm_nw_ = 0:pm_n_w_max-1;
end;%if (pm_n_w_max_use<=0 | pm_n_w_max_use>=pm_n_w_max);
%%%%%%%%;
flag_pm_n_UX_rank_subselect = 0;
if (pm_n_UX_rank_use> 0 & pm_n_UX_rank_use< pm_n_UX_rank);
if (verbose); disp(sprintf(' %% reducing pm_n_UX_rank %d --> %d',pm_n_UX_rank,pm_n_UX_rank_use)); end;
flag_pm_n_UX_rank_subselect = 1;
pm_n_UX_rank_use = max(0,min(pm_n_UX_rank,pm_n_UX_rank_use));
end;%if (pm_n_UX_rank_use> 0 & pm_n_UX_rank_use< pm_n_UX_rank);
if (pm_n_UX_rank_use<=0 | pm_n_UX_rank_use>=pm_n_UX_rank);
pm_n_UX_rank_use = pm_n_UX_rank;
end;%if (pm_n_UX_rank_use<=0 | pm_n_UX_rank_use>=pm_n_UX_rank);
%%%%;
if (flag_pm_n_UX_rank_subselect);
UX_M_k_p_wnM___ = UX_M_k_p_wnM___(:,1:pm_n_UX_rank_use,:); %<-- replace input with reduced version. ;
pm_n_UX_rank = pm_n_UX_rank_use;
parameter.pm_n_UX_rank_use = pm_n_UX_rank;
end;%if (flag_pm_n_UX_rank_subselect);
%%%%%%%%;
flag_n_CTF_rank_use = 0;
if (n_CTF_rank_use> 0 & n_CTF_rank_use< n_CTF_rank);
flag_n_CTF_rank_use = 1;
VSCTF_Mc__ = VSCTF_Mc__(:,1:n_CTF_rank_use); %<-- replace input with reduced version. ;
n_CTF_rank = n_CTF_rank_use; %<-- replace input with reduced version. ;
end;%if (n_CTF_rank_use> 0 & n_CTF_rank_use< n_CTF_rank);

assert(size(UX_M_k_p_wnM___,1)==pm_n_w_max);
assert(size(UX_M_k_p_wnM___,2)==pm_n_UX_rank);
assert(size(UX_M_k_p_wnM___,3)==n_M);
if (flag_pm_n_w_max_subselect | flag_pm_n_UX_rank_subselect);
UX_M_k_p_wnM__ = reshape(UX_M_k_p_wnM___,[pm_n_w_max*pm_n_UX_rank,n_M]);
end;%if (flag_pm_n_w_max_subselect | flag_pm_n_UX_rank_subselect);

assert(size(VSCTF_Mc__,1)==n_M);
assert(size(VSCTF_Mc__,2)==n_CTF_rank);

if (verbose); disp(sprintf(' %% [finished cg_lsq_pm_reduce_1]')); end;
