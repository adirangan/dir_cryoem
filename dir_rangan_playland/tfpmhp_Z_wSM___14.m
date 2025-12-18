function ...
[ ...
 parameter ...
,M_k_q_wkM__ ...
,UX_T_M_l2_dM__ ...
,UX_M_l2_M_ ...
,svd_V_UX_M_lwnM____ ... %<-- or UX_T_M_k_q_dwnM____ ;
,UX_CTF_S_k_q_wnS__ ...
,UX_CTF_S_l2_S_ ...
] = ...
tfpmhp_Z_wSM___14( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,n_M ...
,M_k_p_wkM__ ...
,CTF_k_p_r_k_ ...
,n_S ...
,S_k_p_wkS__ ...
,pm_n_UX_rank ...
,pm_UX_kn__ ...
,pm_X_weight_r_ ...
,FTK ...
,index_nM_ ...
,M_k_q_wkM__ ...
,UX_T_M_l2_dM__ ...
,UX_M_l2_M_ ...
,svd_V_UX_M_lwnM____ ... %<-- or UX_T_M_k_q_dwnM____ ;
,index_nS_ ...
,UX_CTF_S_k_q_wnS__ ...
,UX_CTF_S_l2_S_ ...
);

%%%%%%%%;
% precomputation for tfpmh_Z_wSM___14. ;
% Note that the radial compression is applied here, and is structured to align with the canonical innerproduct. ;
% (i.e., can use radial compression even when the loss is high). ;
% Note that this takes in a single (isotropic) CTF_k_p_r_k_, ;
% Assumes n_w_ = n_w_max*ones(n_k_p_r,1);
%%%%%%%%;

str_thisfunction = 'tfpmhp_Z_wSM___14';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_wk_=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_k_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); pm_n_UX_rank=[]; end; na=na+1;
if (nargin<1+na); pm_UX_kn__=[]; end; na=na+1;
if (nargin<1+na); pm_X_weight_r_=[]; end; na=na+1;
if (nargin<1+na); FTK=[]; end; na=na+1;
if (nargin<1+na); index_nM_=[]; end; na=na+1;
if (nargin<1+na); M_k_q_wkM__=[]; end; na=na+1;
if (nargin<1+na); UX_T_M_l2_dM__=[]; end; na=na+1;
if (nargin<1+na); UX_M_l2_M_=[]; end; na=na+1;
if (nargin<1+na); svd_V_UX_M_lwnM____=[]; end; na=na+1;
if (nargin<1+na); index_nS_=[]; end; na=na+1;
if (nargin<1+na); UX_CTF_S_k_q_wnS__=[]; end; na=na+1;
if (nargin<1+na); UX_CTF_S_l2_S_=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'memory_limit_GB')); parameter.memory_limit_GB = 4.0; end; %<-- parameter_bookmark. ;
memory_limit_GB = parameter.memory_limit_GB;
%%%%%%%%;
if ~isfield(parameter,'flag_precompute_M_k_q_wkM__'); parameter.flag_precompute_M_k_q_wkM__ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_M_k_q_wkM__ = parameter.flag_precompute_M_k_q_wkM__;
if ~isfield(parameter,'flag_precompute_UX_T_M_l2_dM__'); parameter.flag_precompute_UX_T_M_l2_dM__ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_UX_T_M_l2_dM__ = parameter.flag_precompute_UX_T_M_l2_dM__;
if ~isfield(parameter,'flag_precompute_UX_M_l2_M_'); parameter.flag_precompute_UX_M_l2_M_ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_UX_M_l2_M_ = parameter.flag_precompute_UX_M_l2_M_;
if ~isfield(parameter,'flag_precompute_svd_V_UX_M_lwnM____'); parameter.flag_precompute_svd_V_UX_M_lwnM____ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_svd_V_UX_M_lwnM____ = parameter.flag_precompute_svd_V_UX_M_lwnM____;
if ~isfield(parameter,'flag_precompute_UX_CTF_S_k_q_wnS__'); parameter.flag_precompute_UX_CTF_S_k_q_wnS__ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_UX_CTF_S_k_q_wnS__ = parameter.flag_precompute_UX_CTF_S_k_q_wnS__;
if ~isfield(parameter,'flag_precompute_UX_CTF_S_l2_S_'); parameter.flag_precompute_UX_CTF_S_l2_S_ = 0; end; %<-- parameter_bookmark. ;
flag_precompute_UX_CTF_S_l2_S_ = parameter.flag_precompute_UX_CTF_S_l2_S_;
%%%%%%%%;

nf=0;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

n_byte_per_float32 = 4; n_byte_per_float64 = 8;
n_byte_per_complex64 = 8; n_byte_per_complex128 = 16;

n_w_ = n_w_(:);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (n_w_sum~=n_w_max*n_k_p_r); disp(sprintf(' %% Warning, n_w_sum %d ~= n_w_max*n_k_p_r %d*%d in %s',n_w_sum,n_w_max,n_k_p_r,str_thisfunction)); end;
if mod(n_w_max,2)~=0; disp(sprintf(' %% Warning, n_w_max %d in %s',n_w_max,str_thisfunction)); end;
n_delta_v = FTK.n_delta_v;
n_svd_l = FTK.n_svd_l;
tmp_index_d0 = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); assert(numel(tmp_index_d0)>=1); tmp_index_d0=tmp_index_d0(1+0); %<-- should be a single index corresponding to zero-displacement. ;
pm_n_k_p_r = pm_n_UX_rank; pm_n_w_max = n_w_max;
pm_n_w_ = pm_n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_sum = pm_n_k_p_r*pm_n_w_max;
pm_wUX_kn__ = diag(pm_X_weight_r_)*pm_UX_kn__;

%%%%%%%%;
% allocate memory for output. ;
%%%%%%%%;
%%%%%%%%;
if ~isempty(index_nM_);
if flag_precompute_M_k_q_wkM__==1; if isempty(M_k_q_wkM__); M_k_q_wkM__ = zeros(n_w_sum,n_M); end; end;
if flag_precompute_UX_T_M_l2_dM__==1; if isempty(UX_T_M_l2_dM__); UX_T_M_l2_dM__ = zeros(n_delta_v,n_M); end; end;
if flag_precompute_UX_M_l2_M_==1; if isempty(UX_M_l2_M_); UX_M_l2_M_ = zeros(n_M,1); end; end;
if flag_precompute_svd_V_UX_M_lwnM____==1;
if isempty(svd_V_UX_M_lwnM____);
if FTK.flag_tf_vs_bf==1; svd_V_UX_M_lwnM____ = zeros(n_svd_l,n_w_max,pm_n_k_p_r,n_M); end;
if FTK.flag_tf_vs_bf==0; svd_V_UX_M_lwnM____ = zeros(n_delta_v,n_w_max,pm_n_k_p_r,n_M); end; %<-- UX_T_M_k_q_dwnM____. ;
end;%if isempty(svd_V_UX_M_lwnM____);
end;%if flag_precompute_svd_V_UX_M_lwnM____==1;
end;%if ~isempty(index_nM_);
%%%%%%%%;
if ~isempty(index_nS_);
if flag_precompute_UX_CTF_S_k_q_wnS__==1; if isempty(UX_CTF_S_k_q_wnS__); UX_CTF_S_k_q_wnS__ = zeros(pm_n_w_sum,n_S); end; end;
if flag_precompute_UX_CTF_S_l2_S_==1; if isempty(UX_CTF_S_l2_S_); UX_CTF_S_l2_S_ = zeros(n_S,1); end; end;
end;%if ~isempty(index_nS_);
%%%%%%%%;

if (flag_verbose>1); 
if flag_precompute_M_k_q_wkM__==1; tmp_str = 'M_k_q_wkM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if flag_precompute_UX_T_M_l2_dM__==1; tmp_str = 'UX_T_M_l2_dM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if flag_precompute_UX_M_l2_M_==1; tmp_str = 'UX_M_l2_M_'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if flag_precompute_svd_V_UX_M_lwnM____==1; tmp_str = 'svd_V_UX_M_lwnM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if flag_precompute_UX_CTF_S_k_q_wnS__==1; tmp_str = 'UX_CTF_S_k_q_wnS__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if flag_precompute_UX_CTF_S_l2_S_==1; tmp_str = 'UX_CTF_S_l2_S_'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (flag_verbose>1); 

%%%%%%%%%%%%%%%%;
if ~isempty(index_nM_);
%%%%%%%%%%%%%%%%;
n_M_sub = numel(index_nM_);
%%%%;
% Prepare M_k_q_wnM__. ;
%%%%;
if flag_precompute_M_k_q_wkM__==1;
tmp_t = tic();
M_k_q_wkM__(:,1+index_nM_) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+index_nM_));
%[~,UX_M_k_p_wnM__(:,1+index_nM_)] = tfpmhh_pm_wUX_0([],n_k_p_r,pm_n_k_p_r,pm_wUX_kn__,n_w_max,n_M_sub,M_k_p_wkM__(:,1+index_nM_));
%UX_M_k_q_wnM__(:,1+index_nM_) = interp_p_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_M_k_p_wnM__(:,1+index_nM_));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% M_k_q_wnM__: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute M_k_q_wnM__',str_thisfunction),tmp_t);
end;%if flag_precompute_M_k_q_wkM__==1;
%%%%;
% Prepare quasi-images. ;
%%%%;
if flag_precompute_svd_V_UX_M_lwnM____==1;
if FTK.flag_tf_vs_bf==1;
tmp_t = tic();
svd_V_UX_M_lwnM____(:,:,:,1+index_nM_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q_wkM__(:,1+index_nM_),pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_V_UX_M_lwnM____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute svd_V_UX_M_lwnM____',str_thisfunction),tmp_t);
end;%if FTK.flag_tf_vs_bf==1;
if FTK.flag_tf_vs_bf==0;
tmp_t = tic();
svd_V_UX_M_lwnM____(:,:,:,1+index_nM_) = tpmh_UXTM_dwnM____0(FTK,n_k_p_r,k_p_r_,n_w_,n_M_sub,M_k_p_wkM__(:,1+index_nM_),pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_T_M_k_q_dwnM____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute UX_T_M_k_q_dwnM____',str_thisfunction),tmp_t);
end;%if FTK.flag_tf_vs_bf==0;
end;%if flag_precompute_svd_V_UX_M_lwnM____==1;
%%%%;
% Now calculate norms of the translated images. ;
%%%%;
if flag_precompute_UX_T_M_l2_dM__==1;
if FTK.flag_tf_vs_bf==1;
tmp_t = tic();
UX_T_M_l2_dM__(:,1+index_nM_) = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_V_UX_M_lwnM____(:,:,:,1+index_nM_));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__1: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute tfpmh_UX_T_M_l2_dM__1',str_thisfunction),tmp_t);
end;%if FTK.flag_tf_vs_bf==1;
if FTK.flag_tf_vs_bf==0;
tmp_t = tic();
UX_T_M_l2_dM__(:,1+index_nM_) = tfpmh_UX_T_M_l2_dM__0(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_V_UX_M_lwnM____(:,:,:,1+index_nM_));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% tfpmh_UX_T_M_l2_dM__0: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute tfpmh_UX_T_M_l2_dM__0',str_thisfunction),tmp_t);
end;%if FTK.flag_tf_vs_bf==0;
end;%if flag_precompute_UX_T_M_l2_dM__==1;
if flag_precompute_UX_M_l2_M_==1;
UX_M_l2_M_(1+index_nM_) = reshape(UX_T_M_l2_dM__(1+tmp_index_d0,1+index_nM_),[n_M_sub,1]);
end;%if flag_precompute_UX_M_l2_M_==1;
%%%%%%%%%%%%%%%%;
end;%if ~isempty(index_nM_);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if ~isempty(index_nS_);
%%%%%%%%%%%%%%%%;
n_S_sub = numel(index_nS_);
%%%%;
% Prepare UX_CTF_S_k_q_wnS__. ;
%%%%;
if flag_precompute_UX_CTF_S_k_q_wnS__==1;
tmp_t = tic();
CTF_S_sub_k_p_wkS__ = reshape(bsxfun(@times,reshape(CTF_k_p_r_k_,[1,n_k_p_r]),reshape(S_k_p_wkS__(:,1+index_nS_),[n_w_max,n_k_p_r,n_S_sub])),[n_w_sum,n_S_sub]);
[~,UX_CTF_S_sub_k_p_wnS__] = tfpmhh_pm_wUX_0([],n_k_p_r,pm_n_k_p_r,pm_wUX_kn__,n_w_max,n_S_sub,CTF_S_sub_k_p_wkS__);
UX_CTF_S_k_q_wnS__(:,1+index_nS_) = interp_p_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_CTF_S_sub_k_p_wnS__);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_CTF_S_k_q_wnS__: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute UX_CTF_S_k_q_wnS__',str_thisfunction),tmp_t);
end;%if flag_precompute_UX_CTF_S_k_q_wnS__==1;
%%%%;
% Prepare UX_CTF_S_l2_S_. ;
%%%%;
if flag_precompute_UX_CTF_S_l2_S_==1;
tmp_t = tic();
UX_CTF_S_l2_S_(1+index_nS_) = squeeze(sum(bsxfun(@times,reshape(conj(UX_CTF_S_k_q_wnS__(:,1+index_nS_)),[pm_n_w_sum,n_S_sub]),reshape(UX_CTF_S_k_q_wnS__(:,1+index_nS_),[pm_n_w_sum,n_S_sub])),1+0)) / max(1,pm_n_w_max);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_CTF_S_l2_S_: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: precompute UX_CTF_S_l2_S_',str_thisfunction),tmp_t);
end;%if flag_precompute_UX_CTF_S_l2_S_==1;
%%%%%%%%%%%%%%%%;
end;%if ~isempty(index_nS_);
%%%%%%%%%%%%%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

