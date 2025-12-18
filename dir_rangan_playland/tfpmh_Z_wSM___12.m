function ...
[ ...
 parameter ...
,Z_wSM___ ...
,UX_R_CTF_R_S_l2_wS__ ...
,CTF_R_S_l2_wS__ ...
,UX_T_M_l2_dM__ ...
,UX_M_l2_M_ ...
,X_wSM___ ...
,delta_x_wSM___ ...
,delta_y_wSM___ ...
,gamma_z_wSM___ ...
,index_sub_wSM___ ...
,Z_dwSM____ ...
,X_dwSM____ ...
] = ...
tfpmh_Z_wSM___12( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,n_M ...
,M_k_p_wkM__ ...
,CTF_k_p_wk_ ...
,n_S ...
,S_k_p_wkS__ ...
,pm_n_UX_rank ...
,pm_UX_kn__ ...
,pm_X_weight_r_ ...
,FTK ...
);
%%%%%%%%;
% Note that this takes in a single (anisotropic) CTF_k_p_wk_, ;
% assumed to be applied to all the n_M images. ;
%%%%%%%%;

str_thisfunction = 'tfpmh_Z_wSM___12';

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wk_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); pm_n_UX_rank=[]; end; na=na+1;
if (nargin<1+na); pm_UX_kn__=[]; end; na=na+1;
if (nargin<1+na); pm_X_weight_r_=[]; end; na=na+1;
if (nargin<1+na); FTK=[]; end; na=na+1;

if isempty(parameter); parameter = struct('type','parameter'); end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'n_M_per_Mbatch')); parameter.n_M_per_Mbatch = 24; end; %<-- parameter_bookmark. ;
n_M_per_Mbatch = parameter.n_M_per_Mbatch;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 24; end; %<-- parameter_bookmark. ;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
if (~isfield(parameter,'flag_optimize_over_gamma_z')); parameter.flag_optimize_over_gamma_z = 0; end; %<-- parameter_bookmark. ;
flag_optimize_over_gamma_z = parameter.flag_optimize_over_gamma_z;
if (~isfield(parameter,'flag_dwSM')); parameter.flag_dwSM = 0; end; %<-- parameter_bookmark. ;
flag_dwSM = parameter.flag_dwSM;
if (~isfield(parameter,'memory_limit_GB')); parameter.memory_limit_GB = 4.0; end; %<-- parameter_bookmark. ;
memory_limit_GB = parameter.memory_limit_GB;
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
n_delta_v = FTK.n_delta_v;
n_svd_l = FTK.n_svd_l;

%%%%%%%%;
% allocate memory for output. ;
%%%%%%%%;
if flag_optimize_over_gamma_z==0;
Z_wSM___ = zeros(n_w_max,n_S,n_M);
UX_R_CTF_R_S_l2_wS__ = zeros(n_w_max,n_S);
CTF_R_S_l2_wS__ = zeros(n_w_max,n_S);
UX_T_M_l2_dM__ = zeros(n_delta_v,n_M);
UX_M_l2_M_ = zeros(n_M,1);
X_wSM___ = zeros(n_w_max,n_S,n_M);
delta_x_wSM___ = zeros(n_w_max,n_S,n_M);
delta_y_wSM___ = zeros(n_w_max,n_S,n_M);
gamma_z_wSM___ = zeros(n_w_max,n_S,n_M);
index_sub_wSM___ = zeros(n_w_max,n_S,n_M);
end;%if flag_optimize_over_gamma_z==0;
if flag_optimize_over_gamma_z==1;
Z_SM__ = zeros(n_S,n_M);
UX_R_CTF_R_S_l2_wS__ = zeros(n_S);
CTF_R_S_l2_wS__ = zeros(n_S);
UX_T_M_l2_dM__ = zeros(n_delta_v,n_M);
UX_M_l2_M_ = zeros(n_M,1);
X_SM__ = zeros(n_S,n_M);
delta_x_SM__ = zeros(n_S,n_M);
delta_y_SM__ = zeros(n_S,n_M);
gamma_z_SM__ = zeros(n_S,n_M);
index_sub_SM__ = zeros(n_S,n_M);
end;%if flag_optimize_over_gamma_z==1;

Z_dwSM____=[];
X_dwSM____=[];
if flag_dwSM;
n_dwSM = n_delta_v*n_w_max*n_S*n_M;
n_dwSM_GB = n_dwSM*n_byte_per_float64/1e9;
if (flag_verbose); tmp_str = 'X_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,n_dwSM_GB)); end;
if (n_dwSM_GB> memory_limit_GB);
disp(sprintf(' %% Warning, n_dwSM_GB %0.2f in %s',n_dwSM_GB,str_thisfunction));
flag_dwSM = 0;
end;%if (n_dwSM_GB> memory_limit_GB);
if (n_dwSM_GB<=memory_limit_GB);
Z_dwSM____ = zeros(n_delta_v,n_w_max,n_S,n_M);
X_dwSM____ = zeros(n_delta_v,n_w_max,n_S,n_M);
flag_dwSM = 1;
end;%if (n_dwSM_GB<=memory_limit_GB);
end;%if flag_dwSM;

if (flag_verbose>0); 
tmp_str = 'UX_R_CTF_R_S_l2_wS__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'CTF_R_S_l2_wS__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'UX_T_M_l2_dM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'UX_M_l2_M_'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
if flag_optimize_over_gamma_z==0;
tmp_str = 'Z_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'X_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if flag_optimize_over_gamma_z==0;
if flag_optimize_over_gamma_z==1;
tmp_str = 'Z_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'X_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if flag_optimize_over_gamma_z==1;
if (flag_dwSM); tmp_str = 'Z_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
if (flag_dwSM); tmp_str = 'X_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (flag_verbose>0); 

%%%%%%%%;
% Now construct the template-norms: ;
% <(R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_,(R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_> ;
% Note that this does not involve collapsing onto principal-modes. ;
%%%%%%%%;
CTF_R_S_l2_wS__ = zeros(n_w_max,n_S);
SS_k_p_wkS__ = conj(S_k_p_wkS__).*S_k_p_wkS__;
SS_k_q_wkS__ = reshape(interp_p_to_q(n_k_p_r,n_w_,n_w_sum,SS_k_p_wkS__),[n_w_sum,n_S]);
CC_k_p_wk_ = conj(CTF_k_p_wk_).*CTF_k_p_wk_;
CC_k_q_wk_ = reshape(interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CC_k_p_wk_),[n_w_sum,1]);
CTF_R_S_l2_wS__ = ifft(squeeze(sum(bsxfun(@times,reshape(bsxfun(@times,CC_k_q_wk_,conj(SS_k_q_wkS__)),[n_w_max,n_k_p_r,n_S]),reshape(weight_2d_k_p_r_,[1,n_k_p_r,1])),1+1)));
CTF_R_S_l2_wS__ = real(CTF_R_S_l2_wS__);
%%%%%%%%;
% Now re-construct the template-norms, this time limited to radial principal-modes: ;
% <((R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_) * pm_wUX_kn__),((R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_) * pm_wUX_kn__)> ;
% Note that this does yes involve collapsing onto principal-modes. ;
%%%%%%%%;
pm_n_k_p_r = pm_n_UX_rank; pm_n_w_max = n_w_max;
pm_n_w_ = pm_n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_sum = pm_n_k_p_r*pm_n_w_max;
UX_R_CTF_R_S_l2_wS__ = zeros(n_w_max,n_S);
pm_wUX_kn__ = diag(pm_X_weight_r_)*pm_UX_kn__;
UX_SS_k_q_wnS__ = reshape(permute(pagemtimes(transpose(pm_wUX_kn__),permute(reshape(SS_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),1+[1,0,2])),1+[1,0,2]),[pm_n_w_sum,n_S]);
UX_CC_k_q_wn_ = reshape(reshape(CC_k_q_wk_,[n_w_max,n_k_p_r])*pm_wUX_kn__,[pm_n_w_sum,1]);
UX_R_CTF_R_S_l2_wS__ = ifft(squeeze(sum(reshape(bsxfun(@times,UX_CC_k_q_wn_,conj(UX_SS_k_q_wnS__)),[pm_n_w_max,pm_n_k_p_r,n_S]),1+1)));
UX_R_CTF_R_S_l2_wS__ = real(UX_R_CTF_R_S_l2_wS__);
%%%%%%%%

flag_continue=1;
n_M_per_Mbatch = max(1,n_M_per_Mbatch);
n_S_per_Sbatch = max(1,n_S_per_Sbatch);
n_dwSM = n_delta_v*n_w_max*n_S_per_Sbatch*n_M_per_Mbatch;
n_dwSM_GB = n_dwSM*n_byte_per_float64/1e9;
if (n_dwSM_GB> memory_limit_GB); disp(sprintf(' %% Warning, n_dwSM_GB %0.2f > %0.2f in %s',n_dwSM_GB,memory_limit_GB,str_thisfunction)); return; end;

n_Mbatch = ceil(n_M/max(1,n_M_per_Mbatch));
if (flag_verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/max(1,n_S_per_Sbatch));
if (flag_verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nMbatch=0:n_Mbatch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
index_nM_in_Mbatch_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
index_nM_in_Mbatch_ = index_nM_in_Mbatch_(1+efind(index_nM_in_Mbatch_<n_M)); n_M_sub = numel(index_nM_in_Mbatch_);
if (flag_verbose>1); disp(sprintf(' %% nMbatch %d/%d index_nM_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_nM_in_Mbatch_(1+0),index_nM_in_Mbatch_(1+n_M_sub-1))); end;
if (flag_verbose>0 & mod(nMbatch,1)==0); disp(sprintf(' %% nMbatch %d/%d index_nM_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_nM_in_Mbatch_(1+0),index_nM_in_Mbatch_(1+n_M_sub-1))); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (n_M_sub>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t = tic();
M_sub_k_p_wkM__ = M_k_p_wkM__(:,1+index_nM_in_Mbatch_);
M_sub_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_sub_k_p_wkM__);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% M_sub_k_q_wkM__: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: M_sub_k_q_wkM__',str_thisfunction),tmp_t);
UX_M_sub_l2_M_ = zeros(n_M_sub,1);
UX_T_M_sub_l2_dM__ = zeros(n_delta_v,n_M_sub);
tmp_t = tic();
CTF_M_sub_k_p_wkM__ = bsxfun(@times,CTF_k_p_wk_,M_sub_k_p_wkM__);
CTF_M_sub_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_M_sub_k_p_wkM__);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% CTF_M_sub_k_q_wkM__: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: CTF_M_sub_k_q_wkM__',str_thisfunction),tmp_t);
%%%%;
% Prepare quasi-images. ;
%%%%;
tmp_t = tic();
svd_V_UX_CTF_M_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,CTF_M_sub_k_q_wkM__,pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_V_UX_CTF_M_sub_lwnM____: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_V_UX_CTF_M_sub_lwnM____',str_thisfunction),tmp_t);
tmp_t = tic();
svd_V_UX_M_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_sub_k_q_wkM__,pm_n_UX_rank,pm_UX_kn__,pm_X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_V_UX_M_sub_lwnM____: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_V_UX_M_sub_lwnM____',str_thisfunction),tmp_t);
tmp_t = tic();
svd_V_UX_CTF_M_sub_nMwl____ = permute(svd_V_UX_CTF_M_sub_lwnM____,1+[2,3,1,0]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_V_UX_CTF_M_sub_nMwl____: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_V_UX_CTF_M_sub_nMwl____',str_thisfunction),tmp_t);
%%%%;
% Now calculate norms of the translated images. ;
%%%%;
tmp_t = tic();
UX_T_M_sub_l2_dM__ = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_V_UX_M_sub_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% tfpmh_UX_T_M_sub_l2_dm__1: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: tfpmh_UX_T_M_sub_l2_dm__1',str_thisfunction),tmp_t);
tmp_index_d0 = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); assert(numel(tmp_index_d0)>=1); tmp_index_d0=tmp_index_d0(1+0); %<-- should be a single index corresponding to zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_T_M_sub_l2_dM__(1+tmp_index_d0,:),[n_M_sub,1]);
%%%%;
% Store results. ;
%%%%;
UX_M_l2_M_(1+index_nM_in_Mbatch_) = UX_M_sub_l2_M_;
UX_T_M_l2_dM__(:,1+index_nM_in_Mbatch_) = UX_T_M_sub_l2_dM__;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nSbatch=0:n_Sbatch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
index_nS_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_nS_in_Sbatch_ = index_nS_in_Sbatch_(1+efind(index_nS_in_Sbatch_<n_S)); n_S_sub = numel(index_nS_in_Sbatch_);
if (flag_verbose>2); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
if (flag_verbose>1 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_nS_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_nS_in_Sbatch_(1+0),index_nS_in_Sbatch_(1+n_S_sub-1))); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (n_S_sub>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
S_sub_k_p_wkS__ = S_k_p_wkS__(:,1+index_nS_in_Sbatch_);
S_sub_k_q_wkS__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_sub_k_p_wkS__);
S_sub_k_q_wSk___ = permute(reshape(S_sub_k_q_wkS__,[n_w_max,n_k_p_r,n_S_sub]),1+[0,2,1]);
%%%%;
% Prepare UX_S_k_q_wnS__. ;
%%%%;
tmp_t = tic();
UX_S_sub_k_q_wnS__ = reshape(permute(reshape(reshape(S_sub_k_q_wSk___,[n_w_max*n_S_sub,n_k_p_r])*diag(pm_X_weight_r_)*pm_UX_kn__,[n_w_max,n_S_sub,pm_n_UX_rank]),1+[0,2,1]),[n_w_max*pm_n_UX_rank,n_S_sub]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_S_sub_k_q_wnS__: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: UX_S_sub_k_q_wnS__',str_thisfunction),tmp_t);
tmp_t = tic();
UX_S_sub_k_q_nSw___ = permute(reshape(UX_S_sub_k_q_wnS__,[n_w_max,pm_n_UX_rank,n_S_sub]),1+[1,2,0]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% UX_S_sub_k_q_nSw__: %0.2fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: UX_S_sub_k_q_nSw__',str_thisfunction),tmp_t);
%%%%;
% Calculate innerproduct Z_sub_dwSM____. ;
%%%%;
tmp_t = tic();
svd_S_sub_V_UX_CTF_M_sub_SMwl____ = pagemtimes(pagectranspose(UX_S_sub_k_q_nSw___),svd_V_UX_CTF_M_sub_nMwl____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_S_sub_V_UX_CTF_M_sub_SMwl____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_S_sub_V_UX_CTF_M_sub_SMwl____',str_thisfunction),tmp_t);
tmp_t = tic();
svd_S_sub_V_UX_CTF_M_sub_lwSM____ = ifft(permute(svd_S_sub_V_UX_CTF_M_sub_SMwl____,1+[3,2,0,1]),[],1+1)*n_w_max;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_S_sub_V_UX_CTF_M_sub_lwSM____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_S_sub_V_UX_CTF_M_sub_lwSM____',str_thisfunction),tmp_t);
tmp_t = tic();
svd_UES_S_sub_V_UX_CTF_M_sub_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_S_sub_V_UX_CTF_M_sub_lwSM____,[n_svd_l,n_w_max*n_S_sub*n_M_sub]),[n_delta_v,n_w_max,n_S_sub,n_M_sub]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_UES_S_sub_V_UX_CTF_M_sub_dwSM____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: svd_UES_S_sub_V_UX_CTF_M_sub_dwSM____',str_thisfunction),tmp_t);
Z_sub_dwSM____  = real(svd_UES_S_sub_V_UX_CTF_M_sub_dwSM____);
%%%%;
% Calculate correlation. ;
%%%%;
tmp_t = tic();
UX_R_CTF_R_S_sub_l2_wS__ = UX_R_CTF_R_S_l2_wS__(:,1+index_nS_in_Sbatch_);
X_sub_dwSM____ = bsxfun(@rdivide,bsxfun(@rdivide,Z_sub_dwSM____,max(1e-12,reshape(sqrt(UX_R_CTF_R_S_sub_l2_wS__),[1,n_w_max,n_S_sub,1]))),max(1e-12,reshape(sqrt(UX_T_M_sub_l2_dM__),[n_delta_v,1,1,n_M_sub])));
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_sub_dwSM____: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: X_sub_dwSM____',str_thisfunction),tmp_t);
%%%%;
% Store results. ;
%%%%;
if flag_dwSM;
Z_dwSM____(:,:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = Z_sub_dwSM____;
X_dwSM____(:,:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = X_sub_dwSM____;
end;%if flag_dwSM;
%%%%%%%%;
if flag_optimize_over_gamma_z==0;
tmp_t = tic();
n_wSM = n_w_max*n_S_sub*n_M_sub;
[tmp_X_wSM_,tmp_ij_delta_wSM_] = max(reshape(X_sub_dwSM____,[n_delta_v,n_wSM]),[],1+0); %<-- maximize correlation. ;
tmp_index_delta_wSM_ = tmp_ij_delta_wSM_-1;
assert(min(tmp_index_delta_wSM_,[],'all')>=0); assert(max(tmp_index_delta_wSM_,[],'all')<=n_delta_v-1);
tmp_X_wSM___ = reshape(tmp_X_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_index_all_wSM_ = tmp_index_delta_wSM_ + reshape([0:n_wSM-1],size(tmp_index_delta_wSM_))*n_delta_v;
assert(min(tmp_index_all_wSM_,[],'all')>=0); assert(max(tmp_index_all_wSM_,[],'all')<=n_delta_v*n_wSM-1);
tmp_Z_wSM___ = reshape(Z_sub_dwSM____(1+tmp_index_all_wSM_),[n_w_max,n_S_sub,n_M_sub]);
tmp_index_delta_wSM___ = reshape(tmp_index_delta_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_delta_x_wSM___ = FTK.delta_x_(1+tmp_index_delta_wSM___);
tmp_delta_y_wSM___ = FTK.delta_y_(1+tmp_index_delta_wSM___);
tmp_gamma_z_wSM___ = 2*pi*(0:n_w_max-1)/max(1,n_w_max);
Z_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = tmp_Z_wSM___;
X_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = tmp_X_wSM___;
delta_x_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_delta_x_wSM___,[n_w_max,n_S_sub,n_M_sub]);
delta_y_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_delta_y_wSM___,[n_w_max,n_S_sub,n_M_sub]);
gamma_z_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = repmat(tmp_gamma_z_wSM___(:),[1,n_S_sub,n_M_sub]);
index_sub_wSM___(:,1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_index_delta_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_wSM___: time %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: X_wSM___',str_thisfunction),tmp_t);
end;%if flag_optimize_over_gamma_z==0;
%%%%%%%%;
if flag_optimize_over_gamma_z==1;
tmp_t = tic();
n_dw = n_delta_v*n_w_max; n_SM = n_S_sub*n_M_sub;
[tmp_X_SM_,tmp_ij_dw_SM_] = max(reshape(X_sub_dwSM____,[n_dw,n_SM]),[],1+0); %<-- maximize correlation. ;
tmp_index_dw_SM_ = tmp_ij_dw_SM_-1;
assert(min(tmp_index_dw_SM_,[],'all')>=0); assert(max(tmp_index_dw_SM_,[],'all')<=n_dw-1);
tmp_X_SM__ = reshape(tmp_X_SM_,[n_S_sub,n_M_sub]);
tmp_index_all_SM_ = tmp_index_dw_SM_ + reshape([0:n_SM-1],size(tmp_index_dw_SM_))*n_dw;
assert(min(tmp_index_all_SM_,[],'all')>=0); assert(max(tmp_index_all_SM_,[],'all')<=n_dw*n_SM-1);
tmp_Z_SM__ = reshape(Z_sub_dwSM____(1+tmp_index_all_SM_),[n_S_sub,n_M_sub]);
tmp_index_dw_SM__ = reshape(tmp_index_dw_SM_,[n_S_sub,n_M_sub]);
tmp_index_delta_SM__ = mod(tmp_index_dw_SM__,n_delta_v);
tmp_index_gamma_SM__ = (tmp_index_dw_SM__ - tmp_index_delta_SM__)/max(1,n_delta_v);
assert(min(tmp_index_delta_SM__,[],'all')>=0); assert(max(tmp_index_delta_SM__,[],'all')<=n_delta_v-1);
assert(min(tmp_index_gamma_SM__,[],'all')>=0); assert(max(tmp_index_gamma_SM__,[],'all')<=n_w_max-1);
tmp_index_delta_SM__ = reshape(tmp_index_delta_SM__,[n_S_sub,n_M_sub]);
tmp_index_gamma_SM__ = reshape(tmp_index_gamma_SM__,[n_S_sub,n_M_sub]);
tmp_delta_x_SM__ = FTK.delta_x_(1+tmp_index_delta_SM__);
tmp_delta_y_SM__ = FTK.delta_y_(1+tmp_index_delta_SM__);
tmp_gamma_z_SM__ = 2*pi*(tmp_index_gamma_SM__)/max(1,n_w_max);
Z_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = tmp_Z_SM__;
X_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = tmp_X_SM__;
delta_x_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_delta_x_SM__,[n_S_sub,n_M_sub]);
delta_y_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_delta_y_SM__,[n_S_sub,n_M_sub]);
gamma_z_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_gamma_z_SM__,[n_S_sub,n_M_sub]);
index_sub_SM__(1+index_nS_in_Sbatch_,1+index_nM_in_Mbatch_) = reshape(tmp_index_dw_SM_,[n_S_sub,n_M_sub]);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_SM__: time %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,sprintf('%s: X_SM__',str_thisfunction),tmp_t);
end;%if flag_optimize_over_gamma_z==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (n_S_sub>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nSbatch=0:n_Sbatch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if (n_M_sub>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nMbatch=0:n_Mbatch-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if flag_optimize_over_gamma_z==1;
Z_wSM___ = Z_SM__;
X_wSM___ = X_SM__;
delta_x_wSM___ = delta_x_SM__;
delta_y_wSM___ = delta_y_SM__;
gamma_z_wSM___ = gamma_z_SM__;
index_sub_wSM___ = index_sub_SM__;
end;%if flag_optimize_over_gamma_z==1;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

