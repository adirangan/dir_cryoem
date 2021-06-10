function ...
[ ...
 parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_single_neighborhood_SM__9( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,n_CTF_rank ...
,UCTF_UX_S_k_q_wnSc____ ...
,n_M ...
,CTF_UX_S_l2_SM__ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
%%%%%%%%;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_single_neighborhood_SM__9]')); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'n_M_per_Mbatch')); parameter.n_M_per_Mbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_optimize_over_gamma_z')); parameter.flag_optimize_over_gamma_z = 1; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_compute_I_value')); parameter.flag_compute_I_value = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
%%%%%%%%;
n_M_per_Mbatch = parameter.n_M_per_Mbatch;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
flag_optimize_over_gamma_z = parameter.flag_optimize_over_gamma_z;
flag_compute_I_value = parameter.flag_compute_I_value;
tolerance_master = parameter.tolerance_master;

I_value_SM__ = [];
X_SM__ = zeros(n_S,n_M);
delta_x_SM__ = zeros(n_S,n_M);
delta_y_SM__ = zeros(n_S,n_M);
gamma_z_SM__ = zeros(n_S,n_M);
if (flag_compute_I_value); I_value_SM__ = zeros(n_S,n_M); end;
if (verbose>0); 
tmp_str = 'X_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
if (flag_compute_I_value); tmp_str = 'I_value_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (verbose>0); 

tmp_t = tic();
UCTF_UX_S_k_q_nScw____ = permute(UCTF_UX_S_k_q_wnSc____,[2,3,4,1]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UCTF_UX_S_k_q_nScw____: %0.6f',tmp_t)); end;
if (verbose>0);
tmp_str = 'UCTF_UX_S_k_q_nScw____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>0);

n_Mbatch = ceil(n_M/n_M_per_Mbatch);
if (verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nMbatch=0:n_Mbatch-1;
index_M_in_Mbatch_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
index_M_in_Mbatch_ = index_M_in_Mbatch_(find(index_M_in_Mbatch_<n_M)); n_M_sub = numel(index_M_in_Mbatch_);
tmp_VSCTF_Mc__ = VSCTF_Mc__(1+index_M_in_Mbatch_,:);
if (verbose>1); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (verbose>0 & mod(nMbatch,1)==0); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (n_M_sub>0);
tmp_t = tic();
svd_VUXM_nMwl____ = zeros(pm_n_UX_rank,n_M_sub,n_w_max,FTK.n_svd_l);
svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]);
%%%%%%%%;
svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
for nCTF_rank=0:n_CTF_rank-1;
svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = svd_SVUXM_SMwl____(:,:,1+nw,1+nl) + ctranspose(UCTF_UX_S_k_q_nScw____(:,:,1+nCTF_rank,1+nw))*(svd_VUXM_nMwl____(:,:,1+nw,1+nl)*sparse(1:n_M_sub,1:n_M_sub,tmp_VSCTF_Mc__(:,1+nCTF_rank)));
end;%for nCTF_rank=0:n_CTF_rank-1;
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
tmp_t = tic();
svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_lwSM____: %0.6f',tmp_t)); end;
%%%%%%%%;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (verbose>1); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_t = tic();
%%%%%%%%;
svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:);
%%%%%%%%;
svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub]));
%%%%%%%%;
l2_dSM___ = zeros(FTK.n_delta_v,n_S_sub,n_M_sub);
for nM_sub=0:n_M_sub-1;
l2_dSM___(:,:,1+nM_sub) = sqrt(UX_M_l2_dM__(:,1+index_M_in_Mbatch_(1+nM_sub)))*transpose(sqrt(CTF_UX_S_l2_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_(1+nM_sub))));
end;%for nM_sub=0:n_M_sub-1;
n2_dSM___ = 1./max(1e-14,l2_dSM___);
f2_dSM___ = zeros(FTK.n_delta_v,n_S_sub,n_M_sub);
for nM_sub=0:n_M_sub-1;
f2_dSM___(:,:,1+nM_sub) = 1./max(1e-12,sqrt(UX_M_l2_dM__(:,1+index_M_in_Mbatch_(1+nM_sub))))*transpose(sqrt(CTF_UX_S_l2_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_(1+nM_sub))));
end;%for nM_sub=0:n_M_sub-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USESVUXM_dwSM____: %0.6f',tmp_t)); end;
if (nMbatch==0 && nSbatch==0 && verbose>0); 
tmp_str = 'svd_VUXM_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_lwsM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'f2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = tic();
X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ;
%X_sub_dwSM____ = bsxfun(@times,reshape(real(n2_dSM___),[FTK.n_delta_v,1,n_S_sub,n_M_sub]),real(svd_USESVUXM_dwSM____)); %<-- correlation. ;
if (flag_compute_I_value);
I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
end;%if (flag_compute_I_value);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_sub_dwSM____: %0.6f',tmp_t)); end;
%%%%%%%%;
flag_check=0;
if flag_check;
[tmp_X_wSM___,tmp_delta_ij___] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
tmp_X_wSM___ = reshape(tmp_X_wSM___,[n_w_max,n_S_sub,n_M_sub]);
tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
[tmp_X_SM__,tmp_dw_ij__] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
[tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
tmp_X_SM__ = reshape(tmp_X_SM__,[n_S_sub,n_M_sub]);
tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
for nS=0:n_S_sub-1;
for nM=0:n_M_sub-1;
[tmp_X,tmp_gamma_ij] = max(tmp_X_wSM___(:,1+nS,1+nM));
assert(tmp_X_wSM___(tmp_gamma_ij,1+nS,1+nM)==tmp_X_SM__(1+nS,1+nM));
assert(tmp_gamma_ij==tmp_gamma_ij__(1+nS,1+nM));
assert(tmp_delta_ij___(tmp_gamma_ij,1+nS,1+nM)==tmp_delta_ij__(1+nS,1+nM));
end;%for nM=0:n_M_sub-1;
end;%for nS=0:n_S_sub-1;
end;%if flag_check;
%%%%%%%%;
tmp_t = tic();
[tmp_X_SM__,tmp_dw_ij__] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
[tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
assert(min(tmp_delta_ij__)>=1); assert(max(tmp_delta_ij__)<=FTK.n_delta_v);
assert(min(tmp_gamma_ij__)>=1); assert(max(tmp_gamma_ij__)<=n_w_max);
tmp_X_SM__ = reshape(tmp_X_SM__,[n_S_sub,n_M_sub]);
tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
tmp_delta_x__ = FTK.delta_x_(tmp_delta_ij__);
tmp_delta_y__ = FTK.delta_y_(tmp_delta_ij__);
tmp_gamma_z__ = 2*pi*(tmp_gamma_ij__-1)/n_w_max;
X_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_SM__;
delta_x_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x__,[n_S_sub,n_M_sub]);
delta_y_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y__,[n_S_sub,n_M_sub]);
gamma_z_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_gamma_z__,[n_S_sub,n_M_sub]);
if (flag_compute_I_value);
tmp_I_value_use_dwSM___ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max,n_S_sub*n_M_sub]);
tmp_I_value_use_SM_ = zeros(n_S_sub*n_M_sub,1);
tmp_t2=tic();
for nl=0:n_S_sub*n_M_sub-1;
tmp_I_value_use_SM_(1+nl) = tmp_I_value_use_dwSM___(tmp_delta_ij__(1+nl),tmp_gamma_ij__(1+nl),1+nl);
end;%for nl=0:n_S_sub*n_M_sub-1;
I_value_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_SM_,[n_S_sub,n_M_sub]);
tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% I_value_SM__ %0.6fs',tmp_t2)); end;
end;%if (flag_compute_I_value);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
%%%%%%%%;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
end;%if (n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;

if ( (nargout>5) & (isempty(I_value_SM__)) ); I_value_SM__ = ones(n_S,n_M); end;

if (verbose>1); disp(sprintf(' %% [finished ampmh_X_single_neighborhood_SM__9]')); end;


