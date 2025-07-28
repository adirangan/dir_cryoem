function ...
[ ...
 parameter ...
,X_gpu_wSM___ ...
,delta_x_gpu_wSM___ ...
,delta_y_gpu_wSM___ ...
,gamma_z_gpu_wSM___ ...
,I_value_gpu_wSM___ ...
] = ...
ampmh_X_gpu_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_gpu_wnS__ ...
,CTF_UX_S_l2_gpu_ ...
,n_M ...
,svd_VUXM_gpu_lwnM____ ...
,UX_M_l2_gpu_dM__ ...
);
%%%%%%%%;
% based on ampmh_X_wSM___8.m. ;
%%%%%%%%;

str_thisfunction = 'ampmh_X_gpu_wSM___8';

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_M_per_Mbatch')); parameter.n_M_per_Mbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_optimize_over_gamma_z')); parameter.flag_optimize_over_gamma_z = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_compute_I_value')); parameter.flag_compute_I_value = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_compress_S')); parameter.flag_compress_S = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
%%%%%%%%;
flag_verbose = parameter.flag_verbose;
n_M_per_Mbatch = parameter.n_M_per_Mbatch;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
flag_optimize_over_gamma_z = parameter.flag_optimize_over_gamma_z;
flag_compute_I_value = parameter.flag_compute_I_value;
flag_compress_S = parameter.flag_compress_S;
if flag_compress_S~=0; disp(sprintf(' %% Warning, flag_compress_S not implemented in %s',str_thisfunction)); end;
tolerance_master = parameter.tolerance_master;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

f_zero = gpuArray( single(0.0));
%%%%;
if ~strcmp(class(CTF_UX_S_k_q_gpu_wnS__),'gpuArray');
tmp_t = tic();
CTF_UX_S_k_q_gpu_wnS__ = gpuArray( (CTF_UX_S_k_q_gpu_wnS__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_UX_S_k_q_gpu_wnS__ %0.6fs',tmp_t)); end;%if (flag_verbose>0);
end;%if ~strcmp(class(CTF_UX_S_k_q_gpu_wnS__),'gpuArray');
%%%%;
if ~strcmp(class(CTF_UX_S_l2_gpu_),'gpuArray');
tmp_t = tic();
CTF_UX_S_l2_gpu_ = gpuArray( (CTF_UX_S_l2_gpu_));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% CTF_UX_S_l2_gpu_ %0.6fs',tmp_t)); end;%if (flag_verbose>0);
end;%if ~strcmp(class(CTF_UX_S_l2_gpu_),'gpuArray');
%%%%;
if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');
tmp_t = tic();
svd_VUXM_gpu_lwnM____ = gpuArray( (svd_VUXM_gpu_lwnM____));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_VUXM_gpu_lwnM____ %0.6fs',tmp_t)); end;%if (flag_verbose>0);
end;%if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');
%%%%;
if ~strcmp(class(UX_M_l2_gpu_dM__),'gpuArray');
tmp_t = tic();
UX_M_l2_gpu_dM__ = gpuArray( (UX_M_l2_gpu_dM__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% UX_M_l2_gpu_dM__ %0.6fs',tmp_t)); end;%if (flag_verbose>0);
end;%if ~strcmp(class(UX_M_l2_gpu_dM__),'gpuArray');
%%%%;
tmp_t = tic();
svd_U_d_expiw_gpu_s__ = gpuArray( (FTK.svd_U_d_expiw_s__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_U_d_expiw_gpu_s__ %0.6fs',tmp_t)); end;%if (flag_verbose>0);
%%%%;

I_value_gpu_wSM___ = [];
if (flag_optimize_over_gamma_z == 0);
X_gpu_wSM___ = zeros(n_w_max,n_S,n_M,'like',f_zero);
delta_x_gpu_wSM___ = zeros(n_w_max,n_S,n_M,'like',f_zero);
delta_y_gpu_wSM___ = zeros(n_w_max,n_S,n_M,'like',f_zero);
gamma_z_gpu_wSM___ = zeros(n_w_max,n_S,n_M,'like',f_zero);
I_value_gpu_wSM___ = [];
if (flag_compute_I_value); I_value_gpu_wSM___ = zeros(n_w_max,n_S,n_M,'like',f_zero); end;
end;%if (flag_optimize_over_gamma_z == 0);
if (flag_optimize_over_gamma_z == 1);
X_gpu_SM__ = zeros(n_S,n_M,'like',f_zero);
delta_x_gpu_SM__ = zeros(n_S,n_M,'like',f_zero);
delta_y_gpu_SM__ = zeros(n_S,n_M,'like',f_zero);
gamma_z_gpu_SM__ = zeros(n_S,n_M,'like',f_zero);
I_value_gpu_SM__ = [];
if (flag_compute_I_value); I_value_gpu_SM__ = zeros(n_S,n_M,'like',f_zero); end;
end;%if (flag_optimize_over_gamma_z == 1);

%CTF_UX_S_k_q_gpu_wnS___ = reshape(CTF_UX_S_k_q_gpu_wnS__(:,1:n_S),[n_w_max,pm_n_UX_rank,n_S]); %<-- used later. ;
%CTF_UX_S_k_q_gpu_nSw___ = permute(CTF_UX_S_k_q_gpu_wnS___,[2,3,1]);
%if (flag_verbose>0);
%tmp_str = 'CTF_UX_S_k_q_gpu_nSw___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
%end;%if (flag_verbose>0);
conj_CTF_UX_S_k_q_gpu_Snw___ = permute(conj(reshape(CTF_UX_S_k_q_gpu_wnS__(:,1:n_S),[n_w_max,pm_n_UX_rank,n_S])),[3,2,1]);

if (flag_verbose>0); 
if (flag_optimize_over_gamma_z == 0);
tmp_str = 'X_gpu_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_gpu_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_gpu_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_gpu_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
if (flag_compute_I_value); tmp_str = 'I_value_gpu_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (flag_optimize_over_gamma_z == 0);
if (flag_optimize_over_gamma_z == 1);
tmp_str = 'X_gpu_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_gpu_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_gpu_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_gpu_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
if (flag_compute_I_value); tmp_str = 'I_value_gpu_SM__'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (flag_optimize_over_gamma_z == 1);
end;%if (flag_verbose>0); 
n_Mbatch = ceil(n_M/n_M_per_Mbatch);
if (flag_verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nMbatch=0:n_Mbatch-1;
index_M_in_Mbatch_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
index_M_in_Mbatch_ = index_M_in_Mbatch_(find(index_M_in_Mbatch_<n_M)); n_M_sub = numel(index_M_in_Mbatch_);
if (flag_verbose>1); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (flag_verbose>0 & mod(nMbatch,1)==0); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (n_M_sub>0);
tmp_t = tic(); nop=0;
svd_VUXM_gpu_nMwl____ = zeros(pm_n_UX_rank,n_M_sub,n_w_max,FTK.n_svd_l,'like',f_zero);
svd_VUXM_gpu_nMwl____ = permute(svd_VUXM_gpu_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]);
nop = nop + numel(svd_VUXM_gpu_nMwl____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_VUXM_gpu_nMwl____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: svd_VUXM_gpu_nMwl____',tmp_t,1,nop);
%%%%%%%%;
tmp_t = tic(); nop=0;
svd_SVUXM_gpu_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l,'like',f_zero);
svd_SVUXM_gpu_SMwl____ = pagemtimes(conj_CTF_UX_S_k_q_gpu_Snw___,svd_VUXM_gpu_nMwl____);
%%%%;
% Note: Depending on memory constraints and batch size, ;
% the loop below might be faster than the pagemtimes call above. ;
%%%%;
flag_loop=0;
if flag_loop;
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_SVUXM_gpu_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_gpu_nSw___(:,:,1+nw))*svd_VUXM_gpu_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
end;%if flag_loop;
%%%%;
nop = nop + FTK.n_svd_l*n_w_max*n_S*pm_n_UX_rank*n_M_sub;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXM_gpu_SMwl____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: svd_VUXM_SMwl____',tmp_t,1,nop);
tmp_t = tic(); nop=0;
svd_SVUXM_gpu_lwSM____ = permute(ifft(permute(svd_SVUXM_gpu_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
nop = nop + numel(svd_SVUXM_gpu_lwSM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXM_gpu_lwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: svd_VUXM_lwSM____',tmp_t,1,nop);
%%%%%%%%;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (flag_verbose>1); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
%%%%%%%%;
svd_SVUXM_gpu_lwsM____ = svd_SVUXM_gpu_lwSM____(:,:,1+index_S_in_Sbatch_,:);
%%%%%%%%;
tmp_t = tic(); nop=0;
svd_USESVUXM_gpu_dwSM____ = real(reshape(svd_U_d_expiw_gpu_s__*reshape(svd_SVUXM_gpu_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub]));
%%%%%%%%;
l2_gpu_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_gpu_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_gpu_dM__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
n2_gpu_dSM___ = 1./max(1e-14,l2_gpu_dSM___);
f2_gpu_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_gpu_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_gpu_dM__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
ss_S_ = reshape(CTF_UX_S_l2_gpu_(1+index_S_in_Sbatch_),[n_S_sub,1]);
nop = nop + FTK.n_delta_v*FTK.n_svd_l*n_w_max*n_S_sub*n_M_sub + n_S_sub*FTK.n_delta_v*n_M_sub + n_S_sub*FTK.n_delta_v*n_M_sub;
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% svd_USESVUXM_gpu_dwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: svd_USESVUXM_gpu_dwSM____',tmp_t,1,nop);
if (nMbatch==0 && nSbatch==0 && flag_verbose>0); 
tmp_str = 'svd_VUXM_gpu_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_gpu_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_gpu_lwsM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_gpu_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_gpu_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_gpu_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'f2_gpu_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (flag_verbose>1); 
tmp_t = tic(); nop=0;
X_sub_gpu_dwSM____ = repmat(reshape(n2_gpu_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_gpu_dwSM____; %<-- correlation. ;
%X_sub_gpu_dwSM____ = bsxfun(@times,reshape(real(n2_gpu_dSM___),[FTK.n_delta_v,1,n_S_sub,n_M_sub]),real(svd_USESVUXM_gpu_dwSM____)); %<-- correlation. ;
if (flag_compute_I_value);
I_value_sub_dwSM____ = repmat(reshape(f2_gpu_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_gpu_dwSM____; %<-- I_value. ;
I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
end;%if (flag_compute_I_value);
nop = nop + numel(X_sub_gpu_dwSM____);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_sub_gpu_dwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: X_sub_gpu_dwSM____',tmp_t,1,nop);
%%%%%%%%;
flag_check=0;
if flag_check;
[tmp_X_gpu_wSM___,tmp_delta_ij___] = max(reshape(X_sub_gpu_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
tmp_X_gpu_wSM___ = reshape(tmp_X_gpu_wSM___,[n_w_max,n_S_sub,n_M_sub]);
tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
[tmp_X_gpu_SM__,tmp_dw_ij__] = max(reshape(X_sub_gpu_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
[tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
tmp_X_gpu_SM__ = reshape(tmp_X_gpu_SM__,[n_S_sub,n_M_sub]);
tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
for nS=0:n_S_sub-1;
for nM=0:n_M_sub-1;
[tmp_X,tmp_gamma_ij] = max(tmp_X_gpu_wSM___(:,1+nS,1+nM));
assert(tmp_X_gpu_wSM___(tmp_gamma_ij,1+nS,1+nM)==tmp_X_gpu_SM__(1+nS,1+nM));
assert(tmp_gamma_ij==tmp_gamma_ij__(1+nS,1+nM));
assert(tmp_delta_ij___(tmp_gamma_ij,1+nS,1+nM)==tmp_delta_ij__(1+nS,1+nM));
end;%for nM=0:n_M_sub-1;
end;%for nS=0:n_S_sub-1;
end;%if flag_check;
%%%%%%%%;
if (flag_optimize_over_gamma_z == 0);
tmp_t = tic(); nop=0;
[tmp_X_gpu_wSM___,tmp_delta_ij___] = max(reshape(X_sub_gpu_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
assert(min(tmp_delta_ij___)>=1); assert(max(tmp_delta_ij___)<=FTK.n_delta_v);
tmp_X_gpu_wSM___ = reshape(tmp_X_gpu_wSM___,[n_w_max,n_S_sub,n_M_sub]);
tmp_delta_ij___ = reshape(tmp_delta_ij___,[n_w_max,n_S_sub,n_M_sub]);
tmp_delta_x___ = FTK.delta_x_(tmp_delta_ij___);
tmp_delta_y___ = FTK.delta_y_(tmp_delta_ij___);
tmp_gamma_z___ = 2*pi*(0:n_w_max-1)/n_w_max;
X_gpu_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_gpu_wSM___;
delta_x_gpu_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x___,[n_w_max,n_S_sub,n_M_sub]);
delta_y_gpu_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y___,[n_w_max,n_S_sub,n_M_sub]);
gamma_z_gpu_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = repmat(tmp_gamma_z___(:),[1,n_S_sub,n_M_sub]);
if (flag_compute_I_value);
tmp_I_value_use_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
tmp_I_value_use_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1,'like',f_zero);
tmp_t2=tic();
for nl=0:n_w_max*n_S_sub*n_M_sub-1;
tmp_I_value_use_wSM_(1+nl) = tmp_I_value_use_dwSM__(tmp_delta_ij___(1+nl),1+nl);
end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
I_value_gpu_wSM___(:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_t2 = toc(tmp_t2); if (flag_verbose>1); disp(sprintf(' %% I_value_gpu_wSM___ %0.6fs',tmp_t2)); end;
end;%if (flag_compute_I_value);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_gpu_wSM___: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: X_gpu_wSM___',tmp_t,1,nop);
end;%if (flag_optimize_over_gamma_z == 0);
%%%%%%%%;
if (flag_optimize_over_gamma_z == 1);
tmp_t = tic(); nop=0;
[tmp_X_gpu_SM__,tmp_dw_ij__] = max(reshape(X_sub_gpu_dwSM____,[FTK.n_delta_v*n_w_max,n_S_sub*n_M_sub]),[],1); %<-- maximize correlation. ;
[tmp_delta_ij__,tmp_gamma_ij__] = ind2sub([FTK.n_delta_v,n_w_max],tmp_dw_ij__);
assert(min(tmp_delta_ij__)>=1); assert(max(tmp_delta_ij__)<=FTK.n_delta_v);
assert(min(tmp_gamma_ij__)>=1); assert(max(tmp_gamma_ij__)<=n_w_max);
tmp_X_gpu_SM__ = reshape(tmp_X_gpu_SM__,[n_S_sub,n_M_sub]);
tmp_delta_ij__ = reshape(tmp_delta_ij__,[n_S_sub,n_M_sub]);
tmp_gamma_ij__ = reshape(tmp_gamma_ij__,[n_S_sub,n_M_sub]);
tmp_delta_x__ = FTK.delta_x_(tmp_delta_ij__);
tmp_delta_y__ = FTK.delta_y_(tmp_delta_ij__);
tmp_gamma_z__ = 2*pi*(tmp_gamma_ij__-1)/n_w_max;
X_gpu_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = tmp_X_gpu_SM__;
delta_x_gpu_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_x__,[n_S_sub,n_M_sub]);
delta_y_gpu_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_delta_y__,[n_S_sub,n_M_sub]);
gamma_z_gpu_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_gamma_z__,[n_S_sub,n_M_sub]);
if (flag_compute_I_value);
tmp_I_value_use_dwSM___ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max,n_S_sub*n_M_sub]);
tmp_I_value_use_SM_ = zeros(n_S_sub*n_M_sub,1,'like',f_zero);
tmp_t2=tic();
for nl=0:n_S_sub*n_M_sub-1;
tmp_I_value_use_SM_(1+nl) = tmp_I_value_use_dwSM___(tmp_delta_ij__(1+nl),tmp_gamma_ij__(1+nl),1+nl);
end;%for nl=0:n_S_sub*n_M_sub-1;
I_value_gpu_SM__(1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = reshape(tmp_I_value_use_SM_,[n_S_sub,n_M_sub]);
tmp_t2 = toc(tmp_t2); if (flag_verbose>1); disp(sprintf(' %% I_value_gpu_SM__ %0.6fs',tmp_t2)); end;
end;%if (flag_compute_I_value);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% X_gpu_wSM___: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_gpu_wSM___8: X_gpu_wSM___',tmp_t,1,nop);
end;%if (flag_optimize_over_gamma_z == 1);
%%%%%%%%;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
end;%if (n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;

if (flag_optimize_over_gamma_z == 1);
X_gpu_wSM___ = X_gpu_SM__;
delta_x_gpu_wSM___ = delta_x_gpu_SM__;
delta_y_gpu_wSM___ = delta_y_gpu_SM__;
gamma_z_gpu_wSM___ = gamma_z_gpu_SM__;
if (flag_compute_I_value); I_value_gpu_wSM___ = I_value_gpu_SM__; end;
end;%if (flag_optimize_over_gamma_z == 1);

if ( (nargout>5) & (isempty(I_value_gpu_wSM___)) ); I_value_gpu_wSM___ = ones(size(X_gpu_wSM___),'like',f_zero); end;

if (flag_verbose>1); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;


