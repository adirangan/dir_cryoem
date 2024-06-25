function ...
[ ...
 parameter ...
,X_dwSM____ ...
,Y_dwSM____ ...
,I_value_dwSM____ ...
] = ...
ampmh_X_dwSM____8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
%%%%%%%%;
% Uses the templates provided in CTF_UX_S_k_q_wnS__. ;
% Calculates correlations between principal-templates and principal-images across translations. ;
% Generally speaking, the output X_dwSM____(1+ndelta,1+nw,1+nS,1+nM) calculates the correlation between: ;
% Rz(+gamma_z) * S_k_p__(:,1+nS) = S_k_p__(Rz(-gamma_z)*k_p_ , 1+nS), ;
% and ;
% exp(-2*pi*i*dot(k_,delta_)) .* M_k_p__(:,1+nM), ;
% where gamma_z = 2*pi*nw/n_w, ;
% and delta_ = [ delta_x ; delta_y ], ;
% where delta_? = FTK.delta_?_(1+ndelta);
% Put more succinctly, up to normalization we have: ;
% X_dwSM____(1+ndelta,1+nw,1+nS,1+nM) = dot( rotate(S,+gamma_z) , transf(M,+delta_) );
%                                    = sum( conj(rotate(S,+gamma_z)) .* transf(M,+delta_) .* weight_2d_ );
% Note that we project X_dwSM____ onto its real component. ;
% The output Y_dwSM____ is analogous to X_dwSM____, with no normalization. ;
%%%%%%%%;
% Note that ampmh_X_dwSM____8 assumes an isotropic CTF. ;
% I.e., the inputs CTF_UX_S_k_q_wnS__ and CTF_UX_S_l2_ ;
% combine an angularly-independent CTF with the templates S_k_q_wnS__. ;
%%%%%%%%
% If requested, the I_value is calculated as follows: ;
% I_value = <M,S> / <M,M> ;
% X = -<IM-S,IM-S> ;
%   = -I*I*<M,M> + 2*I*<M,S> - <S,S> ;
%   = -<M,S><M,S>/<M,M> + 2*<M,S><M,S>/<M,M> - <S,S> ;
%   = <M,S><M,S>/<M,M> - <S,S> ;
%   = <S,S>*(<M,S><M,S>/<M,M>/<S,S> - 1) ;
%   = <S,S>*(correlation^2 - 1) ;
%%%%%%%%;
% Batches images and templates. ;
% Batch sizes are given by: ;
% parameter.n_M_per_Mbatch and parameter.n_S_per_Sbatch. ;
% Default values are 24 and 24. ;
%%%%%%%%;
% Note that the formulae above allow for negative intensities (i.e., I_value can be < 0). ;
% Consequently, image-template pairs exhibiting a negative correlation can actually be quite likely matches, ;
% since the corresponding (negative) intensity produces a large log-likelihood. ;
% Therefore, these formulae should only be used as long as this 'feature' is acceptable. ;
% In this particular function we threshold I_value to be nonnegative. ;
%%%%%%%%;
% extra fields in the 'parameter' structure can be passed in: ;
%%%%%%%%;
% The user can limit the number of principal-modes by setting pm_n_UX_rank_use>0. ;
% by default, the principal-modes used are those with the largest singular-values (i.e., principal-modes 0:pm_n_UX_rank_use-1). ;
%%%%%%%%;
% The user can limit the number of svd-modes by setting either svd_eps_use>0 or n_svd_l_use>0. ;
% by default, the svd-modes used are those with the largest singular-values (i.e., largest FTK.svd_s_). ;
%%%%%%%%;
% The user can limit the number of delta_v values by setting n_delta_v_use>0. ;
% by default, the delta_x and delta_y used are generated by ampmh_FTK_subselect_1.m (guaranteed to include the zero-point). ;
%%%%%%%%;

str_thisfunction = 'ampmh_X_dwSM____8';

if nargin<1;
%%%%%%%%;
disp(sprintf(' %% testing %s (Warning, not tested yet)',str_thisfunction));
%%%%%%%%;
disp('returning'); return;
end;% if nargin<1;

verbose = 0;
if (verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'n_M_per_Mbatch')); parameter.n_M_per_Mbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'n_S_per_Sbatch')); parameter.n_S_per_Sbatch = 24; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'flag_compress_S')); parameter.flag_compress_S = 0; end; %<-- parameter_bookmark. ;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
%%%%%%%%;
n_M_per_Mbatch = parameter.n_M_per_Mbatch;
n_S_per_Sbatch = parameter.n_S_per_Sbatch;
flag_compress_S = parameter.flag_compress_S;
tolerance_master = parameter.tolerance_master;

flag_compute_Y = (nargout>=1+2);
flag_compute_I_value = (nargout>=1+3);

tmp_t = tic(); nop=0;
[ ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_use_ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM_use__ ...
] = ...
ampmh_X_wSM_reduce_1( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,n_S ...
,CTF_UX_S_k_q_wnS__ ...
,CTF_UX_S_l2_ ...
,n_M ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
);
nop = nop + numel(svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% reduce: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_wSM___8: reduce',tmp_t,1,nop);

I_value_dwSM____ = [];
X_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S,n_M);
delta_x_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S,n_M);
delta_y_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S,n_M);
gamma_z_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S,n_M);
I_value_dwSM____ = [];
if (flag_compute_I_value); I_value_dwSM____ = zeros(FTK.n_delta_v,n_w_max,n_S,n_M); end;

if (flag_compress_S==0);
CTF_UX_S_k_q_wnS___ = reshape(CTF_UX_S_k_q_wnS__(:,1:n_S),[n_w_max,pm_n_UX_rank,n_S]); %<-- used later. ;
CTF_UX_S_k_q_nSw___ = permute(CTF_UX_S_k_q_wnS___,[2,3,1]);
if (verbose>0);
tmp_str = 'CTF_UX_S_k_q_nSw___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>0);
end;%if (flag_compress_S==0);
if (flag_compress_S==1);
SS_k_q_ = svds(CTF_UX_S_k_q_wnS__,min(n_w_max*pm_n_UX_rank,n_S));
n_S_rank = min(efind(cumsum(SS_k_q_)/sum(SS_k_q_)>1-tolerance_master));
if (verbose>0); disp(sprintf(' %% n_S_rank %d/%d',n_S_rank,min(n_w_max*pm_n_UX_rank,n_S))); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(CTF_UX_S_k_q_wnS__,n_S_rank);
US_CTF_UX_S_k_q_nSw___ = permute(reshape(US_k_q__,[n_w_max,pm_n_UX_rank,n_S_rank]),[2,3,1]);
end;%if (flag_compress_S==1);

if (verbose>0); 
tmp_str = 'X_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_x_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_y_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'gamma_z_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
if (flag_compute_I_value); tmp_str = 'I_value_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9)); end;
end;%if (verbose>0); 

n_Mbatch = ceil(n_M/n_M_per_Mbatch);
if (verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nMbatch=0:n_Mbatch-1;
index_M_in_Mbatch_ = nMbatch*n_M_per_Mbatch + (0:n_M_per_Mbatch-1);
index_M_in_Mbatch_ = index_M_in_Mbatch_(find(index_M_in_Mbatch_<n_M)); n_M_sub = numel(index_M_in_Mbatch_);
if (verbose>1); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (verbose>0 & mod(nMbatch,1)==0); disp(sprintf(' %% nMbatch %d/%d index_M_in_Mbatch_ %d-->%d',nMbatch,n_Mbatch,index_M_in_Mbatch_(1+0),index_M_in_Mbatch_(1+n_M_sub-1))); end;
if (n_M_sub>0);
tmp_t = tic(); nop=0;
svd_VUXM_nMwl____ = zeros(pm_n_UX_rank,n_M_sub,n_w_max,FTK.n_svd_l);
svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+index_M_in_Mbatch_),[3,4,2,1]);
nop = nop + numel(svd_VUXM_nMwl____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_nMwl____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_VUXM_nMwl____',tmp_t,1,nop);
%%%%%%%%;
if (flag_compress_S==0);
tmp_t = tic(); nop=0;
svd_SVUXM_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
nop = nop + FTK.n_svd_l*n_w_max*n_S*pm_n_UX_rank*n_M_sub;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_VUXM_SMwl____',tmp_t,1,nop);
tmp_t = tic(); nop=0;
svd_SVUXM_lwSM____ = permute(ifft(permute(svd_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
nop = nop + numel(svd_SVUXM_lwSM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_lwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_VUXM_lwSM____',tmp_t,1,nop);
end;%if (flag_compress_S==0);
%%%%%%%%;
if (flag_compress_S==1);
tmp_t = tic(); nop=0;
svd_USVUXM_SMwl____ = zeros(n_S_rank,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_USVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(US_CTF_UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
nop = nop + FTK.n_svd_l*n_w_max*n_S_rank*pm_n_UX_rank*n_M_sub;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USVUXM_SMwl____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_USVUXM_SMwl____',tmp_t,1,nop);
tmp_t = tic(); nop=0;
svd_USVUXM_SMwl____ = permute(ifft(permute(svd_USVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[3,4,1,2]);
nop = nop + numel(svd_USVUXM_SMwl____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USVUXM_lwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_USVUXM_SMwl____ (permute)',tmp_t,1,nop);
end;%if (flag_compress_S==1);
%%%%%%%%;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (verbose>1); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
%%%%%%%%;
if (flag_compress_S==0);
svd_SVUXM_lwsM____ = svd_SVUXM_lwSM____(:,:,1+index_S_in_Sbatch_,:);
end;%if (flag_compress_S==0);
if (flag_compress_S==1);
tmp_t = tic(); nop=0;
svd_SVUXM_lwsM____ = permute(reshape(VS_k_q__(1+index_S_in_Sbatch_,:)*SS_k_q__*reshape(svd_USVUXM_SMwl____,[n_S_rank,n_M_sub*n_w_max*FTK.n_svd_l]),[n_S_sub,n_M_sub,n_w_max,FTK.n_svd_l]),[4,3,1,2]);
nop = nop + n_S_sub*n_S_rank*n_S_rank + n_S_sub*n_S_rank*n_M_sub*n_w_max*FTK.n_svd_l;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_lwsM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_SVUXM_lwsM____',tmp_t,1,nop);
end;%if (flag_compress_S==1);
%%%%%%%%;
tmp_t = tic(); nop=0;
svd_USESVUXM_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwsM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub]));
%%%%%%%%;
if (flag_compute_Y);
y2_dSM___ = permute(reshape(bsxfun(@plus,reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]),reshape(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_),[1,FTK.n_delta_v*n_M_sub])),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
end;%if (flag_compute_Y);
l2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
n2_dSM___ = 1./max(1e-14,l2_dSM___);
f2_dSM___ = permute(reshape(reshape(sqrt(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM_use__(:,1+index_M_in_Mbatch_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
ss_S_ = reshape(CTF_UX_S_l2_use_(1+index_S_in_Sbatch_),[n_S_sub,1]);
nop = nop + FTK.n_delta_v*FTK.n_svd_l*n_w_max*n_S_sub*n_M_sub + n_S_sub*FTK.n_delta_v*n_M_sub + n_S_sub*FTK.n_delta_v*n_M_sub;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USESVUXM_dwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: svd_USESVUXM_dwSM____',tmp_t,1,nop);
if (nMbatch==0 && nSbatch==0 && verbose>0); 
tmp_str = 'svd_VUXM_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_lwsM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'f2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = tic(); nop=0;
if (flag_compute_Y);
Y_sub_dwSM____ = repmat(reshape(y2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) - 2*svd_USESVUXM_dwSM____; %<-- innerproduct. ;
end;%if (flag_compute_Y);
X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ;
if (flag_compute_I_value);
I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
end;%if (flag_compute_I_value);
nop = nop + numel(X_sub_dwSM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_sub_dwSM____: %0.6f',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_dwSM____8: X_sub_dwSM____',tmp_t,1,nop);
%%%%%%%%;
X_dwSM____(:,:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = X_sub_dwSM____;
if (flag_compute_Y);
Y_dwSM____(:,:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = Y_sub_dwSM____;
end;%if (flag_compute_Y);
if (flag_compute_I_value);
I_value_dwSM____(:,:,1+index_S_in_Sbatch_,1+index_M_in_Mbatch_) = I_value_use_dwSM____;
end;%if (flag_compute_I_value);
%%%%%%%%;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
end;%if (n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;

if (verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;



