function ...
[ ...
 L_wSM___ ...
,delta_j_wSM___ ...
,I_value_wSM___ ...
] = ...
ampmh_X_wSM___3(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_Sbatch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M...
,n_M_Mbatch...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
);
%%%%%%%%;
% Calculates correlations and log-likelihoods between principal-templates and principal-images across translations. ;
% The I_value is calculated as follows: ;
% I_value = <M,S> / <M,M> ;
% X = -<IM-S,IM-S> ;
%   = -I*I*<M,M> + 2*I*<M,S> - <S,S> ;
%   = -<M,S><M,S>/<M,M> + 2*<M,S><M,S>/<M,M> - <S,S> ;
%   = <M,S><M,S>/<M,M> - <S,S> ;
%   = <S,S>*(<M,S><M,S>/<M,M>/<S,S> - 1) ;
%   = <S,S>*(correlation^2 - 1) ;
% Batches images and uses svd for templates. ;
% Also batches templates. ;
%%%%%%%%;
% Note that this allows for negative intensities (i.e., I_value can be < 0). ;
% Consequently, image-template pairs exhibiting a negative correlation can actually be quite likely matches, ;
% since the corresponding (negative) intensity produces a large log-likelihood. ;
% Therefore, this subroutine should only be used as long as this 'feature' is acceptable. ;
%%%%%%%%;
verbose=0;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_wSM___3]')); end;
n_w_max = max(n_w_);
L_wSM___ = zeros(n_w_max,n_S,n_M);
delta_j_wSM___ = zeros(n_w_max,n_S,n_M);
I_value_wSM___ = zeros(n_w_max,n_S,n_M);
UX_US_CTF_k_q_nSw___ = permute(reshape(US_k_q__,[n_w_max,pm_n_UX_rank,n_S_rank]),[2,3,1]);
if (verbose>0); 
tmp_str = 'L_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_j_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'I_value_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'UX_US_CTF_k_q_nSw___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>0); 
n_Mbatch = ceil(n_M/n_M_Mbatch);
if (verbose>1); disp(sprintf(' %% n_Mbatch %d',n_Mbatch)); end;
n_Sbatch = ceil(n_S/n_S_Sbatch);
if (verbose>1); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nMbatch=0:n_Mbatch-1;
M_ij_ = nMbatch*n_M_Mbatch + (0:n_M_Mbatch-1);
M_ij_ = M_ij_(find(M_ij_<n_M)); n_M_sub = numel(M_ij_);
if (verbose>1); disp(sprintf(' %% nMbatch %d/%d M_ij_ %d-->%d',nMbatch,n_Mbatch,M_ij_(1+0),M_ij_(1+n_M_sub-1))); end;
if (n_M_sub>0);
tmp_t = tic();
svd_VUXM_nMwl____ = zeros(pm_n_UX_rank,n_M_sub,n_w_max,FTK.n_svd_l);
svd_VUXM_nMwl____ = permute(svd_VUXM_lwnM____(:,:,:,1+M_ij_),[3,4,2,1]);
svd_USVUXM_SMwl____ = zeros(n_S_rank,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_USVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(UX_US_CTF_k_q_nSw___(:,:,1+nw))*svd_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USVUXM_SMwl____: %0.6f',tmp_t)); end;
tmp_t = tic();
svd_USVUXM_SMwl____ = permute(ifft(permute(svd_USVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[3,4,1,2]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USVUXM_SMwl____: %0.6f',tmp_t)); end;
for nSbatch=0:n_Sbatch-1;
S_ij_ = nSbatch*n_S_Sbatch + (0:n_S_Sbatch-1);
S_ij_ = S_ij_(find(S_ij_<n_S)); n_S_sub = numel(S_ij_);
if (verbose>1); disp(sprintf(' %% nSbatch %d/%d S_ij_ %d-->%d',nSbatch,n_Sbatch,S_ij_(1+0),S_ij_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_t = tic();
svd_SVUXM_lwSM____ = permute(reshape(VS_k_q__(1+S_ij_,:)*SS_k_q__*reshape(svd_USVUXM_SMwl____,[n_S_rank,n_M_sub*n_w_max*FTK.n_svd_l]),[n_S_sub,n_M_sub,n_w_max,FTK.n_svd_l]),[4,3,1,2]);
svd_USESVUXM_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwSM____,[FTK.n_svd_l,n_w_max*n_S_sub*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S_sub,n_M_sub]);
l2_dSM___ = permute(reshape(reshape(sqrt(UX_S_l2_(1+S_ij_)),[n_S_sub,1])*reshape(sqrt(UX_M_l2_dM__(:,1+M_ij_)),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
n2_dSM___ = 1./max(1e-14,l2_dSM___);
f2_dSM___ = permute(reshape(reshape(sqrt(UX_S_l2_(1+S_ij_)),[n_S_sub,1])*reshape(1./max(1e-14,sqrt(UX_M_l2_dM__(:,1+M_ij_))),[1,FTK.n_delta_v*n_M_sub]),[n_S_sub,FTK.n_delta_v,n_M_sub]),[2,1,3]);
ss_S_ = reshape(UX_S_l2_(1+S_ij_),[n_S_sub,1]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USESVUXM_dwSM____: %0.6f',tmp_t)); end;
if (nMbatch==0 && nSbatch==0 && verbose>0); 
tmp_str = 'svd_VUXM_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USVUXM_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_lwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'f2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = tic();
X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____; %<-- correlation. ;
I_value_sub_dwSM____ = repmat(reshape(f2_dSM___,[FTK.n_delta_v,1,n_S_sub,n_M_sub]),[1,n_w_max,1,1]) .* X_sub_dwSM____; %<-- I_value. ;
I_value_use_dwSM____ = max(0,real(I_value_sub_dwSM____));
% L = -I*I*<M,M> + 2*I*<M,S> - <S,S> ;
L_sub_dwSM____ = -abs(I_value_use_dwSM____).^2.*repmat(reshape(real(UX_M_l2_dM__(:,1+M_ij_)),[FTK.n_delta_v,1,1,n_M_sub]),[1,n_w_max,n_S_sub,1]) + 2*I_value_use_dwSM____.*real(svd_USESVUXM_dwSM____) - repmat(reshape(real(ss_S_),[1,1,n_S_sub,1]),[FTK.n_delta_v,n_w_max,1,n_M_sub]); %<-- log-likelihood. ;
if (nMbatch==0 && nSbatch==0 && verbose>0); 
tmp_str = 'L_sub_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% L_sub_dwSM____: %0.6f',tmp_t)); end;
tmp_t = tic();
[tmp_L_,tmp_delta_j_] = max(reshape(L_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]),[],1); %<-- maximize log-likelihood. ;
tmp_delta_j_ = tmp_delta_j_ - 1; assert(min(tmp_delta_j_)>=0); assert(max(tmp_delta_j_)<FTK.n_delta_v);
tmp_I_value_sub_dwSM__ = reshape(I_value_use_dwSM____,[FTK.n_delta_v,n_w_max*n_S_sub*n_M_sub]);
tmp_I_value_sub_wSM_ = zeros(n_w_max*n_S_sub*n_M_sub,1);
tmp_t2=tic();
for nl=0:n_w_max*n_S_sub*n_M_sub-1;
tmp_I_value_sub_wSM_(1+nl) = tmp_I_value_sub_dwSM__(1+tmp_delta_j_(1+nl),1+nl);
end;%for nl=0:n_w_max*n_S_sub*n_M_sub-1;
tmp_t2 = toc(tmp_t2); if (verbose>1); disp(sprintf(' %% tmp_I_value_sub_wSM_ %0.6fs',tmp_t2)); end;
L_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_L_,[n_w_max,n_S_sub,n_M_sub]);
delta_j_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_delta_j_,[n_w_max,n_S_sub,n_M_sub]);
I_value_wSM___(:,1+S_ij_,1+M_ij_) = reshape(tmp_I_value_sub_wSM_,[n_w_max,n_S_sub,n_M_sub]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% L_wSM___: %0.6f',tmp_t)); end;
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
end;%if (n_M_sub>0);
end;%for nMbatch=0:n_Mbatch-1;
if (verbose>1); disp(sprintf(' %% [finished ampmh_X_wSM___3]')); end;


