function ...
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___0(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M...
,n_M_batch...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
);
%%%%%%%%;
% Calculates correlations between principal-templates and principal-images across translations. ;
% Batches images and uses svd for templates. ;
%%%%%%%%;
verbose=1;
if (verbose>0); disp(sprintf(' %% [entering ampmh_X_wSM___0]')); end;
n_w_max = max(n_w_);
X_wSM___ = zeros(n_w_max,n_S,n_M);
delta_j_wSM___ = zeros(n_w_max,n_S,n_M);
UX_US_CTF_k_q_nSw___ = permute(reshape(US_k_q__,[n_w_max,pm_n_UX_rank,n_S_rank]),[2,3,1]);
if (verbose>0); 
tmp_str = 'X_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'delta_j_wSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'UX_US_CTF_k_q_nSw___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>0); 
n_batch = ceil(n_M/n_M_batch);
if (verbose>0); disp(sprintf(' %% n_batch %d',n_batch)); end;
for nbatch=0:n_batch-1;
M_ij_ = nbatch*n_M_batch + (0:n_M_batch-1);
M_ij_ = M_ij_(find(M_ij_<n_M)); n_M_sub = numel(M_ij_);
if (verbose>1); disp(sprintf(' %% nbatch %d/%d M_ij_ %d-->%d',nbatch,n_batch,M_ij_(1+0),M_ij_(1+n_M_sub-1))); end;
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
svd_SVUXM_lwSM____ = permute(reshape(VS_k_q__*SS_k_q__*reshape(svd_USVUXM_SMwl____,[n_S_rank,n_M_sub*n_w_max*FTK.n_svd_l]),[n_S,n_M_sub,n_w_max,FTK.n_svd_l]),[4,3,1,2]);
svd_USESVUXM_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXM_lwSM____,[FTK.n_svd_l,n_w_max*n_S*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]);
l2_dSM___ = permute(reshape(reshape(sqrt(UX_S_l2_),[n_S,1])*reshape(sqrt(UX_M_l2_dM__(:,1+M_ij_)),[1,FTK.n_delta_v*n_M_sub]),[n_S,FTK.n_delta_v,n_M_sub]),[2,1,3]);
n2_dSM___ = 1./max(1e-12,l2_dSM___);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_USESVUXM_dwSM____: %0.6f',tmp_t)); end;
if (nbatch==0 && verbose>1); 
tmp_str = 'svd_VUXM_nMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USVUXM_SMwl____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_SVUXM_lwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'svd_USESVUXM_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'l2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
tmp_str = 'n2_dSM___'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = tic();
X_sub_dwSM____ = repmat(reshape(n2_dSM___,[FTK.n_delta_v,1,n_S,n_M_sub]),[1,n_w_max,1,1]) .* svd_USESVUXM_dwSM____;
if (nbatch==0 && verbose>1); 
tmp_str = 'X_sub_dwSM____'; disp(sprintf(' %% memory: %s --> %0.6f GB',tmp_str,whos(tmp_str).bytes/1e9));
end;%if (verbose>1); 
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_sub_dwSM____: %0.6f',tmp_t)); end;
tmp_t = tic();
[tmp_X_,tmp_delta_j_] = max(reshape(X_sub_dwSM____,[FTK.n_delta_v,n_w_max*n_S*n_M_sub]),[],1);
tmp_delta_j_ = tmp_delta_j_ - 1;
X_wSM___(:,:,1+M_ij_) = reshape(tmp_X_,[n_w_max,n_S,n_M_sub]);
delta_j_wSM___(:,:,1+M_ij_) = reshape(tmp_delta_j_,[n_w_max,n_S,n_M_sub]);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wSM___: %0.6f',tmp_t)); end;
end;%if (n_M_sub>0);
end;%for nbatch=0:n_batch-1;
if (verbose>0); disp(sprintf(' %% [finished ampmh_X_wSM___0]')); end;


