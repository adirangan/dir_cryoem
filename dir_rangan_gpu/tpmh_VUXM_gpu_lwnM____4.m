function VUXM_gpu_lwnM____ = tpmh_VUXM_gpu_lwnM____4(FTK,n_k_p_r,n_w_,n_M,M_k_q_gpu_wrM__,n_UX_rank,UX_gpu_rr__,X_weight_gpu_r_);

flag_verbose=0;

f_zero = gpuArray( single(0.0));
if ~strcmp(class(M_k_q_gpu_wrM__),'gpuArray'); tmp_t = tic(); M_k_q_gpu_wrM__ = gpuArray( (M_k_q_gpu_wrM__)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_gpu_wrM__: time %0.6fs',tmp_t)); end; end;
if ~strcmp(class(UX_gpu_rr__),'gpuArray'); tmp_t = tic(); UX_gpu_rr__ = gpuArray( (UX_gpu_rr__)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% UX_gpu_rr__: time %0.6fs',tmp_t)); end; end;
if ~strcmp(class(X_weight_gpu_r_),'gpuArray'); tmp_t = tic(); X_weight_gpu_r_ = gpuArray( (X_weight_gpu_r_)); tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% X_weight_gpu_r_: time %0.6fs',tmp_t)); end; end;

n_w_max = max(n_w_); n_w_2 = round(n_w_max/2);
l_max = max(abs(FTK.svd_l_));
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
V_gpu_lr__ = gpuArray( (reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r])));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% V_gpu_lr__: time %0.6fs',tmp_t)); end;
tmp_t = tic();
V_UX_gpu_lrn___ = zeros(FTK.n_svd_l,n_k_p_r,n_UX_rank,'like',f_zero);
for nUX_rank=0:n_UX_rank-1;
V_UX_gpu_lrn___(:,:,1+nUX_rank) = V_gpu_lr__*diag(UX_gpu_rr__(:,1+nUX_rank).*X_weight_gpu_r_(:));
end;%for nUX_rank=0:n_UX_rank-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% V_UX_gpu_lrn___: time %0.6fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
assert(mod(n_w_max,2)==0); %<-- assert that n_w_max is even. ;
index_nw_zerobased_from_centered_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place n_w_2 mode first, ;
index_nw_centered_from_zerobased_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place first mode at n_w_2. ;
%%%%;
index_nw_centered_padded_out_start_ = zeros(1+2*l_max,1);
index_nw_centered_padded_out_final_ = zeros(1+2*l_max,1);
index_nw_centered_padded_0in_start_ = zeros(1+2*l_max,1);
index_nw_centered_padded_0in_final_ = zeros(1+2*l_max,1);
for l_shift=-l_max:+l_max;
index_nw_centered_padded_out_start_(1+l_max+l_shift) = 1+l_max;
index_nw_centered_padded_out_final_(1+l_max+l_shift) = n_w_max-1+l_max;
index_nw_centered_padded_0in_start_(1+l_max+l_shift) = 1+l_max+l_shift;
index_nw_centered_padded_0in_final_(1+l_max+l_shift) = n_w_max-1+l_max+l_shift;
end;%for l_shift=-l_max:+l_max;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% index_nw_centered_padded_xxx_xxxxx_: time %0.6fs',tmp_t)); end;
%%%%%%%%;

tmp_t = tic();
V_UX_gpu_nrl___ = permute(V_UX_gpu_lrn___,[3,2,1]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% V_UX_gpu_rnl___: time %0.6fs',tmp_t)); end;
tmp_t = tic();
M_k_q_centered_padded_gpu_rwM___ = cat(1+0,zeros(l_max,n_k_p_r,n_M),circshift(reshape(M_k_q_gpu_wrM__,[n_w_max,n_k_p_r,n_M]),-n_w_2,1),zeros(l_max,n_k_p_r,n_M));
M_k_q_centered_padded_gpu_rwM___(1+l_max,:,:) = 0;
M_k_q_centered_padded_gpu_rwM___ = permute(M_k_q_centered_padded_gpu_rwM___,[2,1,3]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_q_centered_padded_gpu_rwM___: time %0.6fs',tmp_t)); end;
tmp_t = tic();
VUXM_centered_gpu_nwMl____ = zeros(n_UX_rank,n_w_max,n_M,FTK.n_svd_l,'like',f_zero);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
index_nw_centered_padded_out_start = index_nw_centered_padded_out_start_(1+l_max+l_shift);
index_nw_centered_padded_out_final = index_nw_centered_padded_out_final_(1+l_max+l_shift);
index_nw_centered_padded_0in_start = index_nw_centered_padded_0in_start_(1+l_max+l_shift);
index_nw_centered_padded_0in_final = index_nw_centered_padded_0in_final_(1+l_max+l_shift);
index_nw_centered_padded_out_ = [index_nw_centered_padded_out_start:index_nw_centered_padded_out_final];
index_nw_centered_padded_0in_ = [index_nw_centered_padded_0in_start:index_nw_centered_padded_0in_final];
VUXM_centered_gpu_nwMl____(:,2:end,:,1+nl) = pagemtimes(reshape(V_UX_gpu_nrl___(:,:,1+nl),[n_UX_rank,n_k_p_r,1])/n_w_max,reshape(M_k_q_centered_padded_gpu_rwM___(:,1+index_nw_centered_padded_0in_,:),[n_k_p_r,n_w_max-1,n_M]))  ;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% VUXM_centered_padded_lwMn____: time %0.6fs',tmp_t)); end;
tmp_t = tic();
VUXM_gpu_lwnM____ = permute(VUXM_centered_gpu_nwMl____(:,1+index_nw_zerobased_from_centered_,:,:),[4,2,1,3]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% VUXM_gpu_lwnM____: time %0.6fs',tmp_t)); end;



