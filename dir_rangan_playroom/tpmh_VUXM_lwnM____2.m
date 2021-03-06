function VUXM_lwnM____ = tpmh_VUXM_lwnM____2(FTK,n_k_p_r,n_w_,n_M,M_k_q__,n_UX_rank,UX__,X_weight_r_);

verbose=1;
n_w_max = max(n_w_); n_w_2 = round(n_w_max/2);
n_op = n_M*FTK.n_svd_l*n_k_p_r*n_w_max*n_UX_rank;
l_max = max(abs(FTK.svd_l_));
%%%%%%%%;

tmp_t_0 = tic();
VUXM_0_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M,M_k_q__,n_UX_rank,UX__,X_weight_r_);
tmp_t_0 = toc(tmp_t_0); disp(sprintf(' %% strategy_0: %0.6fs',tmp_t_0));
disp(sprintf(' %% n_op %d, Gops %0.2f',n_op,n_op/tmp_t_0/1e9));

%%%%%%%%;
tmp_t = tic();
V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% V_r__ %0.6fs',tmp_t)); end;
V_UX_lrn___ = zeros(FTK.n_svd_l,n_k_p_r,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
V_UX_lrn___(:,:,1+nUX_rank) = V_r__*diag(UX__(:,1+nUX_rank).*X_weight_r_(:));
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
index_nw_out__ = cell(1+2*l_max,1);
index_nw_0in__ = cell(1+2*l_max,1);
n_nw_ = zeros(1+2*l_max,1);
for l_shift=-l_max:+l_max;
index_nw_out_head_ = (0:n_w_2-1-max(0,+l_shift));
index_nw_0in_head_ = [ n_w_max+l_shift:n_w_max-1 , max(0,+l_shift):n_w_2-1-max(0,-l_shift) ];
index_nw_out_tail_ = (n_w_2+1+max(0,-l_shift):n_w_max-1);
index_nw_0in_tail_ = [ n_w_2+1+max(0,+l_shift):n_w_max-1-max(0,-l_shift), 0:+l_shift-1 ];
index_nw_out_ = [ index_nw_out_head_ , index_nw_out_tail_ ];
index_nw_0in_ = [ index_nw_0in_head_ , index_nw_0in_tail_ ];
n_nw = numel(index_nw_out_);
n_nw_(1+l_max+l_shift) = n_nw;
index_nw_out__{1+l_max+l_shift} = index_nw_out_;
index_nw_0in__{1+l_max+l_shift} = index_nw_0in_;
end;%for l_shift=-l_max:+l_max;
%%%%%%%%;
assert(mod(n_w_max,2)==0); %<-- assert that n_w_max is even. ;
index_nw_zerobased_from_centered_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place n_w_2 mode first, ;
index_nw_centered_from_zerobased_ = [ n_w_2:n_w_max-1 , 0:n_w_2-1 ]; %<-- note that we place first mode at n_w_2. ;
index_nw_centered_out_start_ = zeros(1+2*l_max,1);
index_nw_centered_out_final_ = zeros(1+2*l_max,1);
index_nw_centered_0in_start_ = zeros(1+2*l_max,1);
index_nw_centered_0in_final_ = zeros(1+2*l_max,1);
for l_shift=-l_max:+l_max;
index_nw_centered_out_start_(1+l_max+l_shift) = 1+max(0,-l_shift);
index_nw_centered_out_final_(1+l_max+l_shift) = n_w_max-1+min(0,-l_shift);
index_nw_centered_0in_start_(1+l_max+l_shift) = 1+max(0,+l_shift);
index_nw_centered_0in_final_(1+l_max+l_shift) = n_w_max-1+min(0,+l_shift);
end;%for l_shift=-l_max:+l_max;

tmp_t_2 = tic();
M_k_q_rwM___ = permute(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),[2,1,3]);
%%%%%%%%;
VUXM_2_lwnM____ = zeros(FTK.n_svd_l,n_w_max,n_UX_rank,n_M);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
n_nw = n_nw_(1+l_max+l_shift);
index_nw_out_ = index_nw_out__{1+l_max+l_shift};
index_nw_0in_ = index_nw_0in__{1+l_max+l_shift};
for nUX_rank=0:n_UX_rank-1;
VUXM_2_lwnM____(1+nl,1+index_nw_out_,1+nUX_rank,:) = reshape(V_UX_lrn___(1+nl,:,1+nUX_rank)*reshape(M_k_q_rwM___(:,1+index_nw_0in_,:),[n_k_p_r,n_nw*n_M]),[1,n_nw,1,n_M]) / n_w_max ;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t_2 = toc(tmp_t_2); disp(sprintf(' %% strategy_2: %0.6fs',tmp_t_2));
disp(sprintf(' %% n_op %d, Gops %0.2f',n_op,n_op/tmp_t_2/1e9));

tmp_t_3 = tic();
M_k_q_rwM___ = permute(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),[2,1,3]);
%%%%%%%%;
VUXM_3_lwMn____ = zeros(FTK.n_svd_l,n_w_max,n_M,n_UX_rank);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
n_nw = n_nw_(1+l_max+l_shift);
index_nw_out_ = index_nw_out__{1+l_max+l_shift};
index_nw_0in_ = index_nw_0in__{1+l_max+l_shift};
for nUX_rank=0:n_UX_rank-1;
VUXM_3_lwMn____(1+nl,1+index_nw_out_,:,1+nUX_rank) = reshape(V_UX_lrn___(1+nl,:,1+nUX_rank)*reshape(M_k_q_rwM___(:,1+index_nw_0in_,:),[n_k_p_r,n_nw*n_M]),[1,n_nw,n_M,1]) / n_w_max ;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nl=0:FTK.n_svd_l-1;
VUXM_3_lwnM____ = permute(VUXM_3_lwMn____,[1,2,4,3]);
tmp_t_3 = toc(tmp_t_3); disp(sprintf(' %% strategy_3: %0.6fs',tmp_t_3));
disp(sprintf(' %% n_op %d, Gops %0.2f',n_op,n_op/tmp_t_3/1e9));

tmp_t_4 = tic();
tmp_t = tic();
M_k_q_centered_rwM___ = permute(reshape(M_k_q__,[n_w_max,n_k_p_r,n_M]),[2,1,3]);
M_k_q_centered_rwM___ = M_k_q_centered_rwM___(:,1+index_nw_centered_from_zerobased_,:);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q_centered_rwM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
VUXM_centered_4_lwMn____ = zeros(FTK.n_svd_l,n_w_max,n_M,n_UX_rank);
tmp_t = tic();
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
index_nw_centered_out_start = index_nw_centered_out_start_(1+l_max+l_shift);
index_nw_centered_out_final = index_nw_centered_out_final_(1+l_max+l_shift);
index_nw_centered_0in_start = index_nw_centered_0in_start_(1+l_max+l_shift);
index_nw_centered_0in_final = index_nw_centered_0in_final_(1+l_max+l_shift);
index_nw_centered_out_ = [index_nw_centered_out_start:index_nw_centered_out_final];
index_nw_centered_0in_ = [index_nw_centered_0in_start:index_nw_centered_0in_final];
n_nw = numel(index_nw_centered_out_);
for nUX_rank=0:n_UX_rank-1;
for nM=0:n_M-1;
VUXM_centered_4_lwMn____(1+nl,1+index_nw_centered_out_,1+nM,1+nUX_rank) = V_UX_lrn___(1+nl,:,1+nUX_rank)*M_k_q_centered_rwM___(:,1+index_nw_centered_0in_,1+nM) / n_w_max ;
end;%for nM=0:n_M-1;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% VUXM_centered_4_lwnM___: %0.6fs',tmp_t)); end;
tmp_t = tic();
VUXM_4_lwnM____ = permute(VUXM_centered_4_lwMn____(:,1+index_nw_zerobased_from_centered_,:,:),[1,2,4,3]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% VUXM_4_lwnM___: %0.6fs',tmp_t)); end;
tmp_t_4 = toc(tmp_t_4); disp(sprintf(' %% strategy_4: %0.6fs',tmp_t_4));
disp(sprintf(' %% n_op %d, Gops %0.2f',n_op,n_op/tmp_t_4/1e9));

disp(sprintf(' %% 0 vs 2 error: %0.16f',fnorm(VUXM_0_lwnM____ - VUXM_2_lwnM____)/fnorm(VUXM_0_lwnM____)));
disp(sprintf(' %% 0 vs 3 error: %0.16f',fnorm(VUXM_0_lwnM____ - VUXM_3_lwnM____)/fnorm(VUXM_0_lwnM____)));
disp(sprintf(' %% 0 vs 4 error: %0.16f',fnorm(VUXM_0_lwnM____ - VUXM_4_lwnM____)/fnorm(VUXM_0_lwnM____)));

VUXM_lwnM____ = VUXM_0_lwnM____;

return;

