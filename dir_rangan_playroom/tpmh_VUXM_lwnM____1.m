function VUXM_lwnM____ = tpmh_VUXM_lwnM____1(FTK,n_k_p_r,n_w_,n_M,M_k_q__,n_UX_rank,UX__,X_weight_r_);
verbose=1;
n_w_max = max(n_w_);
l_max = max(abs(FTK.svd_l_));
tmp_t = tic();
V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% V_r__ %0.3fs',tmp_t)); end;

disp(sprintf(' %% strategy one: '));
tmp_t_A = tic();
tmp_t = tic();
M_k_q_rwlM____ = innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q__,l_max);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q_rwlM____ %0.3fs',tmp_t)); end;
%%%%%%%%;
VUXM_A_lwnM____ = zeros(FTK.n_svd_l,n_w_max,n_UX_rank,n_M);
tmp_t_VUXM_A_lwnM____ = 0;
for nM=0:n_M-1;
for nUX_rank=0:n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
tmp_t = tic();
VUXM_A_lwnM____(1+nl,:,1+nUX_rank,1+nM) = V_r__(1+nl,:)*diag(UX__(:,1+nUX_rank).*X_weight_r_(:))*M_k_q_rwlM____(:,:,1+l_max+FTK.svd_l_(1+nl),1+nM) / n_w_max ;
tmp_t = toc(tmp_t); tmp_t_VUXM_A_lwnM____ = tmp_t_VUXM_A_lwnM____ + tmp_t;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nM=0:n_M-1;
if (verbose); disp(sprintf(' %% VUXM_A_lwnM____ %0.3fs',tmp_t_VUXM_A_lwnM____)); end;
tmp_t_A = toc(tmp_t_A); disp(sprintf(' %% strategy one: %0.3fs',tmp_t_A));

disp(sprintf(' %% strategy two: '));
tmp_t_B = tic();
tmp_t = tic();
M_k_q_rwMl____ = permute(innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q__,l_max),[1,2,4,3]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q_rwMl____ %0.3fs',tmp_t)); end;
tmp_t_V_r_UX_X__ = 0;
tmp_t_VUXM_wMnl____ = 0;
for nUX_rank=0:n_UX_rank-1;
tmp_t = tic();
tmp_V_r_UX_X__ = V_r__*diag(UX__(:,1+nUX_rank).*X_weight_r_(:));
tmp_t = toc(tmp_t); tmp_t_V_r_UX_X__ = tmp_t_V_r_UX_X__ + tmp_t;
for nl=0:FTK.n_svd_l-1;
tmp_t = tic();
VUXM_wMnl____(:,:,1+nUX_rank,1+nl) = reshape(tmp_V_r_UX_X__(1+nl,:)*reshape(M_k_q_rwMl____(:,:,:,1+l_max+FTK.svd_l_(1+nl)),[n_k_p_r,n_w_max*n_M]),[n_w_max,n_M]) / n_w_max;
tmp_t = toc(tmp_t); tmp_t_VUXM_wMnl____ = tmp_t_VUXM_wMnl____ + tmp_t;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nUX_rank=0:n_UX_rank-1;
tmp_t = tic();
if (verbose); disp(sprintf(' %% V_r_UX_X__ %0.3fs',tmp_t_V_r_UX_X__)); end;
if (verbose); disp(sprintf(' %% VUXM_wMnl____ %0.3fs',tmp_t_VUXM_wMnl____)); end;
VUXM_B_lwnM____ = permute(VUXM_wMnl____,[4,1,3,2]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% VUXM_B_lwnM____ %0.3fs',tmp_t)); end;
tmp_t_B = toc(tmp_t_B); disp(sprintf(' %% strategy two: %0.3fs',tmp_t_B));

disp(sprintf(' %% A vs B error: %0.16f',fnorm(VUXM_A_lwnM____ - VUXM_B_lwnM____)/fnorm(VUXM_A_lwnM____)));

VUXM_lwnM____ = VUXM_A_lwnM____;

return;

