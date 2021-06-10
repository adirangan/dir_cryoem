function VUXM_lw__ = tpmh_VUXM_lw__1(FTK,n_k_p_r,n_w_,M_k_q_,UX_,X_weight_r_);

tmp_t_0 = tic();
VUXM_0_lw__ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,1,M_k_q_,1,UX_,X_weight_r_);
tmp_t_0 = toc(tmp_t_0); disp(sprintf(' %% strategy 0: %0.6fs',tmp_t_0));

verbose=1;
n_w_max = max(n_w_); n_w_2 = round(n_w_max/2);
l_max = max(abs(FTK.svd_l_));
tmp_t = tic();
V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% V_r__ %0.6fs',tmp_t)); end;
V_UX_lr__ = V_r__*diag(UX_(:).*X_weight_r_(:));

tmp_t_1 = tic();
M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,1,M_k_q_,l_max);
%%%%%%%%;
VUXM_1_lw__ = zeros(FTK.n_svd_l,n_w_max);
for nl=0:FTK.n_svd_l-1;
VUXM_1_lw__(1+nl,:) = V_UX_lr__(1+nl,:)*M_k_q_rwl___(:,:,1+l_max+FTK.svd_l_(1+nl)) / n_w_max ;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t_1 = toc(tmp_t_1); disp(sprintf(' %% strategy 1: %0.6fs',tmp_t_1));

tmp_t_2 = tic();
M_k_q_rw__ = transpose(reshape(M_k_q_,[n_w_max,n_k_p_r]));
%%%%%%%%;
VUXM_2_lw__ = zeros(FTK.n_svd_l,n_w_max);
for nl=0:FTK.n_svd_l-1;
l_shift = FTK.svd_l_(1+nl);
index_nw_out_head_ = (0:n_w_2-1-max(0,+l_shift));
index_nw_0in_head_ = [ n_w_max+l_shift:n_w_max-1 , max(0,+l_shift):n_w_2-1-max(0,-l_shift) ];
%disp(sprintf(' %% head %d %d',length(index_nw_out_head_),length(index_nw_0in_head_)));
index_nw_out_tail_ = (n_w_2+1+max(0,-l_shift):n_w_max-1);
index_nw_0in_tail_ = [ n_w_2+1+max(0,+l_shift):n_w_max-1-max(0,-l_shift), 0:+l_shift-1 ];
%disp(sprintf(' %% tail %d %d',length(index_nw_out_tail_),length(index_nw_0in_tail_)));
index_nw_out_ = [ index_nw_out_head_ , index_nw_out_tail_ ];
index_nw_0in_ = [ index_nw_0in_head_ , index_nw_0in_tail_ ];
%disp(sprintf(' %% l_shift %+d',l_shift));
%disp(sprintf(' %% out: %s',num2str(index_nw_out_)))
%disp(sprintf(' %% 0in: %s',num2str(index_nw_0in_)))
VUXM_2_lw__(1+nl,1+index_nw_out_) = V_UX_lr__(1+nl,:)*M_k_q_rw__(:,1+index_nw_0in_) / n_w_max ;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t_2 = toc(tmp_t_2); disp(sprintf(' %% strategy 2: %0.6fs',tmp_t_2));

disp(sprintf(' %% 0 vs 1 error: %0.16f',fnorm(VUXM_0_lw__ - VUXM_1_lw__)/fnorm(VUXM_0_lw__)));
disp(sprintf(' %% 0 vs 2 error: %0.16f',fnorm(VUXM_0_lw__ - VUXM_2_lw__)/fnorm(VUXM_0_lw__)));

VUXM_lw__ = VUXM_0_lw__;

