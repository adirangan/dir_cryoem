function M_resampled_k_p_wkM__ = ...
interp_p_to_p_resample_( ...
 n_k_p_r ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);

n_w_ = reshape(n_w_(1:n_k_p_r),[n_k_p_r,1]);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

n_w_resampled_ = n_w_max*ones(n_k_p_r,1);
n_w_resampled_sum = sum(n_w_resampled_);

for nM=0:n_M-1;
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_resampled_k_q_wk__ = zeros(n_w_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
n_w_2 = round(n_w/2);
index_nw_0in_head_ = n_w_csum_(1+nk_p_r) + [0:n_w_2-1];
index_nw_out_head_ = 0:n_w_2-1;
index_nw_0in_tail_ = n_w_csum_(1+nk_p_r) + n_w + [-n_w_2+1:-1];
index_nw_out_tail_ = n_w_max + [-n_w_2+1:-1];
M_resampled_k_q_wk__(1+index_nw_out_head_,1+nk_p_r) = M_k_q_wk_(1+index_nw_0in_head_);
M_resampled_k_q_wk__(1+index_nw_out_tail_,1+nk_p_r) = M_k_q_wk_(1+index_nw_0in_tail_);
end;%for nk_p_r=0:n_k_p_r-1;
M_resampled_k_q_wk_ = reshape(M_resampled_k_q_wk__,[n_w_resampled_sum,1]);
M_resampled_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_resampled_,n_w_resampled_sum,M_resampled_k_q_wk_);
M_resampled_k_p_wkM__(:,1+nM) = M_resampled_k_p_wk_;
end;%for nM=0:n_M-1;


