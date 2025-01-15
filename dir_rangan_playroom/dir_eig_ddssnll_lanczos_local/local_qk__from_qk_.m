function [tmp_qk__] = local_qk__from_qk_(n_q,n_k_p_r,tmp_qk_);
n_qk = n_q*n_k_p_r;
tmp_qk__ = zeros(n_q,n_k_p_r);
tmp_qk__ = reshape(tmp_qk_,[n_q,n_k_p_r]);

