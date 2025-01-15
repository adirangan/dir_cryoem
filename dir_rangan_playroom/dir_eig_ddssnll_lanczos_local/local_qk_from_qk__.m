function [tmp_qk_] = local_qk_from_qk__(n_q,n_k_p_r,tmp_qk__);
n_qk = n_q*n_k_p_r;
tmp_qk_ = zeros(n_qk,1);
tmp_qk_ = reshape(tmp_qk__,[n_qk,1]);
