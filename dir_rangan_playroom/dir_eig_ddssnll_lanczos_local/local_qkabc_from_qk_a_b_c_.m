function [tmp_qkabc_] = local_qkabc_from_qk_a_b_c_(n_q,n_k_p_r,n_M,tmp_dvol_qk_,tmp_a_M_,tmp_b_M_,tmp_c_M_);
n_qk = n_q*n_k_p_r;
tmp_qkabc_ = zeros(n_qk + 3*n_M,1);
tmp_qkabc_(1:n_qk) = tmp_dvol_qk_(:);
tmp_qkabc_(1*n_qk + 0*n_M + [1:n_M]) = tmp_a_M_(:);
tmp_qkabc_(1*n_qk + 1*n_M + [1:n_M]) = tmp_b_M_(:);
tmp_qkabc_(1*n_qk + 2*n_M + [1:n_M]) = tmp_c_M_(:);
