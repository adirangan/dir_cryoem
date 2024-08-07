function [tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_dvol_yk_ = zeros(n_lm_sum,1);
tmp_a_M_ = zeros(n_M,1);
tmp_b_M_ = zeros(n_M,1);
tmp_c_M_ = zeros(n_M,1);
tmp_dvol_yk_(:) = tmp_ykabc_(1:n_lm_sum);
tmp_a_M_(:) = tmp_ykabc_(1*n_lm_sum + 0*n_M + [1:n_M]);
tmp_b_M_(:) = tmp_ykabc_(1*n_lm_sum + 1*n_M + [1:n_M]);
tmp_c_M_(:) = tmp_ykabc_(1*n_lm_sum + 2*n_M + [1:n_M]);
