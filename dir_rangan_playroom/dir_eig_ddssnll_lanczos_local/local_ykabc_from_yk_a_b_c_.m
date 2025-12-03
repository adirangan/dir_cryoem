function [tmp_ykabc_] = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_);
n_y_ = (l_max_+1).^2; n_y_max = max(n_y_); n_y_sum = sum(n_y_); n_y_csum_ = cumsum([0;n_y_]); l_max_max = max(l_max_);
tmp_ykabc_ = zeros(n_y_sum + 3*n_M,1);
tmp_ykabc_(1:n_y_sum) = tmp_dvol_yk_(:);
tmp_ykabc_(1*n_y_sum + 0*n_M + [1:n_M]) = tmp_a_M_(:);
tmp_ykabc_(1*n_y_sum + 1*n_M + [1:n_M]) = tmp_b_M_(:);
tmp_ykabc_(1*n_y_sum + 2*n_M + [1:n_M]) = tmp_c_M_(:);
