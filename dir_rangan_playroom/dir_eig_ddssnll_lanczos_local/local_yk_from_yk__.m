function [tmp_yk_] = local_yk_from_yk__(n_k_p_r,l_max_,tmp_yk__);
n_y_ = (l_max_+1).^2; n_y_max = max(n_y_); n_y_sum = sum(n_y_); n_y_csum_ = cumsum([0;n_y_]); l_max_max = max(l_max_);
tmp_yk_ = zeros(n_y_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_y = n_y_(1+nk_p_r);
tmp_index_ = n_y_csum_(1+nk_p_r) + (0:n_y-1);
tmp_yk_(1+tmp_index_) = tmp_yk__(1:n_y,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
