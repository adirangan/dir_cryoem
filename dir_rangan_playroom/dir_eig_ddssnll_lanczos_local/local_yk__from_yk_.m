function [tmp_yk__] = local_yk__from_yk_(n_k_p_r,l_max_,tmp_yk_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk__(1:n_lm,1+nk_p_r) = tmp_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
