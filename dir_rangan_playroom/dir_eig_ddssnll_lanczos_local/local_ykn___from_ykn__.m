function [tmp_ykn___] = local_ykn___from_ykn__(n_k_p_r,l_max_,tmp_ykn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykn__,2);
tmp_ykn___ = zeros(n_lm_max,n_k_p_r,n_n);
for nn=0:n_n-1;
tmp_ykn___(:,:,1+nn) = local_yk__from_yk_(n_k_p_r,l_max_,tmp_ykn__(:,1+nn));
end;%for nn=0:n_n-1;
