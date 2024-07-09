function [tmp_ykabcn__] = local_ykabcn__from_ykn__an__bn__cn__(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__,tmp_a_Mn__,tmp_b_Mn__,tmp_c_Mn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykabcn__,2);
tmp_ykabcn__ = zeros(n_lm_sum + 3*n_M,n_n);
for nn=0:n_n-1;
tmp_ykabcn__(:,1+nn) = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__(:,1+nn),tmp_a_Mn__(:,1+nn),tmp_b_Mn__(:,1+nn),tmp_c_Mn__(:,1+nn));
end;%for nn=0:n_n-1;
