function [tmp_f_ykabcn__] = local_weightless_normalize_fn__(n_k_p_r,l_max_,n_M,tmp_f_ykabcn__);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_f_ykabcn__(:,1+nn) = local_weightless_normalize_f_(n_k_p_r,l_max_,n_M,tmp_f_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
