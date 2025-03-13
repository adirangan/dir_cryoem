function [tmp_f_qkabcn__] = local_qkabc_imagecount_normalize_fn__(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabcn__);
n_n = size(tmp_f_qkabcn__,2);
for nn=0:n_n-1;
tmp_f_qkabcn__(:,1+nn) = local_qkabc_imagecount_normalize_f_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabcn__(:,1+nn));
end;%for nn=0:n_n-1;
