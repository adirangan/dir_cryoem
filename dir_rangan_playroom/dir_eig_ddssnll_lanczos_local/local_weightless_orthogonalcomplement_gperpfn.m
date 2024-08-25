function [tmp_g_ykabc_] = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,l_max_,n_M,tmp_f_ykabcn__,tmp_g_ykabc_);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_g_ykabc_ = local_weightless_orthogonalcomplement_gperpf(n_k_p_r,l_max_,n_M,tmp_f_ykabcn__(:,1+nn),tmp_g_ykabc_);
end;%for nn=0:n_n-1;
