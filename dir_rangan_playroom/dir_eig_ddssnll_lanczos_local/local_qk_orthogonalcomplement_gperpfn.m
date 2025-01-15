function [tmp_g_qkabc_] = local_qk_orthogonalcomplement_gperpfn(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabcn__,tmp_g_qkabc_);
n_n = size(tmp_f_qkabcn__,2);
for nn=0:n_n-1;
tmp_g_qkabc_ = local_qk_orthogonalcomplement_gperpf(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabcn__(:,1+nn),tmp_g_qkabc_);
end;%for nn=0:n_n-1;
