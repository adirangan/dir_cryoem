function [tmp_fg_mn__] = local_qkabc_imagecount_fm__bar_dot_gn__(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabcm__,tmp_g_qkabcn__);
n_m = size(tmp_f_qkabcm__,2);
n_n = size(tmp_g_qkabcn__,2);
tmp_fg_mn__ = zeros(n_m,n_n);
for nm=0:n_m-1;
for nn=0:n_n-1;
[tmp_fg_mn__(1+nm,1+nn)] = local_qkabc_imagecount_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabcm__(:,1+nm),tmp_g_qkabcn__(:,1+nn));
end;%for nn=0:n_n-1;
end;%for nm=0:n_m-1;

