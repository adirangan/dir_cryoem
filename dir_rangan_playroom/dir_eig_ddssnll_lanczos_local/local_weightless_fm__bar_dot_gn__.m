function [tmp_fg_mn__] = local_weightless_fm__bar_dot_gn__(n_k_p_r,l_max_,n_M,tmp_f_ykabcm__,tmp_g_ykabcn__);
n_m = size(tmp_f_ykabcm__,2);
n_n = size(tmp_g_ykabcn__,2);
tmp_fg_mn__ = zeros(n_m,n_n);
for nm=0:n_m-1;
for nn=0:n_n-1;
[tmp_fg_mn__(1+nm,1+nn)] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M,tmp_f_ykabcm__(:,1+nm),tmp_g_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
end;%for nm=0:n_m-1;
