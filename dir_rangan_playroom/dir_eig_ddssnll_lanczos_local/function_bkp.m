
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_dvol_yk_ = zeros(n_lm_sum,1);
tmp_a_M_ = zeros(n_M,1);
tmp_b_M_ = zeros(n_M,1);
tmp_c_M_ = zeros(n_M,1);
tmp_dvol_yk_(:) = tmp_ykabc_(1:n_lm_sum);
tmp_a_M_(:) = tmp_ykabc_(1*n_lm_sum + 0*n_M + [1:n_M]);
tmp_b_M_(:) = tmp_ykabc_(1*n_lm_sum + 1*n_M + [1:n_M]);
tmp_c_M_(:) = tmp_ykabc_(1*n_lm_sum + 2*n_M + [1:n_M]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykabc_] = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_ykabc_ = zeros(n_lm_sum + 3*n_M,1);
tmp_ykabc_(1:n_lm_sum) = tmp_dvol_yk_(:);
tmp_ykabc_(1*n_lm_sum + 0*n_M + [1:n_M]) = tmp_a_M_(:);
tmp_ykabc_(1*n_lm_sum + 1*n_M + [1:n_M]) = tmp_b_M_(:);
tmp_ykabc_(1*n_lm_sum + 2*n_M + [1:n_M]) = tmp_c_M_(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_dvol_ykn__,tmp_a_Mn__,tmp_b_Mn__,tmp_c_Mn__] = local_ykn_an__bn__cn__from_ykabcn__(n_k_p_r,l_max_,n_M,tmp_ykabcn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykabcn__,2);
tmp_dvol_ykn__ = zeros(n_lm_sum,n_n);
tmp_a_Mn__ = zeros(n_M,n_n);
tmp_b_Mn__ = zeros(n_M,n_n);
tmp_c_Mn__ = zeros(n_M,n_n);
for nn=0:n_n-1;
[tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_ykabc_);
tmp_dvol_ykn__(:,1+nn) = tmp_dvol_yk_;
tmp_a_Mn__(:,1+nn) = tmp_a_M_;
tmp_b_Mn__(:,1+nn) = tmp_b_M_;
tmp_c_Mn__(:,1+nn) = tmp_c_M_;
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykabcn__] = local_ykabcn__from_ykn__an__bn__cn__(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__,tmp_a_Mn__,tmp_b_Mn__,tmp_c_Mn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykabcn__,2);
tmp_ykabcn__ = zeros(n_lm_sum + 3*n_M,n_n);
for nn=0:n_n-1;
tmp_ykabcn__(:,1+nn) = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__(:,1+nn),tmp_a_Mn__(:,1+nn),tmp_b_Mn__(:,1+nn),tmp_c_Mn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_yk__] = local_yk__from_yk_(n_k_p_r,l_max_,tmp_yk_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk__(1:n_lm,1+nk_p_r) = tmp_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_yk_] = local_yk_from_yk__(n_k_p_r,l_max_,tmp_yk__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk_(1+tmp_index_) = tmp_yk__(1:n_lm,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykn___] = local_ykn___from_ykn__(n_k_p_r,l_max_,tmp_ykn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykn__,2);
tmp_ykn___ = zeros(n_lm_max,n_k_p_r,n_n);
for nn=0:n_n-1;
tmp_ykn___(:,:,1+nn) = local_yk__from_yk_(n_k_p_r,l_max_,tmp_ykn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykn__] = local_ykn__from_ykn___(n_k_p_r,l_max_,tmp_ykn___);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykn__,2);
tmp_ykn__ = zeros(n_lm_sum,n_n);
for nn=0:n_n-1;
tmp_ykn__(:,1+nn) = local_yk_from_yk__(n_k_p_r,l_max_,tmp_ykn___(:,:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_fg] = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
[tmp_f_dvol_yk_,tmp_f_a_M_,tmp_f_b_M_,tmp_f_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_f_ykabc_);
[tmp_g_dvol_yk_,tmp_g_a_M_,tmp_g_b_M_,tmp_g_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_g_ykabc_);
tmp_f_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_f_dvol_yk_);
tmp_g_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_g_dvol_yk_);
tmp_fg_dvol = sum( (conj(tmp_f_dvol_yk__) .* tmp_g_dvol_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) , 'all' ) / scaling_volumetric ;
tmp_fg_a = sum(conj(tmp_f_a_M_) .* tmp_g_a_M_, 'all');
tmp_fg_b = sum(conj(tmp_f_b_M_) .* tmp_g_b_M_, 'all');
tmp_fg_c = sum(conj(tmp_f_c_M_) .* tmp_g_c_M_, 'all');
tmp_fg = tmp_fg_dvol + tmp_fg_a + tmp_fg_b + tmp_fg_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_fg_mn__] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcm__,tmp_g_ykabcn__);
n_m = size(tmp_f_ykabcm__,2);
n_n = size(tmp_g_ykabcn__,2);
tmp_fg_mn__ = zeros(n_m,n_n);
for nm=0:n_m-1;
for nn=0:n_n-1;
[tmp_fg_mn__(1+nm,1+nn)] = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcm__(:,1+nm),tmp_g_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
end;%for nm=0:n_m-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_f_ykabc_] = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_);
tmp_ff = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_f_ykabc_);
tmp_f = sqrt(tmp_ff);
tmp_f_ykabc_ = tmp_f_ykabc_/max(1e-12,tmp_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_f_ykabcn__] = local_normalize_fn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_f_ykabcn__(:,1+nn) = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_g_ykabc_] = local_orthogonalcomplement_gperpf(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_f_ykabc_ = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_);
tmp_fg = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_g_ykabc_ = tmp_g_ykabc_ - tmp_f_ykabc_ * tmp_fg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_g_ykabc_] = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__,tmp_g_ykabc_);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_g_ykabc_ = local_orthogonalcomplement_gperpf(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__(:,1+nn),tmp_g_ykabc_);
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
