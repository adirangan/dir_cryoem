function [tmp_g_ykabc_] = local_orthogonalcomplement_gperpf(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_f_ykabc_ = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_);
tmp_fg = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_g_ykabc_ = tmp_g_ykabc_ - tmp_f_ykabc_ * tmp_fg;
