function [tmp_g_ykabc_] = local_weightless_orthogonalcomplement_gperpf(n_k_p_r,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_f_ykabc_ = local_weightless_normalize_f_(n_k_p_r,l_max_,n_M,tmp_f_ykabc_);
tmp_fg = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_g_ykabc_ = tmp_g_ykabc_ - tmp_f_ykabc_ * tmp_fg;
