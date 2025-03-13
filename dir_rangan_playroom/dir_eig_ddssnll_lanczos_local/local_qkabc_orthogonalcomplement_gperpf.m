function [tmp_g_qkabc_] = local_qkabc_orthogonalcomplement_gperpf(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabc_,tmp_g_qkabc_);
tmp_f_qkabc_ = local_qkabc_normalize_f_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabc_);
tmp_fg = local_qkabc_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabc_,tmp_g_qkabc_);
tmp_g_qkabc_ = tmp_g_qkabc_ - tmp_f_qkabc_ * tmp_fg;
