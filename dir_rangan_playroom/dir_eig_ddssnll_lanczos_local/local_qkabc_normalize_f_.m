function [tmp_f_ykabc_] = local_qkabc_normalize_f_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_ykabc_);
tmp_ff = local_qkabc_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_ykabc_,tmp_f_ykabc_);
tmp_f = sqrt(tmp_ff);
tmp_f_ykabc_ = tmp_f_ykabc_/max(1e-12,tmp_f);
