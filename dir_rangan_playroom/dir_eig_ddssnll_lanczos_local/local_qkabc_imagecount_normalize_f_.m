function [tmp_f_qkabc_] = local_qkabc_imagecount_normalize_f_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabc_);
tmp_ff = local_qkabc_imagecount_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,weight_imagecount_M_,tmp_f_qkabc_,tmp_f_qkabc_);
tmp_f = sqrt(tmp_ff);
tmp_f_qkabc_ = tmp_f_qkabc_/max(1e-12,tmp_f);
