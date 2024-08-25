function [tmp_f_ykabc_] = local_imagecount_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,weight_imagecount_M_,tmp_f_ykabc_);
tmp_ff = local_imagecount_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,weight_imagecount_M_,tmp_f_ykabc_,tmp_f_ykabc_);
tmp_f = sqrt(tmp_ff);
tmp_f_ykabc_ = tmp_f_ykabc_/max(1e-12,tmp_f);
