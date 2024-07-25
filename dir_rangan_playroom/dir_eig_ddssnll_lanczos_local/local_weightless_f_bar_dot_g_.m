function [tmp_fg,tmp_fg_dvol,tmp_fg_a,tmp_fg_b,tmp_fg_c] = local_weightless_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
[tmp_f_dvol_yk_,tmp_f_a_M_,tmp_f_b_M_,tmp_f_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_f_ykabc_);
[tmp_g_dvol_yk_,tmp_g_a_M_,tmp_g_b_M_,tmp_g_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_g_ykabc_);
tmp_fg_dvol = local_weightless_f_dvol_bar_dot_g_dvol_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,tmp_f_dvol_yk_,tmp_g_dvol_yk_);
tmp_fg_a = sum(conj(tmp_f_a_M_) .* tmp_g_a_M_, 'all');
tmp_fg_b = sum(conj(tmp_f_b_M_) .* tmp_g_b_M_, 'all');
tmp_fg_c = sum(conj(tmp_f_c_M_) .* tmp_g_c_M_, 'all');
tmp_fg = tmp_fg_dvol + tmp_fg_a + tmp_fg_b + tmp_fg_c;
