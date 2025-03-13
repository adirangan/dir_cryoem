function [tmp_fg,tmp_fg_dvol,tmp_fg_a,tmp_fg_b,tmp_fg_c] = local_qkabc_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,n_M,tmp_f_qkabc_,tmp_g_qkabc_);
n_qk = n_q*n_k_p_r;
[tmp_f_dvol_qk_,tmp_f_a_M_,tmp_f_b_M_,tmp_f_c_M_] = local_qk_a_b_c_from_qkabc_(n_q,n_k_p_r,n_M,tmp_f_qkabc_);
[tmp_g_dvol_qk_,tmp_g_a_M_,tmp_g_b_M_,tmp_g_c_M_] = local_qk_a_b_c_from_qkabc_(n_q,n_k_p_r,n_M,tmp_g_qkabc_);
tmp_fg_dvol = local_qk_f_dvol_bar_dot_g_dvol_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,tmp_f_dvol_qk_,tmp_g_dvol_qk_);
tmp_fg_a = sum(conj(tmp_f_a_M_) .* tmp_g_a_M_, 'all');
tmp_fg_b = sum(conj(tmp_f_b_M_) .* tmp_g_b_M_, 'all');
tmp_fg_c = sum(conj(tmp_f_c_M_) .* tmp_g_c_M_, 'all');
tmp_fg = tmp_fg_dvol + tmp_fg_a + tmp_fg_b + tmp_fg_c;
