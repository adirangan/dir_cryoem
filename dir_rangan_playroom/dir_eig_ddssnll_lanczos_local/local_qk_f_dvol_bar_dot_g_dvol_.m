function [tmp_fg_dvol] = local_qk_f_dvol_bar_dot_g_dvol_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_,tmp_f_dvol_qk_,tmp_g_dvol_qk_);
n_qk = n_q*n_k_p_r;
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
tmp_f_dvol_qk__ = local_qk__from_qk_(n_q,n_k_p_r,tmp_f_dvol_qk_);
tmp_g_dvol_qk__ = local_qk__from_qk_(n_q,n_k_p_r,tmp_g_dvol_qk_);
weight_3d_riesz_k_p_qk__ = local_qk__from_qk_(n_q,n_k_p_r,weight_3d_riesz_k_p_qk_);
tmp_fg_dvol = sum( (conj(tmp_f_dvol_qk__) .* tmp_g_dvol_qk__) .* weight_3d_riesz_k_p_qk__ , 'all' ) / scaling_volumetric ;

