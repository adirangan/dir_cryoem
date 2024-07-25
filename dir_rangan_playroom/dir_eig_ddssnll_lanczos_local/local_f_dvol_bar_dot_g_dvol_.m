function [tmp_fg_dvol] = local_f_dvol_bar_dot_g_dvol_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,tmp_f_dvol_yk_,tmp_g_dvol_yk_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
tmp_f_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_f_dvol_yk_);
tmp_g_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_g_dvol_yk_);
tmp_fg_dvol = sum( (conj(tmp_f_dvol_yk__) .* tmp_g_dvol_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) , 'all' ) / scaling_volumetric ;

