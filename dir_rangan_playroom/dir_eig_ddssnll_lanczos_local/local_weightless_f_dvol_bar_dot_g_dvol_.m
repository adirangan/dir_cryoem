function [tmp_fg_dvol] = local_weightless_f_dvol_bar_dot_g_dvol_(n_k_p_r,l_max_,tmp_f_dvol_yk_,tmp_g_dvol_yk_);
n_y_ = (l_max_+1).^2; n_y_max = max(n_y_); n_y_sum = sum(n_y_); n_y_csum_ = cumsum([0;n_y_]); l_max_max = max(l_max_);
tmp_f_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_f_dvol_yk_);
tmp_g_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_g_dvol_yk_);
tmp_fg_dvol = sum( conj(tmp_f_dvol_yk__) .* tmp_g_dvol_yk__ , 'all' ) ;

