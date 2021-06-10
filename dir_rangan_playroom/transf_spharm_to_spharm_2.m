function [b_k_Y_] = transf_spharm_to_spharm_2(verbose,n_k_p_r,k_p_r_,l_max_,a_k_Y_,delta_);
% Applies real-space translation by delta_ ;
% to spherical harmonic expansion a_k_Y_, producing b_k_Y_. ;
% Note that a_ is assumed to represent the molecule in fourier-space. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level ;
% n_k_p_r = number of shells. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_value for shell nk_p_r ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_max_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_max_(nk_p_r)+1)^2 coefficients ;
% a_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% delta_ = real array of displacements; [delta_x, delta_y, delta_z] ;
% ;
% outputs: ;
% ;
% b_k_Y_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_k_Y_ corresponds to translated molecule ;

n_lm_ = (l_max_+1).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);

dWtdkd__ = dwignertdkd__(dWtdkd__l_max_max);
delta_p_r = fnorm(+delta_);
Wt___ = expm_dwignertdkd__(dWtdkd__,n_k_p_r,k_p_r_,l_max_,delta_p_r);
Wt_ = wignert_ODE_0(dWtdkd__,Wt___,n_k_p_r,k_p_r_,l_max_,delta_p_r);
delta_z_c_ = transpose(+delta_);
delta_z_p_r = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2 + delta_z_c_(1+2).^2);
delta_z_p_01 = sqrt(delta_z_c_(1+0).^2 + delta_z_c_(1+1).^2);
delta_z_p_azimu_b = atan2(delta_z_c_(1+1),delta_z_c_(1+0));
delta_z_p_polar_a = atan2(delta_z_p_01,delta_z_c_(1+2));
delta_z_p_euler_pos_ = [0,+delta_z_p_polar_a,+delta_z_p_azimu_b];
delta_z_p_euler_neg_ = [-delta_z_p_azimu_b,-delta_z_p_polar_a,0];
delta_z_c_ = [cos(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);...
              sin(delta_z_p_azimu_b)*sin(delta_z_p_polar_a);...
              cos(delta_z_p_polar_a)]...
             *delta_z_p_r;
W_beta_neg__ = wignerd_b(l_max_max,delta_z_p_euler_neg_(1+1));
W_beta_pos__ = wignerd_b(l_max_max,delta_z_p_euler_pos_(1+1));
tmp_Y_form_ = a_k_Y_;
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_neg__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_neg_);
tmp_Y_form_ = Wt_*tmp_Y_form_; 
tmp_Y_form_ = rotate_spharm_to_spharm_2(0,W_beta_pos__,n_k_p_r,k_p_r_,l_max_,tmp_Y_form_,delta_z_p_euler_pos_);
b_k_Y_ = tmp_Y_form_;



