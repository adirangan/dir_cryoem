function [a_k_Y_norm_,a_avg,a_std,u_k_Y_] = spharm_normalize_1(n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_)
% Normalizes the spherical-harmonic expansion to have norm 1, but does not center. ;
n_lm_ = (l_max_+1).^2;
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
weight_Y_row_ = zeros(n_lm_sum,1);
weight_Y_col_ = zeros(n_lm_sum,1);
weight_Y_val_ = zeros(n_lm_sum,1);
e_k_Y_ = zeros(n_lm_sum,1);
na=0;
for nk_p_r=1:n_k_p_r;
tmp_ij_ = na + (1:n_lm_(nk_p_r));
weight_Y_row_(tmp_ij_) = tmp_ij_;
weight_Y_col_(tmp_ij_) = tmp_ij_;
weight_Y_val_(tmp_ij_) = weight_k_p_r_(nk_p_r);
e_k_Y_(1+na) = 1;
na=na+n_lm_(nk_p_r);
end;%for nk_p_r=1:n_k_p_r;
weight_Y_ = sparse(weight_Y_row_,weight_Y_col_,weight_Y_val_,n_lm_sum,n_lm_sum);
e_avg = ctranspose(e_k_Y_)*(weight_Y_*e_k_Y_);
u_k_Y_ = e_k_Y_./max(sqrt(e_avg),1e-12);
a_avg = ctranspose(u_k_Y_)*(weight_Y_*a_k_Y_);
a_k_Y_norm_ = (a_k_Y_ - 0.0d0*a_avg*u_k_Y_); %<-- do not center. ;
a_std = sqrt(ctranspose(a_k_Y_norm_)*(weight_Y_*a_k_Y_norm_));
a_k_Y_norm_ = a_k_Y_norm_./max(a_std,1e-12);
