function [X0] = register_spharm_to_spharm_2(verbose,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_,b_k_Y_);
% registration between molecule_A and molecule_B ;
% ;
% verbose = integer verbosity_level ;
% n_k_p_r = integer maximum k ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk) = k_p_r_value for shell nk ;
% weight_k_p_r_ = real array of length n_k_p_r: weight_k_p_r_(nk) = quadrature weight for shell nk. This should be designed so that 4*pi*sum(weight_k_p_r_) = volume of sphere. ;
% l_max_ = integer array of length n_k_p_r; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_k_Y_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_k_Y_ = complex array of length \sum_{nk} (l_max_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;

n_lm_ = (l_max_+1).^2; n_lm_csum_ = cumsum([0;n_lm_(:)]);
l_max_max = l_max_(end);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);

X0 = 0;
for nk_p_r = 0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
l_max = l_max_(1+nk_p_r); n_lm = n_lm_(1+nk_p_r);
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (1:n_lm);
X0 = X0 + weight_k_p_r_(1+nk_p_r) * sum(conj(a_k_Y_(tmp_ij_)).*(b_k_Y_(tmp_ij_)));
end;%for nk_p_r = 0:n_k_p_r-1;

