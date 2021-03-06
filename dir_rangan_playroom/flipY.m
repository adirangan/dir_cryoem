function a_k_Y_flip_ = flipY(n_k_p_r,l_max_,a_k_Y_orig_);
% inverts the function in real-space. ;
% This sends +x_c_ (a vector in 3-space) to -x_c_. ;
% Correspondingly, this sends +k_c_ (a vector in 3-space) to -k_c_. ;
% and the euler-angles polar_a and azimu_b are sent, respectively, to: ;
% (pi - polar_a) and (pi + azimu_b). ;
% Consequently, the spherical harmonic: ;
% Y_{l}^{m}(+k_c_) = (1/Z_{l}^{m}) P_{l}^{m}(+cos(polar_a)) * exp(+i*m*azimu_b) ;
% is sent to: ;
% Y_{l}^{m}(-k_c_) = (1/Z_{l}^{m}) P_{l}^{m}(+cos(pi - polar_a)) * exp(+i*m*(pi + azimu_b)) ;
%                  = (1/Z_{l}^{m}) P_{l}^{m}(-cos(polar_a)) * exp(+i*m*pi) * exp(+i*m*azimu_b) ;
%                  = (-1)^{l+m} * (-1)^{m} * Y_{l}^{m}(+k_c_) ;
%                  = (-1)^{l} * Y_{l}^{m}(+k_c_). ;
% Thus, to invert x_c_ in real-space we only need to flip the sign ;
% of the odd-order spherical-harmonic coefficients in k_Y_ space. ;
if (numel(l_max_)==1); l_max_ = l_max_(1)*ones(n_k_p_r,1); end;
l_max_max = max(l_max_);
n_lm_ = (1+l_max_).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
a_k_Y_flip_ = a_k_Y_orig_;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
a_k_Y_flip_(1+na) = a_k_Y_orig_(1+na)*(-1)^l_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
