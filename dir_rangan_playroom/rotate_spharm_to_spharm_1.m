function [b_] = rotate_spharm_to_spharm_1(verbose,n_k,k_,l_max_,a_,alpha_);
% rotates spherical harmonic expansion a_ by alpha_, producing b_ ;
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% l_max_ = integer array of length n_k; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% alpha_ equals array of euler angles; [alpha, beta, gamma]. Note that this corresponds to rotation k_ by the inverse [-gamma, -beta, -alpha]. ;
% ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ corresponds to rotated molecule ;

n_lm_ = (l_max_+1).^2;
k_max = k_(end);
l_max = l_max_(end);
m_max_ = -l_max : +l_max;
n_m_max = length(m_max_);

if (verbose>1); disp(sprintf(' %% rotating molecule by: [%0.2f %0.2f %0.2f]',+alpha_)); end;
if (verbose>1); disp(sprintf(' %% rotating coordinate_frame by: [%0.2f %0.2f %0.2f]',-alpha_(3),-alpha_(2),-alpha_(1))); end;
b_ = zeros(size(a_));
for nk=1:n_k;
l_max = l_max_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = zeros(size(a_k_));
W_beta_ = wignerd_b(l_max,+alpha_(2));
for l_val=0:l_max;
W_alpha = diag(exp(+i*[-l_val:+l_val]*-alpha_(3)));
W_gamma = diag(exp(+i*[-l_val:+l_val]*-alpha_(1)));
a_k_tmp = a_k_(1+l_val*(l_val+1) + (-l_val:+l_val));
a_k_tmp = reshape(a_k_tmp,2*l_val+1,1);
b_k_tmp = W_alpha*W_beta_{1+l_val}*W_gamma*a_k_tmp;
b_k_(1+l_val*(l_val+1) + (-l_val:+l_val)) = b_k_tmp;
end;%for l_val=0:l_max;
b_(ix_base + (1:n_lm)) = b_k_;
end;%for nk=1:n_k;

