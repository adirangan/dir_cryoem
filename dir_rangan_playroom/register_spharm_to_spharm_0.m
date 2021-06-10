function [X0] = register_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,b_);
% calculates innerproduct between molecule_A and molecule_B;
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% ;
% X0 = innerproduct ;

n_lm_ = (n_l_+1).^2;
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);

X0 = 0;
verbose_tab = 0;
for nk = 0:n_k-1;
k_val = k_(1+nk);
n_l = n_l_(1+nk); n_lm = (1+n_l)*(1+n_l);
ix_base = sum(n_lm_(1:1+nk-1)); a_k_ = a_(ix_base + (1:n_lm)); b_k_ = b_(ix_base + (1:n_lm));
X0 = X0 + k_val^2 * dot(a_k_,b_k_);
end;%for nk = 0:n_k-1;

%{
X_tmp=0;
for nk = 1:n_k;
k_val = nk;
n_l = n_l_(nk); n_lm = n_lm_(nk); ix_base = sum(n_lm_(1:nk-1));
a_k_ = a_(ix_base + (1:n_lm)); b_k_ = b_(ix_base + (1:n_lm));
for nl = 0:n_l;
m_ = [-nl:+nl];
a_tmp = a_k_(1+nl*(nl+1) + (-nl:+nl));
a_tmp = reshape(a_tmp,2*nl+1,1);
b_tmp = b_k_(1+nl*(nl+1) + (-nl:+nl));
b_tmp = reshape(b_tmp,2*nl+1,1);
X_tmp = X_tmp + k_val^2 * (ctranspose(a_tmp)*b_tmp);
end;%for nl = 0:n_l;
end;%for nk = 1:n_k;
X0 = X_tmp;
 %}
