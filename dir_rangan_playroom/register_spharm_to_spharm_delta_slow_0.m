function [X_slow_] = register_spharm_to_spharm_delta_slow_0(verbose,n_k,k_,l_max_,a_,b_,n_delta,delta__,sample_d);
% tests registration between molecule_A and molecule_B using an array of delta_ (slow only -- via brute force);
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% l_max_ = integer array of length n_k; l_max_(nk) = spherical harmonic order on shell nk; l_max_(nk) corresponds to n_lm_(nk) = (l_max_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% n_delta = integer number of delta_ translation vectors ;
% delta__ = real array of size (n_delta,3), ;
% delta__(ndelta,:) = [delta_x,delta_y,delta_z] for translation number ndelta. ;
% ;
% X_slow_ = complex array of size (n_delta,1) ;
% X_slow_(ndelta) corresponds to the innerproduct between molecule_A and molecule_B, where ;
% the latter has been translated (in real-space) by displacement ;
% delta_ = [delta__(ndelta,:)] = [delta_x,delta_y,delta_z] ;

if (nargin<9); sample_d = 1.0; end; 

n_lm_ = (l_max_+1).^2;
k_max = k_(end);
l_max_max = l_max_(end);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);

[n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k,k_,l_max_,a_,sample_d);
[n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,b_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k,k_,l_max_,b_,sample_d);

X_slow_ = zeros(n_delta,1);

for ndelta = 1:n_delta;
delta_ = delta__(ndelta,:);
[c_all_] = transf_k_p_to_k_p_0(verbose,n_all,k_all_,theta_all_,phi_all_,b_all_,delta_);
X_slow_(ndelta) = dot(a_all_,c_all_.*weight_all_);
end;%for ndelta = 1:n_delta;
