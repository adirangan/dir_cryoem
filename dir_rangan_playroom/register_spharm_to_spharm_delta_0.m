function [X_] = register_spharm_to_spharm_delta_0(verbose,n_k,k_,n_l_,a_,b_,n_delta,delta__,sample_d);
% tests registration between molecule_A and molecule_B using an array of delta_ (uses nufft3d3);
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% a_ corresponds to molecule_a, b_ to molecule_b ;
% n_delta = integer number of delta_ translation vectors ;
% delta__ = real array of size (n_delta,3), ;
% delta__(ndelta,:) = [delta_x,delta_y,delta_z] for translation number ndelta. ;
% ;
% X_ = complex array of dimension (n_delta_x,n_delta_y,n_delta_z). ;
% X_(ndelta_x,ndelta_y,ndelta_z) corresponds to the innerproduct between molecule_A and molecule_B, ; 
% where the latter has been translated (in real-space) by displacement ;
% delta_ = [delta__(ndelta,:)] = [delta_x,delta_y,delta_z] ;

if (nargin<9); sample_d = 1.0; end;

n_lm_ = (n_l_+1).^2;
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);

[n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k,k_,n_l_,a_,sample_d);
[n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,b_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k,k_,n_l_,b_,sample_d);

k_max = k_(n_k);
kx_all_resampled_ = pi*kx_all_/k_max ;
ky_all_resampled_ = pi*ky_all_/k_max ;
kz_all_resampled_ = pi*kz_all_/k_max ;
m_x_ = - delta__(:,1)*(k_max/2)/(pi/2);
m_y_ = - delta__(:,2)*(k_max/2)/(pi/2);
m_z_ = - delta__(:,3)*(k_max/2)/(pi/2);
[AB_x_c_nu3d3_,ier] = nufft3d3(n_all,kx_all_resampled_,ky_all_resampled_,kz_all_resampled_,conj(a_all_).*b_all_.*weight_all_,-1,1e-12,n_delta,m_x_,m_y_,m_z_);

X_ = AB_x_c_nu3d3_;
