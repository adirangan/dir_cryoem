function [b_] = transf_spharm_to_spharm_0(verbose,n_k,k_,n_l_,a_,delta_,sample_d);
% Applies real-space translation by delta_ ;
% to spherical harmonic expansion a_, producing b_. ;
% Note that a_ is assumed to represent the molecule in fourier-space. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level ;
% n_k = integer maximum k ;
% k_ = real array of length n_k; k_(nk) = k_value for shell nk ;
% n_l_ = integer array of length n_k; n_l_(nk) = spherical harmonic order on shell nk; n_l_(nk) corresponds to n_lm_(nk) = (n_l_(nk)+1)^2 coefficients ;
% a_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% delta_ = real array of displacements; [delta_x, delta_y, delta_z] ;
% sample_d = real sampling resolution used to generate points on spheres (passed to sample_shell_0) ;
% ;
% outputs: ;
% ;
% b_ = complex array of length \sum_{nk} (n_lm_(nk)+1)^2 ; coefficients are ordered in a row, with m varying quickly and l varying slowly ;
% b_ corresponds to translated molecule ;

if (nargin<7); sample_d = 1.0; end;

n_lm_ = (n_l_+1).^2;
k_max = k_(end);
n_l_max = n_l_(end);
m_max_ = -n_l_max : +n_l_max;
n_m_max = length(m_max_);

[n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,a_all_,kx_all_,ky_all_,kz_all_] = convert_spharm_to_k_p_0(verbose,n_k,k_,n_l_,a_,sample_d); 

a_all_ = transf_k_p_to_k_p_0(verbose,n_all,k_all_,theta_all_,phi_all_,a_all_,delta_);

[b_] = convert_k_p_to_spharm_0(verbose,n_k,k_,n_l_,n_all,n_sub_,k_all_,theta_all_,phi_all_,weight_all_,a_all_);




