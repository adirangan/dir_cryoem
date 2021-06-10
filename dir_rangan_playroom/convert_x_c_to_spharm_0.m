function [a_] = convert_x_c_to_spharm_0(verbose,n_k_p_r,k_p_r_,weight_k_p_r_,l_val_,X_0_,X_1_,X_2_,a_x_c_,sample_d);
% evaluates spherical-harmonic-expansion a_ (in fourier-space) on a collection of points in real-space. ;
% ;
% inputs: ;
% ;
% verbose = integer verbosity_level. ;
% n_k_p_r = integer maximum k_p_r. ;
% k_p_r_ = real array of length n_k_p_r; k_p_r_(nk_p_r) = k_value for shell nk_p_r. ;
% weight_k_p_r_ = real array of radial quadrature weights. ;
% l_val_ = integer array of length n_k_p_r; l_val_(nk_p_r) = spherical harmonic order on shell nk_p_r; l_val_(nk_p_r) corresponds to n_lm_(nk_p_r) = (l_val_(nk_p_r)+1)^2 coefficients. ;
% X_0_ = real array of length n_x_c ; x_c_0 value to evaluate. ;
% X_1_ = real array of length n_x_c ; x_c_1 value to evaluate. ;
% X_2_ = real array of length n_x_c ; x_c_2 value to evaluate. ;
% a_x_c_ = complex array of length n_x_c ; evaluation of function in real-space. ;
% sample_d = real sampling resolution used to generate points on spheres (passed to sample_shell_0) ;
% ;
% outputs: ;
% ;
% a_ = complex array of length \sum_{nk_p_r} (n_lm_(nk_p_r)+1)^2 ; coefficients are ordered linearly, with m varying quickly and l varying slowly and k varying most slowly. ;

n_lm_ = (l_val_+1).^2;
l_val_max = max(l_val_);
m_max_ = -l_val_max : +l_val_max;
n_m_max = length(m_max_);
n_x_c = numel(X_0_);
assert(numel(X_0_)==numel(X_1_));
assert(numel(X_1_)==numel(X_2_));

n_all_ = zeros(n_k_p_r,1);
for nk_p_r=1:n_k_p_r;
k_p_r = k_p_r_(nk_p_r); 
n_all_(nk_p_r) =  sample_shell_0(k_p_r,sample_d,'L');
end%for nk_p_r=1:n_k_p_r;
n_all = sum(n_all_);
k_p_r_all_ = zeros(n_all,1); azimu_b_all_ = zeros(n_all,1); polar_a_all_ = zeros(n_all,1); weight_all_ = zeros(n_all,1);

n_sub_ = zeros(n_k_p_r,1);
n_sub = 0;
for nk_p_r=1:n_k_p_r;
k_p_r = k_p_r_(nk_p_r); 
n_sub_(nk_p_r) = n_sub;
[length_sub,azimu_b_sub_,polar_a_sub_,weight_sub_] = sample_shell_0(k_p_r,sample_d,'L');
k_p_r_all_(1 + n_sub + (0:length_sub-1)) = k_p_r*ones(length_sub,1);
azimu_b_all_(1 + n_sub + (0:length_sub-1)) = azimu_b_sub_;
polar_a_all_(1 + n_sub + (0:length_sub-1)) = polar_a_sub_;
weight_all_(1 + n_sub + (0:length_sub-1)) = weight_sub_;
n_sub = n_sub + length_sub;
end;%for nk_p_r=1:n_k_p_r;
k_c_0_all_ = k_p_r_all_.*cos(azimu_b_all_).*sin(polar_a_all_);
k_c_1_all_ = k_p_r_all_.*sin(azimu_b_all_).*sin(polar_a_all_);
k_c_2_all_ = k_p_r_all_.*cos(polar_a_all_);

a_all_ = nufft3d3(n_x_c,pi*X_0_(:),pi*X_1_(:),pi*X_2_(:),a_x_c_,-1,1e-12,n_all,2*pi*k_c_0_all_/pi,2*pi*k_c_1_all_/pi,2*pi*k_c_2_all_/pi) / n_x_c / (2*pi);

a_ = convert_k_p_to_spharm_0(verbose,n_k_p_r,k_p_r_,l_val_,n_all,n_sub_,k_p_r_all_,azimu_b_all_,polar_a_all_,weight_all_,a_all_);
