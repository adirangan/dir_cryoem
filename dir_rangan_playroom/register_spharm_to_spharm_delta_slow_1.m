function [X_slow_] = register_spharm_to_spharm_delta_slow_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_quad_,b_k_Y_quad_,n_delta,delta__);
% tests registration between molecule_A and molecule_B using an array of delta_ (slow only -- via brute force);
% delta is applied to b. ;

[a_k_p_all_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_quad_);
[b_k_p_all_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,b_k_Y_quad_);

X_slow_ = zeros(n_delta,1);
for ndelta = 1:n_delta;
delta_ = delta__(ndelta,:);
[c_k_p_all_] = transf_k_p_to_k_p_1(verbose,n_k_all,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,b_k_p_all_,delta_);
X_slow_(ndelta) = sum(conj(a_k_p_all_).*c_k_p_all_.*weight_k_all_);
end;%for ndelta = 1:n_delta;
