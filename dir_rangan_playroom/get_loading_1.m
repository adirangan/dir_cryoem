function M_loading__ = get_loading_1(n_residual_loading,n_residual_iteration,syn_n_order,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,n_w_,syn_n_M,M_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_,image_delta_x_,image_delta_y_);

n_w_sum = sum(n_w_);
T_k_p__ = M_k_p__;
for nM=0:syn_n_M-1;
T_k_p__(:,1+nM) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p__(:,1+nM),+image_delta_x_(1+nM),+image_delta_y_(1+nM));
end;%for nM=0:syn_n_M-1;

M_loading__ = get_loading_0(n_residual_loading,n_residual_iteration,syn_n_order,n_k_p_r,weight_k_p_r_,l_max_,n_w_,syn_n_M,T_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
