function a_k_Y_ = cg_lsq_1(n_order,n_k_p_r,l_max_,n_w_,n_M,M_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
% simply applies conjugate-gradient to the least-squares problem to solve for a_k_Y_. ;

n_lm_ = (1+l_max_).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

a_k_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_Y_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
n_w = n_w_(1+nk_p_r);
tmp_M_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
[k_p_polar_a__,k_p_azimu_b__] = cg_rhs_1(n_M,n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(l_max,cos(linspace(0,pi,n_polar_a)),n_azimu_b);
tensor_to_scatter__ = cg_interpolate_n_1(n_order,n_polar_a,n_azimu_b,n_w*n_M,k_p_polar_a__(:),k_p_azimu_b__(:));
scatter_to_tensor__ = transpose(tensor_to_scatter__);
An__ = @(a_k_Y_) tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]);
At__ = @(a_k_X_) convert_spharm__to_spharm_0(l_max,cg_evaluate_t_1(n_polar_a,n_azimu_b,reshape(scatter_to_tensor__*a_k_X_,[n_polar_a,n_azimu_b]),l_max,legendre_evaluate_mlj___,expil__,expi__));
AtAn__ = @(a_k_Y_) At__(An__(a_k_Y_));
[a_k_Y_(1+tmp_Y_ij_),~] = pcg(AtAn__,At__(reshape(M_k_p__(1+tmp_M_ij_,:),[n_w*n_M,1])));
end;%for nk_p_r=0:n_k_p_r-1;

