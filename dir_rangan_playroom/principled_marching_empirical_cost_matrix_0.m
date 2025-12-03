function ...
[ ...
 X__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p__ ...
);
% empirical (2d) cost. ;

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;

M_k_q_rwM___ = reshape(innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q__,0),[n_k_p_r,n_w_max,n_M]);

X_00__ = zeros(n_k_p_r,n_k_p_r);
for nM=0:n_M-1;
M_k_q_rw__ = M_k_q_rwM___(:,:,1+nM);
tmp_X_00__ = conj(M_k_q_rw__)*transpose(M_k_q_rw__);
X_00__ = X_00__ + tmp_X_00__;
end;%for nM=0:n_M-1;
X_00__ = (2*pi)^2 * X_00__ / n_M ;

X_01__ = zeros(n_k_p_r,n_k_p_r);
M_k_q_r0M_ = sum(reshape(M_k_q_rwM___(:,1+0,:),[n_k_p_r,n_M]),2);
X_01__ = (2*pi)^2 * conj(M_k_q_r0M_) * transpose(M_k_q_r0M_) / n_M^2 ;

X_weight_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
X_weight_r_(1+nk_p_r) = sqrt(weight_2d_k_p_r_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;

X__ = diag(X_weight_r_) * (2*real(X_00__) - 2*real(X_01__)) * diag(X_weight_r_) ;

