function M_k_q_rwlM____ = innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q__,l_max);
n_w_max = max(n_w_);
M_k_q_rwlM____ = zeros(n_k_p_r,n_w_max,1+2*l_max,n_M);
for nM=0:n_M-1;
M_k_q_rwlM____(:,:,:,1+nM) = innerproduct_q_k_stretch_quad_stack___0(n_k_p_r,n_w_,M_k_q__(:,1+nM),l_max);
end;%for nM=0:n_M-1;
