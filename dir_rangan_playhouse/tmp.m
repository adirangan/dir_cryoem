clear;

n_M = 1;
n_k_p_r = 3;
n_w_max = 4;

M_k_q_wkM___ = zeros(n_w_max,n_k_p_r,n_M);
for nM=0:n_M-1;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_max-1;
M_k_q_wkM___(1+nw,1+nk_p_r,1+nM) = (0+nM) + (1+nk_p_r) + (1+nw)*i;
end;%for nw=0:n_w_max-1;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nM=0:n_M-1;
M_k_q_wkM__ = reshape(M_k_q_wkM___,[n_w_max*n_k_p_r,n_M]);

weight_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
weight_k_p_r_(1+nk_p_r) = (1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;

n_w_ = n_w_max*ones(n_k_p_r,1);

M_k_q_rwM___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r, weight_k_p_r_, n_w_, n_M, M_k_q_wkM__);

l_max = 2;
M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack___0(n_k_p_r,n_w_,M_k_q_wkM__(:,1+0),l_max);
