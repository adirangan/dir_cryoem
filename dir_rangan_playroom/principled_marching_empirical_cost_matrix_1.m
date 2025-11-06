function ...
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_1( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
);
% empirical (2d) cost. ;

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

M_k_q_wkM__ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__);
if numel(unique(n_w_))> 1;
M_k_q_rwM___ = reshape(innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,n_M,M_k_q_wkM__,0),[n_k_p_r,n_w_max,n_M]);
end;%if numel(unique(n_w_))> 1;
if numel(unique(n_w_))==1;
M_k_q_rwM___ = permute(reshape(M_k_q_wkM__,[n_w_max,n_k_p_r,n_M]),1+[1,0,2]); M_k_q_rwM___(:,1+n_w_max/2,:) = 0;
end;%if numel(unique(n_w_))==1;

X_00_kk__ = zeros(n_k_p_r,n_k_p_r);
for nM=0:n_M-1;
M_k_q_rw__ = M_k_q_rwM___(:,:,1+nM);
tmp_X_00_kk__ = conj(M_k_q_rw__)*transpose(M_k_q_rw__);
X_00_kk__ = X_00_kk__ + tmp_X_00_kk__;
end;%for nM=0:n_M-1;
X_00_kk__ = (2*pi)^2 * X_00_kk__ / max(1,n_M) ;

X_01_kk__ = zeros(n_k_p_r,n_k_p_r);
M_k_q_r0M_ = sum(reshape(M_k_q_rwM___(:,1+0,:),[n_k_p_r,n_M]),2);
X_01_kk__ = (2*pi)^2 * conj(M_k_q_r0M_) * transpose(M_k_q_r0M_) / max(1,n_M^2) ;

X_weight_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
X_weight_r_(1+nk_p_r) = sqrt(weight_2d_k_p_r_(1+nk_p_r));
end;%for nk_p_r=0:n_k_p_r-1;

X_kk__ = diag(X_weight_r_) * (2*real(X_00_kk__) - 2*real(X_01_kk__)) * diag(X_weight_r_) ;
