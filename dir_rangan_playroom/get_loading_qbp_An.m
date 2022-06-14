function ...
P_k_p_wnM__ ...
= ...
get_loading_qbp_An( ...
 n_k_p_r ...
,l_max_ ...
,n_w_ ...
,a_k_Y_yn_ ...
,legendre_evaluate_ljm___ ...
,tensor_to_scatter__ ...
,n_M ...
,CTF_k_p_rwM__ ...
);

l_max_ = l_max_(1:n_k_p_r);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); 
n_w_max = max(n_w_(1:n_k_p_r));
n_w_sum = sum(n_w_(1:n_k_p_r));
n_w_csum_ = cumsum([0;n_w_(1:n_k_p_r)]);

l_max = max(l_max_); assert(min(l_max_)==l_max);
n_w = max(n_w_); assert(min(n_w_)==n_w);
n_polar_a = ceil(n_w/2);
n_azimu_b = max(1+2*l_max,2*n_polar_a);

for nk_p_r=0:n_k_p_r-1;
An___{1+nk_p_r} = @(a_k_Y_) CTF_k_p_rwM__{1+nk_p_r}.*(tensor_to_scatter__*reshape(cg_evaluate_n_1(l_max,convert_spharm_to_spharm__0(l_max,a_k_Y_),n_polar_a,n_azimu_b,legendre_evaluate_ljm___),[n_polar_a*n_azimu_b,1]));
end;%for nk_p_r=0:n_k_p_r-1;

P_k_p_wnM__ = zeros(n_w_sum,n_M);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = n_lm_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
index_Y_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
index_nw_ = n_w_csum_(1+nk_p_r) + (0:n_w-1);
An__ = An___{1+nk_p_r};
P_k_p_wnM__(1+index_nw_,:) = reshape(An__(a_k_Y_yn_(1+index_Y_)),[n_w,n_M]);
end;%for nk_p_r=0:n_k_p_r-1;

