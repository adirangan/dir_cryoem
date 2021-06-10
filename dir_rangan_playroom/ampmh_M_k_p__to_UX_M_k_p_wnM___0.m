function ...
[UX_M_k_q_wnM___...
,UX_M_k_p_wnM___...
 ] = ...
ampmh_M_k_p__to_UX_M_k_p_wnM___0(...
 n_k_p_r...
,k_p_r_...
,n_w_...
,n_UX_rank...
,UX__...
,X_weight_r_...
,n_M...
,M_k_p__...
,M_k_q__...
,image_delta_x_...
,image_delta_y_...
);

n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

pm_n_k_p_r = n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);

if nargout>1; UX_M_k_p_wnM__ = zeros(n_w_max,n_UX_rank,n_M); end;
UX_M_k_q_wnM__ = zeros(n_w_max,n_UX_rank,n_M);
for nM=0:n_M-1;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p__(:,1+nM),image_delta_x_(1+nM),image_delta_y_(1+nM));
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_M_k_q_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
end;%for nk_p_r=0:n_k_p_r-1;
UX_M_k_q_wnM___(:,1+nUX_rank,1+nM) = tmp_UX_M_k_q_;
if nargout>1; UX_M_k_p_wnM___(:,1+nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_M_k_q_); end;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nM=0:n_M-1;

