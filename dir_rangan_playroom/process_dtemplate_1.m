function ...
process_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_S_k_p__ ...
,da_S_k_p__ ...
,db_S_k_p__ ...
,dc_S_k_p__ ...
,dt_S_k_p__ ...
);
% Scales db_S_k_p__ by sin(polar_a). ;
% This scaling allows da_S_k_p__ and db_S_k_p__ / sin(pm_viewing_polar_a_) to both ;
% correspond to the same absolute speed on the surface of the sphere. ;
pm_n_w_max = max(pm_n_w_);
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);

grad_abc_sigma___ = zeros(3,pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_l2__ = zeros(pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_grad_t___ = zeros(3,pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_grad_t_l2__ = zeros(pm_n_k_p_r,pm_n_viewing_all);
grad_abc_sigma_csum___ = zeros(3,pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_l2_csum__ = zeros(pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_grad_t_csum___ = zeros(3,pm_n_k_p_r,pm_n_viewing_all);
grad_abc_inv_grad_t_l2_csum__ = zeros(pm_n_k_p_r,pm_n_viewing_all);
for nviewing_all=0:pm_n_viewing_all-1;
for nk_p_r=0:pm_n_k_p_r-1;
tmp_index_ = pm_n_w_csum_(1+nk_p_r) + (0:pm_n_w_(1+nk_p_r)-1);
tmp_grad_abc__ = ...
[ ...
 da_S_k_p__(1+tmp_index_,1+nviewing_all) ...
,db_S_k_p__(1+tmp_index_,1+nviewing_all) / sin(pm_viewing_polar_a_all_(1+nviewing_all)) ...
,dc_S_k_p__(1+tmp_index_,1+nviewing_all) ...
];
tmp_grad_t_ = dt_S_k_p__(1+tmp_index_,1+nviewing_all);
tmp_grad_abc_inv_grad_t_ = mldivide(tmp_grad_abc__,tmp_grad_t_);
grad_abc_inv_grad_t___(:,1+nk_p_r,1+nviewing_all) = tmp_grad_abc_inv_grad_t_;
grad_abc_inv_grad_t_l2__(1+nk_p_r,1+nviewing_all) = fnorm(tmp_grad_abc_inv_grad_t_);
tmp_svd_ = svds(tmp_grad_abc__,3);
grad_abc_sigma___(:,1+nk_p_r,1+nviewing_all) = tmp_svd_;
grad_abc_inv_l2__(1+nk_p_r,1+nviewing_all) = sqrt(sum(tmp_svd_.^(-2)));
tmp_index_csum_ = 0:(pm_n_w_csum_(1+nk_p_r)+pm_n_w_(1+nk_p_r)-1);
tmp_grad_abc__ = ...
[ ...
 da_S_k_p__(1+tmp_index_csum_,1+nviewing_all) ...
,db_S_k_p__(1+tmp_index_csum_,1+nviewing_all) / sin(pm_viewing_polar_a_all_(1+nviewing_all)) ...
,dc_S_k_p__(1+tmp_index_csum_,1+nviewing_all) ...
];
tmp_grad_t_ = dt_S_k_p__(1+tmp_index_csum_,1+nviewing_all);
tmp_grad_abc_inv_grad_t_ = mldivide(tmp_grad_abc__,tmp_grad_t_);
grad_abc_inv_grad_t_csum___(:,1+nk_p_r,1+nviewing_all) = tmp_grad_abc_inv_grad_t_;
grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_all) = fnorm(tmp_grad_abc_inv_grad_t_);
tmp_svd_ = svds(tmp_grad_abc__,3);
grad_abc_sigma_csum___(:,1+nk_p_r,1+nviewing_all) = tmp_svd_;
grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_all) = sqrt(sum(tmp_svd_.^(-2)));
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nviewing_all=0:pm_n_viewing_all-1;

grad_abc_inv_l2_ = zeros(pm_n_k_p_r,1);
grad_abc_inv_l2_ = grad_abc_inv_l2__*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]);
grad_abc_inv_l2_csum_ = zeros(pm_n_k_p_r,1);
grad_abc_inv_l2_csum_ = grad_abc_inv_l2_csum__*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]);
grad_abc_inv_grad_t_l2_ = zeros(pm_n_k_p_r,1);
grad_abc_inv_grad_t_l2_ = grad_abc_inv_grad_t_l2__*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]);
grad_abc_inv_grad_t_l2_csum_ = zeros(pm_n_k_p_r,1);
grad_abc_inv_grad_t_l2_csum_ = grad_abc_inv_grad_t_l2_csum__*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]);

process_dtemplate = struct('type','process_dtemplate');
process_dtemplate.grad_abc_inv_l2_ = grad_abc_inv_l2_;
process_dtemplate.grad_abc_sigma___ = grad_abc_sigma___;
process_dtemplate.grad_abc_inv_l2__ = grad_abc_inv_l2__;
process_dtemplate.grad_abc_inv_grad_t___ = grad_abc_inv_grad_t___;
process_dtemplate.grad_abc_inv_grad_t_l2__ = grad_abc_inv_grad_t_l2__;
process_dtemplate.grad_abc_sigma_csum___ = grad_abc_sigma_csum___;
process_dtemplate.grad_abc_inv_l2_csum__ = grad_abc_inv_l2_csum__;
process_dtemplate.grad_abc_inv_grad_t_csum___ = grad_abc_inv_grad_t_csum___;
process_dtemplate.grad_abc_inv_grad_t_l2_csum__ = grad_abc_inv_grad_t_l2_csum__;
process_dtemplate.grad_abc_inv_l2_ = grad_abc_inv_l2_;
process_dtemplate.grad_abc_inv_l2_csum_ = grad_abc_inv_l2_csum_;
process_dtemplate.grad_abc_inv_grad_t_l2_ = grad_abc_inv_grad_t_l2_;
process_dtemplate.grad_abc_inv_grad_t_l2_csum_ = grad_abc_inv_grad_t_l2_csum_;


