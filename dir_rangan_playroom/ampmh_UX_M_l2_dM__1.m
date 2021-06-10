function UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);
n_w_max = max(n_w_);
UX_M_k_q_dwnM____ = reshape(FTK.svd_U_d_expiw_s__(:,:)*reshape(svd_VUXM_lwnM____(:,:,:,:),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M]);
UX_M_l2_dM__(:,:) = sum(reshape(abs(permute(UX_M_k_q_dwnM____,[2,3,1,4])).^2,[n_w_max*pm_n_UX_rank,FTK.n_delta_v,n_M]),1) * n_w_max;
