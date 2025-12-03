function UX_T_M_l2_dM__ = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
% assumes n_w_ = n_w_max*ones(n_k_p_r,1);
str_thisfunction = 'tfpmh_UX_T_M_l2_dM__1';
UX_T_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);
n_w_max = max(n_w_); n_k_p_r = numel(n_w_);
n_w_sum = sum(n_w_); if (n_w_sum~=n_w_max*n_k_p_r); disp(sprintf(' %% Warning, n_w_sum %d ~= n_w_max*n_k_p_r %d*%d in %s',n_w_sum,n_w_max,n_k_p_r,str_thisfunction)); end;
UX_T_M_k_q_dwnM____ = reshape(pagemtimes(FTK.svd_U_d_expiw_s__,svd_VUXM_lwnM____),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M]);
UX_T_M_l2_dM__ = reshape(sum(reshape(permute(abs(UX_T_M_k_q_dwnM____).^2,1+[1,2,0,3]),[n_w_max*pm_n_UX_rank,FTK.n_delta_v,n_M]),1+0) * n_w_max,[FTK.n_delta_v,n_M]);
