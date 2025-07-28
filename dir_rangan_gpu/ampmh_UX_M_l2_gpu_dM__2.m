function UX_M_l2_gpu_dM__ = ampmh_UX_M_l2_gpu_dM__2(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_gpu_lwnM____);
flag_verbose=1;

f_zero = gpuArray( single(0.0));
if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');
tmp_t = tic();
svd_VUXM_gpu_lwnM____ = gpuArray( (svd_VUXM_gpu_lwnM____));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_VUXM_gpu_lwnM____ %0.6fs',tmp_t)); end;
end;%if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');

tmp_t = tic();
svd_U_d_expiw_gpu_s__ = gpuArray( (FTK.svd_U_d_expiw_s__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_U_d_expiw_gpu_s__ %0.6fs',tmp_t)); end;

tmp_t = tic();
UX_M_l2_gpu_dM__ = zeros(FTK.n_delta_v,n_M,'like',f_zero);
n_w_max = max(n_w_);
UX_M_k_q_gpu_dwnM____ = reshape(svd_U_d_expiw_gpu_s__*reshape(svd_VUXM_gpu_lwnM____,[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M]);
UX_M_l2_gpu_dM__(:,:) = sum(reshape(abs(permute(UX_M_k_q_gpu_dwnM____,[2,3,1,4])).^2,[n_w_max*pm_n_UX_rank,FTK.n_delta_v,n_M]),1) * n_w_max;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %%  UX_M_l2_gpu_dM__ %0.6fs',tmp_t)); end;
