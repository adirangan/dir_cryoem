function UX_T_M_l2_gpu_dM__ = tfpmh_UX_T_M_l2_gpu_dM__1(n_delta_v,n_svd_l,svd_U_d_expiw_s_gpu__,n_w_max,n_M,pm_n_UX_rank,svd_VUXM_gpu_lwnM____);
flag_verbose=1;

f_zero = gpuArray( single(0.0));
if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');
tmp_t = tic();
svd_VUXM_gpu_lwnM____ = gpuArray( (svd_VUXM_gpu_lwnM____));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_VUXM_gpu_lwnM____ %0.6fs',tmp_t)); end;
end;%if ~strcmp(class(svd_VUXM_gpu_lwnM____),'gpuArray');
if ~strcmp(class(svd_U_d_expiw_s_gpu__),'gpuArray');
tmp_t = tic();
svd_U_d_expiw_s_gpu__ = gpuArray( (svd_U_d_expiw_s_gpu__));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% svd_U_d_expiw_s_gpu__ %0.6fs',tmp_t)); end;
end;%if ~strcmp(class(svd_U_d_expiw_s_gpu__),'gpuArray');

UX_T_M_k_q_gpu_dwnM____ = reshape(pagemtimes(svd_U_d_expiw_s_gpu__,svd_VUXM_gpu_lwnM____),[n_delta_v,n_w_max,pm_n_UX_rank,n_M]);
UX_T_M_l2_gpu_dM__ = zeros(n_delta_v,n_M);
UX_T_M_l2_gpu_dM__ = reshape(sum(reshape(permute(abs(UX_T_M_k_q_gpu_dwnM____).^2,1+[1,2,0,3]),[n_w_max*pm_n_UX_rank,n_delta_v,n_M]),1+0) * n_w_max,[n_delta_v,n_M]);
