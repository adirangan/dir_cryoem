tmp_t = tic();
tab_a = 0; tab_b = 1; tab_c = 2;
J_ori_gpu_abc012wT____ = gpuArray( zeros(n_3,n_3,n_w_max,n_T) );
pinv_J_ori_gpu_abc012wT____ = gpuArray( zeros(n_3,n_3,n_w_max,n_T) );
J_rot_gpu_abc012wT____ = gpuArray( zeros(n_3,n_3,n_w_max,n_T) );
pinv_J_rot_gpu_abc012wT____ = gpuArray( zeros(n_3,n_3,n_w_max,n_T) );
for nw=0:n_w_max-1;
for nT=0:n_T-1;
J_ori_gpu_abc012__ = [ ...
; dtau_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_a) , dtau_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_a) , dtau_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_a) ...
; dtau_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_b) , dtau_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_b) , dtau_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_b) ...
; dtau_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_c) , dtau_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_c) , dtau_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_c) ...
];
pinv_J_ori_gpu_abc012__ = pinv(J_ori_gpu_abc012__,tolerance_pinv);
J_rot_gpu_abc012__ = [ ...
; dtau_R_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_a) , dtau_R_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_a) , dtau_R_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_a) ...
; dtau_R_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_b) , dtau_R_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_b) , dtau_R_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_b) ...
; dtau_R_k_c_0_gpu_wTabc___(1+nw,1+nT,1+tab_c) , dtau_R_k_c_1_gpu_wTabc___(1+nw,1+nT,1+tab_c) , dtau_R_k_c_2_gpu_wTabc___(1+nw,1+nT,1+tab_c) ...
];
pinv_J_rot_gpu_abc012__ = pinv(J_rot_gpu_abc012__,tolerance_pinv);
J_ori_gpu_abc012wT____(:,:,1+nw,1+nT) = J_ori_gpu_abc012__;
pinv_J_ori_gpu_abc012wT____(:,:,1+nw,1+nT) = pinv_J_ori_gpu_abc012__;
J_rot_gpu_abc012wT____(:,:,1+nw,1+nT) = J_rot_gpu_abc012__;
pinv_J_rot_gpu_abc012wT____(:,:,1+nw,1+nT) = pinv_J_rot_gpu_abc012__;
end;%for nT=0:n_T-1;
end;%for nw=0:n_w_max-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% pinv_J_xxx_gpu_abc012__: time %0.6fs',tmp_t)); end;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Note: compare R_use__*dtemplated012_rot_gpu_ to dtemplated012_ori_gpu_.')); end;
if (flag_verbose>0); disp(sprintf(' %% Copying first-derivatives to dtemplateda_rec_gpu_wkT__.')); end;
tmp_t = tic();
dtemplated0_ori_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated1_ori_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated2_ori_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated0_rot_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated1_rot_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated2_rot_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated0_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated1_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplated2_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplateda_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplatedb_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
dtemplatedc_rec_gpu_wkT__ = gpuArray( zeros(n_w_sum,n_T) );
for nw=0:n_w_max-1;
for nT=0:n_T-1;
J_ori_gpu_abc012__ = J_ori_gpu_abc012wT____(:,:,1+nw,1+nT);
pinv_J_ori_gpu_abc012__ = pinv_J_ori_gpu_abc012wT____(:,:,1+nw,1+nT);
J_rot_gpu_abc012__ = J_rot_gpu_abc012wT____(:,:,1+nw,1+nT);
pinv_J_rot_gpu_abc012__ = pinv_J_rot_gpu_abc012wT____(:,:,1+nw,1+nT);
for nk_p_r=0:n_k_p_r-1;
nwk = nw + nk_p_r*n_w_max;
dtemplatedabc_ori_gpu_ = [ ...
;dtemplateda_ori_gpu_wkT__(1+nwk,1+nT) ...
;dtemplatedb_ori_gpu_wkT__(1+nwk,1+nT) ...
;dtemplatedc_ori_gpu_wkT__(1+nwk,1+nT) ...
];
dtemplated012_ori_gpu_ = pinv_J_ori_gpu_abc012__ * dtemplatedabc_ori_gpu_ ;
dtemplatedabc_rot_gpu_ = [ ...
;dtemplateda_rot_gpu_wkT__(1+nwk,1+nT) ...
;dtemplatedb_rot_gpu_wkT__(1+nwk,1+nT) ...
;dtemplatedc_rot_gpu_wkT__(1+nwk,1+nT) ...
];
dtemplated012_rot_gpu_ = pinv_J_rot_gpu_abc012__ * dtemplatedabc_rot_gpu_ ;
if flag_verbose>2;
fnorm_disp(flag_verbose,'R_use__*dtemplated012_rot_gpu_',R_use__*dtemplated012_rot_gpu_,'dtemplated012_ori_gpu_',dtemplated012_ori_gpu_);
end;%if flag_verbose>2;
dtemplated012_rec_gpu_ = R_use__*dtemplated012_rot_gpu_;
dtemplatedabc_rec_gpu_ = J_ori_gpu_abc012__ * dtemplated012_rec_gpu_;
dtemplated0_ori_gpu_wkT__(1+nwk,1+nT) = dtemplated012_ori_gpu_(1+0);
dtemplated1_ori_gpu_wkT__(1+nwk,1+nT) = dtemplated012_ori_gpu_(1+1);
dtemplated2_ori_gpu_wkT__(1+nwk,1+nT) = dtemplated012_ori_gpu_(1+2);
dtemplated0_rot_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rot_gpu_(1+0);
dtemplated1_rot_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rot_gpu_(1+1);
dtemplated2_rot_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rot_gpu_(1+2);
dtemplated0_rec_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rec_gpu_(1+0);
dtemplated1_rec_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rec_gpu_(1+1);
dtemplated2_rec_gpu_wkT__(1+nwk,1+nT) = dtemplated012_rec_gpu_(1+2);
dtemplateda_rec_gpu_wkT__(1+nwk,1+nT) = dtemplatedabc_rec_gpu_(1+tab_a);
dtemplatedb_rec_gpu_wkT__(1+nwk,1+nT) = dtemplatedabc_rec_gpu_(1+tab_b);
dtemplatedc_rec_gpu_wkT__(1+nwk,1+nT) = dtemplatedabc_rec_gpu_(1+tab_c);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nT=0:n_T-1;
end;%for nw=0:n_w_max-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtemplatedx_rec_gpu_wkT__: time %0.6fs',tmp_t)); end;
