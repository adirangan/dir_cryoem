tab_a = 0; tab_b = 1; tab_c = 2;
tab_0 = 0; tab_1 = 1; tab_2 = 2;
%%%%%%%%;
tmp_t = tic();
J_ori_abc012wT____ = ...
  permute( ...
	   cat( ...
		4 ...
		,dtau_k_c_0_wTabc___ ...
		,dtau_k_c_1_wTabc___ ...
		,dtau_k_c_2_wTabc___ ...
		) ...
	   ,[3,4,1,2] ...
	   );
J_rot_abc012wT____ = ...
  permute( ...
	   cat( ...
		4 ...
		,dtau_R_k_c_0_wTabc___ ...
		,dtau_R_k_c_1_wTabc___ ...
		,dtau_R_k_c_2_wTabc___ ...
		) ...
	   ,[3,4,1,2] ...
	   );
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% J_xxx_abc012wT____: time %0.6fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[~,pinv_J_ori_abc012wT____] = pagepinv_0(struct('tolerance_pinv',tolerance_pinv),J_ori_abc012wT____);
[~,pinv_J_rot_abc012wT____] = pagepinv_0(struct('tolerance_pinv',tolerance_pinv),J_rot_abc012wT____);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% pinv_J_xxx_abc012wT____: time %0.6fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
dtemplatedabc_ori_abc1wkT_____ = reshape(permute(cat(4,dtemplateda_ori_wkT__,dtemplatedb_ori_wkT__,dtemplatedc_ori_wkT__),[3,4,1,2]),[n_3,n_1,n_w_max,n_k_p_r,n_T]);
dtemplated012_ori_0121wkT_____ = pagemtimes(reshape(pinv_J_ori_abc012wT____,[n_3,n_3,n_w_max,1,n_T]),dtemplatedabc_ori_abc1wkT_____);
dtemplatedabc_rot_abc1wkT_____ = reshape(permute(cat(3,dtemplateda_rot_wkT__,dtemplatedb_rot_wkT__,dtemplatedc_rot_wkT__),[3,1,2]),[n_3,n_1,n_w_max,n_k_p_r,n_T]);
dtemplated012_rot_0121wkT_____ = pagemtimes(reshape(pinv_J_rot_abc012wT____,[n_3,n_3,n_w_max,1,n_T]),dtemplatedabc_rot_abc1wkT_____);
dtemplated012_rec_0121wkT_____ = pagemtimes(R_use__,dtemplated012_rot_0121wkT_____);
dtemplatedabc_rec_abc1wkT_____ = pagemtimes(reshape(J_ori_abc012wT____,[n_3,n_3,n_w_max,1,n_T]),dtemplated012_rec_0121wkT_____);
if (flag_verbose>0); disp(sprintf(' %% Note: compare R_use__*dtemplated012_rot_ to dtemplated012_ori_.')); end;
if flag_verbose>2;
fnorm_disp(flag_verbose,'dtemplated012_rec_0121wkT_____',dtemplated012_rec_0121wkT_____,'dtemplated012_ori_0121wkT_____',dtemplated012_ori_0121wkT_____);
end;%if flag_verbose>2;
%%%%;
dtemplated012_ori_wkT012____ = reshape(permute(dtemplated012_ori_0121wkT_____,[3,4,5,1,2]),[n_w_sum,n_T,n_3,n_1]);
dtemplated0_ori_wkT__ = dtemplated012_ori_wkT012____(:,:,1+tab_0);
dtemplated1_ori_wkT__ = dtemplated012_ori_wkT012____(:,:,1+tab_1);
dtemplated2_ori_wkT__ = dtemplated012_ori_wkT012____(:,:,1+tab_2);
%%%%;
dtemplated012_rot_wkT012____ = reshape(permute(dtemplated012_rot_0121wkT_____,[3,4,5,1,2]),[n_w_sum,n_T,n_3,n_1]);
dtemplated0_rot_wkT__ = dtemplated012_rot_wkT012____(:,:,1+tab_0);
dtemplated1_rot_wkT__ = dtemplated012_rot_wkT012____(:,:,1+tab_1);
dtemplated2_rot_wkT__ = dtemplated012_rot_wkT012____(:,:,1+tab_2);
%%%%;
dtemplated012_rec_wkT012____ = reshape(permute(dtemplated012_rec_0121wkT_____,[3,4,5,1,2]),[n_w_sum,n_T,n_3,n_1]);
dtemplated0_rec_wkT__ = dtemplated012_rec_wkT012____(:,:,1+tab_0);
dtemplated1_rec_wkT__ = dtemplated012_rec_wkT012____(:,:,1+tab_1);
dtemplated2_rec_wkT__ = dtemplated012_rec_wkT012____(:,:,1+tab_2);
%%%%;
dtemplatedabc_rec_wkTabc____ = reshape(permute(dtemplatedabc_rec_abc1wkT_____,[3,4,5,1,2]),[n_w_sum,n_T,n_3,n_1]);
dtemplateda_rec_wkT__ = dtemplatedabc_rec_wkTabc____(:,:,1+tab_a);
dtemplatedb_rec_wkT__ = dtemplatedabc_rec_wkTabc____(:,:,1+tab_b);
dtemplatedc_rec_wkT__ = dtemplatedabc_rec_wkTabc____(:,:,1+tab_c);
%%%%;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtemplatedx_rec_wkT__: time %0.6fs',tmp_t)); end;
