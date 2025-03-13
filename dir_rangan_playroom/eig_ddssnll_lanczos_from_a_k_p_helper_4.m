[ ...
 tmp_dvol_a_k_p_qk_ ...
,tmp_dtau_euler_polar_a_M_use_ ...
,tmp_dtau_euler_azimu_b_M_use_ ...
,tmp_dtau_euler_gamma_z_M_use_ ...
] = ...
local_qk_a_b_c_from_qkabc_( ...
 n_q ...
,n_k_p_r ...
,n_M_use ...
,bsxfun( ...
	 @times ...
	 ,denomator_root_weight_3d_riesz_weight_imagecount_qkabc_ ...
	 ,w_tilde_qkabc_ ...
	 ) ...
);
if flag_implicit_dtau;
tmp_dtau_euler_polar_a_M_use_ = []; %<-- set to empty if flag_implicit_dtau. ;
tmp_dtau_euler_azimu_b_M_use_ = []; %<-- set to empty if flag_implicit_dtau. ;
tmp_dtau_euler_gamma_z_M_use_ = []; %<-- set to empty if flag_implicit_dtau. ;
end;%if flag_implicit_dtau;

if flag_disp>2;
disp('%%%%%%%%;');
disp('INPUT');
disp('%%%%%%%%;');
try; disp(sprintf(' %% fnorm(n_qk)=%0.16f',fnorm(n_qk))); catch me; end;
try; disp(sprintf(' %% fnorm(n_qk_csum_)=%0.16f',fnorm(n_qk_csum_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_p_r_qk_)=%0.16f',fnorm(k_p_r_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_p_azimu_b_qk_)=%0.16f',fnorm(k_p_azimu_b_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_p_polar_a_qk_)=%0.16f',fnorm(k_p_polar_a_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_3d_k_p_qk_)=%0.16f',fnorm(weight_3d_k_p_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_shell_qk_)=%0.16f',fnorm(weight_shell_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_k_p_r)=%0.16f',fnorm(n_k_p_r))); catch me; end;
try; disp(sprintf(' %% fnorm(k_p_r_)=%0.16f',fnorm(k_p_r_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_p_r_max)=%0.16f',fnorm(k_p_r_max))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_3d_k_p_r_)=%0.16f',fnorm(weight_3d_k_p_r_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_c_0_qk_)=%0.16f',fnorm(k_c_0_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_c_1_qk_)=%0.16f',fnorm(k_c_1_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(k_c_2_qk_)=%0.16f',fnorm(k_c_2_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_polar_a_k_)=%0.16f',fnorm(n_polar_a_k_))); catch me; end;
try; disp(sprintf(' %% fnorm(polar_a_ka__)=%0.16f',fnorm(polar_a_ka__))); catch me; end;
try; disp(sprintf(' %% fnorm(n_azimu_b_ka__)=%0.16f',fnorm(n_azimu_b_ka__))); catch me; end;
try; disp(sprintf(' %% fnorm(qref_k_eq_d)=%0.16f',fnorm(qref_k_eq_d))); catch me; end;
try; disp(sprintf(' %% fnorm(a_k_p_qk_)=%0.16f',fnorm(a_k_p_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(tmp_dvol_a_k_p_qk_)=%0.16f',fnorm(tmp_dvol_a_k_p_qk_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_w_)=%0.16f',fnorm(n_w_))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_2d_k_p_r_)=%0.16f',fnorm(weight_2d_k_p_r_))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_2d_wk_)=%0.16f',fnorm(weight_2d_wk_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_viewing_S)=%0.16f',fnorm(n_viewing_S))); catch me; end;
try; disp(sprintf(' %% fnorm(S_k_p_q2d_wkS__)=%0.16f',fnorm(S_k_p_q2d_wkS__))); catch me; end;
try; disp(sprintf(' %% fnorm(viewing_polar_a_S_)=%0.16f',fnorm(viewing_polar_a_S_))); catch me; end;
try; disp(sprintf(' %% fnorm(viewing_azimu_b_S_)=%0.16f',fnorm(viewing_azimu_b_S_))); catch me; end;
try; disp(sprintf(' %% fnorm(viewing_weight_S_)=%0.16f',fnorm(viewing_weight_S_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_viewing_polar_a)=%0.16f',fnorm(n_viewing_polar_a))); catch me; end;
try; disp(sprintf(' %% fnorm(viewing_polar_a_)=%0.16f',fnorm(viewing_polar_a_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_viewing_azimu_b_)=%0.16f',fnorm(n_viewing_azimu_b_))); catch me; end;
try; disp(sprintf(' %% fnorm(viewing_gamma_z_S_)=%0.16f',fnorm(viewing_gamma_z_S_))); catch me; end;
try; disp(sprintf(' %% fnorm(n_M)=%0.16f',fnorm(n_M))); catch me; end;
try; disp(sprintf(' %% fnorm(weight_imagecount_M_)=%0.16f',fnorm(weight_imagecount_M_))); catch me; end;
try; disp(sprintf(' %% fnorm(M_k_p_wkM__)=%0.16f',fnorm(M_k_p_wkM__))); catch me; end;
try; disp(sprintf(' %% fnorm(n_CTF)=%0.16f',fnorm(n_CTF))); catch me; end;
try; disp(sprintf(' %% fnorm(index_nCTF_from_nM_)=%0.16f',fnorm(index_nCTF_from_nM_))); catch me; end;
try; disp(sprintf(' %% fnorm(CTF_k_p_wkC__)=%0.16f',fnorm(CTF_k_p_wkC__))); catch me; end;
try; disp(sprintf(' %% fnorm(n_eta)=%0.16f',fnorm(n_eta))); catch me; end;
try; disp(sprintf(' %% fnorm(index_neta_from_nM_)=%0.16f',fnorm(index_neta_from_nM_))); catch me; end;
try; disp(sprintf(' %% fnorm(eta_k_p_wke__)=%0.16f',fnorm(eta_k_p_wke__))); catch me; end;
try; disp(sprintf(' %% fnorm(euler_polar_a_M_)=%0.16f',fnorm(euler_polar_a_M_))); catch me; end;
try; disp(sprintf(' %% fnorm(euler_azimu_b_M_)=%0.16f',fnorm(euler_azimu_b_M_))); catch me; end;
try; disp(sprintf(' %% fnorm(euler_gamma_z_M_)=%0.16f',fnorm(euler_gamma_z_M_))); catch me; end;
try; disp(sprintf(' %% fnorm(tmp_dtau_euler_polar_a_M_use_)=%0.16f',fnorm(tmp_dtau_euler_polar_a_M_use_))); catch me; end;
try; disp(sprintf(' %% fnorm(tmp_dtau_euler_azimu_b_M_use_)=%0.16f',fnorm(tmp_dtau_euler_azimu_b_M_use_))); catch me; end;
try; disp(sprintf(' %% fnorm(tmp_dtau_euler_gamma_z_M_use_)=%0.16f',fnorm(tmp_dtau_euler_gamma_z_M_use_))); catch me; end;
try; disp(sprintf(' %% fnorm(KAPPA)=%0.16f',fnorm(KAPPA))); catch me; end;
try; disp(sprintf(' %% fnorm(wS_from_single_shell_sba__)=%0.16f',fnorm(wS_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(dwSda_from_single_shell_sba__)=%0.16f',fnorm(dwSda_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(dwSdb_from_single_shell_sba__)=%0.16f',fnorm(dwSdb_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdaa_from_single_shell_sba__)=%0.16f',fnorm(ddwSdaa_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdab_from_single_shell_sba__)=%0.16f',fnorm(ddwSdab_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdbb_from_single_shell_sba__)=%0.16f',fnorm(ddwSdbb_from_single_shell_sba__))); catch me; end;
try; disp(sprintf(' %% fnorm(n_R)=%0.16f',fnorm(n_R))); catch me; end;
try; disp(sprintf(' %% fnorm(R_use_R___)=%0.16f',fnorm(R_use_R___))); catch me; end;
try; disp(sprintf(' %% fnorm(a_R_k_p_Rqk__)=%0.16f',fnorm(a_R_k_p_Rqk__))); catch me; end;
try; disp(sprintf(' %% fnorm([])=%0.16f',fnorm([]))); catch me; end;
try; disp(sprintf(' %% fnorm(ba_from_single_shell_Rbaba___)=%0.16f',fnorm(ba_from_single_shell_Rbaba___))); catch me; end;
try; disp(sprintf(' %% fnorm(wS_from_R_single_shell_Rsba___)=%0.16f',fnorm(wS_from_R_single_shell_Rsba___))); catch me; end;
try; disp(sprintf(' %% fnorm(dwSda_from_R_single_shell_Rsba___)=%0.16f',fnorm(dwSda_from_R_single_shell_Rsba___))); catch me; end;
try; disp(sprintf(' %% fnorm(dwSdb_from_R_single_shell_Rsba___)=%0.16f',fnorm(dwSdb_from_R_single_shell_Rsba___))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdaa_from_R_single_shell_Rsba___)=%0.16f',fnorm(ddwSdaa_from_R_single_shell_Rsba___))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdab_from_R_single_shell_Rsba___)=%0.16f',fnorm(ddwSdab_from_R_single_shell_Rsba___))); catch me; end;
try; disp(sprintf(' %% fnorm(ddwSdbb_from_R_single_shell_Rsba___)=%0.16f',fnorm(ddwSdbb_from_R_single_shell_Rsba___))); catch me; end;
end;%if flag_disp>2;

tmp_t = tic();
parameter_ddssnll = parameter;
parameter_ddssnll.flag_verbose = max(0,flag_verbose-1);
parameter_ddssnll.flag_check = 0;
parameter_ddssnll.flag_disp = 0;
[ ...
 parameter_ddssnll ...
,tmp_Hvt_qkabc_ ...
,tmp_Hv_q3d_k_p_qk_ ...
,tmp_Ht_q2d_M3__ ...
,tmp_a_restore_C2M0_k_p_qk_ ...
,tmp_Hvv_q3d_k_p_qk_ ...
,tmp_Hvt_q3d_k_p_qk_ ...
,tmp_Htv_q2d_M3__ ...
,tmp_Htt_q2d_M3__ ...
,tmp_dvol_a_k_p_qk_ ...
,tmp_dtau_euler_polar_a_M_ ...
,tmp_dtau_euler_azimu_b_M_ ...
,tmp_dtau_euler_gamma_z_M_ ...
] = ...
ddssnll_from_a_k_p_4( ...
 parameter_ddssnll ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,qref_k_eq_d ...
,a_k_p_qk_ ...
,tmp_dvol_a_k_p_qk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_viewing_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,n_M ... %<-- original n_M. ;
,weight_imagecount_M_ ... %<-- original weight_imagecount_M_. ;
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,tmp_dtau_euler_polar_a_M_use_ ... %<-- set to empty if flag_implicit_dtau. ;
,tmp_dtau_euler_azimu_b_M_use_ ... %<-- set to empty if flag_implicit_dtau. ;
,tmp_dtau_euler_gamma_z_M_use_ ... %<-- set to empty if flag_implicit_dtau. ;
,KAPPA ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,n_R ...
,R_use_R___ ...
,a_R_k_p_Rqk__ ...
,[] ...
,ba_from_single_shell_Rbaba___ ...
,wS_from_R_single_shell_Rsba___ ...
,dwSda_from_R_single_shell_Rsba___ ...
,dwSdb_from_R_single_shell_Rsba___ ...
,ddwSdaa_from_R_single_shell_Rsba___ ...
,ddwSdab_from_R_single_shell_Rsba___ ...
,ddwSdbb_from_R_single_shell_Rsba___ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ddssnll_from_a_k_p_4 (H): time %0.2fs',tmp_t)); end;

if flag_disp>2;
disp('%%%%%%%%;');
disp('OUTPUT');
disp('%%%%%%%%;');
disp(sprintf(' %% fnorm(tmp_Hvt_qkabc_)=%0.16f',fnorm(tmp_Hvt_qkabc_)));
disp(sprintf(' %% fnorm(tmp_Hv_q3d_k_p_qk_)=%0.16f',fnorm(tmp_Hv_q3d_k_p_qk_)));
disp(sprintf(' %% fnorm(tmp_Ht_q2d_M3__)=%0.16f',fnorm(tmp_Ht_q2d_M3__)));
disp(sprintf(' %% fnorm(tmp_a_restore_C2M0_k_p_qk_)=%0.16f',fnorm(tmp_a_restore_C2M0_k_p_qk_)));
disp(sprintf(' %% fnorm(tmp_Hvv_q3d_k_p_qk_)=%0.16f',fnorm(tmp_Hvv_q3d_k_p_qk_)));
disp(sprintf(' %% fnorm(tmp_Hvt_q3d_k_p_qk_)=%0.16f',fnorm(tmp_Hvt_q3d_k_p_qk_)));
disp(sprintf(' %% fnorm(tmp_Htv_q2d_M3__)=%0.16f',fnorm(tmp_Htv_q2d_M3__)));
disp(sprintf(' %% fnorm(tmp_Htt_q2d_M3__)=%0.16f',fnorm(tmp_Htt_q2d_M3__)));
disp(sprintf(' %% fnorm(tmp_dvol_a_k_p_qk_)=%0.16f',fnorm(tmp_dvol_a_k_p_qk_)));
disp(sprintf(' %% fnorm(tmp_dtau_euler_polar_a_M_)=%0.16f',fnorm(tmp_dtau_euler_polar_a_M_)));
disp(sprintf(' %% fnorm(tmp_dtau_euler_azimu_b_M_)=%0.16f',fnorm(tmp_dtau_euler_azimu_b_M_)));
disp(sprintf(' %% fnorm(tmp_dtau_euler_gamma_z_M_)=%0.16f',fnorm(tmp_dtau_euler_gamma_z_M_)));
end;%if flag_disp>2;

dvol_a_k_p_qk_ = tmp_dvol_a_k_p_qk_;
dtau_euler_polar_a_M_ = tmp_dtau_euler_polar_a_M_;
dtau_euler_azimu_b_M_ = tmp_dtau_euler_azimu_b_M_;
dtau_euler_gamma_z_M_ = tmp_dtau_euler_gamma_z_M_;

w_nottilde_qkabc_ = tmp_Hvt_qkabc_(1:n_qkabc); %<-- ignore alignment component if flag_implicit_dtau. ;
w_tilde_qkabc_ = ...
  bsxfun( ...
	  @times ...
	  ,numerator_root_weight_3d_riesz_weight_imagecount_qkabc_ ...
	  ,w_nottilde_qkabc_ ...
	  );

tmp_dvol_a_k_p_qk_ = [];
tmp_dtau_euler_polar_a_M_use_ = [];
tmp_dtau_euler_azimu_b_M_use_ = [];
tmp_dtau_euler_gamma_z_M_use_ = [];
tmp_dtau_euler_polar_a_M_ = [];
tmp_dtau_euler_azimu_b_M_ = [];
tmp_dtau_euler_gamma_z_M_ = [];

dvol_a_k_p_qki__(:,1+niteration) = dvol_a_k_p_qk_;
dtau_euler_polar_a_Mi__(:,1+niteration) = dtau_euler_polar_a_M_;
dtau_euler_azimu_b_Mi__(:,1+niteration) = dtau_euler_azimu_b_M_;
dtau_euler_gamma_z_Mi__(:,1+niteration) = dtau_euler_gamma_z_M_;
