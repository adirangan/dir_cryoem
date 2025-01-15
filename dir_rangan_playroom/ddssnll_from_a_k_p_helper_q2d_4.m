%%%%%%%%;
% Calculate derivative using ssnll_from_a_k_p_14. ;
%%%%%%%%;
dvol_S_k_p_q2d_wkS__=[];
dtau_S_k_p_q2d_wkS3___=[];
dtau_dvol_S_k_p_q2d_wkS3___=[];
dtau_dtau_S_k_p_q2d_wkS33____=[];
%%%%;
tmp_t = tic();
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('R_use__','var'); R_use__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('dvol_a_R_k_p_qk_','var'); dvol_a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
if ~exist('wS_from_R_single_shell_sba__','var'); wS_from_R_single_shell_sba__=[]; end;
if ~exist('dwSda_from_R_single_shell_sba__','var'); dwSda_from_R_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_R_single_shell_sba__','var'); dwSdb_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_sba__','var'); ddwSdaa_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_R_single_shell_sba__','var'); ddwSdab_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_sba__','var'); ddwSdbb_from_R_single_shell_sba__=[]; end;
parameter_ssnll = struct('type','parameter');
parameter_ssnll.flag_verbose = 0*flag_verbose;
parameter_ssnll.tolerance_master = tolerance_master;
parameter_ssnll.tolerance_pinv = tolerance_pinv;
parameter_ssnll.viewing_k_eq_d = viewing_k_eq_d;
parameter_ssnll.template_k_eq_d = template_k_eq_d;
parameter_ssnll.n_order = n_order;
[ ...
 parameter_ssnll ...
,ssnll_q2d_M_ ...
,ssnll_q2d ...
,S_k_p_q2d_wkS__ ...
,dvol_ssnll_q2d_M_ ...
,dvol_ssnll_q2d ...
,dvol_S_k_p_q2d_wkS__ ...
,dvol_dvol_ssnll_q2d ...
,dtau_ssnll_q2d_M3__ ...
,dtau_ssnll_q2d ...
,dtau_S_k_p_q2d_wkS3___ ...
,dtau_dvol_ssnll_q2d_M3__ ...
,dtau_dvol_ssnll_q2d ...
,dtau_dvol_S_k_p_q2d_wkS3___ ...
,dtau_dtau_ssnll_q2d_M33___ ...
,dtau_dtau_ssnll_q2d ...
,dtau_dtau_S_k_p_q2d_wkS33____ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
ssnll_from_a_k_p_14( ...
 parameter_ssnll ...
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
,a_k_p_qk_ ...
,dvol_a_k_p_qk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,dvol_S_k_p_q2d_wkS__ ...
,dtau_S_k_p_q2d_wkS3___ ...
,dtau_dvol_S_k_p_q2d_wkS3___ ...
,dtau_dtau_S_k_p_q2d_wkS33____ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,n_M ...
,weight_imagecount_M_ ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll_from_a_k_p_14: time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%;
dtau = 1e-4;

ssnll_q2d_mid_M_ = ssnll_q2d_M_;
ssnll_q2d_mid = ssnll_q2d;
S_k_p_q2d_mid_wkS__ = S_k_p_q2d_wkS__;
dvol_ssnll_q2d_mid_M_ = dvol_ssnll_q2d_M_;
dvol_ssnll_q2d_mid = dvol_ssnll_q2d;
dvol_S_k_p_q2d_mid_wkS__ = dvol_S_k_p_q2d_wkS__;
dvol_dvol_ssnll_q2d_mid = dvol_dvol_ssnll_q2d;
dtau_ssnll_q2d_mid_M3__ = dtau_ssnll_q2d_M3__;
dtau_ssnll_q2d_mid = dtau_ssnll_q2d;
dtau_S_k_p_q2d_mid_wkS3___ = dtau_S_k_p_q2d_wkS3___;
dtau_dvol_ssnll_q2d_mid_M3__ = dtau_dvol_ssnll_q2d_M3__;
dtau_dvol_ssnll_q2d_mid = dtau_dvol_ssnll_q2d;
dtau_dvol_S_k_p_q2d_mid_wkS3___ = dtau_dvol_S_k_p_q2d_wkS3___;
dtau_dtau_ssnll_q2d_mid_M33___ = dtau_dtau_ssnll_q2d_M33___;
dtau_dtau_ssnll_q2d_mid = dtau_dtau_ssnll_q2d;
dtau_dtau_S_k_p_q2d_mid_wkS33____ = dtau_dtau_S_k_p_q2d_wkS33____;

%%%%%%%%;
tmp_t = tic();
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('R_use__','var'); R_use__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('dvol_a_R_k_p_qk_','var'); dvol_a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
if ~exist('wS_from_R_single_shell_sba__','var'); wS_from_R_single_shell_sba__=[]; end;
if ~exist('dwSda_from_R_single_shell_sba__','var'); dwSda_from_R_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_R_single_shell_sba__','var'); dwSdb_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_sba__','var'); ddwSdaa_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_R_single_shell_sba__','var'); ddwSdab_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_sba__','var'); ddwSdbb_from_R_single_shell_sba__=[]; end;
tmp_parameter_ssnll = parameter_ssnll;
tmp_parameter_ssnll.flag_verbose = 0*flag_verbose;
tmp_parameter_ssnll.tolerance_master = tolerance_master;
tmp_parameter_ssnll.tolerance_pinv = tolerance_pinv;
tmp_parameter_ssnll.viewing_k_eq_d = viewing_k_eq_d;
tmp_parameter_ssnll.template_k_eq_d = template_k_eq_d;
tmp_parameter_ssnll.n_order = n_order;
[ ...
 tmp_parameter_ssnll ...
,ssnll_q2d_pos_M_ ...
,ssnll_q2d_pos ...
,S_k_p_q2d_pos_wkS__ ...
,dvol_ssnll_q2d_pos_M_ ...
,dvol_ssnll_q2d_pos ...
,dvol_S_k_p_q2d_pos_wkS__ ...
,dvol_dvol_ssnll_q2d_pos ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
ssnll_from_a_k_p_14( ...
 tmp_parameter_ssnll ...
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
,a_k_p_qk_ ...
,dvol_a_k_p_qk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,n_M ...
,weight_imagecount_M_ ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,+(euler_polar_a_M_ + 1.0*dtau*dtau_euler_polar_a_M_) ...
,+(euler_azimu_b_M_ + 1.0*dtau*dtau_euler_azimu_b_M_) ...
,+(euler_gamma_z_M_ + 1.0*dtau*dtau_euler_gamma_z_M_) ...
,[] ...
,[] ...
,[] ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll_from_a_k_p_14: time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('R_use__','var'); R_use__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('dvol_a_R_k_p_qk_','var'); dvol_a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
if ~exist('wS_from_R_single_shell_sba__','var'); wS_from_R_single_shell_sba__=[]; end;
if ~exist('dwSda_from_R_single_shell_sba__','var'); dwSda_from_R_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_R_single_shell_sba__','var'); dwSdb_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_sba__','var'); ddwSdaa_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_R_single_shell_sba__','var'); ddwSdab_from_R_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_sba__','var'); ddwSdbb_from_R_single_shell_sba__=[]; end;
tmp_parameter_ssnll = parameter_ssnll;
tmp_parameter_ssnll.flag_verbose = 0*flag_verbose;
tmp_parameter_ssnll.tolerance_master = tolerance_master;
tmp_parameter_ssnll.tolerance_pinv = tolerance_pinv;
tmp_parameter_ssnll.viewing_k_eq_d = viewing_k_eq_d;
tmp_parameter_ssnll.template_k_eq_d = template_k_eq_d;
tmp_parameter_ssnll.n_order = n_order;
[ ...
 tmp_parameter_ssnll ...
,ssnll_q2d_neg_M_ ...
,ssnll_q2d_neg ...
,S_k_p_q2d_neg_wkS__ ...
,dvol_ssnll_q2d_neg_M_ ...
,dvol_ssnll_q2d_neg ...
,dvol_S_k_p_q2d_neg_wkS__ ...
,dvol_dvol_ssnll_q2d_neg ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
ssnll_from_a_k_p_14( ...
 tmp_parameter_ssnll ...
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
,a_k_p_qk_ ...
,dvol_a_k_p_qk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,n_M ...
,weight_imagecount_M_ ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,+(euler_polar_a_M_ - 1.0*dtau*dtau_euler_polar_a_M_) ...
,+(euler_azimu_b_M_ - 1.0*dtau*dtau_euler_azimu_b_M_) ...
,+(euler_gamma_z_M_ - 1.0*dtau*dtau_euler_gamma_z_M_) ...
,[] ...
,[] ...
,[] ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,dvol_a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll_from_a_k_p_14: time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
dtau_M3__ = [ dtau_euler_polar_a_M_ , dtau_euler_azimu_b_M_ , dtau_euler_gamma_z_M_ ] ;
dtau_dtau_M33___ = bsxfun(@times,reshape(dtau_M3__,[n_M,n_3,n_1]),reshape(dtau_M3__,[n_M,n_1,n_3]));
%%%%%%%%;
dtau_ssnll_q2d_dif_M_ = (ssnll_q2d_pos_M_ - ssnll_q2d_neg_M_)/max(1e-12,2*dtau);
dtau_ssnll_q2d_mid_M_ = sum(dtau_ssnll_q2d_mid_M3__.*dtau_M3__,2);
fnorm_disp(flag_verbose,'dtau_ssnll_q2d_dif_M_',dtau_ssnll_q2d_dif_M_,'dtau_ssnll_q2d_mid_M_',dtau_ssnll_q2d_mid_M_);
dtau_ssnll_q2d_dif = (ssnll_q2d_pos - ssnll_q2d_neg)/max(1e-12,2*dtau);
dtau_ssnll_q2d_mid = dtau_ssnll_q2d_mid;
fnorm_disp(flag_verbose,'dtau_ssnll_q2d_dif',dtau_ssnll_q2d_dif,'dtau_ssnll_q2d_mid',dtau_ssnll_q2d_mid);
dtau_dvol_ssnll_q2d_dif_M_ = (dvol_ssnll_q2d_pos_M_ - dvol_ssnll_q2d_neg_M_)/max(1e-12,2*dtau);
dtau_dvol_ssnll_q2d_mid_M_ = sum(dtau_dvol_ssnll_q2d_mid_M3__.*dtau_M3__,2);
fnorm_disp(flag_verbose,'dtau_dvol_ssnll_q2d_dif_M_',dtau_dvol_ssnll_q2d_dif_M_,'dtau_dvol_ssnll_q2d_mid_M_',dtau_dvol_ssnll_q2d_mid_M_);
dtau_dvol_ssnll_q2d_dif = (dvol_ssnll_q2d_pos - dvol_ssnll_q2d_neg)/max(1e-12,2*dtau);
dtau_dvol_ssnll_q2d_mid = dtau_dvol_ssnll_q2d_mid;
fnorm_disp(flag_verbose,'dtau_dvol_ssnll_q2d_dif',dtau_dvol_ssnll_q2d_dif,'dtau_dvol_ssnll_q2d_mid',dtau_dvol_ssnll_q2d_mid);
dtau_dtau_ssnll_q2d_dif_M_ = ( ssnll_q2d_pos_M_ - 2*ssnll_q2d_mid_M_ + ssnll_q2d_neg_M_ ) / max(1e-12,dtau^2) ;
dtau_dtau_ssnll_q2d_mid_M_ = sum(dtau_dtau_ssnll_q2d_mid_M33___.*dtau_dtau_ssnll_q2d_mid_M33___,[2,3]);
fnorm_disp(flag_verbose,'dtau_dtau_ssnll_q2d_dif_M_',dtau_dtau_ssnll_q2d_dif_M_,'dtau_dtau_ssnll_q2d_mid_M_',dtau_dtau_ssnll_q2d_mid_M_);
dtau_dtau_ssnll_q2d_dif = ( ssnll_q2d_pos - 2*ssnll_q2d_mid + ssnll_q2d_neg ) / max(1e-12,dtau^2) ;
dtau_dtau_ssnll_q2d_mid = dtau_dtau_ssnll_q2d_mid;
fnorm_disp(flag_verbose,'dtau_dtau_ssnll_q2d_dif',dtau_dtau_ssnll_q2d_dif,'dtau_dtau_ssnll_q2d_mid',dtau_dtau_ssnll_q2d_mid);
%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%;

