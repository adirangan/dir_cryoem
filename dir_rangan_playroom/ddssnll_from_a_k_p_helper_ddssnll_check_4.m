%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

dvt = dvt_check;
%%%%%%%%;
% Calculate the gradient: ;
%%%%%%%%;
dvol_ssnll_q3d_k_p_qk_ = + 1.0 * ( conj(a_k_p_qk_) .* (a_restore_C2M0_k_p_qk_) - conj(a_restore_C1M1_k_p_qk_) );
dtau_ssnll_q2d_M3__ = dtau_ssnll_q2d_M3__ ;
%%%%%%%%;
G_qkabc_ = cat(1,dvol_ssnll_q3d_k_p_qk_,dtau_ssnll_q2d_M3__(:,1+0),dtau_ssnll_q2d_M3__(:,1+1),dtau_ssnll_q2d_M3__(:,1+2));
G_bar_qkabc_ = conj(G_qkabc_); %<-- note here we take a conjugate to allow for standard dot-product. ;
%%%%%%%%;

%%%%%%%%;
% Now calculate the difference-quotient and compare. ;
%%%%%%%%;
dvt_qkabc_ = local_qkabc_from_qk_a_b_c_(n_q,n_k_p_r,n_M,dvol_a_k_p_qk_,dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_);
dvt_qkabc_l2 = local_qkabc_imagecount_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_qk_,n_M,weight_imagecount_M_,dvt_qkabc_,dvt_qkabc_);
dssnll_mid_q2d = local_qkabc_imagecount_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_qk_,n_M,weight_imagecount_M_,G_bar_qkabc_,dvt_qkabc_);
%%%%;
dvt_ = transpose(-2:+2);%dvt_ = transpose(-8:+8);
n_dvt = numel(dvt_);
ssnll_tmp_q2d_dvt_ = zeros(n_dvt,1);
if flag_pwe2;
ssnll_tmp_pw2_dvt_ = zeros(n_dvt,1);
end;%if flag_pwe2;
%%%%;
for ndvt=0:n_dvt-1;
if (flag_verbose>0); disp(sprintf(' %% ndvt %d/%d',ndvt,n_dvt)); end;
tmp_euler_polar_a_M_ = euler_polar_a_M_ + 1.0*dvt*dvt_(1+ndvt)*dtau_euler_polar_a_M_;
tmp_euler_azimu_b_M_ = euler_azimu_b_M_ + 1.0*dvt*dvt_(1+ndvt)*dtau_euler_azimu_b_M_;
tmp_euler_gamma_z_M_ = euler_gamma_z_M_ + 1.0*dvt*dvt_(1+ndvt)*dtau_euler_gamma_z_M_;
[tmp_euler_polar_a_M_,tmp_euler_azimu_b_M_] = periodize_polar_a_azimu_b_0(tmp_euler_polar_a_M_,tmp_euler_azimu_b_M_);
tmp_euler_gamma_z_M_ = periodize(tmp_euler_gamma_z_M_,0,2*pi);
tmp_a_k_p_qk_ = a_k_p_qk_ + 1.0*dvt*dvt_(1+ndvt)*dvol_a_k_p_qk_;
tmp_a_R_k_p_Rqk__ = cell(n_R,1);
for nR=0:n_R-1; tmp_a_R_k_p_Rqk__{1+nR} = a_R_k_p_Rqk__{1+nR} + 1.0*dvt*dvt_(1+ndvt)*dvol_a_R_k_p_Rqk__{1+nR}; end;%for nR=0:n_R-1;
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('n_R','var'); n_R = []; end;
if ~exist('R_use_R___','var'); R_use_R___ = []; end;
if ~exist('a_R_k_p_Rqk__','var'); a_R_k_p_Rqk__=[]; end;
if ~exist('dvol_a_R_k_p_Rqk__','var'); dvol_a_R_k_p_Rqk__=[]; end;
if ~exist('ba_from_single_shell_Rbaba___','var'); ba_from_single_shell_Rbaba___=[]; end;
if ~exist('wS_from_R_single_shell_Rsba___','var'); wS_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwSda_from_R_single_shell_Rsba___','var'); dwSda_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwSdb_from_R_single_shell_Rsba___','var'); dwSdb_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_Rsba___','var'); ddwSdaa_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdab_from_R_single_shell_Rsba___','var'); ddwSdab_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_Rsba___','var'); ddwSdbb_from_R_single_shell_Rsba___=[]; end;
[ ...
 ~ ...
,~ ...
,ssnll_tmp_q2d_dvt_(1+ndvt) ...
] = ...
ssnll_from_a_k_p_15( ...
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
,tmp_a_k_p_qk_ ...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_viewing_S ...
,[] ...
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
,tmp_euler_polar_a_M_ ...
,tmp_euler_azimu_b_M_ ...
,tmp_euler_gamma_z_M_ ...
,[] ...
,[] ...
,[] ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,n_R ...
,R_use_R___ ...
,tmp_a_R_k_p_Rqk__ ...
,[] ...
,ba_from_single_shell_Rbaba___ ...
,wS_from_R_single_shell_Rsba___ ...
,dwSda_from_R_single_shell_Rsba___ ...
,dwSdb_from_R_single_shell_Rsba___ ...
,ddwSdaa_from_R_single_shell_Rsba___ ...
,ddwSdab_from_R_single_shell_Rsba___ ...
,ddwSdbb_from_R_single_shell_Rsba___ ...
);
if flag_pwe2;
tmp_n_source_a = n_source_a + n_source_dvol_a;
tmp_v_source_a_ = cat(1,v_source_a_,+1*dvt*dvt_(1+ndvt)*v_source_dvol_a_);
tmp_delta_a_c__ = cat(2,delta_a_c__,delta_dvol_a_c__);
[ ...
 ~ ...
,~ ...
,ssnll_tmp_pw2_dvt_(1+ndvt) ...
] = ...
ssnll_from_plane_wave_expansion_2( ...
 [] ...
,tmp_n_source_a ...
,tmp_v_source_a_ ...
,tmp_delta_a_c__ ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_max ...
,weight_2d_wk_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_phi_C_ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,n_M ...
,weight_imagecount_M_ ...
,tmp_euler_polar_a_M_ ...
,tmp_euler_azimu_b_M_ ...
,tmp_euler_gamma_z_M_ ...
,[] ...
,[] ...
,[] ...
,n_source_b ...
,v_source_b_ ...
,delta_b_c__ ...
,fromb_polar_a_M_ ...
,fromb_azimu_b_M_ ...
,fromb_gamma_z_M_ ...
,M_phi_M_ ...
,n_S ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_gamma_z_S_ ...
);
end;%if flag_pwe2;
end;%for ndvt=0:n_dvt-1;
%%%%;
ssnll_mid_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2);
ssnll_pos_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2+1);
ssnll_neg_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2-1);
dssnll_dif_q2d = (ssnll_pos_q2d - ssnll_neg_q2d)/max(1e-12,2*dvt);
fnorm_disp(flag_verbose,'dssnll_dif_q2d',dssnll_dif_q2d,'dssnll_mid_q2d',dssnll_mid_q2d,' %<-- can be large.');
if flag_pwe2;
ssnll_mid_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2);
ssnll_pos_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2+1);
ssnll_neg_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2-1);
dssnll_dif_pw2 = (ssnll_pos_pw2 - ssnll_neg_pw2)/max(1e-12,2*dvt);
fnorm_disp(flag_verbose,'dssnll_dif_pw2',dssnll_dif_pw2,'dssnll_mid_q2d',dssnll_mid_q2d,' %<-- should be small.');
end;%if flag_pwe2;
%%%%%%%%;

%%%%%%%%;
tmp_x_ = dvt*dvt_;
tmp_y_ = real(ssnll_tmp_q2d_dvt_);
tmp_x4 = sum(tmp_x_.^4);
tmp_x3 = sum(tmp_x_.^3);
tmp_x2 = sum(tmp_x_.^2);
tmp_x1 = sum(tmp_x_.^1);
tmp_x0 = sum(tmp_x_.^0);
tmp_yx2 = sum(tmp_y_.*tmp_x_.^2);
tmp_yx1 = sum(tmp_y_.*tmp_x_.^1);
tmp_yx0 = sum(tmp_y_.*tmp_x_.^0);
tmp_lhs_ = [ ...
 tmp_x4 tmp_x3 tmp_x2 ...
;tmp_x3 tmp_x2 tmp_x1 ...
;tmp_x2 tmp_x1 tmp_x0 ...
];
tmp_rhs_ = [ tmp_yx2 ; tmp_yx1 ; tmp_yx0 ];
tmp_q2d_p_ = tmp_lhs_\tmp_rhs_;
tmp_q2d_p_x_ = polyval(tmp_q2d_p_,tmp_x_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 12;
linewidth_use = 2;
hold on;
plot(tmp_x_,tmp_q2d_p_x_,'k-','LineWidth',linewidth_use);
plot(tmp_x_,ssnll_tmp_q2d_dvt_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
xlabel('dvt_','Interpreter','none');
ylabel('ssnll_tmp_q2d_dvt_','Interpreter','none');
end;%if (flag_disp>1);
%%%%;
dssnll_lsq_q2d = 1.0*tmp_q2d_p_(1+1);
fnorm_disp(flag_verbose,'dssnll_lsq_q2d',dssnll_lsq_q2d,'dssnll_mid_q2d',dssnll_mid_q2d,' %<-- can be large.');
%%%%%%%%;

%%%%%%%%;
if flag_pwe2;
%%%%%%%%;
tmp_y_ = real(ssnll_tmp_pw2_dvt_);
tmp_yx2 = sum(tmp_y_.*tmp_x_.^2);
tmp_yx1 = sum(tmp_y_.*tmp_x_.^1);
tmp_yx0 = sum(tmp_y_.*tmp_x_.^0);
tmp_rhs_ = [ tmp_yx2 ; tmp_yx1 ; tmp_yx0 ];
tmp_pw2_p_ = tmp_lhs_\tmp_rhs_;
tmp_pw2_p_x_ = polyval(tmp_pw2_p_,tmp_x_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 12;
linewidth_use = 2;
hold on;
plot(tmp_x_,tmp_pw2_p_x_,'k-','LineWidth',linewidth_use);
plot(tmp_x_,ssnll_tmp_pw2_dvt_,'ko','MarkerSize',markersize_use,'MarkerFaceColor','c');
xlabel('dvt_','Interpreter','none');
ylabel('ssnll_tmp_pw2_dvt_','Interpreter','none');
end;%if (flag_disp>1);
%%%%;
dssnll_lsq_pw2 = 1.0*tmp_pw2_p_(1+1);
fnorm_disp(flag_verbose,'dssnll_lsq_pw2',dssnll_lsq_pw2,'dssnll_mid_q2d',dssnll_mid_q2d,' %<-- should be small.');
%%%%%%%%;
end;%if flag_pwe2;
%%%%%%%%;

%%%%%%%%;
% Estimate the second-derivative. ;
%%%%%%%%;
ddssnll_mid_q2d = local_qkabc_imagecount_f_bar_dot_g_(n_q,n_k_p_r,weight_3d_riesz_qk_,n_M,weight_imagecount_M_,Hvt_qkabc_,dvt_qkabc_);
%%%%;
ssnll_mid_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2);
ssnll_pos_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2+1);
ssnll_neg_q2d = ssnll_tmp_q2d_dvt_((1+n_dvt)/2-1);
ddssnll_dif_q2d = (ssnll_pos_q2d - 2*ssnll_mid_q2d + ssnll_neg_q2d)/max(1e-12,dvt^2);
fnorm_disp(flag_verbose,'ddssnll_dif_q2d',ddssnll_dif_q2d,'ddssnll_mid_q2d',ddssnll_mid_q2d,' %<-- can be large.');
%%%%;
ddssnll_lsq_q2d = 2.0*tmp_q2d_p_(1+0);
fnorm_disp(flag_verbose,'ddssnll_lsq_q2d',ddssnll_lsq_q2d,'ddssnll_mid_q2d',ddssnll_mid_q2d,' %<-- can be large.');
%%%%;
if (flag_verbose>2);
disp(sprintf(' %% dssnll_mid_q2d:  %+0.16f',dssnll_mid_q2d));
disp(sprintf(' %% dssnll_dif_q2d:  %+0.16f',dssnll_dif_q2d));
disp(sprintf(' %% dssnll_lsq_q2d:  %+0.16f',dssnll_lsq_q2d));
disp(sprintf(' %% ddssnll_mid_q2d: %+0.16f',ddssnll_mid_q2d));
disp(sprintf(' %% ddssnll_dif_q2d: %+0.16f',ddssnll_dif_q2d));
disp(sprintf(' %% ddssnll_lsq_q2d: %+0.16f',ddssnll_lsq_q2d));
end;%if (flag_verbose>2);
%%%%%%%%;

%%%%%%%%;
if flag_pwe2;
ssnll_mid_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2);
ssnll_pos_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2+1);
ssnll_neg_pw2 = ssnll_tmp_pw2_dvt_((1+n_dvt)/2-1);
ddssnll_dif_pw2 = (ssnll_pos_pw2 - 2*ssnll_mid_pw2 + ssnll_neg_pw2)/max(1e-12,dvt^2);
fnorm_disp(flag_verbose,'ddssnll_dif_pw2',ddssnll_dif_pw2,'ddssnll_mid_q2d',ddssnll_mid_q2d,' %<-- should be small.');
%%%%;
ddssnll_lsq_pw2 = 2.0*tmp_pw2_p_(1+0);
fnorm_disp(flag_verbose,'ddssnll_lsq_pw2',ddssnll_lsq_pw2,'ddssnll_mid_q2d',ddssnll_mid_q2d,' %<-- should be small.');
%%%%;
if (flag_verbose>2);
disp(sprintf(' %% dssnll_mid_pw2:  %+0.16f',dssnll_mid_pw2));
disp(sprintf(' %% dssnll_dif_pw2:  %+0.16f',dssnll_dif_pw2));
disp(sprintf(' %% dssnll_lsq_pw2:  %+0.16f',dssnll_lsq_pw2));
disp(sprintf(' %% ddssnll_mid_q2d: %+0.16f',ddssnll_mid_q2d));
disp(sprintf(' %% ddssnll_dif_pw2: %+0.16f',ddssnll_dif_pw2));
disp(sprintf(' %% ddssnll_lsq_pw2: %+0.16f',ddssnll_lsq_pw2));
end;%if (flag_verbose>2);
end;%if flag_pwe2;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
