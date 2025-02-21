%%%%%%%%;
% Calculate volumetric terms. ;
%%%%%%%%;
tmp_t = tic();
parameter_KAPPA = struct('type','KAPPA');
parameter_KAPPA.flag_verbose = 0*flag_verbose;
parameter_KAPPA.kernel_basic_qref_k_eq_d_double = kernel_basic_qref_k_eq_d_double;
parameter_KAPPA.kernel_basic_l_max_use = kernel_basic_l_max_use;
parameter_KAPPA.kernel_basic_l_max_ext = kernel_basic_l_max_ext;
parameter_KAPPA.kernel_basic_l_max_band = kernel_basic_l_max_band;
parameter_KAPPA.flag_recalc_qref_from_data = 1;
parameter_KAPPA.flag_recalc_dtau_qref_from_data = 1;
parameter_KAPPA.flag_recalc_dtau_dtau_qref_from_data = 1;
[ ...
 parameter ...
,KAPPA ...
,a_restore_C2M0_k_p_qk__ ...
,a_restore_C1M1_k_p_qk__ ...
,a_restore_C0M2_k_p_qk__ ...
,dtau_a_restore_C2M0_k_p_qk__ ...
,dtau_a_restore_C1M1_k_p_qk__ ...
,dtau_a_restore_C0M2_k_p_qk__ ...
,dtau_dtau_a_restore_C2M0_k_p_qk__ ...
,dtau_dtau_a_restore_C1M1_k_p_qk__ ...
,dtau_dtau_a_restore_C0M2_k_p_qk__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
kappa_basic_apply_4( ...
 parameter_KAPPA ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,+dtau_euler_polar_a_M_ ...
,+dtau_euler_azimu_b_M_ ...
,+dtau_euler_gamma_z_M_ ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_basic_apply_4 (yes derivatives) time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
if ~isempty(a_restore_C2M0_k_p_qk__); a_restore_C2M0_k_p_qk_ = reshape(a_restore_C2M0_k_p_qk__,[n_qk,1]); end;
if ~isempty(a_restore_C1M1_k_p_qk__); a_restore_C1M1_k_p_qk_ = reshape(a_restore_C1M1_k_p_qk__,[n_qk,1]); end;
if ~isempty(a_restore_C0M2_k_p_qk__); a_restore_C0M2_k_p_qk_ = reshape(a_restore_C0M2_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_a_restore_C2M0_k_p_qk__); dtau_a_restore_C2M0_k_p_qk_ = reshape(dtau_a_restore_C2M0_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_a_restore_C1M1_k_p_qk__); dtau_a_restore_C1M1_k_p_qk_ = reshape(dtau_a_restore_C1M1_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_a_restore_C0M2_k_p_qk__); dtau_a_restore_C0M2_k_p_qk_ = reshape(dtau_a_restore_C0M2_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_dtau_a_restore_C2M0_k_p_qk__); dtau_dtau_a_restore_C2M0_k_p_qk_ = reshape(dtau_dtau_a_restore_C2M0_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_dtau_a_restore_C1M1_k_p_qk__); dtau_dtau_a_restore_C1M1_k_p_qk_ = reshape(dtau_dtau_a_restore_C1M1_k_p_qk__,[n_qk,1]); end;
if ~isempty(dtau_dtau_a_restore_C0M2_k_p_qk__); dtau_dtau_a_restore_C0M2_k_p_qk_ = reshape(dtau_dtau_a_restore_C0M2_k_p_qk__,[n_qk,1]); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%;

a_restore_mid_C2M0_k_p_qk_ = a_restore_C2M0_k_p_qk_;
a_restore_mid_C1M1_k_p_qk_ = a_restore_C1M1_k_p_qk_;
a_restore_mid_C0M2_k_p_qk_ = a_restore_C0M2_k_p_qk_;
dtau_a_restore_mid_C2M0_k_p_qk_ = dtau_a_restore_C2M0_k_p_qk_;
dtau_a_restore_mid_C1M1_k_p_qk_ = dtau_a_restore_C1M1_k_p_qk_;
dtau_a_restore_mid_C0M2_k_p_qk_ = dtau_a_restore_C0M2_k_p_qk_;
dtau_dtau_a_restore_mid_C2M0_k_p_qk_ = dtau_dtau_a_restore_C2M0_k_p_qk_;
dtau_dtau_a_restore_mid_C1M1_k_p_qk_ = dtau_dtau_a_restore_C1M1_k_p_qk_;
dtau_dtau_a_restore_mid_C0M2_k_p_qk_ = dtau_dtau_a_restore_C0M2_k_p_qk_;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing first-derivative (directional derivative only)')); end;
%%%%%%%%;
dtau = 1e-4;
tmp_t = tic();
tmp_parameter_KAPPA = parameter_KAPPA;
tmp_parameter_KAPPA.flag_recalc_qref_from_data = 1;
tmp_parameter_KAPPA.flag_recalc_dtau_qref_from_data = 0;
tmp_parameter_KAPPA.flag_recalc_dtau_dtau_qref_from_data = 0;
[ ...
 ~ ...
,~ ...
,a_restore_pos_C2M0_k_p_qk_ ...
,a_restore_pos_C1M1_k_p_qk_ ...
,a_restore_pos_C0M2_k_p_qk_ ...
] = ...
kappa_basic_apply_4( ...
 tmp_parameter_KAPPA ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,+(euler_polar_a_M_ + dtau*dtau_euler_polar_a_M_) ...
,+(euler_azimu_b_M_ + dtau*dtau_euler_azimu_b_M_) ...
,+(euler_gamma_z_M_ + dtau*dtau_euler_gamma_z_M_) ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_basic_apply_4 (not derivatives) time %0.2fs',tmp_t)); end;
a_restore_pos_C2M0_k_p_qk_ = reshape(a_restore_pos_C2M0_k_p_qk_,[n_qk,1]);
a_restore_pos_C1M1_k_p_qk_ = reshape(a_restore_pos_C1M1_k_p_qk_,[n_qk,1]);
a_restore_pos_C0M2_k_p_qk_ = reshape(a_restore_pos_C0M2_k_p_qk_,[n_qk,1]);
%%%%;
tmp_t = tic();
tmp_parameter_KAPPA.flag_recalc_qref_from_data = 1;
tmp_parameter_KAPPA.flag_recalc_dtau_qref_from_data = 0;
tmp_parameter_KAPPA.flag_recalc_dtau_dtau_qref_from_data = 0;
[ ...
 ~ ...
,~ ...
,a_restore_neg_C2M0_k_p_qk_ ...
,a_restore_neg_C1M1_k_p_qk_ ...
,a_restore_neg_C0M2_k_p_qk_ ...
] = ...
kappa_basic_apply_4( ...
 tmp_parameter_KAPPA ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,+(euler_polar_a_M_ - dtau*dtau_euler_polar_a_M_) ...
,+(euler_azimu_b_M_ - dtau*dtau_euler_azimu_b_M_) ...
,+(euler_gamma_z_M_ - dtau*dtau_euler_gamma_z_M_) ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_basic_apply_4 (not derivatives) time %0.2fs',tmp_t)); end;
a_restore_neg_C2M0_k_p_qk_ = reshape(a_restore_neg_C2M0_k_p_qk_,[n_qk,1]);
a_restore_neg_C1M1_k_p_qk_ = reshape(a_restore_neg_C1M1_k_p_qk_,[n_qk,1]);
a_restore_neg_C0M2_k_p_qk_ = reshape(a_restore_neg_C0M2_k_p_qk_,[n_qk,1]);
%%%%%%%%;

%%%%%%%%;
dtau_a_restore_dif_C2M0_k_p_qk_ = (a_restore_pos_C2M0_k_p_qk_ - a_restore_neg_C2M0_k_p_qk_)/max(1e-12,2*dtau);
dtau_a_restore_dif_C1M1_k_p_qk_ = (a_restore_pos_C1M1_k_p_qk_ - a_restore_neg_C1M1_k_p_qk_)/max(1e-12,2*dtau);
dtau_a_restore_dif_C0M2_k_p_qk_ = (a_restore_pos_C0M2_k_p_qk_ - a_restore_neg_C0M2_k_p_qk_)/max(1e-12,2*dtau);
fnorm_disp(flag_verbose,'dtau_a_restore_dif_C2M0_k_p_qk_',dtau_a_restore_dif_C2M0_k_p_qk_,'dtau_a_restore_mid_C2M0_k_p_qk_',dtau_a_restore_mid_C2M0_k_p_qk_);
fnorm_disp(flag_verbose,'dtau_a_restore_dif_C1M1_k_p_qk_',dtau_a_restore_dif_C1M1_k_p_qk_,'dtau_a_restore_mid_C1M1_k_p_qk_',dtau_a_restore_mid_C1M1_k_p_qk_);
fnorm_disp(flag_verbose,'dtau_a_restore_dif_C0M2_k_p_qk_',dtau_a_restore_dif_C0M2_k_p_qk_,'dtau_a_restore_mid_C0M2_k_p_qk_',dtau_a_restore_mid_C0M2_k_p_qk_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row=2; p_col=3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C2M0_k_p_qk_-dtau_a_restore_mid_C2M0_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau C2M0 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C2M0_k_p_qk_-dtau_a_restore_mid_C2M0_k_p_qk_)./abs(dtau_a_restore_dif_C2M0_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative C2M0 error))','Interpreter','none');
title('log10(abs(dtau C2M0 relative error))','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C1M1_k_p_qk_-dtau_a_restore_mid_C1M1_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau C1M1 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C1M1_k_p_qk_-dtau_a_restore_mid_C1M1_k_p_qk_)./abs(dtau_a_restore_dif_C1M1_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative C1M1 error))','Interpreter','none');
title('log10(abs(dtau C1M1 relative error))','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C0M2_k_p_qk_-dtau_a_restore_mid_C0M2_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau C0M2 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_polar_a_qk_,log10(abs(dtau_a_restore_dif_C0M2_k_p_qk_-dtau_a_restore_mid_C0M2_k_p_qk_)./abs(dtau_a_restore_dif_C0M2_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative C0M2 error))','Interpreter','none');
title('log10(abs(dtau C0M2 relative error))','Interpreter','none');
%%%%;
end;%if (flag_disp>1);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_a_restore_dif_C2M0_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_a_restore_mid_C2M0_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_a_restore_dif_C2M0_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_a_restore_mid_C2M0_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_a_restore_dif_C2M0_k_p_qk_),real(dtau_a_restore_mid_C2M0_k_p_qk_),'ro');
plot(imag(dtau_a_restore_dif_C2M0_k_p_qk_),imag(dtau_a_restore_mid_C2M0_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C2M0');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_a_restore_dif_C1M1_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_a_restore_mid_C1M1_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_a_restore_dif_C1M1_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_a_restore_mid_C1M1_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_a_restore_dif_C1M1_k_p_qk_),real(dtau_a_restore_mid_C1M1_k_p_qk_),'ro');
plot(imag(dtau_a_restore_dif_C1M1_k_p_qk_),imag(dtau_a_restore_mid_C1M1_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C1M1');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_a_restore_dif_C0M2_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_a_restore_mid_C0M2_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_a_restore_dif_C0M2_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_a_restore_mid_C0M2_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_a_restore_dif_C0M2_k_p_qk_),real(dtau_a_restore_mid_C0M2_k_p_qk_),'ro');
plot(imag(dtau_a_restore_dif_C0M2_k_p_qk_),imag(dtau_a_restore_mid_C0M2_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C0M2');
%%%%;
sgtitle(sprintf('dtau_scatterplot'),'Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% testing second-derivative (directional derivative only)')); end;
%%%%%%%%;
dtau_dtau_a_restore_dif_C2M0_k_p_qk_ = (a_restore_pos_C2M0_k_p_qk_ - 2*a_restore_mid_C2M0_k_p_qk_ + a_restore_neg_C2M0_k_p_qk_)/max(1e-12,dtau*dtau);
fnorm_disp(flag_verbose,'dtau_dtau_a_restore_dif_C2M0_k_p_qk_',dtau_dtau_a_restore_dif_C2M0_k_p_qk_,'dtau_dtau_a_restore_mid_C2M0_k_p_qk_',dtau_dtau_a_restore_mid_C2M0_k_p_qk_);
dtau_dtau_a_restore_dif_C1M1_k_p_qk_ = (a_restore_pos_C1M1_k_p_qk_ - 2*a_restore_mid_C1M1_k_p_qk_ + a_restore_neg_C1M1_k_p_qk_)/max(1e-12,dtau*dtau);
fnorm_disp(flag_verbose,'dtau_dtau_a_restore_dif_C1M1_k_p_qk_',dtau_dtau_a_restore_dif_C1M1_k_p_qk_,'dtau_dtau_a_restore_mid_C1M1_k_p_qk_',dtau_dtau_a_restore_mid_C1M1_k_p_qk_);
dtau_dtau_a_restore_dif_C0M2_k_p_qk_ = (a_restore_pos_C0M2_k_p_qk_ - 2*a_restore_mid_C0M2_k_p_qk_ + a_restore_neg_C0M2_k_p_qk_)/max(1e-12,dtau*dtau);
fnorm_disp(flag_verbose,'dtau_dtau_a_restore_dif_C0M2_k_p_qk_',dtau_dtau_a_restore_dif_C0M2_k_p_qk_,'dtau_dtau_a_restore_mid_C0M2_k_p_qk_',dtau_dtau_a_restore_mid_C0M2_k_p_qk_);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row=2; p_col=3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C2M0_k_p_qk_-dtau_dtau_a_restore_mid_C2M0_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau_dtau C2M0 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C2M0_k_p_qk_-dtau_dtau_a_restore_mid_C2M0_k_p_qk_)./abs(dtau_dtau_a_restore_dif_C2M0_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative error))','Interpreter','none');
title('log10(abs(dtau_dtau C2M0 relative error))','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C1M1_k_p_qk_-dtau_dtau_a_restore_mid_C1M1_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau_dtau C1M1 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C1M1_k_p_qk_-dtau_dtau_a_restore_mid_C1M1_k_p_qk_)./abs(dtau_dtau_a_restore_dif_C1M1_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative error))','Interpreter','none');
title('log10(abs(dtau_dtau C1M1 relative error))','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C0M2_k_p_qk_-dtau_dtau_a_restore_mid_C0M2_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(error))','Interpreter','none');
title('log10(abs(dtau_dtau C0M2 error))','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
plot(k_p_polar_a_qk_,log10(abs(dtau_dtau_a_restore_dif_C0M2_k_p_qk_-dtau_dtau_a_restore_mid_C0M2_k_p_qk_)./abs(dtau_dtau_a_restore_dif_C0M2_k_p_qk_)),'.');
xlim([0,pi]); xlabel('k_p_polar_a_qk_','Interpreter','none');
ylabel('log10(abs(relative error))','Interpreter','none');
title('log10(abs(dtau_dtau C0M2 relative error))','Interpreter','none');
%%%%;
end;%if (flag_disp>1);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 1; p_col = 3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_dtau_a_restore_dif_C2M0_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_dtau_a_restore_mid_C2M0_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_dtau_a_restore_dif_C2M0_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_dtau_a_restore_mid_C2M0_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_dtau_a_restore_dif_C2M0_k_p_qk_),real(dtau_dtau_a_restore_mid_C2M0_k_p_qk_),'ro');
plot(imag(dtau_dtau_a_restore_dif_C2M0_k_p_qk_),imag(dtau_dtau_a_restore_mid_C2M0_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C2M0');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_dtau_a_restore_dif_C1M1_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_dtau_a_restore_mid_C1M1_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_dtau_a_restore_dif_C1M1_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_dtau_a_restore_mid_C1M1_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_dtau_a_restore_dif_C1M1_k_p_qk_),real(dtau_dtau_a_restore_mid_C1M1_k_p_qk_),'ro');
plot(imag(dtau_dtau_a_restore_dif_C1M1_k_p_qk_),imag(dtau_dtau_a_restore_mid_C1M1_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C1M1');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
hold on;
rlim_dif_ = prctile(real(dtau_dtau_a_restore_dif_C0M2_k_p_qk_),[  0,100]); rlim_dif_ = mean(rlim_dif_) + 0.5*1.25*diff(rlim_dif_)*[-1,+1]; if diff(rlim_dif_)==0; rlim_dif_=[-1,+1]; end;
rlim_mid_ = prctile(real(dtau_dtau_a_restore_mid_C0M2_k_p_qk_),[  0,100]); rlim_mid_ = mean(rlim_mid_) + 0.5*1.25*diff(rlim_mid_)*[-1,+1]; if diff(rlim_mid_)==0; rlim_mid_=[-1,+1]; end;
rlim_ = [min([rlim_dif_,rlim_mid_]),max([rlim_dif_,rlim_mid_])];
ilim_dif_ = prctile(imag(dtau_dtau_a_restore_dif_C0M2_k_p_qk_),[  0,100]); ilim_dif_ = mean(ilim_dif_) + 0.5*1.25*diff(ilim_dif_)*[-1,+1]; if diff(ilim_dif_)==0; ilim_dif_=[-1,+1]; end;
ilim_mid_ = prctile(imag(dtau_dtau_a_restore_mid_C0M2_k_p_qk_),[  0,100]); ilim_mid_ = mean(ilim_mid_) + 0.5*1.25*diff(ilim_mid_)*[-1,+1]; if diff(ilim_mid_)==0; ilim_mid_=[-1,+1]; end;
ilim_ = [min([ilim_dif_,ilim_mid_]),max([ilim_dif_,ilim_mid_])];
dlim_ = [min([rlim_,ilim_]),max([rlim_,ilim_])];
plot(dlim_,dlim_,'-','Color',0.85*[1,1,1]);
plot(real(dtau_dtau_a_restore_dif_C0M2_k_p_qk_),real(dtau_dtau_a_restore_mid_C0M2_k_p_qk_),'ro');
plot(imag(dtau_dtau_a_restore_dif_C0M2_k_p_qk_),imag(dtau_dtau_a_restore_mid_C0M2_k_p_qk_),'bx');
hold off;
xlim(dlim_); ylim(dlim_);
axisnotick;
xlabel('dif'); ylabel('mid');
title('C0M2');
%%%%;
sgtitle(sprintf('dtau_dtau_scatterplot'),'Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%;
