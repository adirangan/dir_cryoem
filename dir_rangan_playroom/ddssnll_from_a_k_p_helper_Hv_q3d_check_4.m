%%%%%%%%;
Hvv_q3d_k_p_qk_ = dvol_a_k_p_qk_.*a_restore_C2M0_k_p_qk_;
Hvt_q3d_k_p_qk_ = a_k_p_qk_.*dtau_a_restore_C2M0_k_p_qk_ - dtau_a_restore_C1M1_k_p_qk_;
%%%%%%%%;
Hv_q3d_k_p_qk_ = Hvv_q3d_k_p_qk_ + Hvt_q3d_k_p_qk_;
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1; p_col=2*1; np=0;
%%%%;
tmp_k_p_qk_ = Hv_q3d_k_p_qk_;
tmp_str = 'Hv_q3d';
l2_a0 = sum ((conj(tmp_k_p_qk_).*tmp_k_p_qk_).*weight_3d_riesz_qk_ , 'all' )/scaling_volumetric;
disp(sprintf(' %% %s: tmp_k_p_qk_    : l2_a0: %0.16f',tmp_str,l2_a0));
%%%%;
tmp_abs_k_p_qk_ = abs(tmp_k_p_qk_).*sqrt(weight_3d_riesz_qk_);
flag_2d_vs_3d = 0; c_use__ = colormap_81s;
nk_p_r = floor(n_k_p_r/2);
tmp_index_ = n_qk_csum_(1+nk_p_r+0):n_qk_csum_(1+nk_p_r+1)-1;
lim_ = prctile(tmp_abs_k_p_qk_(1+tmp_index_),[ 5,95]);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_qk_(1+tmp_index_) ...
,k_p_azimu_b_qk_(1+tmp_index_) ...
,tmp_abs_k_p_qk_(1+tmp_index_) ... 
,lim_ ... 
,c_use__ ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axis equal; axis vis3d; axisnotick3d;
title(sprintf('abs(%s)',tmp_str),'Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;
