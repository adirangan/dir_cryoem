%%%%%%%%;
Hvv_q3d_k_p_quad_ = dvol_a_k_p_quad_.*a_restore_C2M0_k_p_quad_;
Hvv_q3d_k_Y_quad_yk_ = dvol_a_times_a_restore_C2M0_k_Y_quad_yk_;
Hvv_q3d_k_Y_quad_yk__ = dvol_a_times_a_restore_C2M0_k_Y_quad_yk__;
Hvt_q3d_k_p_quad_ = a_k_p_quad_.*dtau_a_restore_C2M0_k_p_quad_ - dtau_a_restore_C1M1_k_p_quad_;
Hvt_q3d_k_Y_quad_yk_ = a_times_dtau_a_restore_C2M0_k_Y_quad_yk_ - dtau_a_restore_C1M1_k_Y_quad_yk_;
Hvt_q3d_k_Y_quad_yk__ = a_times_dtau_a_restore_C2M0_k_Y_quad_yk__ - dtau_a_restore_C1M1_k_Y_quad_yk__;
%%%%%%%%;
Hv_q3d_k_p_quad_ = Hvv_q3d_k_p_quad_ + Hvt_q3d_k_p_quad_;
Hv_q3d_k_Y_quad_yk__ = Hvv_q3d_k_Y_quad_yk__ + Hvt_q3d_k_Y_quad_yk__;
%%%%%%%%;
Hv_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,Hv_q3d_k_Y_quad_yk__);
%%%%%%%%;
[ ...
 Hv_q3d_k_p_reco_ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,Hv_q3d_k_Y_quad_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
fnorm_disp(flag_verbose,'Hv_q3d_k_p_quad_',Hv_q3d_k_p_quad_,'Hv_q3d_k_p_reco_',Hv_q3d_k_p_reco_);
%%%%%%%%;
if flag_check;
[ ...
 Hv_q3d_k_Y_reco_yk_ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,Hv_q3d_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
fnorm_disp(flag_verbose,'Hv_q3d_k_Y_quad_yk_',Hv_q3d_k_Y_quad_yk_,'Hv_q3d_k_Y_reco_yk_',Hv_q3d_k_Y_reco_yk_);
end;%if flag_check;
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1; p_col=2*1; np=0;
%%%%;
tmp_k_p_quad_ = Hv_q3d_k_p_quad_;
tmp_k_Y_quad_yk__ = Hv_q3d_k_Y_quad_yk__;
tmp_str = 'Hv_q3d';
tmp_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,tmp_k_Y_quad_yk__);
l2_a0 = sum ((conj(tmp_k_Y_quad_yk__).*tmp_k_Y_quad_yk__)*reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) , 'all' )/scaling_volumetric;
disp(sprintf(' %% %s: tmp_k_Y_quad_yk__: l2_a0: %0.16f',tmp_str,l2_a0));
l2_a0 = sum ((conj(tmp_k_p_quad_).*tmp_k_p_quad_).*weight_3d_riesz_k_all_ , 'all' )/scaling_volumetric;
disp(sprintf(' %% %s: tmp_k_p_quad_    : l2_a0: %0.16f',tmp_str,l2_a0));
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,abs(tmp_k_Y_quad_yk_),'k.');
xlabel('Y_l_val_','Interpreter','none');
title(sprintf('abs(%s)',tmp_str),'Interpreter','none');
%%%%;
tmp_abs_k_p_quad_ = abs(tmp_k_p_quad_).*sqrt(weight_3d_riesz_k_all_);
flag_2d_vs_3d = 0; c_use__ = colormap_81s;
nk_p_r = floor(n_k_p_r/2);
tmp_index_ = n_k_all_csum_(1+nk_p_r+0):n_k_all_csum_(1+nk_p_r+1)-1;
lim_ = prctile(tmp_abs_k_p_quad_(1+tmp_index_),[ 5,95]);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_(1+tmp_index_) ...
,k_p_azimu_b_all_(1+tmp_index_) ...
,tmp_abs_k_p_quad_(1+tmp_index_) ... 
,lim_ ... 
,c_use__ ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axis equal; axis vis3d; axisnotick3d;
title(sprintf('abs(%s)',tmp_str),'Interpreter','none');
end;%if (flag_disp>1);
%%%%%%%%;
