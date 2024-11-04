%%%%%%%%;
a_restore_C2M0_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,a_restore_C2M0_k_Y_quad_yk__);
%%%%%%%%;
dtau_a_restore_C2M0_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,dtau_a_restore_C2M0_k_Y_quad_yk__);
%%%%%%%%;
dtau_a_restore_C1M1_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,dtau_a_restore_C1M1_k_Y_quad_yk__);
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row=4; p_col=2*2; np=0;
%%%%;
for ntype=0:8-1;
if ntype==0; tmp_k_p_quad_ = a_k_p_quad_; tmp_k_Y_quad_yk__ = a_k_Y_quad_yk__; tmp_str = 'a_k_Y'; end;
if ntype==1; tmp_k_p_quad_ = dvol_a_k_p_quad_; tmp_k_Y_quad_yk__ = dvol_a_k_Y_quad_yk__; tmp_str = 'dvol_a_k_Y'; end;
if ntype==2; tmp_k_p_quad_ = a_restore_C2M0_k_p_quad_; tmp_k_Y_quad_yk__ = a_restore_C2M0_k_Y_quad_yk__; tmp_str = 'a_restore_C2M0'; end;
if ntype==3; tmp_k_p_quad_ = a_restore_C1M1_k_p_quad_; tmp_k_Y_quad_yk__ = a_restore_C1M1_k_Y_quad_yk__; tmp_str = 'a_restore_C1M1'; end;
if ntype==4; tmp_k_p_quad_ = dtau_a_restore_C2M0_k_p_quad_; tmp_k_Y_quad_yk__ = dtau_a_restore_C2M0_k_Y_quad_yk__; tmp_str = 'dtau_a_restore_C2M0'; end;
if ntype==5; tmp_k_p_quad_ = dtau_a_restore_C1M1_k_p_quad_; tmp_k_Y_quad_yk__ = dtau_a_restore_C1M1_k_Y_quad_yk__; tmp_str = 'dtau_a_restore_C1M1'; end;
if ntype==6; tmp_k_p_quad_ = dvol_a_times_a_restore_C2M0_k_p_quad_; tmp_k_Y_quad_yk__ = dvol_a_times_a_restore_C2M0_k_Y_quad_yk__; tmp_str = 'dvol_a_times_a_restore_C2M0'; end;
if ntype==7; tmp_k_p_quad_ = a_times_dtau_a_restore_C2M0_k_p_quad_; tmp_k_Y_quad_yk__ = a_times_dtau_a_restore_C2M0_k_Y_quad_yk__; tmp_str = 'a_times_dtau_a_restore_C2M0'; end;
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
mlim_ = prctile(real(tmp_k_p_quad_),[  0,100]);
tmp_abs_k_p_quad_ = abs(tmp_k_p_quad_).*sqrt(weight_3d_riesz_k_all_);
flag_2d_vs_3d = 0; c_use__ = colormap_81s;
nk_p_r = floor(n_k_p_r/2);
tmp_index_ = n_k_all_csum_(1+nk_p_r+0):n_k_all_csum_(1+nk_p_r+1)-1;
plim_ = prctile(tmp_abs_k_p_quad_(1+tmp_index_),[ 5,95]);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_(1+tmp_index_) ...
,k_p_azimu_b_all_(1+tmp_index_) ...
,tmp_abs_k_p_quad_(1+tmp_index_) ... 
,plim_ ... 
,c_use__ ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axis equal; axis vis3d; axisnotick3d;
title(sprintf('abs(%s) mlim_[%+0.6f,%+0.6f]',tmp_str,mlim_),'Interpreter','none');
end;%for ntype=0:8-1;
end;%if (flag_disp>1);
%%%%%%%%;
