if flag_check;
%%%%%%%%;
combine_a_restore_C2M0_k_p_qk_ = ...
 a_restore_C2M0_k_p_qk_ ...
+dtau_a_restore_C2M0_k_p_qk_ ...
;
combine_a_restore_C1M1_k_p_qk_ = ...
 a_restore_C1M1_k_p_qk_ ...
+dtau_a_restore_C1M1_k_p_qk_ ...
;
combine_a_restore_C0M2_k_p_qk_ = ...
 a_restore_C0M2_k_p_qk_ ...
+dtau_a_restore_C0M2_k_p_qk_ ...
;
dvol_a_firstorder_times_combine_a_restore_C2M0_k_p_qk_ = ...
 combine_a_restore_C1M1_k_p_qk_ ...
-combine_a_restore_C2M0_k_p_qk_ .* a_k_p_qk_ ...
;
%%%%%%%%;
tmp_H_q3d_k_p_qk = @(dvol_a_k_p_qk_) ...
  ( ...
    + 0.5 * sum( abs(a_k_p_qk_ + dvol_a_k_p_qk_).^2 .* combine_a_restore_C2M0_k_p_qk_ .* weight_3d_riesz_qk_ ) ...
    - 0.5 * 2*real(sum( conj(a_k_p_qk_ + dvol_a_k_p_qk_).^1 .* combine_a_restore_C1M1_k_p_qk_ .* weight_3d_riesz_qk_ )) ...
    + 0.5 * sum( combine_a_restore_C0M2_k_p_qk_ .* weight_3d_riesz_qk_ ) ...
    ) ...
  / scaling_volumetric  ...
  ;
%%%%%%%%;
tmp_dvol_a_k_p_qk_ = dvol_a_k_p_qk_;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,1,1)
dvol_randn_k_p_qk_ = randn(n_qk,1);
dvol_randn_k_p_qk_ = dvol_randn_k_p_qk_ * fnorm(tmp_dvol_a_k_p_qk_) / max(1e-12,fnorm(dvol_randn_k_p_qk_));
hold on;
for ntmp=-20:+20;
plot(ntmp,real(tmp_H_q3d_k_p_qk(tmp_dvol_a_k_p_qk_ + dvol_randn_k_p_qk_*ntmp/20)),'ko');
end;%for ntmp=-20:+20;
xlabel('abitrary step'); ylabel('value of tmp_H_q3d_k_p_qk','Interpreter','none');
title('criticality of tmp_dvol_a_k_p_qk_','Interpreter','none');
sgtitle(sprintf('flag_dtau %d flag_dvol %d (critical only if ~flag_dvol)',flag_dtau,flag_dvol),'Interpreter','none');
drawnow();
end;%if (flag_disp>1);
end;%if flag_check;
%%%%%%%%;
