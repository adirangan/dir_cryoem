if flag_check;
tmp_H_q3d_k_p_quad = @(dvol_a_k_p_quad_) ...
  ( ...
    + 0.5 * sum( abs(a_k_p_quad_ + dvol_a_k_p_quad_).^2 .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    - 0.5 * 2*real(sum( conj(a_k_p_quad_ + dvol_a_k_p_quad_).^1 .* a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ )) ...
    + 0.5 * sum( a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    ) ...
  / scaling_volumetric  ...
  ;
%%%%%%%%;
tmp_dvol_a_k_p_quad_ = dvol_a_k_p_quad_;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,1,1)
dvol_randn_k_p_quad_ = randn(n_k_all,1);
dvol_randn_k_p_quad_ = dvol_randn_k_p_quad_ * fnorm(tmp_dvol_a_k_p_quad_) / max(1e-12,fnorm(dvol_randn_k_p_quad_));
hold on;
for ntmp=-20:+20;
plot(ntmp,real(tmp_H_q3d_k_p_quad(tmp_dvol_a_k_p_quad_ + dvol_randn_k_p_quad_*ntmp/20)),'ko');
end;%for ntmp=-20:+20;
xlabel('abitrary step'); ylabel('value of tmp_H_q3d_k_p_quad','Interpreter','none');
title('criticality of tmp_dvol_a_k_p_quad_','Interpreter','none');
drawnow();
end;%if (flag_disp>1);
end;%if flag_check;
%%%%%%%%;
