%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% comparing first-derivative (reasonably accurate): ')); end;
%%%%%%%%;
dtau_ssnll_q3d = ...
 + 0.5 * sum( (conj(a_times_a_k_Y_quad_yk__) .* dtau_a_restore_C2M0_k_Y_quad_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* dtau_a_restore_C1M1_k_Y_quad_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( dtau_a_restore_C0M2_k_Y_quad_yk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
dtau_ssnll_q3d = dtau_ssnll_q3d / scaling_volumetric ;
fnorm_disp(flag_verbose,'dtau_ssnll_q2d',dtau_ssnll_q2d,'dtau_ssnll_q3d',dtau_ssnll_q3d);
%%%%%%%%;
dtau_ssnll_q3d = ...
 + 0.5 * sum( abs(a_k_p_quad_).^2 .* dtau_a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * 2*real(sum( conj(a_k_p_quad_).^1 .* dtau_a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ )) ...
 + 0.5 * sum( dtau_a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
  ;
dtau_ssnll_q3d = dtau_ssnll_q3d / scaling_volumetric ;
fnorm_disp(flag_verbose,'dtau_ssnll_q2d',dtau_ssnll_q2d,'dtau_ssnll_q3d',dtau_ssnll_q3d);
%%%%%%%%;
dtau_ssnll_q3d = ...
 + 0.5 * sum( (conj(a_k_Y_quad_yk__).^1 .* a_times_dtau_a_restore_C2M0_k_Y_quad_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* dtau_a_restore_C1M1_k_Y_quad_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( dtau_a_restore_C0M2_k_Y_quad_yk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
dtau_ssnll_q3d = dtau_ssnll_q3d / scaling_volumetric ;
fnorm_disp(flag_verbose,'dtau_ssnll_q2d',dtau_ssnll_q2d,'dtau_ssnll_q3d',dtau_ssnll_q3d);
%%%%%%%%;
