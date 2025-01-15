%%%%%%%%;
Hvt_q3d = ...
  + 1.0 * sum( conj(dvol_a_k_p_qk_) .* (a_k_p_qk_) .* (dtau_a_restore_C2M0_k_p_qk_) .* weight_3d_riesz_qk_ ) ...
  - 1.0 * sum( conj(dvol_a_k_p_qk_) .* (dtau_a_restore_C1M1_k_p_qk_) .* weight_3d_riesz_qk_ ) ...
  ;
Hvt_q3d = Hvt_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Hvt_q3d: %0.16f',Hvt_q3d)); end;
Hvt_q2d = conj(Htv_q2d) ; %<-- due to symmetry. ;
fnorm_disp(flag_verbose,'Hvt_q2d',Hvt_q2d,'Hvt_q3d',Hvt_q3d);
%%%%%%%%;
