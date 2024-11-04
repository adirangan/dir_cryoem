%%%%%%%%;
dtau_M33___ = bsxfun(@times,reshape(dtau_M3__,[n_M,3,1]),reshape(dtau_M3__,[n_M,1,3]));
Htt_q2d = sum( bsxfun(@times,dtau_dtau_ssnll_q2d_M33___.*dtau_M33___,weight_imagecount_M_) , 'all' ) ;
if (flag_verbose>0); disp(sprintf(' %% Htt_q2d: %0.16f',Htt_q2d)); end;
Htt_q2d = sum( bsxfun(@times,dtau_M3__,weight_imagecount_M_) .* sum(bsxfun(@times,dtau_dtau_ssnll_q2d_M33___,reshape(dtau_M3__,[n_M,1,3])),[3]) , 'all' ) ;
if (flag_verbose>0); disp(sprintf(' %% Htt_q2d: %0.16f',Htt_q2d)); end;
Htt_q3d = ...
  + 1.0 * sum( conj(a_k_p_quad_) .* (a_k_p_quad_ ) .* (0.5*dtau_dtau_a_restore_C2M0_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(a_k_p_quad_) .* (0.5*dtau_dtau_a_restore_C1M1_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(0.5*dtau_dtau_a_restore_C1M1_k_p_quad_) .* (a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  + 1.0 * sum( (0.5*dtau_dtau_a_restore_C0M2_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  ;
Htt_q3d = Htt_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Htt_q3d: %0.16f',Htt_q3d)); end;
fnorm_disp(flag_verbose,'Htt_q2d',Htt_q2d,'Htt_q3d',Htt_q3d);
%%%%%%%%;