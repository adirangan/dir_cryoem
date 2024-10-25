%%%%%%%%;
Hvv_q3d = ...
 + 1.0 * sum( (conj(dvol_a_times_dvol_a_k_Y_quad_yk__) .* a_restore_C2M0_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
;
Hvv_q3d = Hvv_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Hvv_q3d: %0.16f',Hvv_q3d)); end;
Hvv_q3d = ...
 + 1.0 * sum( (conj(dvol_a_k_Y_quad_yk__) .* dvol_a_times_a_restore_C2M0_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
;
Hvv_q3d = Hvv_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Hvv_q3d: %0.16f',Hvv_q3d)); end;
Hvv_q3d = ...
 + 1.0 * sum( conj(dvol_a_k_Y_quad_yk__) .* bsxfun(@times,dvol_a_times_a_restore_C2M0_k_Y_yk__,reshape(weight_3d_riesz_k_p_r_,[1,n_k_p_r])) , 'all' ) ...
;
Hvv_q3d = Hvv_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Hvv_q3d: %0.16f',Hvv_q3d)); end;
Hvv_q3d = ...
 + 1.0 * sum( conj(dvol_a_k_p_quad_) .* dvol_a_k_p_quad_ .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
;
Hvv_q3d = Hvv_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Hvv_q3d: %0.16f',Hvv_q3d)); end;
%%%%%%%%;
