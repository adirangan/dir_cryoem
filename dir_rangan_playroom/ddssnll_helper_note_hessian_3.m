%%%%%%%%;
% Now construct hessian. ;
%%%%%%%%;
%{
  % As indicated via the above calculation: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_k_p_quad_) .* (a_k_p_quad_) .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    - 0.5 * sum( conj(a_k_p_quad_) .* a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    - 0.5 * sum( conj(a_restore_C1M1_k_p_quad_) .* (a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
    + 0.5 * sum( a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % rewriting this as: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_) .* (a_) .* C2M0_ .* w_ ) ...
    - 0.5 * sum( conj(a_) .* C1M1_ .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    + 0.5 * sum( C0M2_ .* w_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % we have: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_ + dvol_) .* (a_ + dvol_) .* (C2M0_ + dtau_C2M0_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(a_ + dvol_) .* (C1M1_ + dtau_C1M1_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_ + dtau_C1M1_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_ + dvol_) .* w_ ) ...
    + 0.5 * sum( (C0M2_ + dtau_C0M2_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % where we assume: ;
  % 1. dvol_ is a perturbation of the form dvol_k_p_quad_, ;
  % 2. dtau_ is a perturbation of the form dtau_M3__, not accounting for weight_imagecount_M_, ;
  % 3. dtau_CXMX_ is a jacobian (i.e., of the form dtau_CXMX_M3__), and ;
  % 4. dtau_dtau_CXMX__ is a hessian (i.e., of the form dtau_dtau_CXMX_M33__). ;
  % Using this notation, we can expand and retain terms of order 2 and lower: ;
  ssnll_q3d = ...
    ...
    + 0.5 * sum( conj(a_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    ...
    - 0.5 * sum( conj(a_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (a_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    ...
    + 0.5 * sum( (C0M2_) .* w_ ) ...
    + 0.5 * sum( (dtau_C0M2_*dtau_) .* w_ ) ...
    + 0.5 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % Combining terms of like order: ;
  ssnll_q3d = ...
    ...
    + 0.5 * sum( conj(a_) .* (a_) .* (C2M0_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    + 0.5 * sum( (C0M2_) .* w_ ) ...
    ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (a_) .* w_ ) ...
    + 0.5 * sum( (dtau_C0M2_*dtau_) .* w_ ) ...
    ...
    + 0.5 * sum( conj(dvol_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    + 0.5 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % From this expansion we can extract the gradient as a linear-operator acting on [dvol_;dtau_]: ;
  Gradient of ssnll_q3d = ...
      + 1.0 * real(sum( [ conj(a_) .* (C2M0_) - conj(C1M1) ] .* (dvol_) .* w_ )) ...
      - 1.0 * real(sum( [ conj(a_).*dtau_C1M1_ - 0.5*conj(a_).*a_.*dtau_C2M0_ - 0.5*dtau_C0M2_ ] * (dtau_) .* w_ )) ...
    ;
  Gradient of ssnll_q3d = Gradient of ssnll_q3d / scaling_volumetric ;
  % Note that again we have not accounted for weight_imagecount_M_. ;
  % Similarly, we can extract the hessian as a quadratic-kernel acting on [dvol_;dtau_]: ;
  Hessian of ssnll_q3d = ...
    [ctranspose(dvol_) , ctranspose(dtau_)] * [ H ] * [dvol_;dtau_] ...
    ;
  % where the matrix H is a block-matrix with components: ;
  %      [ Hvv | Hvt ]  ;
  %  H = [-----+-----]  ;
  %      [ Htv | Htt ]  ;
  % such that: ;
  Hvv = quadratic-kernel acting on dvol_ = + 1.0 * sum( conj(dvol_) .* (C2M0_) .* (dvol_) .* w_ ) / scaling_volumetric ;
  Htt = quadratic-kernel acting on dtau_ = ...
    + 1.0 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 1.0 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 1.0 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    + 1.0 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ;
  Htt = Htt / scaling_volumetric ;
  % Note that again we have not accounted for weight_imagecount_M_. ;
  % Note that, via ssnll_from_a_k_Y_12, Htt can be rewritten in the simpler form: ;
  Htt = sum( dtau_M3__(1+nM,1+ntau0) * dtau_dtau_ssnll_q2d_M33___(1+nM,1+ntau0,1+ntau1) * dtau_M3__(1+nM,1+ntau1) ). ;
  % The operator Hvt has a domain compatible with dtau_ and a range compatible with dvol_. ;
  Hvt = ...
    + 1.0 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    - 1.0 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    ;
  Hvt = Hvt / scaling_volumetric ;
  % Note that, given the output of kappa_qpro_apply_2, we can rewrite Hvt as:
  Hvt = ...
      + 1.0 * sum( conj(dvol_) .* (a_) .* (|dtau| * dtau_a_restore_C2M0_k_p_quad_) .* w_ ) ...
      - 1.0 * sum( conj(dvol_) .* (|dtau| * dtau_a_restore_C1M1_k_p_quad_) .* w_ ) ...
      ;
  Hvt = Hvt / scaling_volumetric ;
  % The operator Htv has a domain compatible with dvol_ and a range compatible with dtau_. ;
  Htv = ...
    + 1.0 * sum( conj(dtau_C2M0_*dtau_) .* conj(a_) .* (dvol_) .* w_ ) ...
    - 1.0 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    ;
  Htv = Htv / scaling_volumetric ;
  % Note that, via ssnll_from_a_k_Y_12, Htv can be rewritten in the simpler form: ;
  Htv = sum( dtau_M3__(1+nM,1+ntau1) * dtau_dvol_ssnll_q2d_M3__(1+nM,1+ntau1) ). ;
 %}
%%%%%%%%;
