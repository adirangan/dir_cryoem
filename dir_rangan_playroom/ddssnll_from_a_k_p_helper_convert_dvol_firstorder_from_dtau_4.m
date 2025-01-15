%%%%%%%%;
% determine expansion for dvol_a_firstorder. ;
% The likelihood L(F,\tau) can be expressed as: ;
% L(F,t) = L0 + dLdF*dF + dLdt*dt + 0.5*dF*ddLdFF*dF + 0.5*dF*ddLdFdt*dt + 0.5*dt*ddLdtdF*dF + 0.5*dt*ddLdtdt*dt ;
% Or: dL = dLdF*dF + dLdt*dt + 0.5*dF*ddLdFF*dF + 0.5*dF*ddLdFdt*dt + 0.5*dt*ddLdtdF*dF + 0.5*dt*ddLdtdt*dt ;
% Now at a critical-point we assume dL=0 (e.g., dF=dt=0). ;
% However, if dt is nonzero, then dF must be chosen to lie at a minimum of dL. ;
% In this case we have: ;
% 0 = \partial_{dF} dL = dLdF + Re(ddLdFdF*dF) + Re(ddLdFdt*dt) ;
% Alternatively, we can note that (via the volumetric framework) F = F_{bp} = C1M1 / C2M0 ;
% Thus, dF should equal \partial_{\tau}F_{bp} = (dtau_C1M1*C2M0 - C1M1*dtau_C2M0) / C2M0^2 ;
%%%%%%%%;
% dvol_a = (dtau_C1M1*C2M0 - C1M1*dtau_C2M0) / (C2M0^2) ;
%%%%%%%%;
dvol_a_firstorder_k_p_qk_ = ....
( ... 
+dtau_a_restore_C1M1_k_p_qk_.*a_restore_C2M0_k_p_qk_ ...
-a_restore_C1M1_k_p_qk_.*dtau_a_restore_C2M0_k_p_qk_ ...
) ...
./ max(1e-12,abs(a_restore_C2M0_k_p_qk_).^2) ...
;
%%%%%%%%;
dvol_a_firstorder_k_p_qk__ = reshape(dvol_a_firstorder_k_p_qk_,[n_q,n_k_p_r]);
