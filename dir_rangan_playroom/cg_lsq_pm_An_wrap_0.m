function ...
An_a_UX_Y_ycn__ ...
= ...
cg_lsq_pm_An_wrap_0( ...
 n_CTF_rank ...
,n_M ...
,VSCTF_wMc__ ...
,pm_n_k_p_r ...
,pm_n_w_ ...
,pm_l_max_ ...
,a_UX_Y_ycn__ ...
,scatter_from_tensor__ ...
,n_polar_a ...
,n_azimu_b ...
,legendre_evaluate_ljm___ ...
);
pm_n_w_max = max(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_sum = sum(pm_n_w_(1:pm_n_k_p_r));
pm_n_w_csum_ = cumsum([0;pm_n_w_(1:pm_n_k_p_r)]);
pm_l_max_max = max(pm_l_max_);
pm_n_lm_ = (1+pm_l_max_(1:pm_n_k_p_r)).^2;
pm_n_lm_max = max(pm_n_lm_(1:pm_n_k_p_r));
pm_n_lm_sum = sum(pm_n_lm_(1:pm_n_k_p_r));
pm_n_lm_csum_ = cumsum([0;pm_n_lm_(1:pm_n_k_p_r)]);
An_a_UX_Y_ycn__ = zeros(pm_n_w_max*n_M,pm_n_k_p_r);
for pm_nk_p_r=0:pm_n_k_p_r-1;
An_a_UX_Y_ycn__(:,1+pm_nk_p_r) = ...
cg_lsq_pm_An_0( ...
 n_CTF_rank ...
,n_M ...
,VSCTF_wMc__ ...
,pm_n_w_max ...
,pm_l_max_max ...
,a_UX_Y_ycn__(:,1+pm_nk_p_r) ...
,scatter_from_tensor__ ...
,n_polar_a ...
,n_azimu_b ...
,legendre_evaluate_ljm___ ...
);
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;

