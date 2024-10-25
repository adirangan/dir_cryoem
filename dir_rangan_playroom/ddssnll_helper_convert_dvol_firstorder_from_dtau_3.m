%%%%%%%%;
% determine expansion for dvol_a_firstorder. ;
%%%%%%%%;
tmp_k_p_quad_ = (+dtau_a_restore_C1M1_k_p_quad_.*a_restore_C2M0_k_p_quad_ - a_restore_C1M1_k_p_quad_.*dtau_a_restore_C2M0_k_p_quad_)/max(1e-12,abs(a_restore_C2M0_k_p_quad_).^2);
tmp_str = 'dvol_a_firstorder_k_Y_yk__';
tmp_t = tic();
[ ...
 dvol_a_firstorder_k_Y_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %s: convert_k_p_to_spharm_4: time %0.2fs',tmp_str,tmp_t)); end;
%%%%;
dvol_a_firstorder_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,dvol_a_firstorder_k_Y_yk_);
%%%%%%%%;

