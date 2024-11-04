if  isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk__) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_);
dtau_dtau_a_restore_C0M2_k_Y_quad_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,dtau_dtau_a_restore_C0M2_k_Y_quad_yk_);
end;%if  isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk__) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_);
%%%%;
if  isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk__);
dtau_dtau_a_restore_C0M2_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,dtau_dtau_a_restore_C0M2_k_Y_quad_yk__);
end;%if  isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk__);
%%%%;
if  isempty(dtau_dtau_a_restore_C0M2_k_p_quad_) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_);
tmp_yk_ = dtau_dtau_a_restore_C0M2_k_Y_quad_yk_; tmp_str = 'dtau_dtau_a_restore_C0M2_k_Y_quad_yk_';
%%%%;
tmp_t = tic();
%%%%;
[ ...
 tmp_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
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
,tmp_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_a_restore_C0M2: convert_spharm_to_k_p_4: time %0.2fs',tmp_t)); end;
dtau_dtau_a_restore_C0M2_k_p_quad_ = tmp_quad_;
end;%if  isempty(dtau_dtau_a_restore_C0M2_k_p_quad_) & ~isempty(dtau_dtau_a_restore_C0M2_k_Y_quad_yk_);
%%%%;
