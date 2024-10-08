[dvol_a_k_Y_quad_yk_,dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,bsxfun(@times,denomator_root_weight_3d_riesz_ykabc_,w_tilde_ykabc_));

tmp_t = tic();
parameter_ddssnll = parameter;
parameter_ddssnll.flag_verbose = 0;
parameter_ddssnll.flag_check = 0;
parameter_ddssnll.flag_disp = 0;
[ ...
 parameter_ddssnll ...
,w_nottilde_ykabc_ ...
] = ...
ddssnll_1( ...
 parameter_ddssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,a_k_Y_quad_yk__ ...
,dvol_a_k_Y_quad_yk_ ...
,[] ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,n_eta ...
,index_neta_from_nM_ ...
,eta_k_p_r_ke__ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,KAPPA ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ddssnll_1 (H): time %0.2fs',tmp_t)); end;

w_tilde_ykabc_ = bsxfun(@times,numerator_root_weight_3d_riesz_ykabc_,w_nottilde_ykabc_);
