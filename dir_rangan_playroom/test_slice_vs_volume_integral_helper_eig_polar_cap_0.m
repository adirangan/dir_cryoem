n_synth_M = 1024;
polar_cap_width = pi/8;
flag_calc=1;
if flag_calc;
lanczos_n_iteration_max = 32; tmp_lanczos_n_iteration_max = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now try an eigenvalue estimate with noiseless images and polar-cap viewing-angles. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/eig_from_synth_polar_cap.mat',dir_ssnll);
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat); tmp_lanczos_n_iteration_max = numel(tmp_.alph_i_); clear tmp_;
disp(sprintf(' %% %s found, tmp_lanczos_n_iteration_max %.2d',fname_mat,tmp_lanczos_n_iteration_max));
end;%if ( exist(fname_mat,'file'));
if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
disp(sprintf(' %% %s not complete, calculating',fname_mat));
parameter_eig = struct('type','eig');
parameter_eig.flag_verbose = flag_verbose;
parameter_eig.flag_check = 0;
parameter_eig.flag_disp = 0;
parameter_eig.flag_kernel_qpro_d0 = 1;
parameter_eig.flag_kernel_qpro_d1 = 1;
parameter_eig.kernel_qpro_polar_a_pole_north=4.5*pi/24;
parameter_eig.kernel_qpro_polar_a_pole_south=3.5*pi/24;
parameter_eig.kernel_qpro_l_max_use = l_max;
parameter_eig.lanczos_n_iteration_max = lanczos_n_iteration_max;
U_SmallRotation_Delta_ykabc3__ = [];
S_k_p_q2d_wkS__ = S_k_p_wkS__;
%%%%;
if (~exist(fname_mat,'file'));
tmp_lanczos_n_iteration_max = 0;
rng(0);
euler_polar_a_synth_M_ = transpose(pi*mod(1:n_synth_M,2)) + polar_cap_width*0.5*rand(n_synth_M,1);
euler_azimu_b_synth_M_ = 2*pi*rand(n_synth_M,1);
[euler_polar_a_synth_M_,euler_azimu_b_synth_M_] = periodize_polar_a_azimu_b_0(euler_polar_a_synth_M_,euler_azimu_b_synth_M_);
tmp_index_nS_from_nM_ = ...
knnsearch( ...
 [hist2dab_polar_a_S_,hist2dab_azimu_b_S_] ...
,[euler_polar_a_synth_M_,euler_azimu_b_synth_M_] ...
,'K',1)-1;
M_synth_k_p_wkM__ = hist2dab_S_k_p_wkS__(:,1+tmp_index_nS_from_nM_);
euler_polar_a_synth_M_ = hist2dab_polar_a_S_(1+tmp_index_nS_from_nM_);
euler_azimu_b_synth_M_ = hist2dab_azimu_b_S_(1+tmp_index_nS_from_nM_);
euler_gamma_z_synth_M_ = zeros(n_synth_M,1);
v_ykabci__=[];
w_ykabc_=[];
alph_i_=[];
beta_i_=[];
end;%if (~exist(fname_mat,'file'));
%%%%;
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat);
tmp_lanczos_n_iteration_max = numel(tmp_.alph_i_);
tmp_index_nS_from_nM_ = tmp_.tmp_index_nS_from_nM_;
M_synth_k_p_wkM__ = hist2dab_S_k_p_wkS__(:,1+tmp_index_nS_from_nM_);
euler_polar_a_synth_M_ = tmp_.euler_polar_a_synth_M_;
euler_azimu_b_synth_M_ = tmp_.euler_azimu_b_synth_M_;
euler_gamma_z_synth_M_ = tmp_.euler_gamma_z_synth_M_;
v_ykabci__ = tmp_.v_ykabci__;
w_ykabc_ = tmp_.w_ykabc_;
alph_i_ = tmp_.alph_i_;
beta_i_ = tmp_.beta_i_;
clear tmp_;
end;%if ( exist(fname_mat,'file'));
%%%%;
if ~exist('KAPPA','var'); KAPPA=[]; end;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlma____','var'); d0W_betazeta_mlma____=[]; end;
if ~exist('d1W_betazeta_mlma____','var'); d1W_betazeta_mlma____=[]; end;
if ~exist('d2W_betazeta_mlma____','var'); d2W_betazeta_mlma____=[]; end;
if ~exist('U_SmallRotation_Delta_ykabc3__','var'); U_SmallRotation_Delta_ykabc3__=[]; end;
%%%%;
%if ~exist('v_ykabci__ ','var'); v_ykabci__ =[]; end;
%if ~exist('w_ykabc_ ','var'); w_ykabc_ =[]; end;
%if ~exist('alph_i_','var'); alph_i_=[]; end;
%if ~exist('beta_i_ ','var'); beta_i_ =[]; end;
%%%%;
[ ...
 parameter_eig ...
,v_ykabci__  ...
,w_ykabc_  ...
,alph_i_ ...
,beta_i_ ... 
] = ...
eig_ddssnll_lanczos_0( ...
 parameter_eig ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
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
,n_synth_M ...
,M_synth_k_p_wkM__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_synth_M_ ...
,euler_azimu_b_synth_M_ ...
,euler_gamma_z_synth_M_ ...
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
,U_SmallRotation_Delta_ykabc3__ ...
,v_ykabci__  ...
,w_ykabc_  ...
,alph_i_ ...
,beta_i_ ... 
);
save(fname_mat ...
     ,'parameter' ...
     ,'parameter_eig' ...
     ,'n_k_p_r' ...
     ,'k_p_r_' ...
     ,'k_p_r_max' ...
     ,'l_max_' ...
     ,'n_k_all' ...
     ,'n_k_all_csum_' ...
     ,'k_p_r_all_' ...
     ,'k_p_azimu_b_all_' ...
     ,'k_p_polar_a_all_' ...
     ,'weight_3d_k_all_' ...
     ,'weight_shell_k_' ...
     ,'weight_3d_k_p_r_' ...
     ,'weight_3d_riesz_k_p_r_' ...
     ,'weight_3d_riesz_k_all_' ...
     ,'n_w_' ...
     ,'weight_2d_k_p_r_' ...
     ,'weight_2d_wk_' ...
     ,'n_S' ...
     ,'viewing_polar_a_S_' ...
     ,'viewing_azimu_b_S_' ...
     ,'viewing_weight_S_' ...
     ,'n_viewing_polar_a' ...
     ,'viewing_polar_a_' ...
     ,'n_viewing_azimu_b_' ...
     ,'v_ykabci__'  ...
     ,'w_ykabc_' ...
     ,'alph_i_' ...
     ,'beta_i_' ... 
     ,'n_synth_M' ...
     ,'euler_polar_a_synth_M_' ...
     ,'euler_azimu_b_synth_M_' ...
     ,'euler_gamma_z_synth_M_' ...
     ,'tmp_index_nS_from_nM_' ...
     );
%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_calc;
