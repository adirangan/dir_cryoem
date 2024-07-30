flag_calc=1;
if flag_calc;
lanczos_n_iteration_max = 32; tmp_lanczos_n_iteration_max = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_M_use = 1024;
index_nCTF_from_nM_use_ = index_nCTF_from_nM_(1:n_M_use);
str_M = sprintf('M%.4d',n_M_use);
str_eig = sprintf('eig_from_reco_true_%s',str_M);
%%%%%%%%;
% Now try an eigenvalue estimate with empirical images and empirical viewing-angles. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/%s.mat',dir_ssnll,str_eig);
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat); tmp_lanczos_n_iteration_max = numel(tmp_.alph_i_); clear tmp_;
disp(sprintf(' %% %s found, tmp_lanczos_n_iteration_max %.2d',fname_mat,tmp_lanczos_n_iteration_max));
end;%if ( exist(fname_mat,'file'));
if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
disp(sprintf(' %% %s not complete, calculating',fname_mat));
%%%%;
a_k_Y_use_yk_ = a_k_Y_reco_true_yk_;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_p_use_ ...
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
,a_k_Y_use_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
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
tmp_lanczos_n_iteration_max = 0;
euler_polar_a_use_M_ = euler_polar_a_true_(1:n_M_use);
euler_azimu_b_use_M_ = euler_azimu_b_true_(1:n_M_use);
euler_gamma_z_use_M_ = euler_gamma_z_true_(1:n_M_use);
image_delta_x_use_M_ = image_delta_x_true_(1:n_M_use);
image_delta_y_use_M_ = image_delta_y_true_(1:n_M_use);
M_use_k_p_wkM__ = M_k_p_wkM__(:,1:n_M_use);
for nM_use=0:n_M_use-1;
M_use_k_p_wkM__(:,1+nM_use) = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_use_k_p_wkM__(:,1+nM_use),+image_delta_x_use_M_(1+nM_use),+image_delta_y_use_M_(1+nM_use));
end;%for nM_use=0:n_M-1;
if (~exist(fname_mat,'file'));
v_ykabci__=[];
w_ykabc_=[];
alph_i_=[];
beta_i_=[];
end;%if (~exist(fname_mat,'file'));
%%%%;
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat);
tmp_lanczos_n_iteration_max = numel(tmp_.alph_i_);
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
,a_k_Y_use_yk_ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_use_ ...
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
,n_M_use ...
,M_use_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_use_M_ ...
,euler_azimu_b_use_M_ ...
,euler_gamma_z_use_M_ ...
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
     ,'n_M_use' ...
     ,'euler_polar_a_use_M_' ...
     ,'euler_azimu_b_use_M_' ...
     ,'euler_gamma_z_use_M_' ...
     ,'image_delta_x_use_M_' ...
     ,'image_delta_y_use_M_' ...
     );
%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_calc;