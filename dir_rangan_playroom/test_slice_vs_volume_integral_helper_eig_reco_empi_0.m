flag_calc=1;
if flag_calc;
tmp_lanczos_n_iteration_max = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_M_use = 1024;
index_nCTF_from_nM_use_ = index_nCTF_from_nM_(1:n_M_use);
str_M = sprintf('M%.4d',n_M_use);
str_eig = sprintf('eig_from_reco_empi_%s',str_M);
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
a_k_Y_use_yk_ = a_k_Y_reco_empi_yk_;
a_k_Y_use_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_use_yk_);
%%%%;
[ ...
 S_k_p_use_wkS__ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_use_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
S_k_p_use_wkS__ = reshape(S_k_p_use_wkS__,[n_w_sum,n_S]);
%%%%;
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
parameter_eig.kernel_qpro_polar_a_pole_north=KAPPA_pole_north_double;
parameter_eig.kernel_qpro_polar_a_pole_south=KAPPA_pole_south_double;
parameter_eig.kernel_qpro_qref_k_eq_d_double=KAPPA_qref_k_eq_d_double;
parameter_eig.kernel_qpro_l_max_use = l_max;
parameter_eig.lanczos_n_iteration_max = lanczos_n_iteration_max;
U_SmallRotation_Delta_ykabc3__ = []; %<-- construct internally. ;
U_SmallRotation_Delta_ykabc3__ = zeros(n_lm_sum + n_M_use*3,0); %<-- skip entirely. ;
S_k_p_q2d_wkS__ = S_k_p_use_wkS__;
%%%%;
tmp_lanczos_n_iteration_max = 0;
euler_polar_a_use_M_ = euler_polar_a_empi_(1:n_M_use);
euler_azimu_b_use_M_ = euler_azimu_b_empi_(1:n_M_use);
euler_gamma_z_use_M_ = euler_gamma_z_empi_(1:n_M_use);
image_delta_x_use_M_ = image_delta_x_empi_(1:n_M_use);
image_delta_y_use_M_ = image_delta_y_empi_(1:n_M_use);
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

%{
  
dvt = 1e-2;
n_M_use = 1024;
index_nCTF_from_nM_use_ = index_nCTF_from_nM_(1:n_M_use);
a_k_Y_use_yk_ = a_k_Y_reco_empi_yk_;
dtau_euler_polar_a_use_M_ = euler_polar_a_use_M_;
dtau_euler_azimu_b_use_M_ = euler_azimu_b_use_M_;
dtau_euler_gamma_z_use_M_ = euler_gamma_z_use_M_;
dv_ = transpose(-5:+5); n_dv = numel(dv_);
dt_ = transpose(-4:+4); n_dt = numel(dt_);
ssnll_q2d_dvdt__ = zeros(n_dv,n_dt);
for ndv=0:n_dv-1;
for ndt=0:n_dt-1;
if (flag_verbose>1); disp(sprintf(' %% ndv %d/%d ndt %d/%d',ndv,n_dv,ndt,n_dt)); end;
[ ...
 ~ ...
,~ ...
,ssnll_q2d_dvdt__(1+ndv,1+ndt) ...
] = ...
ssnll_from_a_k_Y_12( ...
 [] ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_use_yk__ + 1.0*dvt*dv_(1+ndv)*a_k_Y_use_yk__...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_use_wkS__ + 1.0*dvt*dv_(1+ndv)*S_k_p_use_wkS__...
,[] ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M_use ...
,M_use_k_p_wkM__ ...
,index_nCTF_from_nM_use_ ...
,CTF_k_p_wkC__ ...
,[] ...
,[] ...
,euler_polar_a_use_M_ + 1.0*dvt*dt_(1+ndt)*dtau_euler_polar_a_use_M_ ...
,euler_azimu_b_use_M_ + 1.0*dvt*dt_(1+ndt)*dtau_euler_azimu_b_use_M_ ...
,euler_gamma_z_use_M_ + 1.0*dvt*dt_(1+ndt)*dtau_euler_gamma_z_use_M_ ...
,[] ...
,[] ...
,[] ...
);
if (flag_verbose>0); disp(sprintf(' %% ndv %d/%d ndt %d/%d: ssnll: %+0.6fs',ndv,n_dv,ndt,n_dt,ssnll_q2d_dvdt__(1+ndv,1+ndt))); end;
end;%for ndt=0:n_dt-1;
end;%for ndv=0:n_dv-1;
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
imagesc(ssnll_q2d_dvdt__);colorbar;
xlabel('dt'); ylabel('dv');
subplot(1,2,2);
surfl(ssnll_q2d_dvdt__);colorbar;
xlabel('dt'); ylabel('dv');
%%%%%%%%;

%}

