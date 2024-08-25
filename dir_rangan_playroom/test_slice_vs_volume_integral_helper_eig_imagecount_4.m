%%%%%%%%;
% intended for use with test_slice_vs_volume_integral_shell_10.m ;
%%%%%%%%;

n_3 = 3;
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

lanczos_n_iteration_max = 32;

%%%%%%%%;
% Use a_k_Y_quad_yk_ and associated S_use_k_p_wkS__ (bounded away from poles). ;
% Try an eigenvalue estimate with noiseless images and uniform viewing-angles. ;
%%%%%%%%;
[ ...
 n_viewing_S_use ...
,viewing_azimu_b_S_use_ ...
,viewing_polar_a_S_use_ ...
,viewing_weight_S_use_ ...
,viewing_k_c_0_S_use_ ...
,viewing_k_c_1_S_use_ ...
,viewing_k_c_2_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,0.5/k_p_r_max ...
,'C' ...
,1 ...
) ;
%%%%;
tmp_t = tic();
[ ...
 S_use_k_p_wkS__ ...
,~ ...
,n_S_use ...
,viewing_azimu_b_S_use_ ...
,viewing_polar_a_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_quad_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_viewing_S_use ...
,viewing_azimu_b_S_use_ ...
,viewing_polar_a_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
);
S_use_k_p_wkS__ = reshape(S_use_k_p_wkS__,[n_w_sum,n_S_use]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_use_k_p_wkS__ (pm_template_2): %0.6fs',tmp_t)); end;
%%%%;
tmp_index_ = efind( abs(viewing_polar_a_S_use_-0*pi)>1e-6 & abs(viewing_polar_a_S_use_-1*pi)>1e-6 );
n_M_use = numel(tmp_index_);
viewing_weight_M_use_ = viewing_weight_S_use_(1+tmp_index_);
viewing_weight_M_use_ = viewing_weight_M_use_*sum(viewing_weight_S_use_)/max(1e-12,sum(viewing_weight_M_use_));
M_use_k_p_wkM__ = S_use_k_p_wkS__(:,1+tmp_index_);
euler_polar_a_M_use_ = viewing_polar_a_S_use_(1+tmp_index_);
euler_azimu_b_M_use_ = viewing_azimu_b_S_use_(1+tmp_index_);
euler_gamma_z_M_use_ = zeros(n_M_use,1);

flag_calc=1;
if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
lsigma_dist_ = -3:0.5:+3; n_lsigma_dist = numel(lsigma_dist_);
factor_imagecount_Ms__ = zeros(n_M_use,n_lsigma_dist);
for nlsigma_dist=0:n_lsigma_dist-1;
lsigma_dist = lsigma_dist_(1+nlsigma_dist);
sigma_dist = exp(-abs(lsigma_dist));
if abs(lsigma_dist)<1e-12; factor_imagecount_M_use_ = 1/sqrt(2*pi)^2; end;
if lsigma_dist< 0;
factor_imagecount_M_use_ = 1/sqrt(2*pi)^2/sigma_dist^2*exp(-min((euler_polar_a_M_use_-0).^2,(euler_polar_a_M_use_-pi).^2)/(2*sigma_dist^2));
end;%if lsigma_dist< 0;
if lsigma_dist> 0;
factor_imagecount_M_use_ = 1/sqrt(2*pi)^2/sigma_dist^2*exp(-(euler_polar_a_M_use_-pi/2).^2/(2*sigma_dist^2));
end;%if lsigma_dist> 0;
tmp_f = sum(factor_imagecount_M_use_.*viewing_weight_M_use_);
factor_imagecount_M_use_ = factor_imagecount_M_use_./max(1e-12,tmp_f);
factor_imagecount_Ms__(:,1+nlsigma_dist) = factor_imagecount_M_use_;
end;%for nlsigma_dist=0:n_lsigma_dist-1;
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = ceil(n_lsigma_dist/p_row); np=0;
%flim_ = prctile(factor_imagecount_Ms__,[  0,100],'all');
flim_ = [0,2.0/(4*pi)];
flag_2d_vs_3d=0;
for nlsigma_dist=0:n_lsigma_dist-1;
lsigma_dist = lsigma_dist_(1+nlsigma_dist);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0( ...
 euler_polar_a_M_use_ ... 
,euler_azimu_b_M_use_ ... 
,factor_imagecount_Ms__(:,1+nlsigma_dist) ...
,flim_ ... 
,colormap_beach ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axisnotick3d; axis equal; axis vis3d;
title(sprintf('lsigma %0.2f',lsigma_dist),'Interpreter','none');
end;%for nlsigma_dist=0:n_lsigma_dist-1;
drawnow();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nlsigma_dist=0:n_lsigma_dist-1;
for flag_implicit_dtau = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
lsigma_dist = lsigma_dist_(1+nlsigma_dist);
factor_imagecount_M_use_ = factor_imagecount_Ms__(:,1+nlsigma_dist);
weight_imagecount_M_use_ = viewing_weight_M_use_ .* factor_imagecount_M_use_ ;
if lsigma_dist< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma_dist))); end;
if lsigma_dist>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma_dist))); end;
if (flag_verbose>0); disp(sprintf(' %% nlsigma_dist %.2d/%.2d %+0.2f %s',nlsigma_dist,n_lsigma_dist,lsigma_dist,str_infix)); end;
str_dir_mat = sprintf('%s_mat',dir_ssnll);
str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if flag_implicit_dtau==0; str_fname_nopath_prefix = sprintf('eig_from_synth_%s',str_infix); end;
if flag_implicit_dtau==1; str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s',str_infix); end;
%%%%;
fname_mat = sprintf('%s/%s.mat',str_dir_mat,str_fname_nopath_prefix);
tmp_lanczos_n_iteration_max = 0;
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat);
if isfield(tmp_,'alph_i_'); tmp_lanczos_n_iteration_max = numel(tmp_.alph_i_); end;
if isfield(tmp_,'alph_tilde_i_'); tmp_lanczos_n_iteration_max = numel(tmp_.alph_tilde_i_); end;
clear tmp_;
disp(sprintf(' %% %s found, tmp_lanczos_n_iteration_max %.2d',fname_mat,tmp_lanczos_n_iteration_max));
end;%if ( exist(fname_mat,'file'));
if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
disp(sprintf(' %% %s not complete, calculating',fname_mat));
parameter_eig = struct('type','eig');
parameter_eig.flag_verbose = flag_verbose;
parameter_eig.flag_implicit_dtau = flag_implicit_dtau;
parameter_eig.flag_check = 1;
parameter_eig.flag_disp = 1;
parameter_eig.flag_kernel_qpro_d0 = 1;
parameter_eig.flag_kernel_qpro_d1 = 1;
parameter_eig.kernel_qpro_polar_a_pole_north=KAPPA_pole_north_double;
parameter_eig.kernel_qpro_polar_a_pole_south=KAPPA_pole_south_double;
parameter_eig.kernel_qpro_qref_k_eq_d_double=KAPPA_qref_k_eq_d_double;
parameter_eig.kernel_qpro_l_max_use = l_max;
parameter_eig.lanczos_n_iteration_max = lanczos_n_iteration_max;
U_SmallRotation_Delta_ykabc3__ = []; %<-- construct internally. ;
U_tilde_SmallRotation_Delta_ykabc3__ = []; %<-- construct internally. ;
%%%%;
if ( flag_recalc | ~exist(fname_mat,'file'));
tmp_lanczos_n_iteration_max = 0;
rng(0);
v_tilde_ykabci__=[];
w_tilde_ykabc_=[];
alph_tilde_i_=[];
beta_tilde_i_=[];
end;%if (~exist(fname_mat,'file'));
%%%%;
if (~flag_recalc &  exist(fname_mat,'file'));
tmp_ = load(fname_mat);
tmp_lanczos_n_iteration_max = numel(tmp_.alph_tilde_i_);
U_tilde_SmallRotation_Delta_ykabc3__ = tmp_.U_tilde_SmallRotation_Delta_ykabc3__;
v_tilde_ykabci__ = tmp_.v_tilde_ykabci__;
w_tilde_ykabc_ = tmp_.w_tilde_ykabc_;
alph_tilde_i_ = tmp_.alph_tilde_i_;
beta_tilde_i_ = tmp_.beta_tilde_i_;
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
if ~exist('U_tilde_SmallRotation_Delta_ykabc3__','var'); U_tilde_SmallRotation_Delta_ykabc3__=[]; end;
%%%%;
%if ~exist('v_tilde_ykabci__ ','var'); v_tilde_ykabci__ =[]; end;
%if ~exist('w_tilde_ykabc_ ','var'); w_tilde_ykabc_ =[]; end;
%if ~exist('alph_tilde_i_','var'); alph_tilde_i_=[]; end;
%if ~exist('beta_tilde_i_ ','var'); beta_tilde_i_ =[]; end;
%%%%;
[ ...
 parameter_eig ...
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
] = ...
eig_ddssnll_lanczos_2( ...
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
,a_k_p_form_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S_use ...
,S_use_k_p_wkS__ ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
,n_M_use ...
,weight_imagecount_M_use_ ...
,M_use_k_p_wkM__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,euler_gamma_z_M_use_ ...
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
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
);
save(fname_mat ...
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
     ,'n_S_use' ...
     ,'S_use_k_p_wkS__' ...
     ,'viewing_polar_a_S_use_' ...
     ,'viewing_azimu_b_S_use_' ...
     ,'viewing_weight_S_use_' ...
     ,'n_viewing_polar_a_use' ...
     ,'viewing_polar_a_use_' ...
     ,'n_viewing_azimu_b_use_' ...
     ,'U_tilde_SmallRotation_Delta_ykabc3__' ...
     ,'v_tilde_ykabci__'  ...
     ,'w_tilde_ykabc_' ...
     ,'alph_tilde_i_' ...
     ,'beta_tilde_i_' ... 
     ,'n_M_use' ...
     ,'viewing_weight_M_use_' ...
     ,'weight_imagecount_M_use_' ...
     ,'factor_imagecount_M_use_' ...
     ,'M_use_k_p_wkM__' ...
     ,'euler_polar_a_M_use_' ...
     ,'euler_azimu_b_M_use_' ...
     ,'euler_gamma_z_M_use_' ...
     );
%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file') | tmp_lanczos_n_iteration_max< lanczos_n_iteration_max);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for flag_implicit_dtau = 1;
end;%for nlsigma_dist=0:n_lsigma_dist-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%disp(sprintf(' %% returning before diagnostic')); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nlsigma_dist=0:n_lsigma_dist-1;
for flag_implicit_dtau = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
lsigma_dist = lsigma_dist_(1+nlsigma_dist);
factor_imagecount_M_use_ = factor_imagecount_Ms__(:,1+nlsigma_dist);
weight_imagecount_M_use_ = viewing_weight_M_use_ .* factor_imagecount_M_use_ ;
if lsigma_dist< -1e-12; str_infix = sprintf('lsigma_n%.3d',fix(100*abs(lsigma_dist))); end;
if lsigma_dist>=-1e-12; str_infix = sprintf('lsigma_p%.3d',fix(100*abs(lsigma_dist))); end;
if (flag_verbose>0); disp(sprintf(' %% nlsigma_dist %.2d/%.2d %+0.2f %s',nlsigma_dist,n_lsigma_dist,lsigma_dist,str_infix)); end;
str_dir_mat = sprintf('%s_mat',dir_ssnll);
str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if flag_implicit_dtau==0; str_fname_nopath_prefix = sprintf('eig_from_synth_%s',str_infix); end;
if flag_implicit_dtau==1; str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s',str_infix); end
fname_mat = sprintf('%s/%s.mat',str_dir_mat,str_fname_nopath_prefix);
if ~exist(fname_mat,'file');
disp(sprintf(' %% Warning, %s not found, not running diagnostic',fname_mat));
end;%if ~exist(fname_mat,'file');
if  exist(fname_mat,'file');
tmp_ = load(fname_mat);
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
if ~exist('U_tilde_SmallRotation_Delta_ykabc3__','var'); U_tilde_SmallRotation_Delta_ykabc3__=[]; end;
parameter_eig_diagnostic = struct('type','eig_diagnostic');
parameter_eig_diagnostic.flag_verbose = flag_verbose;
parameter_eig_diagnostic.flag_implicit_dtau = flag_implicit_dtau;
parameter_eig_diagnostic.flag_check = 1;
parameter_eig_diagnostic.flag_disp = 1;
parameter_eig_diagnostic.flag_replot = 1;
parameter_eig_diagnostic.flag_kernel_qpro_d0 = 1;
parameter_eig_diagnostic.flag_kernel_qpro_d1 = 1;
parameter_eig_diagnostic.kernel_qpro_polar_a_pole_north=KAPPA_pole_north_double;
parameter_eig_diagnostic.kernel_qpro_polar_a_pole_south=KAPPA_pole_south_double;
parameter_eig_diagnostic.kernel_qpro_qref_k_eq_d_double=KAPPA_qref_k_eq_d_double;
parameter_eig_diagnostic.kernel_qpro_l_max_use = l_max;
parameter_eig_diagnostic.lanczos_n_iteration_max = lanczos_n_iteration_max;
eig_ddssnll_lanczos_diagnostic_2( ...
 parameter_eig_diagnostic ...
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
,a_k_p_form_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S_use ...
,S_use_k_p_wkS__ ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
,n_M_use ...
,weight_imagecount_M_use_ ...
,M_use_k_p_wkM__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,euler_gamma_z_M_use_ ...
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
,n_x_u_pack ...
,tmp_.U_tilde_SmallRotation_Delta_ykabc3__ ...
,tmp_.v_tilde_ykabci__  ...
,tmp_.w_tilde_ykabc_  ...
,tmp_.alph_tilde_i_ ...
,tmp_.beta_tilde_i_ ... 
,str_dir_mat ...
,str_dir_jpg ...
,str_fname_nopath_prefix ...
);
clear tmp_;
end;%if  exist(fname_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for flag_implicit_dtau = 1;
end;%for nlsigma_dist=0:n_lsigma_dist-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp(sprintf(' %% returning after diagnostic')); return;
