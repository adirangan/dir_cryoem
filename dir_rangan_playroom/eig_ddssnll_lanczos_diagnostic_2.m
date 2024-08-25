function ...
[ ...
 parameter ...
] = ...
eig_ddssnll_lanczos_diagnostic_2( ...
 parameter ...
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
,n_M ...
,weight_imagecount_M_ ...
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
,n_x_u ...
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
,str_dir_mat ...
,str_dir_jpg ...
,str_fname_nopath_prefix ...
);

str_thisfunction = 'eig_ddssnll_lanczos_diagnostic_2';

if nargin<1;
disp(sprintf(' %% testing %s (not implemented)',str_thisfunction));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_yk_=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_all_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_k_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); a_k_p_quad_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_q2d_wkS__=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
if (nargin<1+na); weight_imagecount_M_=[]; end; na=na+1;
if (nargin<1+na); M_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); n_CTF=[]; end; na=na+1;
if (nargin<1+na); index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_r_kC__=[]; end; na=na+1;
if (nargin<1+na); CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); n_eta=[]; end; na=na+1;
if (nargin<1+na); index_neta_from_nM_=[]; end; na=na+1;
if (nargin<1+na); eta_k_p_r_ke__=[]; end; na=na+1;
if (nargin<1+na); eta_k_p_wke__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_=[]; end; na=na+1;
if (nargin<1+na); KAPPA=[]; end; na=na+1;
if (nargin<1+na); Ylm_uklma___=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_sub_uka__=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_sub_uka__=[]; end; na=na+1;
if (nargin<1+na); l_max_uk_=[]; end; na=na+1;
if (nargin<1+na); index_nu_n_k_per_shell_from_nk_p_r_=[]; end; na=na+1;
if (nargin<1+na); index_k_per_shell_uka__=[]; end; na=na+1;
if (nargin<1+na); V_lmm___=[]; end; na=na+1;
if (nargin<1+na); L_lm__=[]; end; na=na+1;
if (nargin<1+na); d0W_betazeta_mlma____=[]; end; na=na+1;
if (nargin<1+na); d1W_betazeta_mlma____=[]; end; na=na+1;
if (nargin<1+na); d2W_betazeta_mlma____=[]; end; na=na+1;
if (nargin<1+na); n_x_u=[]; end; na=na+1;
if (nargin<1+na); U_tilde_SmallRotation_Delta_ykabc3__=[]; end; na=na+1;
if (nargin<1+na); v_tilde_ykabci__=[]; end; na=na+1;
if (nargin<1+na); w_tilde_ykabc_=[]; end; na=na+1;
if (nargin<1+na); alph_tilde_i_=[]; end; na=na+1;
if (nargin<1+na); beta_tilde_i_=[]; end; na=na+1; 
if (nargin<1+na); str_dir_mat=[]; end; na=na+1; 
if (nargin<1+na); str_dir_jpg=[]; end; na=na+1; 
if (nargin<1+na); str_fname_nopath_prefix=[]; end; na=na+1; 

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_surrogate')); parameter.flag_surrogate = 0; end; %<-- parameter_bookmark. ;
flag_surrogate = parameter.flag_surrogate;
if (~isfield(parameter,'flag_ignore_U')); parameter.flag_ignore_U = 0; end; %<-- parameter_bookmark. ;
flag_ignore_U = parameter.flag_ignore_U;
if (~isfield(parameter,'flag_ignore_tau')); parameter.flag_ignore_tau = 0; end; %<-- parameter_bookmark. ;
flag_ignore_tau = parameter.flag_ignore_tau;
if (~isfield(parameter,'flag_implicit_dtau')); parameter.flag_implicit_dtau = 0; end; %<-- parameter_bookmark. ;
flag_implicit_dtau = parameter.flag_implicit_dtau;
if (~isfield(parameter,'flag_check')); parameter.flag_check = 0; end; %<-- parameter_bookmark. ;
flag_check = parameter.flag_check;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
if (~isfield(parameter,'flag_replot')); parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
flag_replot = parameter.flag_replot; nf=0;
if (~isfield(parameter,'randseed')); parameter.randseed = 1; end; %<-- parameter_bookmark. ;
randseed = parameter.randseed; rng(randseed);
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
if (~isfield(parameter,'viewing_k_eq_d')); parameter.viewing_k_eq_d = []; end; %<-- parameter_bookmark. ;
viewing_k_eq_d = parameter.viewing_k_eq_d;
if (~isfield(parameter,'template_k_eq_d')); parameter.template_k_eq_d = -1; end; %<-- parameter_bookmark. ;
template_k_eq_d = parameter.template_k_eq_d;
if (~isfield(parameter,'n_order')); parameter.n_order = 3; end; %<-- parameter_bookmark. ;
n_order = parameter.n_order;
if (~isfield(parameter,'flag_kernel_qpro_d0')); parameter.flag_kernel_qpro_d0 = 1; end; %<-- parameter_bookmark. ;
flag_kernel_qpro_d0 = parameter.flag_kernel_qpro_d0;
if (~isfield(parameter,'flag_kernel_qpro_d1')); parameter.flag_kernel_qpro_d1 = 1; end; %<-- parameter_bookmark. ;
flag_kernel_qpro_d1 = parameter.flag_kernel_qpro_d1;
if (~isfield(parameter,'kernel_qpro_polar_a_pole_north')); parameter.kernel_qpro_polar_a_pole_north = 4.5*pi/24; end; %<-- parameter_bookmark. ;
kernel_qpro_polar_a_pole_north = parameter.kernel_qpro_polar_a_pole_north;
if (~isfield(parameter,'kernel_qpro_polar_a_pole_south')); parameter.kernel_qpro_polar_a_pole_south = 3.5*pi/24; end; %<-- parameter_bookmark. ;
kernel_qpro_polar_a_pole_south = parameter.kernel_qpro_polar_a_pole_south;
if (~isfield(parameter,'kernel_qpro_l_max_use')); parameter.kernel_qpro_l_max_use = max(l_max_); end; %<-- parameter_bookmark. ;
kernel_qpro_l_max_use = parameter.kernel_qpro_l_max_use;
if (~isfield(parameter,'lanczos_n_iteration_max')); parameter.lanczos_n_iteration_max = 8; end; %<-- parameter_bookmark. ;
lanczos_n_iteration_max = parameter.lanczos_n_iteration_max;
lanczos_n_iteration_pre = size(v_tilde_ykabci__,2);
lanczos_n_iteration_cur = max(0,lanczos_n_iteration_max-lanczos_n_iteration_pre);
if (~isfield(parameter,'tolerance_dvol')); parameter.tolerance_dvol = 1e-1; end; %<-- parameter_bookmark. ;
tolerance_dvol = parameter.tolerance_dvol;
if (~isfield(parameter,'tolerance_fish')); parameter.tolerance_fish = 1e-1; end; %<-- parameter_bookmark. ;
tolerance_fish = parameter.tolerance_fish;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

str_dir_mat_sub = sprintf('%s/dir_%s',str_dir_mat,str_fname_nopath_prefix);
if ~exist(str_dir_mat_sub,'dir'); disp(sprintf(' %% mkdir %s',str_dir_mat_sub)); mkdir(str_dir_mat_sub); end;
str_dir_jpg_sub = sprintf('%s/dir_%s',str_dir_jpg,str_fname_nopath_prefix);
if ~exist(str_dir_jpg_sub,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg_sub)); mkdir(str_dir_jpg_sub); end;

if isempty(weight_imagecount_M_); weight_imagecount_M_=ones(n_M,1); end;
n_3 = 3;
n_M_use = n_M;
weight_imagecount_M_use_ = weight_imagecount_M_;
if flag_implicit_dtau;
if (flag_verbose>0); disp(sprintf(' %% flag_implicit_dtau %d, setting n_M_use = 0',flag_implicit_dtau)); end;
n_M_use = 0;
weight_imagecount_M_use_ = [];
end;%if flag_implicit_dtau;

%%%%%%%%;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
%%%%%%%%;
k_c_0_all_ = cos(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_1_all_ = sin(k_p_azimu_b_all_).*sin(k_p_polar_a_all_);
k_c_2_all_ = cos(k_p_polar_a_all_);
%%%%%%%%;
if isempty(n_x_u); n_x_u = 64; end;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u)^3;
%%%%%%%%;

scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
%%%%%%%%;
% Set up weights. ;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
weight_3d_riesz_k_all_ = weight_3d_k_all_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_all_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_all_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_k_all_(1+tmp_index_) = weight_3d_k_all_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_k_all_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
Y_k_val_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
tmp_l_val_ = zeros(n_lm_(1+nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(1+nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_index_) = tmp_l_val_;
Y_m_val_(1+tmp_index_) = tmp_m_val_;
Y_k_val_(1+tmp_index_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
weight_3d_riesz_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
weight_3d_riesz_yk_(1+tmp_index_) = weight_3d_riesz_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
n_ykabc = n_lm_sum + n_M_use*3;
weight_3d_riesz_weight_imagecount_ykabc_ = cat(1,weight_3d_riesz_yk_/scaling_volumetric,repmat(weight_imagecount_M_use_,[3,1]));
numerator_root_weight_3d_riesz_weight_imagecount_ykabc_ = reshape(sqrt(weight_3d_riesz_weight_imagecount_ykabc_),[n_ykabc,1]);
denomator_root_weight_3d_riesz_weight_imagecount_ykabc_ = 1./max(1e-12,reshape(sqrt(weight_3d_riesz_weight_imagecount_ykabc_),[n_ykabc,1]));
%%%%%%%%;

if (flag_disp>0);
%%%%%%%%;
% Visualize a_k_Y_quad_yk_. ;
%%%%%%%%;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_p_quad_ ...
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
,a_k_Y_quad_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_quad_ time %0.2fs',tmp_t));
%%%%;
fname_fig_pre = sprintf('%s/%s_a_k_Y_quad_yk_FIGA',str_dir_jpg_sub,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[05:10:95];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),a_x_u_quad_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
n_iteration = numel(alph_tilde_i_); T_tilde__ = real(spdiags([circshift(beta_tilde_i_,-1),alph_tilde_i_,beta_tilde_i_],[-1,0,+1],n_iteration,n_iteration));
lambda_xi__ = -Inf*ones(n_iteration,n_iteration); %<-- (index_lambda,niteration);
for niteration=0:n_iteration-1;
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
lambda_sub_ = eigs(T_tilde_sub__,[],1+niteration);
lambda_xi__(1:1+niteration,1+niteration) = sort(lambda_sub_,'ascend');
end;%for niteration=0:n_iteration-1;
S_x_ = sort(eigs(T_tilde__,[],n_iteration),'ascend');
n_iteration = n_iteration;
S_x_ = S_x_;
S_x_min = min(S_x_);
S_x_max = max(S_x_);
%%%%;
vv_n4__ = zeros(n_iteration,4);
lambda_ns__ = zeros(n_iteration,n_iteration); %<-- (niteration,index_lamda);
ee_ns4___ = zeros(n_iteration,n_iteration,4);
ff_ns4___ = zeros(n_iteration,n_iteration,4);
for niteration=0:n_iteration-1;
v_tilde_ykabc_ = v_tilde_ykabci__(:,1+niteration);
[v_tilde_dvol_yk_,v_tilde_polar_a_M_use_,v_tilde_azimu_b_M_use_,v_tilde_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_tilde_ykabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M_use,v_tilde_ykabc_,v_tilde_ykabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_vv)); end;
vv_n4__(1+niteration,:) = [tmp_vv_dvol;tmp_vv_a;tmp_vv_b;tmp_vv_c];
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
[TV_tilde_sub__,lambda_sub__] = eigs(T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
for index_lambda=0:1+niteration-1;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
lambda_ns__(1+niteration,1+index_lambda) = lambda_use;
TV_tilde_eig_ = TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ee)); end;
ee_ns4___(1+niteration,1+index_lambda,:) = [tmp_ee_dvol;tmp_ee_a;tmp_ee_b;tmp_ee_c];
v_eig_ykabc_ = bsxfun(@times,denomator_root_weight_3d_riesz_weight_imagecount_ykabc_,v_tilde_eig_ykabc_);
[v_eig_dvol_yk_,v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_eig_ykabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_imagecount_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M_use,weight_imagecount_M_use_,v_eig_ykabc_,v_eig_ykabc_);
str_ff = sprintf('tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ff)); end;
ff_ns4___(1+niteration,1+index_lambda,:) = [tmp_ff_dvol;tmp_ff_a;tmp_ff_b;tmp_ff_c];
end;%for index_lambda=0:1+niteration-1;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ee_ns4___ vs ff_ns4___: %0.16f',fnorm(ee_ns4___ - ff_ns4___)/max(1e-12,fnorm(ee_ns4___)))); end;

%%%%%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s/%s_FIGA',str_dir_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figmed;fig81s;
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
%%%%;
hold on;
plot(repmat([0;n_iteration],[1,n_iteration]),repmat(reshape(S_x_,[1,n_iteration]),[2,1]),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_sml);
ni_xi__ = repmat([1:n_iteration],[n_iteration,1]);
tmp_index_ = efind(isfinite(lambda_xi__));
plot(ni_xi__(1+tmp_index_),lambda_xi__(1+tmp_index_),'r.','MarkerSize',markersize_use);
hold off;
xlabel('iteration'); ylabel('sigma');
xlim([0,1+n_iteration]);
ylim([S_x_min-0.25,S_x_max+0.25]);
%title(str_fname_nopath_prefix,'Interpreter','none');
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
if (flag_disp>2);
fname_fig_pre = sprintf('%s/%s_FIGB',str_dir_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);fig80s();
p_row = 1; p_col = 5; np=0;
fontsize_use = 12;
ilim_ = [-0.125,+1.125];
%%%%;
for pcol=0:4-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(squeeze(ee_ns4___(:,:,1+pcol)),ilim_);
xlabel('index'); ylabel('iteration'); %axisnotick;
title(sprintf('ee_ns4___(:,:,1+%d)',pcol),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
end;%for pcol=0:4-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(vv_n4__,ilim_);
xlabel('1-4'); ylabel('iteration'); %axisnotick;
title('vv_n4__','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
end;%if (flag_disp>2);
%%%%%%%%;

%%%%%%%%;
% Now scan for relevant eigenvectors. ;
%%%%%%%%;
lambda_nrm_ns__ = lambda_ns__./max(1e-12,max(abs(lambda_ns__),[],'all'));
index_dvol_ = efind(ee_ns4___(:,:,1+0)>=tolerance_dvol);
index_fish_ = efind(lambda_nrm_ns__(:,:)<=tolerance_fish);
index_scan_ = intersect(index_dvol_,index_fish_);
n_scan = numel(index_scan_);
if (flag_verbose>1);
for nscan=0:n_scan-1;
index_scan = index_scan_(1+nscan);
nn = mod(index_scan,n_iteration);
ns = floor(index_scan/n_iteration);
disp(sprintf(' %% nscan %d/%d (%d-->(%d,%d)) lambda %0.6f lambda_nrm %0.6f ee_xx1 %0.6f',nscan,n_scan,index_scan,nn,ns,lambda_ns__(1+nn,1+ns),lambda_nrm_ns__(1+nn,1+ns),ee_ns4___(1+nn,1+ns)));
end;%for nscan=0:n_scan-1;
end;%if (flag_verbose>1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Step through scan.; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nscan=0:n_scan-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
index_scan = index_scan_(1+nscan);
nn = mod(index_scan,n_iteration);
ns = floor(index_scan/n_iteration);
niteration = nn;
index_lambda = ns;
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
[TV_tilde_sub__,lambda_sub__] = eigs(T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
assert(abs(lambda_use-lambda_ns__(1+nn,1+ns))<1e-12);
TV_tilde_eig_ = TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('lambda %0.6f: tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',lambda_use,tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ee)); end;
v_eig_ykabc_ = bsxfun(@times,denomator_root_weight_3d_riesz_weight_imagecount_ykabc_,v_tilde_eig_ykabc_);
[v_eig_dvol_yk_,v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_eig_ykabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_imagecount_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M_use,weight_imagecount_M_use_,v_eig_ykabc_,v_eig_ykabc_);
str_ff = sprintf('lambda %0.6f: tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',lambda_use,tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ff)); end;
assert(abs(tmp_ff_dvol-ee_ns4___(1+nn,1+ns,1+0))<1e-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nscan=0:n_scan-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now analyze specific eigenvectors. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
niteration = n_iteration-1;
T_tilde_sub__ = T_tilde__(1:1+niteration,1:1+niteration);
[TV_tilde_sub__,lambda_sub__] = eigs(T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
%%%%%%%%%%%%%%%%;
for index_lambda = 0:niteration; %<-- 0 is most negative, niteration is most positive. ;
%%%%%%%%%%%%%%%%;
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,niteration,index_lambda);
fname_sub_mat = sprintf('%s/%s.mat',str_dir_mat_sub,str_fname_nopath_sub_prefix);
%%%%%%%%;
if ~exist(fname_sub_mat,'file');
%%%%;
% use T_tilde_sub__ to estimate minimum eigenvector. ;
%%%%;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
TV_tilde_eig_ = TV_tilde_sub__(:,ij_use);
v_tilde_eig_ykabc_ = v_tilde_ykabci__(:,1:1+niteration)*TV_tilde_eig_;
[v_tilde_eig_dvol_yk_,v_tilde_eig_polar_a_M_use_,v_tilde_eig_azimu_b_M_use_,v_tilde_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M_use,v_tilde_eig_ykabc_,v_tilde_eig_ykabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;
v_eig_ykabc_ = bsxfun(@times,denomator_root_weight_3d_riesz_weight_imagecount_ykabc_,v_tilde_eig_ykabc_);
[v_eig_dvol_yk_,v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_eig_ykabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_imagecount_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M_use,weight_imagecount_M_use_,v_eig_ykabc_,v_eig_ykabc_);
str_ff = sprintf('tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ff)); end;
%%%%;
% Now evaluate ddssnll and associated rayleigh-quotient. ;
%%%%;
parameter_ddssnll = parameter;
parameter_ddssnll.flag_verbose = flag_verbose;
parameter_ddssnll.flag_implicit_dtau = flag_implicit_dtau;
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
[ ...
 ~ ...
,tmp_Hvt_ykabc_ ...
,tmp_Hv_q3d_k_Y_quad_yk_ ...
,tmp_Hv_q3d_k_Y_quad_yk__ ...
,tmp_Hv_q3d_k_p_quad_ ...
,tmp_Ht_q2d_M3__ ...
,tmp_a_restore_C2M0_k_Y_lmk_ ...
,tmp_a_restore_C2M0_k_p_quad_ ...
,tmp_Hvv_q3d_k_Y_quad_yk_ ...
,tmp_Hvt_q3d_k_Y_quad_yk_ ...
,tmp_Htv_q2d_M3__ ...
,tmp_Htt_q2d_M3__ ...
,tmp_dtau_euler_polar_a_M_ ...
,tmp_dtau_euler_azimu_b_M_ ...
,tmp_dtau_euler_gamma_z_M_ ...
,tmp_n_dvt ... 
,tmp_dvt_ ... 
,tmp_dvt ... 
,tmp_ssnll_tmp_q2d_dvt_ ... 
,tmp_dssnll_mid_q2d ... 
,tmp_dssnll_dif_q2d ... 
,tmp_dssnll_lsq_q2d ... 
,tmp_ddssnll_mid_q2d ... 
,tmp_ddssnll_dif_q2d ... 
,tmp_ddssnll_lsq_q2d ... 
] = ...
ddssnll_2( ...
 parameter_ddssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,[] ...
,v_eig_dvol_yk_ ...
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
,weight_imagecount_M_ ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,CTF_k_p_wkC__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,v_eig_polar_a_M_use_ ...
,v_eig_azimu_b_M_use_ ...
,v_eig_gamma_z_M_use_ ...
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
%%%%;
w_nottilde_ykabc_ = tmp_Hvt_ykabc_(1:n_ykabc); %<-- ignore alignment component if flag_implicit_dtau. ;
w_tilde_ykabc_ = ...
  bsxfun( ...
	  @times ...
	  ,numerator_root_weight_3d_riesz_weight_imagecount_ykabc_ ...
	  ,w_nottilde_ykabc_ ...
	  );
%%%%;
ddssnll_mid_q2d = local_imagecount_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M_use,weight_imagecount_M_use_,w_nottilde_ykabc_,v_eig_ykabc_);
ddssnll_tilde_mid_q2d = local_weightless_f_bar_dot_g_(n_k_p_r,l_max_,n_M_use,w_tilde_ykabc_,v_tilde_eig_ykabc_);
if (flag_verbose>0); disp(sprintf(' %% ddssnll_mid_q2d: %0.16f',ddssnll_mid_q2d)); end;
if (flag_verbose>0); disp(sprintf(' %% ddssnll_tilde_mid_q2d: %0.16f',ddssnll_tilde_mid_q2d)); end;
if (flag_verbose>0); disp(sprintf(' %% ddssnll_mid_q2d vs tmp_ddssnll_mid_q2d: %0.16f',fnorm(ddssnll_mid_q2d - tmp_ddssnll_mid_q2d)/max(1e-12,fnorm(ddssnll_mid_q2d)))); end;
%%%%;
save(fname_sub_mat ...
     ,'n_k_p_r','k_p_r_','weight_3d_k_p_r_','weight_3d_riesz_k_p_r_','l_max_' ...
     ,'n_M','weight_imagecount_M_','n_M_use','weight_imagecount_M_use_' ...
     ,'scaling_volumetric','weight_3d_riesz_weight_imagecount_ykabc_' ...
     ,'niteration','T_tilde_sub__','TV_tilde_sub__','lambda_sub_','lambda_srt_','ij_srt_' ...
     ,'v_tilde_eig_ykabc_','tmp_ee','tmp_ee_dvol','tmp_ee_a','tmp_ee_b','tmp_ee_c' ...
     ,'v_eig_ykabc_','tmp_ff','tmp_ff_dvol','tmp_ff_a','tmp_ff_b','tmp_ff_c' ...
     ,'tmp_Hvt_ykabc_' ...
     ,'tmp_Hv_q3d_k_Y_quad_yk_' ...
     ,'tmp_Ht_q2d_M3__' ...
     ,'tmp_Hvv_q3d_k_Y_quad_yk_' ...
     ,'tmp_Hvt_q3d_k_Y_quad_yk_' ...
     ,'tmp_Htv_q2d_M3__' ...
     ,'tmp_Htt_q2d_M3__' ...
     ,'euler_polar_a_M_' ...
     ,'euler_azimu_b_M_' ...
     ,'euler_gamma_z_M_' ...
     ,'tmp_dtau_euler_polar_a_M_' ...
     ,'tmp_dtau_euler_azimu_b_M_' ...
     ,'tmp_dtau_euler_gamma_z_M_' ...
     ,'tmp_n_dvt' ... 
     ,'tmp_dvt_' ... 
     ,'tmp_dvt' ... 
     ,'tmp_ssnll_tmp_q2d_dvt_' ... 
     ,'tmp_dssnll_mid_q2d' ... 
     ,'tmp_dssnll_dif_q2d' ... 
     ,'tmp_dssnll_lsq_q2d' ... 
     ,'tmp_ddssnll_mid_q2d' ... 
     ,'tmp_ddssnll_dif_q2d' ... 
     ,'tmp_ddssnll_lsq_q2d' ...
     ,'ddssnll_mid_q2d','ddssnll_tilde_mid_q2d' ...
     );
%%%%%%%%;
end;%if ~exist(fname_sub_mat,'file');
%%%%%%%%;
if  exist(fname_sub_mat,'file');
if (flag_verbose>0); disp(sprintf(' %% %s found, not creating',fname_sub_mat)); end;
load(fname_sub_mat);
end;%if  exist(fname_sub_mat,'file');
%%%%%%%%;
str_sgtitle = sprintf('%s niteration %d lambda %+0.2f = exp(%+0.2f) index_lambda %d/%d: \n %s ',str_fname_nopath_prefix,niteration,lambda_use,log(abs(lambda_use)),index_lambda,niteration,str_ee);
%%%%%%%%;
if (flag_disp>0);
%%%%%%%%;
% visualize v_eig_dvol_yk_;
%%%%%%%%;
[v_eig_dvol_yk_,v_eig_polar_a_M_use_,v_eig_azimu_b_M_use_,v_eig_gamma_z_M_use_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M_use,v_eig_ykabc_);
v_k_Y_use_yk_ = v_eig_dvol_yk_;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 v_k_p_use_ ...
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
,v_k_Y_use_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% v_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
v_x_u_use_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,v_k_p_use_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: v_x_u_use_ time %0.2fs',tmp_t));
%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGA',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[05:10:95];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),v_x_u_use_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle(sprintf('v_eig_dvol_yk_: %s niteration %d lambda %+0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_fname_nopath_prefix,niteration,lambda_use,log(abs(lambda_use)),index_lambda,niteration,str_ee),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_sub_jpg));
print('-djpeg',fname_fig_sub_jpg);
end;%if flag_replot | ~exist(fname_fig_sub_jpg,'file');
close(gcf);
end;%if flag_disp;
%%%%%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGB',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
v_x_u_use___ = reshape(v_x_u_use_,[n_x_u,n_x_u,n_x_u]);
subplot(1,3,1);
imagesc(squeeze(real(v_x_u_use___(:,:,end/2))));
axis image; axisnotick; colorbar;
subplot(1,3,2);
imagesc(squeeze(real(v_x_u_use___(:,end/2,:))));
axis image; axisnotick; colorbar;
subplot(1,3,3);
imagesc(squeeze(real(v_x_u_use___(end/2,:,:))));
axis image; axisnotick; colorbar;
sgtitle(str_sgtitle,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_sub_jpg));
print('-djpeg',fname_fig_sub_jpg);
end;%if flag_replot | ~exist(fname_fig_sub_jpg,'file');
close(gcf);
end;%if flag_disp;
%%%%%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGC',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
vmax = max(abs(v_x_u_use_),[],'all');
amax = max(abs(a_x_u_quad_),[],'all');
dvol_ = 1e-1*[-2:+2]; p_col = numel(dvol_);
prct_ = [95,97.5,98.5]; p_row = numel(prct_);
np=0;
for prow=0:p_row-1;
for pcol=0:p_col-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
prct = prct_(1+prow);
dvol = dvol_(1+pcol);
tmp_a_x_u_use_ = a_x_u_quad_/max(1e-12,amax) + dvol*v_x_u_use_/max(1e-12,vmax);
isosurface_f_x_u_1(struct('percent_threshold_',[prct]),tmp_a_x_u_use_);
title(sprintf('prct %0.2f dvol %0.2f',prct,dvol));
end;%for pcol=0:p_col-1;
end;%for prow=0:p_row-1;
sgtitle(str_sgtitle,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_sub_jpg));
print('-djpeg',fname_fig_sub_jpg);
end;%if flag_replot | ~exist(fname_fig_sub_jpg,'file');
close(gcf);
end;%if flag_disp;
%%%%%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGD',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_dtau_euler_polar_a_M_,tmp_dtau_euler_azimu_b_M_,tmp_dtau_euler_gamma_z_M_].^2,weight_imagecount_M_),'all');
tmp_dtau_euler_polar_a_M_ = tmp_dtau_euler_polar_a_M_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_azimu_b_M_ = tmp_dtau_euler_azimu_b_M_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_gamma_z_M_ = tmp_dtau_euler_gamma_z_M_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%{
subplot(1,2,1);
plot_sphere_grid_0;
hold on;
sphere_post__0( ...
 struct('type','sphere_post','post_r_base',1.0/(2*pi)/8) ...
,n_M ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,factor_amplify*tmp_dtau_euler_polar_a_M_ ...
,factor_amplify*tmp_dtau_euler_azimu_b_M_ ...
,factor_amplify*tmp_dtau_euler_gamma_z_M_ ...
);
hold off;
axis equal; axis vis3d; axisnotick3d;
title('v_min_dtau_','Interpreter','none');
%%%%;
subplot(1,2,2);
plot_sphere_grid_0;
hold on;
sphere_compass__0( ...
 struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',0) ...
,n_M ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,tmp_dtau_euler_polar_a_M_ ...
,tmp_dtau_euler_azimu_b_M_ ...
,tmp_dtau_euler_gamma_z_M_ ...
);
hold off;
axis equal; axis vis3d; axisnotick3d;
title('v_min_dtau_','Interpreter','none');
%}
%%%%;
figbig;
subplot(1,1,1);
hold on;
sphere_compass__0( ...
 struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',1) ...
,n_M ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,tmp_dtau_euler_polar_a_M_ ...
,tmp_dtau_euler_azimu_b_M_ ...
,tmp_dtau_euler_gamma_z_M_ ...
);
hold off;
xlim([0,2*pi]); xlabel('azimu_b','Interpreter','none');
ylim([0,1*pi]); ylabel('polar_a','Interpreter','none');
axisnotick;
title('dtau_euler_','Interpreter','none');
%%%%;
sgtitle(str_sgtitle,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_sub_jpg));
print('-djpeg',fname_fig_sub_jpg);
end;%if flag_replot | ~exist(fname_fig_sub_jpg,'file');
close(gcf);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
end;%for index_lambda = 0:niteration; %<-- 0 is most negative, niteration is most positive. ;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Finished analyzing specific eigenvectors. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
end;%if (flag_disp>0);
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

