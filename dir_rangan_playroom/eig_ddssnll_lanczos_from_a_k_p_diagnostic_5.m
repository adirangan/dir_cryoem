function ...
[ ...
 parameter ...
] = ...
eig_ddssnll_lanczos_from_a_k_p_diagnostic_5( ...
 parameter ...
,pm_n_qk ...
,pm_n_qk_csum_ ...
,pm_k_p_r_qk_ ...
,pm_k_p_azimu_b_qk_ ...
,pm_k_p_polar_a_qk_ ...
,pm_weight_3d_k_p_qk_ ...
,pm_weight_shell_qk_ ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_3d_k_p_r_ ...
,pm_k_c_0_qk_ ...
,pm_k_c_1_qk_ ...
,pm_k_c_2_qk_ ...
,pm_n_polar_a_k_ ...
,pm_polar_a_ka__ ...
,pm_n_azimu_b_ka__ ...
,qref_k_eq_d ...
,pm_a_k_p_qk_ ...
,pm_n_w_ ...
,pm_weight_2d_k_p_r_ ...
,pm_weight_2d_wk_ ...
,n_viewing_S_use ...
,pm_S_use_k_p_wkS__ ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
,viewing_gamma_z_S_use_ ...
,n_M_use ...
,weight_imagecount_M_use_ ...
,pm_M_use_k_p_wkM__ ...
,pm_n_CTF ...
,pm_index_nCTF_from_nM_ ...
,pm_CTF_k_p_wkC__ ...
,pm_n_eta ...
,pm_index_neta_from_nM_ ...
,pm_eta_k_p_wke__ ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,euler_gamma_z_M_use_ ...
,pm_KAPPA ...
,pm_wS_from_single_shell_sba__ ...
,pm_dwSda_from_single_shell_sba__ ...
,pm_dwSdb_from_single_shell_sba__ ...
,pm_ddwSdaa_from_single_shell_sba__ ...
,pm_ddwSdab_from_single_shell_sba__ ...
,pm_ddwSdbb_from_single_shell_sba__ ...
,pm_n_R ...
,pm_R_use_R___ ...
,pm_a_R_k_p_Rqk__ ...
,pm_ba_from_single_shell_Rbaba___ ...
,pm_wS_from_R_single_shell_Rsba___ ...
,pm_dwSda_from_R_single_shell_Rsba___ ...
,pm_dwSdb_from_R_single_shell_Rsba___ ...
,pm_ddwSdaa_from_R_single_shell_Rsba___ ...
,pm_ddwSdab_from_R_single_shell_Rsba___ ...
,pm_ddwSdbb_from_R_single_shell_Rsba___ ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,from_pm_UX_kn__ ...
,X_weight_r_ ...
,n_x_u_pack ...
,pm_U_tilde_SmallRotation_Delta_qkabc3__ ...
,pm_v_tilde_qkabci__  ...
,pm_w_tilde_qkabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
,str_dir_mat ...
,str_dir_jpg ...
,str_fname_nopath_prefix ...
);

str_thisfunction = 'eig_ddssnll_lanczos_from_a_k_p_diagnostic_5';

if nargin<1;
disp(sprintf(' %% testing %s (not implemented)',str_thisfunction));
disp('returning'); return;
end;%if nargin<1;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); pm_n_qk=[]; end; na=na+1;
if (nargin<1+na); pm_n_qk_csum_=[]; end; na=na+1;
if (nargin<1+na); pm_k_p_r_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_k_p_azimu_b_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_k_p_polar_a_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_weight_3d_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_weight_shell_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); pm_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); pm_k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); pm_weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); pm_k_c_0_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_k_c_1_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_k_c_2_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_n_polar_a_k_=[]; end; na=na+1;
if (nargin<1+na); pm_polar_a_ka__=[]; end; na=na+1;
if (nargin<1+na); pm_n_azimu_b_ka__=[]; end; na=na+1;
if (nargin<1+na); qref_k_eq_d=[]; end; na=na+1;
if (nargin<1+na); pm_a_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); pm_n_w_=[]; end; na=na+1;
if (nargin<1+na); pm_weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); pm_weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_S_use=[]; end; na=na+1;
if (nargin<1+na); pm_S_use_k_p_wkS__=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_use_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_use_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_use_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a_use=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_use_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_use_=[]; end; na=na+1;
if (nargin<1+na); viewing_gamma_z_S_use_=[]; end; na=na+1;
if (nargin<1+na); n_M_use=[]; end; na=na+1;
if (nargin<1+na); weight_imagecount_M_use_=[]; end; na=na+1;
if (nargin<1+na); pm_M_use_k_p_wkM__=[]; end; na=na+1;
if (nargin<1+na); pm_n_CTF=[]; end; na=na+1;
if (nargin<1+na); pm_index_nCTF_from_nM_=[]; end; na=na+1;
if (nargin<1+na); pm_CTF_k_p_wkC__=[]; end; na=na+1;
if (nargin<1+na); pm_n_eta=[]; end; na=na+1;
if (nargin<1+na); pm_index_neta_from_nM_=[]; end; na=na+1;
if (nargin<1+na); pm_eta_k_p_wke__=[]; end; na=na+1;
if (nargin<1+na); euler_polar_a_M_use_=[]; end; na=na+1;
if (nargin<1+na); euler_azimu_b_M_use_=[]; end; na=na+1;
if (nargin<1+na); euler_gamma_z_M_use_=[]; end; na=na+1;
if (nargin<1+na); pm_KAPPA=[]; end; na=na+1;
if (nargin<1+na); pm_wS_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_dwSda_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_dwSdb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdaa_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdab_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdbb_from_single_shell_sba__=[]; end; na=na+1;
if (nargin<1+na); pm_n_R=[]; end; na=na+1;
if (nargin<1+na); pm_R_use_R___=[]; end; na=na+1;
if (nargin<1+na); pm_a_R_k_p_Rqk__=[]; end; na=na+1;
if (nargin<1+na); pm_ba_from_single_shell_Rbaba___=[]; end; na=na+1;
if (nargin<1+na); pm_wS_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); pm_dwSda_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); pm_dwSdb_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdaa_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdab_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); pm_ddwSdbb_from_R_single_shell_Rsba___=[]; end; na=na+1;
if (nargin<1+na); n_qk=[]; end; na=na+1;
if (nargin<1+na); n_qk_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_qk_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_qk_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_qk_=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); from_pm_UX_kn__=[]; end; na=na+1;
if (nargin<1+na); X_weight_r_=[]; end; na=na+1;
if (nargin<1+na); n_x_u_pack=[]; end; na=na+1;
if (nargin<1+na); pm_U_tilde_SmallRotation_Delta_qkabc3__=[]; end; na=na+1;
if (nargin<1+na); pm_v_tilde_qkabci__ =[]; end; na=na+1;
if (nargin<1+na); pm_w_tilde_qkabc_ =[]; end; na=na+1;
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
if (~isfield(parameter,'flag_recalc')); parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
flag_recalc = parameter.flag_recalc; nf=0;
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
if (~isfield(parameter,'n_order')); parameter.n_order = 5; end; %<-- parameter_bookmark. ;
n_order = parameter.n_order;
%%%%;
if ~isfield(parameter,'kernel_basic_qref_k_eq_d_double'); parameter.kernel_basic_qref_k_eq_d_double=[]; end;
kernel_basic_qref_k_eq_d_double=parameter.kernel_basic_qref_k_eq_d_double;
if (~isfield(parameter,'kernel_basic_l_max_use')); parameter.kernel_basic_l_max_use = round(2*pi*k_p_r_max); end; %<-- parameter_bookmark. ;
kernel_basic_l_max_use = parameter.kernel_basic_l_max_use;
if (~isfield(parameter,'kernel_basic_l_max_ext')); parameter.kernel_basic_l_max_ext = ceil(1.25*kernel_basic_l_max_use); end; %<-- parameter_bookmark. ;
kernel_basic_l_max_ext = parameter.kernel_basic_l_max_ext;
if (~isfield(parameter,'kernel_basic_l_max_band')); parameter.kernel_basic_l_max_band = +Inf; end; %<-- parameter_bookmark. ;
kernel_basic_l_max_band = parameter.kernel_basic_l_max_band;
if (~isfield(parameter,'tolerance_pinv')); parameter.tolerance_pinv = 1e-6; end; %<-- parameter_bookmark. ;
tolerance_pinv = parameter.tolerance_pinv;
%%%%;
if (~isfield(parameter,'lanczos_n_iteration_max')); parameter.lanczos_n_iteration_max = 8; end; %<-- parameter_bookmark. ;
lanczos_n_iteration_max = parameter.lanczos_n_iteration_max;
lanczos_n_iteration_pre = size(pm_v_tilde_qkabci__,2);
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

%%%%%%%%;
%%%%;
if isempty(weight_imagecount_M_use_); weight_imagecount_M_use_=ones(n_M_use,1); end;
n_1 = 1; n_3 = 3; n_9 = 9;
n_M_imp = n_M_use;
weight_imagecount_M_imp_ = weight_imagecount_M_use_;
if flag_implicit_dtau;
if (flag_verbose>0); disp(sprintf(' %% flag_implicit_dtau %d, setting n_M_imp = 0',flag_implicit_dtau)); end;
n_M_imp = 0;
weight_imagecount_M_imp_ = [];
end;%if flag_implicit_dtau;
%%%%%%%%;


pm_n_q = pm_n_qk/max(1,pm_n_k_p_r);
n_q = pm_n_q;
assert(n_qk==n_q*n_k_p_r);
%%%%%%%%;
pm_n_w_max = max(pm_n_w_);
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
%%%%%%%%;
pm_k_c_0_qk_ = cos(pm_k_p_azimu_b_qk_).*sin(pm_k_p_polar_a_qk_);
pm_k_c_1_qk_ = sin(pm_k_p_azimu_b_qk_).*sin(pm_k_p_polar_a_qk_);
pm_k_c_2_qk_ = cos(pm_k_p_polar_a_qk_);
%%%%%%%%;
n_w_ = pm_n_w_max*ones(n_k_p_r,1);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
k_c_0_qk_ = k_p_r_qk_.*cos(k_p_azimu_b_qk_).*sin(k_p_polar_a_qk_);
k_c_1_qk_ = k_p_r_qk_.*sin(k_p_azimu_b_qk_).*sin(k_p_polar_a_qk_);
k_c_2_qk_ = k_p_r_qk_.*cos(k_p_polar_a_qk_);
%%%%%%%%;
if isempty(n_x_u_pack); n_x_u_pack = 64; end;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_0 = mean(diff(x_u_0_));
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_1 = mean(diff(x_u_1_));
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_2 = mean(diff(x_u_2_));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
weight_xxx_u_ = d_x_0 * d_x_1 * d_x_2 ;
%%%%%%%%;

scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
%%%%%%%%;
% Set up weights. ;
%%%%%%%%;
pm_weight_3d_riesz_k_p_r_ = pm_weight_3d_k_p_r_;
pm_weight_3d_riesz_qk_ = pm_weight_3d_k_p_qk_;
%%%%%%%%;
% Now adjust pm_weights according to get_weight_3d_riesz_2. ;
%%%%%%%%;
pm_weight_shell_qk_ = pm_weight_shell_qk_/max(1e-12,sum(pm_weight_shell_qk_))*4*pi*pm_n_k_p_r ;
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)*pi ;
pm_weight_2d_wk_ = ones(pm_n_w_sum,1)/max(1,pm_n_w_sum)/(4*pi) ;
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)/3 ;
pm_weight_3d_k_p_qk_ = ones(pm_n_qk,1)/max(1,pm_n_qk)*((4/3)*pi) ;
pm_weight_3d_riesz_k_p_r_ = ones(pm_n_k_p_r,1)/max(1,pm_n_k_p_r)*pi ;
pm_weight_3d_riesz_qk_ = ones(pm_n_qk,1)/max(1,pm_n_qk)*(4*pi^2) ;
%%%%%%%%;
for pm_nk_p_r=0:pm_n_k_p_r-1;
pm_k_p_r = pm_k_p_r_(1+pm_nk_p_r);
pm_weight_2d_k_p_r = pm_weight_2d_k_p_r_(1+pm_nk_p_r);
pm_weight_3d_k_p_r = pm_weight_3d_k_p_r_(1+pm_nk_p_r);
tmp_index_ = pm_n_qk_csum_(1+pm_nk_p_r):pm_n_qk_csum_(1+pm_nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% pm_k_p_r %0.6f: sum(pm_weight_3d_k_p_qk_(1+tmp_index_))/(pm_weight_3d_k_p_r * 4*pi): %0.16f',pm_k_p_r,sum(pm_weight_3d_k_p_qk_(1+tmp_index_))/(4*pi*pm_weight_3d_k_p_r))); end;
if (flag_verbose> 1); disp(sprintf(' %% pm_k_p_r %0.6f: sum(pm_weight_3d_riesz_qk_(1+tmp_index_))/(pm_weight_2d_k_p_r * 4*pi): %0.16f',pm_k_p_r,sum(pm_weight_3d_riesz_qk_(1+tmp_index_))/(4*pi*pm_weight_2d_k_p_r))); end;
end;%for pm_nk_p_r=0:pm_n_k_p_r-1;
%%%%%%%%;
pm_n_qkabc = pm_n_qk + n_M_imp*3;
pm_weight_3d_riesz_weight_imagecount_qkabc_ = cat(1,pm_weight_3d_riesz_qk_/scaling_volumetric,repmat(weight_imagecount_M_imp_,[3,1]));
pm_numerator_root_weight_3d_riesz_weight_imagecount_qkabc_ = reshape(sqrt(pm_weight_3d_riesz_weight_imagecount_qkabc_),[pm_n_qkabc,1]);
pm_denomator_root_weight_3d_riesz_weight_imagecount_qkabc_ = 1./max(1e-12,reshape(sqrt(pm_weight_3d_riesz_weight_imagecount_qkabc_),[pm_n_qkabc,1]));
%%%%%%%%;

%disp('returning before visualizing volume'); return;

if (flag_disp>0);
%%%%%%%%;
% Visualize a_k_p_qk_. ;
%%%%%%%%;
pm_a_k_p_qk__ = local_qk__from_qk_(pm_n_q,pm_n_k_p_r,pm_a_k_p_qk_);
a_k_p_reco_qk_ = zeros(n_qk,1);
pm_n_UX_rank = pm_n_k_p_r;
for pm_nUX_rank=0:pm_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
if (flag_verbose>1); disp(sprintf(' %% adding pm_nUX_rank %d/%d nk_p_r %d/%d',pm_nUX_rank,pm_n_UX_rank,nk_p_r,n_k_p_r)); end;
tmp_index_ = n_qk_csum_(1+nk_p_r) + (0:n_q-1);
a_k_p_reco_qk_(1+tmp_index_) = a_k_p_reco_qk_(1+tmp_index_) + from_pm_UX_kn__(1+nk_p_r,1+pm_nUX_rank)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_a_k_p_qk__(1:n_q,1+pm_nUX_rank);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,a_k_p_reco_qk_.*(2*pi)^3.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
fname_fig_pre = sprintf('%s/%s_a_k_p_qk_FIGA',str_dir_jpg_sub,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[95:0.5:99.5];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),a_x_u_reco_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[05:10:95];
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%disp('returning before loading lanczos iterates'); return;

%%%%%%%%;
n_iteration = numel(alph_tilde_i_); pm_T_tilde__ = real(spdiags([circshift(beta_tilde_i_,-1),alph_tilde_i_,beta_tilde_i_],[-1,0,+1],n_iteration,n_iteration));
lambda_xi__ = -Inf*ones(n_iteration,n_iteration); %<-- (index_lambda,niteration);
for niteration=0:n_iteration-1;
pm_T_tilde_sub__ = pm_T_tilde__(1:1+niteration,1:1+niteration);
lambda_sub_ = eigs(pm_T_tilde_sub__,[],1+niteration);
lambda_xi__(1:1+niteration,1+niteration) = sort(lambda_sub_,'ascend');
end;%for niteration=0:n_iteration-1;
S_x_ = sort(eigs(pm_T_tilde__,[],n_iteration),'ascend');
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
pm_v_tilde_qkabc_ = pm_v_tilde_qkabci__(:,1+niteration);
[pm_v_tilde_dvol_qk_,pm_v_tilde_polar_a_M_use_,pm_v_tilde_azimu_b_M_use_,pm_v_tilde_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_tilde_qkabc_);
[tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c] = local_qkabc_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,[],n_M_imp,pm_v_tilde_qkabc_,pm_v_tilde_qkabc_);
str_vv = sprintf('tmp_vv %0.2f,tmp_vv_dvol %0.2f,tmp_vv_a %0.2f,tmp_vv_b %0.2f,tmp_vv_c %0.2f',tmp_vv,tmp_vv_dvol,tmp_vv_a,tmp_vv_b,tmp_vv_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_vv)); end;
vv_n4__(1+niteration,:) = [tmp_vv_dvol;tmp_vv_a;tmp_vv_b;tmp_vv_c];
pm_T_tilde_sub__ = pm_T_tilde__(1:1+niteration,1:1+niteration);
[pm_TV_tilde_sub__,lambda_sub__] = eigs(pm_T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
for index_lambda=0:1+niteration-1;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
lambda_ns__(1+niteration,1+index_lambda) = lambda_use;
pm_TV_tilde_eig_ = pm_TV_tilde_sub__(:,ij_use);
pm_v_tilde_eig_qkabc_ = pm_v_tilde_qkabci__(:,1:1+niteration)*pm_TV_tilde_eig_;
[pm_v_tilde_eig_dvol_qk_,pm_v_tilde_eig_polar_a_M_use_,pm_v_tilde_eig_azimu_b_M_use_,pm_v_tilde_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_tilde_eig_qkabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_qkabc_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,[],n_M_imp,pm_v_tilde_eig_qkabc_,pm_v_tilde_eig_qkabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ee)); end;
ee_ns4___(1+niteration,1+index_lambda,:) = [tmp_ee_dvol;tmp_ee_a;tmp_ee_b;tmp_ee_c];
pm_v_eig_qkabc_ = bsxfun(@times,pm_denomator_root_weight_3d_riesz_weight_imagecount_qkabc_,pm_v_tilde_eig_qkabc_);
[pm_v_eig_dvol_qk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_eig_qkabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_qkabc_imagecount_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,pm_weight_3d_riesz_qk_,n_M_imp,weight_imagecount_M_imp_,pm_v_eig_qkabc_,pm_v_eig_qkabc_);
str_ff = sprintf('tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ff)); end;
ff_ns4___(1+niteration,1+index_lambda,:) = [tmp_ff_dvol;tmp_ff_a;tmp_ff_b;tmp_ff_c];
end;%for index_lambda=0:1+niteration-1;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ee_ns4___ vs ff_ns4___: %0.16f',fnorm(ee_ns4___ - ff_ns4___)/max(1e-12,fnorm(ee_ns4___)))); end;

%disp('returning before visualizing lanczos iterates'); return;

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

%disp('returning before scanning for relevant eigenvectors'); return;

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

%disp('returning before stepping through scan'); return;

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
pm_T_tilde_sub__ = pm_T_tilde__(1:1+niteration,1:1+niteration);
[pm_TV_tilde_sub__,lambda_sub__] = eigs(pm_T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
if (abs(lambda_use-lambda_ns__(1+nn,1+ns))> 1e-12); disp(sprintf(' %% Warning, abs(lambda_use-lambda_ns__(1+nn,1+ns)): %0.16f',abs(lambda_use-lambda_ns__(1+nn,1+ns)))); end;
pm_TV_tilde_eig_ = pm_TV_tilde_sub__(:,ij_use);
pm_v_tilde_eig_qkabc_ = pm_v_tilde_qkabci__(:,1:1+niteration)*pm_TV_tilde_eig_;
[pm_v_tilde_eig_dvol_qk_,pm_v_tilde_eig_polar_a_M_use_,pm_v_tilde_eig_azimu_b_M_use_,pm_v_tilde_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_tilde_eig_qkabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_qkabc_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,[],n_M_imp,pm_v_tilde_eig_qkabc_,pm_v_tilde_eig_qkabc_);
str_ee = sprintf('lambda %0.6f: tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',lambda_use,tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ee)); end;
pm_v_eig_qkabc_ = bsxfun(@times,pm_denomator_root_weight_3d_riesz_weight_imagecount_qkabc_,pm_v_tilde_eig_qkabc_);
[pm_v_eig_dvol_qk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_eig_qkabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_qkabc_imagecount_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,pm_weight_3d_riesz_qk_,n_M_imp,weight_imagecount_M_imp_,pm_v_eig_qkabc_,pm_v_eig_qkabc_);
str_ff = sprintf('lambda %0.6f: tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',lambda_use,tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>1); disp(sprintf(' %% %s',str_ff)); end;
if (abs(tmp_ff_dvol-ee_ns4___(1+nn,1+ns,1+0))> 1e-9); disp(sprintf(' %% Warning, abs(tmp_ff_dvol-ee_ns4___(1+nn,1+ns,1+0)): %0.16f',abs(tmp_ff_dvol-ee_ns4___(1+nn,1+ns,1+0)))); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nscan=0:n_scan-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%disp('returning before analyzing specific eigenevectors'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now analyze specific eigenvectors. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
niteration = n_iteration-1;
pm_T_tilde_sub__ = pm_T_tilde__(1:1+niteration,1:1+niteration);
[pm_TV_tilde_sub__,lambda_sub__] = eigs(pm_T_tilde_sub__,[],1+niteration);
lambda_sub_ = diag(lambda_sub__);
[lambda_srt_,ij_srt_] = sort(lambda_sub_,'ascend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for index_lambda = 0:niteration; %<-- 0 is most negative, niteration is most positive. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
ij_use = ij_srt_(1+index_lambda);
lambda_use = lambda_srt_(1+index_lambda);
str_fname_nopath_sub_prefix = sprintf('%s_i%.3dl%.3d',str_fname_nopath_prefix,niteration,index_lambda);
fname_sub_pre = sprintf('%s/%s',str_dir_mat_sub,str_fname_nopath_sub_prefix);
[flag_sub_skip,fname_sub_mat] = open_fname_tmp(fname_sub_pre);
%%%%%%%%%%%%%%%%;
if flag_recalc | ~flag_sub_skip;
%%%%%%%%%%%%%%%%;
% use pm_T_tilde_sub__ to estimate minimum eigenvector. ;
%%%%;
pm_TV_tilde_eig_ = pm_TV_tilde_sub__(:,ij_use);
pm_v_tilde_eig_qkabc_ = pm_v_tilde_qkabci__(:,1:1+niteration)*pm_TV_tilde_eig_;
[pm_v_tilde_eig_dvol_qk_,pm_v_tilde_eig_polar_a_M_use_,pm_v_tilde_eig_azimu_b_M_use_,pm_v_tilde_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_tilde_eig_qkabc_);
[tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c] = local_qkabc_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,[],n_M_imp,pm_v_tilde_eig_qkabc_,pm_v_tilde_eig_qkabc_);
str_ee = sprintf('tmp_ee %0.2f,tmp_ee_dvol %0.2f,tmp_ee_a %0.2f,tmp_ee_b %0.2f,tmp_ee_c %0.2f',tmp_ee,tmp_ee_dvol,tmp_ee_a,tmp_ee_b,tmp_ee_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ee)); end;
pm_v_eig_qkabc_ = bsxfun(@times,pm_denomator_root_weight_3d_riesz_weight_imagecount_qkabc_,pm_v_tilde_eig_qkabc_);
[pm_v_eig_dvol_qk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_eig_qkabc_);
[tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c] = local_qkabc_imagecount_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,pm_weight_3d_riesz_qk_,n_M_imp,weight_imagecount_M_imp_,pm_v_eig_qkabc_,pm_v_eig_qkabc_);
str_ff = sprintf('tmp_ff %0.2f,tmp_ff_dvol %0.2f,tmp_ff_a %0.2f,tmp_ff_b %0.2f,tmp_ff_c %0.2f',tmp_ff,tmp_ff_dvol,tmp_ff_a,tmp_ff_b,tmp_ff_c);
if (flag_verbose>0); disp(sprintf(' %% %s',str_ff)); end;
%%%%;
% Now evaluate ddssnll and associated rayleigh-quotient. ;
%%%%;
parameter_ddssnll = parameter;
parameter_ddssnll.flag_verbose = flag_verbose;
parameter_ddssnll.flag_implicit_dtau = flag_implicit_dtau;
%%%%;
if ~exist('pm_wS_from_single_shell_sba__','var'); pm_wS_from_single_shell_sba__=[]; end;
if ~exist('pm_dwSda_from_single_shell_sba__','var'); pm_dwSda_from_single_shell_sba__=[]; end;
if ~exist('pm_dwSdb_from_single_shell_sba__','var'); pm_dwSdb_from_single_shell_sba__=[]; end;
if ~exist('pm_ddwSdaa_from_single_shell_sba__','var'); pm_ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('pm_ddwSdab_from_single_shell_sba__','var'); pm_ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('pm_ddwSdbb_from_single_shell_sba__','var'); pm_ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('pm_n_R','var'); pm_n_R = []; end;
if ~exist('pm_R_use_R___','var'); pm_R_use_R___ = []; end;
if ~exist('pm_a_R_k_p_Rqk__','var'); pm_a_R_k_p_Rqk__=[]; end;
if ~exist('pm_ba_from_single_shell_Rbaba___','var'); pm_ba_from_single_shell_Rbaba___=[]; end;
if ~exist('pm_wS_from_R_single_shell_Rsba___','var'); pm_wS_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_dwSda_from_R_single_shell_Rsba___','var'); pm_dwSda_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_dwSdb_from_R_single_shell_Rsba___','var'); pm_dwSdb_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_ddwSdaa_from_R_single_shell_Rsba___','var'); pm_ddwSdaa_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_ddwSdab_from_R_single_shell_Rsba___','var'); pm_ddwSdab_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_ddwSdbb_from_R_single_shell_Rsba___','var'); pm_ddwSdbb_from_R_single_shell_Rsba___=[]; end;
if ~exist('pm_KAPPA','var'); pm_KAPPA = []; end;
[ ...
 ~ ...
,tmp_Hvt_qkabc_ ...
,tmp_Hv_q3d_k_p_qk_ ...
,tmp_Ht_q2d_M3__ ...
,tmp_a_restore_C2M0_k_p_qk_ ...
,tmp_Hvv_q3d_k_p_qk_ ...
,tmp_Hvt_q3d_k_p_qk_ ...
,tmp_Htv_q2d_M3__ ...
,tmp_Htt_q2d_M3__ ...
,tmp_dvol_a_k_p_qk_ ...
,tmp_dtau_euler_polar_a_M_use_ ...
,tmp_dtau_euler_azimu_b_M_use_ ...
,tmp_dtau_euler_gamma_z_M_use_ ...
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
ddssnll_from_a_k_p_4( ...
 parameter_ddssnll ...
,pm_n_qk ...
,pm_n_qk_csum_ ...
,pm_k_p_r_qk_ ...
,pm_k_p_azimu_b_qk_ ...
,pm_k_p_polar_a_qk_ ...
,pm_weight_3d_k_p_qk_ ...
,pm_weight_shell_qk_ ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,pm_weight_3d_k_p_r_ ...
,pm_k_c_0_qk_ ...
,pm_k_c_1_qk_ ...
,pm_k_c_2_qk_ ...
,pm_n_polar_a_k_ ...
,pm_polar_a_ka__ ...
,pm_n_azimu_b_ka__ ...
,qref_k_eq_d ...
,pm_a_k_p_qk_ ...
,pm_v_eig_dvol_qk_ ...
,pm_n_w_ ...
,pm_weight_2d_k_p_r_ ...
,pm_weight_2d_wk_ ...
,n_viewing_S_use ...
,pm_S_use_k_p_wkS__ ...
,viewing_polar_a_S_use_ ...
,viewing_azimu_b_S_use_ ...
,viewing_weight_S_use_ ...
,n_viewing_polar_a_use ...
,viewing_polar_a_use_ ...
,n_viewing_azimu_b_use_ ...
,viewing_gamma_z_S_use_ ...
,n_M_use ...
,weight_imagecount_M_use_ ...
,pm_M_use_k_p_wkM__ ...
,pm_n_CTF ...
,pm_index_nCTF_from_nM_ ...
,pm_CTF_k_p_wkC__ ...
,[] ...
,[] ...
,[] ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,euler_gamma_z_M_use_ ...
,pm_v_eig_polar_a_M_use_ ...
,pm_v_eig_azimu_b_M_use_ ...
,pm_v_eig_gamma_z_M_use_ ...
,pm_KAPPA ...
,pm_wS_from_single_shell_sba__ ...
,pm_dwSda_from_single_shell_sba__ ...
,pm_dwSdb_from_single_shell_sba__ ...
,pm_ddwSdaa_from_single_shell_sba__ ...
,pm_ddwSdab_from_single_shell_sba__ ...
,pm_ddwSdbb_from_single_shell_sba__ ...
,pm_n_R ...
,pm_R_use_R___ ...
,pm_a_R_k_p_Rqk__ ...
,pm_ba_from_single_shell_Rbaba___ ...
,pm_wS_from_R_single_shell_Rsba___ ...
,pm_dwSda_from_R_single_shell_Rsba___ ...
,pm_dwSdb_from_R_single_shell_Rsba___ ...
,pm_ddwSdaa_from_R_single_shell_Rsba___ ...
,pm_ddwSdab_from_R_single_shell_Rsba___ ...
,pm_ddwSdbb_from_R_single_shell_Rsba___ ...
);
%%%%;
pm_w_nottilde_qkabc_ = tmp_Hvt_qkabc_(1:pm_n_qkabc); %<-- ignore alignment component if flag_implicit_dtau. ;
pm_w_tilde_qkabc_ = ...
  bsxfun( ...
	  @times ...
	  ,pm_numerator_root_weight_3d_riesz_weight_imagecount_qkabc_ ...
	  ,pm_w_nottilde_qkabc_ ...
	  );
%%%%;
ddssnll_mid_q2d = local_qkabc_imagecount_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,pm_weight_3d_riesz_qk_,n_M_imp,weight_imagecount_M_imp_,pm_w_nottilde_qkabc_,pm_v_eig_qkabc_);
ddssnll_tilde_mid_q2d = local_qkabc_f_bar_dot_g_(pm_n_q,pm_n_k_p_r,[],n_M_imp,pm_w_tilde_qkabc_,pm_v_tilde_eig_qkabc_);
if (flag_verbose>0); disp(sprintf(' %% ddssnll_mid_q2d: %0.16f',ddssnll_mid_q2d)); end;
if (flag_verbose>0); disp(sprintf(' %% ddssnll_tilde_mid_q2d: %0.16f',ddssnll_tilde_mid_q2d)); end;
if (flag_verbose>0); disp(sprintf(' %% ddssnll_mid_q2d vs tmp_ddssnll_mid_q2d: %0.16f',fnorm(ddssnll_mid_q2d - tmp_ddssnll_mid_q2d)/max(1e-12,fnorm(ddssnll_mid_q2d)))); end;
%%%%;
save(fname_sub_mat ...
     ,'pm_n_k_p_r','pm_k_p_r_','pm_k_p_r_max','pm_weight_3d_k_p_r_','pm_weight_3d_riesz_qk_' ...
     ,'n_M_imp','n_M_use','weight_imagecount_M_imp_','weight_imagecount_M_use_' ...
     ,'scaling_volumetric','pm_weight_3d_riesz_weight_imagecount_qkabc_' ...
     ,'niteration','pm_T_tilde_sub__','pm_TV_tilde_sub__','lambda_sub_','lambda_srt_','ij_srt_' ...
     ,'pm_v_tilde_eig_qkabc_','tmp_ee','tmp_ee_dvol','tmp_ee_a','tmp_ee_b','tmp_ee_c' ...
     ,'pm_v_eig_qkabc_','tmp_ff','tmp_ff_dvol','tmp_ff_a','tmp_ff_b','tmp_ff_c' ...
     ,'tmp_Hvt_qkabc_' ...
     ,'tmp_Hv_q3d_k_p_qk_' ...
     ,'tmp_Ht_q2d_M3__' ...
     ,'tmp_a_restore_C2M0_k_p_qk_' ...
     ,'tmp_Hvv_q3d_k_p_qk_' ...
     ,'tmp_Hvt_q3d_k_p_qk_' ...
     ,'tmp_Htv_q2d_M3__' ...
     ,'tmp_Htt_q2d_M3__' ...
     ,'euler_polar_a_M_use_' ...
     ,'euler_azimu_b_M_use_' ...
     ,'euler_gamma_z_M_use_' ...
     ,'tmp_dtau_euler_polar_a_M_use_' ...
     ,'tmp_dtau_euler_azimu_b_M_use_' ...
     ,'tmp_dtau_euler_gamma_z_M_use_' ...
     ,'tmp_n_dvt'... 
     ,'tmp_dvt_'... 
     ,'tmp_dvt'... 
     ,'tmp_ssnll_tmp_q2d_dvt_'... 
     ,'tmp_dssnll_mid_q2d'... 
     ,'tmp_dssnll_dif_q2d'... 
     ,'tmp_dssnll_lsq_q2d'... 
     ,'tmp_ddssnll_mid_q2d'... 
     ,'tmp_ddssnll_dif_q2d'... 
     ,'tmp_ddssnll_lsq_q2d'... 
     ,'ddssnll_mid_q2d','ddssnll_tilde_mid_q2d' ...
    );
close_fname_tmp(fname_sub_pre);
%%%%%%%%%%%%%%%%;
end;%if flag_recalc | ~flag_sub_skip;
%%%%%%%%%%%%%%%%;

%disp('returning before creating subfigures'); return;

fname_fig_sub_A_jpg = sprintf('%s/%s_FIGA.jpg',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_B_jpg = sprintf('%s/%s_FIGB.jpg',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_C_jpg = sprintf('%s/%s_FIGC.jpg',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_D_jpg = sprintf('%s/%s_FIGD.jpg',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_E_jpg = sprintf('%s/%s_FIGE.jpg',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
flag_fig_sub_exist = 1 ...
& exist(fname_fig_sub_A_jpg,'file') ...
& exist(fname_fig_sub_B_jpg,'file') ...
& exist(fname_fig_sub_C_jpg,'file') ...
& exist(fname_fig_sub_D_jpg,'file') ...
& exist(fname_fig_sub_E_jpg,'file') ...
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if  exist(fname_sub_mat,'file');
if flag_replot | ~flag_fig_sub_exist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %s found, not creating',fname_sub_mat)); end;
load(fname_sub_mat);
%%%%%%%%;
str_sgtitle = sprintf('%s niteration %d lambda %+0.2f = exp(%+0.2f) index_lambda %d/%d: \n %s ',str_fname_nopath_prefix,niteration,lambda_use,log(abs(lambda_use)),index_lambda,niteration,str_ee);
%%%%%%%%;
% visualize pm_v_eig_dvol_qk_;
%%%%%%%%;
[pm_v_eig_dvol_qk_,pm_v_eig_polar_a_M_use_,pm_v_eig_azimu_b_M_use_,pm_v_eig_gamma_z_M_use_] = local_qk_a_b_c_from_qkabc_(pm_n_q,pm_n_k_p_r,n_M_imp,pm_v_eig_qkabc_);
pm_v_k_p_qk_ = pm_v_eig_dvol_qk_;
pm_v_k_p_qk__ = local_qk__from_qk_(pm_n_q,pm_n_k_p_r,pm_v_k_p_qk_);
v_k_p_reco_qk_ = zeros(n_qk,1);
pm_n_UX_rank = pm_n_k_p_r;
for nUX_rank=0:pm_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_qk_csum_(1+nk_p_r) + (0:n_q-1);
v_k_p_reco_qk_(1+tmp_index_) = v_k_p_reco_qk_(1+tmp_index_) + from_pm_UX_kn__(1+nk_p_r,1+nUX_rank)/max(1e-12,X_weight_r_(1+nk_p_r))*pm_v_k_p_qk__(1:n_q,1+nUX_rank);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
v_x_u_reco_ = xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,v_k_p_reco_qk_.*(2*pi)^3.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: v_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGA',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
for percent_threshold=[95:0.5:99.5];
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1(struct('percent_threshold_',[percent_threshold]),v_x_u_reco_);
title(sprintf('p %.2f',percent_threshold));
end;%for percent_threshold=[95:0.5:99.5];
sgtitle(sprintf('pm_v_eig_dvol_qk_: %s niteration %d lambda %+0.2f = exp(%+0.2f) index_lambda %d/%d: %s ',str_fname_nopath_prefix,niteration,lambda_use,log(abs(lambda_use)),index_lambda,niteration,str_ee),'Interpreter','none');
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
v_x_u_reco___ = reshape(v_x_u_reco_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
subplot(1,3,1);
imagesc(squeeze(real(v_x_u_reco___(:,:,end/2))));
axis image; axisnotick; colorbar;
subplot(1,3,2);
imagesc(squeeze(real(v_x_u_reco___(:,end/2,:))));
axis image; axisnotick; colorbar;
subplot(1,3,3);
imagesc(squeeze(real(v_x_u_reco___(end/2,:,:))));
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
vmax = max(abs(v_x_u_reco_),[],'all');
amax = max(abs(a_x_u_reco_),[],'all');
dvol_ = 1e-1*[-2:+2]; p_col = numel(dvol_);
prct_ = [97.5,98.0,98.5]; p_row = numel(prct_);
np=0;
for prow=0:p_row-1;
for pcol=0:p_col-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
prct = prct_(1+prow);
dvol = dvol_(1+pcol);
tmp_a_x_u_reco_ = a_x_u_reco_/max(1e-12,amax) + dvol*v_x_u_reco_/max(1e-12,vmax);
isosurface_f_x_u_1(struct('percent_threshold_',[prct]),tmp_a_x_u_reco_);
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
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_dtau_euler_polar_a_M_use_,tmp_dtau_euler_azimu_b_M_use_,tmp_dtau_euler_gamma_z_M_use_].^2,weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
figbig;
subplot(1,1,1);
hold on;
sphere_compass__0( ...
 struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',1) ...
,n_M_use ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
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
%%%%%%%%;
if flag_disp;
fname_fig_sub_pre = sprintf('%s/%s_FIGE',str_dir_jpg_sub,str_fname_nopath_sub_prefix);
fname_fig_sub_jpg = sprintf('%s.jpg',fname_fig_sub_pre);
if flag_replot | ~exist(fname_fig_sub_jpg,'file');
figure(1+nf);nf=nf+1;clf;figmed;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_dtau_euler_polar_a_M_use_,tmp_dtau_euler_azimu_b_M_use_,tmp_dtau_euler_gamma_z_M_use_].^2,weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
figbig;
for np=0:1;
subplot(1,2,1+np);
hold on;
plot_sphere_grid_0(struct('flag_solid',1));
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/4,'flag_2d_vs_3d',0) ...
,n_M_use ...
,euler_polar_a_M_use_ ...
,euler_azimu_b_M_use_ ...
,tmp_dtau_euler_nrm_polar_a_M_ ...
,tmp_dtau_euler_nrm_azimu_b_M_ ...
,tmp_dtau_euler_nrm_gamma_z_M_ ...
);
hold off;
xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d; axisnotick3d;
if np==0; view(0,+90); title('dtau_euler_ view(0,+90)','Interpreter','none'); end;
if np==1; view(0,-90); title('dtau_euler_ view(0,-90)','Interpreter','none'); end;
end;%for np=0:1;
%%%%;
sgtitle(str_sgtitle,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_sub_jpg));
print('-djpeg',fname_fig_sub_jpg);
end;%if flag_replot | ~exist(fname_fig_sub_jpg,'file');
close(gcf);
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_replot | ~flag_fig_sub_exist;
end;%if  exist(fname_sub_mat,'file');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for index_lambda = 0:niteration; %<-- 0 is most negative, niteration is most positive. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Finished analyzing specific eigenvectors. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;


if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

