function ...
ddssnll_0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,a_k_Y_quad_yk__ ...
,aa_k_Y_quad_yk__ ...
,dvol_a_k_Y_quad_yk_ ...
,dvol_a_k_Y_quad_yk__ ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,weight_3d_k_p_r_ ...
,a_k_p_quad_ ...
,dvol_a_k_p_quad_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,dvol_S_k_p_q2d_wkS__ ...
,dtau_S_k_p_q2d_wkS3___ ...
,dtau_dvol_S_k_p_q2d_wkS3___ ...
,dtau_dtau_S_k_p_q2d_wkS33____ ...
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
,index_neta_from_nM_ ...
,n_eta ...
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

str_thisfunction = 'ddssnll_0';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
test_slice_vs_volume_integral_5;
%%%%%%%%;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
end;%if (nargin<1);
%%%%%%%%;

na=0;

if (nargin<1+na); parameter=[]; end; na=na+1;
if (nargin<1+na); n_k_p_r=[]; end; na=na+1;
if (nargin<1+na); k_p_r_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_max=[]; end; na=na+1;
if (nargin<1+na); l_max_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_yk_=[]; end; na=na+1;
if (nargin<1+na); a_k_Y_quad_yk__=[]; end; na=na+1;
if (nargin<1+na); aa_k_Y_quad_yk__=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_Y_quad_yk_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_Y_quad_yk__=[]; end; na=na+1;
if (nargin<1+na); n_k_all=[]; end; na=na+1;
if (nargin<1+na); n_k_all_csum_=[]; end; na=na+1;
if (nargin<1+na); k_p_r_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_azimu_b_all_=[]; end; na=na+1;
if (nargin<1+na); k_p_polar_a_all_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_all_=[]; end; na=na+1;
if (nargin<1+na); weight_shell_k_=[]; end; na=na+1;
if (nargin<1+na); weight_3d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); a_k_p_quad_=[]; end; na=na+1;
if (nargin<1+na); dvol_a_k_p_quad_=[]; end; na=na+1;
if (nargin<1+na); n_w_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_k_p_r_=[]; end; na=na+1;
if (nargin<1+na); weight_2d_wk_=[]; end; na=na+1;
if (nargin<1+na); n_S=[]; end; na=na+1;
if (nargin<1+na); S_k_p_q2d_wkS__=[]; end; na=na+1;
if (nargin<1+na); dvol_S_k_p_q2d_wkS__=[]; end; na=na+1;
if (nargin<1+na); dtau_S_k_p_q2d_wkS3___=[]; end; na=na+1;
if (nargin<1+na); dtau_dvol_S_k_p_q2d_wkS3___=[]; end; na=na+1;
if (nargin<1+na); dtau_dtau_S_k_p_q2d_wkS33____=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_azimu_b_S_=[]; end; na=na+1;
if (nargin<1+na); viewing_weight_S_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_polar_a=[]; end; na=na+1;
if (nargin<1+na); viewing_polar_a_=[]; end; na=na+1;
if (nargin<1+na); n_viewing_azimu_b_=[]; end; na=na+1;
if (nargin<1+na); n_M=[]; end; na=na+1;
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
if (nargin<1+na); dtau_euler_polar_a_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_azimu_b_M_=[]; end; na=na+1;
if (nargin<1+na); dtau_euler_gamma_z_M_=[]; end; na=na+1;
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

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
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
if (~isfield(parameter,'kernel_qpro_l_max_use')); parameter.kernel_qpro_l_max_use = l_max; end; %<-- parameter_bookmark. ;
kernel_qpro_l_max_use = parameter.kernel_qpro_l_max_use;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
%%%%%%%%;

%%%%%%%%;
if isempty(a_k_Y_quad_yk__);
a_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_quad_yk__(1:n_lm,1+nk_p_r) = a_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
end;%if isempty(a_k_Y_quad_yk__);
%%%%%%%%;
if isempty(dvol_a_k_Y_quad_yk__);
dvol_a_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
dvol_a_k_Y_quad_yk__(1:n_lm,1+nk_p_r) = dvol_a_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
end;%if isempty(dvol_a_k_Y_quad_yk__);
%%%%%%%%;
if isempty(a_k_p_quad_);
tmp_t = tic;
[ ...
 a_k_p_quad_ ...
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
,a_k_Y_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t));
end;%if isempty(a_k_p_quad_);
%%%%%%%%;
if isempty(dvol_a_k_p_quad_);
tmp_t = tic;
[ ...
 dvol_a_k_p_quad_ ...
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
,dvol_a_k_Y_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% dvol_a_k_p_quad_ time %0.2fs',tmp_t));
end;%if isempty(dvol_a_k_p_quad_);
%%%%%%%%;

%%%%%%%%;
% Calculate derivative using ssnll_from_a_k_Y_12. ;
%%%%%%%%;
if ~isempty(dvol_S_k_p_q2d_wkS__); if (flag_verbose>0); disp(sprintf(' %% Warning, ~isempty(dvol_S_k_p_q2d_wkS__) in %s',str_thisfunction)); end; end;
if ~isempty(dtau_S_k_p_q2d_wkS3___); if (flag_verbose>0); disp(sprintf(' %% Warning, ~isempty(dtau_S_k_p_q2d_wkS3___) in %s',str_thisfunction)); end; end;
if ~isempty(dtau_dvol_S_k_p_q2d_wkS3___); if (flag_verbose>0); disp(sprintf(' %% Warning, ~isempty(dtau_dvol_S_k_p_q2d_wkS3___) in %s',str_thisfunction)); end; end;
if ~isempty(dtau_dtau_S_k_p_q2d_wkS33____); if (flag_verbose>0); disp(sprintf(' %% Warning, ~isempty(dtau_dtau_S_k_p_q2d_wkS33____) in %s',str_thisfunction)); end; end;
%%%%;
tmp_t = tic();
parameter_ssnll = struct('type','parameter');
parameter_ssnll.flag_verbose = 0*flag_verbose;
parameter_ssnll.viewing_k_eq_d = viewing_k_eq_d;
parameter_ssnll.template_k_eq_d = template_k_eq_d;
parameter_ssnll.n_order = n_order;
[ ...
 parameter_ssnll ...
,ssnll_q2d_M_ ...
,ssnll_q2d ...
,S_k_p_q2d_wkS__ ...
,dvol_ssnll_q2d_M_ ...
,dvol_ssnll_q2d ...
,dvol_S_k_p_q2d_wkS__ ...
,dvol_dvol_ssnll_q2d ...
,dtau_ssnll_q2d_M3__ ...
,dtau_ssnll_q2d ...
,dtau_S_k_p_q2d_wkS3___ ...
,dtau_dvol_ssnll_q2d_M3__ ...
,dtau_dvol_ssnll_q2d ...
,dtau_dvol_S_k_p_q2d_wkS3___ ...
,dtau_dtau_ssnll_q2d_M33___ ...
,dtau_dtau_ssnll_q2d ...
,dtau_dtau_S_k_p_q2d_wkS33____ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
] = ...
ssnll_from_a_k_Y_12( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk__ ...
,dvol_a_k_Y_quad_yk__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,dvol_S_k_p_q2d_wkS__ ...
,dtau_S_k_p_q2d_wkS3___ ...
,dtau_dvol_S_k_p_q2d_wkS3___ ...
,dtau_dtau_S_k_p_q2d_wkS33____ ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nM_ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% ssnll_from_a_k_Y_12: time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% We construct the riesz integration-weights on the sphere. ;
% These are associated with the riesz-potential 1/k^2.5, ;
% or a weighting-function (for the squared-L2-norm) of 1/k. ;
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

%%%%%%%%;
% Calibrate scaling factor. ;
%%%%%%%%;
term_deltafunc = sqrt(2*pi);
term_2 = (pi*k_p_r_max^2)/(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_wk_) vs (pi*k_p_r_max^2)/(4*pi^2): %0.16f',fnorm(sum(weight_2d_wk_) - term_2))); end;
term_3 = (4/3)*pi*k_p_r_max^3;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_all_) vs (4/3)*pi*k_p_r_max^3: %0.16f',fnorm(sum(weight_3d_k_all_) - term_3))); end;
term_3r = (4*pi^2*k_p_r_max^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_riesz__all_) vs 4*pi^2*k_p_r_max^2: %0.16f',fnorm(sum(weight_3d_riesz_k_all_) - term_3r))); end;
scaling_volumetric = term_3r / term_2 / term_deltafunc ;
%%%%%%%%;

%%%%%%%%;
% Calculate volumetric terms. ;
%%%%%%%%;
tmp_t = tic();
parameter_KAPPA = struct('type','KAPPA');
parameter_KAPPA.flag_verbose = 0*flag_verbose;
parameter_KAPPA.flag_kernel_qpro_d0 = flag_kernel_qpro_d0;
parameter_KAPPA.flag_kernel_qpro_d1 = flag_kernel_qpro_d1;
parameter_KAPPA.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
parameter_KAPPA.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
parameter_KAPPA.kernel_qpro_l_max_use = kernel_qpro_l_max_use;
[ ...
 parameter_KAPPA ...
,KAPPA ...
,a_restore_C2M0_k_Y_lmk__ ...
,a_restore_C1M1_k_Y_lmk__ ...
,a_restore_C0M2_k_Y_lmk__ ...
,dtau_a_restore_C2M0_k_Y_lmk__ ...
,dtau_a_restore_C1M1_k_Y_lmk__ ...
,dtau_a_restore_C0M2_k_Y_lmk__ ...
,dtau_dtau_a_restore_C2M0_k_Y_lmk__ ...
,dtau_dtau_a_restore_C1M1_k_Y_lmk__ ...
,dtau_dtau_a_restore_C0M2_k_Y_lmk__ ...
] = ...
kappa_qpro_apply_2( ...
 parameter_KAPPA ...
,KAPPA ...
,n_w_max ...
,n_M ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_qpro_apply_2 (yes derivatives) time %0.2fs',tmp_t)); end;
%%%%%%%%;

if isempty(aa_k_Y_quad_yk__);
%%%%%%%%;
% determine expansion for abs(a).^2. ;
%%%%%%%%;
[ ...
 aa_k_Y_quad_yk_ ...
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
,abs(a_k_p_quad_).^2 ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%;
aa_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
aa_k_Y_quad_yk__(1:n_lm,1+nk_p_r) = aa_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
end;%if isempty(aa_k_Y_quad_yk__);

%%%%%%%%;
% Convert volumetric terms. ;
%%%%%%%%;
for n_type = 0:3*3-1;
ntype = 0;
if n_type==ntype; tmp_yk__ = a_restore_C2M0_k_Y_lmk__; tmp_str = 'a_restore_C2M0_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = a_restore_C1M1_k_Y_lmk__; tmp_str = 'a_restore_C1M1_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = a_restore_C0M2_k_Y_lmk__; tmp_str = 'a_restore_C0M2_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_a_restore_C2M0_k_Y_lmk__; tmp_str = 'dtau_a_restore_C2M0_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_a_restore_C1M1_k_Y_lmk__; tmp_str = 'dtau_a_restore_C1M1_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_a_restore_C0M2_k_Y_lmk__; tmp_str = 'dtau_a_restore_C0M2_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_dtau_a_restore_C2M0_k_Y_lmk__; tmp_str = 'dtau_dtau_a_restore_C2M0_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_dtau_a_restore_C1M1_k_Y_lmk__; tmp_str = 'dtau_dtau_a_restore_C1M1_k_Y_lmk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = dtau_dtau_a_restore_C0M2_k_Y_lmk__; tmp_str = 'dtau_dtau_a_restore_C0M2_k_Y_lmk__'; end; ntype = ntype + 1;
%%%%;
tmp_yk_ = zeros(n_lm_sum,1);
tmp_t = tic();
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk_(1+tmp_index_) = tmp_yk__(1:n_lm,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
[ ...
 tmp_quad_ ...
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
tmp_t = toc(tmp_t); disp(sprintf(' %% %s: convert_spharm_to_k_p_4: time %0.2fs',tmp_str,tmp_t));
ntype = 0;
if n_type==ntype; a_restore_C2M0_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; a_restore_C1M1_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; a_restore_C0M2_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_a_restore_C2M0_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_a_restore_C1M1_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_a_restore_C0M2_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_dtau_a_restore_C2M0_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_dtau_a_restore_C1M1_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; dtau_dtau_a_restore_C0M2_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
end;%for n_type = 0:3*3-1;
%%%%%%%%;

%%%%%%%%;
% Collect terms into ssnll. ;
%%%%%%%%;
ssnll_q3d = ...
 + 0.5 * sum( (conj(aa_k_Y_quad_yk__) .* a_restore_C2M0_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* a_restore_C1M1_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( a_restore_C0M2_k_Y_lmk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
ssnll_q3d = ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% ssnll_q2d vs ssnll_q3d: %0.16f',fnorm(ssnll_q2d-ssnll_q3d)/fnorm(ssnll_q2d))); end;
%%%%%%%%;
ssnll_q3d = ...
 + 0.5 * sum( abs(a_k_p_quad_).^2 .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * 2*real(sum( conj(a_k_p_quad_).^1 .* a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ )) ...
 + 0.5 * sum( a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
;
ssnll_q3d = ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% ssnll_q2d vs ssnll_q3d: %0.16f',fnorm(ssnll_q2d-ssnll_q3d)/fnorm(ssnll_q2d))); end;
%%%%%%%%;
ssnll_q3d = ...
 + 0.5 * sum( conj(a_k_p_quad_) .* (a_k_p_quad_) .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * sum( conj(a_k_p_quad_) .* a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * sum( conj(a_restore_C1M1_k_p_quad_) .* (a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
 + 0.5 * sum( a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
;
ssnll_q3d = ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% ssnll_q2d vs ssnll_q3d: %0.16f',fnorm(ssnll_q2d-ssnll_q3d)/fnorm(ssnll_q2d))); end;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% comparing first-derivative (reasonably accurate): ')); end;
%%%%%%%%;
dtau_ssnll_q3d = ...
 + 0.5 * sum( (conj(aa_k_Y_quad_yk__) .* dtau_a_restore_C2M0_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* dtau_a_restore_C1M1_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( dtau_a_restore_C0M2_k_Y_lmk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
dtau_ssnll_q3d = dtau_ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% dtau_ssnll_q2d vs dtau_ssnll_q3d: %0.16f',fnorm(dtau_ssnll_q2d - dtau_ssnll_q3d)/max(1e-12,fnorm(dtau_ssnll_q2d)))); end;
%%%%%%%%;
dtau_ssnll_q3d = ...
 + 0.5 * sum( abs(a_k_p_quad_).^2 .* dtau_a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * 2*real(sum( conj(a_k_p_quad_).^1 .* dtau_a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ )) ...
 + 0.5 * sum( dtau_a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
  ;
dtau_ssnll_q3d = dtau_ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% dtau_ssnll_q2d vs dtau_ssnll_q3d: %0.16f',fnorm(dtau_ssnll_q2d - dtau_ssnll_q3d)/max(1e-12,fnorm(dtau_ssnll_q2d)))); end;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% comparing second-derivative (not as accurate): ')); end;
%%%%%%%%;
dtau_dtau_ssnll_q3d = ...
 + 0.5 * sum( (conj(aa_k_Y_quad_yk__) .* dtau_dtau_a_restore_C2M0_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* dtau_dtau_a_restore_C1M1_k_Y_lmk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( dtau_dtau_a_restore_C0M2_k_Y_lmk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
dtau_dtau_ssnll_q3d = dtau_dtau_ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_ssnll_q2d vs dtau_dtau_ssnll_q3d: %0.16f',fnorm(dtau_dtau_ssnll_q2d - dtau_dtau_ssnll_q3d)/max(1e-12,fnorm(dtau_dtau_ssnll_q2d)))); end;
%%%%%%%%;
dtau_dtau_ssnll_q3d = ...
 + 0.5 * sum( abs(a_k_p_quad_).^2 .* dtau_dtau_a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
 - 0.5 * 2*real(sum( conj(a_k_p_quad_).^1 .* dtau_dtau_a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ )) ...
 + 0.5 * sum( dtau_dtau_a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
  ;
dtau_dtau_ssnll_q3d = dtau_dtau_ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% dtau_dtau_ssnll_q2d vs dtau_dtau_ssnll_q3d: %0.16f',fnorm(dtau_dtau_ssnll_q2d - dtau_dtau_ssnll_q3d)/max(1e-12,fnorm(dtau_dtau_ssnll_q2d)))); end;
%%%%%%%%;

%%%%%%%%;
% Now construct hessian. ;
%%%%%%%%;
%{
  % As indicated via the above calculation: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_k_p_quad_) .* (a_k_p_quad_) .* a_restore_C2M0_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    - 0.5 * sum( conj(a_k_p_quad_) .* a_restore_C1M1_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    - 0.5 * sum( conj(a_restore_C1M1_k_p_quad_) .* (a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
    + 0.5 * sum( a_restore_C0M2_k_p_quad_ .* weight_3d_riesz_k_all_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % rewriting this as: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_) .* (a_) .* C2M0_ .* w_ ) ...
    - 0.5 * sum( conj(a_) .* C1M1_ .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    + 0.5 * sum( C0M2_ .* w_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % we have: ;
  ssnll_q3d = ...
    + 0.5 * sum( conj(a_ + dvol_) .* (a_ + dvol_) .* (C2M0_ + dtau_C2M0_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(a_ + dvol_) .* (C1M1_ + dtau_C1M1_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_ + dtau_C1M1_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_ + dvol_) .* w_ ) ...
    + 0.5 * sum( (C0M2_ + dtau_C0M2_*dtau_ + 0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % where we assume: ;
  % 1. dvol_ is a perturbation of the form dvol_k_p_quad_, ;
  % 2. dtau_ is a perturbation of the form dtau_M3__, ;
  % 3. dtau_CXMX_ is a jacobian (i.e., of the form dtau_CXMX_M3__), and ;
  % 4. dtau_dtau_CXMX__ is a hessian (i.e., of the form dtau_dtau_CXMX_M33__). ;
  % Using this notation, we can expand and retain terms of order 2 and lower: ;
  ssnll_q3d = ...
    ...
    + 0.5 * sum( conj(a_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    ...
    - 0.5 * sum( conj(a_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (a_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    ...
    + 0.5 * sum( (C0M2_) .* w_ ) ...
    + 0.5 * sum( (dtau_C0M2_*dtau_) .* w_ ) ...
    + 0.5 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % Combining terms of like order: ;
  ssnll_q3d = ...
    ...
    + 0.5 * sum( conj(a_) .* (a_) .* (C2M0_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    + 0.5 * sum( (C0M2_) .* w_ ) ...
    ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (a_) .* w_ ) ...
    + 0.5 * sum( (dtau_C0M2_*dtau_) .* w_ ) ...
    ...
    + 0.5 * sum( conj(dvol_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    + 0.5 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ...
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % From this expansion we can extract the gradient as a linear-operator acting on [dvol_;dtau_]: ;
  Gradient of ssnll_q3d = ...
      + 1.0 * real(sum( [ conj(a_) .* (C2M0_) - conj(C1M1) ] .* (dvol_) .* w_ )) ...
      - 1.0 * real(sum( [ conj(a_).*dtau_C1M1_ - 0.5*conj(a_).*a_.*dtau_C2M0_ - 0.5*dtau_C0M2_ ] * (dtau_) .* w_ )) ...
    ;
  Gradient of ssnll_q3d = Gradient of ssnll_q3d / scaling_volumetric ;
  % Similarly, we can extract the hessian as a quadratic-kernel acting on [dvol_;dtau_]: ;
  Hessian of ssnll_q3d = ...
    [ctranspose(dvol_) , ctranspose(dtau_)] * [ H ] * [dvol_;dtau_] ...
    ;
  % where the matrix H is a block-matrix with components: ;
  %      [ Hvv | Hvt ]  ;
  %  H = [-----+-----]  ;
  %      [ Htv | Htt ]  ;
  % such that: ;
  Hvv = quadratic-kernel acting on dvol_ = + 1.0 * sum( conj(dvol_) .* (C2M0_) .* (dvol_) .* w_ ) / scaling_volumetric ;
  Htt = quadratic-kernel acting on dtau_ = ...
    + 1.0 * sum( conj(a_) .* (a_ ) .* (0.5*ctranspose(dtau_)*dtau_dtau_C2M0__*(dtau_)) .* w_ ) ...
    - 1.0 * sum( conj(a_) .* (0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* w_ ) ...
    - 1.0 * sum( conj(0.5*ctranspose(dtau_)*dtau_dtau_C1M1__*(dtau_)) .* (a_) .* w_ ) ...
    + 1.0 * sum( (0.5*ctranspose(dtau_)*dtau_dtau_C0M2__*(dtau_)) .* w_ ) ...
    ;
  Htt = Htt / scaling_volumetric ;
  % Note that, via ssnll_from_a_k_Y_12, Htt can be rewritten in the simpler form: ;
  Htt = sum( dtau_M3__(1+nM,1+ntau0) * dtau_dtau_ssnll_q2d_M33___(1+nM,1+ntau0,1+ntau1) * dtau_M3__(1+nM,1+ntau1) ). ;
  % The operator Hvt has a domain compatible with dtau_ and a range compatible with dvol_. ;
  Hvt = ...
    + 1.0 * sum( conj(dvol_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    - 1.0 * sum( conj(dvol_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    ;
  Hvt = Hvt / scaling_volumetric ;
  % Note that, given the output of kappa_qpro_apply_2, we can rewrite Hvt as:
  Hvt = ...
      + 1.0 * sum( conj(dvol_) .* (a_) .* (|dtau| * dtau_a_restore_C2M0_k_p_quad_) .* w_ ) ...
      - 1.0 * sum( conj(dvol_) .* (|dtau| * dtau_a_restore_C1M1_k_p_quad_) .* w_ ) ...
      ;
  Hvt = Hvt / scaling_volumetric ;
  % The operator Htv has a domain compatible with dvol_ and a range compatible with dtau_. ;
  Htv = ...
    + 1.0 * sum( conj(dtau_C2M0_*dtau_) .* conj(a_) .* (dvol_) .* w_ ) ...
    - 1.0 * sum( conj(dtau_C1M1_*dtau_) .* (dvol_) .* w_ ) ...
    ;
  Htv = Htv / scaling_volumetric ;
  % Note that, via ssnll_from_a_k_Y_12, Htv can be rewritten in the simpler form: ;
  Htv = sum( dtau_M3__(1+nM,1+ntau1) * dtau_dvol_ssnll_q2d_M3__(1+nM,1+ntau1) ). ;
 %}
%%%%%%%%;
dtau_M3__ = [ ...
,dtau_euler_polar_a_M_ ...
,dtau_euler_azimu_b_M_ ...
,dtau_euler_gamma_z_M_ ...
] ;
dtau_M33___ = bsxfun(@times,reshape(dtau_M3__,[n_M,3,1]),reshape(dtau_M3__,[n_M,1,3]));
Htt_q2d = sum( dtau_dtau_ssnll_q2d_M33___ .* dtau_M33___ , 'all' ) ;
Htt_q3d = ...
  + 1.0 * sum( conj(a_k_p_quad_) .* (a_k_p_quad_ ) .* (0.5*dtau_dtau_a_restore_C2M0_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(a_k_p_quad_) .* (0.5*dtau_dtau_a_restore_C1M1_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(0.5*dtau_dtau_a_restore_C1M1_k_p_quad_) .* (a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  + 1.0 * sum( (0.5*dtau_dtau_a_restore_C0M2_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  ;
Htt_q3d = Htt_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% Htt_q2d vs Htt_q3d: %0.16f',fnorm(Htt_q2d - Htt_q3d)/max(1e-12,fnorm(Htt_q2d)))); end;
%%%%%%%%;
Htv_q3d = ...
  + 1.0 * sum( conj(dtau_a_restore_C2M0_k_p_quad_) .* conj(a_k_p_quad_) .* (dvol_a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(dtau_a_restore_C1M1_k_p_quad_) .* (dvol_a_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  ;
Htv_q3d = Htv_q3d / scaling_volumetric ;
Htv_q2d = sum( dtau_M3__ .* dtau_dvol_ssnll_q2d_M3__ , 'all' );
if (flag_verbose>0); disp(sprintf(' %% Htv_q2d vs Htv_q3d: %0.16f',fnorm(Htv_q2d - Htv_q3d)/max(1e-12,fnorm(Htv_q2d)))); end;
%%%%%%%%;
Hvt_q3d = ...
  + 1.0 * sum( conj(dvol_a_k_p_quad_) .* (a_k_p_quad_) .* (dtau_a_restore_C2M0_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  - 1.0 * sum( conj(dvol_a_k_p_quad_) .* (dtau_a_restore_C1M1_k_p_quad_) .* weight_3d_riesz_k_all_ ) ...
  ;
Hvt_q3d = Hvt_q3d / scaling_volumetric ;
Hvt_q2d = Htv_q2d ; %<-- due to symmetry. ;
if (flag_verbose>0); disp(sprintf(' %% Hvt_q2d vs Hvt_q3d: %0.16f',fnorm(Hvt_q2d - Hvt_q3d)/max(1e-12,fnorm(Hvt_q2d)))); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
