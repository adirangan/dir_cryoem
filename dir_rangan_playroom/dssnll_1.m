function ...
[ ...
 parameter ...
,G_bar_ykabc_ ...
,dvol_ssnll_q3d_k_Y_quad_yk_ ...
,dtau_ssnll_q2d_M3__ ...
,ssnll_q2d ...
] = ...
dssnll_1( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk_ ...
,a_k_Y_quad_yk__ ...
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
);

str_thisfunction = 'dssnll_1';

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

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_check')); parameter.flag_check = 0; end; %<-- parameter_bookmark. ;
flag_check = parameter.flag_check;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
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
kernel_qpro_polar_a_pole_north = min(pi/2,parameter.kernel_qpro_polar_a_pole_north);
parameter.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
if (~isfield(parameter,'kernel_qpro_polar_a_pole_south')); parameter.kernel_qpro_polar_a_pole_south = 3.5*pi/24; end; %<-- parameter_bookmark. ;
kernel_qpro_polar_a_pole_south = min(pi/2,parameter.kernel_qpro_polar_a_pole_south);
parameter.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
if ~isfield(parameter,'kernel_qpro_qref_k_eq_d_double'); parameter.kernel_qpro_qref_k_eq_d_double=0.5; end;
kernel_qpro_qref_k_eq_d_double=parameter.kernel_qpro_qref_k_eq_d_double;
if (~isfield(parameter,'kernel_qpro_l_max_use')); parameter.kernel_qpro_l_max_use = l_max; end; %<-- parameter_bookmark. ;
kernel_qpro_l_max_use = parameter.kernel_qpro_l_max_use;
%%%%;
if ~isfield(parameter,'flag_kernel_full'); parameter.flag_kernel_full=0; end;
flag_kernel_full=parameter.flag_kernel_full;
if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
flag_kernel_full = 1;
parameter.flag_kernel_full = flag_kernel_full;
end;%if (kernel_qpro_polar_a_pole_north + kernel_qpro_polar_a_pole_south > pi-1e-12);
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
a_k_Y_quad_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_quad_yk_);
end;%if isempty(a_k_Y_quad_yk__);
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
,a_k_Y_quad_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t)); end;
end;%if isempty(a_k_p_quad_);
%%%%%%%%;

%%%%%%%%;
% Calculate derivative using ssnll_from_a_k_Y_12. ;
%%%%%%%%;
dtau_S_k_p_q2d_wkS3___=[];
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
,~ ...
,~ ...
,~ ...
,~ ...
,dtau_ssnll_q2d_M3__ ...
,dtau_ssnll_q2d ...
,dtau_S_k_p_q2d_wkS3___ ...
] = ...
ssnll_from_a_k_Y_12( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_quad_yk__ ...
,[] ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_q2d_wkS__ ...
,[] ...
,dtau_S_k_p_q2d_wkS3___ ...
,[] ...
,[] ...
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
,[] ...
,[] ...
,[] ...
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
if (flag_verbose>0); disp(sprintf(' %% scaling_volumetric: %+0.6f',scaling_volumetric)); end;
if (flag_verbose>0); disp(sprintf(' %% (4*pi)^2 * sqrt(pi/2): %+0.6f',(4*pi)^2 * sqrt(pi/2))); end;
%%%%%%%%;

%%%%%%%%;
% Calculate volumetric terms. ;
%%%%%%%%;
tmp_t = tic();
parameter_KAPPA = struct('type','KAPPA');
parameter_KAPPA.flag_verbose = 0*flag_verbose;
parameter_KAPPA.flag_kernel_full = flag_kernel_full;
parameter_KAPPA.flag_kernel_qpro_d0 = flag_kernel_qpro_d0;
parameter_KAPPA.flag_kernel_qpro_d1 = flag_kernel_qpro_d1;
parameter_KAPPA.kernel_qpro_polar_a_pole_north = kernel_qpro_polar_a_pole_north;
parameter_KAPPA.kernel_qpro_polar_a_pole_south = kernel_qpro_polar_a_pole_south;
parameter_KAPPA.kernel_qpro_qref_k_eq_d_double = kernel_qpro_qref_k_eq_d_double;
parameter_KAPPA.kernel_qpro_l_max_use = kernel_qpro_l_max_use;
[ ...
 parameter_KAPPA ...
,KAPPA ...
,a_restore_C2M0_k_Y_yk__ ...
,a_restore_C1M1_k_Y_yk__ ...
,a_restore_C0M2_k_Y_yk__ ...
] = ...
kappa_qpro_apply_2( ...
 parameter_KAPPA ...
,KAPPA ...
,n_w_max ...
,n_M ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,M_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% kappa_qpro_apply_2 (not derivatives) time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Convert volumetric terms. ;
%%%%%%%%;
n_type_ = [0,1];
if flag_check;n_type_ = 0:3*1-1; end;
for n_type = n_type_;
ntype = 0;
if n_type==ntype; tmp_yk__ = a_restore_C2M0_k_Y_yk__; tmp_str = 'a_restore_C2M0_k_Y_yk__'; end; ntype = ntype + 1; %<-- Needed for Hvv_q3d. ;
if n_type==ntype; tmp_yk__ = a_restore_C1M1_k_Y_yk__; tmp_str = 'a_restore_C1M1_k_Y_yk__'; end; ntype = ntype + 1;
if n_type==ntype; tmp_yk__ = a_restore_C0M2_k_Y_yk__; tmp_str = 'a_restore_C0M2_k_Y_yk__'; end; ntype = ntype + 1;
%%%%;
tmp_t = tic();
tmp_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,tmp_yk__);
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %s: convert_spharm_to_k_p_4: time %0.2fs',tmp_str,tmp_t)); end;
ntype = 0;
if n_type==ntype; a_restore_C2M0_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; a_restore_C1M1_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
if n_type==ntype; a_restore_C0M2_k_p_quad_ = tmp_quad_; end; ntype = ntype + 1;
end;%for n_type = n_type_;
%%%%%%%%;

%%%%%%%%;
% determine expansion for a_times_a_restore_C2M0_k_Y_yk__. ;
%%%%%%%%;
[ ...
 a_times_a_restore_C2M0_k_Y_yk_ ...
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
,a_k_p_quad_.*a_restore_C2M0_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
%%%%;
a_times_a_restore_C2M0_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_times_a_restore_C2M0_k_Y_yk_);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

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
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_index_) = weight_3d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% determine expansion for abs(a).^2. ;
%%%%%%%%;
[ ...
 a_times_a_k_Y_quad_yk_ ...
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
a_times_a_k_Y_quad_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_times_a_k_Y_quad_yk_);
%%%%%%%%;

%%%%%%%%;
% Collect terms into ssnll. ;
%%%%%%%%;
ssnll_q3d = ...
 + 0.5 * sum( (conj(a_times_a_k_Y_quad_yk__) .* a_restore_C2M0_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* a_restore_C1M1_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( a_restore_C0M2_k_Y_yk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
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
ssnll_q3d = ...
 + 0.5 * sum( (conj(a_k_Y_quad_yk__) .* a_times_a_restore_C2M0_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) ) ...
 - 0.5 * 2*real(sum( (conj(a_k_Y_quad_yk__).^1 .* a_restore_C1M1_k_Y_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) )) ...
 + 0.5 * sum( a_restore_C0M2_k_Y_yk__(1+0,:) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1])*sqrt(4*pi) ) ...
;
ssnll_q3d = ssnll_q3d / scaling_volumetric ;
if (flag_verbose>0); disp(sprintf(' %% ssnll_q2d vs ssnll_q3d: %0.16f',fnorm(ssnll_q2d-ssnll_q3d)/fnorm(ssnll_q2d))); end;
%%%%%%%%;

%%%%%%%%;
% Now construct gradient. ;
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
  % Using this notation, we can expand and retain terms of order 1 and lower: ;
  ssnll_q3d = ...
    ...
    + 0.5 * sum( conj(a_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(dvol_) .* (a_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (dvol_) .* (C2M0_) .* w_ ) ...
    + 0.5 * sum( conj(a_) .* (a_) .* (dtau_C2M0_*dtau_) .* w_ ) ...
    ...
    - 0.5 * sum( conj(a_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(dvol_) .* (C1M1_) .* w_ ) ...
    - 0.5 * sum( conj(a_) .* (dtau_C1M1_*dtau_) .* w_ ) ...
    ...
    - 0.5 * sum( conj(C1M1_) .* (a_) .* w_ ) ...
    - 0.5 * sum( conj(C1M1_) .* (dvol_) .* w_ ) ...
    - 0.5 * sum( conj(dtau_C1M1_*dtau_) .* (a_) .* w_ ) ...
    ...
    + 0.5 * sum( (C0M2_) .* w_ ) ...
    + 0.5 * sum( (dtau_C0M2_*dtau_) .* w_ ) ...
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
    ;
  ssnll_q3d = ssnll_q3d / scaling_volumetric ;
  % From this expansion we can extract the gradient as a linear-operator acting on [dvol_;dtau_]: ;
  Gradient of ssnll_q3d = ...
      + 1.0 * real(sum( [ conj(a_) .* (C2M0_) - conj(C1M1) ] .* (dvol_) .* w_ )) ...
      - 1.0 * real(sum( [ conj(a_).*dtau_C1M1_ - 0.5*conj(a_).*a_.*dtau_C2M0_ - 0.5*dtau_C0M2_ ] * (dtau_) .* w_ )) ...
    ;
  Gradient of ssnll_q3d = Gradient of ssnll_q3d / scaling_volumetric ;
 %}
%%%%%%%%;

%%%%%%%%;
dvol_ssnll_q3d_k_p_quad_ = + 1.0 * ( conj(a_k_p_quad_) .* (a_restore_C2M0_k_p_quad_) - conj(a_restore_C1M1_k_p_quad_) );
dvol_ssnll_q3d_k_Y_quad_yk__ = + 1.0 * ( conj(a_times_a_restore_C2M0_k_Y_yk__) - conj(a_restore_C1M1_k_Y_yk__) );
dvol_ssnll_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,dvol_ssnll_q3d_k_Y_quad_yk__);
%%%%%%%%
[ ...
 dvol_ssnll_q3d_k_p_reco_ ...
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
,dvol_ssnll_q3d_k_Y_quad_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_verbose>0); disp(sprintf(' %% dvol_ssnll_q3d_k_p_quad_ vs dvol_ssnll_q3d_k_p_reco_: %0.16f',fnorm(dvol_ssnll_q3d_k_p_quad_ - dvol_ssnll_q3d_k_p_reco_)/max(1e-12,fnorm(dvol_ssnll_q3d_k_p_quad_)))); end;
%%%%%%%%;
[ ...
 dvol_ssnll_q3d_k_Y_reco_yk_ ...
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
,dvol_ssnll_q3d_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
if (flag_verbose>0); disp(sprintf(' %% dvol_ssnll_q3d_k_Y_quad_yk_ vs dvol_ssnll_q3d_k_Y_reco_yk_: %0.16f',fnorm(dvol_ssnll_q3d_k_Y_quad_yk_ - dvol_ssnll_q3d_k_Y_reco_yk_)/max(1e-12,fnorm(dvol_ssnll_q3d_k_Y_quad_yk_)))); end;
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
p_row=1; p_col=2*1; np=0;
%%%%;
tmp_k_p_quad_ = dvol_ssnll_q3d_k_p_quad_;
tmp_k_Y_quad_yk__ = dvol_ssnll_q3d_k_Y_quad_yk__;
tmp_str = 'dvol_ssnll_q3d';
tmp_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,tmp_k_Y_quad_yk__);
l2_a0 = sum ((conj(tmp_k_Y_quad_yk__).*tmp_k_Y_quad_yk__)*reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) , 'all' )/scaling_volumetric;
disp(sprintf(' %% %s: tmp_k_Y_quad_yk__: l2_a0: %0.16f',tmp_str,l2_a0));
l2_a0 = sum ((conj(tmp_k_p_quad_).*tmp_k_p_quad_).*weight_3d_riesz_k_all_ , 'all' )/scaling_volumetric;
disp(sprintf(' %% %s: tmp_k_p_quad_    : l2_a0: %0.16f',tmp_str,l2_a0));
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,abs(tmp_k_Y_quad_yk_),'k.');
xlabel('Y_l_val_','Interpreter','none');
title(sprintf('abs(%s)',tmp_str),'Interpreter','none');
%%%%;
tmp_abs_k_p_quad_ = abs(tmp_k_p_quad_).*sqrt(weight_3d_riesz_k_all_);
flag_2d_vs_3d = 0; c_use__ = colormap_81s;
nk_p_r = floor(n_k_p_r/2);
tmp_index_ = n_k_all_csum_(1+nk_p_r+0):n_k_all_csum_(1+nk_p_r+1)-1;
lim_ = prctile(tmp_abs_k_p_quad_(1+tmp_index_),[ 5,95]);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_(1+tmp_index_) ...
,k_p_azimu_b_all_(1+tmp_index_) ...
,tmp_abs_k_p_quad_(1+tmp_index_) ... 
,lim_ ... 
,c_use__ ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
axis equal; axis vis3d; axisnotick3d;
title(sprintf('abs(%s)',tmp_str),'Interpreter','none');
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_check;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now we calculate the gradient: ;
%%%%%%%%;
dvol_ssnll_q3d_k_p_quad_ = + 1.0 * ( conj(a_k_p_quad_) .* (a_restore_C2M0_k_p_quad_) - conj(a_restore_C1M1_k_p_quad_) );
dvol_ssnll_q3d_k_Y_quad_yk__ = + 1.0 * ( conj(a_times_a_restore_C2M0_k_Y_yk__) - conj(a_restore_C1M1_k_Y_yk__) );
dtau_ssnll_q2d_M3__ = dtau_ssnll_q2d_M3__ ;
%%%%%%%%;
dvol_ssnll_q3d_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,dvol_ssnll_q3d_k_Y_quad_yk__);
%%%%%%%%;
G_ykabc_ = cat(1,dvol_ssnll_q3d_k_Y_quad_yk_,dtau_ssnll_q2d_M3__(:,1+0),dtau_ssnll_q2d_M3__(:,1+1),dtau_ssnll_q2d_M3__(:,1+2));
G_bar_ykabc_ = conj(G_ykabc_); %<-- note here we take a conjugate to allow for standard dot-product. ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
