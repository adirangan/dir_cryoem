function ...
[ ...
 parameter ...
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
] = ...
eig_ddssnll_lanczos_1( ...
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
,U_tilde_SmallRotation_Delta_ykabc3__ ...
,v_tilde_ykabci__  ...
,w_tilde_ykabc_  ...
,alph_tilde_i_ ...
,beta_tilde_i_ ... 
);

str_thisfunction = 'eig_ddssnll_lanczos_1';

%%%%%%%%;
if (nargin<1);
%%%%%%%%;
test_eig_ddssnll_lanczos_surrogate_1;
disp(sprintf(' %% returning')); return;
%%%%%%%%;
%test_slice_vs_volume_integral_5;
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
if (nargin<1+na); U_tilde_SmallRotation_Delta_ykabc3__=[]; end; na=na+1;
if (nargin<1+na); v_tilde_ykabci__=[]; end; na=na+1;
if (nargin<1+na); w_tilde_ykabc_=[]; end; na=na+1;
if (nargin<1+na); alph_tilde_i_=[]; end; na=na+1;
if (nargin<1+na); beta_tilde_i_=[]; end; na=na+1; 

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
if (~isfield(parameter,'flag_check')); parameter.flag_check = 0; end; %<-- parameter_bookmark. ;
flag_check = parameter.flag_check;
if (~isfield(parameter,'flag_disp')); parameter.flag_disp = 0; end; %<-- parameter_bookmark. ;
flag_disp = parameter.flag_disp; nf=0;
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
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% setting indices')); end;
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
if isempty(v_tilde_ykabci__); v_tilde_ykabci__ = zeros(n_lm_sum + n_M*3,0); end;
if isempty(w_tilde_ykabc_); w_tilde_ykabc_ = zeros(n_lm_sum + n_M*3,0); end;
if isempty(alph_tilde_i_); alph_tilde_i_ = zeros(0,1); end;
if isempty(beta_tilde_i_); beta_tilde_i_ = zeros(0,1); end;

%%%%%%%%;
if isempty(n_CTF);
n_CTF = 1;
CTF_k_p_r_kC__ = ones(n_k_p_r,n_CTF);
CTF_k_p_r_wkC__ = ones(n_w_sum,n_CTF);
index_nCTF_from_nM_ = zeros(n_M,1);
end;%if isempty(n_CTF);
%%%%;
if isempty(n_eta);
n_eta = 1;
eta_k_p_r_ke__ = ones(n_k_p_r,n_eta);
eta_k_p_r_wke__ = ones(n_w_sum,n_eta);
index_neta_from_nM_ = zeros(n_M,1);
end;%if isempty(n_eta);
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct weight_3d_riesz')); end;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_all_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_all_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set scaling_volumetric')); end;
%%%%%%%%;
term_deltafunc = sqrt(2*pi);
term_2 = (pi*k_p_r_max^2)/(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_wk_) vs (pi*k_p_r_max^2)/(4*pi^2): %0.16f',fnorm(sum(weight_2d_wk_) - term_2))); end;
term_3 = (4/3)*pi*k_p_r_max^3;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_all_) vs (4/3)*pi*k_p_r_max^3: %0.16f',fnorm(sum(weight_3d_k_all_) - term_3))); end;
term_3r = (4*pi^2*k_p_r_max^2);
scaling_volumetric = term_3r / term_2 / term_deltafunc ;
if (flag_verbose>0); disp(sprintf(' %% scaling_volumetric: %+0.6f',scaling_volumetric)); end;
if (flag_verbose>0); disp(sprintf(' %% (4*pi)^2 * sqrt(pi/2): %+0.6f',(4*pi)^2 * sqrt(pi/2))); end;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set weight_3d_riesz_yk_')); end;
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
%%%%%%%%;
n_ykabc = n_lm_sum + n_M*3; index_yk_ = [0:n_lm_sum-1]; index_abc_ = n_lm_sum + [0:n_M*3-1];
weight_3d_riesz_ykabc_ = cat(1,weight_3d_riesz_yk_/scaling_volumetric,ones(n_M*3,1));
numerator_root_weight_3d_riesz_ykabc_ = reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]);
denomator_root_weight_3d_riesz_ykabc_ = 1./max(1e-12,reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]));

%%%%%%%%;
a_k_Y_quad_yk__ = [];
if (size(a_k_Y_quad_yk_,1)==n_lm_max & size(a_k_Y_quad_yk_,2)==n_k_p_r); a_k_Y_quad_yk__ = a_k_Y_quad_yk_; a_k_Y_quad_yk_ = []; end;
if isempty(a_k_Y_quad_yk__); a_k_Y_quad_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_quad_yk_); end;
if isempty(a_k_Y_quad_yk_); a_k_Y_quad_yk_ = local_yk_from_yk__(n_k_p_r,l_max_,a_k_Y_quad_yk__); end;
%%%%%%%%;
if isempty(a_k_p_quad_);
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set a_k_p_quad_')); end;
tmp_t = tic;
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
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_: time %0.2fs',tmp_t));
end;%if isempty(a_k_p_quad_);
%%%%%%%%;

if isempty(U_tilde_SmallRotation_Delta_ykabc3__) & ~flag_ignore_U;
%%%%%%%%;
tmp_t = tic();
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set U_tilde_SmallRotation_Delta_ykabc3__')); end;
parameter_SmallRotation = struct('type','parameter');
parameter_SmallRotation.flag_verbose = flag_verbose;
parameter_SmallRotation.flag_check = 0;
[ ...
 ~ ...
,U_tilde_SmallRotation_Delta_ykabcs__ ...
,U_SmallRotation_Delta_ykabcs__ ...
,S_SmallRotation_Delta_s_ ...
,V_SmallRotation_Delta_ss__ ...
] = ...
U_SmallRotation_1( ...
 parameter_SmallRotation ...
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
);
U_SmallRotation_Delta_ykabc3__ = U_SmallRotation_Delta_ykabcs__(:,1:3);
tmp_t = toc(tmp_t); disp(sprintf(' %% U_SmallRotation_Delta_ykabc3__: time %0.2fs',tmp_t));
U_SmallRotation_Delta_ykabc3__ = local_normalize_fn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__);
[tmp_UU_33__] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,U_SmallRotation_Delta_ykabc3__);
if (flag_verbose>0); disp(sprintf(' %% abs(tmp_UU_33__): \n %s',num2str(transpose(abs(tmp_UU_33__(:))),' %+0.16f %+0.16f %+0.16f\n'))); end;
U_tilde_SmallRotation_Delta_ykabc3__ = U_tilde_SmallRotation_Delta_ykabcs__(:,1:3);
%%%%%%%%;
end;%if isempty(U_tilde_SmallRotation_Delta_ykabc3__) & ~flag_ignore_U;
%%%%%%%%;
if ~flag_ignore_U;
[tmp_UU_33__] = local_weightless_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__,U_tilde_SmallRotation_Delta_ykabc3__);
if (flag_verbose>0); disp(sprintf(' %% abs(tmp_UU_33__): \n %s',num2str(transpose(abs(tmp_UU_33__(:))),' %+0.16f %+0.16f %+0.16f\n'))); end;
end;%if ~flag_ignore_U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% perform lanczos iteration ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
lanczos_n_iteration_sum = lanczos_n_iteration_pre + lanczos_n_iteration_cur;
v_tilde_ykabci__ = cat(2,v_tilde_ykabci__,zeros(n_ykabc,lanczos_n_iteration_cur));
alph_tilde_i_ = cat(1,alph_tilde_i_,zeros(lanczos_n_iteration_cur,1));
beta_tilde_i_ = cat(1,beta_tilde_i_,zeros(lanczos_n_iteration_cur,1));
%%%%%%%%%%%%%%%%;
for niteration=lanczos_n_iteration_pre:lanczos_n_iteration_sum-1;
%%%%%%%%%%%%%%%%;
if niteration==0;
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, setting beta to zero ',niteration)); end;
beta_tilde_i_(1+niteration)=0;
end;%if niteration==0;
if niteration> 0;
beta_tilde_i_(1+niteration) = sqrt(local_weightless_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,w_tilde_ykabc_,w_tilde_ykabc_));
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, beta_tilde_i_(1+niteration) = %0.6f',niteration,beta_tilde_i_(1+niteration))); end;
end;%if niteration> 0;
%%%%%%%%;
if abs(beta_tilde_i_(1+niteration))>=1e-12;
v_tilde_ykabc_ = w_tilde_ykabc_/beta_tilde_i_(1+niteration);
end;%if abs(beta_tilde_i_(1+niteration))>=1e-12;
if abs(beta_tilde_i_(1+niteration))< 1e-12;
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, beta close to zero, generating new v_tilde_ykabc_',niteration)); end;
dvol_yk_ =  rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_quad_yk_,[2*pi,1*pi,2*pi]*rand());
dtau_euler_polar_a_M_ = 1*pi*rand(n_M,1);
dtau_euler_azimu_b_M_ = 2*pi*rand(n_M,1);
dtau_euler_gamma_z_M_ = 2*pi*rand(n_M,1);
v_tilde_ykabc_ = bsxfun(@times,numerator_root_weight_3d_riesz_ykabc_,local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,dvol_yk_,dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_));
if niteration>=1;
v_tilde_ykabc_ = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabci__(:,1+[0:niteration-1]),v_tilde_ykabc_);
v_tilde_ykabc_ = local_weightless_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabc_);
end;%if niteration>=1;
if ~flag_ignore_U; v_tilde_ykabc_ = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__,v_tilde_ykabc_); end;
if flag_ignore_tau; v_tilde_ykabc_(1+index_abc_)=0; end;
v_tilde_ykabc_ = local_weightless_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabc_);
tmp_vv = local_weightless_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabc_,v_tilde_ykabc_);
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, new v_tilde_ykabc_: tmp_vv: %0.16f',niteration,tmp_vv)); end;
if ~flag_ignore_U;
[tmp_Uv_3_] = local_weightless_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__,v_tilde_ykabc_);
if (flag_verbose>0); disp(sprintf(' %% abs(tmp_Uv_3_): %s',num2str(transpose(abs(tmp_Uv_3_)),' %+0.16f'))); end;
end;%if ~flag_ignore_U;
end;%if abs(beta_tilde_i_(1+niteration))< 1e-12;
%%%%%%%%;
v_tilde_ykabci__(:,1+niteration) = v_tilde_ykabc_;
%%%%;
% calculate w_tilde_ykabc_ = PHP * v_tilde_ykabci__(:,1+niteration) ;
%%%%;
w_tilde_ykabc_ = v_tilde_ykabci__(:,1+niteration);
if ~flag_ignore_U; w_tilde_ykabc_ = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__,w_tilde_ykabc_); end; %<-- Projection. ;
if flag_ignore_tau; w_tilde_ykabc_(1+index_abc_)=0; end;
eig_ddssnll_lanczos_helper_1; %<-- apply weighted version of H. ;
if flag_ignore_tau; w_tilde_ykabc_(1+index_abc_)=0; end;
if ~flag_ignore_U; w_tilde_ykabc_ = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_tilde_SmallRotation_Delta_ykabc3__,w_tilde_ykabc_); end; %<-- Projection. ;
%%%%;
alph_tilde_i_(1+niteration) = local_weightless_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,w_tilde_ykabc_,v_tilde_ykabci__(:,1+niteration));
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, alph_tilde_i_(1+niteration) = %0.6f',niteration,alph_tilde_i_(1+niteration))); end;
w_tilde_ykabc_ = w_tilde_ykabc_ - alph_tilde_i_(1+niteration)*v_tilde_ykabci__(:,1+niteration+0);
if niteration>=1; w_tilde_ykabc_ = w_tilde_ykabc_ - beta_tilde_i_(1+niteration)*v_tilde_ykabci__(:,1+niteration-1); end;
%%%%%%%%;
% This next step seems necessary for stability. ;
%%%%%%%%;
w_tilde_ykabc_ = local_weightless_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabci__(:,1:niteration),w_tilde_ykabc_); %<-- should be unnecessary. ;
%%%%%%%%;
if (flag_verbose>0);
[tmp_VV_ii__] = local_weightless_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_tilde_ykabci__(:,1:1+niteration),v_tilde_ykabci__(:,1:1+niteration));
tmp_str_format = ' %% %%'; for ni=0:niteration; tmp_str_format = sprintf('%s %%+0.4f',tmp_str_format); end; tmp_str_format = sprintf('%s\\n',tmp_str_format);
disp(sprintf(' %% abs(tmp_VV_ii__): \n %s',num2str(transpose(abs(tmp_VV_ii__(:))),tmp_str_format)));
end;%if (flag_verbose>0);
%%%%%%%%%%%%%%%%;
end;%for niteration=lanczos_n_iteration_pre:lanczos_n_iteration_sum-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
if (flag_disp>0);
%%%%%%%%;
n_iteration = numel(alph_tilde_i_); T__ = real(spdiags([circshift(beta_tilde_i_,-1),alph_tilde_i_,beta_tilde_i_],[-1,0,+1],n_iteration,n_iteration));
lambda_xi__ = -Inf*ones(n_iteration,n_iteration);
for niteration=0:n_iteration-1;
T_sub__ = T__(1:1+niteration,1:1+niteration);
lambda_sub_ = eigs(T_sub__,[],1+niteration);
lambda_xi__(1:1+niteration,1+niteration) = sort(lambda_sub_,'ascend');
end;%for niteration=0:n_iteration-1;
S_x_ = sort(eigs(T__,[],n_iteration),'ascend');
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fig81s;
markersize_use = 8;
linewidth_sml = 0.5;
linewidth_big = 2;
%%%%;
subplot(1,1,1);
hold on;
plot(repmat([0;n_iteration],[1,n_iteration]),repmat(reshape(S_x_,[1,n_iteration]),[2,1]),'-','Color',0.85*[1,1,1],'LineWidth',linewidth_sml);
ni_xi__ = repmat([1:n_iteration],[n_iteration,1]);
tmp_index_ = efind(isfinite(lambda_xi__));
plot(ni_xi__(1+tmp_index_),lambda_xi__(1+tmp_index_),'r.','MarkerSize',markersize_use);
hold off;
xlabel('iteration'); ylabel('sigma');
xlim([0,1+n_iteration]);
ylim([min(S_x_)-0.25,max(S_x_)+0.25]);
%%%%%%%%;
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;




