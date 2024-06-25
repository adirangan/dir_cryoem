function ...
[ ...
 parameter ...
,v_ykabci__  ...
,w_ykabc_  ...
,alph_i_ ...
,beta_i_ ... 
] = ...
eig_ddssnll_lanczos_0( ...
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
,index_neta_from_nM_ ...
,n_eta ...
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
,U_SmallRotation_Delta_ykabc3__ ...
,v_ykabci__  ...
,w_ykabc_  ...
,alph_i_ ...
,beta_i_ ... 
);

str_thisfunction = 'eig_ddssnll_lanczos_0';

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
if (nargin<1+na); U_SmallRotation_Delta_ykabc3__=[]; end; na=na+1;
if (nargin<1+na); v_ykabci__=[]; end; na=na+1;
if (nargin<1+na); w_ykabc_=[]; end; na=na+1;
if (nargin<1+na); alph_i_=[]; end; na=na+1;
if (nargin<1+na); beta_i_=[]; end; na=na+1; 

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
if (~isfield(parameter,'kernel_qpro_l_max_use')); parameter.kernel_qpro_l_max_use = l_max; end; %<-- parameter_bookmark. ;
kernel_qpro_l_max_use = parameter.kernel_qpro_l_max_use;
if (~isfield(parameter,'lanczos_n_iteration_cur')); parameter.lanczos_n_iteration_cur = 8; end; %<-- parameter_bookmark. ;
lanczos_n_iteration_cur = parameter.lanczos_n_iteration_cur;
lanczos_n_iteration_pre = size(v_ykabci__,2);
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
if isempty(v_ykabci__); v_ykabci__ = zeros(n_lm_sum + 3*n_M,0); end;
if isempty(w_ykabc_); w_ykabc_ = zeros(n_lm_sum + 3*n_M,0); end;
if isempty(alph_i_); alph_i_ = zeros(0,1); end;
if isempty(beta_i_); beta_i_ = zeros(0,1); end;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct weight_3d_riesz')); end;
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
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set scaling_volumetric')); end;
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

if isempty(U_SmallRotation_Delta_ykabc3__);
%%%%%%%%;
tmp_t = tic();
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% set U_SmallRotation_Delta_ykabc3__')); end;
parameter_SmallRotation = struct('type','parameter');
parameter_SmallRotation.flag_verbose = flag_verbose;
parameter_SmallRotation.flag_check = 0;
[ ...
 ~ ...
,U_SmallRotation_Delta_ykabcs__ ...
,S_SmallRotation_Delta_s_ ...
,V_SmallRotation_Delta_ss__ ...
] = ...
U_SmallRotation_0( ...
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
,index_neta_from_nM_ ...
,n_eta ...
,eta_k_p_r_ke__ ...
,eta_k_p_wke__ ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
);
U_SmallRotation_Delta_ykabc3__ = U_SmallRotation_Delta_ykabcs__(:,1:3);
tmp_t = toc(tmp_t); disp(sprintf(' %% U_SmallRotation_Delta_ykabc3__: time %0.2fs',tmp_t));
%%%%%%%%;
end;%if isempty(U_SmallRotation_Delta_ykabc3__);
%%%%%%%%;
U_SmallRotation_Delta_ykabc3__ = local_normalize_fn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__);
[tmp_UU_33__] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,U_SmallRotation_Delta_ykabc3__);
if (flag_verbose>0); disp(sprintf(' %% abs(tmp_UU_33__): \n %s',num2str(transpose(abs(tmp_UU_33__(:))),' %+0.16f %+0.16f %+0.16f\n'))); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% perform lanczos iteration ')); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_ykabc = n_lm_sum + 3*n_M;
lanczos_n_iteration_sum = lanczos_n_iteration_pre + lanczos_n_iteration_cur;
v_ykabci__ = cat(2,v_ykabci__,zeros(n_ykabc,lanczos_n_iteration_cur));
alph_i_ = cat(1,alph_i_,zeros(lanczos_n_iteration_cur,1));
beta_i_ = cat(1,beta_i_,zeros(lanczos_n_iteration_cur,1));
%%%%%%%%%%%%%%%%;
for niteration=lanczos_n_iteration_pre:lanczos_n_iteration_sum-1;
%%%%%%%%%%%%%%%%;
if niteration==0;
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, setting beta to zero ',niteration)); end;
beta_i_(1+niteration)=0;
end;%if niteration==0;
if niteration> 0;
beta_i_(1+niteration) = sqrt(local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,w_ykabc_,w_ykabc_));
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, beta_i_(1+niteration) = %0.6f',niteration,beta_i_(1+niteration))); end;
end;%if niteration> 0;
%%%%%%%%;
if abs(beta_i_(1+niteration))>=1e-12;
v_ykabc_ = w_ykabc_/beta_i_(1+niteration);
end;%if abs(beta_i_(1+niteration))>=1e-12;
if abs(beta_i_(1+niteration))< 1e-12;
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, beta close to zero, generating new v_ykabc_',niteration)); end;
dvol_yk_ =  rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_quad_yk_,[2*pi,1*pi,2*pi]*rand());
dtau_euler_polar_a_M_ = 1*pi*rand(n_M,1);
dtau_euler_azimu_b_M_ = 2*pi*rand(n_M,1);
dtau_euler_gamma_z_M_ = 2*pi*rand(n_M,1);
v_ykabc_ = cat(1,dvol_yk_,dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_);
if niteration>=1;
v_ykabc_ = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_ykabci__(:,1+[0:niteration-1]),v_ykabc_);
v_ykabc_ = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_ykabc_);
end;%if niteration>=1;
v_ykabc_ = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,v_ykabc_);
v_ykabc_ = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_ykabc_);
tmp_vv = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_ykabc_,v_ykabc_);
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, new v_ykabc_: tmp_vv: %0.16f',niteration,tmp_vv)); end;
[tmp_Uv_3_] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,v_ykabc_);
if (flag_verbose>0); disp(sprintf(' %% abs(tmp_Uv_3_): %s',num2str(transpose(abs(tmp_Uv_3_)),' %+0.16f'))); end;
end;%if abs(beta_i_(1+niteration))< 1e-12;
%%%%%%%%;
v_ykabci__(:,1+niteration) = v_ykabc_;
%%;
w_ykabc_ = v_ykabci__(:,1+niteration);
w_ykabc_ = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,w_ykabc_); %<-- Projection. ;
[dvol_a_k_Y_quad_yk_,dtau_euler_polar_a_M_,dtau_euler_azimu_b_M_,dtau_euler_gamma_z_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,w_ykabc_);
eig_ddssnll_lanczos_helper_0; %<-- calculate w_ykabc_ = PHP * v_ykabci__(:,1+niteration) ;
w_ykabc_ = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,U_SmallRotation_Delta_ykabc3__,w_ykabc_); %<-- Projection. ;
%%;
alph_i_(1+niteration) = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,w_ykabc_,v_ykabci__(:,1+niteration));
if (flag_verbose>0); disp(sprintf(' %% %% ni %.2d, alph_i_(1+niteration) = %0.6f',niteration,alph_i_(1+niteration))); end;
w_ykabc_ = w_ykabc_ - alph_i_(1+niteration)*v_ykabci__(:,1+niteration+0);
if niteration>=1; w_ykabc_ = w_ykabc_ - beta_i_(1+niteration)*v_ykabci__(:,1+niteration-1); end;
%%%%;
if (flag_verbose>0);
[tmp_VV_ii__] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,v_ykabci__(:,1:1+niteration),v_ykabci__(:,1:1+niteration));
tmp_str_format = ' % %'; for ni=0:niteration; tmp_str_format = sprintf('%s %%+0.4f',tmp_str_format); end; tmp_str_format = sprintf('%s\\n',tmp_str_format);
disp(sprintf(' %% abs(tmp_VV_ii__): \n %s',num2str(transpose(abs(tmp_VV_ii__(:))),tmp_str_format)));
end;%if (flag_verbose>0);
%%%%%%%%%%%%%%%%;
end;%for niteration=lanczos_n_iteration_pre:lanczos_n_iteration_sum-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_dvol_yk_ = zeros(n_lm_sum,1);
tmp_a_M_ = zeros(n_M,1);
tmp_b_M_ = zeros(n_M,1);
tmp_c_M_ = zeros(n_M,1);
tmp_dvol_yk_(:) = tmp_ykabc_(1:n_lm_sum);
tmp_a_M_(:) = tmp_ykabc_(1*n_lm_sum + 0*n_M + [1:n_M]);
tmp_b_M_(:) = tmp_ykabc_(1*n_lm_sum + 1*n_M + [1:n_M]);
tmp_c_M_(:) = tmp_ykabc_(1*n_lm_sum + 2*n_M + [1:n_M]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykabc_] = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_ykabc_ = zeros(n_lm_sum + 3*n_M,1);
tmp_ykabc_(1:n_lm_sum) = tmp_dvol_yk_(:);
tmp_ykabc_(1*n_lm_sum + 0*n_M + [1:n_M]) = tmp_a_M_(:);
tmp_ykabc_(1*n_lm_sum + 1*n_M + [1:n_M]) = tmp_b_M_(:);
tmp_ykabc_(1*n_lm_sum + 2*n_M + [1:n_M]) = tmp_c_M_(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_dvol_ykn__,tmp_a_Mn__,tmp_b_Mn__,tmp_c_Mn__] = local_ykn_an__bn__cn__from_ykabcn__(n_k_p_r,l_max_,n_M,tmp_ykabcn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykabcn__,2);
tmp_dvol_ykn__ = zeros(n_lm_sum,n_n);
tmp_a_Mn__ = zeros(n_M,n_n);
tmp_b_Mn__ = zeros(n_M,n_n);
tmp_c_Mn__ = zeros(n_M,n_n);
for nn=0:n_n-1;
[tmp_dvol_yk_,tmp_a_M_,tmp_b_M_,tmp_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_ykabc_);
tmp_dvol_ykn__(:,1+nn) = tmp_dvol_yk_;
tmp_a_Mn__(:,1+nn) = tmp_a_M_;
tmp_b_Mn__(:,1+nn) = tmp_b_M_;
tmp_c_Mn__(:,1+nn) = tmp_c_M_;
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykabcn__] = local_ykabcn__from_ykn__an__bn__cn__(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__,tmp_a_Mn__,tmp_b_Mn__,tmp_c_Mn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykabcn__,2);
tmp_ykabcn__ = zeros(n_lm_sum + 3*n_M,n_n);
for nn=0:n_n-1;
tmp_ykabcn__(:,1+nn) = local_ykabc_from_yk_a_b_c_(n_k_p_r,l_max_,n_M,tmp_dvol_ykn__(:,1+nn),tmp_a_Mn__(:,1+nn),tmp_b_Mn__(:,1+nn),tmp_c_Mn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_yk__] = local_yk__from_yk_(n_k_p_r,l_max_,tmp_yk_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk__(1:n_lm,1+nk_p_r) = tmp_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_yk_] = local_yk_from_yk__(n_k_p_r,l_max_,tmp_yk__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
tmp_yk_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
tmp_yk_(1+tmp_index_) = tmp_yk__(1:n_lm,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykn___] = local_ykn___from_ykn__(n_k_p_r,l_max_,tmp_ykn__);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykn__,2);
tmp_ykn___ = zeros(n_lm_max,n_k_p_r,n_n);
for nn=0:n_n-1;
tmp_ykn___(:,:,1+nn) = local_yk__from_yk_(n_k_p_r,l_max_,tmp_ykn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_ykn__] = local_ykn__from_ykn___(n_k_p_r,l_max_,tmp_ykn___);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
n_n = size(tmp_ykn__,2);
tmp_ykn__ = zeros(n_lm_sum,n_n);
for nn=0:n_n-1;
tmp_ykn__(:,1+nn) = local_yk_from_yk__(n_k_p_r,l_max_,tmp_ykn___(:,:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_fg] = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
n_lm_ = (l_max_+1).^2; n_lm_max = max(n_lm_); n_lm_sum = sum(n_lm_); n_lm_csum_ = cumsum([0;n_lm_]); l_max_max = max(l_max_);
scaling_volumetric = (4*pi)^2 * sqrt(pi/2);
[tmp_f_dvol_yk_,tmp_f_a_M_,tmp_f_b_M_,tmp_f_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_f_ykabc_);
[tmp_g_dvol_yk_,tmp_g_a_M_,tmp_g_b_M_,tmp_g_c_M_] = local_yk_a_b_c_from_ykabc_(n_k_p_r,l_max_,n_M,tmp_g_ykabc_);
tmp_f_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_f_dvol_yk_);
tmp_g_dvol_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,tmp_g_dvol_yk_);
tmp_fg_dvol = sum( (conj(tmp_f_dvol_yk__) .* tmp_g_dvol_yk__) * reshape(weight_3d_riesz_k_p_r_,[n_k_p_r,1]) , 'all' ) / scaling_volumetric ;
tmp_fg_a = sum(conj(tmp_f_a_M_) .* tmp_g_a_M_, 'all');
tmp_fg_b = sum(conj(tmp_f_b_M_) .* tmp_g_b_M_, 'all');
tmp_fg_c = sum(conj(tmp_f_c_M_) .* tmp_g_c_M_, 'all');
tmp_fg = tmp_fg_dvol + tmp_fg_a + tmp_fg_b + tmp_fg_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_fg_mn__] = local_fm__bar_dot_gn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcm__,tmp_g_ykabcn__);
n_m = size(tmp_f_ykabcm__,2);
n_n = size(tmp_g_ykabcn__,2);
tmp_fg_mn__ = zeros(n_m,n_n);
for nm=0:n_m-1;
for nn=0:n_n-1;
[tmp_fg_mn__(1+nm,1+nn)] = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcm__(:,1+nm),tmp_g_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
end;%for nm=0:n_m-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_f_ykabc_] = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_);
tmp_ff = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_f_ykabc_);
tmp_f = sqrt(tmp_ff);
tmp_f_ykabc_ = tmp_f_ykabc_/max(1e-12,tmp_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_f_ykabcn__] = local_normalize_fn__(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_f_ykabcn__(:,1+nn) = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__(:,1+nn));
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_g_ykabc_] = local_orthogonalcomplement_gperpf(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_f_ykabc_ = local_normalize_f_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_);
tmp_fg = local_f_bar_dot_g_(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabc_,tmp_g_ykabc_);
tmp_g_ykabc_ = tmp_g_ykabc_ - tmp_f_ykabc_ * tmp_fg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
function [tmp_g_ykabc_] = local_orthogonalcomplement_gperpfn(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__,tmp_g_ykabc_);
n_n = size(tmp_f_ykabcn__,2);
for nn=0:n_n-1;
tmp_g_ykabc_ = local_orthogonalcomplement_gperpf(n_k_p_r,weight_3d_riesz_k_p_r_,l_max_,n_M,tmp_f_ykabcn__(:,1+nn),tmp_g_ykabc_);
end;%for nn=0:n_n-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;




