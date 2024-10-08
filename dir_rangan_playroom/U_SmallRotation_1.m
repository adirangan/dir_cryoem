function ...
[ ...
 parameter ...
,U_tilde_SmallRotation_Delta_ykabcs__ ...
,U_SmallRotation_Delta_ykabcs__ ...
,S_SmallRotation_Delta_s_ ...
,V_SmallRotation_Delta_ss__ ...
] = ...
U_SmallRotation_1( ...
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
);

str_thisfunction = 'U_SmallRotation_1';

if nargin<1;
disp(sprintf(' %% testing %s',str_thisfunction));

disp('returning'); return;
end;%if nargin<1;

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

if isempty(parameter);
parameter = struct('type','parameter');
end;%if isempty(parameter);
%%%%%%%%;
if (~isfield(parameter,'flag_verbose')); parameter.flag_verbose = 0; end; %<-- parameter_bookmark. ;
flag_verbose = parameter.flag_verbose;
if (~isfield(parameter,'flag_check')); parameter.flag_check = 0; end; %<-- parameter_bookmark. ;
flag_check = parameter.flag_check;
if (~isfield(parameter,'tolerance_master')); parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = parameter.tolerance_master;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

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
n_ykabc = n_lm_sum + n_M*3;
weight_3d_riesz_ykabc_ = cat(1,weight_3d_riesz_yk_/scaling_volumetric,ones(3*n_M,1));
numerator_root_weight_3d_riesz_ykabc_ = reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]);
denomator_root_weight_3d_riesz_ykabc_ = 1./max(1e-12,reshape(sqrt(weight_3d_riesz_ykabc_),[n_ykabc,1]));
%%%%%%%%;

%%%%%%%%;
% Now determine space of all such alignments (should have dimension of so3). ;
%%%%%%%%;
n_SmallRotation = 9;
SmallRotation_dvol_a_k_Y_quad_yks__ = zeros(n_lm_sum,n_SmallRotation);
SmallRotation_dtau_euler_polar_a_Ms__ = zeros(n_M,n_SmallRotation);
SmallRotation_dtau_euler_azimu_b_Ms__ = zeros(n_M,n_SmallRotation);
SmallRotation_dtau_euler_gamma_z_Ms__ = zeros(n_M,n_SmallRotation);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nSmallRotation=0:n_SmallRotation-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% nSmallRotation = %d/%d',nSmallRotation,n_SmallRotation)); end;
rng(1+nSmallRotation);
SmallRotation_polar_a = 1*pi*rand();
SmallRotation_azimu_b = 2*pi*rand();
SmallRotation_gamma_z = 2*pi*rand();
SmallRotation_euler_ = [-SmallRotation_azimu_b,-SmallRotation_polar_a,+SmallRotation_gamma_z]; %<-- note the reversal in ordering. ;
SmallRotation_R__ = euler_to_R(SmallRotation_euler_);
euler_azimu_b_one_ = [+1,+0,+0];
tmp_f = +cos(SmallRotation_polar_a);
tmp_g = -sin(SmallRotation_polar_a)*sin(-SmallRotation_gamma_z); %<-- note flipped sign of SmallRotation_gamma_z. ;
tmp_h = -sin(SmallRotation_polar_a)*cos(-SmallRotation_gamma_z);
SmallRotation_Delta_R__ = [ 0 , -tmp_f , -tmp_g ; +tmp_f , 0 , -tmp_h ; +tmp_g , +tmp_h , 0 ];
SmallRotation_dvol_a_k_Y_quad_yk_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,bsxfun(@times,rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_quad_yk_,-flip(SmallRotation_euler_)),+i*Y_m_val_),+SmallRotation_euler_);
%%%%;
tmp_f = +SmallRotation_Delta_R__(1+1,1+0);
tmp_g = +SmallRotation_Delta_R__(1+2,1+0);
tmp_h = +SmallRotation_Delta_R__(1+2,1+1);
SmallRotation_dtau_euler_polar_a_M_ = zeros(n_M,1);
SmallRotation_dtau_euler_azimu_b_M_ = zeros(n_M,1);
SmallRotation_dtau_euler_gamma_z_M_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_polar_a_ori = +euler_polar_a_M_(1+nM);
tmp_azimu_b_ori = +euler_azimu_b_M_(1+nM);
tmp_gamma_z_ori = -euler_gamma_z_M_(1+nM);
%tmp_da_mid = +(sin(tmp_polar_a_ori)*(tmp_g*cos(tmp_azimu_b_ori) + tmp_h*sin(tmp_azimu_b_ori)))/(1 - cos(tmp_polar_a_ori)^2)^(1/2); %<-- cancel out. ;
tmp_da_mid = +(tmp_g*cos(tmp_azimu_b_ori) + tmp_h*sin(tmp_azimu_b_ori));
%tmp_db_mid = -(cos(tmp_azimu_b_ori)*(tmp_f*sin(tmp_polar_a_ori) - tmp_h*cos(tmp_polar_a_ori)*cos(tmp_azimu_b_ori) + tmp_g*cos(tmp_polar_a_ori)*sin(tmp_azimu_b_ori)))/(cos(tmp_azimu_b_ori)*sin(tmp_polar_a_ori)); %<-- remove lower order terms from denominator. ;
tmp_db_mid = -(tmp_f*sin(tmp_polar_a_ori) - tmp_h*cos(tmp_polar_a_ori)*cos(tmp_azimu_b_ori) + tmp_g*cos(tmp_polar_a_ori)*sin(tmp_azimu_b_ori))/sin(tmp_polar_a_ori); %<-- remove lower order terms from denominator. ;
%tmp_dc_mid = -(cos(tmp_gamma_z_ori)*(tmp_h*cos(tmp_azimu_b_ori) - tmp_g*sin(tmp_azimu_b_ori)))/(cos(tmp_gamma_z_ori)*sin(tmp_polar_a_ori)); %<-- remove lower order terms from denominator. ;
tmp_dc_mid = -(tmp_h*cos(tmp_azimu_b_ori) - tmp_g*sin(tmp_azimu_b_ori))/sin(tmp_polar_a_ori); %<-- remove lower order terms from denominator. ;
SmallRotation_dtau_euler_polar_a_M_(1+nM) = +tmp_da_mid;
SmallRotation_dtau_euler_azimu_b_M_(1+nM) = +tmp_db_mid;
SmallRotation_dtau_euler_gamma_z_M_(1+nM) = -tmp_dc_mid;
end;%for nM=0:n_M-1;
%%%%;
SmallRotation_dvol_a_k_Y_quad_yks__(:,1+nSmallRotation) = SmallRotation_dvol_a_k_Y_quad_yk_;
SmallRotation_dtau_euler_polar_a_Ms__(:,1+nSmallRotation) = SmallRotation_dtau_euler_polar_a_M_;
SmallRotation_dtau_euler_azimu_b_Ms__(:,1+nSmallRotation) = SmallRotation_dtau_euler_azimu_b_M_;
SmallRotation_dtau_euler_gamma_z_Ms__(:,1+nSmallRotation) = SmallRotation_dtau_euler_gamma_z_M_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nSmallRotation=0:n_SmallRotation-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
SmallRotation_Delta_ykabcs__ = cat( ...
 1 ...
,SmallRotation_dvol_a_k_Y_quad_yks__ ...
,SmallRotation_dtau_euler_polar_a_Ms__ ...
,SmallRotation_dtau_euler_azimu_b_Ms__ ...
,SmallRotation_dtau_euler_gamma_z_Ms__ ...
);
[U_SmallRotation_Delta_ykabcs__,S_SmallRotation_Delta_ss__,V_SmallRotation_Delta_ss__] = svds(bsxfun(@times,sqrt(weight_3d_riesz_ykabc_),SmallRotation_Delta_ykabcs__),n_SmallRotation);
S_SmallRotation_Delta_s_ = diag(S_SmallRotation_Delta_ss__);
%ctranspose(U_SmallRotation_Delta_ykabcs__)*U_SmallRotation_Delta_ykabcs__,;
U_SmallRotation_Delta_ykabcs__ = bsxfun(@times,1./max(1e-12,sqrt(weight_3d_riesz_ykabc_)),U_SmallRotation_Delta_ykabcs__);
if (flag_verbose>0); disp(sprintf(' %% S_SmallRotation_Delta_s_: %s',num2str(transpose(S_SmallRotation_Delta_s_),' %+0.6f'))); end;
U_tilde_SmallRotation_Delta_ykabcs__ = bsxfun(@times,numerator_root_weight_3d_riesz_ykabc_,U_SmallRotation_Delta_ykabcs__);

if flag_check;
for ns=0:3-1;
tmp_t = tic();
viewing_gamma_z_S_ = zeros(n_S,1);
SmallRotation_dvol_a_k_Y_quad_yk_ = U_SmallRotation_Delta_ykabcs__(1:n_lm_sum,1+ns);
SmallRotation_dtau_euler_polar_a_M_ = U_SmallRotation_Delta_ykabcs__(1*n_lm_sum + 0*n_M + [1:n_M],1+ns);
SmallRotation_dtau_euler_azimu_b_M_ = U_SmallRotation_Delta_ykabcs__(1*n_lm_sum + 1*n_M + [1:n_M],1+ns);
SmallRotation_dtau_euler_gamma_z_M_ = U_SmallRotation_Delta_ykabcs__(1*n_lm_sum + 2*n_M + [1:n_M],1+ns);
parameter_ssnll = struct('type','parameter');
[ ...
 ~ ...
,ssnll_q2d_M_ ...
,ssnll_q2d ...
,S_k_p_q2d_q2d_wkS__ ...
,dvol_ssnll_q2d_M_ ...
,dvol_ssnll_q2d ...
,dvol_S_k_p_q2d_wkS__ ...
,dvol_dvol_ssnll_q2d ...
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
,a_k_Y_quad_yk_ ...
,SmallRotation_dvol_a_k_Y_quad_yk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,[] ...
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
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,[] ...
,[] ...
,euler_polar_a_M_ ...
,euler_azimu_b_M_ ...
,euler_gamma_z_M_ ...
,SmallRotation_dtau_euler_polar_a_M_ ...
,SmallRotation_dtau_euler_azimu_b_M_ ...
,SmallRotation_dtau_euler_gamma_z_M_ ...
);
tmp_t = toc(tmp_t);
if (flag_verbose>0); disp(sprintf(' %% ssnll_from_a_k_Y_12 (yes derivative): %0.6fs',tmp_t)); end;
%%%%%%%%;
SmallRotation_dtau_M3__ = [ ...
,SmallRotation_dtau_euler_polar_a_M_ ...
,SmallRotation_dtau_euler_azimu_b_M_ ...
,SmallRotation_dtau_euler_gamma_z_M_ ...
] ;
dtau_ssnll_q2d_M_ = sum(dtau_ssnll_q2d_M3__.*SmallRotation_dtau_M3__,[2]);
if (flag_verbose>0); disp(sprintf(' %% ns=%d/3 dtau_ssnll_q2d + dvol_ssnll_q2d = %0.6f',ns,dtau_ssnll_q2d + dvol_ssnll_q2d)); end;
end;%for ns=0:3-1;
end;%if flag_check;
%%%%%%%%;

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;

