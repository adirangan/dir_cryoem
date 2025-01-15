flag_verbose=1; flag_disp=1; nf=0;
str_thisfunction = 'test_interpolate_template_6';
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;

%%%%%%%%;
% First form a_x_u_xxx_. ;
%%%%%%%%;
rng(0);
n_x_u = 64*2;
x_p_r_max = 1.0;
x_u_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u));
x_u_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u));
x_u_2_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u)^3;
%%%%;
n_source = 128; n_3 = 3;
delta_3s__ = 0.5*(2*rand(n_3,n_source)-1);
sigma_x = x_p_r_max/16;
g_ = @(x_u_0___,x_u_1___,x_u_2___,delta_3_) 1/sqrt(2*pi)^3 / sigma_x^3 * exp(- ( (x_u_0___ - delta_3_(1+0)).^2 + (x_u_1___ - delta_3_(1+1)).^2 + (x_u_2___ - delta_3_(1+2)).^2 ) / (2*sigma_x^2) ) ;
a_x_u_xxx___ = zeros(n_x_u,n_x_u,n_x_u);
for nsource=0:n_source-1;
a_x_u_xxx___ = a_x_u_xxx___ + g_(x_u_0___,x_u_1___,x_u_2___,delta_3s__(:,1+nsource));
end;%for nsource=0:n_source-1;
a_x_u_xxx_ = reshape(a_x_u_xxx___,[n_xxx_u,1]);
%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_u_1([],a_x_u_xxx___);
title('a_x_u_xxx___','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now set up quadrature on the sphere. ;
%%%%%%%%;
k_int = 48;
k_p_r_max = k_int/(2*pi); 
k_eq_d_double = 1.0/1.0;
k_eq_d = k_eq_d_double/(2*pi);
flag_uniform_over_n_k_p_r = 1;
flag_uniform_over_polar_a = 0;
str_T_vs_L = 'C2';
[ ...
 n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ;
%%%%%%%%;

%%%%%%%%;
% Now form a_k_p_qk_. ;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_qk_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_xxx_(:).*xxx_u_weight_(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_qk_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,a_k_p_qk_.*(2*pi)^3.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
fnorm_disp(flag_verbose,'a_x_u_xxx_',a_x_u_xxx_,'a_x_u_reco_',a_x_u_reco_);
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_reco_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_reco_(:).*xxx_u_weight_(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
fnorm_disp(flag_verbose,'a_k_p_qk_',a_k_p_qk_,'a_k_p_reco_',a_k_p_reco_);
%%%%%%%%;

%%%%%%%%;
% Now form a_k_Y_yk_. ;
%%%%%%%%;
l_max = round(2*pi*k_p_r_max);
n_lm = (l_max+1)^2;
Y_l_val_ = zeros(n_lm,1);
Y_m_val_ = zeros(n_lm,1);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
Y_l_val_(1+na) = l_val;
Y_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
%%%%;
l_max_ = l_max*ones(n_k_p_r,1);
l_max_max = max(l_max_);
m_max_ = -l_max_max:+l_max_max; n_m_max = numel(m_max_);
n_lm_ = n_lm*ones(n_k_p_r,1);
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
Y_l_val__ = repmat(Y_l_val_,[1,n_k_p_r]);
Y_m_val__ = repmat(Y_m_val_,[1,n_k_p_r]);
%%%%%%%%;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
tmp_t = tic();
[ ...
 a_k_Y_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_p_qk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% convert_k_p_to_spharm: %0.6fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ...
 a_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% convert_spharm_to_k_p_4: %0.6fs',tmp_t)); end;
%%%%%%%%;
fnorm_disp(flag_verbose,'a_k_p_qk_',a_k_p_qk_,'a_k_p_reco_',a_k_p_reco_);
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1; clf; figmed; 
subplot(1,2,1); plot(Y_l_val__(:),abs(a_k_Y_yk_(:)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('abs(a_k_Y_yk_)','Interpreter','none');
subplot(1,2,2); plot(Y_m_val__(:),abs(a_k_Y_yk_(:)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('abs(a_k_Y_yk_)','Interpreter','none');
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now generate templates. ;
%%%%%%%%;
viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
template_k_eq_d = -1;
n_w_max = 2*(l_max_max+1); n_w_0in = n_w_max;
n_w = n_w_max; n_w_ = n_w_max*ones(n_k_p_r,1); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
a_k_Y_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_yk_);
%%%%;
flag_check=1;
if flag_check;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max=l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(a_k_Y_yk_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_yk__(1+l_val*(l_val+1)+m_val,1+nk_p_r));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
end;%if flag_check;
%%%%%%%%;
tmp_t = tic();
[ ...
 template_2_wkS___ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_yk__ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
);
n_S = n_viewing_S;
template_2_wkS__ = reshape(template_2_wkS___,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
a_k_Y_lkm___ = [];
[ ...
 template_3_wkS__ ...
] = ...
sph_template_3( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_yk__ ...
,a_k_Y_lkm___ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% sph_template_3: %0.2fs',tmp_t));
disp(sprintf(' %% n_S %d',n_S));
%%%%%%%%;
fnorm_disp(flag_verbose,'template_2_wkS__',template_2_wkS__,'template_3_wkS__',template_3_wkS__);
%%%%%%%%;

%%%%%%%%;
% Define rotation. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;
polar_a_use = +1*pi/5;
azimu_b_use = +1*pi/7;
gamma_z_use = +1*pi/3;
euler_use_ = [gamma_z_use,polar_a_use,azimu_b_use];
R_use__ = Rz(azimu_b_use)*Ry(polar_a_use)*Rz(gamma_z_use);
fnorm_disp(flag_verbose,'R_use__',R_use__,'euler_to_R_0(euler_use_)',euler_to_R_0(euler_use_));
fnorm_disp(flag_verbose,'euler_use_',euler_use_,'R_to_euler_0(R_use__)',R_to_euler_0(R_use__));
%%%%%%%%;
% Now rotate the quadrature-grid. ;
%%%%%%%%;
tmp_qk3__ = cat(2,k_c_0_qk_,k_c_1_qk_,k_c_2_qk_)*transpose(R_use__);
R_k_c_0_qk_ = tmp_qk3__(:,1+0); R_k_c_1_qk_ = tmp_qk3__(:,1+1); R_k_c_2_qk_ = tmp_qk3__(:,1+2);
R_k_p_r01_qk_ = sqrt(R_k_c_0_qk_.^2 + R_k_c_1_qk_.^2);
R_k_p_polar_a_qk_ = atan2(R_k_p_r01_qk_,R_k_c_2_qk_);
R_k_p_azimu_b_qk_ = atan2(R_k_c_1_qk_,R_k_c_0_qk_);
%%%%%%%%;
% Now produce a_R_k_p_qk_. ;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_R_k_p_qk_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_xxx_(:).*xxx_u_weight_(:),-1,1e-12,n_qk,2*pi*R_k_c_0_qk_/eta,2*pi*R_k_c_1_qk_/eta,2*pi*R_k_c_2_qk_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_R_k_p_qk_ time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
% Subselect some templates. ;
% Ensuring that some of the templates have polar_a and azimu_b at multiples of pi. ;
%%%%%%%%;
viewing_gamma_z_S_ = zeros(n_S,1);
n_S_sub = 4;
nS_sub_ = max(0,min(n_S-1,round(n_S*rand(n_S_sub,1))));
tmp_index_a_ = efind(abs(periodize(viewing_polar_a_S_,-pi/4,+pi/4))<1e-2);
tmp_index_b_ = efind(abs(periodize(viewing_azimu_b_S_,-pi/4,+pi/4))<1e-2);
tmp_index_ba_ = intersect(tmp_index_a_,tmp_index_b_);
nS_sub_ = union(tmp_index_ba_(:),nS_sub_);
n_S_sub = numel(nS_sub_);
n_viewing_S_sub = n_S_sub;
viewing_azimu_b_S_sub_ = viewing_azimu_b_S_(1+nS_sub_);
viewing_polar_a_S_sub_ = viewing_polar_a_S_(1+nS_sub_);
viewing_weight_S_sub_ = viewing_weight_S_(1+nS_sub_);
viewing_gamma_z_S_sub_ = viewing_gamma_z_S_(1+nS_sub_);
template_2_wkS_sub__ = template_2_wkS__(:,1+nS_sub_);
template_3_wkS_sub__ = template_3_wkS__(:,1+nS_sub_);
fnorm_disp(flag_verbose,'template_2_wkS_sub__',template_2_wkS_sub__,'template_3_wkS_sub__',template_3_wkS_sub__);
%%%%%%%%;
% Now form templates from a_k_p_quad_. ;
%%%%%%%%;
n_order = 13;
parameter = struct('type','parameter');
parameter.flag_verbose = flag_verbose;
parameter.tolerance_pinv = 1e-6;
tmp_t = tic();
[ ...
 parameter ...
,template_4_wkS_sub__ ...
] = ...
interpolate_template_5( ...
 parameter ...
,n_order ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% interpolate_template S_sub: %0.6fs',tmp_t)); end;
fnorm_disp(flag_verbose,'template_2_wkS_sub__',template_2_wkS_sub__,'template_4_wkS_sub__',template_4_wkS_sub__);
fnorm_disp(flag_verbose,'template_3_wkS_sub__',template_3_wkS_sub__,'template_4_wkS_sub__',template_4_wkS_sub__);
%%%%%%%%;
% Now calculate derivatives. ;
%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose); disp(sprintf(' %% Now testing derivatives: ')); end;
if (flag_verbose); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
%%%%%%%%;
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
[ ...
 parameter ...
,template_mid_x_wkS_sub__ ...
,n_w ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_sub_from_single_shell_sba__ ...
,dtemplateda_mid_wkS_sub__ ...
,dtemplatedb_mid_wkS_sub__ ...
,dtemplatedc_mid_wkS_sub__ ...
,dwS_subda_from_single_shell_sba__ ...
,dwS_subdb_from_single_shell_sba__ ...
,ddtemplatedaa_mid_wkS_sub__ ...
,ddtemplatedab_mid_wkS_sub__ ...
,ddtemplatedac_mid_wkS_sub__ ...
,ddtemplatedbb_mid_wkS_sub__ ...
,ddtemplatedbc_mid_wkS_sub__ ...
,ddtemplatedcc_mid_wkS_sub__ ...
,ddwS_subdaa_from_single_shell_sba__ ...
,ddwS_subdab_from_single_shell_sba__ ...
,ddwS_subdbb_from_single_shell_sba__ ...
] = ...
interpolate_template_5( ...
 parameter ...
,n_order ...
,n_qk ...
,n_qk_csum_ ...
,k_p_r_qk_ ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_3d_k_p_qk_ ...
,weight_shell_qk_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
);
%%%%%%%%;
% Compare with sph_template_3. ;
%%%%%%%%;
a_k_Y_lmk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
tmp_a_k_Y_lm_ = a_k_Y_yk_(1+tmp_index_);
tmp_a_k_Y_lm__ = zeros(1+l_max_max,n_m_max);
l_max = l_max_(1+nk_p_r);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(na==l_val*(l_val+1)+m_val);
tmp_a_k_Y_lm__(1+l_val,1+l_max_max+m_val) = tmp_a_k_Y_lm_(1+na);
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
a_k_Y_lmk___(:,:,1+nk_p_r) = tmp_a_k_Y_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
a_k_Y_lkm___ = permute(a_k_Y_lmk___,[1,3,2]);
%%%%;
tmp_t = tic();
if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlma____','var'); d0W_betazeta_mlma____ = []; end;
if ~exist('d1W_betazeta_mlma____','var'); d1W_betazeta_mlma____ = []; end;
if ~exist('d2W_betazeta_mlma____','var'); d2W_betazeta_mlma____ = []; end;
viewing_gamma_z = 0.0;
[ ...
 template_3_mid_x_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,dtemplate_3da_mid_wkS__ ...
,dtemplate_3db_mid_wkS__ ...
,dtemplate_3dc_mid_wkS__ ...
,d1W_betazeta_mlma____ ...
,ddtemplate_3daa_mid_wkS__ ...
,ddtemplate_3dab_mid_wkS__ ...
,ddtemplate_3dac_mid_wkS__ ...
,ddtemplate_3dbb_mid_wkS__ ...
,ddtemplate_3dbc_mid_wkS__ ...
,ddtemplate_3dcc_mid_wkS__ ...
,d2W_betazeta_mlma____ ...
] = ...
sph_template_3( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_yk__ ...
,a_k_Y_lkm___ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% sph_template_3: %0.2fs',tmp_t));
%%%%;
template_3_mid_x_wkS_sub__ = template_3_mid_x_wkS__(:,1+nS_sub_);
dtemplate_3da_mid_wkS_sub__ = dtemplate_3da_mid_wkS__(:,1+nS_sub_);
dtemplate_3db_mid_wkS_sub__ = dtemplate_3db_mid_wkS__(:,1+nS_sub_);
dtemplate_3dc_mid_wkS_sub__ = dtemplate_3dc_mid_wkS__(:,1+nS_sub_);
ddtemplate_3daa_mid_wkS_sub__ = ddtemplate_3daa_mid_wkS__(:,1+nS_sub_);
ddtemplate_3dab_mid_wkS_sub__ = ddtemplate_3dab_mid_wkS__(:,1+nS_sub_);
ddtemplate_3dac_mid_wkS_sub__ = ddtemplate_3dac_mid_wkS__(:,1+nS_sub_);
ddtemplate_3dbb_mid_wkS_sub__ = ddtemplate_3dbb_mid_wkS__(:,1+nS_sub_);
ddtemplate_3dbc_mid_wkS_sub__ = ddtemplate_3dbc_mid_wkS__(:,1+nS_sub_);
ddtemplate_3dcc_mid_wkS_sub__ = ddtemplate_3dcc_mid_wkS__(:,1+nS_sub_);
%%%%;
%%%%%%%%;
% Display errors. ;
%%%%%%%%;
for nS_sub=0:n_S_sub-1;
disp(sprintf(' %% nS_sub: %d, viewing_polar_a_S_sub_(1+nS_sub)/pi %+0.2f, viewing_azimu_b_S_sub_(1+nS_sub)/pi %+0.2f ',nS_sub,viewing_polar_a_S_sub_(1+nS_sub)/pi,viewing_azimu_b_S_sub_(1+nS_sub)/pi));
%%%%;
fnorm_disp(flag_verbose,'template_2_wkS_sub__(:,1+nS_sub)',template_2_wkS_sub__(:,1+nS_sub),'template_mid_x_wkS_sub__(:,1+nS_sub)',template_mid_x_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'template_2_wkS_sub__(:,1+nS_sub)',template_2_wkS_sub__(:,1+nS_sub),'template_3_mid_x_wkS_sub__(:,1+nS_sub)',template_3_mid_x_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'template_3_wkS_sub__(:,1+nS_sub)',template_3_wkS_sub__(:,1+nS_sub),'template_mid_x_wkS_sub__(:,1+nS_sub)',template_mid_x_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'template_3_wkS_sub__(:,1+nS_sub)',template_3_wkS_sub__(:,1+nS_sub),'template_3_mid_x_wkS_sub__(:,1+nS_sub)',template_3_mid_x_wkS_sub__(:,1+nS_sub));
%%%%;
fnorm_disp(flag_verbose,'+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub),'dtemplateda_mid_wkS_sub__(:,1+nS_sub)',dtemplateda_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub),'dtemplatedb_mid_wkS_sub__(:,1+nS_sub)',dtemplatedb_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub)',-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub),'dtemplatedc_mid_wkS_sub__(:,1+nS_sub)',dtemplatedc_mid_wkS_sub__(:,1+nS_sub));
%%%%;
fnorm_disp(flag_verbose,'+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedaa_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedaa_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedab_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedab_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedac_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedac_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbb_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedbb_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbc_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedbc_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedcc_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedcc_mid_wkS_sub__(:,1+nS_sub));
%%%%;
end;%for nS_sub=0:n_S_sub-1;
%%%%;
%%%%%%%%;

%%%%%%%%;
% Now call interpolate_template_6 using full calculation. ;
%%%%%%%%;
parameter.flag_verbose = 1;
parameter.flag_check = 0;
parameter.flag_parsimonious = 0;
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
[ ...
 parameter ...
,template_mid_x_wkS_sub__ ...
,n_w ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_sub_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dtemplateda_mid_wkS_sub__ ...
,dtemplatedb_mid_wkS_sub__ ...
,dtemplatedc_mid_wkS_sub__ ...
,dwS_subda_from_single_shell_sba__ ...
,dwS_subdb_from_single_shell_sba__ ...
,dtemplateda_rec_wkS_sub__ ...
,dtemplatedb_rec_wkS_sub__ ...
,dtemplatedc_rec_wkS_sub__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddtemplatedaa_mid_wkS_sub__ ...
,ddtemplatedab_mid_wkS_sub__ ...
,ddtemplatedac_mid_wkS_sub__ ...
,ddtemplatedbb_mid_wkS_sub__ ...
,ddtemplatedbc_mid_wkS_sub__ ...
,ddtemplatedcc_mid_wkS_sub__ ...
,ddwS_subdaa_from_single_shell_sba__ ...
,ddwS_subdab_from_single_shell_sba__ ...
,ddwS_subdbb_from_single_shell_sba__ ...
,ddtemplatedaa_rec_wkS_sub__ ...
,ddtemplatedab_rec_wkS_sub__ ...
,ddtemplatedac_rec_wkS_sub__ ...
,ddtemplatedbb_rec_wkS_sub__ ...
,ddtemplatedbc_rec_wkS_sub__ ...
,ddtemplatedcc_rec_wkS_sub__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
interpolate_template_6( ...
 parameter ...
,n_order ...
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
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
);

%%%%%%%%;
% Display errors. ;
%%%%%%%%;
for nS_sub=0:n_S_sub-1;
disp(sprintf(' %% nS_sub: %d, viewing_polar_a_S_sub_(1+nS_sub)/pi %+0.2f, viewing_azimu_b_S_sub_(1+nS_sub)/pi %+0.2f ',nS_sub,viewing_polar_a_S_sub_(1+nS_sub)/pi,viewing_azimu_b_S_sub_(1+nS_sub)/pi));
%%%%;
fnorm_disp(flag_verbose,'+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub),'dtemplateda_rec_wkS_sub__(:,1+nS_sub)',dtemplateda_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub),'dtemplatedb_rec_wkS_sub__(:,1+nS_sub)',dtemplatedb_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub)',-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub),'dtemplatedc_rec_wkS_sub__(:,1+nS_sub)',dtemplatedc_rec_wkS_sub__(:,1+nS_sub));
%%%%;
fnorm_disp(flag_verbose,'+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedaa_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedaa_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedab_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedab_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedac_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedac_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbb_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedbb_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbc_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedbc_rec_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedcc_rec_wkS_sub__(:,1+nS_sub)',ddtemplatedcc_rec_wkS_sub__(:,1+nS_sub));
%%%%;
end;%for nS_sub=0:n_S_sub-1;
%%%%;
disp(sprintf(' %% total:'));
%%%%;
fnorm_disp(flag_verbose,'+dtemplate_3da_mid_wkS_sub__(:,:)',+dtemplate_3da_mid_wkS_sub__(:,:),'dtemplateda_rec_wkS_sub__(:,:)',dtemplateda_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+dtemplate_3db_mid_wkS_sub__(:,:)',+dtemplate_3db_mid_wkS_sub__(:,:),'dtemplatedb_rec_wkS_sub__(:,:)',dtemplatedb_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-dtemplate_3dc_mid_wkS_sub__(:,:)',-dtemplate_3dc_mid_wkS_sub__(:,:),'dtemplatedc_rec_wkS_sub__(:,:)',dtemplatedc_rec_wkS_sub__(:,:));
%%%%;
fnorm_disp(flag_verbose,'+ddtemplate_3daa_mid_wkS_sub__(:,:)',+ddtemplate_3daa_mid_wkS_sub__(:,:),'ddtemplatedaa_rec_wkS_sub__(:,:)',ddtemplatedaa_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dab_mid_wkS_sub__(:,:)',+ddtemplate_3dab_mid_wkS_sub__(:,:),'ddtemplatedab_rec_wkS_sub__(:,:)',ddtemplatedab_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-ddtemplate_3dac_mid_wkS_sub__(:,:)',-ddtemplate_3dac_mid_wkS_sub__(:,:),'ddtemplatedac_rec_wkS_sub__(:,:)',ddtemplatedac_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dbb_mid_wkS_sub__(:,:)',+ddtemplate_3dbb_mid_wkS_sub__(:,:),'ddtemplatedbb_rec_wkS_sub__(:,:)',ddtemplatedbb_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-ddtemplate_3dbc_mid_wkS_sub__(:,:)',-ddtemplate_3dbc_mid_wkS_sub__(:,:),'ddtemplatedbc_rec_wkS_sub__(:,:)',ddtemplatedbc_rec_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dcc_mid_wkS_sub__(:,:)',+ddtemplate_3dcc_mid_wkS_sub__(:,:),'ddtemplatedcc_rec_wkS_sub__(:,:)',ddtemplatedcc_rec_wkS_sub__(:,:));
%%%%;
%%%%%%%%;

%%%%%%%%;
% Now call interpolate_template_6 using parsimonious calculation. ;
%%%%%%%%;
parameter.flag_verbose = 1;
parameter.flag_check = 0;
parameter.flag_parsimonious = 1;
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('a_R_k_p_qk_','var'); a_R_k_p_qk_=[]; end;
if ~exist('ba_from_single_shell_baba__','var'); ba_from_single_shell_baba__=[]; end;
[ ...
 parameter ...
,template_mid_x_wkS_sub__ ...
,n_w ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_sub_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
,wS_from_R_single_shell_sba__ ...
,dtemplateda_mid_wkS_sub__ ...
,dtemplatedb_mid_wkS_sub__ ...
,dtemplatedc_mid_wkS_sub__ ...
,dwS_subda_from_single_shell_sba__ ...
,dwS_subdb_from_single_shell_sba__ ...
,dtemplateda_rec_wkS_sub__ ...
,dtemplatedb_rec_wkS_sub__ ...
,dtemplatedc_rec_wkS_sub__ ...
,dwSda_from_R_single_shell_sba__ ...
,dwSdb_from_R_single_shell_sba__ ...
,ddtemplatedaa_mid_wkS_sub__ ...
,ddtemplatedab_mid_wkS_sub__ ...
,ddtemplatedac_mid_wkS_sub__ ...
,ddtemplatedbb_mid_wkS_sub__ ...
,ddtemplatedbc_mid_wkS_sub__ ...
,ddtemplatedcc_mid_wkS_sub__ ...
,ddwS_subdaa_from_single_shell_sba__ ...
,ddwS_subdab_from_single_shell_sba__ ...
,ddwS_subdbb_from_single_shell_sba__ ...
,ddtemplatedaa_rec_wkS_sub__ ...
,ddtemplatedab_rec_wkS_sub__ ...
,ddtemplatedac_rec_wkS_sub__ ...
,ddtemplatedbb_rec_wkS_sub__ ...
,ddtemplatedbc_rec_wkS_sub__ ...
,ddtemplatedcc_rec_wkS_sub__ ...
,ddwSdaa_from_R_single_shell_sba__ ...
,ddwSdab_from_R_single_shell_sba__ ...
,ddwSdbb_from_R_single_shell_sba__ ...
] = ...
interpolate_template_6( ...
 parameter ...
,n_order ...
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
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
,a_k_p_qk_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in ...
,n_viewing_S_sub ...
,viewing_azimu_b_S_sub_ ...
,viewing_polar_a_S_sub_ ...
,viewing_weight_S_sub_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_sub_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,R_use__ ...
,a_R_k_p_qk_ ...
,ba_from_single_shell_baba__ ...
);

%%%%%%%%;
% Display errors. ;
%%%%%%%%;
for nS_sub=0:n_S_sub-1;
disp(sprintf(' %% nS_sub: %d, viewing_polar_a_S_sub_(1+nS_sub)/pi %+0.2f, viewing_azimu_b_S_sub_(1+nS_sub)/pi %+0.2f ',nS_sub,viewing_polar_a_S_sub_(1+nS_sub)/pi,viewing_azimu_b_S_sub_(1+nS_sub)/pi));
%%%%;
fnorm_disp(flag_verbose,'+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3da_mid_wkS_sub__(:,1+nS_sub),'dtemplateda_mid_wkS_sub__(:,1+nS_sub)',dtemplateda_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub)',+dtemplate_3db_mid_wkS_sub__(:,1+nS_sub),'dtemplatedb_mid_wkS_sub__(:,1+nS_sub)',dtemplatedb_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub)',-dtemplate_3dc_mid_wkS_sub__(:,1+nS_sub),'dtemplatedc_mid_wkS_sub__(:,1+nS_sub)',dtemplatedc_mid_wkS_sub__(:,1+nS_sub));
%%%%;
fnorm_disp(flag_verbose,'+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3daa_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedaa_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedaa_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dab_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedab_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedab_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dac_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedac_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedac_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dbb_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbb_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedbb_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub)',-ddtemplate_3dbc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedbc_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedbc_mid_wkS_sub__(:,1+nS_sub));
fnorm_disp(flag_verbose,'+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub)',+ddtemplate_3dcc_mid_wkS_sub__(:,1+nS_sub),'ddtemplatedcc_mid_wkS_sub__(:,1+nS_sub)',ddtemplatedcc_mid_wkS_sub__(:,1+nS_sub));
%%%%;
end;%for nS_sub=0:n_S_sub-1;
%%%%;
disp(sprintf(' %% total:'));
%%%%;
fnorm_disp(flag_verbose,'+dtemplate_3da_mid_wkS_sub__(:,:)',+dtemplate_3da_mid_wkS_sub__(:,:),'dtemplateda_mid_wkS_sub__(:,:)',dtemplateda_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+dtemplate_3db_mid_wkS_sub__(:,:)',+dtemplate_3db_mid_wkS_sub__(:,:),'dtemplatedb_mid_wkS_sub__(:,:)',dtemplatedb_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-dtemplate_3dc_mid_wkS_sub__(:,:)',-dtemplate_3dc_mid_wkS_sub__(:,:),'dtemplatedc_mid_wkS_sub__(:,:)',dtemplatedc_mid_wkS_sub__(:,:));
%%%%;
fnorm_disp(flag_verbose,'+ddtemplate_3daa_mid_wkS_sub__(:,:)',+ddtemplate_3daa_mid_wkS_sub__(:,:),'ddtemplatedaa_mid_wkS_sub__(:,:)',ddtemplatedaa_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dab_mid_wkS_sub__(:,:)',+ddtemplate_3dab_mid_wkS_sub__(:,:),'ddtemplatedab_mid_wkS_sub__(:,:)',ddtemplatedab_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-ddtemplate_3dac_mid_wkS_sub__(:,:)',-ddtemplate_3dac_mid_wkS_sub__(:,:),'ddtemplatedac_mid_wkS_sub__(:,:)',ddtemplatedac_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dbb_mid_wkS_sub__(:,:)',+ddtemplate_3dbb_mid_wkS_sub__(:,:),'ddtemplatedbb_mid_wkS_sub__(:,:)',ddtemplatedbb_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'-ddtemplate_3dbc_mid_wkS_sub__(:,:)',-ddtemplate_3dbc_mid_wkS_sub__(:,:),'ddtemplatedbc_mid_wkS_sub__(:,:)',ddtemplatedbc_mid_wkS_sub__(:,:));
fnorm_disp(flag_verbose,'+ddtemplate_3dcc_mid_wkS_sub__(:,:)',+ddtemplate_3dcc_mid_wkS_sub__(:,:),'ddtemplatedcc_mid_wkS_sub__(:,:)',ddtemplatedcc_mid_wkS_sub__(:,:));
%%%%%%%%;

