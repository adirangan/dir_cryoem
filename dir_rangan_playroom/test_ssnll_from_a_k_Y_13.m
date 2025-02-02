flag_verbose=1; flag_disp=1; nf=0;
str_thisfunction = 'ssnll_from_a_k_Y_13';
if (flag_verbose>0); disp(sprintf(' %% testing %s',str_thisfunction)); end;
if (flag_verbose>0); disp(sprintf(' %% also see test_ssnll_from_a_k_Y_12')); end;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = 48.0/(2*pi); k_eq_d = 1.0/(2*pi); str_T_vs_L = 'C';
[ ...
 n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
) ;
%%%%;
n_source = 4;
rng(0);
delta_a_c__ = zeros(3,n_source);
delta_b_c__ = zeros(3,n_source);
for nsource=0:n_source-1;
rng(0+nsource);
delta_a_c_ = 2*rand(3,1)-1; delta_a_c_ = delta_a_c_*0.5/k_p_r_max/max(1e-12,fnorm(delta_a_c_)); %<-- ensure small in magnitude. ;
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = 2*rand(3,1)-1; delta_b_c_ = delta_b_c_*0.5/k_p_r_max/max(1e-12,fnorm(delta_b_c_)); %<-- ensure small in magnitude. ;
delta_b_c__(:,1+nsource) = delta_b_c_;
end;%for nsource=0:n_source-1;
a_k_p_form_ = zeros(n_k_all,1);
b_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1) + k_c_2_all_*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = 2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% Now set up spherical-harmonics. ;
%%%%%%%%;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
%%%%;
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
%%%%;
a_k_Y_form_ = zeros(n_lm_sum,1);
b_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c__(:,1+nsource),l_max_);
b_k_Y_form_ = b_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_b_c__(:,1+nsource),l_max_);
end;%for nsource=0:n_source-1;
%%%%;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
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
,a_k_p_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_quad_: %0.16f %%<-- should be <1e-2',fnorm(a_k_Y_form_-a_k_Y_quad_)/fnorm(a_k_Y_form_)));
%%%%;
tmp_t = tic();
[ ...
 b_k_Y_quad_ ...
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
,b_k_p_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% b_k_Y_form_ vs b_k_Y_quad_: %0.16f %%<-- should be <1e-2',fnorm(b_k_Y_form_-b_k_Y_quad_)/fnorm(b_k_Y_form_)));
%%%%;
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
,a_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_p_form_ vs a_k_p_quad_: %0.16f %%<-- should be <1e-2',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
%%%%;
tmp_t = tic;
[ ...
 b_k_p_quad_ ...
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
,b_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% b_k_p_form_ vs b_k_p_quad_: %0.16f %%<-- should be <1e-2',fnorm(b_k_p_form_-b_k_p_quad_)/fnorm(b_k_p_form_)));
%%%%%%%%;
% prepare a_k_Y_form__ and b_k_Y_form__ ;
%%%%%%%%;
a_k_Y_form_yk_ = a_k_Y_form_;
a_k_Y_form_yk__ = zeros(n_lm_max,n_k_p_r);
tmp_t = tic();
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_k_Y_form_yk__(1:n_lm,1+nk_p_r) = a_k_Y_form_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% a_k_Y_form_yk__: %0.6fs',tmp_t)); end;
%%%%;
b_k_Y_form_yk_ = b_k_Y_form_;
b_k_Y_form_yk__ = zeros(n_lm_max,n_k_p_r);
tmp_t = tic();
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
b_k_Y_form_yk__(1:n_lm,1+nk_p_r) = b_k_Y_form_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% b_k_Y_form_yk__: %0.6fs',tmp_t)); end;
%%%%%%%%;
% define rotations. ;
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
% generate templates. ;
%%%%%%%%;
tmp_t = tic();
template_k_eq_d = 1.0/k_p_r_max;
flag_tensor_vs_adap = 1; %<-- tensor grid. ;
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,viewing_k_c_0_all_ ...
,viewing_k_c_1_all_ ...
,viewing_k_c_2_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,template_k_eq_d ...
,str_T_vs_L ...
,flag_tensor_vs_adap ...
) ;
n_S = n_viewing_S;
if (flag_verbose>0); disp(sprintf(' %% n_S %d, n_viewing_polar_a %d, n_viewing_azimu_b [%d,..,%d]',n_S,n_viewing_polar_a,n_viewing_azimu_b_(1+0),n_viewing_azimu_b_(end))); end;
%%%%;
[ ...
 ~ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_ ...
);
%%%%;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,n_w_max ...
,n_S ...
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
,a_k_Y_form_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_k_p_wkS__: %0.6fs',tmp_t)); end;
%%%%%%%%;
% Now get templates for b_k_Y_form. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 T_k_p_wkS__ ...
,n_w_max ...
,n_S ...
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
,b_k_Y_form_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
T_k_p_wkS__ = reshape(T_k_p_wkS__,[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% T_k_p_wkS__: %0.6fs',tmp_t)); end;
%%%%%%%%;
%%%%%%%%;
% Now creating dvol_a_k_Y_form_. ;
%%%%%%%%;
n_source = 4;
rng(1);
delta_dvol_a_c__ = zeros(3,n_source);
for nsource=0:n_source-1;
rng(1+nsource);
delta_dvol_a_c_ = 2*rand(3,1)-1; delta_dvol_a_c_ = delta_dvol_a_c_*0.5/k_p_r_max/max(1e-12,fnorm(delta_dvol_a_c_)); %<-- ensure small in magnitude. ;
delta_dvol_a_c__(:,1+nsource) = delta_dvol_a_c_;
end;%for nsource=0:n_source-1;
dvol_a_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
delta_dvol_a_c_ = delta_dvol_a_c__(:,1+nsource);
dvol_a_k_p_form_ = dvol_a_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_dvol_a_c_(1+0) + k_c_1_all_*delta_dvol_a_c_(1+1) + k_c_2_all_*delta_dvol_a_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%;
dvol_a_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
dvol_a_k_Y_form_ = dvol_a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_dvol_a_c__(:,1+nsource),l_max_);
end;%for nsource=0:n_source-1;
%%%%;
tmp_t = tic;
[ ...
 dvol_a_k_Y_quad_ ...
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
,dvol_a_k_p_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% dvol_a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% dvol_a_k_Y_form_ vs dvol_a_k_Y_quad_: %0.16f %%<-- should be <1e-2',fnorm(dvol_a_k_Y_form_-dvol_a_k_Y_quad_)/fnorm(dvol_a_k_Y_form_)));
%%%%;
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
,dvol_a_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% dvol_a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% dvol_a_k_p_form_ vs dvol_a_k_p_quad_: %0.16f %%<-- should be <1e-2',fnorm(dvol_a_k_p_form_-dvol_a_k_p_quad_)/fnorm(dvol_a_k_p_form_)));
%%%%%%%%;
% prepare dvol_a_k_Y_form__ ;
%%%%%%%%;
dvol_a_k_Y_form_yk_ = dvol_a_k_Y_form_;
dvol_a_k_Y_form_yk__ = zeros(n_lm_max,n_k_p_r);
tmp_t = tic();
for nk_p_r=0:n_k_p_r-1;
n_lm = n_lm_(1+nk_p_r);
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
dvol_a_k_Y_form_yk__(1:n_lm,1+nk_p_r) = dvol_a_k_Y_form_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dvol_a_k_Y_form_yk__: %0.6fs',tmp_t)); end;
%%%%%%%%;
% generate templates. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 dvol_S_k_p_wkS__ ...
,n_w_max ...
,n_S ...
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
,dvol_a_k_Y_form_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
dvol_S_k_p_wkS__ = reshape(dvol_S_k_p_wkS__,[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% dvol_S_k_p_wkS__: %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% set up CTF and eta. ;
% This does not test anisotropic CTF. ;
% Nor does this test anisotropic eta. ;
%%%%%%%%;
n_alpha = 3; CTF_alpha_ = [0.10;0.15;0.30];
n_power = 2; eta_power_ = [0.50;0.75];
n_CTF = n_alpha*n_power;
n_eta = n_alpha*n_power;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
eta_k_p_wke__ = zeros(n_w_sum,n_eta);
for nalpha=0:n_alpha-1;
CTF_alpha = CTF_alpha_(1+nalpha);
for npower=0:n_power-1;
eta_power = eta_power_(1+npower);
nCTF = nalpha+npower*n_alpha;
neta = nalpha+npower*n_alpha;
CTF_k_p_wk_ = (reshape(repmat(reshape(besselj(0,CTF_alpha*k_p_r_),[1,n_k_p_r]),[n_w_max,1]),[n_w_sum,1])).^(1-eta_power);
CTF_k_p_wkC__(:,1+nCTF) = CTF_k_p_wk_;
eta_k_p_wk_ = (reshape(repmat(reshape(besselj(0,CTF_alpha*k_p_r_),[1,n_k_p_r]),[n_w_max,1]),[n_w_sum,1])).^(0+eta_power);
eta_k_p_wke__(:,1+neta) = eta_k_p_wk_;
end;%for npower=0:n_power-1;
end;%for nalpha=0:n_alpha-1;
%%%%;

%%%%%%%%;
% testing full on-grid ssnll. ;
%%%%%%%%;
n_M = n_S;
n_N = floor(n_S/2);
rng(0);
index_nN_from_nM_ = floor(n_N*rand(n_M,1));
index_nM_from_nN__ = cell(n_N,1);
weight_imagecount_N_ = zeros(n_N,1);
for nN=0:n_N-1;
tmp_index_ = efind(index_nN_from_nM_==nN);
index_nM_from_nN__{1+nN} = tmp_index_;
weight_imagecount_N_(1+nN) = numel(tmp_index_);
end;%for nN=0:n_N-1;
weight_imagecount_M_ = ones(n_M,1);
%%%%;
index_nCTF_from_nN_ = zeros(n_N,1);
for nN=0:n_N-1;
index_nCTF_from_nN_(1+nN) = mod(nN,n_CTF);
end;%for nN=0:n_N-1;
index_nCTF_from_nM_ = index_nCTF_from_nN_(1+index_nN_from_nM_);
%%%%;
index_neta_from_nN_ = zeros(n_N,1);
for nN=0:n_N-1;
index_neta_from_nN_(1+nN) = mod(nN,n_eta);
end;%for nN=0:n_N-1;
index_neta_from_nM_ = index_neta_from_nN_(1+index_nN_from_nM_);
%%%%;
euler_polar_a_N_ = viewing_polar_a_S_(1:n_N);
euler_azimu_b_N_ = viewing_azimu_b_S_(1:n_N);
euler_gamma_z_N_ = 2*pi*rand(n_N,1);
euler_polar_a_M_ = euler_polar_a_N_(1+index_nN_from_nM_);
euler_azimu_b_M_ = euler_azimu_b_N_(1+index_nN_from_nM_);
euler_gamma_z_M_ = euler_gamma_z_N_(1+index_nN_from_nM_);
N_k_p_wkN__ = T_k_p_wkS__(:,1:n_N);
for nN=0:n_N-1;
N_k_p_wkN__(:,1+nN) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,N_k_p_wkN__(:,1+nN),+euler_gamma_z_N_(1+nN));
end;%for nN=0:n_N-1;
M_k_p_wkM__ = N_k_p_wkN__(:,1+index_nN_from_nM_);
%%%%;
euler_polar_a_N_ = viewing_polar_a_S_(1:n_N);
euler_azimu_b_N_ = viewing_azimu_b_S_(1:n_N);
euler_gamma_z_N_ = 2*pi*rand(n_N,1);
euler_polar_a_M_ = euler_polar_a_N_(1+index_nN_from_nM_);
euler_azimu_b_M_ = euler_azimu_b_N_(1+index_nN_from_nM_);
euler_gamma_z_M_ = euler_gamma_z_N_(1+index_nN_from_nM_);
dtau_euler_polar_a_N_ = 1*pi*rand(n_N,1);
dtau_euler_azimu_b_N_ = 2*pi*rand(n_N,1);
dtau_euler_gamma_z_N_ = 2*pi*rand(n_N,1);
dtau_euler_polar_a_M_ = dtau_euler_polar_a_N_(1+index_nN_from_nM_);
dtau_euler_azimu_b_M_ = dtau_euler_azimu_b_N_(1+index_nN_from_nM_);
dtau_euler_gamma_z_M_ = dtau_euler_gamma_z_N_(1+index_nN_from_nM_);
%%%%;

if ~exist('V_lmm___','var'); V_lmm___=[]; end;
if ~exist('L_lm__','var'); L_lm__=[]; end;
if ~exist('d0W_betazeta_mlma____','var'); d0W_betazeta_mlma____=[]; end;
if ~exist('d1W_betazeta_mlma____','var'); d1W_betazeta_mlma____=[]; end;
if ~exist('d2W_betazeta_mlma____','var'); d2W_betazeta_mlma____=[]; end;

%%%%;
parameter_ssnll = struct('type','parameter');
[ ...
 parameter ...
,ssnll_M_ ...
,ssnll_M ...
,S_k_p_wkS__ ...
,dvol_ssnll_M_ ...
,dvol_ssnll_M ...
,dvol_S_k_p_wkS__ ...
,dvol_dvol_ssnll_M ...
,dtau_ssnll_M3__ ...
,dtau_ssnll_M ...
,dtau_S_k_p_wkS3___ ...
,dtau_dvol_ssnll_M3__ ...
,dtau_dvol_ssnll_M ...
,dtau_dvol_S_k_p_wkS3___ ...
,dtau_dtau_ssnll_M33___ ...
,dtau_dtau_ssnll_M ...
,dtau_dtau_S_k_p_wkS33____ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
] = ...
ssnll_from_a_k_Y_13( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_form_yk_ ...
,dvol_a_k_Y_form_yk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_wkS__ ...
,dvol_S_k_p_wkS__ ...
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
,weight_imagecount_M_ ...
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
);
%%%%;

%%%%;
parameter_ssnll = struct('type','parameter');
[ ...
 parameter ...
,ssnll_N_ ...
,ssnll_N ...
,S_k_p_wkS__ ...
,dvol_ssnll_N_ ...
,dvol_ssnll_N ...
,dvol_S_k_p_wkS__ ...
,dvol_dvol_ssnll_N ...
,dtau_ssnll_N3__ ...
,dtau_ssnll_N ...
,dtau_S_k_p_wkS3___ ...
,dtau_dvol_ssnll_N3__ ...
,dtau_dvol_ssnll_N ...
,dtau_dvol_S_k_p_wkS3___ ...
,dtau_dtau_ssnll_N33___ ...
,dtau_dtau_ssnll_N ...
,dtau_dtau_S_k_p_wkS33____ ...
,V_lmm___ ...
,L_lm__ ...
,d0W_betazeta_mlma____ ...
,d1W_betazeta_mlma____ ...
,d2W_betazeta_mlma____ ...
] = ...
ssnll_from_a_k_Y_13( ...
 parameter_ssnll ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,l_max_ ...
,a_k_Y_form_yk_ ...
,dvol_a_k_Y_form_yk_ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_S ...
,S_k_p_wkS__ ...
,dvol_S_k_p_wkS__ ...
,[] ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,n_N ...
,weight_imagecount_N_ ...
,N_k_p_wkN__ ...
,index_nCTF_from_nN_ ...
,CTF_k_p_wkC__ ...
,index_neta_from_nN_ ...
,eta_k_p_wke__ ...
,euler_polar_a_N_ ...
,euler_azimu_b_N_ ...
,euler_gamma_z_N_ ...
,dtau_euler_polar_a_N_ ...
,dtau_euler_azimu_b_N_ ...
,dtau_euler_gamma_z_N_ ...
);
%%%%;

disp(sprintf(' %% ssnll_N_(1+index_nN_from_nM_) vs ssnll_M_: %0.16f',fnorm(ssnll_N_(1+index_nN_from_nM_) - ssnll_M_)/max(1e-12,fnorm(ssnll_N_(1+index_nN_from_nM_)))));
disp(sprintf(' %% ssnll_N vs ssnll_M: %0.16f',fnorm(ssnll_N - ssnll_M)/max(1e-12,fnorm(ssnll_N))));
disp(sprintf(' %% dvol_ssnll_N_(1+index_nN_from_nM_) vs dvol_ssnll_M_: %0.16f',fnorm(dvol_ssnll_N_(1+index_nN_from_nM_) - dvol_ssnll_M_)/max(1e-12,fnorm(dvol_ssnll_N_(1+index_nN_from_nM_)))));
disp(sprintf(' %% dvol_ssnll_N vs dvol_ssnll_M: %0.16f',fnorm(dvol_ssnll_N - dvol_ssnll_M)/max(1e-12,fnorm(dvol_ssnll_N))));
disp(sprintf(' %% dvol_dvol_ssnll_N vs dvol_dvol_ssnll_M: %0.16f',fnorm(dvol_dvol_ssnll_N - dvol_dvol_ssnll_M)/max(1e-12,fnorm(dvol_dvol_ssnll_N))));
disp(sprintf(' %% dtau_ssnll_N3__(1+index_nN_from_nM_,:) vs dtau_ssnll_M3__: %0.16f',fnorm(dtau_ssnll_N3__(1+index_nN_from_nM_,:) - dtau_ssnll_M3__)/max(1e-12,fnorm(dtau_ssnll_N3__(1+index_nN_from_nM_,:)))));
disp(sprintf(' %% dtau_ssnll_N vs dtau_ssnll_M: %0.16f',fnorm(dtau_ssnll_N - dtau_ssnll_M)/max(1e-12,fnorm(dtau_ssnll_N))));
disp(sprintf(' %% dtau_dvol_ssnll_N3__(1+index_nN_from_nM_,:) vs dtau_dvol_ssnll_M3__: %0.16f',fnorm(dtau_dvol_ssnll_N3__(1+index_nN_from_nM_,:) - dtau_dvol_ssnll_M3__)/max(1e-12,fnorm(dtau_dvol_ssnll_N3__(1+index_nN_from_nM_,:)))));
disp(sprintf(' %% dtau_dvol_ssnll_N vs dtau_dvol_ssnll_M: %0.16f',fnorm(dtau_dvol_ssnll_N - dtau_dvol_ssnll_M)/max(1e-12,fnorm(dtau_dvol_ssnll_N))));
disp(sprintf(' %% dtau_dtau_ssnll_N33___(1+index_nN_from_nM_,:,:) vs dtau_dtau_ssnll_M33___: %0.16f',fnorm(dtau_dtau_ssnll_N33___(1+index_nN_from_nM_,:,:) - dtau_dtau_ssnll_M33___)/max(1e-12,fnorm(dtau_dtau_ssnll_N33___(1+index_nN_from_nM_,:,:)))));
disp(sprintf(' %% dtau_dtau_ssnll_N vs dtau_dtau_ssnll_M: %0.16f',fnorm(dtau_dtau_ssnll_N - dtau_dtau_ssnll_M)/max(1e-12,fnorm(dtau_dtau_ssnll_N))));

