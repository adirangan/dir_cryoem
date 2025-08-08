function ...
[ ...
 parameter ...
] = ...
test_transforms_8_sub( ...
 parameter ...
);

%%%%%%%%;
% Sets up a simple volume (on a spherical grid), ;
% then generates templates and calculates inner-products. ;
% These calculations consider anisotropic CTF. ;
%%%%%%%%;

str_thisfunction = 'test_transforms_8_sub';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nargin<1;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
string_user = 'rangan';
flag_verbose = 1; %<-- verbosity level. ;
flag_disp=1; nf=0; %<-- display level. ;
flag_replot=1;
k_int = 16; %<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! ;
k_eq_d_double = 1.0; %<-- prefactor for k_eq_d, determines density of sampling in frequency-space. ;
template_k_eq_d_double = 1.0; %<-- prefactor for template_k_eq_d, determines density of viewing-angles on the sphere. ;
n_w_int = 1.0; %<-- prefactor for n_w_max, determines the number of distinct angles (i.e., n_gamma_z) used in frequency-space 2d-polar-grid. ;
parameter = struct('type','parameter');
parameter.string_root = string_root;
parameter.string_user = string_user;
parameter.flag_verbose = flag_verbose;
parameter.flag_disp = flag_disp;
parameter.flag_replot = flag_replot;
parameter.k_int = k_int;
parameter.k_eq_d_double = k_eq_d_double;
parameter.template_k_eq_d_double = template_k_eq_d_double;
parameter.n_w_int = n_w_int;
test_transforms_8_sub(parameter);
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

if isempty(parameter); parameter=struct('type','parameter'); end;
if ~isfield(parameter,'flag_verbose'); parameter.flag_verbose=0; end;
flag_verbose=parameter.flag_verbose;
if ~isfield(parameter,'flag_disp'); parameter.flag_disp=0; end;
flag_disp=parameter.flag_disp; nf=0;
if ~isfield(parameter,'flag_replot'); parameter.flag_replot=0; end;
flag_replot=parameter.flag_replot;
if ~isfield(parameter,'string_root'); parameter.string_root = sprintf('data'); end;
string_root=parameter.string_root;
if ~isfield(parameter,'string_user'); parameter.string_user = sprintf('rangan'); end;
string_user=parameter.string_user;
if ~isfield(parameter,'dir_base'); parameter.dir_base = sprintf('/%s/%s/dir_cryoem/dir_rangan_playroom',string_root,string_user); end;
dir_base=parameter.dir_base;
if ~isfield(parameter,'dir_jpg'); parameter.dir_jpg = sprintf('%s/dir_test_transforms_jpg',dir_base); end;
dir_jpg=parameter.dir_jpg;
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
if ~isfield(parameter,'k_int'); parameter.k_int=32; end;
k_int=parameter.k_int;
if ~isfield(parameter,'k_eq_d_double'); parameter.k_eq_d_double=1.0; end;
k_eq_d_double=parameter.k_eq_d_double;
if ~isfield(parameter,'template_k_eq_d_double'); parameter.template_k_eq_d_double=1.0; end;
template_k_eq_d_double=parameter.template_k_eq_d_double;
if ~isfield(parameter,'n_w_int'); parameter.n_w_int=1.0; end;
n_w_int=parameter.n_w_int;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
dh3d_ = @(kd) 12*pi*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% Set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
flag_uniform_over_n_k_p_r = 1; %<-- we use the same discretization on each shell. ;
flag_uniform_over_polar_a = 0; %<-- however, we allow for different discretizations on each latitude. ;
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
,str_L ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ;
%%%%%%%%;

%%%%%%%%;
% Now define two functions on the sphere, ;
% each a sum of a few plane-waves. ;
%%%%%%%%;
delta_a_c_3s__ = [  ...
 +1.5 , -0.5 , +0.1 , -0.8 ...
;-0.5 , -1.5 , -0.3 , +0.4 ...
;+0.3 , +2.0 , +0.7 , +1.6 ...
] / 2 / k_p_r_max ;
n_source = size(delta_a_c_3s__,2);
a_k_p_form_ = zeros(n_qk,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_a_c_(1+0) + k_c_1_qk_*delta_a_c_(1+1) + k_c_2_qk_*delta_a_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%%%%%;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now test k-quadrature on sphere. ;
%%%%%%%%;
I_a_quad = sum(a_k_p_form_.*weight_3d_k_p_qk_);
I_a_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h3d_(tmp_kd)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'I_a_form',I_a_form,'I_a_quad',I_a_quad,' %%<-- should be <1e-6');
%%%%%%%%;
a_k_p_l2_quad = sum(conj(a_k_p_form_).*a_k_p_form_.*weight_3d_k_p_qk_);
a_k_p_l2_form = 0;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
delta_a_c_0_ = delta_a_c_3s__(:,1+nsource0);
delta_a_c_1_ = delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_0_ - delta_a_c_1_);
tmp_h3d = 4*pi/3; if abs(tmp_kd)>1e-12; tmp_h3d = h3d_(tmp_kd); end;
a_k_p_l2_form = a_k_p_l2_form + tmp_h3d*k_p_r_max^3;
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
fnorm_disp(flag_verbose,'a_k_p_l2_form',a_k_p_l2_form,'a_k_p_l2_quad',a_k_p_l2_quad,' %%<-- should be <1e-6');
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
%%%%%%%%;

%%%%%%%%;
% define rotations in 2d and 3d. ;
%%%%%%%%;
R2 = @(gamma_z) ...
[ +cos(gamma_z) -sin(gamma_z) ; ...
  +sin(gamma_z) +cos(gamma_z) ; ...
] ;
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

%%%%%%%%;
% Sample template viewing-angles. ;
%%%%%%%%;
template_k_eq_d = template_k_eq_d_double/max(1e-12,k_p_r_max);
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,template_k_eq_d ...
,str_L ...
,flag_uniform_over_polar_a ...
);
n_S = n_viewing_S;

%%%%%%%%;
% Now step through and reconstitute the templates themselves. ;
%%%%%%%%;
S_k_p_l2_quad_S_ = zeros(n_S,1);
S_k_p_l2_form_S_ = zeros(n_S,1);
S_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
S_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
S_k_p_wk_ = S_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
S_k_p_l2_quad_S_(1+nS) = sum(conj(S_k_p_wk_).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
tmp_delta_a_c_0_ = tmp_R__*delta_a_c_3s__(:,1+nsource0);
tmp_delta_a_c_1_ = tmp_R__*delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(tmp_delta_a_c_0_(1:2) - tmp_delta_a_c_1_(1:2));
tmp_h2d = (2*pi)^2; if abs(tmp_kd)>1e-12; tmp_h2d = h2d_(tmp_kd); end;
S_k_p_l2_form_S_(1+nS) = S_k_p_l2_form_S_(1+nS) + tmp_h2d/(2*pi)^2 * (pi*k_p_r_max^2);
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
S_k_p_wkS__(:,1+nS) = S_k_p_wk_;
end;%for nS=0:n_S-1;
%%%%%%%%;

%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
S_k_q_wkS__(:,1+nS) = S_k_q_wk_;
end;%for nS=0:n_S-1;
flag_check=0;
if flag_check;
tmp_ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__),[n_w_sum,n_S]); %<-- this accomplishes the same thing. ;
fnorm_disp(flag_verbose,'S_k_q_wkS__',S_k_q_wkS__,'tmp_',tmp_,' %%<-- should be zero');
end;%if flag_check;
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),[1,3,2]);
%%%%%%%%;

flag_delta = 1;
%%%%%%%%;
if ~flag_delta; delta_r_max = 0.0/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested =  1; end;
if  flag_delta; delta_r_max = 0.5/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested = 128; end;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
%%%%%%%%;
n_delta_v = FTK.n_delta_v;

flag_CTF = 1;
%%%%%%%%;
% Now we generate images to test the inner-product calculation. ;
% These images are copied from the templates, and then multiplied by a CTF. ;
% We also shift each image by a fixed displacement. ;
%%%%%%%%;
n_CTF = 3;
CTF_phi_C_ = pi*[+2/3;-1/5;+6/7];
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if ~flag_CTF; CTF_k_p_wkC__(:,1+nCTF) = ones(n_w_sum,1); end;
if  flag_CTF; CTF_k_p_wkC__(:,1+nCTF) = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTF_phi_C_(1+nCTF)); end;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
tmp_t = tic();
n_M = 2*n_S; %<-- pick something bigger than n_S;
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
M_k_p_wkM__ = zeros(n_w_sum,n_M);
index_nCTF_from_nM_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
euler_azimu_b_true_M_ = zeros(n_M,1);
euler_polar_a_true_M_ = zeros(n_M,1);
index_nw_from_nM_ = zeros(n_M,1);
euler_gamma_z_true_M_ = zeros(n_M,1);
index_nd_from_nM_ = zeros(n_M,1);
image_delta_x_true_M_ = zeros(n_M,1);
image_delta_y_true_M_ = zeros(n_M,1);
rng(0);
for nM=0:n_M-1;
nCTF = mod(nM,n_CTF);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
index_nCTF_from_nM_(1+nM) = nCTF;
nS = max(0,min(n_S-1,mod(nM,n_S)));
index_nS_from_nM_(1+nM) = nS;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
azimu_b = viewing_azimu_b_S_(1+nS);
polar_a = viewing_polar_a_S_(1+nS);
index_gamma_z = periodize(nM,0,n_gamma_z);
gamma_z = gamma_z_(1+index_gamma_z);
index_nw_from_nM_(1+nM) = index_gamma_z;
euler_azimu_b_true_M_(1+nM) = azimu_b;
euler_polar_a_true_M_(1+nM) = polar_a;
euler_gamma_z_true_M_(1+nM) = gamma_z;
index_nd_from_nM_(1+nM) = periodize(nM,0,n_delta_v);
delta_x = FTK.delta_x_(1+index_nd_from_nM_(1+nM)); delta_y = FTK.delta_y_(1+index_nd_from_nM_(1+nM));
image_delta_x_true_M_(1+nM) = delta_x;
image_delta_y_true_M_(1+nM) = delta_y;
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z); %<-- note here we rotate S_k_p_wk_ by +gamma_z to form M_k_p_wk_. ;
T_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_wk_,-delta_x,-delta_y); %<-- note here we translate T_k_p_wk_ by -[delta_x,delta_y]. ;
M_k_p_wk_ = CTF_k_p_wk_.*T_k_p_wk_;
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_wkM__(:,1+nM) = M_k_q_wk_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% M_k_p_wkM__: %0.2fs',tmp_t)); end;
clear nCTF CTF_k_p_wk_ nS azimu_b polar_a index_gamma_z gamma_z S_k_p_wk_ T_k_p_wk_ M_k_p_wk_ M_k_q_wk_ ;
%%%%%%%%;

%%%%%%%%;
% Now we pick a single image and a single template, ;
% and test the calculation of the inner-product: ;
% <T(+delta)*M_k_p_wk_,CTF*R(-gamma_z)*S_k_p_wk_>. ;
% For this calculation we will set the translation delta_fix_ to be a fixed 2-vector, ;
% and let gamma_z vary over the on-grid rotations. ;
% TM_CTFRS_form_w_ stores the result of the analytical formula. ;
% TM_CTFRS_quad_w_ stores the result calculated using fourier-bessel qudarature. ;
% TM_CTFRS_brut_w_ stores the direct calculation (using quadrature on the polar-grid). ;
%%%%%%%%;
delta_fix_ = [+0.075;-0.035]; %<-- fixed delta_. ;
nS = max(0,min(n_S-1,round(n_S*2/3)));
nM = max(0,min(n_M-1,round(n_M*1/5)));
index_nS_from_nM = index_nS_from_nM_(1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTFM_k_p_wk_ = CTF_k_p_wk_.*M_k_p_wkM__(:,1+nM);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,CTFM_k_p_wk_,+delta_fix_(1+0),+delta_fix_(1+1));
TM_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,TM_k_p_wk_);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
TM_CTFRS_quad_w_ = ifft(sum(bsxfun(@times,reshape(conj(TM_k_q_wk_).*S_k_q_wk_,[n_w_max,n_k_p_r]),reshape(weight_2d_k_p_r_,[1,n_k_p_r])),2));
TM_CTFRS_brut_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
TM_CTFRS_brut_w_(1+nw) = sum(conj(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,TM_k_p_wk_,+gamma_z)).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'TM_CTFRS_brut_w_',TM_CTFRS_brut_w_,'TM_CTFRS_quad_w_',TM_CTFRS_quad_w_,' %%<-- should be zero');
%%%%%%%%;
azimu_b_S = viewing_azimu_b_S_(1+nS);
polar_a_S = viewing_polar_a_S_(1+nS);
gamma_z_S = 0.0;
R_S__ = Rz(-gamma_z_S)*Ry(-polar_a_S)*Rz(-azimu_b_S);
azimu_b_M = euler_azimu_b_true_M_(1+nM);
polar_a_M = euler_polar_a_true_M_(1+nM);
gamma_z_M = euler_gamma_z_true_M_(1+nM);
index_nd_from_nM = index_nd_from_nM_(1+nM);
delta_x_M = image_delta_x_true_M_(1+nM);
delta_y_M = image_delta_y_true_M_(1+nM);
R_M__ = Rz(+gamma_z_M)*Ry(-polar_a_M)*Rz(-azimu_b_M);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
CTF_phi = CTF_phi_C_(1+nCTF);
TM_CTFRS_form_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = (2*pi*nw)/n_w_max;
TM_CTFRS_form = 0;
for nsource_S=0:n_source-1;
delta_S_3_ = R_S__*delta_a_c_3s__(:,1+nsource_S);
delta_S_2_ = R2(-gamma_z) * delta_S_3_(1:2);
for nsource_M=0:n_source-1;
delta_M_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M);
delta_M_2_ = delta_M_3_(1:2) - [delta_fix_(1+0);delta_fix_(1+1)] + [delta_x_M;delta_y_M];
TM_CTFRS_I = I_xPPx_0(k_p_r_max,CTF_phi,delta_S_2_,CTF_phi,delta_M_2_);
TM_CTFRS_form = TM_CTFRS_form + TM_CTFRS_I;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
TM_CTFRS_form_w_(1+nw) = TM_CTFRS_form;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'TM_CTFRS_form_w_',TM_CTFRS_form_w_,'TM_CTFRS_brut_w_',TM_CTFRS_brut_w_,' %%<-- should be small');
fnorm_disp(flag_verbose,'TM_CTFRS_form_w_',TM_CTFRS_form_w_,'TM_CTFRS_quad_w_',TM_CTFRS_quad_w_,' %%<-- should be small');
%%%%%%%%;

%%%%%%%%;
% Now construct the template-norms: ;
% <R(+gamma_z)*CTF_k_p_wk_.*S_k_p_wk_,R(+gamma_z)*CTF_k_p_wk_.*S_k_p_wk_> ;
% or
% <CTF_k_p_wk_.*R(-gamma_z)*S_k_p_wk_,CTF_k_p_wk_.*R(-gamma_z)*S_k_p_wk_> ;
% Note that this does not involve collapsing onto principal-modes. ;
% CTF_S_l2_wSC_form___ stores the results from the analytical calculation. ;
% CTF_S_l2_wSC_quad___ stores the results using fourier-bessel quadrature. ;
% CTF_S_l2_wSC_brut___ stores the results using a direct-calculation. ;
%%%%%%%%;
CTF_S_l2_wSC_quad___ = zeros(n_w_max,n_S,n_CTF);
SS_k_p_wkS__ = conj(S_k_p_wkS__).*S_k_p_wkS__;
SS_k_q_wkS__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,SS_k_p_wkS__),[n_w_sum,n_S]);
CC_k_p_wkC__ = conj(CTF_k_p_wkC__).*CTF_k_p_wkC__;
CC_k_q_wkC__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,CC_k_p_wkC__),[n_w_sum,n_CTF]);
for nCTF=0:n_CTF-1;
CTF_S_l2_wSC_quad___(:,:,1+nCTF) = ifft(squeeze(sum(bsxfun(@times,reshape(bsxfun(@times,conj(CC_k_q_wkC__(:,1+nCTF)),SS_k_q_wkS__),[n_w_max,n_k_p_r,n_S]),reshape(weight_2d_k_p_r_,[1,n_k_p_r,1])),2)));
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% Now check one template. ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*1/5)));
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R_S__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
nCTF = max(0,min(n_CTF-1,round(n_CTF*1/3)));
tmp_CTF_phi = CTF_phi_C_(1+nCTF);
CTF_S_l2_w_quad_ = CTF_S_l2_wSC_quad___(:,1+nS,1+nCTF);
CTF_S_l2_w_brut_ = zeros(n_w_max,1);
CTF_S_l2_w_form_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
CS_k_p_wk_ = S_k_p_wkS__(:,1+nS).*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),+gamma_z);
CTF_S_l2_w_brut_(1+nw) = sum(conj(CS_k_p_wk_).*CS_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
tmp_S_delta_0_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource0); tmp_S_delta_0_2_ = tmp_S_delta_0_3_(1:2);
for nsource1=0:n_source-1;
tmp_S_delta_1_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource1); tmp_S_delta_1_2_ = tmp_S_delta_1_3_(1:2);
CTF_S_l2_w_form_(1+nw) = CTF_S_l2_w_form_(1+nw) + I_xPPx_0(k_p_r_max,tmp_CTF_phi+gamma_z,tmp_S_delta_0_2_,tmp_CTF_phi+gamma_z,tmp_S_delta_1_2_);
end;%for nsource0=0:n_source-1;
end;%for nsource1=0:n_source-1;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'CTF_S_l2_w_brut_',CTF_S_l2_w_brut_,'CTF_S_l2_w_quad_',CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'CTF_S_l2_w_form_',CTF_S_l2_w_form_,'CTF_S_l2_w_brut_',CTF_S_l2_w_brut_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'CTF_S_l2_w_form_',CTF_S_l2_w_form_,'CTF_S_l2_w_quad_',CTF_S_l2_w_quad_,' %%<-- should be <1e-6');
%%%%%%%%;

flag_pm = 1;
%%%%%%%%%%%%%%%%;
% Now calculate innerproduct Z_dwSM____. ;
% Here Z is the innerproduct: ;
% <T(+delta)*M_k_p_wk_,CTF*R(+gamma_z)*S_k_p_wk_>. ;
% Note the difference in rotation relative to the value calculated above, ;
% this time the rotation should be the same as the one used in CryoLike. 
% Here Z_dwSM_ampm____ is calculated using both: ;
% (1) a compression onto radial-principal-modes, and ;
% (2) a compression of the translation-operator. ;
% (i.e., we use a similar strategy to the driver from ampm). ;
%%%%%%%%%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions.; 
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
Z_dwSM_ampm____ = zeros(n_delta_v,n_w_max,n_S,n_M);
UX_M_l2_M_ = zeros(n_M,1);
UX_M_l2_dM__ = zeros(n_delta_v,n_M);
UX_knC___ = zeros(n_k_p_r,n_UX_rank,n_CTF);
X_weight_rC__ = zeros(n_k_p_r,n_CTF);
%%%%%%%%;
for nCTF=0:n_CTF-1;
if (flag_verbose); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
%%%%;
% Prepare principal-modes. ;
%%%%;
UX_kn__ = eye(n_k_p_r,n_UX_rank);
X_weight_r_ = sqrt(weight_2d_k_p_r_);
if  flag_pm;
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M_sub ...
,M_k_p_wkM__(:,1+index_M_sub_) ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
end;%if  flag_pm;
UX_knC___(:,:,1+nCTF) = UX_kn__;
X_weight_rC__(:,1+nCTF) = X_weight_r_;
%%%%;
% Prepare quasi-images. ;
%%%%;
tmp_t = tic();
CTF_M_sub_k_p_wkM__ = bsxfun(@times,CTF_k_p_wk_,M_k_p_wkM__(:,1+index_M_sub_));
CTF_M_sub_k_q_wkM__ = zeros(n_w_sum,n_M_sub);
for nM_sub=0:n_M_sub-1;
CTF_M_sub_k_q_wkM__(:,1+nM_sub) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,CTF_M_sub_k_p_wkM__(:,1+nM_sub));
end;%for nM_sub=0:n_M_sub-1;
M_sub_k_q_wkM__ = M_k_q_wkM__(:,1+index_M_sub_);
svd_VUXCTFM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,CTF_M_sub_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
svd_VUXM_sub_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M_sub,M_sub_k_q_wkM__,pm_n_UX_rank,UX_kn__,X_weight_r_);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tmpmh_VUXM_lwnM____3: %0.2fs',tmp_t)); end;
%%%%;
% Now calculate norms of the translated images. ;
%%%%;
tmp_t = tic();
UX_M_sub_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_sub_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% ampmh_UX_M_sub_l2_dm__1: %0.2fs',tmp_t)); end;
tmp_index_d0 = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); assert(numel(tmp_index_d0)==1); %<-- should be zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_M_sub_l2_dM__(1+tmp_index_d0,:),[n_M_sub,1]);
%%%%;
% Prepare UX_S_k_q_wnS__. ;
%%%%;
tmp_t = tic();
UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% UX_S_k_q_wnS__: %0.2fs',tmp_t)); end;
%%%%;
% Calculate Z_sub_dwSM____. ;
%%%%;
tmp_t = tic();
UX_S_k_q_nSw___ = permute(reshape(UX_S_k_q_wnS__,[n_w_max,n_UX_rank,n_S]),[2,3,1]);
tmp_t_sub = tic();
svd_VUXCTFM_sub_nMwl____ = permute(svd_VUXCTFM_sub_lwnM____,[3,4,2,1]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose); disp(sprintf(' %% svd_VUXCTFM_sub_nMwl____: %0.2fs',tmp_t_sub)); end;
tmp_t_sub = tic();
svd_SVUXCTFM_sub_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
svd_SVUXCTFM_sub_SMwl____(:,:,1+nw,1+nl) = ctranspose(UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXCTFM_sub_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_SMwl____: %0.6f',tmp_t_sub)); end;
tmp_t_sub = tic();
svd_SVUXCTFM_sub_lwSM____ = permute(ifft(permute(svd_SVUXCTFM_sub_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_sub_lwSM____: %0.6f',tmp_t_sub)); end;
tmp_t_sub = tic(); nop=0;
%svd_USESVUXCTFM_sub_dwSM____ = real(reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXCTFM_sub_lwSM____,[FTK.n_svd_l,n_w_max*n_S*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]));
svd_USESVUXCTFM_sub_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXCTFM_sub_lwSM____,[FTK.n_svd_l,n_w_max*n_S*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_USESVUXCTFM_sub_dwSM____: %0.6f',tmp_t_sub)); end;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% Z_sub_dwSM____: %0.2fs',tmp_t)); end;
%%%%;
% Store results. ;
%%%%;
UX_M_l2_M_(1+index_M_sub_) = UX_M_sub_l2_M_;
UX_M_l2_dM__(:,1+index_M_sub_) = UX_M_sub_l2_dM__;
Z_dwSM_ampm____(:,:,:,1+index_M_sub_) = svd_USESVUXCTFM_sub_dwSM____ ;
%%%%;
clear UX_kn__ X_weight_r_ CTF_k_p_wk_ index_M_sub_ UX_M_sub_l2_dM__ UX_M_sub_l2_M_ ;
clear UX_S_k_q_wnS__ ;
clear svd_VUXCTFM_sub_lwnM____ svd_VUXCTFM_sub_nMwl____ ;
clear svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_lwSM____ ;
clear svd_USESVUXCTFM_sub_dwSM____ ;
%%%%%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
clear nCTF;
%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now estimate landscape of innerproducts across delta_x_ and delta_y_. ;
% Limited to a single image-template pair and a fixed nw. ;
% Z_d_form_ <-- analytical-calculation. ;
% Z_d_ampm_ <-- ampm calculation. ;
% Z_d_brut_ <-- direct-calculation. ;
%%%%%%%%;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
nS = max(0,min(n_S-1,round(n_S*2/5)));
azimu_b_S = viewing_azimu_b_S_(1+nS);
polar_a_S = viewing_polar_a_S_(1+nS);
gamma_z_S = 0.0;
R_S__ = Rz(-gamma_z_S)*Ry(-polar_a_S)*Rz(-azimu_b_S);
nM = max(0,min(n_M-1,round(n_M*4/5)));
index_nS_from_nM = index_nS_from_nM_(1+nM);
azimu_b_M = euler_azimu_b_true_M_(1+nM);
polar_a_M = euler_polar_a_true_M_(1+nM);
gamma_z_M = euler_gamma_z_true_M_(1+nM);
index_nd_from_nM = index_nd_from_nM_(1+nM);
delta_x_M = image_delta_x_true_M_(1+nM);
delta_y_M = image_delta_y_true_M_(1+nM);
R_M__ = Rz(+gamma_z_M)*Ry(-polar_a_M)*Rz(-azimu_b_M);
nw = max(0,min(n_w_max-1,round(n_w_max*3/5)));
Z_d_ampm_ = Z_dwSM_ampm____(:,1+nw,1+nS,1+nM);
Z_d_brut_ = zeros(n_delta_v,1);
Z_d_form_ = zeros(n_delta_v,1);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_phi = CTF_phi_C_(1+nCTF);
%%%%;
gamma_z = (2*pi*nw)/n_w_max;
RS_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
CTF_RS_k_p_wk_ = RS_k_p_wk_.*CTF_k_p_wk_;
CTF_RS_k_p_l2 = sum(conj(CTF_RS_k_p_wk_).*CTF_RS_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
T0M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+0.0,+0.0);
T0M_k_p_l2_brut = sum(conj(T0M_k_p_wk_).*T0M_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
T0M_k_p_l2_form = 0;
for nsource_M0=0:n_source-1;
delta_M0_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M0);
delta_M0_2_ = delta_M0_3_(1:2) - [0*delta_x;0*delta_y] + [delta_x_M;delta_y_M];
for nsource_M1=0:n_source-1;
delta_M1_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M1);
delta_M1_2_ = delta_M1_3_(1:2) - [0*delta_x;0*delta_y] + [delta_x_M;delta_y_M];
T0M_k_p_I = I_xPPx_0(k_p_r_max,CTF_phi,delta_M0_2_,CTF_phi,delta_M1_2_);
T0M_k_p_l2_form = T0M_k_p_l2_form + T0M_k_p_I;
end;%for nsource_M1=0:n_source-1;
end;%for nsource_M0=0:n_source-1;
fnorm_disp(flag_verbose,'T0M_k_p_l2_form',T0M_k_p_l2_form,'T0M_k_p_l2_brut',T0M_k_p_l2_brut,' %%<-- should be <1e-6');
%%;
for ndelta_v=0:n_delta_v-1;
%%;
delta_x = FTK.delta_x_(1+ndelta_v); delta_y = FTK.delta_y_(1+ndelta_v);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
CTF_RS_k_p_TM_k_p = sum(conj(CTF_RS_k_p_wk_).*TM_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
Z_d_brut_(1+ndelta_v) = CTF_RS_k_p_TM_k_p;
%%;
%%%%;
Z_d_form = 0;
for nsource_S=0:n_source-1;
delta_S_3_ = R_S__*delta_a_c_3s__(:,1+nsource_S);
delta_S_2_ = R2(+gamma_z) * delta_S_3_(1:2);
for nsource_M=0:n_source-1;
delta_M_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M);
delta_M_2_ = delta_M_3_(1:2) - [delta_x;delta_y] + [delta_x_M;delta_y_M];
tmp_I = I_xPPx_0(k_p_r_max,CTF_phi,delta_S_2_,CTF_phi,delta_M_2_);
Z_d_form = Z_d_form + tmp_I;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
Z_d_form_(1+ndelta_v) = Z_d_form;
%%;
end;%for ndelta_v=0:n_delta_v-1;
%%;
fnorm_disp(flag_verbose,'Z_d_brut_',Z_d_brut_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_brut_',Z_d_brut_,' %%<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

return;





