function ...
[ ...
 parameter ...
] = ...
test_transforms_8_form_only( ...
 parameter ...
);

%%%%%%%%;
% Sets up a simple volume (on a spherical grid), ;
% then generates templates and calculates inner-products. ;
% These calculations consider anisotropic CTF. ;
%%%%%%%%;

str_thisfunction = 'test_transforms_8_form_only';

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
k_eq_d_double = 1.0; %<-- prefactor for k_eq_d, determines density of sampling in frequency-space (for radial-discretization in graphing only). ;
%%%%;
% volumetric parameters: will influence template and image. ;
%%%%;
k_p_r_max = 32.0; %<-- directly set the highest frequency (both for calculation and for plotting) ;
n_w_max = 512; %<-- directly set the number of inplane angles (both for calculation and for plotting). ;
delta_a_c_3s__ = [  ...
 +1.5 , -0.5 , +0.1 , -0.8 ...
;-0.5 , -1.5 , -0.3 , +0.4 ...
;+0.3 , +2.0 , +0.7 , +1.6 ...
] / 2 / k_p_r_max ; %<-- list of sources: delta_a_c_3s__(1+n3,1+ns) = n3-component of source-index ns. ;
%%%%;
% template euler-orientation. ;
%%%%;
template_polar_a = +2*pi/3 ;
template_azimu_b = +3*pi/5 ;
template_gamma_z = +0*pi/1 ; %<-- typically set to 0*pi when automatically generating templates. ;
%%%%;
% orientation of planar-function-prefactor for template. ;
%%%%;
CTFA_phi = -7*pi/11;
%%%%;
% image orientation and displacement. ;
%%%%;
image_polar_a = +1*pi/3 ;
image_azimu_b = +2*pi/5 ;
image_gamma_z = -3*pi/7 ; %<-- inplane rotation applied to image. Note that, for this image-construction, this has an orientation opposite to the template_gamma_z above. ;
image_delta_x = +0.15;
image_delta_y = -0.05;
%%%%;
% orientation of planar-function-prefactor for image. ;
%%%%;
CTFB_phi = +5*pi/11;
%%%%;
delta_x_external = -1.30/k_p_r_max; %<-- external delta_x. ;
delta_y_external = +2.20/k_p_r_max; %<-- external delta_y. ;
gamma_z_external = -2*pi/3; %<-- external gamma_z. ;
%%%%;
parameter = struct('type','parameter');
parameter.string_root = string_root;
parameter.string_user = string_user;
parameter.flag_verbose = flag_verbose;
parameter.flag_disp = flag_disp;
parameter.flag_replot = flag_replot;
parameter.k_eq_d_double = k_eq_d_double;
parameter.k_p_r_max = k_p_r_max;
parameter.n_w_max = n_w_max;
parameter.delta_a_c_3s__ = delta_a_c_3s__;
parameter.template_polar_a = template_polar_a;
parameter.template_azimu_b = template_azimu_b;
parameter.template_gamma_z = template_gamma_z;
parameter.CTFA_phi = CTFA_phi;
parameter.image_polar_a = image_polar_a;
parameter.image_azimu_b = image_azimu_b;
parameter.image_gamma_z = image_gamma_z;
parameter.image_delta_x = image_delta_x;
parameter.image_delta_y = image_delta_y;
parameter.CTFB_phi = CTFB_phi;
parameter.delta_x_external = delta_x_external;
parameter.delta_y_external = delta_y_external;
parameter.gamma_z_external = gamma_z_external;
test_transforms_8_form_only(parameter);
disp('returning'); return;
end;%if nargin<1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

na=0;
if (nargin<1+na); parameter=[]; end; na=na+1;

%%%%;
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
%%%%;
if ~isfield(parameter,'k_eq_d_double'); parameter.k_eq_d_double=1.0; end;
k_eq_d_double = parameter.k_eq_d_double;
if ~isfield(parameter,'k_p_r_max'); parameter.k_p_r_max=2*pi*32.0; end;
k_p_r_max = parameter.k_p_r_max;
if ~isfield(parameter,'n_w_max'); parameter.n_w_max=512; end;
n_w_max = parameter.n_w_max;
if ~isfield(parameter,'delta_a_c_3s__');
parameter.delta_a_c_3s__ = [  ...
 +1.5 , -0.5 , +0.1 , -0.8 ...
;-0.5 , -1.5 , -0.3 , +0.4 ...
;+0.3 , +2.0 , +0.7 , +1.6 ...
] / 2 / k_p_r_max ; %<-- list of sources: delta_a_c_3s__(1+n3,1+ns) = n3-component of source-index ns. ;
;
end;%if ~isfield(parameter,'delta_a_c_3s__');
delta_a_c_3s__ = parameter.delta_a_c_3s__;
if ~isfield(parameter,'template_polar_a'); parameter.template_polar_a=+2*pi/3; end;
template_polar_a = parameter.template_polar_a;
if ~isfield(parameter,'template_azimu_b'); parameter.template_azimu_b=+3*pi/5; end;
template_azimu_b = parameter.template_azimu_b;
if ~isfield(parameter,'template_gamma_z'); parameter.template_gamma_z=+0*pi/1; end;
template_gamma_z = parameter.template_gamma_z;
if ~isfield(parameter,'CTFA_phi'); parameter.CTFA_phi=-7*pi/11; end;
CTFA_phi = parameter.CTFA_phi;
if ~isfield(parameter,'image_polar_a'); parameter.image_polar_a=+1*pi/3; end;
image_polar_a = parameter.image_polar_a;
if ~isfield(parameter,'image_azimu_b'); parameter.image_azimu_b=+2*pi/5; end;
image_azimu_b = parameter.image_azimu_b;
if ~isfield(parameter,'image_gamma_z'); parameter.image_gamma_z=-3*pi/7; end;
image_gamma_z = parameter.image_gamma_z;
if ~isfield(parameter,'image_delta_x'); parameter.image_delta_x=+0.15; end;
image_delta_x = parameter.image_delta_x;
if ~isfield(parameter,'image_delta_y'); parameter.image_delta_y=-0.05; end;
image_delta_y = parameter.image_delta_y;
if ~isfield(parameter,'CTFB_phi'); parameter.CTFB_phi=+5*pi/11; end;
CTFB_phi = parameter.CTFB_phi;
if ~isfield(parameter,'delta_x_external'); parameter.delta_x_external=-1.30/k_p_r_max; end;
delta_x_external = parameter.delta_x_external;
if ~isfield(parameter,'delta_y_external'); parameter.delta_y_external=+2.20/k_p_r_max; end;
delta_y_external = parameter.delta_y_external;
if ~isfield(parameter,'gamma_z_external'); parameter.gamma_z_external=-2*pi/3; end;
gamma_z_external = parameter.gamma_z_external;
%%%%;

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
% Set up limited version of spherical-quadrature on sphere. ;
%%%%%%%%;
k_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_L ...
);
%%%%%%%%;
% Set up limited version of polar-quadrature on disk. ;
%%%%%%%%;
template_k_eq_d = -1;
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_p_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
,weight_3d_k_p_r_ ...
);
%%%%%%%%;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
%%%%%%%%;

%%%%%%%%;
% read number of sources. ;
%%%%%%%%;
n_source = size(delta_a_c_3s__,2);
if (flag_verbose>0); disp(sprintf(' %% k_p_r_max %0.6f n_k_p_r %d n_w_max %d n_w_sum %d; n_source %d',k_p_r_max,n_k_p_r,n_w_max,n_w_sum,n_source)); end;

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
% Now step through and reconstitute the template. ;
%%%%%%%%;
tmp_azimu_b = template_azimu_b;
tmp_polar_a = template_polar_a;
tmp_gamma_z = template_gamma_z;
tmp_R_S__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b); %<-- note -tmp_gamma_z. ;
S_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource);
S_k_p_wk_ = S_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
S_k_p_l2_brut = sum(conj(S_k_p_wk_).*S_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
S_k_p_l2_form = 0.0;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
tmp_delta_a_c_0_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource0);
tmp_delta_a_c_1_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(tmp_delta_a_c_0_(1:2) - tmp_delta_a_c_1_(1:2));
tmp_h2d = (2*pi)^2; if abs(tmp_kd)>1e-12; tmp_h2d = h2d_(tmp_kd); end;
S_k_p_l2_form = S_k_p_l2_form + tmp_h2d/(2*pi)^2 * (pi*k_p_r_max^2);
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
%%%%%%%%;
fnorm_disp(flag_verbose,'S_k_p_l2_form',S_k_p_l2_form,'S_k_p_l2_brut',S_k_p_l2_brut,' %<-- should be small');

%%%%%%%%;
% Now we generate an image to test the inner-product calculation. ;
% This image itself is a template (i.e. projection of the volume) multiplied by a planar-function (referred to as a CTF). ;
% We also shift the image by a fixed displacement. ;
%%%%%%%%;
CTFB_k_p_wk_ = zeros(n_w_sum,1);
CTFB_k_p_wk_(:,1) = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTFB_phi);
CTFA_k_p_wk_ = zeros(n_w_sum,1);
CTFA_k_p_wk_(:,1) = 2*k_p_r_wk_.*cos(k_p_w_wk_ - CTFA_phi);
%%%%%%%%;
% Now generate image. ;
%%%%%%%%;
tmp_azimu_b = image_azimu_b;
tmp_polar_a = image_polar_a;
tmp_gamma_z = 0.0; %<-- standard inplane-rotation for a template. ;
tmp_R_M__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R_M__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
%%%%%%%%;
T_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_,+image_gamma_z); %<-- note here we rotate T_k_p_wk_ by +gamma_z. ;
T_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_wk_,-image_delta_x,-image_delta_y); %<-- note here we translate T_k_p_wk_ by -[delta_x,delta_y]. ;
M_k_p_wk_ = CTFB_k_p_wk_.*T_k_p_wk_;
clear T_k_p_wk_ ;
%%%%%%%%;

%%%%%%%%;
% Now we test the calculation of the inner-product: ;
% <T(+delta_external_)*M_k_p_wk_,CTFA*R(-gamma_z_external)*S_k_p_wk_>. ;
% TM_CTFRS_form_w_ stores the result of the analytical formula. ;
% TM_CTFRS_brut_w_ stores the direct calculation (using quadrature on the polar-grid). ;
%%%%%%%%;
CTFM_k_p_wk_ = CTFA_k_p_wk_.*M_k_p_wk_;
TCTFM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,CTFM_k_p_wk_,+delta_x_external,+delta_y_external);
TM_CTFRS_brut = sum(conj(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,TCTFM_k_p_wk_,+gamma_z_external)).*S_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
%%%%%%%%;
azimu_b_S = template_azimu_b;
polar_a_S = template_polar_a;
gamma_z_S = template_gamma_z;
R_S__ = Rz(-gamma_z_S)*Ry(-polar_a_S)*Rz(-azimu_b_S);
azimu_b_M = image_azimu_b;
polar_a_M = image_polar_a;
gamma_z_M = image_gamma_z;
delta_x_M = image_delta_x;
delta_y_M = image_delta_y;
R_M__ = Rz(+gamma_z_M)*Ry(-polar_a_M)*Rz(-azimu_b_M); %<-- note sign switch for gamma_z_M. ;
TM_CTFRS_form = 0;
for nsource_S=0:n_source-1;
delta_S_3_ = R_S__*delta_a_c_3s__(:,1+nsource_S);
delta_S_2_ = R2(-gamma_z_external) * delta_S_3_(1:2);
for nsource_M=0:n_source-1;
delta_M_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M);
delta_M_2_ = delta_M_3_(1:2) - [delta_x_external;delta_y_external] + [delta_x_M;delta_y_M];
TM_CTFRS_I = I_xPPx_0(k_p_r_max,CTFA_phi,delta_S_2_,CTFB_phi,delta_M_2_);
TM_CTFRS_form = TM_CTFRS_form + TM_CTFRS_I;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
fnorm_disp(flag_verbose,'TM_CTFRS_form',TM_CTFRS_form,'TM_CTFRS_brut',TM_CTFRS_brut,' %%<-- should be small');
%%%%%%%%;

%%%%%%%%;
% Now construct the template-norms: ;
% <R(+gamma_z_external)*CTFA_k_p_wk_.*S_k_p_wk_,R(+gamma_z_external)*CTFA_k_p_wk_.*S_k_p_wk_> ;
% or
% <CTFA_k_p_wk_.*R(-gamma_z_external)*S_k_p_wk_,CTFA_k_p_wk_.*R(-gamma_z_external)*S_k_p_wk_> ;
% Note that this does not involve collapsing onto principal-modes. ;
% CTFA_S_l2_form stores the results from the analytical calculation. ;
% CTFA_S_l2_brut stores the results using a direct-calculation. ;
%%%%%%%%;
tmp_azimu_b = template_azimu_b;
tmp_polar_a = template_polar_a;
tmp_gamma_z = template_gamma_z;
tmp_R_S__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
CS_k_p_wk_ = S_k_p_wk_.*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTFA_k_p_wk_,+gamma_z_external);
CTFA_S_l2_brut = sum(conj(CS_k_p_wk_).*CS_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
CTFA_S_l2_form = 0;
for nsource0=0:n_source-1;
tmp_S_delta_0_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource0); tmp_S_delta_0_2_ = tmp_S_delta_0_3_(1:2);
for nsource1=0:n_source-1;
tmp_S_delta_1_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource1); tmp_S_delta_1_2_ = tmp_S_delta_1_3_(1:2);
CTFA_S_l2_form = CTFA_S_l2_form + I_xPPx_0(k_p_r_max,CTFA_phi+gamma_z_external,tmp_S_delta_0_2_,CTFA_phi+gamma_z_external,tmp_S_delta_1_2_);
end;%for nsource0=0:n_source-1;
end;%for nsource1=0:n_source-1;
fnorm_disp(flag_verbose,'CTFA_S_l2_form',CTFA_S_l2_form,'CTFA_S_l2_brut',CTFA_S_l2_brut,' %%<-- should be <1e-6');
%%%%%%%%;

%%%%%%%%%%%%%%%%;
% Now calculate innerproduct Z_dwSM____. ;
% Here Z is the innerproduct: ;
% <T(+delta_external_)*M_k_p_wk_,CTFA*R(+gamma_z_external)*S_k_p_wk_>. ;
% Note the difference in rotation relative to the value calculated above, ;
% this time the rotation should be the same as the one used in CryoLike. 
%%%%%%%%%%%%%%%%;
azimu_b_S = template_azimu_b;
polar_a_S = template_polar_a;
gamma_z_S = template_gamma_z;
R_S__ = Rz(-gamma_z_S)*Ry(-polar_a_S)*Rz(-azimu_b_S);
azimu_b_M = image_azimu_b;
polar_a_M = image_polar_a;
gamma_z_M = image_gamma_z;
delta_x_M = image_delta_x;
delta_y_M = image_delta_y;
R_M__ = Rz(+gamma_z_M)*Ry(-polar_a_M)*Rz(-azimu_b_M); %<-- note sign switch for gamma_z_M. ;
Z_d_brut = 0;
Z_d_form = 0;
%%%%;
RS_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z_external);
CTFA_RS_k_p_wk_ = RS_k_p_wk_.*CTFA_k_p_wk_;
CTFA_RS_k_p_l2 = sum(conj(CTFA_RS_k_p_wk_).*CTFA_RS_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
T0M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+0.0,+0.0);
T0M_k_p_l2_brut = sum(conj(T0M_k_p_wk_).*T0M_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
T0M_k_p_l2_form = 0;
for nsource_M0=0:n_source-1;
delta_M0_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M0);
delta_M0_2_ = delta_M0_3_(1:2) - [0*delta_x_external;0*delta_y_external] + [delta_x_M;delta_y_M];
for nsource_M1=0:n_source-1;
delta_M1_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M1);
delta_M1_2_ = delta_M1_3_(1:2) - [0*delta_x_external;0*delta_y_external] + [delta_x_M;delta_y_M];
T0M_k_p_I = I_xPPx_0(k_p_r_max,CTFB_phi,delta_M0_2_,CTFB_phi,delta_M1_2_);
T0M_k_p_l2_form = T0M_k_p_l2_form + T0M_k_p_I;
end;%for nsource_M1=0:n_source-1;
end;%for nsource_M0=0:n_source-1;
fnorm_disp(flag_verbose,'T0M_k_p_l2_form',T0M_k_p_l2_form,'T0M_k_p_l2_brut',T0M_k_p_l2_brut,' %%<-- should be <1e-6');
%%;
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x_external,+delta_y_external);
CTFA_RS_k_p_TM_k_p = sum(conj(CTFA_RS_k_p_wk_).*TM_k_p_wk_.*weight_2d_k_p_wk_)*(2*pi)^2;
Z_d_brut = CTFA_RS_k_p_TM_k_p;
%%;
%%%%;
Z_d_form = 0;
for nsource_S=0:n_source-1;
delta_S_3_ = R_S__*delta_a_c_3s__(:,1+nsource_S);
delta_S_2_ = R2(+gamma_z_external) * delta_S_3_(1:2);
for nsource_M=0:n_source-1;
delta_M_3_ = R_M__*delta_a_c_3s__(:,1+nsource_M);
delta_M_2_ = delta_M_3_(1:2) - [delta_x_external;delta_y_external] + [delta_x_M;delta_y_M];
tmp_I = I_xPPx_0(k_p_r_max,CTFA_phi,delta_S_2_,CTFB_phi,delta_M_2_);
Z_d_form = Z_d_form + tmp_I;
end;%for nsource_M=0:n_source-1;
end;%for nsource_S=0:n_source-1;
%%;
fnorm_disp(flag_verbose,'Z_d_form',Z_d_form,'Z_d_brut',Z_d_brut,' %%<-- should be <1e-2');
%%%%%%%%;

if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 5; np=0;
%%%%;
for ntype=0:5-1;
if ntype==0; tmp_wk_ = M_k_p_wk_; tmp_str = 'M_k_p_wk_'; end;
if ntype==1; tmp_wk_ = TM_k_p_wk_; tmp_str = 'TM_k_p_wk_'; end;
if ntype==2; tmp_wk_ = S_k_p_wk_; tmp_str = 'S_k_p_wk_'; end;
if ntype==3; tmp_wk_ = RS_k_p_wk_; tmp_str = 'RS_k_p_wk_'; end;
if ntype==4; tmp_wk_ = CTFA_RS_k_p_wk_; tmp_str = 'CTFA_RS_k_p_wk_'; end;
tmp_lim_ = max(abs(tmp_wk_))*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_wk_),tmp_lim_,colormap_80s);
axis image; axisnotick; title(sprintf('real(%s)',tmp_str),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(tmp_wk_),tmp_lim_,colormap_80s);
axis image; axisnotick; title(sprintf('imag(%s)',tmp_str),'Interpreter','none');
end;%for ntype=0:4-1;
%%%%;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

return;





