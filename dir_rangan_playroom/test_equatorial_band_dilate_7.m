%%%%%%%%;
% tests equatorial band dilation. ;
% generates figures as well. ;
%%%%%%%%;

clear;
%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_manuscript = sprintf('/%s/rangan/dir_cryoem/dir_spurious_heterogeneity_manuscript',string_root);
if ~exist(dir_manuscript,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript)); mkdir(dir_manuscript); end;
dir_manuscript_jpg = sprintf('%s/dir_M3d_shape_longitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;
flag_verbose = 0;
flag_replot=1;

%%%%%%%%;
% First set up collection of templates. ;
%%%%%%%%;
verbose = 0; flag_disp=1; nf=0;
n_k_p_r = 1; k_p_r_max = 1.0; k_p_r_ = k_p_r_max; weight_3d_k_p_r_ = 1; l_max_ = 0; a_k_Y_ = 0;
viewing_k_eq_d = 0.5*k_p_r_max/(2*pi);;
template_k_eq_d = -1;
n_w_max = 128; n_w_0in_ = n_w_max;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max)); cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
[ ...
 S_k_p__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_w_ ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0_wS__ ...
,template_k_c_1_wS__ ...
,template_k_c_2_wS__ ...
] = ...
get_template_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_S = n_viewing_S;
template_k_r01_wS__ = sqrt(template_k_c_0_wS__.^2 + template_k_c_1_wS__.^2);
template_k_r012_wS__ = sqrt(template_k_c_0_wS__.^2 + template_k_c_1_wS__.^2 + template_k_c_2_wS__.^2);
template_azimu_b_wS__ = atan2(template_k_c_1_wS__,template_k_c_0_wS__);
template_polar_a_wS__ = atan2(template_k_r01_wS__,template_k_c_2_wS__);
viewing_k_c_0_S_ = cos(viewing_azimu_b_S_).*sin(viewing_polar_a_S_);
viewing_k_c_1_S_ = sin(viewing_azimu_b_S_).*sin(viewing_polar_a_S_);
viewing_k_c_2_S_ = cos(viewing_polar_a_S_);

n_test = 16;
for ntest=0:n_test-1;
nS = max(0,min(n_S-1,floor(n_S*ntest/n_test)));
[ ...
 template_k_c_0_w_ ...
,template_k_c_1_w_ ...
,template_k_c_2_w_ ...
,template_azimu_b_w_ ...
,template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,viewing_azimu_b_S_(1+nS) ...
,viewing_polar_a_S_(1+nS) ...
,n_w_max ...
);
tmp_error = ...
 + fnorm(template_k_c_0_wS__(:,1+nS)-template_k_c_0_w_) ...
 + fnorm(template_k_c_1_wS__(:,1+nS)-template_k_c_1_w_) ...
 + fnorm(template_k_c_2_wS__(:,1+nS)-template_k_c_2_w_) ...
 + fnorm(template_azimu_b_wS__(:,1+nS)-template_azimu_b_w_) ...
 + fnorm(template_polar_a_wS__(:,1+nS)-template_polar_a_w_) ...
;
if (verbose>0); disp(sprintf(' %% ntest %d/%d tmp_error %0.16f;',ntest,n_test,tmp_error)); end;
end;%for ntest=0:n_test-1;

Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;

Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];

%%%%%%%%;
% Set parameters for equatorial band dilation. ;
% We assume that each template is rotated (along the equator) by g_dilation, ;
% corresponding to an overall mapping of viewing_azimu_b via f_dilation. ;
% Note that, due to our strategy for constructing templates, ;
% the equatorial band dilation corresponds to incrementing the viewing_azimu_b by g_dilation, ;
% meaning that the template itself will remain indexed by gamma_z.
% This (canonical) indexing via gamma_z will typically result in a shift between the ;
% original gamma_z associated with a particular point (in a particular template) ;
% and the final gamma_z associated with that same point (in the predilated version of that template). ;
%%%%%%%%;
equa_band_dilated_amplitude = 0.015;
g_dilation = @(point_pole_predilated_azimu_b) equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b); %<-- approximation holds well for first nontrivial mode. ;
f_dilation = @(point_pole_predilated_azimu_b) point_pole_predilated_azimu_b + g_dilation(point_pole_predilated_azimu_b);
%%%%%%%%;
% Note that, to first order: ;
% f_dilation_inverse(azimu_b) = azimu_b - g_dilation(azimu_b). ;
%%%%%%%%;
% point_output is indexed by polar_a and azimu_b. ;
% point_pole is indexed by n_w_max. ;
%%%%%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% plot diagram of equatorial band dilation. ;
% Caption: ;
% In this figure we illustrate some properties of the equatorial-band-dilation described in the text. ;
% On the left we show a sphere (grey), representing a single shell within k-space. ;
% The cyan-borded square (with magenta interior) indicates the viewing-angle $\tau$ (in terms of $\alpha$ and $\beta$) corresponding to a particular template. ;
% the slice in k-space associated with the the template $S(\tau)$ corresponds to a great-circle on the spherical-shell, indicated in green. ;
% One point $\vk$ within the template (on the surface of the spherical shell) is indicated with a black-bordered dot (with grey interior). ;
% The yellow-bordered square (again with a magenta interior) indicates the viewing-angle $\tau^{\prime} = f^{-1}(\tau)$, with $f$ defined using an amplitude of $0.15$. ;
% That is to say, the yellow-bordered square indicates the viewing-angle $\tau^{\prime}$ which, after dilation by $f$, produces the original viewing-angle $\tau$ shown in the cyan-bordered square. ;
% We use a magenta dot to indicate the point $\vk^{\prime}$ on the template $S(\tau^{\prime})$ which is mapped to $\vk$ after $S(\tau^{\prime})$ is dilated by $f$.
% On the right-subplot we show not only the cyan-bordered square from the left-subplot (still with a magenta interior), but also several other cyan-bordered squares. ;
% Each of these other cyan-bordered squares indicates the viewing-angle $\tau$ of a template containing the point $\vk$. ;
% These cyan-bordered squares form a great-circle normal to the point $\vk$; the interior color of these cyan-bordered squares varies periodically around this great-circle. ;
% For each of these cyan-bordered squares, we also show a corresponding yellow-bordered square (with matching interior color). ;
% Just as in the right-subplot, each yellow-bordered square indicates the viewing-angle $\tau^{\prime}$ obtained by applying the inverse-dilation to the viewing-angle $\tau$ of the corresponding cyan-bordered square. ;
% For each yellow-bordered square we also use a dot to indicate the point $\vk^{\prime}$ which gets mapped to $\vk$ after dilation; the color of each dot matches the color used for the square. ;
% Note that the colored dots form an orbit supported along the latitudinal line passing through the point $\vk$ itself (i.e., varying along the azimuthal-direction).
% This colored latitudinal-segment in the right-subplot indicates the locations $\{\vk^{\prime}\}$ which contribute to the value of the volume at $\vk$ after dilation of all the templates.
%%%%%%%%;
tmp_equa_band_dilated_amplitude = 0.15;
tmp_g_dilation = @(point_pole_predilated_azimu_b) tmp_equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b);
tmp_f_dilation = @(point_pole_predilated_azimu_b) point_pole_predilated_azimu_b + tmp_g_dilation(point_pole_predilated_azimu_b);
figure(1+nf);nf=nf+1;clf;figmed;
c_hsv__ = colormap('hsv'); n_c_hsv = size(c_hsv__,1);
sym_point = 'o'; sym_pole = 's';
markersize_use = 8; markersize_dot = 4; linewidth_use = 1.5; linewidth_mrk = 2.0;
for np=0:1; subplot(1,2,1+np); plot_sphere_grid_0; hold on; end;
npole_use = floor(0.92*n_w_max);
n_point_a = 1+128; point_output_polar_a_a_ = linspace(0,pi,n_point_a+2); point_output_polar_a_a_ = transpose(point_output_polar_a_a_(2:end-1));
%n_point_a = 1+8; point_output_polar_a_a_ = linspace(0,pi,n_point_a+2); point_output_polar_a_a_ = transpose(point_output_polar_a_a_(2:end-1));
n_point_b = 1+32; point_output_azimu_b_b_ = linspace(0,2*pi,n_point_b+0); point_output_azimu_b_b_ = transpose(point_output_azimu_b_b_(1:end));
%n_point_b = 0+3; point_output_azimu_b_b_ = linspace(0,2*pi,n_point_b+1); point_output_azimu_b_b_ = transpose(point_output_azimu_b_b_(1:n_point_b));
point_output_azimu_b_ab__ = zeros(n_point_a,n_point_b,1);
point_output_polar_a_ab__ = zeros(n_point_a,n_point_b,1);
point_pole_predilated_template_gammax_k_c_0_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_k_c_1_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_k_c_2_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_azimu_b_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_polar_a_abw___ = zeros(n_point_a,n_point_b,n_w_max);
%%%%;
point_pole_azimu_b_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_polar_a_abw___ = zeros(n_point_a,n_point_b,n_w_max);
for npoint_a=0:n_point_a-1;for npoint_b=0:n_point_b-1;
point_output_azimu_b = point_output_azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + point_output_polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_output_azimu_b_ab__(1+npoint_a,1+npoint_b) = point_output_azimu_b;
point_output_polar_a_ab__(1+npoint_a,1+npoint_b) = point_output_polar_a;
for npole=0:n_w_max-1;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
point_pole_azimu_b_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_azimu_b;
point_pole_polar_a_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_polar_a;
end;%for npole=0:n_w_max-1;
end;end;%for npoint_a=0:n_point_a-1;for npoint_b=0:n_point_b-1;
%%%%;
npoint_a = 32; %<-- not too close to the equator. ;
npoint_b = 32;
%%%%%%%%;
% Note that: ;
% point_output_k_c_ = ;
% [+cos(point_output_azimu_b)*cos(point_output_polar_a);+sin(point_output_azimu_b)*cos(point_output_polar_a);-sin(point_output_polar_a)] ;
% which is perpendicular to: ;
% [+cos(point_output_azimu_b)*sin(point_output_polar_a);+sin(point_output_azimu_b)*sin(point_output_polar_a);+cos(point_output_polar_a)] ;
% the latter of which is the first point listed in point_pole_k_c__(:,1+0). ;
%%%%%%%%;
point_output_azimu_b = point_output_azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + point_output_polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_output_azimu_b_ab__(1+npoint_a,1+npoint_b) = point_output_azimu_b;
point_output_polar_a_ab__(1+npoint_a,1+npoint_b) = point_output_polar_a;
%%%%%%%%;
for np=0:1; subplot(1,2,1+np); plot3(point_output_k_c_(1+0),point_output_k_c_(1+1),point_output_k_c_(1+2),sym_point,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0.85*[1,1,1],'MarkerSize',markersize_use,'LineWidth',linewidth_mrk); end;
%%%%%%%%%%%%%%%%;
for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
% Now we step through each of the templates associated with the point point_output. ;
%%%%%%%%;
nc_hsv = max(0,min(n_c_hsv-1,floor(n_c_hsv*(periodize(npole,0,n_w_max/2)/(n_w_max/2))))); %<-- note double winding (i.e., andipodal templates are the same). ;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
[ ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,point_pole_template_azimu_b_w_ ...
,point_pole_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_azimu_b ...
,point_pole_polar_a ...
,n_w_max ...
);
%%%%%%%%;
% Here we determine the predilated template that is associated with the original template mapped to point_output. ;
%%%%%%%%;
tmp_f_error = @(point_pole_predilated_azimu_b) abs(tmp_f_dilation(point_pole_predilated_azimu_b) - point_pole_azimu_b).^2;
point_pole_predilated_azimu_b = fminsearch(tmp_f_error,point_pole_azimu_b);
tmp_point_pole_azimu_b = point_pole_predilated_azimu_b + equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b);
if (verbose>0); disp(sprintf(' %% npole %d/%d, fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b): %0.16f',npole,n_w_max,fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b))); end;
point_pole_predilated_polar_a = point_pole_polar_a;
point_pole_predilated_k_c_ = [ ...
  cos(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; sin(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; cos(point_pole_predilated_polar_a) ...
];
[ ...
 point_pole_predilated_template_k_c_0_w_ ...
,point_pole_predilated_template_k_c_1_w_ ...
,point_pole_predilated_template_k_c_2_w_ ...
,point_pole_predilated_template_azimu_b_w_ ...
,point_pole_predilated_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_predilated_azimu_b ...
,point_pole_predilated_polar_a ...
,n_w_max ...
);
%%%%%%%%;
point_pole_template_gamma0_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;-sin(point_pole_polar_a) ...
];
point_pole_template_sgx_k_c_ = cross(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_template_sgx = dot(point_pole_template_sgx_k_c_,point_pole_k_c_);
point_pole_template_cgx = dot(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_gx = atan2(point_pole_template_sgx,point_pole_template_cgx);
point_pole_template_gammax_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_polar_a)*cos(point_pole_gx) ...
];
% fnorm(point_pole_template_gammax_k_c_ - point_output_k_c_),; %<-- should be 0. ;
point_pole_template_gammax_k_r01 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2);
% fnorm((cos(point_pole_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2 - point_pole_template_gammax_k_r01^2),; %<-- should be 0. ;
point_pole_template_gammax_k_r012 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2 + point_pole_template_gammax_k_c_(1+2).^2);
point_pole_template_gammax_azimu_b = atan2(point_pole_template_gammax_k_c_(1+1),point_pole_template_gammax_k_c_(1+0));
point_pole_template_gammax_polar_a = atan2(point_pole_template_gammax_k_r01,point_pole_template_gammax_k_c_(1+2));
%%%%%%%%;
% Now we determine the gamma_z (denoted gammax) within predilated template that maps to point_output. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_ = [ ...
 +cos(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_predilated_polar_a)*cos(point_pole_gx) ...
];
point_pole_predilated_template_gammax_k_c_0_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+0);
point_pole_predilated_template_gammax_k_c_1_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+1);
point_pole_predilated_template_gammax_k_c_2_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+2);
point_pole_predilated_template_gammax_k_r01 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2);
% fnorm((cos(point_pole_predilated_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2 - point_pole_predilated_template_gammax_k_r01^2),; %<-- should be 0. ;
point_pole_predilated_template_gammax_k_r012 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2 + point_pole_predilated_template_gammax_k_c_(1+2).^2);
point_pole_predilated_template_gammax_azimu_b = atan2(point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+0));
point_pole_predilated_template_gammax_polar_a = atan2(point_pole_predilated_template_gammax_k_r01,point_pole_predilated_template_gammax_k_c_(1+2));
point_pole_predilated_template_gammax_azimu_b_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_azimu_b;
point_pole_predilated_template_gammax_polar_a_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_polar_a;
%{
tmp_da = (point_pole_predilated_azimu_b - point_pole_azimu_b);
tmp_x0 = +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx);
tmp_y0 = +sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx);
tmp_x1 = +cos(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_predilated_azimu_b)*sin(point_pole_gx);
tmp_y1 = +sin(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_predilated_azimu_b)*sin(point_pole_gx);
tmp_dx = -sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - cos(point_pole_azimu_b)*sin(point_pole_gx);
tmp_dy = +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx);
disp(sprintf(' %% fnorm(tmp_dx - (tmp_x1-tmp_x0)/tmp_da) %0.16f',fnorm(tmp_dx - (tmp_x1-tmp_x0)/tmp_da)));
disp(sprintf(' %% fnorm(tmp_dy - (tmp_y1-tmp_y0)/tmp_da) %0.16f',fnorm(tmp_dy - (tmp_y1-tmp_y0)/tmp_da)));
tmp_B0 = point_pole_template_gammax_azimu_b;
tmp_B1 = point_pole_predilated_template_gammax_azimu_b;
tmp_dB = (tmp_dy*tmp_x0 - tmp_y0*tmp_dx)/max(1e-12,(cos(point_pole_predilated_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2);
disp(sprintf(' %% fnorm(tmp_dB - (tmp_B1-tmp_B0)/tmp_da) %0.16f',fnorm(tmp_dB - (tmp_B1-tmp_B0)/tmp_da)));
 %}
% fnorm((point_pole_predilated_template_gammax_azimu_b - point_pole_template_gammax_azimu_b) - (-equa_band_dilated_amplitude*sin(2*point_pole_azimu_b))),; %<-- should be order equa_band_dilated_amplitude.^2. ;
%%%%%%%%;
for np=0:1;
subplot(1,2,1+np);
if (np==0 & npole==npole_use);
plot3(point_pole_k_c_(1+0),point_pole_k_c_(1+1),point_pole_k_c_(1+2),sym_pole,'MarkerEdgeColor','c','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use,'LineWidth',linewidth_mrk);
plot3(point_pole_predilated_k_c_(1+0),point_pole_predilated_k_c_(1+1),point_pole_predilated_k_c_(1+2),sym_pole,'MarkerEdgeColor','y','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use,'LineWidth',linewidth_mrk);
plot3(point_pole_template_k_c_0_w_,point_pole_template_k_c_1_w_,point_pole_template_k_c_2_w_,'g-','LineWidth',linewidth_use);
plot3(point_pole_predilated_template_gammax_k_c_(1+0),point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+2),sym_point,'MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerEdgeColor','none','MarkerSize',markersize_use,'LineWidth',linewidth_mrk);
end;%if (np==0 & npole==npole_use);
if (np==1);
plot3(point_pole_k_c_(1+0),point_pole_k_c_(1+1),point_pole_k_c_(1+2),sym_pole,'MarkerEdgeColor','c','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use,'LineWidth',linewidth_mrk);
plot3(point_pole_predilated_k_c_(1+0),point_pole_predilated_k_c_(1+1),point_pole_predilated_k_c_(1+2),sym_pole,'MarkerEdgeColor','y','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use,'LineWidth',linewidth_mrk);
plot3(point_pole_predilated_template_gammax_k_c_(1+0),point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+2),sym_point,'MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use,'MarkerEdgeColor','none');
end;%if (np==1);
end;%for np=0:1;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
for np=0:1;
subplot(1,2,1+np);
hold off;
view(25,30);
axis vis3d;
xlim([-1,+1]);set(gca,'XTick',[-1,+1],'XTickLabel',{'-','+'}); xlabel('x');
ylim([-1,+1]);set(gca,'YTick',[-1,+1],'YTickLabel',{'-','+'}); ylabel('y');
zlim([-1,+1]);set(gca,'ZTick',[-1,+1],'ZTickLabel',{'-','+'}); zlabel('z');
xlabel('');set(gca,'XTick',[]);
ylabel('');set(gca,'YTick',[]);
zlabel('');set(gca,'ZTick',[]);
end;%for np=0:1;
%%%%%%%%;
drawnow;
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_diagram_FIGA',dir_manuscript_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
sgtitle('','Interpreter','none');
print('-djpeg',sprintf('%s_strip.jpg',fname_fig_pre));
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
end;%if flag_disp;

error('stopping early');

flag_disp=0;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;fig80s; c_hsv__ = colormap('hsv'); n_c_hsv = size(c_hsv__,1);
markersize_big = 8; markersize_use = 6; markersize_dot = 4; linewidth_use = 0.5;
plot_sphere_grid_0;
hold on;
end;%if flag_disp;
%%%%%%%%;
n_point_a = 1+128; point_output_polar_a_a_ = linspace(0,pi,n_point_a+2); point_output_polar_a_a_ = transpose(point_output_polar_a_a_(2:end-1));
%n_point_a = 1+8; point_output_polar_a_a_ = linspace(0,pi,n_point_a+2); point_output_polar_a_a_ = transpose(point_output_polar_a_a_(2:end-1));
n_point_b = 1+32; point_output_azimu_b_b_ = linspace(0,2*pi,n_point_b+0); point_output_azimu_b_b_ = transpose(point_output_azimu_b_b_(1:end));
%n_point_b = 0+3; point_output_azimu_b_b_ = linspace(0,2*pi,n_point_b+1); point_output_azimu_b_b_ = transpose(point_output_azimu_b_b_(1:n_point_b));
point_output_azimu_b_ab__ = zeros(n_point_a,n_point_b,1);
point_output_polar_a_ab__ = zeros(n_point_a,n_point_b,1);
point_pole_predilated_template_gammax_k_c_0_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_k_c_1_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_k_c_2_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_azimu_b_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_predilated_template_gammax_polar_a_abw___ = zeros(n_point_a,n_point_b,n_w_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npoint_a=0:n_point_a-1; 
if (verbose>-1); disp(sprintf(' %% npoint_a %d/%d',npoint_a,n_point_a)); end;
for npoint_b=0:n_point_b-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Note that: ;
% point_output_k_c_ = ;
% [+cos(point_output_azimu_b)*cos(point_output_polar_a);+sin(point_output_azimu_b)*cos(point_output_polar_a);-sin(point_output_polar_a)] ;
% which is perpendicular to: ;
% [+cos(point_output_azimu_b)*sin(point_output_polar_a);+sin(point_output_azimu_b)*sin(point_output_polar_a);+cos(point_output_polar_a)] ;
% the latter of which is the first point listed in point_pole_k_c__(:,1+0). ;
%%%%%%%%;
point_output_azimu_b = point_output_azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + point_output_polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_output_azimu_b_ab__(1+npoint_a,1+npoint_b) = point_output_azimu_b;
point_output_polar_a_ab__(1+npoint_a,1+npoint_b) = point_output_polar_a;
%%%%%%%%;
if flag_disp;
plot3(point_output_k_c_(1+0),point_output_k_c_(1+1),point_output_k_c_(1+2),'mo','MarkerFaceColor',0.85*[1,1,1],'MarkerSize',markersize_use);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
% Now we step through each of the templates associated with the point point_output. ;
%%%%%%%%;
if flag_disp;
nc_hsv = max(0,min(n_c_hsv-1,floor(n_c_hsv*(periodize(npole,0,n_w_max/2)/(n_w_max/2))))); %<-- note double winding (i.e., andipodal templates are the same). ;
end;%if flag_disp;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
[ ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,point_pole_template_azimu_b_w_ ...
,point_pole_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_azimu_b ...
,point_pole_polar_a ...
,n_w_max ...
);
%%%%%%%%;
% Here we determine the predilated template that is associated with the original template mapped to point_output. ;
%%%%%%%%;
tmp_f_error = @(point_pole_predilated_azimu_b) abs(f_dilation(point_pole_predilated_azimu_b) - point_pole_azimu_b).^2;
point_pole_predilated_azimu_b = fminsearch(tmp_f_error,point_pole_azimu_b);
tmp_point_pole_azimu_b = point_pole_predilated_azimu_b + equa_band_dilated_amplitude*sin(2*point_pole_predilated_azimu_b);
if (verbose>0); disp(sprintf(' %% npole %d/%d, fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b): %0.16f',npole,n_w_max,fnorm(point_pole_azimu_b-tmp_point_pole_azimu_b))); end;
point_pole_predilated_polar_a = point_pole_polar_a;
point_pole_predilated_k_c_ = [ ...
  cos(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; sin(point_pole_predilated_azimu_b)*sin(point_pole_predilated_polar_a) ...
; cos(point_pole_predilated_polar_a) ...
];
[ ...
 point_pole_predilated_template_k_c_0_w_ ...
,point_pole_predilated_template_k_c_1_w_ ...
,point_pole_predilated_template_k_c_2_w_ ...
,point_pole_predilated_template_azimu_b_w_ ...
,point_pole_predilated_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 verbose ...
,1 ...
,point_pole_predilated_azimu_b ...
,point_pole_predilated_polar_a ...
,n_w_max ...
);
%%%%%%%%;
point_pole_template_gamma0_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;-sin(point_pole_polar_a) ...
];
point_pole_template_sgx_k_c_ = cross(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_template_sgx = dot(point_pole_template_sgx_k_c_,point_pole_k_c_);
point_pole_template_cgx = dot(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_gx = atan2(point_pole_template_sgx,point_pole_template_cgx);
point_pole_template_gammax_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_polar_a)*cos(point_pole_gx) ...
];
% fnorm(point_pole_template_gammax_k_c_ - point_output_k_c_),; %<-- should be 0. ;
point_pole_template_gammax_k_r01 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2);
% fnorm((cos(point_pole_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2 - point_pole_template_gammax_k_r01^2),; %<-- should be 0. ;
point_pole_template_gammax_k_r012 = sqrt(point_pole_template_gammax_k_c_(1+0).^2 + point_pole_template_gammax_k_c_(1+1).^2 + point_pole_template_gammax_k_c_(1+2).^2);
point_pole_template_gammax_azimu_b = atan2(point_pole_template_gammax_k_c_(1+1),point_pole_template_gammax_k_c_(1+0));
point_pole_template_gammax_polar_a = atan2(point_pole_template_gammax_k_r01,point_pole_template_gammax_k_c_(1+2));
%%%%%%%%;
% Now we determine the gamma_z (denoted gammax) within predilated template that maps to point_output. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_ = [ ...
 +cos(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_predilated_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_predilated_polar_a)*cos(point_pole_gx) ...
];
point_pole_predilated_template_gammax_k_c_0_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+0);
point_pole_predilated_template_gammax_k_c_1_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+1);
point_pole_predilated_template_gammax_k_c_2_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_k_c_(1+2);
point_pole_predilated_template_gammax_k_r01 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2);
% fnorm((cos(point_pole_predilated_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2 - point_pole_predilated_template_gammax_k_r01^2),; %<-- should be 0. ;
point_pole_predilated_template_gammax_k_r012 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2 + point_pole_predilated_template_gammax_k_c_(1+2).^2);
point_pole_predilated_template_gammax_azimu_b = atan2(point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+0));
point_pole_predilated_template_gammax_polar_a = atan2(point_pole_predilated_template_gammax_k_r01,point_pole_predilated_template_gammax_k_c_(1+2));
point_pole_predilated_template_gammax_azimu_b_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_azimu_b;
point_pole_predilated_template_gammax_polar_a_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_predilated_template_gammax_polar_a;
%{
tmp_da = (point_pole_predilated_azimu_b - point_pole_azimu_b);
tmp_x0 = +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx);
tmp_y0 = +sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx);
tmp_x1 = +cos(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_predilated_azimu_b)*sin(point_pole_gx);
tmp_y1 = +sin(point_pole_predilated_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_predilated_azimu_b)*sin(point_pole_gx);
tmp_dx = -sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - cos(point_pole_azimu_b)*sin(point_pole_gx);
tmp_dy = +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx);
disp(sprintf(' %% fnorm(tmp_dx - (tmp_x1-tmp_x0)/tmp_da) %0.16f',fnorm(tmp_dx - (tmp_x1-tmp_x0)/tmp_da)));
disp(sprintf(' %% fnorm(tmp_dy - (tmp_y1-tmp_y0)/tmp_da) %0.16f',fnorm(tmp_dy - (tmp_y1-tmp_y0)/tmp_da)));
tmp_B0 = point_pole_template_gammax_azimu_b;
tmp_B1 = point_pole_predilated_template_gammax_azimu_b;
tmp_dB = (tmp_dy*tmp_x0 - tmp_y0*tmp_dx)/max(1e-12,(cos(point_pole_predilated_polar_a)*cos(point_pole_gx))^2 + (sin(point_pole_gx))^2);
disp(sprintf(' %% fnorm(tmp_dB - (tmp_B1-tmp_B0)/tmp_da) %0.16f',fnorm(tmp_dB - (tmp_B1-tmp_B0)/tmp_da)));
 %}
% fnorm((point_pole_predilated_template_gammax_azimu_b - point_pole_template_gammax_azimu_b) - (-equa_band_dilated_amplitude*sin(2*point_pole_azimu_b))),; %<-- should be order equa_band_dilated_amplitude.^2. ;
%%%%%%%%;
if flag_disp;
if flag_disp> 1;
plot3(point_pole_k_c__(1+0,:),point_pole_k_c__(1+1,:),point_pole_k_c__(1+2,:),'k.','MarkerSize',markersize_dot);
plot3(point_pole_k_c_(1+0),point_pole_k_c_(1+1),point_pole_k_c_(1+2),'co','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use);
plot3(point_pole_predilated_k_c_(1+0),point_pole_predilated_k_c_(1+1),point_pole_predilated_k_c_(1+2),'yo','MarkerFaceColor',c_hsv__(1+nc_hsv,:).^0.5,'MarkerSize',markersize_use);
plot3(point_pole_template_k_c_0_w_,point_pole_template_k_c_1_w_,point_pole_template_k_c_2_w_,'g-','LineWidth',linewidth_use);
plot3(point_pole_template_k_c_0_w_(1+0),point_pole_template_k_c_1_w_(1+0),point_pole_template_k_c_2_w_(1+0),'ko','MarkerSize',markersize_use);
plot3(point_pole_template_gamma0_k_c_(1+0),point_pole_template_gamma0_k_c_(1+1),point_pole_template_gamma0_k_c_(1+2),'kx','MarkerSize',markersize_use);
end;%if flag_disp> 1;
plot3(point_pole_template_gammax_k_c_(1+0),point_pole_template_gammax_k_c_(1+1),point_pole_template_gammax_k_c_(1+2),'ko','MarkerSize',markersize_big);
plot3(point_pole_predilated_template_gammax_k_c_(1+0),point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+2),'.','Color',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_dot);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;end;%for npoint_a=0:n_point_a-1; for npoint_b=0:n_point_b-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
xlim([-1,+1]); ylim([-1,+1]); zlim([-1,+1]); axis equal; axis vis3d;
xlabel('x'); ylabel('y'); zlabel('z');
hold off;
end;%if flag_disp;
%%%%%%%%;

disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_0_abw___ - cos(point_pole_predilated_template_gammax_azimu_b_abw___).*sin(point_pole_predilated_template_gammax_polar_a_abw___))));
disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_1_abw___ - sin(point_pole_predilated_template_gammax_azimu_b_abw___).*sin(point_pole_predilated_template_gammax_polar_a_abw___))));
disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_2_abw___ - cos(point_pole_predilated_template_gammax_polar_a_abw___))));
%%%%%%%%;
% Calculate avg, std and var. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_0_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_0_abw___,3);
point_pole_predilated_template_gammax_k_c_1_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_1_abw___,3);
point_pole_predilated_template_gammax_k_c_2_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_2_abw___,3);
point_pole_predilated_template_gammax_polar_a_avg_ab__ = mean(point_pole_predilated_template_gammax_polar_a_abw___,3);
point_pole_predilated_template_gammax_polar_a_std_ab__ = std(point_pole_predilated_template_gammax_polar_a_abw___,1,3);
point_pole_predilated_template_gammax_polar_a_var_ab__ = var(point_pole_predilated_template_gammax_polar_a_abw___,1,3);
point_pole_predilated_template_gammax_azimu_b_avg_ab__ = mean(point_pole_predilated_template_gammax_azimu_b_abw___,3);
point_pole_predilated_template_gammax_azimu_b_std_ab__ = std(point_pole_predilated_template_gammax_azimu_b_abw___,1,3);
point_pole_predilated_template_gammax_azimu_b_var_ab__ = var(point_pole_predilated_template_gammax_azimu_b_abw___,1,3);
%%%%%%%%;
% k_c_r holds the distance (in k_c_) between each point_pole_predilated_template_gammax_k_c_ and the avg. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_r_abw___ = ...
sqrt( ...
 bsxfun(@minus,point_pole_predilated_template_gammax_k_c_0_abw___,mean(point_pole_predilated_template_gammax_k_c_0_abw___,3)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_1_abw___,mean(point_pole_predilated_template_gammax_k_c_1_abw___,3)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_2_abw___,mean(point_pole_predilated_template_gammax_k_c_2_abw___,3)).^2 ...
);
%point_pole_predilated_template_gammax_k_c_r_avg_ab__ = sqrt(mean(point_pole_predilated_template_gammax_k_c_r_abw___.^2,3)); %<-- take 2-norm, rather than mean norm. ;
point_pole_predilated_template_gammax_k_c_r_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_r_abw___,3); %<-- take mean-norm, rather than 2-norm. ;

%%%%%%%%;
% Notation: ;
% g = g_dilation, f = I + g*epsilon, where we assume g_dilation(b) = epsilon*sin(2*b) ;
% Assume that (a,b) are fixed for a particular point on the sphere: a = polar_a, b = azimu_b. ;
% The value of c = gamma_z will parametrize the poles associated with that point. ;
%%%%;
% Now, all the pole locations for the point (a,b) are given by the same forumula used in get_template:
% (x,y,z) = [ cb*ca*cc - sb*sc , sb*ca*cc + cb*sc , -sa*cc ] ;
% The norm (x^2 + y^2) = (r01)^2 = ( (cb*ca*cc - sb*sc)^2 + (sb*ca*cc + cb*sc)^2 ) ;
% = ( (cb*ca*cc)^2 + (sb*sc)^2 - 2*(cb*ca*cc*sb*sc) + (sb*ca*cc)^2 + (cb*sc)^2 + 2*(sb*ca*cc*cb*sc) ) ;
% = ( (cb*ca*cc)^2 + (sb*sc)^2 + (sb*ca*cc)^2 + (cb*sc)^2 ) ;
% = (ca*cc)^2 + (sc)^2 ;
%%%%;
% The pole-azimu_b at (x,y,z) is: atan2(y,x) = atan((sb*ca*cc+cb*sc)/(cb*ca*cc-sb*sc)) =: phi(a,b,c). ;
%%%%;
% Note that, obviously, as any particular pole is perturbed azimuthally, ;
% the associated point on the sphere (associated with that pole) is perturbed azimuthally by the same amount. ;
% Thus, for any pole, the pole_predilated_azimu_b = pole_original_azimu_b - epsilon*g_dilation(pole_original_azimu_b) ;
% or, more simply: phi_new - phi_old = -epsilon*sin(2*phi_old) ;
% which is the azimuthal change in the point induced by the azimuthal change in the pole. ;
%%%%;
% Now we can estimate the average azimu_b-shift for any point on the sphere by averaging over c (parametrizing the poles associated with that point). ;
% \int_{d\gamma_z} -epsilon * sin(2*phi) = -epsilon * \int_{d\gamma_z} 2*sin(phi)*cos(phi) = -epsilon * \int_{d\gamma_z} * 2 * (x/r01) * (y/r01) ;
% = -epsilon \int_{d\gamma_z} 2 * ( (sb*ca*cc + cb*sc)*(cb*ca*cc - sb*sc) ) / ( (ca*cc)^2 + (sc)^2 ) ;
%%%%;
% Now seeking a separable solution we fix the azimu_b at, say, pi/4 to get:
% = -2*epsilon \int_{d\gamma_z} ( (sqrt(0.5)*ca*cc + sqrt(0.5)*sc)*(sqrt(0.5)*ca*cc - sqrt(0.5)*sc) ) / ( (ca*cc)^2 + (sc)^2 ) ;
% = -epsilon \int_{d\gamma_z} ( (ca*cc)^2 - (sc)^2 ) / ( (ca*cc)^2 + (sc)^2 ) ;
% = -epsilon \int_{d\gamma_z} ( (ca)^2*(cc)^2 - (sc)^2 ) / ( (ca)^2*(cc)^2 + (sc)^2 ) ;
%%%%;
% Now using Gradshteyn & Ryzhik: 3.647 (p402) BI (47)(20) we see: ;
% \int_{0}^{pi/2} dx (cos(x)^p * cos(px)) / (a^2*sin(x)^2 + b^2*cos(x)^2) = (pi/(2*b)) * a^(p-1) / (a+b)^p ;
% (1/(2*pi)) * \int_{0}^{2*pi} dx (cos(x)^2)/(sin(x)^2 + b^2*cos(x)^2) = 1/(b*(1+b)) ;
% +\int_{gamma_z} ( alpha^2 * cos(gamma)^2 )/ ( alpha^2*cos(gamma)^2 + sin(gamma)^2 ) = +alpha^2/(alpha*(1+alpha)) = -alpha/(1+alpha) ;
% +\int_{gamma_z} ( sin(gamma)^2 )/ ( alpha^2*cos(gamma)^2 + sin(gamma)^2 ) = +(1/(alpha*(1+alpha)) - 1/alpha) ;
% +\int_{gamma_z} ( alpha^2*cos(gamma)^2 - sin(gamma)^2 )/ ( alpha^2*cos(gamma)^2 + sin(gamma)^2 ) =  +alpha/(1+alpha) + (1/(alpha*(1+alpha)) - 1/alpha);
% = (alpha^2 - alpha) / (alpha^2 + alpha) = (alpha-1)/(alpha+1) ;
% Do not forget absolute-values: (abs(alpha)-1)/(abs(alpha)+1) ;
%%%%%%%%;
%% Aside: as we perturb the point-polar_a a, the pole-azimu_b phi changes as \partial_{a} phi. ;
%% Note that: datan(y/x) = 1/(1+y^2/x^2) * (y'/x - yx'/x^2) = (x^2)/(x^2+y^2)*(y'x - yx')/(x^2) = (y'x-yx')/(x^2+y^2) ;
%% Thus, \partial_{a} phi = ( (-sb*sa*cc)*(db*ca*cc-sb*sc) - (sb*ca*cc+cb*sc)*(-cb*sa*cc) ) / ( (ca*cc)^2 + (sc)^2 ) ;
%% = ( (cc)^2*(-sb*sa*cb*ca) + (sb)^2*(sa*cc*sc) + (cc)^2*(+sb*ca*cb*sa) + (cb)^2*(sa*cc*sc) ) / ( (ca*cc)^2 + (sc)^2 ) ;
%% = (sa*cc*sc) / ( (ca*cc)^2 + (sc)^2 ) ;
%%%%%%%%;
tmp_x_abw___ = bsxfun(@times,cos(point_output_azimu_b_ab__).*cos(point_output_polar_a_ab__-pi/2),reshape(cos(gamma_z_),[1,1,n_w_max])) - bsxfun(@times,sin(point_output_azimu_b_ab__),reshape(sin(gamma_z_),[1,1,n_w_max]));
tmp_y_abw___ = bsxfun(@times,sin(point_output_azimu_b_ab__).*cos(point_output_polar_a_ab__-pi/2),reshape(cos(gamma_z_),[1,1,n_w_max])) + bsxfun(@times,cos(point_output_azimu_b_ab__),reshape(sin(gamma_z_),[1,1,n_w_max]));
tmp_z_abw___ = sqrt(tmp_x_abw___.^2 + tmp_y_abw___.^2);
tmp_w_abw___ = sqrt(bsxfun(@plus,bsxfun(@times,cos(point_output_polar_a_ab__-pi/2),reshape(cos(gamma_z_),[1,1,n_w_max])).^2,reshape(sin(gamma_z_),[1,1,n_w_max]).^2));
disp(sprintf(' %% tmp_z_abw___ vs tmp_w_abw___: %0.16f',fnorm(tmp_z_abw___ - tmp_w_abw___)/fnorm(tmp_z_abw___)));
tmp_u_ab__ = mean( -equa_band_dilated_amplitude*2*tmp_x_abw___.*tmp_y_abw___./max(1e-12,tmp_w_abw___.^2) , 3 );
point_pole_predilated_template_gammax_azimu_b_av2_ab__ = tmp_u_ab__;
point_pole_predilated_template_gammax_azimu_b_av1_ab__ = periodize(point_pole_predilated_template_gammax_azimu_b_avg_ab__ - point_output_azimu_b_ab__,-pi,+pi);
subplot(1,2,1);imagesc(point_pole_predilated_template_gammax_azimu_b_av1_ab__,equa_band_dilated_amplitude*[-1,+1]);
axisnotick; xlabel('azimu_b','Interpreter','none'); ylabel('polar_a','Interpreter','none')
title('empirical');
subplot(1,2,2);imagesc(point_pole_predilated_template_gammax_azimu_b_av2_ab__,equa_band_dilated_amplitude*[-1,+1]);
axisnotick; xlabel('azimu_b','Interpreter','none'); ylabel('polar_a','Interpreter','none')
title('formula');
%%%%%%%%;

%%%%%%%%;
% fix cb = sb = sqrt(0.5). ;
%%%%%%%%;
tmp = sqrt(0.5); tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a));
tmp_x0_aw__ = bsxfun(@times,tmp.*cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])) - bsxfun(@times,tmp,reshape(sin(gamma_z_),[1,n_w_max]));
tmp_y0_aw__ = bsxfun(@times,tmp.*cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])) + bsxfun(@times,tmp,reshape(sin(gamma_z_),[1,n_w_max]));
tmp_z0_aw__ = sqrt(tmp_x0_aw__.^2 + tmp_y0_aw__.^2);
tmp_w0_aw__ = sqrt(bsxfun(@plus,bsxfun(@times,cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])).^2,reshape(sin(gamma_z_),[1,n_w_max]).^2));
disp(sprintf(' %% tmp_z0_aw__ vs tmp_w0_aw__: %0.16f',fnorm(tmp_z0_aw__ - tmp_w0_aw__)/fnorm(tmp_z0_aw__)));
tmp_u0_a_ = mean( -equa_band_dilated_amplitude*2*tmp_x0_aw__.*tmp_y0_aw__./max(1e-12,tmp_w0_aw__.^2) , 2 );
tmp_u0_ab__ = bsxfun(@times,tmp_u0_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
disp(sprintf(' %% tmp_u_ab__ vs tmp_u0_ab__: %0.16f',fnorm(tmp_u_ab__ - tmp_u0_ab__)/fnorm(tmp_u_ab__))); %<-- should be order equa_band_dilated_amplitude. ;
nb=round((n_point_b-1)/8);
fnorm(tmp_u_ab__(:,1+nb) - tmp_u0_a_)/fnorm(tmp_u_ab__(:,1+nb));
disp(sprintf(' %% tmp_u_ab__(:,1+nb) vs tmp_u0_a_: %0.16f',fnorm(tmp_u_ab__(:,1+nb) - tmp_u0_a_)/fnorm(tmp_u_ab__(:,1+nb)))); %<-- should be order equa_band_dilated_amplitude. ;
%%%%%%%%;
tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a));
tmp_z1_aw__ = bsxfun(@minus,bsxfun(@times,cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])).^2,reshape(sin(gamma_z_),[1,n_w_max]).^2);
tmp_w1_aw__ = bsxfun(@plus ,bsxfun(@times,cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])).^2,reshape(sin(gamma_z_),[1,n_w_max]).^2);
tmp_u1_a_ = mean( -equa_band_dilated_amplitude*tmp_z1_aw__./max(1e-12,tmp_w1_aw__) , 2 );
tmp_u1_ab__ = bsxfun(@times,tmp_u1_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
disp(sprintf(' %% tmp_u_ab__ vs tmp_u1_ab__: %0.16f',fnorm(tmp_u_ab__ - tmp_u1_ab__)/fnorm(tmp_u_ab__))); %<-- should be order equa_band_dilated_amplitude. ;
nb=round((n_point_b-1)/8);
fnorm(tmp_u_ab__(:,1+nb) - tmp_u1_a_)/fnorm(tmp_u_ab__(:,1+nb));
disp(sprintf(' %% tmp_u_ab__(:,1+nb) vs tmp_u1_a_: %0.16f',fnorm(tmp_u_ab__(:,1+nb) - tmp_u1_a_)/fnorm(tmp_u_ab__(:,1+nb)))); %<-- should be order equa_band_dilated_amplitude. ;
%%%%%%%%;
% Gradshteyn & Ryzhik: 3.647 p402: ;
% b0=rand();tmp_z_ = cos(gamma_z_).^2./(sin(gamma_z_).^2 + b0^2*cos(gamma_z_).^2); mean(tmp_z_),; 1/b0/(1+b0),;
% b0=rand();tmp_z_ = 1./(sin(gamma_z_).^2 + b0^2*cos(gamma_z_).^2); mean(tmp_z_),; 1/b0,;
% b0=rand();tmp_z_ = (b0^2*cos(gamma_z_).^2 - sin(gamma_z_).^2)./(b0^2*cos(gamma_z_).^2 + sin(gamma_z_).^2); mean(tmp_z_),; (b0^2+1)/b0/(1+b0) - 1/b0,;
% b0=rand();tmp_z_ = (b0^2*cos(gamma_z_).^2 - sin(gamma_z_).^2)./(b0^2*cos(gamma_z_).^2 + sin(gamma_z_).^2); mean(tmp_z_),; (abs(b0)-1)/(abs(b0)+1),;
%%%%%%%%;
tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a)); tmp_ca_ = cos(tmp_output_polar_a_a_-pi/2);
tmp_u2_a_ = -equa_band_dilated_amplitude*(abs(tmp_ca_)-1)./(abs(tmp_ca_)+1);
tmp_u2_ab__ = bsxfun(@times,tmp_u2_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
disp(sprintf(' %% tmp_u_ab__ vs tmp_u2_ab__: %0.16f',fnorm(tmp_u_ab__ - tmp_u2_ab__)/fnorm(tmp_u_ab__))); %<-- should be order equa_band_dilated_amplitude. ;
nb=round((n_point_b-1)/8);
fnorm(tmp_u_ab__(:,1+nb) - tmp_u2_a_)/fnorm(tmp_u_ab__(:,1+nb));
disp(sprintf(' %% tmp_u_ab__(:,1+nb) vs tmp_u2_a_: %0.16f',fnorm(tmp_u_ab__(:,1+nb) - tmp_u2_a_)/fnorm(tmp_u_ab__(:,1+nb)))); %<-- should be order equa_band_dilated_amplitude. ;
%%%%%%%%;

%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; figbeach;
tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a)); tmp_ca_ = cos(tmp_output_polar_a_a_-pi/2);
tmp_u3_a_ = equa_band_dilated_amplitude .* sqrt(1 - ((abs(tmp_ca_)-1)./(abs(tmp_ca_)+1)).^2);
tmp_u2_a_ = -equa_band_dilated_amplitude*(abs(tmp_ca_)-1)./(abs(tmp_ca_)+1);
tmp_u2_ab__ = bsxfun(@times,tmp_u2_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
subplot(1,3,1);
plot(point_output_polar_a_a_,tmp_u3_a_,'m.-','LineWidth',3);
xlim([0,pi]); xlabel('point_output_polar_a_a_','Interpreter','none'); ylabel('amplitude');
set(gca,'YTick',[],'XTick',pi*[0,0.5,1],'XTickLabel',{'0','\pi/2','\pi'});
title('$\alpha\in[0,2\pi]$','Interpreter','latex');
subplot(1,3,[2,3]);
imagesc(tmp_u2_ab__,max(tmp_u2_ab__,[],'all')*[-1,+1]);
set(gca,'YTick',[1,(n_point_a+1)/2,n_point_a],'YTickLabel',{'0','\pi/2','\pi'});
ylabel('point_output_polar_a_a','Interpreter','none');
set(gca,'XTick',[1,(n_point_b+1)/2,n_point_b],'XTickLabel',{'0','\pi','2\pi'});
xlabel('point_output_azimu_b','Interpreter','none');
title('$\alpha$, $\beta$ varying across sphere','Interpreter','latex');
%%%%;    
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_std_from_viewing_angle_FIGA',dir_manuscript_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% Writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
sgtitle('','Interpreter','none');
print('-djpeg',sprintf('%s_stripped.jpg',fname_fig_pre));
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
%close(gcf);
%%%%%%%%;

%%%%%%%%;
% Now estimate:
% point_pole_predilated_template_gammax_azimu_b_avg_ab__, ;
% point_pole_predilated_template_gammax_azimu_b_std_ab__, ;
% for a variety of viewing-angle distributions. ;
%%%%%%%%;
% We assume the viewing-angle distribution is concentrated around the equatorial-belt: ;
% rho(polar_a) \propto exp( -cos(polar_a)^2 / (2*sigma_a^2) ) ;
%%%%%%%%;
lsigma_a_ = transpose(-4:0.125:0); n_lsigma_a = numel(lsigma_a_);
point_pole_predilated_template_gammax_azimu_b_avg_abs___ = zeros(n_point_a,n_point_b,n_lsigma_a);
point_pole_predilated_template_gammax_azimu_b_std_abs___ = zeros(n_point_a,n_point_b,n_lsigma_a);
for nlsigma_a=0:n_lsigma_a-1;
lsigma_a = lsigma_a_(1+nlsigma_a);
sigma_a = exp(lsigma_a);
%%%%;
tmp_weight_abw___ = zeros(n_point_a,n_point_b,n_w_max);
for npoint_a=0:n_point_a-1;
for npoint_b=0:n_point_b-1;
tmp_weight_w_ = zeros(n_w_max,1);
tmp_c_ = cos(point_pole_polar_a_abw___(1+npoint_a,1+npoint_b,:)); tmp_c_min = min(tmp_c_);
tmp_nlweight_w_ = +(tmp_c_).^2/(2*sigma_a^2); %<-- -log(weight_w_). ;
tmp_w_min = min(tmp_nlweight_w_);
tmp_weight_w_ = exp(-(tmp_nlweight_w_-tmp_w_min)); %<-- at least one positive value. ;
tmp_weight_w_ = tmp_weight_w_/max(1e-12,sum(tmp_weight_w_));
tmp_weight_abw___(1+npoint_a,1+npoint_b,:) = tmp_weight_w_;
end;%for npoint_b=0:n_point_b-1;
end;%for npoint_a=0:n_point_a-1;
%%%%;
tmp_point_pole_predilated_template_gammax_k_c_0_avg_ab__ = sum(point_pole_predilated_template_gammax_k_c_0_abw___.*tmp_weight_abw___,3);
tmp_point_pole_predilated_template_gammax_k_c_1_avg_ab__ = sum(point_pole_predilated_template_gammax_k_c_1_abw___.*tmp_weight_abw___,3);
tmp_point_pole_predilated_template_gammax_k_c_2_avg_ab__ = sum(point_pole_predilated_template_gammax_k_c_2_abw___.*tmp_weight_abw___,3);
tmp_point_pole_predilated_template_gammax_azimu_b_avg_ab__ = sum(point_pole_predilated_template_gammax_azimu_b_abw___.*tmp_weight_abw___,3);
tmp_point_pole_predilated_template_gammax_azimu_b_std_ab__ = sqrt(sum(bsxfun(@minus,point_pole_predilated_template_gammax_azimu_b_abw___,tmp_point_pole_predilated_template_gammax_azimu_b_avg_ab__).^2.*tmp_weight_abw___,3));
point_pole_predilated_template_gammax_azimu_b_avg_abs___(:,:,1+nlsigma_a) = tmp_point_pole_predilated_template_gammax_azimu_b_avg_ab__;
point_pole_predilated_template_gammax_azimu_b_std_abs___(:,:,1+nlsigma_a) = tmp_point_pole_predilated_template_gammax_azimu_b_std_ab__;
end;%for nlsigma_a=0:n_lsigma_a-1;

%%%%%%%%;
% Note that the average azimuthal displacement is not strongly affected by sigma_a. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
p_row = 4;p_col=8;np=0;
for nlsigma_a=0:1:min(n_lsigma_a,p_row*p_col)-1;
subplot(p_row,p_col,1+np);np=np+1;
tmp_avg_ab__ = point_pole_predilated_template_gammax_azimu_b_avg_abs___(:,:,1+nlsigma_a);
tmp_avg_ab__ = tmp_avg_ab__ - point_output_azimu_b_ab__;
tmp_avg_ab__ = periodize(tmp_avg_ab__,-pi,+pi);
tmp_avg_ab__(:,round((1+n_point_b)/2))=0;
imagesc(tmp_avg_ab__,+1.5e-2*[-1,+1]);
title(sprintf('nlsigma_a %d <-- %0.2f',nlsigma_a,exp(lsigma_a_(1+nlsigma_a))),'Interpreter','none'); axisnotick;
end;%for nlsigma_a=0:min(n_lsigma_a,p_row*p_col)-1;
%%%%;
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_avg_vs_sigma_FIGA',dir_manuscript_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% Writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
sgtitle('','Interpreter','none');
print('-djpeg',sprintf('%s_stripped.jpg',fname_fig_pre));
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;

%%%%%%%%;
% However, the average azimuthal diffusion is strongly affected by sigma_a. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
p_row = 4;p_col=8;np=0;
for nlsigma_a=0:1:min(n_lsigma_a,p_row*p_col)-1;
subplot(p_row,p_col,1+np);np=np+1;
tmp_std_ab__ = point_pole_predilated_template_gammax_azimu_b_std_abs___(:,:,1+nlsigma_a);
tmp_std_ab__(:,round((1+n_point_b)/2))=tmp_std_ab__(:,1+0);
imagesc(tmp_std_ab__,[0,+1.5e-2]);
title(sprintf('nlsigma_a %d <-- %0.2f',nlsigma_a,exp(lsigma_a_(1+nlsigma_a))),'Interpreter','none'); axisnotick;
end;%for nlsigma_a=0:1:min(n_lsigma_a,p_row*p_col)-1;
%%%%;
fname_fig_pre = sprintf('%s/M3d_shape_longitudinal_perturbation_std_vs_sigma_FIGA',dir_manuscript_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if ( flag_replot | ~exist(fname_fig_jpg,'file'));
if (flag_verbose> 0); disp(sprintf(' %% writing %s',fname_fig_pre)); end;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% Writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
sgtitle('','Interpreter','none');
print('-djpeg',sprintf('%s_stripped.jpg',fname_fig_pre));
end;%if ( flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;












