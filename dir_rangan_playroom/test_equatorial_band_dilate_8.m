%%%%%%%%;
% tests equatorial band dilation. ;
% focuses on calculation. ;
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

%%%%%%%%;
flag_verbose = 0;
flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;
%%%%;
dir_manuscript = sprintf('/%s/rangan/dir_cryoem/dir_spurious_heterogeneity_manuscript',string_root);
if ~exist(dir_manuscript,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript)); mkdir(dir_manuscript); end;
dir_manuscript_jpg = sprintf('%s/dir_M3d_shape_longitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;
%%%%%%%%;

%%%%%%%%;
% First set up collection of templates. ;
%%%%%%%%;
flag_verbose = 0; flag_disp=1; nf=0;
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
 flag_verbose ...
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
 flag_verbose ...
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
if (flag_verbose>0); disp(sprintf(' %% ntest %d/%d tmp_error %0.16f;',ntest,n_test,tmp_error)); end;
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

%%%%%%%%;
% Now we collect the following arrays: ;
% point_output_*_ab__ : the variables corresponding to points on the surface of the sphere ;
% These are each (n_point_a,n_point_b) arrays. ;
% (note that point_output_polar_a_ begins at -pi/2) ;
% point_pole_*_abw___ : the variables corresponding to poles associated with points on the surface of the sphere. ;
% These are each (n_point_a,n_point_b,n_w_max) arrays. ;
% The point_pole_*_abw___(1+npoint_a,1+npoint_b,:) correspond to point_output_*_ab__(1+npoint_a,1+npoint_b). ;
%%%%%%%%;
n_point_a = 1+ 64; polar_a_a_ = linspace(0,1*pi,n_point_a+2); polar_a_a_ = transpose(polar_a_a_(2:end-1)); %<-- not periodic. ;
n_point_b = 1+128; azimu_b_b_ = linspace(0,2*pi,n_point_b+0); azimu_b_b_ = transpose(azimu_b_b_(1:end-0)); %<-- yes periodic. ;
point_output_azimu_b_ab__ = zeros(n_point_a,n_point_b,1);
point_output_polar_a_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_0_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_1_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_2_ab__ = zeros(n_point_a,n_point_b,1);
point_pole_azimu_b_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_polar_a_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_0_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_1_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_2_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_gamma_z_abw___ = zeros(n_point_a,n_point_b,n_w_max);
tmp_error_0 = 0; tmp_error_1 = 0; tmp_error_2 = 0; tmp_error_3 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npoint_a=0:n_point_a-1; 
if (flag_verbose>-1); disp(sprintf(' %% npoint_a %d/%d ; tmp_error_ [%0.16f,%0.16f,%0.16f,%0.16f]',npoint_a,n_point_a,tmp_error_0,tmp_error_1,tmp_error_2,tmp_error_3)); end;
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
point_output_azimu_b = azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
%point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_output_k_c_ = [+cos(point_output_azimu_b)*cos(point_output_polar_a);+sin(point_output_azimu_b)*cos(point_output_polar_a);-sin(point_output_polar_a)] ;
point_output_k_r01 = sqrt(point_output_k_c_(1+0).^2 + point_output_k_c_(1+1).^2);
point_output_k_c_0_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+0);
point_output_k_c_1_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+1);
point_output_k_c_2_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+2);
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_pole_k_c_0_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+0,:);
point_pole_k_c_1_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+1,:);
point_pole_k_c_2_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+2,:);
point_output_azimu_b_ab__(1+npoint_a,1+npoint_b) = point_output_azimu_b;
point_output_polar_a_ab__(1+npoint_a,1+npoint_b) = point_output_polar_a;
%%%%%%%%%%%%%%%%;
for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
% Now we step through each of the templates associated with the point point_output. ;
%%%%%%%%;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
point_pole_azimu_b_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_azimu_b;
point_pole_polar_a_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_polar_a;
%%%%%%%%;
% Eventually we will determine the predilated template that is associated with the original template mapped to point_output. ;
% Here we determine point_pole_gamma_z, ;
% which is the angle within the point_pole_template_k_c__ ;
% (denoted gamma_z below) ;
% such that the point_pole_template_k_c__ evaluated at that gamma_z ;
% (denoted point_pole_template_gamma_z_k_c_) ;
% is the same as point_output_k_c_. ;
%%%%%%%%;
point_pole_template_gamma0_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;-sin(point_pole_polar_a) ...
];
point_pole_template_sgx_k_c_ = cross(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_template_sgx = dot(point_pole_template_sgx_k_c_,point_pole_k_c_);
point_pole_template_cgx = dot(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_gamma_z = atan2(point_pole_template_sgx,point_pole_template_cgx);
point_pole_gamma_z_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_gamma_z;
point_pole_template_gamma_z_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gamma_z) - sin(point_pole_azimu_b)*sin(point_pole_gamma_z) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gamma_z) + cos(point_pole_azimu_b)*sin(point_pole_gamma_z) ...
;-sin(point_pole_polar_a)*cos(point_pole_gamma_z) ...
];
tmp_error_0 = tmp_error_0 + fnorm(point_pole_template_gamma_z_k_c_ - point_output_k_c_); %<-- should be 0. ;
point_pole_template_gamma_z_k_r01 = sqrt(point_pole_template_gamma_z_k_c_(1+0).^2 + point_pole_template_gamma_z_k_c_(1+1).^2);
tmp_error_1 = tmp_error_1 + fnorm((cos(point_pole_polar_a)*cos(point_pole_gamma_z))^2 + (sin(point_pole_gamma_z))^2 - point_pole_template_gamma_z_k_r01^2); %<-- should be 0. ;
point_pole_template_gamma_z_k_r012 = sqrt(point_pole_template_gamma_z_k_c_(1+0).^2 + point_pole_template_gamma_z_k_c_(1+1).^2 + point_pole_template_gamma_z_k_c_(1+2).^2);
point_pole_template_gamma_z_azimu_b = atan2(point_pole_template_gamma_z_k_c_(1+1),point_pole_template_gamma_z_k_c_(1+0));
tmp_error_2 = tmp_error_2 + fnorm(periodize(point_output_azimu_b - point_pole_template_gamma_z_azimu_b,-pi,+pi)); %<-- should be 0. ;
point_pole_template_gamma_z_polar_a = atan2(point_pole_template_gamma_z_k_r01,point_pole_template_gamma_z_k_c_(1+2));
tmp_error_3 = tmp_error_3 + fnorm(periodize(+pi/2+point_output_polar_a - point_pole_template_gamma_z_polar_a,-pi/2,+pi/2)); %<-- should be 0. ;
%%%%%%%%;
clear point_pole_k_c_;
clear point_pole_k_r01;
clear point_pole_k_r012;
clear point_pole_azimu_b;
clear point_pole_polar_a;
clear point_pole_template_gamma0_k_c_;
clear point_pole_template_sgx_k_c_;
clear point_pole_template_sgx;
clear point_pole_template_cgx;
clear point_pole_gamma_z;
clear point_pole_template_gamma_z_k_c_;
clear point_pole_template_gamma_z_k_r01;
clear point_pole_template_gamma_z_k_r012;
clear point_pole_template_gamma_z_azimu_b;
clear point_pole_template_gamma_z_polar_a;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
clear point_output_azimu_b;
clear point_output_polar_a;
clear point_output_gamma_z = 0;
clear point_output_k_c_;
clear point_output_k_r01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;end;%for npoint_a=0:n_point_a-1; for npoint_b=0:n_point_b-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now if we have a particular (continuous) perturbation field: ;
%%%%;
% g_polar_a = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_polar_a = @(point_pole_polar_a,point_pole_azimu_b) point_pole_polar_a + epsilon_polar_a * g_polar_a(point_pole_polar_a,point_pole_azimu_b) ;
% g_azimu_b = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_azimu_b = @(point_pole_polar_a,point_pole_azimu_b) point_pole_azimu_b + epsilon_azimu_b * g_azimu_b(point_pole_polar_a,point_pole_azimu_b) ;
% g_gamma_z = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_gamma_z = @(point_pole_polar_a,point_pole_azimu_b) point_pole_gamma_z + epsilon_gamma_z * g_gamma_z(point_pole_polar_a,point_pole_azimu_b) ;
%%%%;
% Assuming that epsilon_ = [epsilon_polar_a,epsilon_azimu_b,epsilon_gamma_z] is small, the inverse (to first order) is given by: ;
% f_polar_a_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_polar_a_inverse = point_pole_polar_a - epsilon_polar_a * g_polar_a(point_pole_polar_a,point_pole_azimu_b) ;
% f_azimu_b_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_azimu_b_inverse = point_pole_azimu_b - epsilon_azimu_b * g_azimu_b(point_pole_polar_a,point_pole_azimu_b) ;
% f_gamma_z_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_gamma_z_inverse = point_pole_gamma_z - epsilon_gamma_z * g_gamma_z(point_pole_polar_a,point_pole_azimu_b) ;
% with the final formula (for gamma_z_inverse) being exact. ;
%%%%%%%%;
% Thus, the location on the sphere corresponding to: ;
% point_pole_predilated_template_gamma_z_k_c_ = [...
%  +cos(point_pole_azimu_b_inverse)*cos(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) - sin(point_pole_azimu_b_inverse)*sin(point_pole_gamma_z_inverse) ...
% ;+sin(point_pole_azimu_b_inverse)*cos(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) + cos(point_pole_azimu_b_inverse)*sin(point_pole_gamma_z_inverse) ...
% ;-sin(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) ...
% ];
%%%%%%%%;
% Using the notation: ;
% ca, sa = cos and sin of polar_a, ;
% cb, sb = cos and sin of azimu_b, ;
% cc, sc = cos and sin of gamma_z, ;
% we have the relations: ;
% x0 = +cb*ca*cc - sb*sc ;
% x1 = +sb*ca*cc + cb*sc ;
% x2 = -sa*cc ;
% and the jacobian: ;
% dx0da = -cb*sa*cc ;
% dx0db = -sb*ca*cc - cb*sc ;
% dx0dc = -cb*ca*sc - sb*cc ;
% dx1da = -sb*sa*cc ;
% dx1db = +cb*ca*cc - sb*sc ;
% dx1dc = -sb*ca*sc + cb*cc ;
% dx2da = -ca*cc ;
% dx2db = 0;
% dx2dc = +sa*sc ;
%%%%%%%%;
ca_abw___ = cos(point_pole_polar_a_abw___); sa_abw___ = sin(point_pole_polar_a_abw___);
cb_abw___ = cos(point_pole_azimu_b_abw___); sb_abw___ = sin(point_pole_azimu_b_abw___);
cc_abw___ = cos(point_pole_gamma_z_abw___); sc_abw___ = sin(point_pole_gamma_z_abw___);
dx0da_abw___ = -cb_abw___.*sa_abw___.*cc_abw___ ;
dx0db_abw___ = -sb_abw___.*ca_abw___.*cc_abw___ - cb_abw___.*sc_abw___ ;
dx0dc_abw___ = -cb_abw___.*ca_abw___.*sc_abw___ - sb_abw___.*cc_abw___ ;
dx1da_abw___ = -sb_abw___.*sa_abw___.*cc_abw___ ;
dx1db_abw___ = +cb_abw___.*ca_abw___.*cc_abw___ - sb_abw___.*sc_abw___ ;
dx1dc_abw___ = -sb_abw___.*ca_abw___.*sc_abw___ + cb_abw___.*cc_abw___ ;
dx2da_abw___ = -ca_abw___.*cc_abw___ ;
dx2db_abw___ = +zeros(n_point_a,n_point_b,n_w_max);
dx2dc_abw___ = +sa_abw___.*sc_abw___ ;
%%%%%%%%;
% Now, the accumulated dx0, dx1, dx2 are given by: ;
% g_a_abw___ = bsxfun(@plus,g_polar_a_abwj____,reshape(eps_a_j_,[1,1,1,n_a_j]));
% g_b_abw___ = bsxfun(@plus,g_azimu_b_abwj____,reshape(eps_b_j_,[1,1,1,n_b_j]));
% g_c_abw___ = bsxfun(@plus,g_gamma_z_abwj____,reshape(eps_c_j_,[1,1,1,n_c_j]));
% dx0_abw___ = -dx0da_abw___.*g_a_abw___ - dx0db_abw___.*g_b_abw___ - dx0dc_abw___.*g_c_abw___ ;
% dx1_abw___ = -dx1da_abw___.*g_a_abw___ - dx1db_abw___.*g_b_abw___ - dx1dc_abw___.*g_c_abw___ ;
% dx2_abw___ = -dx2da_abw___.*g_a_abw___ - dx2db_abw___.*g_b_abw___ - dx2dc_abw___.*g_c_abw___ ;
% or, written in terms of the eps_a_j_, eps_b_j_, eps_c_j_, ;
% dx0_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx0da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx0db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx0dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% dx1_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx1da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx1db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx1dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% dx2_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx2da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx2db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx2dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% collectively producing the shift: ;
% dx_abw3____ = cat(4,dx0_abw___,dx1_abw___,dx2_abw___);
% Note that this shift, while defined in terms of a vector in \Real^{3}, ;
% is constrained to lie tangent to the sphere. ;
%%%%%%%%;
% Now we can record viewing_weight_empirical_abw___ as a function of point_pole_polar_a_abw___ and point_pole_azimu_b_abw___. ;
% If, e.g., this is uniform, then viewing_weight_empirical_abw___ = constant. ;
% Note that this does not account for quadrature on the sphere, ;
% as these weights will be interpreted as contributions from each pole on an equispaced ring (i.e., n_pole = n_w_max). ;
% With this measure we can calculate various derivatives. ;
% The viewing_weight associated with each ring is: ;
% viewing_ring_weight_empirical_abw___ = bsxfun(@rdivide,viewing_weight_empirical_abw___,max(1e-12,sum(viewing_weight_empirical_abw___,3))). ;
% Once again, this is interpreted as a function of point_pole_polar_a_abw___ and point_pole_azimu_b_abw___. ;
%%%%%%%%;
% The average shift at each point_output_k_c_ (indexed by ab13) is now determined by: ;
% dx_avg_ab13___ = sum(bsxfun(@times,dx_abw3___,viewing_ring_weight_empirical_abw___),3). ;
% Meanwhile, the average variation at each point_output_k_c_ is now determined by: ;
% dx_var_ab__ = sum(bsxfun(@times,bsxfun(@minus,dx_abw3___,dx_avg_ab13___).^2,viewing_ring_weight_empirical_abw___),[3,4]). ;
%%%%%%%%;
% Note that the dx_var_ab__ can be written as a quadratic-function of eps_: [eps_a_j_;eps_b_j_;eps_c_j_]. ;
% This immediately implies that dx_var_ab__ can be optimized (subject to a norm-constraint on eps_) ;
% by taking the eigenvectors of the appropriate quadratic-kernel. ;
%%%%%%%%;

%%%%%%%%;
% quick check of tmp_b_est. ;
%%%%%%%%;
tmp_a = 1*pi*rand(); ca = cos(tmp_a); sa = sin(tmp_a);
tmp_b = 2*pi*rand(); cb = cos(tmp_b); sb = sin(tmp_b);
tmp_c = 2*pi*rand(); cc = cos(tmp_c); sc = sin(tmp_c);
tmp_x0 = cb*ca*cc - sb*sc ;
tmp_x1 = sb*ca*cc + cb*sc ;
tmp_x2 = -sa*cc ;
tmp_r01 = sqrt(tmp_x0.^2 + tmp_x1.^2);
tmp_b_out = atan2(tmp_x1,tmp_x0);
tmp_a_out = atan2(tmp_r01,tmp_x2);
tmp_b_est = atan2(sc,ca*cc) + tmp_b;
tmp_x_ = Rz(tmp_b)*Ry(tmp_a)*Rz(tmp_c)*[1;0;0];
disp(sprintf(' %% fnorm tmp_x_ - [tmp_x0;tmp_x1;tmp_x2]: %0.16f',fnorm(tmp_x_-[tmp_x0;tmp_x1;tmp_x2])));
tmp_f_a = @(x_) atan2(sqrt(x_(1+0).^2 + x_(1+1).^2),x_(1+2));
disp(sprintf(' %% fnorm tmp_f_a - tmp_a_out: %0.16f',fnorm(tmp_f_a(tmp_x_) - tmp_a_out)));
tmp_f_b = @(x_) atan2(x_(1+1),x_(1+0));
disp(sprintf(' %% fnorm tmp_f_b - tmp_b_out: %0.16f',fnorm(tmp_f_b(tmp_x_) - tmp_b_out)));
disp(sprintf(' %% fnorm tmp_b_est - tmp_b_out: %0.16f',fnorm(periodize(tmp_b_est - tmp_b_out,-pi,+pi))));

%%%%%%%%;
% Note that, typically speaking, an arbitrary alignment-perturbation will result in ;
% a nonzero 'first-order-term' (associated with the average). ;
% More specifically, this average will have a magnitude proportional to fnorm(epsilon_), ;
% and will point essentially arbitrarily at every point_output_k_c_. ;
%%%%;
% The consequence of this haphazard average drift is that any particular template ;
% (originally aligned to the unperturbed volume) ;
% will now no longer align with *any* locally-perturbed template from the perturbed volume. ;
% This mis-alignment will result in a first-order shift in the likelihood ;
% (proportional to fnorm(epsilon_)). ;
%%%%;
% In special cases, such as a longitudinal or latitudinal perturbation, ;
% each template from the original volume will align to first-order with a particular ;
% locally-perturbed template from the perturbed volume. ;
% In this case any mis-alignment will result from the second-order shift in the likelihood ;
% (proportional to fnorm(epsilon_).^2) ;
% associated with the local diffusion in the perturbed volume  ;
% induced by dx_var_ab__. ;
%%%%;
% For these special cases we can bound the minimum-curvature of the likelihood-landscape ;
% by solving the eigenvalue-problem described above. ;
%%%%%%%%;

%%%%%%%%;
% To demonstrate this calculation using 
									   












