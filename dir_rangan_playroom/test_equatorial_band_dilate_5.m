%%%%%%%%;
% tests equatorial band dilation. ;
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
point_pole_predilated_template_gammax_k_c_0_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_0_abw___,3);
point_pole_predilated_template_gammax_k_c_1_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_1_abw___,3);
point_pole_predilated_template_gammax_k_c_2_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_2_abw___,3);
point_pole_predilated_template_gammax_polar_a_avg_ab__ = mean(point_pole_predilated_template_gammax_polar_a_abw___,3);
point_pole_predilated_template_gammax_polar_a_std_ab__ = std(point_pole_predilated_template_gammax_polar_a_abw___,1,3);
point_pole_predilated_template_gammax_polar_a_var_ab__ = var(point_pole_predilated_template_gammax_polar_a_abw___,1,3);
point_pole_predilated_template_gammax_azimu_b_avg_ab__ = mean(point_pole_predilated_template_gammax_azimu_b_abw___,3);
point_pole_predilated_template_gammax_azimu_b_std_ab__ = std(point_pole_predilated_template_gammax_azimu_b_abw___,1,3);
point_pole_predilated_template_gammax_azimu_b_var_ab__ = var(point_pole_predilated_template_gammax_azimu_b_abw___,1,3);
point_pole_predilated_template_gammax_k_c_r_abw___ = ...
sqrt( ...
 bsxfun(@minus,point_pole_predilated_template_gammax_k_c_0_abw___,mean(point_pole_predilated_template_gammax_k_c_0_abw___,3)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_1_abw___,mean(point_pole_predilated_template_gammax_k_c_1_abw___,3)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_2_abw___,mean(point_pole_predilated_template_gammax_k_c_2_abw___,3)).^2 ...
);
%point_pole_predilated_template_gammax_k_c_r_avg_ab__ = sqrt(mean(point_pole_predilated_template_gammax_k_c_r_abw___.^2,3)); %<-- take 2-norm, rather than mean norm. ;
point_pole_predilated_template_gammax_k_c_r_avg_ab__ = mean(point_pole_predilated_template_gammax_k_c_r_abw___,3); %<-- take mean-norm, rather than 2-norm. ;

flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot(1,2,1);
imagesc( ...
	 periodize( ...
		    point_pole_predilated_template_gammax_azimu_b_avg_ab__ ...
		    - point_output_azimu_b_ab__ ...
		    - g_dilation(point_output_azimu_b_ab__).*( ...
							       cos(1*point_output_polar_a_ab__) ...
							       - 0.500*sin(2*point_output_polar_a_ab__).^2 ...
							       - 0.250*sin(4*point_output_polar_a_ab__).^2 ...
							       ) ...
		    ,-pi,+pi) ...
	 ,equa_band_dilated_amplitude*[-1,+1] ...
	 );
colorbar;
%%%%;
subplot(1,2,2);
hold on;
plot(point_output_polar_a_a_,2*point_pole_predilated_template_gammax_azimu_b_var_ab__,'.-');
plot(point_output_polar_a_a_,equa_band_dilated_amplitude.^2*(2/pi)*atan(2*pi*abs(point_output_polar_a_a_-pi/2)),'go');
hold off;
ylim(equa_band_dilated_amplitude.^2*[0,+1]);
colorbar;
%%%%;
end;%if flag_disp;


figure(1+nf);nf=nf+1;clf;figsml;
%tmp_ab__ = g_dilation(point_output_azimu_b_ab__).*(cos(1*point_output_polar_a_ab__) - 0.500*sin(2*point_output_polar_a_ab__).^2	- 0.250*sin(4*point_output_polar_a_ab__).^2);
tmp_ab__ = 0*g_dilation(point_output_azimu_b_ab__) .* ( (pi/2-abs(point_output_polar_a_ab__))/(pi/2) ).^2 ;
plot(point_output_polar_a_a_,point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,1:end-1) - point_output_azimu_b_ab__(:,1:end-1) - tmp_ab__(:,1:end-1),'.-');
ylim(equa_band_dilated_amplitude*[-1,+1]);
xlim([0,pi]);


tmp_ij = 5;
tmp_ab__ = g_dilation(point_output_azimu_b_ab__) .* ( ( (pi/2-abs(point_output_polar_a_ab__))/(pi/2) ).^3 );
tmp2_ab__ = g_dilation(point_output_azimu_b_ab__) .* ( equa_band_dilated_amplitude*(1-cos(4*point_output_polar_a_ab__)) );
tmp_ab__ = 0*tmp_ab__; tmp2_ab__ = 0*tmp2_ab__;
plot(point_output_polar_a_a_,point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,tmp_ij) - point_output_azimu_b_ab__(:,tmp_ij) - tmp_ab__(:,tmp_ij),'r-',point_output_polar_a_a_,tmp2_ab__(:,tmp_ij),'g-');
%plot(point_output_polar_a_a_,point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,tmp_ij) - point_output_azimu_b_ab__(:,tmp_ij),'k-',point_output_polar_a_a_,tmp_ab__(:,tmp_ij),'r-');
%plot(point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,tmp_ij) - point_output_azimu_b_ab__(:,tmp_ij) , tmp_ab__(:,tmp_ij).^2,'.-');
%plot( log(point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,tmp_ij) - point_output_azimu_b_ab__(:,tmp_ij)) , log(tmp_ab__(:,tmp_ij)) ,'.-');
%plot(point_output_polar_a_a_,point_pole_predilated_template_gammax_azimu_b_avg_ab__(:,tmp_ij) - point_output_azimu_b_ab__(:,tmp_ij),'.');

figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
tmp_d_ab__ = bsxfun(@minus ...
		    ,point_pole_predilated_template_gammax_azimu_b_avg_ab__ ...
		    ,point_output_azimu_b_ab__ ...
		    );
tmp_e_ab__ = bsxfun(@rdivide ...
		    ,tmp_d_ab__ ...
		    ,sin(2*point_output_azimu_b_ab__) ...
		    ) ...
  /(equa_band_dilated_amplitude);
tmp_e_ab__ = tmp_e_ab__(:,2:end-1);
hold on;
plot(point_output_polar_a_a_,tmp_e_ab__,'.-');
plot(point_output_polar_a_a_,exp(-abs(point_output_polar_a_a_-pi/2)/((2*pi/12))),'g-');
hold off;
ylim([-0.125,+1.125]);
xlabel('point_output_polar_a_a_','Interpreter','none'); xlim([0,pi]); set(gca,'XTick',pi*[0,0.5,1],'XTickLabel',{'0','\pi/2','\pi'});
title('scaled avg','Interpreter','none');
subplot(1,2,2);
tmp_v_ab__ = point_pole_predilated_template_gammax_azimu_b_var_ab__ ...
  /(0.5*equa_band_dilated_amplitude.^2);
tmp_v_ab__ = tmp_v_ab__(:,2:end-1);
hold on;
plot(point_output_polar_a_a_,tmp_v_ab__,'.-');
plot(point_output_polar_a_a_,1-exp(-abs(point_output_polar_a_a_-pi/2)/((2*pi/24))),'g-');
hold off;
ylim([-0.125,+1.125]);
xlabel('point_output_polar_a_a_','Interpreter','none'); xlim([0,pi]); set(gca,'XTick',pi*[0,0.5,1],'XTickLabel',{'0','\pi/2','\pi'});
title('scaled var','Interpreter','none');
% plot(point_output_polar_a_a_,tmp_e_ab__.^2,point_output_polar_a_a_,1-tmp_v_ab__,'o'); ylim([-0.125,+1.125]); %<-- Note overlap. ;

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

tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a));
tmp_z1_aw__ = bsxfun(@minus,bsxfun(@times,cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])).^2,reshape(sin(gamma_z_),[1,n_w_max]).^2);
tmp_w1_aw__ = bsxfun(@plus ,bsxfun(@times,cos(tmp_output_polar_a_a_-pi/2),reshape(cos(gamma_z_),[1,n_w_max])).^2,reshape(sin(gamma_z_),[1,n_w_max]).^2);
tmp_u1_a_ = mean( -equa_band_dilated_amplitude*tmp_z1_aw__./max(1e-12,tmp_w1_aw__) , 2 );
tmp_u1_ab__ = bsxfun(@times,tmp_u1_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
disp(sprintf(' %% tmp_u_ab__ vs tmp_u1_ab__: %0.16f',fnorm(tmp_u_ab__ - tmp_u1_ab__)/fnorm(tmp_u_ab__))); %<-- should be order equa_band_dilated_amplitude. ;
nb=round((n_point_b-1)/8);
fnorm(tmp_u_ab__(:,1+nb) - tmp_u1_a_)/fnorm(tmp_u_ab__(:,1+nb));
disp(sprintf(' %% tmp_u_ab__(:,1+nb) vs tmp_u1_a_: %0.16f',fnorm(tmp_u_ab__(:,1+nb) - tmp_u1_a_)/fnorm(tmp_u_ab__(:,1+nb)))); %<-- should be order equa_band_dilated_amplitude. ;

% Gradshteyn & Ryzhik: 3.647 p402: ;
% b0=rand();tmp_z_ = cos(gamma_z_).^2./(sin(gamma_z_).^2 + b0^2*cos(gamma_z_).^2); mean(tmp_z_),; 1/b0/(1+b0),;
% b0=rand();tmp_z_ = 1./(sin(gamma_z_).^2 + b0^2*cos(gamma_z_).^2); mean(tmp_z_),; 1/b0,;
% b0=rand();tmp_z_ = (b0^2*cos(gamma_z_).^2 - sin(gamma_z_).^2)./(b0^2*cos(gamma_z_).^2 + sin(gamma_z_).^2); mean(tmp_z_),; (b0^2+1)/b0/(1+b0) - 1/b0,;
% b0=rand();tmp_z_ = (b0^2*cos(gamma_z_).^2 - sin(gamma_z_).^2)./(b0^2*cos(gamma_z_).^2 + sin(gamma_z_).^2); mean(tmp_z_),; (abs(b0)-1)/(abs(b0)+1),;

tmp_output_polar_a_a_ = transpose(linspace(-pi/2,+pi/2,n_point_a)); tmp_ca_ = cos(tmp_output_polar_a_a_-pi/2);
tmp_u2_a_ = -equa_band_dilated_amplitude*(abs(tmp_ca_)-1)./(abs(tmp_ca_)+1);
tmp_u2_ab__ = bsxfun(@times,tmp_u2_a_,reshape(sin(2*point_output_azimu_b_b_),[1,n_point_b]));
disp(sprintf(' %% tmp_u_ab__ vs tmp_u2_ab__: %0.16f',fnorm(tmp_u_ab__ - tmp_u2_ab__)/fnorm(tmp_u_ab__))); %<-- should be order equa_band_dilated_amplitude. ;
nb=round((n_point_b-1)/8);
fnorm(tmp_u_ab__(:,1+nb) - tmp_u2_a_)/fnorm(tmp_u_ab__(:,1+nb));
disp(sprintf(' %% tmp_u_ab__(:,1+nb) vs tmp_u2_a_: %0.16f',fnorm(tmp_u_ab__(:,1+nb) - tmp_u2_a_)/fnorm(tmp_u_ab__(:,1+nb)))); %<-- should be order equa_band_dilated_amplitude. ;











