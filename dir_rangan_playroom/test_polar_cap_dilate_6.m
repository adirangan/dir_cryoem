%%%%%%%%;
% tests polar cap dilation. ;
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
dir_manuscript_jpg = sprintf('%s/dir_M3d_shape_latitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;
flag_replot=1;

%%%%%%%%;
% First set up collection of templates. ;
%%%%%%%%;
flag_verbose = 0; flag_disp=1; nf=0;
n_k_p_r = 1; k_p_r_max = 1.0; k_p_r_ = k_p_r_max; weight_3d_k_p_r_ = 1; l_max_ = 0; a_k_Y_ = 0;
viewing_k_eq_d = 0.5*k_p_r_max/(2*pi);;
template_k_eq_d = -1;
n_w_max = 64; n_w_0in_ = n_w_max;
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

polar_cap_theta = pi/4; polar_cap_dilated_amplitude = 0.1;
index_polar_cap_ = efind(viewing_polar_a_S_<=polar_cap_theta);
n_viewing_polar_cap = numel(index_polar_cap_);
viewing_polar_cap_azimu_b_S_ = viewing_azimu_b_S_(1+index_polar_cap_);
%%%%;
viewing_polar_cap_polar_a_S_ = viewing_polar_a_S_(1+index_polar_cap_);
[ ...
 template_polar_cap_k_c_0_wS__ ...
,template_polar_cap_k_c_1_wS__ ...
,template_polar_cap_k_c_2_wS__ ...
,template_polar_cap_azimu_b_wS__ ...
,template_polar_cap_polar_a_wS__ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,n_viewing_polar_cap ...
,viewing_polar_cap_azimu_b_S_ ...
,viewing_polar_cap_polar_a_S_ ...
,n_w_max ...
);
%%%%;
viewing_polar_cap_dilated_polar_a_S_ = viewing_polar_cap_polar_a_S_ + polar_cap_dilated_amplitude*sin(2*viewing_polar_cap_polar_a_S_);
[ ...
 template_polar_cap_dilated_k_c_0_wS__ ...
,template_polar_cap_dilated_k_c_1_wS__ ...
,template_polar_cap_dilated_k_c_2_wS__ ...
,template_polar_cap_dilated_azimu_b_wS__ ...
,template_polar_cap_dilated_polar_a_wS__ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,n_viewing_polar_cap ...
,viewing_polar_cap_azimu_b_S_ ...
,viewing_polar_cap_dilated_polar_a_S_ ...
,n_w_max ...
);
%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Illustrate polar_cap and polar_cap_dilated. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
%%%%;
subplot(1,2,1);
hold on;
for nS=0:n_viewing_polar_cap-1;
surfline_0( ...
 template_polar_cap_k_c_0_wS__(:,1+nS) ...
,template_polar_cap_k_c_1_wS__(:,1+nS) ...
,template_polar_cap_k_c_2_wS__(:,1+nS) ...
,template_polar_cap_polar_a_wS__(:,1+nS) ...
);
end;%for nS=0:n_viewing_polar_cap-1;
xlim([-1,+1]); ylim([-1,+1]); zlim([-1,+1]); axis equal; axis vis3d;
%%%%;
subplot(1,2,2);
hold on;
for nS=0:n_viewing_polar_cap-1;
surfline_0( ...
 template_polar_cap_dilated_k_c_0_wS__(:,1+nS) ...
,template_polar_cap_dilated_k_c_1_wS__(:,1+nS) ...
,template_polar_cap_dilated_k_c_2_wS__(:,1+nS) ...
,template_polar_cap_dilated_polar_a_wS__(:,1+nS) ...
);
end;%for nS=0:n_viewing_polar_cap_dilated-1;
xlim([-1,+1]); ylim([-1,+1]); zlim([-1,+1]); axis equal; axis vis3d;
%%%%;
end;%if flag_disp;

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
% Set parameters for polar cap dilation. ;
% We assume that each template is dilated (away from the pole) by g_dilation, ;
% corresponding to an overall mapping of viewing_polar_a via f_dilation. ;
% Note that, due to our strategy for constructing templates, ;
% the polar cap dilation corresponds to incrementing the viewing_polar_a by g_dilation, ;
% meaning that the template itself will remain indexed by gamma_z.
% This (canonical) indexing via gamma_z will typically result in a shift between the ;
% original gamma_z associated with a particular point (in a particular template) ;
% and the final gamma_z associated with that same point (in the predilated version of that template). ;
%%%%%%%%;
polar_cap_dilated_amplitude = 0.05;
g_dilation = @(point_pole_predilated_polar_a) polar_cap_dilated_amplitude*sin(2*point_pole_predilated_polar_a); %<-- approximation holds well for first nontrivial mode. ;
f_dilation = @(point_pole_predilated_polar_a) point_pole_predilated_polar_a + g_dilation(point_pole_predilated_polar_a);
%%%%%%%%;
% point_output is indexed by polar_a. ;
% point_pole is indexed by n_w_max. ;
%%%%%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% plot diagram of polar dilation. ;
% Caption: ;
% In this figure we illustrate some properties of the polar-dilation described in the text. ;
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
% Note that the colored dots form a ring; this ring passes through the point $\vk$ itself when the viewing-angle $\tau$ lies on the equator; in this case the viewing-angle $\tau^{\prime}$ is equal to $\tau$ (a special-case of the dilation preserving the location of $\vk$). ;
% The colored ring of dots in the right-subplot indicates the locations $\{\vk^{\prime}\}$ which contribute to the value of the volume at $\vk$ after dilation of all the templates.
%%%%%%%%;
tmp_polar_cap_dilated_amplitude = 0.15;
tmp_g_dilation = @(point_pole_predilated_polar_a) tmp_polar_cap_dilated_amplitude*sin(2*point_pole_predilated_polar_a);
tmp_f_dilation = @(point_pole_predilated_polar_a) point_pole_predilated_polar_a + tmp_g_dilation(point_pole_predilated_polar_a);
figure(1+nf);nf=nf+1;clf;figmed;
c_hsv__ = colormap('hsv'); n_c_hsv = size(c_hsv__,1);
sym_point = 'o'; sym_pole = 's';
markersize_use = 8; markersize_dot = 4; linewidth_use = 1.5; linewidth_mrk = 2.0;
for np=0:1; subplot(1,2,1+np); plot_sphere_grid_0; hold on; end;
npole_use = floor(0.92*n_w_max);
n_point = 1+64;
npoint = 15; %<-- not too close to the pole. ;
point_output_azimu_b = 0;%point_output_azimu_b = 2*pi*npoint/max(1,n_point);
point_output_polar_a = -pi/2 + pi/1*npoint/max(1,n_point-1);
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
%%%%%%%%;
for np=0:1; subplot(1,2,1+np); plot3(point_output_k_c_(1+0),point_output_k_c_(1+1),point_output_k_c_(1+2),sym_point,'MarkerFaceColor',0.85*[1,1,1],'MarkerEdgeColor',[0,0,0],'MarkerSize',markersize_use,'LineWidth',linewidth_mrk); end;
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
point_pole_azimu_b_aw__(1+npoint,1+npole) = point_pole_azimu_b;
point_pole_polar_a_aw__(1+npoint,1+npole) = point_pole_polar_a;
[ ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,point_pole_template_azimu_b_w_ ...
,point_pole_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,1 ...
,point_pole_azimu_b ...
,point_pole_polar_a ...
,n_w_max ...
);
%%%%%%%%;
% Here we determine the predilated template that is associated with the original template mapped to point_output. ;
%%%%%%%%;
tmp_f_error = @(point_pole_predilated_polar_a) abs(tmp_f_dilation(point_pole_predilated_polar_a) - point_pole_polar_a).^2;
point_pole_predilated_polar_a = fminsearch(tmp_f_error,point_pole_polar_a);
tmp_point_pole_polar_a = point_pole_predilated_polar_a + polar_cap_dilated_amplitude*sin(2*point_pole_predilated_polar_a);
if (flag_verbose>1); disp(sprintf(' %% npole %d/%d, fnorm(point_pole_polar_a-tmp_point_pole_polar_a): %0.16f',npole,n_w_max,fnorm(point_pole_polar_a-tmp_point_pole_polar_a))); end;
point_pole_predilated_k_c_ = [ ...
  cos(point_pole_azimu_b)*sin(point_pole_predilated_polar_a) ...
; sin(point_pole_azimu_b)*sin(point_pole_predilated_polar_a) ...
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
 flag_verbose ...
,1 ...
,point_pole_azimu_b ...
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
%%%%%%%%;
% Now we determine the gamma_z (denoted gammax) within predilated template that maps to point_output. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_ = [ ...
 +cos(point_pole_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_predilated_polar_a)*cos(point_pole_gx) ...
];
point_pole_predilated_template_gammax_k_c_0_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+0);
point_pole_predilated_template_gammax_k_c_1_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+1);
point_pole_predilated_template_gammax_k_c_2_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+2);
point_pole_predilated_template_gammax_k_r01 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2);
point_pole_predilated_template_gammax_k_r012 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2 + point_pole_predilated_template_gammax_k_c_(1+2).^2);
point_pole_predilated_template_gammax_azimu_b = atan2(point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+0));
point_pole_predilated_template_gammax_polar_a = atan2(point_pole_predilated_template_gammax_k_r01,point_pole_predilated_template_gammax_k_c_(1+2));
point_pole_predilated_template_gammax_azimu_b_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_azimu_b;
point_pole_predilated_template_gammax_polar_a_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_polar_a;
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
fname_fig_pre = sprintf('%s/M3d_shape_latitudinal_perturbation_diagram_FIGA',dir_manuscript_jpg);
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

%error('stopping early');

flag_disp=1;
%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;fig80s; c_hsv__ = colormap('hsv'); n_c_hsv = size(c_hsv__,1);
markersize_big = 8; markersize_use = 6; markersize_dot = 4; linewidth_use = 0.5;
plot_sphere_grid_0;
hold on;
end;%if flag_disp;
%%%%%%%%;
n_point = 1+64;
point_output_azimu_b_a_ = zeros(n_point,1);
point_output_polar_a_a_ = zeros(n_point,1);
point_pole_azimu_b_aw__ = zeros(n_point,n_w_max);
point_pole_polar_a_aw__ = zeros(n_point,n_w_max);
point_pole_predilated_template_gammax_k_c_0_aw__ = zeros(n_point,n_w_max);
point_pole_predilated_template_gammax_k_c_1_aw__ = zeros(n_point,n_w_max);
point_pole_predilated_template_gammax_k_c_2_aw__ = zeros(n_point,n_w_max);
point_pole_predilated_template_gammax_azimu_b_aw__ = zeros(n_point,n_w_max);
point_pole_predilated_template_gammax_polar_a_aw__ = zeros(n_point,n_w_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npoint=0:n_point-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
point_output_azimu_b = 0;%point_output_azimu_b = 2*pi*npoint/max(1,n_point);
point_output_polar_a = -pi/2 + pi/1*npoint/max(1,n_point-1);
point_output_gamma_z = 0;
point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_output_azimu_b_a_(1+npoint) = point_output_azimu_b;
point_output_polar_a_a_(1+npoint) = point_output_polar_a;
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
point_pole_azimu_b_aw__(1+npoint,1+npole) = point_pole_azimu_b;
point_pole_polar_a_aw__(1+npoint,1+npole) = point_pole_polar_a;
[ ...
 point_pole_template_k_c_0_w_ ...
,point_pole_template_k_c_1_w_ ...
,point_pole_template_k_c_2_w_ ...
,point_pole_template_azimu_b_w_ ...
,point_pole_template_polar_a_w_ ...
] = ...
get_template_single_ring_k_c_0( ...
 flag_verbose ...
,1 ...
,point_pole_azimu_b ...
,point_pole_polar_a ...
,n_w_max ...
);
%%%%%%%%;
% Here we determine the predilated template that is associated with the original template mapped to point_output. ;
%%%%%%%%;
tmp_f_error = @(point_pole_predilated_polar_a) abs(f_dilation(point_pole_predilated_polar_a) - point_pole_polar_a).^2;
point_pole_predilated_polar_a = fminsearch(tmp_f_error,point_pole_polar_a);
tmp_point_pole_polar_a = point_pole_predilated_polar_a + polar_cap_dilated_amplitude*sin(2*point_pole_predilated_polar_a);
if (flag_verbose>1); disp(sprintf(' %% npole %d/%d, fnorm(point_pole_polar_a-tmp_point_pole_polar_a): %0.16f',npole,n_w_max,fnorm(point_pole_polar_a-tmp_point_pole_polar_a))); end;
point_pole_predilated_k_c_ = [ ...
  cos(point_pole_azimu_b)*sin(point_pole_predilated_polar_a) ...
; sin(point_pole_azimu_b)*sin(point_pole_predilated_polar_a) ...
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
 flag_verbose ...
,1 ...
,point_pole_azimu_b ...
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
%%%%%%%%;
% Now we determine the gamma_z (denoted gammax) within predilated template that maps to point_output. ;
%%%%%%%%;
point_pole_predilated_template_gammax_k_c_ = [ ...
 +cos(point_pole_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) - sin(point_pole_azimu_b)*sin(point_pole_gx) ...
;+sin(point_pole_azimu_b)*cos(point_pole_predilated_polar_a)*cos(point_pole_gx) + cos(point_pole_azimu_b)*sin(point_pole_gx) ...
;-sin(point_pole_predilated_polar_a)*cos(point_pole_gx) ...
];
point_pole_predilated_template_gammax_k_c_0_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+0);
point_pole_predilated_template_gammax_k_c_1_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+1);
point_pole_predilated_template_gammax_k_c_2_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_k_c_(1+2);
point_pole_predilated_template_gammax_k_r01 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2);
point_pole_predilated_template_gammax_k_r012 = sqrt(point_pole_predilated_template_gammax_k_c_(1+0).^2 + point_pole_predilated_template_gammax_k_c_(1+1).^2 + point_pole_predilated_template_gammax_k_c_(1+2).^2);
point_pole_predilated_template_gammax_azimu_b = atan2(point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+0));
point_pole_predilated_template_gammax_polar_a = atan2(point_pole_predilated_template_gammax_k_r01,point_pole_predilated_template_gammax_k_c_(1+2));
point_pole_predilated_template_gammax_azimu_b_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_azimu_b;
point_pole_predilated_template_gammax_polar_a_aw__(1+npoint,1+npole) = point_pole_predilated_template_gammax_polar_a;
%%%%%%%%;
if flag_disp;
%plot3(point_pole_k_c__(1+0,:),point_pole_k_c__(1+1,:),point_pole_k_c__(1+2,:),'k.','MarkerSize',markersize_dot);
%plot3(point_pole_k_c_(1+0),point_pole_k_c_(1+1),point_pole_k_c_(1+2),'co','MarkerFaceColor',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_use);
%plot3(point_pole_predilated_k_c_(1+0),point_pole_predilated_k_c_(1+1),point_pole_predilated_k_c_(1+2),'yo','MarkerFaceColor',c_hsv__(1+nc_hsv,:).^0.5,'MarkerSize',markersize_use);
%plot3(point_pole_template_k_c_0_w_,point_pole_template_k_c_1_w_,point_pole_template_k_c_2_w_,'g-','LineWidth',linewidth_use);
%plot3(point_pole_template_k_c_0_w_(1+0),point_pole_template_k_c_1_w_(1+0),point_pole_template_k_c_2_w_(1+0),'ko','MarkerSize',markersize_use);
%plot3(point_pole_template_gamma0_k_c_(1+0),point_pole_template_gamma0_k_c_(1+1),point_pole_template_gamma0_k_c_(1+2),'kx','MarkerSize',markersize_use);
plot3(point_pole_template_gammax_k_c_(1+0),point_pole_template_gammax_k_c_(1+1),point_pole_template_gammax_k_c_(1+2),'ko','MarkerSize',markersize_big);
plot3(point_pole_predilated_template_gammax_k_c_(1+0),point_pole_predilated_template_gammax_k_c_(1+1),point_pole_predilated_template_gammax_k_c_(1+2),'.','Color',c_hsv__(1+nc_hsv,:),'MarkerSize',markersize_dot);
end;%if flag_disp;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for npoint=0:n_point-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
xlim([-1,+1]); ylim([-1,+1]); zlim([-1,+1]); axis equal; axis vis3d;
xlabel('x'); ylabel('y'); zlabel('z');
hold off;
end;%if flag_disp;
%%%%%%%%;

if (flag_verbose> 0);
disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_0_aw__ - cos(point_pole_predilated_template_gammax_azimu_b_aw__).*sin(point_pole_predilated_template_gammax_polar_a_aw__))));
disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_1_aw__ - sin(point_pole_predilated_template_gammax_azimu_b_aw__).*sin(point_pole_predilated_template_gammax_polar_a_aw__))));
disp(sprintf(' %% tmp_error: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_2_aw__ - cos(point_pole_predilated_template_gammax_polar_a_aw__))));
end;%if (flag_verbose> 0);

point_pole_predilated_template_gammax_k_c_0_avg_a_ = mean(point_pole_predilated_template_gammax_k_c_0_aw__,2);
point_pole_predilated_template_gammax_k_c_1_avg_a_ = mean(point_pole_predilated_template_gammax_k_c_1_aw__,2);
point_pole_predilated_template_gammax_k_c_2_avg_a_ = mean(point_pole_predilated_template_gammax_k_c_2_aw__,2);
point_pole_predilated_template_gammax_polar_a_avg_a_ = mean(point_pole_predilated_template_gammax_polar_a_aw__,2);
point_pole_predilated_template_gammax_polar_a_std_a_ = std(point_pole_predilated_template_gammax_polar_a_aw__,1,2);
point_pole_predilated_template_gammax_k_c_r_aw__ = ...
sqrt( ...
 bsxfun(@minus,point_pole_predilated_template_gammax_k_c_0_aw__,mean(point_pole_predilated_template_gammax_k_c_0_aw__,2)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_1_aw__,mean(point_pole_predilated_template_gammax_k_c_1_aw__,2)).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_2_aw__,mean(point_pole_predilated_template_gammax_k_c_2_aw__,2)).^2 ...
);
point_pole_predilated_template_gammax_k_c_r_avg1_a_ = sqrt(mean(point_pole_predilated_template_gammax_k_c_r_aw__.^2,2)); %<-- take 2-norm, rather than mean norm. ;
point_pole_predilated_template_gammax_k_c_r_avg0_a_ = mean(point_pole_predilated_template_gammax_k_c_r_aw__,2); %<-- take mean-norm, rather than 2-norm. ;
disp(sprintf(' %% vs: %0.16f',fnorm(point_pole_predilated_template_gammax_k_c_r_avg1_a_ - point_pole_predilated_template_gammax_k_c_r_avg0_a_)./fnorm(point_pole_predilated_template_gammax_k_c_r_avg1_a_)));

flag_disp=1;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Demonstrate match. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot(1,3,1);
hold on;
plot(point_output_polar_a_a_+pi/2,point_pole_predilated_template_gammax_polar_a_avg_a_-point_output_polar_a_a_-pi/2,'ro-');
plot(point_output_polar_a_a_+pi/2,0.5*g_dilation(point_output_polar_a_a_+pi/2),'gx-');
hold off;
xlim([0,pi]); xlabel('point_output_polar_a_a_+pi/2','Interpreter','none'); ylabel('amplitude'); title('avg polar_a displacement amplitude','Interpreter','none');
%%%%;
subplot(1,3,2);
hold on;
plot(point_output_polar_a_a_+pi/2,+point_pole_predilated_template_gammax_polar_a_std_a_,'ro-');
plot(point_output_polar_a_a_+pi/2,-point_pole_predilated_template_gammax_polar_a_std_a_,'ro-');
plot(point_output_polar_a_a_+pi/2,(0.5/sqrt(2))*g_dilation(point_output_polar_a_a_+pi/2),'gx-');
hold off;
xlim([0,pi]); xlabel('point_output_polar_a_a_+pi/2','Interpreter','none'); ylabel('amplitude'); title('std polar_a displacement amplitude','Interpreter','none');
%%%%;
subplot(1,3,3);
hold on;
plot(point_output_polar_a_a_+pi/2,+point_pole_predilated_template_gammax_k_c_r_avg0_a_,'mo-');
plot(point_output_polar_a_a_+pi/2,-point_pole_predilated_template_gammax_k_c_r_avg0_a_,'mo-');
plot(point_output_polar_a_a_+pi/2,0.5*g_dilation(point_output_polar_a_a_+pi/2),'cx-');
hold off;
xlim([0,pi]); xlabel('point_output_polar_a_a_+pi/2','Interpreter','none'); ylabel('amplitude'); title('avg polar_a displacement radius','Interpreter','none');
%%%%;
fname_fig_pre = sprintf('%s/M3d_shape_latitudinal_perturbation_avg_and_std_FIGA',dir_manuscript_jpg);
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
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Now estimate:
% point_pole_predilated_template_gammax_polar_a_avg_a_, ;
% point_pole_predilated_template_gammax_polar_a_std_a_, ;
% point_pole_predilated_template_gammax_k_c_r_avg0_a_, ;
% for a variety of viewing-angle distributions. ;
%%%%%%%%;
% We assume the viewing-angle distribution is concentrated around the poles: ;
% rho(polar_a) \propto exp( -sin(polar_a)^2 / (2*sigma_a^2) ) ;
%%%%%%%%;
lsigma_a_ = transpose(-4:0.125:0); n_lsigma_a = numel(lsigma_a_);
point_pole_predilated_template_gammax_polar_a_avg_as__ = zeros(n_point,n_lsigma_a);
point_pole_predilated_template_gammax_polar_a_std_as__ = zeros(n_point,n_lsigma_a);
point_pole_predilated_template_gammax_k_c_r_avg0_as__ = zeros(n_point,n_lsigma_a);
for nlsigma_a=0:n_lsigma_a-1;
lsigma_a = lsigma_a_(1+nlsigma_a);
sigma_a = exp(lsigma_a);
%%%%;
tmp_weight_aw__ = zeros(n_point,n_w_max);
for npoint=0:n_point-1;
tmp_weight_w_ = zeros(n_w_max,1);
tmp_s_ = sin(point_pole_polar_a_aw__(1+npoint,:));
tmp_s_min = min(tmp_s_);
tmp_nlweight_w_ = +(tmp_s_-tmp_s_min).^2/(2*sigma_a^2); %<-- -log(weight_w_). ;
tmp_w_min = min(tmp_nlweight_w_);
tmp_weight_w_ = exp(-(tmp_nlweight_w_-tmp_w_min)); %<-- at least one positive value. ;
tmp_weight_w_ = tmp_weight_w_/max(1e-12,sum(tmp_weight_w_));
tmp_weight_aw__(1+npoint,:) = tmp_weight_w_;
end;%for npoint=0:n_point-1;
%%%%;
tmp_point_pole_predilated_template_gammax_k_c_0_avg_a_ = sum(point_pole_predilated_template_gammax_k_c_0_aw__.*tmp_weight_aw__,2);
tmp_point_pole_predilated_template_gammax_k_c_1_avg_a_ = sum(point_pole_predilated_template_gammax_k_c_1_aw__.*tmp_weight_aw__,2);
tmp_point_pole_predilated_template_gammax_k_c_2_avg_a_ = sum(point_pole_predilated_template_gammax_k_c_2_aw__.*tmp_weight_aw__,2);
tmp_point_pole_predilated_template_gammax_polar_a_avg_a_ = sum(point_pole_predilated_template_gammax_polar_a_aw__.*tmp_weight_aw__,2);
tmp_point_pole_predilated_template_gammax_polar_a_std_a_ = sqrt(sum(bsxfun(@minus,point_pole_predilated_template_gammax_polar_a_aw__,tmp_point_pole_predilated_template_gammax_polar_a_avg_a_).^2.*tmp_weight_aw__,2));
tmp_point_pole_predilated_template_gammax_k_c_r_aw__ = ...
sqrt( ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_0_aw__,tmp_point_pole_predilated_template_gammax_k_c_0_avg_a_).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_1_aw__,tmp_point_pole_predilated_template_gammax_k_c_1_avg_a_).^2 ...
+bsxfun(@minus,point_pole_predilated_template_gammax_k_c_2_aw__,tmp_point_pole_predilated_template_gammax_k_c_2_avg_a_).^2 ...
);
tmp_point_pole_predilated_template_gammax_k_c_r_avg0_a_ = sum(tmp_point_pole_predilated_template_gammax_k_c_r_aw__.*tmp_weight_aw__,2); %<-- take mean-norm, rather than 2-norm. ;
point_pole_predilated_template_gammax_polar_a_avg_as__(:,1+nlsigma_a) = tmp_point_pole_predilated_template_gammax_polar_a_avg_a_;
point_pole_predilated_template_gammax_polar_a_std_as__(:,1+nlsigma_a) = tmp_point_pole_predilated_template_gammax_polar_a_std_a_;
point_pole_predilated_template_gammax_k_c_r_avg0_as__(:,1+nlsigma_a) = tmp_point_pole_predilated_template_gammax_k_c_r_avg0_a_;
end;%for nlsigma_a=0:n_lsigma_a-1;

%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed; fig80s;
subplot(1,3,1);
plot(point_output_polar_a_a_+pi/2,abs((0.5/sqrt(2))*g_dilation(point_output_polar_a_a_+pi/2)),'m.-','LineWidth',3);
xlim([0,pi]); xlabel('point_output_polar_a_a_+pi/2','Interpreter','none'); ylabel('amplitude');
set(gca,'YTick',[],'XTick',pi*[0,0.5,1],'XTickLabel',{'0','\pi/2','\pi'});
title('uniform distribution');
subplot(1,3,[2,3]);
imagesc(point_pole_predilated_template_gammax_polar_a_std_as__);
set(gca,'Ytick',[]);
set(gca,'YTick',[1,(n_point+1)/2,n_point],'YTickLabel',{'0','\pi/2','\pi'});
ylabel('point_output_polar_a_a+pi/2','Interpreter','none');
set(gca,'Xtick',[1:4:n_lsigma_a],'XTickLabel',num2str(exp(lsigma_a_(1:4:end)),'%0.2f'));
xlabel('polar cap std sigma (units of sin(polar_a))','Interpreter','none');
title('gaussian-distributed around polar cap with std sigma');
sgtitle('std polar_a displacement amplitude','Interpreter','none');
%%%%;    
fname_fig_pre = sprintf('%s/M3d_shape_latitudinal_perturbation_std_vs_sigma_FIGA',dir_manuscript_jpg);
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
