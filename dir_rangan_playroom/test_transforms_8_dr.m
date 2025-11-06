%%%%%%%%;
% Sets up a simple volume (in spherical-harmonic-coordinates), ;
% then tests the various template-generators. ;
% These calculations consider anisotropic CTF. ;
%%%%%%%%;

str_thisfunction = 'test_transforms_8_dr';

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_verbose = 1; %<-- verbosity level. ;
flag_disp=1; nf=0; %<-- display level. ;
flag_replot = 1;

k_int = 16; %<-- highest frequency (2*pi*k_p_r_max), watch out for aliasing! ;
k_eq_d_double = 1.0; %<-- prefactor for k_eq_d, determines density of sampling in frequency-space. ;
template_k_eq_d_double = 0.5; %<-- prefactor for template_k_eq_d, determines density of viewing-angles on the sphere. ;
n_w_int = 1.0; %<-- prefactor for n_w_max, determines the number of distinct angles (i.e., n_gamma_z) used in frequency-space 2d-polar-grid. ;
n_x_c = max(64,2*k_int); %<-- the number of 'pixels' on a side in the real-space cartesian-grid. ;

dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom',string_root);
dir_jpg = sprintf('%s/dir_test_transforms_jpg',dir_base);
if (~exist(dir_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;

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
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = half_diameter_x_c;
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
weight_xxx_c = (2*x_p_r_max/n_x_c)^3;
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 12;
fontsize_use = 16;
subplot(1,2,1);
plot(x_c_0___(:,:,1+0),x_c_1___(:,:,1+0),'k.','MarkerSize',markersize_use);
xlabel('x_c_0','Interpreter','none');
ylabel('x_c_1','Interpreter','none');
axis equal; axisnotick;
title('2d real-space cartesian-grid');
set(gca,'FontSize',fontsize_use);
subplot(1,2,2);
plot3(x_c_0___(:),x_c_1___(:),x_c_2___(:),'k.','MarkerSize',markersize_use);
xlabel('x_c_0','Interpreter','none');
ylabel('x_c_1','Interpreter','none');
zlabel('x_c_2','Interpreter','none');
axis equal; axis vis3d; axisnotick3d;
title('3d real-space cartesian-grid');
set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/test_transforms_x_c_x___FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
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
% For eventual reconstruction we generate quadrature on shell. ;
%%%%%%%%;
qref_k_eq_d = k_eq_d/max(1e-12,k_p_r_max);
[ ...
 qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,qref_k_eq_d ...
,str_L ...
,flag_uniform_over_polar_a ...
);
tmp_index_ = [n_qk_csum_(end-1):n_qk_csum_(end)-1];
fnorm_disp(flag_verbose,'k_p_azimu_b_qk_(1+tmp_index_)',k_p_azimu_b_qk_(1+tmp_index_),'qref_azimu_b_shell_',qref_azimu_b_shell_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'k_p_polar_a_qk_(1+tmp_index_)',k_p_polar_a_qk_(1+tmp_index_),'qref_polar_a_shell_',qref_polar_a_shell_,' %%<-- should be zero');
%%%%%%%%;
% Visualization. ;
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
markersize_use = 4;
fontsize_use = 12;
p_row = 2; p_col = 4; n_plot = p_row*p_col;
for nplot=0:n_plot-1;
nk_p_r = max(0,min(n_k_p_r-1,round(n_k_p_r*nplot/n_plot)));
tmp_index_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
subplot(p_row,p_col,1+nplot);
plot3(k_c_0_qk_(1+tmp_index_),k_c_1_qk_(1+tmp_index_),k_c_2_qk_(1+tmp_index_),'k.','MarkerSize',markersize_use);
xlabel('k_c_0','Interpreter','none');
ylabel('k_c_1','Interpreter','none');
zlabel('k_c_2','Interpreter','none');
xlim(1.15*k_p_r_max*[-1,+1]);
ylim(1.15*k_p_r_max*[-1,+1]);
zlim(1.15*k_p_r_max*[-1,+1]);
axis vis3d; axisnotick3d;
set(gca,'FontSize',fontsize_use);
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
end;%for nplot=0:n_plot-1;
sgtitle('frequency-space spherical-grid');
fname_fig_pre = sprintf('%s/test_transforms_k_c_x_qk_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now define two functions on the sphere, ;
% each a sum of a few plane-waves. ;
%%%%%%%%;
delta_a_c_3s__ = [  ...
 +1.5 , -0.5 ...
;-0.5 , -1.5 ...
;+0.3 , +2.0 ...
] / 2 / k_p_r_max ;
delta_b_c_3s__ = [  ...
 -0.5 , +0.8 ...
;-1.0 , +0.2 ...
;+1.2 , -0.7 ...
] / 2 / k_p_r_max ;
n_source = size(delta_a_c_3s__,2);
a_k_p_form_ = zeros(n_qk,1);
b_k_p_form_ = zeros(n_qk,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_a_c_(1+0) + k_c_1_qk_*delta_a_c_(1+1) + k_c_2_qk_*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_qk_*delta_b_c_(1+0) + k_c_1_qk_*delta_b_c_(1+1) + k_c_2_qk_*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%%%%%;
% Define frequency-space cartesian-grid. ;
%%%%%%%%;
n_k_c = n_x_c + 2; %<-- just to check dimensions. ;
half_diameter_k_c = k_p_r_max;
diameter_k_c = 2.0d0*half_diameter_k_c;
k_p_r_max = half_diameter_k_c;
k_c_0_ = linspace(-k_p_r_max,+k_p_r_max,n_k_c);
k_c_1_ = linspace(-k_p_r_max,+k_p_r_max,n_k_c);
k_c_2_ = linspace(-k_p_r_max,+k_p_r_max,n_k_c);
[k_c_0___,k_c_1___,k_c_2___] = ndgrid(k_c_0_,k_c_1_,k_c_2_); n_kkk_c = n_k_c^3;
weight_kkk_c = (2*k_p_r_max/n_k_c)^3;
%%%%%%%%;
a_k_c_form___ = zeros(n_k_c,n_k_c,n_k_c);
b_k_c_form___ = zeros(n_k_c,n_k_c,n_k_c);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
a_k_c_form___ = a_k_c_form___ + exp(+i*2*pi*(k_c_0___*delta_a_c_(1+0) + k_c_1___*delta_a_c_(1+1) + k_c_2___*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
b_k_c_form___ = b_k_c_form___ + exp(+i*2*pi*(k_c_0___*delta_b_c_(1+0) + k_c_1___*delta_b_c_(1+1) + k_c_2___*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
%%%%%%%%;

%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 12;
p_row = 2; p_col = 4; np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1; cla;
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_80s();
a_k_c_form_lim_ = prctile(mean(reshape(real(a_k_c_form___),[n_k_c,n_k_c,n_k_c]),3),[ 5,95],'all');
tmp_parameter.vlim_ = 1.0*a_k_c_form_lim_;
tmp_parameter.percent_threshold_ = [ 1,99];
isosurface_f_x_u_1(tmp_parameter,a_k_c_form___);
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d;
title('real(a_k_c_form___)','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
%%%%;
nk_p_r = max(0,min(n_k_p_r-1,floor(1.0*n_k_p_r)));
k_p_r = k_p_r_(1+nk_p_r);
n_qk_csum_pre = n_qk_csum_(1+nk_p_r+0);
n_qk_csum_pos = n_qk_csum_(1+nk_p_r+1);
assert(qref_n_shell==n_qk_csum_pos-n_qk_csum_pre);
a_k_p_sub_fin_ = a_k_p_form_(1+n_qk_csum_pre+[0:qref_n_shell-1]);
a_k_p_sub_fin_lim_ = prctile(abs(a_k_p_sub_fin_),95)*[-1,+1];
a_k_p_ori_sub_fin_lim_ = a_k_p_sub_fin_lim_;
%%%%;
for np=1:n_plot-1;
subplot(p_row,p_col,1+np);
nk_p_r = max(0,min(n_k_p_r-1,floor((np/max(1,n_plot-1))*n_k_p_r)));
k_p_r = k_p_r_(1+nk_p_r);
n_qk_csum_pre = n_qk_csum_(1+nk_p_r+0);
n_qk_csum_pos = n_qk_csum_(1+nk_p_r+1);
assert(qref_n_shell==n_qk_csum_pos-n_qk_csum_pre);
a_k_p_sub_ = a_k_p_form_(1+n_qk_csum_pre+[0:qref_n_shell-1]);
a_k_p_sub_lim_ = prctile(abs(a_k_p_sub_),95)*[-1,+1];
%%%%;
hold on;
imagesc_polar_a_azimu_b_0(qref_polar_a_shell_,qref_azimu_b_shell_,real(a_k_p_sub_),a_k_p_ori_sub_fin_lim_,colormap_80s,0,k_p_r);
n_contour = 16;
tmp_parameter = struct('type','parameter');
tmp_parameter.flag_percent_use=0;
tmp_parameter.vlim_ = 1.0*a_k_p_ori_sub_fin_lim_;
tmp_parameter.vval_ = transpose(linspace(min(tmp_parameter.vlim_),max(tmp_parameter.vlim_),n_contour));
tmp_parameter.flag_k_c_interp = 1;
tmp_parameter.k_p_r_max_use = 1.0*k_p_r_max;
imagesc_S_k_p_3d_2( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
);
hold off;
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d; axis equal;
title(sprintf('nk_p_r %d k_p_r %0.2f',nk_p_r,k_p_r),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
drawnow();
end;%for np=1:n_plot-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_transforms_a_k_p_form_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now test k-quadrature on sphere. ;
%%%%%%%%;
I_a_quad = sum(a_k_p_form_.*weight_3d_k_p_qk_);
I_b_quad = sum(b_k_p_form_.*weight_3d_k_p_qk_);
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h3d_(tmp_kd)*k_p_r_max^3;
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h3d_(tmp_kd)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'I_a_form',I_a_form,'I_a_quad',I_a_quad,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'I_b_form',I_b_form,'I_b_quad',I_b_quad,' %%<-- should be <1e-6');
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
b_k_p_l2_quad = sum(conj(b_k_p_form_).*b_k_p_form_.*weight_3d_k_p_qk_);
b_k_p_l2_form = 0;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
delta_b_c_0_ = delta_b_c_3s__(:,1+nsource0);
delta_b_c_1_ = delta_b_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_0_ - delta_b_c_1_);
tmp_h3d = 4*pi/3; if abs(tmp_kd)>1e-12; tmp_h3d = h3d_(tmp_kd); end;
b_k_p_l2_form = b_k_p_l2_form + tmp_h3d*k_p_r_max^3;
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
fnorm_disp(flag_verbose,'b_k_p_l2_form',b_k_p_l2_form,'b_k_p_l2_quad',b_k_p_l2_quad,' %%<-- should be <1e-6');
%%%%%%%%;

%%%%%%%%;
% Using the k-quadrature on the sphere we can determine the real-space functions a_x_c_form___ and b_x_c_form___ analytically. ;
%%%%%%%%;
a_x_c_slow___ = zeros(n_x_c,n_x_c,n_x_c);
for nx_0=0:n_x_c-1;
for nx_1=0:n_x_c-1;
for nx_2=0:n_x_c-1;
a_x_c_slow = 0.0;
for nsource0=0:n_source-1;
delta_a_c_0_ = delta_a_c_3s__(:,1+nsource0);
delta_a_c_1_ = -[x_c_0_(1+nx_0);x_c_1_(1+nx_1);x_c_2_(1+nx_2)];
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_0_ - delta_a_c_1_);
tmp_h3d = 4*pi/3; if abs(tmp_kd)>1e-12; tmp_h3d = h3d_(tmp_kd); end;
a_x_c_slow = a_x_c_slow + tmp_h3d*k_p_r_max^3;
end;%for nsource0=0:n_source-1;
a_x_c_slow___(1+nx_0,1+nx_1,1+nx_2) = a_x_c_slow;
end;%for nx_2=0:n_x_c-1;
end;%for nx_1=0:n_x_c-1;
end;%for nx_0=0:n_x_c-1;
a_x_c_l2_quad = sum(abs(a_x_c_slow___).^2,'all')*weight_xxx_c;
disp(sprintf(' %% Note l2-loss: a_x_c_l2_quad %+0.6f vs a_k_p_l2_form %+0.6f',a_x_c_l2_quad,a_k_p_l2_form));
%%%%%%%%;
b_x_c_slow___ = zeros(n_x_c,n_x_c,n_x_c);
for nx_0=0:n_x_c-1;
for nx_1=0:n_x_c-1;
for nx_2=0:n_x_c-1;
b_x_c_slow = 0.0;
for nsource0=0:n_source-1;
delta_b_c_0_ = delta_b_c_3s__(:,1+nsource0);
delta_b_c_1_ = -[x_c_0_(1+nx_0);x_c_1_(1+nx_1);x_c_2_(1+nx_2)];
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_0_ - delta_b_c_1_);
tmp_h3d = 4*pi/3; if abs(tmp_kd)>1e-12; tmp_h3d = h3d_(tmp_kd); end;
b_x_c_slow = b_x_c_slow + tmp_h3d*k_p_r_max^3;
end;%for nsource0=0:n_source-1;
b_x_c_slow___(1+nx_0,1+nx_1,1+nx_2) = b_x_c_slow;
end;%for nx_2=0:n_x_c-1;
end;%for nx_1=0:n_x_c-1;
end;%for nx_0=0:n_x_c-1;
b_x_c_l2_quad = sum(abs(b_x_c_slow___).^2,'all')*weight_xxx_c;
disp(sprintf(' %% Note l2-loss: b_x_c_l2_quad %+0.6f vs b_k_p_l2_form %+0.6f',b_x_c_l2_quad,b_k_p_l2_form));
%%%%%%%%;
a_x_c_form___ = zeros(n_x_c,n_x_c,n_x_c);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c_3s__(:,1+nsource);
tmp_kd___ = ...
  2*pi*k_p_r_max ...
  *sqrt( ...
       + (delta_a_c_(1+0) + x_c_0___).^2 ...
       + (delta_a_c_(1+1) + x_c_1___).^2 ...
       + (delta_a_c_(1+2) + x_c_2___).^2 ...
       ) ...
  ;
tmp_h3d___ = h3d_(tmp_kd___);
tmp_index_ = efind(abs(tmp_kd___)<=1e-12);
tmp_h3d___(1+tmp_index_) = 4*pi/3;
a_x_c_form___ = a_x_c_form___ + tmp_h3d___*k_p_r_max.^3 ;
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'a_x_c_slow___',a_x_c_slow___,'a_x_c_form___',a_x_c_form___,' %%<-- should be 0');
a_x_c_l2_quad = sum(abs(a_x_c_form___).^2,'all')*weight_xxx_c;
disp(sprintf(' %% Note l2-loss: a_x_c_l2_quad %+0.6f vs a_k_p_l2_form %+0.6f',a_x_c_l2_quad,a_k_p_l2_form));
%%%%%%%%;
b_x_c_form___ = zeros(n_x_c,n_x_c,n_x_c);
for nsource=0:n_source-1;
delta_b_c_ = delta_b_c_3s__(:,1+nsource);
tmp_kd___ = ...
  2*pi*k_p_r_max ...
  *sqrt( ...
       + (delta_b_c_(1+0) + x_c_0___).^2 ...
       + (delta_b_c_(1+1) + x_c_1___).^2 ...
       + (delta_b_c_(1+2) + x_c_2___).^2 ...
       ) ...
  ;
tmp_h3d___ = h3d_(tmp_kd___);
tmp_index_ = efind(abs(tmp_kd___)<=1e-12);
tmp_h3d___(1+tmp_index_) = 4*pi/3;
b_x_c_form___ = b_x_c_form___ + tmp_h3d___*k_p_r_max.^3 ;
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'b_x_c_slow___',b_x_c_slow___,'b_x_c_form___',b_x_c_form___,' %%<-- should be 0');
b_x_c_l2_quad = sum(abs(b_x_c_form___).^2,'all')*weight_xxx_c;
disp(sprintf(' %% Note l2-loss: b_x_c_l2_quad %+0.6f vs b_k_p_l2_form %+0.6f',b_x_c_l2_quad,b_k_p_l2_form));
%%%%%%%%;
% Note that, due to the l2-loss, ;
% we can not expect to reconstruct a_k_p_form_ from a_x_c_l2_quad_ on the limited real-space cartesian-grid. ;
% (i.e., the high-frequency-components will be lost). ;
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_quad___ = reshape(xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,a_k_p_form_.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta),[n_x_c,n_x_c,n_x_c]);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_quad___: time %0.6fs',tmp_t));
fnorm_disp(flag_verbose,'a_x_c_form___',a_x_c_form___,'a_x_c_quad___',a_x_c_quad___);
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,a_x_c_form___(:).*weight_xxx_c(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_: time %0.6fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_quad_',a_k_p_quad_,' %%<-- can be large (bandlimited)');
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
b_x_c_quad___ = reshape(xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,b_k_p_form_.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta),[n_x_c,n_x_c,n_x_c]);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: b_x_c_quad___: time %0.6fs',tmp_t));
fnorm_disp(flag_verbose,'b_x_c_form___',b_x_c_form___,'b_x_c_quad___',b_x_c_quad___);
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
b_k_p_quad_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,b_x_c_form___(:).*weight_xxx_c(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: b_k_p_quad_: time %0.6fs',tmp_t));
fnorm_disp(flag_verbose,'b_k_p_form_',b_k_p_form_,'b_k_p_quad_',b_k_p_quad_,' %%<-- can be large (bandlimited)');
%%%%%%%%;

%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
%%%%%%%%;
subplot(1,2,1);
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_beach();
a_x_c_form_lim_ = prctile(mean(reshape(real(a_x_c_form___),[n_x_c,n_x_c,n_x_c]),3),[ 5,95],'all');
tmp_parameter.vlim_ = 1.0*a_x_c_form_lim_;
isosurface_f_x_u_1(tmp_parameter,a_x_c_form___);
xlabel('x0'); ylabel('x1'); zlabel('x2'); axisnotick3d;
title('a_x_c_form___','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
subplot(1,2,2);
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_beach();
b_x_c_form_lim_ = prctile(mean(reshape(real(b_x_c_form___),[n_x_c,n_x_c,n_x_c]),3),[ 5,95]);
tmp_parameter.vlim_ = 1.0*b_x_c_form_lim_;
isosurface_f_x_u_1(tmp_parameter,b_x_c_form___);
xlabel('x0'); ylabel('x1'); zlabel('x2'); axisnotick3d;
title('b_x_c_form___','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
drawnow;
fname_fig_pre = sprintf('%s/test_transforms_a_x_c_form___FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now set up and test polar-quadrature-weights on disk. ;
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
,weight_3d_k_p_r_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
n_gamma_z = n_w_max;
gamma_z_ = transpose(linspace(0,2*pi,n_gamma_z+1)); gamma_z_ = gamma_z_(1:n_gamma_z);
%%%%%%%%;
tmp_S_delta_x_c_ = 0.85*[cos(pi/4);sin(pi/4)]/max(1e-12,k_p_r_max);
tmp_S_phi = +1*pi/5;
tmp_T_delta_x_c_ = 1.35*[cos(-pi/3);sin(-pi/3)]/max(1e-12,k_p_r_max);
tmp_T_phi = -3*pi/7;
tmp_S_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_c_1_wk_*tmp_S_delta_x_c_(1+1)));
tmp_plane_S_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - tmp_S_phi);
tmp_T_k_p_wk_ = exp(+2*pi*i*(k_c_0_wk_*tmp_T_delta_x_c_(1+0) + k_c_1_wk_*tmp_T_delta_x_c_(1+1)));
tmp_plane_T_k_p_wk_ = 2*k_p_r_wk_.*cos(k_p_w_wk_ - tmp_T_phi);
tmp_S_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
tmp_T_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_T_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum );
%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; np=0;
Slim_ = max(abs(tmp_S_k_p_wk_),[],'all')*[-1,+1];
Tlim_ = max(abs(tmp_T_k_p_wk_),[],'all')*[-1,+1];
fontsize_use = 12;
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_S_k_p_wk_),Slim_,colormap_80s()); axis image; axisnotick; title('real(tmp_S_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(tmp_S_k_p_wk_),Slim_,colormap_80s()); axis image; axisnotick; title('imag(tmp_S_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_S_x_c_xx_),Slim_,colormap_beach()); axis image; axisnotick; title('real(tmp_S_x_c_xx_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,imag(tmp_S_x_c_xx_),Slim_,colormap_beach()); axis image; axisnotick; title('imag(tmp_S_x_c_xx_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_T_k_p_wk_),Tlim_,colormap_80s()); axis image; axisnotick; title('real(tmp_T_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(tmp_T_k_p_wk_),Tlim_,colormap_80s()); axis image; axisnotick; title('imag(tmp_T_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(tmp_T_x_c_xx_),Tlim_,colormap_beach()); axis image; axisnotick; title('real(tmp_T_x_c_xx_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,imag(tmp_T_x_c_xx_),Tlim_,colormap_beach()); axis image; axisnotick; title('imag(tmp_T_x_c_xx_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
drawnow;
fname_fig_pre = sprintf('%s/test_transforms_test_weight_2d_wk_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%;
I_quad = sum(conj(tmp_plane_S_k_p_wk_.*tmp_S_k_p_wk_).*(tmp_plane_T_k_p_wk_.*tmp_T_k_p_wk_).*weight_2d_wk_)*(2*pi)^2;
I_form = I_xPPx_0(k_p_r_max,tmp_S_phi,tmp_S_delta_x_c_,tmp_T_phi,tmp_T_delta_x_c_);
fnorm_disp(flag_verbose,'I_form',I_form,'I_quad',I_quad,' %%<-- should be <1e-6');
%%%%%%%%;

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
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
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
%%%%%%%%;
a_k_Y_form_ = zeros(n_lm_sum,1);
for nsource=0:n_source-1;
a_k_Y_form_ = a_k_Y_form_ + plane_wave_expansion_1(n_k_p_r,k_p_r_,delta_a_c_3s__(:,1+nsource),l_max_);
end;%for nsource=0:n_source-1;
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
,a_k_p_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_quad_',a_k_Y_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_quad_ ...
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
,a_k_Y_form_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_quad_',a_k_p_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_p_reco_ ...
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
,a_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_reco_',a_k_p_reco_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'a_k_p_quad_',a_k_p_quad_,'a_k_p_reco_',a_k_p_reco_,' %%<-- should be <1e-2');
%%%%%%%%;
tmp_t = tic;
[ ...
 a_k_Y_reco_ ...
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
,a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_reco_',a_k_Y_reco_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'a_k_Y_quad_',a_k_Y_quad_,'a_k_Y_reco_',a_k_Y_reco_,' %%<-- should be <1e-2');
%%%%%%%%;
a_k_Y_form_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_form_);
a_k_Y_form_yk___ = zeros(1+l_max_max,n_m_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
tmp_a_k_Y_form_lm_ = a_k_Y_form_(1+tmp_index_);
tmp_a_k_Y_form_lm__ = zeros(1+l_max_max,n_m_max);
l_max = l_max_(1+nk_p_r);
na=0;
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(na==l_val*(l_val+1)+m_val);
tmp_a_k_Y_form_lm__(1+l_val,1+l_max_max+m_val) = tmp_a_k_Y_form_lm_(1+na);
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
a_k_Y_form_yk___(:,:,1+nk_p_r) = tmp_a_k_Y_form_lm__;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
flag_check=1;
if flag_check;
na=0;
for nk_p_r=0:n_k_p_r-1;
l_max=l_max_(1+nk_p_r);
for l_val=0:l_max;
for m_val=-l_val:+l_val;
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_yk__(1+l_val*(l_val+1)+m_val,1+nk_p_r));
assert(a_k_Y_form_(1+n_lm_csum_(1+nk_p_r)+l_val*(l_val+1)+m_val)==a_k_Y_form_yk___(1+l_val,1+l_max_max+m_val,1+nk_p_r));
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_lm_sum);
end;%if flag_check;
%%%%%%%%;
a_k_Y_form_lkm___ = permute(a_k_Y_form_yk___,1+[0,2,1]);
%%%%%%%%;
a_k_Y_l2_quad = sum((conj(a_k_Y_form_yk__).*a_k_Y_form_yk__)*reshape(weight_3d_k_p_r_,[n_k_p_r,1]));
disp(sprintf(' %% a_k_Y_l2_quad %+0.6f a_k_p_l2_form %+0.6f',a_k_Y_l2_quad,a_k_p_l2_form));

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
% Now generate templates from the volumes in spherical-harmonic coordinates. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_form_yk__,[n_lm_max,n_k_p_r]) ...
,template_k_eq_d_double/max(1e-12,k_p_r_max) ...
,-1 ...
,n_w_max ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t)); end;
n_S = n_viewing_S;
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
S_k_p_l2_quad_S_ = reshape(sum(bsxfun(@times,conj(S_k_p_wkS__).*S_k_p_wkS__,reshape(weight_2d_wk_,[n_w_sum,1])),1),[n_S,1])*(2*pi)^2;
%%%%%%%%;

%%%%%%%%;
% The 3d-frequency-space locations associated with a particular template can be calculated as follows: ;
% (see get_template_1.m);
%%%%%%%%;
nS = max(0,min(n_S-1,floor(0.2*n_S)));
tmp_viewing_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_viewing_polar_a = viewing_polar_a_S_(1+nS);
tmp_viewing_gamma_z = 0.0;
tmp_cc_ = cos(k_p_w_wk_); tmp_sc_ = sin(k_p_w_wk_);
tmp_cb = cos(tmp_viewing_azimu_b); tmp_sb = sin(tmp_viewing_azimu_b);
tmp_ca = cos(tmp_viewing_polar_a); tmp_sa = sin(tmp_viewing_polar_a);
tmp_k_c_0_wk_ = (+tmp_cb*tmp_ca*tmp_cc_ - tmp_sb*tmp_sc_).*k_p_r_wk_;
tmp_k_c_1_wk_ = (+tmp_sb*tmp_ca*tmp_cc_ + tmp_cb*tmp_sc_).*k_p_r_wk_;
tmp_k_c_2_wk_ = (-tmp_sa*tmp_cc_                        ).*k_p_r_wk_;
tmp_0_k_c_wk3__ = cat(2,tmp_k_c_0_wk_,tmp_k_c_1_wk_,tmp_k_c_2_wk_);
%%%%%%%%;
% which can also be calculated via: ;
% (see imagesc_S_k_p_3d_2.m);
%%%%%%%%;
k_c_wk2__ = cat(2,cos(k_p_w_wk_).*k_p_r_wk_,sin(k_p_w_wk_).*k_p_r_wk_);
k_c_wk3__ = cat(2,k_c_wk2__,zeros(n_w_sum,1));
tmp_1_k_c_wk3__ = k_c_wk3__*transpose(Ry(tmp_viewing_polar_a))*transpose(Rz(tmp_viewing_azimu_b));
%%%%%%%%;
fnorm_disp(flag_verbose,'tmp_0_k_c_wk3__',tmp_0_k_c_wk3__,'tmp_1_k_c_wk3__',tmp_1_k_c_wk3__,' %%<-- should be zero');
%%%%%%%%;

%%%%%%%%;
% Now step through and reconstitute the templates themselves. ;
%%%%%%%%;
T_k_p_l2_quad_S_ = zeros(n_S,1);
T_k_p_l2_form_S_ = zeros(n_S,1);
T_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_l2_quad_S_(1+nS) = sum(conj(T_k_p_wk_).*T_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
for nsource1=0:n_source-1;
tmp_delta_a_c_0_ = tmp_R__*delta_a_c_3s__(:,1+nsource0);
tmp_delta_a_c_1_ = tmp_R__*delta_a_c_3s__(:,1+nsource1);
tmp_kd = 2*pi*k_p_r_max*fnorm(tmp_delta_a_c_0_(1:2) - tmp_delta_a_c_1_(1:2));
tmp_h2d = (2*pi)^2; if abs(tmp_kd)>1e-12; tmp_h2d = h2d_(tmp_kd); end;
T_k_p_l2_form_S_(1+nS) = T_k_p_l2_form_S_(1+nS) + tmp_h2d/(2*pi)^2 * (pi*k_p_r_max^2);
end;%for nsource1=0:n_source-1;
end;%for nsource0=0:n_source-1;
if (flag_disp>1);
if mod(nS,128)==0;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_)); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none');
subplot(1,2,2);imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_)); axis image; axisnotick; title('real(T_k_p_wk_)','Interpreter','none');
sgtitle(sprintf('nS %d/%d',nS,n_S));
end;%if mod(nS,128)==0;
end;%if (flag_disp>1);
T_k_p_wkS__(:,1+nS) = T_k_p_wk_;
end;%for nS=0:n_S-1;
fnorm_disp(flag_verbose,'T_k_p_wkS__',T_k_p_wkS__,'S_k_p_wkS__',S_k_p_wkS__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_l2_form_S_',T_k_p_l2_form_S_,'T_k_p_l2_quad_S_',T_k_p_l2_quad_S_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_l2_form_S_',T_k_p_l2_form_S_,'S_k_p_l2_quad_S_',S_k_p_l2_quad_S_,' %%<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 12;
p_row = 3; p_col = 4; n_plot = 4; np=0;
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_beach();
a_x_c_form_lim_ = prctile(mean(reshape(real(a_x_c_form___),[n_x_c,n_x_c,n_x_c]),3),[ 5,95],'all');
tmp_parameter.vlim_ = 1.0*a_x_c_form_lim_;
isosurface_f_x_u_1(tmp_parameter,a_x_c_form___);
xlabel('x0'); ylabel('x1'); zlabel('x2'); axisnotick3d;
title('a_x_c_form___','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
subplot(p_row,p_col,1+np);np=np+1;cla;
tmp_parameter = struct('type','parameter');
tmp_parameter.c_use__ = colormap_80s();
a_k_c_form_lim_ = prctile(mean(reshape(real(a_k_c_form___),[n_k_c,n_k_c,n_k_c]),3),[ 5,95],'all');
tmp_parameter.vlim_ = 1.0*a_k_c_form_lim_;
tmp_parameter.percent_threshold_ = [ 1,99];
isosurface_f_x_u_1(tmp_parameter,a_k_c_form___);
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d;
title('a_k_c_form___','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
nS = max(0,min(n_S-1,floor(0.2*n_S)));
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
nk_p_r = max(0,min(n_k_p_r-1,floor(1.0*n_k_p_r)));
k_p_r = k_p_r_(1+nk_p_r);
n_qk_csum_pre = n_qk_csum_(1+nk_p_r+0);
n_qk_csum_pos = n_qk_csum_(1+nk_p_r+1);
assert(qref_n_shell==n_qk_csum_pos-n_qk_csum_pre);
a_k_p_sub_fin_ = a_k_p_form_(1+n_qk_csum_pre+[0:qref_n_shell-1]);
a_k_p_sub_fin_lim_ = prctile(abs(a_k_p_sub_fin_),95)*[-1,+1];
a_k_p_ori_sub_fin_lim_ = a_k_p_sub_fin_lim_;
%%%%%%%%;
Slim_ = max(abs(S_k_p_wk_),[],'all')*[-1,+1];
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_,colormap_80s()); axis image; axisnotick; title('real(S_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),Slim_,colormap_80s()); axis image; axisnotick; title('imag(S_k_p_wk_)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
%%%%%%%%;
for ntype=0:1;
if ntype==0; tmp_S_k_p_wk_ = real(S_k_p_wkS__(:,1+nS)); tmp_str = 'real'; end;
if ntype==1; tmp_S_k_p_wk_ = imag(S_k_p_wkS__(:,1+nS)); tmp_str = 'imag'; end;
for nplot=0:n_plot-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
nk_p_r = max(0,min(n_k_p_r-1,floor((nplot)/max(1,n_plot)*n_k_p_r)));
k_p_r = k_p_r_(1+nk_p_r);
n_qk_csum_pre = n_qk_csum_(1+nk_p_r+0);
n_qk_csum_pos = n_qk_csum_(1+nk_p_r+1);
assert(qref_n_shell==n_qk_csum_pos-n_qk_csum_pre);
a_k_p_sub_ = a_k_p_form_(1+n_qk_csum_pre+[0:qref_n_shell-1]);
a_k_p_sub_lim_ = prctile(abs(a_k_p_sub_),95)*[-1,+1];
%%%%;
hold on;
if ntype==0; tmp_a_k_p_sub_ = real(a_k_p_sub_); end; if ntype==1; tmp_a_k_p_sub_ = imag(a_k_p_sub_); end;
imagesc_polar_a_azimu_b_0(qref_polar_a_shell_,qref_azimu_b_shell_,tmp_a_k_p_sub_,a_k_p_ori_sub_fin_lim_,colormap_80s,0,k_p_r);
%n_contour = 16;
tmp_parameter = struct('type','parameter');
tmp_parameter.n_contour = 0;
tmp_parameter.flag_percent_use=0;
tmp_parameter.c_use__ = colormap_80s;
tmp_parameter.vlim_ = 1.0*a_k_p_ori_sub_fin_lim_;
%tmp_parameter.vval_ = transpose(linspace(min(tmp_parameter.vlim_),max(tmp_parameter.vlim_),n_contour));
tmp_parameter.flag_k_c_interp = 0;
imagesc_S_k_p_3d_2( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_wk_ ...
,1 ...
,tmp_S_k_p_wk_ ...
,viewing_azimu_b_S_(1+nS) ...
,viewing_polar_a_S_(1+nS) ...
);
hold off;
xlabel('k0'); ylabel('k1'); zlabel('k2'); axisnotick3d;
title(sprintf('nk_p_r %d k_p_r %.2f %s',nk_p_r,k_p_r,tmp_str),'Interpreter','none');
set(gca,'FontSize',fontsize_use);
axis equal;
drawnow;
end;%for nplot=0:n_plot-1;
end;%for ntype=0:1;
%%%%%%%%;
sgtitle(sprintf('azimu_b %+0.2f pi polar_a %+0.2f pi',viewing_azimu_b_S_(1+nS)/pi,viewing_polar_a_S_(1+nS)/pi),'Interpreter','none');
fname_fig_pre = sprintf('%s/test_transforms_S_k_p_wk_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now test out direct spherical-harmonic contruction of S_k_q_wk_. ;
% Note that below we introduce a (unnecessary) tmp_gamma_z to check consistency. ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*2/3)));
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
%tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = pi/12;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'T_k_p_wk_',T_k_p_wk_,' %%<-- should be <1e-2');
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'T_k_q_wk_',T_k_q_wk_,' %%<-- should be <1e-2');
%%%%;
S_k_p_l2_qua2 = sum(conj(S_k_p_wk_).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
S_k_p_l2_quad = S_k_p_l2_quad_S_(1+nS);
fnorm_disp(flag_verbose,'S_k_p_l2_quad',S_k_p_l2_quad,'S_k_p_l2_qua2',S_k_p_l2_qua2,' %%<-- should be <1e-9');
%%%%;
tmp_euler_ = [-tmp_azimu_b,-tmp_polar_a,-tmp_gamma_z];
W_beta__ = wignerd_b(l_max_max,+tmp_euler_(1+1));
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
[ ...
 a_k_Y_rota_ ...
] = ...
rotate_spharm_to_spharm_2( ...
 0*flag_verbose ...
,[] ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,a_k_Y_form_ ...
,tmp_euler_ ...
);
%%%%;
a_k_Y_rota_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_rota_);
a_k_Y_l2_rota = sum(conj(a_k_Y_rota_yk__).*a_k_Y_rota_yk__*reshape(weight_3d_k_p_r_,[n_k_p_r,1]));
%%%%;
Q_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm_(1+nk_p_r)-1];
a_k_Y_sub_ = a_k_Y_rota_(1+tmp_index_);
l_max = l_max_(1+nk_p_r);
for nw=0:n_w_max-1;
nq = periodize(nw,-n_w_max/2,+n_w_max/2); m0_val = nq;
for l_val=0:l_max;
if abs(m0_val)<=l_val;
Q_k_q_wk_(1+nw+nk_p_r*n_w_max) = Q_k_q_wk_(1+nw+nk_p_r*n_w_max) + ...
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
*a_k_Y_sub_(1+(l_val^2 + l_val + m0_val)) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
;
end;%if abs(m0_val)<=l_val;
end;%for l_val=0:l_max;
end;%for nw=0:n_w_max-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + [0:n_lm_(1+nk_p_r)-1];
a_k_Y_sub_ = a_k_Y_form_(1+tmp_index_);
l_max = l_max_(1+nk_p_r);
for nw=0:n_w_max-1;
nq = periodize(nw,-n_w_max/2,+n_w_max/2); m0_val = nq;
for l_val=0:l_max;
for m1_val=-l_val:+l_val;
if abs(m0_val)<=l_val;
R_k_q_wk_(1+nw+nk_p_r*n_w_max) = R_k_q_wk_(1+nw+nk_p_r*n_w_max) + ...
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
*exp(-i*m0_val*(+tmp_euler_(1+2))) ...
*W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*exp(-i*m1_val*(+tmp_euler_(1+0))) ...
*a_k_Y_sub_(1+(l_val^2 + l_val + m1_val)) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
;
end;%if abs(m0_val)<=l_val;
end;%for m1_val=-l_val:+l_val;
end;%for l_val=0:l_max;
end;%for nw=0:n_w_max-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
Q_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,Q_k_q_wk_);
R_k_p_l2_qua2 = sum(conj(R_k_p_wk_).*R_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
Q_k_p_l2_qua2 = sum(conj(Q_k_p_wk_).*Q_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
%%%%;
fnorm_disp(flag_verbose,'R_k_q_wk_',R_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_q_wk_',S_k_q_wk_,'R_k_q_wk_',R_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_q_wk_',T_k_q_wk_,'Q_k_q_wk_',Q_k_q_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_q_wk_',T_k_q_wk_,'R_k_q_wk_',R_k_q_wk_,' %%<-- should be <1e-2');
%%%%%%%%;
fnorm_disp(flag_verbose,'R_k_p_wk_',R_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_p_wk_',S_k_p_wk_,'R_k_p_wk_',R_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_wk_',T_k_p_wk_,'Q_k_p_wk_',Q_k_p_wk_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_wk_',T_k_p_wk_,'R_k_p_wk_',R_k_p_wk_,' %%<-- should be <1e-2');
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 4; np=0;
Slim_ = max(abs(S_k_p_wk_))*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(Q_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('real(Q_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(Q_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('imag(Q_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('real(R_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(R_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('imag(R_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('real(S_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('imag(S_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('real(T_k_p_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_wk_),Slim_,colormap_80s()); axisnotick; axis image; title(sprintf('imag(T_k_p_wk_)'),'Interpreter','none');
sgtitle(sprintf(' nS %d/%d',nS,n_S),'Interpreter','none');
fname_fig_pre = sprintf('%s/test_transforms_S_k_p_wk_nS%.3d_FIGA',dir_jpg,nS);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 4; np=0;
Slim_ = max(abs(S_k_q_wk_))*[-1,+1];
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(Q_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('real(Q_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(Q_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('imag(Q_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('real(R_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(R_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('imag(R_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('real(S_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('imag(S_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('real(T_k_q_wk_)'),'Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla;
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_q_wk_),Slim_,colormap_81s()); axisnotick; axis tight; title(sprintf('imag(T_k_q_wk_)'),'Interpreter','none');
sgtitle(sprintf(' nS %d/%d',nS,n_S),'Interpreter','none');
fname_fig_pre = sprintf('%s/test_transforms_S_k_q_wk_nS%.3d_FIGA',dir_jpg,nS);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now we try and set up a template-operator ;
% for a collection of azimu_b associated with a single polar_a. ;
%%%%%%%%;
polar_a_use = viewing_polar_a_S_(round(n_S/4));
tmp_index_ = efind(abs(viewing_polar_a_S_-polar_a_use)<1e-6);
n_azimu_b_use = numel(tmp_index_);
azimu_b_use_ = viewing_azimu_b_S_(1+tmp_index_);
S_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
S_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
T_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
nS = tmp_index_(1+nazimu_b_use);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
tmp_azimu_b = viewing_azimu_b_S_(1+nS);
tmp_polar_a = viewing_polar_a_S_(1+nS);
tmp_gamma_z = 0.0; %<-- default. ;
tmp_gamma_z = -pi/9;
S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,-tmp_gamma_z);
S_k_p_sub_wkb__(:,1+nazimu_b_use) = S_k_p_wk_;
tmp_R__ = Rz(-tmp_gamma_z)*Ry(-tmp_polar_a)*Rz(-tmp_azimu_b);
T_k_p_wk_ = zeros(n_w_sum,1);
for nsource=0:n_source-1;
tmp_delta_ = tmp_R__*delta_a_c_3s__(:,1+nsource);
T_k_p_wk_ = T_k_p_wk_ + exp(+i*2*pi*(k_c_0_wk_*tmp_delta_(1+0) + k_c_1_wk_*tmp_delta_(1+1)));
end;%for nsource=0:n_source-1;
T_k_p_sub_wkb__(:,1+nazimu_b_use) = T_k_p_wk_;
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
T_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,T_k_p_wk_);
S_k_q_sub_wkb__(:,1+nazimu_b_use) = S_k_q_wk_;
T_k_q_sub_wkb__(:,1+nazimu_b_use) = T_k_q_wk_;
clear S_k_p_wk_ T_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_q_sub_wkb__',S_k_q_sub_wkb__,'T_k_q_sub_wkb__',T_k_q_sub_wkb__,' %%<-- should be <1e-2');
%%%%;
tmp_t = tic();
W_beta__ = wignerd_b(l_max_max,-polar_a_use);
zeta_lm__ = zeros(1+l_max_max,n_m_max);
for l_val=0:l_max_max;
a1=((2*l_val+1)/(4*pi));
Llm__ = legendre(l_val,0,'unnorm');
for m_val=-l_val:+l_val;
if (l_val >0); Llm_ = Llm__(1+abs(m_val),:); end; if (l_val==0); Llm_ = Llm__; end; assert(numel(Llm_)==1);
a2=exp(0.5*lfactorial(l_val-abs(m_val)) - 0.5*lfactorial(l_val+abs(m_val))); c=sqrt(a1)*a2; s=1; % original phase ;
zeta_lm__(1+l_val,1+l_max_max+m_val) = s*c*Llm_(1+0);
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
%%%%;
W_betazeta_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
W_betaones_mlm___ = zeros(n_m_max,1+l_max_max,n_m_max);
for l_val=0:l_max_max;
for m0_val=-l_val:+l_val;
for m1_val=-l_val:+l_val;
W_betazeta_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*zeta_lm__(1+l_val,1+l_max_max+m0_val) ...
 ;
W_betaones_mlm___(1+l_max_max+m0_val,1+l_val,1+l_max_max+m1_val) = ...
 W_beta__{1+l_val}(1+l_val+m0_val,1+l_val+m1_val) ...
*1 ...
 ;
end;%for m1_val=-l_val:+l_val;
end;%for m0_val=-l_val:+l_val;
end;%for l_val=0:l_max_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_betazeta_mlm___: %0.2fs',tmp_t)); end;
%%%%;
flag_check=1;
if flag_check;
tmp_w_ = randn(n_m_max) + i*randn(n_m_max);
tmp_azimu_b_use_ = 2*pi*rand(n_azimu_b_use,1);
tmp_f__ = exp(-i*(-reshape(tmp_azimu_b_use_,[n_azimu_b_use,1]))*reshape(m_max_,[1,n_m_max]));
tmp_fw_0_ = tmp_f__*tmp_w_;
tmp_fw_1_ = xxnufft1d2(n_azimu_b_use,tmp_azimu_b_use_,+1,1e-6,n_m_max,tmp_w_);
fnorm_disp(flag_verbose,'tmp_fw_0_',tmp_fw_0_,'tmp_fw_1_',tmp_fw_1_,' %%<-- should be <1e-6');
clear tmp_azimu_b_use_ tmp_w_ tmp_f__ tmp_fw_0_ tmp_fw_1_ ;
end;%if flag_check;
%%%%;
tmp_t = tic();
W_caza_mkm___ = zeros(n_m_max,n_k_p_r,n_m_max); %<-- diag(exp(+i*m0_val_*tmp_gamma_z))*W_betazeta_ml__*a_k_Y_form_lk__ for each m1_val. ;
for m1_val=-l_max_max:+l_max_max;
W_caza_mkm___(:,:,1+l_max_max+m1_val) = ...
 sqrt(n_w_max) ... %<-- this scaling factor is to ensure that interp_q_to_p is scaled correctly. ;
 *diag(exp(-i*m_max_*(-tmp_gamma_z))) ...
 *reshape(W_betazeta_mlm___(:,:,1+l_max_max+m1_val),[n_m_max,1+l_max_max]) ...
 *reshape(a_k_Y_form_lkm___(:,:,1+l_max_max+m1_val),[1+l_max_max,n_k_p_r]) ...
 ;
end;%for m1_val=-l_max_max:+l_max_max;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_caza_mkm___: %0.2fs',tmp_t)); end;
tmp_t = tic();
W_caza_mmk___ = permute(W_caza_mkm___,1+[2,0,1]);
W_caza_bmk___ = reshape(xxnufft1d2(n_azimu_b_use,azimu_b_use_,+1,1e-6,n_m_max,reshape(W_caza_mmk___,[n_m_max,n_m_max*n_k_p_r])),[n_azimu_b_use,n_m_max,n_k_p_r]);
W_caza_mkb___ = permute(W_caza_bmk___,1+[1,2,0]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% W_caza_mkb___: %0.2fs',tmp_t)); end;
%%%%%%%%;
R_k_q_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
R_k_p_sub_wkb__ = zeros(n_w_sum,n_azimu_b_use);
for nazimu_b_use=0:n_azimu_b_use-1;
R_k_p_wk_ = zeros(n_w_sum,1);
R_k_q_wk_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
for m_val=-l_max_max:+l_max_max;
nq = m_val; if (nq<0); nq=nq+n_w_max; end;
R_k_q_wk_(1+nq+nk_p_r*n_w_max) = W_caza_mkb___(1+l_max_max+m_val,1+nk_p_r,1+nazimu_b_use);
end;%for m_val=-l_max_max:+l_max_max;
end;%for nk_p_r=0:n_k_p_r-1;
R_k_q_sub_wkb__(:,1+nazimu_b_use) = R_k_q_wk_;
R_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,R_k_q_wk_);
R_k_p_sub_wkb__(:,1+nazimu_b_use) = R_k_p_wk_;
end;%for nazimu_b_use=0:n_azimu_b_use-1;
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,' %%<-- should be <1e-2');
%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
Slim_ = max(abs(S_k_p_sub_wkb__),[],'all')*[-1,+1];
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_sub_wkb__(:,1+nazimu_b_use)),Slim_,colormap_80s); axis image; axisnotick;
%title(sprintf('real(S_k_p_sub_wkb__(:,1+%d))',nazimu_b_use),'Interpreter','none');
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
fname_fig_pre = sprintf('%s/test_transforms_S_k_p_sub_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>0);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 6; p_col = ceil(n_azimu_b_use/p_row); np=0;
Slim_ = max(abs(S_k_p_sub_wkb__),[],'all')*[-1,+1];
for nazimu_b_use=0:n_azimu_b_use-1;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(R_k_p_sub_wkb__(:,1+nazimu_b_use)),Slim_,colormap_80s); axis image; axisnotick;
%title(sprintf('real(R_k_p_sub_wkb__(:,1+%d))',nazimu_b_use),'Interpreter','none');
title(sprintf('nazimu_b_use %d',nazimu_b_use),'Interpreter','none');
end;%for nazimu_b_use=0:n_azimu_b_use-1;
fname_fig_pre = sprintf('%s/test_transforms_R_k_p_sub_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
end;%if (flag_disp>1);
%%%%%%%%;
[ ...
 U_k_p_sub_wkb__ ...
,W_betazeta_mlm___ ...
] = ...
sph_template_single_polar_a_3( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_form_lkm___ ...
,n_w_max ...
,polar_a_use ...
,n_azimu_b_use ...
,azimu_b_use_ ...
,tmp_gamma_z ...
,W_betazeta_mlm___ ...
);
fnorm_disp(flag_verbose,'S_k_p_sub_wkb__',S_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'T_k_p_sub_wkb__',T_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'R_k_p_sub_wkb__',R_k_p_sub_wkb__,'U_k_p_sub_wkb__',U_k_p_sub_wkb__,' %%<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%;
% We check the following orthogonality relation: ;
% \sum_{nq} conj(W_beta__{1+l0}(1+nq,1+m0)) * W_beta__{1+l1}(1+nq,1+m1) ;
% Verdict: should be delta_{m0,m1} when l0==l1. ;
%%%%%%%%;
flag_check=1;
if flag_check;
sum0_lmlm____ = zeros(1+l_max_max,n_m_max,1+l_max_max,n_m_max);
for l0_val=0:l_max_max;
disp(sprintf(' %% slow calculation: l0_val %d/%d',l0_val,l_max_max));
for l1_val=0:l_max_max;
lm_val = min(l0_val,l1_val);
for m0_val=-l0_val:+l0_val;
for m1_val=-l1_val:+l1_val;
for q_val=-lm_val:+lm_val;
sum0_lmlm____(1+l0_val,1+l_max_max+m0_val,1+l1_val,1+l_max_max+m1_val) = ...
 sum0_lmlm____(1+l0_val,1+l_max_max+m0_val,1+l1_val,1+l_max_max+m1_val) ...
+ 1 ...
*conj(W_beta__{1+l0_val}(1+l0_val+q_val,1+l0_val+m0_val)) ...
*W_beta__{1+l1_val}(1+l1_val+q_val,1+l1_val+m1_val) ...
;
end;%for q_val=-lm_val:+lm_val;
end;%for m1_val=-l1_val:+l1_val;
end;%for m0_val=-l0_val:+l0_val;
end;%for l1_val=0:l_max_max;
end;%for l0_val=0:l_max_max;
%%%%;
%figure(1);clf;figsml;imagesc(reshape(log10(abs(sum0_lmlm____)),[(1+l_max_max)*n_m_max,(1+l_max_max)*n_m_max]),[-5,0]);colorbar; axis image;
%figure(2);clf;figsml;imagesc(reshape(log10(abs(permute(sum0_lmlm____,1+[0,2,1,3]))),[(1+l_max_max)^2,n_m_max^2]),[-5,0]);colorbar; axis image;
%%%%;
disp(sprintf(' %% nnz sum0_lmlm____: %d / %d --> %0.2f',nnz(sum0_lmlm____),numel(sum0_lmlm____),nnz(sum0_lmlm____)/numel(sum0_lmlm____)));
%%%%;
sum1_lmlm____ = zeros(1+l_max_max,n_m_max,1+l_max_max,n_m_max);
tmp_mlm__ = reshape(W_betaones_mlm___,[n_m_max,(1+l_max_max)*n_m_max]);
sum1_lmlm____ = reshape(ctranspose(tmp_mlm__)*tmp_mlm__,[1+l_max_max,n_m_max,1+l_max_max,n_m_max]);
disp(sprintf(' %% nnz sum1_lmlm____: %d / %d --> %0.2f',nnz(sum1_lmlm____),numel(sum1_lmlm____),nnz(sum1_lmlm____)/numel(sum1_lmlm____)));
fnorm_disp(flag_verbose,'sum0_lmlm____',sum0_lmlm____,'sum1_lmlm____',sum1_lmlm____,' %%<-- should be close to zero');
%%%%;
error_kroneker = 0;
for l_val=0:l_max_max;
tmp_mm__ = squeeze(sum0_lmlm____(1+l_val,1+l_max_max+[-l_val:+l_val],1+l_val,1+l_max_max+[-l_val:+l_val]));
error_kroneker = error_kroneker + fnorm(tmp_mm__ - eye(1+2*l_val));
end;%for l_val=0:l_max_max;
disp(sprintf(' %% error_kroneker: %0.16f %<-- should be close to zero',error_kroneker));
%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
S_k_q_wkS__(:,1+nS) = S_k_q_wk_;
end;%for nS=0:n_S-1;
tmp_ = reshape(interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__),[n_w_sum,n_S]); %<-- this accomplishes the same thing. ;
fnorm_disp(flag_verbose,'S_k_q_wkS__',S_k_q_wkS__,'tmp_',tmp_,' %%<-- should be zero');
S_k_q_wSk___ = permute(reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),1+[0,2,1]);
%%%%%%%%;

flag_delta = 1;
%%%%%%%%;
if ~flag_delta; delta_r_max = 0.0/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested =   1; end;
if  flag_delta; delta_r_max = 0.5/max(1e-12,k_p_r_max); svd_eps = 1e-6; n_delta_v_requested = 128; end;
FTK = ampmh_FTK_2(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
%%%%%%%%;
n_delta_v = FTK.n_delta_v;

flag_CTF = 1;
%%%%%%%%;
% Now we generate images for the likelihood-calculation. ;
% This test accounts for anisotropic CTF. ;
% We also shift each image by an 'on-grid' displacement. ;
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
% quickly test vectorized construction. ;
%%%%%%%%;
N_k_p_wkM__ = zeros(n_w_sum,n_M);
N_k_p_wkM__ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_).*reshape(transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+index_nS_from_nM_),gamma_z_(1+index_nw_from_nM_)),-FTK.delta_x_(1+index_nd_from_nM_),-FTK.delta_y_(1+index_nd_from_nM_)),[n_w_sum,n_M]);
N_k_q_wkM__ = reshape(interp_p_to_q(n_k_p_r,n_w_,n_w_sum,N_k_p_wkM__),[n_w_sum,n_M]);
fnorm_disp(flag_verbose,'M_k_p_wkM__',M_k_p_wkM__,'N_k_p_wkM__',N_k_p_wkM__,' %%<-- should be zero');
fnorm_disp(flag_verbose,'M_k_q_wkM__',M_k_q_wkM__,'N_k_q_wkM__',N_k_q_wkM__,' %%<-- should be zero');
%%%%%%%%;

%%%%%%%%;
% quickly test inner-product: <R(+gamma_z)*M_k_p_wk_,S_k_p_wk_> ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*2/3)));
nM = max(0,min(n_M-1,round(n_M*1/5)));
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM));
M_k_p_wk_ = CTF_k_p_wk_.*M_k_p_wkM__(:,1+nM);
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_);
tmp_1_w_ = ifft(sum(bsxfun(@times,reshape(conj(M_k_q_wk_).*S_k_q_wk_,[n_w_max,n_k_p_r]),reshape(weight_2d_k_p_r_,[1,n_k_p_r])),2));
tmp_0_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
tmp_0_w_(1+nw) = sum(conj(rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_,+gamma_z)).*S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'tmp_0_w_',tmp_0_w_,'tmp_1_w_',tmp_1_w_,' %%<-- should be zero');
%%%%%%%%;

%%%%%%%%;
% Now construct the template-norms: ;
% <(R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_,(R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_> ;
% Note that this does not involve collapsing onto principal-modes. ;
%%%%%%%%;
R_CTF_S_l2_wSC_quad___ = zeros(n_w_max,n_S,n_CTF);
SS_k_p_wkS__ = conj(S_k_p_wkS__).*S_k_p_wkS__;
SS_k_q_wkS__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,SS_k_p_wkS__),[n_w_sum,n_S]);
CC_k_p_wkC__ = conj(CTF_k_p_wkC__).*CTF_k_p_wkC__;
CC_k_q_wkC__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,CC_k_p_wkC__),[n_w_sum,n_CTF]);
for nCTF=0:n_CTF-1;
R_CTF_S_l2_wSC_quad___(:,:,1+nCTF) = ifft(squeeze(sum(bsxfun(@times,reshape(bsxfun(@times,conj(CC_k_q_wkC__(:,1+nCTF)),SS_k_q_wkS__),[n_w_max,n_k_p_r,n_S]),reshape(weight_2d_k_p_r_,[1,n_k_p_r,1])),1+1)));
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
R_CTF_S_l2_w_quad_ = R_CTF_S_l2_wSC_quad___(:,1+nS,1+nCTF);
R_CTF_S_l2_w_qua2_ = zeros(n_w_max,1);
R_CTF_S_l2_w_form_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
R_CTF_S_k_p_wk_ = S_k_p_wkS__(:,1+nS).*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),+gamma_z);
R_CTF_S_l2_w_qua2_(1+nw) = sum(conj(R_CTF_S_k_p_wk_).*R_CTF_S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
for nsource0=0:n_source-1;
tmp_S_delta_0_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource0); tmp_S_delta_0_2_ = tmp_S_delta_0_3_(1:2);
for nsource1=0:n_source-1;
tmp_S_delta_1_3_ = tmp_R_S__*delta_a_c_3s__(:,1+nsource1); tmp_S_delta_1_2_ = tmp_S_delta_1_3_(1:2);
R_CTF_S_l2_w_form_(1+nw) = R_CTF_S_l2_w_form_(1+nw) + I_xPPx_0(k_p_r_max,tmp_CTF_phi+gamma_z,tmp_S_delta_0_2_,tmp_CTF_phi+gamma_z,tmp_S_delta_1_2_);
end;%for nsource0=0:n_source-1;
end;%for nsource1=0:n_source-1;
end;%for nw=0:n_w_max-1;
fnorm_disp(flag_verbose,'R_CTF_S_l2_w_qua2_',R_CTF_S_l2_w_qua2_,'R_CTF_S_l2_w_quad_',R_CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'R_CTF_S_l2_w_qua2_',R_CTF_S_l2_w_qua2_,'R_CTF_S_l2_w_form_',R_CTF_S_l2_w_form_,' %%<-- should be <1e-6');
fnorm_disp(flag_verbose,'R_CTF_S_l2_w_quad_',R_CTF_S_l2_w_quad_,'R_CTF_S_l2_w_form_',R_CTF_S_l2_w_form_,' %%<-- should be <1e-6');
%%%%%%%%;
% check against dir_ascii. ;
%%%%%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_ascii = 3;
for nascii=0:n_ascii-1;
na=0;
if nascii==na; tmp_str_ascii = 'R_CTF_S_l2_w_quad_'; tmp_local = R_CTF_S_l2_w_quad_; tab_p=0; end; na=na+1;
if nascii==na; tmp_str_ascii = 'R_CTF_S_l2_w_qua2_'; tmp_local = R_CTF_S_l2_w_qua2_; tab_p=0; end; na=na+1;
if nascii==na; tmp_str_ascii = 'R_CTF_S_l2_w_form_'; tmp_local = R_CTF_S_l2_w_form_; tab_p=0; end; na=na+1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_ascii,tmp_local,'ascii',tmp_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
end;%for nascii=0:n_ascii-1;
%%%%%%%%;

flag_pm = 1;
%%%%%%%%%%%%%%%%;
% Now calculate innerproduct Z_dwSM____. ;
% Here we account for anisotropic CTF. ;
%%%%%%%%%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. Note that, in principle, this could be different for different nCTF. ; 
pm_n_UX_rank = n_UX_rank; %<-- just to check dimension. ;
Z_dwSM_ampm____ = zeros(n_delta_v,n_w_max,n_S,n_M);
UX_M_l2_M_ = zeros(n_M,1);
UX_T_M_l2_dM__ = zeros(n_delta_v,n_M);
UX_knC___ = zeros(n_k_p_r,n_UX_rank,n_CTF);
X_weight_rC__ = zeros(n_k_p_r,n_CTF);
%%%%%%%%;
for nCTF=0:n_CTF-1;
if (flag_verbose); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: CTF_k_p_wk_: %0.16f',nCTF,n_CTF,fnorm(CTF_k_p_wk_))); end;
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
M_sub_k_p_wkM__ = M_k_p_wkM__(:,1+index_M_sub_);
%%%%;
% Prepare principal-modes. ;
%%%%;
UX_kn__ = eye(n_k_p_r,n_UX_rank);
X_weight_r_ = sqrt(weight_2d_k_p_r_);
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: X_weight_r_: %0.16f',nCTF,n_CTF,fnorm(X_weight_r_))); end;
if  flag_pm;
UX_kn__ = zeros(n_k_p_r,n_UX_rank);
X_weight_r_ = zeros(n_k_p_r,1);
[ ...
 X_kk__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_1( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M_sub ...
,M_sub_k_p_wkM__ ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
end;%if  flag_pm;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: UX_kn__: %0.16f',nCTF,n_CTF,fnorm(UX_kn__))); end;
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
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: M_sub_k_p_wkM__: %0.16f',nCTF,n_CTF,fnorm(M_sub_k_p_wkM__))); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: M_sub_k_q_wkM__: %0.16f',nCTF,n_CTF,fnorm(M_sub_k_q_wkM__))); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: CTF_M_sub_k_p_wkM__: %0.16f',nCTF,n_CTF,fnorm(CTF_M_sub_k_p_wkM__))); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: CTF_M_sub_k_q_wkM__: %0.16f',nCTF,n_CTF,fnorm(CTF_M_sub_k_q_wkM__))); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_VUXCTFM_sub_lwnM____: %0.16f',nCTF,n_CTF,fnorm(svd_VUXCTFM_sub_lwnM____))); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_VUXM_sub_lwnM____: %0.16f',nCTF,n_CTF,fnorm(svd_VUXM_sub_lwnM____))); end;
%%%%;
% Now calculate norms of the translated images. ;
%%%%;
tmp_t = tic();
UX_T_M_sub_l2_dM__ = tfpmh_UX_T_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_sub_lwnM____);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% tfpmh_UX_T_M_sub_l2_dm__1: %0.2fs',tmp_t)); end;
tmp_index_d0 = intersect(efind(FTK.delta_x_==0),efind(FTK.delta_y_==0)); assert(numel(tmp_index_d0)==1); %<-- should be zero-displacement. ;
UX_M_sub_l2_M_ = reshape(UX_T_M_sub_l2_dM__(1+tmp_index_d0,:),[n_M_sub,1]);
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: UX_M_sub_l2_M_: %0.16f',nCTF,n_CTF,fnorm(UX_M_sub_l2_M_))); end;
%%%%;
% Prepare UX_S_k_q_wnS__. ;
%%%%;
tmp_t = tic();
UX_S_k_q_wnS__ = reshape(permute(reshape(reshape(S_k_q_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),1+[0,2,1]),[n_w_max*pm_n_UX_rank,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% UX_S_k_q_wnS__: %0.2fs',tmp_t)); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: UX_S_k_q_wnS__: %0.16f',nCTF,n_CTF,fnorm(UX_S_k_q_wnS__))); end;
%%%%;
% Calculate Z_sub_dwSM____. ;
%%%%;
tmp_t = tic();
UX_S_k_q_nSw___ = permute(reshape(UX_S_k_q_wnS__,[n_w_max,n_UX_rank,n_S]),1+[1,2,0]);
tmp_t_sub = tic();
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: UX_S_k_q_nSw___: %0.16f',nCTF,n_CTF,fnorm(UX_S_k_q_nSw___))); end;
svd_VUXCTFM_sub_nMwl____ = permute(svd_VUXCTFM_sub_lwnM____,1+[2,3,1,0]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose); disp(sprintf(' %% svd_VUXCTFM_sub_nMwl____: %0.2fs',tmp_t_sub)); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_VUXCTFM_sub_nMwl____: %0.16f',nCTF,n_CTF,fnorm(svd_VUXCTFM_sub_nMwl____))); end;
tmp_t_sub = tic();
%svd_SVUXCTFM_sub_SMwl____ = zeros(n_S,n_M_sub,n_w_max,FTK.n_svd_l);
%for nl=0:FTK.n_svd_l-1; for nw=0:n_w_max-1;
%svd_SVUXCTFM_sub_SMwl____(:,:,1+nw,1+nl) = ctranspose(UX_S_k_q_nSw___(:,:,1+nw))*svd_VUXCTFM_sub_nMwl____(:,:,1+nw,1+nl);
%end;end;%for nw=0:n_w_max-1;for nl=0:FTK.n_svd_l-1;
svd_SVUXCTFM_sub_SMwl____ = pagemtimes(pagectranspose(UX_S_k_q_nSw___),svd_VUXCTFM_sub_nMwl____);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_sub_SMwl____: %0.6fs',tmp_t_sub)); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_SVUXCTFM_sub_SMwl____: %0.16f',nCTF,n_CTF,fnorm(svd_SVUXCTFM_sub_SMwl____))); end;
tmp_t_sub = tic();
%svd_SVUXCTFM_sub_lwSM____ = permute(ifft(permute(svd_SVUXCTFM_sub_SMwl____,1+[2,3,0,1]),[],1)*n_w_max,1+[1,0,2,3]);
svd_SVUXCTFM_sub_lwSM____ = ifft(permute(svd_SVUXCTFM_sub_SMwl____,1+[3,2,0,1]),[],1+1)*n_w_max;
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_SVUXCTFM_sub_lwSM____: %0.6fs',tmp_t_sub)); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_SVUXCTFM_sub_lwSM____: %0.16f',nCTF,n_CTF,fnorm(svd_SVUXCTFM_sub_lwSM____))); end;
tmp_t_sub = tic(); nop=0;
svd_USESVUXCTFM_sub_dwSM____ = reshape(FTK.svd_U_d_expiw_s__*reshape(svd_SVUXCTFM_sub_lwSM____,[FTK.n_svd_l,n_w_max*n_S*n_M_sub]),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]);
%svd_USESVUXCTFM_sub_dwSM____ = reshape(pagemtimes(FTK.svd_U_d_expiw_s__,svd_SVUXCTFM_sub_lwSM____),[FTK.n_delta_v,n_w_max,n_S,n_M_sub]);
tmp_t_sub = toc(tmp_t_sub); if (flag_verbose>1); disp(sprintf(' %% svd_USESVUXCTFM_sub_dwSM____: %0.6fs',tmp_t_sub)); end;
if (flag_verbose>1); disp(sprintf(' %% nCTF %d/%d: svd_USESVUXCTFM_sub_dwSM____: %0.16f',nCTF,n_CTF,fnorm(svd_USESVUXCTFM_sub_dwSM____))); end;
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% Z_sub_dwSM____: %0.2fs',tmp_t)); end;
%%%%;
% Store results. ;
%%%%;
UX_M_l2_M_(1+index_M_sub_) = UX_M_sub_l2_M_;
UX_T_M_l2_dM__(:,1+index_M_sub_) = UX_T_M_sub_l2_dM__;
Z_dwSM_ampm____(:,:,:,1+index_M_sub_) = svd_USESVUXCTFM_sub_dwSM____ ;
%%%%;
clear UX_kn__ X_weight_r_ CTF_k_p_wk_ index_M_sub_ UX_T_M_sub_l2_dM__ UX_M_sub_l2_M_ ;
clear UX_S_k_q_wnS__ ;
clear svd_VUXCTFM_sub_lwnM____ svd_VUXCTFM_sub_nMwl____ ;
clear svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_SMwl____ svd_SVUXCTFM_sub_lwSM____ ;
clear svd_USESVUXCTFM_sub_dwSM____ ;
%%%%%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
clear nCTF;
%%%%%%%%;
% check against dir_ascii. ;
%%%%%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_ascii = 2;
for nascii=0:n_ascii-1;
na=0;
if nascii==na; tmp_str_ascii = 'UX_T_M_l2_dM__'; tmp_local = UX_T_M_l2_dM__; tab_p=0; end; na=na+1;
if nascii==na; tmp_str_ascii = 'UX_M_l2_M_'; tmp_local = UX_M_l2_M_; tab_p=0; end; na=na+1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_ascii,tmp_local,'ascii',tmp_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
end;%for nascii=0:n_ascii-1;
%%%%%%%%;

%%%%%%%%;
% Now re-construct the template-norms, this time limited to radial principal-modes: ;
% <((R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_) * wUX_kn__),((R(+gamma_z)*CTF_k_p_wk_).*S_k_p_wk_) * wUX_kn__)> ;
% Note that this does yes involve collapsing onto principal-modes. ;
%%%%%%%%;
pm_n_k_p_r = pm_n_UX_rank; pm_n_w_max = n_w_max;
pm_n_w_ = pm_n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_sum = pm_n_k_p_r*pm_n_w_max;
UX_R_CTF_S_l2_wSC_quad___ = zeros(n_w_max,n_S,n_CTF);
SS_k_p_wkS__ = conj(S_k_p_wkS__).*S_k_p_wkS__;
SS_k_q_wkS__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,SS_k_p_wkS__),[n_w_sum,n_S]);
CC_k_p_wkC__ = conj(CTF_k_p_wkC__).*CTF_k_p_wkC__;
CC_k_q_wkC__ = reshape(interp_p_to_q_block_0(n_k_p_r,n_w_,n_w_sum,CC_k_p_wkC__),[n_w_sum,n_CTF]);
for nCTF=0:n_CTF-1;
UX_kn__ = UX_knC___(:,:,1+nCTF);
X_weight_r_ = X_weight_rC__(:,1+nCTF);
wUX_kn__ = diag(X_weight_r_)*UX_kn__;
UX_SS_k_q_wnS__ = reshape(permute(pagemtimes(transpose(wUX_kn__),permute(reshape(SS_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),1+[1,0,2])),1+[1,0,2]),[pm_n_w_sum,n_S]);
UX_CC_k_q_wn_ = reshape(reshape(CC_k_q_wkC__(:,1+nCTF),[n_w_max,n_k_p_r])*wUX_kn__,[pm_n_w_sum,1]);
UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF) = ifft(squeeze(sum(reshape(bsxfun(@times,conj(UX_CC_k_q_wn_),UX_SS_k_q_wnS__),[pm_n_w_max,pm_n_k_p_r,n_S]),1+1)));
clear UX_kn__ X_weight_r_ wUX_kn__ ;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% Now check one template and one CTF. ;
%%%%%%%%;
nS = max(0,min(n_S-1,round(n_S*1/5)));
nCTF = max(0,min(n_CTF-1,round(n_CTF*2/3)));
UX_kn__ = UX_knC___(:,:,1+nCTF); X_weight_r_ = X_weight_rC__(:,1+nCTF); wUX_kn__ = diag(X_weight_r_)*UX_kn__;
R_CTF_S_l2_w_quad_ = R_CTF_S_l2_wSC_quad___(:,1+nS,1+nCTF);
UX_R_CTF_S_l2_w_quad_ = UX_R_CTF_S_l2_wSC_quad___(:,1+nS,1+nCTF);
R_CTF_S_l2_w_qua2_ = zeros(n_w_max,1);
UX_R_CTF_S_l2_w_qua2_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
R_CTF_S_k_p_wk_ = S_k_p_wkS__(:,1+nS).*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),+gamma_z);
R_CTF_S_l2_w_qua2_(1+nw) = sum(conj(R_CTF_S_k_p_wk_).*R_CTF_S_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
UX_R_CTF_S_k_p_wn_ = reshape(reshape(S_k_p_wkS__(:,1+nS).*rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),+gamma_z),[n_w_max,n_k_p_r])*wUX_kn__,[pm_n_w_sum,1]);
UX_R_CTF_S_l2_w_qua2_(1+nw) = sum(conj(UX_R_CTF_S_k_p_wn_).*UX_R_CTF_S_k_p_wn_)/max(1,n_w_max);
end;%for nw=0:n_w_max-1;
clear UX_kn__ X_weight_r_ wUX_kn__ UX_R_CTF_S_k_p_wn_ ;
fnorm_disp(flag_verbose,'R_CTF_S_l2_w_qua2_',R_CTF_S_l2_w_qua2_,'R_CTF_S_l2_w_quad_',R_CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'R_CTF_S_l2_w_qua2_',R_CTF_S_l2_w_qua2_,'UX_R_CTF_S_l2_w_quad_',UX_R_CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'UX_R_CTF_S_l2_w_qua2_',UX_R_CTF_S_l2_w_qua2_,'R_CTF_S_l2_w_quad_',R_CTF_S_l2_w_quad_,' %%<-- should be zero');
fnorm_disp(flag_verbose,'UX_R_CTF_S_l2_w_qua2_',UX_R_CTF_S_l2_w_qua2_,'UX_R_CTF_S_l2_w_quad_',UX_R_CTF_S_l2_w_quad_,' %%<-- should be zero');
%%%%%%%%;

%%%%%%%%;
% Now estimate landscape of innerproducts across delta_x_ and delta_y_. ;
% Limited to a single image-template pair and a fixed nw. ;
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
Z_d_quad_ = zeros(n_delta_v,1);
Z_d_form_ = zeros(n_delta_v,1);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_phi = CTF_phi_C_(1+nCTF);
%%%%;
gamma_z = (2*pi*nw)/max(1,n_w_max);
RS_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_wk_,+gamma_z);
CTF_RS_k_p_wk_ = RS_k_p_wk_.*CTF_k_p_wk_;
CTF_RS_k_p_l2 = sum(conj(CTF_RS_k_p_wk_).*CTF_RS_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
T0M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+0.0,+0.0);
T0M_k_p_l2 = sum(conj(T0M_k_p_wk_).*T0M_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
%%;
for ndelta_v=0:n_delta_v-1;
%%;
delta_x = FTK.delta_x_(1+ndelta_v); delta_y = FTK.delta_y_(1+ndelta_v);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
CTF_RS_k_p_TM_k_p = sum(conj(CTF_RS_k_p_wk_).*TM_k_p_wk_.*weight_2d_wk_)*(2*pi)^2;
Z_d_quad_(1+ndelta_v) = CTF_RS_k_p_TM_k_p;
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
fnorm_disp(flag_verbose,'Z_d_quad_',Z_d_quad_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_ampm_',Z_d_ampm_,' %%<-- should be <1e-2');
fnorm_disp(flag_verbose,'Z_d_form_',Z_d_form_,'Z_d_quad_',Z_d_quad_,' %%<-- should be <1e-2');
%%%%%%%%;
% check against dir_ascii. ;
%%%%%%%%;
dir_ascii = '../dir_rangan_python/dir_ascii';
n_ascii = 3;
for nascii=0:n_ascii-1;
na=0;
if nascii==na; tmp_str_ascii = 'Z_d_quad_'; tmp_local = Z_d_quad_; tab_p=1; end; na=na+1;
if nascii==na; tmp_str_ascii = 'Z_d_ampm_'; tmp_local = Z_d_ampm_; tab_p=1; end; na=na+1;
if nascii==na; tmp_str_ascii = 'Z_d_form_'; tmp_local = Z_d_form_; tab_p=1; end; na=na+1;
fname_ascii = sprintf('%s/%s.ascii',dir_ascii,tmp_str_ascii);
if ~exist(fname_ascii,'file'); disp(sprintf(' %% Warning, %s not found',fname_ascii)); end;
if  exist(fname_ascii,'file');
if (flag_verbose>1); disp(sprintf(' %% %s found, loading',fname_ascii)); end;
fid = fopen(fname_ascii,'r');
if tab_p==0; tmp_ascii = reshape(cell2mat(textscan(fid,'%f')),size(tmp_local)); end;
if tab_p==1; tmp_ascii = reshape(cell2mat(textscan(fid,'(%f)')),size(tmp_local)); end;
fclose(fid);
fnorm_disp(flag_verbose,tmp_str_ascii,tmp_local,'ascii',tmp_ascii,' %% <-- should be <1e-6');
end;%if  exist(fname_ascii,'file');
end;%for nascii=0:n_ascii-1;
%%%%%%%%;
if (flag_disp>0);
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 3; np=0;
Zlim_ = prctile(real(Z_d_form_),[ 5,95],'all');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_form_),Zlim_); axisnotick; axis image; title('real(Z_d_form_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,imag(Z_d_form_),Zlim_); axisnotick; axis image; title('imag(Z_d_form_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_quad_),Zlim_); axisnotick; axis image; title('real(Z_d_quad_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,imag(Z_d_quad_),Zlim_); axisnotick; axis image; title('imag(Z_d_quad_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,real(Z_d_ampm_),Zlim_); axisnotick; axis image; title('real(Z_d_ampm_)','Interpreter','none');
subplot_t(p_row,p_col,1+np);np=np+1;cla; imagesc_d_FTK([],FTK.delta_x_,FTK.delta_y_,imag(Z_d_ampm_),Zlim_); axisnotick; axis image; title('imag(Z_d_ampm_)','Interpreter','none');
fname_fig_pre = sprintf('%s/test_transforms_Z_d_form_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%close(gcf);
end;%if (flag_disp>0);
%%%%%%%%;

%%%%%%%%;
% Now try and pick out optimal poses from displacements. ;
%%%%%%%%;
euler_azimu_b_opti_M_ = zeros(n_M,1);
euler_polar_a_opti_M_ = zeros(n_M,1);
euler_gamma_z_opti_M_ = zeros(n_M,1);
image_delta_x_opti_M_ = zeros(n_M,1);
image_delta_y_opti_M_ = zeros(n_M,1);
for nM=0:n_M-1;
if (flag_verbose>0); if mod(nM,32)==0; disp(sprintf(' %% nM %.3d/%.3d',nM,n_M)); end; end;
nCTF = index_nCTF_from_nM_(1+nM);
tmp_UX_T_M_l2_d_ = UX_T_M_l2_dM__(:,1+nM);
tmp_UX_R_CTF_S_l2_wS__ = flipud(UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
tmp_Z_dwS_ampm___ = Z_dwSM_ampm____(:,:,:,1+nM);
tmp_X_dwS_ampm___ = bsxfun(@rdivide,bsxfun(@rdivide,tmp_Z_dwS_ampm___,max(1e-12,reshape(sqrt(tmp_UX_R_CTF_S_l2_wS__),[1,n_w_max,n_S]))),max(1e-12,reshape(sqrt(tmp_UX_T_M_l2_d_),[n_delta_v,1,1])));
[~,tmp_ij] = max(tmp_X_dwS_ampm___,[],'all','linear'); tmp_index = tmp_ij-1;
tmp_index_tmp = tmp_index;
tmp_nd = mod(tmp_index_tmp,n_delta_v); tmp_index_tmp = (tmp_index_tmp - tmp_nd)/max(1,n_delta_v);
tmp_nw = mod(tmp_index_tmp,  n_w_max); tmp_index_tmp = (tmp_index_tmp - tmp_nw)/max(1,  n_w_max);
tmp_nS = mod(tmp_index_tmp,      n_S); tmp_index_tmp = (tmp_index_tmp - tmp_nS)/max(1,      n_S);
assert(fnorm(tmp_index_tmp)==0);
euler_azimu_b_opti_M_(1+nM) = viewing_azimu_b_S_(1+tmp_nS);
euler_polar_a_opti_M_(1+nM) = viewing_polar_a_S_(1+tmp_nS);
euler_gamma_z_opti_M_(1+nM) = (2*pi*tmp_nw)/max(1,n_w_max);
image_delta_x_opti_M_(1+nM) = FTK.delta_x_(1+tmp_nd);
image_delta_y_opti_M_(1+nM) = FTK.delta_y_(1+tmp_nd);
end;%for nM=0:n_M-1;
%%%%%%%%;
flag_check=0;
%%%%%%%%;
if flag_check;
euler_azimu_b_pose_M_ = zeros(n_M,1);
euler_polar_a_pose_M_ = zeros(n_M,1);
euler_gamma_z_pose_M_ = zeros(n_M,1);
image_delta_x_pose_M_ = zeros(n_M,1);
image_delta_y_pose_M_ = zeros(n_M,1);
index_true_pose_prct_M_ = zeros(n_M,1);
for nM=0:n_M-1;
if (flag_verbose>0); if mod(nM,32)==0; disp(sprintf(' %% nM %.3d/%.3d',nM,n_M)); end; end;
nCTF = index_nCTF_from_nM_(1+nM);
tmp_UX_T_M_l2_d_ = UX_T_M_l2_dM__(:,1+nM);
tmp_UX_R_CTF_S_l2_wS__ = flipud(UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
tmp_Z_dwS_ampm___ = Z_dwSM_ampm____(:,:,:,1+nM);
tmp_X_dwS_ampm___ = bsxfun(@rdivide,bsxfun(@rdivide,tmp_Z_dwS_ampm___,max(1e-12,reshape(sqrt(tmp_UX_R_CTF_S_l2_wS__),[1,n_w_max,n_S]))),max(1e-12,reshape(sqrt(tmp_UX_T_M_l2_d_),[n_delta_v,1,1])));
[tmp_X_dwS_sort_,tmp_ij_] = sort(tmp_X_dwS_ampm___(:),'ascend'); [~,tmp_ji_] = sort(tmp_ij_,'ascend'); tmp_index_ = tmp_ij_-1;
tmp_index_tmp_ = tmp_index_;
tmp_nd_ = mod(tmp_index_tmp_,n_delta_v); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nd_)/max(1,n_delta_v);
tmp_nw_ = mod(tmp_index_tmp_,  n_w_max); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nw_)/max(1,  n_w_max);
tmp_nS_ = mod(tmp_index_tmp_,      n_S); tmp_index_tmp_ = (tmp_index_tmp_ - tmp_nS_)/max(1,      n_S);
assert(fnorm(tmp_index_tmp_)==0);
opt_nd = tmp_nd_(end); opt_delta_x = FTK.delta_x_(1+opt_nd); opt_delta_y = FTK.delta_y_(1+opt_nd);
opt_nw = tmp_nw_(end); opt_gamma_z = gamma_z_(1+opt_nw);
opt_nS = tmp_nS_(end); opt_polar_a = viewing_polar_a_S_(1+opt_nS); opt_azimu_b = viewing_azimu_b_S_(1+opt_nS);
euler_azimu_b_pose_M_(1+nM) = opt_azimu_b;
euler_polar_a_pose_M_(1+nM) = opt_polar_a;
euler_gamma_z_pose_M_(1+nM) = opt_gamma_z;
image_delta_x_pose_M_(1+nM) = opt_delta_x;
image_delta_y_pose_M_(1+nM) = opt_delta_y;
tru_nd = index_nd_from_nM_(1+nM); tru_delta_x = image_delta_x_true_M_(1+nM); tru_delta_y = image_delta_y_true_M_(1+nM);
tru_nw = index_nw_from_nM_(1+nM); tru_gamma_z = gamma_z_(1+tru_nw);
tru_nS = index_nS_from_nM_(1+nM); tru_polar_a = euler_polar_a_true_M_(1+nM); tru_azimu_b = euler_azimu_b_true_M_(1+nM);
index_true_pose = tru_nd + (tru_nw + (tru_nS)*n_w_max)*n_delta_v ;
index_true_pose_prct = tmp_ji_(1+index_true_pose)/max(1,n_delta_v*n_w_max*n_S);
index_true_pose_prct_M_(1+nM) = index_true_pose_prct;
if (flag_verbose>1);
disp(sprintf(' %% nM %.3d/%.3d: ',nM,n_M));
disp(sprintf(' %% opt_nS %.3d tru_nS %.3d opt_nw %.3d tru_nw %.3d opt_nd %.3d tru_nd %.3d',opt_nS,tru_nS,opt_nw,tru_nw,opt_nd,tru_nd));
disp(sprintf(' %% opt_azimu_b %+0.2f opt_polar_a %+0.2f opt_gamma_z %+0.2f opt_delta_x %+0.2f opt_delta_y %+0.2f',opt_azimu_b,opt_polar_a,opt_gamma_z,opt_delta_x,opt_delta_y));
disp(sprintf(' %% tru_azimu_b %+0.2f tru_polar_a %+0.2f tru_gamma_z %+0.2f tru_delta_x %+0.2f tru_delta_y %+0.2f',tru_azimu_b,tru_polar_a,tru_gamma_z,tru_delta_x,tru_delta_y));
end;%if (flag_verbose>0);
end;%for nM=0:n_M-1;
%%%%%%%%;
fnorm_disp(flag_verbose,'euler_azimu_b_opti_M_',euler_azimu_b_opti_M_,'euler_azimu_b_pose_M_',euler_azimu_b_pose_M_);
fnorm_disp(flag_verbose,'euler_polar_a_opti_M_',euler_polar_a_opti_M_,'euler_polar_a_pose_M_',euler_polar_a_pose_M_);
fnorm_disp(flag_verbose,'euler_gamma_z_opti_M_',euler_gamma_z_opti_M_,'euler_gamma_z_pose_M_',euler_gamma_z_pose_M_);
fnorm_disp(flag_verbose,'image_delta_x_opti_M_',image_delta_x_opti_M_,'image_delta_x_pose_M_',image_delta_x_pose_M_);
fnorm_disp(flag_verbose,'image_delta_y_opti_M_',image_delta_y_opti_M_,'image_delta_y_pose_M_',image_delta_y_pose_M_);
%%%%%%%%;
end;%if flag_check;
%%%%%%%%;

%%%%%%%%;
% Now compare with results from tfpmh_Z_wSM___12. ;
%%%%%%%%;
nCTF = max(0,min(n_CTF-1,round(n_CTF*2/3)));
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
pm_n_UX_rank = n_UX_rank;
pm_UX_kn__ = UX_knC___(:,:,1+nCTF);
pm_X_weight_r_ = X_weight_rC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
M_sub_k_p_wkM__ = M_k_p_wkM__(:,1+index_M_sub_);
tmp_parameter=struct('type','parameter');
tmp_parameter.flag_verbose = 1;
tmp_parameter.flag_optimize_over_gamma_z = 0;
[ ...
 tmp_parameter ...
,tmp_Z_wSM___ ...
,tmp_UX_R_CTF_S_l2_wS__ ...
,tmp_R_CTF_S_l2_wS__ ...
,tmp_UX_T_M_l2_dM__ ...
,tmp_UX_M_l2_M_ ...
,tmp_X_wSM___ ...
,tmp_delta_x_wSM___ ...
,tmp_delta_y_wSM___ ...
,tmp_gamma_z_wSM___ ...
,tmp_index_sub_wSM___ ...
,tmp_Z_dwSM____ ...
,tmp_X_dwSM____ ...
] = ...
tfpmh_Z_wSM___12( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_M_sub ...
,M_sub_k_p_wkM__ ...
,CTF_k_p_wk_ ...
,n_S ...
,S_k_p_wkS__ ...
,pm_n_UX_rank ...
,pm_UX_kn__ ...
,pm_X_weight_r_ ...
,FTK ...
);
fnorm_disp(flag_verbose,'tmp_R_CTF_S_l2_wS__',tmp_R_CTF_S_l2_wS__,'R_CTF_S_l2_wSC_quad___(:,:,1+nCTF)',R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
fnorm_disp(flag_verbose,'tmp_UX_R_CTF_S_l2_wS__',tmp_UX_R_CTF_S_l2_wS__,'UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF)',UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
fnorm_disp(flag_verbose,'tmp_UX_T_M_l2_dM__',tmp_UX_T_M_l2_dM__,'UX_T_M_l2_dM__(:,1+index_M_sub_)',UX_T_M_l2_dM__(:,1+index_M_sub_));
fnorm_disp(flag_verbose,'tmp_UX_M_l2_M_',tmp_UX_M_l2_M_,'UX_M_l2_M_(1+index_M_sub_)',UX_M_l2_M_(1+index_M_sub_));
tmp_Z_dwSM_ampm____ = Z_dwSM_ampm____(:,:,:,1+index_M_sub_);
n_wSM = n_w_max*n_S*n_M_sub;
tmp_index_all_wSM_ = tmp_index_sub_wSM___(:) + reshape([0:n_wSM-1],size(tmp_index_sub_wSM___(:)))*n_delta_v;
fnorm_disp(flag_verbose,'tmp_Z_wSM___(:)',tmp_Z_wSM___(:),'tmp_Z_dwSM_ampm____(1+tmp_index_all_wSM_)',tmp_Z_dwSM_ampm____(1+tmp_index_all_wSM_));
tmp_Z_errrel = 0.0;
tmp_X_errrel = 0.0;
for nM_sub=0:n_M_sub-1;
nM = index_M_sub_(1+nM_sub);
nw = max(0,min(n_w_max-1,floor(n_w_max*rand())));
nS = max(0,min(n_S-1,floor(n_S*rand())));
ndelta = tmp_index_sub_wSM___(1+nw,1+nS,1+nM_sub);
tmp_Z_ampm = tmp_Z_dwSM_ampm____(1+ndelta,1+nw,1+nS,1+nM_sub);
tmp_X_ampm = tmp_Z_ampm/max(1e-12,sqrt(UX_R_CTF_S_l2_wSC_quad___(1+nw,1+nS,1+nCTF)))/max(1e-12,sqrt(tmp_UX_T_M_l2_dM__(1+ndelta,1+nM_sub)));
tmp_Z_tfpm = tmp_Z_wSM___(1+nw,1+nS,1+nM_sub);
tmp_X_tfpm = tmp_X_wSM___(1+nw,1+nS,1+nM_sub);
tmp_Z_errrel = tmp_Z_errrel + fnorm(tmp_Z_ampm-tmp_Z_tfpm)/max(1e-12,fnorm(tmp_Z_ampm));
tmp_X_errrel = tmp_X_errrel + fnorm(tmp_X_ampm-tmp_X_tfpm)/max(1e-12,fnorm(tmp_X_ampm));
end;%for nM_sub=0:n_M_sub-1;
disp(sprintf(' %% tmp_Z_errrel: %0.16f',tmp_Z_errrel));
disp(sprintf(' %% tmp_X_errrel: %0.16f',tmp_X_errrel));
%%%%%%%%;
nCTF = max(0,min(n_CTF-1,round(n_CTF*1/3)));
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
pm_n_UX_rank = n_UX_rank;
pm_UX_kn__ = UX_knC___(:,:,1+nCTF);
pm_X_weight_r_ = X_weight_rC__(:,1+nCTF);
index_M_sub_ = efind(index_nCTF_from_nM_==nCTF); n_M_sub = numel(index_M_sub_);
M_sub_k_p_wkM__ = M_k_p_wkM__(:,1+index_M_sub_);
tmp_parameter=struct('type','parameter');
tmp_parameter.flag_verbose = 1;
tmp_parameter.flag_optimize_over_gamma_z = 1;
[ ...
 tmp_parameter ...
,tmp_Z_SM__ ...
,tmp_UX_R_CTF_S_l2_wS__ ...
,tmp_R_CTF_S_l2_wS__ ...
,tmp_UX_T_M_l2_dM__ ...
,tmp_UX_M_l2_M_ ...
,tmp_X_SM__ ...
,tmp_delta_x_SM__ ...
,tmp_delta_y_SM__ ...
,tmp_gamma_z_SM__ ...
,tmp_index_sub_SM__ ...
,tmp_Z_dwSM____ ...
,tmp_X_dwSM____ ...
] = ...
tfpmh_Z_wSM___12( ...
 tmp_parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,n_M_sub ...
,M_sub_k_p_wkM__ ...
,CTF_k_p_wk_ ...
,n_S ...
,S_k_p_wkS__ ...
,pm_n_UX_rank ...
,pm_UX_kn__ ...
,pm_X_weight_r_ ...
,FTK ...
);
fnorm_disp(flag_verbose,'tmp_R_CTF_S_l2_wS__',tmp_R_CTF_S_l2_wS__,'R_CTF_S_l2_wSC_quad___(:,:,1+nCTF)',R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
fnorm_disp(flag_verbose,'tmp_UX_R_CTF_S_l2_wS__',tmp_UX_R_CTF_S_l2_wS__,'UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF)',UX_R_CTF_S_l2_wSC_quad___(:,:,1+nCTF));
fnorm_disp(flag_verbose,'tmp_UX_T_M_l2_dM__',tmp_UX_T_M_l2_dM__,'UX_T_M_l2_dM__(:,1+index_M_sub_)',UX_T_M_l2_dM__(:,1+index_M_sub_));
fnorm_disp(flag_verbose,'tmp_UX_M_l2_M_',tmp_UX_M_l2_M_,'UX_M_l2_M_(1+index_M_sub_)',UX_M_l2_M_(1+index_M_sub_));
tmp_Z_dwSM_ampm____ = Z_dwSM_ampm____(:,:,:,1+index_M_sub_);
n_dw = n_delta_v*n_w_max; n_SM = n_S*n_M_sub;
tmp_index_all_SM_ = tmp_index_sub_SM__(:) + reshape([0:n_SM-1],size(tmp_index_sub_SM__(:)))*n_dw;
fnorm_disp(flag_verbose,'tmp_Z_SM__(:)',tmp_Z_SM__(:),'tmp_Z_dwSM_ampm____(1+tmp_index_all_SM_)',tmp_Z_dwSM_ampm____(1+tmp_index_all_SM_));
tmp_Z_errrel = 0.0;
tmp_X_errrel = 0.0;
for nM_sub=0:n_M_sub-1;
nM = index_M_sub_(1+nM_sub);
nS = max(0,min(n_S-1,floor(n_S*rand())));
ndw = tmp_index_sub_SM__(1+nS,1+nM_sub);
ndelta = mod(ndw,n_delta_v);
nw = (ndw-ndelta)/max(1,n_delta_v);
tmp_Z_ampm = tmp_Z_dwSM_ampm____(1+ndelta,1+nw,1+nS,1+nM_sub);
tmp_X_ampm = tmp_Z_ampm/max(1e-12,sqrt(UX_R_CTF_S_l2_wSC_quad___(1+nw,1+nS,1+nCTF)))/max(1e-12,sqrt(tmp_UX_T_M_l2_dM__(1+ndelta,1+nM_sub)));
tmp_Z_tfpm = tmp_Z_SM__(1+nS,1+nM_sub);
tmp_X_tfpm = tmp_X_SM__(1+nS,1+nM_sub);
tmp_Z_errrel = tmp_Z_errrel + fnorm(tmp_Z_ampm-tmp_Z_tfpm)/max(1e-12,fnorm(tmp_Z_ampm));
tmp_X_errrel = tmp_X_errrel + fnorm(tmp_X_ampm-tmp_X_tfpm)/max(1e-12,fnorm(tmp_X_ampm));
end;%for nM_sub=0:n_M_sub-1;
disp(sprintf(' %% tmp_Z_errrel: %0.16f',tmp_Z_errrel));
disp(sprintf(' %% tmp_X_errrel: %0.16f',tmp_X_errrel));
%%%%%%%%;

disp('returning before 3d-reconstructions'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now perform various 3d-reconstructions. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
a_k_p_relerror = @(a_k_p_test_) sqrt(sum(abs(a_k_p_form_ - a_k_p_test_).^2.*weight_3d_k_p_qk_)./max(1e-12,sum(abs(a_k_p_form_).^2.*weight_3d_k_p_qk_)));
a_k_Y_relerror = @(a_k_Y_test_) sqrt(sum(abs(a_k_Y_form_ - a_k_Y_test_).^2.*weight_Y_)./max(1e-12,sum(abs(a_k_Y_form_).^2.*weight_Y_))); 

%%%%%%%%;
% First we center the images via the appropriate image-specific translation: ;
%%%%%%%%;
TM_k_p_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
delta_x = image_delta_x_true_M_(1+nM); delta_y = image_delta_y_true_M_(1+nM);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
TM_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,+delta_x,+delta_y);
TM_k_p_wkM__(:,1+nM) = TM_k_p_wk_;
end;%for nM=0:n_M-1;
%%%%%%%%;
parameter_kappa = struct('type','parameter');
parameter_kappa.kernel_basic_l_max_use = l_max_max;
parameter_kappa.kernel_basic_qref_k_eq_d_double = k_eq_d/max(1e-12,k_p_r_max);
parameter_kappa.flag_recalc_qref_from_data = 1;
parameter_kappa.flag_recalc_dtau_qref_from_data = 0;
parameter_kappa.flag_recalc_dtau_dtau_qref_from_data = 0;
KAPPA = [];
weight_imagecount_M_ = ones(n_M,1);
tmp_t = tic();
[ ...
 parameter_kappa ...
,KAPPA ...
,a_from_M_C2M0_k_p_qk__ ...
,a_from_M_C1M1_k_p_qk__ ...
,a_from_M_C0M2_k_p_qk__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,KAPPA ...
,n_w_max ...
,n_M ...
,weight_imagecount_M_ ...
,+euler_polar_a_true_M_ ...
,+euler_azimu_b_true_M_ ...
,+euler_gamma_z_true_M_ ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,TM_k_p_wkM__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
a_k_p_from_M_ = reshape(a_from_M_C1M1_k_p_qk__./max(1e-12,a_from_M_C2M0_k_p_qk__),[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_M_',a_k_p_from_M_);
if (flag_verbose>0); disp(sprintf(' %% a_k_p_form_ vs a_k_p_from_M_: volumetric-relative-error: %0.16f',a_k_p_relerror(a_k_p_from_M_))); end;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter_kappa ...
,~ ...
,a_from_S_C2M0_k_p_qk__ ...
,a_from_S_C1M1_k_p_qk__ ...
,a_from_S_C0M2_k_p_qk__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,[] ...
,n_w_max ...
,n_S ...
,ones(n_S,1) ...
,+viewing_polar_a_S_ ...
,+viewing_azimu_b_S_ ...
,+zeros(n_S,1) ...
,[] ...
,[] ...
,[] ...
,n_k_p_r ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,[] ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
a_k_p_from_S_ = reshape(a_from_S_C1M1_k_p_qk__./max(1e-12,a_from_S_C2M0_k_p_qk__),[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_S_',a_k_p_from_S_);
if (flag_verbose>0); disp(sprintf(' %% a_k_p_form_ vs a_k_p_from_S_: volumetric-relative-error: %0.16f',a_k_p_relerror(a_k_p_from_S_))); end;
%%%%%%%%;
qbp_eps = 1e-6;
[ ...
 a_k_Y_qbpr_M_ ...
,~ ...
,a_k_p_qbpr_M_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_true_M_ ...
,euler_azimu_b_true_M_ ...
,euler_gamma_z_true_M_ ...
,+image_delta_x_true_M_ ...
,+image_delta_y_true_M_ ...
);
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_qbpr_M_',a_k_Y_qbpr_M_);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_qbpr_M_: volumetric-relative-error: %0.16f',a_k_Y_relerror(a_k_Y_qbpr_M_))); end;
%%%%%%%%;
qbp_eps = 1e-6;
[ ...
 a_k_Y_qbpr_S_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,zeros(n_S,1) ...
);
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_qbpr_S_',a_k_Y_qbpr_S_);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_qbpr_S_: volumetric-relative-error: %0.16f',a_k_Y_relerror(a_k_Y_qbpr_S_))); end;
%%%%%%%%;
n_order = 11;
tmp_t = tic();
[ ...
 a_k_Y_lsqr_S_ ...
] = ...
cg_lsq_6( ...
 n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,zeros(n_S,1) ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% cg_lsq_6: %0.2fs',tmp_t)); end;
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_lsqr_S_',a_k_Y_lsqr_S_);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_lsqr_S_: volumetric-relative-error: %0.16f',a_k_Y_relerror(a_k_Y_lsqr_S_))); end;
%%%%%%%%;
n_order = 11;
tmp_t = tic();
[ ...
 a_k_Y_lsqr_M_ ...
] = ...
cg_lsq_6( ...
 n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_M ...
,M_k_p_wkM__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,euler_polar_a_true_M_ ...
,euler_azimu_b_true_M_ ...
,euler_gamma_z_true_M_ ...
,+image_delta_x_true_M_ ...
,+image_delta_y_true_M_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% cg_lsq_6: %0.2fs',tmp_t)); end;
fnorm_disp(flag_verbose,'a_k_Y_form_',a_k_Y_form_,'a_k_Y_lsqr_M_',a_k_Y_lsqr_M_);
if (flag_verbose>0); disp(sprintf(' %% a_k_Y_form_ vs a_k_Y_lsqr_M_: volumetric-relative-error: %0.16f',a_k_Y_relerror(a_k_Y_lsqr_M_))); end;
%%%%%%%%;

%%%%%%%%;
% Now perform a 3d-reconstruction using kappa_basic_apply_4, ;
% but limited to a fixed number of principal-modes. ;
% Note that this procedure is only efficient if: ;
%  (i) the CTFs are isotropic, and ;
% (ii) the number of eigen-CTFs is low. ;
% so as a first demonstration we will use only sythetic-images (i.e., the templates). ;
%%%%%%%%;

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
,n_S ...
,S_k_p_wkS__ ...
);
[UX_kn__,SX_k__,VX_kn__] = svds(X_kk__,n_UX_rank);
end;%if  flag_pm;
%%%%;

%%%%;
% Prepare UX_S_k_p_wkS__. ;
%%%%;
tmp_t = tic();
S_k_p_wSk___ = reshape(permute(reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),1+[0,2,1]),[n_w_max,n_S,n_k_p_r]);
UX_S_k_p_wnS__ = reshape(permute(reshape(reshape(S_k_p_wSk___,[n_w_max*n_S,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max,n_S,pm_n_UX_rank]),1+[0,2,1]),[n_w_max*pm_n_UX_rank,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% UX_S_k_p_wnS__: %0.2fs',tmp_t)); end;
%%%%;
parameter_kappa = struct('type','parameter');
parameter_kappa.kernel_basic_l_max_use = l_max_max;
parameter_kappa.kernel_basic_qref_k_eq_d_double = k_eq_d/max(1e-12,k_p_r_max);
parameter_kappa.flag_recalc_qref_from_data = 1;
parameter_kappa.flag_recalc_dtau_qref_from_data = 0;
parameter_kappa.flag_recalc_dtau_dtau_qref_from_data = 0;
KAPPA = [];
weight_imagecount_S_ = ones(n_S,1);
tmp_t = tic();
[ ...
 parameter_kappa ...
,KAPPA ...
,UX_a_from_UX_S_C2M0_k_p_qn__ ...
,UX_a_from_UX_S_C1M1_k_p_qn__ ...
,UX_a_from_UX_S_C0M2_k_p_qn__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,KAPPA ...
,n_w_max ...
,n_S ...
,ones(n_S,1) ...
,+viewing_polar_a_S_ ...
,+viewing_azimu_b_S_ ...
,+zeros(n_S,1) ...
,[] ...
,[] ...
,[] ...
,n_UX_rank ...
,UX_S_k_p_wnS__ ...
,[] ...
,[] ...
,[] ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
%%%%;
UX_a_from_UX_S_k_p_qn__ = reshape(UX_a_from_UX_S_C1M1_k_p_qn__./max(1e-12,UX_a_from_UX_S_C2M0_k_p_qn__),[qref_n_shell,n_UX_rank]);
a_from_UX_S_k_p_qk__ = zeros(qref_n_shell,n_k_p_r);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
a_from_UX_S_k_p_qk__(:,1+nk_p_r) = a_from_UX_S_k_p_qk__(:,1+nk_p_r) + UX_kn__(1+nk_p_r,1+nUX_rank)/max(1e-12,X_weight_r_(1+nk_p_r))*UX_a_from_UX_S_k_p_qn__(:,1+nUX_rank);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
a_from_UX_S_k_p_qk__ = UX_a_from_UX_S_k_p_qn__*transpose(UX_kn__)*diag(1./max(1e-12,X_weight_r_));
a_k_p_from_UX_S_ = reshape(a_from_UX_S_k_p_qk__,[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_UX_S_',a_k_p_from_UX_S_);
a_from_UX_S_k_p_qk__ = UX_a_from_UX_S_k_p_qn__*transpose(UX_kn__)*diag(1./max(1e-12,X_weight_r_));
a_k_p_from_UX_S_ = reshape(a_from_UX_S_k_p_qk__,[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_UX_S_',a_k_p_from_UX_S_);
if (flag_verbose>0); disp(sprintf(' %% a_k_p_form_ vs a_k_p_from_UX_S_: volumetric-relative-error: %0.16f',a_k_p_relerror(a_k_p_from_UX_S_))); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now we will perform something similar, ;
% except this time we will use synthetic-images ;
% which are associated with one of three isotropic ICTFs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
n_ICTF = 4;
ICTF_alpha_C_ = transpose(linspace(0.05,2*2.4048/k_p_r_max,n_ICTF)); %<-- note that first root of besselj(0,.) is ~2.4048 ;
ICTF_k_p_r_kC__ = zeros(n_k_p_r,n_ICTF);
for nICTF=0:n_ICTF-1;
if ~flag_CTF; ICTF_k_p_r_kC__(:,1+nICTF) = ones(n_k_p_r,1); end;
if  flag_CTF; ICTF_k_p_r_kC__(:,1+nICTF) = besselj(0,ICTF_alpha_C_(1+nICTF)*k_p_r_); end;
end;%for nICTF=0:n_ICTF-1;
n_N = n_ICTF*n_S;
N_k_p_wkN__ = zeros(n_w_sum,n_N);
index_nICTF_from_nN_ = zeros(n_N,1);
euler_azimu_b_true_N_ = zeros(n_N,1);
euler_polar_a_true_N_ = zeros(n_N,1);
euler_gamma_z_true_N_ = zeros(n_N,1);
for nICTF=0:n_ICTF-1;
ICTF_k_p_r_k_ = ICTF_k_p_r_kC__(:,1+nICTF);
N_k_p_wkN__(:,1+nICTF*n_S+[0:n_S-1]) = reshape(bsxfun(@times,reshape(S_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(ICTF_k_p_r_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
index_nICTF_from_nN_(1+nICTF*n_S+[0:n_S-1]) = nICTF;
euler_azimu_b_true_N_(1+nICTF*n_S+[0:n_S-1]) = viewing_azimu_b_S_;
euler_polar_a_true_N_(1+nICTF*n_S+[0:n_S-1]) = viewing_polar_a_S_;
euler_gamma_z_true_N_(1+nICTF*n_S+[0:n_S-1]) = zeros(n_S,1);
end;%for nICTF=0:n_ICTF-1;
%%%%%%%%;
a_from_UX_N_C1M1_k_p_qkC___ = zeros(qref_n_shell,n_k_p_r,n_CTF);
a_from_UX_N_C2M0_k_p_qkC___ = zeros(qref_n_shell,n_k_p_r,n_CTF);
%%%%%%%%%%%%%%%%;
for nICTF=0:n_ICTF-1;
%%%%%%%%%%%%%%%%;
ICTF_k_p_r_k_ = ICTF_k_p_r_kC__(:,1+nICTF);
index_N_sub_ = efind(index_nICTF_from_nN_==nICTF);
n_N_sub = numel(index_N_sub_);
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
,n_N_sub ...
,N_k_p_wkN__(:,1+index_N_sub_) ...
);
[tmp_UX_kn__,tmp_SX_k__,tmp_VX_kn__] = svds(X_kk__,n_UX_rank);
UX_kn__(:,1:n_UX_rank) = tmp_UX_kn__(:,1:n_UX_rank);
end;%if  flag_pm;
%%%%;
% Prepare UX_N_k_p_wkN__. ;
%%%%;
tmp_t = tic();
N_sub_k_p_wNk___ = reshape(permute(reshape(N_k_p_wkN__(:,1+index_N_sub_),[n_w_max,n_k_p_r,n_N_sub]),1+[0,2,1]),[n_w_max,n_N_sub,n_k_p_r]);
UX_N_sub_k_p_wnN__ = reshape(permute(reshape(reshape(N_sub_k_p_wNk___,[n_w_max*n_N_sub,n_k_p_r])*diag(X_weight_r_)*UX_kn__,[n_w_max,n_N_sub,pm_n_UX_rank]),1+[0,2,1]),[n_w_max*pm_n_UX_rank,n_N_sub]);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% UX_N_sub_k_p_wnN__: %0.2fs',tmp_t)); end;
%%%%;
parameter_kappa = struct('type','parameter');
parameter_kappa.kernel_basic_l_max_use = l_max_max;
parameter_kappa.kernel_basic_qref_k_eq_d_double = k_eq_d/max(1e-12,k_p_r_max);
parameter_kappa.flag_recalc_qref_from_data = 1;
parameter_kappa.flag_recalc_dtau_qref_from_data = 0;
parameter_kappa.flag_recalc_dtau_dtau_qref_from_data = 0;
KAPPA = [];
weight_imagecount_N_sub_ = ones(n_N_sub,1);
tmp_t = tic();
[ ...
 parameter_kappa ...
,KAPPA ...
,UX_a_from_UX_N_sub_C2M0_k_p_qn__ ...
,UX_a_from_UX_N_sub_C1M1_k_p_qn__ ...
,UX_a_from_UX_N_sub_C0M2_k_p_qn__ ...
] = ...
kappa_basic_apply_4( ...
 parameter_kappa ...
,KAPPA ...
,n_w_max ...
,n_N_sub ...
,weight_imagecount_N_sub_ ...
,+euler_polar_a_true_N_(1+index_N_sub_) ...
,+euler_azimu_b_true_N_(1+index_N_sub_) ...
,+euler_gamma_z_true_N_(1+index_N_sub_) ...
,[] ...
,[] ...
,[] ...
,n_UX_rank ...
,UX_N_sub_k_p_wnN__ ...
,[] ...
,[] ...
,[] ...
,qref_k_eq_d ...
,qref_n_shell ...
,qref_azimu_b_shell_ ...
,qref_polar_a_shell_ ...
,qref_weight_shell_ ...
,qref_k_c_0_shell_ ...
,qref_k_c_1_shell_ ...
,qref_k_c_2_shell_ ...
,qref_n_polar_a ...
,qref_polar_a_ ...
,qref_n_azimu_b_ ...
);
tmp_t = toc(tmp_t); if (flag_verbose); disp(sprintf(' %% kappa_basic_apply_4: %0.2fs',tmp_t)); end;
%%%%;
a_from_UX_N_C1M1_k_p_qk__ = UX_a_from_UX_N_sub_C1M1_k_p_qn__*transpose(UX_kn__)*diag(1./max(1e-12,X_weight_r_))*diag(ICTF_k_p_r_k_.^1);
a_from_UX_N_C2M0_k_p_qk__ = repmat(UX_a_from_UX_N_sub_C2M0_k_p_qn__(:,1),[1,n_k_p_r])*diag(ICTF_k_p_r_k_.^2);
a_from_UX_N_C1M1_k_p_qkC___(:,:,1+nICTF) = a_from_UX_N_C1M1_k_p_qk__;
a_from_UX_N_C2M0_k_p_qkC___(:,:,1+nICTF) = a_from_UX_N_C2M0_k_p_qk__;
clear UX_kn__ tmp_UX_kn__ tmp_SX_k__ tmp_VX_kn__ X_kk__ ;
clear N_sub_k_p_wNk___ UX_N_sub_k_p_wnN__ ;
clear KAPPA UX_a_from_UX_N_sub_C2M0_k_p_qn__ UX_a_from_UX_N_sub_C1M1_k_p_qn__ UX_a_from_UX_N_sub_C0M2_k_p_qn__ ;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%for nICTF=0:n_ICTF-1;
%%%%%%%%%%%%%%%%;
a_from_UX_N_k_p_qk__ = sum(a_from_UX_N_C1M1_k_p_qkC___,3)./max(1e-12,sum(max(1e-12,a_from_UX_N_C2M0_k_p_qkC___),3));
a_k_p_from_UX_N_ = reshape(a_from_UX_N_k_p_qk__,[qref_n_shell*n_k_p_r,1]);
fnorm_disp(flag_verbose,'a_k_p_form_',a_k_p_form_,'a_k_p_from_UX_N_',a_k_p_from_UX_N_);
if (flag_verbose>0); disp(sprintf(' %% a_k_p_form_ vs a_k_p_from_UX_N_: volumetric-relative-error: %0.16f',a_k_p_relerror(a_k_p_from_UX_N_))); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

return;





