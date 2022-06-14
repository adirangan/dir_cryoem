%%%%%%%%;
% test 1d-version of cryo alignment. ;
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
nf=0;
flag_recalc = 1;
flag_plot = 0;
flag_replot = 1;
flag_center_volume = 1;
flag_center_image = 0;
flag_invert = 0;
tolerance_master = 1e-2;
global_n_iteration = 16;
global_rseed = 0;
rng(global_rseed);
Pixel_A = 4.0;
%%%%%%%%;

global_flag_volume_type = 4;
snr_per_A = 0.01; %<-- decrease to increase noise. ;
f_rand_0in = 0.0;
n_M_factor_0in = 4; %<-- n_M = n_w_max * n_M_factor_0in. ;
k_p_r_max_factor_0in = 2; %<-- k_p_r_max = 48 * k_p_r_max_factor_0in. ;
n_x_u_factor_0in = 4; %<-- n_x_u = 64 * n_x_u_factor_0in. ;
n_pick_0in = 12; %<-- strong overlap at n_pick_0in = ceil(sqrt(n_w_max)). ;
local_rseed = 0;

%%%%%%%%;
fname_prefix = 'SR0A_0';
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
%%%%%%%%;

str_vol = sprintf('v%.2d',global_flag_volume_type);
str_snr = sprintf('s%.2d',round(10*-log10(snr_per_A)));
str_fra = sprintf('f%.2d',round(100*f_rand_0in));
str_nMf = sprintf('M%.2d',n_M_factor_0in);
str_npi = sprintf('p%.2d',n_pick_0in);
str_rng = sprintf('r%.2d',local_rseed);
string_infix = sprintf('%s%s%s%s%s%s' ...
,str_vol ...
,str_snr ...
,str_nMf ...
,str_npi ...
,str_fra ...
,str_rng ...
);

verbose=1;
%%%%%%%%;
% First create consensus volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u = n_x_u_factor_0in*64;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
[x_u_0__,x_u_1__] = ndgrid(x_u_0_,x_u_1_); n_xx_u = n_x_u^2;
%%%%%%%%;
if global_flag_volume_type==0;
x_sigma = 2*1/n_x_u;
a_x_u_ = zeros(n_x_u,n_x_u);
x_s_0 = 0; x_s_1 = 0;
x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
a_x_u_ = a_x_u_ +  1/(sqrt(2*pi)*x_sigma).^2 * exp(-x_u_r__.^2/(2*x_sigma^2));
end;%if global_flag_volume_type==0;
%%%%%%%%;
if global_flag_volume_type==1;
a_x_u_ = zeros(n_x_u,n_x_u);
tmp_ = [ ...
	      blockletter_0('F',[n_x_u_factor_0in*6,n_x_u_factor_0in*1,1]) ...
             ,zeros(n_x_u_factor_0in*6*7,n_x_u_factor_0in*3) ...
	     ,blockletter_0('I',[n_x_u_factor_0in*6,n_x_u_factor_0in*1,1]) ...
             ,zeros(n_x_u_factor_0in*6*7,n_x_u_factor_0in*3) ...
	     ,blockletter_0('N',[n_x_u_factor_0in*6,n_x_u_factor_0in*1,1]) ...
             ,zeros(n_x_u_factor_0in*6*7,n_x_u_factor_0in*3) ...
	     ,blockletter_0('Y',[n_x_u_factor_0in*6,n_x_u_factor_0in*1,1]) ...
             ,zeros(n_x_u_factor_0in*6*7,n_x_u_factor_0in*3) ...
	     ,blockletter_0('U',[n_x_u_factor_0in*6,n_x_u_factor_0in*1,1]) ...
	 ];
a_x_u_( ceil(n_x_u/2 - size(tmp_,1)/2) + [0:size(tmp_,1)-1] ,  ceil(n_x_u/2 - size(tmp_,2)/2) + [0:size(tmp_,2)-1] ) = tmp_;
end;%if global_flag_volume_type==1;
%%%%%%%%;
if global_flag_volume_type==2;
x_sigma = 4*1/n_x_u;
a_x_u_ = zeros(n_x_u,n_x_u);
n_source = 128;
for nsource=0:n_source-1;
x_s_0 = (nsource/n_source).^(0.5)*0.75/sqrt(2)*sin(7.5*pi*nsource/n_source);
x_s_1 = (nsource/n_source).^(0.5)*0.75/sqrt(2)*cos(7.7*pi*nsource/n_source);
x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
a_x_u_ = a_x_u_ +  1/(sqrt(2*pi)*x_sigma).^2 * exp(-x_u_r__.^2/(2*x_sigma^2));
end;%for nsource=0:n_source-1;
end;%if global_flag_volume_type==2;
%%%%%%%%;
if global_flag_volume_type==3;
x_sigma = 1/n_x_u;
a_x_u_ = zeros(n_x_u,n_x_u);
n_source = 256;
for nsource=0:n_source-1;
x_s_0 = (nsource/n_source).^(0.5)*0.75/sqrt(2)*sin(4.7*pi*nsource/n_source);
x_s_1 = (nsource/n_source).^(0.5)*0.75/sqrt(2)*cos(4.3*pi*nsource/n_source);
x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_x_sigma = x_sigma * ( 8 + (n_source-nsource)/n_source * 8 );
a_x_u_ = a_x_u_ +  1/(sqrt(2*pi)*tmp_x_sigma).^2 * exp(-x_u_r__.^2/(2*tmp_x_sigma^2));
end;%for nsource=0:n_source-1;
end;%if global_flag_volume_type==3;
%%%%%%%%;
if global_flag_volume_type==4;
x_sigma = 1/n_x_u;
a_x_u_ = zeros(n_x_u,n_x_u);
n_source = 512;
for nsource=0:n_source-1;
x_s_0 = (0.15 + 0.85*nsource/n_source).^(0.5)*0.85/sqrt(2)*sin(40.7*pi*nsource/n_source);
x_s_1 = (0.15 + 0.85*nsource/n_source).^(0.5)*0.85/sqrt(2)*cos(31.5*pi*nsource/n_source);
x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_x_sigma = x_sigma * ( 6 - (n_source-nsource)/n_source * 4 );
a_x_u_ = a_x_u_ +  1/(sqrt(2*pi)*tmp_x_sigma).^2 * exp(-x_u_r__.^2/(2*tmp_x_sigma^2));
end;%for nsource=0:n_source-1;
x_s_0 = 0; x_s_1 = 0;
x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_x_sigma = 0.25;
a_x_u_ = a_x_u_ .*  (1/(sqrt(2*pi)*tmp_x_sigma).^2 * exp(-x_u_r__.^2/(2*tmp_x_sigma^2)));
a_x_u_ = reshape(atan(a_x_u_(:)/prctile(a_x_u_(:),99)),[n_x_u,n_x_u]);
end;%if global_flag_volume_type==4;
%%%%%%%%;
xx_u_weight_ = (2*x_p_r_max/n_x_u)^2;
dx = diameter_x_c / n_x_u;
a_x_u_l1 = sum(abs(a_x_u_).^1,'all')*dx^2;
a_x_u_l2 = sum(abs(a_x_u_).^2,'all')*dx^2;
if (verbose); disp(sprintf(' %% a_x_u_l1 = %0.16f',a_x_u_l1)); end;
if (verbose); disp(sprintf(' %% a_x_u_l2 = %0.16f',a_x_u_l2)); end;
%%%%%%%%;
% Calculate moments. ;
%%%%%%%%;
a_rho_x_u_ = a_x_u_ + min(a_x_u_,[],'all');
a_rho_x_c_0_avg = sum(x_u_0__.^1.*a_rho_x_u_/sum(a_rho_x_u_,'all'),'all');
a_rho_x_c_1_avg = sum(x_u_1__.^1.*a_rho_x_u_/sum(a_rho_x_u_,'all'),'all');
a_rho_x_c_0_std = sum((x_u_0__ - a_rho_x_c_0_avg).^2.*a_rho_x_u_/sum(a_rho_x_u_,'all'),'all');
a_rho_x_c_1_std = sum((x_u_1__ - a_rho_x_c_1_avg).^2.*a_rho_x_u_/sum(a_rho_x_u_,'all'),'all');
a_rho_x_c_avg_ = [a_rho_x_c_0_avg ; a_rho_x_c_1_avg ];
a_rho_x_c_std_ = [a_rho_x_c_0_std ; a_rho_x_c_1_std ];
disp(sprintf(' %% a_rho_x_c_std_ vs a_rho_x_c_avg_: %0.2f',fnorm(a_rho_x_c_std_)/fnorm(a_rho_x_c_avg_)));
if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u);
disp(sprintf(' %% Warning, molecule may not be well centered. Consider recentering.'));
end;%if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u);
%%%%%%%%;
% Possible to re-center. ;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u-1,-n_x_u/2,+n_x_u/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u-1,-n_x_u/2,+n_x_u/2)/2; %<-- box has diameter 2. ;
[k_u_0__,k_u_1__] = ndgrid(k_u_0_,k_u_1_); n_kk_u = n_x_u^2;
b_rho_x_u_ = real(ifftn(fftn(a_rho_x_u_).*exp(+i*2*pi*(k_u_0__*a_rho_x_c_0_avg + k_u_1__*a_rho_x_c_1_avg))));
b_rho_x_c_0_avg = sum(x_u_0__.^1.*b_rho_x_u_/sum(b_rho_x_u_,'all'),'all');
b_rho_x_c_1_avg = sum(x_u_1__.^1.*b_rho_x_u_/sum(b_rho_x_u_,'all'),'all');
b_rho_x_c_0_std = sum((x_u_0__ - b_rho_x_c_0_avg).^2.*b_rho_x_u_/sum(b_rho_x_u_,'all'),'all');
b_rho_x_c_1_std = sum((x_u_1__ - b_rho_x_c_1_avg).^2.*b_rho_x_u_/sum(b_rho_x_u_,'all'),'all');
b_rho_x_c_avg_ = [b_rho_x_c_0_avg ; b_rho_x_c_1_avg];
b_rho_x_c_std_ = [b_rho_x_c_0_std ; b_rho_x_c_1_std];
disp(sprintf(' %% b_rho_x_c_std_ vs b_rho_x_c_avg_: %0.2f',fnorm(b_rho_x_c_std_)/fnorm(b_rho_x_c_avg_)));
flag_plot=0;
if flag_plot;
figure(1);clf;figmed;figbeach();
alim_ = [0,prctile(a_x_u_,99,'all')];
subplot(1,2,1);imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,a_rho_x_u_,alim_); axis image; axisnotick;
subplot(1,2,2);imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,b_rho_x_u_,alim_); axis image; axisnotick;
if flag_plot; disp('returning'); return; end;
end;%if flag_plot;
%%%%%%%%;
a_x_u_base_ = a_x_u_;
if flag_center_volume;
disp(sprintf(' %% centering volume'));
a_x_u_base_ = b_rho_x_u_;
end;%if flag_center_volume;
%%%%%%%%;
if flag_plot;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_base_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,2,1); imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,a_x_u_,alim_); title('a'); axis image; axisnotick;
subplot(1,2,2); imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,b_rho_x_u_,alim_); title('b'); axis image; axisnotick;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;
%%%%%%%%;
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','n_x_u' ...
     ,'x_u_0_','x_u_1_','x_u_0__','x_u_1__','n_x_u','n_xx_u','xx_u_weight_' ...
     ,'a_x_u_base_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

k_p_r_max = k_p_r_max_factor_0in*48/(2*pi); k_eq_d = 1.0/(2*pi);
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
,J_node_ ...
,J_weight_ ...
,J_chebfun_ ...
,J_polyval_ ...
] = ...
sample_sphere_7( ...
 verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;

n_w_max = 2*(n_k_p_r);
template_k_eq_d = -1; n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%template_k_eq_d = 1/(2*pi); n_w_0in_ = [];
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_weight_2d_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);

dx = diameter_x_c / n_x_u;
a_k_p_quad_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,a_x_u_base_ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
) ...
* sqrt(n_x_u^2)*dx^2;
a_x_u_reco_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,a_k_p_quad_.*weight_2d_k_all_*(2*pi)^2 ...
) ...
* sqrt(n_x_u^2) * n_w_sum;
a_k_q_quad_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,a_k_p_quad_);
a_k_q_flip_ = image_q_flip_0(n_k_p_r,n_w_,a_k_q_quad_);
a_k_p_flip_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,a_k_q_flip_);
a_x_u_flip_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,a_k_p_flip_.*weight_2d_k_all_*(2*pi)^2 ...
) ...
* sqrt(n_x_u^2) * n_w_sum;

a_x_u_base_l2 = sum(abs(a_x_u_base_(:)).^2.*dx^2);
a_k_p_quad_l2 = sum(abs(a_k_p_quad_(:)).^2.*weight_2d_k_all_) * (2*pi)^2;
a_k_q_quad_l2 = sum(abs(a_k_q_quad_(:)).^2.*weight_2d_k_all_) * (2*pi)^2;
a_x_u_reco_l2 = sum(abs(a_x_u_reco_).^2,'all').*dx^2;
a_x_u_flip_l2 = sum(abs(a_x_u_flip_).^2,'all').*dx^2;
disp(sprintf(' %% a_x_u_base_l2 = %0.16f',a_x_u_base_l2));
disp(sprintf(' %% a_k_p_quad_l2 = %0.16f',a_k_p_quad_l2));
disp(sprintf(' %% a_k_q_quad_l2 = %0.16f',a_k_q_quad_l2));
disp(sprintf(' %% a_x_u_reco_l2 = %0.16f',a_x_u_reco_l2));
disp(sprintf(' %% a_x_u_flip_l2 = %0.16f',a_x_u_flip_l2));

if flag_plot;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_reco_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,3,1);imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,a_x_u_base_,[0,1]); axis image; axisnotick;
title('a_x_u_base_','Interpreter','none');
subplot(1,3,2);imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,real(a_x_u_reco_),[0,1]); axis image; axisnotick;
title('a_x_u_reco_','Interpreter','none');
subplot(1,3,3);imagesc_c_mask(n_x_u,x_u_0,n_x_u,x_u_1_,real(a_x_u_flip_),[0,1]); axis image; axisnotick;
title('a_x_u_flip_','Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

if flag_plot;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_p_quad_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,2,1); 
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(a_k_p_quad_)); axis image; axisnotick; title('real');
subplot(1,2,2); 
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(a_k_p_quad_)); axis image; axisnotick; title('imag');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

flag_check=0;
if flag_check;
%%%%%%%%;
% test rotations. ;
%%%%%%%%;
gamma_z_ = 2*pi*[0:n_w_max-1]/n_w_max; 
%%%%;
X_0_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = gamma_z_(1+nw);
b_k_p_quad_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,a_k_p_quad_,gamma_z);
X_0 = sum(conj(b_k_p_quad_).*a_k_p_flip_.*weight_2d_k_all_) * (2*pi)^2;
X_0_w_(1+nw) = X_0;
end;%for nw=0:n_w_max-1;
%%%%;
X_1_q_ = innerproduct_q_k_stretch_0(n_k_p_r,k_p_r_,n_w_,n_w_sum,a_k_q_quad_,a_k_q_flip_);
X_1_w_ = ifft(X_1_q_)*n_w_max;
%%%%;
X_2_q_ = sum(conj(reshape(a_k_q_quad_,[n_w_max,n_k_p_r])).*reshape(a_k_q_flip_,[n_w_max,n_k_p_r]).*reshape(weight_2d_k_all_,[n_w_max,n_k_p_r]),2) * (2*pi)^2;
X_2_w_ = ifft(X_2_q_)*n_w_max;
%%%%;
figure(1);clf;figmed;
subplot(1,3,1);
plot(gamma_z_,real(X_0_w_),'r.-',gamma_z_,imag(X_0_w_),'b.-');
xlabel('gamma');ylabel('X'); title('X_0','Interpreter','none');
subplot(1,3,2);
plot(gamma_z_,real(X_1_w_),'r.-',gamma_z_,imag(X_1_w_),'b.-');
xlabel('gamma');ylabel('X'); title('X_1','Interpreter','none');
subplot(1,3,3);
plot(gamma_z_,real(X_2_w_),'r.-',gamma_z_,imag(X_2_w_),'b.-');
xlabel('gamma');ylabel('X'); title('X_2','Interpreter','none');
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% ignore CTF for now. ;
%%%%%%%%;
CTF_k_p_r_ = ones(n_k_p_r,1);

tolerance_pm = tolerance_master;
%%%%%%%%;
% Find principal-modes. ;
%%%%%%%%;
[ ...
 X_kk__ ...
,X_weight_r_k_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,1 ...
,a_k_p_quad_ ...
);
n_UX_rank = n_k_p_r - 1; %<-- just to check dimensions. ;
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank_0in = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_kn__ = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_k_ = tmp_SX_(1+[0:n_UX_rank-1]);

%%%%%%%%;
% set up principal-volume. ;
%%%%%%%%;
UX_a_k_p_quad_wn__ = zeros(n_w_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
UX_a_k_p_quad_wn__(:,1+nUX_rank) = reshape(a_k_p_quad_,[n_w_max,n_k_p_r])*(UX_kn__(:,1+nUX_rank).*X_weight_r_k_);
end;%for nUX_rank=0:n_UX_rank-1;
UX_a_k_p_quad_l2 = sum(abs(UX_a_k_p_quad_wn__(:)).^2) / n_w_max;
disp(sprintf(' %% a_k_p_quad_l2 = %0.16f',a_k_p_quad_l2));
disp(sprintf(' %% UX_a_k_p_quad_l2 = %0.16f',UX_a_k_p_quad_l2));

%%%%%%%%;
% set up scattershot-images. ;
%%%%%%%%;
sigma_noise_per_Pixel = 0;
if (snr_per_A>0); sigma_noise_per_Pixel = a_k_p_quad_l2/max(1e-12,snr_per_A)/sqrt(n_x_u)/sqrt(Pixel_A); end;
n_M = 1024; %<-- number of images. ;
n_M = n_w_max*n_M_factor_0in; %<-- number of images. ;
n_pick = ceil(sqrt(n_w_max)); %<-- controls overlap. ;
n_pick = n_pick_0in; %<--controls overlap. ;
M_k_p_pkM___ = zeros(n_pick,n_k_p_r,n_M);
UX_M_k_p_pnM__ = zeros(n_pick,n_UX_rank,n_M);
euler_gamma_z_true_pM__ = zeros(n_pick,n_M);
index_nw_from_nM_true_pM__ = zeros(n_pick,n_M);
na=0;
for nM=0:n_M-1;
if (verbose); if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
tmp_p_ = randperm(n_w_max);
nw_ = tmp_p_(1+[0:n_pick-1]);
nw_ = periodize(nw_-nw_(1+0),0,n_w_max);
nw_ = periodize(nw_+mod(nM,n_w_max),0,n_w_max);
index_nw_from_nM_true_pM__(:,1+nM) = nw_;
euler_gamma_z_true_pM__(:,1+nM) = 2*pi*nw_/n_w_max;
%%%%;
b_x_u_temp_ = a_x_u_base_ + sigma_noise_per_Pixel*rand(n_x_u,n_x_u);
b_k_p_temp_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,b_x_u_temp_ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
) ...
* sqrt(n_x_u^2)*dx^2;
b_k_p_temp_wk__ = reshape(b_k_p_temp_,[n_w_max,n_k_p_r]);
UX_b_k_p_temp_wn__ = zeros(n_w_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
UX_b_k_p_temp_wn__(:,1+nUX_rank) = reshape(b_k_p_temp_,[n_w_max,n_k_p_r])*(UX_kn__(:,1+nUX_rank).*X_weight_r_k_);
end;%for nUX_rank=0:n_UX_rank-1;
%%%%;
M_k_p_pkM___(:,:,1+nM) = b_k_p_temp_wk__(1+nw_,:);
UX_M_k_p_pnM___(:,:,1+nM) = UX_b_k_p_temp_wn__(1+nw_,:);
clear b_k_p_temp_wk__ UX_b_k_p_temp_wn__;
end;%for nM=0:n_M-1;

euler_gamma_z_true_M_ = reshape(euler_gamma_z_true_pM__(1+0,:),[n_M,1]);
index_nw_from_nM_true_M_ = reshape(index_nw_from_nM_true_pM__(1+0,:),[n_M,1]);
dgamma_z_pM__ = bsxfun(@minus,euler_gamma_z_true_pM__,euler_gamma_z_true_pM__(1+0,:));
dgamma_z_pM__ = periodize(dgamma_z_pM__,0,2*pi);
dnw_pM__ = bsxfun(@minus,index_nw_from_nM_true_pM__,index_nw_from_nM_true_pM__(1+0,:));

if flag_plot;
%%%%%%%%;
% plot principal-images. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/UX_M_k_p_pnM__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;colormap('hsv');
gamma_z_ = 2*pi*transpose([0:n_w_max-1])/n_w_max;
c_hsv__ = colormap('hsv'); n_c_hsv = size(c_hsv__,1);
markersize_use = 8;
zlim_ = 0 + 1.5*std(abs(UX_a_k_p_quad_wn__(:,1)),1,'all')*[-1,+1];
p_top = min(6,pm_n_UX_rank_0in); p_row = 2; p_col = ceil(p_top/p_row);
for nUX_rank=0:p_top-1;
subplot(p_row,p_col,1+nUX_rank);
hold on;
tmp_z_ = UX_a_k_p_quad_wn__(:,1+nUX_rank);
surfline_0(real(tmp_z_),imag(tmp_z_),gamma_z_);
tmp_z__ = UX_M_k_p_pnM___(:,1+nUX_rank,:);
tmp_c__ = euler_gamma_z_true_pM__(:,:);
scatter(real(tmp_z__(:)),imag(tmp_z__(:)),markersize_use,tmp_c__(:),'filled');
hold off;
xlim(zlim_); ylim(zlim_); axis square;
xlabel('real'); ylabel('imag');
title(sprintf('nUX_rank %d',nUX_rank),'Interpreter','none');
end;%for nUX_rank=0:p_top-1;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

%%%%%%%%;
% Now set up simple model linking the observations to a model a_k_q_quad_. ;
%%%%%%%%;
gamma_z_M_ = zeros(n_M,1);
gamma_z_M_ = euler_gamma_z_true_pM__(1+0,:);
[ M_from_a_pMq___ ] = ...
M_from_a_SR0A_0( ...
 n_w_max ...
,gamma_z_M_ ...
,n_pick ...
,dgamma_z_pM__ ...
);
M_from_a_pMq__ = reshape(M_from_a_pMq___,[n_pick*n_M,n_w_max]);
a_from_M_qpM__ = pinv(M_from_a_pMq__,1e-6);
%%%%%%%%;
% build best possible full reconstruction. ;
%%%%%%%%;
a_k_q_reco_wk__ = zeros(n_w_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_M_pM_ = reshape(M_k_p_pkM___(:,1+nk_p_r,:),[n_pick*n_M,1]);
a_k_q_reco_wk__(:,1+nk_p_r) = a_from_M_qpM__ * tmp_M_pM_;
end;%for nk_p_r=0:n_k_p_r-1;
a_k_p_reco_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,a_k_q_reco_wk__(:));
a_k_p_reco_l2 = sum(abs(a_k_p_reco_wk_(:)).^2.*weight_2d_k_all_) * (2*pi)^2;
%%%%%%%%;
a_x_u_reco_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,a_k_p_reco_wk_.*weight_2d_k_all_*(2*pi)^2 ...
) ...
* sqrt(n_x_u^2) * n_w_sum;
%%%%%%%%;
% test best possible principal-mode reconstruction. ;
%%%%%%%%;
UX_a_k_q_quad_wn_ = interp_p_to_q(n_UX_rank,n_w_,n_w_max*n_UX_rank,UX_a_k_p_quad_wn__);
UX_a_k_q_quad_wn__ = reshape(UX_a_k_q_quad_wn_,[n_w_max,n_UX_rank]);
UX_a_k_q_reco_wn__ = zeros(n_w_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_M_pM_ = reshape(UX_M_k_p_pnM___(:,1+nUX_rank,:),[n_pick*n_M,1]);
UX_a_k_q_reco_wn__(:,1+nUX_rank) = a_from_M_qpM__ * tmp_UX_M_pM_;
end;%for nUX_rank=0:n_UX_rank-1;
l2_best_n_ = zeros(n_UX_rank,1);
X_best_n_ = zeros(n_UX_rank,1);
for nUX_rank=0:n_UX_rank-1;
tmp_0_ = UX_a_k_q_quad_wn__(:,1+[0:nUX_rank]); tmp_0_ = tmp_0_(:);
tmp_1_ = UX_a_k_q_reco_wn__(:,1+[0:nUX_rank]); tmp_1_ = tmp_1_(:);
l2_best_n_(1+nUX_rank) = fnorm(tmp_0_-tmp_1_)/fnorm(tmp_0_);
X_best_n_(1+nUX_rank) = real(corr(tmp_0_,tmp_1_));
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
% plot best recovery as a function of pm_n_UX_rank. ;
%%%%%%%%;
if flag_plot;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_best_n_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(1:n_UX_rank,l2_best_n_,'k.-','LineWidth',2,'MarkerSize',16);
xlabel('n_UX_rank','Interpreter','none');
ylabel('l2'); xlim([0,n_UX_rank]); ylim([-0.25,1.5]);
grid on;
subplot(1,2,2);
plot(1:n_UX_rank,X_best_n_,'k.-','LineWidth',2,'MarkerSize',16);
xlabel('n_UX_rank','Interpreter','none');
ylabel('X'); xlim([0,n_UX_rank]); ylim([-0.25,1.5]);
grid on;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now set up simple alignment operation. ;
%%%%%%%%;
euler_gamma_z_Ma_l2_n_ = zeros(n_UX_rank,1);
euler_gamma_z_aM_l2_n_ = zeros(n_UX_rank,1);
for nUX_rank=0:n_UX_rank-1;
parameter = struct('type','parameter');
parameter.f_rand = f_rand_0in;
[ ...
 parameter ...
,gamma_z_aM_M_ ...
,M_vs_a_l2_aM_M_ ...
,gamma_z_Ma_M_ ...
,M_vs_a_l2_Ma_M_ ...
,M_vs_a_l2_wM__ ...
] = ...
gamma_z_from_a_M_SR0A_0( ...
 parameter ...
,n_w_max ...
,1+nUX_rank ...
,UX_a_k_p_quad_wn__(:,1+[0:nUX_rank]) ...
,n_pick ...
,n_M ...
,UX_M_k_p_pnM___(:,1+[0:nUX_rank],:) ...
,dnw_pM__ ...
);
euler_gamma_z_Ma_l2_n_(1+nUX_rank) = sum(abs(periodize(euler_gamma_z_true_M_-gamma_z_Ma_M_,-pi,+pi)).^2);
euler_gamma_z_aM_l2_n_(1+nUX_rank) = sum(abs(periodize(euler_gamma_z_true_M_-gamma_z_aM_M_,-pi,+pi)).^2);
if (verbose); disp(sprintf(' %% nUX_rank %d/%d: Ma_l2 %0.3f aM_l2 %0.3f',nUX_rank,n_UX_rank,euler_gamma_z_Ma_l2_n_(1+nUX_rank),euler_gamma_z_aM_l2_n_(1+nUX_rank))); end;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
% plot alignment as a function of pm_n_UX_rank. ;
%%%%%%%%;
if flag_plot;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/euler_gamma_z_l2_n_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);
plot(1:n_UX_rank,euler_gamma_z_Ma_l2_n_,'k.-','LineWidth',2,'MarkerSize',16);
xlabel('n_UX_rank','Interpreter','none');
ylabel('l2'); xlim([0,n_UX_rank]);
title('euler_gamma_z_Ma_l2_n_','Interpreter','none');
grid on;
subplot(1,2,2);
plot(1:n_UX_rank,euler_gamma_z_Ma_l2_n_,'k.-','LineWidth',2,'MarkerSize',16);
xlabel('n_UX_rank','Interpreter','none');
ylabel('X'); xlim([0,n_UX_rank]);
title('euler_gamma_z_aM_l2_n_','Interpreter','none');
grid on;
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if flag_plot;
%%%%%%%%;
end;%if flag_check;

%%%%%%%%;
% Now try alternating minimization. ;
%%%%%%%%;
UX_a_k_q_quad_wn_ = interp_p_to_q(n_UX_rank,n_w_,n_w_max*n_UX_rank,UX_a_k_p_quad_wn__);
UX_a_k_q_quad_wn__ = reshape(UX_a_k_q_quad_wn_,[n_w_max,n_UX_rank]); %<-- ground truth. ;
%%%%;

pm_n_UX_rank_max = min(n_UX_rank,2*pm_n_UX_rank_0in);
X_inf___ = zeros(global_n_iteration,pm_n_UX_rank_max,2);
kld_inf___ = zeros(global_n_iteration,pm_n_UX_rank_max,2);
flag_vis = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for local_flag_Ma_vs_aM=[0:1];
%for pm_n_UX_rank=[2,pm_n_UX_rank_0in,2*pm_n_UX_rank_0in];%for pm_n_UX_rank=1:pm_n_UX_rank_max;
for pm_n_UX_rank=pm_n_UX_rank_0in;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

string_infix_use = string_infix;
if (local_flag_Ma_vs_aM);
string_infix_use = sprintf('%sa1',string_infix);
end;%if (local_flag_Ma_vs_aM);

if flag_vis;
figure(1+nf);nf=nf+1;clf;figbig; figbeach();
d_r = 2*pi/(n_M/2);
p_row = 4; p_col = ceil((5+global_n_iteration)/p_row); np=0;
subplot(p_row,p_col,1+np);np=np+1;
x_s_0 = 0; x_s_1 = 0; x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_index_ = efind(x_u_r__(:)<half_diameter_x_c/sqrt(2));
aavg = mean(real(a_x_u_base_(1+tmp_index_))); astd = std(real(a_x_u_base_(1+tmp_index_)),1);
alim_ = aavg + 2.5*astd*[-1,+1];
alim_ = prctile(real(a_x_u_base_(1+tmp_index_)),[ 1,99]);
imagesc_c_mask(n_x_u,x_u_0_,n_x_u,x_u_1_,real(a_x_u_base_),alim_);
axis image; axisnotick; title('orig');
subplot(p_row,p_col,1+np);np=np+1;
x_s_0 = 0; x_s_1 = 0; x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_index_ = efind(x_u_r__(:)<half_diameter_x_c/sqrt(2));
aavg = mean(real(a_x_u_reco_xx__(1+tmp_index_))); astd = std(real(a_x_u_reco_xx__(1+tmp_index_)),1);
alim_ = aavg + 2.5*astd*[-1,+1];
alim_ = prctile(real(a_x_u_reco_xx__(1+tmp_index_)),[ 1,99]);
hold on;
imagesc_c_mask(n_x_u,x_u_0_,n_x_u,x_u_1_,real(a_x_u_reco_xx__),alim_);
axis image; axisnotick; title('reco');
image_SR0A_gamma_z_0( ...
 n_M ...
,n_w_max ...
,index_nw_from_nM_true_M_ ...
,index_nw_from_nM_true_M_ ...
,0 ...
,0 ...
,1 ...
,d_r ...
);
hold off;
drawnow();
end;%if flag_vis;
if (verbose>0); disp(sprintf(' %% pm_n_UX_rank %d/%d',pm_n_UX_rank,pm_n_UX_rank_max)); end;
UX_a_k_p_quad_sub_l2 = sum(abs(UX_a_k_p_quad_wn__(:,1+[0:pm_n_UX_rank-1])).^2,'all') / n_w_max;
rng(local_rseed);
nw_M_ = max(0,min(n_w_max-1,floor(n_w_max*rand(n_M,1))));
gamma_z_M_ = 2*pi*nw_M_/n_w_max;
parameter = struct('type','parameter');
parameter.f_rand = f_rand_0in;
parameter.flag_Ma_vs_aM = local_flag_Ma_vs_aM;
%%%%%%%%;
if (~isfield(parameter,'n_iteration')); parameter.n_iteration = global_n_iteration; end;
if (~isfield(parameter,'flag_Ma_vs_aM')); parameter.flag_Ma_vs_aM = 1; end;
n_iteration = parameter.n_iteration;
flag_Ma_vs_aM = parameter.flag_Ma_vs_aM;
X_i_ = zeros(n_iteration,1);
Y_i_ = zeros(n_iteration,1);
kld_i_ = zeros(n_iteration,1);
UX_b_k_q_reco_wni___ = zeros(n_w_max,pm_n_UX_rank,n_iteration);
UX_b_k_p_reco_wni___ = zeros(n_w_max,pm_n_UX_rank,n_iteration);
nw_Mi__ = zeros(n_M,n_iteration);
gamma_z_Mi__ = zeros(n_M,n_iteration);
%%%%%%%%;
for niteration=0:n_iteration-1;
%%%%%%%%;
if (verbose>1); disp(sprintf(' %% niteration %d/%d',niteration,n_iteration)); end;
[ M_from_a_pMq___ ] = ...
M_from_a_SR0A_0( ...
 n_w_max ...
,gamma_z_M_ ...
,n_pick ...
,dgamma_z_pM__ ...
);
M_from_a_pMq__ = reshape(M_from_a_pMq___,[n_pick*n_M,n_w_max]);
a_from_M_qpM__ = pinv(M_from_a_pMq__,1e-6);
UX_b_k_q_reco_wn__ = zeros(n_w_max,pm_n_UX_rank);
for nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_M_pM_ = reshape(UX_M_k_p_pnM___(:,1+nUX_rank,:),[n_pick*n_M,1]);
UX_b_k_q_reco_wn__(:,1+nUX_rank) = a_from_M_qpM__ * tmp_UX_M_pM_;
end;%for nUX_rank=0:pm_n_UX_rank-1;
%%%%;
b_k_q_reco_wk__ = zeros(n_w_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_M_pM_ = reshape(M_k_p_pkM___(:,1+nk_p_r,:),[n_pick*n_M,1]);
b_k_q_reco_wk__(:,1+nk_p_r) = a_from_M_qpM__ * tmp_M_pM_;
end;%for nk_p_r=0:n_k_p_r-1;
b_k_p_reco_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,b_k_q_reco_wk__(:));
b_k_p_reco_l2 = sum(abs(b_k_p_reco_wk_(:)).^2.*weight_2d_k_all_) * (2*pi)^2;
%%%%%%%%;
if flag_vis;
b_x_u_reco_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_u ...
,diameter_x_c ...
,n_x_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,b_k_p_reco_wk_.*weight_2d_k_all_*(2*pi)^2 ...
) ...
* sqrt(n_x_u^2) * n_w_sum;
subplot(p_row,p_col,1+np);np=np+1;
x_s_0 = 0; x_s_1 = 0; x_u_r__ = sqrt((x_u_0__ - x_s_0).^2 + (x_u_1__-x_s_1).^2);
tmp_index_ = efind(x_u_r__(:)<half_diameter_x_c/sqrt(2));
aavg = mean(real(b_x_u_reco_xx__(1+tmp_index_))); astd = std(real(b_x_u_reco_xx__(1+tmp_index_)),1);
alim_ = aavg + 2.5*astd*[-1,+1];
alim_ = prctile(real(b_x_u_reco_xx__(1+tmp_index_)),[ 1,99]);
hold on;
imagesc_c_mask(n_x_u,x_u_0_,n_x_u,x_u_1_,real(b_x_u_reco_xx__),alim_);
axis image; axisnotick; title(sprintf('ni %d',niteration));
image_SR0A_gamma_z_0( ...
 n_M ...
,n_w_max ...
,index_nw_from_nM_true_M_ ...
,nw_M_ ...
,0 ...
,0 ...
,1 ...
,d_r ...
);
hold off;
drawnow();
end;%if flag_vis;
%%%%;
UX_b_k_p_reco_wn__ = reshape(interp_q_to_p(pm_n_UX_rank,n_w_max*ones(pm_n_UX_rank,1),n_w_max*pm_n_UX_rank,UX_b_k_q_reco_wn__(:)),[n_w_max,pm_n_UX_rank]);
%%%%;
[ ...
 parameter ...
,nw_aM_M_ ...
,gamma_z_aM_M_ ...
,M_vs_a_l2_aM_M_ ...
,nw_Ma_M_ ...
,gamma_z_Ma_M_ ...
,M_vs_a_l2_Ma_M_ ...
,M_vs_a_l2_wM__ ...
] = ...
gamma_z_from_a_M_SR0A_0( ...
 parameter ...
,n_w_max ...
,pm_n_UX_rank ...
,UX_b_k_p_reco_wn__(:,1+[0:pm_n_UX_rank-1]) ...
,n_pick ...
,n_M ...
,UX_M_k_p_pnM___(:,1+[0:pm_n_UX_rank-1],:) ...
,dnw_pM__ ...
);
%%%%;
if (flag_Ma_vs_aM==0); nw_M_ = nw_aM_M_; gamma_z_M_ = gamma_z_aM_M_; end;
if (flag_Ma_vs_aM==1); nw_M_ = nw_Ma_M_; gamma_z_M_ = gamma_z_Ma_M_; end;
%%%%;
UX_b_k_p_reco_sub_l2 = sum(abs(UX_b_k_p_reco_wn__(:,1+[0:pm_n_UX_rank-1])).^2,'all') / n_w_max;
tmp_X_0_q_ = sum(conj(UX_a_k_q_quad_wn__(:,1+[0:pm_n_UX_rank-1])).*UX_b_k_q_reco_wn__,2)/n_w_max;
tmp_X_0_w_ = ifft(tmp_X_0_q_)*n_w_max;
tmp_X_1_q_ = sum(conj(UX_a_k_q_quad_wn__(:,1+[0:pm_n_UX_rank-1])).*reshape(image_q_flip_0(pm_n_UX_rank,n_w_max,UX_b_k_q_reco_wn__),[n_w_max,pm_n_UX_rank]),2)/n_w_max;
tmp_X_1_w_ = ifft(tmp_X_1_q_)*n_w_max;
tmp_X = max(max(real(tmp_X_0_w_)),max(real(tmp_X_1_w_)));
%%%%;
tmp_Y_0_q_ = sum(conj(reshape(a_k_q_quad_,[n_w_max,n_k_p_r])).*reshape(b_k_q_reco_wk__(:),[n_w_max,n_k_p_r]).*reshape(weight_2d_k_all_,[n_w_max,n_k_p_r]),2) * (2*pi)^2;
tmp_Y_0_w_ = ifft(tmp_Y_0_q_)*n_w_max;
tmp_Y_1_q_ = sum(conj(reshape(a_k_q_flip_,[n_w_max,n_k_p_r])).*reshape(b_k_q_reco_wk__(:),[n_w_max,n_k_p_r]).*reshape(weight_2d_k_all_,[n_w_max,n_k_p_r]),2) * (2*pi)^2;
tmp_Y_1_w_ = ifft(tmp_Y_1_q_)*n_w_max;
tmp_Y = max(max(real(tmp_Y_0_w_)),max(real(tmp_Y_1_w_)));
%%%%;
p_nw_w_ = sparse(1+nw_M_,1,1,n_w_max,1)/n_M;
plp_nw_w_ = p_nw_w_.*log(p_nw_w_);
tmp_index_ = efind(isfinite(plp_nw_w_));
kld = sum(plp_nw_w_(1+tmp_index_)) + log(n_w_max);
%%%%;
X_i_(1+niteration) = tmp_X / sqrt(UX_a_k_p_quad_sub_l2*UX_b_k_p_reco_sub_l2);
Y_i_(1+niteration) = tmp_Y / sqrt(a_k_p_quad_l2*b_k_p_reco_l2);
UX_b_k_q_reco_wni___(:,:,1+niteration) = UX_b_k_q_reco_wn__;
UX_b_k_p_reco_wni___(:,:,1+niteration) = UX_b_k_p_reco_wn__;
gamma_z_Mi__(:,1+niteration) = gamma_z_M_;
nw_Mi__(:,1+niteration) = nw_M_;
kld_i_(1+niteration) = kld;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
if flag_vis;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(1:n_iteration,X_i_,'k.-','LineWidth',2,'MarkerSize',16);
plot(1:n_iteration,Y_i_,'r.-','LineWidth',2,'MarkerSize',16);
hold off;
legend({'X','Y'},'Location','SouthWest');
xlim([1,n_iteration]);xlabel('niteration');
ylim([-0.1,+1.1]);ylabel('X,Y'); grid on;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot(1:n_iteration,kld_i_,'b.-','LineWidth',2,'MarkerSize',16);
xlim([1,n_iteration]);xlabel('niteration');
ylim([-0.1,+1.1]);ylabel('kld'); grid on;
%%%%;
h__ = zeros(n_w_max,n_iteration);
for niteration=0:n_iteration-1;
h__(:,1+niteration) = full(sparse(1+nw_Mi__(:,1+niteration),1,1,n_w_max,1));
end;%for niteration=0:n_iteration-1;
subplot(p_row,p_col,1+np);np=np+1;
imagesc(h__,[0,2*n_M/n_w_max]); xlabel('iteration'); ylabel('nw','Interpreter','none');
%%%%;
sgtitle(sprintf('%s: pm_n_UX_rank %d',string_infix_use,pm_n_UX_rank),'Interpreter','none');
%%%%;
end;%if flag_vis;

X_inf___(:,pm_n_UX_rank,1+local_flag_Ma_vs_aM) = X_i_;
kld_inf___(:,pm_n_UX_rank,1+local_flag_Ma_vs_aM) = kld_i_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for pm_n_UX_rank=1:pm_n_UX_rank_max;
end;%for local_flag_Ma_vs_aM=[0:1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



