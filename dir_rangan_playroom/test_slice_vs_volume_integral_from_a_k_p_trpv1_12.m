%%%%%%%%;
% Applying eig_ddssnll_lanczos_from_a_k_p_4 to trpv1. ;
% Limiting to principal modes. ;
%%%%%%%%;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

k_int = 8;

str_thisfunction = 'test_slice_vs_volume_integral_from_a_k_p_trpv1_12';
flag_recalc=0; flag_replot=0;
flag_verbose=1; flag_disp=1; nf=0;

[~,str_hostname] = system('hostname');
flag_256G = 0 ...
| ~isempty(strfind(str_hostname,'crunchy')) ...
| ~isempty(strfind(str_hostname,'linserv')) ...
;
flag_128G = 0 ...
| ~isempty(strfind(str_hostname,'crunchy')) ...
| ~isempty(strfind(str_hostname,'linserv')) ...
| ~isempty(strfind(str_hostname,'xcalibr8')) ...
| ~isempty(strfind(str_hostname,'snappy')) ...
;
%%;
if ~exist('k_int','var'); k_int = 32; if flag_256G; k_int = 48; end; end;
%%;
if k_int==8;
k_eq_d_double = 0.50;
t_eq_d_double = 0.50;
n_w_int = 2;
pm_n_UX_rank_max = 3;
n_order = 5;
lanczos_n_iteration_max = pm_n_UX_rank_max*(k_int+1)^2; %<-- This should be close to a fully resolved calculation with 3 principal-modes. ;
end;%if k_int==8;
%%;
if k_int==16;
k_eq_d_double = 0.50;
t_eq_d_double = 0.50;
n_w_int = 2;
n_order = 5;
pm_n_UX_rank_max = +Inf;
lanczos_n_iteration_max = 256;
end;%if k_int==16;
%%;
if k_int==32;
k_eq_d_double = 0.50;
t_eq_d_double = 0.50;
n_w_int = 2;
n_order = 9;
pm_n_UX_rank_max = +Inf;
lanczos_n_iteration_max = 64;
end;%if k_int==32;
%%;
if k_int==48;
k_eq_d_double = 1.00;
t_eq_d_double = 1.00;
n_w_int = 1;
n_order = 9;
pm_n_UX_rank_max = +Inf;
lanczos_n_iteration_max = 32;
end;%if k_int==48;
%%;
KAPPA_flag_kernel_full = 1;
KAPPA_basic_l_max_use = k_int;
KAPPA_basic_l_max_ext = ceil(1.25*KAPPA_basic_l_max_use);
KAPPA_basic_l_max_band = floor(KAPPA_basic_l_max_use/2);
KAPPA_qref_k_eq_d_double = k_eq_d_double;
%%;
flag_calc = flag_256G | (k_int<=32 & flag_128G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
global_parameter=struct('type','parameter');
global_parameter.flag_recalc = flag_recalc;
global_parameter.flag_replot = flag_replot;
fname_prefix='trpv1_x0';
dir_nopath_data_star='trpv1';
Pixel_Spacing=1.2156;
fname_nopath_volume='emd_5778.mrc';
fname_nopath_star='tv1_relion_data.star';
%%%%%%%%;
if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_volume')); global_parameter.flag_center_volume = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_center_image')); global_parameter.flag_center_image = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_invert')); global_parameter.flag_invert = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
flag_center_volume = global_parameter.flag_center_volume;
flag_center_image = global_parameter.flag_center_image;
flag_invert = global_parameter.flag_invert;
tolerance_master = global_parameter.tolerance_master;
nf=0;
%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_ssnll_ori = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_ssnll_k%d',string_root,fname_prefix_xfix,k_int);
dir_ssnll = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_ssnll_from_a_k_p_k%d',string_root,fname_prefix_xfix,k_int);
if (~exist(sprintf('%s_mat',dir_ssnll),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_ssnll)); mkdir(sprintf('%s_mat',dir_ssnll)); end;
if (~exist(sprintf('%s_jpg',dir_ssnll),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_ssnll)); mkdir(sprintf('%s_jpg',dir_ssnll)); end;
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
%%%%%%%%;

%%%%%%%%;
% We can examine the volume avg and std as follows: ;
%%%%%%%%;
fname_emd = sprintf('%s/%s',dir_data_star,fname_nopath_volume);
a_x_u_load_ = cast(ReadMRC(fname_emd),'double');
n_x_u = size(a_x_u_load_,1);
if (flag_verbose>0); disp(sprintf(' %% sqrt(sum(max(0,a_x_u_load_.^2))): %+0.6f',sqrt(sum(max(0,a_x_u_load_.^2),'all')))); end;
if (flag_verbose>0); disp(sprintf(' %% sqrt(sum(max(0,a_x_u_load_).^2)): %+0.6f',sqrt(sum(max(0,a_x_u_load_).^2,'all')))); end;
if (flag_verbose>0); disp(sprintf(' %% sum(max(0,a_x_u_load_)): %+0.6f',sum(max(0,a_x_u_load_),'all'))); end;
if (flag_verbose>0); disp(sprintf(' %% sum(a_x_u_load_): %+0.6f',sum(a_x_u_load_,'all'))); end;
clear a_x_u_load_;
%%%%%%%%;

%%%%%%%%;
% First create consensus volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_ssnll);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_emd = sprintf('%s/%s',dir_data_star,fname_nopath_volume);
a_x_u_load_ = cast(ReadMRC(fname_emd),'double');
n_x_u = size(a_x_u_load_,1);
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64;
n_pack = n_x_u/n_x_u_pack;
pack_row_ij_ = zeros(n_x_u_pack,1);
pack_col_ij_ = zeros(n_x_u_pack,1);
pack_val_ij_ = zeros(n_x_u_pack,1);
na=0;
for nx_u=0:n_x_u-1;
pack_row_ij_(1+na) = 1+nx_u;
pack_col_ij_(1+na) = 1+floor(nx_u/n_pack);
pack_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nx_u=0:n_x_u-1;
x_u_pack_ = sparse(pack_row_ij_,pack_col_ij_,pack_val_ij_,n_x_u,n_x_u_pack);
a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u,n_x_u_pack]),[3,1,2]),[n_x_u_pack*n_x_u_pack,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),[3,1,2]);
clear a_x_u_load_;
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_0 = mean(diff(x_u_0_));
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_1 = mean(diff(x_u_1_));
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack); d_x_2 = mean(diff(x_u_2_));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
weight_xxx_u_ = d_x_0 * d_x_1 * d_x_2 ;
%%%%%%%%;
% Calculate moments. ;
%%%%%%%%;
a_rho_x_u_pack_ = a_x_u_pack_ + min(a_x_u_pack_,[],'all');
a_rho_x_c_0_avg = sum(x_u_0___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_avg = sum(x_u_1___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_avg = sum(x_u_2___.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_0_std = sum((x_u_0___ - a_rho_x_c_0_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_std = sum((x_u_1___ - a_rho_x_c_1_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_std = sum((x_u_2___ - a_rho_x_c_2_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_avg_ = [a_rho_x_c_0_avg ; a_rho_x_c_1_avg ; a_rho_x_c_2_avg];
a_rho_x_c_std_ = [a_rho_x_c_0_std ; a_rho_x_c_1_std ; a_rho_x_c_2_std];
disp(sprintf(' %% a_rho_x_c_std_ vs a_rho_x_c_avg_: %0.2f',fnorm(a_rho_x_c_std_)/fnorm(a_rho_x_c_avg_)));
if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
disp(sprintf(' %% Warning, molecule may not be well centered. Consider recentering.'));
end;%if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
%%%%%%%%;
% Possible to re-center. ;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_2_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
[K_u_0_,K_u_1_,K_u_2_] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_K_u = n_x_u_pack^3;
b_rho_x_u_pack_ = real(ifftn(fftn(a_rho_x_u_pack_).*exp(+i*2*pi*(K_u_0_*a_rho_x_c_0_avg + K_u_1_*a_rho_x_c_1_avg + K_u_2_*a_rho_x_c_2_avg))));
b_rho_x_c_0_avg = sum(x_u_0___.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_1_avg = sum(x_u_1___.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_2_avg = sum(x_u_2___.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_0_std = sum((x_u_0___ - b_rho_x_c_0_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_1_std = sum((x_u_1___ - b_rho_x_c_1_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_2_std = sum((x_u_2___ - b_rho_x_c_2_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_avg_ = [b_rho_x_c_0_avg ; b_rho_x_c_1_avg ; b_rho_x_c_2_avg];
b_rho_x_c_std_ = [b_rho_x_c_0_std ; b_rho_x_c_1_std ; b_rho_x_c_2_std];
disp(sprintf(' %% b_rho_x_c_std_ vs b_rho_x_c_avg_: %0.2f',fnorm(b_rho_x_c_std_)/fnorm(b_rho_x_c_avg_)));
%%%%%%%%;
a_x_u_base_ = a_x_u_pack_;
if flag_center_volume;
disp(sprintf(' %% centering volume'));
a_x_u_base_ = b_rho_x_u_pack_;
end;%if flag_center_volume;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/a_x_u_base_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); isosurface_f_x_c_0(a_x_u_pack_,98.5); title('a packed');
subplot(2,2,2); isosurface_f_x_c_0(a_x_u_pack_,[97.5,98.5,99.5]); title('a packed');
subplot(2,2,3); isosurface_f_x_c_0(b_rho_x_u_pack_,98.5); title('b packed');
subplot(2,2,4); isosurface_f_x_c_0(b_rho_x_u_pack_,[97.5,98.5,99.5]); title('b packed');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','n_x_u_pack','n_pack','pack_row_ij_','pack_col_ij_','pack_val_ij_','x_u_pack_' ...
     ,'x_u_0_','x_u_1_','x_u_2_','x_u_0___','x_u_1___','x_u_2___','n_x_u','n_xxx_u','weight_xxx_u_','a_x_u_base_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
dx_pack = diameter_x_c/n_x_u_pack; %dx_pack = mean(diff(x_u_0_));
[x_u_0__,x_u_1__] = ndgrid(x_u_0_,x_u_1_);
dxx_pack = dx_pack*dx_pack;%dxx_pack = mean(diff(x_u_0_))*mean(diff(x_u_1_));
dxxx_pack = dx_pack*dx_pack*dx_pack;%dxxx_pack = mean(diff(x_u_0_))*mean(diff(x_u_1_))*mean(diff(x_u_2_));

%%%%%%%%;
% simple visualization of a_base. ;
%%%%%%%%
flag_plot=1;
if flag_plot;
fname_fig_pre = sprintf('%s_jpg/a_x_u_base_prct_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
prctile_ = [94.0:0.5:99.5]; n_prctile = numel(prctile_);
figure(1);clf;figbig;
p_row = 3; p_col = 4; ns=0;
for nprctile=0:n_prctile-1;
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_p = prctile_(1+nprctile);
isosurface_f_x_c_0(a_x_u_base_,tmp_p);
title(sprintf('a base %.1f',tmp_p));
end;%for nprctile=0:n_prctile-1;
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_C2 = 'C2';
flag_uniform_over_n_k_p_r = 1;
flag_uniform_over_polar_a = 0; %<-- This is set to match test_ssnll_from_a_k_p_15 ;
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
,str_C2 ...
,flag_uniform_over_n_k_p_r ...
,flag_uniform_over_polar_a ...
) ;
%%%%;
% and generate consistant quadrature-grid on shell. ;
%%%%;
qref_k_eq_d = k_eq_d_double/(2*pi)/max(1e-12,k_p_r_max);
tmp_t = tic();
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
,str_C2 ...
,flag_uniform_over_polar_a ...
) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% sample_shell_6: %0.2fs',tmp_t)); end;
assert(fnorm(qref_azimu_b_shell_-k_p_azimu_b_qk_(1:qref_n_shell))<1e-12);
assert(fnorm(qref_polar_a_shell_-k_p_polar_a_qk_(1:qref_n_shell))<1e-12);
assert(fnorm(sum(qref_weight_shell_)-4*pi)<1e-12);
n_q = qref_n_shell;
weight_3d_k_p_qk__ = reshape(weight_3d_k_p_qk_,[n_q,n_k_p_r]);
weight_shell_qk__ = reshape(weight_shell_qk_,[n_q,n_k_p_r]);
assert(fnorm(4*pi*weight_3d_k_p_r_ - reshape(sum(weight_3d_k_p_qk__,1),[1,n_k_p_r]))<1e-12);
%%%%%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; n_plot = p_row*p_col;
for nplot=0:n_plot-1;
nk_p_r = max(0,min(n_k_p_r-1,round(n_k_p_r*nplot/n_plot)));
tmp_index_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
subplot(p_row,p_col,1+nplot);
plot3(k_c_0_qk_(1+tmp_index_),k_c_1_qk_(1+tmp_index_),k_c_2_qk_(1+tmp_index_),'.');
axis equal; axis vis3d; axisnotick3d;
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
end;%for nplot=0:n_plot-1;
end;%if flag_disp;
%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in = n_w_max; n_w_0in_ = n_w_max*ones(n_k_p_r,1);
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

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_p_quad_.mat',dir_ssnll);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_base_(:).*weight_xxx_u_(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_qk,2*pi*k_c_0_qk_*eta,2*pi*k_c_1_qk_*eta,2*pi*k_c_2_qk_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_p_qk_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_reco_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_reco_(:).*weight_xxx_u_(:),-1,1e-12,n_qk,2*pi*k_c_0_qk_/eta,2*pi*k_c_1_qk_/eta,2*pi*k_c_2_qk_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
a_x_u_base_l2 = sum(abs(a_x_u_base_).^2.*weight_xxx_u_,'all');
if (flag_verbose>0); disp(sprintf(' %% a_x_u_base_l2: %+0.8f',a_x_u_base_l2)); end;
a_k_p_quad_l2 = sum(abs(a_k_p_quad_).^2.*weight_3d_k_p_qk_,'all');
if (flag_verbose>0); disp(sprintf(' %% a_k_p_quad_l2: %+0.8f',a_k_p_quad_l2)); end;
a_x_u_reco_l2 = sum(abs(a_x_u_reco_).^2.*weight_xxx_u_,'all');
if (flag_verbose>0); disp(sprintf(' %% a_x_u_reco_l2: %+0.8f',a_x_u_reco_l2)); end;
a_k_p_reco_l2 = sum(abs(a_k_p_reco_).^2.*weight_3d_k_p_qk_,'all');
if (flag_verbose>0); disp(sprintf(' %% a_k_p_reco_l2: %+0.8f',a_k_p_reco_l2)); end;
%%%%;
disp(sprintf(' %% xxnufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_base_(:)-a_x_u_reco_)/fnorm(a_x_u_base_(:))));
disp(sprintf(' %% at this point one should ensure that a_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
save(fname_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_qk','n_qk_csum_' ...
     ,'k_p_r_qk_','k_p_azimu_b_qk_','k_p_polar_a_qk_' ...
     ,'weight_3d_k_p_qk_','weight_shell_qk_' ...
     ,'n_k_p_r','k_p_r_' ...
     ,'weight_3d_k_p_r_' ...
     ,'k_c_0_qk_','k_c_1_qk_','k_c_2_qk_' ...
     ,'a_k_p_quad_' ...
     ,'a_x_u_reco_' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/a_k_p_quad_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;
clf;
plot(k_p_r_qk_,log10(abs(a_k_p_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
figbig;
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

a_k_p_qk_ = a_k_p_quad_ ;

%%%%%%%%;
% generate templates S_k_p_wk_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/S_k_p_wkS__.mat',dir_ssnll);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic();
viewing_k_eq_d = t_eq_d_double/k_p_r_max;
template_k_eq_d = t_eq_d_double/k_p_r_max;
str_L = 'L'; %<-- Use C2 to allow for off-grid interpolation. ;
flag_tensor_vs_adap = 0; %<-- use adaptive grid for ddssnll_from_a_k_p, but need tensor grid for ddsssnll_3. ;
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,viewing_k_c_0_S_ ...
,viewing_k_c_1_S_ ...
,viewing_k_c_2_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,template_k_eq_d ...
,str_C2 ...
,flag_tensor_vs_adap ...
) ;
n_S = n_viewing_S;
viewing_gamma_z_S_ = zeros(n_S,1);
if (flag_verbose>0); disp(sprintf(' %% n_S %d, n_viewing_polar_a %d, n_viewing_azimu_b [%d,..,%d]',n_S,n_viewing_polar_a,n_viewing_azimu_b_(1+0),n_viewing_azimu_b_(end))); end;
%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% generate templates using it7. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%;
tmp_t = tic();
if ~exist('wS_from_single_shell_sba__','var'); wS_from_single_shell_sba__=[]; end;
if ~exist('dwSda_from_single_shell_sba__','var'); dwSda_from_single_shell_sba__=[]; end;
if ~exist('dwSdb_from_single_shell_sba__','var'); dwSdb_from_single_shell_sba__=[]; end;
if ~exist('ddwSdaa_from_single_shell_sba__','var'); ddwSdaa_from_single_shell_sba__=[]; end;
if ~exist('ddwSdab_from_single_shell_sba__','var'); ddwSdab_from_single_shell_sba__=[]; end;
if ~exist('ddwSdbb_from_single_shell_sba__','var'); ddwSdbb_from_single_shell_sba__=[]; end;
if ~exist('n_R','var'); n_R = []; end;
if ~exist('R_use_R___','var'); R_use_R___ = []; end;
if ~exist('a_R_k_p_Rqk__','var'); a_R_k_p_Rqk__=[]; end;
if ~exist('ba_from_single_shell_Rbaba___','var'); ba_from_single_shell_Rbaba___=[]; end;
if ~exist('wS_from_R_single_shell_Rsba___','var'); wS_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwSda_from_R_single_shell_Rsba___','var'); dwSda_from_R_single_shell_Rsba___=[]; end;
if ~exist('dwSdb_from_R_single_shell_Rsba___','var'); dwSdb_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdaa_from_R_single_shell_Rsba___','var'); ddwSdaa_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdab_from_R_single_shell_Rsba___','var'); ddwSdab_from_R_single_shell_Rsba___=[]; end;
if ~exist('ddwSdbb_from_R_single_shell_Rsba___','var'); ddwSdbb_from_R_single_shell_Rsba___=[]; end;
parameter_interpolate = struct('type','parameter');
parameter_interpolate.flag_verbose = 0;
parameter_interpolate.flag_check = 0;
[ ...
 parameter_interpolate ...
,template_from_a_k_p_wkS__ ...
,n_w ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,n_S_sub ...
,index_nS_from_nS_sub_ ...
,n_R ...
,R_use_R___ ...
,viewing_R_polar_a_mod_S_subR__ ...
,viewing_R_azimu_b_mod_S_subR__ ...
,nR_from_nS_sub_ ...
,a_R_k_p_Rqk__ ...
,ba_from_single_shell_Rbaba___ ...
,wS_from_R_single_shell_Rsba___ ...
,dtemplateda_from_a_k_p_wkS__ ...
,dtemplatedb_from_a_k_p_wkS__ ...
,dtemplatedc_from_a_k_p_wkS__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,dwSda_from_R_single_shell_Rsba___ ...
,dwSdb_from_R_single_shell_Rsba___ ...
,ddtemplatedaa_from_a_k_p_wkS__ ...
,ddtemplatedab_from_a_k_p_wkS__ ...
,ddtemplatedac_from_a_k_p_wkS__ ...
,ddtemplatedbb_from_a_k_p_wkS__ ...
,ddtemplatedbc_from_a_k_p_wkS__ ...
,ddtemplatedcc_from_a_k_p_wkS__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,ddwSdaa_from_R_single_shell_Rsba___ ...
,ddwSdab_from_R_single_shell_Rsba___ ...
,ddwSdbb_from_R_single_shell_Rsba___ ...
] = ...
interpolate_template_7( ...
 parameter_interpolate ...
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
,-1 ...
,n_w_0in ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,viewing_gamma_z_S_ ...
,wS_from_single_shell_sba__ ...
,dwSda_from_single_shell_sba__ ...
,dwSdb_from_single_shell_sba__ ...
,ddwSdaa_from_single_shell_sba__ ...
,ddwSdab_from_single_shell_sba__ ...
,ddwSdbb_from_single_shell_sba__ ...
,n_R ...
,R_use_R___ ...
,a_R_k_p_Rqk__ ...
,ba_from_single_shell_Rbaba___ ...
,wS_from_R_single_shell_Rsba___ ...
,dwSda_from_R_single_shell_Rsba___ ...
,dwSdb_from_R_single_shell_Rsba___ ...
,ddwSdaa_from_R_single_shell_Rsba___ ...
,ddwSdab_from_R_single_shell_Rsba___ ...
,ddwSdbb_from_R_single_shell_Rsba___ ...
);
S_k_p_wkS__ = template_from_a_k_p_wkS__;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_it7_k_p_wkS__ (interpolate_template_7): time %0.6fs',tmp_t)); end;
save(fname_mat ...
     ,'n_w_max','template_k_eq_d' ...
     ,'S_k_p_wkS__' ...
     ,'n_w_' ...
     ,'weight_2d_k_p_r_' ...
     ,'weight_2d_wk_' ...
     ,'n_viewing_S' ...
     ,'viewing_azimu_b_S_' ...
     ,'viewing_polar_a_S_' ...
     ,'viewing_weight_S_' ...
     ,'viewing_k_c_0_S_' ...
     ,'viewing_k_c_1_S_' ...
     ,'viewing_k_c_2_S_' ...
     ,'n_viewing_polar_a' ...
     ,'viewing_polar_a_' ...
     ,'n_viewing_azimu_b_' ...
     ,'n_S','n_w_sum','n_w_csum_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/S_k_p_wkS__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p_wkS__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wkS__(:,1+nS)),Slim_,colormap_80s);
axis equal; axisnotick; title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('Sample S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
% Now convert templates to S_k_q_wkS__. ;
%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/S_k_q_wkS__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p_wkS__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(S_k_q_wkS__(:,1+nS)),Slim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_q_wkS__(:,1+nS))/max(abs(S_k_q_wkS__(:)))),[-4,0],colormap_80s);
title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/M_k_p_wkM__.mat',dir_ssnll_ori);
if ~exist(fname_mat,'file'); fname_mat = sprintf('%s_mat/M_k_p_wkM__.mat',dir_ssnll); end;
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_M = 1024*8;
[ ...
 M_x_c_xxM___ ...
,index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_1( ...
 dir_data_star ...
,fname_nopath_star ...
,n_M ...
);
if (fnorm(Voltage_CTF_)< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_ = 300*ones(n_M,1); end;
if flag_invert; M_x_c_xxM___ = -M_x_c_xxM___; end;
%%%%%%%%;
% Remove any edge artefacts, mean center and normalize each image. ;
%%%%%%%%;
disp(sprintf(' %% Removing edge-artefacts'));
n_M_ext_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[M_x_c_xxM___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c_xxM___(:,:,1+nM),4,0.5,2,0);
end;%for nM=0:n_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
n_x_M_u = size(M_x_c_xxM___,1);
assert(n_x_M_u==size(M_x_c_xxM___,2));
disp(sprintf(' %% typical edge-artefact covers %0.6f = (%0.6f)^2 of image.',median(n_M_ext_/n_x_M_u^2),median(sqrt(n_M_ext_/n_x_M_u^2))));
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_M_u-1]/n_x_M_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
x_p_r__ = sqrt(x_c_0__.^2 + x_c_1__.^2);
x_c_mask__ = x_p_r__<=half_diameter_x_c;
M_abs_x_c_0_avg_ = zeros(n_M,1);
M_abs_x_c_1_avg_ = zeros(n_M,1);
M_mask_abs_x_c_0_avg_ = zeros(n_M,1);
M_mask_abs_x_c_1_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
M_abs_x_c_ = abs(squeeze(M_x_c_xxM___(:,:,1+nM))); %<-- no mask. ;
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nM) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nM) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
M_mask_abs_x_c_ = abs(squeeze(M_x_c_xxM___(:,:,1+nM)).*x_c_mask__); %<-- radial mask. ;
M_mask_abs_avg = mean(M_mask_abs_x_c_,'all');
M_mask_abs_x_c_0_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_0__,'all');
M_mask_abs_x_c_1_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_1__,'all');
M_mask_abs_x_c_0_avg_(1+nM) = M_mask_abs_x_c_0_avg;
M_mask_abs_x_c_1_avg_(1+nM) = M_mask_abs_x_c_1_avg;
clear M_mask_abs_x_c_;
end;%for nM=0:n_M-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/M_abs_x_c_avg_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figmed;
%%%%;
subplot(1,2,1);
plot(M_abs_x_c_0_avg_,M_abs_x_c_1_avg_,'.');
xlabel('M_abs_x_c_0_avg_','Interpreter','none');
ylabel('M_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('M_abs_x_c_ (no radial mask)','Interpreter','none');
%%%%;
subplot(1,2,2);
plot(M_mask_abs_x_c_0_avg_,M_mask_abs_x_c_1_avg_,'.');
xlabel('M_mask_abs_x_c_0_avg_','Interpreter','none');
ylabel('M_mask_abs_x_c_1_avg_','Interpreter','none');
axis equal; grid on;
title('M_mask_abs_x_c_ (max radial mask)','Interpreter','none');
%%%%;
tmp_corr_ = corr([M_abs_x_c_0_avg_,M_abs_x_c_1_avg_],[M_mask_abs_x_c_0_avg_,M_mask_abs_x_c_1_avg_]);
sgtitle(sprintf('correlation x,y = ( %0.4f , %0.4f )',tmp_corr_(1,1),tmp_corr_(2,2)));
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
% Now convert images to M_k_p_wkM__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
M_k_p_wkM__ = zeros(n_w_sum,n_M);
N_k_p_wkM__ = zeros(n_w_sum,n_M);
image_center_delta_x_c_0_ = zeros(n_M,1);
image_center_delta_x_c_1_ = zeros(n_M,1);
n_x_M_center = max(n_w_max,n_x_M_u/4);
if (flag_verbose>0); disp(sprintf(' %% n_x_M_center %d',n_x_M_center)); end;
M_x_c_avg_xx__ = zeros(n_x_M_center,n_x_M_center);
M_x_c_std_xx__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_avg_xx__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_std_xx__ = zeros(n_x_M_center,n_x_M_center);
for nM=0:n_M-1;
if (mod(nM,128)==0); if (flag_verbose>0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
M_x_c_xx__ = squeeze(M_x_c_xxM___(:,:,1+nM));
M_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_xx__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_wk_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
% Now *ALSO* center images after filtering/masking/thresholding: ;
parameter = struct('type','parameter');
[ ...
 parameter ...
,N_k_p_ ...
,M_x_c_0in_ ...
,M_x_c_out_ ...
,tmp_delta_x_c_0 ...
,tmp_delta_x_c_1 ...
] = ...
image_center_1( ...
 parameter ...
,n_x_M_center ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_ ...
,weight_2d_wk_ ...
);
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
N_k_p_wkM__(:,1+nM) = N_k_p_;
image_center_delta_x_c_0_(1+nM) = tmp_delta_x_c_0;
image_center_delta_x_c_1_(1+nM) = tmp_delta_x_c_1;
M_x_c_avg_xx__ = M_x_c_avg_xx__ + reshape(real(M_x_c_0in_),[n_x_M_center,n_x_M_center]);
M_x_c_std_xx__ = M_x_c_std_xx__ + reshape(real(M_x_c_0in_).^2,[n_x_M_center,n_x_M_center]);
N_x_c_avg_xx__ = N_x_c_avg_xx__ + reshape(real(M_x_c_out_),[n_x_M_center,n_x_M_center]);
N_x_c_std_xx__ = N_x_c_std_xx__ + reshape(real(M_x_c_out_).^2,[n_x_M_center,n_x_M_center]);
end;%for nM=0:n_M-1;
M_x_c_avg_xx__ = M_x_c_avg_xx__/n_M; M_x_c_std_xx__ = M_x_c_std_xx__/n_M; M_x_c_std_xx__ = sqrt(M_x_c_std_xx__ - M_x_c_avg_xx__.^2);
N_x_c_avg_xx__ = N_x_c_avg_xx__/n_M; N_x_c_std_xx__ = N_x_c_std_xx__/n_M; N_x_c_std_xx__ = sqrt(N_x_c_std_xx__ - N_x_c_avg_xx__.^2);
%%%%%%%%;
save(fname_mat ...
     ,'n_M','n_x_M_u','Pixel_Spacing','x_c_0_','x_c_1_','x_c_0__','x_c_1__' ...
     ,'n_M_ext_' ...
     ,'M_abs_x_c_0_avg_','M_abs_x_c_1_avg_' ...
     ,'M_mask_abs_x_c_0_avg_','M_mask_abs_x_c_1_avg_' ...
     ,'image_center_delta_x_c_0_','image_center_delta_x_c_1_' ...
     ,'M_k_p_wkM__','N_k_p_wkM__' ...
     ,'M_x_c_avg_xx__','M_x_c_std_xx__' ...
     ,'N_x_c_avg_xx__','N_x_c_std_xx__' ...
     ,'index_nCTF_from_nM_' ...
     ,'index_nM_from_nCTF_' ...
     ,'Voltage_CTF_' ...
     ,'DefocusU_CTF_' ...
     ,'DefocusV_CTF_' ...
     ,'DefocusAngle_CTF_' ...
     ,'SphericalAberration_CTF_' ...
     ,'AmplitudeContrast_CTF_' ...
     );
clear M_x_c_xxM___;
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/M_all_center_FIGB',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;figbig;figbeach();
subplot(1,2,1); imagesc(M_x_c_std_xx__); axis image; axisnotick;
title('M_x_c_std_xx__','Interpreter','none');
subplot(1,2,2); imagesc(N_x_c_std_xx__); axis image; axisnotick;
title('N_x_c_std_xx__','Interpreter','none');
sgtitle(sprintf(' %% centering filtered images'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));

%%%%%%%%;
% Now convert images to M_k_q_wkM__. ;
%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_wkM__ time %0.2fs',tmp_t)); end;
tmp_t = tic();
N_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
N_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,N_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_wkM__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/M_k_q_wkM__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
Mlim_ = std(abs(M_k_p_wkM__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nM = max(0,min(n_M-1,floor(n_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(M_k_q_wkM__(:,1+nM)),Mlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(M_k_q_wkM__(:,1+nM))/max(abs(M_k_q_wkM__(:)))),[-4,0],colormap_80s);
title(sprintf('nM %d',nM));
end;%for nl=0:15-1;
sgtitle(sprintf('M_k_q_wkM__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/M_k_q_norm_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
M_k_q_norm_ = sqrt(sum(abs(M_k_q_wkM__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(M_k_q_norm_/max(M_k_q_norm_)),[-4,0],colormap_80s);
title('M_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
M_k_q_norm_ = sqrt(sum(abs(M_k_q_wkM__).^2,2)); tmp_eps = M_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,M_k_q_norm_-tmp_eps)/max(M_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('M_k_q_norm_ - tmp_eps','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/N_k_q_wkM__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
Nlim_ = std(abs(N_k_p_wkM__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nM = max(0,min(n_M-1,floor(n_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(N_k_q_wkM__(:,1+nM)),Nlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(N_k_q_wkM__(:,1+nM))/max(abs(N_k_q_wkM__(:)))),[-4,0],colormap_80s);
title(sprintf('nM %d',nM));
end;%for nl=0:15-1;
sgtitle(sprintf('N_k_q_wkM__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/N_k_q_norm_',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
N_k_q_norm_ = sqrt(sum(abs(N_k_q_wkM__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(N_k_q_norm_/max(N_k_q_norm_)),[-4,0],colormap_80s);
title('N_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
N_k_q_norm_ = sqrt(sum(abs(N_k_q_wkM__).^2,2)); tmp_eps = N_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,N_k_q_norm_-tmp_eps)/max(N_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('N_k_q_norm_ - tmp_eps','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate CTF functions. ;
%{
Compare niko_ctf with section 28 from ; ;
https://link.springer.com/book/10.1007/978-0-387-76501-3 ;
The relevant formulae are 28.3 and 28.4
(which define the spatial-frequency 'u' in terms of inverse-angstroms -- without a 2*pi built in) ;
and then later on in 28.33 and 28.34,
where the magnitude of the spatial-frequency 'u' is multiplied by lambda to provide a non-dimensional 'angle' (see 28.32). ;
With this definition it seems as though 'u' is the wavenumber -- i.e., the number of full wavelengths per unit-length. ;
So, in the matlab-code below, a wave-number of k_p_r_max = 48/(2*pi) then corresponds to a wave-number of 48, ;
with the correction that -- since our box is of length 2 -- 
there are actually k_p_r_max*2 = 96 full waves in the box when k_p_r = 48/(2*pi). ;
 ;
Now when we call niko_ctf below, ;
we set tmp_k_c_1 and tmp_k_c_2 to equal (2*pi)*k_p_r*cos(tmp_theta) and (2*pi)*k_p_r*sin(tmp_theta), respectively. ;
This corresponds to passing in essentially (2*pi)*k_p_r; ;
i.e., when k_p_r = k_p_r_max we pass in the number 48. ;
Now within niko_ctf this (2*pi)*k_p_r gets multiplied by 'thetatr': ;
thetatr = ' CTF_lambda / Box_size_in_angstroms / pi ' (see line 941 and 949) ;
ultimately producing: ;
(2*pi) * k_p_r * thetatr = 2 * k_p_r * CTF_lambda / Box_size_in_angstroms . ;
Now CTF_Lambda is the electron-wavelength in angstroms (see line 931) . ;
Thus, within niko_ctf, the product ;
angle = rad*thetatr = u*lambda, where ;
lambda = electron-wavelength in angstroms ;
u = 2 * k_p_r / Box_size_in_angstroms = wavenumber in inverse-angstroms ;
Where, as we discussed above, ;
the extra factor of 2 is because we define our wavenumbers assuming that the box is side-length 2 (rather than 1). ;
%}
%%%%%%%%;
fname_mat = sprintf('%s_mat/CTF_k_p_wkC__.mat',dir_ssnll);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_CTF = numel(index_nM_from_nCTF_);
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if (mod(nCTF,100)==0); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_Spherical_Aberration = SphericalAberration_CTF_(1+nCTF);% spherical aberration of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = Voltage_CTF_(1+nCTF);% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = DefocusU_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_V = DefocusV_CTF_(1+nCTF);% defocus values (in Angstroms) ;
CTF_Defocus_Angle = DefocusAngle_CTF_(1+nCTF);% angle of astigmatism ;
%CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
CTF_Amplitude_Contrast = AmplitudeContrast_CTF_(1+nCTF);% CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in CTF ;
tmp_w2=CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in CTF ;
%  CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification;
CTF_Object_Pixel_Size = Pixel_Spacing;% pixel size of the scanner in physical space in Angstroms ;
CTF_lambda_per_box = CTF_lambda/(n_x_M_u*CTF_Object_Pixel_Size);% n_x_M_u*CTF_Object_Pixel_Size is the box size in Angstroms ;
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/n_w_(1+nk);
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle,CTF_lambda_per_box/pi,tmp_k_c_1,tmp_k_c_2);
CTF_k_p_wkC__(1+na,1+nCTF) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r_kC__(1+nk_p_r,1+nCTF) = mean(CTF_k_p_wkC__(1+tmp_index_,1+nCTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_wk_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_k_(1+nk_p_r) = mean(CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_kM__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1:n_M));
SCTF_ = svd(CTF_k_p_r_kM__);
n_CTF_rank = max(find(SCTF_/max(SCTF_)>tolerance_master));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r_kM__,n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;
% Now plot out some of the CTF-functions for varying anisotropy. ;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/CTF_k_p_sample__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;
[tmp_anisotropy_,CTF_anisotropy_index_] = sort(DefocusU_CTF_ - DefocusV_CTF_,'ascend'); CTF_anisotropy_index_ = CTF_anisotropy_index_ - 1;
for nl=0:15-1;
subplot(3,5,1+nl);
tmp_nCTF = max(0,min(n_CTF-1,floor(n_CTF*nl/(15-1)))); nCTF = CTF_anisotropy_index_(1+tmp_nCTF);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,CTF_k_p_wkC__(:,1+nCTF),[-1,+1],colormap_beach());
axis image; axisnotick;
title(sprintf('nCTF %d anisotropy %0.2f',nCTF,tmp_anisotropy_(1+tmp_nCTF)));
end;%for nl=0:15-1;
sgtitle(sprintf('CTF_k_p_wkC__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
% Now determine the CTF cross correlation. ;
% This depends  on index_nCTF_from_nM_. ;
%%%%%%%%;
tmp_CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))),2);
tmp_CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_k_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xavg__ = tmp_CTF_avg_k_p_r_k_ * transpose(tmp_CTF_avg_k_p_r_k_);
CTF_k_p_r_xcor__ = CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1))) * transpose(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+(0:n_M-1)))) / n_M;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/CTF_k_p_xcor__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figmed;
figbeach();
subplot(1,2,1); imagesc(CTF_k_p_r_xavg__); axis image; axisnotick; title('CTF_k_p_r_xavg__','Interpreter','none');
subplot(1,2,2); imagesc(CTF_k_p_r_xcor__); axis image; axisnotick; title('CTF_k_p_r_xcor__','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
save(fname_mat ...
     ,'n_CTF' ...
     ,'index_nCTF_from_nM_' ...
     ,'CTF_k_p_wkC__' ...
     ,'CTF_k_p_r_kC__' ...
     ,'CTF_k_p_r_kM__' ...
     ,'n_CTF_rank' ...
     ,'SCTF_','UCTF_kc__','VSCTF_Mc__' ...
     ,'CTF_avg_k_p_wk_' ...
     ,'CTF_avg_k_p_r_k_' ...
     ,'CTF_k_p_r_xavg__' ...
     ,'CTF_k_p_r_xcor__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/CTF_k_p_r_xxxx__',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;
clf;
colormap(colormap_beach());
subplot(1,2,1);
imagesc(CTF_k_p_r_xavg__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('CTF_k_p_r_xavg__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
subplot(1,2,2);
imagesc(CTF_k_p_r_xcor__,[0,1]);
axis image;
set(gca,'XTick',1:n_k_p_r,'XTickLabel',k_p_r_); xtickangle(90);
set(gca,'YTick',1:n_k_p_r,'YTickLabel',k_p_r_);
title('CTF_k_p_r_xcor__ cross correlation','Interpreter','none');
xlabel('k_p_r','Interpreter','none');
ylabel('k_p_r','Interpreter','none');
figbig;
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% First cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
%%%%%%%%;
n_cluster = 1+max(index_ncluster_from_nCTF_);
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/VSCTF_Mc_scatter_FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;
clf;figsml;figbeach();
markersize_use = 32;
scatter(VSCTF_Mc__(:,1+0),VSCTF_Mc__(:,1+min(1,n_CTF_rank-1)),markersize_use,index_ncluster_from_nM_,'filled');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate the empirical principal-modes for each of the nCTF. ;
%%%%%%%%;
X_2d_Nemp_d1_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_Nemp_d1_weight_rc__ = zeros(n_k_p_r,n_cluster);
X_2d_Memp_d1_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_Memp_d1_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
[X_2d_Nemp_d1_kk__,X_2d_Nemp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_index_nM_from_ncluster,N_k_p_wkM__(:,1+index_nM_from_ncluster_));
X_2d_Nemp_d1_kkc___(:,:,1+ncluster) = X_2d_Nemp_d1_kk__;
X_2d_Nemp_d1_weight_rc__(:,1+ncluster) = X_2d_Nemp_d1_weight_r_;
clear X_2d_Nemp_d1_kk__ X_2d_Nemp_d1_weight_r_ ;
[X_2d_Memp_d1_kk__,X_2d_Memp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_index_nM_from_ncluster,M_k_p_wkM__(:,1+index_nM_from_ncluster_));
X_2d_Memp_d1_kkc___(:,:,1+ncluster) = X_2d_Memp_d1_kk__;
X_2d_Memp_d1_weight_rc__(:,1+ncluster) = X_2d_Memp_d1_weight_r_;
clear X_2d_Memp_d1_kk__ X_2d_Memp_d1_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/X_2d_Memp_d1_kkc___FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(real(X_2d_Memp_d1_kkc___(:,:,1+ncluster))); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/X_2d_Nemp_d1_kkc___FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(real(X_2d_Nemp_d1_kkc___(:,:,1+ncluster))); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;


%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_Y_quad_yk_.mat',dir_ssnll_ori);
if (~exist(fname_mat,'file')); disp(sprintf(' %% Warning, %s not found',fname_mat)); end;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Now calculate the idealized principal-modes for each nCTF. ;
%%%%%%%%;
X_2d_xavg_d0_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_d0_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
tmp_CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
tmp_CTF_k_p_r_xavg_kk__ = tmp_CTF_k_p_r_xavg_k_*transpose(tmp_CTF_k_p_r_xavg_k_);
tmp_delta_sigma = 0;
[X_2d_xavg_d0_kk__,X_2d_xavg_d0_weight_r_] = principled_marching_cost_matrix_6(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_yk_,tmp_CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
X_2d_xavg_d0_kkc___(:,:,1+ncluster) = X_2d_xavg_d0_kk__;
X_2d_xavg_d0_weight_rc__(:,1+ncluster) = X_2d_xavg_d0_weight_r_;
clear X_2d_xavg_d0_kk__ X_2d_xavg_d0_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/X_2d_xavg_d0_kkc___FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(real(X_2d_xavg_d0_kkc___(:,:,1+ncluster))); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate the idealized principal-modes for unit CTF. ;
%%%%%%%%;
X_2d_x1_d0_kk__ = zeros(n_k_p_r,n_k_p_r);
X_2d_x1_d0_weight_r_ = zeros(n_k_p_r,1);
[X_2d_x1_d0_kk__,X_2d_x1_d0_weight_r_] = principled_marching_cost_matrix_6(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_yk_);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/X_2d_x1_d0_kk___FIGA',dir_ssnll);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
figure(1+nf);nf=nf+1;clf;figsml;figbeach();
subplot(1,1,1);
imagesc(real(X_2d_x1_d0_kk__)); axis image; axisnotick;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now use ideal principal-modes to align images to volume. ;
%%%%%%%%;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
flag_N_vs_M = flag_center_image;
if flag_N_vs_M==0; tmp_N_k_p_wkM__ = M_k_p_wkM__; tmp_str = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_N_k_p_wkM__ = N_k_p_wkM__; tmp_str = 'N'; end;%if flag_N_vs_M==1;
%%%%;
if flag_center_image==1; 
M_k_p_wkM__ = N_k_p_wkM__;
M_k_q_wkM__ = N_k_q_wkM__;
end;%if flag_center_image==1; 
if flag_center_image==0;
N_k_p_wkM__ = M_k_p_wkM__;
N_k_q_wkM__ = M_k_q_wkM__;
end;%if flag_center_image==0;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_Y_reco_from_%s__.mat',dir_ssnll_ori,tmp_str);
if ~exist(fname_mat,'file'); disp(sprintf(' %% Warning, %s not found',fname_mat)); end;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
euler_polar_a_empi_tail_ = tmp_.euler_polar_a_Mi__(:,2:end);
euler_azimu_b_empi_tail_ = tmp_.euler_azimu_b_Mi__(:,2:end);
euler_gamma_z_empi_tail_ = tmp_.euler_gamma_z_Mi__(:,2:end);
euler_polar_a_empi_ = tmp_.euler_polar_a_Mi__(:,end);
euler_polar_a_empi_dif_ = periodize(bsxfun(@minus,tmp_.euler_polar_a_Mi__(:,2:end),euler_polar_a_empi_),-pi/2,+pi/2); %<-- not symmetry. ;
euler_polar_a_tavg_ = periodize( euler_polar_a_empi_ + mean(euler_polar_a_empi_dif_(:,end-3:end-0),2) ,0,1*pi);
euler_azimu_b_empi_ = tmp_.euler_azimu_b_Mi__(:,end);
euler_azimu_b_empi_dif_ = periodize(bsxfun(@minus,tmp_.euler_azimu_b_Mi__(:,2:end),euler_azimu_b_empi_),-pi/4,+pi/4); %<-- 4-fold symmetry. ;
euler_azimu_b_tavg_ = periodize( euler_azimu_b_empi_ + mean(euler_azimu_b_empi_dif_(:,end-3:end-0),2) ,0,2*pi);
euler_gamma_z_empi_ = tmp_.euler_gamma_z_Mi__(:,end);
euler_gamma_z_empi_dif_ = periodize(bsxfun(@minus,tmp_.euler_gamma_z_Mi__(:,2:end),euler_gamma_z_empi_),-pi/2,+pi/2); %<-- not symmetry. ;
euler_gamma_z_tavg_ = periodize( euler_gamma_z_empi_ + mean(euler_gamma_z_empi_dif_(:,end-3:end-0),2) ,0,2*pi);
image_delta_x_empi_ = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_empi_ = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_empi = tmp_.corr_a_k_Y_i_(end);
a_k_Y_reco_empi_yk_ = tmp_.a_k_Y_reco_yki__(:,end);
a_k_Y_reco_empi_yk__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_reco_empi_yk_);
image_X_value_empi_ = tmp_.image_X_value_Mi__(:,end-1);
clear tmp_;
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s_mat/X_TM_A_TM_.mat',dir_ssnll_ori);
if ~exist(fname_mat,'file'); disp(sprintf(' %% Warning, %s not found',fname_mat)); end;
if  exist(fname_mat,'file');
if (flag_verbose>0); disp(sprintf(' %% %s found, not creating',fname_mat)); end;
load(fname_mat);
end;%if  exist(fname_mat,'file');
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% mean(A_TM_M_): %0.6f',mean(A_TM_M_))); end;

%%%%%%%%;
% Simple low-accuracy check of transformations (note imperfectly matched spatial grid). ;
%%%%%%%%;
flag_check=1;
if flag_check;
T_x_c_ = exp(-(x_u_0__.^2 + x_u_1__.^2)/(2*(half_diameter_x_c/6).^2)); T_x_c_ = T_x_c_/max(1e-12,sqrt(sum(abs(T_x_c_).^2,'all')*dxx_pack));
T_x_c_l2 = sum(abs(T_x_c_).^2,'all')*dxx_pack;
T_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,T_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dxx_pack ;
T_k_p_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_T_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
tmp_T_x_c_l2 = sum(abs(tmp_T_x_c_).^2,'all')*dxx_pack;
tmp_T_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,tmp_T_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack^2)*dxx_pack ;
innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_T_k_p_,tmp_T_k_p_);
tmp_T_k_p_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_T_k_p_,tmp_T_k_p_);
if (flag_verbose>0); disp(sprintf(' %% T_x_c_l2: %0.6f',T_x_c_l2)); end;
if (flag_verbose>0); disp(sprintf(' %% T_k_p_l2: %0.6f',T_k_p_l2)); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_T_x_c_l2: %0.6f',tmp_T_x_c_l2)); end;
if (flag_verbose>0); disp(sprintf(' %% tmp_T_k_p_l2: %0.6f',tmp_T_k_p_l2)); end;
end;%if flag_check;
%%%%%%%%;
% Now the reconstruction a_k_Y_reco_empi_yk_ is designed so that the average 'dosage' a:=TM/TT is 1.0. ;
% That is to say, there is (on average) one aligned-template hidden inside each image. ;
% If we assume that the dosage a is TM/TT for each image, then the residual-norm is: ;
% RR = <M-aT,M-aT> = MM - 2*a*MT + a*a*TT = MM - 2*MT*MT/TT + MT*MT/TT = MM-MT*MT/TT. ;
% If, on the other hand, we assume that a = abar = 1.0, then the residual-norm is: ;
% RR = MM - 2*MT + TT. ;
%%%%%%%%;
A_avg = mean(A_TM_M_,'all');
RR_opt_M_ = MM_l2_M_ - (TM_M_).^2./max(1e-12,TT_l2_M_);
RR_opt = mean(RR_opt_M_,'all');
MM_bar = mean(MM_l2_M_,'all');
RR_bar_M_ = MM_l2_M_ - 2*TM_M_ + TT_l2_M_;
RR_bar = mean(RR_bar_M_,'all');
R_opt = mean(sqrt(RR_opt_M_),'all');
R_bar = mean(sqrt(RR_bar_M_),'all');
TT_bar = mean(TT_l2_M_);
T_bar = mean(sqrt(TT_l2_M_));
if (flag_verbose>0); disp(sprintf(' %% A_avg: %0.6f',A_avg)); end;
if (flag_verbose>0); disp(sprintf(' %% RR_opt: %0.6f',RR_opt)); end;
if (flag_verbose>0); disp(sprintf(' %% R_opt: %0.6f',R_opt)); end;
if (flag_verbose>0); disp(sprintf(' %% RR_bar: %0.6f',RR_bar)); end;
if (flag_verbose>0); disp(sprintf(' %% R_bar: %0.6f',R_bar)); end;
if (flag_verbose>0); disp(sprintf(' %% TT_bar: %0.6f',TT_bar)); end;
if (flag_verbose>0); disp(sprintf(' %% T_bar: %0.6f',T_bar)); end;
%%%%%%%%;
% This RR_bar is the typical l2-norm for the residual image. ;
% Let us assume that the overall molecular mass is 400kDa. ;
% Also assume that the molecule is roughly encased within a box of ~1/2 side-length. ;
% This means that the number of pixels holding the signal is ~n_x_u_pack^2/4 = 1000. ;
% Thus, there should be roughly 400kDa/1000 = 400Da of 'signal' packed into each of 1000 pixels in a synthetic image. ;
% Assuming each pixel is roughly 1AA, or 1 square-angstrom, we have 400 Da/AA packed into each pixel of the sythetic image support. ;
% Thus, the norm-squared of a typical template is roughly 1000AA*[400Da/AA]^2 = 16e7 DD/AA %<-- Daltons^2 per Angstrom^2. ;
% Now TT_bar=0.0042 corresponds to roughly 16e7 DD/AA. ;
% For the noise, we expect individual pixel-values of sigma_D/AA (with sigma_D ~ sigma*N(0,1) in Daltons). ;
% Now the norm-squared for the noise should be n_x_u_pack^2 AA * [sigma_D/AA]^2 = n_x_u_pack^2 * sigma^2 DD/AA. ;
% Thus RR_bar corresponds to 4e3 sigma^2 DD/AA. ;
% In this case RR_bar=0.054~13*TT_bar corresponds to 4e3 sigma^2 DD/AA. ;
% Thus, 13*16e7 DD/AA = 4e3 * sigma^2 DD/AA, implying that sigma^2 is about 52e4, or sigma_D is about 700 Daltons. ;
% This corresponds to a per-pixel snr of about 0.5 or so at this resolution. ;
% In our 'image-driven' units, sigma^2 is about RR_bar/n_x_u_pack^2 = 0.000013 = 1.3e-5 --> 1/sigma^2 = 7.7e4. ;
%%%%%%%%;
% Thus, in our 'image-driven' units (which have determined the magnitude of a_k_Y_reco_empi_yk_), ;
% we have a typical template with l2-norm of TT_bar=0.0042, corresponding to 16e7 DD/AA. ;
% Meanwhile,a  typical noise-image has l2-norm of RR_bar=0.054, corresponding to 13*16e7 DD/AA. ;
% The different between the l2-norm of a poorly-aligned image (say, MM_bar) and a well-aligned image (say, RR_bar) is ~0.058-0.054 = 0.004, ;
% which is obviously the l2-norm of a typical template (i.e., TT_bar~0.004). ;
% So a typical noise-level (sigma-squared) of 52e4 implies that the log-relative-likelihood between an aligned- and unaligned-image is: ;
% 16e7/52e4 ~ 300, which is tremendous. ;
%%%%%%%%;
% So now if we multiply any particular dvol (e.g., v_eig_dvol_k_Y_yk_) by a constant "A", ;
% then the formula: ;
% ddssnll(a + delta*dvol) = ddssnll(a) + 0.5* H(a;dvol) * delta^2 ;
% becomes: ;
% ddssnll(a + (delta/A)*A*dvol) = ddssnll(a) + 0.5* H(a;dvol)*A^2 * (delta/A)^2 ;
% implying that the hessian H is multiplied by A^2. ;
% So if we shrink dvol by multiplying by a small A, then the hessian decreases by A^2. ;
% I.e.: H(a;A*dvol) = A^2*H(a;dvol). ;
%%%%%%%%;
% Now, as an example for trpv1, we can examine nlsigma=4; index_lambda=2. ;
% In this case we have: ;
% lambda_dif = 0.36, lambda_lsq = 0.27. ;
% [~,~,tmp_v_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_v_eig_dvol_yk_) --> 14.06. ;
% [~,~,tmp_a_std] = spharm_normalize_2(pm_n_k_p_r,pm_k_p_r_,pm_weight_3d_k_p_r_,pm_l_max_,pm_a_k_Y_quad_yk_) --> 0.156. ;
% tmp_v_std/tmp_a_std --> 90. ;
% similarly: ;
% fnorm(v_x_u_reco_)/fnorm(tmp_a_x_u_reco_) --> 72. ;
% similarly: ;
% [pm_a_k_Y_quad_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_a_k_Y_quad_yk_,pm_a_k_Y_quad_yk_)) --> 0.0144 ;
% [pm_v_tilde_eig_dvol_lr] = sqrt(local_weightless_f_bar_dot_g_(pm_n_k_p_r,pm_l_max_,0,pm_v_tilde_eig_dvol_yk_,pm_v_tilde_eig_dvol_yk_)) --> 1.00 ;
% [pm_v_eig_dvol_lr] = sqrt(local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_v_eig_dvol_yk_,pm_v_eig_dvol_yk_)) --> 1.00 ;
% pm_v_eig_dvol_lr/pm_a_k_Y_quad_lr --> 69. ;
% So we can imagine multiplying pm_v_eig_dvol_yk_ by (1/90) to produce something with equivalent norm to pm_a_k_Y_quad_yk_. ;
% This then multiplies the value of H by (1/90^2) --> lambda_lsq*(tmp_a_std/tmp_v_std)^2 --> 3.4e-5 (in image-driven units). ;
% Note that sum(weight_imagecount_M_use_) --> 1.0, ;
% so our likelihood takes the form of:
% ssnll(a+delta*dvol) - ssnll(a) = N_image * (0.5 * 3.4e-5 * delta^2) ; %<-- here a and dvol have comparable l2-norm. ;
% nll(a+delta*dvol) - nll(a) = sigma^2 * N_image * (0.5 * 3.4e-5 * delta^2) ; %<-- delta could be thought of as dimensionless. ;
% nll(a+delta*dvol) - nll(a) = 7.7e4 * N_image * (0.5 * 3.4e-5 * delta^2) ;
% nll(a+delta*dvol) - nll(a) = 1.3 [bits/image] * N_image * delta^2 ;
% log(20) [bits] = 1.3 [bits/image] * [1000 images] * delta^2 --> 0.0023 = delta^2 --> 0.048 = delta. ;
% Meaning that about 0.048 (or 1/21) of the volume dvol can be added onto a before the new hypothesis can be rejected in favor of the data. ;
% Now, noting that a and dvol are basically perpendicular: ;
% [pm_a_dot_v_lr] = local_f_bar_dot_g_(pm_n_k_p_r,pm_weight_3d_riesz_k_p_r_,pm_l_max_,0,pm_a_k_Y_quad_yk_,pm_v_eig_dvol_yk_) --> 1.4e-5. ;
% pm_a_dot_v_lr/max(1e-12,pm_a_k_Y_quad_lr*pm_v_eig_dvol_lr) --> 1e-3. ;
% we can conclude that delta = 0.048 corresponds to a shift (on a unit-sphere encompassing pm_a_k_Y_quad_yk_ and pm_v_eig_dvol_yk_) of sin(theta) = 0.048. ;
% Thus, cos(theta) is roughly 1-0.5*delta^2 = 1-0.0023/2 = 0.9989 --> a loss of about 0.001 of the projection onto the full volume. ;
% By contrast, the minimal eigenvector in the 'empirically-sampled' or 'uniformly-sampled' case is about 10, ;
% meaning that an analogous calculation results in: ;
% nll(a+delta*dvol) - nll(a) = 7.7e4 * N_image * (0.5 * lambda_lsq / 90^2 * delta^2). ;
% 3.0 [bits] = 4 [bits/image/lambda] * [1000 images] * lambda_lsq * delta^2. ;
% sqrt(3/4e3/lambda_lsq) = delta --> 0.05 for lambda_lsq=0.27, but 0.0087 for lambda_lsq=10 ;
% Similarly, we get delta ~ 0.11 for lambda_lsq=0.06. ;
%%%%%%%%;
 
%%%%%%%%;
hist2dab_k_eq_d = 0.5/k_p_r_max;
[ ...
 n_hist2dab_S ...
,hist2dab_azimu_b_S_ ...
,hist2dab_polar_a_S_ ...
,hist2dab_weight_S_ ...
,hist2dab_k_c_0_S_ ...
,hist2dab_k_c_1_S_ ...
,hist2dab_k_c_2_S_ ...
,n_hist2dab_polar_a ...
,hist2dab_polar_a_ ...
,n_hist2dab_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,hist2dab_k_eq_d ...
,'L' ... %<-- exclude pole. ;
,0 ... %<-- adaptive grid. ;
) ;
%%%%;
if flag_disp;
fname_fig_pre = sprintf('%s_jpg/euler_tavg_FIGA',dir_ssnll_ori);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024+512,768]);
fontsize_use = 12;
subplot(1,3,1);
hlim_ = [];
c_use__ = colormap_80s;
flag_2d_vs_3d = 0;
flag_loghist_vs_hist = 0;
[ ...
 h_raw_ab_ ...
 h_w3d_ab_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 hist2dab_polar_a_S_ ...
,hist2dab_azimu_b_S_ ...
,hist2dab_weight_S_ ...
,euler_polar_a_tavg_(:) ...
,euler_azimu_b_tavg_(:) ...
,hlim_ ...
,c_use__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
xlabel('x');ylabel('y');zlabel('z');axisnotick3d;
title('empirical $\mu(\tau)$','Interpreter','latex');
colorbar(gca,'off');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,[2,3]);
hlim_ = [];
c_use__ = colormap_80s;
flag_2d_vs_3d = 1;
flag_loghist_vs_hist = 0;
hist2d_polar_a_azimu_b_0( ...
 hist2dab_polar_a_S_ ...
,hist2dab_azimu_b_S_ ...
,hist2dab_weight_S_ ...
,euler_polar_a_tavg_(:) ...
,euler_azimu_b_tavg_(:) ...
,hlim_ ...
,c_use__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
xlabel('azimu_b','Interpreter','none');
ylabel('polar_a','Interpreter','none');
axisnotick;
title('empirical $\mu(\tau)$','Interpreter','latex');
colorbar(gca);
set(gca,'FontSize',fontsize_use);
%%%%;
sgtitle(fname_fig_pre,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_pre));
end;%if ( exist(fname_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% generate refined templates hist2dab_S_k_p_wk_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
template_k_eq_d = -1;
tmp_t = tic();
[ ...
 hist2dab_S_k_p_wkS__ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max_max ...
,n_k_p_r ...
,a_k_Y_quad_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_hist2dab_S ...
,hist2dab_azimu_b_S_ ...
,hist2dab_polar_a_S_ ...
,hist2dab_weight_S_ ...
,n_hist2dab_polar_a ...
,hist2dab_polar_a_ ...
,n_hist2dab_azimu_b_ ...
);
hist2dab_S_k_p_wkS__ = reshape(hist2dab_S_k_p_wkS__,[n_w_sum,n_hist2dab_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% hist2dab_S_k_p_wkS__ (pm_template_2): %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct weight_3d_riesz')); end;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
weight_3d_riesz_qk_ = weight_3d_k_p_qk_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_qk_csum_(1+nk_p_r):n_qk_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_p_qk_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_p_qk_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_qk_(1+tmp_index_) = weight_3d_k_p_qk_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_qk_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_qk_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

if flag_calc;
%%%%%%%%;
% Use noiseless templates and synthetic viewing-angle distribution. ;
%%%%%%%%;
%test_slice_vs_volume_integral_helper_eig_imagecount_5;
%%%%%%%%;
% Use noiseless templates and synthetic viewing-angle distribution and principal-modes. ;
%%%%%%%%;
test_slice_vs_volume_integral_from_a_k_p_helper_eig_imagecount_7;
%%%%%%%%;
end;%if flag_calc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
