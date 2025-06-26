%%%%%%%%;
% Testing Hemaglu molecule from Pilar (20250614). ;
% updated (20250619) to use batch_cat_0, rather than the original star-file from particles_mismatch. ;
%%%%%%%%;

platform = 'rusty';%platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data/rangan'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home/rangan'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home/rangan'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home/rangan'; end;
if (strcmp(platform,'ceph')); setup_ceph; string_root = 'mnt/home/rangan/ceph'; end;

str_thisfunction = 'test_Hemaglu_1';
flag_verbose=1; flag_disp=1; nf=0;

k_int = 48;
k_eq_d_double = 1.00;
t_eq_d_double = 1.00;
n_w_int = 1;
f_rand = 0;
pm_n_UX_rank_max = +Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

global_parameter = struct('type','parameter');
global_parameter.flag_recalc = 0;
global_parameter.flag_replot = 0;
global_parameter.flag_center_volume_via_trans = 0;
global_parameter.flag_center_volume_via_value = 1;
global_parameter.flag_center_image_via_trans = 1;
global_parameter.flag_center_image_via_value = 1;
global_parameter.flag_invert_image_via_value = 1;
global_parameter.tolerance_master = 1e-2;
fname_prefix='Hemaglu_Smallset';
Pixel_Spacing = 1.03; %<-- Angstroms per pixel. ;
sigma_angstrom = 3.00; %<-- std of atomic gaussian in angstroms. ;
dir_noroot_volume='Hemaglu_Smallset/models';
fname_noroot_volume_={'6wxb_CA.pdb','6wxb_CA_def.pdb'};
dir_noroot_star='Hemaglu_Smallset/particles_mismatch';
fname_noroot_star='batch_cat_0.star';

%%%%%%%%;
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_pm = sprintf('/%s/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
dir_pm_txt = sprintf('%s_txt',dir_pm); if (~exist(dir_pm_txt,'dir')); disp(sprintf(' %% mkdir %s',dir_pm_txt)); mkdir(sprintf('%s',dir_pm_txt)); end;
dir_pm_mat = sprintf('%s_mat',dir_pm); if (~exist(dir_pm_mat,'dir')); disp(sprintf(' %% mkdir %s',dir_pm_mat)); mkdir(sprintf('%s',dir_pm_mat)); end;
dir_pm_jpg = sprintf('%s_jpg',dir_pm); if (~exist(dir_pm_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_pm_jpg)); mkdir(sprintf('%s',dir_pm_jpg)); end;
dir_pm_individual_jpg = sprintf('%s_individual_jpg',dir_pm); if (~exist(dir_pm_individual_jpg,'dir')); disp(sprintf(' %% mkdir %s',dir_pm_individual_jpg)); mkdir(sprintf('%s',dir_pm_individual_jpg)); end;
dir_star = sprintf('/%s/dir_cryoem/dir_%s',string_root,dir_noroot_star);
n_volume = numel(fname_noroot_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
flag_recalc = global_parameter.flag_recalc;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
flag_replot = global_parameter.flag_replot;
if (~isfield(global_parameter,'flag_center_volume_via_trans')); global_parameter.flag_center_volume_via_trans = 0; end; %<-- parameter_bookmark. ;
flag_center_volume_via_trans = global_parameter.flag_center_volume_via_trans;
if (~isfield(global_parameter,'flag_center_volume_via_value')); global_parameter.flag_center_volume_via_value = 0; end; %<-- parameter_bookmark. ;
flag_center_volume_via_value = global_parameter.flag_center_volume_via_value;
if (~isfield(global_parameter,'flag_center_image_via_trans')); global_parameter.flag_center_image_via_trans = 0; end; %<-- parameter_bookmark. ;
flag_center_image_via_trans = global_parameter.flag_center_image_via_trans;
if (~isfield(global_parameter,'flag_center_image_via_value')); global_parameter.flag_center_image_via_value = 0; end; %<-- parameter_bookmark. ;
flag_center_image_via_value = global_parameter.flag_center_image_via_value;
if (~isfield(global_parameter,'flag_invert_image_via_value')); global_parameter.flag_invert_image_via_value = 0; end; %<-- parameter_bookmark. ;
flag_invert_image_via_value = global_parameter.flag_invert_image_via_value;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
tolerance_master = global_parameter.tolerance_master;
nf=0;

%%%%%%%%;
% Load or create downsampled volume. ;
%%%%%%%%;
%%%%%%%%;
fname_mat = sprintf('%s/a_x_u_base_vxxx__.mat',dir_pm_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_x_u_file = 256; %<-- box diameter in pixels. ;
n_x_u_pack = 64;
n_file_pack = n_x_u_file/n_x_u_pack; assert(n_file_pack==round(n_file_pack));
a_x_u_pack_vxxx__ = cell(n_volume,1);
a_x_u_base_vxxx__ = cell(n_volume,1);
%%%%%%%%;
for nvolume=0:n_volume-1;
fname_a_pdb = sprintf('/%s/dir_cryoem/dir_%s/%s',string_root,dir_noroot_volume,fname_noroot_volume_{1+nvolume});
parameter = struct('type','parameter');
parameter.flag_verbose = flag_verbose;
parameter.flag_disp = 1; parameter.nf=nf; nf=nf+1;
parameter.flag_center = 1;
parameter.sigma_angstrom = sigma_angstrom;
[ ...
 parameter ...
,a_x_u_load_ ...
] = ...
a_x_u___from_pdb_2( ...
 parameter ...
,fname_a_pdb ...
,n_x_u_file ...
,Pixel_Spacing ...
,n_file_pack ...
);
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
a_x_u_pack_vxxx__{1+nvolume} = a_x_u_pack_;
clear a_x_u_load_;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
% Calculate moments. ;
%%%%%%%%;
for nvolume=0:n_volume-1;
%%%%%%%%;
a_x_u_pack_= a_x_u_pack_vxxx__{1+nvolume};
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
if flag_center_volume_via_trans;
disp(sprintf(' %% centering volume by translating in real-space'));
a_x_u_base_ = b_rho_x_u_pack_;
end;%if flag_center_volume_via_trans;
%%%%%%%%;
if flag_center_volume_via_value; a_x_u_base_ = a_x_u_base_ - mean(a_x_u_base_,'all'); end;
a_x_u_base_vxxx__{1+nvolume} = a_x_u_base_;
%%%%%%%%;
fname_fig = sprintf('%s/a_x_u_base_%d_',dir_pm_jpg,nvolume);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); isosurface_f_x_c_0(a_x_u_pack_,98.5); title('a packed');
subplot(2,2,2); isosurface_f_x_c_0(a_x_u_pack_,[97.5,98.5,99.5]); title('a packed');
subplot(2,2,3); isosurface_f_x_c_0(b_rho_x_u_pack_,98.5); title('b packed');
subplot(2,2,4); isosurface_f_x_c_0(b_rho_x_u_pack_,[97.5,98.5,99.5]); title('b packed');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','n_x_u_pack','n_pack','pack_row_ij_','pack_col_ij_','pack_val_ij_','x_u_pack_' ...
     ,'x_u_0_','x_u_1_','x_u_2_','x_u_0___','x_u_1___','x_u_2___','n_x_u','n_xxx_u','xxx_u_weight_','a_x_u_base_vxxx__' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% simple visualization of a_base. ;
%%%%%%%%
flag_plot=1;
if flag_plot;
for nvolume=0:n_volume-1;
a_x_u_base_ = a_x_u_base_vxxx__{1+nvolume};
prctile_ = [94.0:0.5:99.5]; n_prctile = numel(prctile_);
figure(1);clf;figbig;
p_row = 3; p_col = 4; ns=0;
for nprctile=0:n_prctile-1;
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_p = prctile_(1+nprctile);
isosurface_f_x_c_0(a_x_u_base_,tmp_p);
title(sprintf('a base %.1f',tmp_p));
end;%for nprctile=0:n_prctile-1;
fname_fig = sprintf('%s/a_x_u_base_%d_',dir_pm_jpg,nvolume);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
close(gcf);
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for nvolume=0:n_volume-1;
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
fname_mat = sprintf('%s/a_k_p_quad_vqk__.mat',dir_pm_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_k_p_quad_vqk__ = cell(n_volume,1);
a_x_u_reco_vxxx__ = cell(n_volume,1);
for nvolume=0:n_volume-1;
a_x_u_base_ = a_x_u_base_vxxx__{1+nvolume};
tmp_t = tic;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi);
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
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_base_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_x_u_base_(:)',a_x_u_base_(:),'a_x_u_reco_',a_x_u_reco_);
disp(sprintf(' %% at this point one should ensure that a_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
a_k_p_quad_vqk__{1+nvolume} = a_k_p_quad_;
a_x_u_reco_vxxx__{1+nvolume} = a_x_u_reco_;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
save(fname_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_' ...
     ,'k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_' ...
     ,'weight_3d_k_all_','weight_shell_k_' ...
     ,'n_k_p_r','k_p_r_' ...
     ,'weight_3d_k_p_r_' ...
     ,'k_c_0_all_','k_c_1_all_','k_c_2_all_' ...
     ,'J_node_','J_weight_','J_chebfun_','J_polyval_' ...
     ,'a_k_p_quad_vqk__' ...
     ,'a_x_u_reco_vxxx__' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s/a_k_p_quad_',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row=1;p_col=n_volume;
for nvolume=0:n_volume-1;
a_k_p_quad_ = a_k_p_quad_vqk__{1+nvolume};
subplot(p_row,p_col,1+nvolume);
plot(k_p_r_all_,log10(abs(a_k_p_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
end;%for nvolume=0:n_volume-1;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
fname_mat = sprintf('%s/a_k_Y_quad_vy__.mat',dir_pm_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
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
a_k_Y_quad_vy__ = cell(n_volume,1);
a_k_p_reco_vqk__ = cell(n_volume,1);
for nvolume=0:n_volume-1;
a_k_p_quad_ = a_k_p_quad_vqk__{1+nvolume};
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(0*flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(0*flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_p_quad_',a_k_p_quad_,'a_k_p_reco_',a_k_p_reco_,' %<-- should be <1e-2');
%%%%%%%%;
% if necessary can call: a_k_Y_quad__ = local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_quad_); 
%%%%%%%%;
a_k_Y_quad_vy__{1+nvolume} = a_k_Y_quad_;
a_k_p_reco_vqk__{1+nvolume} = a_k_p_reco_;
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
save(fname_mat ...
     ,'l_max_' ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'Y_l_val_','Y_m_val_','Y_k_val_','weight_Y_' ...
     ,'a_k_Y_quad_vy__' ...
     ,'a_k_p_reco_vqk__' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s/a_k_Y_quad_A',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row = n_volume; p_col = 3; np=0;
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_vy__{1+nvolume};
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(p_row,p_col,1+np);np=np+1;
plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
end;%for nvolume=0:n_volume-1;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/a_k_Y_quad_',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
p_row=1;p_col=n_volume;
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_vy__{1+nvolume};
subplot(p_row,p_col,1+nvolume);
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_)),[-10,0],colormap_beach());
xlabel('l_val','Interpreter','none');
ylabel('m_val','Interpreter','none');
zlabel('k_p_r','Interpreter','none');
axis vis3d; view([-65,+20]);
title(sprintf('a_k_Y_quad_%d_',nvolume),'Interpreter','none');
end;%for nvolume=0:n_volume-1;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% generate templates S_k_p_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
fname_mat = sprintf('%s/S_k_p_vwkS___.mat',dir_pm_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
S_k_p_vwkS___ = cell(n_volume,1);
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_vy__{1+nvolume};
n_w_max = 2*(l_max_max+1);
template_k_eq_d = -1;
viewing_k_eq_d = t_eq_d_double*1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_k_p__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
,template_k_c_0__ ...
,template_k_c_1__ ...
,template_k_c_2__ ...
] = ...
get_template_1( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max*ones(n_k_p_r,1) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (flag_verbose>1); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
S_k_p_vwkS___{1+nvolume} = S_k_p__;
end;%for nvolume=0:n_volume-1;
save(fname_mat ...
     ,'n_w_max','template_k_eq_d','viewing_k_eq_d' ...
     ,'S_k_p_vwkS___' ...
     ,'n_w_' ...
     ,'weight_2d_k_p_r_' ...
     ,'weight_2d_k_all_' ...
     ,'n_viewing_all' ...
     ,'viewing_azimu_b_all_' ...
     ,'viewing_polar_a_all_' ...
     ,'viewing_weight_all_' ...
     ,'n_viewing_polar_a' ...
     ,'viewing_polar_a_' ...
     ,'n_viewing_azimu_b_' ...
     ,'template_k_c_0__' ...
     ,'template_k_c_1__' ...
     ,'template_k_c_2__' ...
     ,'n_S','n_w_sum','n_w_csum_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
for nvolume=0:n_volume-1;
S_k_p__ = S_k_p_vwkS___{1+nvolume};
fname_fig = sprintf('%s/S_k_p_%d__',dir_pm_jpg,nvolume);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p__(:,1+nS)),Slim_,colormap_80s);
axis equal; axisnotick; title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('Sample S_k_p__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for nvolume=0:n_volume-1;
%%%%%%%%;
flag_check=1;
if flag_check;
%%%%%%%%;
% test pm_template_2 vs get_template_1. ;
%%%%%%%%;
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_vy__{1+nvolume};
S_k_p__ = S_k_p_vwkS___{1+nvolume};
tmp_t = tic();
[ ...
 tmp_S_k_p__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(local_yk__from_yk_(n_k_p_r,l_max_,a_k_Y_quad_),[n_lm_max,n_k_p_r]) ...
,t_eq_d_double/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
tmp_S_k_p__ = reshape(tmp_S_k_p__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
fnorm_disp(flag_verbose,'S_k_p__',S_k_p__,'tmp_S_k_p__',tmp_S_k_p__,' %<-- should be <1e-12');
end;%for nvolume=0:n_volume-1;
end;%if flag_check;

%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
fname_mat = sprintf('%s/M_k_p__.mat',dir_pm_mat);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_M = 1024;
[ ...
 M_x_c___ ...
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
 dir_star ...
,fname_noroot_star ...
,n_M ...
);
n_M = size(M_x_c___,3);
if (fnorm(Voltage_CTF_)< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_ = 300*ones(n_M,1); end;
if flag_invert_image_via_value; M_x_c___ = -M_x_c___; end;
%%%%%%%%;
% Remove any edge artefacts, mean center and normalize each image. ;
%%%%%%%%;
disp(sprintf(' %% Removing edge-artefacts'));
n_M_ext_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[M_x_c___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c___(:,:,1+nM),4,0.5,2,0);
end;%for nM=0:n_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
n_x_M_u = size(M_x_c___,1);
assert(n_x_M_u==size(M_x_c___,2));
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
M_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nM))); %<-- no mask. ;
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nM) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nM) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
M_mask_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nM)).*x_c_mask__); %<-- radial mask. ;
M_mask_abs_avg = mean(M_mask_abs_x_c_,'all');
M_mask_abs_x_c_0_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_0__,'all');
M_mask_abs_x_c_1_avg = mean(M_mask_abs_x_c_/M_mask_abs_avg.*x_c_1__,'all');
M_mask_abs_x_c_0_avg_(1+nM) = M_mask_abs_x_c_0_avg;
M_mask_abs_x_c_1_avg_(1+nM) = M_mask_abs_x_c_1_avg;
clear M_mask_abs_x_c_;
end;%for nM=0:n_M-1;
%%%%%%%%;
fname_fig = sprintf('%s/M_abs_x_c_avg_',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now convert images to M_k_p__. ;
%%%%%%%%;
dx = diameter_x_c/n_x_M_u;
M_k_p__ = zeros(n_w_sum,n_M);
N_k_p__ = zeros(n_w_sum,n_M);
image_center_delta_x_c_0_ = zeros(n_M,1);
image_center_delta_x_c_1_ = zeros(n_M,1);
n_x_M_center = max(n_w_max,n_x_M_u/4);
if (flag_verbose>1); disp(sprintf(' %% n_x_M_center %d',n_x_M_center)); end;
M_x_c_avg__ = zeros(n_x_M_center,n_x_M_center);
M_x_c_std__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_avg__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_std__ = zeros(n_x_M_center,n_x_M_center);
for nM=0:n_M-1;
if (mod(nM,128)==0); if (flag_verbose>1); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
if flag_center_image_via_value; M_x_c_ = M_x_c_ - mean(M_x_c_,'all'); end;
M_k_p_ = interp_x_c_to_k_p_xxnufft(n_x_M_u,diameter_x_c,n_x_M_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_M_u^2)*dx^2 ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
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
,M_k_p_ ...
,weight_2d_k_all_ ...
);
M_k_p__(:,1+nM) = M_k_p_;
N_k_p__(:,1+nM) = N_k_p_;
image_center_delta_x_c_0_(1+nM) = tmp_delta_x_c_0;
image_center_delta_x_c_1_(1+nM) = tmp_delta_x_c_1;
M_x_c_avg__ = M_x_c_avg__ + reshape(real(M_x_c_0in_),[n_x_M_center,n_x_M_center]);
M_x_c_std__ = M_x_c_std__ + reshape(real(M_x_c_0in_).^2,[n_x_M_center,n_x_M_center]);
N_x_c_avg__ = N_x_c_avg__ + reshape(real(M_x_c_out_),[n_x_M_center,n_x_M_center]);
N_x_c_std__ = N_x_c_std__ + reshape(real(M_x_c_out_).^2,[n_x_M_center,n_x_M_center]);
end;%for nM=0:n_M-1;
M_x_c_avg__ = M_x_c_avg__/n_M; M_x_c_std__ = M_x_c_std__/n_M; M_x_c_std__ = sqrt(M_x_c_std__ - M_x_c_avg__.^2);
N_x_c_avg__ = N_x_c_avg__/n_M; N_x_c_std__ = N_x_c_std__/n_M; N_x_c_std__ = sqrt(N_x_c_std__ - N_x_c_avg__.^2);
%%%%%%%%;
save(fname_mat ...
     ,'n_M','n_x_M_u','Pixel_Spacing','x_c_0_','x_c_1_','x_c_0__','x_c_1__' ...
     ,'n_M_ext_' ...
     ,'M_abs_x_c_0_avg_','M_abs_x_c_1_avg_' ...
     ,'M_mask_abs_x_c_0_avg_','M_mask_abs_x_c_1_avg_' ...
     ,'image_center_delta_x_c_0_','image_center_delta_x_c_1_' ...
     ,'M_k_p__','N_k_p__' ...
     ,'M_x_c_avg__','M_x_c_std__' ...
     ,'N_x_c_avg__','N_x_c_std__' ...
     ,'index_nCTF_from_nM_' ...
     ,'index_nM_from_nCTF_' ...
     ,'Voltage_CTF_' ...
     ,'DefocusU_CTF_' ...
     ,'DefocusV_CTF_' ...
     ,'DefocusAngle_CTF_' ...
     ,'SphericalAberration_CTF_' ...
     ,'AmplitudeContrast_CTF_' ...
     );
clear M_x_c___;
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s/M_all_center_FIGB',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figbig;figbeach();
subplot(1,2,1); imagesc(M_x_c_std__); axis image; axisnotick;
title('M_x_c_std__','Interpreter','none');
subplot(1,2,2); imagesc(N_x_c_std__); axis image; axisnotick;
title('N_x_c_std__','Interpreter','none');
sgtitle(sprintf(' %% centering filtered images'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now convert images to M_k_q__. ;
%%%%%%%%;
tmp_t = tic();
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q__ time %0.2fs',tmp_t)); end;
tmp_t = tic();
N_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
N_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,N_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
fname_fig = sprintf('%s/M_k_q__',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Mlim_ = std(abs(M_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nM = max(0,min(n_M-1,floor(n_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(M_k_q__(:,1+nM)),Mlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(M_k_q__(:,1+nM))/max(abs(M_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nM %d',nM));
end;%for nl=0:15-1;
sgtitle(sprintf('M_k_q__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/N_k_q__',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Nlim_ = std(abs(N_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nM = max(0,min(n_M-1,floor(n_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(N_k_q__(:,1+nM)),Nlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(N_k_q__(:,1+nM))/max(abs(N_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nM %d',nM));
end;%for nl=0:15-1;
sgtitle(sprintf('N_k_q__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s/CTF_k_p_wkC__.mat',dir_pm_mat);
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
fname_fig = sprintf('%s/CTF_k_p_sample__',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s/CTF_k_p_xcor__',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
figbeach();
subplot(1,2,1); imagesc(CTF_k_p_r_xavg__); axis image; axisnotick; title('CTF_k_p_r_xavg__','Interpreter','none');
subplot(1,2,2); imagesc(CTF_k_p_r_xcor__); axis image; axisnotick; title('CTF_k_p_r_xcor__','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s/CTF_k_p_r_xxxx__',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now cluster the CTF based on tolerance_cluster. ;
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
fname_fig = sprintf('%s/VSCTF_Mc_scatter_FIGA',dir_pm_jpg);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;figsml;figbeach();
markersize_use = 32;
scatter(VSCTF_Mc__(:,1+0),VSCTF_Mc__(:,1+min(1,n_CTF_rank-1)),markersize_use,index_ncluster_from_nM_,'filled');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now use ideal principal-modes to align images to volume. ;
%%%%%%%%;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
for flag_N_vs_M = 0:1;
if flag_N_vs_M==0; tmp_N_k_p__ = M_k_p__; tmp_str_N_vs_M = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_N_k_p__ = N_k_p__; tmp_str_N_vs_M = 'N'; end;%if flag_N_vs_M==1;
for nvolume=0:n_volume-1;
a_k_Y_quad_ = a_k_Y_quad_vy__{1+nvolume};
fname_mat = sprintf('%s/a_k_Y_reco_%d_from_%s__.mat',dir_pm_mat,nvolume,tmp_str_N_vs_M);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.tolerance_pm = 1e-6;
parameter.t_eq_d_double = t_eq_d_double;
parameter.n_iteration = 8;
parameter.delta_r_max = tmp_delta_r_max;
parameter.n_delta_v_requested = 24;
parameter.delta_r_upb = tmp_delta_r_upb;
parameter.svd_eps = 1e-6;
parameter.qbp_eps = 1e-2;
parameter.flag_MS_vs_SM = 0;
parameter.f_rand = f_rand;
parameter.fname_align_a_k_Y_pre = sprintf('%s/a_k_Y_reco_%d_from_%s__',dir_pm_mat,nvolume,tmp_str_N_vs_M);
[ ...
 parameter ...
,a_k_Y_reco_yki__ ...
,corr_a_k_Y_i_ ...
,euler_polar_a_Mi__ ...
,euler_azimu_b_Mi__ ...
,euler_gamma_z_Mi__ ...
,image_delta_x_acc_Mi__ ...
,image_delta_y_acc_Mi__ ...
,image_delta_x_upd_Mi__ ...
,image_delta_y_upd_Mi__ ...
,flag_image_delta_upd_Mi__ ...
,image_I_value_Mi__ ...
,image_X_value_Mi__ ...
,image_S_index_Mi__ ...
,tmp_n_S ...
,tmp_template_viewing_azimu_b_all_ ...
,tmp_template_viewing_polar_a_all_ ...
,X_SMi___ ...
,delta_x_SMi___ ...
,delta_y_SMi___ ...
,gamma_z_SMi___ ...
,I_value_SMi___ ...
,CTF_UX_S_l2_SMi___ ...
,UX_M_l2_SMi___ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,n_M ...
,tmp_N_k_p__ ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
assert(tmp_n_S==n_S);
assert(fnorm(tmp_template_viewing_azimu_b_all_-viewing_azimu_b_all_)<1e-6);
assert(fnorm(tmp_template_viewing_polar_a_all_-viewing_polar_a_all_)<1e-6);
if (~exist(fname_mat,'file'));
save(fname_mat ...
     ,'parameter' ...
     ,'a_k_Y_reco_yki__' ...
     ,'corr_a_k_Y_i_' ...
     ,'euler_polar_a_Mi__' ...
     ,'euler_azimu_b_Mi__' ...
     ,'euler_gamma_z_Mi__' ...
     ,'image_delta_x_acc_Mi__' ...
     ,'image_delta_y_acc_Mi__' ...
     ,'image_delta_x_upd_Mi__' ...
     ,'image_delta_y_upd_Mi__' ...
     ,'flag_image_delta_upd_Mi__' ...
     ,'image_I_value_Mi__' ...
     ,'image_X_value_Mi__' ...
     ,'image_S_index_Mi__' ...
     ,'X_SMi___' ...
     ,'delta_x_SMi___' ...
     ,'delta_y_SMi___' ...
     ,'gamma_z_SMi___' ...
     ,'I_value_SMi___' ...
     ,'CTF_UX_S_l2_SMi___' ...
     ,'UX_M_l2_SMi___' ...
     );
end;%if (~exist(fname_mat,'file'));
end;%if (~exist(fname_mat,'file'));
end;%for nvolume=0:n_volume-1;
clear tmp_N_k_p__;
end;%for flag_N_vs_M = 0:1;
%%%%%%%%;

%%%%%%%%;
euler_polar_a_rfit_Mv__ = zeros(n_M,n_volume);
euler_azimu_b_rfit_Mv__ = zeros(n_M,n_volume);
euler_gamma_z_rfit_Mv__ = zeros(n_M,n_volume);
image_delta_x_rfit_Mv__ = zeros(n_M,n_volume);
image_delta_y_rfit_Mv__ = zeros(n_M,n_volume);
corr_a_k_Y_rfit_v_ = zeros(n_volume,1);
image_X_value_rfit_Mv__ = zeros(n_M,n_volume);
image_X_SMv___ = zeros(n_S,n_M,n_volume);
a_k_Y_reco_ykv__ = zeros(n_lm_sum,n_volume);
flag_N_vs_M = flag_center_image_via_trans;
if flag_N_vs_M==0; tmp_N_k_p__ = M_k_p__; tmp_str_N_vs_M = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_N_k_p__ = N_k_p__; tmp_str_N_vs_M = 'N'; end;%if flag_N_vs_M==1;
for nvolume=0:n_volume-1;
fname_mat = sprintf('%s/a_k_Y_reco_%d_from_%s__.mat',dir_pm_mat,nvolume,tmp_str_N_vs_M);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
if isfield(tmp_,'parameter'); f_rand = tmp_.parameter.f_rand; end;
euler_polar_a_rfit_Mv__(:,1+nvolume) = tmp_.euler_polar_a_Mi__(:,end);
euler_azimu_b_rfit_Mv__(:,1+nvolume) = tmp_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_rfit_Mv__(:,1+nvolume) = tmp_.euler_gamma_z_Mi__(:,end);
image_delta_x_rfit_Mv__(:,1+nvolume) = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_rfit_Mv__(:,1+nvolume) = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_rfit_v_(1+nvolume) = tmp_.corr_a_k_Y_i_(end);
image_X_value_rfit_Mv__(:,1+nvolume) = tmp_.image_X_value_Mi__(:,end-1);
image_X_SMv___(:,:,1+nvolume) = tmp_.X_SMi___(:,:,end-1);
CTF_UX_S_l2_SMv___(:,:,1+nvolume) = tmp_.CTF_UX_S_l2_SMi___(:,:,end-1);
UX_M_l2_SMv___(:,:,1+nvolume) = tmp_.UX_M_l2_SMi___(:,:,end-1);
a_k_Y_reco_ykv__(:,1+nvolume) = tmp_.a_k_Y_reco_yki__(:,end);
clear tmp_;
end;%if ( exist(fname_mat,'file'));
end;%for nvolume=0:n_volume-1;
%%%%%%%%;

if flag_disp;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 2; np=0;
for ntype=0:4-1;
if ntype==0; nvolume=0; tmp_a_k_yk_ =    a_k_Y_quad_vy__{1+nvolume}; tmp_str_title = sprintf('%s quad',fname_noroot_volume_{1+nvolume}); end;
if ntype==1; nvolume=0; tmp_a_k_yk_ = a_k_Y_reco_ykv__(:,1+nvolume); tmp_str_title = sprintf('%s reco',fname_noroot_volume_{1+nvolume}); end;
if ntype==2; nvolume=1; tmp_a_k_yk_ =    a_k_Y_quad_vy__{1+nvolume}; tmp_str_title = sprintf('%s quad',fname_noroot_volume_{1+nvolume}); end;
if ntype==3; nvolume=1; tmp_a_k_yk_ = a_k_Y_reco_ykv__(:,1+nvolume); tmp_str_title = sprintf('%s reco',fname_noroot_volume_{1+nvolume}); end;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
[ ... 
  tmp_a_x_u_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,tmp_a_k_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
subplot(p_row,p_col,1+np);np=np+1;cla;
isosurface_f_x_u_1([],tmp_a_x_u_); title(tmp_str_title,'Interpreter','none');
drawnow();
end;%for ntype=0:4-1;
%%%%;
fname_fig_pre = sprintf('%s/a_x_u_reco_x_from_%s_FIGA',dir_pm_jpg,tmp_str_N_vs_M);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%close(gcf);
%%%%%%%%;
end;%if flag_disp;


%%%%%%%%;
% Create diagnostic figures for each image. ;
%%%%%%%%;
if flag_N_vs_M==0; tmp_N_k_p__ = M_k_p__; tmp_str_N_vs_M = 'M'; end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1; tmp_N_k_p__ = N_k_p__; tmp_str_N_vs_M = 'N'; end;%if flag_N_vs_M==1;
gamma_z_ = transpose(linspace(0,2*pi,n_w_max+1)); gamma_z_ = gamma_z_(1:n_w_max);
delta_max = 0.1; n_delta_x = 1+8; n_delta_y = 1+8;
tmp_delta_x_ = transpose(linspace(-delta_max,+delta_max,n_delta_x));
tmp_delta_y_ = transpose(linspace(-delta_max,+delta_max,n_delta_y));
nvolume = 0;
[~,tmp_ij_sort_] = sort(image_X_value_rfit_Mv__(:,1+nvolume),'descend'); tmp_index_sort_ = tmp_ij_sort_-1;
%[~,tmp_ij_sort_] = sort(abs(image_delta_y_rfit_Mv__(:,1+nvolume)),'descend'); tmp_index_sort_ = tmp_ij_sort_-1;
%[~,tmp_ij_sort_] = sort(abs(image_delta_x_rfit_Mv__(:,1+nvolume)) + abs(image_delta_y_rfit_Mv__(:,1+nvolume)),'ascend'); tmp_index_sort_ = tmp_ij_sort_-1;
for nM_srt = 0:n_M-1;
nM = tmp_index_sort_(1+nM_srt);
euler_polar_a_rfit = euler_polar_a_rfit_Mv__(1+nM,1+nvolume);
euler_azimu_b_rfit = euler_azimu_b_rfit_Mv__(1+nM,1+nvolume);
euler_gamma_z_rfit = euler_gamma_z_rfit_Mv__(1+nM,1+nvolume);
image_delta_x_rfit = image_delta_x_rfit_Mv__(1+nM,1+nvolume);
image_delta_y_rfit = image_delta_y_rfit_Mv__(1+nM,1+nvolume);
image_X_value_rfit = image_X_value_rfit_Mv__(1+nM,1+nvolume);
image_index_S_rfit = efind( (viewing_polar_a_all_==euler_polar_a_rfit) & (viewing_azimu_b_all_==euler_azimu_b_rfit) );
assert(numel(image_index_S_rfit)==1);
nS = image_index_S_rfit;
%assert(image_X_value_rfit==image_X_SMv___(1+nS,1+nM,1+nvolume)); %<-- note, these will not match when f_rand==0.05 in pm_align_M_k_p_to_a_k_Y_3. ;
tmp_image_X_S_sort_ = sort(image_X_SMv___(:,1+nM,1+nvolume),'descend');
tmp_index_S_sort = min(efind(tmp_image_X_S_sort_==image_X_value_rfit));
assert(tmp_index_S_sort<=ceil(f_rand*n_S));
fname_fig_pre = sprintf('%s/diagnostic_nM_%.3d_FIGA',dir_pm_individual_jpg,nM);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if  exist(fname_fig_jpg,'file'); if (flag_verbose>0); disp(sprintf(' %% %s found, not creating',fname_fig_jpg)); end; end;
if flag_replot | ~exist(fname_fig_jpg,'file');
if (flag_verbose>0); disp(sprintf(' %% %s not found, creating',fname_fig_jpg)); end;
S_k_p_wk_ = S_k_p_vwkS___{1+nvolume}(:,1+nS);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_S_k_p_wk_ = CTF_k_p_wk_.*S_k_p_wk_;
R_CTF_S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_,+euler_gamma_z_rfit);
R_CTF_S_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,R_CTF_S_k_p_wk_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum );
tmp_N_k_p_wk_ = tmp_N_k_p__(:,1+nM);
tmp_T_N_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_N_k_p_wk_,+image_delta_x_rfit,+image_delta_y_rfit);
tmp_T_N_x_c_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum );
tmp_SS = real(sum(conj(R_CTF_S_k_p_wk_).*R_CTF_S_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
tmp_NN = real(sum(conj(tmp_T_N_k_p_wk_).*tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
tmp_SN = real(sum(conj(R_CTF_S_k_p_wk_).*tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
image_X_value_calc = tmp_SN/max(1e-12,sqrt(tmp_SS*tmp_NN));
%%%%;
image_X_value_calc_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
tmp_gamma_z = gamma_z_(1+nw);
tmp_R_CTF_S_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,CTF_S_k_p_wk_,+tmp_gamma_z);
tmp_SS = real(sum(conj(tmp_R_CTF_S_k_p_wk_).*tmp_R_CTF_S_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
tmp_SN = real(sum(conj(tmp_R_CTF_S_k_p_wk_).*tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
image_X_value_calc_w_(1+nw) = tmp_SN/max(1e-12,sqrt(tmp_SS*tmp_NN));
end;%for nw=0:n_w_max-1;
%%%%;
image_X_value_calc_d__ = zeros(n_delta_x,n_delta_y);
for ndelta_x=0:n_delta_x-1;for ndelta_y=0:n_delta_y-1;
tmp_delta_x=tmp_delta_x_(1+ndelta_x); tmp_delta_y=tmp_delta_y_(1+ndelta_y);
tmp_tmp_T_N_k_p_wk_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_T_N_k_p_wk_,+tmp_delta_x,+tmp_delta_y);
tmp_NN = real(sum(conj(tmp_tmp_T_N_k_p_wk_).*tmp_tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
tmp_SN = real(sum(conj(R_CTF_S_k_p_wk_).*tmp_tmp_T_N_k_p_wk_.*weight_2d_k_all_*(2*pi)^2));
image_X_value_calc_d__(1+ndelta_x,1+ndelta_y) = tmp_SN/max(1e-12,sqrt(tmp_SS*tmp_NN));
end;end;%for ndelta_x=0:n_delta_x-1;for ndelta_y=0:n_delta_y-1;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
%figure(1);clf;figmed;
linewidth_use = 3;
markersize_use = 12;
p_row=2; p_col=2; np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,R_CTF_S_x_c_xx_,[],colormap_beach); axis image; axisnotick;
title('R_CTF_S_x_c_xx_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,tmp_T_N_x_c_xx_,[],colormap_beach); axis image; axisnotick;
title('tmp_T_N_x_c_xx_','Interpreter','none');
sgtitle(sprintf('nM_srt %.3d/%.3d nM %.3d/%.3d X %+0.2f (%+0.2f) gamma/pi %0.2f [dx0,dx1] [%+0.2f,%+0.2f]',nM_srt,n_M,nM,n_M,image_X_value_rfit,image_X_value_calc,euler_gamma_z_rfit/pi,image_delta_x_rfit,image_delta_y_rfit),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
hold on;
plot(gamma_z_/pi,image_X_value_calc_w_,'r-','LineWidth',linewidth_use);
plot(euler_gamma_z_rfit/pi,image_X_value_calc,'ro','MarkerSize',markersize_use)
plot(euler_gamma_z_rfit/pi,image_X_value_rfit,'go','MarkerSize',markersize_use)
xlim([0,2]); ylim([-1,+1]); grid on; xlabel('gamma_z','Interpreter','none'); ylabel('X'); title('image_X_value_calc_w_','Interpreter','none');
tmp_subplot = subplot(p_row,p_col,1+np);np=np+1;
imagesc(image_X_value_calc_d__,[-1,+1]); colorbar; axis image; axisnotick;
title('image_X_value_calc_d__','Interpreter','none');
xlabel('delta_x','Interpreter','none');
xlabel('delta_y','Interpreter','none');
colormap(tmp_subplot,colormap_81s);
drawnow();
%%%%;
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
%%%%%%%%;
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
end;%for nM_srt = 0:n_M-1;

%%%%%%%%;
% Now calculate the likelihoods. ;
% Note that here we use image_X_SMv___, rather than image_X_value_rfit_Mv__ (the latter taken from the top f_rand templates). ;
% plot(reshape(image_X_SMv___(:,:,1+0),[n_S*n_M,1]),reshape(image_X_SMv___(:,:,1+1),[n_S*n_M,1]),'.',[0,1],[0,1],'k-');
% plot(reshape(max(image_X_SMv___(:,:,1+0),[],1),[n_M,1]),reshape(max(image_X_SMv___(:,:,1+1),[],1),[n_M,1]),'.',[0,1],[0,1],'k-');
% plot(image_X_value_rfit_Mv__(:,1+0),image_X_value_rfit_Mv__(:,1+1),'.',[0,1],[0,1],'k-');
%%%%%%%%;
I_costheta_SMv___ = real(image_X_SMv___); %<-- similarity measure. ;
I_sinthetasquared_SMv___ = 1 - I_costheta_SMv___.^2; %<-- dis-similarity measure. ;
I_Tcen_Tcen_SMv___ = CTF_UX_S_l2_SMv___;
I_Mcen_Mcen_SMv___ = UX_M_l2_SMv___;
I_E_E = n_x_u_pack*n_x_u_pack; %<-- assuming that E_x_c_xx__ = ones(n_x_u_pack,n_x_u_pack). ;
n_pixel = I_E_E; %<-- It is not unreasonable to expect that the number of relevant pixels is given by the l2-norm of the mask. ;
%%%%%%%%;

%%%%%%%%;
% [TAG_NEGATIVE_LOG_LIKELIHOOD_MU_COMPLEX] ;
%%%%%%%%;
ssnll_SMv___ = 0.5*I_Tcen_Tcen_SMv___.*I_sinthetasquared_SMv___;
negative_log_likelihood_SMv___ = ...
-(3/2 - n_pixel/2) * log(2*pi) ...
+(1/2) * log(max(1e-12,I_Mcen_Mcen_SMv___)) ...
+(1.0) * log(2) ...
-(2 - n_pixel/2) * log(max(1e-12,ssnll_SMv___)) ...
- gammaln(n_pixel/2 - 2) ...
+ log(max(1e-12,I_E_E)) ...
;
negative_log_likelihood_optimal_Mv__ = reshape(min(negative_log_likelihood_SMv___,[],1),[n_M,n_volume]);
negative_log_likelihood_average_Mv__ = ...
negative_log_likelihood_optimal_Mv__ ...
-log( ...
      reshape( ...
	       sum( ...
		    bsxfun(@times ...
			   ,exp( ...
				 -bsxfun(@minus ...
					 ,negative_log_likelihood_SMv___ ...
					 ,reshape(negative_log_likelihood_optimal_Mv__,[1,n_M,n_volume]) ...
					 ) ...
				 ) ...
			   ,reshape(viewing_weight_all_,[n_S,1,1])/(4*pi*k_p_r_max^2) ...
			   ) ...
		    ,1 ...
		    ) ...
	       ,[n_M,n_volume] ...
	       ) ...
      ) ...
;
%%%%%%%%;
fname_fig_pre = sprintf('%s/nll_FIGA',dir_pm_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
fontsize_use = 14;
subplot(1,1,1);
hold on;
plot(negative_log_likelihood_optimal_Mv__(:,1+0),negative_log_likelihood_optimal_Mv__(:,1+1),'ko','MarkerSize',markersize_use);
plot(negative_log_likelihood_average_Mv__(:,1+0),negative_log_likelihood_average_Mv__(:,1+1),'rx','MarkerSize',markersize_use);
hold off;
xlabel(sprintf('nll %s',fname_noroot_volume_{1+0}),'Interpreter','none');
ylabel(sprintf('nll %s',fname_noroot_volume_{1+1}),'Interpreter','none');
legend({'optimal','average'},'Location','NorthWest');
title('negative-log-likelihood');
set(gca,'FontSize',fontsize_use);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;

%%%%%%%%;
% [TAG_NEGATIVE_LOG_LIKELIHOOD_MU_REAL] ;
%%%%%%%%;
ssnll_SMv___ = 0.5*I_Tcen_Tcen_SMv___.*I_sinthetasquared_SMv___;
negative_log_likelihood_SMv___ = ...
-(2/2 - n_pixel/2) * log(2*pi) ...
+(1/2) * log(max(1e-12,I_Mcen_Mcen_SMv___)) ...
+(1.0) * log(2) ...
-(3/2 - n_pixel/2) * log(max(1e-12,ssnll_SMv___)) ...
- gammaln(n_pixel/2 - 3/2) ...
+ log(max(1e-12,sqrt(I_E_E))) ...
;
negative_log_likelihood_optimal_Mv__ = reshape(min(negative_log_likelihood_SMv___,[],1),[n_M,n_volume]);
negative_log_likelihood_average_Mv__ = ...
negative_log_likelihood_optimal_Mv__ ...
-log( ...
      reshape( ...
	       sum( ...
		    bsxfun(@times ...
			   ,exp( ...
				 -bsxfun(@minus ...
					 ,negative_log_likelihood_SMv___ ...
					 ,reshape(negative_log_likelihood_optimal_Mv__,[1,n_M,n_volume]) ...
					 ) ...
				 ) ...
			   ,reshape(viewing_weight_all_,[n_S,1,1])/(4*pi*k_p_r_max^2) ...
			   ) ...
		    ,1 ...
		    ) ...
	       ,[n_M,n_volume] ...
	       ) ...
      ) ...
;
nll_opt_Mv__ = negative_log_likelihood_optimal_Mv__;
nll_avg_Mv__ = negative_log_likelihood_average_Mv__;
%%%%%%%%;
fname_fig_pre = sprintf('%s/nll_FIGB',dir_pm_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;clf;figsml;
markersize_use = 8;
fontsize_use = 14;
subplot(1,1,1);
hold on;
plot(negative_log_likelihood_optimal_Mv__(:,1+0),negative_log_likelihood_optimal_Mv__(:,1+1),'ko','MarkerSize',markersize_use);
plot(negative_log_likelihood_average_Mv__(:,1+0),negative_log_likelihood_average_Mv__(:,1+1),'rx','MarkerSize',markersize_use);
hold off;
xlabel(sprintf('nll %s',fname_noroot_volume_{1+0}),'Interpreter','none');
ylabel(sprintf('nll %s',fname_noroot_volume_{1+1}),'Interpreter','none');
legend({'optimal','average'},'Location','NorthWest');
title('negative-log-likelihood');
set(gca,'FontSize',fontsize_use);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);
%%%%%%%%;

nvolume=0;
T_6wxb_CA_ori_nvolume_ = nvolume*ones(n_M,1);
T_6wxb_CA_ori_nM_ = transpose([0:n_M-1]);
T_6wxb_CA_ori_nll_opt_ = nll_opt_Mv__(:,1+nvolume);
T_6wxb_CA_ori_nll_avg_ = nll_avg_Mv__(:,1+nvolume);
T_6wxb_CA_ori_X_ = image_X_value_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_euler_polar_a_rfit_ = euler_polar_a_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_euler_azimu_b_rfit_ = euler_azimu_b_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_euler_gamma_z_rfit_ = euler_gamma_z_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_image_delta_x_rfit_ = image_delta_x_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_image_delta_y_rfit_ = image_delta_y_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_ori_ = table( ...
 T_6wxb_CA_ori_nvolume_ ...
,T_6wxb_CA_ori_nM_ ...
,T_6wxb_CA_ori_nll_opt_ ...
,T_6wxb_CA_ori_nll_avg_ ...
,T_6wxb_CA_ori_X_ ...
,T_6wxb_CA_ori_euler_polar_a_rfit_ ...
,T_6wxb_CA_ori_euler_azimu_b_rfit_ ...
,T_6wxb_CA_ori_euler_gamma_z_rfit_ ...
,T_6wxb_CA_ori_image_delta_x_rfit_ ...
,T_6wxb_CA_ori_image_delta_y_rfit_ ...
);
fname_csv = sprintf('%s/T_6wxb_CA_ori_.csv',dir_pm_txt);
if flag_recalc | ~exist(fname_csv,'file');
disp(sprintf(' %% %s not found, creating',fname_csv));
writetable(T_6wxb_CA_ori_,fname_csv);
end;%if flag_recalc | ~exist(fname_csv,'file');
fname_tsv = sprintf('%s/Ttab_6wxb_CA_ori_.csv',dir_pm_txt);
if flag_recalc | ~exist(fname_tsv,'file');
disp(sprintf(' %% %s not found, creating',fname_tsv));
writetable(T_6wxb_CA_ori_,fname_tsv,'Delimiter','\t');
end;%if flag_recalc | ~exist(fname_tsv,'file');

nvolume=1;
T_6wxb_CA_def_nvolume_ = nvolume*ones(n_M,1);
T_6wxb_CA_def_nM_ = transpose([0:n_M-1]);
T_6wxb_CA_def_nll_opt_ = nll_opt_Mv__(:,1+nvolume);
T_6wxb_CA_def_nll_avg_ = nll_avg_Mv__(:,1+nvolume);
T_6wxb_CA_def_X_ = image_X_value_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_euler_polar_a_rfit_ = euler_polar_a_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_euler_azimu_b_rfit_ = euler_azimu_b_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_euler_gamma_z_rfit_ = euler_gamma_z_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_image_delta_x_rfit_ = image_delta_x_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_image_delta_y_rfit_ = image_delta_y_rfit_Mv__(:,1+nvolume);
T_6wxb_CA_def_ = table( ...
 T_6wxb_CA_def_nvolume_ ...
,T_6wxb_CA_def_nM_ ...
,T_6wxb_CA_def_nll_opt_ ...
,T_6wxb_CA_def_nll_avg_ ...
,T_6wxb_CA_def_X_ ...
,T_6wxb_CA_def_euler_polar_a_rfit_ ...
,T_6wxb_CA_def_euler_azimu_b_rfit_ ...
,T_6wxb_CA_def_euler_gamma_z_rfit_ ...
,T_6wxb_CA_def_image_delta_x_rfit_ ...
,T_6wxb_CA_def_image_delta_y_rfit_ ...
);
fname_csv = sprintf('%s/T_6wxb_CA_def_.csv',dir_pm_txt);
if flag_recalc | ~exist(fname_csv,'file');
disp(sprintf(' %% %s not found, creating',fname_csv));
writetable(T_6wxb_CA_def_,fname_csv);
end;%if flag_recalc | ~exist(fname_csv,'file');
fname_tsv = sprintf('%s/Ttab_6wxb_CA_def_.csv',dir_pm_txt);
if flag_recalc | ~exist(fname_tsv,'file');
disp(sprintf(' %% %s not found, creating',fname_tsv));
writetable(T_6wxb_CA_def_,fname_tsv,'Delimiter','\t');
end;%if flag_recalc | ~exist(fname_tsv,'file');

dir_tmp = '/data/rangan/dir_cryoem/dir_Hemaglu_Smallset';
fname_tmp_txt = sprintf('%s/Adi/ll-bioem-original',dir_tmp);
ll_bioem_M_ = textread(fname_tmp_txt);
fname_tmp_txt = sprintf('%s/Adi/ll-cryoLike-int-original',dir_tmp);
ll_cryol_M_ = textread(fname_tmp_txt);
figure(1+nf);nf=nf+1;clf;figmed;
markersize_use = 8;
fontsize_use = 12;
p_row = 1; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1; cla
plot(ll_bioem_M_,ll_cryol_M_,'ko','MarkerFaceColor','r','MarkerSize',markersize_use);
xlabel('bioem');ylabel('cryol'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla
plot(ll_bioem_M_,-T_6wxb_CA_ori_nll_opt_,'ko','MarkerFaceColor','r','MarkerSize',markersize_use);
xlabel('bioem');ylabel('matlab'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1; cla
plot(ll_cryol_M_,-T_6wxb_CA_ori_nll_opt_,'ko','MarkerFaceColor','r','MarkerSize',markersize_use);
xlabel('cryol');ylabel('matlab'); set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/scatterplot_ll_FIGA',dir_pm_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
close(gcf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning');return;

