function ...
[ ...
 global_parameter ...
] = ...
test_pm_riby_2( ...
 global_parameter ...
,fname_prefix ...
,fname_infix ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_micrograph ...
);

if (nargin<1);
table_data__ = { ...
'ribosome_yeast_Micrograph' , 'riby' , 1.06 , 'parsed_6Q8Y_whole_LSU_match3.pdb' , '147_Mar12_12.21.27_159_0.mrc' ; ...
};
n_experiment = size(table_data__,1);
%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
for nexperiment=8;%for nexperiment=(randperm(n_experiment)-1);
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
fname_infix = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_micrograph = table_data__{1+nexperiment,1+na}; na=na+1;
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,fname_infix,Pixel_Spacing,fname_nopath_volume,fname_nopath_micrograph ...
	     ));
global_parameter = struct('type','parameter');
global_parameter = ...
test_pm_riby_2( ...
 global_parameter ...
,fname_prefix ...
,fname_infix ...
,Pixel_Spacing ...
,fname_nopath_volume ...
,fname_nopath_micrograph ...
);
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%;
disp('returning');return;
end;%if (nargin<1);

% try: ;
% global_parameter=[];fname_prefix='ribosome_yeast_Micrograph';fname_infix='riby';Pixel_Spacing=1.06;fname_nopath_volume='parsed_6Q8Y_whole_LSU_match3_sim_60.mrc';fname_nopath_micrograph='147_Mar12_12.21.27_159_0.mrc';

verbose=1;
if (verbose); disp(sprintf(' %% [entering test_pm_riby_2]')); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
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
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,fname_prefix);
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
 fname_nopath_volume ... %<-- class 0. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;

%%%%%%%%;
% First create consensus volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_pm);
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
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
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
fname_fig = sprintf('%s_jpg/a_x_u_base_',dir_pm);
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
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','n_x_u_pack','n_pack','pack_row_ij_','pack_col_ij_','pack_val_ij_','x_u_pack_' ...
     ,'x_u_0_','x_u_1_','x_u_2_','x_u_0___','x_u_1___','x_u_2___','n_x_u','n_xxx_u','xxx_u_weight_','a_x_u_base_' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
dx_pack = diameter_x_c/max(1,n_x_u_pack);

%%%%%%%%;
% simple visualization of a_base. ;
%%%%%%%%
flag_plot=1;
if flag_plot;
fname_fig = sprintf('%s_jpg/a_x_u_base_FIGB_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
prctile_ = [94.0:0.5:99.5]; n_prctile = numel(prctile_);
figure(1);clf;figbig;
p_row = 3; p_col = 4; ns=0;
for nprctile=0:n_prctile-1;
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_p = prctile_(1+nprctile);
isosurface_f_x_c_0(a_x_u_base_,tmp_p);
title(sprintf('a base %.1f',tmp_p));
end;%for nprctile=0:n_prctile-1;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_plot;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_p_quad_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic;
k_p_r_max = 1*48/(2*pi); k_eq_d = 1.0/(2*pi);
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
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_base_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_base_(:)-a_x_u_reco_)/fnorm(a_x_u_base_(:))));
disp(sprintf(' %% at this point one should ensure that a_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
save(fname_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_' ...
     ,'k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_' ...
     ,'weight_3d_k_all_','weight_shell_k_' ...
     ,'n_k_p_r','k_p_r_' ...
     ,'weight_3d_k_p_r_' ...
     ,'k_c_0_all_','k_c_1_all_','k_c_2_all_' ...
     ,'J_node_','J_weight_','J_chebfun_','J_polyval_' ...
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
fname_fig = sprintf('%s_jpg/a_k_p_quad_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
plot(k_p_r_all_,log10(abs(a_k_p_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
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
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_Y_quad_.mat',dir_pm);
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
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_))); %<-- this should be 2-3 digits. ;
%%%%%%%%;
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
save(fname_mat ...
     ,'l_max_' ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'Y_l_val_','Y_m_val_','Y_k_val_','weight_Y_' ...
     ,'a_k_Y_quad_' ...
     ,'a_k_p_reco_' ...
     ,'a_k_Y_quad__' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_Y_quad_A',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
subplot(1,3,1); plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,2); plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,3); plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
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
fname_fig = sprintf('%s_jpg/a_k_Y_quad_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_)),[-10,0],colormap_beach());
xlabel('l_val','Interpreter','none');
ylabel('m_val','Interpreter','none');
zlabel('k_p_r','Interpreter','none');
axis vis3d; view([-65,+20]);
title('a_k_Y_quad_','Interpreter','none');
%figbig;
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
fname_mat = sprintf('%s_mat/S_k_p_wkS__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_w_max = 2*(l_max_max+1);
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
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
 0*verbose ...
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
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% get_template_1: %0.2fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%;
tmp_t = tic();
S_x_u_xxS___ = zeros(n_x_u_pack,n_x_u_pack,n_S);
for nS=0:n_S-1;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_x_u_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum );
S_x_u_xx__ = reshape(S_x_u_xx_,[n_x_u_pack,n_x_u_pack]);
S_x_u_xxS___(:,:,1+nS) = S_x_u_xx__;
end;%for nS=0:n_S-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% S_x_u_xxS___: %0.2fs',tmp_t)); end;
%%%%;
S_mask_threshold = tolerance_master;
M_x_u_xxS___ = zeros(n_x_u_pack,n_x_u_pack,n_S); %<-- template masks. ;
M_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_x_u_xx__ = S_x_u_xxS___(:,:,1+nS);
S_x_u_xx_ = reshape(S_x_u_xx__,[n_x_u_pack^2,1]);
[tmp_val_srt_,tmp_ij_ori_from_srt_] = sort(abs(S_x_u_xx_).^2,'descend');
tmp_val_srt_cdf_ = cumsum(tmp_val_srt_)/max(1e-12,sum(tmp_val_srt_));
tmp_ij_cut = min(find(sqrt(tmp_val_srt_cdf_)>1-S_mask_threshold));
tmp_ij_ori_from_srt_sub_ = tmp_ij_ori_from_srt_(1:tmp_ij_cut);
M_x_u_xx_ = zeros(n_x_u_pack^2,1); M_x_u_xx_(tmp_ij_ori_from_srt_sub_) = 1.0;
M_x_u_xx__ = reshape(M_x_u_xx_,[n_x_u_pack,n_x_u_pack]);
M_k_p_wk_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,M_x_u_xx__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack*n_x_u_pack)*dx_pack*dx_pack;
%N_x_u_xx_ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,M_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum );
%N_x_u_xx__ = reshape(N_x_u_xx_,[n_x_u_pack,n_x_u_pack]);
M_x_u_xxS___(:,:,1+nS) = M_x_u_xx__;
M_k_p_wkS__(:,1+nS) = M_k_p_wk_;
end;%for nS=0:n_S-1;
%%%%;
save(fname_mat ...
     ,'n_w_max','template_k_eq_d','viewing_k_eq_d' ...
     ,'S_k_p_wkS__' ...
     ,'S_x_u_xxS___' ...
     ,'S_mask_threshold' ...
     ,'M_k_p_wkS__' ...
     ,'M_x_u_xxS___' ...
     ,'n_w_' ...
     ,'weight_2d_k_p_r_' ...
     ,'weight_2d_wk_' ...
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
disp(sprintf(' %% %s found, not creating',fname_fig));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_x_u_xxS_vs_M_x_u_xxS_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;figbeach;
p_row = 4; p_col = 6; np=0; n_l = floor(p_row*p_col/2);
for nl=0:n_l-1;
nS = floor(n_S*nl/n_l);
subplot(p_row,p_col,1+np);np=np+1;
imagesc(S_x_u_xxS___(:,:,1+nS)); axis image; axisnotick; title(sprintf('S %.3d',nS));
subplot(p_row,p_col,1+np);np=np+1;
imagesc(M_x_u_xxS___(:,:,1+nS)); axis image; axisnotick; title(sprintf('M %.3d',nS));
end;%for nl=0:n_l-1;
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_p_wkS__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p_wkS__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wkS__(:,1+nS)),Slim_,colormap_80s);
axis equal; axisnotick; title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('Sample S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now convert templates to S_k_q__. ;
%%%%%%%%;
S_k_q__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_q__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p_wkS__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(S_k_q__(:,1+nS)),Slim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_q__(:,1+nS))/max(abs(S_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
% test pm_template_2 vs get_template_1. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 tmp_S_k_p_wkS__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_k_Y_quad__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
tmp_S_k_p_wkS__ = reshape(tmp_S_k_p_wkS__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
disp(sprintf(' %% S_k_p_wkS__ vs tmp_S_k_p_wkS__: %0.16f',fnorm(S_k_p_wkS__-tmp_S_k_p_wkS__)/fnorm(S_k_p_wkS__)));
end;%if flag_check;

%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/O_x_u_pack_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
fname_mrc = sprintf('%s/%s',dir_data_star,fname_nopath_micrograph);
O_x_u_load_ = cast(ReadMRC(fname_mrc),'double');
n_O_x_u_ = size(O_x_u_load_);
n_O_x_u_pack_ = n_O_x_u_/n_pack;
pack_0_row_ij_ = zeros(n_O_x_u_pack_(1+0),1);
pack_0_col_ij_ = zeros(n_O_x_u_pack_(1+0),1);
pack_0_val_ij_ = zeros(n_O_x_u_pack_(1+0),1);
na=0;
for nO_x_u=0:n_O_x_u_(1+0)-1;
pack_0_row_ij_(1+na) = 1+nO_x_u;
pack_0_col_ij_(1+na) = 1+floor(nO_x_u/n_pack);
pack_0_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nO_x_u=0:n_O_x_u_(1+0)-1;
pack_1_row_ij_ = zeros(n_O_x_u_pack_(1+1),1);
pack_1_col_ij_ = zeros(n_O_x_u_pack_(1+1),1);
pack_1_val_ij_ = zeros(n_O_x_u_pack_(1+1),1);
na=0;
for nO_x_u=0:n_O_x_u_(1+1)-1;
pack_1_row_ij_(1+na) = 1+nO_x_u;
pack_1_col_ij_(1+na) = 1+floor(nO_x_u/n_pack);
pack_1_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nO_x_u=0:n_O_x_u_(1+1)-1;
x_u_pack_0_ = sparse(pack_0_row_ij_,pack_0_col_ij_,pack_0_val_ij_,n_O_x_u_(1+0),n_O_x_u_pack_(1+0));
x_u_pack_1_ = sparse(pack_1_row_ij_,pack_1_col_ij_,pack_1_val_ij_,n_O_x_u_(1+1),n_O_x_u_pack_(1+1));
O_x_u_pack_ = transpose(x_u_pack_0_)*O_x_u_load_*x_u_pack_1_;
if flag_invert; O_x_u_pack_ = -O_x_u_pack_; end;
n_M = numel(O_x_u_pack_);
%%%%%%%%;
% Now subtract off a parabolic fit. ;
%%%%%%%%;
tmp_j0_ = transpose(1:n_O_x_u_pack_(1+0));
tmp_j1_ = transpose(1:n_O_x_u_pack_(1+1));
[tmp_j0__,tmp_j1__] = ndgrid(tmp_j0_,tmp_j1_);
tmp_RHS_00 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_01 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^1) , 'all' );
tmp_RHS_10 = sum( O_x_u_pack_ .* (tmp_j0__.^1) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_11 = sum( O_x_u_pack_ .* (tmp_j0__.^1) .* (tmp_j1__.^1) , 'all' );
tmp_RHS_02 = sum( O_x_u_pack_ .* (tmp_j0__.^0) .* (tmp_j1__.^2) , 'all' );
tmp_RHS_20 = sum( O_x_u_pack_ .* (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_RHS_ = [ ...
 tmp_RHS_00 ...
;tmp_RHS_01 ...
;tmp_RHS_10 ...
;tmp_RHS_11 ...
;tmp_RHS_02 ...
;tmp_RHS_20 ...
];
tmp_LHS_00 = sum( (tmp_j0__.^0) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_01 = sum( (tmp_j0__.^0) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_10 = sum( (tmp_j0__.^1) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_11 = sum( (tmp_j0__.^1) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_02 = sum( (tmp_j0__.^0) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_20 = sum( (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_12 = sum( (tmp_j0__.^1) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_21 = sum( (tmp_j0__.^2) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_03 = sum( (tmp_j0__.^0) .* (tmp_j1__.^3) , 'all' );
tmp_LHS_30 = sum( (tmp_j0__.^3) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_20 = sum( (tmp_j0__.^2) .* (tmp_j1__.^0) , 'all' );
tmp_LHS_22 = sum( (tmp_j0__.^2) .* (tmp_j1__.^2) , 'all' );
tmp_LHS_13 = sum( (tmp_j0__.^1) .* (tmp_j1__.^3) , 'all' );
tmp_LHS_31 = sum( (tmp_j0__.^3) .* (tmp_j1__.^1) , 'all' );
tmp_LHS_04 = sum( (tmp_j0__.^0) .* (tmp_j1__.^4) , 'all' );
tmp_LHS_40 = sum( (tmp_j0__.^4) .* (tmp_j1__.^0) , 'all' );
tmp_LHS__ = [ ...
  tmp_LHS_00 , tmp_LHS_01 , tmp_LHS_10 , tmp_LHS_11 , tmp_LHS_02 , tmp_LHS_20 ...
; tmp_LHS_01 , tmp_LHS_02 , tmp_LHS_11 , tmp_LHS_12 , tmp_LHS_03 , tmp_LHS_21 ...
; tmp_LHS_10 , tmp_LHS_11 , tmp_LHS_20 , tmp_LHS_21 , tmp_LHS_12 , tmp_LHS_30 ...
; tmp_LHS_11 , tmp_LHS_12 , tmp_LHS_21 , tmp_LHS_22 , tmp_LHS_13 , tmp_LHS_31 ...
; tmp_LHS_02 , tmp_LHS_03 , tmp_LHS_12 , tmp_LHS_13 , tmp_LHS_04 , tmp_LHS_22 ...
; tmp_LHS_20 , tmp_LHS_21 , tmp_LHS_30 , tmp_LHS_31 , tmp_LHS_22 , tmp_LHS_40 ...
];
tmp_axx_ = tmp_LHS__ \ tmp_RHS_;
a00 = tmp_axx_(1+0);
a01 = tmp_axx_(1+1);
a10 = tmp_axx_(1+2);
a11 = tmp_axx_(1+3);
a02 = tmp_axx_(1+4);
a20 = tmp_axx_(1+5);
tmp_p = @(j0,j1) ...
  a00.*(j0.^0).*(j1.^0) ...
+ a01.*(j0.^0).*(j1.^1) ...
+ a10.*(j0.^1).*(j1.^0) ...
+ a11.*(j0.^1).*(j1.^1) ...
+ a02.*(j0.^0).*(j1.^2) ...
+ a20.*(j0.^2).*(j1.^0) ...
;
p_x_u_pack_ = tmp_p(tmp_j0__,tmp_j1__);
P_x_u_pack_ = O_x_u_pack_ - p_x_u_pack_;
%%%%;
save(fname_mat ...
,'n_O_x_u_','n_O_x_u_pack_','n_pack','O_x_u_pack_','n_M' ...
,'p_x_u_pack_','P_x_u_pack_' ...
);
%%%%%%%%;
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/O_x_u_pack_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,3,1); imagesc(O_x_u_pack_); axis image; axisnotick; title('O_x_u_pack_','Interpreter','none');
subplot(1,3,2); imagesc(p_x_u_pack_); axis image; axisnotick; title('p_x_u_pack_','Interpreter','none');
subplot(1,3,3); imagesc(P_x_u_pack_); axis image; axisnotick; title('O_x_u_pack_-p_x_u_pack_','Interpreter','none');
sgtitle(fname_fig,'Interpreter','none');
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
%%%%%%%%;
fname_mat = sprintf('%s_mat/CTF_k_p_wkC__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_CTF = 1;
index_nCTF_from_nM_ = sparse(n_M,1);
Voltage_CTF = 300; %<-- set to 300kV: see https://www.ebi.ac.uk/empiar/EMPIAR-10998/ ; 
DefocusU_CTF = 4284; %<-- see kexin email: defocus 1 = 4284A ;
DefocusV_CTF = 4089; %<-- see kexin email: defocus 1 = 4089A ;
DefocusAngle_CTF = 18.39; %<-- see kexin email: defocus angle = 18.39A. ;
SphericalAberration_CTF = 2.7; %<-- Spherical aberration is usually set to 2.7mm. ;
AmplitudeContrast_CTF = 0.07; %<-- and amplitude contrast to 0.07. ;
%%%%%%%%;
n_x_M_u = n_x_u_pack;
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
if (mod(nCTF,100)==0); disp(sprintf(' %% nCTF %d/%d',nCTF,n_CTF)); end;
CTF_Spherical_Aberration = SphericalAberration_CTF;% spherical aberration of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = Voltage_CTF;% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = DefocusU_CTF;% defocus values (in Angstroms) ;
CTF_Defocus_V = DefocusV_CTF;% defocus values (in Angstroms) ;
CTF_Defocus_Angle = DefocusAngle_CTF;% angle of astigmatism ;
CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
CTF_Amplitude_Contrast = AmplitudeContrast_CTF;% CTF_Amplitude Contrast ;
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
tmp_ctf_value = ...
niko_ctf( ...
 CTF_Spherical_Aberration ... 
,CTF_lambda ... 
,tmp_w1 ... 
,tmp_w2 ... 
,CTF_Defocus_U ... 
,CTF_Defocus_V ... 
,CTF_Defocus_Angle ... 
,CTF_lambda_per_box/pi ... 
,tmp_k_c_1 ... 
,tmp_k_c_2 ...
);
CTF_k_p_wkC__(1+na,1+nCTF) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% normalize each CTF. ;
%%%%%%%%;
for nCTF=0:n_CTF-1;
CTF_k_p_wk_ = CTF_k_p_wkC__(:,1+nCTF);
CTF_k_p_avg = sum(CTF_k_p_wk_.*weight_2d_wk_)/(pi*k_p_r_max^2)*(4*pi^2);
%CTF_k_p_wk_ = CTF_k_p_wk_ - CTF_k_p_avg; %<-- do not subtract off the average. ;
CTF_k_p_std = sqrt( sum(abs(CTF_k_p_wk_).^2.*weight_2d_wk_)/(pi*k_p_r_max^2)*(4*pi^2) );
CTF_k_p_wk_ = CTF_k_p_wk_/max(1e-12,CTF_k_p_std);
CTF_k_p_wkC__(:,1+nCTF) = CTF_k_p_wk_;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r_kC__(1+nk_p_r,1+nCTF) = mean(CTF_k_p_wkC__(1+tmp_index_,1+nCTF));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+0);
CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__,2);
PSF_avg_x_u_xx_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,CTF_avg_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum;
PSF_avg_x_u_xx__ = real(reshape(PSF_avg_x_u_xx_,[n_x_u_pack,n_x_u_pack]));
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_wk_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_k_(1+nk_p_r) = mean(CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_kM__ = CTF_k_p_r_kC__(:,1+0);
SCTF_ = svd(CTF_k_p_r_kM__);
n_CTF_rank = max(find(SCTF_/max(SCTF_)>tolerance_master));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r_kM__,n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;
% Now determine the CTF cross correlation. ;
%%%%%%%%;
tmp_CTF_avg_k_p_wk_ = mean(CTF_k_p_wkC__(:,:),2);
tmp_CTF_avg_k_p_r_k_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_k_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_wk_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xavg__ = tmp_CTF_avg_k_p_r_k_ * transpose(tmp_CTF_avg_k_p_r_k_);
CTF_k_p_r_xcor__ = CTF_k_p_r_kC__(:,1) * transpose(CTF_k_p_r_kC__(:,1));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/CTF_k_p_xcor__',dir_pm);
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
     ,'CTF_k_p_r_k_' ...
     ,'CTF_k_p_r_kC__' ...
     ,'CTF_k_p_r_kM__' ...
     ,'n_CTF_rank' ...
     ,'SCTF_','UCTF_kc__','VSCTF_Mc__' ...
     ,'CTF_avg_k_p_wk_' ...
     ,'PSF_avg_x_u_xx__' ...
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
fname_fig = sprintf('%s_jpg/CTF_k_p_r_xxxx__',dir_pm);
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
% Now calculate the idealized principal-modes for each nCTF. ;
%%%%%%%%;
X_2d_xavg_d0_kkC___ = zeros(n_k_p_r,n_k_p_r,n_CTF);
X_2d_xavg_d0_weight_rC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
tmp_CTF_k_p_r_xavg_k_ = CTF_k_p_r_kC__(:,1+nCTF);
tmp_CTF_k_p_r_xavg_kk__ = tmp_CTF_k_p_r_xavg_k_*transpose(tmp_CTF_k_p_r_xavg_k_);
tmp_delta_sigma = 0;
[X_2d_xavg_d0_kk__,X_2d_xavg_d0_weight_r_] = principled_marching_cost_matrix_6(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,tmp_CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
X_2d_xavg_d0_kkC___(:,:,1+nCTF) = X_2d_xavg_d0_kk__;
X_2d_xavg_d0_weight_rC__(:,1+nCTF) = X_2d_xavg_d0_weight_r_;
clear X_2d_xavg_d0_kk__ X_2d_xavg_d0_weight_r_ ;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_2d_xavg_d0_kkC___FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figsml;figbeach();
nCTF=0;
imagesc(X_2d_xavg_d0_kkC___(:,:,1+nCTF)); axis image; axisnotick;
sgtitle(fname_fig,'Interpreter','none');
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
% Now calculate masked correlations. ;
% Q_sub_x_u_pack_ is obtained by: ;
% (i) first extracting a smaller subset of the full micrograph. ;
% (ii) then convolving it with the PSF. ;
% Now, for each location, we normalize S and then compute the following: ;
% <Q1,S> = innerproduct between Q^1 and S. ;
% < 1,S> = S_sum = sum of template. ;
% <Q1,M> = innerproduct between Q^1 and M = sum of Q^1 over support of mask. ;
% <Q2,M> = innerproduct between Q^2 and M = sum of Q^2 over support of mask. ;
% < 1,M> = M_sum = sum of mask (i.e., area of support of template). ;
% With these terms we can compute: ;
% Q1_avg = <Q1,M>/M_sum ;
% Q1_var = <Q2,M>/M_sum - Q1_avg^2 ;
% Q1_std = sqrt(Q1_var) ;
% <Q1_nrm,S> = <(Q1 - Q1_avg)/Q1_std,S> = <Q1,S>/Q1_std - <1,S>*Q1_avg/Q1_std;
%%%%%%%%;
n_Q_sub_x_u_pack_ = floor(n_O_x_u_pack_/4);
Q_sub_ij_0_ = floor(n_O_x_u_pack_(1+0)*1.0/4) + [1:n_Q_sub_x_u_pack_(1+0)];
Q_sub_ij_1_ = floor(n_O_x_u_pack_(1+1)*2.5/4) + [1:n_Q_sub_x_u_pack_(1+1)];
Q_x_u_pack_ = conv2(P_x_u_pack_,PSF_avg_x_u_xx__,'same');
Q_sub_x_u_pack_ = Q_x_u_pack_(Q_sub_ij_0_,Q_sub_ij_1_);
flag_disp=1;
if flag_disp;
tmp_Q_x_u_pack_ = Q_x_u_pack_;
tmp_p = prctile(tmp_Q_x_u_pack_,100,'all');
tmp_Q_x_u_pack_(min(Q_sub_ij_0_),Q_sub_ij_1_) = tmp_p;
tmp_Q_x_u_pack_(max(Q_sub_ij_0_),Q_sub_ij_1_) = tmp_p;
tmp_Q_x_u_pack_(Q_sub_ij_0_,min(Q_sub_ij_1_)) = tmp_p;
tmp_Q_x_u_pack_(Q_sub_ij_0_,max(Q_sub_ij_1_)) = tmp_p;
fname_fig = sprintf('%s_jpg/Q_sub_x_u_pack_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;colormap('cool');
subplot(1,2,1); imagesc(tmp_Q_x_u_pack_); axis image; axisnotick; title('Q','Interpreter','none');
subplot(1,2,2); imagesc(Q_sub_x_u_pack_); axis image; axisnotick; title('Q_sub','Interpreter','none');
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_disp;
%%%%%%%%;
k_ref_p_r_max = 2*k_p_r_max; k_ref_eq_d = k_eq_d;
[ ...
 n_k_ref_p_r ...
,k_ref_p_r_ ...
,weight_3d_k_ref_p_r_ ...
] = ...
get_weight_3d_1( ...
 0*verbose ...
,k_ref_p_r_max ...
,k_ref_eq_d ...
);
%%%%;
l_ref_max_upb = round(2*pi*k_ref_p_r_max); l_ref_max_max = l_ref_max_upb;
n_w_ref_0in_ = 2*(l_ref_max_max+1)*ones(n_k_ref_p_r,1);
template_k_eq_d = -1;
[ ...
 n_w_ref_ ...
,weight_2d_ref_k_p_r_ ...
,weight_2d_ref_wk_ ...
] = ...
get_weight_2d_2( ...
 0*verbose ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,template_k_eq_d ...
,n_w_ref_0in_ ...
);
%%%%%%%%;
if ~exist('scatter_from_tensor_order','var'); scatter_from_tensor_order = []; end;
if ~exist('scatter_from_tensor_zswk___','var'); scatter_from_tensor_zswk___ = []; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_t_start = tic();
nw_stride = 2;
n_gamma_z = floor(max(n_w_ref_)/nw_stride);
X_QS_xxzS____ = zeros(n_Q_sub_x_u_pack_(1+0),n_Q_sub_x_u_pack_(1+1),n_gamma_z,n_S);
for nS=0:n_S-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (verbose>-1);
if mod(nS,8)==0;
tmp_t_final = toc(tmp_t_start);
disp(sprintf(' %% nS %d/%d t %0.2fs/%0.2fs',nS,n_S,tmp_t_final,tmp_t_final*n_S/max(1,nS)));
end;%if mod(nS,8)==0;
end;%if (verbose>-1);
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
S_x_u_xx__ = real( interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_u_pack*n_x_u_pack) * n_w_sum );
M_k_p_wk_ = M_k_p_wkS__(:,1+nS); %<-- but do not normalize the mask M_k_p_wk_. ;
M_x_u_xx__ = M_x_u_xxS___(:,:,1+nS);
M_num = sum(M_x_u_xx__,'all');
S_x_u_avg = sum(S_x_u_xx__.*M_x_u_xx__,'all')/max(1,M_num);
S_x_u_nrm_xx__ = (S_x_u_xx__ - S_x_u_avg).*M_x_u_xx__;
S_x_u_var = sum(S_x_u_nrm_xx__.^2,'all')/max(1,M_num);
S_x_u_std = sqrt(S_x_u_var);
S_x_u_nrm_xx__ = S_x_u_nrm_xx__/max(1e-12,S_x_u_std);
S_k_p_nrm_wk_ = interp_x_c_to_k_p_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,S_x_u_nrm_xx__,n_k_p_r,k_p_r_,n_w_)*sqrt(n_x_u_pack*n_x_u_pack)*dx_pack*dx_pack;
S_x_u_nrm_num = sum(S_x_u_nrm_xx__,'all');
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
% compare effects of increasing k_ref_p_r_max. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 ~ ...
,Q1M_A_x_u_xxz___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 [] ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,M_x_u_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
tmp_t = tic();
[ ...
 ~ ...
,Q1M_B_x_u_xxz___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 [] ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,M_x_u_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,[] ...
,weight_2d_wk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
Q1M_C_x_u_xx__ = conv2(Q_sub_x_u_pack_.^1,M_x_u_xx__,'same');
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1); imagesc(Q1M_C_x_u_xx__); axis image; axisnotick; title('C'); colorbar;
subplot(1,3,2); imagesc(Q1M_A_x_u_xxz___(:,:,1)); axis image; axisnotick; title('A'); colorbar;
subplot(1,3,3); imagesc(Q1M_B_x_u_xxz___(:,:,1)); axis image; axisnotick; title('B'); colorbar;
figure(1+nf);nf=nf+1;clf;figmed;
n_h = 128; np=0;
hlim_CA_ = [ ...
 min(min(Q1M_C_x_u_xx__,[],'all'),min(Q1M_A_x_u_xxz___(:,:,1),[],'all')) ...
 max(max(Q1M_C_x_u_xx__,[],'all'),max(Q1M_A_x_u_xxz___(:,:,1),[],'all')) ...
];
h2_CA__ = hist2d_0(Q1M_C_x_u_xx__,Q1M_A_x_u_xxz___(:,:,1),n_h,n_h,hlim_CA_,hlim_CA_);
subplot(2,2,1+np);np=np+1;
imagesc(log2(1+h2_CA__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('Q1M_C_x_u_xx__ vs Q1M_A_x_u_xxz___(:,:,1)','Interpreter','none');
subplot(2,2,1+np);np=np+1;
plot(Q1M_C_x_u_xx__(:),reshape(Q1M_A_x_u_xxz___(:,:,1),[prod(n_Q_sub_x_u_pack_),1]),'.');
axisnotick; xlim(hlim_CA_);ylim(hlim_CA_); axis square;
title('Q1M_C_x_u_xx__ vs Q1M_A_x_u_xxz___(:,:,1)','Interpreter','none');
hlim_CB_ = [ ...
 min(min(Q1M_C_x_u_xx__,[],'all'),min(Q1M_B_x_u_xxz___(:,:,1),[],'all')) ...
 max(max(Q1M_C_x_u_xx__,[],'all'),max(Q1M_B_x_u_xxz___(:,:,1),[],'all')) ...
];
h2_CB__ = hist2d_0(Q1M_C_x_u_xx__,Q1M_B_x_u_xxz___(:,:,1),n_h,n_h,hlim_CB_,hlim_CB_);
subplot(2,2,1+np);np=np+1;
imagesc(log2(1+h2_CB__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('Q1M_C_x_u_xx__ vs Q1M_B_x_u_xxz___(:,:,1)','Interpreter','none');
subplot(2,2,1+np);np=np+1;
plot(Q1M_C_x_u_xx__(:),reshape(Q1M_B_x_u_xxz___(:,:,1),[prod(n_Q_sub_x_u_pack_),1]),'.');
axisnotick; xlim(hlim_CB_);ylim(hlim_CB_); axis square;
title('Q1M_C_x_u_xx__ vs Q1M_B_x_u_xxz___(:,:,1)','Interpreter','none');
%%%%%%%%;
tmp_t = tic();
[ ...
 ~ ...
,Q1S_A_x_u_xxz___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 [] ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,S_x_u_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
tmp_t = tic();
[ ...
 ~ ...
,Q1S_B_x_u_xxz___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 [] ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,S_x_u_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,[] ...
,weight_2d_wk_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
Q1S_C_x_u_xx__ = conv2(Q_sub_x_u_pack_.^1,S_x_u_xx__,'same');
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,3,1); imagesc(Q1S_C_x_u_xx__); axis image; axisnotick; title('C'); colorbar;
subplot(1,3,2); imagesc(Q1S_A_x_u_xxz___(:,:,1)); axis image; axisnotick; title('A'); colorbar;
subplot(1,3,3); imagesc(Q1S_B_x_u_xxz___(:,:,1)); axis image; axisnotick; title('B'); colorbar;
figure(1+nf);nf=nf+1;clf;figmed;
n_h = 128; np=0;
hlim_CA_ = [ ...
 min(min(Q1S_C_x_u_xx__,[],'all'),min(Q1S_A_x_u_xxz___(:,:,1),[],'all')) ...
 max(max(Q1S_C_x_u_xx__,[],'all'),max(Q1S_A_x_u_xxz___(:,:,1),[],'all')) ...
];
h2_CA__ = hist2d_0(Q1S_C_x_u_xx__,Q1S_A_x_u_xxz___(:,:,1),n_h,n_h,hlim_CA_,hlim_CA_);
subplot(2,2,1+np);np=np+1;
imagesc(log2(1+h2_CA__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('Q1S_C_x_u_xx__ vs Q1S_A_x_u_xxz___(:,:,1)','Interpreter','none');
subplot(2,2,1+np);np=np+1;
plot(Q1S_C_x_u_xx__(:),reshape(Q1S_A_x_u_xxz___(:,:,1),[prod(n_Q_sub_x_u_pack_),1]),'.');
axisnotick; xlim(hlim_CA_);ylim(hlim_CA_); axis square;
title('Q1S_C_x_u_xx__ vs Q1S_A_x_u_xxz___(:,:,1)','Interpreter','none');
hlim_CB_ = [ ...
 min(min(Q1S_C_x_u_xx__,[],'all'),min(Q1S_B_x_u_xxz___(:,:,1),[],'all')) ...
 max(max(Q1S_C_x_u_xx__,[],'all'),max(Q1S_B_x_u_xxz___(:,:,1),[],'all')) ...
];
h2_CB__ = hist2d_0(Q1S_C_x_u_xx__,Q1S_B_x_u_xxz___(:,:,1),n_h,n_h,hlim_CB_,hlim_CB_);
subplot(2,2,1+np);np=np+1;
imagesc(log2(1+h2_CB__)); axis image; axisnotick; colorbar; set(gca,'ydir','normal');
title('Q1S_C_x_u_xx__ vs Q1S_B_x_u_xxz___(:,:,1)','Interpreter','none');
subplot(2,2,1+np);np=np+1;
plot(Q1S_C_x_u_xx__(:),reshape(Q1S_B_x_u_xxz___(:,:,1),[prod(n_Q_sub_x_u_pack_),1]),'.');
axisnotick; xlim(hlim_CB_);ylim(hlim_CB_); axis square;
title('Q1S_C_x_u_xx__ vs Q1S_B_x_u_xxz___(:,:,1)','Interpreter','none');
%%%%%%%%;
end;%if flag_check;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
if flag_check;
%%%%%%%%;
% calculate correlations. ;
%%%%%%%%;
Q1M_C_x_u_xx__ = conv2(Q_sub_x_u_pack_.^1,M_x_u_xx__,'same');
tmp_ij_r = 3; tmp_ij_c = 8;
tmp_Q_xx__ = Q_sub_x_u_pack_(tmp_ij_r + [1:n_x_u_pack],tmp_ij_c + [1:n_x_u_pack]);
tmp_Q1M_A = sum(tmp_Q_xx__.^1.*fliplr(flipud(M_x_u_xx__)),'all');
tmp_Q1M_B = Q1M_C_x_u_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_Q1M_A vs tmp_Q1M_B: %0.16f',fnorm(tmp_Q1M_A-tmp_Q1M_B)/fnorm(tmp_Q1M_A)));
tmp_ij = find(abs(Q1M_C_x_u_xx__-tmp_Q1M_A)<1e-12);
tmp_index = tmp_ij-1;
tmp_index_r = floor(mod(tmp_index,n_Q_sub_x_u_pack_(1+0)));
tmp_index_c = floor(tmp_index/n_Q_sub_x_u_pack_(1+0));
assert(tmp_ij_r==tmp_index_r+1-floor(n_x_u_pack/2));
assert(tmp_ij_c==tmp_index_c+1-floor(n_x_u_pack/2));
Q2M_C_x_u_xx__ = conv2(Q_sub_x_u_pack_.^2,M_x_u_xx__,'same');
tmp_Q2M_A = sum(tmp_Q_xx__.^2.*fliplr(flipud(M_x_u_xx__)),'all');
tmp_M_sum = sum(M_x_u_xx__,'all');
tmp_std_A = sqrt(tmp_Q2M_A/tmp_M_sum - (tmp_Q1M_A/tmp_M_sum).^2);
tmp_ = tmp_Q_xx__.^1.*fliplr(flipud(M_x_u_xx__)); tmp_ = tmp_(find(tmp_~=0));
tmp_std_B = std(tmp_,1);
disp(sprintf(' %% tmp_std_A vs tmp_std_B: %0.16f',fnorm(tmp_std_A-tmp_std_B)/fnorm(tmp_std_A)));
tmp_ij_flip_ = find(fliplr(flipud(M_x_u_xx__)));
tmp_S_nrm_ = fliplr(flipud(S_x_u_nrm_xx__)); tmp_S_nrm_ = tmp_S_nrm_(tmp_ij_flip_);
tmp_Q_ = tmp_Q_xx__(tmp_ij_flip_);
tmp_Q_nrm_ = (tmp_Q_-mean(tmp_Q_))/max(1e-12,std(tmp_Q_,1));
tmp_X_QS_B = dot(tmp_Q_nrm_,tmp_S_nrm_)/max(1,tmp_M_sum);
tmp_Q1S_A = sum(tmp_Q_xx__.^1.*fliplr(flipud(S_x_u_nrm_xx__)),'all');
tmp_Q_avg_A = tmp_Q1M_A/tmp_M_sum;
tmp_Q_var_A = tmp_Q2M_A/tmp_M_sum - tmp_Q_avg_A.^2;
tmp_Q_std_A = sqrt(tmp_Q_var_A);
tmp_X_QS_A = tmp_Q1S_A/tmp_Q_std_A/max(1,tmp_M_sum);
disp(sprintf(' %% tmp_X_QS_A vs tmp_X_QS_B: %0.16f',fnorm(tmp_X_QS_A-tmp_X_QS_B)/fnorm(tmp_X_QS_A)));
Q1S_C_x_u_xx__ = conv2(Q_sub_x_u_pack_.^1,S_x_u_xx__.*M_x_u_xx__,'same');
tmp_Q_xx__ = Q_sub_x_u_pack_(tmp_ij_r + [1:n_x_u_pack],tmp_ij_c + [1:n_x_u_pack]);
tmp_Q1S_A = sum(tmp_Q_xx__.^1.*fliplr(flipud(S_x_u_xx__.*M_x_u_xx__)),'all');
tmp_Q1S_B = Q1S_C_x_u_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_Q1S_A vs tmp_Q1S_B: %0.16f',fnorm(tmp_Q1S_A-tmp_Q1S_B)/fnorm(tmp_Q1S_A)));
tmp_Q_avg = sum(tmp_Q_xx__.^1.*fliplr(flipud(M_x_u_xx__)),'all')/tmp_M_sum;
tmp_Q_var = sum(tmp_Q_xx__.^2.*fliplr(flipud(M_x_u_xx__)),'all')/tmp_M_sum - tmp_Q_avg.^2;
tmp_Q_std = sqrt(tmp_Q_var);
tmp_S1_sum = sum(S_x_u_xx__.^1.*M_x_u_xx__,'all');
tmp_S2_sum = sum(S_x_u_xx__.^2.*M_x_u_xx__,'all');
tmp_S_avg = tmp_S1_sum/tmp_M_sum;
tmp_S_var = tmp_S2_sum/tmp_M_sum - tmp_S_avg.^2;
tmp_S_std = sqrt(tmp_S_var);
tmp_X_QS_A = sum( ((tmp_Q_xx__.^1-tmp_Q_avg)/tmp_Q_std) .* fliplr(flipud(((S_x_u_xx__-tmp_S_avg).*M_x_u_xx__)/tmp_S_std)) , 'all')/tmp_M_sum;
%tmp_X_QS_A = sum( ((tmp_Q_xx__.^1-1*tmp_Q_avg)/tmp_Q_std) .* fliplr(flipud(((S_x_u_xx__-1*tmp_S_avg).*M_x_u_xx__)/tmp_S_std)) , 'all')/tmp_M_sum;
Q_avg_xx__ = Q1M_C_x_u_xx__/tmp_M_sum;
tmp_Q_avg_A = tmp_Q_avg;
tmp_Q_avg_B = Q_avg_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_Q_avg_A vs tmp_Q_avg_B: %0.16f',fnorm(tmp_Q_avg_A-tmp_Q_avg_B)/fnorm(tmp_Q_avg_A)));
Q_var_xx__ = Q2M_C_x_u_xx__/tmp_M_sum - Q_avg_xx__.^2;
tmp_Q_var_A = tmp_Q_var;
tmp_Q_var_B = Q_var_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_Q_var_A vs tmp_Q_var_B: %0.16f',fnorm(tmp_Q_var_A-tmp_Q_var_B)/fnorm(tmp_Q_var_A)));
Q_std_xx__ = sqrt(Q_var_xx__);
tmp_Q_std_A = tmp_Q_std;
tmp_Q_std_B = Q_std_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_Q_std_A vs tmp_Q_std_B: %0.16f',fnorm(tmp_Q_std_A-tmp_Q_std_B)/fnorm(tmp_Q_std_A)));
X_QS_xx__ = (Q1S_C_x_u_xx__ - Q_avg_xx__*tmp_S_avg*tmp_M_sum)./max(1e-12,Q_std_xx__*tmp_S_std)/tmp_M_sum;
%X_QS_xx__ = (Q1S_C_x_u_xx__ - Q_avg_xx__*tmp_S1_sum)./max(1e-12,Q_std_xx__*tmp_S_std)/tmp_M_sum;
tmp_X_QS_B = X_QS_xx__(tmp_ij_r+floor(n_x_u_pack/2),tmp_ij_c+floor(n_x_u_pack/2));
disp(sprintf(' %% tmp_X_QS_A vs tmp_X_QS_B: %0.16f',fnorm(tmp_X_QS_A-tmp_X_QS_B)/fnorm(tmp_X_QS_A)));
figure(1+nf);nf=nf+1;clf;figmed;fig80s;
subplot(1,2,1); imagesc(Q_sub_x_u_pack_); 
axis image; axisnotick; title('Q_sub_x_u_pack_','Interpreter','none'); colorbar;
subplot(1,2,2); imagesc(X_QS_xx__);
axis image; axisnotick; title('X_QS_xx__','Interpreter','none'); colorbar;
end;%if flag_check;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
%%%%%%%%;
parameter_innerproduct = struct('type','parameter','nw_stride',2);
tmp_t = tic();
[ ...
 ~ ...
,Q1S_x_u_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,S_x_u_nrm_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ ...
 ~ ...
,Q1M_x_u_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^1 ...
,n_x_u_pack ...
,n_x_u_pack ...
,M_x_u_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
tmp_t = tic();
[ ...
 ~ ...
,Q2M_x_u_xxz___ ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
] = ...
innerproduct_O_x_u_vs_S_x_u_2( ...
 parameter_innerproduct ...
,n_Q_sub_x_u_pack_(1+0) ...
,n_Q_sub_x_u_pack_(1+1) ...
,Q_sub_x_u_pack_.^2 ...
,n_x_u_pack ...
,n_x_u_pack ...
,M_x_u_xx__ ...
,n_k_ref_p_r ...
,k_ref_p_r_ ...
,k_ref_p_r_max ...
,n_w_ref_ ...
,[] ...
,weight_2d_ref_wk_ ...
,[] ...
,scatter_from_tensor_order ...
,scatter_from_tensor_zswk___ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% innerproduct_O_x_u_vs_S_x_u_2: %0.2fs',tmp_t)); end;
%%%%;
Q_avg_xxz___ = Q1M_x_u_xxz___/M_num;
Q_var_xxz___ = Q2M_x_u_xxz___/M_num - Q_avg_xxz___.^2;
Q_std_xxz___ = sqrt(Q_var_xxz___);
X_QS_xxz___ = Q1S_x_u_xxz___./max(1e-12,Q_std_xxz___)/max(1,M_num);
X_QS_xxzS____(:,:,:,1+nS) = X_QS_xxz___;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nS=0:n_S-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% quick snapshot of subset of micrograph. ;
% comparing original (right) with CTF-modulated (left) with correlation (middle). ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_QS_max_xx__FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;colormap('cool');
subplot(1,3,1); imagesc(Q_sub_x_u_pack_);
colorbar;axis image;axisnotick; title('Q_sub_x_u_pack_','Interpreter','none');
set(gca,'ydir','normal');
subplot(1,3,2); imagesc(max(X_QS_xxzS____,[],[3,4]));
colorbar;axis image;axisnotick; title('max(X_QS_xxzS____,[],[3,4])','Interpreter','none');
set(gca,'ydir','normal');
subplot(1,3,3); imagesc(P_x_u_pack_(Q_sub_ij_0_,Q_sub_ij_1_));
colorbar; axis image; axisnotick; title('P_x_u_pack_(Q_sub_ij_0_,Q_sub_ij_1_)','Interpreter','none');
set(gca,'ydir','normal');
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_disp;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% quick snapshot of subset of micrograph. ;
% uses different colormaps for various viewing-angles. ;
%%%%%%%%;
[X_QS_max_xx__,X_QS_ij_max_xx__] = max(X_QS_xxzS____,[],[3,4],'linear');
X_QS_index_max_xx__ = X_QS_ij_max_xx__-1;
X_QS_index_max_xx__ = floor(X_QS_index_max_xx__/prod(n_Q_sub_x_u_pack_));
X_QS_index_z_xx__ = mod(X_QS_index_max_xx__,n_gamma_z);
X_QS_index_S_xx__ = floor(X_QS_index_max_xx__/n_gamma_z);
X_QS_gamma_z_xx__ = (2*pi*X_QS_index_z_xx__)/max(1,n_gamma_z);
%%%%;
X_QS_cz_xx3___ = zeros(prod(n_Q_sub_x_u_pack_),3);
tmp_cz__ = colormap_phase; n_cz = size(tmp_cz__,1);
for nxx=0:prod(n_Q_sub_x_u_pack_)-1;
ncz = max(0,min(n_cz-1,floor(n_cz*(X_QS_gamma_z_xx__(1+nxx)-0)/(2*pi))));
X_QS_cz_xx3___(1+nxx,:) = tmp_cz__(1+ncz,:);
end;%for nxx=0:prod(n_Q_sub_x_u_pack_)-1;
X_QS_cz_xx3___ = reshape(X_QS_cz_xx3___,[n_Q_sub_x_u_pack_(1+0),n_Q_sub_x_u_pack_(1+1),3]);
%%%%;
X_QS_azimu_b_xx__ = viewing_azimu_b_all_(1+X_QS_index_S_xx__);
X_QS_polar_a_xx__ = viewing_polar_a_all_(1+X_QS_index_S_xx__);
X_QS_cab_xx3___ = reshape(colormap_polar_a_azimu_b_2d(X_QS_polar_a_xx__,X_QS_azimu_b_xx__),[n_Q_sub_x_u_pack_(1+0),n_Q_sub_x_u_pack_(1+1),3]);
%%%%;
tmp_0_ = transpose(1:n_Q_sub_x_u_pack_(1+0));
tmp_1_ = transpose(1:n_Q_sub_x_u_pack_(1+1));
[tmp_0__,tmp_1__]=ndgrid(tmp_0_,tmp_1_);
tmp_x__ = bsxfun(@plus,repmat(transpose(tmp_1__(:)),[4,1]),[-0.5;+0.5;+0.5;-0.5]);
tmp_y__ = bsxfun(@plus,repmat(transpose(tmp_0__(:)),[4,1]),[-0.5;-0.5;+0.5;+0.5]);
%%%%;
fname_fig = sprintf('%s_jpg/X_QS_max_xx__FIGB',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;colormap('cool');
ax(1) = subplot(3,3,[7]); imagesc(1:64); axisnotick; title('correlation','Interpreter','none');
ax(2) = subplot(3,3,[8]); colormap_polar_a_azimu_b_2d; title('azimu_b,polar_a','Interpreter','none');
ax(3) = subplot(3,3,[9]); imagesc(1:64); axisnotick; title('gamma_z','Interpreter','none');
subplot(3,3,[3,6]);
tmp_p = patch(tmp_x__,tmp_y__,reshape(X_QS_cz_xx3___,[1,prod(n_Q_sub_x_u_pack_),3]));
set(tmp_p,'EdgeColor','none');
axis image;axisnotick; title('X_QS_cz_xx3___','Interpreter','none');
subplot(3,3,[2,5]);
tmp_p = patch(tmp_x__,tmp_y__,reshape(X_QS_cab_xx3___,[1,prod(n_Q_sub_x_u_pack_),3]));
set(tmp_p,'EdgeColor','none');
axis image;axisnotick; title('X_QS_cab_xx3___','Interpreter','none');
subplot(3,3,[1,4]);
imagesc(X_QS_max_xx__);set(gca,'ydir','normal');
axis image;axisnotick; title('X_QS_max_xx__','Interpreter','none');
colormap(ax(1),colormap('cool'));
colormap(ax(3),colormap_phase);
sgtitle(fname_fig,'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%if flag_disp;



verbose=1;
if (verbose); disp(sprintf(' %% [finished test_pm_riby_2]')); end;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
