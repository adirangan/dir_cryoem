%%%%%%%%;
% makes vfigs (i.e., volume-figs) for Principled_Marching_2?.tex. ;
%%%%%%%%;

%%%%%%%%;
% first recapitulates test_pm_trpv1_2. ;
%%%%%%%%;
table_data__ = { ...
'p28hRPT1_x0' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
};
n_experiment = 1; nexperiment=0;
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
global_parameter = struct('type','parameter');
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 1; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;

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
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,dir_nopath_data_star);
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
b_rho_x_u_pack_ = real(ifftn(fftn(a_rho_x_u_pack_).*exp(-i*2*pi*(K_u_0_*a_rho_x_c_0_avg + K_u_1_*a_rho_x_c_1_avg + K_u_2_*a_rho_x_c_2_avg))));
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

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/a_k_p_quad_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
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

flag_invert=0;

%%%%%%%%;
% Now generate figures for paper. ;
%%%%%%%%;

fname_mat = sprintf('%s_mat/pm_vfig_b_x_u_xi__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_b_k_Y_ = { ...
 'X_2d_xcor_d0_a1t0000p10r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p10r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p10r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p15r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p15r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p15r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p20r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p20r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p20r2_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p25r0_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p25r1_a_k_Y_.mat' ...
,'X_2d_xcor_d0_a1t0000p25r2_a_k_Y_.mat' ...
};
n_fname_b_k_Y = numel(fname_b_k_Y_);
b_k_Y_yi__ = zeros(n_lm_sum,n_fname_b_k_Y);
for nfname_b_k_Y=0:n_fname_b_k_Y-1;
fname_b_k_Y = fname_b_k_Y_{1+nfname_b_k_Y};
tmp_ = load(sprintf('%s_mat/%s',dir_pm,fname_b_k_Y));
b_k_Y_yi__(:,1+nfname_b_k_Y) = tmp_.a_k_Y_reco_;
clear tmp_;
end;%for nfname_b_k_Y=0:n_fname_b_k_Y-1;
%%%%%%%%;
b_k_p_ki__ = zeros(n_k_all,n_fname_b_k_Y);
b_x_u_xi__ = zeros(n_x_u_pack^3,n_fname_b_k_Y);
for nfname_b_k_Y=0:n_fname_b_k_Y-1;
tmp_b_k_Y_ = b_k_Y_yi__(:,1+nfname_b_k_Y);
tmp_t = tic;
[tmp_b_k_p_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,tmp_b_k_Y_);
tmp_t = toc(tmp_t); disp(sprintf(' %% tmp_b_k_Y_ --> tmp_b_k_p_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
tmp_b_x_u_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,tmp_b_k_p_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: tmp_b_x_u_ time %0.2fs',tmp_t));
tmp_b_x_u_ = reshape(tmp_b_x_u_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
b_k_p_ki__(:,1+nfname_b_k_Y) = tmp_b_k_p_(:);
b_x_u_xi__(:,1+nfname_b_k_Y) = tmp_b_x_u_(:);
end;%for nfname_b_k_Y=0:n_fname_b_k_Y-1;
save(fname_mat ...
     ,'fname_b_k_Y_','n_fname_b_k_Y' ...
     ,'b_k_Y_yi__' ...
     ,'b_k_p_ki__' ...
     ,'b_x_u_xi__' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_fig = sprintf('%s_jpg/pm_vfig_volume_FIGA_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
p_row = 3; p_col = 4; ns=0;
tmp_n_x = round(0.5352*n_x_u_pack/2);
tmp_x_ij = n_x_u_pack/2 + [-tmp_n_x:+tmp_n_x];
tmp_prct_0 = 75;
tmp_prct_1 = 90;
tmp_prct_2 = 97.5;
alim_ = mean(real(a_x_u_base_(tmp_x_ij,tmp_x_ij,tmp_x_ij)),'all') + std(real(a_x_u_base_(tmp_x_ij,tmp_x_ij,tmp_x_ij)),1,'all')*2.5*[-1,+1];
for pcol=0:p_col-1;
tmp_a_x_u_ = reshape(a_x_u_base_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
tmp_str = 'ground truth';
if (pcol> 0);
tmp_a_x_u_ = reshape(b_x_u_xi__(:,1+pcol-1),[n_x_u_pack,n_x_u_pack,n_x_u_pack]);
tmp_str = sprintf('sample %d',1+pcol-1);
end;%if (pcol> 0);
subplot(p_row,p_col,1+0*p_col+pcol);cla;
isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct_0);
axisnotick; xlabel([]); ylabel([]); zlabel([]); title(sprintf('%s %0.1f%%',tmp_str,tmp_prct_0))
subplot(p_row,p_col,1+1*p_col+pcol);cla;
isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct_1);
axisnotick; xlabel([]); ylabel([]); zlabel([]); title(sprintf('%s %0.1f%%',tmp_str,tmp_prct_1))
subplot(p_row,p_col,1+2*p_col+pcol);cla;
isosurface_f_x_u_0(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,tmp_x_ij),tmp_prct_2);
axisnotick; xlabel([]); ylabel([]); zlabel([]); title(sprintf('%s %0.1f%%',tmp_str,tmp_prct_2))
%subplot(p_row,p_col,1+1*p_col+pcol);cla;
%imagesc(squeeze(real(tmp_a_x_u_(tmp_x_ij,tmp_x_ij,round(n_x_u_pack/2)))),alim_);
%subplot(p_row,p_col,1+2*p_col+pcol);cla;
%imagesc(squeeze(real(tmp_a_x_u_(tmp_x_ij,round(n_x_u_pack/2),tmp_x_ij))),alim_);
end;%for pcol=0:p_col-1;
set(gcf,'Position',1+[0,0,1024+512,768]);
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now find principal-modes for ground-truth volume. ;
%%%%%%%%;
a_k_Y_mlk_ = a_k_Y_quad_;
a_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,a_k_Y_mlk_);
n_l_max = 1+l_max_max;
%%%%%%%%;
% calculate radial principal-modes. ;
%%%%%%%%;
tmp_t = tic();
X_kk__ = zeros(n_k_p_r,n_k_p_r);
X_weight_r_ = sqrt(weight_3d_k_p_r_);
for nk_p_r_0=0:n_k_p_r-1;
a_k_Y_ml0__ = a_k_Y_mlk___(:,:,1+nk_p_r_0);
X_weight_0 = X_weight_r_(1+nk_p_r_0);
for nk_p_r_1=0:n_k_p_r-1;
a_k_Y_ml1__ = a_k_Y_mlk___(:,:,1+nk_p_r_1);
X_weight_1 = X_weight_r_(1+nk_p_r_1);
X_weight_01 = X_weight_0*X_weight_1;
tmp_X = sum(conj(a_k_Y_ml0__(:,2:end)).*a_k_Y_ml1__(:,2:end),'all');
X_kk__(1+nk_p_r_0,1+nk_p_r_1) = real(tmp_X)*X_weight_01*(4*pi)^2;
end;%for nk_p_r_1=0:n_k_p_r-1;
end;%for nk_p_r_0=0:n_k_p_r-1;
tmp_X_kk__ = principled_marching_cost_matrix_2(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__);
disp(sprintf(' %% X_kk__ vs tmp_X_kk__: %0.16f',fnorm(X_kk__-tmp_X_kk__)/fnorm(X_kk__)));
clear tmp_X_kk__;
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% X_kk__: %0.3fs',tmp_t)); end;
tot_t_X_kk = tmp_t;
%%%%%%%%;
% Now compress k. ;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX__,SX__,VX__] = svds(X_kk__,n_UX_rank); SX_ = diag(SX__);
a_UX_Y_mln___ = zeros(n_m_max,n_l_max,n_UX_rank);
a_UX_Y_mln___ = reshape(reshape(a_k_Y_mlk___,[n_m_max*n_l_max,n_k_p_r])*diag(X_weight_r_)*UX__,[n_m_max,n_l_max,n_UX_rank]);
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% UX__: %0.3fs',tmp_t)); end;
tot_t_UX_a = tmp_t;
%%%%%%%%;
% Now calculate order principal-modes. ;
%%%%%%%%;
tmp_t = tic();
pm_n_UX_rank = n_UX_rank;
Z_ll__ = zeros(n_l_max,n_l_max);
for nl_max_0=1:n_l_max-1;
a_UX_Y_m0n__ = a_UX_Y_mln___(:,1+nl_max_0,1:pm_n_UX_rank);
for nl_max_1=1:n_l_max-1;
a_UX_Y_m1n__ = a_UX_Y_mln___(:,1+nl_max_1,1:pm_n_UX_rank);
tmp_Z = sum(conj(a_UX_Y_m0n__).*a_UX_Y_m1n__,'all');
Z_ll__(1+nl_max_0,1+nl_max_1) = real(tmp_Z)*(4*pi)^2;
end;%for nl_max_1=1:n_l_max-1;
end;%for nl_max_0=1:n_l_max-1;
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% Z_ll__: %0.3fs',tmp_t)); end;
tot_t_Z_ll = tmp_t;
%%%%%%%%;
tmp_t = tic();
tmp_Z_ll__ = zeros(n_l_max,n_l_max);
for nl_max_0=1:n_l_max-1;
a_k_Y_m0k__ = bsxfun(@times,a_k_Y_mlk___(:,1+nl_max_0,1:n_k_p_r),reshape(X_weight_r_,[1,1,n_k_p_r]));
for nl_max_1=1:n_l_max-1;
a_k_Y_m1k__ = bsxfun(@times,a_k_Y_mlk___(:,1+nl_max_1,1:n_k_p_r),reshape(X_weight_r_,[1,1,n_k_p_r]));
tmp_Z = sum(conj(a_k_Y_m0k__).*a_k_Y_m1k__,'all');
tmp_Z_ll__(1+nl_max_0,1+nl_max_1) = real(tmp_Z)*(4*pi)^2;
end;%for nl_max_1=1:n_l_max-1;
end;%for nl_max_0=1:n_l_max-1;
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% tmp_Z_ll__: %0.3fs',tmp_t)); end;
tot_t_Z_ll = tmp_t;
disp(sprintf(' %% Z_ll__ vs tmp_Z_ll__: %0.16f',fnorm(tmp_Z_ll__ - Z_ll__)/fnorm(Z_ll__)));
%%%%%%%%;
tmp_t = tic();
n_UZ_rank = n_l_max-0; %<-- to match norm. ;
[UZ__,SZ__,VZ__] = svds(Z_ll__,n_UZ_rank); SZ_ = diag(SZ__);
tmp_t = toc(tmp_t);
if (verbose); disp(sprintf(' %% UZ__: %0.3fs',tmp_t)); end;
tot_t_UZ = tmp_t;
%%%%%%%%;
flag_plot=0;
if flag_plot;
llim_ = [-20,-5];
%%%%;
figure(1);clf;figbig;figbeach;
for nk_p_r=0:n_k_p_r-1;
subplot(7,7,1+nk_p_r);
imagesc(log(abs(a_k_Y_mlk___(:,:,1+nk_p_r))),llim_)
end;%for nk_p_r=0:16-1;
title('a_k_Y_','Interpreter','none');
%%%%;
llim_ = [-20,-5];
figure(2);clf;figbig;figbeach;
for nUX_rank=0:n_UX_rank-1;
subplot(7,7,1+nUX_rank);
imagesc(log(abs(a_UX_Y_mln___(:,:,1+nUX_rank))),llim_);
end;%for nUX_rank=0:16-1;
title('a_UX_Y_','Interpreter','none');
%%%%;
end;%if flag_plot;

%%%%%%%%;
% Now perform stripped down calculation of innerproducts. ;
%%%%%%%%;
%%%%%%%%;
% calculate d_. ;
%%%%%%%%;
n_w_max = 2*(l_max_max+1);
n_beta = n_w_max; beta_ = transpose(linspace(-pi,+pi,n_beta+1)); beta_ = beta_(1:end-1);
t_0in = tic;
d_mmlb____ = zeros(n_m_max,n_m_max,n_l_max,n_beta);
d_mmzb____ = zeros(n_m_max,n_m_max,n_UZ_rank,n_beta);
for nbeta=0:n_beta-1;
if (verbose); if (mod(nbeta,16)==0); disp(sprintf(' %% nbeta %d/%d',nbeta,n_beta)); end; end;
beta = beta_(1+nbeta);
W_ = wignerd_b(l_max_max,-beta);
n_W_ = zeros(1,1+l_max_max); for (l_val=0:l_max_max); n_W_(1+l_val) = numel(W_{1+l_val}); end;
d_mml___ = zeros(n_m_max,n_m_max,1+l_max_max);
for l_val=0:l_max_max;
d_mml___(1+l_max_max + [-l_val:+l_val],1+l_max_max + [-l_val:+l_val],1+l_val) = W_{1+l_val};
end;%for l_val=0:l_max_max;
d_mmz___ = reshape(reshape(d_mml___,[n_m_max*n_m_max,n_l_max])*UZ__,[n_m_max,n_m_max,n_UZ_rank]);
d_mmlb____(:,:,:,1+nbeta) = d_mml___;
d_mmzb____(:,:,:,1+nbeta) = d_mmz___;
end;%for nbeta=0:n_beta-1;
t_out = toc(t_0in);
if (verbose>1); disp(sprintf(' %% calculate d_: t %0.6f',t_out)); end;

fname_mat = sprintf('%s_mat/pm_vfig_timing_FIGB__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));

%%%%%%%%;
% Now time. ;
%%%%%%%%;
t_btob_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
t_sumn_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
t_abab_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
t_sumz_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
t_fft2_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
n_mult_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
t_sumn_b_ = zeros(n_fname_b_k_Y,1);
t_abab_b_ = zeros(n_fname_b_k_Y,1);
t_sumz_b_ = zeros(n_fname_b_k_Y,1);
t_fft2_b_ = zeros(n_fname_b_k_Y,1);
n_mult_b_ = zeros(n_fname_b_k_Y,1);
n_mult_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
l2er_X_X_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
corr_X_X_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
prct_X_X_zb__ = zeros(n_UZ_rank,n_fname_b_k_Y);
%%%%%%%%;
a_k_Y_mlk_ = a_k_Y_quad_;
a_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,a_k_Y_mlk_);
a_k_Y_mkl___ = permute(a_k_Y_mlk___,[1,3,2]);
a_UX_Y_mln___ = reshape(reshape(a_k_Y_mlk___,[n_m_max*n_l_max,n_k_p_r])*diag(X_weight_r_)*UX__,[n_m_max,n_l_max,n_UX_rank]);
a_UX_Y_mnl___ = permute(a_UX_Y_mln___,[1,3,2]);
for nfname_b_k_Y=0:n_fname_b_k_Y-1;
if (verbose>0); disp(sprintf(' %% nfname_b_k_Y %d/%d',nfname_b_k_Y,n_fname_b_k_Y)); end;
b_k_Y_mlk_ = b_k_Y_yi__(:,1+nfname_b_k_Y);
b_k_Y_mlk___ = convert_spharm_to_spharm___0(n_k_p_r,l_max_,b_k_Y_mlk_);
b_k_Y_mkl___ = permute(b_k_Y_mlk___,[1,3,2]);
[ ...
 X0_mmb___ ...
 t_sumn_0 ...
 t_abab_0 ...
 t_sumz_0 ...
 t_fft2_0 ...
 n_mult_0 ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 verbose ...
,n_m_max ...
,n_k_p_r ...
,n_l_max ...
,weight_3d_k_p_r_ ...
,a_k_Y_mkl___ ...
,b_k_Y_mkl___ ...
,[] ...
,[] ...
,n_beta ...
,d_mmlb____ ...
);
X0_mmb___ = real(X0_mmb___);
t_sumn_b_(1+nfname_b_k_Y) = t_sumn_0;
t_abab_b_(1+nfname_b_k_Y) = t_abab_0;
t_sumz_b_(1+nfname_b_k_Y) = t_sumz_0;
t_fft2_b_(1+nfname_b_k_Y) = t_fft2_0;
n_mult_b_(1+nfname_b_k_Y) = n_mult_0;
%%%%%%%%;
drank = 1;
for nUZ_rank=n_UZ_rank-1:-drank:0;
if (verbose>1); disp(sprintf(' %% nUZ_rank %d/%d',nUZ_rank,n_UZ_rank)); end;
nUX_rank = min(n_UX_rank-1,nUZ_rank);
%%%%%%%%;
% compress b_. ;
%%%%%%%%;
b_UX_Y_mln___ = zeros(n_m_max,n_l_max,1+nUX_rank);
tmp_t = tic();
b_UX_Y_mln___ = reshape(reshape(b_k_Y_mlk___,[n_m_max*n_l_max,n_k_p_r])*diag(X_weight_r_)*UX__(:,1:1+nUX_rank),[n_m_max,n_l_max,1+nUX_rank]);
b_UX_Y_mnl___ = permute(b_UX_Y_mln___,[1,3,2]);
tmp_t = toc(tmp_t);
if (verbose>1); disp(sprintf(' %% b_UX_Y_mln___: %0.3fs',tmp_t)); end;
t_btob_1 = tmp_t;
%%%%%%%%;
[ ...
 X1_mmb___ ...
 t_sumn_1 ...
 t_abab_1 ...
 t_sumz_1 ...
 t_fft2_1 ...
 n_mult_1 ...
] = ...
register_spharm_to_spharm_single_beta_3_stripped_0( ...
 0*verbose ...
,n_m_max ...
,1+nUX_rank ...
,n_l_max ...
,[] ...
,a_UX_Y_mnl___ ...
,b_UX_Y_mnl___ ...
,1+nUZ_rank ...
,UZ__ ...
,n_beta ...
,d_mmzb____ ...
);
X1_mmb___ = real(X1_mmb___);
t_btob_zb__(1+nUZ_rank,1+nfname_b_k_Y) = t_btob_1;
t_sumn_zb__(1+nUZ_rank,1+nfname_b_k_Y) = t_sumn_1;
t_abab_zb__(1+nUZ_rank,1+nfname_b_k_Y) = t_abab_1;
t_sumz_zb__(1+nUZ_rank,1+nfname_b_k_Y) = t_sumz_1;
t_fft2_zb__(1+nUZ_rank,1+nfname_b_k_Y) = t_fft2_1;
n_mult_zb__(1+nUZ_rank,1+nfname_b_k_Y) = n_mult_1;
%%%%;
l2er_X_X = fnorm(X0_mmb___ - X1_mmb___)/fnorm(X0_mmb___);
corr_X_X = corr(X0_mmb___(:),X1_mmb___(:));
[~,tmp_index_nmmb] = max(X1_mmb___,[],'all','linear'); tmp_index_nmmb = tmp_index_nmmb-1;
tmp_X0 = X0_mmb___(1+tmp_index_nmmb);
prct_X_X = sum(X0_mmb___(:)>tmp_X0)/numel(X0_mmb___);
%%%%;
l2er_X_X_zb__(1+nUZ_rank,1+nfname_b_k_Y) = l2er_X_X;
corr_X_X_zb__(1+nUZ_rank,1+nfname_b_k_Y) = corr_X_X;
prct_X_X_zb__(1+nUZ_rank,1+nfname_b_k_Y) = prct_X_X;
%%%%%%%%;
end;%for nUZ_rank=n_UZ_rank-1:-1:0;
end;%for nfname_b_k_Y=0:n_fname_b_k_Y-1;
%%%%%%%%;
save(fname_mat ...
     ,'n_UX_rank','n_UZ_rank' ...
     ,'tot_t_X_kk','tot_t_UX_a','tot_t_Z_ll','tot_t_UZ' ...
     ,'t_sumn_b_','t_abab_b_','t_sumz_b_','t_fft2_b_','n_mult_b_' ...
     ,'t_btob_zb__' ...
     ,'t_sumn_zb__' ...
     ,'t_abab_zb__' ...
     ,'t_sumz_zb__' ...
     ,'t_fft2_zb__' ...
     ,'n_mult_zb__' ...
     ,'l2er_X_X_zb__' ...
     ,'corr_X_X_zb__' ...
     ,'prct_X_X_zb__' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

t_pre0 = tot_t_X_kk + tot_t_UX_a + tot_t_Z_ll + tot_t_UZ;
t_comp_z_ = sum(t_btob_zb__ + t_sumn_zb__ + t_sumz_zb__ + t_fft2_zb__,2);
t_comp_0 = sum(t_sumn_b_ + t_abab_b_ + t_sumz_b_ + t_fft2_b_);
f_comp_z_ = (8*t_comp_0)./(t_pre0 + 8*t_comp_z_);
n_mult_0 = sum(n_mult_b_);
n_mult_z_ = sum(n_mult_zb__,2);
f_mult_z_ = (n_mult_0)./n_mult_z_;
  
l2er_X_X_z_ = mean(l2er_X_X_zb__,2);
corr_X_X_z_ = mean(corr_X_X_zb__,2);
prct_X_X_z_ = mean(prct_X_X_zb__,2);

fname_fig = sprintf('%s_jpg/pm_vfig_timing_FIGB_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+256,512]);
tmp_ij_ = 1:n_UZ_rank;
subplot(1,3,1);
plot(tmp_ij_,l2er_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','k');
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0.0,1.0]); ylabel('error'); set(gca,'YTick',0.0:0.05:1.0); grid on;
legend({'frob'},'Location','East');
title('relative error');
%%%%;
subplot(1,3,2);
hold on;
plot(tmp_ij_,0+corr_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,1-prct_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','c');
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'},'Location','East');
title('correlation');
%%%%;
subplot(1,3,3);
hold on;
plot(tmp_ij_,f_comp_z_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,f_mult_z_(tmp_ij_),'ko-','MarkerFaceColor','c');
plot(tmp_ij_,ones(1,numel(tmp_ij_)),'k-','LineWidth',0.5);
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0,16]); ylabel('factor'); set(gca,'YTick',0:1:16); grid on;
legend({'time','ops'});
title('speedup');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_vfig_timing_FIGC_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024+256,512]);
tmp_ij_ = 1:14;
subplot(1,3,1);
plot(tmp_ij_,l2er_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','k');
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',1:2:max(tmp_ij_)); xlabel('rank $H$','Interpreter','latex');
ylim([0.0,1.0]); ylabel('error'); set(gca,'YTick',0.0:0.05:1.0); grid on;
legend({'frob'},'Location','East');
title('relative error');
%%%%;
subplot(1,3,2);
hold on;
plot(tmp_ij_,0+corr_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,1-prct_X_X_z_(tmp_ij_),'ko-','MarkerFaceColor','c');
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',1:2:max(tmp_ij_)); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'},'Location','East');
title('correlation');
%%%%;
subplot(1,3,3);
hold on;
plot(tmp_ij_,f_comp_z_(tmp_ij_),'ko-','MarkerFaceColor','r');
plot(tmp_ij_,f_mult_z_(tmp_ij_),'ko-','MarkerFaceColor','c');
plot(tmp_ij_,ones(1,numel(tmp_ij_)),'k-','LineWidth',0.5);
hold off;
xlim([0,max(tmp_ij_)+1]); set(gca,'XTick',1:2:max(tmp_ij_)); xlabel('rank $H$','Interpreter','latex');
ylim([0,16]); ylabel('factor'); set(gca,'YTick',0:1:16); grid on;
legend({'time','ops'});
title('speedup');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_vfig_spectrum_FIGD_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
figure(1+nf);nf=nf+1;figmed;
n_svd = 16;
subplot(1,2,1);
plot(SX_,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
set(gca,'XTick',1:n_svd); xlim([0.5,n_svd+0.5]);
xlabel('rank');
ylabel('singular-value','Interpreter','latex');
grid on;
title(sprintf('$C(\\vec{u})$ for the first %d ranks',n_svd),'Interpreter','latex');
subplot(1,2,2);
plot(SZ_,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
set(gca,'XTick',1:n_svd); xlim([0.5,n_svd+0.5]);
xlabel('rank');
ylabel('singular-value','Interpreter','latex');
grid on;
title(sprintf('$D(\\vec{v})$ for the first %d ranks',n_svd),'Interpreter','latex');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;











