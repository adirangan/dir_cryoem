%%%%%%%%;
clear;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

flag_recalc = 0;
flag_replot = 0;
flag_center_volume = 1;
flag_center_image = 1;
flag_invert = 0;
tolerance_master = 1e-2;
nf=0;

%%%%%%%%;
% parameters for synthetic data. ;
%%%%%%%%;
n_x_u_0in = 64; %<-- increase to increase number of pixels (also decrease Pixel_A_0in). ;
delta_sigma_true_0in = 0.00; %<-- increase to increase translation-spread of images. ;
snr_per_AA_0in = 0.0001; %<-- decrease to increase noise. ;
n_M_0in = 1024; %<-- number of images. ;
n_CTF_0in = 1; %<-- number of distinct CTFs. ;
Pixel_A_0in = 4.0; %<-- Pixel_Spacing in Angstroms, if n_x_u increases, this should decrease. ;
Voltage_0in = 300; %<-- no need to change. ;
Defocus_0in = 2.5e4; %<-- increase to increase number of CTF zero-crossings. ;
Defocus_relative_spread_0in = 0.50; %<-- increase to increase n_CTF_rank. ;
Amplitude_Contrast_0in  = 0.1; %<-- no need to change. ;
%%%%%%%%;
parameter = struct('type','parameter');
[ ...
 parameter ...
,fname_xfix ...
] = ...
spharm_to_mrcs_xfix_0( ...
 parameter ...
,n_x_u_0in ...
,delta_sigma_true_0in ...
,snr_per_AA_0in ...
,n_M_0in ...
,n_CTF_0in ...
,Pixel_A_0in ...
,Voltage_0in ...
,Defocus_0in ...
,Defocus_relative_spread_0in ...
,Amplitude_Contrast_0in ...
);
%%%%%%%%;
fname_prefix = 'FINYU';
fname_prefix_xfix = sprintf('%s_%s',fname_prefix,fname_xfix);

dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_pm',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_relion',string_root,fname_prefix_xfix);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_%s',string_root,fname_prefix_xfix);
Pixel_Spacing = Pixel_A_0in; %<-- in angstroms, just an estimation. ;
fname_nopath_volume = sprintf('%s.mat',fname_prefix_xfix); %<-- synthetic. ;
%%%%%%%%;
% all classes and subclasses. ;
%%%%%%%%;
fname_nopath_volume_ = ...
{ ...
  sprintf('%s.mat',fname_prefix_xfix) ... %<-- class 0. ;
};
n_volume = numel(fname_nopath_volume_);
flag_het = 0; if (n_volume> 1); flag_het = 1; end;
fname_nopath_star = sprintf('%s.star',fname_prefix_xfix);

fname_mat = sprintf('%s/%s',dir_data_star,fname_nopath_volume);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_x_u = n_x_u_0in;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
a_x_u___ = zeros(n_x_u,n_x_u,n_x_u);
scale_blockletter_ = [3,3,3];
tmp_blockletter___ = blockletter_0('',scale_blockletter_);
size_blockletter_ = size(tmp_blockletter___);
z_gap = size_blockletter_(1+2);
y_gap = floor(size_blockletter_(1+1)/4);
x_gap = floor(size_blockletter_(1+0)/8);
nz = floor(n_x_u/2 - size_blockletter_(1+2)*2.5 - z_gap*2);
ny = floor(n_x_u/2 - size_blockletter_(1+1)*0.5 - y_gap*2.5);
nx = floor(n_x_u/2 - size_blockletter_(1+0)*0.5 - x_gap*2.5);
a_x_u___( ...
 1+nx+[0:size_blockletter_(1+0)-1] ...
,1+ny+[0:size_blockletter_(1+1)-1] ...
,1+nz+[0:size_blockletter_(1+2)-1] ...
	  ) = ...
blockletter_0('F',scale_blockletter_);
nz = nz + z_gap + size_blockletter_(1+2);
ny = ny + y_gap;
nx = nx + x_gap;
a_x_u___( ...
 1+nx+[0:size_blockletter_(1+0)-1] ...
,1+ny+[0:size_blockletter_(1+1)-1] ...
,1+nz+[0:size_blockletter_(1+2)-1] ...
	  ) = ...
blockletter_0('I',scale_blockletter_);
nz = nz + z_gap + size_blockletter_(1+2);
ny = ny + y_gap;
nx = nx + x_gap;
a_x_u___( ...
 1+nx+[0:size_blockletter_(1+0)-1] ...
,1+ny+[0:size_blockletter_(1+1)-1] ...
,1+nz+[0:size_blockletter_(1+2)-1] ...
	  ) = ...
blockletter_0('N',scale_blockletter_);
nz = nz + z_gap + size_blockletter_(1+2);
ny = ny + y_gap;
nx = nx + x_gap;
a_x_u___( ...
 1+nx+[0:size_blockletter_(1+0)-1] ...
,1+ny+[0:size_blockletter_(1+1)-1] ...
,1+nz+[0:size_blockletter_(1+2)-1] ...
	  ) = ...
blockletter_0('Y',scale_blockletter_);
nz = nz + z_gap + size_blockletter_(1+2);
ny = ny + y_gap;
nx = nx + x_gap;
a_x_u___( ...
 1+nx+[0:size_blockletter_(1+0)-1] ...
,1+ny+[0:size_blockletter_(1+1)-1] ...
,1+nz+[0:size_blockletter_(1+2)-1] ...
	  ) = ...
blockletter_0('U',scale_blockletter_);
nz = nz + z_gap + size_blockletter_(1+2);
ny = ny + y_gap;
nx = nx + x_gap;
a_x_u___ = flip(flip(permute(a_x_u___,[2,3,1]),3),1);
a_x_u_load_ = a_x_u___;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_load_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figsml;
isosurface_f_x_c_0(a_x_u_load_,98.5); title('a load');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
save(fname_mat,'a_x_u_load_');
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% First create consensus volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_pack_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%fname_emd = sprintf('%s/%s',dir_data_star,fname_nopath_volume); a_x_u_load_ = cast(ReadMRC(fname_emd),'double');
load(sprintf('%s/%s',dir_data_star,fname_nopath_volume),'a_x_u_load_');
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
fname_fig = sprintf('%s_jpg/a_x_u_pack_',dir_pm);
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
     ,'x_u_0_','x_u_1_','x_u_2_','x_u_0___','x_u_1___','x_u_2___','n_x_u','n_xxx_u','xxx_u_weight_','a_x_u_pack_','b_rho_x_u_pack_','a_x_u_base_' ...
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
     ,'n_k_all','n_k_all_csum_','k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_','weight_3d_k_all_','weight_shell_k_','n_k_p_r','k_p_r_','weight_3d_k_p_r_','k_c_0_all_','k_c_1_all_','k_c_2_all_','J_node_','J_weight_','J_chebfun_','J_polyval_' ...
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
% generate templates S_k_p_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/S_k_p__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
n_w_max = 2*(l_max_max+1);
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
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
save(fname_mat ...
     ,'n_w_max','template_k_eq_d','viewing_k_eq_d' ...
     ,'S_k_p__' ...
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
disp(sprintf(' %% %s found, not creating',fname_fig));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_p__',dir_pm);
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
%%%%%%%%;
% Now convert templates to S_k_q__. ;
%%%%%%%%;
S_k_q__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_q__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(S_k_q__(:,1+nS)),Slim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_q__(:,1+nS))/max(abs(S_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('S_k_p__'),'Interpreter','none');
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
 tmp_S_k_p__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_k_Y_quad__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
tmp_S_k_p__ = reshape(tmp_S_k_p__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
disp(sprintf(' %% S_k_p__ vs tmp_S_k_p__: %0.16f',fnorm(S_k_p__-tmp_S_k_p__)/fnorm(S_k_p__)));
end;%if flag_check;

%%%%%%%%;
% Now generate mrcs- and star-file if they do not exist. ;
%%%%%%%%;
fname_mrcs = sprintf('%s/%s.mrcs',dir_data_star,fname_prefix_xfix);
fname_star = sprintf('%s/%s.star',dir_data_star,fname_prefix_xfix);
fname_mat = sprintf('%s_mat/euler_form.mat',dir_pm);
flag_exist = exist(fname_mrcs,'file') & exist(fname_star,'file') & exist(fname_mat,'file');
if ~flag_exist;
disp(sprintf(' %% %s, %s and %s not found, creating',fname_mrcs,fname_star,fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
[ ...
 parameter ...
,fname_mrcs_nopath ...
,fname_star_nopath ...
,fname_xfix ...
,M_x_c___ ...
,euler_polar_a_form_ ...
,euler_azimu_b_form_ ...
,euler_gamma_z_form_ ...
,image_delta_x_form_ ...
,image_delta_y_form_ ...
] = ...
spharm_to_mrcs_0( ...
 parameter ...
,k_p_r_max ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,half_diameter_x_c ...
,n_x_u_0in ...
,delta_sigma_true_0in ...
,snr_per_AA_0in ...
,n_M_0in ...
,n_CTF_0in ...
,Pixel_A_0in ...
,Voltage_0in ...
,Defocus_0in ...
,Defocus_relative_spread_0in ...
,Amplitude_Contrast_0in ...
,dir_data_star ...
,fname_prefix ...
);
save(fname_mat ...
     ,'parameter' ...
     ,'euler_polar_a_form_' ...
     ,'euler_azimu_b_form_' ...
     ,'euler_gamma_z_form_' ...
     ,'image_delta_x_form_' ...
     ,'image_delta_y_form_' ...
     );
end;%if ~flag_exist;
%%%%%%%%;
if flag_exist;
disp(sprintf(' %% %s, %s and %s found, not creating',fname_mrcs,fname_star,fname_mat));
load(fname_mat);
end;%if flag_exist;

%%%%%%%%;
% Now load images and CTF parameters from the star-file. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/M_k_p__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Note: can run image_center_1_wrap_0.m as a diagnostic. ;
%%%%%%%%;
image_center_1_wrap_1;
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
 dir_data_star ...
,fname_nopath_star ...
,n_M ...
);
if (fnorm(Voltage_CTF_)< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_ = 300*ones(n_M,1); end;
if flag_invert; M_x_c___ = -M_x_c___; end;
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
fname_fig = sprintf('%s_jpg/M_abs_x_c_avg_',dir_pm);
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
if (verbose); disp(sprintf(' %% n_x_M_center %d',n_x_M_center)); end;
M_x_c_avg__ = zeros(n_x_M_center,n_x_M_center);
M_x_c_std__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_avg__ = zeros(n_x_M_center,n_x_M_center);
N_x_c_std__ = zeros(n_x_M_center,n_x_M_center);
for nM=0:n_M-1;
if (mod(nM,128)==0); if (verbose); disp(sprintf(' %% nM %d/%d',nM,n_M)); end; end;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
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
     ,'n_M','n_x_M_u','x_c_0_','x_c_1_','x_c_0__','x_c_1__' ...
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
fname_fig = sprintf('%s_jpg/M_all_center_FIGB',dir_pm);
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
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q__ time %0.2fs',tmp_t)); end;
tmp_t = tic();
N_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
N_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,N_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% M_k_q__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_k_q__',dir_pm);
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
fname_fig = sprintf('%s_jpg/M_k_q_norm_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
M_k_q_norm_ = sqrt(sum(abs(M_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(M_k_q_norm_/max(M_k_q_norm_)),[-4,0],colormap_80s);
title('M_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
M_k_q_norm_ = sqrt(sum(abs(M_k_q__).^2,2)); tmp_eps = M_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,M_k_q_norm_-tmp_eps)/max(M_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('M_k_q_norm_ - tmp_eps','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/N_k_q__',dir_pm);
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
fname_fig = sprintf('%s_jpg/N_k_q_norm_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
N_k_q_norm_ = sqrt(sum(abs(N_k_q__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(N_k_q_norm_/max(N_k_q_norm_)),[-4,0],colormap_80s);
title('N_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
N_k_q_norm_ = sqrt(sum(abs(N_k_q__).^2,2)); tmp_eps = N_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,N_k_q_norm_-tmp_eps)/max(N_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('N_k_q_norm_ - tmp_eps','Interpreter','none');
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
fname_fig = sprintf('%s_jpg/CTF_k_p_sample__',dir_pm);
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
fname_fig = sprintf('%s_jpg/VSCTF_Mc_scatter_FIGA',dir_pm);
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
[X_2d_Nemp_d1_kk__,X_2d_Nemp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_index_nM_from_ncluster,N_k_p__(:,1+index_nM_from_ncluster_));
X_2d_Nemp_d1_kkc___(:,:,1+ncluster) = X_2d_Nemp_d1_kk__;
X_2d_Nemp_d1_weight_rc__(:,1+ncluster) = X_2d_Nemp_d1_weight_r_;
clear X_2d_Nemp_d1_kk__ X_2d_Nemp_d1_weight_r_ ;
[X_2d_Memp_d1_kk__,X_2d_Memp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_index_nM_from_ncluster,M_k_p__(:,1+index_nM_from_ncluster_));
X_2d_Memp_d1_kkc___(:,:,1+ncluster) = X_2d_Memp_d1_kk__;
X_2d_Memp_d1_weight_rc__(:,1+ncluster) = X_2d_Memp_d1_weight_r_;
clear X_2d_Memp_d1_kk__ X_2d_Memp_d1_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_2d_Memp_d1_kkc___FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(X_2d_Memp_d1_kkc___(:,:,1+ncluster)); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
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
fname_fig = sprintf('%s_jpg/X_2d_Nemp_d1_kkc___FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(X_2d_Nemp_d1_kkc___(:,:,1+ncluster)); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
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
[X_2d_xavg_d0_kk__,X_2d_xavg_d0_weight_r_] = principled_marching_cost_matrix_6(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,tmp_CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
X_2d_xavg_d0_kkc___(:,:,1+ncluster) = X_2d_xavg_d0_kk__;
X_2d_xavg_d0_weight_rc__(:,1+ncluster) = X_2d_xavg_d0_weight_r_;
clear X_2d_xavg_d0_kk__ X_2d_xavg_d0_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_2d_xavg_d0_kkc___FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;figbeach();
ns=0;
for ncluster=0:min(20,n_cluster)-1;
subplot(4,5,1+ns);ns=ns+1;
imagesc(X_2d_xavg_d0_kkc___(:,:,1+ncluster)); axis image; axisnotick;
end;%for ncluster=0:n_cluster-1;
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
% Now use ideal principal-modes to align images to volume. ;
%%%%%%%%;
for flag_N_vs_M = 0:1;
if flag_N_vs_M==0; 
tmp_N_k_p__ = M_k_p__; tmp_str = 'M';
tmp_image_delta_x_form_ = image_delta_x_form_;
tmp_image_delta_y_form_ = image_delta_y_form_;
end;%if flag_N_vs_M==0; 
if flag_N_vs_M==1;
tmp_N_k_p__ = N_k_p__; tmp_str = 'N';
tmp_image_delta_x_form_ = image_delta_x_form_ + image_center_delta_x_c_0_;
tmp_image_delta_y_form_ = image_delta_y_form_ + image_center_delta_x_c_1_;
end;%if flag_N_vs_M==1;
fname_mat = sprintf('%s_mat/a_k_Y_reco_from_%s__.mat',dir_pm,tmp_str);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.n_iteration = 8;
parameter.delta_r_max = 0.03;
parameter.n_delta_v_requested = 24;
parameter.delta_r_upb = 0.25;
parameter.fname_align_a_k_Y_pre = sprintf('%s_mat/a_k_Y_reco_from_%s__',dir_pm,tmp_str);
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
] = ...
pm_align_M_k_p_to_a_k_Y_2( ...
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
,euler_polar_a_form_ ...
,euler_azimu_b_form_ ...
,euler_gamma_z_form_ ...
,tmp_image_delta_x_form_ ...
,tmp_image_delta_y_form_ ...
);
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
     );
clear tmp_N_k_p__;
end;%if (~exist(fname_mat,'file'));
end;%for flag_N_vs_M = 0:1;
%%%%%%%%;
flag_N_vs_M = flag_center_image;
fname_mat = sprintf('%s_mat/a_k_Y_reco_from_%s__.mat',dir_pm,tmp_str);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
tmp_ = load(fname_mat);
euler_polar_a_true_ = tmp_.euler_polar_a_Mi__(:,end);
euler_azimu_b_true_ = tmp_.euler_azimu_b_Mi__(:,end);
euler_gamma_z_true_ = tmp_.euler_gamma_z_Mi__(:,end);
image_delta_x_true_ = tmp_.image_delta_x_acc_Mi__(:,end) + tmp_.image_delta_x_upd_Mi__(:,end);
image_delta_y_true_ = tmp_.image_delta_y_acc_Mi__(:,end) + tmp_.image_delta_y_upd_Mi__(:,end);
corr_a_k_Y_true = tmp_.corr_a_k_Y_i_(end);
image_X_value_true_ = tmp_.image_X_value_Mi__(:,end-1);
clear tmp_;
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
% Check to see if M_k_p__ or N_k_p__ aligns more accurately with the consensus volume. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_vs_N_corr_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
tmp_fname_mat = sprintf('%s_mat/a_k_Y_reco_from_%s__.mat',dir_pm,'M');
tmp_M_ = load(tmp_fname_mat);
tmp_fname_mat = sprintf('%s_mat/a_k_Y_reco_from_%s__.mat',dir_pm,'N');
tmp_N_ = load(tmp_fname_mat);
figure(1+nf);nf=nf+1;figsml;
clf;ns=0;
subplot(1,2,1+ns);ns=ns+1;
plot(real(tmp_M_.corr_a_k_Y_i_),real(tmp_N_.corr_a_k_Y_i_),'o-',[0,1],[0,1],'k-'); 
xlim([0,1]);ylim([0,1]);
axis square; xlabel('M corr'); ylabel('N corr'); title('corr_a_k_Y_i_','Interpreter','none');
subplot(1,2,1+ns);ns=ns+1;
plot(tmp_M_.image_X_value_Mi__(:,end-1),tmp_N_.image_X_value_Mi__(:,end-1),'o',[0,1],[0,1],'k-');
xlim([0,1]);ylim([0,1]);
axis square; xlabel('M X_value','Interpreter','none'); ylabel('N X_value','Interpreter','none');
title('X_value','Interpreter','none');
% Note that image 1+nM=1+525 has a particularly high correlation (0.7). ;
clear tmp_M_ tmp_N_;
sgtitle(sprintf('M vs N for consensus volume'),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/euler_angle_true_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
dscale = 2.5;
delta_sigma_true = std([image_delta_x_true_;image_delta_y_true_],1);
figure(1+nf);nf=nf+1;clf; set(gcf,'Position',1+[0,0,512*3,512]);
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_(1:n_M),+euler_azimu_b_true_(1:n_M),0.35);
markersize_euler = 50;
c2d_delta__ = colormap_gaussian_2d(+image_delta_x_true_(1:n_M),+image_delta_y_true_(1:n_M),dscale*delta_sigma_true,0.35);
markersize_delta = 50;
%%%%%%%%;
subplot(1,3,1+[0,1]);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_true_(1:n_M),1*pi-euler_polar_a_true_(1:n_M),markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title(sprintf('true view (corr %0.4f)',corr_a_k_Y_true));
%%%%%%%%;
subplot(1,3,1+2);
hold on;
patch(dscale*delta_sigma_true*[-1;+1;+1;-1],dscale*delta_sigma_true*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_true_(1:n_M),+image_delta_y_true_(1:n_M),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma_true*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

if flag_center_image==1; 
% do nothing. ;
end;%if flag_center_image==1; 
if flag_center_image==0;
N_k_p__ = M_k_p__;
N_k_q__ = M_k_q__;
end;%if flag_center_image==0;

%%%%%%%%;
% Now test out principled marching. ;
% First define delta_sigma (for translations). ;
%%%%%%%%;
delta_sigma = 1.0 * std([image_delta_x_true_;image_delta_y_true_]); %<-- no reduction. ;
%%%%%%%%;
% Then measure X_TM_. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/X_TM_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
X_TM_ = zeros(n_M,1);
T_k_p__ = zeros(n_w_sum,n_M);
flag_plot=0;
if (flag_plot); 
figure(1+nf);nf=nf+1;clf;figbig;
end;%if (flag_plot); 
%%%%%%%%;
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
for nM=0:n_M-1;
tmp_euler_polar_a = +euler_polar_a_true_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_true_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
N_k_p_ = N_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,N_k_p_,N_k_p_);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,N_k_p_);
X_TM_(1+nM) = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
T_k_p__(:,1+nM) = T_k_p_;
if flag_plot;
figbeach();
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(N_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(M(k))');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(N_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(M(k))');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(T(k))');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(T(k))');
drawnow;
end;%if flag_plot;
end;%for nM=0:n_M-1;
%%%%%%%%;
clear T_k_p__ T_k_p_ ;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma' ...
     ,'X_TM_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% Now plot correlation X_TM_ vs viewing-angle and translation. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_TM_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
markersize_use = 32;
e_c2d__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_true_(1:n_M),+euler_azimu_b_true_(1:n_M),0.35);
d_c2d__ = colormap_gaussian_2d(image_delta_x_true_,image_delta_y_true_,delta_sigma,0.35);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
nc_beach_ = max(0,min(n_c_beach-1,floor(n_c_beach*X_TM_)));
n_h = 64; lh_lim_ = [0,10];
figure(1+nf);nf=nf+1;
figbeach();
figbig;
subplot(2,2,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_true_(1:n_M),1*pi-euler_polar_a_true_(1:n_M),markersize_use,e_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
subplot(2,2,2);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_true_(1:n_M),1*pi-euler_polar_a_true_(1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view (correlation)');
subplot(2,2,3);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(image_delta_x_true_(1:n_M),image_delta_y_true_(1:n_M),markersize_use,d_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
subplot(2,2,4);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(image_delta_x_true_(1:n_M),image_delta_y_true_(1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta (correlation)');
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
% Now determine best possible correlation with ground truth, ;
% using the image-subset given. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_Y_0lsq_reco_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_n_order = 5; tmp_n_M = n_M;
tmp_euler_polar_a_ = +euler_polar_a_true_(1+(0:n_M-1));
tmp_euler_azimu_b_ = +euler_azimu_b_true_(1+(0:n_M-1));
tmp_euler_gamma_z_ = +euler_gamma_z_true_(1+(0:n_M-1));
tmp_image_delta_x_ = +image_delta_x_true_(1+(0:n_M-1));
tmp_image_delta_y_ = +image_delta_y_true_(1+(0:n_M-1));
%%%%%%%%;
tmp_t = tic;
a_k_Y_0lsq_reco_ = ...
cg_lsq_4( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_M ...
,N_k_p__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% N_k_p__ --> a_k_Y_0lsq_reco_ time %0.2fs',tmp_t));
[ ...
 X_0lsq_reco ...
,X_0lsq_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,a_k_Y_0lsq_reco_ ...
);
disp(sprintf(' %% X_0lsq_reco %0.3f flag_flip %d',X_0lsq_reco,X_0lsq_flag_flip));
%%%%%%%%;
tmp_t = tic;
[a_k_p_0lsq_reco] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_0lsq_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_0lsq_reco_ --> a_k_p_0lsq_reco time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; 
a_x_u_0lsq_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_0lsq_reco.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_0lsq_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
save(fname_mat ...
     ,'a_k_Y_0lsq_reco_' ...
     ,'a_k_p_0lsq_reco' ...
     ,'a_x_u_0lsq_reco_' ...
     ,'X_0lsq_reco' ...
     ,'X_0lsq_flag_flip' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_0lsq_reco_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
isosurface_f_x_u_0(reshape(real(a_x_u_0lsq_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('a_x_u_0lsq_reco_: corr %0.4f',X_0lsq_reco),'Interpreter','none');
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
% Now determine best possible correlation with ground truth, ;
% using the image-subset given. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_k_Y_0qbp_reco_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
qbp_eps = tolerance_master; tmp_n_M = n_M;
tmp_euler_polar_a_ = +euler_polar_a_true_(1+(0:n_M-1));
tmp_euler_azimu_b_ = +euler_azimu_b_true_(1+(0:n_M-1));
tmp_euler_gamma_z_ = +euler_gamma_z_true_(1+(0:n_M-1));
tmp_image_delta_x_ = +image_delta_x_true_(1+(0:n_M-1));
tmp_image_delta_y_ = +image_delta_y_true_(1+(0:n_M-1));
%%%%%%%%;
tmp_t = tic;
a_k_Y_0qbp_reco_ = ...
qbp_6( ...
 qbp_eps ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_M ...
,N_k_p__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p_wkC__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% N_k_p__ --> a_k_Y_0qbp_reco_ time %0.2fs',tmp_t));
[ ...
 X_0qbp_reco ...
,X_0qbp_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,a_k_Y_0qbp_reco_ ...
);
disp(sprintf(' %% X_0qbp_reco %0.3f flag_flip %d',X_0qbp_reco,X_0qbp_flag_flip));
%%%%%%%%;
tmp_t = tic;
[a_k_p_0qbp_reco] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_0qbp_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_0qbp_reco_ --> a_k_p_0qbp_reco time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; 
a_x_u_0qbp_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_0qbp_reco.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_0qbp_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
save(fname_mat ...
     ,'a_k_Y_0qbp_reco_' ...
     ,'a_k_p_0qbp_reco' ...
     ,'a_x_u_0qbp_reco_' ...
     ,'X_0qbp_reco' ...
     ,'X_0qbp_flag_flip' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_0qbp_reco_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
isosurface_f_x_u_0(reshape(real(a_x_u_0qbp_reco_),n_x_u_pack,n_x_u_pack,n_x_u_pack),[90,95,99]);
title(sprintf('a_x_u_0qbp_reco_: corr %0.4f',X_0qbp_reco),'Interpreter','none');
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
date_diff_threshold = 0.25;
flag_force_create_mat=0;flag_force_create_tmp=0;
delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
tolerance_pm_ = 0.1.^[1:0.5:3];
n_tolerance_pm = numel(tolerance_pm_);
delta_r_max_factor_ = [0.00,0.25,0.50];
n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
flag_compute = 1 & ( strcmp(platform,'access1') | 0*strcmp(platform,'rusty') | 0*strcmp(platform,'eval1') );
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
%%%%%%%%;
a_k_Y_true_ = a_k_Y_quad_;
euler_polar_a_true_ = +euler_polar_a_true_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_true_(1+(0:n_M-1));
euler_gamma_z_true_ = +euler_gamma_z_true_(1+(0:n_M-1));
image_delta_x_true_ = +image_delta_x_true_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_true_(1+(0:n_M-1));
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
for flag_alternate_MS_vs_SM = [0:1];
% test with: ;
% ndat_rseed=0; dat_rseed = dat_rseed_(1+ndat_rseed); ndelta_r_max_factor = 0; delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor); delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2)); delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); ntolerance_pm=0; tolerance_pm = tolerance_pm_(1+ntolerance_pm); flag_alternate_MS_vs_SM=1;

parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.tolerance_pm = tolerance_pm;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
parameter.flag_alternate_MS_vs_SM = flag_alternate_MS_vs_SM;
parameter = ...
ampmut_wrap_wrap_5( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_max ...
,l_max_ ...
,n_M ...
,N_k_p__ ...
,n_CTF ...
,index_nCTF_from_nM_ ...
,CTF_k_p_r_kC__ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);

%%%%;
if (parameter.n_complete_calculation> 1);
if strcmp(platform,'rusty');
disp(sprintf(' %% parameter.n_complete_calculation %d, returning',parameter.n_complete_calculation));
return; %<-- halt after one calculation on rusty. ;
end;%if strcmp(platform,'rusty');
end;%if (parameter.n_complete_calculation> 1);
%%%%;
end;%for flag_alternate_MS_vs_SM = [0:1];
end;%for ntolerance_pm=0:n_tolerance_pm-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
