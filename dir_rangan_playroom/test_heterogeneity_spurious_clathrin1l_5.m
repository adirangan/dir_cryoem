%%%%%%%%;
% tests heterogeneity_spurious on 6sct_one_leg (clathrin). ;
% Adding calculation of bayesian likelihood. ;
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

flag_recalc = 0;
flag_replot = 0;
flag_center = 0;
tolerance_master = 1e-2;
nf=0;

dir_data = sprintf('/%s/rangan/dir_cryoem/dir_clathrin_montage_Micrograph',string_root);
dir_pm = sprintf('%s/dir_pm',dir_data);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
Pixel_Spacing = 1.00; %<-- in angstroms, unsure. ;
fname_nopath_volume = '6sct_one_leg_bf60_center.mrc';
flag_center_volume = 0;

%%%%%%%%;
% First create consensus volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_base_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_emd = sprintf('%s/%s',dir_data,fname_nopath_volume);
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
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); isosurface_f_x_c_0(a_x_u_pack_,98.5); title('a packed');
subplot(2,2,2); isosurface_f_x_c_0(a_x_u_pack_,[97.5,98.5,99.5]); title('a packed');
subplot(2,2,3); isosurface_f_x_c_0(b_rho_x_u_pack_,98.5); title('b packed');
subplot(2,2,4); isosurface_f_x_c_0(b_rho_x_u_pack_,[97.5,98.5,99.5]); title('b packed');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
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

flag_plot=1;
if flag_plot;
%%%%%%%%;
% simple visualization of a_base. ;
%%%%%%%%
prctile_ = [94.0:0.5:99.5]; n_prctile = numel(prctile_);
figure(1);clf;figbig;
p_row = 3; p_col = 4; ns=0;
for nprctile=0:n_prctile-1;
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_p = prctile_(1+nprctile);
isosurface_f_x_c_0(a_x_u_base_,tmp_p);
title(sprintf('a base %.1f',tmp_p));
end;%for nprctile=0:n_prctile-1;
%%%%;
fname_fig_pre = sprintf('%s_jpg/a_x_u_base_FIGA_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
close(gcf);
%%%%;
end;%if flag_plot;
%%%%%%%%;

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
fname_fig_pre = sprintf('%s_jpg/a_k_p_quad_FIGA_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1; clf; figbig;
plot(k_p_r_all_,log10(abs(a_k_p_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
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
fname_fig_pre = sprintf('%s_jpg/a_k_Y_quad_FIGA_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(1,3,1); plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,2); plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,3); plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

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
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/S_k_p__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p__(:,1+nS)),Slim_,colormap_80s);
axis equal; axisnotick; title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
% Now convert templates to S_k_q__. ;
%%%%%%%%;
S_k_q__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/S_k_q__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = std(abs(S_k_p__(:)))*2.5*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(S_k_q__(:,1+nS)),Slim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_q__(:,1+nS))/max(abs(S_k_q__(:)))),[-4,0],colormap_80s);
title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
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
% Now rotate the volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_p231_x_u_base_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_p231_x_u_base_ = permute(a_x_u_base_,[2,3,1]);
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_p231_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_p231_x_u_base_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_p231_k_p_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_p231_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_p231_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_p231_x_u_reco_ time %0.2fs',tmp_t));
%%%%;
disp(sprintf(' %% xxnufft3d3: a_p231_x_u_reco error: %0.16f',fnorm(a_p231_x_u_base_(:)-a_p231_x_u_reco_)/fnorm(a_p231_x_u_base_(:))));
disp(sprintf(' %% at this point one should ensure that a_p231_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
%%%%%%%%;
tmp_t = tic;
[a_p231_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p231_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_p231_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_p231_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_p231_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_p231_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p231_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_p231_k_Y_quad_ --> a_p231_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% xxnufft3d3: a_p231_k_p_reco error: %0.16f',fnorm(a_p231_k_p_quad_-a_p231_k_p_reco_)/fnorm(a_p231_k_p_quad_))); %<-- this should be 2-3 digits. ;
%%%%;
a_p231_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p231_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p231_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
[ ...
 S_p231_k_p__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(a_p231_k_Y_quad__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
S_p231_k_p__ = reshape(S_p231_k_p__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
save(fname_mat ...
     ,'a_p231_x_u_base_','a_p231_k_p_quad_','a_p231_x_u_reco_','a_p231_k_Y_quad_','a_p231_k_p_reco_','a_p231_k_Y_quad__','S_p231_k_p__' ...
     );
end;%if (flag_recalc | ~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% calculate norm of templates. ;
%%%%%%%%;
S_l2_S_ = zeros(n_S,1);
S_l2_S_ = sum(abs(S_k_p__).^2 .* weight_2d_k_all_,1);
[~,tmp_ij] = min(S_l2_S_); tmp_nS = tmp_ij-1;
S_k_p_ = S_k_p__(:,1+tmp_nS);
S_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
%%%%%%%%;
S_p231_l2_S_ = zeros(n_S,1);
S_p231_l2_S_ = sum(abs(S_p231_k_p__).^2 .* weight_2d_k_all_,1);
[~,tmp_ij] = min(S_p231_l2_S_); tmp_nS = tmp_ij-1;
S_p231_k_p_ = S_p231_k_p__(:,1+tmp_nS);
S_p231_x_c_ = interp_k_p_to_x_c_xxnufft(n_x_u_pack,diameter_x_c,n_x_u_pack,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_p231_k_p_.*weight_2d_k_all_*(2*pi)^2)*sqrt(n_x_u_pack^2) * n_w_sum;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/a_p231_x_u_base_FIGB_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
p_row = 2; p_col = 3; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_parameter = struct('type','parameter','percent_threshold_',99.5);
isosurface_f_x_u_1(tmp_parameter,a_x_u_base_); title('a_x_u_base_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,S_l2_S_ ... 
,[] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;
title('S_l2_S_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,fliplr(transpose(fliplr(real(S_x_c_)))),[],colormap_beach);
axisnotick; axis image;
title('S_x_c_ (min)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_parameter = struct('type','parameter','percent_threshold_',99.5);
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_base_); title('a_p231_x_u_base_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,S_p231_l2_S_ ... 
,[] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;
title('S_p231_l2_S_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_u_pack,x_u_0_,n_x_u_pack,x_u_1_,fliplr(transpose(fliplr(real(S_p231_x_c_)))),[],colormap_beach);
axisnotick; axis image;
title('S_p231_x_c_ (min)','Interpreter','none');
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

n_M = 1024;
weight_2d_wk_ = weight_2d_k_all_;
delta_sigma = 0.1;
delta_r_max_factor = 0.25;
delta_sigma_use = delta_sigma;
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max = delta_r_max_use;
tolerance_pm = tolerance_master;
svd_eps = min(1e-2,tolerance_master);
n_delta_v_requested = 24;

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
DefocusU_CTF = 5000; %<-- see kexin email: defocus 1 = 5000A ;
DefocusV_CTF = 5000; %<-- see kexin email: defocus 1 = 5000A ;
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
clear tmp_k_c_1 tmp_k_c_2 ;
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
clear tmp_CTF_avg_k_p_wk_ tmp_CTF_avg_k_p_r_k_ ;
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/CTF_k_p_xcor__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figmed;
figbeach();
subplot(1,2,1); imagesc(CTF_k_p_r_xavg__); axis image; axisnotick; title('CTF_k_p_r_xavg__','Interpreter','none');
subplot(1,2,2); imagesc(CTF_k_p_r_xcor__); axis image; axisnotick; title('CTF_k_p_r_xcor__','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
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
fname_fig_pre = sprintf('%s_jpg/CTF_k_p_r_xxxx__',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
figure(1+nf);nf=nf+1;clf;figbig;
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
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
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
% Now calculate average CTFs for each cluster. ;
%%%%%%%%;
CTF_k_p_r_xavg_kc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_k_p_r_xavg_kc__(:,1+ncluster) = CTF_k_p_r_xavg_k_;
end;%for ncluster=0:n_cluster-1;

%%%%%%%%;
% Recalculate templates using unit norm volume. ;
%%%%%%%%;
[ ...
 X0_quad ...
,C0_quad ...
] = ...
register_spharm_to_spharm_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_quad_ ...
,a_p231_k_Y_quad_ ...
);
%%%%%%%%;
a_p231_k_Y_norm_yk_ = a_p231_k_Y_quad_/max(1e-12,sqrt(X0_quad));
a_p231_k_Y_norm_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p231_k_Y_norm_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p231_k_Y_norm_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
[ ...
 X0_norm ...
,C0_norm ...
] = ...
register_spharm_to_spharm_3( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_yk_ ...
,a_p231_k_Y_norm_yk_ ...
);
%%%%%%%%;

%%%%%%%%;
% Now calculate idealized principal-modes for each nCTF. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
a_p231_k_Y_norm_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p231_k_Y_norm_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p231_k_Y_norm_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic();
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
CTF_k_p_r_xavg_kk__ = CTF_k_p_r_xavg_k_*transpose(CTF_k_p_r_xavg_k_);
tmp_delta_sigma = delta_r_max;
if (tmp_delta_sigma> 0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_p231_k_Y_norm_yk_,CTF_k_p_r_xavg_kk__,tmp_delta_sigma);
end;%if (tmp_delta_sigma> 0);
if (tmp_delta_sigma==0);
[X_2d_xavg_dx_kk__,X_2d_xavg_dx_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_p231_k_Y_norm_yk__,CTF_k_p_r_xavg_kk__);
end;%if (tmp_delta_sigma==0);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_2d_xavg_dx_kkc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
tmp_X__ = X_2d_xavg_dx_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(tmp_X__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'X__',tmp_t);
%%%%%%%%;

pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
if (verbose);
for ncluster=0:n_cluster-1;
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
disp(sprintf(' %% ncluster %.2d/%.2d --> pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,pm_n_UX_rank_max));
end;%for ncluster=0:n_cluster-1;
end;%if (verbose);

FTK = [];
if isempty(FTK);
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK);
assert(FTK.svd_d_max>=delta_r_max);
assert(FTK.n_delta_v>=n_delta_v_requested);

FTK_0 = [];
if isempty(FTK_0);
tmp_t = tic();
FTK_0 = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,0,svd_eps,0);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK_0: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_FTK_1',tmp_t);
end;%if isempty(FTK_0);
assert(FTK_0.svd_d_max==0);
assert(FTK_0.n_delta_v==1);

%%%%%%%%;
% Now estimate distances between templates. ;
% Use one of the CTF clusters for now. ;
%%%%%%%%;
S_p231_k_p_wkS__ = S_p231_k_p__;
S_p231_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_p231_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_p231_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_p_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_q_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_t = tic();
UX_T_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_l2_dS__: %0.3fs',tmp_t)); end;
if (verbose); disp(sprintf(' %% average l2-norm of templates: %0.16f',mean(UX_T_l2_dS__(:))/(pi*k_p_r_max^2))); end;
tmp_TT_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_wkS__(:,1+nS),T_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
tmp_index = efind((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
UX_T_l2_S_ = transpose(UX_T_l2_dS__(1+tmp_index,:));
if (verbose); disp(sprintf(' %% tmp_TT_S_ vs UX_T_l2_S_: %0.16f',fnorm(tmp_TT_S_ - UX_T_l2_S_)/fnorm(tmp_TT_S_))); end;
flag_plot=0;
if flag_plot;
tmp_index = efind((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
subplot(1,2,1); hold on; 
plot(0:n_S-1,UX_T_l2_S_/(pi*k_p_r_max^2),'rx'); xlabel('nS'); ylabel('l2');
plot(0:n_S-1,tmp_TT_S_/(pi*k_p_r_max^2),'bo'); xlabel('nS'); ylabel('l2');
hold off;
subplot(1,2,2); plot(UX_T_l2_S_/(pi*k_p_r_max^2),tmp_TT_S_/(pi*k_p_r_max^2),'g.'); xlabel('l2_A'); ylabel('l2_B');
end;%if flag_plot;
%%%%%%%%;
tmp_t = tic();
[UX_T_k_q_wnS___,UX_T_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_k_p_wnS__ = reshape(UX_T_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_k_q_wnS__ = reshape(UX_T_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Visualize: ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
figure(3);clf;figbig;fig80s;
subplot(2,2,1);imagesc(reshape(permute(log10(abs(svd_VUXT_lwnS____)),[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]));axisnotick; colorbar;
subplot(2,2,2);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(abs(svd_VUXT_lwnS____).^2,[1,3,4,2]),[FTK.n_svd_l*n_UX_rank*n_S,n_w_max]),1))));
subplot(2,2,3);imagesc(reshape(permute(reshape(log10(abs(UX_T_k_q_wnS__)),[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]));axisnotick;colorbar;
subplot(2,2,4);plot(0:n_w_max-1,log10(sqrt(mean(reshape(permute(reshape(abs(UX_T_k_q_wnS__).^2,[n_w_max,n_UX_rank,n_S]),[2,3,1]),[n_UX_rank*n_S,n_w_max]),1))));
end;%if flag_plot;
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
tmp_t = tic();
TT_k_q_ = svd(UX_T_k_q_wnS__);
n_T_rank = min(find(TT_k_q_/TT_k_q_(1)<1e-3)); if isempty(n_T_rank); n_T_rank = min(size(UX_T_k_q_wnS__)); end;
[UT_k_q__,ST_k_q__,VT_k_q__] = svds(UX_T_k_q_wnS__,n_T_rank);
if (verbose>1); disp(sprintf(' %% n_S %d --> n_T_rank %d',n_S,n_T_rank)); end;
%%%%%%%%;
tmp_index = efind((FTK.delta_x_.^2 + FTK.delta_y_.^2)<1e-6);
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_TT__ ...
,delta_x_TT__ ...
,delta_y_TT__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_k_q_wnS__ ...
,UX_T_l2_S_ ...
,n_S ...
,svd_VUXT_lwnS____ ...
,UX_T_l2_dS__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wTT___: %0.3fs',tmp_t)); end;
%%%%%%%%;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now load the Micrograph. ;
% Taken from test_pm_clathrin_8.m ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_clathrin_montage_Micrograph',string_root);
flag_invert = 1;
flag_center_volume = 0;
fname_prefix='clathrin_montage_Micrograph';
fname_infix='clathrin';
Pixel_Spacing=1.06;
fname_nopath_volume='6sct_one_leg_bf60_center.mrc';
fname_nopath_micrograph='simulated_clathrin_t1000A_d5000A_1.06Apix_p1_384.mrc';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now begin the distortion. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First use 0lsq to reconstruct the volume. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_reco_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_p231_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_all_ ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_reco_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_reco_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_reco_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_reco_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only templates with equatorial viewing-angles. ;
%%%%%%%%;
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_); tmp_n_S = numel(tmp_nS_);
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_equa_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_S ...
,S_p231_k_p_wkS__(:,1+tmp_nS_) ...
,[] ...
,[] ...
,viewing_polar_a_all_(1+tmp_nS_) ...
,viewing_azimu_b_all_(1+tmp_nS_) ...
,zeros(tmp_n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_equa_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_equa_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_equa_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_equa_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only templates with equatorial viewing-angles. ;
% This time we distort/squish the azimuthal viewing angles. ;
%%%%%%%%;
tmp_squi = @(azimu_b) periodize(atan2(2*sin(azimu_b),cos(azimu_b)),0,2*pi);
viewing_polar_a_mode = mode(viewing_polar_a_all_);
tmp_nS_ = efind(abs(viewing_polar_a_all_-viewing_polar_a_mode)<1e-6); n_l = numel(tmp_nS_); tmp_n_S = numel(tmp_nS_);
flag_disp=0;
if flag_disp;
tmp_azimu_b_ = sort(viewing_azimu_b_all_(1+tmp_nS_),'ascend'); 
figure(1+nf);nf=nf+1;figsml; 
plot(tmp_azimu_b_,tmp_squi(tmp_azimu_b_),'.');
xlabel('tmp_azimu_b_','Interpreter','none');
ylabel('tmp_squi(tmp_azimu_b_)','Interpreter','none');
end;%if flag_disp;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_squi_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_S ...
,S_p231_k_p_wkS__(:,1+tmp_nS_) ...
,[] ...
,[] ...
,viewing_polar_a_all_(1+tmp_nS_) ...
,tmp_squi(viewing_azimu_b_all_(1+tmp_nS_)) ...
,zeros(tmp_n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_squi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_squi_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_squi_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_squi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now repeat with a more refined set of templates. ;
%%%%%%%%;
n_S_refine = tmp_n_S * 8;
viewing_polar_a_refine_all_ = (pi/2)*ones(n_S_refine,1);
viewing_azimu_b_refine_all_ = linspace(0,2*pi,1+n_S_refine); viewing_azimu_b_refine_all_ = transpose(viewing_azimu_b_refine_all_(1:n_S_refine));
tmp_t = tic();
[ ...
 S_p231_refine_k_p_wkS__ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_p231_k_Y_norm_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
,n_S_refine ...
,viewing_azimu_b_refine_all_ ...
,viewing_polar_a_refine_all_ ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% pm_template_2: %0.6fs',tmp_t)); end;
S_p231_refine_k_p_wkS__ = reshape(S_p231_refine_k_p_wkS__,[n_w_max*n_k_p_r,n_S_refine]);
%%%%%%%%;
% Now use 0lsq to reconstruct the volume, ;
% but using only the refined templates with equatorial viewing-angles. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_refi_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S_refine ...
,S_p231_refine_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_refine_all_ ...
,viewing_azimu_b_refine_all_ ...
,zeros(n_S_refine,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_refi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_refi_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_refi_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_refi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
% This time we distort/squish the refined azimuthal viewing angles. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_resq_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S_refine ...
,S_p231_refine_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_refine_all_ ...
,tmp_squi(viewing_azimu_b_refine_all_) ...
,zeros(n_S_refine,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_resq_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_resq_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_resq_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_resq_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

if (verbose);
disp(sprintf(' %% a_p231_k_Y_norm_yk_ vs a_p231_k_Y_norm_reco_yk_: %0.16f',fnorm(a_p231_k_Y_norm_yk_ - a_p231_k_Y_norm_reco_yk_)/fnorm(a_p231_k_Y_norm_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_yk_ vs a_p231_k_Y_norm_equa_yk_: %0.16f',fnorm(a_p231_k_Y_norm_yk_ - a_p231_k_Y_norm_equa_yk_)/fnorm(a_p231_k_Y_norm_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_yk_ vs a_p231_k_Y_norm_refi_yk_: %0.16f',fnorm(a_p231_k_Y_norm_yk_ - a_p231_k_Y_norm_refi_yk_)/fnorm(a_p231_k_Y_norm_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_reco_yk_ vs a_p231_k_Y_norm_equa_yk_: %0.16f',fnorm(a_p231_k_Y_norm_reco_yk_ - a_p231_k_Y_norm_equa_yk_)/fnorm(a_p231_k_Y_norm_reco_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_equa_yk_ vs a_p231_k_Y_norm_squi_yk_: %0.16f',fnorm(a_p231_k_Y_norm_equa_yk_ - a_p231_k_Y_norm_squi_yk_)/fnorm(a_p231_k_Y_norm_equa_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_equa_yk_ vs a_p231_k_Y_norm_resq_yk_: %0.16f',fnorm(a_p231_k_Y_norm_equa_yk_ - a_p231_k_Y_norm_resq_yk_)/fnorm(a_p231_k_Y_norm_equa_yk_)));
disp(sprintf(' %% a_p231_k_Y_norm_squi_yk_ vs a_p231_k_Y_norm_resq_yk_: %0.16f',fnorm(a_p231_k_Y_norm_squi_yk_ - a_p231_k_Y_norm_resq_yk_)/fnorm(a_p231_k_Y_norm_squi_yk_)));
end;%if (verbose);

flag_disp=1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
p_row = 2; p_col = 3; np=0;
tmp_parameter = struct('type','parameter','percent_threshold_',99.75);
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_reco_xxx_);
title('a_p231_x_u_norm_reco_xxx_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_equa_xxx_);
title('a_p231_x_u_norm_equa_xxx_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_squi_xxx_);
title('a_p231_x_u_norm_squi_xxx_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_refi_xxx_);
title('a_p231_x_u_norm_refi_xxx_','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_resq_xxx_);
title('a_p231_x_u_norm_resq_xxx_','Interpreter','none');
close(gcf);
end;%if flag_disp;

%%%%%%%%;
% Now use the templates with polar viewing-angles to form a new volume. ;
% This uses template_k_c_2__ from get_template_1. ;
%%%%%%%%;
equ_theta = pi*0.25 ; %<-- some band around the equator. ;
cap_theta = equ_theta; %<-- cap_theta = pi*0.25;
cap_sigma = cap_theta/3; %<-- so now 3 stds captures cap_theta. ;
tmp_index_pole_ = efind(min(viewing_polar_a_all_,pi-viewing_polar_a_all_)<=cap_theta);
tmp_index_equa_ = setdiff(0:n_S-1,tmp_index_pole_);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_equ_vs_pol_FIGA_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
view_ab_ = [00,05];
%%%%;
%subplot(1,3,1);
%plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
%xlabel('k0'); ylabel(''); zlabel('k2');
%xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
%set(gca,'XTick',[],'YTick',[],'ZTick',[]);
%view(view_ab_);
%%%%;
subplot(1,3,3);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,n_S - numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_pole_) ...
,viewing_polar_a_all_(1+tmp_index_pole_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.05*[-1,+1]); ylim(1.05*[-1,+1]); zlim(1.25*[0.0,+1]); axis equal; %axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view([00,0]);%view(view_ab_);
title('zoom in on polar cap');
%%%%;
subplot(1,3,2);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,n_S - numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_pole_) ...
,viewing_polar_a_all_(1+tmp_index_pole_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view(view_ab_);
title('non-equatorial templates');
%%%%;
subplot(1,3,1);
plot_sphere_grid_0(struct('linecolor_a',0.35*[1,1,1],'linecolor_b',0.65*[1,1,1]));
imagesc_S_k_p_3d_belt_3( ...
 struct('type','parameter','k_p_r_max',0.975,'flag_fill',0,'linewidth_use',2,'c_line_use_',[1,0,1]) ...
,numel(tmp_index_equa_) ...
,viewing_azimu_b_all_(1+tmp_index_equa_) ...
,viewing_polar_a_all_(1+tmp_index_equa_) ...
);
xlabel('k0'); ylabel(''); zlabel('k2');
xlim(1.5*[-1,+1]); ylim(1.5*[-1,+1]); zlim(1.5*[-1,+1]); axis vis3d;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
view(view_ab_);
title('equatorial templates');
%%%%;
sgtitle('template [$\alpha$,$\beta$] in cyan, data in magenta','Interpreter','latex');
set(gcf,'Position',1+[0,0,1024*1.5,512]);
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;

%%%%%%%%;
% Now modify the polar cap. ;
%%%%%%%%;
polar_cap_dilated_amplitude = 0.125;
template_k_r01_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2);
template_k_r012_wkS__ = sqrt(template_k_c_0__.^2 + template_k_c_1__.^2 + template_k_c_2__.^2);
template_azimu_b_wkS__ = atan2(template_k_c_1__,template_k_c_0__);
template_polar_a_wkS__ = atan2(template_k_r01_wkS__,template_k_c_2__);
S_ampl_k_p_wkS__ = S_p231_k_p_wkS__;
g_cap_wkS__ = exp(-template_polar_a_wkS__.^2/(2*cap_sigma^2));
ind_0in_wkS__ = (abs(template_polar_a_wkS__-pi/2)<=pi/2-cap_theta); %<-- near equator. ;
ind_out_wkS__ = (abs(template_polar_a_wkS__-pi/2)> pi/2-cap_theta); %<-- near polar cap. ;
S_cap_avg = real(mean(S_p231_k_p_wkS__(find(ind_out_wkS__)))); %<-- average of (unused) polar cap values. ;
S_ampl_k_p_wkS__ = S_ampl_k_p_wkS__.*ind_0in_wkS__ + S_cap_avg.*(1+1*g_cap_wkS__).*ind_out_wkS__; %<-- modify the poles. ;
%%%%%%%%;
% Now we can actually use all the templates in S_ampl_k_p_wkS__, since we have set the polar data to be consistent. ;
%%%%%%%%;
tmp_t = tic();
tmp_n_order = 5;
a_p231_k_Y_norm_rema_yk_ = ...
cg_lsq_6( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,n_S ...
,S_ampl_k_p_wkS__ ...
,[] ...
,[] ...
,viewing_polar_a_all_ + polar_cap_dilated_amplitude*sin(2*viewing_polar_a_all_) ...
,viewing_azimu_b_all_ ...
,zeros(n_S,1) ...
,[] ...
,[] ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_k_Y_norm_refi_yk_ time %0.2fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[ ... 
  a_p231_x_u_norm_rema_xxx_ ...
] = ...
convert_spharm_to_x_c_4( ...
 k_eq_d ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_p231_k_Y_norm_rema_yk_ ...
,half_diameter_x_c ...
,n_x_u_pack ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% a_p231_x_u_norm_refi_xxx_ time %0.2fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now calculate templates. ;
%%%%%%%%;
a_p231_k_Y_norm_refi_yk_ = a_p231_k_Y_norm_refi_yk_/max(1e-12,sqrt(register_spharm_to_spharm_3(verbose,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p231_k_Y_norm_refi_yk_,a_p231_k_Y_norm_refi_yk_)));
a_p231_k_Y_norm_refi_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p231_k_Y_norm_refi_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p231_k_Y_norm_refi_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = tic();
[ ...
 S_p231_k_p_norm_refi_wkS__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_p231_k_Y_norm_refi_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
S_p231_k_p_norm_refi_wkS__ = reshape(S_p231_k_p_norm_refi_wkS__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;
a_p231_k_Y_norm_rema_yk_ = a_p231_k_Y_norm_rema_yk_/max(1e-12,sqrt(register_spharm_to_spharm_3(verbose,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_p231_k_Y_norm_rema_yk_,a_p231_k_Y_norm_rema_yk_)));
a_p231_k_Y_norm_rema_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_p231_k_Y_norm_rema_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_p231_k_Y_norm_rema_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
tmp_t = tic();
[ ...
 S_p231_k_p_norm_rema_wkS__ ...
,~ ...
,tmp_n_S ...
,~ ...
,~ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max_max ...
,n_k_p_r ...
,reshape(a_p231_k_Y_norm_rema_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
S_p231_k_p_norm_rema_wkS__ = reshape(S_p231_k_p_norm_rema_wkS__,[n_w_max*n_k_p_r,tmp_n_S]);
%%%%%%%%;

%%%%%%%%;
% Now estimate distances between templates. ;
% Use one of the CTF clusters for now. ;
%%%%%%%%;
S_p231_k_q_norm_refi_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_p231_k_q_norm_refi_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_p231_k_p_norm_refi_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
S_p231_k_q_norm_rema_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_p231_k_q_norm_rema_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_p231_k_p_norm_rema_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_refi_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_p_norm_refi_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_refi_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_q_norm_refi_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_refi_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_refi_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_refi_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_t = tic();
tmp_TT_refi_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_refi_k_p_wkS__(:,1+nS),T_refi_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_refi_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
%%%%%%%%;
tmp_t = tic();
UX_T_refi_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_refi_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_refi_l2_dS__: %0.3fs',tmp_t)); end;
tmp_t = tic();
[UX_T_refi_k_q_wnS___,UX_T_refi_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_refi_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_refi_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_refi_k_p_wnS__ = reshape(UX_T_refi_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_refi_k_q_wnS__ = reshape(UX_T_refi_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
ncluster = 0; %<-- pick one of the CTF clusters. ;
UX_kn__ = UX_knc___(:,:,1+ncluster);
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
T_rema_k_p_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_p_norm_rema_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
T_rema_k_q_wkS__ = reshape(bsxfun(@times,reshape(S_p231_k_q_norm_rema_wkS__,[n_w_max,n_k_p_r,n_S]),reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1])),[n_w_sum,n_S]);
tmp_t = tic();
svd_VUXT_rema_lwnS____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_S,T_rema_k_q_wkS__,n_UX_rank,UX_kn__,sqrt(weight_2d_k_p_r_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXT_rema_lwnS____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the templates. ;
%%%%%%%%;
tmp_t = tic();
UX_T_rema_l2_dS__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_S,n_UX_rank,svd_VUXT_rema_lwnS____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_rema_l2_dS__: %0.3fs',tmp_t)); end;
tmp_t = tic();
tmp_TT_rema_S_ = zeros(n_S,1);
for nS=0:n_S-1;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_rema_k_p_wkS__(:,1+nS),T_rema_k_p_wkS__(:,1+nS))/(2*pi);
tmp_TT_rema_S_(1+nS) = tmp_TT;
end;%for nS=0:n_S-1;
%%%%%%%%;
tmp_t = tic();
[UX_T_rema_k_q_wnS___,UX_T_rema_k_p_wnS___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,n_UX_rank,n_S,svd_VUXT_rema_lwnS____,zeros(n_S,1),zeros(n_S,1));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_T_rema_k_q_wnS___: %0.6fs',tmp_t)); end;
UX_T_rema_k_p_wnS__ = reshape(UX_T_rema_k_p_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
UX_T_rema_k_q_wnS__ = reshape(UX_T_rema_k_q_wnS___(:,1:n_UX_rank,:),[n_w_max*n_UX_rank,n_S]);
%%%%%%%%;
% Calculate ampmh_X_wSM___8. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z = 1;
parameter.flag_compute_I_value = 0;
parameter.tolerance_master = tolerance_master;
parameter.pm_n_UX_rank_use = n_UX_rank;
tmp_t = tic();
[ ...
 parameter ...
,X_T_rema_T_refi__ ...
,delta_x_T_rema_T_refi__ ...
,delta_y_T_rema_T_refi__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,n_UX_rank ...
,n_S ...
,UX_T_rema_k_q_wnS__ ...
,tmp_TT_rema_S_ ...
,n_S ...
,svd_VUXT_refi_lwnS____ ...
,UX_T_refi_l2_dS__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_wTT___: %0.3fs',tmp_t)); end;
%%%%%%%%;

sigma_bayesian_ = transpose([0,2.^[-12:0.125:0]]); n_sigma_bayesian = numel(sigma_bayesian_);
[ ...
 ~ ...
,ssnll_rema_refi_Ms__ ...
,D2_rema_refi_SM__ ...
] = ...
ssnll_from_X_0( ...
 [] ...
,n_S ...
,n_S ...
,X_T_rema_T_refi__ ...
,tmp_TT_rema_S_*tmp_TT_refi_S_(1+0)/max(1e-12,tmp_TT_rema_S_(1+0)) ...
,tmp_TT_refi_S_ ...
,viewing_weight_all_ ...
,n_sigma_bayesian ...
,sigma_bayesian_ ...
);
%%%%;

%%%%%%%%;
% 20231104: fix later. ;
% unsure how to scale the template-magnitudes to calculate ssnll. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1);cla;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,ssnll_rema_refi_Ms__(:,1+0) ... 
,[] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;
subplot(1,2,2);cla;
flag_2d_vs_3d = 0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_all_ ... 
,viewing_azimu_b_all_ ... 
,max(X_T_rema_T_refi__,[],2) ... 
,[0.90,1.00] ... 
,colormap_81s ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
set(gca,'XTick',[],'YTick',[],'ZTick',[]); xlabel('x0'); ylabel('x1'); zlabel('x2'); axis equal; axis vis3d;

flag_disp=1;
if flag_disp;
prct_ = [99.55,99.65,99.75,99.85,99.95]; n_prct = numel(prct_);
for nprct=0:n_prct-1;
tmp_prct = prct_(1+nprct);
tmp_parameter = struct('type','parameter','percent_threshold_',[tmp_prct]);
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_p%.4d_FIGA',dir_pm,round(100*tmp_prct));
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
figure(1+nf);nf=nf+1;clf;figmed;fontsize_use = 12;
p_row = 1; p_col = 2; pcol=0;
subplot(p_row,p_col,1+pcol);pcol=pcol+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_refi_xxx_);
%title(sprintf('a_p231_x_u_norm_refi_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
title(sprintf('A: %.1f%%',tmp_prct),'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+pcol);pcol=pcol+1;
isosurface_f_x_u_1(tmp_parameter,a_p231_x_u_norm_rema_xxx_);
%title(sprintf('a_p231_x_u_norm_rema_xxx_: %.1f%%',tmp_prct),'Interpreter','none');
title(sprintf('B: %.1f%%',tmp_prct),'Interpreter','none'); set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1,512]); drawnow();
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
%close(gcf);
end;%if (~exist(fname_fig_jpg,'file'));
if ( exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if ( exist(fname_fig_jpg,'file'));
%%%%%%%%;
end;%for nprct=0:n_prct-1;
end;%if flag_disp;
%%%%%%%%



