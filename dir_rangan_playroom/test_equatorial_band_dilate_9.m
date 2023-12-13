%%%%%%%%;
% Tests equatorial band dilation. ;
% Uses LetB1 as a case-study for longitudinal perturbation. ;
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

flag_verbose = 0;
flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

%%%%%%%%;
dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_LetB1/dir_pm',string_root);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_LetB1/dir_relion',string_root);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_LetB1',string_root);
Pixel_Spacing = 1.31; %<-- in angstroms, from https://www.ebi.ac.uk/pdbe/emdb/empiar/entry/10350/ ;
fname_nopath_volume = 'emd_20993.map';
fname_nopath_star = 'job_569_model_1.star';
%%%%;
dir_manuscript = sprintf('/%s/rangan/dir_cryoem/dir_spurious_heterogeneity_manuscript',string_root);
if ~exist(dir_manuscript,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript)); mkdir(dir_manuscript); end;
dir_manuscript_jpg = sprintf('%s/dir_M3d_shape_longitudinal_perturbation_jpg',dir_manuscript);
if ~exist(dir_manuscript_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_manuscript_jpg)); mkdir(dir_manuscript_jpg); end;
%%%%%%%%;

%%%%%%%%;
% First load volume. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/a_x_u_pack_.mat',dir_pm);
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
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = n_x_u_pack^3;
X_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
% Calculate moments. ;
%%%%%%%%;
a_rho_x_u_pack_ = a_x_u_pack_ + min(a_x_u_pack_,[],'all');
a_rho_x_c_0_avg = sum(X_u_0_.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_avg = sum(X_u_1_.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_avg = sum(X_u_2_.^1.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_0_std = sum((X_u_0_ - a_rho_x_c_0_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_1_std = sum((X_u_1_ - a_rho_x_c_1_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_2_std = sum((X_u_2_ - a_rho_x_c_2_avg).^2.*a_rho_x_u_pack_/sum(a_rho_x_u_pack_,'all'),'all');
a_rho_x_c_avg_ = [a_rho_x_c_0_avg ; a_rho_x_c_1_avg ; a_rho_x_c_2_avg];
a_rho_x_c_std_ = [a_rho_x_c_0_std ; a_rho_x_c_1_std ; a_rho_x_c_2_std];
disp(sprintf(' %% a_rho_x_c_std_ vs a_rho_x_c_avg_: %0.2f',fnorm(a_rho_x_c_std_)/fnorm(a_rho_x_c_avg_)));
if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
disp(sprintf(' %% Warning! molecule may not be well centered. Consider recentering.'));
end;%if (max(abs(a_rho_x_c_avg_))> diameter_x_c/n_x_u_pack);
%%%%%%%%;
% Possible to re-center. ;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_2_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
[K_u_0_,K_u_1_,K_u_2_] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_K_u = n_x_u_pack^3;
b_rho_x_u_pack_ = real(ifftn(fftn(a_rho_x_u_pack_).*exp(-i*2*pi*(K_u_0_*a_rho_x_c_0_avg + K_u_1_*a_rho_x_c_1_avg + K_u_2_*a_rho_x_c_2_avg))));
b_rho_x_c_0_avg = sum(X_u_0_.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_1_avg = sum(X_u_1_.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_2_avg = sum(X_u_2_.^1.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_0_std = sum((X_u_0_ - b_rho_x_c_0_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_1_std = sum((X_u_1_ - b_rho_x_c_1_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_2_std = sum((X_u_2_ - b_rho_x_c_2_avg).^2.*b_rho_x_u_pack_/sum(b_rho_x_u_pack_,'all'),'all');
b_rho_x_c_avg_ = [b_rho_x_c_0_avg ; b_rho_x_c_1_avg ; b_rho_x_c_2_avg];
b_rho_x_c_std_ = [b_rho_x_c_0_std ; b_rho_x_c_1_std ; b_rho_x_c_2_std];
disp(sprintf(' %% b_rho_x_c_std_ vs b_rho_x_c_avg_: %0.2f',fnorm(b_rho_x_c_std_)/fnorm(b_rho_x_c_avg_)));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_pack_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(1,2,1); isosurface_f_x_c_0(a_x_u_pack_,98.5); title('a packed');
subplot(1,2,2); isosurface_f_x_c_0(a_x_u_pack_,[97.5,98.5,99.5]); title('a packed');
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
     ,'x_u_0_','x_u_1_','x_u_2_','X_u_0_','X_u_1_','X_u_2_','n_x_u','n_X_u','X_u_weight_','a_x_u_pack_' ...
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
flag_verbose=0;
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
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = nufft3d3(n_X_u,X_u_0_(:)*eta,X_u_1_(:)*eta,X_u_2_(:)*eta,a_x_u_pack_(:).*X_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_pack_(:)-a_x_u_reco_)/fnorm(a_x_u_pack_(:))));
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
flag_verbose=0;
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
[a_k_Y_quad_] = convert_k_p_to_spharm_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_))); %<-- this should be 2-3 digits. ;
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
if (flag_verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_w_max %d',n_viewing_all,n_viewing_polar_a,max(n_w_))); end;
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

flag_verbose=0;
fname_mat = sprintf('%s_mat/a_k_p_uni_dbdb_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now generate uniform (i.e., non-adaptive) k_p_r_ grid on sphere. ;
%%%%%%%%;
tmp_t = tic;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
flag_unif_vs_adap = 1;
flag_tensor_vs_adap = 1;
[ ...
 n_k_uni_all ...
,n_k_uni_all_csum_ ...
,k_p_r_uni_all_ ...
,k_p_azimu_b_uni_all_ ...
,k_p_polar_a_uni_all_ ...
,weight_3d_k_uni_all_ ...
,weight_shell_k_uni_all_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,k_c_0_uni_all_ ...
,k_c_1_uni_all_ ...
,k_c_2_uni_all_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,n_polar_a_k_ ...
,polar_a_ka__ ...
,n_azimu_b_ka__ ...
] = ...
sample_sphere_7( ...
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
% Now, given a_k_Y_quad_, construct a_k_Y_dbdb_, ;
% which is the second-derivative of a_k_Y_quad_ with respect to azimu_b. ;
%%%%%%%%;
a_k_Y_dbdb_ = - Y_m_val_.^2 .* a_k_Y_quad_ ;
tmp_t = tic;
[a_k_p_dbdb_] = convert_spharm_to_k_p_1(flag_verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_dbdb_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_dbdb_ --> a_k_p_dbdb_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_p_uni_dbdb_] = convert_spharm_to_k_p_1(flag_verbose,n_k_uni_all,n_k_uni_all_csum_,k_p_r_uni_all_,k_p_azimu_b_uni_all_,k_p_polar_a_uni_all_,weight_3d_k_uni_all_,weight_shell_k_uni_all_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_dbdb_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_dbdb_ --> a_k_p_uni_dbdb_ time %0.2fs',tmp_t));
%%%%%%%%;
% Note that a_k_p_dbdb_ will need to be scaled by the inverse of: ;
% a_k_p_dbdb_fact_ = k_p_r_all_.^2 .* sin(k_p_polar_a_all_).^2 ;
%%%%%%%%;
a_k_p_dbdb_fact_ = k_p_r_all_.^2 .* sin(k_p_polar_a_all_).^2 ; %<-- factor relating a_k_p_dbdb_ to the azimuthal-component of the laplacian. ;
a_k_p_uni_dbdb_fact_ = k_p_r_uni_all_.^2 .* sin(k_p_polar_a_uni_all_).^2 ; %<-- factor relating a_k_p_uni_dbdb_ to the azimuthal-component of the laplacian. ;
%%%%%%%%;
save(fname_mat ...
,'k_p_r_max','k_eq_d' ...
,'n_k_uni_all' ...
,'n_k_uni_all_csum_' ...
,'k_p_r_uni_all_' ...
,'k_p_azimu_b_uni_all_' ...
,'k_p_polar_a_uni_all_' ...
,'weight_3d_k_uni_all_' ...
,'weight_shell_k_uni_all_' ...
,'n_k_p_r' ...
,'k_p_r_' ...
,'weight_3d_k_p_r_' ...
,'k_c_0_uni_all_' ...
,'k_c_1_uni_all_' ...
,'k_c_2_uni_all_' ...
,'n_polar_a_k_' ...
,'polar_a_ka__' ...
,'n_azimu_b_ka__' ...
,'a_k_Y_dbdb_','a_k_p_dbdb_','a_k_p_uni_dbdb_','a_k_p_dbdb_fact_','a_k_p_uni_dbdb_fact_' ...
);
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
% Now reform a_k_p_uni_dbdb_ and a_k_p_uni_dbdb_fact_ into a uniform tensor grid. ;
%%%%%%%%;
n_k_p_uni_polar_a = n_polar_a_k_(1+0);
k_p_uni_polar_a_ = polar_a_ka__{1+0};
n_k_p_uni_azimu_b = n_azimu_b_ka__{1+0}(1+0);
k_p_uni_azimu_b_ = linspace(0,2*pi,1+n_k_p_uni_azimu_b); k_p_uni_azimu_b_ = transpose(k_p_uni_azimu_b_(1:n_k_p_uni_azimu_b));
k_p_azimu_b_uni_bak___ = reshape(k_p_azimu_b_uni_all_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
k_p_polar_a_uni_bak___ = reshape(k_p_polar_a_uni_all_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
k_p_r_uni_bak___ = reshape(k_p_r_uni_all_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
weight_3d_k_uni_bak___ = reshape(weight_3d_k_uni_all_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
a_k_p_uni_dbdb_bak___ = reshape(a_k_p_uni_dbdb_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
a_k_p_uni_dbdb_fact_bak___ = reshape(a_k_p_uni_dbdb_fact_,[n_k_p_uni_azimu_b,n_k_p_uni_polar_a,n_k_p_r]);
%%%%%%%%;
% Now we can examine: ;
% plot(weight_3d_k_uni_all_.*abs(a_k_p_uni_dbdb_./max(1e-12,a_k_p_uni_dbdb_fact_)),'.') ;
% to confirm that the factor in the denominator does not produce spurious large values. ;
%%%%%%%%;

%%%%%%%%;
% Define rotations. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;

flag_verbose=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat = sprintf('%s_mat/test_hessian_0_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% In the following calculation we will focus on a longitudinal (i.e., equatorial-band) perturbation. ;
% We assume that each template is rotated (along the longitudinal-direction) by a perturbation g_azimu_b. ;
% Note that, due to our strategy for constructing templates, ;
% the longitudinal-dilation corresponds to incrementing the viewing_azimu_b by g_dilation, ;
% meaning that the template itself will remain indexed by gamma_z.
%%%%%%%%;
% Note that, to first order: ;
% f_azimu_b_inverse(azimu_b) = azimu_b - g_azimu_b(azimu_b). ;
%%%%%%%%;
% point_output is indexed by polar_a and azimu_b. ;
% point_pole is indexed by n_w_max. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now we collect the following arrays: ;
% point_output_*_ab__ : the variables corresponding to points on the surface of the sphere ;
% These are each (n_point_a,n_point_b) arrays. ;
% (note that point_output_polar_a_ begins at -pi/2) ;
% point_pole_*_abw___ : the variables corresponding to poles associated with points on the surface of the sphere. ;
% These are each (n_point_a,n_point_b,n_w_max) arrays. ;
% The point_pole_*_abw___(1+npoint_a,1+npoint_b,:) correspond to point_output_*_ab__(1+npoint_a,1+npoint_b). ;
%%%%%%%%;
n_w_max = 128; n_w_0in_ = n_w_max;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = transpose(gamma_z_(1:n_w_max)); cc_ = cos(gamma_z_); sc_ = sin(gamma_z_);
n_point_a = 1+ 64; polar_a_a_ = linspace(0,1*pi,n_point_a+2); polar_a_a_ = transpose(polar_a_a_(2:end-1)); %<-- not periodic. ;
n_point_b = 1+128; azimu_b_b_ = linspace(0,2*pi,n_point_b+0); azimu_b_b_ = transpose(azimu_b_b_(1:end-0)); %<-- yes periodic. ;
point_output_azimu_b_ab__ = zeros(n_point_a,n_point_b,1);
point_output_polar_a_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_0_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_1_ab__ = zeros(n_point_a,n_point_b,1);
point_output_k_c_2_ab__ = zeros(n_point_a,n_point_b,1);
point_pole_azimu_b_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_polar_a_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_0_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_1_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_k_c_2_abw___ = zeros(n_point_a,n_point_b,n_w_max);
point_pole_gamma_z_abw___ = zeros(n_point_a,n_point_b,n_w_max);
tmp_error_0 = 0; tmp_error_1 = 0; tmp_error_2 = 0; tmp_error_3 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for npoint_a=0:n_point_a-1; 
if (flag_verbose>-1); disp(sprintf(' %% npoint_a %d/%d ; tmp_error_ [%0.16f,%0.16f,%0.16f,%0.16f]',npoint_a,n_point_a,tmp_error_0,tmp_error_1,tmp_error_2,tmp_error_3)); end;
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
point_output_azimu_b = azimu_b_b_(1+npoint_b); %<-- yes periodic. ;
point_output_polar_a = -pi/2 + polar_a_a_(1+npoint_a); %<-- not periodic. ;
point_output_gamma_z = 0;
%point_output_k_c_ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [1;0;0];
point_output_k_c_ = [+cos(point_output_azimu_b)*cos(point_output_polar_a);+sin(point_output_azimu_b)*cos(point_output_polar_a);-sin(point_output_polar_a)] ;
point_output_k_r01 = sqrt(point_output_k_c_(1+0).^2 + point_output_k_c_(1+1).^2);
point_output_k_c_0_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+0);
point_output_k_c_1_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+1);
point_output_k_c_2_ab__(1+npoint_a,1+npoint_b) = point_output_k_c_(1+2);
point_pole_k_c__ = Rz(point_output_azimu_b) * Ry(point_output_polar_a) * Rz(point_output_gamma_z) * [zeros(1,n_w_max);transpose(sc_);transpose(cc_)];
point_pole_k_c_0_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+0,:);
point_pole_k_c_1_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+1,:);
point_pole_k_c_2_abw___(1+npoint_a,1+npoint_b,:) = point_pole_k_c__(1+2,:);
point_output_azimu_b_ab__(1+npoint_a,1+npoint_b) = point_output_azimu_b;
point_output_polar_a_ab__(1+npoint_a,1+npoint_b) = point_output_polar_a;
%%%%%%%%%%%%%%%%;
for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
% Now we step through each of the templates associated with the point point_output. ;
%%%%%%%%;
point_pole_k_c_ = point_pole_k_c__(:,1+npole);
point_pole_k_r01 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2);
point_pole_k_r012 = sqrt(point_pole_k_c_(1+0).^2 + point_pole_k_c_(1+1).^2 + point_pole_k_c_(1+2).^2);
point_pole_azimu_b = atan2(point_pole_k_c_(1+1),point_pole_k_c_(1+0));
point_pole_polar_a = atan2(point_pole_k_r01,point_pole_k_c_(1+2));
point_pole_azimu_b_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_azimu_b;
point_pole_polar_a_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_polar_a;
%%%%%%%%;
% Eventually we will determine the predilated template that is associated with the original template mapped to point_output. ;
% Here we determine point_pole_gamma_z, ;
% which is the angle within the point_pole_template_k_c__ ;
% (denoted gamma_z below) ;
% such that the point_pole_template_k_c__ evaluated at that gamma_z ;
% (denoted point_pole_template_gamma_z_k_c_) ;
% is the same as point_output_k_c_. ;
%%%%%%%%;
point_pole_template_gamma0_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a) ...
;-sin(point_pole_polar_a) ...
];
point_pole_template_sgx_k_c_ = cross(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_template_sgx = dot(point_pole_template_sgx_k_c_,point_pole_k_c_);
point_pole_template_cgx = dot(point_pole_template_gamma0_k_c_,point_output_k_c_);
point_pole_gamma_z = atan2(point_pole_template_sgx,point_pole_template_cgx);
point_pole_gamma_z_abw___(1+npoint_a,1+npoint_b,1+npole) = point_pole_gamma_z;
point_pole_template_gamma_z_k_c_ = [...
 +cos(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gamma_z) - sin(point_pole_azimu_b)*sin(point_pole_gamma_z) ...
;+sin(point_pole_azimu_b)*cos(point_pole_polar_a)*cos(point_pole_gamma_z) + cos(point_pole_azimu_b)*sin(point_pole_gamma_z) ...
;-sin(point_pole_polar_a)*cos(point_pole_gamma_z) ...
];
tmp_error_0 = tmp_error_0 + fnorm(point_pole_template_gamma_z_k_c_ - point_output_k_c_); %<-- should be 0. ;
point_pole_template_gamma_z_k_r01 = sqrt(point_pole_template_gamma_z_k_c_(1+0).^2 + point_pole_template_gamma_z_k_c_(1+1).^2);
tmp_error_1 = tmp_error_1 + fnorm((cos(point_pole_polar_a)*cos(point_pole_gamma_z))^2 + (sin(point_pole_gamma_z))^2 - point_pole_template_gamma_z_k_r01^2); %<-- should be 0. ;
point_pole_template_gamma_z_k_r012 = sqrt(point_pole_template_gamma_z_k_c_(1+0).^2 + point_pole_template_gamma_z_k_c_(1+1).^2 + point_pole_template_gamma_z_k_c_(1+2).^2);
point_pole_template_gamma_z_azimu_b = atan2(point_pole_template_gamma_z_k_c_(1+1),point_pole_template_gamma_z_k_c_(1+0));
tmp_error_2 = tmp_error_2 + fnorm(periodize(point_output_azimu_b - point_pole_template_gamma_z_azimu_b,-pi,+pi)); %<-- should be 0. ;
point_pole_template_gamma_z_polar_a = atan2(point_pole_template_gamma_z_k_r01,point_pole_template_gamma_z_k_c_(1+2));
tmp_error_3 = tmp_error_3 + fnorm(periodize(+pi/2+point_output_polar_a - point_pole_template_gamma_z_polar_a,-pi/2,+pi/2)); %<-- should be 0. ;
%%%%%%%%;
clear point_pole_k_c_;
clear point_pole_k_r01;
clear point_pole_k_r012;
clear point_pole_azimu_b;
clear point_pole_polar_a;
clear point_pole_template_gamma0_k_c_;
clear point_pole_template_sgx_k_c_;
clear point_pole_template_sgx;
clear point_pole_template_cgx;
clear point_pole_gamma_z;
clear point_pole_template_gamma_z_k_c_;
clear point_pole_template_gamma_z_k_r01;
clear point_pole_template_gamma_z_k_r012;
clear point_pole_template_gamma_z_azimu_b;
clear point_pole_template_gamma_z_polar_a;
%%%%%%%%%%%%%%%%;
end;%for npole=0:n_w_max-1;
%%%%%%%%%%%%%%%%;
clear point_output_azimu_b;
clear point_output_polar_a;
clear point_output_gamma_z = 0;
clear point_output_k_c_;
clear point_output_k_r01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;end;%for npoint_a=0:n_point_a-1; for npoint_b=0:n_point_b-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
save(fname_mat ...
,'n_w_max' ...
,'n_w_0in_' ...
,'gamma_z_' ...
,'cc_' ...
,'sc_' ...
,'n_point_a' ...
,'polar_a_a_' ...
,'n_point_b' ...
,'azimu_b_b_' ...
,'point_output_azimu_b_ab__' ...
,'point_output_polar_a_ab__' ...
,'point_output_k_c_0_ab__' ...
,'point_output_k_c_1_ab__' ...
,'point_output_k_c_2_ab__' ...
,'point_pole_azimu_b_abw___' ...
,'point_pole_polar_a_abw___' ...
,'point_pole_k_c_0_abw___' ...
,'point_pole_k_c_1_abw___' ...
,'point_pole_k_c_2_abw___' ...
,'point_pole_gamma_z_abw___' ...
);
end;%if (~exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% set up interpolator. ;
%%%%%%%%;
tmp_t = tic();
n_k_p_uni_shell = n_k_p_uni_azimu_b * n_k_p_uni_polar_a ;
[k_p_uni_azimu_b_ba_,k_p_uni_polar_a_ba_] = ndgrid(k_p_uni_azimu_b_,k_p_uni_polar_a_);
n_order = 5;
flag_polar_a_ascend_vs_descend = (polar_a_a_(end)>=polar_a_a_(1));
[ ...
 shell_scatter_from_tensor_sba__ ...
] = ...
shell_k_p_scatter_from_tensor_interpolate_n_5( ...
 n_order ...
,n_point_b ...
,n_point_a ...
,polar_a_a_ ...
,n_k_p_uni_shell ...
,k_p_uni_azimu_b_ba_ ...
,k_p_uni_polar_a_ba_ ...
,flag_polar_a_ascend_vs_descend ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% shell_scatter_from_tensor_sba__ time %0.2fs',tmp_t));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now if we have a particular (continuous) perturbation field: ;
%%%%;
% g_polar_a = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_polar_a = @(point_pole_polar_a,point_pole_azimu_b) point_pole_polar_a + epsilon_polar_a * g_polar_a(point_pole_polar_a,point_pole_azimu_b) ;
% g_azimu_b = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_azimu_b = @(point_pole_polar_a,point_pole_azimu_b) point_pole_azimu_b + epsilon_azimu_b * g_azimu_b(point_pole_polar_a,point_pole_azimu_b) ;
% g_gamma_z = @(point_pole_polar_a,point_pole_azimu_b) < insert function here > ;
% f_gamma_z = @(point_pole_polar_a,point_pole_azimu_b) point_pole_gamma_z + epsilon_gamma_z * g_gamma_z(point_pole_polar_a,point_pole_azimu_b) ;
%%%%;
% Assuming that epsilon_ = [epsilon_polar_a,epsilon_azimu_b,epsilon_gamma_z] is small, the inverse (to first order) is given by: ;
% f_polar_a_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_polar_a_inverse = point_pole_polar_a - epsilon_polar_a * g_polar_a(point_pole_polar_a,point_pole_azimu_b) ;
% f_azimu_b_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_azimu_b_inverse = point_pole_azimu_b - epsilon_azimu_b * g_azimu_b(point_pole_polar_a,point_pole_azimu_b) ;
% f_gamma_z_inverse(point_pole_polar_a,point_pole_azimu_b) = ... ;
% point_pole_gamma_z_inverse = point_pole_gamma_z - epsilon_gamma_z * g_gamma_z(point_pole_polar_a,point_pole_azimu_b) ;
% with the final formula (for gamma_z_inverse) being exact. ;
%%%%%%%%;
% Thus, the location on the sphere corresponding to: ;
% point_pole_predilated_template_gamma_z_k_c_ = [...
%  +cos(point_pole_azimu_b_inverse)*cos(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) - sin(point_pole_azimu_b_inverse)*sin(point_pole_gamma_z_inverse) ...
% ;+sin(point_pole_azimu_b_inverse)*cos(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) + cos(point_pole_azimu_b_inverse)*sin(point_pole_gamma_z_inverse) ...
% ;-sin(point_pole_polar_a_inverse)*cos(point_pole_gamma_z_inverse) ...
% ];
%%%%%%%%;
% Using the notation: ;
% ca, sa = cos and sin of polar_a, ;
% cb, sb = cos and sin of azimu_b, ;
% cc, sc = cos and sin of gamma_z, ;
% we have the relations: ;
% x0 = +cb*ca*cc - sb*sc ;
% x1 = +sb*ca*cc + cb*sc ;
% x2 = -sa*cc ;
% and the jacobian: ;
% dx0da = -cb*sa*cc ;
% dx0db = -sb*ca*cc - cb*sc ;
% dx0dc = -cb*ca*sc - sb*cc ;
% dx1da = -sb*sa*cc ;
% dx1db = +cb*ca*cc - sb*sc ;
% dx1dc = -sb*ca*sc + cb*cc ;
% dx2da = -ca*cc ;
% dx2db = 0;
% dx2dc = +sa*sc ;
%%%%%%%%;
ca_abw___ = cos(point_pole_polar_a_abw___); sa_abw___ = sin(point_pole_polar_a_abw___);
cb_abw___ = cos(point_pole_azimu_b_abw___); sb_abw___ = sin(point_pole_azimu_b_abw___);
cc_abw___ = cos(point_pole_gamma_z_abw___); sc_abw___ = sin(point_pole_gamma_z_abw___);
x0_abw___ = +cb_abw___.*ca_abw___.*cc_abw___ - sb_abw___.*sc_abw___ ;
x1_abw___ = +sb_abw___.*ca_abw___.*cc_abw___ + cb_abw___.*sc_abw___ ;
x2_abw___ = -sa_abw___.*cc_abw___ ;
dx0da_abw___ = -cb_abw___.*sa_abw___.*cc_abw___ ;
dx0db_abw___ = -sb_abw___.*ca_abw___.*cc_abw___ - cb_abw___.*sc_abw___ ;
dx0dc_abw___ = -cb_abw___.*ca_abw___.*sc_abw___ - sb_abw___.*cc_abw___ ;
dx1da_abw___ = -sb_abw___.*sa_abw___.*cc_abw___ ;
dx1db_abw___ = +cb_abw___.*ca_abw___.*cc_abw___ - sb_abw___.*sc_abw___ ;
dx1dc_abw___ = -sb_abw___.*ca_abw___.*sc_abw___ + cb_abw___.*cc_abw___ ;
dx2da_abw___ = -ca_abw___.*cc_abw___ ;
dx2db_abw___ = +zeros(n_point_a,n_point_b,n_w_max);
dx2dc_abw___ = +sa_abw___.*sc_abw___ ;
%%%%%%%%;
% Now, the accumulated dx0, dx1, dx2 are given by: ;
% g_a_abw___ = bsxfun(@plus,g_polar_a_abwj____,reshape(eps_a_j_,[1,1,1,n_a_j]));
% g_b_abw___ = bsxfun(@plus,g_azimu_b_abwj____,reshape(eps_b_j_,[1,1,1,n_b_j]));
% g_c_abw___ = bsxfun(@plus,g_gamma_z_abwj____,reshape(eps_c_j_,[1,1,1,n_c_j]));
% dx0_abw___ = -dx0da_abw___.*g_a_abw___ - dx0db_abw___.*g_b_abw___ - dx0dc_abw___.*g_c_abw___ ;
% dx1_abw___ = -dx1da_abw___.*g_a_abw___ - dx1db_abw___.*g_b_abw___ - dx1dc_abw___.*g_c_abw___ ;
% dx2_abw___ = -dx2da_abw___.*g_a_abw___ - dx2db_abw___.*g_b_abw___ - dx2dc_abw___.*g_c_abw___ ;
% or, written in terms of the eps_a_j_, eps_b_j_, eps_c_j_, ;
% dx0_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx0da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx0db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx0dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% dx1_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx1da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx1db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx1dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% dx2_abw___ = ... ;
% + bsxfun(@times,bsxfun(@times,-dx2da_abw___,g_polar_a_abwj____),reshape(eps_a_j_,[1,1,1,n_a_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx2db_abw___,g_azimu_b_abwj____),reshape(eps_b_j_,[1,1,1,n_b_j])) ... ;
% + bsxfun(@times,bsxfun(@times,-dx2dc_abw___,g_gamma_z_abwj____),reshape(eps_c_j_,[1,1,1,n_c_j])) ... ;
% ;
% collectively producing the shift: ;
% dx_abw3____ = cat(4,dx0_abw___,dx1_abw___,dx2_abw___);
% Note that this shift, while defined in terms of a vector in \Real^{3}, ;
% is constrained to lie tangent to the sphere. ;
%%%%%%%%;
% Now we can record viewing_weight_empirical_abw___ as a function of point_pole_polar_a_abw___ and point_pole_azimu_b_abw___. ;
% If, e.g., this is uniform, then viewing_weight_empirical_abw___ = constant. ;
% Note that this does not account for quadrature on the sphere, ;
% as these weights will be interpreted as contributions from each pole on an equispaced ring (i.e., n_pole = n_w_max). ;
% With this measure we can calculate various derivatives. ;
% The viewing_weight associated with each ring is: ;
% viewing_ring_weight_empirical_abw___ = bsxfun(@rdivide,viewing_weight_empirical_abw___,max(1e-12,sum(viewing_weight_empirical_abw___,3))). ;
% Once again, this is interpreted as a function of point_pole_polar_a_abw___ and point_pole_azimu_b_abw___. ;
%%%%%%%%;
% The average shift at each point_output_k_c_ (indexed by ab13) is now determined by: ;
% dx_avg_ab13___ = sum(bsxfun(@times,dx_abw3___,viewing_ring_weight_empirical_abw___),3). ;
% Meanwhile, the average variation at each point_output_k_c_ is now determined by: ;
% dx_var_ab__ = sum(bsxfun(@times,bsxfun(@minus,dx_abw3___,dx_avg_ab13___).^2,viewing_ring_weight_empirical_abw___),[3,4]). ;
%%%%%%%%;
% Note that the dx_var_ab__ can be written as a quadratic-function of eps_: [eps_a_j_;eps_b_j_;eps_c_j_]. ;
% This immediately implies that dx_var_ab__ can be optimized (subject to a norm-constraint on eps_) ;
% by taking the eigenvectors of the appropriate quadratic-kernel. ;
%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% quick check of tmp_b_est. ;
%%%%%%%%;
tmp_a = 1*pi*rand(); ca = cos(tmp_a); sa = sin(tmp_a);
tmp_b = 2*pi*rand(); cb = cos(tmp_b); sb = sin(tmp_b);
tmp_c = 2*pi*rand(); cc = cos(tmp_c); sc = sin(tmp_c);
tmp_x0 = cb*ca*cc - sb*sc ;
tmp_x1 = sb*ca*cc + cb*sc ;
tmp_x2 = -sa*cc ;
tmp_r01 = sqrt(tmp_x0.^2 + tmp_x1.^2);
tmp_b_out = atan2(tmp_x1,tmp_x0);
tmp_a_out = atan2(tmp_r01,tmp_x2);
tmp_b_est = atan2(sc,ca*cc) + tmp_b;
tmp_x_ = Rz(tmp_b)*Ry(tmp_a)*Rz(tmp_c)*[1;0;0];
disp(sprintf(' %% fnorm tmp_x_ - [tmp_x0;tmp_x1;tmp_x2]: %0.16f',fnorm(tmp_x_-[tmp_x0;tmp_x1;tmp_x2])));
tmp_f_a = @(x_) atan2(sqrt(x_(1+0).^2 + x_(1+1).^2),x_(1+2));
disp(sprintf(' %% fnorm tmp_f_a - tmp_a_out: %0.16f',fnorm(tmp_f_a(tmp_x_) - tmp_a_out)));
tmp_f_b = @(x_) atan2(x_(1+1),x_(1+0));
disp(sprintf(' %% fnorm tmp_f_b - tmp_b_out: %0.16f',fnorm(tmp_f_b(tmp_x_) - tmp_b_out)));
disp(sprintf(' %% fnorm tmp_b_est - tmp_b_out: %0.16f',fnorm(periodize(tmp_b_est - tmp_b_out,-pi,+pi))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Note that, typically speaking, an arbitrary alignment-perturbation will result in ;
% a nonzero 'first-order-term' (associated with the average). ;
% More specifically, this average will have a magnitude proportional to fnorm(epsilon_), ;
% and will point essentially arbitrarily at every point_output_k_c_. ;
%%%%;
% The consequence of this haphazard average drift is that any particular template ;
% (originally aligned to the unperturbed volume) ;
% will now no longer align with *any* locally-perturbed template from the perturbed volume. ;
% This mis-alignment will result in a first-order shift in the likelihood ;
% (proportional to fnorm(epsilon_)). ;
%%%%;
% In special cases, such as a longitudinal or latitudinal perturbation, ;
% each template from the original volume will align to first-order with a particular ;
% locally-perturbed template from the perturbed volume. ;
% In this case any mis-alignment will result from the second-order shift in the likelihood ;
% (proportional to fnorm(epsilon_).^2) ;
% associated with the local diffusion in the perturbed volume  ;
% induced by dx_var_ab__. ;
%%%%;
% For these special cases we can bound the minimum-curvature of the likelihood-landscape ;
% by solving the eigenvalue-problem described above. ;
%%%%%%%%;

%%%%%%%%;
% Visualize mean and variance for particular g_azimu_b. ;
%%%%%%%%;
equa_sigma = pi/2 * (1/14);%equa_sigma = pi/2 * (1/14);
viewing_weight_empirical_abw___ = ones(n_point_a,n_point_b,n_w_max);
if (equa_sigma> 0);
viewing_weight_empirical_abw___ = exp(-(periodize(point_pole_polar_a_abw___,0,pi)-pi/2).^2/(2*equa_sigma.^2)); %<-- roughly 90% of mass within pi/14 of equator. ;
end;%if (equa_sigma> 0);
viewing_ring_weight_empirical_abw___ = bsxfun(@rdivide,viewing_weight_empirical_abw___,max(1e-12,sum(viewing_weight_empirical_abw___,3))) ;
%%%%%%%%;
% Test quadratic representation for a single mode: ;
%%%%%%%%;
g_freq = +2;
g_azimu_b_abw___ = sin(g_freq*point_pole_azimu_b_abw___);
dx0db_avg_ab__ = sum(dx0db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx1db_avg_ab__ = sum(dx1db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx2db_avg_ab__ = sum(dx2db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx0db_var_ab__ = sum(bsxfun(@minus,dx0db_abw___.*g_azimu_b_abw___,dx0db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dx1db_var_ab__ = sum(bsxfun(@minus,dx1db_abw___.*g_azimu_b_abw___,dx1db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dx2db_var_ab__ = sum(bsxfun(@minus,dx2db_abw___.*g_azimu_b_abw___,dx2db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dxdb_avg_ab__ = sqrt(dx0db_avg_ab__.^2 + dx1db_avg_ab__.^2 + dx2db_avg_ab__.^2) ;
dxdb_var_ab__ =     (dx0db_var_ab__ + dx1db_var_ab__ + dx2db_var_ab__) ;
dxdb_std_ab__ = sqrt(dx0db_var_ab__ + dx1db_var_ab__ + dx2db_var_ab__) ;
P0_ab__ = +sum(dx0db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P1_ab__ = +sum(dx1db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P2_ab__ = +sum(dx2db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
Q_ab__ = ...
+sum((dx0db_abw___.*g_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g_azimu_b_abw___),3) ...
+sum((dx1db_abw___.*g_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g_azimu_b_abw___),3) ...
+sum((dx2db_abw___.*g_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g_azimu_b_abw___),3) ...
-P0_ab__.*P0_ab__ -P1_ab__.*P1_ab__ -P2_ab__.*P2_ab__ ...
;
disp(sprintf(' %% Q_ab__ vs dxdb_var_ab__: %0.16f',fnorm(Q_ab__ - dxdb_var_ab__)/fnorm(Q_ab__)));
%%%%%%%%;
% Test quadratic representation for two modes. ;
%%%%%%%%
g0_freq = +2; g1_freq = +4;
g0_azimu_b_abw___ = sin(g0_freq*point_pole_azimu_b_abw___);
g1_azimu_b_abw___ = cos(g1_freq*point_pole_azimu_b_abw___);
eps0_azimu_b = -0.5; eps1_azimu_b = +sqrt(3)/2; 
g_azimu_b_abw___ = eps0_azimu_b*g0_azimu_b_abw___ + eps1_azimu_b*g1_azimu_b_abw___ ;
dx0db_avg_ab__ = sum(dx0db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx1db_avg_ab__ = sum(dx1db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx2db_avg_ab__ = sum(dx2db_abw___.*g_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3);
dx0db_var_ab__ = sum(bsxfun(@minus,dx0db_abw___.*g_azimu_b_abw___,dx0db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dx1db_var_ab__ = sum(bsxfun(@minus,dx1db_abw___.*g_azimu_b_abw___,dx1db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dx2db_var_ab__ = sum(bsxfun(@minus,dx2db_abw___.*g_azimu_b_abw___,dx2db_avg_ab__).^2.*viewing_ring_weight_empirical_abw___,3);
dxdb_avg_ab__ = sqrt(dx0db_avg_ab__.^2 + dx1db_avg_ab__.^2 + dx2db_avg_ab__.^2) ;
dxdb_var_ab__ =     (dx0db_var_ab__ + dx1db_var_ab__ + dx2db_var_ab__) ;
dxdb_std_ab__ = sqrt(dx0db_var_ab__ + dx1db_var_ab__ + dx2db_var_ab__) ;
P0_abg___ = zeros(n_point_a,n_point_b,2);
P1_abg___ = zeros(n_point_a,n_point_b,2);
P2_abg___ = zeros(n_point_a,n_point_b,2);
P0_abg___(:,:,1+0) = +sum(dx0db_abw___.*g0_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P1_abg___(:,:,1+0) = +sum(dx1db_abw___.*g0_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P2_abg___(:,:,1+0) = +sum(dx2db_abw___.*g0_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P0_abg___(:,:,1+1) = +sum(dx0db_abw___.*g1_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P1_abg___(:,:,1+1) = +sum(dx1db_abw___.*g1_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P2_abg___(:,:,1+1) = +sum(dx2db_abw___.*g1_azimu_b_abw___.*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
Q_abgg____ = zeros(n_point_a,n_point_b,2,2);
Q_abgg____(:,:,1+0,1+0) = ...
+sum((dx0db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g0_azimu_b_abw___),3) ...
+sum((dx1db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g0_azimu_b_abw___),3) ...
+sum((dx2db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g0_azimu_b_abw___),3) ...
-P0_abg___(:,:,1+0).*P0_abg___(:,:,1+0) -P1_abg___(:,:,1+0).*P1_abg___(:,:,1+0) -P2_abg___(:,:,1+0).*P2_abg___(:,:,1+0) ...
;
Q_abgg____(:,:,1+0,1+1) = ...
+sum((dx0db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g1_azimu_b_abw___),3) ...
+sum((dx1db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g1_azimu_b_abw___),3) ...
+sum((dx2db_abw___.*g0_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g1_azimu_b_abw___),3) ...
-P0_abg___(:,:,1+0).*P0_abg___(:,:,1+1) -P1_abg___(:,:,1+0).*P1_abg___(:,:,1+1) -P2_abg___(:,:,1+0).*P2_abg___(:,:,1+1) ...
;
Q_abgg____(:,:,1+1,1+0) = ...
+sum((dx0db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g0_azimu_b_abw___),3) ...
+sum((dx1db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g0_azimu_b_abw___),3) ...
+sum((dx2db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g0_azimu_b_abw___),3) ...
-P0_abg___(:,:,1+1).*P0_abg___(:,:,1+0) -P1_abg___(:,:,1+1).*P1_abg___(:,:,1+0) -P2_abg___(:,:,1+1).*P2_abg___(:,:,1+0) ...
;
Q_abgg____(:,:,1+1,1+1) = ...
+sum((dx0db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g1_azimu_b_abw___),3) ...
+sum((dx1db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g1_azimu_b_abw___),3) ...
+sum((dx2db_abw___.*g1_azimu_b_abw___).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g1_azimu_b_abw___),3) ...
-P0_abg___(:,:,1+1).*P0_abg___(:,:,1+1) -P1_abg___(:,:,1+1).*P1_abg___(:,:,1+1) -P2_abg___(:,:,1+1).*P2_abg___(:,:,1+1) ...
;
Q_ab__ = ...
+eps0_azimu_b*Q_abgg____(:,:,1+0,1+0)*eps0_azimu_b ...
+eps0_azimu_b*Q_abgg____(:,:,1+0,1+1)*eps1_azimu_b ...
+eps1_azimu_b*Q_abgg____(:,:,1+1,1+0)*eps0_azimu_b ...
+eps1_azimu_b*Q_abgg____(:,:,1+1,1+1)*eps1_azimu_b ...
;
disp(sprintf(' %% Q_ab__ vs dxdb_var_ab__: %0.16f',fnorm(Q_ab__ - dxdb_var_ab__)/fnorm(Q_ab__)));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024*2,768]);
markersize_use = 8;
%%%%;
subplot(1,4,[1,2]);
hold on;
plot(periodize(point_output_azimu_b_ab__(:),-pi,+pi),pi/2 + point_output_polar_a_ab__(:),'k.','MarkerSize',1);
for npoint_a=0:4:n_point_a-1;
for npoint_b=0:8:n_point_b-1;
plot(periodize(point_output_azimu_b_ab__(1+npoint_a,1+npoint_b),-pi,+pi),pi/2 + point_output_polar_a_ab__(1+npoint_a,1+npoint_b),'go','MarkerSize',markersize_use);
tmp_x0_w_ = squeeze(x0_abw___(1+npoint_a,1+npoint_b,:));
tmp_x1_w_ = squeeze(x1_abw___(1+npoint_a,1+npoint_b,:));
tmp_x2_w_ = squeeze(x2_abw___(1+npoint_a,1+npoint_b,:));
tmp_r01_w_ = sqrt(tmp_x0_w_.^2 + tmp_x1_w_.^2);
tmp_azimu_b_w_ = atan2(tmp_x1_w_,tmp_x0_w_);
tmp_polar_a_w_ = atan2(tmp_r01_w_,tmp_x2_w_);
plot(tmp_azimu_b_w_,tmp_polar_a_w_,'bs','MarkerSize',markersize_use);
tmp_g_w_ = squeeze(g_azimu_b_abw___(1+npoint_a,1+npoint_b,:));
tmp_viewing_ring_weight_empirical_w_ = squeeze(viewing_ring_weight_empirical_abw___(1+npoint_a,1+npoint_b,:));
tmp_eps_b = 1e-2;
tmp_dx0db_w_ = tmp_x0_w_ + squeeze(dx0db_abw___(1+npoint_a,1+npoint_b,:)).*tmp_g_w_*tmp_eps_b;
tmp_dx1db_w_ = tmp_x1_w_ + squeeze(dx1db_abw___(1+npoint_a,1+npoint_b,:)).*tmp_g_w_*tmp_eps_b;
tmp_dx2db_w_ = tmp_x2_w_ + squeeze(dx2db_abw___(1+npoint_a,1+npoint_b,:)).*tmp_g_w_*tmp_eps_b;
tmp_dx0db_avg = sum(tmp_dx0db_w_.*tmp_viewing_ring_weight_empirical_w_);
tmp_dx1db_avg = sum(tmp_dx1db_w_.*tmp_viewing_ring_weight_empirical_w_);
tmp_dx2db_avg = sum(tmp_dx2db_w_.*tmp_viewing_ring_weight_empirical_w_);
tmp_dr01_avg = sqrt(tmp_dx0db_avg.^2 + tmp_dx1db_avg.^2);
tmp_azimu_b_avg = atan2(tmp_dx1db_avg,tmp_dx0db_avg);
tmp_polar_a_avg = atan2(tmp_dr01_avg,tmp_dx2db_avg);
plot(tmp_azimu_b_avg,tmp_polar_a_avg,'mh','MarkerSize',markersize_use);
tmp_dx0db_var = sum((tmp_dx0db_w_-tmp_dx0db_avg).^2.*tmp_viewing_ring_weight_empirical_w_);
assert(abs(tmp_dx0db_var - tmp_eps_b.^2*dx0db_var_ab__(1+npoint_a,1+npoint_b))<1e-12);
tmp_dx1db_var = sum((tmp_dx1db_w_-tmp_dx1db_avg).^2.*tmp_viewing_ring_weight_empirical_w_);
assert(abs(tmp_dx1db_var - tmp_eps_b.^2*dx1db_var_ab__(1+npoint_a,1+npoint_b))<1e-12);
tmp_dx2db_var = sum((tmp_dx2db_w_-tmp_dx2db_avg).^2.*tmp_viewing_ring_weight_empirical_w_);
assert(abs(tmp_dx2db_var - tmp_eps_b.^2*dx2db_var_ab__(1+npoint_a,1+npoint_b))<1e-12);
tmp_dxdb_var = tmp_dx0db_var + tmp_dx1db_var + tmp_dx2db_var ;
assert(abs(tmp_dxdb_var - tmp_eps_b.^2*dxdb_var_ab__(1+npoint_a,1+npoint_b))<1e-12);
tmp_dr01_w_ = sqrt(tmp_dx0db_w_.^2 + tmp_dx1db_w_.^2);
tmp_dazimu_b_w_ = atan2(tmp_dx1db_w_,tmp_dx0db_w_);
tmp_dpolar_a_w_ = atan2(tmp_dr01_w_,tmp_dx2db_w_);
plot(tmp_dazimu_b_w_,tmp_dpolar_a_w_,'rx','MarkerSize',markersize_use);
end;%for npoint_a=0:n_point_a-1;
end;%for npoint_b=0:n_point_b-1;
hold off;
xlabel('azimu_b','Interpreter','none'); ylabel('polar_a','Interpreter','none');
xlim([-1*pi,+1*pi]); ylim([0,1*pi]);
set(gca,'XTick',pi*[0:0.25:2],'XTickLabel',{'0\pi/4','1\pi/4','2\pi/4','3\pi/4','4\pi/4','5\pi/4','6\pi/4','7\pi/4','8\pi/4'});
set(gca,'YTick',pi*[0:0.25:1],'YTickLabel',{'0\pi/4','1\pi/4','2\pi/4','3\pi/4','4\pi/4'});
%%%%;
subplot(1,4,[3,4]); fig81s;
tmp_ab__ = dxdb_var_ab__./max(1e-12,sin(pi/2 + point_output_polar_a_ab__).^2);
imagesc(tmp_ab__,[0,max(tmp_ab__,[],'all')]);colorbar;
title('azimuthal variance (scaled)');
axisnotick;
xlabel('azimu_b','Interpreter','none'); ylabel('polar_a','Interpreter','none');
%%%%;
clear tmp_x0_w_ ;
clear tmp_x1_w_ ;
clear tmp_x2_w_ ;
clear tmp_r01_w_ ;
clear tmp_azimu_b_w_ ;
clear tmp_polar_a_w_ ;
clear tmp_g_w_ ;
clear tmp_viewing_ring_weight_empirical_w_ ;
clear tmp_eps_b ;
clear tmp_dx0db_w_ ;
clear tmp_dx1db_w_ ;
clear tmp_dx2db_w_ ;
clear tmp_dx0db_avg ;
clear tmp_dx1db_avg ;
clear tmp_dx2db_avg ;
clear tmp_dr01_avg ;
clear tmp_azimu_b_avg ;
clear tmp_polar_a_avg ;
clear tmp_dx0db_var ;
clear tmp_dx1db_var ;
clear tmp_dx2db_var ;
clear tmp_dxdb_var ;
clear tmp_dr01_w_ ;
clear tmp_dazimu_b_w_ ;
clear tmp_dpolar_a_w_ ;
clear tmp_ab__;
clear equa_sigma viewing_weight_empirical_abw___ viewing_ring_weight_empirical_abw___ g_freq g_azimu_b_abw___ ;
clear dx0db_avg_ab__ dx1db_avg_ab__ dx2db_avg_ab__ dx0db_var_ab__ dx1db_var_ab__ dx2db_var_ab__ ;
%%%%%%%%;

%%%%%%%%;
% We demonstrate this calculation using LetB1, ;
% Note that the true viewing_weight_empirical_abw___ corresponds to equa_sigma close to pi/14/2. ;
%%%%%%%%;
n_g2 = 17;
g_freq_g2_ = 2:2:+(n_g2-1)/2; g_freq_g2_ = setdiff(g_freq_g2_,0); n_g2 = numel(g_freq_g2_); %<-- frequency 0 is a global rotation. ;
g_freq_g_ = repmat(g_freq_g2_,[2,1]); g_freq_g_ = g_freq_g_(:); n_g = numel(g_freq_g_);
g_azimu_b_abwg____ = zeros(n_point_a,n_point_b,n_w_max,n_g);
for ng=0:2:n_g-1;
g_freq = g_freq_g_(1+ng+0);
g_azimu_b_abwg____(:,:,:,1+ng+0) = cos(g_freq*point_pole_azimu_b_abw___);
g_freq = g_freq_g_(1+ng+1);
g_azimu_b_abwg____(:,:,:,1+ng+1) = sin(g_freq*point_pole_azimu_b_abw___);
end;%for ng=0:2:n_g-1;
%%%%%%%%;
equa_sigma_ = pi/2 * exp(linspace(log2(4),log2(1/16),25)); n_equa_sigma = numel(equa_sigma_);
Q_ggs___ = zeros(n_g,n_g,n_equa_sigma);
%%%%%%%%%%%%%%%%;
for nequa_sigma=0:n_equa_sigma-1;
%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% nequa_sigma %d/%d',nequa_sigma,n_equa_sigma)); end;
equa_sigma = equa_sigma_(1+nequa_sigma);
viewing_weight_empirical_abw___ = exp(-(periodize(point_pole_polar_a_abw___,0,pi)-pi/2).^2/(2*equa_sigma.^2)); %<-- roughly 90% of mass within pi/14 of equator. ;
viewing_ring_weight_empirical_abw___ = bsxfun(@rdivide,viewing_weight_empirical_abw___,max(1e-12,sum(viewing_weight_empirical_abw___,3))) ;
%%%%%%%%;
% Now the second-order contribution to the variance has a quadratic-kernel of the form: ;
%%%%%%%%;
tmp_t = tic();
P0_abg___ = zeros(n_point_a,n_point_b,n_g);
P1_abg___ = zeros(n_point_a,n_point_b,n_g);
P2_abg___ = zeros(n_point_a,n_point_b,n_g);
for ng0=0:n_g-1;
P0_abg___(:,:,1+ng0) = +sum((dx0db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P1_abg___(:,:,1+ng0) = +sum((dx1db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
P2_abg___(:,:,1+ng0) = +sum((dx2db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___,3) ; %<-- average. ;
end;%for ng0=0:n_g-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% PX_abg___ time %0.2fs',tmp_t));
%%%%;
tmp_t = tic();
Q_abgg____ = zeros(n_point_a,n_point_b,n_g,n_g);
for ng0=0:n_g-1;for ng1=0:n_g-1;
tmp_term0_ab__ = ...
+sum((dx0db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___.*(dx0db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng1)),3) ...
+sum((dx1db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___.*(dx1db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng1)),3) ...
+sum((dx2db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng0)).*viewing_ring_weight_empirical_abw___.*(dx2db_abw___.*g_azimu_b_abwg____(:,:,:,1+ng1)),3) ...
;
tmp_term1_ab__ = P0_abg___(:,:,1+ng0).*P0_abg___(:,:,1+ng1) + P1_abg___(:,:,1+ng0).*P1_abg___(:,:,1+ng1) + P2_abg___(:,:,1+ng0).*P2_abg___(:,:,1+ng1) ;
Q_abgg____(:,:,1+ng0,1+ng1) = tmp_term0_ab__ - tmp_term1_ab__ ;
clear tmp_term0_ab__ tmp_term1_ab__ ;
end;end;%for ng1=0:n_g-1;for ng0=0:n_g-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% Q_abgg____ time %0.2fs',tmp_t));
%%%%;
tmp_t = tic();
Q_gg__ = zeros(n_g,n_g);
tmp_wa_bak___ = real(weight_3d_k_uni_bak___ .* abs(a_k_p_uni_dbdb_bak___) ./ max(1e-12,a_k_p_uni_dbdb_fact_bak___)) ; %<-- take abs of concavity. ;
for ng0=0:n_g-1;for ng1=0:n_g-1;
tmp_Q_ba__ = reshape(shell_scatter_from_tensor_sba__ * reshape(permute(Q_abgg____(:,:,1+ng0,1+ng1),[2,1]),[n_point_b*n_point_a,1]),[n_k_p_uni_azimu_b,n_k_p_uni_polar_a]);
Q_gg__(1+ng0,1+ng1) = sum(bsxfun(@times,tmp_Q_ba__,tmp_wa_bak___),'all');
clear tmp_Q_ba__ ;
end;end;%for ng0=0:n_g-1;for ng1=0:n_g-1;
clear tmp_wa_bak___ ;
tmp_t = toc(tmp_t); disp(sprintf(' %% Q_gg__ time %0.2fs',tmp_t));
Q_ggs___(:,:,1+nequa_sigma) = Q_gg__;
clear viewing_weight_empirical_abw___ viewing_ring_weight_empirical_abw___ P0_abg___ P1_abg___ P2_abg___ Q_abgg____ Q_gg__ ;
%%%%%%%%;
end;%for nequa_sigma=0:n_equa_sigma-1;
%%%%%%%%%%%%%%%%;
S_sg__ = zeros(n_equa_sigma,n_g);
for nequa_sigma=0:n_equa_sigma-1;
S_sg__(1+nequa_sigma,:) = eigs(Q_ggs___(:,:,1+nequa_sigma),n_g);
end;%for nequa_sigma=0:n_equa_sigma-1;
[tmp_V_,tmp_S] = eigs(Q_ggs___(:,:,end),1,'smallestabs');
%%%%%%%%;

fname_fig_pre = sprintf('%s_jpg/test_equatorial_band_dilate_LetB1_FisherInformation_FIGA',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_jpg = sprintf('%s_stripped.jpg',fname_fig_pre);
fname_fig_stripped_eps = sprintf('%s_stripped.eps',fname_fig_pre);
%%%%%%%%;
if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
%%%%;
g_azimu_b_ba = @(azimu_b,polar_a) ...
+0.125*tmp_V_(1+0)*sin(g_freq_g_(1+0)*azimu_b) ...
+0.125*tmp_V_(1+1)*cos(g_freq_g_(1+1)*azimu_b) ...
+0.125*tmp_V_(1+2)*sin(g_freq_g_(1+2)*azimu_b) ...
+0.125*tmp_V_(1+3)*cos(g_freq_g_(1+3)*azimu_b) ...
+0.125*tmp_V_(1+4)*sin(g_freq_g_(1+4)*azimu_b) ...
+0.125*tmp_V_(1+5)*cos(g_freq_g_(1+5)*azimu_b) ...
+0.125*tmp_V_(1+6)*sin(g_freq_g_(1+6)*azimu_b) ...
+0.125*tmp_V_(1+7)*cos(g_freq_g_(1+7)*azimu_b) ...
;
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
linewidth_use = 3;
p_row = 2; p_col = 2; np=0;
fontsize_use = 12;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
isosurface_f_x_u_1([],a_x_u_reco_); title('LetB1');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot_sphere_grid_0(struct('flag_solid',1));
hold on;
sphere_fun_post__0( ...
 [] ...
,[] ...
,[] ...
,[] ...
,[] ...
,g_azimu_b_ba ...
,[] ...
);
hold off;
xlabel('x0');
ylabel('x1');
zlabel('x2');
axis equal; %axis vis3d;
view(0,0); title('perturbation (side)');
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
plot_sphere_grid_0(struct('flag_solid',1));
hold on;
sphere_fun_post__0( ...
 [] ...
,[] ...
,[] ...
,[] ...
,[] ...
,g_azimu_b_ba ...
,[] ...
);
hold off;
xlabel('x0');
ylabel('x1');
zlabel('x2');
axis equal; %axis vis3d;
view(0,90); title('perturbation (top)');
set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
tmp_S_s_ = S_sg__(:,end);
plot(log(equa_sigma_*2/pi),tmp_S_s_/max(1e-12,max(tmp_S_s_)),'k-','LineWidth',linewidth_use);
set(gca,'XTick',flip(log(equa_sigma_*2/pi))); xtickangle(90); grid on;
xlim([min(log(equa_sigma_*2/pi))-0.5,max(log(equa_sigma_*2/pi))+0.5]);
xlabel('log(2*equa_sigma/pi)','Interpreter','none');
ylabel('scaled curvature');
set(gca,'Ytick',[0:0.1:1],'YTickLabel',{'0',[],[],[],[],[],[],[],[],[],'max'});
title('minimum curvature');
set(gca,'FontSize',fontsize_use);
%%%%;
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_stripped_jpg);
%print('-depsc',fname_fig_stripped_eps);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
%print('-depsc',fname_fig_eps);
%%%%;
close(gcf);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
if  exist(fname_fig_jpg,'file');
disp(sprintf(' %% %s found, not creating',fname_fig_jpg));
end;%if  exist(fname_fig_jpg,'file');
%%%%%%%%;











