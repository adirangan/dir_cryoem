clear;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

flag_recalc = 0;
flag_replot = 0;
tolerance_master = 1e-2;
nf=0;

dir_pm = sprintf('/%s/rangan/dir_cryoem/dir_trpv1/dir_pm',string_root);
if (~exist(sprintf('%s_mat',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_pm)); mkdir(sprintf('%s_mat',dir_pm)); end;
if (~exist(sprintf('%s_jpg',dir_pm),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_pm)); mkdir(sprintf('%s_jpg',dir_pm)); end;
dir_relion = sprintf('/%s/rangan/dir_cryoem/dir_trpv1/dir_relion',string_root);
if (~exist(sprintf('%s_mat',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_relion)); mkdir(sprintf('%s_mat',dir_relion)); end;
if (~exist(sprintf('%s_jpg',dir_relion),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_relion)); mkdir(sprintf('%s_jpg',dir_relion)); end;
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);
dir_data_marina = sprintf('/%s/rangan/dir_cryoem/dir_trpv1/data_nosym/',string_root);
dir_data_star = sprintf('/%s/rangan/dir_cryoem/dir_trpv1',string_root);
Pixel_Spacing = 1.2156; %<-- in angstroms, from https://www.ebi.ac.uk/pdbe/emdb/empiar/entry/10005/ ;

h2d_ = @(kd) (2*pi)^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^2;
dh2d_ = @(kd) (2*pi)^2*0.5*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) (2*pi)^3*3*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^3;
dh3d_ = @(kd) (2*pi)^3*9*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% First load trpv1c molecule on x_u grid. ;
%%%%%%%%;
fname_dims = sprintf('%s/dims',dir_data_marina);
tmp_ = textread(fname_dims); n_x_u = tmp_(1); n_image = tmp_(2); clear tmp_;
fname_mat = sprintf('%s_mat/a_x_u_pack_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_density = sprintf('%s/density_clean',dir_data_marina);
b_x_u_load_ = textread(fname_density); b_x_u_load_ = reshape(b_x_u_load_(:,1),n_x_u,n_x_u,n_x_u);
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_res = 64;
n_pack = n_x_u/x_u_res;
pack_row_ij_ = zeros(x_u_res,1);
pack_col_ij_ = zeros(x_u_res,1);
pack_val_ij_ = zeros(x_u_res,1);
na=0;
for nx_u=0:n_x_u-1;
pack_row_ij_(1+na) = 1+nx_u;
pack_col_ij_(1+na) = 1+floor(nx_u/n_pack);
pack_val_ij_(1+na) = 1/n_pack;
na=na+1;
end;%for nx_u=0:n_x_u-1;
x_u_pack_ = sparse(pack_row_ij_,pack_col_ij_,pack_val_ij_,n_x_u,x_u_res);
b_x_u_pack_ = reshape(b_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
b_x_u_pack_ = reshape(permute(reshape(b_x_u_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
b_x_u_pack_ = reshape(permute(reshape(b_x_u_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
b_x_u_pack_ = permute(reshape(b_x_u_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = x_u_res^3;
X_u_weight_ = (2*x_p_r_max/x_u_res)^3;
%%%%%%%%;
% Now compare to the emd_5778.mrc volume. ;
%%%%%%%%;
fname_emd = sprintf('%s/emd_5778.mrc',dir_data_star);
a_x_u_load_ = cast(ReadMRC(fname_emd),'double');
a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_x_u_load_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1); isosurface_f_x_c_0(a_x_u_load_,98.5); title('a original');
subplot(2,2,2); isosurface_f_x_c_0(a_x_u_pack_,98.5); title('a packed');
subplot(2,2,3); isosurface_f_x_c_0(b_x_u_load_,98.5); title('b original');
subplot(2,2,4); isosurface_f_x_c_0(b_x_u_pack_,98.5); title('b packed');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Note that the correlation is not that high!. ;
%%%%%%%%;
disp(sprintf(' %% corr(a_x_u_load_,b_x_u_load_) %0.4f',corr(a_x_u_load_(:),b_x_u_load_(:))));
disp(sprintf(' %% corr(a_x_u_pack_,b_x_u_pack_) %0.4f',corr(a_x_u_pack_(:),b_x_u_pack_(:))));
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','x_u_res','n_pack','pack_row_ij_','pack_col_ij_','pack_val_ij_','x_u_pack_' ...
     ,'x_u_0_','x_u_1_','x_u_2_','X_u_0_','X_u_1_','X_u_2_','n_X_u','X_u_weight_','a_x_u_pack_' ...
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


%%%%%%%%;
% Now load images. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/M_k_p__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
grid_x_u_ = linspace(-1,+1,n_x_u+1); grid_x_u_ = x_p_r_max*grid_x_u_(1:end-1);
fname_image_mda = sprintf('%s/images_mda',dir_data_marina);
N_x_c___ = MDA_read_r8(fname_image_mda);
%%%%%%%%;
% Now read images from the mrc-stack. ;
%%%%%%%%;
fname_mrc = sprintf('%s/tv1.mrcs',dir_data_star);
M_x_c___ = cast(ReadMRC(fname_mrc,4081,1024),'double');
%%%%%%%%;
% They seem to match reasonably well ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_x_u__vs_N_x_u__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figsml;
Mlim_ = [min(M_x_c___(:)),max(M_x_c___(:))];
plot(N_x_c___(:),M_x_c___(:),'.',Mlim_,Mlim_,'k-');
xlim(Mlim_);ylim(Mlim_);grid on;
xlabel('images');ylabel('star');
title(sprintf('corr %0.4f',corr(N_x_c___(:),M_x_c___(:))));
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
disp(sprintf(' %% corr(N_x_c___,M_x_c___) %0.4f',corr(N_x_c___(:),M_x_c___(:))));
clear N_x_c___ ;
n_M = size(M_x_c___,3);
%%%%%%%%;
% Now examine image-centroids. ;
%%%%%%%%;
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
M_abs_x_c_0_avg_ = zeros(n_M,1);
M_abs_x_c_1_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
M_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nM)));
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nM) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nM) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
end;%for nM=0:n_M-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_abs_x_c_avg_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figmed;
plot(M_abs_x_c_0_avg_,M_abs_x_c_1_avg_,'.');
xlabel('M_abs_x_c_0_avg_','Interpreter','none');
ylabel('M_abs_x_c_1_avg_','Interpreter','none');
axis equal;
title('M_abs_x_c_','Interpreter','none');
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
M_k_p__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_x_c_ = squeeze(M_x_c___(:,:,1+nM));
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*M_abs_x_c_0_avg_(1+nM),-1*M_abs_x_c_1_avg_(1+nM));
M_k_p__(:,1+nM) = M_k_p_;
end;%for nM=0:n_M-1;
%%%%%%%%;
save(fname_mat ...
     ,'n_M','x_c_0_','x_c_1_','x_c_0__','x_c_1__' ...
     ,'M_abs_x_c_0_avg_','M_abs_x_c_1_avg_' ...
     ,'M_k_p__' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% Now convert images to M_k_q__. ;
%%%%%%%%;
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
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

%{
%%%%%%%%;
% Now load euler angles from original star file. ;
%%%%%%%%;
fname_euler_angle = sprintf('%s/euler_angles',dir_data_marina);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
euler_angle_load__ = cell2mat(euler_angle_load_);
fclose(fp);
%%%%%%%%;
fname_star = sprintf('%s/tv1_relion_data.star',dir_data_star);
[blockNames,blockData,ok]=ReadStarFile(fname_star);
tmp_ij_ = 4081 + [0:n_M-1];
euler_angle_star__ = [ ...
 blockData{1}.rlnAngleRot(tmp_ij_) ...
,blockData{1}.rlnAngleTilt(tmp_ij_) ...
,blockData{1}.rlnAnglePsi(tmp_ij_) ...
,blockData{1}.rlnOriginX(tmp_ij_) ...
,blockData{1}.rlnOriginY(tmp_ij_) ...
];
%%%%%%%%;
fname_fig = sprintf('%s_jpg/euler_angle_original_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,3,1);
plot(euler_angle_star__(:,1),euler_angle_load__(1:n_M,1),'.');
xlabel('star');ylabel('load');title('AngleRot');
axis equal; grid on;
subplot(2,3,2);
plot(euler_angle_star__(:,2),euler_angle_load__(1:n_M,2),'.');
xlabel('star');ylabel('load');title('AngleTilt');
axis equal; grid on;
subplot(2,3,3);
plot(euler_angle_star__(:,3),euler_angle_load__(1:n_M,3),'.');
xlabel('star');ylabel('load');title('AnglePsi');
axis equal; grid on;
subplot(2,3,4);
plot(euler_angle_star__(:,4),euler_angle_load__(1:n_M,4),'.');
xlabel('star');ylabel('load');title('OriginX');
axis equal; grid on;
subplot(2,3,5);
plot(euler_angle_star__(:,5),euler_angle_load__(1:n_M,5),'.');
xlabel('star');ylabel('load');title('OriginY');
axis equal; grid on;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Appears as though 4-fold symmetry was used to calculate star. ;
%%%%%%%%;
 %}

%%%%%%%%;
% Now load euler angles from marina star file. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/star_blockData.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_star = sprintf('%s/tv1_run002_data.star',dir_data_star);
[star_blockNames,star_blockData,ok]=ReadStarFile(fname_star);
save(fname_mat ...
     ,'fname_star','star_blockNames','star_blockData' ...
     );
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_euler_angle = sprintf('%s/euler_angles',dir_data_marina);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
euler_angle_load__ = cell2mat(euler_angle_load_);
fclose(fp);
%%%%%%%%;
tmp_ij_ = 1 + [0:n_M-1];
euler_angle_star__ = [ ...
 star_blockData{1}.rlnAngleRot(tmp_ij_) ...
,star_blockData{1}.rlnAngleTilt(tmp_ij_) ...
,star_blockData{1}.rlnAnglePsi(tmp_ij_) ...
,star_blockData{1}.rlnOriginX(tmp_ij_) ...
,star_blockData{1}.rlnOriginY(tmp_ij_) ...
];
%%%%%%%%;
fname_fig = sprintf('%s_jpg/euler_angle_marina_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,3,1);
plot(euler_angle_star__(:,1),euler_angle_load__(1:n_M,1),'.');
xlabel('star');ylabel('load');title('AngleRot');
axis equal; grid on;
subplot(2,3,2);
plot(euler_angle_star__(:,2),euler_angle_load__(1:n_M,2),'.');
xlabel('star');ylabel('load');title('AngleTilt');
axis equal; grid on;
subplot(2,3,3);
plot(euler_angle_star__(:,3),euler_angle_load__(1:n_M,3),'.');
xlabel('star');ylabel('load');title('AnglePsi');
axis equal; grid on;
subplot(2,3,4);
plot(euler_angle_star__(:,4),euler_angle_load__(1:n_M,4),'.');
xlabel('star');ylabel('load');title('OriginX');
axis equal; grid on;
subplot(2,3,5);
plot(euler_angle_star__(:,5),euler_angle_load__(1:n_M,5),'.');
xlabel('star');ylabel('load');title('OriginY');
axis equal; grid on;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% These match. ;
%%%%%%%%;
disp(sprintf(' %% euler_angle_star__ vs euler_angle_load__: %0.16f',fnorm(euler_angle_star__ - euler_angle_load__(1:n_M,:))/fnorm(euler_angle_star__)));
%%%%%%%%;
euler_angle_convert__ = zeros(3,n_M); %<-- [polar_a,azimu_b,gamma_z]. ;
for nimage=0:n_M-1;
euler_angle_convert__(:,1+nimage) = convert_euler_relion_to_marina([euler_angle_star__(1+nimage,1+0),euler_angle_star__(1+nimage,1+1),euler_angle_star__(1+nimage,1+2)]);
end;%for nimage=0:n_M-1;
euler_polar_a_star_ = transpose(euler_angle_convert__(1+0,1:n_M));
euler_azimu_b_star_ = transpose(euler_angle_convert__(1+1,1:n_M));
euler_gamma_z_star_ = transpose(euler_angle_convert__(1+2,1:n_M));
image_delta_x_star_ = euler_angle_star__(1:n_M,1+3)*(2/n_x_u);
image_delta_y_star_ = euler_angle_star__(1:n_M,1+4)*(2/n_x_u);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/euler_angle_star_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
dscale = 2.5;
delta_sigma_star = std([image_delta_x_star_;image_delta_y_star_],1);
figure(1+nf);nf=nf+1;clf; set(gcf,'Position',1+[0,0,512*3,512]);
c2d_euler__ = colormap_polar_a_azimu_b_2d(+euler_polar_a_star_(1:n_M),+euler_azimu_b_star_(1:n_M),0.35);
markersize_euler = 50;
c2d_delta__ = colormap_gaussian_2d(+image_delta_x_star_(1:n_M),+image_delta_y_star_(1:n_M),dscale*delta_sigma_star,0.35);
markersize_delta = 50;
%%%%%%%%;
subplot(1,3,1+[0,1]);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_azimu_b_star_(1:n_M),1*pi-euler_polar_a_star_(1:n_M),markersize_euler,c2d_euler__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('star view');
%%%%%%%%;
subplot(1,3,1+2);
hold on;
patch(dscale*delta_sigma_star*[-1;+1;+1;-1],dscale*delta_sigma_star*[-1;-1;+1;+1],'k');
scatter(+image_delta_x_star_(1:n_M),+image_delta_y_star_(1:n_M),markersize_delta,c2d_delta__(1:n_M,:),'filled');
hold off;
axisnotick;
axis(dscale*delta_sigma_star*[-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('star delta');
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
% Now load ctf parameters. ;
%%%%%%%%;
fname_num_ctf = sprintf('%s/num_ctf',dir_data_marina);
n_ctf = textread(fname_num_ctf);
fname_ctf_idx = sprintf('%s/ctf_idx',dir_data_marina);
CTF_idx_ = textread(fname_ctf_idx)-1;
fname_mscope_params = sprintf('%s/mscope_params',dir_data_marina);
mscope_params = textread(fname_mscope_params);
CTF_Voltage_kV_ = mscope_params(1)*ones(n_ctf,1);
CTF_Spherical_Aberration_ = mscope_params(2)*ones(n_ctf,1);
CTF_Detector_Pixel_Size_ = mscope_params(3)*ones(n_ctf,1);
CTF_Amplitude_Contrast_ = mscope_params(4)*ones(n_ctf,1);
fname_ctf_params = sprintf('%s/ctf_params',dir_data_marina);
ctf_params_ = textread(fname_ctf_params);
CTF_Defocus_U_ = ctf_params_(:,1);
CTF_Defocus_V_ = ctf_params_(:,2);
CTF_Defocus_Angle_ = ctf_params_(:,3);
%%%%%%%%;
disp(sprintf(' %% CTF_Voltage_kV_ vs star_blockData{1}.rlnVoltage: %0.16f',fnorm(CTF_Voltage_kV_(1+CTF_idx_) - star_blockData{1}.rlnVoltage)));
disp(sprintf(' %% CTF_Spherical_Aberration_ vs star_blockData{1}.rlnSphericalAberration: %0.16f',fnorm(CTF_Spherical_Aberration_(1+CTF_idx_) - star_blockData{1}.rlnSphericalAberration)));
disp(sprintf(' %% CTF_Amplitude_Contrast_ vs star_blockData{1}.rlnAmplitudeContrast: %0.16f',fnorm(CTF_Amplitude_Contrast_(1+CTF_idx_) - star_blockData{1}.rlnAmplitudeContrast)));
disp(sprintf(' %% CTF_Defocus_U_ vs star_blockData{1}.rlnDefocusU: %0.16f',fnorm(CTF_Defocus_U_(1+CTF_idx_) - star_blockData{1}.rlnDefocusU)));
disp(sprintf(' %% CTF_Defocus_V_ vs star_blockData{1}.rlnDefocusV: %0.16f',fnorm(CTF_Defocus_V_(1+CTF_idx_) - star_blockData{1}.rlnDefocusV)));
disp(sprintf(' %% CTF_Defocus_Angle_*180/pi vs star_blockData{1}.rlnDefocusAngle: %0.16f',fnorm(CTF_Defocus_Angle_(1+CTF_idx_)*180/pi - star_blockData{1}.rlnDefocusAngle)/fnorm(CTF_Defocus_Angle_(1+CTF_idx_)*180/pi)));
%%%%%%%%;
% Note the one-to-one correspondence between CTF_idx_ and MicrographName. ;
%%%%%%%%;
u_MicrographName_ = unique(star_blockData{1}.rlnMicrographName); n_u_MicrographName = numel(u_MicrographName_);
for nu_MicrographName=0:n_u_MicrographName-1;
tmp_index_ = efind(star_blockData{1}.rlnMicrographName==u_MicrographName_(1+nu_MicrographName));
assert(numel(unique(CTF_idx_(1+tmp_index_)))==1);
end;%for nu_MicrographName=0:n_u_MicrographName-1;
assert(numel(unique(CTF_idx_))==n_u_MicrographName);
%%%%%%%%;
n_ctf = numel(unique(star_blockData{1}.rlnMicrographName));
CTF_Detector_Pixel_Size_ = mscope_params(3)*ones(n_ctf,1);
CTF_Voltage_kV_ = zeros(n_ctf,1);
CTF_Spherical_Aberration_ = zeros(n_ctf,1);
CTF_Amplitude_Contrast_ = zeros(n_ctf,1);
CTF_Defocus_U_ = zeros(n_ctf,1);
CTF_Defocus_V_ = zeros(n_ctf,1);
CTF_Defocus_Angle_ = zeros(n_ctf,1);
CTF_index_ = zeros(n_M,1);
for nctf=0:n_ctf-1;
MicrographName = u_MicrographName_(1+nctf);
tmp_index_ = efind(star_blockData{1}.rlnMicrographName==MicrographName);
CTF_index_(1+tmp_index_) = nctf;
tmp_index = min(tmp_index_);
CTF_Voltage_kV_(1+nctf) = star_blockData{1}.rlnVoltage(1+tmp_index);
CTF_Spherical_Aberration_(1+nctf) = star_blockData{1}.rlnSphericalAberration(1+tmp_index);
CTF_Amplitude_Contrast_(1+nctf) = star_blockData{1}.rlnAmplitudeContrast(1+tmp_index);
CTF_Defocus_U_(1+nctf) = star_blockData{1}.rlnDefocusU(1+tmp_index);
CTF_Defocus_V_(1+nctf) = star_blockData{1}.rlnDefocusV(1+tmp_index);
CTF_Defocus_Angle_(1+nctf) = star_blockData{1}.rlnDefocusAngle(1+tmp_index);
end;%for nctf=0:n_ctf-1;
%%%%%%%%;

%%%%%%%%;
fname_mat = sprintf('%s_mat/CTF_k_p__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
CTF_k_p__ = zeros(n_w_sum,n_ctf);
CTF_uni_k_p__ = zeros(n_w_max*n_k_p_r,n_ctf);
for nctf=0:n_ctf-1;
if (mod(nctf,100)==0); disp(sprintf(' %% nctf %d/%d',nctf,n_ctf)); end;
CTF_Spherical_Aberration = CTF_Spherical_Aberration_(1+nctf);% spherical aberation of the lens in mm ;
CTF_Spherical_Aberration=CTF_Spherical_Aberration*(10.0d0^7.0d0);% convert into Angstroms ;
CTF_Voltage_kV = CTF_Voltage_kV_(1+nctf);% voltage in kVolts ;
CTF_Voltage_1V=CTF_Voltage_kV*1000.0 ;% convert into Volts ;
CTF_lambda = 12.2643247/sqrt(CTF_Voltage_1V+CTF_Voltage_1V^2*0.978466d-6);% electron wavelength in Angstroms ;
CTF_Defocus_U = CTF_Defocus_U_(1+nctf);% defocus values (in Angstroms) ;
CTF_Defocus_V = CTF_Defocus_V_(1+nctf);% defocus values (in Angstroms) ;
CTF_Defocus_Angle = CTF_Defocus_Angle_(1+nctf);% angle of astigmatism ;
%CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ; %<-- already in radians! make sure not to convert twice!;
CTF_Amplitude_Contrast = CTF_Amplitude_Contrast_(1+nctf);% CTF_Amplitude Contrast ;
tmp_w1=sqrt(1.0d0-CTF_Amplitude_Contrast^2);% weights for the amplitude and phase contrasts in CTF ;
tmp_w2=CTF_Amplitude_Contrast;% weights for the amplitude and phase contrasts in CTF ;
%  CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size/CTF_Magnification;
CTF_Object_Pixel_Size = CTF_Detector_Pixel_Size_(1+nctf);% pixel size of the scanner in physical space (not magnified) in Angstroms ;
% Note that here we use n_x_u instead of n_x_u_ori, since the magnification has been increased accordingly. ;
CTF_lambda_per_box = CTF_lambda/(n_x_u*CTF_Object_Pixel_Size);% n_x_u_max*CTF_Object_Pixel_Size is the box size in Angstroms ;
%  call envelope_fxn(ngridr,xnodesr/pi,D,envelope);
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_(1+nk)-1;
tmp_theta = (2.0d0*pi*nw)/n_w_(1+nk);
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle,CTF_lambda_per_box,tmp_k_c_1/pi,tmp_k_c_2/pi);
CTF_k_p__(1+na,1+nctf) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_(1+nk)-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
na=0;
for nk = 0:n_k_p_r-1;
for nw=0:n_w_max-1;
tmp_theta = (2.0d0*pi*nw)/n_w_max;
tmp_k_c_1 = (2.0d0*pi)*k_p_r_(1+nk)*cos(tmp_theta);
tmp_k_c_2 = (2.0d0*pi)*k_p_r_(1+nk)*sin(tmp_theta);
tmp_ctf_value = niko_ctf(CTF_Spherical_Aberration,CTF_lambda,tmp_w1,tmp_w2,CTF_Defocus_U,CTF_Defocus_V,CTF_Defocus_Angle,CTF_lambda_per_box,tmp_k_c_1/pi,tmp_k_c_2/pi);
CTF_uni_k_p__(1+na,1+nctf) = -tmp_ctf_value;
na = na+1;
end;%for nw=0:n_w_max-1;
end;%for nk = 0:n_k_p_r-1;
%%%%;
end;%for nctf=0:n_ctf-1;
%%%%%%%%;
CTF_k_p_r__ = zeros(n_k_p_r,n_ctf);
for nctf=0:n_ctf-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_k_p__(1+tmp_index_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_ctf-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
SCTF_ = svd(CTF_k_p_r__(:,1+CTF_index_(1:n_M)));
n_CTF_rank = min(efind(cumsum(SCTF_,'reverse')/sum(SCTF_)<1e-2));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_index_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;
% Now plot out some of the CTF-functions for varying anisotropy. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/CTF_k_p_sample__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
[tmp_anisotropy_,CTF_anisotropy_index_] = sort(CTF_Defocus_U_ - CTF_Defocus_V_,'ascend'); CTF_anisotropy_index_ = CTF_anisotropy_index_ - 1;
for nl=0:15-1;
subplot(3,5,1+nl);
tmp_nctf = max(0,min(n_ctf-1,floor(n_ctf*nl/(15-1)))); nctf = CTF_anisotropy_index_(1+tmp_nctf);
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,CTF_k_p__(:,1+nctf),[-1,+1],colormap_beach());
axis image; axisnotick;
title(sprintf('nctf %d anisotropy %0.2f',nctf,tmp_anisotropy_(1+tmp_nctf)));
end;%for nl=0:15-1;
sgtitle(sprintf('CTF_k_p__'),'Interpreter','none');
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
% This depends  on CTF_index_. ;
%%%%%%%%;
tmp_CTF_avg_k_p_ = mean(CTF_k_p__(:,1+CTF_index_(1+(0:n_M-1))),2);
tmp_CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xavg__ = tmp_CTF_avg_k_p_r_ * transpose(tmp_CTF_avg_k_p_r_);
CTF_k_p_r_xcor__ = CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1))) * transpose(CTF_k_p_r__(:,1+CTF_index_(1+(0:n_M-1)))) / n_M;
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
     ,'n_ctf' ...
     ,'CTF_index_' ...
     ,'CTF_Voltage_kV_' ...
     ,'CTF_Spherical_Aberration_' ...
     ,'CTF_Detector_Pixel_Size_' ...
     ,'CTF_Amplitude_Contrast_' ...
     ,'ctf_params_' ...
     ,'CTF_Defocus_U_' ...
     ,'CTF_Defocus_V_' ...
     ,'CTF_Defocus_Angle_' ...
     ,'CTF_k_p__' ...
     ,'CTF_k_p_r__' ...
     ,'n_CTF_rank' ...
     ,'SCTF_','UCTF_kc__','VSCTF_Mc__' ...
     ,'CTF_avg_k_p_' ...
     ,'CTF_avg_k_p_r_' ...
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

flag_sup=0;
if flag_sup;
%%%%%%%%;
% Now convert a larger set of images from the star-file. ;
%%%%%%%%;
n_sup = 8; n_M_sup = n_M*n_sup;
%%%%%%%%;
M_sup_k_p__ = zeros(n_w_sum,n_M_sup);
M_sup_abs_x_c_0_avg_ = zeros(n_M_sup,1);
M_sup_abs_x_c_1_avg_ = zeros(n_M_sup,1);
for nsup=0:n_sup-1;
disp(sprintf(' %% nsup %d/%d',nsup,n_sup));
fname_mrc = sprintf('%s/tv1.mrcs',dir_data_star);
tmp_M_x_c___ = cast(ReadMRC(fname_mrc,4081 + n_M*nsup,n_M),'double');
tmp_M_abs_x_c_0_avg_ = zeros(n_M,1);
tmp_M_abs_x_c_1_avg_ = zeros(n_M,1);
for nM=0:n_M-1;
tmp_M_abs_x_c_ = abs(squeeze(tmp_M_x_c___(:,:,1+nM)));
tmp_M_abs_avg = mean(tmp_M_abs_x_c_,'all');
tmp_M_abs_x_c_0_avg = mean(tmp_M_abs_x_c_/tmp_M_abs_avg.*x_c_0__,'all');
tmp_M_abs_x_c_1_avg = mean(tmp_M_abs_x_c_/tmp_M_abs_avg.*x_c_1__,'all');
tmp_M_abs_x_c_0_avg_(1+nM) = tmp_M_abs_x_c_0_avg;
tmp_M_abs_x_c_1_avg_(1+nM) = tmp_M_abs_x_c_1_avg;
clear tmp_M_abs_x_c_;
end;%for nM=0:n_M-1;
tmp_M_k_p__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
tmp_M_x_c_ = squeeze(tmp_M_x_c___(:,:,1+nM));
tmp_M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
% Now *DO* translate according to tmp_M_abs_x_c_0_avg_ and tmp_M_abs_x_c_1_avg_. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,-1*tmp_M_abs_x_c_0_avg_(1+nM),-1*tmp_M_abs_x_c_1_avg_(1+nM));
tmp_M_k_p__(:,1+nM) = tmp_M_k_p_;
end;%for nM=0:n_M-1;
M_sup_k_p__(:,1+n_M*nsup+[0:n_M-1]) = tmp_M_k_p__;
M_sup_abs_x_c_0_avg_(1+n_M*nsup+[0:n_M-1]) = tmp_M_abs_x_c_0_avg_;
M_sup_abs_x_c_1_avg_(1+n_M*nsup+[0:n_M-1]) = tmp_M_abs_x_c_1_avg_;
clear tmp_*;
end;%for nsup=0:n_sup-1;
%%%%%%%%;
end;%if flag_sup;

%%%%%%%%;
% Now test out principled marching. ;
% First define delta_sigma (for translations). ;
%%%%%%%%;
delta_sigma = 1.0 * std([image_delta_x_star_;image_delta_y_star_]); %<-- no reduction. ;
image_delta_x_star_plus_M_abs_x_c_0_avg_ = image_delta_x_star_ + M_abs_x_c_0_avg_;
image_delta_y_star_plus_M_abs_x_c_1_avg_ = image_delta_y_star_ + M_abs_x_c_1_avg_;
%%%%%%%%;
% Then set up the cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/X__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
[X_3d_xcor_d0__,X_3d_xcor_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xcor__);
[X_3d_xavg_d0__,X_3d_xavg_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xavg__); 
[X_2d_xcor_d0__,X_2d_xcor_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xcor__);
[X_2d_xavg_d0__,X_2d_xavg_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xavg__); 
[X_2d_xcor_d1__,X_2d_xcor_d1_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,CTF_k_p_r_xcor__,delta_sigma);
[X_2d_xavg_d1__,X_2d_xavg_d1_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,CTF_k_p_r_xavg__,delta_sigma);
[X_2d_Suni_d0__,X_2d_Suni_d0_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_S,sparse(1:n_w_sum,1:n_w_sum,CTF_avg_k_p_,n_w_sum,n_w_sum)*S_k_p__);
[X_2d_Memp_d1__,X_2d_Memp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,M_k_p__);
if flag_sup;
[X_2d_Msup_d1__,X_2d_Msup_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M_sup,M_sup_k_p__);
end;%if flag_sup;
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
tmp_euler_polar_a = +euler_angle_convert__(1,1+nM);
tmp_euler_azimu_b = +euler_angle_convert__(2,1+nM);
tmp_euler_gamma_z = -euler_angle_convert__(3,1+nM);
tmp_image_delta_x = +1.0*image_delta_x_star_plus_M_abs_x_c_0_avg_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_star_plus_M_abs_x_c_1_avg_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_index_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_);
X_TM_(1+nM) = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
T_k_p__(:,1+nM) = T_k_p_;
if flag_plot;
figbeach();
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(M(k))');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(M(k))');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(T(k))');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(T(k))');
drawnow;
end;%if flag_plot;
end;%for nM=0:n_M-1;
%%%%%%%%;
[X_2d_Semp_d1__,X_2d_Semp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,T_k_p__);
clear T_k_p__ T_k_p_ ;
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX_3d_xcor_d0__,SX_3d_xcor_d0__,VX_3d_xcor_d0__] = svds(X_3d_xcor_d0__,n_UX_rank); SX_3d_xcor_d0_ = diag(SX_3d_xcor_d0__);
[UX_3d_xavg_d0__,SX_3d_xavg_d0__,VX_3d_xavg_d0__] = svds(X_3d_xavg_d0__,n_UX_rank); SX_3d_xavg_d0_ = diag(SX_3d_xavg_d0__);
[UX_2d_xcor_d0__,SX_2d_xcor_d0__,VX_2d_xcor_d0__] = svds(X_2d_xcor_d0__,n_UX_rank); SX_2d_xcor_d0_ = diag(SX_2d_xcor_d0__);
[UX_2d_xavg_d0__,SX_2d_xavg_d0__,VX_2d_xavg_d0__] = svds(X_2d_xavg_d0__,n_UX_rank); SX_2d_xavg_d0_ = diag(SX_2d_xavg_d0__);
[UX_2d_xcor_d1__,SX_2d_xcor_d1__,VX_2d_xcor_d1__] = svds(X_2d_xcor_d1__,n_UX_rank); SX_2d_xcor_d1_ = diag(SX_2d_xcor_d1__);
[UX_2d_xavg_d1__,SX_2d_xavg_d1__,VX_2d_xavg_d1__] = svds(X_2d_xavg_d1__,n_UX_rank); SX_2d_xavg_d1_ = diag(SX_2d_xavg_d1__);
[UX_2d_Suni_d0__,SX_2d_Suni_d0__,VX_2d_Suni_d0__] = svds(X_2d_Suni_d0__,n_UX_rank); SX_2d_Suni_d0_ = diag(SX_2d_Suni_d0__);
[UX_2d_Semp_d1__,SX_2d_Semp_d1__,VX_2d_Semp_d1__] = svds(X_2d_Semp_d1__,n_UX_rank); SX_2d_Semp_d1_ = diag(SX_2d_Semp_d1__);
[UX_2d_Memp_d1__,SX_2d_Memp_d1__,VX_2d_Memp_d1__] = svds(X_2d_Memp_d1__,n_UX_rank); SX_2d_Memp_d1_ = diag(SX_2d_Memp_d1__);
if flag_sup;
[UX_2d_Msup_d1__,SX_2d_Msup_d1__,VX_2d_Msup_d1__] = svds(X_2d_Msup_d1__,n_UX_rank); SX_2d_Msup_d1_ = diag(SX_2d_Msup_d1__);
end;%if flag_sup;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','n_UX_rank' ...
     ,'X_3d_xcor_d0__','X_3d_xcor_d0_weight_r_','UX_3d_xcor_d0__','SX_3d_xcor_d0_','VX_3d_xcor_d0__' ...
     ,'X_3d_xavg_d0__','X_3d_xavg_d0_weight_r_','UX_3d_xavg_d0__','SX_3d_xavg_d0_','VX_3d_xavg_d0__' ...
     ,'X_2d_xcor_d0__','X_2d_xcor_d0_weight_r_','UX_2d_xcor_d0__','SX_2d_xcor_d0_','VX_2d_xcor_d0__' ...
     ,'X_2d_xavg_d0__','X_2d_xavg_d0_weight_r_','UX_2d_xavg_d0__','SX_2d_xavg_d0_','VX_2d_xavg_d0__' ...
     ,'X_2d_xcor_d1__','X_2d_xcor_d1_weight_r_','UX_2d_xcor_d1__','SX_2d_xcor_d1_','VX_2d_xcor_d1__' ...
     ,'X_2d_xavg_d1__','X_2d_xavg_d1_weight_r_','UX_2d_xavg_d1__','SX_2d_xavg_d1_','VX_2d_xavg_d1__' ...
     ,'X_2d_Suni_d0__','X_2d_Suni_d0_weight_r_','UX_2d_Suni_d0__','SX_2d_Suni_d0_','VX_2d_Suni_d0__' ...
     ,'X_2d_Semp_d1__','X_2d_Semp_d1_weight_r_','UX_2d_Semp_d1__','SX_2d_Semp_d1_','VX_2d_Semp_d1__' ...
     ,'X_2d_Memp_d1__','X_2d_Memp_d1_weight_r_','UX_2d_Memp_d1__','SX_2d_Memp_d1_','VX_2d_Memp_d1__' ...
     ,'X_TM_' ...
     );
if flag_sup;
save(fname_mat ...
     ,'X_2d_Msup_d1__','X_2d_Msup_d1_weight_r_','UX_2d_Msup_d1__','SX_2d_Msup_d1_','VX_2d_Msup_d1__' ...
     ,'-append' ...
     );
end;%if flag_sup;
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
e_c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_convert__(1,1:n_M),+euler_angle_convert__(2,1:n_M),0.35);
d_c2d__ = colormap_gaussian_2d(image_delta_x_star_plus_M_abs_x_c_0_avg_,image_delta_y_star_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
nc_beach_ = max(0,min(n_c_beach-1,floor(n_c_beach*X_TM_)));
n_h = 64; lh_lim_ = [0,10];
figure(1+nf);nf=nf+1;
figbeach();
figbig;
subplot(2,2,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_convert__(2,1:n_M),1*pi-euler_angle_convert__(1,1:n_M),markersize_use,e_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
subplot(2,2,2);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_convert__(2,1:n_M),1*pi-euler_angle_convert__(1,1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view (correlation)');
subplot(2,2,3);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(image_delta_x_star_plus_M_abs_x_c_0_avg_(1:n_M),image_delta_y_star_plus_M_abs_x_c_1_avg_(1:n_M),markersize_use,d_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
subplot(2,2,4);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(image_delta_x_star_plus_M_abs_x_c_0_avg_(1:n_M),image_delta_y_star_plus_M_abs_x_c_1_avg_(1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
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
% Now see how much of the variance associated with one kernel is captured by the principal-modes of the other. ;
% Given two ranks r0 and r1, we can project the rank r1 approximation of X1__ ;
% given by Y1__ = \approx U1__(:,1:r1)*S1__(1:r1,1:r1)*transpose(V1__(:,1:r1)) ;
% onto the first r0 principal-components of X0__, producing: ;
% X0_from_Y1__ = transpose(U0__(:,1:r0)) * Y1 * V0__(:,1:r0) ;
% Then measure the spectrum of this projection: ;
% svd_X0_from_Y1_ = svds(X0_vrom_X1__,r1), ;
% and then finally compare this spectrum to the original spectrum of Y1. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_2d_xcor_d0_vs_2d_Memp_d1_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf; hold on;
%%%%;
total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__ = zeros(n_UX_rank,n_UX_rank);
total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__ = zeros(n_UX_rank,n_UX_rank);
for tmp_r0=0:n_UX_rank-1;
for tmp_r1=0:n_UX_rank-1;
total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_Memp_d1__(:,1:1+tmp_r1)) * UX_2d_xcor_d0__(:,1:1+tmp_r0) * diag(SX_2d_xcor_d0_(1:1+tmp_r0)) * transpose(VX_2d_xcor_d0__(:,1:1+tmp_r0)) * VX_2d_Memp_d1__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_xcor_d0_(1:1+tmp_r0));
total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_xcor_d0__(:,1:1+tmp_r1)) * UX_2d_Memp_d1__(:,1:1+tmp_r0) * diag(SX_2d_Memp_d1_(1:1+tmp_r0)) * transpose(VX_2d_Memp_d1__(:,1:1+tmp_r0)) * VX_2d_xcor_d0__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_Memp_d1_(1:1+tmp_r0));
end;%for tmp_r1=0:n_UX_rank-1;
end;%for tmp_r0=0:n_UX_rank-1;
%%%%;
plot(1:n_UX_rank,diag(total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__),'ro-','LineWidth',2);
plot(1:n_UX_rank,diag(total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__),'bo-','LineWidth',2);
hold off;
legend({'total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__','total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__'},'Location','SouthEast');
xlabel('rank'); ylabel('shared variance');
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
if flag_sup;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_2d_xcor_d0_vs_2d_Msup_d1_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf; hold on;
%%%%;
total_variance_X_2d_Msup_d1_from_X_2d_xcor_d0__ = zeros(n_UX_rank,n_UX_rank);
total_variance_X_2d_xcor_d0_from_X_2d_Msup_d1__ = zeros(n_UX_rank,n_UX_rank);
for tmp_r0=0:n_UX_rank-1;
for tmp_r1=0:n_UX_rank-1;
total_variance_X_2d_Msup_d1_from_X_2d_xcor_d0__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_Msup_d1__(:,1:1+tmp_r1)) * UX_2d_xcor_d0__(:,1:1+tmp_r0) * diag(SX_2d_xcor_d0_(1:1+tmp_r0)) * transpose(VX_2d_xcor_d0__(:,1:1+tmp_r0)) * VX_2d_Msup_d1__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_xcor_d0_(1:1+tmp_r0));
total_variance_X_2d_xcor_d0_from_X_2d_Msup_d1__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_xcor_d0__(:,1:1+tmp_r1)) * UX_2d_Msup_d1__(:,1:1+tmp_r0) * diag(SX_2d_Msup_d1_(1:1+tmp_r0)) * transpose(VX_2d_Msup_d1__(:,1:1+tmp_r0)) * VX_2d_xcor_d0__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_Msup_d1_(1:1+tmp_r0));
end;%for tmp_r1=0:n_UX_rank-1;
end;%for tmp_r0=0:n_UX_rank-1;
%%%%;
plot(1:n_UX_rank,diag(total_variance_X_2d_Msup_d1_from_X_2d_xcor_d0__),'ro-','LineWidth',2);
plot(1:n_UX_rank,diag(total_variance_X_2d_xcor_d0_from_X_2d_Msup_d1__),'bo-','LineWidth',2);
hold off;
legend({'total_variance_X_2d_Msup_d1_from_X_2d_xcor_d0__','total_variance_X_2d_xcor_d0_from_X_2d_Msup_d1__'},'Location','SouthEast');
xlabel('rank'); ylabel('shared variance');
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
fname_fig = sprintf('%s_jpg/X_2d_Memp_d1_vs_2d_Msup_d1_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf; hold on;
%%%%;
total_variance_X_2d_Msup_d1_from_X_2d_Memp_d1__ = zeros(n_UX_rank,n_UX_rank);
total_variance_X_2d_Memp_d1_from_X_2d_Msup_d1__ = zeros(n_UX_rank,n_UX_rank);
for tmp_r0=0:n_UX_rank-1;
for tmp_r1=0:n_UX_rank-1;
total_variance_X_2d_Msup_d1_from_X_2d_Memp_d1__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_Msup_d1__(:,1:1+tmp_r1)) * UX_2d_Memp_d1__(:,1:1+tmp_r0) * diag(SX_2d_Memp_d1_(1:1+tmp_r0)) * transpose(VX_2d_Memp_d1__(:,1:1+tmp_r0)) * VX_2d_Msup_d1__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_Memp_d1_(1:1+tmp_r0));
total_variance_X_2d_Memp_d1_from_X_2d_Msup_d1__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_Memp_d1__(:,1:1+tmp_r1)) * UX_2d_Msup_d1__(:,1:1+tmp_r0) * diag(SX_2d_Msup_d1_(1:1+tmp_r0)) * transpose(VX_2d_Msup_d1__(:,1:1+tmp_r0)) * VX_2d_Memp_d1__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_Msup_d1_(1:1+tmp_r0));
end;%for tmp_r1=0:n_UX_rank-1;
end;%for tmp_r0=0:n_UX_rank-1;
%%%%;
plot(1:n_UX_rank,diag(total_variance_X_2d_Msup_d1_from_X_2d_Memp_d1__),'ro-','LineWidth',2);
plot(1:n_UX_rank,diag(total_variance_X_2d_Memp_d1_from_X_2d_Msup_d1__),'bo-','LineWidth',2);
hold off;
legend({'total_variance_X_2d_Msup_d1_from_X_2d_Memp_d1__','total_variance_X_2d_Memp_d1_from_X_2d_Msup_d1__'},'Location','SouthEast');
xlabel('rank'); ylabel('shared variance');
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
end;%if flag_sup;

%%%%%%%%;
% Verdict: Msup is basically identical to Memp. ;
%%%%%%%%;

%%%%%%%%;
% Now visualize principal-volumes. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/vis_UX_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
vis_n_UX_rank = 8;
vis_UX_2d_xcor_d0_a_k_Y_quad__ = zeros(n_lm_sum,vis_n_UX_rank);
vis_UX_2d_xcor_d1_a_k_Y_quad__ = zeros(n_lm_sum,vis_n_UX_rank);
vis_UX_2d_Memp_d1_a_k_Y_quad__ = zeros(n_lm_sum,vis_n_UX_rank);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = [n_lm_csum_(1+nk_p_r):n_lm_csum_(1+nk_p_r+1)-1];
tmp_a_k_Y_quad_ = a_k_Y_quad_(1+tmp_index_);
for nUX_rank=0:vis_n_UX_rank-1;
vis_UX_2d_xcor_d0_a_k_Y_quad__(1+tmp_index_,1+nUX_rank) = tmp_a_k_Y_quad_ * UX_2d_xcor_d0__(1+nk_p_r,1+nUX_rank)*X_2d_xcor_d0_weight_r_(1+nk_p_r);
vis_UX_2d_xcor_d1_a_k_Y_quad__(1+tmp_index_,1+nUX_rank) = tmp_a_k_Y_quad_ * UX_2d_xcor_d1__(1+nk_p_r,1+nUX_rank)*X_2d_xcor_d1_weight_r_(1+nk_p_r);
vis_UX_2d_Memp_d1_a_k_Y_quad__(1+tmp_index_,1+nUX_rank) = tmp_a_k_Y_quad_ * UX_2d_Memp_d1__(1+nk_p_r,1+nUX_rank)*X_2d_Memp_d1_weight_r_(1+nk_p_r);
end;%for nUX_rank=0:vis_n_UX_rank-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
vis_UX_2d_xcor_d0_a_k_p_reco__ = zeros(n_k_all,vis_n_UX_rank);
vis_UX_2d_xcor_d1_a_k_p_reco__ = zeros(n_k_all,vis_n_UX_rank);
vis_UX_2d_Memp_d1_a_k_p_reco__ = zeros(n_k_all,vis_n_UX_rank);
for nUX_rank=0:vis_n_UX_rank-1;
tmp_t = tic;
vis_UX_2d_xcor_d0_a_k_p_reco__(:,1+nUX_rank) = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,vis_UX_2d_xcor_d0_a_k_Y_quad__(:,1+nUX_rank));
tmp_t = toc(tmp_t); disp(sprintf(' %% vis_UX_2d_xcor_d0_a_k_Y_quad_ --> vis_UX_2d_xcor_d0_a_k_p_reco_ time %0.2fs',tmp_t));
tmp_t = tic;
vis_UX_2d_xcor_d1_a_k_p_reco__(:,1+nUX_rank) = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,vis_UX_2d_xcor_d1_a_k_Y_quad__(:,1+nUX_rank));
tmp_t = toc(tmp_t); disp(sprintf(' %% vis_UX_2d_xcor_d1_a_k_Y_quad_ --> vis_UX_2d_xcor_d1_a_k_p_reco_ time %0.2fs',tmp_t));
tmp_t = tic;
vis_UX_2d_Memp_d1_a_k_p_reco__(:,1+nUX_rank) = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,vis_UX_2d_Memp_d1_a_k_Y_quad__(:,1+nUX_rank));
tmp_t = toc(tmp_t); disp(sprintf(' %% vis_UX_2d_Memp_d1_a_k_Y_quad_ --> vis_UX_2d_Memp_d1_a_k_p_reco_ time %0.2fs',tmp_t));
end;%for nUX_rank=0:vis_n_UX_rank-1;
%%%%%%%%;
vis_UX_2d_xcor_d0_a_x_u_reco__ = zeros(n_X_u,vis_n_UX_rank);
vis_UX_2d_xcor_d1_a_x_u_reco__ = zeros(n_X_u,vis_n_UX_rank);
vis_UX_2d_Memp_d1_a_x_u_reco__ = zeros(n_X_u,vis_n_UX_rank);
for nUX_rank=0:vis_n_UX_rank-1;
eta = pi/k_p_r_max; tmp_t = tic;
vis_UX_2d_xcor_d0_a_x_u_reco__(:,1+nUX_rank) = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,vis_UX_2d_xcor_d0_a_k_p_reco__(:,1+nUX_rank).*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: vis_UX_2d_xcor_d0_a_x_u_reco_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
vis_UX_2d_xcor_d1_a_x_u_reco__(:,1+nUX_rank) = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,vis_UX_2d_xcor_d1_a_k_p_reco__(:,1+nUX_rank).*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: vis_UX_2d_xcor_d1_a_x_u_reco_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
vis_UX_2d_Memp_d1_a_x_u_reco__(:,1+nUX_rank) = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,vis_UX_2d_Memp_d1_a_k_p_reco__(:,1+nUX_rank).*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: vis_UX_2d_Memp_d1_a_x_u_reco_ time %0.2fs',tmp_t));
end;%for nUX_rank=0:vis_n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'vis_n_UX_rank' ...
     ,'vis_UX_2d_xcor_d0_a_k_Y_quad__','vis_UX_2d_xcor_d0_a_k_p_reco__','vis_UX_2d_xcor_d0_a_x_u_reco__' ...
     ,'vis_UX_2d_xcor_d1_a_k_Y_quad__','vis_UX_2d_xcor_d1_a_k_p_reco__','vis_UX_2d_xcor_d1_a_x_u_reco__' ...
     ,'vis_UX_2d_Memp_d1_a_k_Y_quad__','vis_UX_2d_Memp_d1_a_k_p_reco__','vis_UX_2d_Memp_d1_a_x_u_reco__' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/vis_UX_2d_xcor_d0_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
tmp_ = reshape(real(vis_UX_2d_xcor_d0_a_x_u_reco__(:,vis_n_UX_rank)),[x_u_res,x_u_res,x_u_res]);
tmp_ = tmp_(:,:,1+floor(x_u_res/2));
vis_avg = mean(tmp_(:)); vis_std = std(tmp_(:)); 
vis_min = min(tmp_(:)); vis_max = max(tmp_(:));
vis_lim_ = prctile(tmp_(:),[ 1,99]); vis_lim_ = max(abs(vis_lim_))*[-1,+1]; vis_lim_ = 1.5*vis_lim_;
%vis_val_ = prctile(tmp_(:),[90,95,99]);
vis_val_ = prctile(tmp_(:),[99.5]);
%%%%;
figure(1+nf);nf=nf+1;clf;
for nUX_rank=0:vis_n_UX_rank-1;
subplot(4,vis_n_UX_rank,1+nUX_rank+0*vis_n_UX_rank);
tmp_ = reshape(real(vis_UX_2d_xcor_d0_a_x_u_reco__(:,1+nUX_rank)),[x_u_res,x_u_res,x_u_res]);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_xcor_d0_(1+nUX_rank)/sum(SX_2d_xcor_d0_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+1*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_xcor_d0_(1+nUX_rank)/sum(SX_2d_xcor_d0_)));
tmp_ = reshape(real(sum(vis_UX_2d_xcor_d0_a_x_u_reco__(:,1+[0:nUX_rank]),2)),[x_u_res,x_u_res,x_u_res]);
subplot(4,vis_n_UX_rank,1+nUX_rank+2*vis_n_UX_rank);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_xcor_d0_(1+(0:nUX_rank)))/sum(SX_2d_xcor_d0_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+3*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_xcor_d0_(1+(0:nUX_rank)))/sum(SX_2d_xcor_d0_)));
end;%for nUX_rank=0:vis_n_UX_rank-1;
figbig;
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/vis_UX_2d_xcor_d1_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
tmp_ = reshape(real(vis_UX_2d_xcor_d1_a_x_u_reco__(:,vis_n_UX_rank)),[x_u_res,x_u_res,x_u_res]);
tmp_ = tmp_(:,:,1+floor(x_u_res/2));
vis_avg = mean(tmp_(:)); vis_std = std(tmp_(:)); 
vis_min = min(tmp_(:)); vis_max = max(tmp_(:));
vis_lim_ = prctile(tmp_(:),[ 1,99]); vis_lim_ = max(abs(vis_lim_))*[-1,+1]; vis_lim_ = 1.5*vis_lim_;
%vis_val_ = prctile(tmp_(:),[90,95,99]);
vis_val_ = prctile(tmp_(:),[99.5]);
%%%%;
figure(1+nf);nf=nf+1;clf;
for nUX_rank=0:vis_n_UX_rank-1;
subplot(4,vis_n_UX_rank,1+nUX_rank+0*vis_n_UX_rank);
tmp_ = reshape(real(vis_UX_2d_xcor_d1_a_x_u_reco__(:,1+nUX_rank)),[x_u_res,x_u_res,x_u_res]);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_xcor_d1_(1+nUX_rank)/sum(SX_2d_xcor_d1_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+1*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_xcor_d1_(1+nUX_rank)/sum(SX_2d_xcor_d1_)));
tmp_ = reshape(real(sum(vis_UX_2d_xcor_d1_a_x_u_reco__(:,1+[0:nUX_rank]),2)),[x_u_res,x_u_res,x_u_res]);
subplot(4,vis_n_UX_rank,1+nUX_rank+2*vis_n_UX_rank);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_xcor_d1_(1+(0:nUX_rank)))/sum(SX_2d_xcor_d1_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+3*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_xcor_d1_(1+(0:nUX_rank)))/sum(SX_2d_xcor_d1_)));
end;%for nUX_rank=0:vis_n_UX_rank-1;
figbig;
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/vis_UX_2d_Memp_d1_FIGA',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
tmp_ = reshape(real(vis_UX_2d_Memp_d1_a_x_u_reco__(:,vis_n_UX_rank)),[x_u_res,x_u_res,x_u_res]);
tmp_ = tmp_(:,:,1+floor(x_u_res/2));
vis_avg = mean(tmp_(:)); vis_std = std(tmp_(:)); 
vis_min = min(tmp_(:)); vis_max = max(tmp_(:));
vis_lim_ = prctile(tmp_(:),[ 1,99]); vis_lim_ = max(abs(vis_lim_))*[-1,+1]; vis_lim_ = 1.5*vis_lim_;
%vis_val_ = prctile(tmp_(:),[90,95,99]);
vis_val_ = prctile(tmp_(:),[99.5]);
%%%%;
figure(1+nf);nf=nf+1;clf;
for nUX_rank=0:vis_n_UX_rank-1;
subplot(4,vis_n_UX_rank,1+nUX_rank+0*vis_n_UX_rank);
tmp_ = reshape(real(vis_UX_2d_Memp_d1_a_x_u_reco__(:,1+nUX_rank)),[x_u_res,x_u_res,x_u_res]);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_Memp_d1_(1+nUX_rank)/sum(SX_2d_Memp_d1_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+1*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: r %d (\\sigma %0.2f)',1+nUX_rank,SX_2d_Memp_d1_(1+nUX_rank)/sum(SX_2d_Memp_d1_)));
tmp_ = reshape(real(sum(vis_UX_2d_Memp_d1_a_x_u_reco__(:,1+[0:nUX_rank]),2)),[x_u_res,x_u_res,x_u_res]);
subplot(4,vis_n_UX_rank,1+nUX_rank+2*vis_n_UX_rank);
isosurface_f_x_u_0(tmp_,[],vis_val_);
title(sprintf('R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_Memp_d1_(1+(0:nUX_rank)))/sum(SX_2d_Memp_d1_)));
subplot(4,vis_n_UX_rank,1+nUX_rank+3*vis_n_UX_rank);
imagesc(tmp_(:,:,32),vis_lim_); fig80s; axis image; axisnotick; xlabel('x'); ylabel('y');
title(sprintf('z=0: R %d (\\Sigma %0.2f)',1+nUX_rank,sum(SX_2d_Memp_d1_(1+(0:nUX_rank)))/sum(SX_2d_Memp_d1_)));
end;%for nUX_rank=0:vis_n_UX_rank-1;
figbig;
%%%%;
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
% Decide which score to use. ;
% This will be based on the str_cost. ;
%%%%%%%%;
str_cost = sprintf('Memp');
%%%%%%%%;
if strfind(str_cost,'2d_xcor');
X__ = X_2d_xcor_d1__ ;
X_weight_r_ = X_2d_xcor_d1_weight_r_ ;
UX__ = UX_2d_xcor_d1__ ;
SX__ = diag(SX_2d_xcor_d1_) ;
SX_ = SX_2d_xcor_d1_ ;
VX__ = VX_2d_xcor_d1__ ;
end;%if strfind(str_cost,'2d_xcor');
%%%%%%%%;
if strfind(str_cost,'Memp');
X__ = X_2d_Memp_d1__ ;
X_weight_r_ = X_2d_Memp_d1_weight_r_ ;
UX__ = UX_2d_Memp_d1__ ;
SX__ = diag(SX_2d_Memp_d1_) ;
SX_ = SX_2d_Memp_d1_ ;
VX__ = VX_2d_Memp_d1__ ;
end;%if strfind(str_cost,'Memp');
%%%%%%%%;
fname_mat = sprintf('%s_mat/X_%s__.mat',dir_pm,str_cost);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now form a_CTF_avg_UX_Y_quad__ ;
%%%%%%%%;
a_CTF_avg_UX_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_avg_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) = a_CTF_avg_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','n_UX_rank' ...
     ,'str_cost' ...
     ,'X__','X_weight_r_' ...
     ,'UX__','SX__','SX_','VX__','a_CTF_avg_UX_Y_quad__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/X_%s_A',dir_pm,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;
colormap(colormap_beach());
subplot(1,2,1); imagesc(log10(abs(UX__)),[-3,0]); xlabel('rank'); ylabel('shell'); title('log10(abs(UX)) [-3,0]'); 
subplot(1,2,2); plot(log10(abs(diag(SX__))),'ko'); xlabel('rank'); ylabel('log10(\sigma)'); title('log10(SX)');
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
fname_fig = sprintf('%s_jpg/X_%s_B',dir_pm,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;
%%%%%%%%;
% First set up a tensor-product spherical grid (in k_p_ space). ;
%%%%%%%%;
k_u_res = 64;
k_u_polar_a_ = linspace(0,pi,k_u_res);
k_u_azimu_b_ = linspace(0,2*pi,2*k_u_res);
[K_u_polar_a_,K_u_azimu_b_] = ndgrid(k_u_polar_a_,k_u_azimu_b_); n_K_u = k_u_res*2*k_u_res;
K_u_weight_ = sin(K_u_polar_a_);
%%%%%%%%;
% Now look at the functions on each shell associated with these 'principal-modes'. ;
%%%%%%%%;
n_plot = 6;
%plot_nk_p_r_ = max(1,min(n_k_p_r,round(linspace(1,n_k_p_r,n_plot))));
plot_nk_p_r_ = 0:n_plot-1;
quad_lim_ = 0.5 * abs(a_CTF_avg_UX_Y_quad__(1,1)) * [-1,+1];
for nplot=0:n_plot-1;
nk_p_r = plot_nk_p_r_(1+nplot);
[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_CTF_avg_UX_Y_quad__(:,1+nk_p_r)),k_u_res,2*k_u_res);
subplot(3,n_plot,1 + nplot + 0*n_plot); imagesc(real(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('real nk_p_r: %d, real(quad)',nk_p_r),'Interpreter','none');
subplot(3,n_plot,1 + nplot + 1*n_plot); imagesc(imag(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('imag nk_p_r: %d, imag(quad)',nk_p_r),'Interpreter','none');
subplot(3,n_plot,1 + nplot + 2*n_plot); imagesc( abs(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf(' abs nk_p_r: %d,  abs(quad)',nk_p_r),'Interpreter','none');
end;%for nplot=0:n_plot-1;
colormap(colormap_beach());
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
% Now examine principled-images. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/UX_%s_M_k_p_wnM___.mat',dir_pm,str_cost);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
UX_M_k_p_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
UX_M_k_q_wnM___ = zeros(n_w_max,n_UX_rank,n_M);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_ = UX__(:,1+nUX_rank);
for nM=0:n_M-1;
tmp_M_k_p__ = reshape(M_k_p__(:,1+nM),[n_w_max,n_k_p_r]);
tmp_M_k_q__ = reshape(M_k_q__(:,1+nM),[n_w_max,n_k_p_r]);
UX_M_k_p_wnM___(:,1+nUX_rank,1+nM) = tmp_M_k_p__*tmp_UX_;
UX_M_k_q_wnM___(:,1+nUX_rank,1+nM) = tmp_M_k_q__*tmp_UX_;
end;%for nM=0:n_M-1;
end;%for nUX_rank=0:n_UX_rank-1;
save(fname_mat ...
     ,'UX_M_k_p_wnM___' ...
     ,'UX_M_k_q_wnM___' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/UX_%s_M_k_p_wnM___',dir_pm,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
gamma_z_ = linspace(0,2*pi,n_w_max+1); gamma_z_ = gamma_z_(1:end-1);
UX_M_k_p_wMn___ = permute(UX_M_k_p_wnM___,[1,3,2]);
clim_ = max(abs(UX_M_k_p_wMn___(:,:,1)),[],'all')*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl);
nUX_rank = max(0,min(n_UX_rank-1,floor(n_UX_rank*nl/(15-1))));
plot(gamma_z_,real(UX_M_k_p_wMn___(:,:,1+nUX_rank)),'k-','LineWidth',0.25);
xlim([0,2*pi]);ylim(clim_);
xlabel('gamma');ylabel('UX_M_k_p_','Interpreter','none'); title(sprintf('nUX_rank %d',nUX_rank),'Interpreter','none');
end;%for nl=0:15-1;
clear UX_M_k_p_wMn___;
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
fname_mat = sprintf('%s_mat/c_k_Y_.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_n_order = 5; tmp_n_M = n_M;
tmp_euler_polar_a_ = +euler_polar_a_star_(1+(0:n_M-1));
tmp_euler_azimu_b_ = +euler_azimu_b_star_(1+(0:n_M-1));
tmp_euler_gamma_z_ = -euler_gamma_z_star_(1+(0:n_M-1)); %<-- note sign change. ;
tmp_image_delta_x_ = +image_delta_x_star_plus_M_abs_x_c_0_avg_(1+(0:n_M-1));
tmp_image_delta_y_ = +image_delta_y_star_plus_M_abs_x_c_1_avg_(1+(0:n_M-1));
%%%%%%%%;
tmp_t = tic;
c_k_Y_reco_ = ...
cg_lsq_4( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_M ...
,M_k_p__ ...
,CTF_index_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> c_k_Y_reco_ time %0.2fs',tmp_t));
[ ...
 X_best_reco ...
,X_best_flag_flip ...
] = ...
register_spharm_to_spharm_wigner_wrap_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_k_Y_quad_ ...
,c_k_Y_reco_ ...
);
disp(sprintf(' %% X_best_reco %0.3f flag_flip %d',X_best_reco,X_best_flag_flip));
%%%%%%%%;
tmp_t = tic;
[c_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,c_k_Y_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_reco_ --> c_k_p_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; 
c_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,c_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: c_x_u_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
save(fname_mat ...
     ,'c_k_Y_reco_' ...
     ,'c_k_p_reco_' ...
     ,'c_x_u_reco_' ...
     ,'X_best_reco' ...
     ,'X_best_flag_flip' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/c_x_u_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
isosurface_f_x_u_0(reshape(real(c_x_u_reco_),x_u_res,x_u_res,x_u_res),[90,95,99]);
title(sprintf('c_x_u_: corr %0.4f',X_best_reco),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

flag_compute = 1 & strcmp(platform,'access1');
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
%%%%%%%%;
date_diff_threshold = 0.25;
flag_force_create_mat=0;flag_force_create_tmp=0;
delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [1,2,4,6,8,10,12,14,16];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
% testing: ;
%{

a_k_Y_true_ = a_k_Y_quad_;
euler_polar_a_true_ = +euler_polar_a_star_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_star_(1+(0:n_M-1));
euler_gamma_z_true_ = -euler_gamma_z_star_(1+(0:n_M-1)); %<-- note sign change. ;
image_delta_x_true_ = +image_delta_x_star_plus_M_abs_x_c_0_avg_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_star_plus_M_abs_x_c_1_avg_(1+(0:n_M-1));
%%%%;
ndat_rseed=0;
dat_rseed = dat_rseed_(1+ndat_rseed);
ndelta_r_max_factor=2;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
ndat_n_UX_rank=3;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
ampmut_wrap_wrap_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,dat_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_ctf ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);
%%%%;

%}

a_k_Y_true_ = a_k_Y_quad_;
euler_polar_a_true_ = +euler_polar_a_star_(1+(0:n_M-1));
euler_azimu_b_true_ = +euler_azimu_b_star_(1+(0:n_M-1));
euler_gamma_z_true_ = -euler_gamma_z_star_(1+(0:n_M-1)); %<-- note sign change. ;
image_delta_x_true_ = +image_delta_x_star_plus_M_abs_x_c_0_avg_(1+(0:n_M-1));
image_delta_y_true_ = +image_delta_y_star_plus_M_abs_x_c_1_avg_(1+(0:n_M-1));
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
parameter = struct('type','parameter');
parameter.rseed = dat_rseed;
parameter.delta_r_max = delta_r_max_use;
parameter.delta_r_upb = delta_r_max_upb;
parameter.dir_pm = dir_pm;
ampmut_wrap_wrap_3( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,dat_n_UX_rank ...
,UX__ ...
,X_weight_r_ ...
,n_M ...
,M_k_p__ ...
,n_ctf ...
,n_CTF_rank ...
,CTF_index_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_true_ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_ ...
,euler_gamma_z_true_ ...
,image_delta_x_true_ ...
,image_delta_y_true_ ...
);
%%%%;
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

disp('returning'); return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Collect ampmut runs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

n_iteration = 2*16;
delta_sigma_use = delta_sigma;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
dat_n_UX_rank_ = [1,2,4,6,8,10,12,14,16];  n_dat_n_UX_rank = numel(dat_n_UX_rank_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
ut_align_a_CTF_avg_UX_Y____ = zeros(n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
ut_align_a_k_Y____ = zeros(n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
vt_align_a_CTF_avg_UX_Y____ = zeros(n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
vt_align_a_k_Y____ = zeros(n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
xt_align_a_CTF_avg_UX_Y____ = zeros(2*n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
xt_align_a_k_Y____ = zeros(2*n_iteration,n_dat_rseed,n_delta_r_max_factor,n_dat_n_UX_rank);
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
%%%%;
ut_fname_pre = sprintf('%s_mat/ut%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_ut_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',ut_fname_pre);
tmp_ = load(fname_ut_align_a_CTF_avg_UX_Y_mat);
ut_align_a_CTF_avg_UX_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
xt_align_a_CTF_avg_UX_Y____(0*n_iteration + [1:n_iteration],1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
fname_ut_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',ut_fname_pre);
tmp_ = load(fname_ut_align_a_k_Y_mat);
ut_align_a_k_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
xt_align_a_k_Y____(0*n_iteration + [1:n_iteration],1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
vt_fname_pre = sprintf('%s_mat/vt%.4dn%.2dr%d',dir_pm,floor(1000*delta_r_max_use),dat_n_UX_rank,dat_rseed);
fname_vt_align_a_CTF_avg_UX_Y_mat = sprintf('%s_align_a_CTF_avg_UX_Y_.mat',vt_fname_pre);
tmp_ = load(fname_vt_align_a_CTF_avg_UX_Y_mat);
vt_align_a_CTF_avg_UX_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
xt_align_a_CTF_avg_UX_Y____(1*n_iteration + [1:n_iteration],1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
fname_vt_align_a_k_Y_mat = sprintf('%s_align_a_k_Y_.mat',vt_fname_pre);
tmp_ = load(fname_vt_align_a_k_Y_mat);
vt_align_a_k_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
xt_align_a_k_Y____(1*n_iteration + [1:n_iteration],1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank) = tmp_.X_best_;
%%%%;
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;

fname_ue1t_align_a_CTF_avg_UX_Y_mat = sprintf('%s_mat/ue1t0056n16r0_align_a_CTF_avg_UX_Y_.mat',dir_pm);
tmp_ = load(fname_ue1t_align_a_CTF_avg_UX_Y_mat);
ue1t_align_a_CTF_avg_UX_Y_ = tmp_.X_best_;
fname_ue1t_align_a_k_Y_mat = sprintf('%s_mat/ue1t0056n16r0_align_a_k_Y_.mat',dir_pm);
tmp_ = load(fname_ue1t_align_a_k_Y_mat);
ue1t_align_a_k_Y_ = tmp_.X_best_;

figure(1);clf;figbig;
markersize_use = 8;
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
str_symbol_ = {'o-','^-','s-','p-'};
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
subplot(3,3,1+ndat_n_UX_rank);
hold on;
plot([1,2*n_iteration],X_best_reco + [0,0],'k-');
plot(0*16+0.5+[0,0],[0,1],'k-');
plot(1*16+0.5+[0,0],[0,1],'k-');
plot(2*16+0.5+[0,0],[0,1],'k-');
plot(3*16+0.5+[0,0],[0,1],'k-');
plot(4*16+0.5+[0,0],[0,1],'k-');
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
nc = max(0,min(n_c_80s-1,floor(n_c_80s*ndat_n_UX_rank/(n_dat_n_UX_rank-1))));
str_symbol = str_symbol_{1+ndelta_r_max_factor};
plot(1:2*n_iteration,xt_align_a_k_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank),str_symbol,'MarkerSize',markersize_use,'Color',c_80s__(1+nc,:),'MarkerFaceColor',c_80s__(1+nc,:),'MarkerEdgeColor','k','LineWidth',2);
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;
if (dat_n_UX_rank==16);
plot(1:n_iteration,ue1t_align_a_k_Y_,'ko','MarkerSize',markersize_use,'Color','k');
end;%if (dat_n_UX_rank==16);
hold off;
xlim([0,2*n_iteration+1]);
ylim([0.50,0.90]);
xlabel('iteration');
ylabel('correlation');
grid on;
title(sprintf('n %.2d',dat_n_UX_rank));
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
sgtitle('xt_align_a_k_Y____','Interpreter','none');

figure(2);clf;figbig;
markersize_use = 8;
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
str_symbol_ = {'o-','^-','s-','p-'};
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
subplot(3,3,1+ndat_n_UX_rank);
hold on;
plot(0*16+0.5+[0,0],[0,1],'k-');
plot(1*16+0.5+[0,0],[0,1],'k-');
plot(2*16+0.5+[0,0],[0,1],'k-');
plot(3*16+0.5+[0,0],[0,1],'k-');
plot(4*16+0.5+[0,0],[0,1],'k-');
%%%%%%%%;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
nc = max(0,min(n_c_80s-1,floor(n_c_80s*ndat_n_UX_rank/(n_dat_n_UX_rank-1))));
str_symbol = str_symbol_{1+ndelta_r_max_factor};
plot(1:2*n_iteration,xt_align_a_CTF_avg_UX_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank),str_symbol,'MarkerSize',markersize_use,'Color',c_80s__(1+nc,:),'MarkerFaceColor',c_80s__(1+nc,:),'MarkerEdgeColor','k','LineWidth',2);
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;
if (dat_n_UX_rank==16);
plot(1:n_iteration,ue1t_align_a_CTF_avg_UX_Y_,'ko','MarkerSize',markersize_use,'Color','k');
end;%if (dat_n_UX_rank==16);
hold off;
xlim([0,2*n_iteration+1]);
ylim([0.10,0.75]);
xlabel('iteration');
ylabel('correlation');
grid on;
title(sprintf('n %.2d',dat_n_UX_rank));
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
sgtitle('xt_align_a_CTF_avg_UX_Y____','Interpreter','none');

figure(3);clf;figbig;
markersize_use = 8;
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
str_symbol_ = {'o-','^-','s-','p-'};
hold on;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2)); %<-- allow large accumulated translations. ;
for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
dat_n_UX_rank = dat_n_UX_rank_(1+ndat_n_UX_rank);
%%%%;
nc = max(0,min(n_c_80s-1,floor(n_c_80s*ndat_n_UX_rank/(n_dat_n_UX_rank-1))));
str_symbol = str_symbol_{1+ndelta_r_max_factor};
tmp_x_ = xt_align_a_CTF_avg_UX_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank);
tmp_y_ = xt_align_a_k_Y____(:,1+ndat_rseed,1+ndelta_r_max_factor,1+ndat_n_UX_rank);
plot(tmp_x_,tmp_y_,str_symbol,'MarkerSize',markersize_use,'Color',c_80s__(1+nc,:),'MarkerFaceColor',c_80s__(1+nc,:),'MarkerEdgeColor','k','LineWidth',2);
%%%%;
end;%for ndat_n_UX_rank=0:n_dat_n_UX_rank-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;
hold off;
sgtitle('xt_align_a_CTF_avg_UX_Y____ vs xt_align_a_k_Y____','Interpreter','none');
xlim([0.10,0.75]);
ylim([0.50,0.90]);
xlabel('xt_align_a_CTF_avg_UX_Y____','Iteration','none');
ylabel('xt_align_a_k_Y____','Iteration','none');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Relion runs. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First identify mrcs file containing image stack. ;
%%%%%%%%;
fname_mscope_params = sprintf('%s/mscope_params',dir_data_marina);
mscope_params = textread(fname_mscope_params);
Voltage_kV = mscope_params(1);
SphericalAberration = mscope_params(2);
PixelSize = mscope_params(3);
AmplitudeContrast = mscope_params(4);
%%%%%%%%;
fname_image_mrcs = sprintf('%s/tv1.mrcs',dir_data_star);
%%%%%%%%;

n_M_relion_ = 1024*[1,2,4,8]; n_n_M_relion = numel(n_M_relion_);
for nn_M_relion=0:n_n_M_relion-1;
%%%%%%%%%%%%%%%%;
n_M_relion = n_M_relion_(1+nn_M_relion);
jname = sprintf('job_%d',n_M_relion);
tmp_slash = '\';
dir_job = sprintf('%s_mat/%s',dir_relion,jname);
if (~exist(dir_job,'dir')); sprintf(' %% mkdir %s',dir_job); mkdir(dir_job); end;
%%%%%%%%;
% Now write star file linking to the mrcs image stack. ;
%%%%%%%%;
tmp_ImageName_list_ = cell(n_M_relion,1);
tmp_OpticsGroup_list_ = cell(n_M_relion,1);
tmp_DefocusU_list_ = cell(n_M_relion,1);
tmp_DefocusV_list_ = cell(n_M_relion,1);
tmp_DefocusAngle_list_ = cell(n_M_relion,1);
for nM=0:n_M_relion-1;
nctf = CTF_index_(1+nM);
tmp_ImageName_list_{1+nM} = sprintf('%d@../tv1.mrcs',4081+nM); %<-- start at frame 4081. ;
tmp_OpticsGroup_list_{1+nM} = sprintf('1');
tmp_DefocusU_list_{1+nM} = sprintf('%0.16f',CTF_Defocus_U_(1+nctf));
tmp_DefocusV_list_{1+nM} = sprintf('%0.16f',CTF_Defocus_V_(1+nctf));
tmp_DefocusAngle_list_{1+nM} = sprintf('%0.16f',CTF_Defocus_Angle_(1+nctf));
end;%for nM=0:n_M_relion-1;
fname_star = sprintf('%s_mat/tv1_relion_%s.star',dir_relion,jname);
fp = fopen(fname_star,'w');
fprintf(fp,'# version 30001\n');
fprintf(fp,'\n');
fprintf(fp,'data_optics\n');
fprintf(fp,'\n');
fprintf(fp,'loop_\n');
fprintf(fp,'_rlnOpticsGroupName #1\n');
fprintf(fp,'_rlnOpticsGroup #2\n');
fprintf(fp,'_rlnMtfFileName #3\n');
fprintf(fp,'_rlnVoltage #4\n');
fprintf(fp,'_rlnSphericalAberration #5\n');
fprintf(fp,'_rlnAmplitudeContrast #6\n');
fprintf(fp,'_rlnImagePixelSize #7\n');
fprintf(fp,'_rlnImageSize #8\n');
fprintf(fp,'_rlnImageDimensionality #9\n');
fprintf(fp,'opticsGroup1 1 mtf_k2_300kV.star %f %f %f %f %d 2\n',Voltage_kV,SphericalAberration,AmplitudeContrast,PixelSize,n_x_u);
fprintf(fp,'\n');
fprintf(fp,'# version 30001\n');
fprintf(fp,'\n');
fprintf(fp,'data_particles\n');
fprintf(fp,'\n');
fprintf(fp,'loop_\n');
fprintf(fp,'_rlnImageName #3\n');
fprintf(fp,'_rlnOpticsGroup #4\n');
fprintf(fp,'_rlnDefocusU #6\n');
fprintf(fp,'_rlnDefocusV #7\n');
fprintf(fp,'_rlnDefocusAngle #8\n');
for nM=0:n_M_relion-1;
fprintf(fp,'%s %s %s %s %s\n',tmp_ImageName_list_{1+nM},tmp_OpticsGroup_list_{1+nM},tmp_DefocusU_list_{1+nM},tmp_DefocusV_list_{1+nM},tmp_DefocusAngle_list_{1+nM});
end;%for nM=0:n_M_relion-1;
fclose(fp);
%%%%%%%%;
% Now write shell script to call command. ;
%%%%%%%%;
fname_command = sprintf('%s_mat/tv1_relion_%s.sh',dir_relion,jname);
fp = fopen(fname_command,'w');
fprintf(fp,'%s/relion_refine %s\n',dir_relion_bin,tmp_slash);
fprintf(fp,'--o %s/run %s\n',jname,tmp_slash);
fprintf(fp,'--sgd_ini_iter 50 %s\n',tmp_slash);
fprintf(fp,'--sgd_inbetween_iter 200 %s\n',tmp_slash);
fprintf(fp,'--sgd_fin_iter 50 %s\n',tmp_slash);
fprintf(fp,'--sgd_write_iter 20 %s\n',tmp_slash);
fprintf(fp,'--sgd_ini_resol 35 %s\n',tmp_slash);
fprintf(fp,'--sgd_fin_resol 15 %s\n',tmp_slash);
fprintf(fp,'--sgd_ini_subset 100 %s\n',tmp_slash);
fprintf(fp,'--sgd_fin_subset 500 %s\n',tmp_slash);
fprintf(fp,'--sgd %s\n',tmp_slash);
fprintf(fp,'--denovo_3dref %s\n',tmp_slash);
fprintf(fp,'--i %s %s\n',sprintf('tv1_relion_%s.star',jname),tmp_slash);
fprintf(fp,'--ctf %s\n',tmp_slash);
fprintf(fp,'--K 1 %s\n',tmp_slash);
fprintf(fp,'--sym C1 %s\n',tmp_slash);
fprintf(fp,'--flatten_solvent %s\n',tmp_slash);
fprintf(fp,'--zero_mask %s\n',tmp_slash);
fprintf(fp,'--dont_combine_weights_via_disc %s\n',tmp_slash);
fprintf(fp,'--pool 3 %s\n',tmp_slash);
fprintf(fp,'--pad 1 %s\n',tmp_slash);
fprintf(fp,'--skip_gridding %s\n',tmp_slash);
fprintf(fp,'--particle_diameter 400 %s\n',tmp_slash);
fprintf(fp,'--oversampling 1 %s\n',tmp_slash);
fprintf(fp,'--healpix_order 1 %s\n',tmp_slash);
fprintf(fp,'--offset_range 8 %s\n',tmp_slash);
fprintf(fp,'--offset_step 2 %s\n',tmp_slash);
fprintf(fp,'--j 6 %s\n',tmp_slash);
fprintf(fp,'--pipeline_control %s %s\n',jname,tmp_slash);
fprintf(fp,';\n');
fprintf(fp,'\n');
fclose(fp);
disp(sprintf(' %% fname_command: %s',fname_command));
type(fname_command);
%%%%%%%%%%%%%%%%;
end;%for nn_M_relion=0:n_n_M_relion-1;

for jname = {'job_4096','job_2048','job_0'};
%%%%%%%%;
dir_job = sprintf('%s_mat/%s',dir_relion,jname{1});
%%%%%%%%;
% Now examine output. ;
%%%%%%%%;
relion_iter_ = 0:20:300; n_relion_iter = numel(relion_iter_);
a_x_u_relion_pack__ = cell(n_relion_iter,1);
%a_k_p_relion_quad__ = cell(n_relion_iter,1);
for nrelion_iter=0:n_relion_iter-1;
relion_iter = relion_iter_(1+nrelion_iter);
fname_mrc = sprintf('%s/run_it%.3d_class001.mrc',dir_job,relion_iter);
a_x_u_relion_ = cast(ReadMRC(fname_mrc),'double');
a_x_u_relion_pack_ = reshape(a_x_u_relion_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = reshape(permute(reshape(a_x_u_relion_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = reshape(permute(reshape(a_x_u_relion_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = permute(reshape(a_x_u_relion_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
eta = pi/k_p_r_max; 
%a_k_p_relion_quad_ = nufft3d3(n_X_u,X_u_0_(:)*eta,X_u_1_(:)*eta,X_u_2_(:)*eta,a_x_u_relion_pack_(:).*X_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
a_x_u_relion_pack__{1+nrelion_iter} = a_x_u_relion_pack_;
%a_k_p_relion_quad__{1+nrelion_iter} = a_k_p_relion_quad_;
end;%for nrelion_iter=0:n_relion_iter-1;
%%%%%%%%;
fname_fig = sprintf('%s/a_x_u_relion_pack__',dir_job);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig)));
figure(1);clf;figbig;
prows = 3; pcols = 6; np=0;
subplot(prows,pcols,1+np);np=np+1;
isosurface_f_x_c_0(a_x_u_pack_,98.5); axis vis3d; title('ori');
for nrelion_iter=0:n_relion_iter-1;
subplot(prows,pcols,1+np);np=np+1;
relion_iter = relion_iter_(1+nrelion_iter);
isosurface_f_x_c_0(a_x_u_relion_pack__{1+nrelion_iter},98.5); axis vis3d; title(sprintf('relion it%.3d',relion_iter));
end;%for nrelion_iter=0:n_relion_iter-1;
disp(sprintf(' %% writing %s',sprintf('%s.jpg',fname_fig)));
print('-djpeg',sprintf('%s.jpg',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s/a_x_u_relion_pack_',dir_job);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',sprintf('%s.jpg',fname_fig)));
figure(2);clf;figbig;
tmp_prctile_ = [95,98.5,99];
subplot(1,2,1);
isosurface_f_x_c_0(a_x_u_pack_,tmp_prctile_); axis vis3d; title('ori');
subplot(1,2,2);
nrelion_iter = n_relion_iter-1; relion_iter = relion_iter_(1+nrelion_iter);
isosurface_f_x_c_0(a_x_u_relion_pack__{1+nrelion_iter},tmp_prctile_); axis vis3d; title(sprintf('relion it%.3d',relion_iter));
disp(sprintf(' %% writing %s',sprintf('%s.jpg',fname_fig)));
print('-djpeg',sprintf('%s.jpg',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now check each iteration for correlation. ;
%%%%%%%%;
verbose=1;
fname_mat = sprintf('%s/X_relion_.mat',dir_job);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
X_relion_ = zeros(n_relion_iter,1);
for nrelion_iter=n_relion_iter-1:-1:0;
relion_iter = relion_iter_(1+nrelion_iter);
fname_mrc = sprintf('%s/run_it%.3d_class001.mrc',dir_job,relion_iter);
tmp_t = tic();
n_X_u = x_u_res^3; eta = pi/k_p_r_max; 
a_k_p_relion_quad_ = nufft3d3(n_X_u,X_u_0_(:)*eta,X_u_1_(:)*eta,X_u_2_(:)*eta,a_x_u_relion_pack__{1+nrelion_iter}(:).*X_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_k_p_relion_quad_: %0.3fs',tmp_t));
tmp_t = tic();
[a_k_Y_relion_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_relion_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_spharm_1: a_k_Y_relion_quad_: %0.3fs',tmp_t));
[tmp_X_relion_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,a_k_Y_relion_quad_);
[tmp_X_relion_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,a_k_Y_relion_quad_));
if (tmp_X_relion_orig>=tmp_X_relion_flip); flag_flip = 0; end;
if (tmp_X_relion_orig< tmp_X_relion_flip); flag_flip = 1; end;
tmp_X_relion_reco = max(tmp_X_relion_orig,tmp_X_relion_flip);
if (verbose); disp(sprintf(' %% nrelion_iter %.2d/%.2d: tmp_X_relion_reco %0.3f <-- flag_flip %d',nrelion_iter,n_relion_iter,tmp_X_relion_reco,flag_flip)); end;
X_relion_(1+nrelion_iter) = tmp_X_relion_reco;
end;%for nrelion_iter=0:n_relion_iter-1;
save(fname_mat,'relion_iter_','X_relion_');
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
%%%%%%%%;
fname_fig = sprintf('%s/X_relion_',dir_job);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
figure(3);clf;
plot(relion_iter_,X_relion_,'ko-');
xlabel('relion iter');ylabel('correlation');
xlim([0,max(relion_iter_)+1]);ylim([0,1]);
title('X_true','interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for jname = {'job_4096','job_2048','job_0'};


