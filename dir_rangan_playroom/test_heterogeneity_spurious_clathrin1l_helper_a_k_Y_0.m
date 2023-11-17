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

