% Here str_cost = '3d' corrects for volume element, ;
% whereas str_cost = '2d' corrects for area element. ;
% Here str_CTF_factor = 'rescaled' for CTF-rescaled images, ;
% whereas str_CTF_factor = 'original' for original CTF used in images. ;
% 1. the images are *NOT* pretranslated, and ;
% 2. the delta_sigma used to calculate the cost is *NOT* set to 0. ;
% Moreover, the n_w_ used to calculate M_k_p__ is not adaptive. ;
% Now check to see if the experimental images can themselves be used to recapitulate the principal-modes. ;
% Now we investigate the role played by f_rand. ;
% Now we investigate the stability / convergence of alternating minimization. ;
clear;


%platform = 'access1';
platform = 'OptiPlex';
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

flag_recalc = 0;
flag_replot = 0;
str_CTF_factor = 'original';

dir_trunk = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching_rib80s',string_root);

h2d_ = @(kd) (2*pi)^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^2;
dh2d_ = @(kd) (2*pi)^2*0.5*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) (2*pi)^3*3*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^3;
dh3d_ = @(kd) (2*pi)^3*9*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;

%%%%%%%%;
% Quick check of nufft1d3, using k=frequency ;
%%%%%%%%;
gx = @(x_,mu,sg) 1/sqrt(2*pi)/sg * exp(-(x_-mu).^2/2/sg.^2);
gk = @(k_,mu,sg) gx(k_,0,1/sg).*exp(-i*k_*mu)/sg;
x_r_max = 10;
[x_all_,x_weight_all_] = chebpts(128,x_r_max*[-1,+1]);
n_x_all = numel(x_all_);
k_r_max = 15;
[k_all_,k_weight_all_] = chebpts(256,k_r_max*[-1,+1]);
n_k_all = numel(k_all_);
sg = 1.5; mu = 2.5;
f_x_c_form_ = gx(x_all_,mu,sg);
f_k_c_form_ = gk(k_all_,mu,sg);
eta = pi/x_r_max;
f_k_c_quad_ = nufft1d3(n_x_all,x_all_*eta,f_x_c_form_.*x_weight_all_(:),-1,1e-12,n_k_all,k_all_/eta)/sqrt(2*pi);
eta = pi/k_r_max;
f_x_c_quad_ = nufft1d3(n_k_all,k_all_*eta,f_k_c_form_.*k_weight_all_(:),+1,1e-12,n_x_all,x_all_/eta)/sqrt(2*pi);
flag_plot=0;
if flag_plot;
subplot(1,2,1);plot(k_all_,f_k_c_form_,'k.-',k_all_,f_k_c_quad_,'ro-'); xlabel('k'); ylabel('f_k');
subplot(1,2,2);plot(x_all_,f_x_c_form_,'k.-',x_all_,f_x_c_quad_,'ro-'); xlabel('x'); ylabel('f_x');
end;%if flag_plot;
disp(sprintf(' %% nufft1d3: f_k_c_quad error: %0.16f',fnorm(f_k_c_form_-f_k_c_quad_)/fnorm(f_k_c_form_)));
disp(sprintf(' %% nufft1d3: f_x_c_quad error: %0.16f',fnorm(f_x_c_form_-f_x_c_quad_)/fnorm(f_x_c_form_)));
%%%%%%%%;
% Quick check of nufft1d3, using k=wavenumber ;
%%%%%%%%;
gx = @(x_,mu,sg) 1/sqrt(2*pi)/sg * exp(-(x_-mu).^2/2/sg.^2);
gk = @(k_,mu,sg) gx(k_,0,1/sg).*exp(-i*k_*mu)/sg;
x_r_max = 10;
[x_all_,x_weight_all_] = chebpts(128,x_r_max*[-1,+1]);
n_x_all = numel(x_all_);
k_r_max = 15/(2*pi);
[k_all_,k_weight_all_] = chebpts(256,k_r_max*[-1,+1]);
n_k_all = numel(k_all_);
sg = 1.5; mu = 2.5;
f_x_c_form_ = gx(x_all_,mu,sg);
f_k_c_form_ = gk(2*pi*k_all_,mu,sg);
eta = pi/x_r_max;
f_k_c_quad_ = nufft1d3(n_x_all,x_all_*eta,f_x_c_form_.*x_weight_all_(:),-1,1e-12,n_k_all,2*pi*k_all_/eta)/sqrt(2*pi);
eta = pi/k_r_max;
f_x_c_quad_ = nufft1d3(n_k_all,2*pi*k_all_*eta,f_k_c_form_.*2*pi.*k_weight_all_(:),+1,1e-12,n_x_all,x_all_/eta)/sqrt(2*pi);
flag_plot=0;
if flag_plot;
subplot(1,2,1);plot(k_all_,f_k_c_form_,'k.-',k_all_,f_k_c_quad_,'ro-'); xlabel('k'); ylabel('f_k');
subplot(1,2,2);plot(x_all_,f_x_c_form_,'k.-',x_all_,f_x_c_quad_,'ro-'); xlabel('x'); ylabel('f_x');
end;%if flag_plot;
disp(sprintf(' %% nufft1d3: f_k_c_quad error: %0.16f',fnorm(f_k_c_form_-f_k_c_quad_)/fnorm(f_k_c_form_)));
disp(sprintf(' %% nufft1d3: f_x_c_quad error: %0.16f',fnorm(f_x_c_form_-f_x_c_quad_)/fnorm(f_x_c_form_)));

%%%%%%%%;
% First load rib80s molecule on x_u grid. ;
%%%%%%%%;
dir_data = sprintf('/%s/rangan/dir_cryoem/dir_rib80s/new_joakim_relion_run/subs128/data',string_root);
%dir_data = sprintf('/%s/rangan/dir_cryoem/dir_rib80s/new_joakim_relion_run/data1',string_root);
fname_dims = sprintf('%s/dims',dir_data);
n_x_u_ori = 360;
tmp_ = textread(fname_dims); n_x_u = tmp_(1); n_image = tmp_(2); clear tmp_;
fname_density = sprintf('%s/density_run002_subs128',dir_data);
a_x_u_load_ = textread(fname_density); a_x_u_load_ = reshape(a_x_u_load_(:,1),n_x_u,n_x_u,n_x_u);
flag_plot=0; if flag_plot; figure(1); isosurface_f_x_c_0(a_x_u_load_,98.5); end;
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
a_x_u_pack_ = reshape(a_x_u_load_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = reshape(permute(reshape(a_x_u_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
a_x_u_pack_ = permute(reshape(a_x_u_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
flag_plot=0; if flag_plot; figure(1); subplot(1,2,1); isosurface_f_x_c_0(a_x_u_load_,98.5); subplot(1,2,2); isosurface_f_x_c_0(a_x_u_pack_,98.5); end;
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = x_u_res^3;
X_u_weight_ = (2*x_p_r_max/x_u_res)^3;

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_a_k_p_quad_.mat',dir_trunk);
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
,weight_k_all_ ...
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
a_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_quad_.*(2*pi)^3.*weight_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_x_u_reco error: %0.16f',fnorm(a_x_u_pack_(:)-a_x_u_reco_)/fnorm(a_x_u_pack_(:))));
disp(sprintf(' %% at this point one should ensure that a_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
save(fname_mat ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_','k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_','weight_k_all_','weight_shell_k_','n_k_p_r','k_p_r_','weight_3d_k_p_r_','k_c_0_all_','k_c_1_all_','k_c_2_all_','J_node_','J_weight_','J_chebfun_','J_polyval_' ...
     ,'a_k_p_quad_' ...
     ,'a_x_u_reco_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_a_k_p_quad_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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

%{
%%%%%%%%;
% Now see what the function looks like on a uniform x_c_ grid. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_b_x_u_reco__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_x_u_reco__ = zeros(n_X_u,n_k_p_r);
eta = pi/k_p_r_max; 
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+1+nk_p_r)-1;
tmp_n_k_all = numel(1+tmp_index_);
tmp_t = tic;
a_x_u_reco__(:,1+nk_p_r) = nufft3d3(tmp_n_k_all,2*pi*k_c_0_all_(1+tmp_index_)*eta,2*pi*k_c_1_all_(1+tmp_index_)*eta,2*pi*k_c_2_all_(1+tmp_index_)*eta,a_k_p_quad_(1+tmp_index_).*(2*pi)^3.*weight_k_all_(1+tmp_index_),+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nk_p_r %.3d/%.3d nufft3d3: a_x_u_reco__(:,%.3d) time %0.2fs',nk_p_r,n_k_p_r,nk_p_r,tmp_t));
end;%for nk_p_r=0:n_k_p_r-1;
b_x_u_reco__ = cumsum(a_x_u_reco__,2);
save(fname_mat ...
     ,'a_x_u_reco__','b_x_u_reco__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_b_x_u_reco__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
subplot(2,2,1); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,16)),x_u_res,x_u_res,x_u_res),[90,95,99]); title(sprintf('K<=%0.2f',k_p_r_(16)));
subplot(2,2,2); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,24)),x_u_res,x_u_res,x_u_res),[90,95,99]); title(sprintf('K<=%0.2f',k_p_r_(24)));
subplot(2,2,3); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,32)),x_u_res,x_u_res,x_u_res),[90,95,99]); title(sprintf('K<=%0.2f',k_p_r_(32)));
subplot(2,2,4); plot(k_p_r_,corr(a_x_u_pack_(:),real(b_x_u_reco__(:,:))),'ko'); xlabel('max K'); ylabel('correlation'); title('correlation of csum');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
 %}

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_a_k_Y_quad_.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
l_max_upb = 36;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
%l_max_(1+nk_p_r) = 1+ceil(2*pi*k_p_r_(1+nk_p_r));
%l_max_(1+nk_p_r) = 1+ceil(2*pi*k_p_r_(1+nk_p_r));
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
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ --> a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_)));
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
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_a_k_Y_quad_A',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_a_k_Y_quad_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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

%%%%%%%%
% Now generate templates. ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_S_k_p__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
template_k_eq_d = k_eq_d*2;
viewing_k_eq_d = k_eq_d*8;
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
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,a_k_Y_quad_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
);
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
save(fname_mat ...
     ,'template_k_eq_d' ...
     ,'viewing_k_eq_d' ...
     ,'S_k_p__','n_S' ...
     ,'n_w_','weight_2d_k_p_r_','weight_2d_k_all_','n_viewing_all','viewing_azimu_b_all_','viewing_polar_a_all_','n_viewing_polar_a','viewing_polar_a_','n_viewing_azimu_b_','template_k_c_0__','template_k_c_1__','template_k_c_2__' ...
     ,'n_w_max','n_w_sum','n_w_csum_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_S_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
subplot(2,3,1+nplot);
nviewing_all = max(0,min(n_viewing_all-1,round(n_viewing_all*nplot/n_plot)));
imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(S_k_p__(:,1+nviewing_all)),[],colormap_beach()); 
axis image; axisnotick; title(sprintf('S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
end;%for nplot=0:n_plot-1;
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
% Ensure that the templates associated with a single shell correspond (exactly) to the rings within the full templates. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_S_k_p__A',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
nk_p_r = 12;
tmp_S_k_p_0__ = S_k_p__(1+n_w_csum_(1+nk_p_r)+[0:n_w_(1+nk_p_r)-1],:);
colormap(colormap_beach());
subplot(1,2,1); 
imagesc(abs(tmp_S_k_p_0__)); xlabel('nviewing_all','Interpreter','none'); ylabel('ngamma_z','Interpreter','none'); title(sprintf('template ring %d',nk_p_r));
tmp_a_k_Y_quad_ = a_k_Y_quad_(1+n_lm_csum_(1+nk_p_r)+[0:n_lm_(1+nk_p_r)-1]);
tmp_S_k_p_1__ = get_template_0(verbose,1,k_p_r_(1+nk_p_r),k_p_r_max,weight_3d_k_p_r_(1+nk_p_r),l_max_(1+nk_p_r),tmp_a_k_Y_quad_,viewing_k_eq_d,template_k_eq_d);
subplot(1,2,2); 
imagesc(abs(tmp_S_k_p_1__)); xlabel('nviewing_all','Interpreter','none'); ylabel('ngamma_z','Interpreter','none'); title(sprintf('ring template %d',nk_p_r));
disp(sprintf(' %% nk_p_r %d/%d tmp_S_k_p_0__ vs tmp_S_k_p_1__: %0.16f',nk_p_r,n_k_p_r,fnorm(tmp_S_k_p_0__-tmp_S_k_p_1__)/fnorm(tmp_S_k_p_0__)));
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
% For images we will use the same n_w_ as that used by templates. ;
% We will also use the same quadrature weights for integration in 2d. ;
% In addition to the adaptively defined polar grid, ;
% we also generate M_uni_k_p__ on a uniform grid. ;
%%%%%%%%;
n_image_sub = 1024; disp(sprintf(' %% Warning! setting n_image_sub = 1024'));

%%%%%%%%;
% extract ctf function. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_CTF_k_p__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
fname_num_ctf = sprintf('%s/num_ctf',dir_data);
n_ctf = textread(fname_num_ctf);
fname_ctf_idx = sprintf('%s/ctf_idx',dir_data);
CTF_idx_ = textread(fname_ctf_idx)-1;
fname_mscope_params = sprintf('%s/mscope_params',dir_data);
mscope_params = textread(fname_mscope_params);
CTF_Voltage_kV_ = mscope_params(1)*ones(n_ctf,1);
CTF_Spherical_Aberration_ = mscope_params(2)*ones(n_ctf,1);
CTF_Detector_Pixel_Size_ = mscope_params(3)*ones(n_ctf,1);
CTF_Amplitude_Contrast_ = mscope_params(4)*ones(n_ctf,1);
fname_ctf_params = sprintf('%s/ctf_params',dir_data);
ctf_params_ = textread(fname_ctf_params);
CTF_Defocus_U_ = ctf_params_(:,1);
CTF_Defocus_V_ = ctf_params_(:,2);
CTF_Defocus_Angle_ = ctf_params_(:,3);
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
CTF_Defocus_Angle = CTF_Defocus_Angle*pi/180.0d0;% convert into radians ;
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
SCTF_ = svd(CTF_k_p_r__(:,1+CTF_idx_(1:n_image_sub)));
n_CTF_rank = min(efind(cumsum(SCTF_,'reverse')/sum(SCTF_)<1e-2));
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_idx_(1:n_image_sub)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
%%%%%%%%;
CTF_uni_k_p_r__ = zeros(n_k_p_r,n_ctf);
for nctf=0:n_ctf-1;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_max*nk_p_r + (0:n_w_max-1);
CTF_uni_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_uni_k_p__(1+tmp_index_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_ctf-1;
CTF_uni_avg_k_p_ = mean(CTF_uni_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_uni_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_uni_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_max*nk_p_r + (0:n_w_max-1);
CTF_uni_avg_k_p_r_(1+nk_p_r) = mean(CTF_uni_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
[UCTF_uni_kc__,SCTF_uni_c__,VCTF_uni_Mc__] = svds(CTF_uni_k_p_r__(:,1+CTF_idx_(1:n_image_sub)),n_CTF_rank);
VSCTF_uni_Mc__ = VCTF_uni_Mc__*SCTF_uni_c__;
%%%%%%%%;
% Now determine the CTF cross correlation. ;
% This depends  on CTF_idx_. ;
%%%%%%%%;
tmp_CTF_avg_k_p_ = mean(CTF_k_p__(:,1+CTF_idx_(1+(0:n_image_sub-1))),2);
tmp_CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_(1+tmp_index_));
end;%for nk_p_r=0:n_k_p_r-1;
CTF_k_p_r_xavg__ = tmp_CTF_avg_k_p_r_ * transpose(tmp_CTF_avg_k_p_r_);
CTF_k_p_r_xcor__ = CTF_k_p_r__(:,1+CTF_idx_(1+(0:n_image_sub-1))) * transpose(CTF_k_p_r__(:,1+CTF_idx_(1+(0:n_image_sub-1)))) / n_image_sub;
%%%%%%%%;
save(fname_mat ...
     ,'n_ctf' ...
     ,'CTF_idx_' ...
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
     ,'CTF_uni_k_p__' ...
     ,'CTF_uni_k_p_r__' ...
     ,'CTF_uni_avg_k_p_' ...
     ,'CTF_uni_avg_k_p_r_' ...
     ,'UCTF_uni_kc__','VSCTF_uni_Mc__' ...
     ,'CTF_k_p_r_xavg__' ...
     ,'CTF_k_p_r_xcor__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
subplot(2,3,1+nplot);
imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_k_p__(:,1+nplot)),[-1,+1],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title(sprintf('CTF(k) nctf %d',nplot));
end;%for nplot=0:n_plot-1;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_k_p_r_xxxx__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
fname_euler_angle = sprintf('%s/euler_angles',dir_data);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
fclose(fp);
euler_angle_marina_ = zeros(3,n_image); %<-- [polar_a,azimu_b,gamma_z]. ;
for nimage=0:n_image-1;
euler_angle_marina_(:,1+nimage) = convert_euler_relion_to_marina([euler_angle_load_{1}(1+nimage),euler_angle_load_{2}(1+nimage),euler_angle_load_{3}(1+nimage)]);
end;%for nimage=0:n_image-1;
delta_read_x_ = euler_angle_load_{4}*(2/n_x_u);
delta_read_y_ = euler_angle_load_{5}*(2/n_x_u);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_euler_angle_marina_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
subplot(2,5,1+0);hist(euler_angle_marina_(1,:),linspace(0,pi,128)); axis tight; title('polar_a');
subplot(2,5,2+0);hist(euler_angle_marina_(2,:),linspace(0,2*pi,128)); axis tight; title('azimu_b');
subplot(2,5,3+0);hist(euler_angle_marina_(3,:),linspace(0,2*pi,128)); axis tight; title('gamma_z');
subplot(2,5,4+0);hist(delta_read_x_,linspace(-0.2,+0.2,128)); axis tight; title('delta_x');
subplot(2,5,5+0);hist(delta_read_y_,linspace(-0.2,+0.2,128)); axis tight; title('delta_y');
tmp_index_ = 0:1024-1;
subplot(2,5,1+5);hist(euler_angle_marina_(1,1+tmp_index_),linspace(0,pi,128)); axis tight; title('polar_a sub');
subplot(2,5,2+5);hist(euler_angle_marina_(2,1+tmp_index_),linspace(0,2*pi,128)); axis tight; title('azimu_b sub');
subplot(2,5,3+5);hist(euler_angle_marina_(3,1+tmp_index_),linspace(0,2*pi,128)); axis tight; title('gamma_z sub');
subplot(2,5,4+5);hist(delta_read_x_(1+tmp_index_),linspace(-0.2,+0.2,128)); axis tight; title('delta_x sub');
subplot(2,5,5+5);hist(delta_read_y_(1+tmp_index_),linspace(-0.2,+0.2,128)); axis tight; title('delta_y sub');
figbig;
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
% Now look at some of the experimental images associated with rib80s. ;
% Note that these will *NOT* be corrected/centered according to the (presumed) displacements. ;
% Moreover, to save space, we will not store the M_x_c___. ;
%%%%%%%%;
fname_M_x_c_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_M_x_c___.mat',dir_trunk);
if (flag_recalc | ~exist(fname_M_x_c_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_M_x_c_mat));
fname_euler_angle = sprintf('%s/euler_angles',dir_data);
fp = fopen(fname_euler_angle,'r');
euler_angle_load_ = textscan(fp,'%f%f%f\n%f%f\n');
fclose(fp);
euler_angle_marina_ = zeros(3,n_image); %<-- [polar_a,azimu_b,gamma_z]. ;
for nimage=0:n_image-1;
euler_angle_marina_(:,1+nimage) = convert_euler_relion_to_marina([euler_angle_load_{1}(1+nimage),euler_angle_load_{2}(1+nimage),euler_angle_load_{3}(1+nimage)]);
end;%for nimage=0:n_image-1;
delta_read_x_ = euler_angle_load_{4}*(2/n_x_u);
delta_read_y_ = euler_angle_load_{5}*(2/n_x_u);
% Note that delta_read_r_ is not too large. ;
% prctile(sqrt(delta_read_x_(:).^2 + delta_read_y_(:).^2),95) ; %<-- 0.1182. ;
grid_x_u_ = linspace(-1,+1,n_x_u+1); grid_x_u_ = x_p_r_max*grid_x_u_(1:end-1);
fname_image_mda = sprintf('%s/images_mda',dir_data);
M_x_c___ = MDA_read_r8(fname_image_mda);
n_image_sub = min(n_image_sub,size(M_x_c___,3));
save(fname_M_x_c_mat ...
     ,'euler_angle_load_' ...
     ,'euler_angle_marina_' ...
     ,'delta_read_x_' ...
     ,'delta_read_y_' ...
     ,'grid_x_u_' ...
     ,'n_image_sub' ...
     );
end;%if (~exist(fname_M_x_c_mat,'file'));
if ( exist(fname_M_x_c_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_M_x_c_mat));
load(fname_M_x_c_mat);
end;%if ( exist(fname_M_x_c_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_euler_angle_marina_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
subplot(2,2,3); bar(linspace(0,1*pi,32),hist(euler_angle_marina_(1,:),linspace(0,1*pi,32))); xlabel('polar_a','Interpreter','none'); ylabel('#');
subplot(2,2,2); bar(linspace(0,2*pi,64),hist(euler_angle_marina_(2,:),linspace(0,2*pi,64))); xlabel('azimu_b','Interpreter','none'); ylabel('#');
subplot(2,2,1); bar(linspace(0,2*pi,64),hist(euler_angle_marina_(3,:),linspace(0,2*pi,64))); xlabel('gamma_z','Interpreter','none'); ylabel('#');
subplot(2,2,4); imagesc(log10(1+hist2d_0(euler_angle_marina_(2,:),euler_angle_marina_(1,:),64,32,[0,2*pi],[0,1*pi])),[0,3]); colormap(colormap_beach()); xlabel('azimu_b','Interpreter','none'); ylabel('polar_a','Interpreter','none'); axisnotick;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_delta_read_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
imagesc(log10(1+hist2d_0(delta_read_x_,delta_read_y_,n_x_u,n_x_u,[-1,+1],[-1,+1])),[0,3]); 
axis image; xlabel('x'); ylabel('y'); title('delta read');
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
% Compare images with shifted images. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_M_x_c___center_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
N_x_c___ = M_x_c___;
O_x_c___ = M_x_c___;
P_x_c___ = M_x_c___;
for nimage_sub=0:n_image_sub-1;
M_x_c_ = squeeze(M_x_c___(:,:,1+nimage_sub));
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
N_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,-1*delta_read_x_(1+nimage_sub),-1*delta_read_y_(1+nimage_sub));
N_x_c___(:,:,1+nimage_sub) = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,N_k_p_.*weight_2d_k_all_);
O_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
O_x_c___(:,:,1+nimage_sub) = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,O_k_p_.*weight_2d_k_all_);
P_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,0*delta_read_x_(1+nimage_sub),0*delta_read_y_(1+nimage_sub));
P_x_c___(:,:,1+nimage_sub) = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,P_k_p_.*weight_2d_k_all_);
end;%for nimage_sub=0:n_image_sub-1;
M_x_c_avg__ = real(mean(M_x_c___,3));
N_x_c_avg__ = real(mean(N_x_c___,3));
O_x_c_avg__ = real(mean(O_x_c___,3));
P_x_c_avg__ = real(mean(P_x_c___,3));
subplot(2,2,1); imagesc(M_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('M (orig)');
subplot(2,2,2); imagesc(N_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('N (-delta)');
subplot(2,2,3); imagesc(O_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('O (+delta)');
subplot(2,2,4); imagesc(P_x_c_avg__); axis image; axisnotick; colormap(colormap_beach()); title('P (none)');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
clear N_x_c___ O_x_c___ P_x_c___ ;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_M_x_c___sample',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
for nimage_sub=1;%for nimage_sub=1:n_image_sub;
M_x_c_ = squeeze(M_x_c___(:,:,nimage_sub));
M_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
N_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1));
%figure(1+nimage_sub); 
figure(1); 
colormap(colormap_beach()); 
subplot(2,2,1); imagesc_c(n_x_u,grid_x_u_,n_x_u,grid_x_u_,real(M_x_c_),[],colormap_beach()); 
set(gca,'XTick',[],'YTick',[]); axis image; title('M(x)');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(M(k))');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(M(k))');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,real(N_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(N(k))');
drawnow();
end;%for nimage_sub=1:n_image_sub;
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
% Now generate M_k_p__ and M_uni_k_p__. ;
% if str_CTF_factor contains 'rescale', then we rescale each image by CTF_avg_k_p_r_ ./ CTF_k_p_r__(:,1+CTF_idx_(1+nM)). ;
% Note that this will also rescale the noise for each image. ;
% if str_CTF_factor contains 'original', then we do not rescale (and use the original CTF associated with each image). ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_M_k_p___.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
O_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds not centered images (i.e., original). ;
M_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds yes centered images (except, actually not centered). ;
M_uni_k_p__ = zeros(n_w_max*n_k_p_r,n_image_sub); %<-- holds yes centered images (except, actually not centered) on uniform grid. ;
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,100)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
tmp_M_x_c_ = squeeze(M_x_c___(:,:,1+nimage_sub));
%%%%%%%%;
% First convert tmp_M_x_c_ on the (nonuniform) polar grid used for the templates. ;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
% Now do *NOT* translate according to delta_read__. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
% Now *DO* translate according to delta_read__. ;
%tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
O_k_p__(:,1+nimage_sub) = tmp_O_k_p_; %<-- not centered (i.e., original). ;
M_k_p__(:,1+nimage_sub) = tmp_M_k_p_; %<-- yes centered (except, actually not centered). ;
%%%%%%%%;
% Second convert tmp_M_x_c_ on a uniform polar grid to create principled-images. ;
%%%%%%%%;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
% Now do *NOT* translate according to delta_read__. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
% Now *DO* translate according to delta_read__. ;
%tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
tmp_M_k_p__ = reshape(tmp_M_k_p_,n_w_max,n_k_p_r);
M_uni_k_p__(:,1+nimage_sub) = tmp_M_k_p__(:);
end;%for nimage_sub=0:n_image_sub-1;
save(fname_mat ...
     ,'M_k_p__' ...
     ,'M_uni_k_p__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

%%%%%%%%;
delta_sigma = 1.0 * std([delta_read_x_;delta_read_y_]); %<-- no reduction. ;
str_cost = sprintf('2d_xcor_d%.4d',floor(1000*delta_sigma));
%str_cost = sprintf('2d_emp_N%d',n_image_sub);
%%%%%%%%;
if ( strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s',str_cost);
end;%if ( strcmp(str_CTF_factor,'original'));
if (~strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s_%s',str_cost,str_CTF_factor);
end;%if (~strcmp(str_CTF_factor,'original'));

flag_check=0;
if flag_check;
%%%%%%%%;
% Now check ctf aggregated over all images. ;
%%%%%%%%;
n_M = n_image_sub;
M_k_p_l2__ = zeros(n_k_p_r,n_M);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
M_k_p_l2__(1+nk_p_r,:) = sqrt(mean(abs(M_k_p__(1+tmp_index_,:)).^2,1));
end;%for nk_p_r=0:n_k_p_r-1;
n_S = n_viewing_all;
S_k_p_l2__ = zeros(n_k_p_r,n_S);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
S_k_p_l2__(1+nk_p_r,:) = sqrt(mean(abs(S_k_p__(1+tmp_index_,:)).^2,1));
end;%for nk_p_r=0:n_k_p_r-1;
figure(1);clf;
colormap(colormap_beach());
subplot(2,2,1); imagesc(M_k_p_l2__); axisnotick; title('M');
subplot(2,2,2); imagesc(S_k_p_l2__); axisnotick; title('S');
subplot(2,2,3); plot(k_p_r_,mean(M_k_p_l2__,2)./mean(S_k_p_l2__,2));
subplot(2,2,4); plot(k_p_r_,CTF_avg_k_p_r_);
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now check ctf per image. ;
% Seems like the ctf is actually a reasonable fit. ;
%%%%%%%%;
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
%%%%%%%%;
n_M = n_image_sub;
for nM=0:n_M-1;
tmp_euler_polar_a = +euler_angle_marina_(1,1+nM);
tmp_euler_azimu_b = +euler_angle_marina_(2,1+nM);
tmp_euler_gamma_z = -euler_angle_marina_(3,1+nM);
tmp_image_delta_x = +1.0*delta_read_x_(1+nM);
tmp_image_delta_y = +1.0*delta_read_y_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p__(:,1+nS);
%%%%%%%%;
M_k_p_l2_ = zeros(n_k_p_r,1);
S_k_p_l2_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
M_k_p_l2_(1+nk_p_r) = sqrt(mean(abs(M_k_p_(1+tmp_index_)).^2));
S_k_p_l2_(1+nk_p_r) = sqrt(mean(abs(S_k_p_(1+tmp_index_)).^2));
end;%for nk_p_r=0:n_k_p_r-1;
hold on; plot(CTF_k_p_r__(:,1+CTF_idx_(1+nM)),real(M_k_p_l2_./S_k_p_l2_),'-'); hold off;
%subplot(1,2,1); hold on; plot(k_p_r_,CTF_k_p_r__(:,1+CTF_idx_(1+nM)),'k-'); hold off;
%subplot(1,2,2); hold on; plot(k_p_r_,real(M_k_p_l2_./S_k_p_l2_),'r-'); hold off;
%%%%%%%%;
end;%for nM=0:n_M-1;
%%%%%%%%;
end;%if flag_check;

verbose=1;
%%%%%%%%;
% Now test out principled marching. ;
% First set up cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_X_%s.mat',dir_trunk,str_cost);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
if strfind(str_cost,'3d'); 
tmp_weight_k_p_r_ = weight_3d_k_p_r_;
end;%if strfind(str_cost,'3d'); 
if strfind(str_cost,'2d'); 
tmp_weight_k_p_r_ = weight_2d_k_p_r_;
end;%if strfind(str_cost,'2d'); 
tmp_CTF__ = CTF_k_p_r_xcor__; %<-- default. ;
if strfind(str_cost,'xcor');
tmp_CTF__ = CTF_k_p_r_xcor__; %<-- use CTF_k_p_r_xcor__ here to account for average cross-correlation of CTF functions. ;
end;%if strfind(str_cost,'xcor');
if strfind(str_cost,'xavg');
tmp_CTF__ = CTF_k_p_r_xavg__; %<-- use CTF_k_p_r_xavg__ here to account for mean CTF function. ;
end;%if strfind(str_cost,'xavg');
[X_d0__,X_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,tmp_weight_k_p_r_,l_max_max,a_k_Y_quad__,tmp_CTF__); 
if (strfind(str_cost,'xcor') | strfind(str_cost,'xavg'));
if (delta_sigma>0); [X__,X_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,tmp_weight_k_p_r_,l_max_,1,1,a_k_Y_quad_,tmp_CTF__,delta_sigma); end;
if (delta_sigma==0); [X__,X_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,tmp_weight_k_p_r_,l_max_max,a_k_Y_quad__,tmp_CTF__); end;
end;%if (strfind(str_cost,'xcor') | strfind(str_cost,'xavg'));
if (strfind(str_cost,'emp'));
[X__,X_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_S,S_k_p__);
end;%if (strfind(str_cost,'emp'));
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX__,SX__,VX__] = svds(X__,n_UX_rank);
n_UX_rank = numel(find(diag(SX__)/SX__(1,1)>1e-3));
[UX__,SX__,VX__] = svds(X__,n_UX_rank);
[UX_d0__,SX_d0__,VX_d0__] = svds(X_d0__,n_UX_rank);
a_CTF_UX_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) = a_CTF_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','str_cost','str_combine'...
     ,'X_d0__','X_d0_weight_r_' ...
     ,'X__','X_weight_r_' ...
     ,'n_UX_rank'...
     ,'UX__','SX__','VX__','a_CTF_UX_Y_quad__' ...
     ,'UX_d0__','SX_d0__','VX_d0__'...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
% Check cost for each principal-mode. ;
%%%%%%%%;
flag_check=0;
if flag_check;
for nUX_rank=0:min(36,n_UX_rank-1);%for nUX_rank=0:n_UX_rank-1;
[tmp_X,tmp_X_ori,tmp_X_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_CTF_UX_Y_quad__(:,1+nUX_rank),a_CTF_UX_Y_quad__(:,1+nUX_rank));
tmp_Z = transpose(UX__(:,1+nUX_rank))*X__*(UX__(:,1+nUX_rank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %+0.16f tmp_X %+0.16f tmp_X_ori*tmp_weight_so3 %+0.16f tmp_X_tau %+0.16f ratio %+0.16f',nUX_rank,n_UX_rank,tmp_Z,tmp_X,tmp_X_ori*tmp_weight_so3,tmp_X_tau,(tmp_X_ori*tmp_weight_so3)/tmp_X_tau));
end;%for nUX_rank=0:n_UX_rank-1;
end;%if flag_check;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_X_%s_A',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_X_%s_B',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
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
quad_lim_ = 0.5 * abs(a_CTF_UX_Y_quad__(1,1)) * [-1,+1];
for nplot=0:n_plot-1;
nk_p_r = plot_nk_p_r_(1+nplot);
[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_CTF_UX_Y_quad__(:,1+nk_p_r)),k_u_res,2*k_u_res);
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
% Now generate principled-[templates*CTF_avg]. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_CTF_UX_%s_S_k_p___.mat',dir_trunk,str_cost);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
CTF_UX_S_k_p_wSn___ = zeros(n_w_max,n_S,n_UX_rank);
CTF_UX_S_k_q_wSn___ = zeros(n_w_max,n_S,n_UX_rank);
for nS=0:n_S-1;
if (mod(nS,100)==0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
tmp_CTF_S_k_p_ = S_k_p__(:,1+nS).*CTF_avg_k_p_; %<-- use average templates here under assumption that templates are used alone. ;
tmp_CTF_S_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_CTF_S_k_p_);
tmp_CTF_S_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_CTF_S_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
tmp_CTF_S_k_p__ = reshape(tmp_CTF_S_k_p_,n_w_max,n_k_p_r);
for nUX_rank=0:n_UX_rank-1;
tmp_CTF_UX_S_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_CTF_UX_S_k_p_ = tmp_CTF_UX_S_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_CTF_S_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
CTF_UX_S_k_p_wSn___(:,1+nS,1+nUX_rank) = tmp_CTF_UX_S_k_p_;
CTF_UX_S_k_q_wSn___(:,1+nS,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_CTF_UX_S_k_p_);
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nS=0:n_S-1;
save(fname_mat ...
     ,'CTF_UX_S_k_p_wSn___' ...
     ,'CTF_UX_S_k_q_wSn___' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_UX_%s_S_k_p___',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(real(squeeze(CTF_UX_S_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_UX_%s_S_k_q___',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(log10(abs(squeeze(CTF_UX_S_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank))
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_UX_%s_S_k_q___spectrum',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
n_plot=6;
for nplot=0:n_plot-1;
tmp_ = log10(abs(svds(reshape(permute(CTF_UX_S_k_q_wSn___(:,:,1+(0:nplot)),[1,3,2]),[n_w_max*(1+nplot),n_S]),n_w_max)));
subplot(2,3,1+nplot);
plot(tmp_,'o'); xlim([1,n_w_max]); ylim([-7+tmp_(1),tmp_(1)+1]);
xlabel('pc'); ylabel('log10(sigma)');
title(sprintf('spectrum of CTF_UX_S_k_q_wSn___(:,:,nUX_rank==1:%d)',1+nplot),'Interpreter','none');
end;%for nplot=0:n_plot-1;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_CTF_UX_%s_S_k_p___A',dir_trunk,str_cost);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1);
for nUX_rank=0:min(n_UX_rank-1,15-1);
subplot(3,5,1+nUX_rank); hold on;
for nS=0:32:n_S-1;
nc = max(0,min(n_c-1,floor(n_c*nS/n_S)));
plot(2*pi*(0:n_w_max-1)/n_w_max,real(CTF_UX_S_k_p_wSn___(:,1+nS,1+nUX_rank)),'.','Color',c_(1+nc,:)); 
end;%for nS=0:n_S-1;
xlim([0,2*pi]);
ylim(4e-6*[-1,+1]); 
xlabel('gamma_z','Interpreter','none');
ylabel('real(CTF_UX_S_k_p)','Interpreter','none');
title(sprintf('rank %d',nUX_rank));
end; %for nUX_rank=0:15-1;
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
% Now generate principled-images. ;
% if str_CTF_factor contains 'rescale', then we rescale each image by CTF_avg_k_p_r_ ./ CTF_k_p_r__(:,1+CTF_idx_(1+nM)). ;
% Note that this will also rescale the noise for each image. ;
% if str_CTF_factor contains 'original', then we do not rescale (and use the original CTF associated with each image). ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_UX_%s_M_k_p___.mat',dir_trunk,str_combine);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
O_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds not centered images (i.e., original). ;
M_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds yes centered images (except, actually not centered). ;
M_uni_k_p__ = zeros(n_w_max*n_k_p_r,n_image_sub); %<-- holds yes centered images (except, actually not centered) on uniform grid. ;
UX_M_k_p___ = zeros(n_w_max,n_image_sub,n_UX_rank);
UX_M_k_q___ = zeros(n_w_max,n_image_sub,n_UX_rank);
fname_image_mda = sprintf('%s/images_mda',dir_data);
tmp_M_x_c___ = MDA_read_r8(fname_image_mda);
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,100)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
tmp_M_x_c_ = squeeze(tmp_M_x_c___(:,:,1+nimage_sub));
%%%%%%%%;
% First convert tmp_M_x_c_ on the (nonuniform) polar grid used for the templates. ;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
% Now do *NOT* translate according to delta_read__. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
% Now *DO* translate according to delta_read__. ;
%tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
O_k_p__(:,1+nimage_sub) = tmp_O_k_p_; %<-- not centered (i.e., original). ;
M_k_p__(:,1+nimage_sub) = tmp_M_k_p_; %<-- yes centered (except, actually not centered). ;
%%%%%%%%;
% Second convert tmp_M_x_c_ on a uniform polar grid to create principled-images. ;
%%%%%%%%;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
% Now do *NOT* translate according to delta_read__. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
% Now *DO* translate according to delta_read__. ;
%tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+1*delta_read_x_(1+nimage_sub),+1*delta_read_y_(1+nimage_sub));
tmp_M_k_p__ = reshape(tmp_M_k_p_,n_w_max,n_k_p_r);
M_uni_k_p__(:,1+nimage_sub) = tmp_M_k_p__(:);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_M_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+nimage_sub);
if strfind(str_CTF_factor,'rescale'); CTF_factor = CTF_avg_k_p_r_(1+nk_p_r) / max(1e-6,CTF_k_p_r__(1+nk_p_r,1+nctf)); end;
if strfind(str_CTF_factor,'original'); CTF_factor = 1; end;
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*CTF_factor*tmp_M_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
UX_M_k_p___(:,1+nimage_sub,1+nUX_rank) = tmp_UX_M_k_p_;
UX_M_k_q___(:,1+nimage_sub,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_M_k_p_);
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nimage_sub=0:n_image_sub-1;
clear tmp_M_x_c___;
save(fname_mat ...
     ,'O_k_p__' ...
     ,'M_k_p__' ...
     ,'M_uni_k_p__' ...
     ,'UX_M_k_p___' ...
     ,'UX_M_k_q___' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_UX_%s_M_k_p___',dir_trunk,str_combine);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(UX_M_k_p___(:,:,1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(real(squeeze(UX_M_k_p___(:,:,1+nUX_rank))),clim_); xlabel('nimage'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_UX_%s_M_k_q___',dir_trunk,str_combine);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(log10(abs(squeeze(UX_M_k_q___(:,:,1+nUX_rank)))),[-5,0]); xlabel('nimage'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-5,0]',1+nUX_rank))
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_UX_%s_M_k_q___spectrum',dir_trunk,str_combine);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=min(n_UX_rank-1,6);
for nplot=0:n_plot-1;
tmp_ = log10(abs(svds(reshape(permute(UX_M_k_q___(:,:,1+(0:nplot)),[1,3,2]),[n_w_max*(1+nplot),n_image_sub]),n_w_max)));
subplot(2,3,1+nplot);
plot(tmp_,'o'); xlim([1,n_w_max]); ylim([-7+tmp_(1),tmp_(1)+1]);
xlabel('pc'); ylabel('log10(sigma)');
title(sprintf('spectrum of UX_M_k_q___(:,:,nUX_rank==1:%d)',1+nplot),'Interpreter','none');
end;%for nplot=0:n_plot-1;
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_UX_%s_M_k_p___.mat',dir_trunk,str_combine);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

flag_check=0;
if flag_check;
%%%%%%%%;
% check cg_lsq_3. ;
%%%%%%%%;
tmp_verbose=1;
tmp_pm_n_UX_rank = 2;
tmp_pm_n_order = 5;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_CTF_index_ = 1; tmp_pm_CTF_k_p__ = ones(tmp_pm_n_w_sum,1);
tmp_t = toc(tmp_t); if (tmp_verbose); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
tmp_euler_polar_a_ = +euler_angle_marina_(1,1+(0:n_image_sub-1));
tmp_euler_azimu_b_ = +euler_angle_marina_(2,1+(0:n_image_sub-1));
tmp_euler_gamma_z_ = -euler_angle_marina_(3,1+(0:n_image_sub-1));
tmp_n_M = n_image_sub;
tmp_UX_M_k_p_wnM__ = reshape(permute(UX_M_k_p___(:,1:tmp_n_M,1:tmp_pm_n_UX_rank),[1,3,2]),[tmp_pm_n_w_sum,tmp_n_M]);
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
tmp_a_CTF_UX_Y_0lsq_ = ...
cg_lsq_3( ...
 tmp_pm_n_order ...
,tmp_pm_n_k_p_r ...
,tmp_pm_l_max_ ...
,tmp_pm_n_w_ ...
,tmp_n_M ...
,tmp_UX_M_k_p_wnM__ ...
,tmp_pm_CTF_index_ ...
,tmp_pm_CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
);
tmp_t = toc(tmp_t); if (tmp_verbose); disp(sprintf(' %% cg_lsq_3 for a_CTF_UX_Y_0lsq_: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
[ ...
 tmp_CTF_UX_S_k_p__ ...
,~ ...
,~ ...
,~ ...
,tmp_n_viewing_all ...
,tmp_viewing_azimu_b_all_ ...
,tmp_viewing_polar_a_all_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_template_0( ...
 tmp_verbose ...
,tmp_pm_n_k_p_r ...
,tmp_pm_k_p_r_ ...
,tmp_pm_k_p_r_max ...
,tmp_pm_weight_k_p_r_ ...
,tmp_pm_l_max_ ...
,tmp_a_CTF_UX_Y_0lsq_ ...
,tmp_viewing_k_eq_d ...
,-1 ...
,tmp_pm_n_w_ ...
);
tmp_t = toc(tmp_t); if (tmp_verbose); disp(sprintf(' %% get_template_0 for tmp_CTF_UX_S_k_p__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% check cg_lsq_pm_0. ;
%%%%%%%%;
tmp_image_delta_x_ = zeros(tmp_n_M,1);
tmp_image_delta_y_ = zeros(tmp_n_M,1);
tmp_image_I_value_ = ones(tmp_n_M,1);
tmp_n_CTF_rank = 3;
[tmp_UCTF_kc__,tmp_SCTF_c__,tmp_VCTF_Mc__] = svds(CTF_k_p_r__(:,1+CTF_idx_(1:n_image_sub)),tmp_n_CTF_rank); tmp_VSCTF_Mc__ = tmp_VCTF_Mc__*tmp_SCTF_c__;
tmp_a_UCTF_UX_Y_0lsq_ync__ = ...
cg_lsq_pm_0( ...
 tmp_pm_n_order ...
,tmp_pm_n_k_p_r ...
,tmp_pm_k_p_r_ ...
,tmp_pm_l_max_ ...
,tmp_pm_n_w_ ...
,tmp_n_M ...
,tmp_UX_M_k_p_wnM__ ...
,tmp_n_CTF_rank ...
,tmp_VSCTF_Mc__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
,tmp_image_I_value_...
);
%%%%%%%%;
tmp_t = tic();
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
tmp_UX_UCTF_S_k_p_wSc___ = zeros(tmp_pm_n_w_sum,tmp_n_viewing_all,tmp_n_CTF_rank);
for tmp_nCTF_rank=0:tmp_n_CTF_rank-1;
[ ...
 tmp_UX_UCTF_S_k_p_wSc___(:,:,1+tmp_nCTF_rank) ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
get_template_0( ...
 tmp_verbose ...
,tmp_pm_n_k_p_r ...
,tmp_pm_k_p_r_ ...
,tmp_pm_k_p_r_max ...
,tmp_pm_weight_k_p_r_ ...
,tmp_pm_l_max_ ...
,tmp_a_UCTF_UX_Y_0lsq_ync__(:,1+tmp_nCTF_rank) ...
,tmp_viewing_k_eq_d ...
,-1 ...
,tmp_pm_n_w_ ...
);
end;%for tmp_nCTF_rank=0:tmp_n_CTF_rank-1;
tmp_t = toc(tmp_t); if (tmp_verbose); disp(sprintf(' %% get_template_0 for tmp_CTF_UX_S_k_p__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_error_ = zeros(tmp_n_M,1);
for tmp_nM=0:tmp_n_M-1;
tmp_UX_VSCTF_UCTF_S_k_p__ = zeros(tmp_pm_n_w_sum,tmp_n_viewing_all);
for tmp_nCTF_rank=0:tmp_n_CTF_rank-1;
tmp_UX_VSCTF_UCTF_S_k_p__ = tmp_UX_VSCTF_UCTF_S_k_p__ + tmp_VSCTF_Mc__(1+tmp_nM,1+tmp_nCTF_rank)*tmp_UX_UCTF_S_k_p_wSc___(:,:,1+tmp_nCTF_rank);
end;%for tmp_nCTF_rank=0:tmp_n_CTF_rank-1;
tmp_error_(1+tmp_nM) = fnorm(tmp_CTF_UX_S_k_p__ - tmp_UX_VSCTF_UCTF_S_k_p__)/fnorm(tmp_CTF_UX_S_k_p__);
end;%for tmp_nM=0:tmp_n_M-1;
figure(1);clf;
subplot(2,1,1);
plot(1:tmp_n_M,tmp_error_,'k.'); xlabel('nM'); ylabel('S error'); grid on;
xlim([1,tmp_n_M]);
title('tmp_CTF_UX_S_k_p__ versus tmp_UX_VSCTF_UCTF_S_k_p__','Interpreter','none');
subplot(2,1,2);
figbeach; imagesc(CTF_k_p_r__(:,1+CTF_idx_(1:tmp_n_M)),[0.3,1.0]);
axisnotick;
xlabel('nM'); ylabel('k');
%%%%%%%%;
end;%if flag_check;

fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_c_k_Y_.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_n_order = 5; tmp_n_M = n_image_sub;
tmp_euler_polar_a_ = +euler_angle_marina_(1,1+(0:n_image_sub-1));
tmp_euler_azimu_b_ = +euler_angle_marina_(2,1+(0:n_image_sub-1));
tmp_euler_gamma_z_ = -euler_angle_marina_(3,1+(0:n_image_sub-1));
tmp_image_delta_x_ = +1.0*delta_read_x_(1+(0:n_image_sub-1));
tmp_image_delta_y_ = +1.0*delta_read_y_(1+(0:n_image_sub-1));
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
,CTF_idx_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> c_k_Y_reco_ time %0.2fs',tmp_t));
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,c_k_Y_reco_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,c_k_Y_reco_));
X_best_reco = max(tmp_X_best_orig,tmp_X_best_flip);
disp(sprintf(' %% X_best_reco %0.3f',X_best_reco));
%%%%%%%%;
tmp_t = tic;
[c_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,c_k_Y_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_reco_ --> c_k_p_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
eta = pi/k_p_r_max; 
c_x_u_reco_ = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,c_k_p_reco_.*(2*pi)^3.*weight_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: c_x_u_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
save(fname_mat ...
     ,'c_k_Y_reco_' ...
     ,'c_k_p_reco_' ...
     ,'c_x_u_reco_' ...
     ,'X_best_reco' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_c_x_u_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
isosurface_f_x_u_0(reshape(real(c_x_u_reco_),x_u_res,x_u_res,x_u_res),[90,95,99]);
title('c_x_u_','Interpreter','none');
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

flag_check=0;
if flag_check;
%%%%%%%%;
% check advantage gained for specialized svd of UX__(:,~).*X_weight_r_.*besselj(~,~) ;
% Very little advantage at all, especially for large N_pixel. ;
%%%%%%%%;
tmp_eps = 1e-3; 
tmp_n_UX_rank = 15; tmp_N =  0.75;
for tmp_N = [0.5,0.75,1.00,1.50,2.00];
tmp_FTK_0 = gen_UXJsvd_FTK_8(n_k_p_r,k_p_r_,0*ones(n_k_p_r,1) + 1*max(abs(UX__(:,1:tmp_n_UX_rank)),[],2).*X_weight_r_,k_p_r_max,tmp_N,tmp_eps,25,32,33);
tmp_FTK_1 = gen_UXJsvd_FTK_8(n_k_p_r,k_p_r_,1*ones(n_k_p_r,1) + 0*max(abs(UX__(:,1:tmp_n_UX_rank)),[],2).*X_weight_r_,k_p_r_max,tmp_N,tmp_eps,25,32,33);
disp(sprintf(' %% tmp_N %0.2f tmp_n_UX_rank %d: tmp_FTK_0.n_svd_l/tmp_FTK_1.n_svd_l: %0.2f',tmp_N,tmp_n_UX_rank,tmp_FTK_0.n_svd_l/tmp_FTK_1.n_svd_l));
end;%for tmp_N = [0.5,0.75,1.00,1.50,2.00];
end;%flag_check=0;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now calculate the typical directional-derivative (with respect to viewing-angle) of the template-operator applied to the true molecule. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_d0_CTF_UXI_S_k_p__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
pm_n_UX_rank = n_k_p_r;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
n_UXI_rank = pm_n_k_p_r;
UXI__=eye(pm_n_k_p_r);
a_CTF_UXI_Y_quad__ = zeros(pm_n_lm_max,n_UXI_rank);
for nUXI_rank=0:n_UXI_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) = a_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) + UXI__(1+nk_p_r,1+nUXI_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UXI_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUXI_rank=0:n_UXI_rank-1;
%%%%%%%%;
pm_template_k_eq_d = k_eq_d;
pm_viewing_k_eq_d = k_eq_d;
[ ...
 d0_CTF_UXI_S_k_p__ ...
,da_CTF_UXI_S_k_p__ ...
,db_CTF_UXI_S_k_p__ ...
,dc_CTF_UXI_S_k_p__ ...
,dt_CTF_UXI_S_k_p__ ...
,pm_n_w_ ...
,pm_weight_2d_k_p_r_ ...
,pm_weight_2d_k_all_ ...
,pm_n_viewing_all ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,pm_n_viewing_polar_a ...
,pm_viewing_polar_a_ ...
,pm_n_viewing_azimu_b_ ...
,pm_template_k_c_0__ ...
,pm_template_k_c_1__ ...
,pm_template_k_c_2__ ...
] = ...
get_dtemplate_0( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,[] ...
,pm_l_max_ ...
,a_CTF_UXI_Y_quad__(:) ...
,pm_viewing_k_eq_d ...
,0*pm_template_k_eq_d ...
,pm_n_w_ ...
);
%%%%%%%%;
% Now determine scaling factor for signal-to-noise-ratio snr. ;
%%%%%%%%;
% We assume that X_weight_r_ is given by sqrt(weight_2d_k_p_r_);
% This means that sum(X_weight_r_.^2) = pi*k_p_r_max.^2. ;
% Similarly, it means that each image- and template-ring is normalized to have the same variance. ;
% Assuming that each image- and template-ring are subdivided into the same number of arc-elements (i.e., that n_w_ is uniform), ;
% we can estimate the average 'l2-signal' per arc-element simply by taking the l2-norm: ;
% l2_signal = sqrt(mean(abs(d0_CTF_UXI_S_k_p__).^2 * reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]) / (4*pi))). ;
% With this baseline, the variance l2_sigma^2 associated with a given l2_snr is: ;
% l2_sigma = l2_signal/l2_snr. ;
%%%%%%%%;
l2_signal = sqrt(mean(abs(d0_CTF_UXI_S_k_p__).^2 * reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]) / (4*pi))) ;
%%%%%%%%;
% Now we can estimate the critical viewing-angle-error \eta associated with each snr and surface-multiplicity J_{sub}. ;
% To be more explicit: ;
% First we assume that we have viewing-angle errors \eta_{j} for each image M_{j}. ;
% These \eta_{j} will be drawn from a gaussian-distribution in [viewing_azimu_b,viewing_polar_a,viewing_gamma_z] with variance \sigma_{\eta}^{2}/. ;
% We also assume that each image is a CTF*template polluted by errors \eps_{j}, ;
% where the \eps_{j} are iid gaussian with variance \sigma_{\eps}^{2} per entry (i.e., per degree of freedom in image-space). ;
% Now, assuming that the number of images is large, these two sources of noise will induce two errors in the lsq-step. ;
% 1. The first of the lsq-errors, due to the image-errors \eps_{j} will be random. ;
%    If we represent the lsq-problem as: ;
%    P*a_k_Y_0lsq_ = image-data, with P representing the evaluation-operator, ;
%    then the image-data will be the exact template-data + an additional term stacking the various \eps_{j}. ;
%    The isotropic gaussian error \eps will induce a gaussian error in P*a_k_Y_0lsq_ that is also isotropic ;
%    (with the same variance \sigma_{\eps}^{2}), albeit limited to the range of P. ;
%    This will then induce an error in the coefficients of a_k_Y_0lsq_ (which will not exactly equal those of a_k_Y_diff_). ;
%    Assuming that there are multiple images, and that the viewing-angles are uniformly distributed, ;
%    The difference a_k_Y_0lsq_ - a_k_Y_diff_ =: \chi will be drawn from an isotropic gaussian with variance ;
%    \sigma_{\chi}^{2} = \sigma_{eta}^{2}/J_{sub}, ;
%    where J_{sub} represents the multiplicity of each surface point (on the sphere) as represented across all the images. ;
% 2. The second of these lsq-errors, due to the viewing-angle-errors \eta_{j} will be systematic. ;
%    Assuming the viewing-angles are distributed uniformly, the \eta_{j} will induce an error in ;
%    calculating a_k_Y_0lsq_ that is similar to the effect of surface-diffusion: ;
%    a_k_Y_diff___(nk_p_r,l_val,m_val) \approx exp(-l_val*(l_val+1)*\sigma_{eta}^{2}) a_k_Y_true___(nk_p_r,l_val,m_val), ;
%    implying that: ;
%    a_k_Y_diff___(nk_p_r,l_val,m_val) - a_k_Y_true___(nk_p_r,l_val,m_val) \approx -l_val*(l_val+1)*\sigma_{eta}^{2}. ;
%%%%%%%%;
% The subsequent alignment-step will then produce several errors to each image viewing-angle. ;
% Assuming that both \sigma_{\eps} and \sigma_{\eta} are small: ;
% 0. The isotropic eps_{j} for each image will induce an alignment-error E0_{j} corresponding to the solution ;
%    to the least-squares problem: \grad_{\tau}(S\tau_{j}^{true} a_k_Y_0lsq_) * E0_{j} = eps_{j}. ;
%    The resulting E0_{j} will be drawn from an anisotropic gaussian with variances given by sigma_{\eps}^{2} multiplied by ;
%    the inverse-squares of the singular values of \grad_{\tau}(S\tau_{j}^{true} a_k_Y_0lsq_). ;
%    Note that, to lowest order, \grad_{\tau}(S\tau_{j}^{true} a_k_Y_0lsq_) \approx \grad_{\tau}(S\tau_{j}^{true} a_k_Y_true_)
% 1. The random lsq-error will correspond to a random alignment-error E1_{j}. ;
%    Much like E0_{j}, E1_{j} will (approximately) solve the least-squares problem: ;
%    \grad_{\tau}(S\tau_{j}^{true} a_k_Y_true_) * E0_{j} = \chi_{j}. ;
%    Because the variance of \chi is 1/J_{sub} times the variance of \eps, ;
%    the magnitude of E1_{j} will be correspondingly lower than that of E0_{j}. ;
% 2. The systematic lsq-error will correspond to a systematic alignment-error E2_{j} solving: ;
%    \grad_{\tau}(S\tau_{j}^{true} a_k_Y_0lsq_) * E2_{j} = S\tau_{j}^{true} (a_k_Y_diff_ - a_k_Y_true_). ;
%    Note that the systematic image/template-error S\tau_{j}^{true} (a_k_Y_diff_ - a_k_Y_true_)
%    is proportional to \sigma_{\eta}^{2}, as will be E2_{j}. ;
%    Each E2_{j} can then be written as: ;
%    E2_{j} = Pinv(\grad_{\tau}^{\inv}(S\tau_{j}^{true} a_k_Y_0lsq_)) * \partial_{t} (S\tau_{j}^{true} a_k_Y_true_) * t, ;
%    where t = \sigma_{\eta}^{2}. ;
%    Below we will denote this expression as: ;
%    |E2_{j}| = E2_bar_{j} * t = E2_bar_{j} * \sigma_{\eta}^{2} ;
% Note that we typically expect each of these errors to be uncorrelated (i.e., at a random orientation to one another). ;
% Consequently, the total variance of \eps_{j} + chi_{j} is \sigma_{\eps}^{2}*(1+1/J_{sub}) \approx \sigma_{\eps}^{2}. ;
% Using similar logic, we expect the typical magnitude of the alignment-error to be: ;
% |\eta_{j}^{new}| \approx E2_bar_{j} * \sigma_{\eta}^{2} + |\grad_{\tau}^{\inv}(S\tau_{j}^{true} a_k_Y_true_)| * \sigma_{\eps}*sqrt(1+1/J_{sub}). ;
%%%%%%%%;
% If the typical |\eta_{j}^{new}| is less than \sigma_{\eta}, then we expect the alignment-errors to shrink as alternating-minimization proceeds. ;
% On the other hand, if the typical |\eta_{j}^{new}| is larger than \sigma_{\eta}, then we expect that alignment-errors will grow. ;
% Notably, we expect to see the alignment-errors for a given viewing-angle converge to a size given roughly by \eta^{equi}_{j}: ;
% \eta^{equi}_{j} = E2_bar_{j} * [\eta^{equi}_{j}]^2 + |\grad_{\tau}^{\inv}(S\tau_{j}^{true} a_k_Y_true_)| * \sigma_{\eps}*sqrt(1+1/J_{sub}). ;
% We solve for this critical \eta^{\equi}_{j} later on. ;
save(fname_mat ...
     ,'l2_signal' ...
     ,'pm_template_k_eq_d' ...
     ,'pm_viewing_k_eq_d' ...
     ,'d0_CTF_UXI_S_k_p__' ...
     ,'da_CTF_UXI_S_k_p__' ...
     ,'db_CTF_UXI_S_k_p__' ...
     ,'dc_CTF_UXI_S_k_p__' ...
     ,'dt_CTF_UXI_S_k_p__' ...
     ,'pm_n_w_' ...
     ,'pm_weight_2d_k_p_r_' ...
     ,'pm_weight_2d_k_all_' ...
     ,'pm_n_viewing_all' ...
     ,'pm_viewing_azimu_b_all_' ...
     ,'pm_viewing_polar_a_all_' ...
     ,'pm_viewing_weight_all_' ...
     ,'pm_n_viewing_polar_a' ...
     ,'pm_viewing_polar_a_' ...
     ,'pm_n_viewing_azimu_b_' ...
     ,'pm_template_k_c_0__' ...
     ,'pm_template_k_c_1__' ...
     ,'pm_template_k_c_2__' ...
     ,'pm_n_w_max' ...
     ,'pm_n_w_sum' ...
     ,'pm_n_w_csum_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_d0_CTF_UXI_S_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
nviewing_all = max(0,min(pm_n_viewing_all-1,round(pm_n_viewing_all*nplot/n_plot)));
subplot(6,5,1+5*nplot+0);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(d0_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('d0_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+1);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(da_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('da_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+2);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(db_CTF_UXI_S_k_p__(:,1+nviewing_all)/sin(pm_viewing_polar_a_all_(1+nviewing_all))),[],colormap_beach());
axis image; axisnotick; title(sprintf('db_CTF_UXI_S_k_p__(:,1+nv) / sin(polar_a) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+3);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(dc_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('dc_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+4);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(dt_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('dt_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
end;%for nplot=0:n_plot-1;
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
% Now calculate the typical variance induced by the gradient. ;
%%%%%%%%;
process_UXI_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_CTF_UXI_S_k_p__ ...
,da_CTF_UXI_S_k_p__ ...
,db_CTF_UXI_S_k_p__ ...
,dc_CTF_UXI_S_k_p__ ...
,dt_CTF_UXI_S_k_p__ ...
);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_d0_CTF_UXI_S_grad_FIGA',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
l2_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_l2_snr=numel(l2_snr_);
eta_equi_upb___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
eta_equi_lob___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
l2_sigma = 0; if l2_snr>0; l2_sigma = l2_signal/l2_snr; end;
for nk_p_r=0:pm_n_k_p_r-1;
for nviewing_all=0:pm_n_viewing_all-1;
tmp_a = process_UXI_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_all);
tmp_b = -1.0d0;
tmp_c = process_UXI_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_all)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb = +Inf; eta_equi_lob = +Inf; end;
if (tmp_d>=0); eta_equi_upb = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
eta_equi_upb___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_upb;
eta_equi_lob___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_lob;
end;%for nviewing_all=0:pm_n_viewing_all-1;
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
eta_equi_upb__ = reshape(reshape(eta_equi_upb___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
eta_equi_lob__ = reshape(reshape(eta_equi_lob___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
%%%%;
c__ = colormap_beach();
c__(end,:) = [0,0,0];
pm_nk_p_r_ = round(linspace(0,pm_n_k_p_r-1,9));
np=0;
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
for nl=0:numel(pm_nk_p_r_)-1;%for nk_p_r=0:pm_n_k_p_r-1;
nk_p_r = pm_nk_p_r_(1+nl);
subplot(n_l2_snr,9,1+np);
imagesc_polar_a_azimu_b_0(pm_viewing_polar_a_all_,pm_viewing_azimu_b_all_,log10(abs(eta_equi_lob___(1+nl2_snr,1+nk_p_r,:))),[-2,0],c__);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('snr %0.2f R %d',l2_snr,nk_p_r));
np=np+1;
end;%for nl=0:numel(pm_nk_p_r_)-1;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
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
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now calculate the typical directional-derivative (with respect to viewing-angle) of the template-operator applied to the true molecule. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80s_0_d0_CTF_UX_S_k_p__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
pm_n_UX_rank = 9;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_); pm_n_lm_max = max(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
pm_template_k_eq_d = k_eq_d;
pm_viewing_k_eq_d = k_eq_d;
[ ...
 d0_CTF_UX_S_k_p__ ...
,da_CTF_UX_S_k_p__ ...
,db_CTF_UX_S_k_p__ ...
,dc_CTF_UX_S_k_p__ ...
,dt_CTF_UX_S_k_p__ ...
,pm_n_w_ ...
,pm_weight_2d_k_p_r_ ...
,pm_weight_2d_k_all_ ...
,pm_n_viewing_all ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,pm_n_viewing_polar_a ...
,pm_viewing_polar_a_ ...
,pm_n_viewing_azimu_b_ ...
,pm_template_k_c_0__ ...
,pm_template_k_c_1__ ...
,pm_template_k_c_2__ ...
] = ...
get_dtemplate_0( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,[] ...
,pm_l_max_ ...
,a_CTF_UX_Y_quad__(:) ...
,pm_viewing_k_eq_d ...
,0*pm_template_k_eq_d ...
,pm_n_w_ ...
);
save(fname_mat ...
     ,'pm_template_k_eq_d' ...
     ,'pm_viewing_k_eq_d' ...
     ,'d0_CTF_UX_S_k_p__' ...
     ,'da_CTF_UX_S_k_p__' ...
     ,'db_CTF_UX_S_k_p__' ...
     ,'dc_CTF_UX_S_k_p__' ...
     ,'dt_CTF_UX_S_k_p__' ...
     ,'pm_n_w_' ...
     ,'pm_weight_2d_k_p_r_' ...
     ,'pm_weight_2d_k_all_' ...
     ,'pm_n_viewing_all' ...
     ,'pm_viewing_azimu_b_all_' ...
     ,'pm_viewing_polar_a_all_' ...
     ,'pm_viewing_weight_all_' ...
     ,'pm_n_viewing_polar_a' ...
     ,'pm_viewing_polar_a_' ...
     ,'pm_n_viewing_azimu_b_' ...
     ,'pm_template_k_c_0__' ...
     ,'pm_template_k_c_1__' ...
     ,'pm_template_k_c_2__' ...
     ,'pm_n_w_max' ...
     ,'pm_n_w_sum' ...
     ,'pm_n_w_csum_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_d0_CTF_UX_S_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
nviewing_all = max(0,min(pm_n_viewing_all-1,round(pm_n_viewing_all*nplot/n_plot)));
subplot(6,5,1+5*nplot+0);
imagesc(real(reshape(d0_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('d0_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+1);
imagesc(real(reshape(da_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('da_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+2);
imagesc(real(reshape(db_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('db_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+3);
imagesc(real(reshape(dc_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('dc_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+4);
imagesc(real(reshape(dt_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('dt_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
end;%for nplot=0:n_plot-1;
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
% Now calculate the typical variance induced by the gradient. ;
%%%%%%%%;
process_UX_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_CTF_UX_S_k_p__ ...
,da_CTF_UX_S_k_p__ ...
,db_CTF_UX_S_k_p__ ...
,dc_CTF_UX_S_k_p__ ...
,dt_CTF_UX_S_k_p__ ...
);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_d0_CTF_UX_S_grad_FIGA',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
l2_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_l2_snr=numel(l2_snr_);
eta_equi_upb___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
eta_equi_lob___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
l2_sigma = 0; if l2_snr>0; l2_sigma = l2_signal/l2_snr; end;
for nk_p_r=0:pm_n_k_p_r-1;
for nviewing_all=0:pm_n_viewing_all-1;
tmp_a = process_UX_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_all);
tmp_b = -1.0d0;
tmp_c = process_UX_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_all)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb = +Inf; eta_equi_lob = +Inf; end;
if (tmp_d>=0); eta_equi_upb = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
eta_equi_upb___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_upb;
eta_equi_lob___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_lob;
end;%for nviewing_all=0:pm_n_viewing_all-1;
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
eta_equi_upb__ = reshape(reshape(eta_equi_upb___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
eta_equi_lob__ = reshape(reshape(eta_equi_lob___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
%%%%;
c__ = colormap_beach();
c__(end,:) = [0,0,0];
np=0;
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
for nk_p_r=0:pm_n_k_p_r-1;
subplot(n_l2_snr,pm_n_k_p_r,1+np);
imagesc_polar_a_azimu_b_0(pm_viewing_polar_a_all_,pm_viewing_azimu_b_all_,log10(abs(eta_equi_lob___(1+nl2_snr,1+nk_p_r,:))),[-2,0],c__);
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick;
title(sprintf('snr %0.2f R %d',l2_snr,nk_p_r));
np=np+1;
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
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
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% check to see if translations can be estimated using only a few principal-modes. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-2;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 1;
dat_n_iteration = 4;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_UX_ = sum(UX__(:,1:dat_n_UX_rank).^2,2);
tmp_UX_csum_ = cumsum(tmp_UX_,'reverse')/sum(tmp_UX_);
index_k_p_r_use = min(efind(tmp_UX_csum_<1e-2));
tmp_n_k_p_r = 1+index_k_p_r_use;
n_w_uni_ = n_w_max*ones(tmp_n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(n_w_uni_),:);
tmp_l_max_ = l_max_(1:tmp_n_k_p_r);
tmp_n_lm_ = (tmp_l_max_+1).^2;
tmp_n_lm_max = max(tmp_n_lm_);
tmp_n_lm_sum = sum(tmp_n_lm_);
tmp_n_lm_csum_ = cumsum([0;tmp_n_lm_]);
tmp_a_CTF_UX_Y_quad__ = zeros(tmp_n_lm_max,dat_n_UX_rank);
for dat_nUX_rank=0:dat_n_UX_rank-1;
for tmp_nk_p_r=0:tmp_n_k_p_r-1;
tmp_l_max = tmp_l_max_(1+tmp_nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = tmp_n_lm_csum_(1+tmp_nk_p_r) + (0:tmp_n_lm-1);
tmp_a_CTF_UX_Y_quad__(1:tmp_n_lm,1+dat_nUX_rank) = ...
tmp_a_CTF_UX_Y_quad__(1:tmp_n_lm,1+dat_nUX_rank) ...
+ UX__(1+tmp_nk_p_r,1+dat_nUX_rank)*X_weight_r_(1+tmp_nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+tmp_nk_p_r);
%<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for tmp_nk_p_r=0:tmp_n_k_p_r-1;
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
dat_f_rand = 0.05; dat_flag_plot=0;
tpmutam_experimental_0( ...
 dat_infix ...
,dat_rseed ...
,dat_n_M ...
,dat_M_k_p__ ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,dir_trunk ...
,tmp_n_k_p_r ...
,k_p_r_(1:tmp_n_k_p_r) ...
,k_p_r_(tmp_n_k_p_r) ...
,n_w_uni_ ...
,weight_3d_k_p_r_(1:tmp_n_k_p_r) ...
,weight_2d_k_p_r_(1:tmp_n_k_p_r) ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__(1:tmp_n_k_p_r,:) ...
,l_max_(1:tmp_n_k_p_r) ...
,UX__(1:tmp_n_k_p_r,1:dat_n_UX_rank) ...
,X_weight_r_(1:tmp_n_k_p_r) ...
,tmp_a_CTF_UX_Y_quad__ ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,[] ...
,dat_f_rand ...
,dat_flag_plot ...
);
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%;
% Now check to see if the translations discovered correlate with the true translations. ;
%%%%%%%%;
c_ = colormap_beach(); n_c = size(c_,1);
%hsv_ = colormap('hsv'); n_hsv = size(hsv_,1);
%nhsv_ = max(0,min(n_hsv-1,floor(n_hsv*(pi+atan2(delta_read_y_(1:n_image_sub),delta_read_x_(1:n_image_sub)))/(2*pi))));
%scatter(delta_read_x_(1:n_image_sub),delta_read_y_(1:n_image_sub),15,hsv_(1+nhsv_,:),'filled'); axis([-1,+1,-1,+1]); axis square;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_delta_nd%d_FIGA',dir_trunk,ndelta_r_max_factor);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%;
figure(1+ndelta_r_max_factor); clf;
%%%%%%%%;
subplot(2,2,2);
c2d_ = colormap_gaussian_2d(delta_read_x_,delta_read_y_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_x_(1:n_image_sub),delta_read_y_(1:n_image_sub),64,c2d_(1:n_image_sub,:),'filled');
%scatter(delta_read_x_(1:end),delta_read_y_(1:end),25,c2d_(1:end,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
figbig;
%%%%%%%%;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
for dat_nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*(ndat_rseed + dat_nUX_rank*n_dat_rseed)/(n_dat_rseed*dat_n_UX_rank))));
dat_rseed = dat_rseed_(1+ndat_rseed);
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
%%%%%%%%;
if ( exist(MS_fname_mat,'file'));
tmp_MS_ = load(MS_fname_mat);
subplot(2,2,1);
hold on;
plot((0:dat_n_iteration-1) + 0*dat_n_iteration,tmp_MS_.X_best_MS_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
if (dat_nUX_rank==dat_n_UX_rank-1); if (dat_rseed==0);
for niteration=0:dat_n_iteration-1;
subplot(4,dat_n_iteration,2*dat_n_iteration + 1+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration),32,c2d_(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('MS %d',niteration));
end;%for niteration=0:dat_n_iteration-1;
end;end;%if (dat_nUX_rank==dat_n_UX_rank-1); if (dat_rseed==0);
clear tmp_MS_;
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
tmp_SM_ = load(SM_fname_mat);
subplot(2,2,1);
hold on;
subplot(2,2,1);
plot((0:dat_n_iteration-1) + 1*dat_n_iteration,tmp_SM_.X_best_SM_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
if (dat_nUX_rank==dat_n_UX_rank-1); if (dat_rseed==0);
for niteration=0:dat_n_iteration-1;
subplot(4,dat_n_iteration,3*dat_n_iteration + 1+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration),32,c2d_(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('SM %d',niteration));
end;%for niteration=0:dat_n_iteration-1;
end;end;%if (dat_nUX_rank==dat_n_UX_rank-1); if (dat_rseed==0);
clear tmp_SM_;
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%;
end;%if flag_check;

flag_compute=0;
if flag_compute;
%%%%%%%%;
% Now set up alternating minimization using experimental principled-images, ;
% but with a variable range of translations. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 5;
dat_n_iteration = 32;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.0625,0.125,0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_pre_ = load(sprintf('%s_mat/tpmutameux_2d_xcor_d0196_ut0120le5v064_n1024_SM_nUX000rng%.3d.mat',dir_trunk,dat_rseed));
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__;
a_CTF_UX_Y_quad__ = zeros(n_lm_max,dat_n_UX_rank);
for dat_nUX_rank=0:dat_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
l_max = l_max_(1+nk_p_r);
n_lm = (l_max+1).^2;
index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm-1);
a_CTF_UX_Y_quad__(1:n_lm,1+dat_nUX_rank) = ...
a_CTF_UX_Y_quad__(1:n_lm,1+dat_nUX_rank) ...
+ UX__(1+nk_p_r,1+dat_nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+index_)*CTF_avg_k_p_r_(1+nk_p_r);
%<-- use average CTF here, under the assumption that a_CTF_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
dat_f_rand = 0.05; dat_flag_plot=0;
tpmutam_experimental_0( ...
 dat_infix ...
,dat_rseed ...
,dat_n_M ...
,dat_M_k_p__ ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,dir_trunk ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_uni_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,delta_r_max_use ...
,delta_r_max_upb ...
,svd_eps ...
,n_delta_v_requested ...
,n_CTF_rank ...
,CTF_idx_ ...
,CTF_k_p_r__ ...
,l_max_ ...
,UX__(:,1:dat_n_UX_rank) ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,tmp_pre_.euler_polar_a_SM__(:,end) ...
,tmp_pre_.euler_azimu_b_SM__(:,end) ...
,tmp_pre_.euler_gamma_z_SM__(:,end) ...
,tmp_pre_.image_delta_x_SM__(:,end) ...
,tmp_pre_.image_delta_y_SM__(:,end) ...
,tmp_pre_.image_I_value_SM__(:,end) ...
,dat_f_rand ...
,dat_flag_plot ...
);
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%;
% experimental: Collect MS and SM output. ;
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_MS____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_SM____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
tmp_ = load(MS_fname_mat);
dat_X_best_MS___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_MS_(:);
flag_dat_X_best_MS___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_MS____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_MS__(:,:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_SM_(:);
flag_dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_SM____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_SM__(:,:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_dat_X_best_MS___,'all'),numel(flag_dat_X_best_MS___)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_dat_X_best_SM___,'all'),numel(flag_dat_X_best_SM___)));
if (sum(flag_dat_X_best_MS___,'all') | sum(flag_dat_X_best_SM___,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_best_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_X_best_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_best_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_X_best_SM___(:,1+nUX_rank,:));
tmp_j_SM_ = find(sum(tmp_F_SM__,1));
if (numel(tmp_j_SM_)>0);
plot(1*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_SM__(:,tmp_j_SM_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_SM_)>0);
end;%for nUX_rank=0:dat_n_UX_rank-1;
plot((dat_n_iteration-0.5)*[1,1],[-1,1],'k-');
xlim([0,2*dat_n_iteration-1]); ylim([-0.125,1.000]);
xlabel('Phase-1 <--> Phase-2');
ylabel(sprintf('correlation'),'Interpreter','none');
grid on;
plot([0,2*dat_n_iteration-1],X_best_reco*[1,1],'k-');
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_combine,delta_r_max_use,str_delta),'Interpreter','none');
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
end;%if (sum(flag_dat_X_best_MS___,'all') | sum(flag_dat_X_best_SM___,'all'));
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%if flag_compute;

flag_compute=0;
if flag_compute;
%%%%%%%%;
% Now set up alternating minimization using synthetic principled-images. ;
%%%%%%%%;
syn_n_order = 5;
syn_n_M = 1024;
syn_n_UX_rank = 5;
syn_n_iteration = 32;
%%%%%%%%;
svd_eps_ = 10.^[-3:0.5:-1]; n_svd_eps = numel(svd_eps_);
%n_delta_v_requested_ = [1,16,32,64,128]; n_n_delta_v_requested = numel(n_delta_v_requested_);
%delta_r_max_use_ = sqrt(2)*delta_sigma*[0.00,0.05,0.15,0.50,1.00]*sqrt(log(20^2)); %<-- assumes percentile of 1/20. ;
n_delta_v_requested_ = [1,16,32]; n_n_delta_v_requested = numel(n_delta_v_requested_);
delta_r_max_use_ = sqrt(2)*delta_sigma*[0.00,0.05,0.15]*sqrt(log(20^2)); %<-- assumes percentile of 1/20. ;
syn_rseed_ = [0:2]; n_syn_rseed = numel(syn_rseed_);
syn_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_syn_snr=numel(syn_snr_);
for nsvd_eps=0:n_svd_eps-1;
svd_eps = svd_eps_(1+nsvd_eps);
for nn_delta_v_requested=0:n_n_delta_v_requested-1;
n_delta_v_requested = n_delta_v_requested_(1+nn_delta_v_requested);
delta_r_max_use = delta_r_max_use_(1+nn_delta_v_requested);
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
syn_infix = sprintf('%s_%s',str_combine,str_delta);
for nsyn_rseed=0:n_syn_rseed-1;
syn_rseed = syn_rseed_(1+nsyn_rseed);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
n_molecule = 1;
molecule_density_ = 1;
a_k_Y_quad__ = a_k_Y_quad_;
a_CTF_UX_Y_quad_mavg__ = a_CTF_UX_Y_quad__(:,1:syn_n_UX_rank);
syn_f_rand = 0.05; syn_flag_plot=0;
tpmhtam_synthetic_12(...
 syn_infix...
,syn_rseed...
,syn_n_M...
,syn_snr...
,syn_n_UX_rank...
,syn_n_iteration...
,syn_n_order...
,dir_trunk...
,n_x_u...
,diameter_x_c...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,delta_r_max_use...
,svd_eps...
,n_delta_v_requested...
,n_w_max...
,CTF_uni_avg_k_p_...
,l_max_...
,n_molecule...
,molecule_density_...
,a_k_Y_quad__...
,UX__...
,X_weight_r_...
,a_CTF_UX_Y_quad_mavg__...
,syn_f_rand...
,syn_flag_plot...
);
end;%for nsyn_snr=0:n_syn_snr-1;
end;%for nsyn_rseed=0:n_syn_rseed-1;
end;%for nn_delta_v_requested=0:n_n_delta_v_requested-1;
end;%for nsvd_eps=0:n_svd_eps-1;
%%%%%%%%;
% synthetic: Collect MS and SM output. ;
%%%%%%%%;
for nsvd_eps=0:n_svd_eps-1;
svd_eps = svd_eps_(1+nsvd_eps);
for nn_delta_v_requested=0:n_n_delta_v_requested-1;
n_delta_v_requested = n_delta_v_requested_(1+nn_delta_v_requested);
delta_r_max_use = delta_r_max_use_(1+nn_delta_v_requested);
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
syn_infix = sprintf('%s_%s',str_combine,str_delta);
%%%%%%%%;
syn_a_k_Y_0lsq_X_best__ = zeros(n_syn_snr,syn_n_UX_rank);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
for nUX_rank=0:syn_n_UX_rank-1;
fname_mat = sprintf('%s_mat/tpmhtamsux_%s_M_k_p___n%.3ds%.4dr%.3d.mat',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, skipping',fname_mat));
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, loading',fname_mat));
load(fname_mat,'syn_a_k_Y_0lsq_X_best');
syn_a_k_Y_0lsq_X_best__(1+nsyn_snr,1+nUX_rank) = syn_a_k_Y_0lsq_X_best;
clear syn_a_k_Y_0lsq_X_best;
end;%if ( exist(fname_mat,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;
end;%for nsyn_snr=0:n_syn_snr-1;
%%%%%%%%;
syn_X_best_MS____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
flag_syn_X_best_MS____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
syn_M_loading_MS_____ = zeros(3,syn_n_M,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
syn_M_molecule_type____ = zeros(syn_n_M,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
syn_X_best_SM____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
flag_syn_X_best_SM____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
syn_M_loading_SM_____ = zeros(3,syn_n_M,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
syn_M_molecule_type____ = zeros(syn_n_M,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
for nsyn_rseed=0:n_syn_rseed-1;
syn_rseed = syn_rseed_(1+nsyn_rseed);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
for nUX_rank=0:syn_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmhtamsux_%s_n%.3ds%.4d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr));
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,syn_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
tmp_ = load(MS_fname_mat);
syn_X_best_MS____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.X_best_MS_(:);
flag_syn_X_best_MS____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = 1;
syn_M_loading_MS_____(:,:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.syn_M_loading_MS__(:,:);
syn_M_molecule_type____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.syn_M_molecule_type_(:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
syn_X_best_SM____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.X_best_SM_(:);
flag_syn_X_best_SM____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = 1;
syn_M_loading_SM_____(:,:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.syn_M_loading_SM__(:,:);
syn_M_molecule_type____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.syn_M_molecule_type_(:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;
end;%for nsyn_snr=0:n_syn_snr-1;
end;%for nsyn_rseed=0:n_syn_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_syn_X_best_MS____,'all'),numel(flag_syn_X_best_MS____)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_syn_X_best_SM____,'all'),numel(flag_syn_X_best_SM____)));
if (sum(flag_syn_X_best_MS____,'all') | sum(flag_syn_X_best_SM____,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/tpmhtamsux_%s_n%.3dsxxxx_X_best_A',dir_trunk,syn_infix,syn_n_M);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
for nsyn_snr=0:n_syn_snr-1;
subplot(2,3,1+nsyn_snr);
hold on;
for nUX_rank=0:syn_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/syn_n_UX_rank)));
tmp_X_MS__ = squeeze(syn_X_best_MS____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_F_MS__ = squeeze(flag_syn_X_best_MS____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*syn_n_iteration + (0:syn_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(syn_X_best_SM____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_F_SM__ = squeeze(flag_syn_X_best_SM____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_j_SM_ = find(sum(tmp_F_SM__,1));
if (numel(tmp_j_SM_)>0);
plot(1*syn_n_iteration + (0:syn_n_iteration-1),tmp_X_SM__(:,tmp_j_SM_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_SM_)>0);
end;%for nUX_rank=0:syn_n_UX_rank-1;
plot(0:2*syn_n_iteration-1,syn_a_k_Y_0lsq_X_best__(1+nsyn_snr,1+0)*ones(1,2*syn_n_iteration),'k-');
plot((syn_n_iteration-0.5)*[1,1],[-1,1],'k-');
xlim([0,2*syn_n_iteration-1]); ylim([-0.125,1.000]);
xlabel('Phase-1 <--> Phase-2');
ylabel(sprintf('correlation'),'Interpreter','none');
grid on;
title(sprintf('syn_snr %0.2f',syn_snr_(1+nsyn_snr)),'Interpreter','none');
end;%for nsyn_snr=0:n_syn_snr-1;
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
end;%if (sum(flag_syn_X_best_MS____,'all') | sum(flag_syn_X_best_SM____,'all'));
end;%for nn_delta_v_requested=0:n_n_delta_v_requested-1;
end;%for nsvd_eps=0:n_svd_eps-1;
end;%if flag_compute;

disp('returning'); return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Appendix functions start here. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now check to see if the experimental images can themselves be used to recapitulate the principal-modes. ;
%%%%%%%%;
n_M = n_image_sub;
[emp_X_1M__] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,M_k_p__);
[emp_X_1M_U__,emp_X_1M_S__,emp_X_1M_V__] = svd(emp_X_1M__); emp_X_1M_S_ = diag(emp_X_1M_S__);
n_S = n_viewing_all;
[emp_X_CS__] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_S,S_k_p__.*repmat(CTF_avg_k_p_,[1,n_S]));
[emp_X_CS_U__,emp_X_CS_S__,emp_X_CS_V__] = svd(emp_X_CS__); emp_X_CS_S_ = diag(emp_X_CS_S__);
[tru_X_d1__] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,1,1,a_k_Y_quad_,CTF_k_p_r_xcor__,delta_sigma);
[tru_X_d1_U__,tru_X_d1_S__,tru_X_d1_V__] = svd(tru_X_d1__); tru_X_d1_S_ = diag(tru_X_d1_S__);
[tru_X_d0__] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xcor__); 
[tru_X_d0_U__,tru_X_d0_S__,tru_X_d0_V__] = svd(tru_X_d0__); tru_X_d0_S_ = diag(tru_X_d0_S__);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_X_compare_A',dir_trunk);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
figure(1);clf;
colormap(colormap_beach());
subplot(2,2,1); imagesc(tru_X_d0__); axis image; axisnotick; title('tru_X_d0__','Interpreter','none');
subplot(2,2,2); imagesc(tru_X_d1__); axis image; axisnotick; title('tru_X_d1__','Interpreter','none');
subplot(2,2,3); imagesc(emp_X_CS__); axis image; axisnotick; title('emp_X_CS__','Interpreter','none');
subplot(2,2,4); imagesc(emp_X_1M__); axis image; axisnotick; title('emp_X_1M__','Interpreter','none');
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
%%%%%%%%;
tmp_X_1M_from_d1__ = transpose(emp_X_1M_U__)*tru_X_d1__*emp_X_1M_V__;
X_1M_from_d1_sigma__ = zeros(n_k_p_r,n_k_p_r);
tmp_X_d1_from_1M__ = transpose(tru_X_d1_U__)*emp_X_1M__*tru_X_d1_V__;
X_d1_from_1M_sigma__ = zeros(n_k_p_r,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_ = svds(tmp_X_1M_from_d1__(1:1+nk_p_r,1:1+nk_p_r),1+nk_p_r);
X_1M_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r) = tmp_;
tmp_ = svds(tmp_X_d1_from_1M__(1:1+nk_p_r,1:1+nk_p_r),1+nk_p_r);
X_d1_from_1M_sigma__(1:1+nk_p_r,1+nk_p_r) = tmp_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_X_compare_B',dir_trunk);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1);
subplot(2,2,1); 
hold on;
stairs(-0.5+(1:n_k_p_r),tru_X_d1_S_,'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),X_1M_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('sigma');
hold off;
title('X_1M_from_d1_sigma__','Interpreter','none');
subplot(2,2,2); 
hold on;
stairs(-0.5+(1:n_k_p_r),emp_X_1M_S_,'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),X_d1_from_1M_sigma__(1:1+nk_p_r,1+nk_p_r),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('sigma');
hold off;
title('X_d1_from_1M_sigma__','Interpreter','none');
subplot(2,2,3); 
hold on;
stairs(-0.5+(1:n_k_p_r),log10(tru_X_d1_S_),'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),log10(X_1M_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r)),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('log10(sigma)');
hold off;
title('X_1M_from_d1_sigma__','Interpreter','none');
subplot(2,2,4); 
hold on;
stairs(-0.5+(1:n_k_p_r),log10(emp_X_1M_S_),'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),log10(X_d1_from_1M_sigma__(1:1+nk_p_r,1+nk_p_r)),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('log10(sigma)');
hold off;
title('X_d1_from_1M_sigma__','Interpreter','none');
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
% Also check to see if alternating-minimization using an empirical-cost ;
% can then be used to produce a newer (less noisy) cost .;
%%%%%%%%;
tmp_rseed=0;
fname_mat = sprintf('%s_mat/tpmhtameux_2d_emp_N1024_d0116le7v064_n1024_SM_nUX004rng%.3d.mat',dir_trunk,tmp_rseed);
tmp_emp_SM_ = load(fname_mat);
%%%%%%%%;
% calculate molecule using final image parameters. ;
%%%%%%%%;
tmp_euler_polar_a_ = tmp_emp_SM_.euler_polar_a_SM__(:,end);
tmp_euler_azimu_b_ = tmp_emp_SM_.euler_azimu_b_SM__(:,end);
tmp_euler_gamma_z_ = tmp_emp_SM_.euler_gamma_z_SM__(:,end);
tmp_image_delta_x_ = tmp_emp_SM_.image_delta_x_SM__(:,end);
tmp_image_delta_y_ = tmp_emp_SM_.image_delta_y_SM__(:,end);
tmp_n_order = 5; tmp_n_M = 1024;
tmp_t = tic;
tmp_e_k_Y_reco_ = ...
cg_lsq_4( ...
 tmp_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,tmp_n_M ...
,M_k_p__ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> tmp_e_k_Y_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
[tmp_e_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,tmp_e_k_Y_reco_);
[tmp_e_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,tmp_e_k_Y_reco_));
tmp_e_X_best_reco = max(tmp_e_X_best_orig,tmp_e_X_best_flip);
disp(sprintf(' %% tmp_e_X_best_reco %0.3f',tmp_e_X_best_reco));
%%%%%%%%;
% now calculate principal-modes of estimated molecule. ;
%%%%%%%%;
tmp_str_cost = sprintf('2d_xcor_d%.4d',floor(1000*delta_sigma));
if strfind(tmp_str_cost,'3d'); 
tmp_weight_k_p_r_ = weight_3d_k_p_r_;
end;%if strfind(tmp_str_cost,'3d'); 
if strfind(tmp_str_cost,'2d'); 
tmp_weight_k_p_r_ = weight_2d_k_p_r_;
end;%if strfind(tmp_str_cost,'2d'); 
if strfind(tmp_str_cost,'xcor');
tmp_CTF__ = CTF_k_p_r_xcor__; %<-- use CTF_k_p_r_xcor__ here to account for average cross-correlation of CTF functions. ;
end;%if strfind(tmp_str_cost,'xcor');
if strfind(tmp_str_cost,'xavg');
tmp_CTF__ = CTF_k_p_r_xavg__; %<-- use CTF_k_p_r_xavg__ here to account for mean CTF function. ;
end;%if strfind(tmp_str_cost,'xavg');
if (delta_sigma>0); [emf_X_e1__] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,tmp_weight_k_p_r_,l_max_,1,1,tmp_e_k_Y_reco_,tmp_CTF__,delta_sigma); end;
if (delta_sigma==0); [emf_X_e0__] = principled_marching_cost_matrix_3(n_k_p_r,tmp_weight_k_p_r_,l_max_max,tmp_e_k_Y_reco__,tmp_CTF__); end;
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
[emf_X_e1_U__,emf_X_e1_S__,emf_X_e1_V__] = svds(emf_X_e1__,n_k_p_r); emf_X_e1_S_ = diag(emf_X_e1_S__);

tmp_X_e1_from_d1__ = transpose(emf_X_e1_U__)*tru_X_d1__*emf_X_e1_V__;
X_e1_from_d1_sigma__ = zeros(n_k_p_r,n_k_p_r);
tmp_X_d1_from_e1__ = transpose(tru_X_d1_U__)*emf_X_e1__*tru_X_d1_V__;
X_d1_from_e1_sigma__ = zeros(n_k_p_r,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_ = svds(tmp_X_e1_from_d1__(1:1+nk_p_r,1:1+nk_p_r),1+nk_p_r);
X_e1_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r) = tmp_;
tmp_ = svds(tmp_X_d1_from_e1__(1:1+nk_p_r,1:1+nk_p_r),1+nk_p_r);
X_d1_from_e1_sigma__(1:1+nk_p_r,1+nk_p_r) = tmp_;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80s_0_X_compare_C',dir_trunk);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1);
subplot(2,2,1); 
hold on;
stairs(-0.5+(1:n_k_p_r),tru_X_d1_S_,'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),X_e1_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('sigma');
hold off;
title('X_e1_from_d1_sigma__','Interpreter','none');
subplot(2,2,2); 
hold on;
stairs(-0.5+(1:n_k_p_r),emf_X_e1_S_,'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),X_d1_from_e1_sigma__(1:1+nk_p_r,1+nk_p_r),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('sigma');
hold off;
title('X_d1_from_e1_sigma__','Interpreter','none');
subplot(2,2,3); 
hold on;
stairs(-0.5+(1:n_k_p_r),log10(tru_X_d1_S_),'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),log10(X_e1_from_d1_sigma__(1:1+nk_p_r,1+nk_p_r)),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('log10(sigma)');
hold off;
title('X_e1_from_d1_sigma__','Interpreter','none');
subplot(2,2,4); 
hold on;
stairs(-0.5+(1:n_k_p_r),log10(emf_X_e1_S_),'k-','LineWidth',2);
for nk_p_r=0:n_k_p_r-1;
nc = max(0,min(n_c-1,floor(n_c*nk_p_r/n_k_p_r)));
stairs(-0.5+(1:1+nk_p_r),log10(X_d1_from_e1_sigma__(1:1+nk_p_r,1+nk_p_r)),'-','Color',c_(1+nc,:),'LineWidth',2);
end;%for nk_p_r=0:n_k_p_r-1;
xlim([0.5,0.5+n_k_p_r]);
xlabel('rank'); ylabel('log10(sigma)');
hold off;
title('X_d1_from_e1_sigma__','Interpreter','none');
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
% Now load one of the files from the experimental runs. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 5;
dat_n_iteration = 32;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
ndelta_r_max_factor=4; %ndelta_r_max_factor=5; %<-- the best choice seems to be slightly larger than the actual set of translations. ;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
ndat_rseed=0;
dat_rseed = dat_rseed_(1+ndat_rseed);
%%%%%%%%;
dat_X_best_MS__ = zeros(dat_n_iteration,dat_n_UX_rank);
flag_dat_X_best_MS__ = zeros(dat_n_iteration,dat_n_UX_rank);
dat_X_best_SM__ = zeros(dat_n_iteration,dat_n_UX_rank);
flag_dat_X_best_SM__ = zeros(dat_n_iteration,dat_n_UX_rank);
for nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
tmp_ = load(MS_fname_mat);
dat_X_best_MS__(:,1+nUX_rank) = tmp_.X_best_MS_(:);
flag_dat_X_best_MS__(:,1+nUX_rank) = 1;
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM__(:,1+nUX_rank) = tmp_.X_best_SM_(:);
flag_dat_X_best_SM__(:,1+nUX_rank) = 1;
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
%%%%%%%%;
% Now rerun with different values for f_rand. ;
%%%%%%%%;
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__;
% Here we copy over the MS file. ;
fname_mat_pre = sprintf('%s_mat/tpmhtameux_%s_M_k_p___n%.3dr%.3d.mat',dir_trunk,dat_infix,dat_n_M,dat_n_UX_rank);
f_rand_ = [0.00,0.01,0.025,0.05,0.10,0.15]; n_f_rand = numel(f_rand_);
for nf_rand=0:n_f_rand-1;
f_rand = f_rand_(1+nf_rand);
dat_infix_pos = sprintf('%s_f%.3d',dat_infix,floor(100*f_rand));
fname_mat_pos = sprintf('%s_mat/tpmhtameux_%s_M_k_p___n%.3dr%.3d.mat',dir_trunk,dat_infix_pos,dat_n_M,dat_n_UX_rank);
string_command = sprintf(' ln %s %s;\n',fname_mat_pre,fname_mat_pos);
disp(string_command); system(string_command);
fname_0 = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
for nUX_rank=0:dat_n_UX_rank-1;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat_pre = sprintf('%s.mat',MS_fname_pre);
fname_0_pos = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix_pos,dat_n_M);
MS_fname_pre_pos = sprintf('%s_MS_%s',fname_0_pos,fname_2);
MS_fname_mat_pos = sprintf('%s.mat',MS_fname_pre_pos);
string_command = sprintf(' ln %s %s;\n',MS_fname_mat_pre,MS_fname_mat_pos);
disp(string_command); system(string_command);
end;%for nUX_rank=0:dat_n_UX_rank-1;
disp(sprintf(' %% calling tmpmhtam_experimental_10 with dat_infix_pos %s',dat_infix_pos));
flag_plot=0;
tpmhtam_experimental_11(...
 dat_infix_pos...
,dat_rseed...
,dat_n_M...
,dat_M_k_p__...
,dat_n_UX_rank...
,dat_n_iteration...
,dat_n_order...
,dir_trunk...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,n_w_uni_...
,weight_3d_k_p_r_...
,weight_2d_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,delta_r_max_use...
,svd_eps...
,n_delta_v_requested...
,CTF_uni_avg_k_p_...
,CTF_uni_avg_k_p_r_...
,l_max_...
,UX__...
,X_weight_r_...
,a_CTF_UX_Y_quad__(:,1:dat_n_UX_rank)...
,f_rand...
,flag_plot...
);
end;%for nf_rand=0:n_f_rand-1;
%%%%%%%%;
% Now collect the results. ;
%%%%%%%%;
X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_f_rand);
flag_found__ = zeros(dat_n_UX_rank,n_f_rand);
for nf_rand=0:n_f_rand-1;
f_rand = f_rand_(1+nf_rand);
dat_infix_pos = sprintf('%s_f%.3d',dat_infix,floor(100*f_rand));
for nUX_rank=0:dat_n_UX_rank-1;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
fname_0_pos = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix_pos,dat_n_M);
SM_fname_pre_pos = sprintf('%s_SM_%s',fname_0_pos,fname_2);
SM_fname_mat_pos = sprintf('%s.mat',SM_fname_pre_pos);
if ( exist(SM_fname_mat_pos,'file'));
tmp_ = load(SM_fname_mat_pos);
X_best_SM___(:,1+nUX_rank,1+nf_rand) = tmp_.X_best_SM_;
flag_found__(1+nUX_rank,1+nf_rand) = 1;
clear tmp_;
end;%if ( exist(SM_fname_mat_pos,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for nf_rand=0:n_f_rand-1;
%%%%%%%%;
% plot fig by f_rand. ;
%%%%%%%%;
dat_infix_pos = sprintf('%s_fxxx',dat_infix);
fname_fig = sprintf('%s_jpg/tpmhtameux_%s_n%.3d_X_best_SM_FIGA',dir_trunk,dat_infix_pos,dat_n_M);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;
c_ = colormap_beach(); n_c=size(c_,1);
for nf_rand=0:n_f_rand-1;
f_rand = f_rand_(1+nf_rand);
subplot(2,3,1+nf_rand); hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(nUX_rank/dat_n_UX_rank*n_c)));
if (flag_found__(1+nUX_rank,1+nf_rand));
plot(1:dat_n_iteration,X_best_SM___(:,1+nUX_rank,1+nf_rand),'x-','Color',c_(1+nc,:),'LineWidth',2);
plot(1:dat_n_iteration,X_best_SM___(:,1+nUX_rank,1+nf_rand),'ko');
xlabel('SM iteration'); ylabel('correlation');
xlim([1,dat_n_iteration]); ylim([-0.1,1.0]); grid on;
title(sprintf('f_rand %0.6f',f_rand),'Interpreter','none');
end;%if (flag_found__(1+nUX_rank,1+nf_rand));
end;%for nUX_rank=0:dat_n_UX_rank-1;
hold off;
end;%for nf_rand=0:n_f_rand-1;
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
% plot fig by dat_n_UX_rank. ;
%%%%%%%%;
dat_infix_pos = sprintf('%s_fxxx',dat_infix);
fname_fig = sprintf('%s_jpg/tpmhtameux_%s_n%.3d_X_best_SM_FIGB',dir_trunk,dat_infix_pos,dat_n_M);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(2);clf;
c_ = colormap_beach(); n_c=size(c_,1);
for nUX_rank=0:dat_n_UX_rank-1;
subplot(2,3,1+nUX_rank); hold on;
for nf_rand=0:n_f_rand-1;
f_rand = f_rand_(1+nf_rand);
nc = max(0,min(n_c-1,floor(nf_rand/n_f_rand*n_c)));
if (flag_found__(1+nUX_rank,1+nf_rand));
plot(1:dat_n_iteration,X_best_SM___(:,1+nUX_rank,1+nf_rand),'x-','Color',c_(1+nc,:),'LineWidth',2);
plot(1:dat_n_iteration,X_best_SM___(:,1+nUX_rank,1+nf_rand),'ko');
xlabel('SM iteration'); ylabel('correlation');
xlim([1,dat_n_iteration]); ylim([-0.1,1.0]); grid on;
title(sprintf('nUX_rank %d',nUX_rank),'Interpreter','none');
end;%if (flag_found__(1+nUX_rank,1+nf_rand));
end;%for nUX_rank=0:dat_n_UX_rank-1;
hold off;
end;%for nf_rand=0:n_f_rand-1;
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
% Now visualize X_wSM___ for a few of the images. ;
%%%%%%%%;
nf_rand=0; %<-- this will not be used here. ;
f_rand = f_rand_(1+nf_rand);
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__;
dat_infix_pos = sprintf('%s_f%.3d',dat_infix,floor(100*f_rand));
fname_mat_pos = sprintf('%s_mat/tpmhtameux_%s_M_k_p___n%.3dr%.3d.mat',dir_trunk,dat_infix_pos,dat_n_M,dat_n_UX_rank);
tmp_0_ = load(fname_mat_pos);
nUX_rank = 3; %for nUX_rank=0:dat_n_UX_rank-1;

fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
fname_0_pos = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix_pos,dat_n_M);
MS_fname_pre_pos = sprintf('%s_MS_%s',fname_0_pos,fname_2);
MS_fname_mat_pos = sprintf('%s.mat',MS_fname_pre_pos);
tmp_1_ = load(MS_fname_mat_pos);
%%%%%%%%;
tmp_pm_n_UX_rank = 1+nUX_rank;
tmp_pm_n_order = dat_n_order;
tmp_pm_n_k_p_r = tmp_pm_n_UX_rank;
tmp_pm_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_k_p_r_max = 1;
tmp_pm_n_w_ = n_w_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_w_max = n_w_max;
tmp_pm_n_w_sum = sum(tmp_pm_n_w_);
tmp_pm_n_w_csum_ = cumsum([0;tmp_pm_n_w_]);
tmp_pm_l_max_ = l_max_max*ones(tmp_pm_n_k_p_r,1);
tmp_pm_n_lm_ = (1+tmp_pm_l_max_).^2; tmp_pm_n_lm_sum = sum(tmp_pm_n_lm_);
tmp_pm_weight_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
tmp_pm_weight_2d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
%tmp_pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ;
tmp_pm_CTF_index_ = 1; tmp_pm_CTF_k_p__ = ones(tmp_pm_n_w_sum,1);
%%%%%%%%;
% Find out which on-grid displacements correspond most closely to current displacements. ;
% Then use current displacements to form principal-images. ;
%%%%%%%%;
tmp_t = tic();
[tmp_UX_M_k_q_wnM___,tmp_UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(tmp_0_.FTK,n_w_uni_,tmp_pm_n_UX_rank,dat_n_M,tmp_0_.svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:),tmp_1_.image_delta_x_MS__(:,end),tmp_1_.image_delta_y_MS__(:,end));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
tmp_a_CTF_UX_Y_0lsq_ = ...
cg_lsq_3( ...
 tmp_pm_n_order ...
,tmp_pm_n_k_p_r ...
,tmp_pm_l_max_ ...
,tmp_pm_n_w_ ...
,dat_n_M ...
,reshape(tmp_UX_M_k_p_wnM___,[tmp_pm_n_w_max*tmp_pm_n_k_p_r,dat_n_M]) ...
,tmp_pm_CTF_index_ ...
,tmp_pm_CTF_k_p__ ...
,tmp_1_.euler_polar_a_MS__(:,end) ...
,tmp_1_.euler_azimu_b_MS__(:,end) ...
,tmp_1_.euler_gamma_z_MS__(:,end) ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% cg_lsq_3 for a_CTF_UX_Y_0lsq_: %0.3fs',tmp_t)); end;
disp(sprintf(' %% fnorm(tmp_a_CTF_UX_Y_0lsq_ - tmp_1_.a_CTF_UX_Y_0lsq_MS__(:,end))/fnorm(tmp_a_CTF_UX_Y_0lsq_) = %0.16f',fnorm(tmp_a_CTF_UX_Y_0lsq_ - tmp_1_.a_CTF_UX_Y_0lsq_MS__(:,end))/fnorm(tmp_a_CTF_UX_Y_0lsq_)));
%%%%%%%%;
% use current model to generate current principal-templates. ;
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
[tmp_S_k_p__,~,~,~,tmp_n_viewing_all,tmp_viewing_azimu_b_all_,tmp_viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,tmp_pm_n_k_p_r,tmp_pm_k_p_r_,tmp_pm_k_p_r_max,tmp_pm_weight_k_p_r_,tmp_pm_l_max_,tmp_a_CTF_UX_Y_0lsq_,tmp_viewing_k_eq_d,-1,tmp_pm_n_w_);
if (tmp_verbose>0); disp(sprintf(' %% tmp_viewing_k_eq_d %0.3f, tmp_n_viewing_all %d',tmp_viewing_k_eq_d,tmp_n_viewing_all)); end;
tmp_n_S = tmp_n_viewing_all;
tmp_UX_S_l2_ = zeros(tmp_n_S,1);
for nS=0:tmp_n_S-1;
tmp_UX_S_l2_(1+nS) = innerproduct_p_quad(tmp_pm_n_k_p_r,tmp_pm_k_p_r_,tmp_pm_weight_2d_k_p_r_/(2*pi),tmp_pm_n_w_,tmp_pm_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:tmp_n_S-1;
tmp_S_k_q__ = zeros(tmp_pm_n_w_sum,tmp_n_S);
for nS=0:tmp_n_S-1;
tmp_S_k_q__(:,1+nS) = interp_p_to_q(tmp_pm_n_k_p_r,tmp_pm_n_w_,tmp_pm_n_w_sum,tmp_S_k_p__(:,1+nS)); 
end;%for nS=0:tmp_n_S-1; 
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_template_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
tmp_t = tic();
tmp_SS_k_q_ = svd(tmp_S_k_q__);
tmp_n_S_rank = min(find(tmp_SS_k_q_/tmp_SS_k_q_(1)<1e-3)); if isempty(tmp_n_S_rank); tmp_n_S_rank = min(size(tmp_S_k_q__)); end;
[tmp_US_k_q__,tmp_SS_k_q__,tmp_VS_k_q__] = svds(tmp_S_k_q__,tmp_n_S_rank);
if (verbose); disp(sprintf(' %% tmp_n_S %d --> tmp_n_S_rank %d',tmp_n_S,tmp_n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
% Batches images into batches of size n_M_Mbatch. ;
% Batches templates into batches of size n_S_Sbatch. ;
% Only stores the optimal translation for each image. ;
%%%%%%%%;
tmp_t = tic();
tmp_UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(tmp_0_.FTK,n_w_uni_,dat_n_M,tmp_pm_n_UX_rank,tmp_0_.svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
tmp_t = tic();
tmp_n_M_Mbatch = 24;
tmp_n_S_Sbatch = 24;
[tmp_X_wSM___,tmp_delta_j_wSM___] = ...
ampmh_X_wSM___1(...
 tmp_0_.FTK...
,n_w_uni_...
,tmp_pm_n_UX_rank...
,tmp_n_S...
,tmp_n_S_rank...
,tmp_n_S_Sbatch...
,tmp_US_k_q__...
,tmp_SS_k_q__...
,tmp_VS_k_q__...
,tmp_UX_S_l2_...
,dat_n_M...
,tmp_n_M_Mbatch...
,tmp_0_.svd_VUXM_lwnM____(:,:,1:tmp_pm_n_UX_rank,:)...
,tmp_UX_M_l2_dM__...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% tmp_X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% We could save tmp_X_wSM___, but it is almost 2GB in size. ;
% Instead just look at a few pictures. ;
%%%%%%%%;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
fname_fig = sprintf('%s_jpg/tpmhtameux_%s_n%.3d_MS_%s_final_X_S_',dir_trunk,dat_infix,dat_n_M,fname_2);
if (1 | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;
%tmp_lim_X_wSM__ = [min(real(tmp_X_wSM___(:))),max(real(tmp_X_wSM___(:)))];
n_plot=15;
for nplot=0:n_plot-1;
subplot(3,5,1+nplot);
nM = max(0,min(dat_n_M-1,floor(nplot/n_plot*dat_n_M)));
tmp_X_S_ = max(real(tmp_X_wSM___(:,:,1+nM)),[],1);
%tmp_lim_ = tmp_lim_X_wSM__; %tmp_lim_ = [-1,1];
tmp_lim_ = [min(tmp_X_S_(:)),max(tmp_X_S_(:))];
hold on; 
imagesc_polar_a_azimu_b_0(tmp_viewing_polar_a_all_,tmp_viewing_azimu_b_all_,tmp_X_S_,tmp_lim_);
[~,tmp_index] = max(tmp_X_S_); tmp_index = tmp_index - 1;
plot(tmp_viewing_azimu_b_all_(1+tmp_index),tmp_viewing_polar_a_all_(1+tmp_index),'kx','MarkerSize',10); 
plot(tmp_viewing_azimu_b_all_(1+tmp_index),tmp_viewing_polar_a_all_(1+tmp_index),'ko','MarkerSize',10); 
hold off;
xlim([0,2*pi]); ylim([0,1*pi]); axisnotick; 
title(sprintf('nM %.4d X %0.3f',nM,max(tmp_X_S_(:))));
end;%for nplot=0:n_plot-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%end;%for nUX_rank=0:dat_n_UX_rank-1;

%%%%%%%%;
% Now try and determine an appropriate (viewing-angle-specific) estimate for f_rand. ;
%%%%%%%%;
nUX_rank = 3;
fname_2 = sprintf('nUX%.3drng%.3d',nUX_rank,dat_rseed);
fname_0 = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
tmp_1_ = load(MS_fname_mat);
%%%%%%%%;
% First use estimated viewing-angles after MS_phase to reconstruct full molecule in k_p. ;
%%%%%%%%;
tmp_euler_polar_a_ = tmp_1_.euler_polar_a_MS__(:,end);
tmp_euler_azimu_b_ = tmp_1_.euler_azimu_b_MS__(:,end);
tmp_euler_gamma_z_ = tmp_1_.euler_gamma_z_MS__(:,end);
tmp_image_delta_x_ = tmp_1_.image_delta_x_MS__(:,end);
tmp_image_delta_y_ = tmp_1_.image_delta_y_MS__(:,end);
tmp_t = tic;
tmp_d_k_Y_reco_ = ...
cg_lsq_4( ...
 dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,l_max_ ...
,n_w_ ...
,dat_n_M ...
,M_k_p__ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_ ...
,tmp_euler_azimu_b_ ...
,tmp_euler_gamma_z_ ...
,tmp_image_delta_x_ ...
,tmp_image_delta_y_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> tmp_d_k_Y_reco_ time %0.2fs',tmp_t));
%%%%%%%%;
[tmp_d_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,tmp_d_k_Y_reco_);
[tmp_d_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,tmp_d_k_Y_reco_));
tmp_d_X_best_reco = max(tmp_d_X_best_orig,tmp_d_X_best_flip);
disp(sprintf(' %% tmp_d_X_best_reco %0.3f',tmp_d_X_best_reco));

%%%%%%%%;
% now calculate principal-modes of estimated molecule. ;
%%%%%%%%;

if strfind(str_cost,'3d'); 
tmp_weight_k_p_r_ = weight_3d_k_p_r_;
end;%if strfind(str_cost,'3d'); 
if strfind(str_cost,'2d'); 
tmp_weight_k_p_r_ = weight_2d_k_p_r_;
end;%if strfind(str_cost,'2d'); 
if strfind(str_cost,'xcor');
tmp_CTF__ = CTF_k_p_r_xcor__; %<-- use CTF_k_p_r_xcor__ here to account for average cross-correlation of CTF functions. ;
end;%if strfind(str_cost,'xcor');
if strfind(str_cost,'xavg');
tmp_CTF__ = CTF_k_p_r_xavg__; %<-- use CTF_k_p_r_xavg__ here to account for mean CTF function. ;
end;%if strfind(str_cost,'xavg');
if (delta_sigma>0); [tmp_d_X__,tmp_d_X_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,tmp_weight_k_p_r_,l_max_,1,1,tmp_d_k_Y_reco_,tmp_CTF__,delta_sigma); end;
if (delta_sigma==0); [tmp_d_X__,tmp_d_X_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,tmp_weight_k_p_r_,l_max_max,tmp_d_k_Y_reco__,tmp_CTF__); end;
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
[tmp_d_UX__,tmp_d_SX__,tmp_d_VX__] = svds(tmp_d_X__,est_n_UX_rank);
%%%%%%%%;
% Now pass estimated molecule to tpmhtamsf_0. ;
%%%%%%%%;

est_rseed = 0;
est_snr = 0.10;
est_n_M = dat_n_M;
est_n_UX_rank = dat_n_UX_rank;
est_delta_r_max = 1*delta_r_max_use;%est_delta_r_max = delta_r_max_use;
[...
 est_prctile__...
,est_vdist__...
,tmp_n_viewing_all...
,tmp_viewing_azimu_b_all_...
,tmp_viewing_polar_a_all_...
,true_viewing_polar_a_...
,true_viewing_azimu_b_...
,true_viewing_gamma_z_...
,true_viewing_delta_x_...
,true_viewing_delta_y_...
,est_M_S_index_...
] = ...
tpmhtamsf_0(...
 est_rseed...
,est_n_M...
,est_snr...
,est_n_UX_rank...
,dir_trunk...
,n_x_u...
,diameter_x_c...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_3d_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,est_delta_r_max...
,svd_eps...
,n_delta_v_requested...
,n_w_max...
,CTF_uni_avg_k_p_...
,l_max_...
,tmp_d_k_Y_reco_...
,tmp_d_UX__...
,tmp_d_X_weight_r_...
);


disp('returning'); return;
