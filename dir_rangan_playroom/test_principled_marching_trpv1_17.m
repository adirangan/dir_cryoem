% Here str_cost = '3d' corrects for volume element, ;
% whereas str_cost = '2d' corrects for area element. ;
% Here str_CTF_factor = 'rescaled' for CTF-rescaled images, ;
% whereas str_CTF_factor = 'original' for original CTF used in images. ;
% Similar to test_principled_marching_trpv1_13, except that: ;
% 1. the images are *NOT* pretranslated, and ;
% 2. the delta_sigma used to calculate the cost is *NOT* set to 0. ;
% Moreover, the n_w_ used to calculate M_k_p__ is not adaptive. ;
% Now we investigate the role played by f_rand. ;
clear; setup;

str_CTF_factor = 'original';

dir_trunk = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching');

h2d_ = @(kd) (2*pi)^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^2;
dh2d_ = @(kd) (2*pi)^2*0.5*(besselj(-1,kd) - besselj(+3,kd));
h3d_ = @(kd) (2*pi)^3*3*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = (2*pi)^3;
dh3d_ = @(kd) (2*pi)^3*9*( (kd.^2/3 - 1) .* sin(kd) + (kd).*cos(kd) ) ./ kd.^4 ;
flag_check=0;
if flag_check;
x_ = transpose(linspace(0,5,1024)); y_ = 0.5*(x_(2:end)+x_(1:end-1));
subplot(1,2,1); plot(y_,dh2d_(y_),'k-',y_,diff(h2d_(x_))/mean(diff(x_)),'ro'); title('dh2d'); xlabel('kd');
subplot(1,2,2); plot(y_,dh3d_(y_),'k-',y_,diff(h3d_(x_))/mean(diff(x_)),'ro'); title('dh3d'); xlabel('kd');
end;%if flag_check;

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
% First load trpv1 molecule on x_u grid. ;
%%%%%%%%;
dir_data = '/data/rangan/dir_cryoem/dir_trpv1/data_nosym';
fname_dims = sprintf('%s/dims',dir_data);
tmp_ = textread(fname_dims); n_x_u = tmp_(1); n_image = tmp_(2); clear tmp_;
fname_density = sprintf('%s/density_clean',dir_data);
a_x_u_load_ = textread(fname_density); a_x_u_load_ = reshape(a_x_u_load_,n_x_u,n_x_u,n_x_u);
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_a_k_p_quad_.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic;
k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_,J_node_,J_weight_,J_chebfun_,J_polyval_] = sample_sphere_7(verbose,k_p_r_max,k_eq_d,'L') ;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_a_k_p_quad_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_b_x_u_reco__.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_x_u_reco__ = zeros(n_X_u,n_k_p_r);
eta = pi/k_p_r_max; 
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+1+nk_p_r)-1;
tmp_n_k_all = numel(1+tmp_ij_);
tmp_t = tic;
a_x_u_reco__(:,1+nk_p_r) = nufft3d3(tmp_n_k_all,2*pi*k_c_0_all_(1+tmp_ij_)*eta,2*pi*k_c_1_all_(1+tmp_ij_)*eta,2*pi*k_c_2_all_(1+tmp_ij_)*eta,a_k_p_quad_(1+tmp_ij_).*(2*pi)^3.*weight_k_all_(1+tmp_ij_),+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_b_x_u_reco__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_a_k_Y_quad_.mat',dir_trunk);
if (~exist(fname_mat,'file'));
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
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
Y_l_val_(1+tmp_ij_) = tmp_l_val_;
Y_m_val_(1+tmp_ij_) = tmp_m_val_;
Y_k_val_(1+tmp_ij_) = k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_(1+tmp_ij_) = weight_3d_k_p_r_(1+nk_p_r);
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
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_(1+tmp_ij_);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_a_k_Y_quad_A',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_a_k_Y_quad_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_S_k_p__.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
template_k_eq_d = k_eq_d*2;
viewing_k_eq_d = k_eq_d*8;
[S_k_p__,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,l_max_,a_k_Y_quad_,viewing_k_eq_d,template_k_eq_d);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_S_k_p__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_S_k_p__A',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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

flag_check=0;
if flag_check;
%%%%%%%%;
% Now test least-squares solver for some of the shells. ;
%%%%%%%%;
for nk_p_r=0:8:n_k_p_r-1;%for nk_p_r=0:n_k_p_r-1;
tmp_S_k_p__ = S_k_p__(1+n_w_csum_(1+nk_p_r)+[0:n_w_(1+nk_p_r)-1],:);
viewing_gamma_z_all_ = 2*pi*rand(n_viewing_all,1); tmp_M_k_p__ = zeros(n_w_(1+nk_p_r),n_viewing_all);
for nviewing_all=0:n_viewing_all-1; tmp_M_k_p__(:,1+nviewing_all) = rotate_p2p_fx(1,n_w_(1+nk_p_r),n_w_(1+nk_p_r),tmp_S_k_p__(:,1+nviewing_all),+viewing_gamma_z_all_(1+nviewing_all)); end; %<-- Note that the sign of viewing_gamma_z is the same as the next line. ;
[tmp_k_p_polar_a__,tmp_k_p_azimu_b__] = cg_rhs_1(n_viewing_all,n_w_(1+nk_p_r),viewing_polar_a_all_,viewing_azimu_b_all_,+viewing_gamma_z_all_); %<-- Note that the sign of viewing_gamma_z is the same as the previous line. ;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_polar_a = max(15,1+2*tmp_l_max); tmp_n_azimu_b = 1+2*tmp_n_polar_a;
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(tmp_l_max,cos(linspace(0,pi,tmp_n_polar_a)),tmp_n_azimu_b);
tmp_n_order = 5;
tensor_to_scatter__ = cg_interpolate_n_1(tmp_n_order,tmp_n_polar_a,tmp_n_azimu_b,n_w_(1+nk_p_r)*n_viewing_all,tmp_k_p_polar_a__(:),tmp_k_p_azimu_b__(:));
scatter_to_tensor__ = transpose(tensor_to_scatter__);
tmp_An__ = @(a_k_Y_) tensor_to_scatter__*reshape(cg_evaluate_n_1(tmp_l_max,convert_spharm_to_spharm__0(tmp_l_max,a_k_Y_),tmp_n_polar_a,tmp_n_azimu_b,legendre_evaluate_ljm___),[tmp_n_polar_a*tmp_n_azimu_b,1]);
tmp_At__ = @(a_k_X_) convert_spharm__to_spharm_0(tmp_l_max,cg_evaluate_t_1(tmp_n_polar_a,tmp_n_azimu_b,reshape(scatter_to_tensor__*a_k_X_,[tmp_n_polar_a,tmp_n_azimu_b]),tmp_l_max,legendre_evaluate_mlj___,expil__,expi__));
tmp_AtAn__ = @(a_k_Y_) tmp_At__(tmp_An__(a_k_Y_));
tmp_a_k_Y_0lsq_ = pcg(tmp_AtAn__,tmp_At__(tmp_M_k_p__(:)));
tmp_a_k_Y_quad_ = a_k_Y_quad_(1+n_lm_csum_(1+nk_p_r)+[0:n_lm_(1+nk_p_r)-1]);
disp(sprintf(' %% lsq: nk_p_r %d/%d, tmp_a_k_Y_quad_ vs tmp_a_k_Y_0lsq_: %0.16f',nk_p_r,n_k_p_r,fnorm(tmp_a_k_Y_quad_-tmp_a_k_Y_0lsq_)/fnorm(tmp_a_k_Y_quad_)));
end;%for nk_p_r=0:n_k_p_r-1;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now test out alternating minimization using images copied from templates ;
% generated using just one of the shells (with no added noise). ;
%%%%%%%%;
tmp_verbose=0;
nk_p_r = 16; tmp_n_k_p_r = 1;
tmp_a_k_Y_quad_ = a_k_Y_quad_(1+n_lm_csum_(1+nk_p_r)+[0:n_lm_(1+nk_p_r)-1]); %<-- use as ground truth for this shell. ;
tmp_k_p_r = k_p_r_(1+nk_p_r); tmp_n_w = n_w_(1+nk_p_r); tmp_n_w_max = tmp_n_w; tmp_n_w_sum = tmp_n_w; 
tmp_weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r); tmp_weight_2d_k_p_r = 1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_polar_a = max(15,1+2*tmp_l_max); tmp_n_azimu_b = 1+2*tmp_n_polar_a;
[legendre_evaluate_ljm___,legendre_evaluate_mlj___,expil__,expi__] = legendre_evaluate_ljm___0(tmp_l_max,cos(linspace(0,pi,tmp_n_polar_a)),tmp_n_azimu_b);
tmp_n_order = 5;
tmp_S_k_p__ = S_k_p__(1+n_w_csum_(1+nk_p_r)+[0:tmp_n_w-1],:);
viewing_gamma_z_all_ = 2*pi*rand(n_viewing_all,1); 
M_k_p__ = zeros(n_w_sum,n_viewing_all);
tmp_M_k_p__ = zeros(tmp_n_w,n_viewing_all);
tmp_M_k_q__ = zeros(tmp_n_w,n_viewing_all);
for nviewing_all=0:n_viewing_all-1; 
M_k_p__(:,1+nviewing_all) = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,S_k_p__(:,1+nviewing_all),+viewing_gamma_z_all_(1+nviewing_all));
tmp_M_k_p__(:,1+nviewing_all) = rotate_p2p_fx(tmp_n_k_p_r,tmp_n_w,tmp_n_w_sum,S_k_p__(1+n_w_csum_(1+nk_p_r)+[0:tmp_n_w-1],1+nviewing_all),+viewing_gamma_z_all_(1+nviewing_all));
tmp_M_k_q__(:,1+nviewing_all) = interp_p_to_q(tmp_n_k_p_r,tmp_n_w,tmp_n_w_sum,tmp_M_k_p__(:,1+nviewing_all));
end;%for nviewing_all=0:n_viewing_all-1; 
n_M = n_viewing_all; 
tmp_n_M = n_viewing_all;
tmp_MM_ = zeros(tmp_n_M,1);
for nM=0:tmp_n_M-1;
tmp_MM_(1+nM) = innerproduct_p_quad(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_M_k_p__(:,1+nM),tmp_M_k_p__(:,1+nM));
end;%for nM=0:tmp_n_M-1;
%%%%%%%%;
% initialize current euler-angles randomly. ;
%%%%%%%%;
euler_polar_a_ = pi*rand(tmp_n_M,1); euler_azimu_b_ = 2*pi*rand(tmp_n_M,1); euler_gamma_z_ = 2*pi*rand(tmp_n_M,1);
n_iteration = 32;
for niteration=0:n_iteration-1;
%%%%%%%%;
% use current euler-angles to solve for current model (on single shell). ;
%%%%%%%%;
[tmp_k_p_polar_a__,tmp_k_p_azimu_b__] = cg_rhs_1(tmp_n_M,tmp_n_w,euler_polar_a_,euler_azimu_b_,+euler_gamma_z_);
tensor_to_scatter__ = cg_interpolate_n_1(tmp_n_order,tmp_n_polar_a,tmp_n_azimu_b,tmp_n_w*tmp_n_M,tmp_k_p_polar_a__(:),tmp_k_p_azimu_b__(:));
scatter_to_tensor__ = transpose(tensor_to_scatter__);
tmp_An__ = @(a_k_Y_) tensor_to_scatter__*reshape(cg_evaluate_n_1(tmp_l_max,convert_spharm_to_spharm__0(tmp_l_max,a_k_Y_),tmp_n_polar_a,tmp_n_azimu_b,legendre_evaluate_ljm___),[tmp_n_polar_a*tmp_n_azimu_b,1]);
tmp_At__ = @(a_k_X_) convert_spharm__to_spharm_0(tmp_l_max,cg_evaluate_t_1(tmp_n_polar_a,tmp_n_azimu_b,reshape(scatter_to_tensor__*a_k_X_,[tmp_n_polar_a,tmp_n_azimu_b]),tmp_l_max,legendre_evaluate_mlj___,expil__,expi__));
tmp_AtAn__ = @(a_k_Y_) tmp_At__(tmp_An__(a_k_Y_));
[tmp_a_k_Y_0lsq_,~] = pcg(tmp_AtAn__,tmp_At__(tmp_M_k_p__(:)));
%%%%%%%%;
% Compare current model (on single shell) to tmp_a_k_Y_quad_ (also on single shell). ;
%%%%%%%%;
[X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(tmp_n_k_p_r,tmp_k_p_r,k_p_r_max,tmp_weight_3d_k_p_r,0,tmp_l_max,tmp_a_k_Y_quad_,tmp_a_k_Y_0lsq_);
disp(sprintf(' %% [single shell]: tmp_a_k_Y_quad_ vs tmp_a_k_Y_lsq0_: correlation %+0.6f',X_best));
if (mod(niteration,8)==7);
%%%%%%%%;
% use current euler-angles to solve for current model (across all shells). ;
%%%%%%%%;
a_k_Y_0lsq_ = cg_lsq_1(tmp_n_order,n_k_p_r,l_max_,n_w_,n_M,M_k_p__,euler_polar_a_,euler_azimu_b_,euler_gamma_z_);
%%%%%%%%;
% Compare current model (across all shells) to tmp_a_k_Y_quad_ (across all shells). ;
%%%%%%%%;
[X_best,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,a_k_Y_0lsq_);
disp(sprintf(' %% [across all shells] a_k_Y_quad_ vs a_k_Y_0lsq_: correlation %+0.6f',X_best));
end;%if (mod(niteration,8)==7);
%%%%%%%%;
% use current model (on single shell) to generate current templates. ;
%%%%%%%%;
tmp_viewing_k_eq_d = k_eq_d*9; %<-- make this slightly different from viewing_k_eq_d to test out code. ;
[tmp_S_k_p__,~,~,~,tmp_n_viewing_all,tmp_viewing_azimu_b_all_,tmp_viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,tmp_n_k_p_r,tmp_k_p_r,k_p_r_max,tmp_weight_3d_k_p_r,tmp_l_max,tmp_a_k_Y_0lsq_,tmp_viewing_k_eq_d,-1,tmp_n_w);
tmp_n_S = tmp_n_viewing_all;
tmp_SS_ = zeros(tmp_n_S,1);
for nS=0:tmp_n_S-1;
tmp_SS_(1+nS) = innerproduct_p_quad(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:tmp_n_S-1;
tmp_S_k_q__ = zeros(tmp_n_w_sum,tmp_n_S);
for nS=0:tmp_n_S-1; tmp_S_k_q__(:,1+nS) = interp_p_to_q(tmp_n_k_p_r,tmp_n_w,tmp_n_w_sum,tmp_S_k_p__(:,1+nS)); end;%for nS=0:tmp_n_S-1; 
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
%%%%%%%%;
tmp_X___ = zeros(tmp_n_w_max,tmp_n_M,tmp_n_S);
for nS=0:tmp_n_S-1;
for nM=0:tmp_n_M-1;
tmp_X___(:,1+nM,1+nS) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r,tmp_weight_2d_k_p_r/(2*pi),tmp_n_w,tmp_n_w_sum,tmp_S_k_q__(:,1+nS),tmp_M_k_q__(:,1+nM)))*tmp_n_w_max/sqrt(tmp_MM_(1+nM))/sqrt(tmp_SS_(1+nS)) ; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
end;%for nM=0:tmp_n_M-1;
end;%for nS=0:tmp_n_S-1;
%%%%%%%%;
% Use current correlations to udate current euler-angles. ;
%%%%%%%%;
flag_M_used_ = zeros(tmp_n_M,1);
tmp_permutation_ = randperm(tmp_n_S)-1;
nS=0; nM_sum=0;
while (sum(flag_M_used_)<tmp_n_M);
index_M_unused_ = find(flag_M_used_==0)-1;
[~,index_wM_best] = max(real(tmp_X___(:,1+index_M_unused_,1+tmp_permutation_(1+nS))),[],'all','linear'); index_wM_best = index_wM_best-1;
[nw_best,index_M_best] = ind2sub([tmp_n_w_max,numel(index_M_unused_)],1+index_wM_best); 
nw_best = nw_best-1; index_M_best = index_M_best-1;
nM_best = index_M_unused_(1+index_M_best);
flag_M_used_(1+nM_best)=1;
euler_polar_a_(1+nM_best) = tmp_viewing_polar_a_all_(1+tmp_permutation_(1+nS));
euler_azimu_b_(1+nM_best) = tmp_viewing_azimu_b_all_(1+tmp_permutation_(1+nS));
euler_gamma_z_(1+nM_best) = 2*pi*nw_best/tmp_n_w_max;
nS = nS+1; if (nS>=tmp_n_S); nS=0; end;
nM_sum = nM_sum+1;
end;%while (sum(flag_M_used_)<tmp_n_M);
assert(nM_sum==tmp_n_M);
%%%%%%%%;
% Now go back to beginning of loop. ;
%%%%%%%%;
end;%for niteration=0:n_iteration-1;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Test innerproducts. ;
%%%%%%%%;
template_2d_k_c_0_ = zeros(n_w_sum,1);
template_2d_k_c_1_ = zeros(n_w_sum,1);
template_2d_k_p_r_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
for nw=0:n_w_(1+nk_p_r)-1;
gamma_z = 2*pi*nw/max(1,n_w_(1+nk_p_r));
template_2d_k_c_0_(1+na) = k_p_r*cos(gamma_z);
template_2d_k_c_1_(1+na) = k_p_r*sin(gamma_z);
template_2d_k_p_r_(1+na) = k_p_r;
na=na+1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
assert(na==n_w_sum);
tmp_delta_x_c_ = [+0.08;-0.12];
S_k_p_form_ = exp(+i*2*pi*(template_2d_k_c_0_*tmp_delta_x_c_(1+0) + template_2d_k_c_1_*tmp_delta_x_c_(1+1)));
I_quad = sum(S_k_p_form_.*weight_2d_k_all_)*(2*pi)^2;
I_form = h2d_(2*pi*k_p_r_max*fnorm(tmp_delta_x_c_))/(2*pi)^2 * (pi*k_p_r_max^2);
disp(sprintf(' %% testing template quadrature: plane-wave: I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));
I_quad = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,ones(n_w_sum,1),S_k_p_form_);
disp(sprintf(' %% testing template quadrature: plane-wave: I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));
S_k_p_form_ = 0.3*template_2d_k_p_r_.^3 - 0.9*template_2d_k_p_r_.^5;
I_quad = sum(S_k_p_form_.*weight_2d_k_all_)*(2*pi)^2;
I_form = 0.3*2*pi*k_p_r_max^5/5 - 0.9*2*pi*k_p_r_max^7/7;
disp(sprintf(' %% testing template quadrature: polynomial: I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));
I_quad = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,ones(n_w_sum,1),S_k_p_form_);
disp(sprintf(' %% testing template quadrature: polynomial: I_form vs I_quad: %0.16f',fnorm(I_form-I_quad)/fnorm(I_form)));
%%%%%%%%;
tmp_M_k_p_ = S_k_p__(:,max(1,min(n_viewing_all,round(n_viewing_all*0.25))));
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_S_k_p_ = S_k_p__(:,max(1,min(n_viewing_all,round(n_viewing_all*0.75))));
tmp_S_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_);
X0_= zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
tmp_T_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_,gamma_z); %<-- T = R(+gamma)*S = S(R(-gamma)*k) ;
X0_(1+nw) = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_T_k_p_,tmp_M_k_p_); %<-- <R(+gamma)*S,M> = <S,R(-gamma)*M> ;
end;%for nw=0:n_w_max-1;
X1_ = ifft(innerproduct_q_k_stretch_quad_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_S_k_q_,tmp_M_k_q_))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
%%%%%%%%;
disp(sprintf(' %% fnorm(X0_-X1_)/fnorm(X0_) = %0.16f',fnorm(X0_-X1_)/fnorm(X0_)));
flag_plot=0;
if flag_plot;
tmp_T_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_,+pi/6);
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(tmp_S_k_p_),[],colormap_beach()); 
axisnotick; axis image; title('real S');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(tmp_T_k_p_),[],colormap_beach()); 
axisnotick; axis image; title('real T=R(+pi/6)*S=S(R(-pi/6)*k)');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),imag(tmp_S_k_p_),[],colormap_beach()); 
axisnotick; axis image; title('imag S');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),imag(tmp_T_k_p_),[],colormap_beach()); 
axisnotick; axis image; title('imag T=R(+pi/6)*S=S(R(-pi/6)*k)');
figbig;
end;%if flag_plot;
end;%if flag_check;

%%%%%%%%;
% For images we will use the same n_w_ as that used by templates. ;
% We will also use the same quadrature weights for integration in 2d. ;
% In addition to the adaptively defined polar grid, ;
% we also generate M_uni_k_p__ on a uniform grid. ;
%%%%%%%%;

%%%%%%%%;
% extract ctf function. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_CTF_k_p__.mat',dir_trunk);
if (~exist(fname_mat,'file'));
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
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_k_p__(1+tmp_ij_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_ctf-1;
CTF_avg_k_p_ = mean(CTF_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
CTF_avg_k_p_r_(1+nk_p_r) = mean(CTF_avg_k_p_(1+tmp_ij_));
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
CTF_uni_k_p_r__ = zeros(n_k_p_r,n_ctf);
for nctf=0:n_ctf-1;
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_max*nk_p_r + (0:n_w_max-1);
CTF_uni_k_p_r__(1+nk_p_r,1+nctf) = mean(CTF_uni_k_p__(1+tmp_ij_,1+nctf));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nctf=0:n_ctf-1;
CTF_uni_avg_k_p_ = mean(CTF_uni_k_p__,2);
%imagesc_p(n_k_p_r,k_p_r_,n_w_,sum(n_w_),real(CTF_uni_avg_k_p_(:)),[-1,+1],colormap_beach());
CTF_uni_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_max*nk_p_r + (0:n_w_max-1);
CTF_uni_avg_k_p_r_(1+nk_p_r) = mean(CTF_uni_avg_k_p_(1+tmp_ij_));
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% Now determine the CTF cross correlation. ;
% This depends  on CTF_idx_. ;
%%%%%%%%;
n_image_sub = 1024; disp(sprintf(' %% Warning! setting n_image_sub = 1024')); 
tmp_CTF_avg_k_p_ = mean(CTF_k_p__(:,1+CTF_idx_(1+(0:n_image_sub-1))),2);
tmp_CTF_avg_k_p_r_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_CTF_avg_k_p_r_(1+nk_p_r) = mean(tmp_CTF_avg_k_p_(1+tmp_ij_));
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
     ,'CTF_avg_k_p_' ...
     ,'CTF_avg_k_p_r_' ...
     ,'CTF_uni_k_p__' ...
     ,'CTF_uni_k_p_r__' ...
     ,'CTF_uni_avg_k_p_' ...
     ,'CTF_uni_avg_k_p_r_' ...
     ,'CTF_k_p_r_xavg__' ...
     ,'CTF_k_p_r_xcor__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_CTF_k_p__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_CTF_k_p_r_xxxx__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
% Now look at some of the experimental images associated with trpv1. ;
% Note that these will *NOT* be corrected/centered according to the (presumed) displacements. ;
% Moreover, to save space, we will not store the M_x_c___. ;
%%%%%%%%;
fname_M_x_c_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_M_x_c___.mat',dir_trunk);
if (~exist(fname_M_x_c_mat,'file'));
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
n_image_sub = size(M_x_c___,3);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_euler_angle_marina_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_delta_read_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_M_x_c___center_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_M_x_c___sample',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_M_k_p___.mat',dir_trunk);
if (~exist(fname_mat,'file'));
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
%str_cost = sprintf('2d_xcor_d%.4d',floor(1000*delta_sigma));
str_cost = sprintf('2d_emp_N%d',n_image_sub);
%%%%%%%%;
if ( strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s',str_cost);
end;%if ( strcmp(str_CTF_factor,'original'));
if (~strcmp(str_CTF_factor,'original'));
str_combine = sprintf('%s_%s',str_cost,str_CTF_factor);
end;%if (~strcmp(str_CTF_factor,'original'));

verbose=1;
%%%%%%%%;
% Now test out principled marching. ;
% First set up cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_X_%s.mat',dir_trunk,str_cost);
if (~exist(fname_mat,'file'));
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
a_UX_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) = a_UX_Y_quad__(1:tmp_n_lm,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_ij_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','str_cost','str_combine'...
     ,'X_d0__','X_d0_weight_r_' ...
     ,'X__','X_weight_r_' ...
     ,'n_UX_rank'...
     ,'UX__','SX__','VX__','a_UX_Y_quad__' ...
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
[tmp_X,tmp_X_ori,tmp_X_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_UX_Y_quad__(:,1+nUX_rank),a_UX_Y_quad__(:,1+nUX_rank));
tmp_Z = transpose(UX__(:,1+nUX_rank))*X__*(UX__(:,1+nUX_rank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %+0.16f tmp_X %+0.16f tmp_X_ori*tmp_weight_so3 %+0.16f tmp_X_tau %+0.16f ratio %+0.16f',nUX_rank,n_UX_rank,tmp_Z,tmp_X,tmp_X_ori*tmp_weight_so3,tmp_X_tau,(tmp_X_ori*tmp_weight_so3)/tmp_X_tau));
end;%for nUX_rank=0:n_UX_rank-1;
end;%if flag_check;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_X_%s_A',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_X_%s_B',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
quad_lim_ = 0.5 * abs(a_UX_Y_quad__(1,1)) * [-1,+1];
for nplot=0:n_plot-1;
nk_p_r = plot_nk_p_r_(1+nplot);
[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_UX_Y_quad__(:,1+nk_p_r)),k_u_res,2*k_u_res);
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_UX_%s_S_CTF_k_p___.mat',dir_trunk,str_cost);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
UX_S_CTF_k_p_wSn___ = zeros(n_w_max,n_S,n_UX_rank);
UX_S_CTF_k_q_wSn___ = zeros(n_w_max,n_S,n_UX_rank);
for nS=0:n_S-1;
if (mod(nS,100)==0); disp(sprintf(' %% nS %d/%d',nS,n_S)); end;
tmp_S_CTF_k_p_ = S_k_p__(:,1+nS).*CTF_avg_k_p_; %<-- use average templates here under assumption that templates are used alone. ;
tmp_S_CTF_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_S_CTF_k_p_);
tmp_S_CTF_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_S_CTF_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
tmp_S_CTF_k_p__ = reshape(tmp_S_CTF_k_p_,n_w_max,n_k_p_r);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_S_CTF_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_UX_S_CTF_k_p_ = tmp_UX_S_CTF_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_S_CTF_k_p__(:,1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
UX_S_CTF_k_p_wSn___(:,1+nS,1+nUX_rank) = tmp_UX_S_CTF_k_p_;
UX_S_CTF_k_q_wSn___(:,1+nS,1+nUX_rank) = interp_p_to_q(1,n_w_max,n_w_max,tmp_UX_S_CTF_k_p_);
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nS=0:n_S-1;
save(fname_mat ...
     ,'UX_S_CTF_k_p_wSn___' ...
     ,'UX_S_CTF_k_q_wSn___' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_S_CTF_k_p___',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=1-1; clim_ = 1.5*std(real(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank)),1,'all')*[-1,+1];
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(real(squeeze(UX_S_CTF_k_p_wSn___(:,:,1+nUX_rank))),clim_); xlabel('nS'); ylabel('gamma'); title(sprintf('real(rank==%d) [%0.2f,%0.2f]',1+nUX_rank,clim_));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_S_CTF_k_q___',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
nUX_rank=min(n_UX_rank-1, 1-1); subplot(2,3,1); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 2-1); subplot(2,3,2); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 4-1); subplot(2,3,3); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1, 8-1); subplot(2,3,4); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,16-1); subplot(2,3,5); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank));
nUX_rank=min(n_UX_rank-1,32-1); subplot(2,3,6); imagesc(log10(abs(squeeze(UX_S_CTF_k_q_wSn___(:,:,1+nUX_rank)))),[-15,-5]); xlabel('nS'); ylabel('q'); title(sprintf('log10(abs(rank %d)) [-15,-5]',1+nUX_rank))
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_S_CTF_k_q___spectrum',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
n_plot=6;
for nplot=0:n_plot-1;
tmp_ = log10(abs(svds(reshape(permute(UX_S_CTF_k_q_wSn___(:,:,1+(0:nplot)),[1,3,2]),[n_w_max*(1+nplot),n_S]),n_w_max)));
subplot(2,3,1+nplot);
plot(tmp_,'o'); xlim([1,n_w_max]); ylim([-7+tmp_(1),tmp_(1)+1]);
xlabel('pc'); ylabel('log10(sigma)');
title(sprintf('spectrum of UX_S_CTF_k_q_wSn___(:,:,nUX_rank==1:%d)',1+nplot),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_S_CTF_k_p___A',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1);
for nUX_rank=0:min(n_UX_rank-1,15-1);
subplot(3,5,1+nUX_rank); hold on;
for nS=0:32:n_S-1;
nc = max(0,min(n_c-1,floor(n_c*nS/n_S)));
plot(2*pi*(0:n_w_max-1)/n_w_max,real(UX_S_CTF_k_p_wSn___(:,1+nS,1+nUX_rank)),'.','Color',c_(1+nc,:)); 
end;%for nS=0:n_S-1;
xlim([0,2*pi]);
ylim(4e-8*[-1,+1]); 
xlabel('gamma_z','Interpreter','none');
ylabel('real(UX_S_CTF_k_p)','Interpreter','none');
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_UX_%s_M_k_p___.mat',dir_trunk,str_combine);
if (~exist(fname_mat,'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_M_k_p___',dir_trunk,str_combine);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_M_k_q___',dir_trunk,str_combine);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_UX_%s_M_k_q___spectrum',dir_trunk,str_combine);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_UX_%s_M_k_p___.mat',dir_trunk,str_combine);
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

fname_mat = sprintf('%s_mat/test_principled_marching_trpv1_17_c_k_Y_.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_n_order = 5; tmp_n_M = n_image_sub;
tmp_euler_polar_a_ = +euler_angle_marina_(1,1+(0:n_image_sub-1));
tmp_euler_azimu_b_ = +euler_angle_marina_(2,1+(0:n_image_sub-1));
tmp_euler_gamma_z_ = -euler_angle_marina_(3,1+(0:n_image_sub-1));
tmp_image_delta_x_ = +1.0*delta_read_x_(1+(0:n_image_sub-1));
tmp_image_delta_y_ = +1.0*delta_read_y_(1+(0:n_image_sub-1));
%%%%%%%%;
tmp_t = tic;
c_k_Y_reco_ = cg_lsq_4(tmp_n_order,n_k_p_r,k_p_r_,l_max_,n_w_,tmp_n_M,M_k_p__,CTF_idx_,CTF_k_p__,tmp_euler_polar_a_,tmp_euler_azimu_b_,tmp_euler_gamma_z_,tmp_image_delta_x_,tmp_image_delta_y_);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_c_x_u_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
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
% Set up a trial run to test alternating minimization using experimental principled-images. ;
%%%%%%%%;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_iteration = 32;
dat_rseed = 0;
rng(dat_rseed);
nUX_rank = 4;
tmp_n_k_p_r = 1+nUX_rank;
tmp_k_p_r_ = ones(1+nUX_rank,1); tmp_k_p_r_max = 1;
tmp_n_w_ = n_w_max*ones(1+nUX_rank,1); tmp_n_w_max = max(tmp_n_w_); tmp_n_w_sum = sum(tmp_n_w_); 
tmp_weight_3d_k_p_r_ = ones(1+nUX_rank,1); tmp_weight_2d_k_p_r_ = ones(1+nUX_rank,1);
tmp_l_max_ = l_max_max*ones(1+nUX_rank,1);
tmp_n_lm_ = (1+tmp_l_max_).^2; tmp_n_lm_sum = sum(tmp_n_lm_);
tmp_a_k_Y_quad_ = reshape(a_UX_Y_quad__(:,1+(0:nUX_rank)),[tmp_n_lm_sum,1]); %<-- use as ground truth for this set of principled-image-rings. ;
tmp_M_k_p__ = reshape(permute(UX_M_k_p___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
tmp_M_k_q__ = reshape(permute(UX_M_k_q___(:,1+(0:dat_n_M-1),1+(0:nUX_rank)),[1,3,2]),[tmp_n_w_sum,dat_n_M]);
%%%%%%%%;
tmp_rseed=dat_rseed;tmp_n_iteration=dat_n_iteration;tmp_n_iteration_register=1;tmp_viewing_k_eq_d = 1/(2*pi);tmp_n_order = dat_n_order; 
flag_MS_vs_SM = 1; f_rand = [];
[dat_X_best_allshell_MS_,~,euler_polar_a_MS__,euler_azimu_b_MS__,euler_gamma_z_MS__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r_,tmp_k_p_r_max,tmp_weight_3d_k_p_r_,tmp_weight_2d_k_p_r_,tmp_n_w_,dat_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max_,tmp_a_k_Y_quad_,[],[],[],flag_MS_vs_SM,f_rand);
%%%%%%%%;
tmp_n_order = 5; tmp_n_M = n_image_sub;
tmp_euler_polar_a_ = euler_polar_a_MS__(:,end);
tmp_euler_azimu_b_ = euler_azimu_b_MS__(:,end);
tmp_euler_gamma_z_ = euler_gamma_z_MS__(:,end);
tmp_t = tic;
tmp_d_k_Y_reco_ = cg_lsq_3(tmp_n_order,n_k_p_r,l_max_,n_w_,tmp_n_M,M_k_p__,CTF_idx_,CTF_k_p__,tmp_euler_polar_a_,tmp_euler_azimu_b_,tmp_euler_gamma_z_);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> tmp_d_k_Y_reco_ time %0.2fs',tmp_t));
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,tmp_d_k_Y_reco_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,tmp_d_k_Y_reco_));
X_best_MS = max(tmp_X_best_orig,tmp_X_best_flip);
%%%%%%%%;
flag_MS_vs_SM = 0; f_rand = 0.05;
[dat_X_best_allshell_SM_,~,euler_polar_a_SM__,euler_azimu_b_SM__,euler_gamma_z_SM__] = am_2(tmp_rseed,tmp_n_iteration,tmp_n_iteration_register,tmp_viewing_k_eq_d,tmp_n_order,tmp_n_k_p_r,tmp_k_p_r_,tmp_k_p_r_max,tmp_weight_3d_k_p_r_,tmp_weight_2d_k_p_r_,tmp_n_w_,dat_n_M,tmp_M_k_p__,tmp_M_k_q__,1,ones(tmp_n_w_sum,1),tmp_l_max_,tmp_a_k_Y_quad_,euler_polar_a_MS__(:,end),euler_azimu_b_MS__(:,end),euler_gamma_z_MS__(:,end),flag_MS_vs_SM,f_rand);
%%%%%%%%;
tmp_n_order = 5; tmp_n_M = n_image_sub;
tmp_euler_polar_a_ = euler_polar_a_SM__(:,end);
tmp_euler_azimu_b_ = euler_azimu_b_SM__(:,end);
tmp_euler_gamma_z_ = euler_gamma_z_SM__(:,end);
tmp_t = tic;
tmp_d_k_Y_reco_ = cg_lsq_3(tmp_n_order,n_k_p_r,l_max_,n_w_,tmp_n_M,M_k_p__,CTF_idx_,CTF_k_p__,tmp_euler_polar_a_,tmp_euler_azimu_b_,tmp_euler_gamma_z_);
tmp_t = toc(tmp_t); disp(sprintf(' %% M_k_p__ --> tmp_d_k_Y_reco_ time %0.2fs',tmp_t));
[tmp_X_best_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,tmp_d_k_Y_reco_);
[tmp_X_best_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,tmp_d_k_Y_reco_));
X_best_SM = max(tmp_X_best_orig,tmp_X_best_flip);
%%%%%%%%;
disp(sprintf(' %% dat_rseed %d nUX_rank %d: first phase X_best_MS %0.3f --> %0.3f ; second phase X_best_SM %0.3f --> %0.3f ',dat_rseed,nUX_rank,dat_X_best_allshell_MS_(end),X_best_MS,dat_X_best_allshell_SM_(end),X_best_SM));
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Test out a single translation using one of the templates. ;
%%%%%%%%;
tmp_ns = 0; %<-- pick one of the templates. ;
tmp_S_k_p_ = S_k_p__(:,1+tmp_ns);
tmp_S_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_);
tmp_S_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_S_k_p_.*weight_2d_k_all_*(2*pi)^2);
tmp_N_pixel = 1.5; tmp_delta_max = tmp_N_pixel/(2*pi*k_p_r_max)*pi*sqrt(2);
tmp_delta_ = tmp_delta_max*0.95*[cos(pi/3);sin(pi/3)];
tmp_T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_S_k_p_,tmp_delta_(1+0),tmp_delta_(1+1));
tmp_T_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_T_k_p_);
tmp_T_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_T_k_p_.*weight_2d_k_all_*(2*pi)^2);
%%%%%%%%;
%dir_svd = '/data/rangan/dir_cryoem/dir_rangan_playground/dir_gen_Jsvd_6';
%FTK = get_svd_FTK_2(1e-3,k_p_r_,n_k_p_r,1,tmp_delta_max,0,1,dir_svd);
FTK = gen_Jsvd_FTK_7(k_p_r_max,tmp_N_pixel,1e-3);
FTK.n_delta_v = 1;
FTK.delta_x_ = tmp_delta_(1+0);
FTK.delta_y_ = tmp_delta_(1+1);
FTK.svd_d_max = tmp_delta_max;
FTK.svd_r_max = 2*pi*k_p_r_max;
%%%%%%%%;
tmp_R_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_S_k_q_,tmp_delta_(1+0),tmp_delta_(1+1));
tmp_R_k_p_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,tmp_R_k_q_);
tmp_R_x_c_ = interp_k_p_to_x_c_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,tmp_R_k_p_.*weight_2d_k_all_*(2*pi)^2);
flag_plot=0;
if flag_plot;
clf;
subplot(2,3,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_S_k_p_),[],colormap_beach()); title('S_k_p_','Interpreter','none'); axisnotick; axis image;
subplot(2,3,4); imagesc(real(tmp_S_x_c_)); title('S_x_c_','Interpreter','none'); axisnotick; axis image; set(gca,'Ydir','normal');
subplot(2,3,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_T_k_p_),[],colormap_beach()); title('T_k_p_','Interpreter','none'); axisnotick; axis image;
subplot(2,3,5); imagesc(real(tmp_T_x_c_)); title('T_x_c_','Interpreter','none'); axisnotick; axis image; set(gca,'Ydir','normal');
subplot(2,3,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_R_k_p_),[],colormap_beach()); title('R_k_p_','Interpreter','none'); axisnotick; axis image;
subplot(2,3,6); imagesc(real(tmp_R_x_c_)); title('R_x_c_','Interpreter','none'); axisnotick; axis image; set(gca,'Ydir','normal');
figbig;
end;%if flag_plot;
disp(sprintf(' %% tmp_T_k_q_ vs tmp_R_k_q_: %0.16f',fnorm(tmp_T_k_q_-tmp_R_k_q_)/fnorm(tmp_T_k_q_)));
%%%%%%%%;
% Now test out multiple translations across one image and one template. ;
%%%%%%%%;
tmp_N_pixel = 1.5; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,tmp_N_pixel,1e-3);
tmp_delta_max = tmp_N_pixel/(2*pi*k_p_r_max)*pi*sqrt(2);
tmp_n_delta_v = 8;
tmp_delta__ = tmp_delta_max/sqrt(2)*(2*rand(2,tmp_n_delta_v)-1);
FTK.n_delta_v = tmp_n_delta_v;
FTK.delta_x_ = transpose(tmp_delta__(1+0,:));
FTK.delta_y_ = transpose(tmp_delta__(1+1,:));
FTK.svd_d_max = tmp_delta_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
%%%%%%%%;
tmp_ns = 0; %<-- pick one of the templates. ;
tmp_S_k_p_ = S_k_p__(:,1+tmp_ns);
tmp_S_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_);
tmp_nm = 0; %<-- pick one of the images. ;
tmp_M_k_p_ = M_k_p__(:,1+tmp_nm);
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_X0__ = zeros(tmp_n_delta_v,n_w_max); %<-- direct calculation of tmp_X0__(1+ndelta_v,1+nw) = <R(+gamma)*S,T(+delta)M> ;
tmp_X1__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations, brute-force for translations. ;
tmp_X2__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations, direct svd for translations. ;
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations and svd for translations, along with blas. ;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_T_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_T_k_p_);
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
tmp_R_k_p_ = rotate_p2p_fx(n_k_p_r,n_w_,n_w_sum,tmp_S_k_p_,gamma_z); %<-- R = R(+gamma)*S = S(R(-gamma)*k) ;
tmp_X0__(1+ndelta_v,1+nw) = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_R_k_p_,tmp_T_k_p_); %<-- <R(+gamma)*S,M> = <S,R(-gamma)*M> ;
end;%for nw=0:n_w_max-1;
tmp_X1__(1+ndelta_v,:) = ifft(innerproduct_q_k_stretch_quad_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_S_k_q_,tmp_T_k_q_))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_X2__(1+ndelta_v,:) = ifft(innerproduct_q_k_stretch_quad_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_S_k_q_,tmp_T_k_q_))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
tmp_delta_r = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_S_k_q_rw__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,weight_2d_k_p_r_/(2*pi),n_w_,1,tmp_S_k_q_);
%%%%%%%%;
tmp_VSM__ = zeros(FTK.n_svd_l,n_w_max);
for nl=0:FTK.n_svd_l-1;
tmp_VSM__(1+nl,:) = ifft(tmp_V_r__(1+nl,:)*(conj(tmp_S_k_q_rw__(:,:)).*tmp_M_k_q_rwl___(:,:,1+tmp_l_max+FTK.svd_l_(1+nl))))*n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
tmp_USEVSM__ = (tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*tmp_VSM__;
%%%%%%%%;
tmp_X3__ = tmp_USEVSM__;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
colormap(colormap_beach());
subplot(2,2,1+0); imagesc(real(tmp_X0__)); axisnotick; title('X0');
subplot(2,2,1+1); imagesc(real(tmp_X1__)); axisnotick; title('X1');
subplot(2,2,1+2); imagesc(real(tmp_X2__)); axisnotick; title('X2');
subplot(2,2,1+3); imagesc(real(tmp_X3__)); axisnotick; title('X3');
figbig;
end;%if flag_plot;
%%%%%%%%;
% Now test out multiple translations across one principal-image and one principal-template. ;
%%%%%%%%;
tmp_N_pixel = 1.5; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,tmp_N_pixel,1e-3);
tmp_delta_max = tmp_N_pixel/(2*pi*k_p_r_max)*pi*sqrt(2);
tmp_n_delta_v = 8;
tmp_delta__ = tmp_delta_max/sqrt(2)*(2*rand(2,tmp_n_delta_v)-1);
FTK.n_delta_v = tmp_n_delta_v;
FTK.delta_x_ = transpose(tmp_delta__(1+0,:));
FTK.delta_y_ = transpose(tmp_delta__(1+1,:));
FTK.svd_d_max = tmp_delta_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
%%%%%%%%;
tmp_ns = 0; %<-- pick one of the principal-templates. ;
tmp_n_UX_rank = 1+1; %<-- pick maximum principal-mode. ;
tmp_n_k_p_r = tmp_n_UX_rank;
tmp_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_n_w_ = n_w_max*ones(tmp_n_k_p_r,1);
tmp_n_w_sum = sum(tmp_n_w_);
tmp_n_w_max = max(tmp_n_w_);
tmp_weight_2d_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_UX_S_CTF_k_p__ = squeeze(UX_S_CTF_k_p_wSn___(:,1+tmp_ns,1+(0:tmp_n_UX_rank-1)));
tmp_UX_S_CTF_k_q__ = squeeze(UX_S_CTF_k_q_wSn___(:,1+tmp_ns,1+(0:tmp_n_UX_rank-1)));
flag_plot=0;
if flag_plot; clf; hold on; plot(1:n_w_max,real(tmp_UX_S_CTF_k_p__),'LineWidth',3); hold off; end;
tmp_nm = 0; %<-- pick one of the images. ;
tmp_M_k_p_ = M_k_p__(:,1+tmp_nm); %<-- subplot(1,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(tmp_M_k_p_),[],colormap_beach()); subplot(1,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(tmp_M_k_p_),[],colormap_beach());
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_M_k_q__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1,tmp_M_k_q_);
tmp_UX_M_k_q__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
tmp_UX_M_k_p__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+tmp_nm);
if strfind(str_CTF_factor,'rescale'); CTF_factor = CTF_avg_k_p_r_(1+nk_p_r) / max(1e-6,CTF_k_p_r__(1+nk_p_r,1+nctf)); end;
if strfind(str_CTF_factor,'original'); CTF_factor = 1; end;
tmp_UX_M_k_q__(:,1+nUX_rank) = tmp_UX_M_k_q__(:,1+nUX_rank) + transpose(UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*CTF_factor*tmp_M_k_q__(1+nk_p_r,:));
end;%for nk_p_r=0:n_k_p_r-1;
tmp_UX_M_k_p__(:,1+nUX_rank) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_M_k_q__(:,1+nUX_rank));
end;%for nUX_rank=0:tmp_n_UX_rank-1;
if flag_plot; hold on; plot(1:n_w_max,real(tmp_UX_M_k_p__),'LineWidth',2); hold off; end;
tmp_X0__ = zeros(tmp_n_delta_v,n_w_max); %<-- direct calculation of tmp_X0__(1+ndelta_v,1+nw) = <R(+gamma)*S,T(+delta)M> ;
tmp_X1__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations, brute-force for translations. ;
tmp_X2__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations, direct svd for translations. ;
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max); %<-- use bessel-coordinates for rotations and svd for translations, along with blas. ;
for ndelta_v=0:tmp_n_delta_v-1;
tmp_T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_p_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_T_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_T_k_p_);
tmp_T_k_q__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1,tmp_T_k_q_);
tmp_UX_T_k_q__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
tmp_UX_T_k_p__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
tmp_UX_T_k_q__(:,1+nUX_rank)=0;
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+tmp_nm);
if strfind(str_CTF_factor,'rescale'); CTF_factor = CTF_avg_k_p_r_(1+nk_p_r) / max(1e-6,CTF_k_p_r__(1+nk_p_r,1+nctf)); end;
if strfind(str_CTF_factor,'original'); CTF_factor = 1; end;
tmp_UX_T_k_q__(:,1+nUX_rank) = tmp_UX_T_k_q__(:,1+nUX_rank) + transpose(UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*CTF_factor*tmp_T_k_q__(1+nk_p_r,:));
end;%for nk_p_r=0:n_k_p_r-1;
tmp_UX_T_k_p__(:,1+nUX_rank) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_T_k_q__(:,1+nUX_rank));
end;%for nUX_rank=0:tmp_n_UX_rank-1;
if flag_plot; hold on; plot(1:n_w_max,real(tmp_UX_T_k_p__),'--','LineWidth',1); hold off; end;
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
tmp_UX_R_CTF_k_p__ = reshape(rotate_p2p_fx(tmp_n_k_p_r,tmp_n_w_,tmp_n_w_sum,tmp_UX_S_CTF_k_p__(:),gamma_z),[tmp_n_w_max,tmp_n_UX_rank]) ; %<-- R = R(+gamma)*S = S(R(-gamma)*k),
tmp_X0__(1+ndelta_v,1+nw) = innerproduct_p_quad(tmp_n_k_p_r,tmp_k_p_r_,tmp_weight_2d_k_p_r_/(2*pi),tmp_n_w_,tmp_n_w_sum,tmp_UX_R_CTF_k_p__(:),tmp_UX_T_k_p__(:)); %<-- <R(+gamma)*S,M> = <S,R(-gamma)*M> ;
end;%for nw=0:n_w_max-1;
tmp_X1__(1+ndelta_v,:) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r_,tmp_weight_2d_k_p_r_/(2*pi),tmp_n_w_,tmp_n_w_sum,tmp_UX_S_CTF_k_q__(:),tmp_UX_T_k_q__(:)))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_T_k_q__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1,tmp_T_k_q_);
tmp_UX_T_k_q__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
tmp_UX_T_k_p__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
tmp_UX_T_k_q__(:,1+nUX_rank)=0;
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+tmp_nm);
if strfind(str_CTF_factor,'rescale'); CTF_factor = CTF_avg_k_p_r_(1+nk_p_r) / max(1e-6,CTF_k_p_r__(1+nk_p_r,1+nctf)); end;
if strfind(str_CTF_factor,'original'); CTF_factor = 1; end;
tmp_UX_T_k_q__(:,1+nUX_rank) = tmp_UX_T_k_q__(:,1+nUX_rank) + transpose(UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*CTF_factor*tmp_T_k_q__(1+nk_p_r,:));
end;%for nk_p_r=0:n_k_p_r-1;
tmp_UX_T_k_p__(:,1+nUX_rank) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_T_k_q__(:,1+nUX_rank));
end;%for nUX_rank=0:tmp_n_UX_rank-1;
tmp_X2__(1+ndelta_v,:) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r_,tmp_weight_2d_k_p_r_/(2*pi),tmp_n_w_,tmp_n_w_sum,tmp_UX_S_CTF_k_q__(:),tmp_UX_T_k_q__(:)))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
end;%for ndelta_v=0:tmp_n_delta_v-1;
disp(sprintf(' %% tmp_X0__ vs tmp_X1__: %0.16f',fnorm(tmp_X0__-tmp_X1__)/fnorm(tmp_X0__)));
disp(sprintf(' %% tmp_X0__ vs tmp_X2__: %0.16f',fnorm(tmp_X0__-tmp_X2__)/fnorm(tmp_X0__)));
%%%%%%%%;
tmp_X3__ = zeros(tmp_n_delta_v,n_w_max);
%%%%%%%%;
tmp_U_d__ = zeros(tmp_n_delta_v,FTK.n_svd_l);
tmp_delta_r_ = zeros(tmp_n_delta_v,1);
for ndelta_v=0:tmp_n_delta_v-1;
tmp_U_d_ = zeros(FTK.n_svd_l,1);
tmp_delta_r_(1+ndelta_v) = fnorm(tmp_delta__(:,1+ndelta_v));
tmp_delta_r = fnorm(tmp_delta__(:,1+ndelta_v));
if (tmp_delta_r>FTK.svd_d_max);
disp(sprintf(' %% Warning, tmp_delta_r %0.6f > svd_d_max %0.6f',tmp_delta_r,FTK.svd_d_max));
end;%if;
tmp_svd_d = (tmp_delta_r_(1+ndelta_v) - FTK.svd_d_m)/FTK.svd_d_c;
for nl=0:FTK.n_svd_l-1;
tmp_U_d__(1+ndelta_v,1+nl) = polyval_r8_reverse_0(FTK.n_svd_d,FTK.svd_U_d_(1+0+nl*FTK.n_svd_d+(0:FTK.n_svd_d-1)),1,tmp_svd_d);
end;%for nl=0:FTK.n_svd_l-1;
end;%for ndelta_v=0:tmp_n_delta_v-1;
%%%%%%%%;
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_V_r__ = zeros(FTK.n_svd_l,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
if (2*pi*k_p_r_(1+nk_p_r)>FTK.svd_r_max);
disp(sprintf(' %% Warning, 2*pi*k_p_r_(1+nk_p_r) %0.6f > FTK.svd_r_max %0.6f',2*pi*k_p_r_(1+nk_p_r),FTK.svd_r_max));
end;%if;
tmp_svd_r = (2*pi*k_p_r_(1+nk_p_r) - FTK.svd_r_m)/FTK.svd_r_c;
for nl=0:FTK.n_svd_l-1;
tmp_V_r__(1+nl,1+nk_p_r) = polyval_r8_reverse_0(FTK.n_svd_r,FTK.svd_V_r_(1+0+nl*FTK.n_svd_r+(0:FTK.n_svd_r-1)),1,tmp_svd_r);
end;%for nl=0:FTK.n_svd_l-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% Note that, at this point these two should match: ;
% hold on; for ndelta_v=0:tmp_n_delta_v-1; plot(2*pi*k_p_r_,real(tmp_U_d__(1+ndelta_v,:).*tmp_expiw__(1+ndelta_v,:))*diag(FTK.svd_s_)*tmp_V_r__(:,:),'ro',2*pi*k_p_r_,real(exp(-i*2*pi*k_p_r_*tmp_delta_r_(1+ndelta_v)*cos(tmp_delta_w_(1+ndelta_v)))),'r-',2*pi*k_p_r_,imag(tmp_U_d__(1+ndelta_v,:).*tmp_expiw__(1+ndelta_v,:))*diag(FTK.svd_s_)*tmp_V_r__(:,:),'bo',2*pi*k_p_r_,imag(exp(-i*2*pi*k_p_r_*tmp_delta_r_(1+ndelta_v)*cos(tmp_delta_w_(1+ndelta_v)))),'b-'); end; hold off;
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q__ = zeros(n_w_sum,1+2*tmp_l_max);
ic=0;
for nk_p_r=0:n_k_p_r-1;
for nw=0:n_w_(1+nk_p_r)-1;
tmp_n_w_t = floor(1.0d0*n_w_(1+nk_p_r)/2.0d0);
for tmp_I_l=-tmp_l_max:+tmp_l_max;
nwc = nw;
if (nwc>=tmp_n_w_t); nwc = nwc - n_w_(1+nk_p_r); end;%if;
flag_ict_overflow = 0;
nwd = nwc + tmp_I_l;
if (abs(nwd)<tmp_n_w_t); nwt = periodize(nwd,0,n_w_(1+nk_p_r));  else; nwt = 0; flag_ict_overflow = 1; end;%if;
ict = ic-nw+nwt;
if (flag_ict_overflow==0); tmp_M_k_q__(1+ic,1+tmp_l_max+tmp_I_l) = tmp_M_k_q_(1+ict); end;%if;
end;%for tmp_I_l=-tmp_l_max:+tmp_l_max;
ic = ic + 1;
end;%for nw=0:n_w_(1+nk_p_r)-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% Now at this point, these two should match (original inner product). ;
%%%%%%%%;
ndelta_v=0; gamma_z = 0;
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_X2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_S_k_q_,tmp_T_k_q_);
tmp_X3 = 0;
for nl=0:FTK.n_svd_l-1;
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_w_csum_(1+nk_p_r) + (0:n_w_(1+nk_p_r)-1);
tmp_X3 = tmp_X3 + tmp_U_d__(1+ndelta_v,1+nl)*tmp_expiw__(1+ndelta_v,1+nl)*FTK.svd_s_(1+nl)*tmp_V_r__(1+nl,1+nk_p_r)*sum(conj(tmp_S_k_q_(1+tmp_ij_)).*tmp_M_k_q__(1+tmp_ij_,1+tmp_l_max+FTK.svd_l_(1+nl)))*weight_2d_k_p_r_(1+nk_p_r)/n_w_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nl=0:FTK.n_svd_l-1;
disp(sprintf(' %% tmp_X2 vs tmp_X3: %0.16f+i%0.16f %0.16f+i%0.16f',real(tmp_X2),imag(tmp_X2),real(tmp_X3),imag(tmp_X3)));
%%%%%%%%;
% These should also match (original inner product). ;
%%%%%%%%;
ndelta_v=0; gamma_z = 0;
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_X2_ = zeros(n_w_max,1); tmp_X2_(:) = ifft(innerproduct_q_k_stretch_quad_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_/(2*pi),n_w_,n_w_sum,tmp_S_k_q_,tmp_T_k_q_))*n_w_max;
tmp_X3_ = zeros(n_w_max,1);
for nl=0:FTK.n_svd_l-1;
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_ij_get_ = n_w_csum_(1+nk_p_r) + (0:tmp_n_w-1);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_X3_(1+tmp_ij_set_) = tmp_X3_(1+tmp_ij_set_) + tmp_U_d__(1+ndelta_v,1+nl)*tmp_expiw__(1+ndelta_v,1+nl)*FTK.svd_s_(1+nl)*tmp_V_r__(1+nl,1+nk_p_r)*conj(tmp_S_k_q_(1+tmp_ij_get_)).*tmp_M_k_q__(1+tmp_ij_get_,1+tmp_l_max+FTK.svd_l_(1+nl))*weight_2d_k_p_r_(1+nk_p_r)/tmp_n_w;
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_X3_(1+tmp_ij_set_) = tmp_X3_(1+tmp_ij_set_) + tmp_U_d__(1+ndelta_v,1+nl)*tmp_expiw__(1+ndelta_v,1+nl)*FTK.svd_s_(1+nl)*tmp_V_r__(1+nl,1+nk_p_r)*conj(tmp_S_k_q_(1+tmp_ij_get_)).*tmp_M_k_q__(1+tmp_ij_get_,1+tmp_l_max+FTK.svd_l_(1+nl))*weight_2d_k_p_r_(1+nk_p_r)/tmp_n_w;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_X3_ = ifft(tmp_X3_)*n_w_max;
disp(sprintf(' %% tmp_X2_ vs tmp_X3_: %0.16f',fnorm(tmp_X2_-tmp_X3_)/fnorm(tmp_X2_)));
%%%%%%%%;
% These should also match (principal inner product). ;
%%%%%%%%;
ndelta_v=0; gamma_z = 0;
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q_,tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_T_k_q__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1,tmp_T_k_q_);
tmp_UX_T_k_q__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
tmp_UX_T_k_p__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_UX_T_k_q__(:,1+nUX_rank) = tmp_UX_T_k_q__(:,1+nUX_rank) + transpose(UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_T_k_q__(1+nk_p_r,:));
end;%for nk_p_r=0:n_k_p_r-1;
tmp_UX_T_k_p__(:,1+nUX_rank) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_T_k_q__(:,1+nUX_rank));
end;%for nUX_rank=0:tmp_n_UX_rank-1;
tmp_X2_ = zeros(n_w_max,1); tmp_X2_(:) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r_,tmp_weight_2d_k_p_r_/(2*pi),tmp_n_w_,tmp_n_w_sum,tmp_UX_S_CTF_k_q__(:),tmp_UX_T_k_q__(:)))*tmp_n_w_max;
tmp_X3_ = zeros(n_w_max,1);
for nl=0:FTK.n_svd_l-1;
for nUX_rank=0:tmp_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_ij_get_ = n_w_csum_(1+nk_p_r) + (0:tmp_n_w-1);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_X3_(1+tmp_ij_set_) = tmp_X3_(1+tmp_ij_set_) + tmp_U_d__(1+ndelta_v,1+nl)*tmp_expiw__(1+ndelta_v,1+nl)*FTK.svd_s_(1+nl)*tmp_V_r__(1+nl,1+nk_p_r)*conj(tmp_UX_S_CTF_k_q__(1+tmp_ij_set_,1+nUX_rank)).*tmp_M_k_q__(1+tmp_ij_get_,1+tmp_l_max+FTK.svd_l_(1+nl))*UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r) / n_w_max;
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_X3_(1+tmp_ij_set_) = tmp_X3_(1+tmp_ij_set_) + tmp_U_d__(1+ndelta_v,1+nl)*tmp_expiw__(1+ndelta_v,1+nl)*FTK.svd_s_(1+nl)*tmp_V_r__(1+nl,1+nk_p_r)*conj(tmp_UX_S_CTF_k_q__(1+tmp_ij_set_,1+nUX_rank)).*tmp_M_k_q__(1+tmp_ij_get_,1+tmp_l_max+FTK.svd_l_(1+nl))*UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r) / n_w_max;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:tmp_n_UX_rank-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_X3_ = ifft(tmp_X3_)*n_w_max;
disp(sprintf(' %% tmp_X2_ vs tmp_X3_: %0.16f',fnorm(tmp_X2_-tmp_X3_)/fnorm(tmp_X2_)));
%%%%%%%%;
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1+2*tmp_l_max,tmp_M_k_q__);
tmp_UX_S_CTF_k_q_1wn___ = innerproduct_q_k_stretch_quad_stack__0(1,tmp_n_w_max/(2*pi),tmp_n_w_max,tmp_n_UX_rank,tmp_UX_S_CTF_k_q__);
%%%%%%%%;
tmp_VUXM___ = zeros(FTK.n_svd_l,n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
tmp_VUXM___(1+nl,:,1+nUX_rank) = tmp_V_r__(1+nl,:)*diag(UX__(:,1+nUX_rank).*X_weight_r_(:))*tmp_M_k_q_rwl___(:,:,1+tmp_l_max+FTK.svd_l_(1+nl)) / n_w_max ;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nUX_rank=0:tmp_n_UX_rank-1;
%tmp_SVUXM__ = zeros(FTK.n_svd_l,n_w_max);
%for nl=0:FTK.n_svd_l-1;
%for nUX_rank=0:tmp_n_UX_rank-1;
%tmp_SVUXM__(1+nl,:) = tmp_SVUXM__(1+nl,:) + conj(tmp_UX_S_CTF_k_q_1wn___(1+0,:,1+nUX_rank)).*tmp_VUXM___(1+nl,:,1+nUX_rank);
%end;%for nUX_rank=0:tmp_n_UX_rank-1;
%tmp_SVUXM__(1+nl,:) = ifft(tmp_SVUXM__(1+nl,:))*n_w_max;
%end;%for nl=0:FTK.n_svd_l-1;
tmp_SVUXM__ = ifft(sum(repmat(conj(reshape(tmp_UX_S_CTF_k_q__,[1,tmp_n_w_max,tmp_n_UX_rank])),[FTK.n_svd_l,1,1]).*tmp_VUXM___,3),[],2)*tmp_n_w_max;
tmp_USESVUXM__ = (tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*tmp_SVUXM__;
%%%%%%%%;
tmp_X3__ = tmp_USESVUXM__;
disp(sprintf(' %% tmp_X0__ vs tmp_X3__: %0.16f',fnorm(tmp_X0__-tmp_X3__)/fnorm(tmp_X0__)));
disp(sprintf(' %% tmp_X2__ vs tmp_X3__: %0.16f',fnorm(tmp_X2__-tmp_X3__)/fnorm(tmp_X2__)));
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
colormap(colormap_beach());
subplot(2,2,1+0); imagesc(real(tmp_X0__)); axisnotick; title('X0');
subplot(2,2,1+1); imagesc(real(tmp_X1__)); axisnotick; title('X1');
subplot(2,2,1+2); imagesc(real(tmp_X2__)); axisnotick; title('X2');
subplot(2,2,1+3); imagesc(real(tmp_X3__)); axisnotick; title('X3');
figbig;
end;%if flag_plot;
%%%%%%%%;
tmp_X4__ = zeros(tmp_n_delta_v,n_w_max); 
tmp_U_d__ = permute(reshape(FTK.svd_polyval_U_d_,[FTK.n_svd_l,tmp_n_delta_v]),[2,1]);
tmp_V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q_rwl___ = innerproduct_q_k_stretch_quad_stack___0(n_k_p_r,n_w_,tmp_M_k_q_,tmp_l_max);
tmp_VUXM___ = zeros(FTK.n_svd_l,n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
tmp_VUXM___(1+nl,:,1+nUX_rank) = tmp_V_r__(1+nl,:)*diag(UX__(:,1+nUX_rank).*X_weight_r_(:))*tmp_M_k_q_rwl___(:,:,1+tmp_l_max+FTK.svd_l_(1+nl)) / n_w_max ;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nUX_rank=0:tmp_n_UX_rank-1;
tmp_SVUXM__ = ifft(sum(repmat(conj(reshape(tmp_UX_S_CTF_k_q__,[1,tmp_n_w_max,tmp_n_UX_rank])),[FTK.n_svd_l,1,1]).*tmp_VUXM___,3),[],2)*tmp_n_w_max;
tmp_USESVUXM__ = (tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*tmp_SVUXM__;
tmp_X4__ = tmp_USESVUXM__;
disp(sprintf(' %% tmp_X0__ vs tmp_X4__: %0.16f',fnorm(tmp_X0__-tmp_X4__)/fnorm(tmp_X0__)));
disp(sprintf(' %% tmp_X2__ vs tmp_X4__: %0.16f',fnorm(tmp_X2__-tmp_X4__)/fnorm(tmp_X2__)));
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now test out timing for multiple translations across multiple principal-images and multiple principal-templates. ;
%%%%%%%%;
tmp_N_pixel = 1.00; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,tmp_N_pixel,1e-3);
tmp_delta_max = tmp_N_pixel/(2*pi*k_p_r_max)*pi*sqrt(2);
tmp_n_delta_v = 64;
tmp_delta__ = tmp_delta_max/sqrt(2)*(2*rand(2,tmp_n_delta_v)-1);
FTK.n_delta_v = tmp_n_delta_v;
FTK.delta_x_ = transpose(tmp_delta__(1+0,:));
FTK.delta_y_ = transpose(tmp_delta__(1+1,:));
FTK.svd_d_max = tmp_delta_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
%%%%%%%%;
tmp_n_S = 63; tmp_nS_ = floor(n_viewing_all*rand(tmp_n_S,1));
tmp_n_UX_rank = 1+3; %<-- pick maximum principal-mode. ;
tmp_n_k_p_r = tmp_n_UX_rank;
tmp_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_n_w_ = n_w_max*ones(tmp_n_k_p_r,1);
tmp_n_w_sum = sum(tmp_n_w_);
tmp_n_w_max = max(tmp_n_w_);
tmp_weight_2d_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_UX_S_CTF_k_p_wnS___ = permute(squeeze(UX_S_CTF_k_p_wSn___(:,1+tmp_nS_,1+(0:tmp_n_UX_rank-1))),[1,3,2]);
tmp_UX_S_CTF_k_q_wnS___ = permute(squeeze(UX_S_CTF_k_q_wSn___(:,1+tmp_nS_,1+(0:tmp_n_UX_rank-1))),[1,3,2]);
%%%%%%%%;
tmp_n_M = 65; tmp_nM_ = floor(n_image_sub*rand(tmp_n_M,1));
tmp_M_k_p__ = M_k_p__(:,1+tmp_nM_);
tmp_M_k_q__ = zeros(n_w_sum,tmp_n_M); for nM=0:tmp_n_M-1; tmp_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p__(:,1+nM)); end;
%%%%%%%%;
tmp_X3____ = zeros(tmp_n_delta_v,n_w_max,tmp_n_S,tmp_n_M);
%%%%%%%%;
tmp_U_d__ = permute(reshape(FTK.svd_polyval_U_d_,[FTK.n_svd_l,tmp_n_delta_v]),[2,1]);
tmp_V_r__ = reshape(FTK.svd_polyval_V_r_,[FTK.n_svd_l,n_k_p_r]);
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
%%%%%%%%;
tmp_l_max = max(abs(FTK.svd_l_));
tmp_M_k_q_rwlM____ = innerproduct_q_k_stretch_quad_stack____1(n_k_p_r,n_w_,tmp_n_M,tmp_M_k_q__,tmp_l_max);
%%%%%%%%;
tmp_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,tmp_n_UX_rank,tmp_n_M);
for nM=0:tmp_n_M-1;
for nUX_rank=0:tmp_n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
tmp_VUXM_lwnM____(1+nl,:,1+nUX_rank,1+nM) = tmp_V_r__(1+nl,:)*diag(UX__(:,1+nUX_rank).*X_weight_r_(:))*tmp_M_k_q_rwlM____(:,:,1+tmp_l_max+FTK.svd_l_(1+nl),1+nM) / n_w_max ;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nUX_rank=0:tmp_n_UX_rank-1;
end;%for nM=0:tmp_n_M-1;
tmp_t=tic();
tmp_SVUXM_lwMS____ = zeros(FTK.n_svd_l,n_w_max,tmp_n_M,tmp_n_S);
for nS=0:tmp_n_S-1;
for nM=0:tmp_n_M-1;
for nl=0:FTK.n_svd_l-1;
for nUX_rank=0:tmp_n_UX_rank-1;
tmp_SVUXM_lwMS____(1+nl,:,1+nM,1+nS) = tmp_SVUXM_lwMS____(1+nl,:,1+nM,1+nS) + reshape(conj(tmp_UX_S_CTF_k_q_wnS___(:,1+nUX_rank,1+nS)),[1,n_w_max]).*reshape(tmp_VUXM_lwnM____(1+nl,:,1+nUX_rank,1+nM),[1,n_w_max]);
end;%for nUX_rank=0:tmp_n_UX_rank-1;
tmp_SVUXM_lwMS____(1+nl,:,1+nM,1+nS) = ifft(tmp_SVUXM_lwMS____(1+nl,:,1+nM,1+nS))*n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
end;%for nM=0:tmp_n_M-1;
end;%for nS=0:tmp_n_S-1;
tmp_t=toc(tmp_t); disp(sprintf(' %% tmp_SVUXM_lwMS____: SMlr nested loops: %0.4f',tmp_t));
%%%%%%%%;
tmp_t=tic();
for nS=0:tmp_n_S-1;
for nM=0:tmp_n_M-1;
tmp_SVUXM_lwMS____(:,:,1+nM,1+nS) = ifft(sum(repmat(conj(reshape(tmp_UX_S_CTF_k_q_wnS___(:,:,1+nS),[1,tmp_n_w_max,tmp_n_UX_rank])),[FTK.n_svd_l,1,1]).*tmp_VUXM_lwnM____(:,:,:,1+nM),3),[],2)*tmp_n_w_max;
end;%for nM=0:tmp_n_M-1;
end;%for nS=0:tmp_n_S-1;
tmp_t=toc(tmp_t); disp(sprintf(' %% tmp_SVUXM_lwMS____: SM nested loops: %0.4f',tmp_t));
%%%%%%%%;
tmp_UX_S_CTF_k_q_nSw___ = permute(tmp_UX_S_CTF_k_q_wnS___,[2,3,1]);
tmp_VUXM_nMwl____ = permute(tmp_VUXM_lwnM____,[3,4,2,1]);
tmp_t=tic();
tmp_SVUXM_SMwl____ = zeros(tmp_n_S,tmp_n_M,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
tmp_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(tmp_UX_S_CTF_k_q_nSw___(:,:,1+nw))*tmp_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_SVUXM_lwSM____ = permute(ifft(permute(tmp_SVUXM_SMwl____,[3,4,1,2]),[],1)*n_w_max,[2,1,3,4]);
tmp_t=toc(tmp_t); disp(sprintf(' %% tmp_SVUXM_lwMS____: lw nested loops ifft 1: %0.4f',tmp_t));
tmp_t=tic();
tmp_SVUXM_SMwl____ = zeros(tmp_n_S,tmp_n_M,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
tmp_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(tmp_UX_S_CTF_k_q_nSw___(:,:,1+nw))*tmp_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_SVUXM_lwSM____ = ifft(permute(tmp_SVUXM_SMwl____,[4,3,1,2]),[],2)*n_w_max;
tmp_t=toc(tmp_t); disp(sprintf(' %% tmp_SVUXM_lwMS____: lw nested loops ifft 2: %0.4f',tmp_t));
%%%%%%%%;
tmp_USESVUXM_dwSM____ = reshape((tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*reshape(tmp_SVUXM_lwSM____,[FTK.n_svd_l,n_w_max*tmp_n_S*tmp_n_M]),[tmp_n_delta_v,n_w_max,tmp_n_S,tmp_n_M]);
%%%%%%%%;
tmp_X3____ = tmp_USESVUXM_dwSM____;
end;%if flag_check;

flag_check = 0;
if flag_check;
%%%%%%%%;
% Now form tmp_X3____ from scratch: ;
%%%%%%%%;
tmp_N_pixel = 1.00; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,tmp_N_pixel,1e-3);
tmp_delta_max = tmp_N_pixel/(2*pi*k_p_r_max)*pi*sqrt(2);
tmp_n_delta_v = 64;
tmp_delta__ = tmp_delta_max/sqrt(2)*(2*rand(2,tmp_n_delta_v)-1);
FTK.n_delta_v = tmp_n_delta_v;
FTK.delta_x_ = transpose(tmp_delta__(1+0,:));
FTK.delta_y_ = transpose(tmp_delta__(1+1,:));
FTK.svd_d_max = tmp_delta_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
%%%%%%%%;
tmp_n_S = 63; tmp_nS_ = floor(n_viewing_all*rand(tmp_n_S,1));
tmp_n_UX_rank = 1+2; %<-- pick maximum principal-mode. ;
tmp_n_k_p_r = tmp_n_UX_rank;
tmp_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_n_w_ = n_w_max*ones(tmp_n_k_p_r,1);
tmp_n_w_sum = sum(tmp_n_w_);
tmp_n_w_max = max(tmp_n_w_);
tmp_weight_2d_k_p_r_ = ones(tmp_n_k_p_r,1);
tmp_UX_S_CTF_k_q_wnS___ = permute(squeeze(UX_S_CTF_k_q_wSn___(:,1+tmp_nS_,1+(0:tmp_n_UX_rank-1))),[1,3,2]);
tmp_UX_S_CTF_k_q_nSw___ = permute(squeeze(UX_S_CTF_k_q_wSn___(:,1+tmp_nS_,1+(0:tmp_n_UX_rank-1))),[3,2,1]);
%%%%%%%%;
tmp_n_M = 65; tmp_nM_ = floor(n_image_sub*rand(tmp_n_M,1));
tmp_M_k_p__ = M_k_p__(:,1+tmp_nM_);
tmp_M_k_q__ = zeros(n_w_sum,tmp_n_M); for nM=0:tmp_n_M-1; tmp_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p__(:,1+nM)); end;
%%%%%%%%;
tmp_X3_dwSM____ = zeros(tmp_n_delta_v,n_w_max,tmp_n_S,tmp_n_M);
%%%%%%%%;
tmp_U_d__ = permute(reshape(FTK.svd_polyval_U_d_,[FTK.n_svd_l,tmp_n_delta_v]),[2,1]);
tmp_delta_w_ = atan2(tmp_delta__(1+1,:),tmp_delta__(1+0,:));
tmp_expiw__ = exp(-i*(pi/2 - reshape(tmp_delta_w_,[tmp_n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
tmp_VUXM_nMwl____ = tpmh_VUXM_nMwl____0(FTK,n_k_p_r,n_w_,tmp_n_M,tmp_M_k_q__,tmp_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
tmp_t=tic();
tmp_SVUXM_SMwl____ = zeros(tmp_n_S,tmp_n_M,n_w_max,FTK.n_svd_l);
for nl=0:FTK.n_svd_l-1;
for nw=0:n_w_max-1;
tmp_SVUXM_SMwl____(:,:,1+nw,1+nl) = ctranspose(tmp_UX_S_CTF_k_q_nSw___(:,:,1+nw))*tmp_VUXM_nMwl____(:,:,1+nw,1+nl);
end;%for nw=0:n_w_max-1;
end;%for nl=0:FTK.n_svd_l-1;
tmp_SVUXM_lwSM____ = ifft(permute(tmp_SVUXM_SMwl____,[4,3,1,2]),[],2)*n_w_max;
tmp_t=toc(tmp_t); disp(sprintf(' %% tmp_SVUXM_lwMS____: lw nested loops ifft 2: %0.4f',tmp_t));
%%%%%%%%;
tmp_USESVUXM_dwSM____ = reshape((tmp_U_d__.*tmp_expiw__)*diag(FTK.svd_s_)*reshape(tmp_SVUXM_lwSM____,[FTK.n_svd_l,n_w_max*tmp_n_S*tmp_n_M]),[tmp_n_delta_v,n_w_max,tmp_n_S,tmp_n_M]);
tmp_X3_dwSM____ = tmp_USESVUXM_dwSM____;
%%%%%%%%;
% now check a few of these. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
ndelta_v = floor(tmp_n_delta_v*rand());
nM = floor(tmp_n_M*rand());
nS = floor(tmp_n_S*rand());
tmp_X2_ = zeros(n_w_max,1);
tmp_T_k_q_ = transf_svd_q_to_q_FTK_6(FTK,n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_M_k_q__(:,1+nM),tmp_delta__(1+0,1+ndelta_v),tmp_delta__(1+1,1+ndelta_v));
tmp_T_k_q__ = innerproduct_q_k_stretch_quad_stack__0(n_k_p_r,n_w_/(2*pi),n_w_,1,tmp_T_k_q_);
tmp_UX_T_k_q__ = zeros(tmp_n_w_max,tmp_n_UX_rank);
for nUX_rank=0:tmp_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+nM);
if strfind(str_CTF_factor,'rescale'); CTF_factor = CTF_avg_k_p_r_(1+nk_p_r) / max(1e-6,CTF_k_p_r__(1+nk_p_r,1+nctf)); end;
if strfind(str_CTF_factor,'original'); CTF_factor = 1; end;
tmp_UX_T_k_q__(:,1+nUX_rank) = tmp_UX_T_k_q__(:,1+nUX_rank) + transpose(UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*CTF_factor*tmp_T_k_q__(1+nk_p_r,:));
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:tmp_n_UX_rank-1;
tmp_X2_(:) = ifft(innerproduct_q_k_stretch_quad_0(tmp_n_k_p_r,tmp_k_p_r_,tmp_weight_2d_k_p_r_/(2*pi),tmp_n_w_,tmp_n_w_sum,reshape(tmp_UX_S_CTF_k_q_wnS___(:,:,1+nS),[tmp_n_w_max*tmp_n_UX_rank,1]),tmp_UX_T_k_q__(:)))*n_w_max; %<-- multiplication by n_w_max not needed in fortran fftw_plan_back. ;
disp(sprintf(' %% ni %.3d ndelta_v %.3d nS %.3d nM %.3d : tmp_X2_ vs tmp_X3_dwSM____: %0.16f',niteration,ndelta_v,nS,nM,fnorm(tmp_X2_ - reshape(tmp_X3_dwSM____(1+ndelta_v,:,1+nS,1+nM),[n_w_max,1]))/fnorm(tmp_X2_)));
end;%for niteration=0:n_iteration-1;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% check timing of ampmh_X_wSM___1.m ;
% Seems as though it does indeed take roughly FTK.n_svd_l times the cost of a single translation. ;
%%%%%%%%;
%delta_r_max = 0*delta_sigma*sqrt(log(20^2));
delta_r_max = delta_sigma*sqrt(log(20^2));
delta_r_p = 0.05;
delta_r_s = delta_r_max/sqrt(2*log(1/delta_r_p));
delta_r_N = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2));
svd_eps = 1e-3;
%n_delta_v_requested = 1;
n_delta_v_requested = 128;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
disp(sprintf(' %% p-val %0.4f delta_r_max %0.6f sigma %0.4f N_pixel %0.4f --> FTK.n_svd_l %d',delta_r_p,delta_r_max,delta_r_s,delta_r_N,FTK.n_svd_l));
%%%%%%%%;
% prepare precomputation for ampm. ;
%%%%%%%%;
dat_n_M = n_image_sub/4;
dat_n_UX_rank = 1+3;
dat_M_k_p__ = M_k_p__(:,1:dat_n_M);
pm_n_UX_rank = dat_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ;
pm_CTF_index_ = 1; pm_CTF_k_p__ = ones(pm_n_w_sum,1);
%%%%%%%%;
tmp_t = tic();
dat_M_k_q__ = zeros(n_w_sum,dat_n_M);
for nM=0:dat_n_M-1;
dat_M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,dat_M_k_p__(:,1+nM));
end;%for nM=0:dat_n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% dat_M_k_q__: %0.3fs',tmp_t)); end;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,dat_n_M,dat_M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,dat_n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
[dat_UX_M_k_q_wnM_d0___,dat_UX_M_k_p_wnM_d0___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,dat_n_M,svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,:),zeros(dat_n_M,1),zeros(dat_n_M,1));
tmp_t = tic();
%%%%%%%%;
tmp_verbose=0;
tmp_viewing_k_eq_d = 1/(2*pi);
[tmp_S_k_p__,~,~,~,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,pm_n_k_p_r,pm_k_p_r_,pm_k_p_r_max,pm_weight_k_p_r_,pm_l_max_,a_UX_Y_quad__(:,1:pm_n_UX_rank),tmp_viewing_k_eq_d,-1,pm_n_w_);
if (tmp_verbose>0); disp(sprintf(' %% viewing_k_eq_d %0.3f, n_viewing_all %d',viewing_k_eq_d,n_viewing_all)); end;
n_S = n_viewing_all;
UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
UX_S_l2_(1+nS) = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
tmp_S_k_q__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
tmp_S_k_q__(:,1+nS) = interp_p_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS)); 
end;%for nS=0:n_S-1; 
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_template_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
tmp_t = tic();
SS_k_q_ = svd(tmp_S_k_q__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(tmp_S_k_q__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(tmp_S_k_q__,n_S_rank);
if (verbose); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
SS_k_q_ = svd(tmp_S_k_q__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(tmp_S_k_q__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(tmp_S_k_q__,n_S_rank);
if (verbose); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Use current templates to calculate current innerproducts/correlations. ;
% Batches images into batches of size n_M_Mbatch. ;
% Batches templates into batches of size n_S_Sbatch. ;
% Only stores the optimal translation for each image. ;
%%%%%%%%;
tmp_t = tic();
n_M_Mbatch = 48;
n_S_Sbatch = 48;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___1(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_Sbatch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,dat_n_M...
,n_M_Mbatch...
,svd_VUXM_lwnM____...
,UX_M_l2_dM__...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% X_wsM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% check advantage gained for specialized svd of UX__(:,~).*X_weight_r_.*besselj(~,~) ;
%%%%%%%%;
tmp_eps = 1e-3; tmp_N =  0.5; 
tmp_FTK_0 = gen_UXJsvd_FTK_8(n_k_p_r,k_p_r_,0*ones(n_k_p_r,1) + 1*max(abs(UX__(:,1:3)),[],2).*X_weight_r_,k_p_r_max,tmp_N,tmp_eps,25,32,33);
tmp_FTK_1 = gen_UXJsvd_FTK_8(n_k_p_r,k_p_r_,1*ones(n_k_p_r,1) + 0*max(abs(UX__(:,1:3)),[],2).*X_weight_r_,k_p_r_max,tmp_N,tmp_eps,25,32,33);
disp(tmp_FTK_0.n_svd_l/tmp_FTK_1.n_svd_l);
end;%flag_check=0;

flag_check=0;
if flag_check;
%%%%%%%%;
% check timing for translating principal-images. ;
%%%%%%%%;
tmp_n_UX_rank = 2;
pm_n_k_p_r = tmp_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
delta_r_max = delta_sigma*sqrt(log(20^2));
pm_N_pixel = delta_r_max * (2*pi*k_p_r_max) / (pi*sqrt(2)) ; 
FTK = gen_Jsvd_FTK_7(k_p_r_max,pm_N_pixel,1e-3);
if delta_r_max==0; [FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_] = get_delta_2(delta_r_max,  1); end;
if delta_r_max> 0; [FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_] = get_delta_2(delta_r_max,128); end;
FTK.svd_d_max = delta_r_max;
FTK.svd_polyval_U_d_ = get_svd_polyval_U_d_0(FTK.svd_d_max,FTK.n_svd_d,FTK.svd_d_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_U_d_,FTK.n_delta_v,FTK.delta_x_,FTK.delta_y_);
FTK.svd_r_max = 2*pi*k_p_r_max;
FTK.svd_polyval_V_r_ = get_svd_polyval_V_r_0(FTK.svd_r_max,FTK.n_svd_r,FTK.svd_r_,FTK.n_svd_l,FTK.svd_l_,FTK.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
FTK.svd_expiw__ = exp(-i*(pi/2 - reshape(atan2(FTK.delta_y_,FTK.delta_x_),[FTK.n_delta_v,1]))*reshape(FTK.svd_l_,[1,FTK.n_svd_l]));
FTK.svd_U_d_expiw_s__ = (permute(reshape(FTK.svd_polyval_U_d_,[FTK.n_svd_l,FTK.n_delta_v]),[2,1]).*FTK.svd_expiw__)*diag(FTK.svd_s_);
%%%%%%%%;
n_M = n_image_sub;
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
%%%%%%%%;
tmp_image_delta_j_ = max(0,min(FTK.n_delta_v-1,periodize(0:n_M-1,0,FTK.n_delta_v)));
tmp_image_delta_x_ = FTK.delta_x_(1+tmp_image_delta_j_);
tmp_image_delta_y_ = FTK.delta_y_(1+tmp_image_delta_j_);
%%%%%%%%;
UX_M_k_q_wnM_0___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- brute. ;
UX_M_k_p_wnM_0___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- brute. ;
UX_M_k_q_wnM_1___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- blas. ;
UX_M_k_p_wnM_1___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- blas. ;
tmp_t = tic();
for nM=0:n_M-1;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p__(:,1+nM),tmp_image_delta_x_(1+nM),tmp_image_delta_y_(1+nM));
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
for nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_M_k_q_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_UX_M_k_q_(1+tmp_ij_set_) = tmp_UX_M_k_q_(1+tmp_ij_set_) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
end;%for nk_p_r=0:n_k_p_r-1;
UX_M_k_q_wnM_0___(:,1+nUX_rank,1+nM) = tmp_UX_M_k_q_;
UX_M_k_p_wnM_0___(:,1+nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,tmp_UX_M_k_q_);
end;%for nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_0___: %0.6fs',tmp_t));
tmp_t = tic();
[UX_M_k_q_wnM_1___,UX_M_k_p_wnM_1___] = ampmh_M_k_p__to_UX_M_k_p_wnM___0(n_k_p_r,k_p_r_,n_w_,pm_n_UX_rank,UX__,X_weight_r_,n_M,M_k_p__,M_k_q__,tmp_image_delta_x_,tmp_image_delta_y_);
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_1___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_1___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_1___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_1___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_1___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
tmp_t = tic();
svd_VUXM_nMwl____ = tpmh_VUXM_nMwl____0(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_nMwl____: %0.6fs',tmp_t));
UX_M_k_q_wnM_2___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- svd. ;
UX_M_k_p_wnM_2___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- svd. ;
tmp_t = tic();
for nM=0:n_M-1;
for nUX_rank=0:pm_n_UX_rank-1;
for nl=0:FTK.n_svd_l-1;
UX_M_k_q_wnM_2___(:,1+nUX_rank,1+nM) = UX_M_k_q_wnM_2___(:,1+nUX_rank,1+nM) + FTK.svd_U_d_expiw_s__(1+tmp_image_delta_j_(1+nM),1+nl)*reshape(svd_VUXM_nMwl____(1+nUX_rank,1+nM,:,1+nl),[n_w_max,1]) * n_w_max;
end;%for nl=0:FTK.n_svd_l-1;
UX_M_k_p_wnM_2___(:,1+nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,UX_M_k_q_wnM_2___(:,1+nUX_rank,1+nM));
end;%for nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_2___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_2___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_2___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_2___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_2___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
UX_M_k_q_wnM_3___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- svd with blas. ;
UX_M_k_p_wnM_3___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- svd with blas. ;
tmp_t = tic();
%svd_VUXM_lwnM____ = permute(svd_VUXM_nMwl____,[4,3,1,2]);
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_lwnM____: %0.6fs',tmp_t));
tmp_t = tic();
for nM=0:n_M-1;
UX_M_k_q_wnM_3___(:,:,1+nM) = reshape(FTK.svd_U_d_expiw_s__(1+tmp_image_delta_j_(1+nM),:)*reshape(svd_VUXM_lwnM____(:,:,:,1+nM),[FTK.n_svd_l,n_w_max*pm_n_UX_rank]),[n_w_max,pm_n_UX_rank]) * n_w_max;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_q_wnM_3___: %0.6fs',tmp_t));
tmp_t = tic();
for nM=0:n_M-1;
for nUX_rank=0:pm_n_UX_rank-1;
UX_M_k_p_wnM_3___(:,1+nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,UX_M_k_q_wnM_3___(:,1+nUX_rank,1+nM));
end;%for nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_3___: %0.6fs',tmp_t));
tmp_t = tic();
UX_M_k_p_wnM_3___ = ifft(UX_M_k_q_wnM_3___,[],1)*sqrt(n_w_max);
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_3___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_3___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_3___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_3___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_3___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
svd_eps = 1e-3;
FTK = ampmh_FTK_0(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps);
[UX_M_k_q_wnM_5___,UX_M_k_p_wnM_5___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M,svd_VUXM_lwnM____,tmp_image_delta_x_,tmp_image_delta_y_);
disp(sprintf(' %% UX_M_k_q_wnM_3___ vs UX_M_k_q_wnM_5___: %0.16f',fnorm(UX_M_k_q_wnM_3___-UX_M_k_q_wnM_5___)/fnorm(UX_M_k_q_wnM_3___)));
disp(sprintf(' %% UX_M_k_p_wnM_3___ vs UX_M_k_p_wnM_5___: %0.16f',fnorm(UX_M_k_p_wnM_3___-UX_M_k_p_wnM_5___)/fnorm(UX_M_k_p_wnM_3___)));
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% check UXJsvd. ;
%%%%%%%%;
UX_M_k_q_wnM_4___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- UXJ svd with blas. ;
UX_M_k_p_wnM_4___ = zeros(n_w_max,pm_n_UX_rank,n_M); %<-- UXJ svd with blas. ;
UXFTK_ = cell(pm_n_UX_rank,1);
for nUX_rank=0:pm_n_UX_rank-1;
UXFTK_{1+nUX_rank} = gen_UXJsvd_FTK_8(n_k_p_r,k_p_r_,UX__(:,1+nUX_rank).*X_weight_r_,k_p_r_max,pm_N_pixel,1e-3);
if delta_r_max==0; [UXFTK_{1+nUX_rank}.n_delta_v,UXFTK_{1+nUX_rank}.delta_x_,UXFTK_{1+nUX_rank}.delta_y_] = get_delta_2(delta_r_max,  1); end;
if delta_r_max> 0; [UXFTK_{1+nUX_rank}.n_delta_v,UXFTK_{1+nUX_rank}.delta_x_,UXFTK_{1+nUX_rank}.delta_y_] = get_delta_2(delta_r_max,128); end;
UXFTK_{1+nUX_rank}.svd_d_max = delta_r_max;
UXFTK_{1+nUX_rank}.svd_polyval_U_d_ = get_svd_polyval_U_d_0(UXFTK_{1+nUX_rank}.svd_d_max,UXFTK_{1+nUX_rank}.n_svd_d,UXFTK_{1+nUX_rank}.svd_d_,UXFTK_{1+nUX_rank}.n_svd_l,UXFTK_{1+nUX_rank}.svd_l_,UXFTK_{1+nUX_rank}.svd_U_d_,UXFTK_{1+nUX_rank}.n_delta_v,UXFTK_{1+nUX_rank}.delta_x_,UXFTK_{1+nUX_rank}.delta_y_);
UXFTK_{1+nUX_rank}.svd_r_max = 2*pi*k_p_r_max;
UXFTK_{1+nUX_rank}.svd_polyval_V_r_ = get_svd_polyval_V_r_0(UXFTK_{1+nUX_rank}.svd_r_max,UXFTK_{1+nUX_rank}.n_svd_r,UXFTK_{1+nUX_rank}.svd_r_,UXFTK_{1+nUX_rank}.n_svd_l,UXFTK_{1+nUX_rank}.svd_l_,UXFTK_{1+nUX_rank}.svd_V_r_,n_k_p_r,2*pi*k_p_r_);
UXFTK_{1+nUX_rank}.svd_expiw__ = exp(-i*(pi/2 - reshape(atan2(UXFTK_{1+nUX_rank}.delta_y_,UXFTK_{1+nUX_rank}.delta_x_),[UXFTK_{1+nUX_rank}.n_delta_v,1]))*reshape(UXFTK_{1+nUX_rank}.svd_l_,[1,UXFTK_{1+nUX_rank}.n_svd_l]));
UXFTK_{1+nUX_rank}.svd_U_d_expiw_s__ = (permute(reshape(UXFTK_{1+nUX_rank}.svd_polyval_U_d_,[UXFTK_{1+nUX_rank}.n_svd_l,UXFTK_{1+nUX_rank}.n_delta_v]),[2,1]).*UXFTK_{1+nUX_rank}.svd_expiw__)*diag(UXFTK_{1+nUX_rank}.svd_s_);
end;%for nUX_rank=0:pm_n_UX_rank-1;
svd_l_max = 0; for nUX_rank=0:pm_n_UX_rank-1; svd_l_max = max(svd_l_max,UXFTK_{1+nUX_rank}.n_svd_l); end;
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM_____ = cell(pm_n_UX_rank,1);
for nUX_rank=0:pm_n_UX_rank-1;
svd_VUXM_lwnM_____{1+nUX_rank} = zeros(UXFTK_{1+nUX_rank}.n_svd_l,n_w_max,1,n_M);
svd_VUXM_lwnM_____{1+nUX_rank}(:,:,1,:) = permute(reshape(tpmh_VUXM_nMwl____0(UXFTK_{1+nUX_rank},n_k_p_r,n_w_,n_M,M_k_q__,1,ones(n_k_p_r,1),ones(n_k_p_r,1)),[1,n_M,n_w_max,UXFTK_{1+nUX_rank}.n_svd_l]),[4,3,1,2]);
end;%for nUX_rank=0:pm_n_UX_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_lwmN____: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM_____ = cell(pm_n_UX_rank,1);
for nUX_rank=0:pm_n_UX_rank-1;
svd_VUXM_lwnM_____{1+nUX_rank} = zeros(UXFTK_{1+nUX_rank}.n_svd_l,n_w_max,1,n_M);
svd_VUXM_lwnM_____{1+nUX_rank}(:,:,1,:) = tpmh_VUXM_lwnM____0(UXFTK_{1+nUX_rank},n_k_p_r,n_w_,n_M,M_k_q__,1,ones(n_k_p_r,1),ones(n_k_p_r,1));
end;%for nUX_rank=0:pm_n_UX_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% svd_VUXM_lwmN____: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
for nUX_rank=0:pm_n_UX_rank-1;
for nM=0:n_M-1;
UX_M_k_q_wnM_4___(:,1+nUX_rank,1+nM) = reshape(UXFTK_{1+nUX_rank}.svd_U_d_expiw_s__(1+tmp_image_delta_j_(1+nM),:)*reshape(svd_VUXM_lwnM_____{1+nUX_rank}(:,:,1,1+nM),[UXFTK_{1+nUX_rank}.n_svd_l,n_w_max*1]),[n_w_max,1]) * n_w_max;
end;%for nM=0:n_M-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_q_wnM_4___: %0.6fs',tmp_t));
tmp_t = tic();
for nM=0:n_M-1;
for nUX_rank=0:pm_n_UX_rank-1;
UX_M_k_p_wnM_4___(:,1+nUX_rank,1+nM) = interp_q_to_p(1,n_w_max,n_w_max,UX_M_k_q_wnM_4___(:,1+nUX_rank,1+nM));
end;%for nUX_rank=0:pm_n_UX_rank-1;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% UX_M_k_p_wnM_4___: %0.6fs',tmp_t));
disp(sprintf(' %% UX_M_k_q_wnM_0___ vs UX_M_k_q_wnM_4___: %0.16f',fnorm(UX_M_k_q_wnM_0___-UX_M_k_q_wnM_4___)/fnorm(UX_M_k_q_wnM_0___)));
disp(sprintf(' %% UX_M_k_p_wnM_0___ vs UX_M_k_p_wnM_4___: %0.16f',fnorm(UX_M_k_p_wnM_0___-UX_M_k_p_wnM_4___)/fnorm(UX_M_k_p_wnM_0___)));
%%%%%%%%;
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% Now check timing for calculating the norms (across all translations) of the principal-images. ;
%%%%%%%%;
tmp_n_UX_rank = 2;
pm_n_UX_rank = tmp_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
delta_r_max = delta_sigma*sqrt(log(20^2)); svd_eps = 1e-3;
FTK = ampmh_FTK_0(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps);
%%%%%%%%;
n_M = n_image_sub;
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
%%%%%%%%;
tmp_t = tic();
n_M_sub = 32;
UX_M_l2_dM_0__ = zeros(FTK.n_delta_v,n_M_sub); %<-- brute. ;
for nM=0:n_M_sub-1;
for ndelta_v=0:FTK.n_delta_v-1;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p__(:,1+nM),FTK.delta_x_(1+ndelta_v),FTK.delta_y_(1+ndelta_v));
tmp_M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,tmp_M_k_p_);
tmp_UX_M_k_q__ = zeros(n_w_max,pm_n_UX_rank);
for nUX_rank=0:pm_n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_n_w = n_w_(1+nk_p_r);
tmp_n_w_2 = round(tmp_n_w/2);
tmp_ij_set_ = (0:tmp_n_w_2-1); tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(0:tmp_n_w_2-1);
tmp_UX_M_k_q__(1+tmp_ij_set_,1+nUX_rank) = tmp_UX_M_k_q__(1+tmp_ij_set_,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
tmp_ij_set_ = (tmp_n_w_2+1:tmp_n_w-1)-tmp_n_w+n_w_max; tmp_ij_get_ = n_w_csum_(1+nk_p_r)+(tmp_n_w_2+1:tmp_n_w-1);
tmp_UX_M_k_q__(1+tmp_ij_set_,1+nUX_rank) = tmp_UX_M_k_q__(1+tmp_ij_set_,1+nUX_rank) + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_q_(1+tmp_ij_get_);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:pm_n_UX_rank-1;
tmp_UX_M_k_p__ = interp_q_to_p(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,tmp_UX_M_k_q__);
UX_M_l2_dM_0__(1+ndelta_v,1+nM) = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,tmp_UX_M_k_p__,tmp_UX_M_k_p__);
end;%for ndelta_v=0:FTK.n_delta_v-1;
end;%for nM=0:n_M_sub-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% n_M_sub %d UX_M_l2_dM_0__: %0.6fs',n_M_sub,tmp_t));
tmp_t = tic();
n_M_sub = 32;
UX_M_l2_dM_1__ = zeros(FTK.n_delta_v,n_M_sub); %<-- svd. ;
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q__(:,1:n_M_sub),pm_n_UX_rank,UX__,X_weight_r_);
UX_M_k_q_dwnM____ = reshape(FTK.svd_U_d_expiw_s__(:,:)*reshape(svd_VUXM_lwnM____(:,:,:,:),[FTK.n_svd_l,n_w_max*pm_n_UX_rank*n_M_sub]),[FTK.n_delta_v,n_w_max,pm_n_UX_rank,n_M_sub]);
UX_M_l2_dM_1__(:,:) = sum(reshape(abs(permute(UX_M_k_q_dwnM____,[2,3,1,4])).^2,[n_w_max*pm_n_UX_rank,FTK.n_delta_v,n_M_sub]),1) * n_w_max;
tmp_t = toc(tmp_t); disp(sprintf(' %% n_M_sub %d UX_M_l2_dM_1__: %0.6fs',n_M_sub,tmp_t));
disp(sprintf(' %% UX_M_l2_dM_0__ vs UX_M_l2_dM_1__: %0.16f',fnorm(UX_M_l2_dM_0__-UX_M_l2_dM_1__)/fnorm(UX_M_l2_dM_0__)));
tmp_t = tic();
UX_M_l2_dM_2__ = ampmh_UX_M_l2_dM__0(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q__(:,1:n_M_sub),pm_n_UX_rank,UX__,X_weight_r_,[]);
tmp_t = toc(tmp_t); disp(sprintf(' %% n_M_sub %d UX_M_l2_dM_2__: %0.6fs',n_M_sub,tmp_t));
disp(sprintf(' %% UX_M_l2_dM_1__ vs UX_M_l2_dM_2__: %0.16f',fnorm(UX_M_l2_dM_1__-UX_M_l2_dM_2__)/fnorm(UX_M_l2_dM_1__)));
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M_sub,M_k_q__(:,1:n_M_sub),pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = tic();
UX_M_l2_dM_2__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M_sub,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); disp(sprintf(' %% n_M_sub %d UX_M_l2_dM_2__: %0.6fs',n_M_sub,tmp_t));
disp(sprintf(' %% UX_M_l2_dM_1__ vs UX_M_l2_dM_2__: %0.16f',fnorm(UX_M_l2_dM_1__-UX_M_l2_dM_2__)/fnorm(UX_M_l2_dM_1__)));
end;%if flag_check;

flag_check=0;
if flag_check;
%%%%%%%%;
% check innerproducts between principal-templates and principal-images. ;
%%%%%%%%;
tmp_n_UX_rank = 2;
pm_n_UX_rank = tmp_n_UX_rank;
pm_n_k_p_r = pm_n_UX_rank;
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_k_p_r_max = 1;
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_n_w_max = n_w_max;
pm_n_w_sum = sum(pm_n_w_);
pm_n_w_csum_ = cumsum([0;pm_n_w_]);
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_lm_ = (1+pm_l_max_).^2; pm_n_lm_sum = sum(pm_n_lm_);
pm_weight_3d_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_2d_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
delta_r_max = delta_sigma*sqrt(log(20^2)); svd_eps = 1e-3; n_delta_v_requested = 128;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
%%%%%%%%;
n_M = n_image_sub;
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p__(:,1+nM));
end;%for nM=0:n_M-1;
%%%%%%%%;
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____0(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
%%%%%%%%;
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
%%%%%%%%;
tmp_S_k_p__ = reshape(permute(UX_S_CTF_k_p_wSn___(:,:,1:pm_n_UX_rank),[1,3,2]),[n_w_max*pm_n_UX_rank,n_S]);
tmp_t = tic();
n_S = n_viewing_all;
UX_S_l2_ = zeros(n_S,1);
for nS=0:n_S-1;
UX_S_l2_(1+nS) = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS),tmp_S_k_p__(:,1+nS));
end;%for nS=0:n_S-1;
tmp_S_k_q__ = zeros(pm_n_w_sum,n_S);
for nS=0:n_S-1;
tmp_S_k_q__(:,1+nS) = interp_p_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,tmp_S_k_p__(:,1+nS)); 
end;%for nS=0:n_S-1; 
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% get_template_0: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now take svd of principal-templates. ;
%%%%%%%%;
tmp_t = tic();
SS_k_q_ = svd(tmp_S_k_q__);
n_S_rank = min(find(SS_k_q_/SS_k_q_(1)<1e-3)); if isempty(n_S_rank); n_S_rank = min(size(tmp_S_k_q__)); end;
[US_k_q__,SS_k_q__,VS_k_q__] = svds(tmp_S_k_q__,n_S_rank);
if (verbose); disp(sprintf(' %% n_S %d --> n_S_rank %d',n_S,n_S_rank)); end;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% svd of templates: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
n_M_sub = 128*1;
n_M_batch = 32/2;
n_S_batch = 32/2;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___2(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_batch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M_sub...
,n_M_batch...
,svd_VUXM_lwnM____(:,:,:,1:n_M_sub)...
,UX_M_l2_dM__(:,1:n_M_sub)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_M_sub %d X_wSM___: %0.3fs',n_M_sub,tmp_t)); end;
%%%%%%%%;
% Now check a few. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
nM = max(0,min(n_M_sub-1,floor(n_M_sub*rand())));
nS = max(0,min(n_S-1,floor(n_S*rand())));
nw = max(0,min(n_w_max-1,floor(n_w_max*rand())));
X_1 = X_wSM___(1+nw,1+nS,1+nM);
delta_j = delta_j_wSM___(1+nw,1+nS,1+nM);
tmp_image_delta_x = FTK.delta_x_(1+delta_j); tmp_image_delta_y = FTK.delta_y_(1+delta_j);
[UX_M_k_q_wn__,UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,svd_VUXM_lwnM____(:,:,:,1+nM),tmp_image_delta_x,tmp_image_delta_y);
UX_S_CTF_k_q_ = US_k_q__*SS_k_q__*ctranspose(VS_k_q__(1+nS,:));
X_0 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,rotate_q_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_S_CTF_k_q_,2*pi*nw/n_w_max),UX_M_k_q_wn__(:)) / sqrt(UX_S_l2_(1+nS) * UX_M_l2_dM__(1+delta_j,1+nM));
disp(sprintf(' %% X_0 vs X_1: %0.16f',fnorm(X_0-X_1)/fnorm(X_0)));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
tmp_t = tic();
n_M_sub = 128*1;
n_M_batch = 6;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___0(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M_sub...
,n_M_batch...
,svd_VUXM_lwnM____(:,:,:,1:n_M_sub)...
,UX_M_l2_dM__(:,1:n_M_sub)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_M_sub %d X_wSM___: %0.3fs',n_M_sub,tmp_t)); end;
%%%%%%%%;
% Now check a few. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
nM = max(0,min(n_M_sub-1,floor(n_M_sub*rand())));
nS = max(0,min(n_S-1,floor(n_S*rand())));
nw = max(0,min(n_w_max-1,floor(n_w_max*rand())));
X_1 = X_wSM___(1+nw,1+nS,1+nM);
delta_j = delta_j_wSM___(1+nw,1+nS,1+nM);
tmp_image_delta_x = FTK.delta_x_(1+delta_j); tmp_image_delta_y = FTK.delta_y_(1+delta_j);
[UX_M_k_q_wn__,UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,svd_VUXM_lwnM____(:,:,:,1+nM),tmp_image_delta_x,tmp_image_delta_y);
UX_S_CTF_k_q_ = US_k_q__*SS_k_q__*ctranspose(VS_k_q__(1+nS,:));
X_0 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,rotate_q_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_S_CTF_k_q_,2*pi*nw/n_w_max),UX_M_k_q_wn__(:)) / sqrt(UX_S_l2_(1+nS) * UX_M_l2_dM__(1+delta_j,1+nM));
disp(sprintf(' %% X_0 vs X_1: %0.16f',fnorm(X_0-X_1)/fnorm(X_0)));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
tmp_t = tic();
n_M_sub = 128*1;
n_M_batch = 32;
n_S_batch = 32;
[X_wSM___,delta_j_wSM___] = ...
ampmh_X_wSM___1(...
 FTK...
,n_w_...
,pm_n_UX_rank...
,n_S...
,n_S_rank...
,n_S_batch...
,US_k_q__...
,SS_k_q__...
,VS_k_q__...
,UX_S_l2_...
,n_M_sub...
,n_M_batch...
,svd_VUXM_lwnM____(:,:,:,1:n_M_sub)...
,UX_M_l2_dM__(:,1:n_M_sub)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% n_M_sub %d X_wSM___: %0.3fs',n_M_sub,tmp_t)); end;
%%%%%%%%;
% Now check a few. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
nM = max(0,min(n_M_sub-1,floor(n_M_sub*rand())));
nS = max(0,min(n_S-1,floor(n_S*rand())));
nw = max(0,min(n_w_max-1,floor(n_w_max*rand())));
X_1 = X_wSM___(1+nw,1+nS,1+nM);
delta_j = delta_j_wSM___(1+nw,1+nS,1+nM);
tmp_image_delta_x = FTK.delta_x_(1+delta_j); tmp_image_delta_y = FTK.delta_y_(1+delta_j);
[UX_M_k_q_wn__,UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,svd_VUXM_lwnM____(:,:,:,1+nM),tmp_image_delta_x,tmp_image_delta_y);
UX_S_CTF_k_q_ = US_k_q__*SS_k_q__*ctranspose(VS_k_q__(1+nS,:));
X_0 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,rotate_q_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_S_CTF_k_q_,2*pi*nw/n_w_max),UX_M_k_q_wn__(:)) / sqrt(UX_S_l2_(1+nS) * UX_M_l2_dM__(1+delta_j,1+nM));
disp(sprintf(' %% X_0 vs X_1: %0.16f',fnorm(X_0-X_1)/fnorm(X_0)));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
%%%%%%%%;
% Now check MS and SM. ;
%%%%%%%%;
tmp_t = tic();
[...
 tmp_euler_polar_a_MS_...
,tmp_euler_azimu_b_MS_...
,tmp_euler_gamma_z_MS_...
,tmp_image_delta_x_MS_...
,tmp_image_delta_y_MS_...
,tmp_image_X_value_MS_...
,tmp_image_S_index_MS_...
] = ...
ampmh_MS_0(...
 n_w_...
,n_S...
,viewing_polar_a_all_...
,viewing_azimu_b_all_...
,n_M_sub...
,X_wSM___...
,FTK.delta_x_(1+delta_j_wSM___)...
,FTK.delta_y_(1+delta_j_wSM___)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% MS : %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now check a few. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
nM = max(0,min(n_M_sub-1,floor(n_M_sub*rand())));
image_S_index = tmp_image_S_index_MS_(1+nM);
tmp_euler_gamma_z = tmp_euler_gamma_z_MS_(1+nM);
tmp_image_delta_x = tmp_image_delta_x_MS_(1+nM);
tmp_image_delta_y = tmp_image_delta_y_MS_(1+nM);
X_1 = tmp_image_X_value_MS_(1+nM);
[UX_M_k_q_wn__,UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,svd_VUXM_lwnM____(:,:,:,1+nM),tmp_image_delta_x,tmp_image_delta_y);
UX_M_l2 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,UX_M_k_q_wn__(:),UX_M_k_q_wn__(:));
UX_S_CTF_k_q_ = US_k_q__*SS_k_q__*ctranspose(VS_k_q__(1+image_S_index,:));
X_0 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,rotate_q_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_S_CTF_k_q_,tmp_euler_gamma_z),UX_M_k_q_wn__(:)) / sqrt(UX_S_l2_(1+image_S_index) * UX_M_l2);
disp(sprintf(' %% X_0 vs X_1: %0.16f',fnorm(X_0-X_1)/fnorm(X_0)));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
f_rand = 0.00;
tmp_t = tic();
[...
 tmp_euler_polar_a_SM_...
,tmp_euler_azimu_b_SM_...
,tmp_euler_gamma_z_SM_...
,tmp_image_delta_x_SM_...
,tmp_image_delta_y_SM_...
,tmp_image_X_value_SM_...
,tmp_image_S_index_SM_...
] = ...
ampmh_SM_0(...
 f_rand...
,n_w_...
,n_S...
,viewing_polar_a_all_...
,viewing_azimu_b_all_...
,n_M_sub...
,X_wSM___...
,FTK.delta_x_(1+delta_j_wSM___)...
,FTK.delta_y_(1+delta_j_wSM___)...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% SM : %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now check a few. ;
%%%%%%%%;
n_iteration = 6;
for niteration=0:n_iteration-1;
nM = max(0,min(n_M_sub-1,floor(n_M_sub*rand())));
image_S_index = tmp_image_S_index_SM_(1+nM);
tmp_euler_gamma_z = tmp_euler_gamma_z_SM_(1+nM);
tmp_image_delta_x = tmp_image_delta_x_SM_(1+nM);
tmp_image_delta_y = tmp_image_delta_y_SM_(1+nM);
X_1 = tmp_image_X_value_SM_(1+nM);
[UX_M_k_q_wn__,UX_M_k_p_wn__] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,1,svd_VUXM_lwnM____(:,:,:,1+nM),tmp_image_delta_x,tmp_image_delta_y);
UX_M_l2 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,UX_M_k_q_wn__(:),UX_M_k_q_wn__(:));
UX_S_CTF_k_q_ = US_k_q__*SS_k_q__*ctranspose(VS_k_q__(1+image_S_index,:));
X_0 = innerproduct_p_quad(pm_n_k_p_r,pm_k_p_r_,pm_weight_2d_k_p_r_/(2*pi),pm_n_w_,pm_n_w_sum,rotate_q_to_q(pm_n_k_p_r,pm_n_w_,pm_n_w_sum,UX_S_CTF_k_q_,tmp_euler_gamma_z),UX_M_k_q_wn__(:)) / sqrt(UX_S_l2_(1+image_S_index) * UX_M_l2);
disp(sprintf(' %% X_0 vs X_1: %0.16f',fnorm(X_0-X_1)/fnorm(X_0)));
end;%for niteration=0:n_iteration-1;
%%%%%%%%;
end;%if flag_check;

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
delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_combine,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
n_w_uni_ = n_w_max*ones(n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__;
tpmhtam_experimental_11(...
 dat_infix...
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
,a_UX_Y_quad__(:,1:dat_n_UX_rank)...
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
fname_0 = sprintf('%s_mat/tpmhtameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
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
fname_fig = sprintf('%s_jpg/tpmhtameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
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

%{
%%%%%%%%;
% Now set up alternating minimization using synthetic principled-images. ;
%%%%%%%%;
delta_r_max_use = sqrt(2)*delta_sigma*sqrt(log(20^2)); 
syn_n_order = 5;
syn_n_M = 1024;
syn_n_UX_rank = 5;
syn_n_iteration = 32;
%%%%%%%%;
svd_eps_ = 10.^[-3:0.5:-1]; n_svd_eps = numel(svd_eps_);
n_delta_v_requested_ = [128,64,32,16]; n_n_delta_v_requested = numel(n_delta_v_requested_);
syn_rseed_ = [0:2]; n_syn_rseed = numel(syn_rseed_);
syn_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_syn_snr=numel(syn_snr_);
for nsvd_eps=0:n_svd_eps-1;
svd_eps = svd_eps_(1+nsvd_eps);
for nn_delta_v_requested=0:n_n_delta_v_requested-1;
n_delta_v_requested = n_delta_v_requested_(1+nn_delta_v_requested);
str_delta = sprintf('d%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
syn_infix = sprintf('%s_%s',str_combine,str_delta);
for nsyn_rseed=0:n_syn_rseed-1;
syn_rseed = syn_rseed_(1+nsyn_rseed);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
n_molecule = 1;
molecule_density_ = 1;
a_k_Y_quad__ = a_k_Y_quad_;
a_UX_Y_quad_mavg__ = a_UX_Y_quad__(:,1:syn_n_UX_rank);
tpmhtam_synthetic_11(...
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
,a_UX_Y_quad_mavg__...
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
 %}

disp('returning'); return;

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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_X_compare_A',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_trpv1_17_X_compare_B',dir_trunk);
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
,a_UX_Y_quad__(:,1:dat_n_UX_rank)...
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
tmp_a_UX_Y_0lsq_ = cg_lsq_3(tmp_pm_n_order,tmp_pm_n_k_p_r,tmp_pm_l_max_,tmp_pm_n_w_,dat_n_M,reshape(tmp_UX_M_k_p_wnM___,[tmp_pm_n_w_max*tmp_pm_n_k_p_r,dat_n_M]),tmp_pm_CTF_index_,tmp_pm_CTF_k_p__,tmp_1_.euler_polar_a_MS__(:,end),tmp_1_.euler_azimu_b_MS__(:,end),tmp_1_.euler_gamma_z_MS__(:,end));
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% cg_lsq_3 for a_UX_Y_0lsq_: %0.3fs',tmp_t)); end;
disp(sprintf(' %% fnorm(tmp_a_UX_Y_0lsq_ - tmp_1_.a_UX_Y_0lsq_MS__(:,end))/fnorm(tmp_a_UX_Y_0lsq_) = %0.16f',fnorm(tmp_a_UX_Y_0lsq_ - tmp_1_.a_UX_Y_0lsq_MS__(:,end))/fnorm(tmp_a_UX_Y_0lsq_)));
%%%%%%%%;
% use current model to generate current principal-templates. ;
%%%%%%%%;
tmp_t = tic();
tmp_verbose=0;
tmp_viewing_k_eq_d = 1/(2*pi)/sqrt(1);
[tmp_S_k_p__,~,~,~,tmp_n_viewing_all,tmp_viewing_azimu_b_all_,tmp_viewing_polar_a_all_,~,~,~,~,~,~] = get_template_0(tmp_verbose,tmp_pm_n_k_p_r,tmp_pm_k_p_r_,tmp_pm_k_p_r_max,tmp_pm_weight_k_p_r_,tmp_pm_l_max_,tmp_a_UX_Y_0lsq_,tmp_viewing_k_eq_d,-1,tmp_pm_n_w_);
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
tmp_d_k_Y_reco_ = cg_lsq_4(dat_n_order,n_k_p_r,k_p_r_,l_max_,n_w_,dat_n_M,M_k_p__,CTF_idx_,CTF_k_p__,tmp_euler_polar_a_,tmp_euler_azimu_b_,tmp_euler_gamma_z_,tmp_image_delta_x_,tmp_image_delta_y_);
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
