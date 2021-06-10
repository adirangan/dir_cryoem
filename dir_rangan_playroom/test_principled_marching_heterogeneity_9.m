% Here str_cost = '3d' corrects for volume element, ;
% whereas str_cost = '2d' corrects for area element. ;
% Here str_CTF_factor = 'rescaled' for CTF-rescaled images, ;
% whereas str_CTF_factor = 'original' for original CTF used in images. ;
% str_CTF_factor = 'one' means CTF=1. ;
clear; setup;

str_cost = '2d_xcor';
str_CTF_factor = 'one';
if (~strcmp(str_CTF_factor,'one')); str_combine = sprintf('%s_%s',str_cost,str_CTF_factor); end;
if ( strcmp(str_CTF_factor,'one')); str_combine = sprintf('%s',str_cost); end;

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
% Define a family of molecules in k-space. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
x_u_res = 64;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = x_u_res^3;
X_u_weight_ = (2*x_p_r_max/x_u_res)^3;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_heterogeneity_9_a_k_p_form__.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_t = tic;
n_molecule = 4; molecule_density_ = [0.35;0.25;0.25;0.15];
%k_p_r_max = 48/(2*pi); k_eq_d = 1.0/(2*pi);
k_p_r_max = 36/(2*pi); k_eq_d = 1.0/(2*pi);
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_,J_node_,J_weight_,J_chebfun_,J_polyval_] = sample_sphere_7(verbose,k_p_r_max,k_eq_d,'L') ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;
n_source = 128; sg = 1/18; mu_source__ = zeros(3,n_source,n_molecule); mu_r__ = zeros(n_source,n_molecule);
for nmolecule=0:n_molecule-1;
for nsource=0:n_source-1;
t = nsource/(n_source-1);
mu_r__(1+nsource,1+nmolecule) = t + (t>0.75)*((t-0.75)/0.5)^2*(nmolecule-(n_molecule-1)/2);
mu_source__(1+0,1+nsource,1+nmolecule) = 1/3*mu_r__(1+nsource,1+nmolecule) * sin(2*pi*2*t);
mu_source__(1+1,1+nsource,1+nmolecule) = 1/3*mu_r__(1+nsource,1+nmolecule) * cos(2*pi*2*t);
mu_source__(1+2,1+nsource,1+nmolecule) = -1/3 + 2/3*t ;
end;%for nsource=0:n_source-1;
end;%for nmolecule=0:n_molecule-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_mu_source___',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
hold on;
c_ = colormap_beach(); n_c = size(c_,1);
for nmolecule=0:n_molecule-1;
nc = max(0,min(n_c-1,floor(n_c*nmolecule/n_molecule)));
plot3(mu_source__(1+0,:,1+nmolecule),mu_source__(1+1,:,1+nmolecule),mu_source__(1+2,:,1+nmolecule),'o','MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%for nmolecule=0:n_molecule-1;
axis vis3d;
view([+25,+25]);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
a_k_p_form__ = zeros(n_k_all,n_molecule);
a_x_c_form__ = zeros(n_X_u,n_molecule);
for nmolecule=0:n_molecule-1;
for nsource=0:n_source-1;
a_k_p_form__(:,1+nmolecule) = a_k_p_form__(:,1+nmolecule) + ...
  gk(2*pi*k_c_0_all_(:),mu_source__(1+0,1+nsource,1+nmolecule),sg) .* ...
  gk(2*pi*k_c_1_all_(:),mu_source__(1+1,1+nsource,1+nmolecule),sg) .* ...
  gk(2*pi*k_c_2_all_(:),mu_source__(1+2,1+nsource,1+nmolecule),sg) ;
end;%for nsource=0:n_source-1;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_reco__(:,1+nmolecule) = nufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_form__(:,1+nmolecule).*(2*pi)^3.*weight_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_u_reco_ time %0.2fs',tmp_t));
end;%for nmolecule=0:n_molecule-1;
save(fname_mat ...
     ,'half_diameter_x_c','diameter_x_c','x_p_r_max','x_u_res','x_u_0_','x_u_1_','x_u_2_' ...
     ,'X_u_0_','X_u_1_','X_u_2_','n_X_u','X_u_weight_' ...
     ,'k_p_r_max','k_eq_d' ...
     ,'n_k_all','n_k_all_csum_','k_p_r_all_','k_p_azimu_b_all_','k_p_polar_a_all_','weight_k_all_','weight_shell_k_','n_k_p_r','k_p_r_','weight_k_p_r_','k_c_0_all_','k_c_1_all_','k_c_2_all_','J_node_','J_weight_','J_chebfun_','J_polyval_' ...
     ,'n_molecule','molecule_density_','n_source','sg','mu_r__','mu_source__' ...
     ,'a_k_p_form__' ...
     ,'a_x_u_reco__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_a_x_u_reco__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
for nmolecule=0:n_molecule-1;
subplot(2,2,1+nmolecule); isosurface_f_x_u_0(reshape(real(a_x_u_reco__(:,1+nmolecule)),x_u_res,x_u_res,x_u_res),[90,95,99]); title(sprintf('%d',nmolecule));
end;%for nmolecule=0:n_molecule-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
clf;
for nmolecule=0:n_molecule-1;
subplot(2,2,1+nmolecule); isosurface_f_x_u_0(reshape(real(a_x_u_reco__(:,1+nmolecule)),x_u_res,x_u_res,x_u_res),[99]); title(sprintf('%d',nmolecule));
end;%for nmolecule=0:n_molecule-1;
figbig;
disp(sprintf(' %% writing %s_B',fname_fig));
print('-djpeg',sprintf('%s_B.jpg',fname_fig));
print('-depsc',sprintf('%s_B.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now see what the function looks like on a uniform x_c_ grid. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_heterogeneity_9_b_x_u_reco_0__.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
nmolecule=0;
b_x_u_reco__ = zeros(n_X_u,n_k_p_r);
eta = pi/k_p_r_max; 
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+1+nk_p_r)-1;
tmp_n_k_all = numel(1+tmp_ij_);
tmp_t = tic;
b_x_u_reco__(:,1+nk_p_r) = nufft3d3(tmp_n_k_all,2*pi*k_c_0_all_(1+tmp_ij_)*eta,2*pi*k_c_1_all_(1+tmp_ij_)*eta,2*pi*k_c_2_all_(1+tmp_ij_)*eta,a_k_p_form__(1+tmp_ij_,1+nmolecule).*(2*pi)^3.*weight_k_all_(1+tmp_ij_),+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nk_p_r %.3d/%.3d nufft3d3: b_x_u_reco__(:,%.3d) time %0.2fs',nk_p_r,n_k_p_r,nk_p_r,tmp_t));
end;%for nk_p_r=0:n_k_p_r-1;
b_x_u_reco__ = cumsum(b_x_u_reco__,2);
save(fname_mat ...
     ,'b_x_u_reco__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if (~exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_b_x_u_reco_0__',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
subplot(2,2,1); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,16)),x_u_res,x_u_res,x_u_res),[99]); title(sprintf('K<=%0.2f',k_p_r_(16)));
subplot(2,2,2); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,24)),x_u_res,x_u_res,x_u_res),[99]); title(sprintf('K<=%0.2f',k_p_r_(24)));
subplot(2,2,3); isosurface_f_x_u_0(reshape(real(b_x_u_reco__(:,32)),x_u_res,x_u_res,x_u_res),[99]); title(sprintf('K<=%0.2f',k_p_r_(32)));
subplot(2,2,4); plot(k_p_r_,corr(real(b_x_u_reco__(:,end)),real(b_x_u_reco__(:,:))),'ko'); xlabel('max K'); ylabel('correlation'); title('correlation of csum');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now convert to a_k_Y__ ; 
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_heterogeneity_9_a_k_Y_quad__.mat',dir_trunk);
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
weight_Y_(1+tmp_ij_) = weight_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
tmp_t = tic;
a_k_Y_quad__ = zeros(n_lm_sum,n_molecule);
for nmolecule=0:n_molecule-1;
[a_k_Y_quad__(:,1+nmolecule)] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_p_form__(:,1+nmolecule));
end;%for nmolecule=0:n_molecule-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
for nmolecule=0:n_molecule-1;
[a_k_p_reco__(:,1+nmolecule)] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_quad__(:,1+nmolecule));
end;%for nmolecule=0:n_molecule-1;
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad__ --> a_k_p_reco__ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_reco error: %0.16f',fnorm(a_k_p_form__-a_k_p_reco__)/fnorm(a_k_p_form__)));
%%%%%%%%;
a_k_Y_quad___ = zeros(n_lm_max,n_k_p_r,n_molecule);
for nmolecule=0:n_molecule-1;
for nk_p_r=0:n_k_p_r-1;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad___(1:n_lm_(1+nk_p_r),1+nk_p_r,1+nmolecule) = a_k_Y_quad__(1+tmp_ij_,1+nmolecule);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nmolecule=0:n_molecule-1;
save(fname_mat ...
     ,'l_max_' ...
     ,'n_lm_','n_lm_max','n_lm_sum','n_lm_csum_','l_max_max','m_max_','n_m_max' ...
     ,'Y_l_val_','Y_m_val_','Y_k_val_','weight_Y_' ...
     ,'a_k_Y_quad__' ...
     ,'a_k_p_reco__' ...
     ,'a_k_Y_quad___' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_a_k_Y_quad_0_A',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
nmolecule=0;
a_k_Y_quad_ = a_k_Y_quad__(:,1+nmolecule);
subplot(1,3,1); plot(Y_l_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,2); plot(Y_m_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,3); plot(Y_k_val_,log10(abs(a_k_Y_quad_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_a_k_Y_quad_0_',dir_trunk);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
nmolecule=0;
a_k_Y_quad_ = a_k_Y_quad__(:,1+nmolecule);
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_)),[-10,0],colormap_beach());
clear a_k_Y_quad_;
xlabel('l_val','Interpreter','none');
ylabel('m_val','Interpreter','none');
zlabel('k_p_r','Interpreter','none');
axis vis3d; view([-65,+20]);
title('a_k_Y_quad_','Interpreter','none');
%figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%
% Now generate template grid. ;
%%%%%%%%;
verbose=0;
fname_mat = sprintf('%s_mat/test_principled_marching_heterogeneity_9_weight_2d_k_p_r_.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
template_k_eq_d = k_eq_d*2;
viewing_k_eq_d = k_eq_d*8;
[~,n_w_,weight_2d_k_p_r_,weight_2d_k_all_,n_viewing_all,viewing_azimu_b_all_,viewing_polar_a_all_,n_viewing_polar_a,viewing_polar_a_,n_viewing_azimu_b_,template_k_c_0__,template_k_c_1__,template_k_c_2__] = get_template_0(verbose,n_k_p_r,k_p_r_,k_p_r_max,weight_k_p_r_,l_max_,[],viewing_k_eq_d,template_k_eq_d);
n_S = n_viewing_all; n_w_max = max(n_w_); n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
save(fname_mat ...
     ,'template_k_eq_d','viewing_k_eq_d','n_S' ...
     ,'n_w_','weight_2d_k_p_r_','weight_2d_k_all_','n_viewing_all','viewing_azimu_b_all_','viewing_polar_a_all_','n_viewing_polar_a','viewing_polar_a_','n_viewing_azimu_b_','template_k_c_0__','template_k_c_1__','template_k_c_2__' ...
     ,'n_w_max','n_w_sum','n_w_csum_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% For images we will use the same n_w_ as that used by templates. ;
% We will also use the same quadrature weights for integration in 2d. ;
%%%%%%%%;
%%%%%%%%;
% Set up dummy ctf. ;
%%%%%%%%;
CTF_k_p_r_xcor__ = ones(n_k_p_r,n_k_p_r);
CTF_k_p_r_xavg__ = ones(n_k_p_r,n_k_p_r);
CTF_avg_k_p_ = ones(n_w_sum,1);
CTF_avg_k_p_r_ = ones(n_k_p_r,1);

verbose=1;
%%%%%%%%;
% Now test out principled marching. ;
% First set up cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_heterogeneity_9_X_%s.mat',dir_trunk,str_cost);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
if strfind(str_cost,'3d'); 
tmp_weight_k_p_r_ = weight_k_p_r_;
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
[X_,X_weight_r_] = principled_marching_cost_matrix_4(n_k_p_r,tmp_weight_k_p_r_,l_max_,n_molecule,molecule_density_,a_k_Y_quad__,tmp_CTF__); 
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[UX_,SX_,VX_] = svds(X_,n_UX_rank);
a_UX_Y_quad___ = zeros(n_lm_max,n_UX_rank,n_molecule);
for nmolecule=0:n_molecule-1;
for nUX_rank=0:n_UX_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_ij_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_UX_Y_quad___(1:tmp_n_lm,1+nUX_rank,1+nmolecule) = a_UX_Y_quad___(1:tmp_n_lm,1+nUX_rank,1+nmolecule) + UX_(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad__(1+tmp_ij_,1+nmolecule)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_UX_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUX_rank=0:n_UX_rank-1;
end;%for nmolecule=0:n_molecule-1;
%%%%%%%%;
save(fname_mat ...
     ,'X_','X_weight_r_' ...
     ,'n_UX_rank','UX_','SX_','VX_','a_UX_Y_quad___' ...
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
nmolecule=0;
for nUX_rank=0:min(36,n_UX_rank-1);%for nUX_rank=0:n_UX_rank-1;
[tmp_X,tmp_X_ori,tmp_X_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_UX_Y_quad___(:,1+nUX_rank,1+nmolecule),a_UX_Y_quad___(:,1+nUX_rank,1+nmolecule));
tmp_Z = transpose(UX_(:,1+nUX_rank))*X_*(UX_(:,1+nUX_rank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %+0.16f tmp_X %+0.16f tmp_X_ori*tmp_weight_so3 %+0.16f tmp_X_tau %+0.16f ratio %+0.16f',nUX_rank,n_UX_rank,tmp_Z,tmp_X,tmp_X_ori*tmp_weight_so3,tmp_X_tau,(tmp_X_ori*tmp_weight_so3)/tmp_X_tau));
end;%for nUX_rank=0:n_UX_rank-1;
end;%if flag_check;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_X_%s_A',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
colormap(colormap_beach());
subplot(1,2,1); imagesc(log10(abs(UX_)),[-3,0]); xlabel('rank'); ylabel('shell'); title('log10(abs(UX)) [-3,0]'); 
subplot(1,2,2); plot(log10(abs(diag(SX_))),'ko'); xlabel('rank'); ylabel('log10(\sigma)'); title('log10(SX)');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_heterogeneity_9_X_%s_B',dir_trunk,str_cost);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
nmolecule=0;
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
quad_lim_ = 0.5 * abs(a_UX_Y_quad___(1,1,1)) * [-1,+1];
for nplot=0:n_plot-1;
nk_p_r = plot_nk_p_r_(1+nplot);
[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_UX_Y_quad___(:,1+nk_p_r,1+nmolecule)),k_u_res,2*k_u_res);
subplot(3,n_plot,1 + nplot + 0*n_plot); imagesc(real(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('real nk_p_r: %d, real(quad)',nk_p_r),'Interpreter','none');
subplot(3,n_plot,1 + nplot + 1*n_plot); imagesc(imag(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('imag nk_p_r: %d, imag(quad)',nk_p_r),'Interpreter','none');
subplot(3,n_plot,1 + nplot + 2*n_plot); imagesc( abs(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf(' abs nk_p_r: %d,  abs(quad)',nk_p_r),'Interpreter','none');
clear b_k_p_quad_;
end;%for nplot=0:n_plot-1;
colormap(colormap_beach());
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

disp('returning'); return; 

a_UX_Y_quad_mavg__ = zeros(n_lm_max,n_UX_rank);
for nmolecule=0:n_molecule-1;
a_UX_Y_quad_mavg__ = a_UX_Y_quad_mavg__ + molecule_density_(1+nmolecule)*a_UX_Y_quad___(:,:,1+nmolecule);
end;%for nmolecule=0:n_molecule-1;
%%%%%%%%;
% Now set up alternating minimization using synthetic principled-images. ;
%%%%%%%%;
syn_infix = str_combine;
syn_n_order = 5;
syn_n_M = 1024*2;
syn_n_UX_rank = 16;
syn_n_iteration = 32;
syn_rseed_ = [0:1]; n_syn_rseed = numel(syn_rseed_);
syn_snr_ = [0 , 2.^[-1:-0.5:-5]]; n_syn_snr=numel(syn_snr_);
for nsyn_rseed=0:n_syn_rseed-1;
syn_rseed = syn_rseed_(1+nsyn_rseed);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
tpmhham_synthetic_5(...
 syn_infix...
,syn_rseed...
,syn_n_M...
,syn_snr...
,syn_n_UX_rank...
,syn_n_iteration...
,syn_n_order...
,dir_trunk...
,x_u_res...
,diameter_x_c...
,n_k_p_r...
,k_p_r_...
,k_p_r_max...
,weight_k_p_r_...
,viewing_k_eq_d...
,template_k_eq_d...
,CTF_avg_k_p_...
,CTF_avg_k_p_r_...
,l_max_...
,n_molecule...
,molecule_density_...
,a_k_Y_quad__...
,UX_...
,X_weight_r_...
,a_UX_Y_quad_mavg__...
);
end;%for nsyn_snr=0:n_syn_snr-1;
end;%for nsyn_rseed=0:n_syn_rseed-1;
%%%%%%%%;
% Collect output. ;
%%%%%%%%;
syn_X_best_allshell____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
flag_syn_X_best_allshell____ = zeros(syn_n_iteration,syn_n_UX_rank,n_syn_snr,n_syn_rseed);
for nsyn_rseed=0:n_syn_rseed-1;
syn_rseed = syn_rseed_(1+nsyn_rseed);
for nsyn_snr=0:n_syn_snr-1;
syn_snr = syn_snr_(1+nsyn_snr);
for nUX_rank=0:syn_n_UX_rank-1;
fname_pre = sprintf('%s_mat/tpmhham_synthetic_5_UX_%s_Y_n%.3ds%.4dr%.3d_allshell_rng%.3dnUX%.3d',dir_trunk,syn_infix,syn_n_M,floor(1000*syn_snr),syn_n_UX_rank,syn_rseed,nUX_rank);
fname_mat = sprintf('%s.mat',fname_pre);
if ( exist(fname_mat,'file'));
tmp_ = load(fname_mat);
syn_X_best_allshell____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = tmp_.syn_X_best_allshell_(:);
flag_syn_X_best_allshell____(:,1+nUX_rank,1+nsyn_snr,1+nsyn_rseed) = 1;
clear tmp_;
end;%if ( exist(fname_mat,'file'));
end;%for nUX_rank=0:syn_n_UX_rank-1;
end;%for nsyn_snr=0:n_syn_snr-1;
end;%for nsyn_rseed=0:n_syn_rseed-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/tpmhham_synthetic_5_UX_%s_Y_n%.3dsxxxrxxx_syn_X_best_allshell____',dir_trunk,syn_infix,syn_n_M);
if (~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
for nsyn_snr=0:n_syn_snr-1;
subplot(2,5,1+nsyn_snr);
hold on;
for nUX_rank=0:syn_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/syn_n_UX_rank)));
tmp_X_ = squeeze(syn_X_best_allshell____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_F_ = squeeze(flag_syn_X_best_allshell____(:,1+nUX_rank,1+nsyn_snr,:));
tmp_ij_ = find(sum(tmp_F_,1));
if (numel(tmp_ij_)>0);
plot(0:syn_n_iteration-1,tmp_X_(:,tmp_ij_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_ij_)>0);
end;%for nUX_rank=0:syn_n_UX_rank-1;
%legend(num2str(0:syn_n_UX_rank-1),'Location','NorthWest');
xlim([0,syn_n_iteration-1]); ylim([-0.125,1.000]);
xlabel('niteration');
ylabel(sprintf('correlation with a_UX_Y_quad__(:,:)'),'Interpreter','none');
grid on;
title(sprintf('syn_snr %0.2f',syn_snr_(1+nsyn_snr)),'Interpreter','none');
end;%for nsyn_snr=0:n_syn_snr-1;
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

disp('returning'); return;

