% Here str_cost = '3d' corrects for volume element, ;
% whereas str_cost = '2d' corrects for area element. ;
% 1. the images are *NOT* pretranslated, and ;
% 2. the delta_sigma used to calculate the cost is *NOT* set to 0. ;
% Moreover, the n_w_ used to calculate M_k_p__ is not adaptive. ;
% Now check to see if the experimental images can themselves be used to recapitulate the principal-modes. ;
% Now we investigate the role played by f_rand. ;
% Now we investigate the stability / convergence of alternating minimization. ;
% 20201028: note that registration can be done with N_wavelength==0 (i.e., zero displacement), since the lsq-step should be close to centered. ;
% 20201029: need to update tpmutamf_0 to visualize viewing angles. ; %<-- seems as though viewing angles are reasonably accurate. ;
% 20201029: need to update score to take into account (highly) nonuniform distribution of viewing angles. ;
% 20201029: need to try ampm with more principal modes (and tpmutam_experimental with smaller displacement radius). ;
% 20201029: then need to test on synthetic rib80s data set. ;
% 20201029: and also perhaps on the trpv1 data-set just to be sure. ;
% 20201030: update methods to allow for nonuniform pm_grid. ; 
%<-- on second thought, this seems like a waste, given how most of the principal-modes use most of the frequencies to some degree. ;
% 20201101: seems as though 2d_xcor cost is quite good, while Memp is not so good. ;
% 20201106: trying Memp cost again. ;
% 20201110: seems to work reasonably well actually. Moving back to trpv1 to double check. ;
% 20201222: mean centering of images. ;
% 20210301: attempting local search. ;
clear;

platform = 'access1';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;

flag_recalc = 0;
flag_replot = 0;

dir_trunk = sprintf('/%s/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching_rib80sc',string_root);
%dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_root);
string_rusty_root = 'mnt/home';
dir_relion_bin = sprintf('/%s/rangan/relion/build/bin',string_rusty_root);

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
% First load rib80sc molecule on x_u grid. ;
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_a_k_p_quad_.mat',dir_trunk);
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
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_a_k_p_quad_',dir_trunk);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_b_x_u_reco__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
a_x_u_reco__ = zeros(n_X_u,n_k_p_r);
eta = pi/k_p_r_max; 
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+1+nk_p_r)-1;
tmp_n_k_all = numel(1+tmp_index_);
tmp_t = tic;
a_x_u_reco__(:,1+nk_p_r) = nufft3d3(tmp_n_k_all,2*pi*k_c_0_all_(1+tmp_index_)*eta,2*pi*k_c_1_all_(1+tmp_index_)*eta,2*pi*k_c_2_all_(1+tmp_index_)*eta,a_k_p_quad_(1+tmp_index_).*(2*pi)^3.*weight_3d_k_all_(1+tmp_index_),+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_b_x_u_reco__',dir_trunk);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_a_k_Y_quad_.mat',dir_trunk);
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
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
tmp_t = tic;
[a_k_p_reco_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_quad_);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_a_k_Y_quad_A',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_a_k_Y_quad_',dir_trunk);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_S_k_p__.mat',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_S_k_p__',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_S_k_p__A',dir_trunk);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_CTF_k_p__.mat',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_k_p__',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_k_p_r_xxxx__',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_euler_angle_marina_',dir_trunk);
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
% Now look at some of the experimental images associated with rib80sc. ;
% Note that these will *NOT* be corrected/centered according to the (presumed) displacements. ;
% Moreover, to save space, we will not store the M_x_c___. ;
%%%%%%%%;
fname_M_x_c_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_M_x_c___.mat',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_euler_angle_marina_',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_delta_read_',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_M_x_c___center_',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_M_x_c___sample',dir_trunk);
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
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_M_k_p___.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
fname_image_mda = sprintf('%s/images_mda',dir_data);
M_x_c___ = MDA_read_r8(fname_image_mda);
x_c_0_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
x_c_1_ = -half_diameter_x_c + transpose([0:n_x_u-1]/n_x_u)*diameter_x_c;
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_);
M_abs_x_c_0_avg_ = zeros(n_image_sub,1);
M_abs_x_c_1_avg_ = zeros(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
M_abs_x_c_ = abs(squeeze(M_x_c___(:,:,1+nimage_sub)));
M_abs_avg = mean(M_abs_x_c_,'all');
M_abs_x_c_0_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_0__,'all');
M_abs_x_c_1_avg = mean(M_abs_x_c_/M_abs_avg.*x_c_1__,'all');
M_abs_x_c_0_avg_(1+nimage_sub) = M_abs_x_c_0_avg;
M_abs_x_c_1_avg_(1+nimage_sub) = M_abs_x_c_1_avg;
clear M_abs_x_c_;
end;%for nimage_sub=0:n_image_sub-1;
%%%%%%%%;
O_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds not centered images (i.e., original). ;
M_k_p__ = zeros(n_w_sum,n_image_sub); %<-- holds yes centered images (except, actually not centered). ;
M_uni_k_p__ = zeros(n_w_max*n_k_p_r,n_image_sub); %<-- holds yes centered images (except, actually not centered) on uniform grid. ;
for nimage_sub=0:n_image_sub-1;
if (mod(nimage_sub,100)==0); disp(sprintf(' %% nimage_sub %d/%d',nimage_sub,n_image_sub)); end;
tmp_M_x_c_ = squeeze(M_x_c___(:,:,1+nimage_sub));
%%%%%%%%;
% First convert tmp_M_x_c_ on the (nonuniform) polar grid used for the templates. ;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_) ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
% Now do *NOT* translate according to delta_read__. ;
% tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
O_k_p__(:,1+nimage_sub) = tmp_O_k_p_; %<-- not centered (i.e., original). ;
M_k_p__(:,1+nimage_sub) = tmp_M_k_p_; %<-- yes centered (according to ). ;
%%%%%%%%;
% Second convert tmp_M_x_c_ on a uniform polar grid to create principled-images. ;
%%%%%%%%;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
% Now do *NOT* translate according to delta_read__. ;
% tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
tmp_M_k_p__ = reshape(tmp_M_k_p_,n_w_max,n_k_p_r);
M_uni_k_p__(:,1+nimage_sub) = tmp_M_k_p__(:);
end;%for nimage_sub=0:n_image_sub-1;
delta_read_plus_M_abs_x_c_0_avg_ = delta_read_x_(1:n_image_sub) + M_abs_x_c_0_avg_;
delta_read_plus_M_abs_x_c_1_avg_ = delta_read_y_(1:n_image_sub) + M_abs_x_c_1_avg_;
save(fname_mat ...
     ,'M_k_p__' ...
     ,'M_uni_k_p__' ...
     ,'x_c_0_','x_c_1_' ...
     ,'M_abs_x_c_0_avg_','M_abs_x_c_1_avg_','delta_read_plus_M_abs_x_c_0_avg_','delta_read_plus_M_abs_x_c_1_avg_' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));

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
tmp_image_delta_x = +1.0*delta_read_plus_M_abs_x_c_0_avg_(1+nM);
tmp_image_delta_y = +1.0*delta_read_plus_M_abs_x_c_1_avg_(1+nM);
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
% First define delta_sigma (for translations). ;
%%%%%%%%;
delta_sigma = 1.0 * std([delta_read_plus_M_abs_x_c_0_avg_;delta_read_plus_M_abs_x_c_1_avg_]); %<-- no reduction. ;
%%%%%%%%;
% Then set up the cost matrix. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_X__.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
n_M = n_image_sub;
[X_3d_xcor_d0__,X_3d_xcor_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xcor__);
[X_3d_xavg_d0__,X_3d_xavg_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_3d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xavg__); 
[X_2d_xcor_d0__,X_2d_xcor_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xcor__);
[X_2d_xavg_d0__,X_2d_xavg_d0_weight_r_] = principled_marching_cost_matrix_3(n_k_p_r,weight_2d_k_p_r_,l_max_max,a_k_Y_quad__,CTF_k_p_r_xavg__); 
[X_2d_xcor_d1__,X_2d_xcor_d1_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,CTF_k_p_r_xcor__,delta_sigma);
[X_2d_xavg_d1__,X_2d_xavg_d1_weight_r_] = principled_marching_cost_matrix_5(n_k_p_r,k_p_r_,weight_2d_k_p_r_,l_max_,[],[],a_k_Y_quad_,CTF_k_p_r_xavg__,delta_sigma);
[X_2d_Suni_d0__,X_2d_Suni_d0_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_S,sparse(1:n_w_sum,1:n_w_sum,CTF_avg_k_p_,n_w_sum,n_w_sum)*S_k_p__);
[X_2d_Memp_d1__,X_2d_Memp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,M_k_p__);
%%%%%%%%;
X_TM_ = zeros(n_M,1);
T_k_p__ = zeros(n_w_sum,n_M);
flag_plot=0;
%%%%%%%%;
n_S = n_viewing_all;
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
for nM=0:n_M-1;
tmp_euler_polar_a = +euler_angle_marina_(1,1+nM);
tmp_euler_azimu_b = +euler_angle_marina_(2,1+nM);
tmp_euler_gamma_z = -euler_angle_marina_(3,1+nM);
tmp_image_delta_x = +1.0*delta_read_plus_M_abs_x_c_0_avg_(1+nM);
tmp_image_delta_y = +1.0*delta_read_plus_M_abs_x_c_1_avg_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
S_k_p_ = S_k_p__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_idx_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_MM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,M_k_p_,M_k_p_);
tmp_TM = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,M_k_p_);
X_TM_(1+nM) = real(tmp_TM)/sqrt(tmp_TT*tmp_MM);
T_k_p__(:,1+nM) = T_k_p_;
if flag_plot;
figure(1); clf;
figbeach();
subplot(2,2,1); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(M(k))');
subplot(2,2,2); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(M(k))');
subplot(2,2,3); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('real(T(k))');
subplot(2,2,4); imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(T_k_p_),[],colormap_beach());
set(gca,'XTick',[],'YTick',[]); axis image; title('imag(T(k))');
figbig;
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
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
% Now plot correlation X_TM_ vs viewing-angle and translation. ;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_X_TM_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
markersize_use = 32;
e_c2d__ = colormap_polar_a_azimu_b_2d(+euler_angle_marina_(1,1:n_M),+euler_angle_marina_(2,1:n_M),0.35);
d_c2d__ = colormap_gaussian_2d(delta_read_plus_M_abs_x_c_0_avg_,delta_read_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
nc_beach_ = max(0,min(n_c_beach-1,floor(n_c_beach*X_TM_)));
n_h = 64; lh_lim_ = [0,10];
figure(3);
figbeach();
figbig;
subplot(2,2,1);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),markersize_use,e_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view');
subplot(2,2,2);
hold on;
patch(2*pi*[+0;+1;+1;+0],1*pi*[+0;+0;+1;+1],0.65*[1,1,1],'EdgeColor','none');
scatter(0*pi+euler_angle_marina_(2,1:n_M),1*pi-euler_angle_marina_(1,1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
hold off;
axisnotick;
axis([0,2*pi,0,1*pi]); axis equal; xlabel('azimu_b'); ylabel('polar_a'); title('true view (correlation)');
subplot(2,2,3);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_M),delta_read_plus_M_abs_x_c_1_avg_(1:n_M),markersize_use,d_c2d__(1:n_M,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
subplot(2,2,4);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_M),delta_read_plus_M_abs_x_c_1_avg_(1:n_M),markersize_use,c_beach__(1+nc_beach_,:),'filled');
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
total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__ = zeros(n_UX_rank,n_UX_rank);
total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__ = zeros(n_UX_rank,n_UX_rank);
for tmp_r0=0:n_UX_rank-1;
for tmp_r1=0:n_UX_rank-1;
total_variance_X_2d_Memp_d1_from_X_2d_xcor_d0__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_Memp_d1__(:,1:1+tmp_r1)) * UX_2d_xcor_d0__(:,1:1+tmp_r0) * diag(SX_2d_xcor_d0_(1:1+tmp_r0)) * transpose(VX_2d_xcor_d0__(:,1:1+tmp_r0)) * VX_2d_Memp_d1__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_xcor_d0_(1:1+tmp_r0));
total_variance_X_2d_xcor_d0_from_X_2d_Memp_d1__(1+tmp_r1,1+tmp_r0) = sum(svds( transpose(UX_2d_xcor_d0__(:,1:1+tmp_r1)) * UX_2d_Memp_d1__(:,1:1+tmp_r0) * diag(SX_2d_Memp_d1_(1:1+tmp_r0)) * transpose(VX_2d_Memp_d1__(:,1:1+tmp_r0)) * VX_2d_xcor_d0__(:,1:1+tmp_r1) , tmp_r1 ))/sum(SX_2d_Memp_d1_(1:1+tmp_r0));
end;%for tmp_r1=0:n_UX_rank-1;
end;%for tmp_r0=0:n_UX_rank-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_X_2d_xcor_d0_vs_2d_Memp_d1_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf; hold on;
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

%%%%%%%%;
% Now visualize principal-volumes. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_vis_UX_.mat',dir_trunk);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_vis_UX_2d_xcor_d0_FIGA',dir_trunk);
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
figure(1);clf;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_vis_UX_2d_xcor_d1_FIGA',dir_trunk);
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
figure(1);clf;
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_vis_UX_2d_Memp_d1_FIGA',dir_trunk);
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
figure(1);clf;
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
%str_cost = sprintf('2d_xcor_d%.4d',floor(1000*delta_sigma));
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_X_%s__.mat',dir_trunk,str_cost);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now form a_CTF_UX_Y_quad__ ;
%%%%%%%%;
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
     ,'delta_sigma','n_UX_rank' ...
     ,'str_cost' ...
     ,'X__','X_weight_r_' ...
     ,'UX__','SX__','SX_','VX__','a_CTF_UX_Y_quad__' ...
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_X_%s_A',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_X_%s_B',dir_trunk,str_cost);
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

flag_check=0;
if flag_check;
%%%%%%%%;
% Now generate principled-[templates*CTF_avg]. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_CTF_UX_%s_S_k_p___.mat',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_UX_%s_S_k_p___',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_UX_%s_S_k_q___',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_UX_%s_S_k_q___spectrum',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_CTF_UX_%s_S_k_p___A',dir_trunk,str_cost);
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
end;%if flag_check;

%%%%%%%%;
% Now generate principled-images. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_UX_%s_M_k_p___.mat',dir_trunk,str_cost);
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
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
% Now do *NOT* translate according to delta_read__. ;
% tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
O_k_p__(:,1+nimage_sub) = tmp_O_k_p_; %<-- not centered (i.e., original). ;
M_k_p__(:,1+nimage_sub) = tmp_M_k_p_; %<-- yes centered (except, actually not centered). ;
%%%%%%%%;
% Second convert tmp_M_x_c_ on a uniform polar grid to create principled-images. ;
%%%%%%%%;
tmp_O_k_p_ = interp_x_c_to_k_p_nufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,tmp_M_x_c_,n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1)) ;
% Now *DO* translate according to M_abs_x_c_0_avg_ and M_abs_x_c_1_avg_. ;
tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,-1*M_abs_x_c_0_avg_(1+nimage_sub),-1*M_abs_x_c_1_avg_(1+nimage_sub));
% Now do *NOT* translate according to delta_read__. ;
% tmp_M_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_max*ones(n_k_p_r,1),n_w_max*n_k_p_r,tmp_O_k_p_,+0*delta_read_x_(1+nimage_sub),+0*delta_read_y_(1+nimage_sub));
tmp_M_k_p__ = reshape(tmp_M_k_p_,n_w_max,n_k_p_r);
M_uni_k_p__(:,1+nimage_sub) = tmp_M_k_p__(:);
for nUX_rank=0:n_UX_rank-1;
tmp_UX_M_k_p_ = zeros(n_w_max,1);
for nk_p_r=0:n_k_p_r-1;
nctf = CTF_idx_(1+nimage_sub);
tmp_UX_M_k_p_ = tmp_UX_M_k_p_ + UX__(1+nk_p_r,1+nUX_rank)*X_weight_r_(1+nk_p_r)*tmp_M_k_p__(:,1+nk_p_r);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_UX_%s_M_k_p___',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_UX_%s_M_k_q___',dir_trunk,str_cost);
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
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_UX_%s_M_k_q___spectrum',dir_trunk,str_cost);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_UX_%s_M_k_p___.mat',dir_trunk,str_cost);
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
tmp_pm_weight_3d_k_p_r_ = ones(tmp_pm_n_k_p_r,1);
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
,tmp_pm_weight_3d_k_p_r_ ...
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
,tmp_pm_weight_3d_k_p_r_ ...
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

fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_c_k_Y_.mat',dir_trunk);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
tmp_n_order = 5; tmp_n_M = n_image_sub;
tmp_euler_polar_a_ = +euler_angle_marina_(1,1+(0:n_image_sub-1));
tmp_euler_azimu_b_ = +euler_angle_marina_(2,1+(0:n_image_sub-1));
tmp_euler_gamma_z_ = -euler_angle_marina_(3,1+(0:n_image_sub-1));
tmp_image_delta_x_ = +1.0*delta_read_plus_M_abs_x_c_0_avg_(1+(0:n_image_sub-1));
tmp_image_delta_y_ = +1.0*delta_read_plus_M_abs_x_c_1_avg_(1+(0:n_image_sub-1));
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
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_c_x_u_',dir_trunk);
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_d0_CTF_UXI_S_k_p__.mat',dir_trunk);
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
b_k_Y_quad_ = rotate_spharm_to_spharm_2(0,[],n_k_p_r,k_p_r_,l_max_,a_k_Y_quad_,[0;+pi/2;0]);
%%%%%%%%;
n_UXI_rank = pm_n_k_p_r;
UXI__=eye(pm_n_k_p_r);
a_CTF_UXI_Y_quad__ = zeros(pm_n_lm_max,n_UXI_rank);
b_CTF_UXI_Y_quad__ = zeros(pm_n_lm_max,n_UXI_rank);
for nUXI_rank=0:n_UXI_rank-1;
for nk_p_r=0:n_k_p_r-1;
tmp_l_max = l_max_(1+nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:tmp_n_lm-1);
a_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) = a_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) + UXI__(1+nk_p_r,1+nUXI_rank)*X_weight_r_(1+nk_p_r)*a_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that a_CTF_UXI_Y_quad_ will be used alone. ;
b_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) = b_CTF_UXI_Y_quad__(1:tmp_n_lm,1+nUXI_rank) + UXI__(1+nk_p_r,1+nUXI_rank)*X_weight_r_(1+nk_p_r)*b_k_Y_quad_(1+tmp_index_)*CTF_avg_k_p_r_(1+nk_p_r); %<-- use average CTF here, under the assumption that b_CTF_UXI_Y_quad_ will be used alone. ;
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nUXI_rank=0:n_UXI_rank-1;
%%%%%%%%;
pm_template_k_eq_d = k_eq_d;
pm_viewing_k_eq_d = k_eq_d;
[ ...
 d0_a_CTF_UXI_S_k_p__ ...
,da_a_CTF_UXI_S_k_p__ ...
,db_a_CTF_UXI_S_k_p__ ...
,dc_a_CTF_UXI_S_k_p__ ...
,dt_a_CTF_UXI_S_k_p__ ...
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
[ ...
 d0_b_CTF_UXI_S_k_p__ ...
,da_b_CTF_UXI_S_k_p__ ...
,db_b_CTF_UXI_S_k_p__ ...
,dc_b_CTF_UXI_S_k_p__ ...
,dt_b_CTF_UXI_S_k_p__ ...
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
,~ ...
] = ...
get_dtemplate_0( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,[] ...
,pm_l_max_ ...
,b_CTF_UXI_Y_quad__(:) ...
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
l2_a_signal = sqrt(mean(abs(d0_a_CTF_UXI_S_k_p__).^2 * reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]) / (4*pi))) ;
l2_b_signal = sqrt(mean(abs(d0_b_CTF_UXI_S_k_p__).^2 * reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1]) / (4*pi))) ;
l2_signal = 0.5*(l2_a_signal + l2_b_signal);
%%%%%%%%;
save(fname_mat ...
     ,'l2_a_signal','l2_b_signal','l2_signal' ...
     ,'pm_template_k_eq_d' ...
     ,'pm_viewing_k_eq_d' ...
     ,'d0_a_CTF_UXI_S_k_p__' ...
     ,'da_a_CTF_UXI_S_k_p__' ...
     ,'db_a_CTF_UXI_S_k_p__' ...
     ,'dc_a_CTF_UXI_S_k_p__' ...
     ,'dt_a_CTF_UXI_S_k_p__' ...
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
     ,'d0_b_CTF_UXI_S_k_p__' ...
     ,'da_b_CTF_UXI_S_k_p__' ...
     ,'db_b_CTF_UXI_S_k_p__' ...
     ,'dc_b_CTF_UXI_S_k_p__' ...
     ,'dt_b_CTF_UXI_S_k_p__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_d0_CTF_UXI_S_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
nviewing_all = max(0,min(pm_n_viewing_all-1,round(pm_n_viewing_all*nplot/n_plot)));
subplot(6,5,1+5*nplot+0);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(d0_a_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('d0_a_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+1);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(da_a_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('da_a_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+2);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(db_a_CTF_UXI_S_k_p__(:,1+nviewing_all)/sin(pm_viewing_polar_a_all_(1+nviewing_all))),[],colormap_beach());
axis image; axisnotick; title(sprintf('db_a_CTF_UXI_S_k_p__(:,1+nv) / sin(polar_a) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+3);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(dc_a_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('dc_a_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+4);
imagesc_p(pm_n_k_p_r,1:pm_n_k_p_r,pm_n_w_,sum(pm_n_w_),real(dt_a_CTF_UXI_S_k_p__(:,1+nviewing_all)),[],colormap_beach());
axis image; axisnotick; title(sprintf('dt_a_CTF_UXI_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
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
%%%%%%%%;
% Now calculate the typical variance induced by the gradient. ;
%%%%%%%%;
process_UXI_a_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_a_CTF_UXI_S_k_p__ ...
,da_a_CTF_UXI_S_k_p__ ...
,db_a_CTF_UXI_S_k_p__ ...
,dc_a_CTF_UXI_S_k_p__ ...
,dt_a_CTF_UXI_S_k_p__ ...
);
process_UXI_b_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_b_CTF_UXI_S_k_p__ ...
,da_b_CTF_UXI_S_k_p__ ...
,db_b_CTF_UXI_S_k_p__ ...
,dc_b_CTF_UXI_S_k_p__ ...
,dt_b_CTF_UXI_S_k_p__ ...
);
%%%%%%%%;
pm_viewing_a_k_c_0_ = zeros(n_viewing_all,1);
pm_viewing_a_k_c_1_ = zeros(n_viewing_all,1);
pm_viewing_a_k_c_2_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_0_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_1_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_2_ = zeros(n_viewing_all,1);
for nviewing_all=0:pm_n_viewing_all-1;
pm_viewing_azimu_b = pm_viewing_azimu_b_all_(1+nviewing_all);
pm_viewing_polar_a = pm_viewing_polar_a_all_(1+nviewing_all);
pm_viewing_a_k_c_0_(1+nviewing_all) = cos(pm_viewing_azimu_b)*sin(pm_viewing_polar_a);
pm_viewing_a_k_c_1_(1+nviewing_all) = sin(pm_viewing_azimu_b)*sin(pm_viewing_polar_a);
pm_viewing_a_k_c_2_(1+nviewing_all) = cos(pm_viewing_polar_a);
tmp_a_ = [ cos(pm_viewing_azimu_b)*sin(pm_viewing_polar_a) ; sin(pm_viewing_azimu_b)*sin(pm_viewing_polar_a) ; cos(pm_viewing_polar_a) ];
tmp_b_ = [ 0 , 0 , 1 ; 0 , 1 , 0 ; 1 , 0 , 0 ] * tmp_a_; %<-- rotation of viewing-angle by -pi/2 about positive y-axis. ;
pm_viewing_b_k_c_0_(1+nviewing_all) = tmp_b_(1+0);
pm_viewing_b_k_c_1_(1+nviewing_all) = tmp_b_(1+1);
pm_viewing_b_k_c_2_(1+nviewing_all) = tmp_b_(1+2);
end;%for nviewing_all=0:pm_n_viewing_all-1;
index_knn_b_from_a_ = knnsearch( [ pm_viewing_b_k_c_0_ , pm_viewing_b_k_c_1_ , pm_viewing_b_k_c_2_ ] , [ pm_viewing_a_k_c_0_ , pm_viewing_a_k_c_1_ , pm_viewing_a_k_c_2_ ] ) - 1;
%%%%%%%%;
l2_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_l2_snr=numel(l2_snr_);
eta_equi_upb___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
eta_equi_lob___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
l2_sigma = 0; if l2_snr>0; l2_sigma = l2_signal/l2_snr; end;
for nk_p_r=0:pm_n_k_p_r-1;
for nviewing_all=0:pm_n_viewing_all-1;
tmp_a = process_UXI_a_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_all);
tmp_b = -1.0d0;
tmp_c = process_UXI_a_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_all)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb_a = +Inf; eta_equi_lob_a = +Inf; end;
if (tmp_d>=0); eta_equi_upb_a = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob_a = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
nviewing_b_from_a = index_knn_b_from_a_(1+nviewing_all);
tmp_a = process_UXI_b_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_b_from_a);
tmp_b = -1.0d0;
tmp_c = process_UXI_b_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_b_from_a)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb_b = +Inf; eta_equi_lob_b = +Inf; end;
if (tmp_d>=0); eta_equi_upb_b = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob_b = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
eta_equi_upb = min(eta_equi_upb_a,eta_equi_upb_b*pm_viewing_weight_all_(1+nviewing_b_from_a)/pm_viewing_weight_all_(1+nviewing_all));
eta_equi_lob = min(eta_equi_lob_a,eta_equi_lob_b*pm_viewing_weight_all_(1+nviewing_b_from_a)/pm_viewing_weight_all_(1+nviewing_all));
eta_equi_upb___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_upb;
eta_equi_lob___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_lob;
end;%for nviewing_all=0:pm_n_viewing_all-1;
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
eta_equi_upb__ = reshape(reshape(eta_equi_upb___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
eta_equi_lob__ = reshape(reshape(eta_equi_lob___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_d0_CTF_UXI_S_grad_FIGA',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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
fname_mat = sprintf('%s_mat/test_principled_marching_rib80sc_0_d0_CTF_UX_S_k_p__.mat',dir_trunk);
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
b_CTF_UX_Y_quad__ = zeros(n_lm_max,n_UX_rank);
for nUX_rank=0:n_UX_rank-1;
b_CTF_UX_Y_quad__(:,1+nUX_rank) = rotate_spharm_to_spharm_2(0,[],1,1,l_max_max,a_CTF_UX_Y_quad__(:,1+nUX_rank),[0;+pi/2;0]);
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
pm_template_k_eq_d = k_eq_d;
pm_viewing_k_eq_d = k_eq_d;
[ ...
 d0_a_CTF_UX_S_k_p__ ...
,da_a_CTF_UX_S_k_p__ ...
,db_a_CTF_UX_S_k_p__ ...
,dc_a_CTF_UX_S_k_p__ ...
,dt_a_CTF_UX_S_k_p__ ...
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
[ ...
 d0_b_CTF_UX_S_k_p__ ...
,da_b_CTF_UX_S_k_p__ ...
,db_b_CTF_UX_S_k_p__ ...
,dc_b_CTF_UX_S_k_p__ ...
,dt_b_CTF_UX_S_k_p__ ...
] = ...
get_dtemplate_0( ...
 verbose ...
,pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_k_p_r_max ...
,[] ...
,pm_l_max_ ...
,b_CTF_UX_Y_quad__(:) ...
,pm_viewing_k_eq_d ...
,0*pm_template_k_eq_d ...
,pm_n_w_ ...
);
save(fname_mat ...
     ,'pm_template_k_eq_d' ...
     ,'pm_viewing_k_eq_d' ...
     ,'d0_a_CTF_UX_S_k_p__' ...
     ,'da_a_CTF_UX_S_k_p__' ...
     ,'db_a_CTF_UX_S_k_p__' ...
     ,'dc_a_CTF_UX_S_k_p__' ...
     ,'dt_a_CTF_UX_S_k_p__' ...
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
     ,'d0_b_CTF_UX_S_k_p__' ...
     ,'da_b_CTF_UX_S_k_p__' ...
     ,'db_b_CTF_UX_S_k_p__' ...
     ,'dc_b_CTF_UX_S_k_p__' ...
     ,'dt_b_CTF_UX_S_k_p__' ...
     );
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_d0_CTF_UX_S_k_p__',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
n_plot=6;
for nplot=0:n_plot-1;
nviewing_all = max(0,min(pm_n_viewing_all-1,round(pm_n_viewing_all*nplot/n_plot)));
subplot(6,5,1+5*nplot+0);
imagesc(real(reshape(d0_a_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('d0_a_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+1);
imagesc(real(reshape(da_a_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('da_a_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+2);
imagesc(real(reshape(db_a_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('db_a_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+3);
imagesc(real(reshape(dc_a_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('dc_a_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
subplot(6,5,1+5*nplot+4);
imagesc(real(reshape(dt_a_CTF_UX_S_k_p__(:,1+nviewing_all),[pm_n_w_max,pm_n_k_p_r]))); colormap(colormap_beach());
axisnotick; title(sprintf('dt_a_CTF_UX_S_k_p__(:,1+nv) %d',nviewing_all),'Interpreter','none');
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
process_UX_a_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_a_CTF_UX_S_k_p__ ...
,da_a_CTF_UX_S_k_p__ ...
,db_a_CTF_UX_S_k_p__ ...
,dc_a_CTF_UX_S_k_p__ ...
,dt_a_CTF_UX_S_k_p__ ...
);
process_UX_b_dtemplate = ...
process_dtemplate_1( ...
 pm_n_k_p_r ...
,pm_k_p_r_ ...
,pm_n_viewing_all ...
,pm_n_w_ ...
,pm_viewing_azimu_b_all_ ...
,pm_viewing_polar_a_all_ ...
,pm_viewing_weight_all_ ...
,d0_b_CTF_UX_S_k_p__ ...
,da_b_CTF_UX_S_k_p__ ...
,db_b_CTF_UX_S_k_p__ ...
,dc_b_CTF_UX_S_k_p__ ...
,dt_b_CTF_UX_S_k_p__ ...
);
%%%%%%%%;
pm_viewing_a_k_c_0_ = zeros(n_viewing_all,1);
pm_viewing_a_k_c_1_ = zeros(n_viewing_all,1);
pm_viewing_a_k_c_2_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_0_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_1_ = zeros(n_viewing_all,1);
pm_viewing_b_k_c_2_ = zeros(n_viewing_all,1);
for nviewing_all=0:pm_n_viewing_all-1;
pm_viewing_azimu_b = pm_viewing_azimu_b_all_(1+nviewing_all);
pm_viewing_polar_a = pm_viewing_polar_a_all_(1+nviewing_all);
pm_viewing_a_k_c_0_(1+nviewing_all) = cos(pm_viewing_azimu_b)*sin(pm_viewing_polar_a);
pm_viewing_a_k_c_1_(1+nviewing_all) = sin(pm_viewing_azimu_b)*sin(pm_viewing_polar_a);
pm_viewing_a_k_c_2_(1+nviewing_all) = cos(pm_viewing_polar_a);
tmp_a_ = [ cos(pm_viewing_azimu_b)*sin(pm_viewing_polar_a) ; sin(pm_viewing_azimu_b)*sin(pm_viewing_polar_a) ; cos(pm_viewing_polar_a) ];
tmp_b_ = [ 0 , 0 , 1 ; 0 , 1 , 0 ; 1 , 0 , 0 ] * tmp_a_; %<-- rotation of viewing-angle by -pi/2 about positive y-axis. ;
pm_viewing_b_k_c_0_(1+nviewing_all) = tmp_b_(1+0);
pm_viewing_b_k_c_1_(1+nviewing_all) = tmp_b_(1+1);
pm_viewing_b_k_c_2_(1+nviewing_all) = tmp_b_(1+2);
end;%for nviewing_all=0:pm_n_viewing_all-1;
index_knn_b_from_a_ = knnsearch( [ pm_viewing_b_k_c_0_ , pm_viewing_b_k_c_1_ , pm_viewing_b_k_c_2_ ] , [ pm_viewing_a_k_c_0_ , pm_viewing_a_k_c_1_ , pm_viewing_a_k_c_2_ ] ) - 1;
%%%%%%%%;
l2_snr_ = [0 , 2.^[-2:-0.5:-4]]; n_l2_snr=numel(l2_snr_);
eta_equi_upb___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
eta_equi_lob___ = zeros(n_l2_snr,pm_n_k_p_r,pm_n_viewing_all);
for nl2_snr=0:n_l2_snr-1;
l2_snr = l2_snr_(1+nl2_snr);
l2_sigma = 0; if l2_snr>0; l2_sigma = l2_signal/l2_snr; end;
for nk_p_r=0:pm_n_k_p_r-1;
for nviewing_all=0:pm_n_viewing_all-1;
tmp_a = process_UX_a_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_all);
tmp_b = -1.0d0;
tmp_c = process_UX_a_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_all)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb_a = +Inf; eta_equi_lob_a = +Inf; end;
if (tmp_d>=0); eta_equi_upb_a = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob_a = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
nviewing_b_from_a = index_knn_b_from_a_(1+nviewing_all);
tmp_a = process_UX_b_dtemplate.grad_abc_inv_grad_t_l2_csum__(1+nk_p_r,1+nviewing_b_from_a);
tmp_b = -1.0d0;
tmp_c = process_UX_b_dtemplate.grad_abc_inv_l2_csum__(1+nk_p_r,1+nviewing_b_from_a)*l2_sigma;
tmp_d = tmp_b^2 - 4*tmp_a*tmp_c;
if (tmp_d< 0);  eta_equi_upb_b = +Inf; eta_equi_lob_b = +Inf; end;
if (tmp_d>=0); eta_equi_upb_b = (-tmp_b + sqrt(tmp_d))/(2*tmp_a); eta_equi_lob_b = (-tmp_b - sqrt(tmp_d))/(2*tmp_a); end;
eta_equi_upb = min(eta_equi_upb_a,eta_equi_upb_b*pm_viewing_weight_all_(1+nviewing_b_from_a)/pm_viewing_weight_all_(1+nviewing_all));
eta_equi_lob = min(eta_equi_lob_a,eta_equi_lob_b*pm_viewing_weight_all_(1+nviewing_b_from_a)/pm_viewing_weight_all_(1+nviewing_all));
eta_equi_upb___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_upb;
eta_equi_lob___(1+nl2_snr,1+nk_p_r,1+nviewing_all) = eta_equi_lob;
end;%for nviewing_all=0:pm_n_viewing_all-1;
end;%for nk_p_r=0:pm_n_k_p_r-1;
end;%for nl2_snr=0:n_l2_snr-1;
eta_equi_upb__ = reshape(reshape(eta_equi_upb___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
eta_equi_lob__ = reshape(reshape(eta_equi_lob___,[n_l2_snr*pm_n_k_p_r,pm_n_viewing_all])*reshape(pm_viewing_weight_all_,[pm_n_viewing_all,1])/(4*pi),[n_l2_snr,pm_n_k_p_r]);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/test_principled_marching_rib80sc_0_d0_CTF_UX_S_grad_FIGA',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
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

flag_compute = 0 & strcmp(platform,'access1');
if flag_compute;
test_mp_init_0; %<-- test initialization using pca. ;
end;%if flag_compute;

flag_compute = 0 & strcmp(platform,'access1');
if flag_compute;
test_am_init_2; %<-- test initialization using displacements. ;
end;%if flag_compute;

flag_compute = 0 & strcmp(platform,'access1');
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_UX_ = sum(UX__(:,1:dat_n_UX_rank).^2,2);
tmp_UX_csum_ = cumsum(tmp_UX_,'reverse')/sum(tmp_UX_);
index_k_p_r_use = min(efind(tmp_UX_csum_<1e-2));
tmp_n_k_p_r = 1+index_k_p_r_use;
tmp_n_w_uni_ = n_w_max*ones(tmp_n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(tmp_n_w_uni_),:);
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
dat_f_rand = 0.05; dat_flag_plot=0; dat_order_limit_MS = -1; dat_flag_mp_init = 0;
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
,tmp_n_w_uni_ ...
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
,dat_order_limit_MS ...
,dat_flag_mp_init ...
);
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now check to see if the translations discovered correlate with the true translations. ;
%%%%%%%%;
flag_replot = 0 & strcmp(platform,'access1');
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
tmp_dir_jpg = sprintf('%s_jpg/dir_tpmutameux_%s',dir_trunk,dat_infix);
if (~exist(tmp_dir_jpg,'dir')); disp(sprintf(' %% %s not found, creating',tmp_dir_jpg)); mkdir(tmp_dir_jpg); end;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_fig = sprintf('%s/tpmutameux_%s_nUX%.3drng%.3d_FIGd',tmp_dir_jpg,dat_infix,dat_nUX_rank,dat_rseed);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); 
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1); markersize_use = 16;
%%%%%%%%;
subplot(5,7,2);
c2d__ = colormap_gaussian_2d(delta_read_plus_M_abs_x_c_0_avg_,delta_read_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_image_sub),delta_read_plus_M_abs_x_c_1_avg_(1:n_image_sub),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
nc = max(0,min(n_c-1,floor(n_c*(ndat_rseed + dat_nUX_rank*n_dat_rseed)/(n_dat_rseed*dat_n_UX_rank))));
tmp_flag_exist = 1; if (~exist(MS_fname_mat,'file') & ~exist(SM_fname_mat,'file')); tmp_flag_exist = 0; end;
%%%%%%%%;
if ( exist(MS_fname_mat,'file'));
tmp_MS_ = load(MS_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 0*dat_n_iteration,tmp_MS_.X_best_MS_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_MS00__ = [tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_MS00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('MS %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_MS_;
end;%if ( exist(MS_fname_mat,'file'));
%%%%%%%%;
if ( exist(SM_fname_mat,'file'));
tmp_SM_ = load(SM_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 1*dat_n_iteration,tmp_SM_.X_best_SM_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+dat_n_iteration+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_SM00__ = [tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_SM00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('SM %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_SM_;
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
figbig;
if tmp_flag_exist;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if tmp_flag_exist;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file')); 
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now run through each of the runs and update mat-file with the true correlation. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
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
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
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
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% experimental: Collect MS and SM output. ;
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('ut%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_MS____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
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
dat_X_true_MS___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_MS_(:);
flag_dat_exist_MS___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_MS____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_MS__(:,:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_SM_(:);
dat_X_true_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_SM_(:);
flag_dat_exist_SM___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_SM____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_SM__(:,:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_dat_exist_MS___,'all'),numel(flag_dat_exist_MS___)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_dat_exist_SM___,'all'),numel(flag_dat_exist_SM___)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now plot typical figures illustrating correlations as a function of iteration. ;
%%%%%%%%;
if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
flag_replot = 1;
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_best_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_best_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_true_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_true_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_true_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
end;%if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

flag_compute = 1 & strcmp(platform,'access1');
flag_skip_to_plot = 0; %<-- skip over main computation and plot figures (as best as possible). ;
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
% adding extra step where the initial MS phase limits the order of the test-function. ;
% adding extra step where initial data is determined by mp_init_0;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
order_limit_MS_ = [0,1,2,4,8]; n_order_limit_MS = numel(order_limit_MS_);
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mpolut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_UX_ = sum(UX__(:,1:dat_n_UX_rank).^2,2);
tmp_UX_csum_ = cumsum(tmp_UX_,'reverse')/sum(tmp_UX_);
index_k_p_r_use = min(efind(tmp_UX_csum_<1e-2));
tmp_n_k_p_r = 1+index_k_p_r_use;
tmp_n_w_uni_ = n_w_max*ones(tmp_n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(tmp_n_w_uni_),:);
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
dat_f_rand = 0.05; dat_flag_plot=0; dat_flag_mp_init=1;
if ~flag_skip_to_plot;
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
,tmp_n_w_uni_ ...
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
,dat_order_limit_MS ...
,dat_flag_mp_init ...
);
end;%if ~flag_skip_to_plot;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now check to see if the translations discovered correlate with the true translations. ;
%%%%%%%%;
flag_replot = 0 & strcmp(platform,'access1');
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mpolut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
tmp_dir_jpg = sprintf('%s_jpg/dir_tpmutameux_%s',dir_trunk,dat_infix);
if (~exist(tmp_dir_jpg,'dir')); disp(sprintf(' %% %s not found, creating',tmp_dir_jpg)); mkdir(tmp_dir_jpg); end;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_fig = sprintf('%s/tpmutameux_%s_nUX%.3drng%.3d_FIGd',tmp_dir_jpg,dat_infix,dat_nUX_rank,dat_rseed);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); 
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1); markersize_use = 16;
%%%%%%%%;
subplot(5,7,2);
c2d__ = colormap_gaussian_2d(delta_read_plus_M_abs_x_c_0_avg_,delta_read_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_image_sub),delta_read_plus_M_abs_x_c_1_avg_(1:n_image_sub),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
nc = max(0,min(n_c-1,floor(n_c*(ndat_rseed + dat_nUX_rank*n_dat_rseed)/(n_dat_rseed*dat_n_UX_rank))));
tmp_flag_exist = 1; if (~exist(MS_fname_mat,'file') & ~exist(SM_fname_mat,'file')); tmp_flag_exist = 0; end;
%%%%%%%%;
if ( exist(MS_fname_mat,'file'));
tmp_MS_ = load(MS_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 0*dat_n_iteration,tmp_MS_.X_best_MS_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_MS00__ = [tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_MS00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('MS %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_MS_;
end;%if ( exist(MS_fname_mat,'file'));
%%%%%%%%;
if ( exist(SM_fname_mat,'file'));
tmp_SM_ = load(SM_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 1*dat_n_iteration,tmp_SM_.X_best_SM_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+dat_n_iteration+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_SM00__ = [tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_SM00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('SM %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_SM_;
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
figbig;
if tmp_flag_exist;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if tmp_flag_exist;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file')); 
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now run through each of the runs and update mat-file with the true correlation. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
order_limit_MS_ = [0,1,2,4,8]; n_order_limit_MS = numel(order_limit_MS_);
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mpolut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
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
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutamf_1( ...
 sprintf('%s_mat',dir_trunk) ...
,sprintf('%s_jpg',dir_trunk) ...
,fname_0 ...
,fname_1 ...
,fname_2 ...
,n_image_sub ...
,M_k_p__ ...
,dat_nUX_rank ...
,dat_n_UX_rank ...
,dat_n_iteration ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
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
,CTF_k_p__ ...
,l_max_ ...
,UX__ ...
,X_weight_r_ ...
,a_CTF_UX_Y_quad__ ...
,[] ...
,a_k_Y_quad_ ...
,euler_angle_marina_ ...
,delta_read_x_ ...
,delta_read_y_ ...
,delta_sigma ...
);
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% experimental: Collect MS and SM output. ;
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mpolut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_MS____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
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
dat_X_true_MS___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_MS_(:);
flag_dat_exist_MS___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_MS____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_MS__(:,:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_SM_(:);
dat_X_true_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_SM_(:);
flag_dat_exist_SM___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_SM____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_SM__(:,:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_dat_exist_MS___,'all'),numel(flag_dat_exist_MS___)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_dat_exist_SM___,'all'),numel(flag_dat_exist_SM___)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now plot typical figures illustrating correlations as a function of iteration. ;
%%%%%%%%;
if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
flag_replot = 1;
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_best_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_best_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_true_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_true_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_true_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
end;%if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

disp('returning'); return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Appendix functions start here. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

flag_compute = 0 & strcmp(platform,'access1');
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
% adding extra step where initial data is determined by mp_init_0;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mput%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_UX_ = sum(UX__(:,1:dat_n_UX_rank).^2,2);
tmp_UX_csum_ = cumsum(tmp_UX_,'reverse')/sum(tmp_UX_);
index_k_p_r_use = min(efind(tmp_UX_csum_<1e-2));
tmp_n_k_p_r = 1+index_k_p_r_use;
tmp_n_w_uni_ = n_w_max*ones(tmp_n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(tmp_n_w_uni_),:);
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
,tmp_n_w_uni_ ...
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
,-1 ...
);
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now check to see if the translations discovered correlate with the true translations. ;
%%%%%%%%;
flag_replot = 0 & strcmp(platform,'access1');
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mput%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
tmp_dir_jpg = sprintf('%s_jpg/dir_tpmutameux_%s',dir_trunk,dat_infix);
if (~exist(tmp_dir_jpg,'dir')); disp(sprintf(' %% %s not found, creating',tmp_dir_jpg)); mkdir(tmp_dir_jpg); end;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_fig = sprintf('%s/tpmutameux_%s_nUX%.3drng%.3d_FIGd',tmp_dir_jpg,dat_infix,dat_nUX_rank,dat_rseed);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); 
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1); markersize_use = 16;
%%%%%%%%;
subplot(5,7,2);
c2d__ = colormap_gaussian_2d(delta_read_plus_M_abs_x_c_0_avg_,delta_read_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_image_sub),delta_read_plus_M_abs_x_c_1_avg_(1:n_image_sub),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
nc = max(0,min(n_c-1,floor(n_c*(ndat_rseed + dat_nUX_rank*n_dat_rseed)/(n_dat_rseed*dat_n_UX_rank))));
tmp_flag_exist = 1; if (~exist(MS_fname_mat,'file') & ~exist(SM_fname_mat,'file')); tmp_flag_exist = 0; end;
%%%%%%%%;
if ( exist(MS_fname_mat,'file'));
tmp_MS_ = load(MS_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 0*dat_n_iteration,tmp_MS_.X_best_MS_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_MS00__ = [tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_MS00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('MS %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_MS_;
end;%if ( exist(MS_fname_mat,'file'));
%%%%%%%%;
if ( exist(SM_fname_mat,'file'));
tmp_SM_ = load(SM_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 1*dat_n_iteration,tmp_SM_.X_best_SM_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+dat_n_iteration+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_SM00__ = [tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_SM00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('SM %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_SM_;
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
figbig;
if tmp_flag_exist;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if tmp_flag_exist;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file')); 
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now run through each of the runs and update mat-file with the true correlation. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mput%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% experimental: Collect MS and SM output. ;
%%%%%%%%;
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('mput%.4dle%dv%.3d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_MS____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
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
dat_X_true_MS___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_MS_(:);
flag_dat_exist_MS___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_MS____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_MS__(:,:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_SM_(:);
dat_X_true_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_SM_(:);
flag_dat_exist_SM___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_SM____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_SM__(:,:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_dat_exist_MS___,'all'),numel(flag_dat_exist_MS___)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_dat_exist_SM___,'all'),numel(flag_dat_exist_SM___)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now plot typical figures illustrating correlations as a function of iteration. ;
%%%%%%%%;
if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
flag_replot = 1;
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_best_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_best_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_true_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_true_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_true_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
end;%if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

flag_compute = 0 & strcmp(platform,'access1');
if flag_compute;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% check to see if translations can be estimated (and appropriately updated) using only a few principal-modes. ;
% adding extra step where the initial MS phase limits the order of the test-function. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
order_limit_MS_ = [0,1,2,4,8]; n_order_limit_MS = numel(order_limit_MS_);
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('olut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
tmp_UX_ = sum(UX__(:,1:dat_n_UX_rank).^2,2);
tmp_UX_csum_ = cumsum(tmp_UX_,'reverse')/sum(tmp_UX_);
index_k_p_r_use = min(efind(tmp_UX_csum_<1e-2));
tmp_n_k_p_r = 1+index_k_p_r_use;
tmp_n_w_uni_ = n_w_max*ones(tmp_n_k_p_r,1);
dat_M_k_p__ = M_uni_k_p__(1:sum(tmp_n_w_uni_),:);
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
,tmp_n_w_uni_ ...
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
,dat_order_limit_MS ...
);
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now check to see if the translations discovered correlate with the true translations. ;
%%%%%%%%;
flag_replot = 0 & strcmp(platform,'access1');
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('olut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
tmp_dir_jpg = sprintf('%s_jpg/dir_tpmutameux_%s',dir_trunk,dat_infix);
if (~exist(tmp_dir_jpg,'dir')); disp(sprintf(' %% %s not found, creating',tmp_dir_jpg)); mkdir(tmp_dir_jpg); end;
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_fig = sprintf('%s/tpmutameux_%s_nUX%.3drng%.3d_FIGd',tmp_dir_jpg,dat_infix,dat_nUX_rank,dat_rseed);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig)); 
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_pre = sprintf('%s_MS_%s',fname_0,fname_2);
MS_fname_mat = sprintf('%s.mat',MS_fname_pre);
SM_fname_pre = sprintf('%s_SM_%s',fname_0,fname_2);
SM_fname_mat = sprintf('%s.mat',SM_fname_pre);
figure(1);clf;
c_ = colormap_beach(); n_c = size(c_,1); markersize_use = 16;
%%%%%%%%;
subplot(5,7,2);
c2d__ = colormap_gaussian_2d(delta_read_plus_M_abs_x_c_0_avg_,delta_read_plus_M_abs_x_c_1_avg_,delta_sigma,0.35);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(delta_read_plus_M_abs_x_c_0_avg_(1:n_image_sub),delta_read_plus_M_abs_x_c_1_avg_(1:n_image_sub),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title('true delta');
%%%%%%%%;
nc = max(0,min(n_c-1,floor(n_c*(ndat_rseed + dat_nUX_rank*n_dat_rseed)/(n_dat_rseed*dat_n_UX_rank))));
tmp_flag_exist = 1; if (~exist(MS_fname_mat,'file') & ~exist(SM_fname_mat,'file')); tmp_flag_exist = 0; end;
%%%%%%%%;
if ( exist(MS_fname_mat,'file'));
tmp_MS_ = load(MS_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 0*dat_n_iteration,tmp_MS_.X_best_MS_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_MS00__ = [tmp_MS_.image_delta_x_MS__(:,1+niteration),tmp_MS_.image_delta_y_MS__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_MS00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('MS %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_MS_;
end;%if ( exist(MS_fname_mat,'file'));
%%%%%%%%;
if ( exist(SM_fname_mat,'file'));
tmp_SM_ = load(SM_fname_mat);
subplot(5,7,1);
hold on;
plot((0:dat_n_iteration-1) + 1*dat_n_iteration,tmp_SM_.X_best_SM_,'o-','Color',c_(1+nc,:),'LineWidth',2);
xlim([0,2*dat_n_iteration-1]);
ylim([-0.2,1.0]); grid on;
xlabel('iteration'); ylabel('correlation');
hold off;
for niteration=0:dat_n_iteration-1;
subplot(5,7,3+dat_n_iteration+niteration);
hold on;
patch([-1;+1;+1;-1],[-1;-1;+1;+1],'k');
scatter(tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration),markersize_use,c2d__(1:n_image_sub,:),'filled');
hold off;
axisnotick;
tmp_delta_true__ = [delta_read_plus_M_abs_x_c_0_avg_(1:dat_n_M),delta_read_plus_M_abs_x_c_1_avg_(1:dat_n_M)];
tmp_delta_SM00__ = [tmp_SM_.image_delta_x_SM__(:,1+niteration),tmp_SM_.image_delta_y_SM__(:,1+niteration)];
tmp_error_ = sqrt(sum((tmp_delta_true__ - tmp_delta_SM00__).^2,2));
axis([-1,+1,-1,+1]); axis square; xlabel('dx'); ylabel('dy'); title(sprintf('SM %d e %0.2f',niteration,mean(tmp_error_)));
end;%for niteration=0:dat_n_iteration-1;
clear tmp_SM_;
end;%if ( exist(SM_fname_mat,'file'));
%%%%%%%%;
figbig;
if tmp_flag_exist;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
end;%if tmp_flag_exist;
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file')); 
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now run through each of the runs and update mat-file with the true correlation. ;
%%%%%%%%;
delta_sigma_use = delta_sigma;
svd_eps = 1e-3;
dat_n_order = 5;
dat_n_M = n_image_sub;
dat_n_UX_rank = 8;
dat_n_iteration = 16;
dat_rseed_ = [0:2]; n_dat_rseed = numel(dat_rseed_);
%delta_r_max_factor_ = [0,0.125,0.25,0.5,1,sqrt(2),2]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0]; n_delta_r_max_factor = numel(delta_r_max_factor_);
%delta_r_max_factor_ = [0.125,0.25,0.50,0.75]; n_delta_r_max_factor = numel(delta_r_max_factor_);
delta_r_max_factor_ = [0.25]; n_delta_r_max_factor = numel(delta_r_max_factor_);
order_limit_MS_ = [0,1,2,4,8]; n_order_limit_MS = numel(order_limit_MS_);
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('olut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
for ndat_rseed=0:n_dat_rseed-1;
dat_rseed = dat_rseed_(1+ndat_rseed);
for dat_nUX_rank=0:dat_n_UX_rank-1;
fname_0 = sprintf('%s_mat/tpmutameux_%s_n%.3d',dir_trunk,dat_infix,dat_n_M);
fname_2 = sprintf('nUX%.3drng%.3d',dat_nUX_rank,dat_rseed);
MS_fname_mat = sprintf('%s_MS_%s.mat',fname_0,fname_2);
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(MS_fname_mat,'file'));
disp(sprintf(' %% %s found',MS_fname_mat));
fname_1 = 'MS';
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(MS_fname_mat,'file'));
if ( exist(SM_fname_mat,'file'));
disp(sprintf(' %% %s found',SM_fname_mat));
fname_1 = 'SM';
tpmutam_X_true_0( ...
 fname_0 ...
,fname_1 ...
,fname_2 ...
,dat_n_M ...
,M_k_p__ ...
,dat_n_order ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,n_w_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,CTF_idx_ ...
,CTF_k_p__ ...
,l_max_ ...
,a_k_Y_quad_ ...
);
end;%if ( exist(SM_fname_mat,'file'));
end;%for dat_nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% experimental: Collect MS and SM output. ;
%%%%%%%%;
for norder_limit_MS=0:n_order_limit_MS-1;
dat_order_limit_MS = order_limit_MS_(1+norder_limit_MS);
for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
delta_r_max_factor = delta_r_max_factor_(1+ndelta_r_max_factor);
delta_r_max_use = delta_r_max_factor * delta_sigma_use*sqrt(log(20^2));
delta_r_max_upb = 2.0 * delta_sigma_use*sqrt(log(20^2));
if (delta_r_max_factor<=0); n_delta_v_requested=1; end;
if (delta_r_max_factor> 0) & (delta_r_max_factor<=1); n_delta_v_requested=64; end;
if (delta_r_max_factor> 1) & (delta_r_max_factor<=2); n_delta_v_requested=128; end;
str_delta = sprintf('olut%.4dle%dv%.3dol%d',floor(1000*delta_r_max_use),round(-log(svd_eps)),n_delta_v_requested,dat_order_limit_MS);
dat_infix = sprintf('%s_%s',str_cost,str_delta);
dat_X_best_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_MS___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_M_loading_MS____ = zeros(3,dat_n_M,dat_n_UX_rank,n_dat_rseed);
dat_X_best_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
dat_X_true_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
flag_dat_exist_SM___ = zeros(dat_n_iteration,dat_n_UX_rank,n_dat_rseed);
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
dat_X_true_MS___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_MS_(:);
flag_dat_exist_MS___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_MS____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_MS__(:,:);
clear tmp_;
end;%if ( exist(MS_fname_mat,'file'));
SM_fname_mat = sprintf('%s_SM_%s.mat',fname_0,fname_2);
if ( exist(SM_fname_mat,'file'));
tmp_ = load(SM_fname_mat);
dat_X_best_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_best_SM_(:);
dat_X_true_SM___(:,1+nUX_rank,1+ndat_rseed) = tmp_.X_true_SM_(:);
flag_dat_exist_SM___(:,1+nUX_rank,1+ndat_rseed) = 1;
dat_M_loading_SM____(:,:,1+nUX_rank,1+ndat_rseed) = tmp_.dat_M_loading_SM__(:,:);
clear tmp_;
end;%if ( exist(SM_fname_mat,'file'));
end;%for nUX_rank=0:dat_n_UX_rank-1;
end;%for ndat_rseed=0:n_dat_rseed-1;
disp(sprintf(' %% MS: found %d/%d',sum(flag_dat_exist_MS___,'all'),numel(flag_dat_exist_MS___)));
disp(sprintf(' %% SM: found %d/%d',sum(flag_dat_exist_SM___,'all'),numel(flag_dat_exist_SM___)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%;
% Now plot typical figures illustrating correlations as a function of iteration. ;
%%%%%%%%;
if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
%%%%%%%%;
% Generate figures. ;
%%%%%%%%;
flag_replot = 1;
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_best_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_best_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_best_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
fname_fig = sprintf('%s_jpg/tpmutameux_%s_n%.3d_X_true_A',dir_trunk,dat_infix,dat_n_M);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
clf;
c_ = colormap_beach(); n_c = size(c_,1); %c_ = diag(0.85 + 0.15*transpose(linspace(-1,+1,n_c)).^2) * c_;
subplot(1,1,1);
hold on;
for nUX_rank=0:dat_n_UX_rank-1;
nc = max(0,min(n_c-1,floor(n_c*nUX_rank/dat_n_UX_rank)));
tmp_X_MS__ = squeeze(dat_X_true_MS___(:,1+nUX_rank,:));
tmp_F_MS__ = squeeze(flag_dat_exist_MS___(:,1+nUX_rank,:));
tmp_j_MS_ = find(sum(tmp_F_MS__,1));
if (numel(tmp_j_MS_)>0);
plot(0*dat_n_iteration + (0:dat_n_iteration-1),tmp_X_MS__(:,tmp_j_MS_),'o-','Color',c_(1+nc,:),'LineWidth',2,'MarkerFaceColor',c_(1+nc,:),'MarkerEdgeColor','k');
end;%if (numel(tmp_j_MS_)>0);
tmp_X_SM__ = squeeze(dat_X_true_SM___(:,1+nUX_rank,:));
tmp_F_SM__ = squeeze(flag_dat_exist_SM___(:,1+nUX_rank,:));
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
title(sprintf('%s delta_r_max_use %0.4f (%s)',str_cost,delta_r_max_use,str_delta),'Interpreter','none');
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
end;%if (sum(flag_dat_exist_MS___,'all') | sum(flag_dat_exist_SM___,'all'));
end;%for ndelta_r_max_factor=0:n_delta_r_max_factor-1;
end;%for norder_limit_MS=0:n_order_limit_MS-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_compute;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% debug functions start here. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
flag_compute = 1 & strcmp(platform,'OptiPlex');
flag_skip_to_plot = 0; %<-- skip over main computation and plot figures (as best as possible). ;
if flag_compute;
test_principled_marching_rib80sc_2_sub_debug;
end;%if flag_compute;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% relion comparison start here. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First write mrcs file containing image stack. ;
%%%%%%%%;
fname_image_mda = sprintf('%s/images_mda',dir_data);
M_x_c___ = MDA_read_r8(fname_image_mda);
n_image_sub = min(n_image_sub,size(M_x_c___,3));
M_x_c___ = M_x_c___(:,:,1:n_image_sub);
%%%%%%%%;
fname_mscope_params = sprintf('%s/mscope_params',dir_data);
mscope_params = textread(fname_mscope_params);
Voltage_kV = mscope_params(1);
SphericalAberration = mscope_params(2);
PixelSize = mscope_params(3);
AmplitudeContrast = mscope_params(4);
%%%%%%%%;
fname_image_mrcs = sprintf('%s_relion/test_pm_rib80sc_2.mrcs',dir_trunk);
fp = WriteMRCHeader(M_x_c___,PixelSize,fname_image_mrcs); fwrite(fp,M_x_c___,'float32'); fclose(fp);
%%%%%%%%;

%%%%%%%%;
% Now write star file linking to the mrcs image stack. ;
%%%%%%%%;
tmp_ImageName_list_ = cell(n_image_sub,1);
tmp_OpticsGroup_list_ = cell(n_image_sub,1);
tmp_DefocusU_list_ = cell(n_image_sub,1);
tmp_DefocusV_list_ = cell(n_image_sub,1);
tmp_DefocusAngle_list_ = cell(n_image_sub,1);
for nimage_sub=0:n_image_sub-1;
nctf = CTF_idx_(1+nimage_sub);
tmp_ImageName_list_{1+nimage_sub} = sprintf('%d@test_pm_rib80sc_2.mrcs',1+nimage_sub);
tmp_OpticsGroup_list_{1+nimage_sub} = sprintf('1');
tmp_DefocusU_list_{1+nimage_sub} = sprintf('%0.16f',CTF_Defocus_U_(1+nctf));
tmp_DefocusV_list_{1+nimage_sub} = sprintf('%0.16f',CTF_Defocus_V_(1+nctf));
tmp_DefocusAngle_list_{1+nimage_sub} = sprintf('%0.16f',CTF_Defocus_Angle_(1+nctf));
end;%for nimage_sub=0:n_image_sub-1;
fname_star = sprintf('%s_relion/test_pm_rib80sc_2.star',dir_trunk);
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
for nimage_sub=0:n_image_sub-1;
fprintf(fp,'%s %s %s %s %s\n',tmp_ImageName_list_{1+nimage_sub},tmp_OpticsGroup_list_{1+nimage_sub},tmp_DefocusU_list_{1+nimage_sub},tmp_DefocusV_list_{1+nimage_sub},tmp_DefocusAngle_list_{1+nimage_sub});
end;%for nimage_sub=0:n_image_sub-1;
fclose(fp);

%%%%%%%%;
% Now write shell script to call command. ;
%%%%%%%%;
jname = 'job_0';
tmp_slash = '\';
dir_job = sprintf('%s_relion/%s',dir_trunk,jname);
if (~exist(dir_job,'dir')); sprintf(' %% mkdir %s',dir_job); mkdir(dir_job); end;
fname_command = sprintf('%s_relion/test_pm_rib80sc_2.sh',dir_trunk);
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
fprintf(fp,'--i %s %s\n','test_pm_rib80sc_2.star',tmp_slash);
fprintf(fp,'--ctf %s\n',tmp_slash);
fprintf(fp,'--K 1 %s\n',tmp_slash);
fprintf(fp,'--sym C1 %s\n',tmp_slash);
fprintf(fp,'--flatten_solvent %s\n',tmp_slash);
fprintf(fp,'--zero_mask %s\n',tmp_slash);
fprintf(fp,'--dont_combine_weights_via_disc %s\n',tmp_slash);
fprintf(fp,'--pool 3 %s\n',tmp_slash);
fprintf(fp,'--pad 2 %s\n',tmp_slash);
fprintf(fp,'--skip_gridding %s\n',tmp_slash);
fprintf(fp,'--particle_diameter 400 %s\n',tmp_slash);
fprintf(fp,'--oversampling 1 %s\n',tmp_slash);
fprintf(fp,'--healpix_order 1 %s\n',tmp_slash);
fprintf(fp,'--offset_range 6 %s\n',tmp_slash);
fprintf(fp,'--offset_step 4 %s\n',tmp_slash);
fprintf(fp,'--j 6 %s\n',tmp_slash);
fprintf(fp,'--pipeline_control %s %s\n',jname,tmp_slash);
fprintf(fp,';\n');
fprintf(fp,'\n');
fclose(fp);
disp(sprintf(' %% fname_command: %s',fname_command));
type(fname_command);

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
%a_k_p_relion_quad_ = nufft3d3(n_X_u,X_u_0_(:)*eta,X_u_1_(:)*eta,X_u_2_(:)*eta,a_x_u_relion_pack_(:).*X_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
a_x_u_relion_pack__{1+nrelion_iter} = a_x_u_relion_pack_;
%a_k_p_relion_quad__{1+nrelion_iter} = a_k_p_relion_quad_;
end;%for nrelion_iter=0:n_relion_iter-1;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(1);clf;figbig;
prows = 3; pcols = 6; np=0;
subplot(prows,pcols,1+np);np=np+1;
isosurface_f_x_c_0(a_x_u_load_,98.5); axis vis3d; title('ori');
for nrelion_iter=0:n_relion_iter-1;
subplot(prows,pcols,1+np);np=np+1;
relion_iter = relion_iter_(1+nrelion_iter);
isosurface_f_x_c_0(a_x_u_relion_pack__{1+nrelion_iter},98.5); axis vis3d; title(sprintf('relion it%.3d',relion_iter));
end;%for nrelion_iter=0:n_relion_iter-1;
end;%if flag_plot;
%%%%%%%%;
flag_plot=1;
if flag_plot;
figure(2);clf;figbig;
tmp_prctile_ = [95,98.5,99];
subplot(1,2,1);
isosurface_f_x_c_0(a_x_u_load_,tmp_prctile_); axis vis3d; title('ori');
subplot(1,2,2);
nrelion_iter = n_relion_iter-1; relion_iter = relion_iter_(1+nrelion_iter);
isosurface_f_x_c_0(a_x_u_relion_pack__{1+nrelion_iter},tmp_prctile_); axis vis3d; title(sprintf('relion it%.3d',relion_iter));
end;%if flag_plot;

%%%%%%%%;
% Now check each iteration for correlation. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/X_relion_.mat',dir_trunk);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
X_relion_ = zeros(n_relion_iter,1);
for nrelion_iter=0:n_relion_iter-1;
relion_iter = relion_iter_(1+nrelion_iter);
fname_mrc = sprintf('%s/run_it%.3d_class001.mrc',dir_job,relion_iter);
tmp_t = tic();
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

fname_fig = sprintf('%s_jpg/X_relion_',dir_trunk);
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
% Now repeat for the different runs job_rib80s_1024_0 through job_rib80s_8192_0. ;
%%%%%%%%;
jname_ = {'job_rib80s_1024_0','job_rib80s_2048_0','job_rib80s_4096_0','job_rib80s_8192_0'};
n_jname = numel(jname_);
for njname=0:n_jname-1;
jname = jname_{1+njname};
fname_mat = sprintf('%s_mat/%s_X_relion_.mat',dir_trunk,jname);
if (~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_slash = '\';
dir_job = sprintf('%s_relion/%s',dir_trunk,jname);
relion_iter_ = 0:20:300; n_relion_iter = numel(relion_iter_);
%%%%%%%%;
a_x_u_relion_pack__ = cell(n_relion_iter,1);
for nrelion_iter=0:n_relion_iter-1;
relion_iter = relion_iter_(1+nrelion_iter);
fname_mrc = sprintf('%s/run_it%.3d_class001.mrc',dir_job,relion_iter);
a_x_u_relion_ = cast(ReadMRC(fname_mrc),'double');
a_x_u_relion_pack_ = reshape(a_x_u_relion_,[n_x_u*n_x_u,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = reshape(permute(reshape(a_x_u_relion_pack_,[n_x_u,n_x_u,x_u_res]),[3,1,2]),[n_x_u*x_u_res,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = reshape(permute(reshape(a_x_u_relion_pack_,[x_u_res,n_x_u,x_u_res]),[3,1,2]),[x_u_res*x_u_res,n_x_u])*x_u_pack_;
a_x_u_relion_pack_ = permute(reshape(a_x_u_relion_pack_,[x_u_res,x_u_res,x_u_res]),[3,1,2]);
a_x_u_relion_pack__{1+nrelion_iter} = a_x_u_relion_pack_;
end;%for nrelion_iter=0:n_relion_iter-1;
%%%%%%%%;
X_relion_ = zeros(n_relion_iter,1);
a_k_Y_relion_quad__ = cell(n_relion_iter,1);
for nrelion_iter=0:n_relion_iter-1;
relion_iter = relion_iter_(1+nrelion_iter);
fname_mrc = sprintf('%s/run_it%.3d_class001.mrc',dir_job,relion_iter);
tmp_t = tic();
a_k_p_relion_quad_ = nufft3d3(n_X_u,X_u_0_(:)*eta,X_u_1_(:)*eta,X_u_2_(:)*eta,a_x_u_relion_pack__{1+nrelion_iter}(:).*X_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_k_p_relion_quad_: %0.3fs',tmp_t));
tmp_t = tic();
[a_k_Y_relion_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_relion_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_spharm_1: a_k_Y_relion_quad_: %0.3fs',tmp_t));
a_k_Y_relion_quad__{1+nrelion_iter} = a_k_Y_relion_quad_;
[tmp_X_relion_orig,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,a_k_Y_relion_quad_);
[tmp_X_relion_flip,~,~,~,~,~,~] = register_spharm_to_spharm_wigner_0(n_k_p_r,k_p_r_,k_p_r_max,weight_3d_k_p_r_,0,l_max_,a_k_Y_quad_,flipY(n_k_p_r,l_max_,a_k_Y_relion_quad_));
if (tmp_X_relion_orig>=tmp_X_relion_flip); flag_flip = 0; end;
if (tmp_X_relion_orig< tmp_X_relion_flip); flag_flip = 1; end;
tmp_X_relion_reco = max(tmp_X_relion_orig,tmp_X_relion_flip);
if (verbose); disp(sprintf(' %% nrelion_iter %.2d/%.2d: tmp_X_relion_reco %0.3f <-- flag_flip %d',nrelion_iter,n_relion_iter,tmp_X_relion_reco,flag_flip)); end;
X_relion_(1+nrelion_iter) = tmp_X_relion_reco;
end;%for nrelion_iter=0:n_relion_iter-1;
save(fname_mat ...
     ,'relion_iter_' ...
     ,'a_k_Y_relion_quad__' ...
     ,'X_relion_' ...
     );
end;%if (~exist(fname_mat,'file'));
load(fname_mat);
end;%for njname=0:n_jname-1;
%%%%%%%%;
% produce figure. ;
%%%%%%%%;

fname_fig = sprintf('%s_jpg/job_rib80s_xxxx_X_relion_',dir_trunk);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
%%%%%%%%;
figure(4);clf;figmed;
c_80s__ = colormap_80s(); n_c_80s = size(c_80s__,1);
hold on; 
n_jname = numel(jname_);
for njname=0:n_jname-1;
jname = jname_{1+njname};
fname_mat = sprintf('%s_mat/%s_X_relion_.mat',dir_trunk,jname);
tmp_ = load(fname_mat); X_relion_ = tmp_.X_relion_;
nc = max(0,min(n_c_80s-1,floor(n_c_80s*njname/n_jname)));
plot(relion_iter_,X_relion_,'o-','Color',c_80s__(1+nc,:),'MarkerFaceColor',c_80s__(1+nc,:),'LineWidth',2);
end;%for njname=0:n_jname-1;
hold off;
set(gca,'Xtick',relion_iter_,'XTickLabel',relion_iter_);
xlim([0,max(relion_iter_)+1]); xlabel('relion iteration');
ylim([0,1]); ylabel('correlation'); grid on;
title('X_true','interpreter','none');
legend({'1024','2048','4096','8192'},'Location','NorthWest');
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning'); return;
