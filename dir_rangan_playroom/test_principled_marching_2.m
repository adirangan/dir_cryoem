clear; setup;
%%%%%%%%;
% Quick check of nufft1d3 ;
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
% Now define a molecule in real-space using a tensor-product grid. ;
%%%%%%%%;
x_p_r_max = 1.0;
n_x_c_0 = 84; [x_c_0_,x_c_0_weight_] = chebpts(n_x_c_0,x_p_r_max*[-1,+1]);
n_x_c_1 = 85; [x_c_1_,x_c_1_weight_] = chebpts(n_x_c_1,x_p_r_max*[-1,+1]);
n_x_c_2 = 86; [x_c_2_,x_c_2_weight_] = chebpts(n_x_c_2,x_p_r_max*[-1,+1]);
[X_c_0_,X_c_1_,X_c_2_] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_X_c = n_x_c_0*n_x_c_1*n_x_c_2;
[X_c_0_weight_,X_c_1_weight_,X_c_2_weight_] = ndgrid(x_c_0_weight_,x_c_1_weight_,x_c_2_weight_);
X_c_weight_ = X_c_0_weight_.*X_c_1_weight_.*X_c_2_weight_;
clear X_c_0_weight_ X_c_1_weight_ X_c_2_weight_;
k_p_r_max = 24.0;
n_k_c_0 = 94; [k_c_0_,k_c_0_weight_] = chebpts(n_k_c_0,k_p_r_max*[-1,+1]);
n_k_c_1 = 95; [k_c_1_,k_c_1_weight_] = chebpts(n_k_c_1,k_p_r_max*[-1,+1]);
n_k_c_2 = 96; [k_c_2_,k_c_2_weight_] = chebpts(n_k_c_2,k_p_r_max*[-1,+1]);
[K_c_0_,K_c_1_,K_c_2_] = ndgrid(k_c_0_,k_c_1_,k_c_2_); n_K_c = n_k_c_0*n_k_c_1*n_k_c_2;
[K_c_0_weight_,K_c_1_weight_,K_c_2_weight_] = ndgrid(k_c_0_weight_,k_c_1_weight_,k_c_2_weight_);
K_c_weight_ = K_c_0_weight_.*K_c_1_weight_.*K_c_2_weight_;
clear K_c_0_weight_ K_c_1_weight_ K_c_2_weight_;
n_source = 128; sg = 1/12; mu_source_ = zeros(n_source,3);
tmp_t = tic;
for nsource=1:n_source;
t = (nsource-1)/(n_source-1);
mu_source_(nsource,1+0) = 1/3*t * sin(2*pi*2*t);
mu_source_(nsource,1+1) = 1/3*t * cos(2*pi*2*t);
mu_source_(nsource,1+2) = -1/3 + 2/3*t ;
end;%for nsource=1:n_source;
a_x_c_form_ = zeros(n_X_c,1);
a_k_c_form_ = zeros(n_K_c,1);
for nsource=1:n_source;
a_x_c_form_ = a_x_c_form_ + gx(X_c_0_(:),mu_source_(nsource,1+0),sg).*gx(X_c_1_(:),mu_source_(nsource,1+1),sg).*gx(X_c_2_(:),mu_source_(nsource,1+2),sg);
a_k_c_form_ = a_k_c_form_ + gk(K_c_0_(:),mu_source_(nsource,1+0),sg).*gk(K_c_1_(:),mu_source_(nsource,1+1),sg).*gk(K_c_2_(:),mu_source_(nsource,1+2),sg);
end;%for nsource=1:n_source;
tmp_t = toc(tmp_t); disp(sprintf(' %% %d sources: time %0.2fs',n_source,tmp_t));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_c_quad_ = nufft3d3(n_X_c,X_c_0_(:)*eta,X_c_1_(:)*eta,X_c_2_(:)*eta,a_x_c_form_(:).*X_c_weight_(:),-1,1e-12,n_K_c,K_c_0_(:)/eta,K_c_1_(:)/eta,K_c_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_k_c_quad_ time %0.2fs',tmp_t));
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_quad_ = nufft3d3(n_K_c,K_c_0_(:)*eta,K_c_1_(:)*eta,K_c_2_(:)*eta,a_k_c_form_(:).*K_c_weight_(:),+1,1e-12,n_X_c,X_c_0_(:)/eta,X_c_1_(:)/eta,X_c_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_c_quad_ time %0.2fs',tmp_t));
flag_plot=0;
if flag_plot;
subplot(1,2,1); plot(sqrt(X_c_0_(:).^2+X_c_1_(:).^2+X_c_2_(:).^2),log10(abs(a_x_c_form_)),'.'); xlabel('|x_c|'); ylabel('log10(|a(x_c)|)');
subplot(1,2,2); plot(sqrt(K_c_0_(:).^2+K_c_1_(:).^2+K_c_2_(:).^2),log10(abs(a_k_c_form_)),'.'); xlabel('|k_c|'); ylabel('log10(|a(k_c)|)');
end;%flag_plot=0;
disp(sprintf(' %% nufft3d3: a_k_c_quad error: %0.16f',fnorm(a_k_c_form_-a_k_c_quad_)/fnorm(a_k_c_form_)));
disp(sprintf(' %% nufft3d3: a_x_c_quad error: %0.16f',fnorm(a_x_c_form_-a_x_c_quad_)/fnorm(a_x_c_form_)));

%%%%%%%%;
% Now convert to a_k_p_ ;
%%%%%%%%;
verbose=0;
tmp_t = tic;
k_eq_d = 1.0;
[n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,k_c_0_all_,k_c_1_all_,k_c_2_all_,J_node_,J_weight_,J_chebfun_,J_polyval_] = sample_sphere_7(verbose,k_p_r_max,k_eq_d,'L') ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
a_k_p_form_ = zeros(n_k_all,1);
tmp_t = tic;
for nsource=1:n_source;
a_k_p_form_ = a_k_p_form_ + gk(k_c_0_all_,mu_source_(nsource,1+0),sg).*gk(k_c_1_all_,mu_source_(nsource,1+1),sg).*gk(k_c_2_all_,mu_source_(nsource,1+2),sg);
end;%for nsource=1:n_source;
tmp_t = toc(tmp_t); disp(sprintf(' %% %d sources: time %0.2fs',n_source,tmp_t));
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_quad_ = nufft3d3(n_X_c,X_c_0_(:)*eta,X_c_1_(:)*eta,X_c_2_(:)*eta,a_x_c_form_(:).*X_c_weight_(:),-1,1e-12,n_k_all,k_c_0_all_/eta,k_c_1_all_/eta,k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_quad error: %0.16f',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_quad_ = nufft3d3(n_k_all,k_c_0_all_*eta,k_c_1_all_*eta,k_c_2_all_*eta,a_k_p_form_.*weight_k_all_,+1,1e-12,n_X_c,X_c_0_(:)/eta,X_c_1_(:)/eta,X_c_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_c_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_x_c_quad error: %0.16f',fnorm(a_x_c_form_-a_x_c_quad_)/fnorm(a_x_c_form_)));
disp(sprintf(' %% at this point one should ensure that a_k_p_quad_ on the outer shells (i.e., near k_p_r_max) has decayed to the desired precision.'));
flag_plot=0;
if flag_plot;
plot(k_p_r_all_,log10(abs(a_k_p_quad_)),'.'); xlabel('k'); ylabel('log10(|a(k)|)');
end;%if flag_plot;

%%%%%%%%;
% Now see what the function looks like on a uniform x_c_ grid. ;
%%%%%%%%;
x_u_res = 64;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,x_u_res);
[X_u_0_,X_u_1_,X_u_2_] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_X_u = x_u_res^3;
X_u_weight_ = (2*x_p_r_max/x_u_res)^3;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_quad_ = nufft3d3(n_k_all,k_c_0_all_*eta,k_c_1_all_*eta,k_c_2_all_*eta,a_k_p_form_.*weight_k_all_,+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nufft3d3: a_x_u_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
a_x_u_quad__ = zeros(n_X_u,n_k_p_r);
eta = pi/k_p_r_max; 
for nk_p_r=1:n_k_p_r;
tmp_ij_ = (1+n_k_all_csum_(nk_p_r)):n_k_all_csum_(1+nk_p_r);
tmp_n_k_all = numel(tmp_ij_);
tmp_t = tic;
a_x_u_quad__(:,nk_p_r) = nufft3d3(tmp_n_k_all,k_c_0_all_(tmp_ij_)*eta,k_c_1_all_(tmp_ij_)*eta,k_c_2_all_(tmp_ij_)*eta,a_k_p_form_(tmp_ij_).*weight_k_all_(tmp_ij_),+1,1e-12,n_X_u,X_u_0_(:)/eta,X_u_1_(:)/eta,X_u_2_(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% nk_p_r %.3d nufft3d3: a_x_u_quad__(:,%.3d) time %0.2fs',nk_p_r,nk_p_r,tmp_t));
end;%for nk_p_r=1:n_k_p_r;
b_x_u_quad__ = cumsum(a_x_u_quad__,2);
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
subplot(2,3,1); isosurface_f_x_u_0(reshape(real(a_x_u_quad_),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('full');
subplot(2,3,2); isosurface_f_x_u_0(reshape(real(b_x_u_quad__(:, 5)),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('K<= 5');
subplot(2,3,3); isosurface_f_x_u_0(reshape(real(b_x_u_quad__(:,10)),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('K<=10');
subplot(2,3,4); isosurface_f_x_u_0(reshape(real(b_x_u_quad__(:,15)),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('K<=15');
subplot(2,3,5); isosurface_f_x_u_0(reshape(real(b_x_u_quad__(:,20)),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('K<=20');
subplot(2,3,6); isosurface_f_x_u_0(reshape(real(b_x_u_quad__(:,25)),x_u_res,x_u_res,x_u_res),[90,95,99]); axisnotick3d; title('K<=25');
figbig;
fname_pre = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching_jpg/test_principled_marching_2_FIGB');
disp(sprintf(' %% writing %s',fname_pre));
print('-djpeg',sprintf('%s.jpg',fname_pre));
print('-depsc',sprintf('%s.eps',fname_pre));
end;%if flag_plot;

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
verbose=0;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=1:n_k_p_r;
%l_max_(nk_p_r) = 1+ceil(2*pi*k_p_r_(nk_p_r));
l_max_(nk_p_r) = 1+ceil(k_p_r_(nk_p_r));
end;%for nk_p_r=1:n_k_p_r;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_); dWtdkd__l_max_max = 2*l_max_max;
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
Y_l_val_ = zeros(n_lm_sum,1);
Y_m_val_ = zeros(n_lm_sum,1);
for nk_p_r=1:n_k_p_r;
l_max = l_max_(nk_p_r);
tmp_l_val_ = zeros(n_lm_(nk_p_r),1);
tmp_m_val_ = zeros(n_lm_(nk_p_r),1);
na=0; 
for l_val=0:l_max;
for m_val=-l_val:+l_val;
tmp_l_val_(1+na) = l_val;
tmp_m_val_(1+na) = m_val;
na=na+1;
end;%for m_val=-l_val:+l_val;
end;%for l_val=0:l_max;
tmp_ij_ = n_lm_csum_(nk_p_r) + (1:n_lm_(nk_p_r));
Y_l_val_(tmp_ij_) = tmp_l_val_;
Y_m_val_(tmp_ij_) = tmp_m_val_;
end;%for nk_p_r=1:n_k_p_r;
weight_Y_ = zeros(n_lm_sum,1);
for nk_p_r=1:n_k_p_r;
tmp_ij_ = n_lm_csum_(nk_p_r) + (1:n_lm_(nk_p_r));
weight_Y_(tmp_ij_) = weight_k_p_r_(nk_p_r);
end;%for nk_p_r=1:n_k_p_r;
%%%%%%%%;
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_p_form_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% Here we should ensure that a_k_Y_ on the outer shells has decayed to the desired precision.'));
disp(sprintf(' %% Moreover, we should ensure that a_k_Y_(l,m) has decayed for large l,m at each shell.'));
flag_plot=0;
if flag_plot;
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_)),[-10,0],colormap_beach());
end;%if flag_plot;
tmp_t = tic;
[a_k_p_quad_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_k_p_r_,l_max_,a_k_Y_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad --> a_k_p_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% nufft3d3: a_k_p_quad error: %0.16f',fnorm(a_k_p_form_-a_k_p_quad_)/fnorm(a_k_p_form_)));
%%%%%%%%;
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=1:n_k_p_r;
tmp_ij_ = n_lm_csum_(nk_p_r) + (1:n_lm_(nk_p_r));
a_k_Y_quad__(1:n_lm_(nk_p_r),nk_p_r) = a_k_Y_quad_(tmp_ij_);
end;%for nk_p_r=1:n_k_p_r;

verbose=1;
%%%%%%%%;
% Now test out principled marching. ;
% First set up cost matrix. ;
%%%%%%%%;
n_polar_a = n_m_max; polar_a_ = linspace(-pi,pi,n_polar_a+1); polar_a_ = polar_a_(1:end-1);
n_azimu_b = n_m_max; azimu_b_ = linspace(0,2*pi,n_azimu_b+1); azimu_b_ = azimu_b_(1:end-1);
n_gamma_z = n_m_max; gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = gamma_z_(1:end-1);
weight_so3 = (2*pi)*(2*pi)*4; %<-- total volume of so3. ;
weight_sub = ((2*pi)/n_m_max)^3; %<-- abs(sin(polar_a))*weight_sub is used for each summand on so3. ;
%%%%%%%%;
X_ori_ = zeros(n_k_p_r,n_k_p_r);
for nk_p_r_1=1:n_k_p_r;
for nk_p_r_2=nk_p_r_1:n_k_p_r;
tmp_l_max = min(l_max_(nk_p_r_1),l_max_(nk_p_r_2));
tmp_n_lm = (tmp_l_max+1).^2;
tmp_ij_1_ = n_lm_csum_(nk_p_r_1) + (1:tmp_n_lm);
tmp_ij_2_ = n_lm_csum_(nk_p_r_2) + (1:tmp_n_lm);
polar_a = 0;
X_ori_(nk_p_r_1,nk_p_r_2) = register_spharm_to_spharm_2(verbose,1,1,1,tmp_l_max,a_k_Y_quad_(tmp_ij_1_),a_k_Y_quad_(tmp_ij_2_));
X_ori_(nk_p_r_2,nk_p_r_1) = conj(X_ori_(nk_p_r_1,nk_p_r_2));
tmp_sum = 0;
for npolar_a=1:n_polar_a;
polar_a = polar_a_(npolar_a);
[tmp_X_tau__] = register_spharm_to_spharm_single_beta_2(verbose,1,1,1,tmp_l_max,a_k_Y_quad_(tmp_ij_1_),a_k_Y_quad_(tmp_ij_2_),polar_a,n_m_max,azimu_b_,n_m_max,gamma_z_,[],[],[],[]);
tmp_sum = tmp_sum + sum(tmp_X_tau__,'all')*abs(sin(polar_a))*weight_sub; %<-- need quadrature weight to ensure uniform measure over SO3. ;
end;%for npolar_a=1:n_polar_a;
if (verbose>1); disp(sprintf(' %% nk_p_r_1 %d nk_p_r_2 %d l_max %d,%d-->%d n_lm %d,%d-->%d --> sum %0.16f',nk_p_r_1,nk_p_r_2,l_max_(nk_p_r_1),l_max_(nk_p_r_2),tmp_l_max,n_lm_(nk_p_r_1),n_lm_(nk_p_r_2),tmp_n_lm,real(tmp_sum))); end;
X_tau_(nk_p_r_1,nk_p_r_2) = tmp_sum;
X_tau_(nk_p_r_2,nk_p_r_1) = conj(tmp_sum);
end;%for nk_p_r_2=nk_p_r_1:n_k_p_r;
if (verbose>0); disp(sprintf(' %% nk_p_r_1 %d nk_p_r_2 %d l_max %d,%d-->%d n_lm %d,%d-->%d --> sum %0.16f',nk_p_r_1,nk_p_r_2,l_max_(nk_p_r_1),l_max_(nk_p_r_2),tmp_l_max,n_lm_(nk_p_r_1),n_lm_(nk_p_r_2),tmp_n_lm,real(tmp_sum))); end;
end;%for nk_p_r_1=1:n_k_p_r;
%%%%%%%%;
X_ = real(X_ori_)*weight_so3 - real(X_tau_);
%%%%%%%%;
% Compare to cost matrix calculated using a_k_Y_quad__: ;
%%%%%%%%;
flag_check=0;
if flag_check;
Z_ = principled_marching_cost_matrix_0(n_k_p_r,l_max_max,a_k_Y_quad__);
disp(sprintf(' %% X_ vs Z_: %0.16f',fnorm(Z_-X_)/fnorm(X_)));
end;%if flag_check;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
[tmp_UX_,tmp_SX_,tmp_VX_] = svds(X_,n_k_p_r);
a_UX_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nrank=1:n_k_p_r;
for nk_p_r=1:n_k_p_r;
tmp_l_max = l_max_(nk_p_r);
tmp_n_lm = (tmp_l_max+1).^2;
tmp_ij_ = n_lm_csum_(nk_p_r) + (1:tmp_n_lm);
a_UX_Y_quad__(1:tmp_n_lm,nrank) = a_UX_Y_quad__(1:tmp_n_lm,nrank) + tmp_UX_(nk_p_r,nrank)*a_k_Y_quad_(tmp_ij_);
end;%for nk_p_r=1:n_k_p_r;
end;%for nrank=1:n_k_p_r;
%%%%%%%%;
% Check cost for each principal-mode. ;
%%%%%%%%;
%%%%%%%%;
% First set up a tensor-product spherical grid (in k_p_ space). ;
%%%%%%%%;
k_u_res = 64;
k_u_polar_a_ = linspace(0,pi,k_u_res);
k_u_azimu_b_ = linspace(0,2*pi,2*k_u_res);
[K_u_polar_a_,K_u_azimu_b_] = ndgrid(k_u_polar_a_,k_u_azimu_b_); n_K_u = k_u_res*2*k_u_res;
K_u_weight_ = sin(K_u_polar_a_);
%%%%%%%%;
for nrank=1:n_k_p_r;
%[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_UX_Y_quad__(:,nrank)),k_u_res,2*k_u_res);
[tmp_X,tmp_X_ori,tmp_X_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_UX_Y_quad__(:,nrank),a_UX_Y_quad__(:,nrank));
tmp_Z = transpose(tmp_UX_(:,nrank))*X_*(tmp_UX_(:,nrank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %0.2f tmp_X %0.2f tmp_X_ori*tmp_weight_so3 %0.2f tmp_X_tau %0.2f ratio %0.2f',nrank,n_k_p_r,tmp_Z,tmp_X,tmp_X_ori*tmp_weight_so3,tmp_X_tau,(tmp_X_ori*tmp_weight_so3)/tmp_X_tau));
end;%for nrank=1:n_k_p_r;

%%%%%%%%;
% Now do the same, but for the continuous version. ;
%%%%%%%%;
verbose=1;
%%%%%%%%;
% First check orthogonality/orthonormality using chebfun. ;
%%%%%%%%;
I_form_ = eye(1+n_k_p_r);
I_cheb_ = zeros(1+n_k_p_r,1+n_k_p_r);
for nk_p_r_1=0:n_k_p_r;
for nk_p_r_2=nk_p_r_1:n_k_p_r;
I_cheb_(1+nk_p_r_1,1+nk_p_r_2) = integral(J_chebfun_{1+nk_p_r_1}*J_chebfun_{1+nk_p_r_2}*chebfun(@(x) (x+1)^2));
I_cheb_(1+nk_p_r_2,1+nk_p_r_1) = I_cheb_(1+nk_p_r_1,1+nk_p_r_2);
end;%for nk_p_r_2=nk_p_r_1:n_k_p_r;
end;%for nk_p_r_1=0:n_k_p_r;
disp(sprintf(' %% I_form_ vs I_cheb_: %0.16f',fnorm(I_form_-I_cheb_)/fnorm(I_form_)));
%%%%%%%%;
% Now check orthogonality/orthonormality using quadrature. ;
%%%%%%%%;
I_form_ = eye(n_k_p_r);
I_quad_ = zeros(n_k_p_r,n_k_p_r);
for nk_p_r_1=0:n_k_p_r-1;
for nk_p_r_2=nk_p_r_1:n_k_p_r-1;
I_quad_(1+nk_p_r_1,1+nk_p_r_2) = sum(transpose(J_polyval_(1+nk_p_r_1,:)).*transpose(J_polyval_(1+nk_p_r_2,:)).*J_weight_(:));
I_quad_(1+nk_p_r_2,1+nk_p_r_1) = I_quad_(1+nk_p_r_1,1+nk_p_r_2);
end;%for nk_p_r_2=nk_p_r_1:n_k_p_r-1;
end;%for nk_p_r_1=0:n_k_p_r-1;
disp(sprintf(' %% I_form_ vs I_quad_: %0.16f',fnorm(I_form_-I_quad_)/fnorm(I_form_)));
%%%%%%%%;
% Now use quadrature to determine the jacobi-coefficients of each of the spherical-harmonic terms. ;
%%%%%%%%;
a_J_Y_quad__ = a_k_Y_quad__*diag(J_weight_)*transpose(J_polyval_(1:n_k_p_r,:));
%%%%%%%%;
% Now check to see if the jacobi-coefficients accurately reproduce the original spherical-harmonic terms. ;
%%%%%%%%;
a_k_Y_reco__ = a_J_Y_quad__*J_polyval_(1:n_k_p_r,:);
flag_plot=0;
if flag_plot;
n_plot = 12;
tmp_ij_ = max(1,min(n_k_p_r,round(linspace(1,n_k_p_r,n_plot+1))));
for nplot=1:n_plot;
subplot(3,4,nplot); cla; hold on;
tmp_ij_sub_ = tmp_ij_(nplot):(tmp_ij_(nplot+1)-1);
plot(k_p_r_,transpose(real(a_k_Y_quad__(tmp_ij_sub_,:))),'r.-');
plot(k_p_r_,transpose(imag(a_k_Y_quad__(tmp_ij_sub_,:))),'b.-');
plot(k_p_r_,transpose(real(a_k_Y_reco__(tmp_ij_sub_,:))),'ro');
plot(k_p_r_,transpose(imag(a_k_Y_reco__(tmp_ij_sub_,:))),'bo');
hold off; 
xlabel('k_p_r_','Interpreter','none'); title(sprintf('%d-->%d',tmp_ij_sub_(1),tmp_ij_sub_(end)));
end;%for nplot=1:n_plot;
figbig;
end;%if flag_plot;
disp(sprintf(' %% a_k_Y_quad__ vs a_k_Y_reco__: %0.16f',fnorm(a_k_Y_quad__ - a_k_Y_reco__)/fnorm(a_k_Y_quad__)));
%%%%%%%%;
% Now build cost-matrix for a_J_Y_quad_. ;
%%%%%%%%;
[Y_,Y_ori_,Y_tau_] = principled_marching_cost_matrix_0(n_k_p_r,l_max_max,a_J_Y_quad__);
%%%%%%%%;
% Now perform singular-value-decomposition: ;
%%%%%%%%;
[tmp_UY_,tmp_SY_,tmp_VY_] = svds(Y_,n_k_p_r);
a_UY_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nrank=1:n_k_p_r;
for nk_p_r=1:n_k_p_r;
a_UY_Y_quad__(:,nrank) = a_UY_Y_quad__(:,nrank) + tmp_UY_(nk_p_r,nrank)*a_J_Y_quad__(:,nk_p_r);
end;%for nk_p_r=1:n_k_p_r;
end;%for nrank=1:n_k_p_r;
a_UX_Y_reco__ = a_UY_Y_quad__;
%%%%%%%%;
% check cost for each princpled-mode. ;
%%%%%%%%;
for nrank=1:n_k_p_r;
[tmp_Y,tmp_Y_ori,tmp_Y_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_UY_Y_quad__(:,nrank),a_UY_Y_quad__(:,nrank));
tmp_Z = transpose(tmp_UY_(:,nrank))*Y_*(tmp_UY_(:,nrank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %0.2f tmp_Y %0.2f tmp_Y_ori*tmp_weight_so3 %0.2f tmp_Y_tau %0.2f ratio %0.2f',nrank,n_k_p_r,tmp_Z,tmp_Y,tmp_Y_ori*tmp_weight_so3,tmp_Y_tau,(tmp_Y_ori*tmp_weight_so3)/tmp_Y_tau));
[tmp_X,tmp_X_ori,tmp_X_tau,tmp_weight_so3] = principled_marching_cost_0(verbose,n_m_max,l_max_max,a_UX_Y_reco__(:,nrank),a_UX_Y_reco__(:,nrank));
tmp_Z = transpose(tmp_UX_(:,nrank))*X_*(tmp_UX_(:,nrank));
disp(sprintf(' %% mode %.3d/%.3d: tmp_Z %0.2f tmp_X %0.2f tmp_X_ori*tmp_weight_so3 %0.2f tmp_X_tau %0.2f ratio %0.2f',nrank,n_k_p_r,tmp_Z,tmp_X,tmp_X_ori*tmp_weight_so3,tmp_X_tau,(tmp_X_ori*tmp_weight_so3)/tmp_X_tau));
end;%for nrank=1:n_k_p_r;
%%%%%%%%;
% Note that there is a very different weighting to the two cost functions, ;
% potentially explaining the large discrepancy between these two calculations. ;
%%%%%%%%;
flag_plot=0;
if flag_plot;
clf;
colormap(colormap_beach());
subplot(2,2,1); imagesc(log10(abs(tmp_UX_)),[-3,0]); xlabel('rank'); ylabel('shell'); title('log10(abs(UX)) [-3,0]'); 
subplot(2,2,2); plot(log10(abs(diag(tmp_SX_))),'ko'); xlabel('rank'); ylabel('log10(\sigma)'); title('log10(SX)');
subplot(2,2,3); imagesc(log10(abs(tmp_UY_)),[-3,0]); xlabel('rank'); ylabel('order'); title('log10(abs(UJ)) [-3,0]'); 
subplot(2,2,4); plot(log10(abs(diag(tmp_SY_))),'ko'); xlabel('rank'); ylabel('log10(\sigma)'); title('log10(SJ)');
figbig;
fname_pre = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching_jpg/test_principled_marching_2_FIGC');
disp(sprintf(' %% writing %s',fname_pre));
print('-djpeg',sprintf('%s.jpg',fname_pre));
print('-depsc',sprintf('%s.eps',fname_pre));
end;%if flag_plot;

%%%%%%%%;
% Now look at the functions on each shell associated with these 'principal-modes'. ;
%%%%%%%%;
flag_plot=1;
% Now convert the a_UX_Y_reco__ onto this spherical grid (in k_p space). ;
if flag_plot;
n_plot = 6;
%plot_nk_p_r_ = max(1,min(n_k_p_r,round(linspace(1,n_k_p_r,n_plot))));
plot_nk_p_r_ = 1:n_plot;
quad_lim_ = 0.5 * abs(a_UX_Y_quad__(1,1)) * [-1,+1];
reco_lim_ = 0.5 * abs(a_UX_Y_reco__(1,1)) * [-1,+1];
for nplot=1:n_plot;
nk_p_r = plot_nk_p_r_(nplot);
[b_k_p_quad_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_UX_Y_quad__(:,nk_p_r)),k_u_res,2*k_u_res);
subplot(6,n_plot,1 + (nplot-1) + 0*n_plot); imagesc(real(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('real nk_p_r: %d, real(quad)',nk_p_r),'Interpreter','none');
subplot(6,n_plot,1 + (nplot-1) + 1*n_plot); imagesc(imag(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf('imag nk_p_r: %d, imag(quad)',nk_p_r),'Interpreter','none');
subplot(6,n_plot,1 + (nplot-1) + 2*n_plot); imagesc( abs(b_k_p_quad_),quad_lim_); axisnotick; title(sprintf(' abs nk_p_r: %d,  abs(quad)',nk_p_r),'Interpreter','none');
[b_k_p_reco_] = reshape(convert_spharm_to_k_p_1(verbose,n_K_u,0,ones(n_K_u,1),K_u_azimu_b_(:),K_u_polar_a_(:),K_u_weight_(:),K_u_weight_(:),1,1,1,l_max_max,a_UX_Y_reco__(:,nk_p_r)),k_u_res,2*k_u_res);
subplot(6,n_plot,1 + (nplot-1) + 3*n_plot); imagesc(real(b_k_p_reco_),reco_lim_); axisnotick; title(sprintf('real nk_p_r: %d, real(reco)',nk_p_r),'Interpreter','none');
subplot(6,n_plot,1 + (nplot-1) + 4*n_plot); imagesc(imag(b_k_p_reco_),reco_lim_); axisnotick; title(sprintf('imag nk_p_r: %d, imag(reco)',nk_p_r),'Interpreter','none');
subplot(6,n_plot,1 + (nplot-1) + 5*n_plot); imagesc( abs(b_k_p_reco_),reco_lim_); axisnotick; title(sprintf(' abs nk_p_r: %d,  abs(reco)',nk_p_r),'Interpreter','none');
end;%for nplot=1:n_plot;
colormap(colormap_beach());
figbig;
fname_pre = sprintf('/data/rangan/dir_cryoem/dir_rangan_playroom/dir_principled_marching_jpg/test_principled_marching_2_FIGD');
disp(sprintf(' %% writing %s',fname_pre));
print('-djpeg',sprintf('%s.jpg',fname_pre));
print('-depsc',sprintf('%s.eps',fname_pre));
end;%if flag_plot;

