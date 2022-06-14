global_parameter=struct('type','parameter');
global_parameter.tolerance_master = 1e-6;
str_thisfunction = 'test_slice_vs_volume_integral_0';
verbose=1;
if (verbose); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

if isempty(global_parameter); global_parameter = struct('type','parameter'); end;
if (~isfield(global_parameter,'flag_recalc')); global_parameter.flag_recalc = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'flag_replot')); global_parameter.flag_replot = 0; end; %<-- parameter_bookmark. ;
if (~isfield(global_parameter,'tolerance_master')); global_parameter.tolerance_master = 1e-2; end; %<-- parameter_bookmark. ;
flag_recalc = global_parameter.flag_recalc;
flag_replot = global_parameter.flag_replot;
tolerance_master = global_parameter.tolerance_master;
nf=0;

%%%%%%%%;
n_x_u = 256;
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
%%%%%%%%;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;
%%%%%%%%;
k_u_0_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_1_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
k_u_2_ = periodize(0:n_x_u_pack-1,-n_x_u_pack/2,+n_x_u_pack/2)/2; %<-- box has diameter 2. ;
[k_u_0___,k_u_1___,k_u_2___] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_kkk_u = n_x_u_pack^3;

%%%%%%%%;
% Now set up grids for k_p_ ;
%%%%%%%%;
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

%%%%%%%%;
% Now set up grids for k_Y_ ;
%%%%%%%%;
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
Y_l_max_val_ = zeros(n_lm_sum,1);
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
Y_l_max_val_(1+tmp_index_) = l_max;
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
rng(0);
a_k_Y_base_ = ...
 ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ) ...
.* exp(-(Y_l_val_ + abs(Y_m_val_)).^2./(2*(Y_l_max_val_/4).^2)) ...
.* exp(-(Y_k_val_-k_p_r_max/2).^2./(2*(k_p_r_max/8)).^2) ...
;
flag_check=0;
if flag_check;
tmp_t = tic;
[a_k_p_quad_] = convert_spharm_to_k_p_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_base_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_base_ --> a_k_p_quad_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_Y_quad_] = convert_k_p_to_spharm_1(verbose,n_k_all,n_k_all_csum_,k_p_r_all_,k_p_azimu_b_all_,k_p_polar_a_all_,weight_3d_k_all_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_quad_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_quad --> a_k_Y_quad_ time %0.2fs',tmp_t));
disp(sprintf(' %% a_k_Y_base_ vs a_k_Y_quad_: %0.16f',fnorm(a_k_Y_base_-a_k_Y_quad_)/fnorm(a_k_Y_quad_))); %<-- this should be 2-3 digits. ;
end;%if flag_check;
a_k_Y_quad_ = a_k_Y_base_;
%%%%%%%%;
a_k_Y_quad__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% generate templates S_k_p_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
n_w_max = 2*(l_max_max+1);
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
template_k_eq_d = -1;
viewing_k_eq_d = 1.0; %<-- subsample just a few templates for visualization. ;
tmp_t = tic();
[ ...
 n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,~ ...
,~ ...
,~ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_5( ...
 k_p_r_max ...
,viewing_k_eq_d ...
,'L' ...
) ; %<-- obtain viewing angles on outer shell. ;
n_viewing_azimu_b_sum = sum(n_viewing_azimu_b_);
n_viewing_azimu_b_csum_ = cumsum([0;n_viewing_azimu_b_]);
if (verbose); disp(sprintf(' %% n_viewing_all %d n_viewing_polar_a %d n_viewing_azimu_b_sum %d',n_viewing_all,n_viewing_polar_a,n_viewing_azimu_b_sum)); end;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
] = ...
get_weight_2d_1( ...
 verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
n_w_max = max(n_w_);
n_w_sum = sum(n_w_);
n_w_csum_ = cumsum([0;n_w_]);
if (verbose); disp(sprintf(' %% n_w_max %d n_w_sum %d',n_w_max,n_w_sum)); end;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,~ ...
,n_S ...
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
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
n_S = n_viewing_all; n_w_max = max(n_w_);
%%%%%%%%;
flag_check=0;
if flag_check;
%%%%%%%%;
% test pm_template_2 vs get_template_1. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 tmp_S_k_p_wkS__ ...
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
disp(sprintf(' %% S_k_p_wkS__ vs tmp_S_k_p_wkS__: %0.16f',fnorm(S_k_p_wkS__-tmp_S_k_p_wkS__)/fnorm(S_k_p_wkS__)));
end;%if flag_check;
%%%%%%%%;

flag_disp=0;
if flag_disp;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wkS__(:,1)),[],colormap_beach());
end;%if flag_disp;

weight_Y_2d_ = zeros(n_lm_sum,1);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
weight_Y_2d_(1+tmp_index_) = weight_2d_k_p_r_(1+nk_p_r);
end;%for nk_p_r=0:n_k_p_r-1;

weight_2d_k_p_wkS__ = reshape(weight_2d_k_all_*4*pi^2,[n_w_sum,1])*reshape(viewing_weight_all_/(k_p_r_max^2),[1,n_S]);

E_3d_k_Y = sum(weight_Y_);
E_3d_k_p = sum(weight_3d_k_all_);
E_2d_k_p = sum(weight_2d_k_p_wkS__, 'all' );
E_2d_k_Y = sum(weight_Y_2d_);

L2_3d_k_Y = sum(abs(a_k_Y_base_).^2 .* weight_Y_);
flag_check=0;
if flag_check;
L2_3d_k_p = sum(abs(a_k_p_quad_).^2 .* weight_3d_k_all_);
end;%if flag_check;

L2_2d_k_p = sum(abs(S_k_p_wkS__).^2 .* weight_2d_k_p_wkS__,'all');
L2_2d_k_Y = sum(abs(a_k_Y_base_).^2 .* weight_Y_2d_);
if (verbose); disp(sprintf(' %% L2_2d_k_p %0.6f vs L2_2d_k_Y %0.6f: %0.16f',L2_2d_k_p,L2_2d_k_Y,fnorm(L2_2d_k_p - L2_2d_k_Y)/fnorm(L2_2d_k_Y))); end;

%%%%%%%%;
% Now throw in a set of images. ;
%%%%%%%%;
n_M = 5;
b_k_Y_quad_ykM__ = zeros(n_lm_sum,n_M);
b_k_Y_quad_ykM___ = zeros(n_lm_max,n_k_p_r,n_M);
for nM=0:n_M-1;
rng(1+nM);
b_k_Y_quad_ykM__(:,1+nM) = ...
 ( randn(n_lm_sum,1) + i*randn(n_lm_sum,1) ) ...
.* exp(-(Y_l_val_ + abs(Y_m_val_)).^2./(2*(Y_l_max_val_/4).^2)) ...
.* exp(-(Y_k_val_-k_p_r_max/2).^2./(2*(k_p_r_max/8)).^2) ...
;
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_k_Y_quad_ykM___(1:n_lm_(1+nk_p_r),1+nk_p_r,1+nM) = b_k_Y_quad_ykM__(1+tmp_index_,1+nM);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nM=0:n_M-1;
%%%%%%%%;
b_avg_k_Y_quad_yk_ = mean(b_k_Y_quad_ykM__,2);
b_avg_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
b_avg_k_Y_quad_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = b_avg_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
% generate templates M_k_p_wkSM___ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
M_k_p_wkSM___ = zeros(n_w_sum,n_S,n_M);
for nM=0:n_M-1;
if (verbose); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
tmp_t = tic();
[ ...
 M_k_p_wkS__ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(b_k_Y_quad_ykM___(:,:,1+nM),[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
M_k_p_wkS__ = reshape(M_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
M_k_p_wkSM___(:,:,1+nM) = M_k_p_wkS__; clear M_k_p_wkS__;
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));
end;%for nM=0:n_M-1;
if (verbose); disp(sprintf(' %% M_avg')); end;
tmp_t = tic();
[ ...
 M_avg_k_p_wkS__ ...
] = ...
pm_template_2( ...
 verbose ...
,l_max ...
,n_k_p_r ...
,reshape(b_avg_k_Y_quad_yk__,[n_lm_max,n_k_p_r]) ...
,viewing_k_eq_d/k_p_r_max ...
,template_k_eq_d ...
,n_w_max ...
);
M_avg_k_p_wkS__ = reshape(M_avg_k_p_wkS__,[n_w_max*n_k_p_r,n_S]);
tmp_t = toc(tmp_t); disp(sprintf(' %% pm_template_2: %0.2fs',tmp_t));

L2_2d_alM_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_alM_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_k_Y_quad_ykM__,a_k_Y_base_)).^2,weight_Y_2d_) , 'all' );
if (verbose); disp(sprintf(' %% L2_2d_alM_k_p %0.6f vs L2_2d_alM_k_Y %0.6f: %0.16f',L2_2d_alM_k_p,L2_2d_alM_k_Y,fnorm(L2_2d_alM_k_p - L2_2d_alM_k_Y)/fnorm(L2_2d_alM_k_Y))); end;

L2_2d_avg_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_avg_k_p_wkS__,S_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_avg_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_avg_k_Y_quad_yk_,a_k_Y_base_)).^2,weight_Y_2d_) , 'all' );
if (verbose); disp(sprintf(' %% L2_2d_avg_k_p %0.6f vs L2_2d_avg_k_Y %0.6f: %0.16f',L2_2d_avg_k_p,L2_2d_avg_k_Y,fnorm(L2_2d_avg_k_p - L2_2d_avg_k_Y)/fnorm(L2_2d_avg_k_Y))); end;

L2_2d_var_k_p = sum( bsxfun(@times,abs(bsxfun(@minus,M_k_p_wkSM___,M_avg_k_p_wkS__)).^2,weight_2d_k_p_wkS__) , 'all' );
L2_2d_var_k_Y = sum( bsxfun(@times,abs(bsxfun(@minus,b_k_Y_quad_ykM__,b_avg_k_Y_quad_yk_)).^2,weight_Y_2d_) , 'all' );
if (verbose); disp(sprintf(' %% L2_2d_var_k_p %0.6f vs L2_2d_var_k_Y %0.6f: %0.16f',L2_2d_var_k_p,L2_2d_var_k_Y,fnorm(L2_2d_var_k_p - L2_2d_var_k_Y)/fnorm(L2_2d_var_k_Y))); end;

disp(sprintf(' %% L2_2d_alM_k_Y vs (L2_2d_avg_k_Y*n_M + L2_2d_var_k_Y): %0.16f',fnorm(L2_2d_alM_k_Y-(L2_2d_avg_k_Y*n_M + L2_2d_var_k_Y))/fnorm(L2_2d_alM_k_Y)));



flag_check=0;
if flag_check;
%%%%%%%%;
% Now test out approximation to sum of exponentials: ;
%%%%%%%%;
nf=0;
rng(0);
sigma_ = exp(transpose([-4:0.125:0]));
lsigma_ = log(sigma_); 
n_sigma = numel(sigma_); 
n_a = 10; 
a_ = sort(rand(n_a,1),'ascend');
e_sa__ = bsxfun(@rdivide,-reshape(a_.^2,[1,n_a]),2*reshape(sigma_.^2,[n_sigma,1]));
f_sa__ = exp(e_sa__); 
f_s_ = sum(f_sa__,2); 
lf_s_ = log(f_s_); 
g_s_ = f_sa__(:,1);
lg_s_ = log(g_s_);
d_sa__ = bsxfun(@minus,e_sa__,e_sa__(:,1));
%lh_s_ = log(g_s_.*sum(exp(d_sa__),2));
lh_s_ = lg_s_ + sum(exp(d_sa__(:,2:end)),2);
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
plot(lsigma_,lf_s_,'k.-');
plot(lsigma_,lg_s_,'rx-');
plot(lsigma_,lh_s_,'go-');
hold off;
xlabel('lsigma'); ylabel('lsum');
end;%if flag_check;
