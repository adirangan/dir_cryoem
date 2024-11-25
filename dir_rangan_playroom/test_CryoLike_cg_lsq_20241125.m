function nf = test_CryoLike_cg_lsq_20241125(k_eq_d_double,k_int,n_x_u,nf);
%%%%%%%%;
% Tests 3d-integration. ;
%%%%%%%%;
% Input: ;
% k_eq_d_double: double resolution-parameter. ;
%     k_eq_d_double is defined as (2*pi)*k_eq_d. ;
%     k_eq_d is the distance (along the equator) between sample-point in k-space on the largest k-shell. ;
% k_int: int largest wavenumber. ;
%     k_p_r_max is defined as k_int/(2*pi). ;
% n_x_u: int number of spatial-points alongside box of side-length diameter_x_c:=2.0. ;
% nf: int figure number for plotting. ;
%%%%%%%%;
% Additional internal parameters. ;
% sigma_k: This is the standard-deviation of the gaussian-envelope in k-space. ;
%     I set it to be k_p_r_max/4.0, but you can increase it. ;
%     Doing so will narrow the support of the spatial-gaussian, ;
%     increasing the resolution at the cost of integration-accuracy. ;
% delta_x_c_max: This is the maximum translation allowed for generating the synthetic volume. ;
%     I set it to be half-diameter_max/2.0, but you can increase it. ;
%     Doing so will reduce the integration and reconstruction-accuracy. ;
%%%%%%%%;
% defaults: k_eq_d_double = 1.0; k_int = 65; n_x_u = 128; nf=0;

str_thisfunction = 'test_CryoLike_cg_lsq_20241125';

if (nargin<1);
k_eq_d_double_ = [1.0,4.0];
n_k_eq_d_double = numel(k_eq_d_double_);
nf=0;
for nk_eq_d_double=0:n_k_eq_d_double-1;
k_eq_d_double = k_eq_d_double_(1+nk_eq_d_double);
nf = test_CryoLike_cg_lsq_20241125(k_eq_d_double,[],[],nf);
end;%for nk_eq_d_double=0:n_k_eq_d_double-1;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if nargin<1+na; k_eq_d_double = []; end; na=na+1;
if nargin<1+na; k_int = []; end; na=na+1;
if nargin<1+na; n_x_u = []; end; na=na+1;
if nargin<1+na; nf = []; end; na=na+1;

% defaults for testing: ;
% clear; str_thisfunction = 'test_CryoLike_cg_lsq_20241125'; k_eq_d_double = []; k_int = []; n_x_u = []; nf = [];

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

rng(0);
flag_recalc=0; flag_replot=0;
flag_verbose=1; flag_disp=1;

if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

if isempty(k_eq_d_double); k_eq_d_double = 1.0; end;
if isempty(k_int); k_int = 65; end;
if isempty(n_x_u); n_x_u = 128; end;
if isempty(nf); nf = 0; end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% setting k_eq_d_double = %0.6f %%<-- should be roughly 1.0 or less for accurate integration. ',k_eq_d_double)); end;
if (flag_verbose>0); disp(sprintf(' %% setting k_int = %d ',k_int)); end;
if (flag_verbose>0); disp(sprintf(' %% setting n_x_u = %d ',n_x_u)); end;

half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = half_diameter_x_c;
%x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_0_ = transpose(x_u_0_(1:n_x_u));
%x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_1_ = transpose(x_u_1_(1:n_x_u));
%x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u+1); x_u_2_ = transpose(x_u_2_(1:n_x_u));
x_u_0_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_0 = mean(diff(x_u_0_));
x_u_1_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_1 = mean(diff(x_u_1_));
x_u_2_ = transpose(linspace(-x_p_r_max,+x_p_r_max,n_x_u)); d_x_2 = mean(diff(x_u_2_));
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u^3;
weight_xxx_u_ = d_x_0 * d_x_1 * d_x_2 ;

%%%%%%%%;
% Now set up k-quadrature on sphere. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_T_vs_L = 'C';
flag_unif_vs_adap = 0; flag_tensor_vs_adap = 0; %<-- This is set to match test_ssnll_from_a_k_Y_12 ;
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
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 4; n_plot = p_row*p_col;
for nplot=0:n_plot-1;
nk_p_r = max(0,min(n_k_p_r-1,round(n_k_p_r*nplot/n_plot)));
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
subplot(p_row,p_col,1+nplot);
plot3(k_c_0_all_(1+tmp_index_),k_c_1_all_(1+tmp_index_),k_c_2_all_(1+tmp_index_),'.');
axis equal; axis vis3d; axisnotick3d;
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
end;%for nplot=0:n_plot-1;
end;%if flag_disp;
%%%%;

%%%%%%%%;
% test quadrature. ;
%%%%%%%%
h3d_ = @(kd) 4*pi*( sin(kd) - (kd).*cos(kd) ) ./ kd.^3 ; % calculates <f_j,f_k>, normalized so that <f,f> = 4*pi/3;
n_source = 53;
delta_a_c__ = zeros(3,n_source);
delta_b_c__ = zeros(3,n_source);
for nsource=0:n_source-1;
rng(1+nsource);
delta_a_c_ = half_diameter_x_c*(2*rand(3,1)-1);%delta_a_c_ = 2*randn(3,1)/k_p_r_max;
delta_a_c__(:,1+nsource) = delta_a_c_;
delta_b_c_ = half_diameter_x_c*(2*rand(3,1)-1);%delta_b_c_ = 2*randn(3,1)/k_p_r_max;
delta_b_c__(:,1+nsource) = delta_b_c_;
end;%for nsource=0:n_source-1;
a_k_p_form_ = zeros(n_k_all,1);
b_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
a_k_p_form_ = a_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_a_c_(1+0) + k_c_1_all_*delta_a_c_(1+1) + k_c_2_all_*delta_a_c_(1+2)));
delta_b_c_ = delta_b_c__(:,1+nsource);
b_k_p_form_ = b_k_p_form_ + exp(+i*2*pi*(k_c_0_all_*delta_b_c_(1+0) + k_c_1_all_*delta_b_c_(1+1) + k_c_2_all_*delta_b_c_(1+2)));
end;%for nsource=0:n_source-1;
I_a_quad = sum(a_k_p_form_.*weight_3d_k_all_);
I_b_quad = sum(b_k_p_form_.*weight_3d_k_all_);
I_a_form = 0;
I_b_form = 0;
for nsource=0:n_source-1;
delta_a_c_ = delta_a_c__(:,1+nsource);
delta_b_c_ = delta_b_c__(:,1+nsource);
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_a_c_);
I_a_form = I_a_form + h3d_(tmp_kd)*k_p_r_max^3;
tmp_kd = 2*pi*k_p_r_max*fnorm(delta_b_c_);
I_b_form = I_b_form + h3d_(tmp_kd)*k_p_r_max^3;
end;%for nsource=0:n_source-1;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Testing quadrature:')); end;
if (flag_verbose>0); disp(sprintf(' %% %% I_a_form vs I_a_quad %0.16f %%<-- should be <1e-6',fnorm(I_a_form-I_a_quad)/fnorm(I_a_form))); end;
if (flag_verbose>0); disp(sprintf(' %% %% I_b_form vs I_b_quad %0.16f %%<-- should be <1e-6',fnorm(I_b_form-I_b_quad)/fnorm(I_b_form))); end;
%%%%%%%%;
clear n_source delta_a_c__ delta_b_c__ a_k_p_form_ b_k_p_form_ ;
%%%%%%%%;

%%%%%%%%;
% Now integrate gaussian. ;
% The width of the gaussian is specifically chosen so that it is well-integrated in both x_u and k_p. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Transforming centered gaussian:')); end;
sigma_k = k_p_r_max/4.0; sigma_x = 1.0/max(1e-12,sigma_k)/(2*pi); %<-- fix sigma_k. ;
%sigma_x = 1/16; sigma_k = 1.0/max(1e-12,sigma_x)/(2*pi); %<-- fix sigma_x. ;
nrm_x_from_k = sigma_x^3 * sqrt((2*pi)^3);
g_k_p_form_ = 1/sqrt(2*pi)^3 .* 1/max(1e-12,sigma_k^3) .* exp(-(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2)/max(1e-12,2*sigma_k^2)) ;
g_k_p_form_l1 = sum(g_k_p_form_.*weight_3d_k_all_,'all');
g_k_p_full_l1 = 1.0;
g_k_p_form_l2 = sqrt(sum(abs(g_k_p_form_).^2.*weight_3d_k_all_,'all'));
g_k_p_full_l2 = sqrt( (sqrt(2*pi) * sigma_k * sqrt(2)).^(-3) );
g_x_u_form_ = nrm_x_from_k * 1/sqrt(2*pi)^3 .* 1/max(1e-12,sigma_x^3) .* exp(-(x_u_0___.^2 + x_u_1___.^2 + x_u_2___.^2)/max(1e-12,2*sigma_x^2)) ;
g_x_u_form_l1 = sum(g_x_u_form_.*weight_xxx_u_,'all');
g_x_u_full_l1 = nrm_x_from_k*1.0;
g_x_u_form_l2 = sqrt(sum(abs(g_x_u_form_).^2.*weight_xxx_u_,'all'));
g_x_u_full_l2 = nrm_x_from_k * sqrt( (sqrt(2*pi) * sigma_x * sqrt(2)).^(-3) );
if (flag_verbose>0); disp(sprintf(' %% %% g_k_p_full_l1 vs g_k_p_form_l1: %0.16f',fnorm(g_k_p_full_l1 - g_k_p_form_l1)/max(1e-12,fnorm(g_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_k_p_full_l2 vs g_k_p_form_l2: %0.16f',fnorm(g_k_p_full_l2 - g_k_p_form_l2)/max(1e-12,fnorm(g_k_p_full_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_x_u_full_l1 vs g_x_u_form_l1: %0.16f',fnorm(g_x_u_full_l1 - g_x_u_form_l1)/max(1e-12,fnorm(g_x_u_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_x_u_full_l2 vs g_x_u_form_l2: %0.16f',fnorm(g_x_u_full_l2 - g_x_u_form_l2)/max(1e-12,fnorm(g_x_u_full_l2)))); end;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
g_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,g_x_u_form_(:).*weight_xxx_u_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: g_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
g_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,g_k_p_form_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: g_x_u_quad_ time %0.2fs',tmp_t)); end;
g_x_u_quad_ = reshape(g_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
g_k_p_quad_l1 = sum(g_k_p_quad_.*weight_3d_k_all_,'all');
g_k_p_quad_l2 = sqrt(sum(abs(g_k_p_quad_).^2.*weight_3d_k_all_,'all'));
g_x_u_quad_l1 = sum(g_x_u_quad_.*weight_xxx_u_,'all');
g_x_u_quad_l2 = sqrt(sum(abs(g_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% g_k_p_full_l1 vs g_k_p_quad_l1: %0.16f',fnorm(g_k_p_full_l1 - g_k_p_quad_l1)/max(1e-12,fnorm(g_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_k_p_full_l2 vs g_k_p_quad_l2: %0.16f',fnorm(g_k_p_full_l2 - g_k_p_quad_l2)/max(1e-12,fnorm(g_k_p_full_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_x_u_full_l1 vs g_x_u_quad_l1: %0.16f',fnorm(g_x_u_full_l1 - g_x_u_quad_l1)/max(1e-12,fnorm(g_x_u_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% g_x_u_full_l2 vs g_x_u_quad_l2: %0.16f',fnorm(g_x_u_full_l2 - g_x_u_quad_l2)/max(1e-12,fnorm(g_x_u_full_l2)))); end;

%%%%%%%%;
% Now we can translate the gaussian in x_u, ;
% but we make sure to keep the support within the box. ;
% In other words, we make sure we can still integrate it to about 4-6 digits of accuracy. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Transforming shifted gaussian:')); end;
delta_x_c_max = half_diameter_x_c/2.0;
delta_x_c_ = randn(3,1); delta_x_c_ = delta_x_c_max*delta_x_c_/max(1e-12,fnorm(delta_x_c_));
h_k_p_form_ = g_k_p_form_ .* exp(-i*2*pi*( k_c_0_all_.*delta_x_c_(1+0) + k_c_1_all_.*delta_x_c_(1+1) + k_c_2_all_.*delta_x_c_(1+2) ));
h_k_p_form_l1 = sum(h_k_p_form_.*weight_3d_k_all_,'all');
h_k_p_full_l1 = exp(-(delta_x_c_(1+0).^2 + delta_x_c_(1+1).^2 + delta_x_c_(1+2).^2)/max(1e-12,2*sigma_x^2)) ;
h_k_p_form_l2 = sqrt(sum(abs(h_k_p_form_).^2.*weight_3d_k_all_,'all'));
h_k_p_full_l2 = sqrt( (sqrt(2*pi) * sigma_k * sqrt(2)).^(-3) );
h_x_u_form_ = nrm_x_from_k * 1/sqrt(2*pi)^3 .* 1/max(1e-12,sigma_x^3) .* exp(-((x_u_0___-delta_x_c_(1+0)).^2 + (x_u_1___-delta_x_c_(1+1)).^2 + (x_u_2___-delta_x_c_(1+2)).^2)/max(1e-12,2*sigma_x^2)) ;
h_x_u_form_l1 = sum(h_x_u_form_.*weight_xxx_u_,'all');
h_x_u_full_l1 = nrm_x_from_k*1.0;
h_x_u_form_l2 = sqrt(sum(abs(h_x_u_form_).^2.*weight_xxx_u_,'all'));
h_x_u_full_l2 = nrm_x_from_k * sqrt( (sqrt(2*pi) * sigma_x * sqrt(2)).^(-3) );
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l1 vs h_k_p_form_l1: %0.16f %%<-- only accurate if shifted gaussian decays within k-domain. ',fnorm(h_k_p_full_l1 - h_k_p_form_l1)/max(1e-12,fnorm(h_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l2 vs h_k_p_form_l2: %0.16f %%<-- should be accurate. ',fnorm(h_k_p_full_l2 - h_k_p_form_l2)/max(1e-12,fnorm(h_k_p_full_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_x_u_full_l1 vs h_x_u_form_l1: %0.16f',fnorm(h_x_u_full_l1 - h_x_u_form_l1)/max(1e-12,fnorm(h_x_u_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_x_u_full_l2 vs h_x_u_form_l2: %0.16f',fnorm(h_x_u_full_l2 - h_x_u_form_l2)/max(1e-12,fnorm(h_x_u_full_l2)))); end;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
h_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,h_x_u_form_(:).*weight_xxx_u_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: h_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
h_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,h_k_p_form_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: h_x_u_quad_ time %0.2fs',tmp_t)); end;
h_x_u_quad_ = reshape(h_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
h_k_p_quad_l1 = sum(h_k_p_quad_.*weight_3d_k_all_,'all');
h_k_p_quad_l2 = sqrt(sum(abs(h_k_p_quad_).^2.*weight_3d_k_all_,'all'));
h_x_u_quad_l1 = sum(h_x_u_quad_.*weight_xxx_u_,'all');
h_x_u_quad_l2 = sqrt(sum(abs(h_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l1 vs h_k_p_quad_l1: %0.16f %%<-- only accurate if shifted gaussian decays within k-domain. ',fnorm(h_k_p_full_l1 - h_k_p_quad_l1)/max(1e-12,fnorm(h_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l2 vs h_k_p_quad_l2: %0.16f %%<-- should be accurate. ',fnorm(h_k_p_full_l2 - h_k_p_quad_l2)/max(1e-12,fnorm(h_k_p_full_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_x_u_full_l1 vs h_x_u_quad_l1: %0.16f',fnorm(h_x_u_full_l1 - h_x_u_quad_l1)/max(1e-12,fnorm(h_x_u_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_x_u_full_l2 vs h_x_u_quad_l2: %0.16f',fnorm(h_x_u_full_l2 - h_x_u_quad_l2)/max(1e-12,fnorm(h_x_u_full_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;

if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
linewidth_use = 2;
markersize_use = 12;
p_row = 2; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
g_k_p_lim_ = max(abs(g_k_p_form_(:)))*[-1,+1];
plot(g_k_p_lim_,g_k_p_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(g_k_p_form_(:)),real(g_k_p_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(g_k_p_lim_); ylim(g_k_p_lim_); axisnotick; axis square;
xlabel('g_k_p_form_','Interpreter','none');
ylabel('g_k_p_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
g_x_u_lim_ = max(abs(g_x_u_form_(:)))*[-1,+1];
plot(g_x_u_lim_,g_x_u_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(g_x_u_form_(:)),real(g_x_u_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(g_x_u_lim_); ylim(g_x_u_lim_); axisnotick; axis square;
xlabel('g_x_u_form_','Interpreter','none');
ylabel('g_x_u_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
h_k_p_lim_ = max(abs(h_k_p_form_(:)))*[-1,+1];
plot(h_k_p_lim_,h_k_p_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(h_k_p_form_(:)),real(h_k_p_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(h_k_p_lim_); ylim(h_k_p_lim_); axisnotick; axis square;
xlabel('h_k_p_form_','Interpreter','none');
ylabel('h_k_p_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
h_x_u_lim_ = max(abs(h_x_u_form_(:)))*[-1,+1];
plot(h_x_u_lim_,h_x_u_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(h_x_u_form_(:)),real(h_x_u_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(h_x_u_lim_); ylim(h_x_u_lim_); axisnotick; axis square;
xlabel('h_x_u_form_','Interpreter','none');
ylabel('h_x_u_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
str_sgtitle = sprintf('k_eq_d_double %0.2f k_int %d n_x_u %d',k_eq_d_double,k_int,n_x_u);
sgtitle(str_sgtitle,'Interpreter','none');
drawnow();
end;%if flag_disp;

%%%%%%%%;
% Now generate a volume. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% generating arbitrary volume:')); end;
n_source = 128;
delta_x_c_3s__ = zeros(3,n_source);
for nsource=0:n_source-1;
%delta_x_c_ = randn(3,1); delta_x_c_ = delta_x_c_max*delta_x_c_/max(1e-12,fnorm(delta_x_c_)); %<-- random 
tmp_s = 2*pi*nsource/max(1,n_source-1); tmp_r = 2.0*nsource/max(1,n_source-1) - 1.0; tmp_r = 1.0;
tmp_delta_x_0 = delta_x_c_max*tmp_r*cos(2*tmp_s);
tmp_delta_x_1 = delta_x_c_max*tmp_r*sin(3*tmp_s);
tmp_delta_x_2 = delta_x_c_max*tmp_r*0.5*(cos(5*tmp_s) + sin(7*tmp_s));
delta_x_c_ = [tmp_delta_x_0;tmp_delta_x_1;tmp_delta_x_2];
delta_x_c_3s__(:,1+nsource) = delta_x_c_;
end;%for nsource=0:n_source-1;
j_x_u_form_ = zeros(n_x_u,n_x_u,n_x_u);
for nsource=0:n_source-1;
delta_x_c_ = delta_x_c_3s__(:,1+nsource);
j_x_u_form_ = j_x_u_form_ + nrm_x_from_k * 1/sqrt(2*pi)^3 .* 1/max(1e-12,sigma_x^3) .* exp(-((x_u_0___-delta_x_c_(1+0)).^2 + (x_u_1___-delta_x_c_(1+1)).^2 + (x_u_2___-delta_x_c_(1+2)).^2)/max(1e-12,2*sigma_x^2)) ;
end;%for nsource=0:n_source-1;
j_x_u_form_l1 = sum(j_x_u_form_.*weight_xxx_u_,'all');
j_x_u_form_l2 = sqrt(sum(abs(j_x_u_form_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
j_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,j_x_u_form_(:).*weight_xxx_u_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: j_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
j_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,j_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;;
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% xxnufft3d3: j_x_u_quad_ time %0.2fs',tmp_t)); end;
j_x_u_quad_ = reshape(j_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
j_x_u_quad_l1 = sum(j_x_u_quad_.*weight_xxx_u_,'all');
j_x_u_quad_l2 = sqrt(sum(abs(j_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% j_x_u_form_l1 vs j_x_u_quad_l1: %0.16f',fnorm(j_x_u_form_l1 - j_x_u_quad_l1)/max(1e-12,fnorm(j_x_u_form_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% j_x_u_form_l2 vs j_x_u_quad_l2: %0.16f',fnorm(j_x_u_form_l2 - j_x_u_quad_l2)/max(1e-12,fnorm(j_x_u_form_l2)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;

if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 16;
val_zoom = sqrt(2);
prctile_use = 98;
p_row = 1; p_col = 2; np=0;
j_x_u_bnd = max(abs(j_x_u_form_(:)));
parameter = struct('type','parameter');
parameter.vval_ = prctile(abs(j_x_u_form_),prctile_use,'all');
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
isosurface_f_x_u_1(parameter,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),j_x_u_form_));
axisnotick3d; axis equal; axis vis3d;
title('original');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
isosurface_f_x_u_1(parameter,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),j_x_u_quad_));
axisnotick3d; axis equal; axis vis3d;
title('reconstructed');
set(gca,'FontSize',fontsize_use);
%%%%;
str_sgtitle = sprintf('k_eq_d_double %0.2f k_int %d n_x_u %d',k_eq_d_double,k_int,n_x_u);
sgtitle(str_sgtitle,'Interpreter','none');
drawnow();
end;%if flag_disp;

%%%%%%%%;
% Note that, in the nufft1d3, we have: ;
% x_u_0___(:)*eta <-- ranges from pi*[-1:+1]. ;
% 2*pi*k_c_0_all_/eta <-- ranges from 2*x_p_r_max*k_p_r_max*[-1:+1]. ;
% The product [x_u_0___(:)*eta]*[2*pi*k_c_0_all_/eta] yields 2*pi*x_u_0*k_c_0, as desired. ;
%%%%;
% For the nufft1d2, we have: ;
% k0 <-- -n_x_u/2:+(n_x_u-1)/2. ;
% Thus, we need x[j] to range from: ;
% -2*pi*k_p_r_max*x_p_r_max/(n_x_u/2) up to +2*pi*k_p_r_max*x_p_r_max/((n_x_u-1)/2). ;
% This is impossible with an affine transformation, ;
% so instead we scale x[j] by 2*x_p_r_max/max(1,n_x_u-1) to produce x_tilde, ;
% and then multiply the result by the appropriate phase-factor (see below). ;
%%%%%%%%;
tmp_x_p_r = 1.25;
tmp_x_ = x_u_0_*tmp_x_p_r;
tmp_f_x_ = randn(n_x_u,1) + i*randn(n_x_u,1);
tmp_k_ = randn(8,1);
tmp_isgn = -1;
tmp_f_k_1d3_ = finufft1d3(tmp_x_,tmp_f_x_,tmp_isgn,1e-12,tmp_k_);
tmp_k_tilde_ = tmp_k_.*(2*tmp_x_p_r/max(1,n_x_u-1));
tmp_f_k_1d2_ = finufft1d2(tmp_k_tilde_,tmp_isgn,1e-12,tmp_f_x_).*exp(+i*tmp_isgn*tmp_k_tilde_/2); %<-- note the phase-factor. ;
fnorm_disp(flag_verbose,'tmp_f_k_1d3_',tmp_f_k_1d3_,'tmp_f_k_1d2_',tmp_f_k_1d2_);
%%%%%%%%;
% And to go backwards we do the reverse. ;
%%%%%%%%;
tmp_x_p_r = 1.25;
tmp_x_ = x_u_0_*tmp_x_p_r;
tmp_k_ = randn(8,1);
tmp_f_k_ = randn(8,1) + i*randn(8,1);
tmp_isgn = +1;
tmp_f_x_1d3_ = finufft1d3(tmp_k_,tmp_f_k_,tmp_isgn,1e-12,tmp_x_);
tmp_k_tilde_ = tmp_k_.*(2*tmp_x_p_r/max(1,n_x_u-1));
tmp_f_x_1d1_ = finufft1d1(tmp_k_tilde_,tmp_f_k_.*exp(+i*tmp_isgn*tmp_k_tilde_/2),tmp_isgn,1e-12,n_x_u); %<-- note the phase-factor. ;
fnorm_disp(flag_verbose,'tmp_f_x_1d3_',tmp_f_x_1d3_,'tmp_f_x_1d1_',tmp_f_x_1d1_);

%%%%%%%%;
% Now demonstrate that we can recover the same volume in fourier-space using: ;
% (i) a combination of the nufft3d2 and nufft3d1, instead of ;
% (ii) the nufft3d3. ;
%%%%;
% Note that, in the nufft3d3, we have: ;
% x_u_0___(:)*eta <-- ranges from pi*[-1:+1]. ;
% 2*pi*k_c_0_all_/eta <-- ranges from 2*x_p_r_max*k_p_r_max*[-1:+1]. ;
% The product [x_u_0___(:)*eta]*[2*pi*k_c_0_all_/eta] yields 2*pi*x_u_0*k_c_0, as desired. ;
%%%%;
% For the nufft3d2, we have: ;
% k0 <-- -n_x_u/2:+(n_x_u-1)/2 <-- -64:+63. ;
% Thus, we need x[j] to range from -2*pi*k_p_r_max*x_p_r_max/(n_x_u/2) up to +2*pi*k_p_r_max*x_p_r_max/((n_x_u-1)/2). 
%%%%%%%%;
k_tilde_c_0_all_ = 2*pi*k_c_0_all_.*(2*x_p_r_max/max(1,n_x_u-1));
k_tilde_c_1_all_ = 2*pi*k_c_1_all_.*(2*x_p_r_max/max(1,n_x_u-1));
k_tilde_c_2_all_ = 2*pi*k_c_2_all_.*(2*x_p_r_max/max(1,n_x_u-1));
%%%%;
tmp_isgn = -1;
tmp_t = tic;
j_k_p_q3d2_ = finufft3d2(k_tilde_c_0_all_,k_tilde_c_1_all_,k_tilde_c_2_all_,tmp_isgn,1e-12,reshape(j_x_u_form_(:).*weight_xxx_u_(:),[n_x_u,n_x_u,n_x_u])).*exp(+i*tmp_isgn*(k_tilde_c_0_all_+k_tilde_c_1_all_+k_tilde_c_2_all_)/2)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3);
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% finufft3d2: j_k_p_q3d2_ time %0.2fs',tmp_t)); end;
fnorm_disp(flag_verbose,'j_k_p_quad_',j_k_p_quad_,'j_k_p_q3d2_',j_k_p_q3d2_);
%%%%;
tmp_isgn = +1;
tmp_t = tic;
j_x_u_q3d1_ = finufft3d1(k_tilde_c_0_all_,k_tilde_c_1_all_,k_tilde_c_2_all_,j_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_.*exp(+i*tmp_isgn*(k_tilde_c_0_all_+k_tilde_c_1_all_+k_tilde_c_2_all_)/2),tmp_isgn,1e-12,n_x_u,n_x_u,n_x_u)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3);
if (flag_verbose>0); tmp_t = toc(tmp_t); disp(sprintf(' %% %% finufft3d1: j_x_u_q3d1_ time %0.2fs',tmp_t)); end;
j_x_u_q3d1_ = reshape(j_x_u_q3d1_,[n_x_u,n_x_u,n_x_u]);
fnorm_disp(flag_verbose,'j_x_u_quad_',j_x_u_quad_,'j_x_u_q3d1_',j_x_u_q3d1_);
%%%%;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;



