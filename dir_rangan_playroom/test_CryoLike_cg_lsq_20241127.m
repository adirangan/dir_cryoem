function nf = test_CryoLike_cg_lsq_20241127(k_eq_d_double,k_int,n_x_u,nf);
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
% defaults: k_eq_d_double = 1.0; k_int = 48; n_x_u = 64; nf=0;

str_thisfunction = 'test_CryoLike_cg_lsq_20241127';

if (nargin<1);
k_eq_d_double_ = [1.0,4.0];
n_k_eq_d_double = numel(k_eq_d_double_);
nf=0;
for nk_eq_d_double=0:n_k_eq_d_double-1;
k_eq_d_double = k_eq_d_double_(1+nk_eq_d_double);
nf = test_CryoLike_cg_lsq_20241127(k_eq_d_double,[],[],nf);
end;%for nk_eq_d_double=0:n_k_eq_d_double-1;
disp(sprintf(' %% returning')); return;
end;%if (nargin<1);

na=0;
if nargin<1+na; k_eq_d_double = []; end; na=na+1;
if nargin<1+na; k_int = []; end; na=na+1;
if nargin<1+na; n_x_u = []; end; na=na+1;
if nargin<1+na; nf = []; end; na=na+1;

% defaults for testing: ;
% clear; str_thisfunction = 'test_CryoLike_cg_lsq_20241127'; k_eq_d_double = []; k_int = []; n_x_u = []; nf = [];

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
if isempty(k_int); k_int = 48; end;
if isempty(n_x_u); n_x_u = 64; end;
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
fnorm_disp(flag_verbose,'I_a_form',I_a_form,'I_a_quad',I_a_quad);
fnorm_disp(flag_verbose,'I_b_form',I_b_form,'I_b_quad',I_b_quad);
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
fnorm_disp(flag_verbose,'g_k_p_full_l1',g_k_p_full_l1,'g_k_p_form_l1',g_k_p_form_l1);
fnorm_disp(flag_verbose,'g_k_p_full_l2',g_k_p_full_l2,'g_k_p_form_l2',g_k_p_form_l2);
fnorm_disp(flag_verbose,'g_x_u_full_l1',g_x_u_full_l1,'g_x_u_form_l1',g_x_u_form_l1);
fnorm_disp(flag_verbose,'g_x_u_full_l2',g_x_u_full_l2,'g_x_u_form_l2',g_x_u_form_l2);
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
g_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,g_x_u_form_(:).*weight_xxx_u_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: g_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
g_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,g_k_p_form_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: g_x_u_quad_ time %0.2fs',tmp_t)); end;
g_x_u_quad_ = reshape(g_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
g_k_p_quad_l1 = sum(g_k_p_quad_.*weight_3d_k_all_,'all');
g_k_p_quad_l2 = sqrt(sum(abs(g_k_p_quad_).^2.*weight_3d_k_all_,'all'));
g_x_u_quad_l1 = sum(g_x_u_quad_.*weight_xxx_u_,'all');
g_x_u_quad_l2 = sqrt(sum(abs(g_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
fnorm_disp(flag_verbose,'g_k_p_full_l1',g_k_p_full_l1,'g_k_p_quad_l1',g_k_p_quad_l1);
fnorm_disp(flag_verbose,'g_k_p_full_l2',g_k_p_full_l2,'g_k_p_quad_l2',g_k_p_quad_l2);
fnorm_disp(flag_verbose,'g_x_u_full_l1',g_x_u_full_l1,'g_x_u_quad_l1',g_x_u_quad_l1);
fnorm_disp(flag_verbose,'g_x_u_full_l2',g_x_u_full_l2,'g_x_u_quad_l2',g_x_u_quad_l2);

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
fnorm_disp(flag_verbose,'h_k_p_full_l1',h_k_p_full_l1,'h_k_p_form_l1',h_k_p_form_l1);
fnorm_disp(flag_verbose,'h_k_p_full_l2',h_k_p_full_l2,'h_k_p_form_l2',h_k_p_form_l2);
fnorm_disp(flag_verbose,'h_x_u_full_l1',h_x_u_full_l1,'h_x_u_form_l1',h_x_u_form_l1);
fnorm_disp(flag_verbose,'h_x_u_full_l2',h_x_u_full_l2,'h_x_u_form_l2',h_x_u_form_l2);
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
h_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,h_x_u_form_(:).*weight_xxx_u_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: h_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
h_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,h_k_p_form_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: h_x_u_quad_ time %0.2fs',tmp_t)); end;
h_x_u_quad_ = reshape(h_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
h_k_p_quad_l1 = sum(h_k_p_quad_.*weight_3d_k_all_,'all');
h_k_p_quad_l2 = sqrt(sum(abs(h_k_p_quad_).^2.*weight_3d_k_all_,'all'));
h_x_u_quad_l1 = sum(h_x_u_quad_.*weight_xxx_u_,'all');
h_x_u_quad_l2 = sqrt(sum(abs(h_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l1 vs h_k_p_quad_l1: %0.16f %%<-- only accurate if shifted gaussian decays within k-domain. ',fnorm(h_k_p_full_l1 - h_k_p_quad_l1)/max(1e-12,fnorm(h_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l2 vs h_k_p_quad_l2: %0.16f %%<-- should be accurate. ',fnorm(h_k_p_full_l2 - h_k_p_quad_l2)/max(1e-12,fnorm(h_k_p_full_l2)))); end;
fnorm_disp(flag_verbose,'h_k_p_full_l1',h_k_p_full_l1,'h_k_p_quad_l1',h_k_p_quad_l1);
fnorm_disp(flag_verbose,'h_k_p_full_l2',h_k_p_full_l2,'h_k_p_quad_l2',h_k_p_quad_l2);
fnorm_disp(flag_verbose,'h_x_u_full_l1',h_x_u_full_l1,'h_x_u_quad_l1',h_x_u_quad_l1);
fnorm_disp(flag_verbose,'h_x_u_full_l2',h_x_u_full_l2,'h_x_u_quad_l2',h_x_u_quad_l2);
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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: j_k_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
j_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,j_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: j_x_u_quad_ time %0.2fs',tmp_t)); end;
j_x_u_quad_ = reshape(j_x_u_quad_,[n_x_u,n_x_u,n_x_u]);
%%%%;
j_x_u_quad_l1 = sum(j_x_u_quad_.*weight_xxx_u_,'all');
j_x_u_quad_l2 = sqrt(sum(abs(j_x_u_quad_).^2.*weight_xxx_u_,'all'));
%%%%%%%%;
fnorm_disp(flag_verbose,'j_x_u_form_l1',j_x_u_form_l1,'j_x_u_quad_l1',j_x_u_quad_l1);
fnorm_disp(flag_verbose,'j_x_u_form_l2',j_x_u_form_l2,'j_x_u_quad_l2',j_x_u_quad_l2);

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
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% finufft3d2: j_k_p_q3d2_ time %0.2fs',tmp_t)); end;
fnorm_disp(flag_verbose,'j_k_p_quad_',j_k_p_quad_,'j_k_p_q3d2_',j_k_p_q3d2_);
%%%%;
tmp_isgn = +1;
tmp_t = tic;
j_x_u_q3d1_ = finufft3d1(k_tilde_c_0_all_,k_tilde_c_1_all_,k_tilde_c_2_all_,j_k_p_quad_.*(2*pi)^3.*weight_3d_k_all_.*exp(+i*tmp_isgn*(k_tilde_c_0_all_+k_tilde_c_1_all_+k_tilde_c_2_all_)/2),tmp_isgn,1e-12,n_x_u,n_x_u,n_x_u)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% finufft3d1: j_x_u_q3d1_ time %0.2fs',tmp_t)); end;
j_x_u_q3d1_ = reshape(j_x_u_q3d1_,[n_x_u,n_x_u,n_x_u]);
fnorm_disp(flag_verbose,'j_x_u_quad_',j_x_u_quad_,'j_x_u_q3d1_',j_x_u_q3d1_);
%%%%%%%%;

%%%%%%%%;
% Now generate some templates. ;
%%%%%%%%;
n_S = 1024;
n_w_max = 2+2*round(k_p_r_max/max(1e-12,k_eq_d));
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_sum = n_w_max*n_k_p_r;
n_w_csum_ = cumsum([0;n_w_]);
viewing_polar_a_S_ = 1*pi*rand(n_S,1);
viewing_azimu_b_S_ = 2*pi*rand(n_S,1);
viewing_gamma_z_S_ = 2*pi*rand(n_S,1);
%%%%%%%%;
% Note that the convention in cg_rhs_2 is to ;
% subtract viewing_gamma_z from inplane_gamma_z. ;
% This may not be consistent with _fourier_circles. ;
%%%%%%%%;
[ ...
 k_p_polar_a_wS__ ...
,k_p_azimu_b_wS__ ...
,k_c_0_wS__ ...
,k_c_1_wS__ ...
,k_c_2_wS__ ...
] = ...
cg_rhs_2( ...
 n_S ...
,n_w_max ...
,viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_gamma_z_S_ ...
);
%%%%%%%%;
k_c_0_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_0_wS__,[n_w_max,1,n_S]));
k_c_1_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_1_wS__,[n_w_max,1,n_S]));
k_c_2_wkS___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_2_wS__,[n_w_max,1,n_S]));
%%%%%%%%;
% Now we construct the forward- and adjoint-evaluation-operator. ;
% Note that these produce the 'unscaled' templates. ;
% That is, they do not yet account for the 2d-weight. ;
% As such, the variance per-degree-of-freedom is not uniform, ;
% and these are not immediately appropriate for a least-squares problem. ;
%%%%%%%%;
A_forward__ = @(a_x_u_test_) ...
reshape( xxnufft3d3(n_xxx_u,x_u_0___(:)*(pi/x_p_r_max),x_u_1___(:)*(pi/x_p_r_max),x_u_2___(:)*(pi/x_p_r_max),a_x_u_test_(:).*weight_xxx_u_(:),-1,1e-12,n_w_max*n_k_p_r*n_S,2*pi*k_c_0_wkS___(:)/(pi/x_p_r_max),2*pi*k_c_1_wkS___(:)/(pi/x_p_r_max),2*pi*k_c_2_wkS___(:)/(pi/x_p_r_max))/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) , [n_w_sum,n_S] ) ;
%%%%;
A_adjoint__ = @(a_k_p_test_) ...
reshape( xxnufft3d3(n_k_all,2*pi*k_c_0_wkS___(:)*(pi/k_p_r_max),2*pi*k_c_1_wkS___(:)*(pi/k_p_r_max),2*pi*k_c_2_wkS___(:)*(pi/k_p_r_max),reshape(a_k_p_test_,[n_w_sum*n_S,1]),+1,1e-12,n_xxx_u,x_u_0___(:)/(pi/k_p_r_max),x_u_1___(:)/(pi/k_p_r_max),x_u_2___(:)/(pi/k_p_r_max))/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) * weight_xxx_u_(:) , [n_x_u,n_x_u,n_x_u] );
%%%%;
% Now we generate the templates: ;
%%%%;
tmp_t = tic;
S_k_p_wkS__ = A_forward__(j_x_u_form_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% finufft3d3: S_k_p_wkS__ time %0.2fs',tmp_t)); end;
%%%%;
% Now we generate a selection of synthetic-images, copied from the templates, ;
% as well as a collection of CTF-functions. ;
%%%%;
n_M = n_S;
n_CTF = 8; %<-- some number of CTF-functions. ;
CTF_k_p_r_kC__ = zeros(n_k_p_r,n_CTF);
for nCTF=0:n_CTF-1;
CTF_k_p_r_k_ = cos(nCTF*2*pi*k_p_r_/max(1e-12,k_p_r_max));
CTF_k_p_r_kC__(:,1+nCTF) = CTF_k_p_r_k_;
end;%for nCTF=0:n_CTF-1;
index_nCTF_from_nM_ = mod(transpose([0:n_M-1]),n_CTF);
M_k_p_wkM__ = S_k_p_wkS__;
for nM=0:n_M-1;
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
M_k_p_wk__ = reshape(M_k_p_wk_,[n_w_max,n_k_p_r]);
M_k_p_wk__ = bsxfun(@times,M_k_p_wk__,reshape(CTF_k_p_r_k_,[1,n_k_p_r]));
M_k_p_wk_ = reshape(M_k_p_wk__,[n_w_sum,1]);
M_k_p_wkM__(:,1+nM) = M_k_p_wk_;
end;%for nM=0:n_M-1;
%%%%;

if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 6; np=0;
for nM=0:min(p_row*p_col/3,min(n_M,n_CTF))-1;
nS = nM;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS);
Slim_ = prctile(real(S_k_p_wk_),75,'all'); Slim_ = 1.25*Slim_*[-1,+1];
M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
nCTF = index_nCTF_from_nM_(1+nM);
CTF_k_p_r_k_ = CTF_k_p_r_kC__(:,1+nCTF);
subplot(p_row,p_col,1+np);np=np+1;cla;
plot(k_p_r_,CTF_k_p_r_k_);
xlim([0,k_p_r_max]); xlabel('k_p_r_','Interpreter','none');
ylim([-1,+1]); ylabel('CTF_k_p_r_','Interpreter','none');
title(sprintf('nM %d: CTF',nM),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),Slim_);
axis image; axisnotick;
title(sprintf('nM %d: S_k_p_wk_',nM),'Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_wk_),Slim_);
axis image; axisnotick;
title(sprintf('nM %d: M_k_p_wk_',nM),'Interpreter','none');
end;%for nM=0:min(n_M,p_row)-1;
end;%if flag_disp;

%%%%%%%%;
% Now we perform a quick and dirty back-propagation. ;
% This will serve as an initial-guess for the least-squares below. ;
% Note that the implementation below uses a uniform cartesian-grid in frequency-space ;
% (i.e., k_u_0___, k_u_1___ and k_u_2___). ;
% We could certainly apply the same strategy for a spherical-grid ;
% as long as we define the quadrature-weights appropriately. ;
%%%%%%%%;
n_k_u = n_x_u;
half_diameter_k_c = k_p_r_max;
diameter_k_c = 2.0d0*half_diameter_k_c;
k_u_0_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_u)); d_k_0 = mean(diff(k_u_0_));
k_u_1_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_u)); d_k_1 = mean(diff(k_u_1_));
k_u_2_ = transpose(linspace(-k_p_r_max,+k_p_r_max,n_k_u)); d_k_2 = mean(diff(k_u_2_));
[k_u_0___,k_u_1___,k_u_2___] = ndgrid(k_u_0_,k_u_1_,k_u_2_); n_kkk_u = n_k_u^3;
weight_kkk_u_ = d_k_0 * d_k_1 * d_k_2 ;
index_k_u_0_wkS___ = max(0,min(n_k_u-1,round((k_c_0_wkS___+half_diameter_k_c)/max(1e-12,d_k_0))));
index_k_u_1_wkS___ = max(0,min(n_k_u-1,round((k_c_1_wkS___+half_diameter_k_c)/max(1e-12,d_k_1))));
index_k_u_2_wkS___ = max(0,min(n_k_u-1,round((k_c_2_wkS___+half_diameter_k_c)/max(1e-12,d_k_2))));
index_k_u_wkS___ = index_k_u_0_wkS___ + (index_k_u_1_wkS___ + index_k_u_2_wkS___*n_k_u)*n_k_u ;
CTF_k_p_wkM__ = reshape(repmat(reshape(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),[1,n_k_p_r,n_M]),[n_w_max,1,1]),[n_w_sum,n_M]);
M0C2_k_u_ = sparse(1+index_k_u_wkS___(:),1+0,abs(CTF_k_p_wkM__(:)).^2,n_kkk_u,1);
M_upb = max(abs(M_k_p_wkM__),[],'all'); %<-- temporarily rescale so that division threshold applies. ;
M1C1_k_u_ = sparse(1+index_k_u_wkS___(:),1+0,CTF_k_p_wkM__(:).*M_k_p_wkM__(:)/max(1e-12,M_upb),n_kkk_u,1);
j_k_u_qbpu_ = M_upb*full(M1C1_k_u_./max(1e-12,M0C2_k_u_));
eta = pi/k_p_r_max; tmp_t = tic;
j_x_u_qbpu_ = xxnufft3d3(n_kkk_u,2*pi*k_u_0___(:)*eta,2*pi*k_u_1___(:)*eta,2*pi*k_u_2___(:)*eta,j_k_u_qbpu_.*weight_kkk_u_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: j_x_u_qbpu_ time %0.2fs',tmp_t)); end;
j_x_u_qbpu_ = reshape(j_x_u_qbpu_,[n_x_u,n_x_u,n_x_u]);

%%%%%%%%;
% Here we quickly test the least-squares algorithm. ;
%%%%%%%%;
tmp_n = 64;
tmp_A_n__ = randn(tmp_n+0,tmp_n+0);
tmp_A_t__ = transpose(tmp_A_n__);
lhs_tru_ = rand(tmp_n,1);
b_tru_ = tmp_A_n__*lhs_tru_;
LSQ_forward__ = @(x_) tmp_A_n__*x_;
LSQ_adjoint__ = @(b_) tmp_A_t__*b_;
LSQ__ = @(x_) LSQ_adjoint__(LSQ_forward__(x_));
RHS_form_ = LSQ__(lhs_tru_);
RHS_ = LSQ_adjoint__(b_tru_);
fnorm_disp(flag_verbose,'RHS_form_',RHS_form_,'RHS_',RHS_);
flag_verbose = 1;
flag_disp = 1;
tolerance_cg_lsq = 1e-4;
n_iteration = tmp_n;
niteration=0;
lhs_tmp_ = zeros(tmp_n,1);
res_ = RHS_;
pcg_ = res_;
beta_num = sum(res_.^2);
beta_den = 1.0;
beta = beta_num/max(1e-12,beta_den);
flag_continue = 1;
while flag_continue
if flag_verbose; fprintf(' %% niteration = %.3d beta_num = %0.6f\n', niteration, beta_num); end;
An_pcg_ = LSQ_forward__(pcg_);
zeta = sum(abs(An_pcg_).^2,'all');
alph = beta_num / max(1e-12, zeta);
lhs_tmp_ = lhs_tmp_ + alph*pcg_ ;
res_ = res_ - alph * LSQ_adjoint__(An_pcg_);
beta_den = beta_num;
beta_num = sum(res_.^2);
beta = beta_num / max(1e-12,beta_den);
pcg_ = res_ + beta * pcg_ ;
niteration = niteration + 1;
if niteration >= n_iteration; flag_continue=0; end;
if sqrt(beta) < tolerance_cg_lsq; flag_continue=0; end;
end;%while;
fnorm_disp(flag_verbose,'lhs_tru_',lhs_tru_,'lhs_tmp_',lhs_tmp_);
lhs_pcg_ = pcg(LSQ__,RHS_,tolerance_cg_lsq,n_iteration); %<-- this is the prepackaged matlab version. ;
fnorm_disp(flag_verbose,'lhs_tru_',lhs_tru_,'lhs_pcg_',lhs_pcg_);

%%%%%%%%;
% Note that, typically, the distribution of points within each template can be nonuniform. ;
% To account for this we need the template-weights. ;
% Here the convention is that: ;
% sum(weight_2d_k_p_r_) = pi*k_p_r_max^2. ;
% sum(weight_2d_wk_)*(2*pi)^2 = pi*k_p_r_max^2. ;
% We will use these further below. ;
%%%%%%%%;
[ ...
 ~ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_ ...
);
%%%%%%%%;

%%%%%%%%;
% Now we incorporate the weights into the least-squares problem. ;
% Now the forward-operator is designed so that: ;
% LSQ_forward__ * a_x_u_test_ \approx bsxfun(@times,M_k_p_wkM___,reshape(sqrt(weight_2d_k_p_r_),[1,n_k_p_r,1])), ;
% precisely so that the residual is: ;
% R_k_p_wkM__ = bsxfun(@times,T_k_p_wkM___ - M_k_p_wkM___,reshape(sqrt(weight_2d_k_p_r_),[1,n_k_p_r,1])), ;
% and the standard l2-norm of the residual is: ;
% sum(abs(R_k_p_wkM__).^2,'all') = sum(bsxfun(@times,abs(T_k_p_wkM___ - M_k_p_wkM___).^2,reshape(weight_2d_k_p_r_,[1,n_k_p_r,1])),'all'), ;
% which is the sigma-squared-negative-log-likelihood (ssnll) in the low-temperature limit. ;
%%%%%%%%;
LSQ_forward__ = @(a_x_u_test_) reshape(bsxfun(@times,bsxfun(@times,reshape(A_forward__(reshape(a_x_u_test_,[n_x_u,n_x_u,n_x_u])),[n_w_max,n_k_p_r,n_M]),reshape(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),[1,n_k_p_r,n_M])),reshape(sqrt(weight_2d_k_p_r_),[1,n_k_p_r,1])),[n_w_sum*n_M,1]);
LSQ_adjoint__ = @(a_k_p_test_) reshape(A_adjoint__(reshape(bsxfun(@times,bsxfun(@times,reshape(a_k_p_test_,[n_w_max,n_k_p_r,n_M]),reshape(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),[1,n_k_p_r,n_M])),reshape(sqrt(weight_2d_k_p_r_),[1,n_k_p_r,1])),[n_w_sum,n_M])),[n_x_u*n_x_u*n_x_u,1]);
%%%%%%%%;
% Now we set up the normal-equations for the least-squares problem above: ;
%%%%%%%%;
LSQ__ = @(a_x_u_test_) LSQ_adjoint__(LSQ_forward__(a_x_u_test_));
tmp_t = tic();
RHS_form_ = LSQ__(j_x_u_form_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% RHS_ time %0.2fs',tmp_t)); end;
tmp_t = tic();
RHS_ = LSQ_adjoint__(reshape(bsxfun(@times,reshape(M_k_p_wkM__,[n_w_max,n_k_p_r,n_M]),reshape(sqrt(weight_2d_k_p_r_),[1,n_k_p_r,1])),[n_w_sum*n_M,1]));
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% RHS_ time %0.2fs',tmp_t)); end;
fnorm_disp(flag_verbose,'RHS_form_',RHS_form_,'RHS_',RHS_);
%%%%%%%%;
% Now we finally run the conjugate-gradient algorithm. ;
%%%%%%%%;
% Conjugate Gradient Loop
flag_verbose = 1;
flag_disp = 1;
tolerance_cg_lsq = 1e-4;
n_iteration = 20;
niteration=0;
lhs_tmp_ = zeros(n_x_u*n_x_u*n_x_u,1); %<-- default. ;
lhs_tmp_ = j_x_u_qbpu_(:); %<-- initial guess from quadrature-back-propagation. ;
tmp_t = tic();
res_ = RHS_;
if fnorm(lhs_tmp_)>tolerance_cg_lsq; res_ = RHS_ - LSQ__(lhs_tmp_); end;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% precomputation: time %0.2fs',tmp_t)); end;
pcg_ = res_;
beta_num = sum(res_.^2);
beta_den = 1.0;
beta = beta_num/max(1e-12,beta_den);
%%%%;
if flag_disp;
figure(1+nf);clf;figbig;
val_zoom = sqrt(2);
prctile_use = 98;
j_x_u_bnd = max(abs(j_x_u_form_(:)));
p_row = 1; p_col = 2;
parameter = struct('type','parameter');
parameter.vval_ = prctile(abs(j_x_u_form_),prctile_use,'all');
subplot(p_row,p_col,1); cla;
isosurface_f_x_u_1(parameter,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),j_x_u_form_));
axisnotick3d; axis equal; axis vis3d;
title('original');
end;%if flag_disp;
%%%%;
flag_continue = 1;
while flag_continue
if flag_verbose; fprintf(' %% niteration = %.3d beta_num = %0.6f\n', niteration, beta_num); end;
tmp_t = tic();
An_pcg_ = LSQ_forward__(pcg_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% An_pcg_: time %0.2fs',tmp_t)); end;
zeta = sum(abs(An_pcg_).^2,'all');
alph = beta_num / max(1e-12, zeta);
lhs_tmp_ = lhs_tmp_ + alph*pcg_ ;
tmp_t = tic();
res_ = res_ - alph * LSQ_adjoint__(An_pcg_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% At_An_pcg_: time %0.2fs',tmp_t)); end;
beta_den = beta_num;
beta_num = sum(res_.^2);
beta = beta_num / max(1e-12,beta_den);
pcg_ = res_ + beta * pcg_ ;
niteration = niteration + 1;
if niteration >= n_iteration; flag_continue=0; end;
if sqrt(beta) < tolerance_cg_lsq; flag_continue=0; end;
%%%%;
if flag_disp;
figure(1+nf);
subplot(p_row,p_col,2); cla;
parameter_tmp = struct('type','parameter');
parameter_tmp.vval_ = prctile(abs(lhs_tmp_),prctile_use,'all');
isosurface_f_x_u_1(parameter_tmp,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),reshape(real(lhs_tmp_),[n_x_u,n_x_u,n_x_u])));
axisnotick3d; axis equal; axis vis3d;
title(sprintf('niteration %.3d beta_num %0.6f',niteration,beta_num),'Interpreter','none');
drawnow();
end;%if flag_disp;
%%%%;
end;%while;

if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;




