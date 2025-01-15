function nf = test_CryoLike_fft_template_20241210(k_eq_d_double,k_int,n_x_u,nf);
%%%%%%%%;
% Tests template-generation from a_x_u_. ;
%%%%%%%%;
% defaults: k_eq_d_double = 1.0; k_int = 48; n_x_u = 64; nf=0;

str_thisfunction = 'test_CryoLike_fft_template_20241210';

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
% clear; str_thisfunction = 'test_CryoLike_fft_template_20241210'; k_eq_d_double = []; k_int = []; n_x_u = []; nf = [];

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
title('reconstructed (quad)');
set(gca,'FontSize',fontsize_use);
%%%%;
str_sgtitle = sprintf('k_eq_d_double %0.2f k_int %d n_x_u %d',k_eq_d_double,k_int,n_x_u);
sgtitle(str_sgtitle,'Interpreter','none');
drawnow();
end;%if flag_disp;

%%%%%%%%;
% Now convert to j_k_Y_ ; 
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
if ~exist('Ylm_uklma___','var'); Ylm_uklma___=[]; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__=[]; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__=[]; end;
if ~exist('l_max_uk_','var'); l_max_uk_=[]; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_=[]; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__=[]; end;
[ ...
 j_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_k_p_to_spharm_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,j_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% j_k_p_quad_ --> j_k_Y_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
tmp_t = tic;
[ ...
 j_k_p_reco_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
] = ...
convert_spharm_to_k_p_4( ...
 0*flag_verbose ...
,n_k_all ...
,n_k_all_csum_ ...
,k_p_r_all_ ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_3d_k_all_ ...
,weight_shell_k_ ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,j_k_Y_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% j_k_Y_quad_ --> j_k_p_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'j_k_p_quad_',j_k_p_quad_,'j_k_p_reco_',j_k_p_reco_);
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
j_x_u_reco_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,j_k_p_reco_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^3) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft3d3: j_x_u_reco_ time %0.2fs',tmp_t)); end;
j_x_u_reco_ = reshape(j_x_u_reco_,[n_x_u,n_x_u,n_x_u]);
%%%%%%%%;
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
isosurface_f_x_u_1(parameter,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),j_x_u_reco_));
axisnotick3d; axis equal; axis vis3d;
title('reconstructed (reco)');
set(gca,'FontSize',fontsize_use);
%%%%;
str_sgtitle = sprintf('k_eq_d_double %0.2f k_int %d n_x_u %d',k_eq_d_double,k_int,n_x_u);
sgtitle(str_sgtitle,'Interpreter','none');
drawnow();
end;%if flag_disp;
%%%%%%%%;

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
n_w_max = 2+2*round(k_p_r_max/max(1e-12,k_eq_d));
viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
template_k_eq_d = -1;
tmp_t = tic;
[ ...
 S_0_k_p_wkS___ ...
,n_w_max ...
,n_viewing_all ...
,viewing_azimu_b_all_ ...
,viewing_polar_a_all_ ...
,viewing_weight_all_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,local_yk__from_yk_(n_k_p_r,l_max_,j_k_Y_quad_) ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max ...
);
n_S = n_viewing_all;
n_w_ = n_w_max*ones(n_k_p_r,1);
n_w_max = max(n_w_);
n_w_sum = n_w_max*n_k_p_r;
n_w_csum_ = cumsum([0;n_w_]);
S_0_k_p_wkS__ = reshape(S_0_k_p_wkS___,[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% pm_template_2: S_0_k_p_wkS__ time %0.2fs',tmp_t)); end;
%%%%;
% copy parameters. ;
%%%%;
viewing_polar_a_S_ = viewing_polar_a_all_;
viewing_azimu_b_S_ = viewing_azimu_b_all_;
viewing_gamma_z_S_ = 2*pi*rand(n_S,1);
%%%%%%%%;
% Note that the convention in cg_rhs_2 is to ;
% subtract viewing_gamma_z from inplane_gamma_z. ;
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
% Now we construct the forward-evaluation-operator. ;
%%%%%%%%;
A_forward__ = @(a_x_u_test_) ...
reshape( xxnufft3d3(n_xxx_u,x_u_0___(:)*(pi/x_p_r_max),x_u_1___(:)*(pi/x_p_r_max),x_u_2___(:)*(pi/x_p_r_max),a_x_u_test_(:).*weight_xxx_u_(:),-1,1e-12,n_w_max*n_k_p_r*n_S,2*pi*k_c_0_wkS___(:)/(pi/x_p_r_max),2*pi*k_c_1_wkS___(:)/(pi/x_p_r_max),2*pi*k_c_2_wkS___(:)/(pi/x_p_r_max))/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^3) , [n_w_sum,n_S] ) ;
%%%%;
tmp_t = tic;
S_1_k_p_wkS__ = A_forward__(j_x_u_form_);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% finufft3d3: S_1_k_p_wkS__ time %0.2fs',tmp_t)); end;
%%%%;
S_2_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_2_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_1_k_p_wkS__(:,1+nS),-viewing_gamma_z_S_(1+nS));
end;%for nS=0:n_S-1;
%%%%;
fnorm_disp(flag_verbose,'S_0_k_p_wkS__',S_0_k_p_wkS__,'S_2_k_p_wkS__',S_2_k_p_wkS__);

if flag_disp;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 4; p_col = 6; np=0;
for np=0:min(p_row*p_col/2,n_S)-1;
nS = max(0,min(n_S-1,round(n_S*np/max(1,p_row*p_col/2))));
S_0_k_p_wk_ = S_0_k_p_wkS__(:,1+nS);
Slim_ = prctile(real(S_0_k_p_wk_),75,'all'); Slim_ = 1.25*Slim_*[-1,+1];
S_2_k_p_wk_ = S_2_k_p_wkS__(:,1+nS);
subplot(p_row,p_col,1+2*np+0);cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_0_k_p_wk_),Slim_);
axis image; axisnotick;
title(sprintf('nS %d: S_0_k_p_wk_',nS),'Interpreter','none');
subplot(p_row,p_col,1+2*np+1);cla;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_2_k_p_wk_),Slim_);
axis image; axisnotick;
title(sprintf('nS %d: S_2_k_p_wk_',nS),'Interpreter','none');
np=np+1;
end;%for nS=0:min(p_row*p_col,n_S)-1;
end;%if flag_disp;

%%%%%%%%;
% Note that, for the standard 'padded' fft, we have: ;
% f_x_ori_ = array of size N_ori. ;
% x_ori_ ranges from -tmp_x_p_r to +tmp_x_p_r --> x_ori_(1+nj) = -tmp_x_p_r + nj*2*tmp_x_p_r/(N_ori-1). ;
% 2pik_pad_ = 2*pi*[0:N_pad-1]/(2*tmp_x_p_r)*(N_ori-1)/(N_pad);
% f_k_pad_ = xfft(f_x_ori_,N_pad).*exp(-i*tmp_isgn*2pik_pad_*tmp_x_p_r) = array of size N_pad. ;
% Note the extra phase-factor for f_k_pad_. ;
% Note also that 2pik_pad_ can be periodized. ;
% With this convention d2pik_pad = (2*pi)/(2*tmp_x_p_r)*(n_x_u-1)/(n_x_u_pad). ;
% Thus, 2pik_pad_ ranges from 2*pi/(4*tmp_x_p_r)*(N_ori-1)*[-1,+1]. ;
% This implies that, in order for this strategy to work, ;
% 2*pi*k_p_r_max must be less than 2*pi/(4*x_p_r_max)*(n_x_u-1). ;
% i.e., k_p_r_max must be less than (n_x_u-1)/(4*x_p_r_max). ;
%%%%%%%%;
tmp_x_p_r = 1.25;
tmp_x_ = x_u_0_*tmp_x_p_r;
tmp_f_x_ = randn(n_x_u,1) + i*randn(n_x_u,1);
n_x_u_pad = 3*n_x_u;
tmp_2pik_fft1d0_ = transpose(periodize([0:n_x_u_pad-1],-n_x_u_pad/2,+n_x_u_pad/2))*(2*pi)/(2*tmp_x_p_r)*(n_x_u-1)/(n_x_u_pad);
for tmp_isgn = [-1,+1];
if tmp_isgn == -1; tmp_f_k_fft_ =       1.0* fft(tmp_f_x_,n_x_u_pad).*exp(-i*tmp_isgn*tmp_2pik_fft1d0_*tmp_x_p_r); end;
if tmp_isgn == +1; tmp_f_k_fft_ = n_x_u_pad*ifft(tmp_f_x_,n_x_u_pad).*exp(-i*tmp_isgn*tmp_2pik_fft1d0_*tmp_x_p_r); end;
tmp_2pik_fft1d3_ = tmp_2pik_fft1d0_;
tmp_f_k_1d3_ = finufft1d3(tmp_x_,tmp_f_x_,tmp_isgn,1e-12,tmp_2pik_fft1d3_);
tmp_2pik_tilde_ = tmp_2pik_fft1d3_.*(2*tmp_x_p_r/max(1,n_x_u-1));
tmp_f_k_1d2_ = finufft1d2(tmp_2pik_tilde_,tmp_isgn,1e-12,tmp_f_x_).*exp(+i*tmp_isgn*tmp_2pik_tilde_/2);
fnorm_disp(flag_verbose,'tmp_f_k_fft_',tmp_f_k_fft_,'tmp_f_k_1d3_',tmp_f_k_1d3_);
fnorm_disp(flag_verbose,'tmp_f_k_1d3_',tmp_f_k_1d3_,'tmp_f_k_1d2_',tmp_f_k_1d2_);
end;%for tmp_isgn = [-1,+1];

%%%%%%%%;
% With this in mind we check k_p_r_max. ;
%%%%%%%%;
flag_check_k = (k_p_r_max< (n_x_u-1)/(4*x_p_r_max)) ;
if (flag_verbose>0); disp(sprintf(' %% flag_check_k %d: k_p_r_max %0.6f (n_x_u-1)/(4*x_p_r_max) %0.6f',flag_check_k,k_p_r_max,(n_x_u-1)/(4*x_p_r_max))); end;
if ~flag_check_k; disp(sprintf(' %% Warning: flag_check_k %d: k_p_r_max %0.6f >= (n_x_u-1)/(4*x_p_r_max) %0.6f',flag_check_k,k_p_r_max,(n_x_u-1)/(4*x_p_r_max))); end;

%%%%%%%%;
% Now we create j_k_u_fft3d0___. ;
%%%%%%%%;
n_x_u_pad = 3*n_x_u;
tmp_2pik_fft1d0_ = transpose(periodize([0:n_x_u_pad-1],-n_x_u_pad/2,+n_x_u_pad/2))*(2*pi)/(2*x_p_r_max)*(n_x_u-1)/(n_x_u_pad);
tmp_2pik_fft1d0_lim_ = [min(tmp_2pik_fft1d0_),max(tmp_2pik_fft1d0_)];
[tmp_2pik_0_fft3d0___,tmp_2pik_1_fft3d0___,tmp_2pik_2_fft3d0___] = ndgrid(tmp_2pik_fft1d0_,tmp_2pik_fft1d0_,tmp_2pik_fft1d0_);
tmp_isgn = -1;
if tmp_isgn == -1; j_k_u_fft3d0___ =       1.0  *fftshift( fftn(reshape(j_x_u_form_,[n_x_u,n_x_u,n_x_u]),[n_x_u_pad,n_x_u_pad,n_x_u_pad]).*exp(-i*tmp_isgn*(tmp_2pik_0_fft3d0___+tmp_2pik_1_fft3d0___+tmp_2pik_2_fft3d0___)*x_p_r_max)); end;
if tmp_isgn == +1; j_k_u_fft3d0___ = n_x_u_pad^3*fftshift(ifftn(reshape(j_x_u_form_,[n_x_u,n_x_u,n_x_u]),[n_x_u_pad,n_x_u_pad,n_x_u_pad]).*exp(-i*tmp_isgn*(tmp_2pik_0_fft3d0___+tmp_2pik_1_fft3d0___+tmp_2_2pik_fft3d0___)*x_p_r_max)); end;
tmp_2pik_0_fft3d3___ = tmp_2pik_0_fft3d0___;
tmp_2pik_1_fft3d3___ = tmp_2pik_1_fft3d0___;
tmp_2pik_2_fft3d3___ = tmp_2pik_2_fft3d0___;
j_k_u_fft3d3___ = fftshift(reshape(finufft3d3(x_u_0___(:),x_u_1___(:),x_u_2___(:),j_x_u_form_(:),tmp_isgn,1e-9,tmp_2pik_0_fft3d3___(:),tmp_2pik_1_fft3d3___(:),tmp_2pik_2_fft3d3___(:)),[n_x_u_pad,n_x_u_pad,n_x_u_pad]));
fnorm_disp(flag_verbose,'j_k_u_fft3d0___',j_k_u_fft3d0___,'j_k_u_fft3d3___',j_k_u_fft3d3___);

%%%%%%%%;
% Now we create the interpolation operator. ;
%%%%%%%%;
GB_max = 32;
n_order = 3;
n_S_per_Sbatch = 128;
GB_per_Sbatch = (1+n_order)^4 * n_w_max*n_k_p_r*n_S_per_Sbatch*16/1e9; %<-- n_nodex^4*n_scatter used for flag_dd. ;
if (flag_verbose>0); disp(sprintf(' %% GB_per_Sbatch %0.2f',GB_per_Sbatch)); end;
if (GB_per_Sbatch>=GB_max); disp(sprintf(' %% Warning, GB_per_Sbatch %d',GB_per_Sbatch)); end;
scatter_from_tensor_s012__ = sparse([],[],[],n_w_max*n_k_p_r*0,n_x_u_pad^3,0);
n_Sbatch = ceil(n_S/n_S_per_Sbatch);
if (flag_verbose>0); disp(sprintf(' %% n_Sbatch %d',n_Sbatch)); end;
for nSbatch=0:n_Sbatch-1;
index_S_in_Sbatch_ = nSbatch*n_S_per_Sbatch + (0:n_S_per_Sbatch-1);
index_S_in_Sbatch_ = index_S_in_Sbatch_(find(index_S_in_Sbatch_<n_S)); n_S_sub = numel(index_S_in_Sbatch_);
if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (flag_verbose>0 & mod(nSbatch,32)==0); disp(sprintf(' %% nSbatch %d/%d index_S_in_Sbatch_ %d-->%d',nSbatch,n_Sbatch,index_S_in_Sbatch_(1+0),index_S_in_Sbatch_(1+n_S_sub-1))); end;
if (n_S_sub>0);
tmp_t = tic();
k_c_0_wkS_sub___ = k_c_0_wkS___(:,:,1+index_S_in_Sbatch_);
k_c_1_wkS_sub___ = k_c_1_wkS___(:,:,1+index_S_in_Sbatch_);
k_c_2_wkS_sub___ = k_c_2_wkS___(:,:,1+index_S_in_Sbatch_);
[ ...
 scatter_from_tensor_sub_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_x_u_pad ...
,n_x_u_pad ...
,n_x_u_pad ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,n_w_max*n_k_p_r*n_S_sub ...
,2*pi*k_c_0_wkS_sub___(:) ...
,2*pi*k_c_1_wkS_sub___(:) ...
,2*pi*k_c_2_wkS_sub___(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% nSbatch %d/%d: scatter_from_tensor_sub_s012__: %0.6fs',nSbatch,n_Sbatch,tmp_t)); end;
scatter_from_tensor_s012__ = cat(1,scatter_from_tensor_s012__,scatter_from_tensor_sub_s012__);
end;%if (n_S_sub>0);
end;%for nSbatch=0:n_Sbatch-1;
%%%%%%%%;
% Now use the interpolation operator to build the templates. ;
%%%%%%%%;
tmp_t = tic;
S_3_k_p_wkS__ = reshape(scatter_from_tensor_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% scatter_from_tensor_s012__: S_3_k_p_wkS__ time %0.2fs',tmp_t)); end;
%%%%;
S_4_k_p_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_4_k_p_wkS__(:,1+nS) = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_3_k_p_wkS__(:,1+nS),-viewing_gamma_z_S_(1+nS));
end;%for nS=0:n_S-1;
%%%%;
fnorm_disp(flag_verbose,'S_0_k_p_wkS__',S_0_k_p_wkS__,'S_4_k_p_wkS__',S_4_k_p_wkS__);

%%%%%%%%;
% Compare results to fft_template_4. ;
%%%%%%%%;
parameter_fft = struct('type','parameter');
parameter_fft.flag_verbose=1;
parameter_fft.GB_max = 32;
parameter_fft.n_order = 3;
parameter_fft.n_S_per_Sbatch = 128;
parameter_fft.n_x_u_pad_factor = 3.0;
[ ...
 parameter ...
,S_5_k_p_wkS__ ...
,b_k_u_fft3d0___ ...
] = ...
fft_template_4( ...
 parameter_fft ...
,n_x_u ...
,x_p_r_max ...
,j_x_u_form_ ...
,k_p_r_max ...
,k_eq_d ...
,viewing_k_eq_d ...
,template_k_eq_d ...
,n_w_max ...
);
fnorm_disp(flag_verbose,'S_4_k_p_wkS__',S_4_k_p_wkS__,'S_5_k_p_wkS__',S_5_k_p_wkS__);

%%%%%%%%;
% Now test out differentiation for a small subset of templates. ;
%%%%%%%%;
n_3 = 3;
n_S_sub = 5;
viewing_polar_a_S_sub_ = viewing_polar_a_all_(1:n_S_sub);
viewing_azimu_b_S_sub_ = viewing_azimu_b_all_(1:n_S_sub);
viewing_gamma_z_S_sub_ = 2*pi*rand(n_S_sub,1);
dtau = 1e-5;
dviewing_polar_a_S_sub_ = randn(n_S_sub,1);
dviewing_azimu_b_S_sub_ = randn(n_S_sub,1);
dviewing_gamma_z_S_sub_ = randn(n_S_sub,1);
dtau_S_sub3__ = [ dviewing_polar_a_S_sub_ , dviewing_azimu_b_S_sub_ , dviewing_gamma_z_S_sub_ ];
dtau_S_sub33___ = bsxfun(@times,reshape(dtau_S_sub3__,[n_S_sub,3,1]),reshape(dtau_S_sub3__,[n_S_sub,1,3]));
%%%%%%%%;
% Note that the convention in cg_rhs_2 is to ;
% subtract viewing_gamma_z from inplane_gamma_z. ;
%%%%%%%%;
[ ...
 k_p_polar_a_wS_sub__ ...
,k_p_azimu_b_wS_sub__ ...
,k_c_0_wS_sub__ ...
,k_c_1_wS_sub__ ...
,k_c_2_wS_sub__ ...
,k_p_r01_wS_sub__ ...
,dtau_k_p_polar_a_wS_sub3___ ...
,dtau_k_p_azimu_b_wS_sub3___ ...
,dtau_k_c_0_wS_sub3___ ...
,dtau_k_c_1_wS_sub3___ ...
,dtau_k_c_2_wS_sub3___ ...
,dtau_k_p_r01_wS_sub3___ ...
,dtau_dtau_k_p_polar_a_wS_sub33____ ...
,dtau_dtau_k_p_azimu_b_wS_sub33____ ...
,dtau_dtau_k_c_0_wS_sub33____ ...
,dtau_dtau_k_c_1_wS_sub33____ ...
,dtau_dtau_k_c_2_wS_sub33____ ...
,dtau_dtau_k_p_r01_wS_sub33____ ...
] = ...
cg_rhs_2( ...
 n_S_sub ...
,n_w_max ...
,viewing_polar_a_S_sub_ ...
,viewing_azimu_b_S_sub_ ...
,viewing_gamma_z_S_sub_ ...
);
%%%%;
[ ...
 ~ ...
,~ ...
,k_c_0_wS_sub_pos__ ...
,k_c_1_wS_sub_pos__ ...
,k_c_2_wS_sub_pos__ ...
] = ...
cg_rhs_2( ...
 n_S_sub ...
,n_w_max ...
,viewing_polar_a_S_sub_ + 1.0*dtau*dviewing_polar_a_S_sub_ ...
,viewing_azimu_b_S_sub_ + 1.0*dtau*dviewing_azimu_b_S_sub_ ...
,viewing_gamma_z_S_sub_ + 1.0*dtau*dviewing_gamma_z_S_sub_ ...
);
[ ...
 ~ ...
,~ ...
,k_c_0_wS_sub_neg__ ...
,k_c_1_wS_sub_neg__ ...
,k_c_2_wS_sub_neg__ ...
] = ...
cg_rhs_2( ...
 n_S_sub ...
,n_w_max ...
,viewing_polar_a_S_sub_ - 1.0*dtau*dviewing_polar_a_S_sub_ ...
,viewing_azimu_b_S_sub_ - 1.0*dtau*dviewing_azimu_b_S_sub_ ...
,viewing_gamma_z_S_sub_ - 1.0*dtau*dviewing_gamma_z_S_sub_ ...
);
%%%%%%%%;
k_p_polar_a_wkS_sub___ = bsxfun(@times,ones(1,n_k_p_r),reshape(k_p_polar_a_wS_sub__,[n_w_max,1,n_S_sub]));
k_p_azimu_b_wkS_sub___ = bsxfun(@times,ones(1,n_k_p_r),reshape(k_p_azimu_b_wS_sub__,[n_w_max,1,n_S_sub]));
k_c_0_wkS_sub___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_0_wS_sub__,[n_w_max,1,n_S_sub]));
k_c_1_wkS_sub___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_1_wS_sub__,[n_w_max,1,n_S_sub]));
k_c_2_wkS_sub___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_2_wS_sub__,[n_w_max,1,n_S_sub]));
dtau_k_c_0_wkS_sub3____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1]),reshape(dtau_k_c_0_wS_sub3___,[n_w_max,1,n_S_sub,n_3]));
dtau_k_c_1_wkS_sub3____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1]),reshape(dtau_k_c_1_wS_sub3___,[n_w_max,1,n_S_sub,n_3]));
dtau_k_c_2_wkS_sub3____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1]),reshape(dtau_k_c_2_wS_sub3___,[n_w_max,1,n_S_sub,n_3]));
dtau_dtau_k_c_0_wkS_sub33_____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1,1]),reshape(dtau_dtau_k_c_0_wS_sub33____,[n_w_max,1,n_S_sub,n_3,n_3]));
dtau_dtau_k_c_1_wkS_sub33_____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1,1]),reshape(dtau_dtau_k_c_1_wS_sub33____,[n_w_max,1,n_S_sub,n_3,n_3]));
dtau_dtau_k_c_2_wkS_sub33_____ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1,1,1]),reshape(dtau_dtau_k_c_2_wS_sub33____,[n_w_max,1,n_S_sub,n_3,n_3]));
%%%%;
k_c_0_wkS_sub_pos___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_0_wS_sub_pos__,[n_w_max,1,n_S_sub]));
k_c_1_wkS_sub_pos___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_1_wS_sub_pos__,[n_w_max,1,n_S_sub]));
k_c_2_wkS_sub_pos___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_2_wS_sub_pos__,[n_w_max,1,n_S_sub]));
k_c_0_wkS_sub_neg___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_0_wS_sub_neg__,[n_w_max,1,n_S_sub]));
k_c_1_wkS_sub_neg___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_1_wS_sub_neg__,[n_w_max,1,n_S_sub]));
k_c_2_wkS_sub_neg___ = bsxfun(@times,reshape(k_p_r_,[1,n_k_p_r,1]),reshape(k_c_2_wS_sub_neg__,[n_w_max,1,n_S_sub]));
%%%%;
tmp_t = tic();
[ ...
 scatter_from_tensor_sub_s012__ ...
 dscatterd0_from_tensor_sub_s012__ ...
 dscatterd1_from_tensor_sub_s012__ ...
 dscatterd2_from_tensor_sub_s012__ ...
 ddscatterd00_from_tensor_sub_s012__ ...
 ddscatterd01_from_tensor_sub_s012__ ...
 ddscatterd02_from_tensor_sub_s012__ ...
 ddscatterd11_from_tensor_sub_s012__ ...
 ddscatterd12_from_tensor_sub_s012__ ...
 ddscatterd22_from_tensor_sub_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_x_u_pad ...
,n_x_u_pad ...
,n_x_u_pad ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,n_w_max*n_k_p_r*n_S_sub ...
,2*pi*k_c_0_wkS_sub___(:) ...
,2*pi*k_c_1_wkS_sub___(:) ...
,2*pi*k_c_2_wkS_sub___(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% scatter_from_tensor_sub_s012__: %0.6fs',tmp_t)); end;
tmp_t = tic();
[ ...
 scatter_from_tensor_sub_pos_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_x_u_pad ...
,n_x_u_pad ...
,n_x_u_pad ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,n_w_max*n_k_p_r*n_S_sub ...
,2*pi*k_c_0_wkS_sub_pos___(:) ...
,2*pi*k_c_1_wkS_sub_pos___(:) ...
,2*pi*k_c_2_wkS_sub_pos___(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% scatter_from_tensor_sub_pos_s012__: %0.6fs',tmp_t)); end;
tmp_t = tic();
[ ...
 scatter_from_tensor_sub_neg_s012__ ...
] = ...
volume_k_c_scatter_from_tensor_interpolate_n_8( ...
 n_order ...
,n_x_u_pad ...
,n_x_u_pad ...
,n_x_u_pad ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,tmp_2pik_fft1d0_lim_ ...
,n_w_max*n_k_p_r*n_S_sub ...
,2*pi*k_c_0_wkS_sub_neg___(:) ...
,2*pi*k_c_1_wkS_sub_neg___(:) ...
,2*pi*k_c_2_wkS_sub_neg___(:) ...
);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% scatter_from_tensor_sub_neg_s012__: %0.6fs',tmp_t)); end;
%%%%;
   S_sub_mid_k_p_wkS_sub___ = reshape(   scatter_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
dSd0_sub_mid_k_p_wkS_sub___ = reshape(dscatterd0_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
dSd1_sub_mid_k_p_wkS_sub___ = reshape(dscatterd1_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
dSd2_sub_mid_k_p_wkS_sub___ = reshape(dscatterd2_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd00_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd00_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd01_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd01_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd02_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd02_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd11_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd11_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd12_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd12_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
ddSd22_sub_mid_k_p_wkS_sub___ = reshape(ddscatterd22_from_tensor_sub_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
dtau_S_sub_mid_k_p_wkS_sub___ = ...
+ (2*pi)*dSd0_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_k_c_0_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4) ...
+ (2*pi)*dSd1_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_k_c_1_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4) ...
+ (2*pi)*dSd2_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_k_c_2_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4) ...
;
dtau_k_c_0_wkS_sub___ = sum(bsxfun(@times,dtau_k_c_0_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4);
dtau_k_c_1_wkS_sub___ = sum(bsxfun(@times,dtau_k_c_1_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4);
dtau_k_c_2_wkS_sub___ = sum(bsxfun(@times,dtau_k_c_2_wkS_sub3____,reshape(dtau_S_sub3__,[1,1,n_S_sub,n_3])),4);
dtau_dtau_S_sub_mid_k_p_wkS_sub___ = ...
+ (2*pi)*dSd0_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_dtau_k_c_0_wkS_sub33_____,reshape(dtau_S_sub33___,[1,1,n_S_sub,n_3,n_3])),[4,5]) ...
+ (2*pi)*dSd1_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_dtau_k_c_1_wkS_sub33_____,reshape(dtau_S_sub33___,[1,1,n_S_sub,n_3,n_3])),[4,5]) ...
+ (2*pi)*dSd2_sub_mid_k_p_wkS_sub___ .* sum(bsxfun(@times,dtau_dtau_k_c_2_wkS_sub33_____,reshape(dtau_S_sub33___,[1,1,n_S_sub,n_3,n_3])),[4,5]) ...
+ 1.0*(2*pi)*(2*pi)*ddSd00_sub_mid_k_p_wkS_sub___ .* dtau_k_c_0_wkS_sub___ .* dtau_k_c_0_wkS_sub___ ...
+ 2.0*(2*pi)*(2*pi)*ddSd01_sub_mid_k_p_wkS_sub___ .* dtau_k_c_0_wkS_sub___ .* dtau_k_c_1_wkS_sub___ ...
+ 2.0*(2*pi)*(2*pi)*ddSd02_sub_mid_k_p_wkS_sub___ .* dtau_k_c_0_wkS_sub___ .* dtau_k_c_2_wkS_sub___ ...
+ 1.0*(2*pi)*(2*pi)*ddSd11_sub_mid_k_p_wkS_sub___ .* dtau_k_c_1_wkS_sub___ .* dtau_k_c_1_wkS_sub___ ...
+ 2.0*(2*pi)*(2*pi)*ddSd12_sub_mid_k_p_wkS_sub___ .* dtau_k_c_1_wkS_sub___ .* dtau_k_c_2_wkS_sub___ ...
+ 1.0*(2*pi)*(2*pi)*ddSd22_sub_mid_k_p_wkS_sub___ .* dtau_k_c_2_wkS_sub___ .* dtau_k_c_2_wkS_sub___ ...
;
%%%%;
S_sub_pos_k_p_wkS_sub___ = reshape(   scatter_from_tensor_sub_pos_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
S_sub_neg_k_p_wkS_sub___ = reshape(   scatter_from_tensor_sub_neg_s012__*j_k_u_fft3d0___(:).*weight_xxx_u_(:),[n_w_max,n_k_p_r,n_S_sub]);
dtau_S_sub_dif_k_p_wkS_sub___ = (S_sub_pos_k_p_wkS_sub___ - S_sub_neg_k_p_wkS_sub___)/max(1e-12,2*dtau);
dtau_dtau_S_sub_dif_k_p_wkS_sub___ = (S_sub_pos_k_p_wkS_sub___ - 2*S_sub_mid_k_p_wkS_sub___ + S_sub_neg_k_p_wkS_sub___)/max(1e-12,dtau^2);
%%%%;
fnorm_disp(flag_verbose,'dtau_S_sub_dif_k_p_wkS_sub___',dtau_S_sub_dif_k_p_wkS_sub___,'dtau_S_sub_mid_k_p_wkS_sub___',dtau_S_sub_mid_k_p_wkS_sub___);
tmp_index_ = efind(abs(dtau_dtau_S_sub_dif_k_p_wkS_sub___)<max(abs(dtau_dtau_S_sub_mid_k_p_wkS_sub___(:))));
if (flag_verbose>0); disp(sprintf(' %% Note: some of dtau_dtau_S_sub_dif_k_p_wkS_sub___ are inaccurate: using numel(tmp_index_) %d/%d',numel(tmp_index_),(n_w_sum*n_S_sub))); end;
fnorm_disp(flag_verbose,'dtau_dtau_S_sub_dif_k_p_wkS_sub___(1+tmp_index_)',dtau_dtau_S_sub_dif_k_p_wkS_sub___(1+tmp_index_),'dtau_dtau_S_sub_mid_k_p_wkS_sub___(1+tmp_index_)',dtau_dtau_S_sub_mid_k_p_wkS_sub___(1+tmp_index_));
flag_check=0;
if flag_check;
figure(1+nf);nf=nf+1;clf;figsml;
plot(sort(log(abs(dtau_dtau_S_sub_dif_k_p_wkS_sub___(:)-dtau_dtau_S_sub_mid_k_p_wkS_sub___(:)))),'.');
xlabel('index','Interpreter','none');
ylabel('log(abs(diff))','Interpreter','none');
title('sort(abs(dtau_dtau_S_sub_dif_k_p_wkS_sub___(:)-dtau_dtau_S_sub_mid_k_p_wkS_sub___(:)))','Interpreter','none');
end;%if flag_check;

%%%%%%%%;


if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;




