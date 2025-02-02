%%%%%%%%;
% Test integration in fourier and physical space. ;
% Specifically looks at functions which are zero on fourier-polar grid. ;
%%%%%%%%;

platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;

dir_jpg = sprintf('/%s/rangan/dir_cryoem/dir_CryoBIFE_MD/dir_CryoLike_fig',string_root);

flag_verbose = 1;
flag_plot = 1; nf=0;
n_x_c = 256; %<-- this is the high-resolution spatial-grid. ;
k_hi_int = ceil(1*n_x_c/1); %<-- this is the hi-frequency fourier-grid. ;
k_lo_int = ceil(1*n_x_c/2); %<-- this is the lo-frequency fourier-grid. ;
k_eq_d_double = 1.0; %<-- 0.5 means increased radial-resolution in fourier-space. ;
n_w_int = 1; %<-- 2 means increased angular-resolution in fourier-space. ;

%%%%%%%%;
% First define integral of <f,f>. ;
%%%%%%%%;
h2d_ = @(kd) 4*pi^2*(besselj(0,kd) + besselj(2,kd)); % calculates <f_j,f_k>, normalized so that <f,f> = (4*pi^2);
dh2d_ = @(kd) 4*pi^3*(besselj(-1,kd) - besselj(+3,kd));

%%%%%%%%;
% Define spatial grid. ;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = half_diameter_x_c;
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c); dx_0 = mean(diff(x_c_0_));
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c); dx_1 = mean(diff(x_c_1_));
[x_c_0__,x_c_1__] = ndgrid(x_c_0_,x_c_1_); n_xx_c = n_x_c^2;
dx = dx_0; weight_xx_c_ = dx^2;

%%%%%%%%;
% Now set up quadrature for k_hi. ;
%%%%%%%%;
k_hi_p_r_max = k_hi_int/(2*pi); k_hi_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_hi_p_r ...
,k_hi_p_r_ ...
,weight_3d_k_hi_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_hi_p_r_max ...
,k_hi_eq_d ...
,str_L ...
);
n_w_hi_max = n_w_int*2*(k_hi_int+1);
n_w_hi_0in_ = n_w_hi_max*ones(n_k_hi_p_r,1);
[ ...
 n_w_hi_ ...
,weight_2d_k_hi_p_r_ ...
,weight_2d_hi_wk_ ...
,k_hi_p_r_wk_ ...
,k_hi_p_w_wk_ ...
,k_hi_c_0_wk_ ...
,k_hi_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_hi_p_r ...
,k_hi_p_r_ ...
,k_hi_p_r_max ...
,-1 ...
,n_w_hi_0in_ ...
);
n_w_hi_sum = sum(n_w_hi_);
n_w_hi_csum_ = cumsum([0;n_w_hi_]);
%%%%%%%%;
% Now set up quadrature for k_lo. ;
%%%%%%%%;
k_lo_p_r_max = k_lo_int/(2*pi); k_lo_eq_d = k_eq_d_double/(2*pi); str_L = 'L';
[ ...
 n_k_lo_p_r ...
,k_lo_p_r_ ...
,weight_3d_k_lo_p_r_ ...
] = ...
get_weight_3d_1( ...
 flag_verbose ...
,k_lo_p_r_max ...
,k_lo_eq_d ...
,str_L ...
);
n_w_lo_max = n_w_int*2*(k_lo_int+1);
n_w_lo_0in_ = n_w_lo_max*ones(n_k_lo_p_r,1);
[ ...
 n_w_lo_ ...
,weight_2d_k_lo_p_r_ ...
,weight_2d_lo_wk_ ...
,k_lo_p_r_wk_ ...
,k_lo_p_w_wk_ ...
,k_lo_c_0_wk_ ...
,k_lo_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_lo_p_r ...
,k_lo_p_r_ ...
,k_lo_p_r_max ...
,-1 ...
,n_w_lo_0in_ ...
);
n_w_lo_sum = sum(n_w_lo_);
n_w_lo_csum_ = cumsum([0;n_w_lo_]);
%%%%%%%%;
% And set default k to k_lo. ;
%%%%%%%%;
k_int = k_lo_int;
k_p_r_max = k_lo_p_r_max;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now define and integrate a gaussian. ;
% The width of the gaussian is specifically chosen so that it is well-integrated in both x_c and k_p. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Transforming centered gaussian:')); end;
sigma_k = k_p_r_max/8.0; sigma_x = 1.0/max(1e-12,sigma_k)/(2*pi); %<-- fix sigma_k. ;
%sigma_x = 1/16; sigma_k = 1.0/max(1e-12,sigma_x)/(2*pi); %<-- fix sigma_x. ;
nrm_x_from_k = sigma_x^2 * sqrt((2*pi)^2);
g_k_hi_p_form_ = 1/sqrt(2*pi)^2 .* 1/max(1e-12,sigma_k^2) .* exp(-(k_hi_c_0_wk_.^2 + k_hi_c_1_wk_.^2)/max(1e-12,2*sigma_k^2)) ;
g_k_hi_p_form_l1 = sum(g_k_hi_p_form_.*weight_2d_hi_wk_*(2*pi)^2,'all');
g_k_p_full_l1 = 1.0;
g_k_hi_p_form_l2 = sqrt(sum(abs(g_k_hi_p_form_).^2.*weight_2d_hi_wk_*(2*pi)^2,'all'));
g_k_p_full_l2 = sqrt( (sqrt(2*pi) * sigma_k * sqrt(2)).^(-2) );
g_x_c_form_ = nrm_x_from_k * 1/sqrt(2*pi)^2 .* 1/max(1e-12,sigma_x^2) .* exp(-(x_c_0__.^2 + x_c_1__.^2)/max(1e-12,2*sigma_x^2)) ;
g_x_c_form_l1 = sum(g_x_c_form_.*weight_xx_c_,'all');
g_x_c_full_l1 = nrm_x_from_k*1.0;
g_x_c_form_l2 = sqrt(sum(abs(g_x_c_form_).^2.*weight_xx_c_,'all'));
g_x_c_full_l2 = nrm_x_from_k * sqrt( (sqrt(2*pi) * sigma_x * sqrt(2)).^(-2) );
fnorm_disp(flag_verbose,'g_k_p_full_l1',g_k_p_full_l1,'g_k_hi_p_form_l1',g_k_hi_p_form_l1);
fnorm_disp(flag_verbose,'g_k_p_full_l2',g_k_p_full_l2,'g_k_hi_p_form_l2',g_k_hi_p_form_l2);
fnorm_disp(flag_verbose,'g_x_c_full_l1',g_x_c_full_l1,'g_x_c_form_l1',g_x_c_form_l1);
fnorm_disp(flag_verbose,'g_x_c_full_l2',g_x_c_full_l2,'g_x_c_form_l2',g_x_c_form_l2);
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
g_k_hi_p_quad_ = xxnufft2d3(n_xx_c,x_c_0__(:)*eta,x_c_1__(:)*eta,g_x_c_form_(:).*weight_xx_c_(:),-1,1e-12,n_w_hi_sum,2*pi*k_hi_c_0_wk_/eta,2*pi*k_hi_c_1_wk_/eta)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: g_k_hi_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_hi_p_r_max; tmp_t = tic;
g_x_hi_c_quad_ = xxnufft2d3(n_w_hi_sum,2*pi*k_hi_c_0_wk_*eta,2*pi*k_hi_c_1_wk_*eta,g_k_hi_p_form_.*(2*pi)^4.*weight_2d_hi_wk_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: g_x_hi_c_quad_ time %0.2fs',tmp_t)); end;
g_x_hi_c_quad_ = reshape(g_x_hi_c_quad_,[n_x_c,n_x_c]);
%%%%;
g_k_hi_p_quad_l1 = sum(g_k_hi_p_quad_.*weight_2d_hi_wk_*(2*pi)^2,'all');
g_k_hi_p_quad_l2 = sqrt(sum(abs(g_k_hi_p_quad_).^2.*weight_2d_hi_wk_*(2*pi)^2,'all'));
g_x_hi_c_quad_l1 = sum(g_x_hi_c_quad_.*weight_xx_c_,'all');
g_x_hi_c_quad_l2 = sqrt(sum(abs(g_x_hi_c_quad_).^2.*weight_xx_c_,'all'));
%%%%%%%%;
fnorm_disp(flag_verbose,'g_k_p_full_l1',g_k_p_full_l1,'g_k_hi_p_quad_l1',g_k_hi_p_quad_l1);
fnorm_disp(flag_verbose,'g_k_hi_p_form_l1',g_k_hi_p_form_l1,'g_k_hi_p_quad_l1',g_k_hi_p_quad_l1);
fnorm_disp(flag_verbose,'g_k_p_full_l2',g_k_p_full_l2,'g_k_hi_p_quad_l2',g_k_hi_p_quad_l2);
fnorm_disp(flag_verbose,'g_k_hi_p_form_l2',g_k_hi_p_form_l2,'g_k_hi_p_quad_l2',g_k_hi_p_quad_l2);
fnorm_disp(flag_verbose,'g_x_c_full_l1',g_x_c_full_l1,'g_x_hi_c_quad_l1',g_x_hi_c_quad_l1);
fnorm_disp(flag_verbose,'g_x_c_form_l1',g_x_c_form_l1,'g_x_hi_c_quad_l1',g_x_hi_c_quad_l1);
fnorm_disp(flag_verbose,'g_x_c_full_l2',g_x_c_full_l2,'g_x_hi_c_quad_l2',g_x_hi_c_quad_l2);
fnorm_disp(flag_verbose,'g_x_c_form_l2',g_x_c_form_l2,'g_x_hi_c_quad_l2',g_x_hi_c_quad_l2);

%%%%%%%%;
% Now we can translate the gaussian in x_c, ;
% but we make sure to keep the support within the box. ;
% In other words, we make sure we can still integrate it to about 4-6 digits of accuracy. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Transforming shifted gaussian:')); end;
delta_x_c_max = half_diameter_x_c/4.0;
delta_x_c_ = randn(2,1); delta_x_c_ = delta_x_c_max*delta_x_c_/max(1e-12,fnorm(delta_x_c_));
h_k_hi_p_form_ = g_k_hi_p_form_ .* exp(-i*2*pi*( k_hi_c_0_wk_.*delta_x_c_(1+0) + k_hi_c_1_wk_.*delta_x_c_(1+1) ));
h_k_hi_p_form_l1 = sum(h_k_hi_p_form_.*weight_2d_hi_wk_*(2*pi)^2,'all');
h_k_p_full_l1 = exp(-(delta_x_c_(1+0).^2 + delta_x_c_(1+1).^2)/max(1e-12,2*sigma_x^2)) ;
h_k_hi_p_form_l2 = sqrt(sum(abs(h_k_hi_p_form_).^2.*weight_2d_hi_wk_*(2*pi)^2,'all'));
h_k_p_full_l2 = sqrt( (sqrt(2*pi) * sigma_k * sqrt(2)).^(-2) );
h_x_c_form_ = nrm_x_from_k * 1/sqrt(2*pi)^2 .* 1/max(1e-12,sigma_x^2) .* exp(-((x_c_0__-delta_x_c_(1+0)).^2 + (x_c_1__-delta_x_c_(1+1)).^2)/max(1e-12,2*sigma_x^2)) ;
h_x_c_form_l1 = sum(h_x_c_form_.*weight_xx_c_,'all');
h_x_c_full_l1 = nrm_x_from_k*1.0;
h_x_c_form_l2 = sqrt(sum(abs(h_x_c_form_).^2.*weight_xx_c_,'all'));
h_x_c_full_l2 = nrm_x_from_k * sqrt( (sqrt(2*pi) * sigma_x * sqrt(2)).^(-2) );
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l1 vs h_k_hi_p_form_l1: %0.16f %%<-- only accurate if shifted gaussian decays within k-domain. ',fnorm(h_k_p_full_l1 - h_k_hi_p_form_l1)/max(1e-12,fnorm(h_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l2 vs h_k_hi_p_form_l2: %0.16f %%<-- should be accurate. ',fnorm(h_k_p_full_l2 - h_k_hi_p_form_l2)/max(1e-12,fnorm(h_k_p_full_l2)))); end;
fnorm_disp(flag_verbose,'h_k_p_full_l1',h_k_p_full_l1,'h_k_hi_p_form_l1',h_k_hi_p_form_l1);
fnorm_disp(flag_verbose,'h_k_p_full_l2',h_k_p_full_l2,'h_k_hi_p_form_l2',h_k_hi_p_form_l2);
fnorm_disp(flag_verbose,'h_x_c_full_l1',h_x_c_full_l1,'h_x_c_form_l1',h_x_c_form_l1);
fnorm_disp(flag_verbose,'h_x_c_full_l2',h_x_c_full_l2,'h_x_c_form_l2',h_x_c_form_l2);
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
%%%%%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
h_k_hi_p_quad_ = xxnufft2d3(n_xx_c,x_c_0__(:)*eta,x_c_1__(:)*eta,h_x_c_form_(:).*weight_xx_c_(:),-1,1e-12,n_w_hi_sum,2*pi*k_hi_c_0_wk_/eta,2*pi*k_hi_c_1_wk_/eta)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: h_k_hi_p_quad_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_hi_p_r_max; tmp_t = tic;
h_x_hi_c_quad_ = xxnufft2d3(n_w_hi_sum,2*pi*k_hi_c_0_wk_*eta,2*pi*k_hi_c_1_wk_*eta,h_k_hi_p_form_.*(2*pi)^4.*weight_2d_hi_wk_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: h_x_hi_c_quad_ time %0.2fs',tmp_t)); end;
h_x_hi_c_quad_ = reshape(h_x_hi_c_quad_,[n_x_c,n_x_c]);
%%%%;
h_k_hi_p_quad_l1 = sum(h_k_hi_p_quad_.*weight_2d_hi_wk_*(2*pi)^2,'all');
h_k_hi_p_quad_l2 = sqrt(sum(abs(h_k_hi_p_quad_).^2.*weight_2d_hi_wk_*(2*pi)^2,'all'));
h_x_hi_c_quad_l1 = sum(h_x_hi_c_quad_.*weight_xx_c_,'all');
h_x_hi_c_quad_l2 = sqrt(sum(abs(h_x_hi_c_quad_).^2.*weight_xx_c_,'all'));
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l1 vs h_k_hi_p_quad_l1: %0.16f %%<-- only accurate if shifted gaussian decays within k-domain. ',fnorm(h_k_p_full_l1 - h_k_hi_p_quad_l1)/max(1e-12,fnorm(h_k_p_full_l1)))); end;
if (flag_verbose>0); disp(sprintf(' %% %% h_k_p_full_l2 vs h_k_hi_p_quad_l2: %0.16f %%<-- should be accurate. ',fnorm(h_k_p_full_l2 - h_k_hi_p_quad_l2)/max(1e-12,fnorm(h_k_p_full_l2)))); end;
fnorm_disp(flag_verbose,'h_k_p_full_l1',h_k_p_full_l1,'h_k_hi_p_quad_l1',h_k_hi_p_quad_l1);
fnorm_disp(flag_verbose,'h_k_p_full_l2',h_k_p_full_l2,'h_k_hi_p_quad_l2',h_k_hi_p_quad_l2);
fnorm_disp(flag_verbose,'h_x_c_full_l1',h_x_c_full_l1,'h_x_hi_c_quad_l1',h_x_hi_c_quad_l1);
fnorm_disp(flag_verbose,'h_x_c_full_l2',h_x_c_full_l2,'h_x_hi_c_quad_l2',h_x_hi_c_quad_l2);
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;

if flag_plot>1;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
linewidth_use = 2;
markersize_use = 12;
p_row = 2; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
g_k_hi_p_lim_ = max(abs(g_k_hi_p_form_(:)))*[-1,+1];
plot(g_k_hi_p_lim_,g_k_hi_p_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(g_k_hi_p_form_(:)),real(g_k_hi_p_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(g_k_hi_p_lim_); ylim(g_k_hi_p_lim_); axisnotick; axis square;
xlabel('g_k_hi_p_form_','Interpreter','none');
ylabel('g_k_hi_p_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
g_x_c_lim_ = max(abs(g_x_c_form_(:)))*[-1,+1];
plot(g_x_c_lim_,g_x_c_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(g_x_c_form_(:)),real(g_x_hi_c_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(g_x_c_lim_); ylim(g_x_c_lim_); axisnotick; axis square;
xlabel('g_x_c_form_','Interpreter','none');
ylabel('g_x_hi_c_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
h_k_hi_p_lim_ = max(abs(h_k_hi_p_form_(:)))*[-1,+1];
plot(h_k_hi_p_lim_,h_k_hi_p_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(h_k_hi_p_form_(:)),real(h_k_hi_p_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(h_k_hi_p_lim_); ylim(h_k_hi_p_lim_); axisnotick; axis square;
xlabel('h_k_hi_p_form_','Interpreter','none');
ylabel('h_k_hi_p_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1; hold on;
h_x_c_lim_ = max(abs(h_x_c_form_(:)))*[-1,+1];
plot(h_x_c_lim_,h_x_c_lim_,'-','LineWidth',linewidth_use,'Color',0.85*[1,1,1]);
plot(real(h_x_c_form_(:)),real(h_x_hi_c_quad_(:)),'c.','MarkerSize',markersize_use);
xlim(h_x_c_lim_); ylim(h_x_c_lim_); axisnotick; axis square;
xlabel('h_x_c_form_','Interpreter','none');
ylabel('h_x_hi_c_quad_','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
str_sgtitle = sprintf('k_eq_d_double %0.2f k_int %d n_x_c %d',k_eq_d_double,k_int,n_x_c);
sgtitle(str_sgtitle,'Interpreter','none');
drawnow();
end;%if flag_plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now define and integrate a bandlimited plane-wave. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% %% %% %% %% %% %% %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% Transforming bandlimited plane-wave:')); end;
tmp_S_delta_x_c_ = [cos(pi/3);sin(pi/3)]/k_hi_p_r_max;
S_k_hi_p_form_wk_ = exp(+2*pi*i*(k_hi_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_hi_c_1_wk_*tmp_S_delta_x_c_(1+1)));
S_k_lo_p_form_wk_ = exp(+2*pi*i*(k_lo_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_lo_c_1_wk_*tmp_S_delta_x_c_(1+1)));
tmp_kd_hi__ = 2*pi*k_hi_p_r_max*sqrt( (tmp_S_delta_x_c_(1+0) + x_c_0__).^2 + (tmp_S_delta_x_c_(1+1) + x_c_1__).^2 );
S_x_hi_c_form_xx_ = reshape(h2d_(tmp_kd_hi__(:))/(2*pi)^2 * (pi*k_hi_p_r_max^2),[n_x_c,n_x_c]);
tmp_kd_lo__ = 2*pi*k_lo_p_r_max*sqrt( (tmp_S_delta_x_c_(1+0) + x_c_0__).^2 + (tmp_S_delta_x_c_(1+1) + x_c_1__).^2 );
S_x_lo_c_form_xx_ = reshape(h2d_(tmp_kd_lo__(:))/(2*pi)^2 * (pi*k_lo_p_r_max^2),[n_x_c,n_x_c]);
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
S_k_hi_p_quad_wk_ = xxnufft2d3(n_xx_c,x_c_0__(:)*eta,x_c_1__(:)*eta,S_x_hi_c_form_xx_(:).*weight_xx_c_(:),-1,1e-12,n_w_hi_sum,2*pi*k_hi_c_0_wk_/eta,2*pi*k_hi_c_1_wk_/eta)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: S_k_hi_p_quad_wk_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_hi_p_r_max; tmp_t = tic;
S_x_hi_c_quad_xx_ = xxnufft2d3(n_w_hi_sum,2*pi*k_hi_c_0_wk_*eta,2*pi*k_hi_c_1_wk_*eta,S_k_hi_p_form_wk_.*(2*pi)^4.*weight_2d_hi_wk_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: S_x_hi_c_quad_xx_ time %0.2fs',tmp_t)); end;
S_x_hi_c_quad_xx_ = reshape(S_x_hi_c_quad_xx_,[n_x_c,n_x_c]);
%%%%;
R_x_hi_c_quad_xx_ = S_x_hi_c_quad_xx_ - S_x_hi_c_form_xx_;
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
S_k_lo_p_quad_wk_ = xxnufft2d3(n_xx_c,x_c_0__(:)*eta,x_c_1__(:)*eta,S_x_lo_c_form_xx_(:).*weight_xx_c_(:),-1,1e-12,n_w_lo_sum,2*pi*k_lo_c_0_wk_/eta,2*pi*k_lo_c_1_wk_/eta)/sqrt(2*pi)/sqrt(2*pi) * (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: S_k_lo_p_quad_wk_ time %0.2fs',tmp_t)); end;
%%%%;
eta = pi/k_lo_p_r_max; tmp_t = tic;
S_x_lo_c_quad_xx_ = xxnufft2d3(n_w_lo_sum,2*pi*k_lo_c_0_wk_*eta,2*pi*k_lo_c_1_wk_*eta,S_k_lo_p_form_wk_.*(2*pi)^4.*weight_2d_lo_wk_,+1,1e-12,n_xx_c,x_c_0__(:)/eta,x_c_1__(:)/eta)/sqrt(2*pi)/sqrt(2*pi) / (sqrt(2*pi)^2) ;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% %% xxnufft2d3: S_x_lo_c_quad_xx_ time %0.2fs',tmp_t)); end;
S_x_lo_c_quad_xx_ = reshape(S_x_lo_c_quad_xx_,[n_x_c,n_x_c]);
%%%%;
R_x_lo_c_quad_xx_ = S_x_lo_c_quad_xx_ - S_x_lo_c_form_xx_;
%%%%;
if flag_plot>1;
S_x_lim_ = max(abs(S_x_lo_c_quad_xx_),[],'all')*0.75*[-1,+1];
p_row = 2; p_col = 3; np=0;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_hi_c_form_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_hi_c_form_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_hi_c_quad_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_hi_c_quad_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(R_x_hi_c_quad_xx_),S_x_lim_); axis image; axisnotick; title('real(R_x_hi_c_quad_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_hi_c_form_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_hi_c_form_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_lo_c_quad_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_lo_c_quad_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(R_x_lo_c_quad_xx_),S_x_lim_); axis image; axisnotick; title('real(R_x_lo_c_quad_xx_)','Interpreter','none');
end;%if flag_plot;
%%%%%%%%;
R_x_hi_c_quad_l2 = sum(conj(R_x_hi_c_quad_xx_).*R_x_hi_c_quad_xx_,'all')*weight_xx_c_;
S_x_hi_c_quad_l2 = sum(conj(S_x_hi_c_quad_xx_).*S_x_hi_c_quad_xx_,'all')*weight_xx_c_;
S_x_hi_c_form_l2 = sum(conj(S_x_hi_c_form_xx_).*S_x_hi_c_form_xx_,'all')*weight_xx_c_;
if (flag_verbose>0); disp(sprintf(' %% S_x_hi_c_form_l2: %0.6f S_x_hi_c_quad_l2: %0.6f R_x_hi_c_quad_l2: %0.6f ratio: %0.16f',S_x_hi_c_form_l2,S_x_hi_c_quad_l2,R_x_hi_c_quad_l2,R_x_hi_c_quad_l2/max(1e-12,S_x_hi_c_form_l2))); end;
%%%%;
R_x_lo_c_quad_l2 = sum(conj(R_x_lo_c_quad_xx_).*R_x_lo_c_quad_xx_,'all')*weight_xx_c_;
S_x_lo_c_quad_l2 = sum(conj(S_x_lo_c_quad_xx_).*S_x_lo_c_quad_xx_,'all')*weight_xx_c_;
S_x_lo_c_form_l2 = sum(conj(S_x_lo_c_form_xx_).*S_x_lo_c_form_xx_,'all')*weight_xx_c_;
if (flag_verbose>0); disp(sprintf(' %% S_x_lo_c_form_l2: %0.6f S_x_lo_c_quad_l2: %0.6f R_x_lo_c_quad_l2: %0.6f ratio: %0.16f',S_x_lo_c_form_l2,S_x_lo_c_quad_l2,R_x_lo_c_quad_l2,R_x_lo_c_quad_l2/max(1e-12,S_x_lo_c_form_l2))); end;
%%%%%%%%;
% test interpolation-function on symmetric x_c grid. ;
%%%%%%%%;
flag_sym = 1;
S_x_hi_c_qua2_xx_ = ...
interp_k_p_to_x_c_sym_xxnufft( ...
 n_x_c ...
,diameter_x_c ...
,n_x_c ...
,diameter_x_c ...
,flag_sym ...
,n_k_hi_p_r ...
,k_hi_p_r_ ...
,k_hi_p_r_max ...
,n_w_hi_ ...
,weight_2d_hi_wk_ ...
,S_k_hi_p_form_wk_ ...
);
fnorm_disp(flag_verbose,'S_x_hi_c_qua2_xx_',S_x_hi_c_qua2_xx_,'S_x_hi_c_quad_xx_',S_x_hi_c_quad_xx_,' %<-- should be <1e-6');
flag_sym = 1;
S_k_hi_p_qua2_wk_ = ...
interp_x_c_sym_to_k_p_xxnufft( ...
 n_x_c ...
,diameter_x_c ...
,n_x_c ...
,diameter_x_c ...
,flag_sym ...
,S_x_hi_c_form_xx_ ...
,n_k_hi_p_r ...
,k_hi_p_r_ ...
,k_hi_p_r_max ...
,n_w_hi_ ...
,weight_2d_hi_wk_ ...
);
fnorm_disp(flag_verbose,'S_k_hi_p_qua2_wk_',S_k_hi_p_qua2_wk_,'S_k_hi_p_quad_wk_',S_k_hi_p_quad_wk_,' %<-- should be <1e-6');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Generate two bandlimited plane-waves. ;
%%%%%%%%;
tmp_S_delta_x_c_ = [cos(pi/4);sin(pi/4)]/k_hi_p_r_max;
tmp_T_delta_x_c_ = [cos(pi/3);sin(pi/3)]/k_hi_p_r_max;
S_k_hi_p_wk_ = exp(+2*pi*i*(k_hi_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_hi_c_1_wk_*tmp_S_delta_x_c_(1+1)));
T_k_hi_p_wk_ = exp(+2*pi*i*(k_hi_c_0_wk_*tmp_T_delta_x_c_(1+0) + k_hi_c_1_wk_*tmp_T_delta_x_c_(1+1)));
S_k_lo_p_wk_ = exp(+2*pi*i*(k_lo_c_0_wk_*tmp_S_delta_x_c_(1+0) + k_lo_c_1_wk_*tmp_S_delta_x_c_(1+1)));
T_k_lo_p_wk_ = exp(+2*pi*i*(k_lo_c_0_wk_*tmp_T_delta_x_c_(1+0) + k_lo_c_1_wk_*tmp_T_delta_x_c_(1+1)));
S_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,S_k_hi_p_wk_) );
T_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,T_k_hi_p_wk_) );
S_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,S_k_lo_p_wk_) );
T_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,T_k_lo_p_wk_) );
%%%%%%%%;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;
p_row = 2; p_col = 4; np=0;
k_p_axis_ = k_hi_p_r_max*[-1,+1,-1,+1];
S_x_lim_ = max(abs(S_x_from_lo_c_xx_),[],'all')*0.75*[-1,+1];
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(S_k_hi_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(S_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_from_hi_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(T_k_hi_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(T_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_x_from_hi_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(S_k_lo_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(S_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_x_from_lo_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(T_k_lo_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(T_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_x_from_lo_c_xx_)','Interpreter','none');
fname_fig_pre = sprintf('%s/test_CryoLike_physical_vs_fourier_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now test integration of one wave against the other. ;
%%%%%%%%;
SdotT_x_hi = sum(conj(S_x_from_hi_c_xx_).*T_x_from_hi_c_xx_,'all')*weight_xx_c_;
SdotT_k_hi = sum(conj(S_k_hi_p_wk_).*T_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
tmp_kd_hi = 2*pi*k_hi_p_r_max*fnorm(tmp_S_delta_x_c_ - tmp_T_delta_x_c_);
SdotT_form_hi = h2d_(tmp_kd_hi)/(2*pi)^2 * (pi*k_hi_p_r_max^2);
SdotT_x_lo = sum(conj(S_x_from_lo_c_xx_).*T_x_from_lo_c_xx_,'all')*weight_xx_c_;
SdotT_k_lo = sum(conj(S_k_lo_p_wk_).*T_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
tmp_kd_lo = 2*pi*k_lo_p_r_max*fnorm(tmp_S_delta_x_c_ - tmp_T_delta_x_c_);
SdotT_form_lo = h2d_(tmp_kd_lo)/(2*pi)^2 * (pi*k_lo_p_r_max^2);
fnorm_disp(flag_verbose,'SdotT_form_hi',SdotT_form_hi,'SdotT_k_hi',SdotT_k_hi,' %<-- should be <1e-6');
fnorm_disp(flag_verbose,'SdotT_form_hi',SdotT_form_hi,'SdotT_x_hi',SdotT_x_hi,' %<-- should be <1e-2');
fnorm_disp(flag_verbose,'SdotT_form_lo',SdotT_form_lo,'SdotT_k_lo',SdotT_k_lo,' %<-- should be <1e-6');
fnorm_disp(flag_verbose,'SdotT_form_lo',SdotT_form_lo,'SdotT_x_lo',SdotT_x_lo,' %<-- should be <1e-2');
%%%%%%%%;

%%%%%%%%;
% Now test integration of norms. ;
%%%%%%%%;
S_x_from_hi_l2 = sum(conj(S_x_from_hi_c_xx_).*S_x_from_hi_c_xx_,'all')*weight_xx_c_;
S_k_hi_l2 = sum(conj(S_k_hi_p_wk_).*S_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'pi*k_hi_p_r_max^2',pi*k_hi_p_r_max^2,'S_k_hi_l2',S_k_hi_l2,' %<-- should be <1e-6');
fnorm_disp(flag_verbose,'pi*k_hi_p_r_max^2',pi*k_hi_p_r_max^2,'S_x_from_hi_l2',S_x_from_hi_l2,' %<-- should be <1e-2');
S_x_from_lo_l2 = sum(conj(S_x_from_lo_c_xx_).*S_x_from_lo_c_xx_,'all')*weight_xx_c_;
S_k_lo_l2 = sum(conj(S_k_lo_p_wk_).*S_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'pi*k_lo_p_r_max^2',pi*k_lo_p_r_max^2,'S_k_lo_l2',S_k_lo_l2,' %<-- should be <1e-6');
fnorm_disp(flag_verbose,'pi*k_lo_p_r_max^2',pi*k_lo_p_r_max^2,'S_x_from_lo_l2',S_x_from_lo_l2,' %<-- should be <1e-2');
%%%%%%%%;
% test inverse transformation. ;
%%%%%%%%;
S_k_hi_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_x_from_hi_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_) ;
R_k_hi_p_wk_ = S_k_hi_p_wk_ - S_k_hi_reco_p_wk_;
R_k_hi_l2 = sum(conj(R_k_hi_p_wk_).*R_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% S_k_hi_l2 %0.6f, R_k_hi_l2 %0.6f, ratio %0.6f',S_k_hi_l2,R_k_hi_l2,R_k_hi_l2./max(1e-12,S_k_hi_l2))); end;
S_k_lo_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_x_from_lo_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_) ;
R_k_lo_p_wk_ = S_k_lo_p_wk_ - S_k_lo_reco_p_wk_;
R_k_lo_l2 = sum(conj(R_k_lo_p_wk_).*R_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% S_k_lo_l2 %0.6f, R_k_lo_l2 %0.6f, ratio %0.6f',S_k_lo_l2,R_k_lo_l2,R_k_lo_l2./max(1e-12,S_k_lo_l2))); end;
%%%%%%%%;

%%%%%%%%;
% Now restrict to annulus in k-space. ;
%%%%%%%%;
S_annulus_x_c_xx_ = S_x_from_hi_c_xx_ - S_x_from_lo_c_xx_;
T_annulus_x_c_xx_ = T_x_from_hi_c_xx_ - T_x_from_lo_c_xx_;
S_annulus_k_hi_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_annulus_x_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
T_annulus_k_hi_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_annulus_x_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
S_annulus_k_lo_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_annulus_x_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
T_annulus_k_lo_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_annulus_x_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
S_annulus_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,S_annulus_k_hi_p_wk_) );
T_annulus_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,T_annulus_k_hi_p_wk_) );
S_annulus_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,S_annulus_k_lo_p_wk_) );
T_annulus_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,T_annulus_k_lo_p_wk_) );
%%%%%%%%;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;
p_row = 2; p_col = 6; np=0;
k_p_axis_ = k_hi_p_r_max*[-1,+1,-1,+1];
S_x_lim_ = max(abs(S_x_from_lo_c_xx_),[],'all')*0.75*[-1,+1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_annulus_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_annulus_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(S_annulus_k_hi_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(S_annulus_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_annulus_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_annulus_x_from_hi_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_annulus_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_annulus_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(T_annulus_k_hi_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(T_annulus_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_annulus_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_annulus_x_from_hi_c_xx_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_annulus_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_annulus_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(S_annulus_k_lo_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(S_annulus_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_annulus_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_annulus_x_from_lo_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_annulus_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_annulus_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(T_annulus_k_lo_p_wk_),[-1,+1]); axis(k_p_axis_); axis square; axisnotick; title('real(T_annulus_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_annulus_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_annulus_x_from_lo_c_xx_)','Interpreter','none');
%%%%;
fname_fig_pre = sprintf('%s/test_CryoLike_physical_vs_fourier_FIGB',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now test integration of one annulus against the other. ;
%%%%%%%%;
SdotT_annulus_x_hi = sum(conj(S_annulus_x_from_hi_c_xx_).*T_annulus_x_from_hi_c_xx_,'all')*weight_xx_c_;
SdotT_annulus_k_hi = sum(conj(S_annulus_k_hi_p_wk_).*T_annulus_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
tmp_kd_hi = 2*pi*k_hi_p_r_max*fnorm(tmp_S_delta_x_c_ - tmp_T_delta_x_c_);
tmp_kd_lo = 2*pi*k_lo_p_r_max*fnorm(tmp_S_delta_x_c_ - tmp_T_delta_x_c_);
SdotT_annulus_form_hi = h2d_(tmp_kd_hi)/(2*pi)^2 * (pi*k_hi_p_r_max^2) - h2d_(tmp_kd_lo)/(2*pi)^2 * (pi*k_lo_p_r_max^2) ;
SdotT_annulus_x_lo = sum(conj(S_annulus_x_from_lo_c_xx_).*T_annulus_x_from_lo_c_xx_,'all')*weight_xx_c_;
SdotT_annulus_k_lo = sum(conj(S_annulus_k_lo_p_wk_).*T_annulus_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'SdotT_annulus_form_hi',SdotT_annulus_form_hi,'SdotT_annulus_k_hi',SdotT_annulus_k_hi,' %<-- should be <1e-1');
fnorm_disp(flag_verbose,'SdotT_annulus_form_hi',SdotT_annulus_form_hi,'SdotT_annulus_x_hi',SdotT_annulus_x_hi,' %<-- should be <1e-1');
fnorm_disp(flag_verbose,'SdotT_annulus_form_hi',SdotT_annulus_form_hi,'SdotT_annulus_k_lo',SdotT_annulus_k_lo,' %<-- should be ~1.0');
fnorm_disp(flag_verbose,'SdotT_annulus_form_hi',SdotT_annulus_form_hi,'SdotT_annulus_x_lo',SdotT_annulus_x_lo,' %<-- should be ~1.0');
%%%%%%%%;

%%%%%%%%;
% Now set up a noisy pair of test-functions. ;
%%%%%%%%;
rng(0);
%%%%%%%%;
S_noise_x_c_xx_ = randn(n_x_c,n_x_c);
S_noise_k_hi_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_noise_x_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
S_noise_k_hi_l2 = sum(conj(S_noise_k_hi_p_wk_).*S_noise_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
S_noise_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,S_noise_k_hi_p_wk_) );
S_noise_x_from_hi_l2 = sum(conj(S_noise_x_from_hi_c_xx_).*S_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
S_noise_k_hi_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_noise_x_from_hi_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
R_noise_k_hi_p_wk_ = S_noise_k_hi_p_wk_ - S_noise_k_hi_reco_p_wk_;
R_noise_k_hi_l2 = sum(conj(R_noise_k_hi_p_wk_).*R_noise_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% S_noise_k_hi_l2 %0.6f, R_noise_k_hi_l2 %0.6f, ratio %0.6f',S_noise_k_hi_l2,R_noise_k_hi_l2,R_noise_k_hi_l2./max(1e-12,S_noise_k_hi_l2))); end;
S_noise_x_from_hi_reco_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,S_noise_k_hi_reco_p_wk_) );
R_noise_x_from_hi_c_xx_ = S_noise_x_from_hi_c_xx_ - S_noise_x_from_hi_reco_c_xx_;
R_noise_x_from_hi_l2 = sum(conj(R_noise_x_from_hi_c_xx_).*R_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
if (flag_verbose); disp(sprintf(' %% S_noise_x_from_hi_l2 %0.6f, R_noise_x_from_hi_l2 %0.6f, ratio %0.6f',S_noise_x_from_hi_l2,R_noise_x_from_hi_l2,R_noise_x_from_hi_l2./max(1e-12,S_noise_x_from_hi_l2))); end;
%%%%%%%%;
S_noise_k_lo_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_noise_x_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
S_noise_k_lo_l2 = sum(conj(S_noise_k_lo_p_wk_).*S_noise_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
S_noise_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,S_noise_k_lo_p_wk_) );
S_noise_x_from_lo_l2 = sum(conj(S_noise_x_from_lo_c_xx_).*S_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
S_noise_k_lo_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,S_noise_x_from_lo_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
R_noise_k_lo_p_wk_ = S_noise_k_lo_p_wk_ - S_noise_k_lo_reco_p_wk_;
R_noise_k_lo_l2 = sum(conj(R_noise_k_lo_p_wk_).*R_noise_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% S_noise_k_lo_l2 %0.6f, R_noise_k_lo_l2 %0.6f, ratio %0.6f',S_noise_k_lo_l2,R_noise_k_lo_l2,R_noise_k_lo_l2./max(1e-12,S_noise_k_lo_l2))); end;
S_noise_x_from_lo_reco_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,S_noise_k_lo_reco_p_wk_) );
R_noise_x_from_lo_c_xx_ = S_noise_x_from_lo_c_xx_ - S_noise_x_from_lo_reco_c_xx_;
R_noise_x_from_lo_l2 = sum(conj(R_noise_x_from_lo_c_xx_).*R_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
if (flag_verbose); disp(sprintf(' %% S_noise_x_from_lo_l2 %0.6f, R_noise_x_from_lo_l2 %0.6f, ratio %0.6f',S_noise_x_from_lo_l2,R_noise_x_from_lo_l2,R_noise_x_from_lo_l2./max(1e-12,S_noise_x_from_lo_l2))); end;
%%%%%%%%;
rng(1);
T_noise_x_c_xx_ = 1/sqrt(2)*S_noise_x_c_xx_ + 1/sqrt(2)*randn(n_x_c,n_x_c); %<-- close to pi/4 angle with S_noise_x_c_xx_. ;
T_noise_k_hi_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_noise_x_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
T_noise_k_hi_l2 = sum(conj(T_noise_k_hi_p_wk_).*T_noise_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
T_noise_x_from_hi_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,T_noise_k_hi_p_wk_) );
T_noise_x_from_hi_l2 = sum(conj(T_noise_x_from_hi_c_xx_).*T_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
T_noise_k_hi_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_noise_x_from_hi_c_xx_,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_);
R_noise_k_hi_p_wk_ = T_noise_k_hi_p_wk_ - T_noise_k_hi_reco_p_wk_;
R_noise_k_hi_l2 = sum(conj(R_noise_k_hi_p_wk_).*R_noise_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% T_noise_k_hi_l2 %0.6f, R_noise_k_hi_l2 %0.6f, ratio %0.6f',T_noise_k_hi_l2,R_noise_k_hi_l2,R_noise_k_hi_l2./max(1e-12,T_noise_k_hi_l2))); end;
T_noise_x_from_hi_reco_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_hi_p_r,k_hi_p_r_,k_hi_p_r_max,n_w_hi_,weight_2d_hi_wk_,T_noise_k_hi_reco_p_wk_) );
R_noise_x_from_hi_c_xx_ = T_noise_x_from_hi_c_xx_ - T_noise_x_from_hi_reco_c_xx_;
R_noise_x_from_hi_l2 = sum(conj(R_noise_x_from_hi_c_xx_).*R_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
if (flag_verbose); disp(sprintf(' %% T_noise_x_from_hi_l2 %0.6f, R_noise_x_from_hi_l2 %0.6f, ratio %0.6f',T_noise_x_from_hi_l2,R_noise_x_from_hi_l2,R_noise_x_from_hi_l2./max(1e-12,T_noise_x_from_hi_l2))); end;
%%%%%%%%;
T_noise_k_lo_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_noise_x_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
T_noise_k_lo_l2 = sum(conj(T_noise_k_lo_p_wk_).*T_noise_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
T_noise_x_from_lo_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,T_noise_k_lo_p_wk_) );
T_noise_x_from_lo_l2 = sum(conj(T_noise_x_from_lo_c_xx_).*T_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
T_noise_k_lo_reco_p_wk_ = interp_x_c_sym_to_k_p_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,T_noise_x_from_lo_c_xx_,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_);
R_noise_k_lo_p_wk_ = T_noise_k_lo_p_wk_ - T_noise_k_lo_reco_p_wk_;
R_noise_k_lo_l2 = sum(conj(R_noise_k_lo_p_wk_).*R_noise_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
if (flag_verbose); disp(sprintf(' %% T_noise_k_lo_l2 %0.6f, R_noise_k_lo_l2 %0.6f, ratio %0.6f',T_noise_k_lo_l2,R_noise_k_lo_l2,R_noise_k_lo_l2./max(1e-12,T_noise_k_lo_l2))); end;
T_noise_x_from_lo_reco_c_xx_ = real( interp_k_p_to_x_c_sym_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,flag_sym,n_k_lo_p_r,k_lo_p_r_,k_lo_p_r_max,n_w_lo_,weight_2d_lo_wk_,T_noise_k_lo_reco_p_wk_) );
R_noise_x_from_lo_c_xx_ = T_noise_x_from_lo_c_xx_ - T_noise_x_from_lo_reco_c_xx_;
R_noise_x_from_lo_l2 = sum(conj(R_noise_x_from_lo_c_xx_).*R_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
if (flag_verbose); disp(sprintf(' %% T_noise_x_from_lo_l2 %0.6f, R_noise_x_from_lo_l2 %0.6f, ratio %0.6f',T_noise_x_from_lo_l2,R_noise_x_from_lo_l2,R_noise_x_from_lo_l2./max(1e-12,T_noise_x_from_lo_l2))); end;
%%%%%%%%;
if flag_plot;
figure(1+nf);nf=nf+1;figbig;
p_row = 4; p_col = 4; np=0;
k_p_axis_ = k_hi_p_r_max*[-1,+1,-1,+1];
S_x_lim_ = max(abs(S_noise_x_c_xx_),[],'all')*0.75*[-1,+1];
S_k_lim_ = max(abs(S_noise_k_hi_p_wk_),[],'all')*0.75*[-1,+1];
T_x_lim_ = max(abs(T_noise_x_c_xx_),[],'all')*0.75*[-1,+1];
T_k_lim_ = max(abs(T_noise_k_hi_p_wk_),[],'all')*0.75*[-1,+1];
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(S_noise_k_hi_p_wk_),S_k_lim_); axis(k_p_axis_); axis square; axisnotick; title('real(S_noise_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_from_hi_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_from_hi_reco_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_from_hi_reco_c_xx_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_noise_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_hi_p_r,k_hi_p_r_,n_w_hi_,n_w_hi_sum,real(T_noise_k_hi_p_wk_),S_k_lim_); axis(k_p_axis_); axis square; axisnotick; title('real(T_noise_k_hi_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_from_hi_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_noise_x_from_hi_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_from_hi_reco_c_xx_),T_x_lim_); axis image; axisnotick; title('real(T_noise_x_from_hi_reco_c_xx_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(S_noise_k_lo_p_wk_),S_k_lim_); axis(k_p_axis_); axis square; axisnotick; title('real(S_noise_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_from_lo_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_noise_x_from_lo_reco_c_xx_),S_x_lim_); axis image; axisnotick; title('real(S_noise_x_from_lo_reco_c_xx_)','Interpreter','none');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_noise_x_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_p(n_k_lo_p_r,k_lo_p_r_,n_w_lo_,n_w_lo_sum,real(T_noise_k_lo_p_wk_),S_k_lim_); axis(k_p_axis_); axis square; axisnotick; title('real(T_noise_k_lo_p_wk_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_from_lo_c_xx_),S_x_lim_); axis image; axisnotick; title('real(T_noise_x_from_lo_c_xx_)','Interpreter','none');
subplot(p_row,p_col,1+np);np=np+1;imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(T_noise_x_from_lo_reco_c_xx_),T_x_lim_); axis image; axisnotick; title('real(T_noise_x_from_lo_reco_c_xx_)','Interpreter','none');
%%%%;
fname_fig_pre = sprintf('%s/test_CryoLike_physical_vs_fourier_FIGC',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if flag_plot;
%%%%%%%%;

%%%%%%%%;
% Now test integration of one noisy-function against the other. ;
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Note difference between not filtered and yes filtered (k_hi_p_r_max): ')); end;
S_raw_dot_T_raw_noise_x_hi = sum(conj(S_noise_x_c_xx_).*T_noise_x_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_raw_noise_x_hi = sum(conj(S_noise_x_from_hi_c_xx_).*T_noise_x_c_xx_,'all')*weight_xx_c_;
S_raw_dot_T_fil_noise_x_hi = sum(conj(S_noise_x_c_xx_).*T_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_fil_noise_x_hi = sum(conj(S_noise_x_from_hi_c_xx_).*T_noise_x_from_hi_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_fil_noise_k_hi = sum(conj(S_noise_k_hi_p_wk_).*T_noise_k_hi_p_wk_.*weight_2d_hi_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_hi',S_fil_dot_T_fil_noise_k_hi,'S_raw_dot_T_raw_noise_x_hi',S_raw_dot_T_raw_noise_x_hi,' %<-- could be very large, especially if n_x_c is large');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_hi',S_fil_dot_T_fil_noise_k_hi,'S_fil_dot_T_raw_noise_x_hi',S_fil_dot_T_raw_noise_x_hi,' %<-- should be small?');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_hi',S_fil_dot_T_fil_noise_k_hi,'S_raw_dot_T_fil_noise_x_hi',S_raw_dot_T_fil_noise_x_hi,' %<-- should be small?');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_hi',S_fil_dot_T_fil_noise_k_hi,'S_fil_dot_T_fil_noise_x_hi',S_fil_dot_T_fil_noise_x_hi,' %<-- may be large if k_hi_p_r_max is close to bandlimit');
%%%%;
d2_S_raw_T_raw_noise_x_hi = sum(conj(S_noise_x_c_xx_-T_noise_x_c_xx_).*(S_noise_x_c_xx_-T_noise_x_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_raw_noise_x_hi = sum(conj(S_noise_x_from_hi_c_xx_-T_noise_x_c_xx_).*(S_noise_x_from_hi_c_xx_-T_noise_x_c_xx_),'all')*weight_xx_c_;
d2_S_raw_T_fil_noise_x_hi = sum(conj(S_noise_x_c_xx_-T_noise_x_from_hi_c_xx_).*(S_noise_x_c_xx_-T_noise_x_from_hi_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_fil_noise_x_hi = sum(conj(S_noise_x_from_hi_c_xx_-T_noise_x_from_hi_c_xx_).*(S_noise_x_from_hi_c_xx_-T_noise_x_from_hi_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_fil_noise_k_hi = sum(conj(S_noise_k_hi_p_wk_-T_noise_k_hi_p_wk_).*(S_noise_k_hi_p_wk_-T_noise_k_hi_p_wk_).*weight_2d_hi_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_hi',d2_S_fil_T_fil_noise_k_hi,'d2_S_raw_T_raw_noise_x_hi',d2_S_raw_T_raw_noise_x_hi,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_hi',d2_S_fil_T_fil_noise_k_hi,'d2_S_fil_T_raw_noise_x_hi',d2_S_fil_T_raw_noise_x_hi,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_hi',d2_S_fil_T_fil_noise_k_hi,'d2_S_raw_T_fil_noise_x_hi',d2_S_raw_T_fil_noise_x_hi,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_hi',d2_S_fil_T_fil_noise_k_hi,'d2_S_fil_T_fil_noise_x_hi',d2_S_fil_T_fil_noise_x_hi,' %<-- should be <1e-2');
%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% Note difference between not filtered and yes filtered (k_lo_p_r_max): ')); end;
S_raw_dot_T_raw_noise_x_lo = sum(conj(S_noise_x_c_xx_).*T_noise_x_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_raw_noise_x_lo = sum(conj(S_noise_x_from_lo_c_xx_).*T_noise_x_c_xx_,'all')*weight_xx_c_;
S_raw_dot_T_fil_noise_x_lo = sum(conj(S_noise_x_c_xx_).*T_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_fil_noise_x_lo = sum(conj(S_noise_x_from_lo_c_xx_).*T_noise_x_from_lo_c_xx_,'all')*weight_xx_c_;
S_fil_dot_T_fil_noise_k_lo = sum(conj(S_noise_k_lo_p_wk_).*T_noise_k_lo_p_wk_.*weight_2d_lo_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_lo',S_fil_dot_T_fil_noise_k_lo,'S_raw_dot_T_raw_noise_x_lo',S_raw_dot_T_raw_noise_x_lo,' %<-- could be very large, especially if n_x_c is large');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_lo',S_fil_dot_T_fil_noise_k_lo,'S_fil_dot_T_raw_noise_x_lo',S_fil_dot_T_raw_noise_x_lo,' %<-- should be small?');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_lo',S_fil_dot_T_fil_noise_k_lo,'S_raw_dot_T_fil_noise_x_lo',S_raw_dot_T_fil_noise_x_lo,' %<-- should be small?');
fnorm_disp(flag_verbose,'S_fil_dot_T_fil_noise_k_lo',S_fil_dot_T_fil_noise_k_lo,'S_fil_dot_T_fil_noise_x_lo',S_fil_dot_T_fil_noise_x_lo,' %<-- may be large if k_lo_p_r_max is close to bandlimit');
%%%%;
d2_S_raw_T_raw_noise_x_lo = sum(conj(S_noise_x_c_xx_-T_noise_x_c_xx_).*(S_noise_x_c_xx_-T_noise_x_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_raw_noise_x_lo = sum(conj(S_noise_x_from_lo_c_xx_-T_noise_x_c_xx_).*(S_noise_x_from_lo_c_xx_-T_noise_x_c_xx_),'all')*weight_xx_c_;
d2_S_raw_T_fil_noise_x_lo = sum(conj(S_noise_x_c_xx_-T_noise_x_from_lo_c_xx_).*(S_noise_x_c_xx_-T_noise_x_from_lo_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_fil_noise_x_lo = sum(conj(S_noise_x_from_lo_c_xx_-T_noise_x_from_lo_c_xx_).*(S_noise_x_from_lo_c_xx_-T_noise_x_from_lo_c_xx_),'all')*weight_xx_c_;
d2_S_fil_T_fil_noise_k_lo = sum(conj(S_noise_k_lo_p_wk_-T_noise_k_lo_p_wk_).*(S_noise_k_lo_p_wk_-T_noise_k_lo_p_wk_).*weight_2d_lo_wk_*(2*pi)^2);
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_lo',d2_S_fil_T_fil_noise_k_lo,'d2_S_raw_T_raw_noise_x_lo',d2_S_raw_T_raw_noise_x_lo,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_lo',d2_S_fil_T_fil_noise_k_lo,'d2_S_fil_T_raw_noise_x_lo',d2_S_fil_T_raw_noise_x_lo,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_lo',d2_S_fil_T_fil_noise_k_lo,'d2_S_raw_T_fil_noise_x_lo',d2_S_raw_T_fil_noise_x_lo,' %<-- could be large');
fnorm_disp(flag_verbose,'d2_S_fil_T_fil_noise_k_lo',d2_S_fil_T_fil_noise_k_lo,'d2_S_fil_T_fil_noise_x_lo',d2_S_fil_T_fil_noise_x_lo,' %<-- should be <1e-2');
%%%%%%%%;



