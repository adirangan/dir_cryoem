clear;

flag_verbose=1; nf=0;
flag_disp=1;
flag_replot = 0;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

%%%%%%%%;
% Visualizing the different basis functions. ;
%%%%%%%%;
f_res = 1; flag_tensor_vs_adap = 0;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64*f_res;
n_x_M_u = n_x_u_pack;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u);
[x_u_0__,x_u_1__] = ndgrid(x_u_0_,x_u_1_);
x_p_r__ = sqrt(x_u_0__.^2 + x_u_1__.^2);
x_p_w__ = atan2(x_u_1__,x_u_0__);
dx = diameter_x_c/n_x_M_u;
%%%%%%%%;
k_p_r_max = 48/(2*pi)*f_res; k_eq_d = 0.5/(2*pi);
[ ...
 n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
] = ...
get_weight_3d_1( ...
 max(0,flag_verbose-1) ...
,k_p_r_max ...
,k_eq_d ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_3d_k_p_r_)*4*pi - volume: %0.16f',sum(weight_3d_k_p_r_)*4*pi - 4/3*pi*k_p_r_max^3)); end;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max); %<-- typically sufficient for 2-3 digits of precision. ;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
l_max_max = max(l_max_);
n_w_max = fix(1/k_eq_d/(2*pi))*2*(l_max_max+1);
if flag_tensor_vs_adap==0; template_k_eq_d = k_eq_d; n_w_0in_ = []; end;
if flag_tensor_vs_adap==1; n_w_0in_ = n_w_max*ones(n_k_p_r,1); template_k_eq_d = -1; end; %<-- default value of n_w_max. ;
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_all_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 max(0,flag_verbose-1) ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_k_p_r_) - area %0.16f',sum(weight_2d_k_p_r_) - pi*k_p_r_max^2)); end;
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_k_all_)*(4*pi^2) - area %0.16f',sum(weight_2d_k_all_)*(4*pi^2) - pi*k_p_r_max^2)); end;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% disp(sprintf(' %% returning before visualizing the different noise-distributions')); return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

str_type = 'random';%str_type = 'gaussian';
L_x_c_ = zeros(n_x_M_u,n_x_M_u);
if strcmp(str_type,'gaussian');
sigma_x = 1/32;
delta_x_ = [0.3;0.12];
y_u_0__ = x_u_0__ - delta_x_(1+0);
y_u_1__ = x_u_1__ - delta_x_(1+1);
y_p_r__ = sqrt(y_u_0__.^2 + y_u_1__.^2);
y_p_w__ = atan2(y_u_1__,y_u_0__);
L_x_c_ = 1/(2*pi)/sigma_x^2 * exp(-y_p_r__.^2/(2*sigma_x^2));
end;%if strcmp(str_type,'gaussian');
if strcmp(str_type,'random');
sigma_x = 10.00; %<-- variance per unit area. ;
L_x_c_ = (sigma_x/diameter_x_c)*randn(n_x_M_u,n_x_M_u);
end;%if strcmp(str_type,'random');
%%%%;
L_k_p_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,L_x_c_ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
)*sqrt(n_x_M_u^2)*dx^2 ;
M_x_c_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,L_k_p_.*weight_2d_k_all_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*2 ;
M_k_p_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,M_x_c_ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
)*sqrt(n_x_M_u^2)*dx^2 ;
L_x_c_l2 = sum(abs(L_x_c_).^2 .* dx.^2,'all');
M_x_c_l2 = sum(abs(M_x_c_).^2 .* dx.^2,'all');
L_k_p_l2 = sum(abs(L_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
M_k_p_l2 = sum(abs(M_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% L_x_c_l2 vs M_x_c_l2: %0.16f',fnorm(L_x_c_l2 - M_x_c_l2)/fnorm(L_x_c_l2))); end;
if (flag_verbose>0); disp(sprintf(' %% L_k_p_l2 vs M_k_p_l2: %0.16f',fnorm(L_k_p_l2 - M_k_p_l2)/fnorm(L_k_p_l2))); end;
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_l2 - M_k_p_l2)/fnorm(M_x_c_l2))); end;
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_k_q_l2 = sum(abs(M_k_q_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_k_p_l2 vs M_k_q_l2: %0.16f',fnorm(M_k_p_l2 - M_k_q_l2)/fnorm(M_k_p_l2))); end;
N_k_p_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_);
N_k_p_l2 = sum(abs(N_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_k_p_l2 vs N_k_p_l2: %0.16f',fnorm(M_k_p_l2 - N_k_p_l2)/fnorm(M_k_p_l2))); end;
%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
p_row = 2; p_col = 3; np=0;
fontsize_use = 18;
k_p_r_max_use = k_p_r_max;
M_x_clim_ = [-1,+1]*prctile(abs(real(M_x_c_)),100,'all');
M_k_clim_ = [-1,+1]*prctile(abs(real(M_k_p_)),100,'all');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(L_x_c_),M_x_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$ original','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(L_k_p_),M_k_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Re\{\hat{A}(\vec{k})\}$ filtered$\times 1$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(L_k_p_),M_k_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Im\{\hat{A}(\vec{k})\}$ filtered$\times 1$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(M_x_c_),M_x_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$ filtered$\times 2$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),M_k_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Re\{\hat{A}(\vec{k})\}$ filtered$\times 3$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(M_k_p_),M_k_clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Im\{\hat{A}(\vec{k})\}$ filtered$\times 3$','Interpreter','latex');
%%%%;
fname_fig_pre = sprintf('/%s/rangan/dir_cryoem/dir_pm_manuscript/dir_pm_fig/pm_fig_image_basis_noise_FIGA',string_root);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if ~exist(fname_fig_jpg);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml
p_row = 1; p_col = 1; np=0;
fontsize_use = 18;
k_p_r_max_use = k_p_r_max;
E_k_clim_ = [-1,+1]*prctile(abs(L_k_p_-M_k_p_),100,'all');
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(L_k_p_-M_k_p_),M_k_clim_,colormap_81s); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$|\hat{A}(\vec{k})|$ filtered $\times 1$ minus $\times 3$','Interpreter','latex');
%%%%;
fname_fig_pre = sprintf('/%s/rangan/dir_cryoem/dir_pm_manuscript/dir_pm_fig/pm_fig_image_basis_noise_FIGB',string_root);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if ~exist(fname_fig_jpg);
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
disp(sprintf(' %% returning before visualizing the different choice of basis')); return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%k_q_k_lim_ = [0,k_p_r_max/2];
k_p_r_max_use = 32;
k_q_k_lim_ = [0,k_p_r_max_use/2];
k_q_q_lim_ = [-n_w_max,+n_w_max]/4;

for nl=0:3-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nl==0;
M_x_c_ = zeros(n_x_M_u,n_x_M_u);
sigma_x = 1/32;
delta_x_ = [0.3;0.12];
y_u_0__ = x_u_0__ - delta_x_(1+0);
y_u_1__ = x_u_1__ - delta_x_(1+1);
y_p_r__ = sqrt(y_u_0__.^2 + y_u_1__.^2);
y_p_w__ = atan2(y_u_1__,y_u_0__);
M_x_c_ = 1/(2*pi)/sigma_x^2 * exp(-y_p_r__.^2/(2*sigma_x^2));
M_k_p_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,M_x_c_ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
)*sqrt(n_x_M_u^2)*dx^2 ;
N_x_c_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_.*weight_2d_k_all_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*2 ;
M_x_c_l2 = sum(abs(M_x_c_).^2 .* dx.^2,'all');
N_x_c_l2 = sum(abs(N_x_c_).^2 .* dx.^2,'all');
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs N_x_c_l2: %0.16f',fnorm(M_x_c_l2 - N_x_c_l2)/fnorm(M_x_c_l2))); end;
M_k_p_l2 = sum(abs(M_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_l2 - M_k_p_l2)/fnorm(M_x_c_l2))); end;
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_k_q_l2 = sum(abs(M_k_q_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_l2 - M_k_q_l2)/fnorm(M_x_c_l2))); end;
N_k_p_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_);
N_k_p_l2 = sum(abs(N_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs N_k_p_l2: %0.16f',fnorm(M_x_c_l2 - N_k_p_l2)/fnorm(M_x_c_l2))); end;
end;%if nl==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nl==1;
k_c_0_ = zeros(n_w_sum,1);
k_c_1_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
tmp_gamma_ = transpose(linspace(0,2*pi,1+n_w)); tmp_gamma_ = tmp_gamma_(1:end-1);
k_c_0_(1+tmp_index_) = k_p_r*cos(tmp_gamma_);
k_c_1_(1+tmp_index_) = k_p_r*sin(tmp_gamma_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
sigma_k = k_p_r_max/64;
delta_k_ = [4,3]*2;
j_c_0_ = k_c_0_ - delta_k_(1+0);
j_c_1_ = k_c_1_ - delta_k_(1+1);
j_p_r_ = sqrt(j_c_0_.^2 + j_c_1_.^2);
j_p_w_ = atan2(j_c_1_,j_c_0_);
M_k_p_ = 1/(2*pi)/sigma_k^2 * exp(-j_p_r_.^2/(2*sigma_k^2));
M_x_c_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_.*weight_2d_k_all_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*2 ;
M_x_c_l2 = sum(abs(M_x_c_).^2 .* dx.^2,'all');
M_k_p_l2 = sum(abs(M_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_l2 - M_k_p_l2)/fnorm(M_x_c_l2))); end;
M_k_q_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_);
M_k_q_l2 = sum(abs(M_k_q_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_l2 - M_k_q_l2)/fnorm(M_x_c_l2))); end;
end;%if nl==1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nl==2;
k_q_q_ = zeros(n_w_sum,1);
k_q_k_ = zeros(n_w_sum,1);
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
tmp_q_ = periodize([0:n_w-1],-n_w/2,+n_w/2);
k_q_q_(1+tmp_index_) = tmp_q_;
k_q_k_(1+tmp_index_) = k_p_r;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
sigma_k = k_p_r_max*2/256;
sigma_q = n_w_max/2/256;
delta_q_ = [5;8];
j_q_q_ = k_q_q_ - delta_q_(1+0);
j_q_k_ = k_q_k_ - delta_q_(1+1);
M_k_q_ = 1/(2*pi)/(sigma_k*sigma_q) * exp(-j_q_q_.^2/(2*sigma_q^2)) .* exp(-j_q_k_.^2/(2*sigma_k^2));
M_k_p_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_);
M_x_c_ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_.*weight_2d_k_all_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*2 ;
M_x_c_l2 = sum(abs(M_x_c_).^2 .* dx.^2,'all');
M_k_p_l2 = sum(abs(M_k_p_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_l2 - M_k_p_l2)/fnorm(M_x_c_l2))); end;
M_k_q_l2 = sum(abs(M_k_q_).^2 .* weight_2d_k_all_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_l2 - M_k_q_l2)/fnorm(M_x_c_l2))); end;
end;%if nl==2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

if flag_disp;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024*2,1024*0.60]);figpm;
p_row = 1; p_col = 5;
fontsize_use = 18;
%%%%;
subplot(p_row,p_col,[1,2]);
clim_ = [-1,+1]*prctile(abs(real(M_x_c_)),100,'all');
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(M_x_c_),clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,[3,4]);
clim_ = [-1,+1]*prctile(abs(real(M_k_p_)),100,'all');
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_),clim_,colormap_pm); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_use*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_use*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Re\{\hat{A}(\vec{k})\}$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,[5]);
clim_ = [-1,+1]*prctile(abs(real(M_k_q_)),100,'all');
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_q_),clim_,colormap_pm); axisnotick;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),max(k_q_k_lim_)*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),max(k_q_q_lim_)*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
xlim(k_q_k_lim_);ylim(k_q_q_lim_);
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(max(k_q_k_lim_)*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(max(k_q_q_lim_)*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k$','Interpreter','latex');ylabel('$q$','Interpreter','latex');title('$\Re\{a(k;q)\}$','Interpreter','latex');
%%%%;
fname_fig_pre = sprintf('/%s/rangan/dir_cryoem/dir_pm_manuscript/dir_pm_fig/pm_fig_image_basis_nl%d_FIGA',string_root,nl);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
if flag_replot | ~exist(fname_fig_jpg);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
print('-depsc',fname_fig_eps);
end;%if ~exist(fname_fig_jpg);
%%%%;
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nl=0:3-1;





