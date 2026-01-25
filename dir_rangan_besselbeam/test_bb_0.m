%%%%%%%%;
% Paraxial approximation to shrodinger-equation. ;
% Beam headed in z-direction. ;
% Using polar (i.e., radial) description of x-y-plane for now. ;
%%%%%%%%;

flag_verbose=1; nf=0;
flag_disp = 1; nf=0; flag_replot = 1;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem/dir_rangan_besselbeam/dir_jpg',string_root);

%%%%%%%%;
% This is basically the same as test_basis_3.m ;
% (using the same scale in fourier-space). ;
%%%%%%%%;

%%%%%%%%;
% Visualizing the different basis functions. ;
%%%%%%%%;
f_res = 1;
k_int = 48*f_res;
k_eq_d_double = 0.5;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u_pack = 64*f_res;
n_x_M_u = n_x_u_pack;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1);
x_u_0_ = transpose(x_u_0_(1:n_x_M_u));
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_M_u+1);
x_u_1_ = transpose(x_u_1_(1:n_x_M_u));
[x_u_0_xx__,x_u_1_xx__] = ndgrid(x_u_0_,x_u_1_);
x_u_r_xx__ = sqrt(x_u_0_xx__.^2 + x_u_1_xx__.^2);
x_u_w_xx__ = atan2(x_u_1_xx__,x_u_0_xx__);
dx = diameter_x_c/n_x_M_u;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi);
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
n_w_max = 2*(l_max_max+1);
template_k_eq_d = -1; %<-- default value of n_w_max. ;
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
%%%%%%%%;
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_k_wk_ ...
] = ...
get_weight_2d_1( ...
 max(0,flag_verbose-1) ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,template_k_eq_d ...
,n_w_0in_ ...
);
if (flag_verbose>0); disp(sprintf(' %% sum(weight_2d_k_p_r_) - area %0.16f',sum(weight_2d_k_p_r_) - pi*k_p_r_max^2)); end;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
k_p_r_wk_ = reshape(transpose(repmat(k_p_r_,[1,n_w_max])),[n_w_sum,1]);

k_p_r_max_disp = max(32,ceil(k_p_r_max));
k_q_k_lim_ = [0,k_p_r_max_disp/2];
k_q_q_lim_ = [-n_w_max,+n_w_max]/4;

flag_disp=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
for nl=0:3-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if nl==0;
M_x_c_xx__ = zeros(n_x_M_u,n_x_M_u);
sigma_x = 0.5/k_p_r_max;
delta_x_ = (0.125)*[cos(0*pi/4);sin(0*pi/4)];
y_u_0_xx__ = x_u_0_xx__ - delta_x_(1+0);
y_u_1_xx__ = x_u_1_xx__ - delta_x_(1+1);
y_p_r__ = sqrt(y_u_0_xx__.^2 + y_u_1_xx__.^2);
y_p_w__ = atan2(y_u_1_xx__,y_u_0_xx__);
M_x_c_xx__ = 1/(2*pi)/sigma_x^2 * exp(-y_p_r__.^2/(2*sigma_x^2));
M_k_p_wk_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,M_x_c_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
)*sqrt(n_x_M_u^2)*dx^2 ;
N_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
M_x_c_xx__l2 = sum(abs(M_x_c_xx__).^2 .* dx.^2,'all');
N_x_c_l2 = sum(abs(N_x_c_xx__).^2 .* dx.^2,'all');
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs N_x_c_l2: %0.16f',fnorm(M_x_c_xx__l2 - N_x_c_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_p_l2 = sum(abs(M_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_l2 = sum(abs(M_k_q_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_q_l2)/fnorm(M_x_c_xx__l2))); end;
N_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_wk_);
N_k_p_l2 = sum(abs(N_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs N_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - N_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
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
sigma_k = 0.5;
%delta_k_ = +8.0*[cos(3*pi/4);sin(3*pi/4)];
delta_k_ = +8.0*[cos(0*pi/4);sin(0*pi/4)];
j_c_0_ = k_c_0_ - delta_k_(1+0);
j_c_1_ = k_c_1_ - delta_k_(1+1);
j_p_r_ = sqrt(j_c_0_.^2 + j_c_1_.^2);
j_p_w_ = atan2(j_c_1_,j_c_0_);
M_k_p_wk_ = 1/(2*pi)/sigma_k^2 * exp(-j_p_r_.^2/(2*sigma_k^2));
M_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
M_x_c_xx__l2 = sum(abs(M_x_c_xx__).^2 .* dx.^2,'all');
M_k_p_l2 = sum(abs(M_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_l2 = sum(abs(M_k_q_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_q_l2)/fnorm(M_x_c_xx__l2))); end;
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
delta_q_ = [6;8];
j_q_q_ = k_q_q_ - delta_q_(1+0);
j_q_k_ = k_q_k_ - delta_q_(1+1);
M_k_q_wk_ = 1/(2*pi)/(sigma_k*sigma_q) * exp(-j_q_q_.^2/(2*sigma_q^2)) .* exp(-j_q_k_.^2/(2*sigma_k^2));
M_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_wk_);
M_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
M_x_c_xx__l2 = sum(abs(M_x_c_xx__).^2 .* dx.^2,'all');
M_k_p_l2 = sum(abs(M_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_q_l2 = sum(abs(M_k_q_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_q_l2)/fnorm(M_x_c_xx__l2))); end;
end;%if nl==2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_disp;
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024*2,1024*0.60]);
gamma_use = 2.0;
colormap(colormap_pm().^gamma_use);
p_row = 1; p_col = 5;
fontsize_use = 18;
%%%%;
subplot(p_row,p_col,[1,2]);
clim_ = [-1,+1]*prctile(abs(real(M_x_c_xx__)),100,'all');
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(M_x_c_xx__),clim_,colormap_pm.^gamma_use); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,[3,4]);
clim_ = [-1,+1]*prctile(abs(real(M_k_p_wk_)),100,'all');
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_p_wk_),clim_,colormap_pm.^gamma_use); axis image;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_disp*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),k_p_r_max_disp*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(k_p_r_max_disp*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(k_p_r_max_disp*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k_{1}$','Interpreter','latex');ylabel('$k_{2}$','Interpreter','latex');title('$\Re\{\hat{A}(\vec{k})\}$','Interpreter','latex');
%%%%;
subplot(p_row,p_col,[5]);
clim_ = [-1,+1]*prctile(abs(real(M_k_q_wk_)),100,'all');
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(M_k_q_wk_),clim_,colormap_pm.^gamma_use); axisnotick;
xl = arrayfun(@(x)xline(x,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),max(k_q_k_lim_)*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
yl = arrayfun(@(y)yline(y,'-','Color',0.85*[1,1,1],'LineWidth',2,'Alpha',0.25),max(k_q_q_lim_)*linspace(-half_diameter_x_c,+half_diameter_x_c,1+8));
xlim(k_q_k_lim_);ylim(k_q_q_lim_);
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',fix(max(k_q_k_lim_)*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'YTick',fix(max(k_q_q_lim_)*[-half_diameter_x_c,0,+half_diameter_x_c]));
set(gca,'TickLength',[0,0]);
xlabel('$k$','Interpreter','latex');ylabel('$q$','Interpreter','latex');title('$\Re\{a(k;q)\}$','Interpreter','latex');
%%%%;
fname_fig_pre = sprintf('%s/pm_fig_image_basis_nl%d_FIGB',dir_base,nl);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_disp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% simple plane-wave. ;
%%%%%%%%;
str_type = 'ringlens';
%%%%;
M_x_c_xx__ = zeros(n_x_M_u,n_x_M_u);
if strcmp(str_type,'gaussian');
sigma_x = 4*0.5/k_p_r_max;
delta_x_ = (0.125)*[cos(0*pi/4);sin(0*pi/4)];
y_u_0_xx__ = x_u_0_xx__ - delta_x_(1+0);
y_u_1_xx__ = x_u_1_xx__ - delta_x_(1+1);
y_p_r__ = sqrt(y_u_0_xx__.^2 + y_u_1_xx__.^2);
y_p_w__ = atan2(y_u_1_xx__,y_u_0_xx__);
M_x_c_xx__ = 1/(2*pi)/sigma_x^2 * exp(-y_p_r__.^2/(2*sigma_x^2));
end;%if strcmp(str_type,'gaussian');
if strcmp(str_type,'planewave');
M_x_c_xx__ = exp(i*2*pi*0.125*(x_u_0_xx__ + x_u_1_xx__));
end;%if strcmp(str_type,'planewave');
if strcmp(str_type,'quadlens');
M_x_c_xx__ = exp(i*2*pi*(0 - 1.5*(x_u_r_xx__.^2)));
end;%if strcmp(str_type,'quadlens');
if strcmp(str_type,'linelens');
M_x_c_xx__ = exp(i*2*pi*(0 - 0.5*(x_u_r_xx__.^1)));
end;%if strcmp(str_type,'linelens');
if strcmp(str_type,'ringlens');
M_x_c_xx__ = exp(i*2*pi*(0 - 1.5*(x_u_r_xx__.^2))).*(x_u_r_xx__>=0.5).*(x_u_r_xx__<=0.6);
end;%if strcmp(str_type,'ringlens');
%%%%;
M_k_p_wk_ = ...
interp_x_c_to_k_p_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,M_x_c_xx__ ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
)*sqrt(n_x_M_u^2)*dx^2 ;
%%%%;
N_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
%%%%;
M_x_c_xx__l2 = sum(abs(M_x_c_xx__).^2 .* dx.^2,'all');
N_x_c_l2 = sum(abs(N_x_c_xx__).^2 .* dx.^2,'all');
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs N_x_c_l2: %0.16f',fnorm(M_x_c_xx__l2 - N_x_c_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_p_l2 = sum(abs(M_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
M_k_q_wk_ = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_);
M_k_q_l2 = sum(abs(M_k_q_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs M_k_q_l2: %0.16f',fnorm(M_x_c_xx__l2 - M_k_q_l2)/fnorm(M_x_c_xx__l2))); end;
N_k_p_wk_ = interp_q_to_p(n_k_p_r,n_w_,n_w_sum,M_k_q_wk_);
N_k_p_l2 = sum(abs(N_k_p_wk_).^2 .* weight_2d_k_wk_,'all')*(4*pi^2);
if (flag_verbose>0); disp(sprintf(' %% M_x_c_xx__l2 vs N_k_p_l2: %0.16f',fnorm(M_x_c_xx__l2 - N_k_p_l2)/fnorm(M_x_c_xx__l2))); end;
%%%%%%%%;

%%%%%%%%;
lambda = 1.0; %<-- here measure z in terms of lambda. ;
dzl = 1.0/64.0;
zl_max = 2.0;
n_zl = zl_max/dzl;
tmp_M_k_p_wk_ = M_k_p_wk_;
%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
fontsize_use = 18;
gamma_use = 2.0;
colormap(colormap_pm().^gamma_use);
clim_ = [];
%%%%;
for nzl=0:n_zl-1;
tmp_M_k_p_wk_ = tmp_M_k_p_wk_.*exp(-i*pi*lambda*dzl*k_p_r_wk_.^2); %<-- paraxial fourier propagator ;
tmp_N_x_c_xx__ = ...
interp_k_p_to_x_c_xxnufft( ...
 n_x_M_u ...
,diameter_x_c ...
,n_x_M_u ...
,diameter_x_c ...
,n_k_p_r ...
,k_p_r_ ...
,n_w_ ...
,tmp_M_k_p_wk_.*weight_2d_k_wk_*(4*pi^2) ...
)/(sqrt(n_x_M_u^2)*dx^2)*n_w_max^2*4 ;
if isempty(clim_); clim_ = [-1,+1]*prctile(abs(real(tmp_N_x_c_xx__)),100,'all'); end;
subplot(1,1,1); cla;
imagesc_c(n_x_M_u,x_u_0_,n_x_M_u,x_u_1_,real(tmp_N_x_c_xx__),clim_,colormap_pm.^gamma_use); axis image;
set(gca,'FontSize',fontsize_use);
set(gca,'XTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'YTick',[-half_diameter_x_c,0,+half_diameter_x_c]);
set(gca,'TickLength',[0,0]);
xlabel('$x_{1}$','Interpreter','latex');ylabel('$x_{2}$','Interpreter','latex');title('$\Re\{A(\vec{x})\}$','Interpreter','latex');
drawnow();
end;%for nzl=0:n_zl-1;
%%%%;



