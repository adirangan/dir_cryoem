function test_heterogeneity_spurious_0();
%%%%%%%%;
% Sets up idealized single shell examples. ;
%%%%%%%%;

flag_verbose=1;
flag_disp=1; nf=0;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_c = 64;
x_c_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
x_c_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_c);
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c^3;
xxx_c_weight_ = (2*x_p_r_max/n_x_c)^3;
k_c_0_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
k_c_1_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
k_c_2_ = periodize(0:n_x_c-1,-n_x_c/2,+n_x_c/2)/2; %<-- box has diameter 2. ;
[k_c_0___,k_c_1___,k_c_2___] = ndgrid(k_c_0_,k_c_1_,k_c_2_); n_K_c = n_x_c^3;
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
 flag_verbose ...
,k_p_r_max ...
,k_eq_d ...
,'L' ...
) ; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
tmp_t = toc(tmp_t); disp(sprintf(' %% sample_sphere_7: time %0.2fs',tmp_t));
%%%%%%%%;

parameter = struct('type','parameter');
parameter.sigma_g = 0.125;
parameter.str_sigma_g = sprintf('sg%d',round(parameter.sigma_g*1000));
parameter.mu_g = 0.35;
parameter.str_mu_g = sprintf('sm%d',round(parameter.mu_g*1000));
parameter.backgr_v = +0.125 + i*0.125;
parameter.str_backgr_v = sprintf('bv%di%d',8+round(real(parameter.backgr_v)*8),8+round(imag(parameter.backgr_v)*8));
parameter.square_r = 0.75*1/sqrt(2);
parameter.str_square_r = sprintf('sr%d',round(parameter.square_r*sqrt(2)*8));
parameter.square_v = +1.0 + i*0.500;
parameter.str_square_v = sprintf('sv%di%d',8+round(real(parameter.square_v)*8),8+round(imag(parameter.square_v)*8));
parameter.circle_r = 0.25;
parameter.str_circle_r = sprintf('cr%d',round(parameter.circle_r*8));
parameter.circle_v = 1.0 + i*0.500;
parameter.str_circle_v = sprintf('cv%di%d',8+round(real(parameter.circle_v)*8),8+round(imag(parameter.circle_v)*8));
parameter.circle_c_0 = parameter.square_r-1.25*parameter.circle_r;
parameter.circle_c_1 = parameter.square_r-1.25*parameter.circle_r;
parameter.e_u = parameter.square_r*sqrt(2);
parameter.e_v = 2*parameter.circle_r;
parameter.x_shift = 2.0;
parameter.view_az = +60;
parameter.view_el = +30;
parameter.infix = sprintf( ...
 '%s%s%s%s%s%s%s' ...
,parameter.str_sigma_g ...
,parameter.str_mu_g ...
,parameter.str_backgr_v ...
,parameter.str_square_r ...
,parameter.str_square_v ...
,parameter.str_circle_r ...
,parameter.str_circle_v ...
);

%{
tmp_index_ = efind(k_p_r_all_==max(k_p_r_all_)); n_k_sub = numel(tmp_index_);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+tmp_index_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+tmp_index_);
k_c_0_sub_ = k_c_0_all_(1+tmp_index_);
k_c_1_sub_ = k_c_1_all_(1+tmp_index_);
k_c_2_sub_ = k_c_2_all_(1+tmp_index_);
a_k_p_sub_ = shell_val_A(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
b_k_p_sub_ = shell_val_B(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
c_k_p_sub_ = shell_val_C(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
figure(1+nf);nf=nf+1;figmed;
subplot(1,3,1); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(a_k_p_sub_),[-1,+1],colormap_beach(),0);
subplot(1,3,2); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(b_k_p_sub_),[-1,+1],colormap_beach(),0);
subplot(1,3,3); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(c_k_p_sub_),[-1,+1],colormap_beach(),0);
%}

a_k_p_all_ = sphere_val_A(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
b_k_p_all_ = sphere_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
c_k_p_all_ = sphere_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
%%%%%%%%;
%{
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_all_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,a_x_c_base_(:).*xxx_c_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_all_ time %0.2fs',tmp_t));
%}
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_c_xxx_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_all_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_c_xxx_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
b_x_c_xxx_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,b_k_p_all_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: b_x_c_xxx_ time %0.2fs',tmp_t));
%%%%%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
c_x_c_xxx_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,c_k_p_all_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_c,x_c_0___(:)/eta,x_c_1___(:)/eta,x_c_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: c_x_c_xxx_ time %0.2fs',tmp_t));
%%%%%%%%;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGA',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;figmed;
fontsize_use = 24;
prctile_ = [97.5];%prctile_ = [97.5, 98.5, 99.5];
subplot(1,3,1); isosurface_f_x_u_0(a_x_c_xxx_ + parameter.x_shift,prctile_); title('A'); view(parameter.view_az,parameter.view_el);
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,2); isosurface_f_x_u_0(b_x_c_xxx_ + parameter.x_shift,prctile_); title('B'); view(parameter.view_az,parameter.view_el);
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,3); isosurface_f_x_u_0(c_x_c_xxx_ + parameter.x_shift,prctile_); title('C'); view(parameter.view_az,parameter.view_el);
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function val_A_all_ = shell_val_A(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
circle_c_0 = parameter.circle_c_0;
circle_c_1 = parameter.circle_c_1;
backgr_v_all_ = real(backgr_v) + i*imag(backgr_v)*(k_c_2_all_> 0) - i*imag(backgr_v)*(k_c_2_all_< 0);
square_v_all_ = real(square_v) + i*imag(square_v)*(k_c_2_all_> 0) - i*imag(square_v)*(k_c_2_all_< 0);
circle_v_all_ = real(circle_v) + i*imag(circle_v)*(k_c_2_all_> 0) - i*imag(circle_v)*(k_c_2_all_< 0);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
k_n_0_all_ = k_c_0_all_./max(1e-12,k_p_r_all_);
k_n_1_all_ = k_c_1_all_./max(1e-12,k_p_r_all_);
k_n_r01_p_all_ = sqrt((k_n_0_all_ - circle_c_0).^2 + (k_n_1_all_ - circle_c_1).^2);
k_n_r01_n_all_ = sqrt((k_n_0_all_ + circle_c_0).^2 + (k_n_1_all_ + circle_c_1).^2);
%k_n_s_all_ = (k_n_0_all_.^8 + k_n_1_all_.^8).^(1/8);
k_n_s_all_ = max(abs(k_n_0_all_),abs(k_n_1_all_));
val_A_all_ = backgr_v_all_ ...
 + square_v_all_.*(k_n_s_all_<=square_r) ...
 + circle_v_all_.*(k_n_r01_p_all_<=circle_r).*(k_c_2_all_> 0) ...
 + circle_v_all_.*(k_n_r01_n_all_<=circle_r).*(k_c_2_all_< 0) ...
  ;

function val_B_all_ = shell_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
circle_c_0 = parameter.circle_c_0;
circle_c_1 = parameter.circle_c_1;
circle_r = square_r*sqrt(2);
circle_c_0 = 0;
circle_c_1 = 0;
backgr_v_all_ = real(backgr_v) + i*imag(backgr_v)*(k_c_2_all_> 0) - i*imag(backgr_v)*(k_c_2_all_< 0);
circle_v_all_ = real(circle_v) + i*imag(circle_v)*(k_c_2_all_> 0) - i*imag(circle_v)*(k_c_2_all_< 0);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
k_n_0_all_ = k_c_0_all_./max(1e-12,k_p_r_all_);
k_n_1_all_ = k_c_1_all_./max(1e-12,k_p_r_all_);
k_n_r01_p_all_ = sqrt((k_n_0_all_ - circle_c_0).^2 + (k_n_1_all_ - circle_c_1).^2);
k_n_r01_n_all_ = sqrt((k_n_0_all_ + circle_c_0).^2 + (k_n_1_all_ + circle_c_1).^2);
val_B_all_ = backgr_v_all_ ...
 + circle_v_all_.*(k_n_r01_p_all_<=circle_r).*(k_c_2_all_> 0) ...
 + circle_v_all_.*(k_n_r01_n_all_<=circle_r).*(k_c_2_all_< 0) ...
  ;

function val_C_all_ = shell_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
circle_c_0 = parameter.circle_c_0;
circle_c_1 = parameter.circle_c_1;
e_u = parameter.e_u;
e_v = parameter.e_v;
backgr_v_all_ = real(backgr_v) + i*imag(backgr_v)*(k_c_2_all_> 0) - i*imag(backgr_v)*(k_c_2_all_< 0);
square_v_all_ = real(square_v) + i*imag(square_v)*(k_c_2_all_> 0) - i*imag(square_v)*(k_c_2_all_< 0);
circle_v_all_ = real(circle_v) + i*imag(circle_v)*(k_c_2_all_> 0) - i*imag(circle_v)*(k_c_2_all_< 0);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
k_n_0_all_ = k_c_0_all_./max(1e-12,k_p_r_all_);
k_n_1_all_ = k_c_1_all_./max(1e-12,k_p_r_all_);
k_n_u_all_ = (k_n_0_all_ + k_n_1_all_)/sqrt(2);
k_n_v_all_ = (k_n_0_all_ - k_n_1_all_)/sqrt(2);
k_n_e_all_ = sqrt(k_n_u_all_.^2./e_u.^2 + k_n_v_all_.^2./e_v.^2);
k_n_r01_p_all_ = sqrt((k_n_0_all_ - circle_c_0).^2 + (k_n_1_all_ - circle_c_1).^2);
k_n_r01_n_all_ = sqrt((k_n_0_all_ + circle_c_0).^2 + (k_n_1_all_ + circle_c_1).^2);
val_C_all_ = backgr_v_all_ ...
 + square_v_all_.*(k_n_e_all_<=1) ...
 + circle_v_all_.*(k_n_r01_p_all_<=circle_r).*(k_c_2_all_> 0) ...
 + circle_v_all_.*(k_n_r01_n_all_<=circle_r).*(k_c_2_all_< 0) ...
  ;

function val_A_all_ = sphere_val_A(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
sigma_g = parameter.sigma_g; mu_g = parameter.mu_g;
g_r_all_ = exp(-(k_p_r_all_-mu_g*k_p_r_max).^2/(2*sigma_g.^2.*k_p_r_max.^2));
val_A_all_ = g_r_all_.*shell_val_A(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);

function val_B_all_ = sphere_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
sigma_g = parameter.sigma_g; mu_g = parameter.mu_g;
g_r_all_ = exp(-(k_p_r_all_-mu_g*k_p_r_max).^2/(2*sigma_g.^2.*k_p_r_max.^2));
val_B_all_ = g_r_all_.*shell_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);

function val_C_all_ = sphere_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
sigma_g = parameter.sigma_g; mu_g = parameter.mu_g;
g_r_all_ = exp(-(k_p_r_all_-mu_g*k_p_r_max).^2/(2*sigma_g.^2.*k_p_r_max.^2));
val_C_all_ = g_r_all_.*shell_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);


