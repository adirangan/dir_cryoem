function test_heterogeneity_spurious_1();
%%%%%%%%;
% Sets up idealized single shell examples. ;
%%%%%%%%;

str_thisfunction = 'test_heterogeneity_spurious_1';

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
k_p_r_max = (1.0)*48/(2*pi); k_eq_d = sqrt(1.0)*1.0/(2*pi);
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
parameter.str_shape = 's2';
parameter.sigma_g = 0.125/2;
parameter.str_sigma_g = sprintf('sg%d',round(parameter.sigma_g*1000));
parameter.mu_g = 0.450;
parameter.str_mu_g = sprintf('sm%d',round(parameter.mu_g*1000));
parameter.backgr_v = +0.0 + i*0.0;
parameter.str_backgr_v = sprintf('bv%di%d',8+round(real(parameter.backgr_v)*8),8+round(imag(parameter.backgr_v)*8));
parameter.square_r = 0.75*1/sqrt(2);
parameter.str_square_r = sprintf('sr%d',round(parameter.square_r*sqrt(2)*8));
parameter.square_v = +1.0 + i*0.500;
parameter.str_square_v = sprintf('sv%di%d',8+round(real(parameter.square_v)*8),8+round(imag(parameter.square_v)*8));
parameter.circle_r = 0.85;
parameter.str_circle_r = sprintf('cr%d',round(parameter.circle_r*8));
parameter.circle_v = +0.5 + i*0.250;
parameter.str_circle_v = sprintf('cv%di%d',8+round(real(parameter.circle_v)*8),8+round(imag(parameter.circle_v)*8));
parameter.circle_c_0 = sqrt(2.0);
parameter.circle_c_1 = sqrt(0.5);
parameter.infix = sprintf( ...
 '%s%s%s%s%s%s%s%s' ...
,parameter.str_shape ...
,parameter.str_sigma_g ...
,parameter.str_mu_g ...
,parameter.str_backgr_v ...
,parameter.str_square_r ...
,parameter.str_square_v ...
,parameter.str_circle_r ...
,parameter.str_circle_v ...
);

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

flag_disp=1;
if flag_disp;
tmp_index_ = efind(k_p_r_all_==max(k_p_r_all_)); n_k_sub = numel(tmp_index_);
k_p_azimu_b_sub_ = k_p_azimu_b_all_(1+tmp_index_);
k_p_polar_a_sub_ = k_p_polar_a_all_(1+tmp_index_);
k_c_0_sub_ = k_c_0_all_(1+tmp_index_);
k_c_1_sub_ = k_c_1_all_(1+tmp_index_);
k_c_2_sub_ = k_c_2_all_(1+tmp_index_);
a_k_p_sub_ = shell_val_A(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
b_k_p_sub_ = shell_val_B(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
c_k_p_sub_ = shell_val_C(parameter,k_p_r_max,n_k_sub,k_c_0_sub_,k_c_1_sub_,k_c_2_sub_);
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGB',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;figbig;
p_row = 2; p_col = 6; np=0;
a_k_p_sub_lim_ = max(abs(a_k_p_sub_),[],'all')*[-1,+1];
%%%%;
for nl=0:2;
if nl==0; tmp_d_k_p_sub_ = a_k_p_sub_; tmp_d_k_p_sub_lim_ = a_k_p_sub_lim_; tmp_str = 'A'; end;
if nl==1; tmp_d_k_p_sub_ = b_k_p_sub_; tmp_d_k_p_sub_lim_ = a_k_p_sub_lim_; tmp_str = 'B'; end;
if nl==2; tmp_d_k_p_sub_ = c_k_p_sub_; tmp_d_k_p_sub_lim_ = a_k_p_sub_lim_; tmp_str = 'C'; end;
subplot(p_row,p_col,1+np);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('real $F_{%s}(k)$',tmp_str),'Interpreter','latex'); set(gca,'XTick',[],'YTick',[],'ZTick',[])
subplot(p_row,p_col,1+np+p_col);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('imag $F_{%s}(k)$',tmp_str),'Interpreter','latex'); set(gca,'XTick',[],'YTick',[],'ZTick',[])
np=np+1;
subplot(p_row,p_col,1+np);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('real $F_{%s}(k)$ (pole)',tmp_str),'Interpreter','latex'); view(0,90); set(gca,'XTick',[],'YTick',[],'ZTick',[])
subplot(p_row,p_col,1+np+p_col);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('imag $F_{%s}(k)$ (pole)',tmp_str),'Interpreter','latex'); view(0,90); set(gca,'XTick',[],'YTick',[],'ZTick',[])
np=np+1;
end;%for nl=0:2;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.75,512]);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
close(gcf);
%disp('returning');return;
end;%if flag_disp;

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
tmp_parameter = struct('type','parameter');
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),95,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),[ 1,99],'all');
subplot(1,3,1); isosurface_f_x_u_1(tmp_parameter,a_x_c_xxx_); title('$F_{A}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,2); isosurface_f_x_u_1(tmp_parameter,b_x_c_xxx_); title('$F_{B}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,3); isosurface_f_x_u_1(tmp_parameter,c_x_c_xxx_); title('$F_{C}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
close(gcf);

if (flag_verbose> 0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

function ...
[ ...
 k_p_r_all_ ...
,k_n_0_all_ ...
,k_n_1_all_ ...
,k_n_2_all_ ...
,k_n_azimu_b_all_ ...
,k_n_polar_a_all_ ...
,k_n_r01_all_ ...
] = ...
val_helper_(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
k_p_r_all_ = sqrt(k_c_0_all_.^2 + k_c_1_all_.^2 + k_c_2_all_.^2);
k_n_0_all_ = k_c_0_all_./max(1e-12,k_p_r_all_);
k_n_1_all_ = k_c_1_all_./max(1e-12,k_p_r_all_);
k_n_2_all_ = k_c_2_all_./max(1e-12,k_p_r_all_);
k_n_r01_all_ = sqrt(k_n_0_all_.^2 + k_n_1_all_.^2);
k_n_azimu_b_all_ = atan2(k_n_1_all_,k_n_0_all_);
k_n_polar_a_all_ = atan2(k_n_r01_all_,k_n_2_all_);

function val_A_all_ = shell_val_A(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
[k_p_r_all_,k_n_0_all_,k_n_1_all_,k_n_2_all_,k_n_azimu_b_all_,k_n_polar_a_all_,k_n_r01_all_] = val_helper_(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
%%%%;
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
val_A_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'s2'));
dp_n_ = [0;0;1];
dp_n_c_all_ = k_n_0_all_*dp_n_(1+0)+k_n_1_all_*dp_n_(1+1)+k_n_2_all_*dp_n_(1+2);
ds_n_ = [0;1;0];
ds_n_c_all_ = k_n_0_all_*ds_n_(1+0)+k_n_1_all_*ds_n_(1+1)+k_n_2_all_*ds_n_(1+2);
val_A_all_ = ...
+1.0*(abs(dp_n_c_all_)<=circle_r) .* ( ...
+(real(backgr_v)+i*imag(backgr_v)).*( (k_n_2_all_> 0) ) ...
+(real(backgr_v)-i*imag(backgr_v)).*( (k_n_2_all_<=0) ) ...
) ...
+1.0*(abs(dp_n_c_all_)> circle_r) .* ( ...
+(real(square_v)+i*imag(square_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_c_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s2'));
%%%%%%%%;

function val_B_all_ = shell_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
[k_p_r_all_,k_n_0_all_,k_n_1_all_,k_n_2_all_,k_n_azimu_b_all_,k_n_polar_a_all_,k_n_r01_all_] = val_helper_(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
%%%%;
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
circle_c_0 = parameter.circle_c_0;
circle_c_1 = parameter.circle_c_1;
val_B_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'s2'));
dp_n_ = [0;0;1];
dp_n_c_all_ = k_n_0_all_*dp_n_(1+0)+k_n_1_all_*dp_n_(1+1)+k_n_2_all_*dp_n_(1+2);
k_n_rx_all_ = sqrt((k_n_0_all_/circle_c_0).^2 + (k_n_1_all_/circle_c_1).^2 + k_n_2_all_.^2);
dp_n_cx_all_ = ( (k_n_0_all_/circle_c_0)*dp_n_(1+0) + (k_n_1_all_/circle_c_1)*dp_n_(1+1) + k_n_2_all_*dp_n_(1+2) )./max(1e-12,k_n_rx_all_);
ds_n_ = [0;1;0];
ds_n_c_all_ = k_n_0_all_*ds_n_(1+0)+k_n_1_all_*ds_n_(1+1)+k_n_2_all_*ds_n_(1+2);
val_B_all_ = ...
+1.0*(abs(dp_n_cx_all_)<=circle_r) .* ( ...
+(real(backgr_v)+i*imag(backgr_v)).*( (k_n_2_all_> 0) ) ...
+(real(backgr_v)-i*imag(backgr_v)).*( (k_n_2_all_<=0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(square_v)+i*imag(square_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s2'));
%%%%%%%%;

function val_C_all_ = shell_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
[k_p_r_all_,k_n_0_all_,k_n_1_all_,k_n_2_all_,k_n_azimu_b_all_,k_n_polar_a_all_,k_n_r01_all_] = val_helper_(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
%%%%;
backgr_v = parameter.backgr_v;
square_r = parameter.square_r; square_v = parameter.square_v;
circle_r = parameter.circle_r; circle_v = parameter.circle_v;
circle_c_0 = parameter.circle_c_0;
circle_c_1 = parameter.circle_c_1;
val_C_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'s2'));
dp_n_ = [0;0;1];
dp_n_c_all_ = k_n_0_all_*dp_n_(1+0)+k_n_1_all_*dp_n_(1+1)+k_n_2_all_*dp_n_(1+2);
k_n_rx_all_ = sqrt((k_n_0_all_/circle_c_0).^2 + (k_n_1_all_/circle_c_1).^2 + k_n_2_all_.^2);
dp_n_cx_all_ = ( (k_n_0_all_/circle_c_0)*dp_n_(1+0) + (k_n_1_all_/circle_c_1)*dp_n_(1+1) + k_n_2_all_*dp_n_(1+2) )./max(1e-12,k_n_rx_all_);
ds_n_ = [0;1;-0.15]; ds_n_ = ds_n_/fnorm(ds_n_);
ds_n_c_all_ = k_n_0_all_*ds_n_(1+0)+k_n_1_all_*ds_n_(1+1)+k_n_2_all_*ds_n_(1+2);
val_C_all_ = ...
+1.0*(abs(dp_n_cx_all_)<=circle_r) .* ( ...
+(real(backgr_v)+i*imag(backgr_v)).*( (k_n_2_all_> 0) ) ...
+(real(backgr_v)-i*imag(backgr_v)).*( (k_n_2_all_<=0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(square_v)+i*imag(square_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (ds_n_c_all_< 0) & (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (ds_n_c_all_> 0) & (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s2'));
%%%%%%%%;

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


