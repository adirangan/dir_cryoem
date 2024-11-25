function test_heterogeneity_spurious_3();
%%%%%%%%;
% Sets up idealized single shell examples. ;
% Beach-ball example for niko and kexin. ;
% Equatorial projections are interchangable, ;
% but polar projection is unique. ;
%%%%%%%%;

str_thisfunction = 'test_heterogeneity_spurious_3';

flag_verbose=1;
flag_disp=1; nf=0;
%%%%%%%%;
half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_c = 64*2;
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
k_p_r_max = (0.5)*48/(2*pi); k_eq_d = sqrt(0.5)*1.0/(2*pi);
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
%parameter.str_shape = 's1';
parameter.sigma_g = 0.125/1.0;
parameter.str_sigma_g = sprintf('sg%d',round(parameter.sigma_g*1000));
parameter.mu_g = 0.650;
parameter.str_mu_g = sprintf('sm%d',round(parameter.mu_g*1000));
parameter.backgr_v = +0.0 + i*0.0;
parameter.str_backgr_v = sprintf('bv%di%d',8+round(real(parameter.backgr_v)*8),8+round(imag(parameter.backgr_v)*8));
parameter.square_r = 0.75*1/sqrt(2);
parameter.str_square_r = sprintf('sr%d',round(parameter.square_r*sqrt(2)*8));
parameter.square_v = +1.0 + i*0.500;
parameter.str_square_v = sprintf('sv%di%d',8+round(real(parameter.square_v)*8),8+round(imag(parameter.square_v)*8));
parameter.circle_r = 0.65;
parameter.str_circle_r = sprintf('cr%d',round(parameter.circle_r*8));
parameter.circle_v = +0.5 + i*0.250;
parameter.str_circle_v = sprintf('cv%di%d',8+round(real(parameter.circle_v)*8),8+round(imag(parameter.circle_v)*8));
parameter.circle_c_0 = 0.5;
parameter.circle_c_1 = 1.0;
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
%prctile_use = 95;
prctile_use = 90;

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
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('real $F_{%s}(k)$',tmp_str),'Interpreter','latex'); set(gca,'XTick',[],'YTick',[],'ZTick',[]); view(-60,15); 
subplot(p_row,p_col,1+np+p_col);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('imag $F_{%s}(k)$',tmp_str),'Interpreter','latex'); set(gca,'XTick',[],'YTick',[],'ZTick',[]); view(-60,15); 
np=np+1;
subplot(p_row,p_col,1+np);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('real $F_{%s}(k)$ (pole)',tmp_str),'Interpreter','latex'); view(0,90); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(p_row,p_col,1+np+p_col);
imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(tmp_d_k_p_sub_),tmp_d_k_p_sub_lim_,colormap_beach(),0); title(sprintf('imag $F_{%s}(k)$ (pole)',tmp_str),'Interpreter','latex'); view(0,90); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
np=np+1;
end;%for nl=0:2;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.75,512]);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%close(gcf);
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
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),prctile_use,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),[ 1,99],'all');
subplot(1,3,1); isosurface_f_x_u_1(tmp_parameter,a_x_c_xxx_); title('$F_{A}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,2); isosurface_f_x_u_1(tmp_parameter,b_x_c_xxx_); title('$F_{B}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,3); isosurface_f_x_u_1(tmp_parameter,c_x_c_xxx_); title('$F_{C}(x)$','Interpreter','latex');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%close(gcf);

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
if (strcmp(parameter.str_shape,'s1'));
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
+(real(square_v)+i*imag(square_v)).*( (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_c_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s2'));
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
circle_c_1 = 1.0*parameter.circle_c_1;
val_B_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'s1'));
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
+(real(square_v)+i*imag(square_v)).*( (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s1'));
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
circle_c_0 = 0.5*parameter.circle_c_0;
circle_c_1 = 1.0*parameter.circle_c_1;
val_C_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'s1'));
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
+(real(square_v)+i*imag(square_v)).*( (k_n_2_all_> 0) ) ...
+(real(square_v)-i*imag(square_v)).*( (k_n_2_all_< 0) ) ...
) ...
+1.0*(abs(dp_n_cx_all_)> circle_r) .* ( ...
+(real(circle_v)+i*imag(circle_v)).*( (k_n_2_all_> 0) ) ...
+(real(circle_v)-i*imag(circle_v)).*( (k_n_2_all_< 0) ) ...
) ...
;
return;
end;%if (strcmp(parameter.str_shape,'s1'));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%{
parameter = struct('type','parameter');
parameter.str_shape = 'beachball';
parameter.sigma_g = 0.125/2;
parameter.str_sigma_g = sprintf('sg%d',round(parameter.sigma_g*1000));
parameter.mu_g = 0.450/1;
parameter.str_mu_g = sprintf('sm%d',round(parameter.mu_g*1000));
parameter.infix = sprintf( ...
 '%s%s%s' ...
,parameter.str_shape ...
,parameter.str_sigma_g ...
,parameter.str_mu_g ...
);

flag_disp=1;
if flag_disp;
%%%%%%%%;
% illustrate single shell. ;
%%%%%%%%;
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
a_k_p_sub_lim_ = prctile(abs(a_k_p_sub_),95,'all')*[-1,+1];
subplot(2,3,1); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(a_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('real(F(k)) ori','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,2); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(b_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('real(F(k)) rot+','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,3); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(c_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('real(F(k)) rot-','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,4); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(a_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('imag(F(k)) ori','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,5); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(b_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('imag(F(k)) rot+','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,6); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(c_k_p_sub_),a_k_p_sub_lim_,colormap_beach(),0); title('imag(F(k)) rot-','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGK',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
disp('returning'); return;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% illustrate equatorial templates. ;
%%%%%%%%;
n_w_max = ceil(2*(2*pi*k_p_r_max));
n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,E_k_p_r_wk_ ...
,E_k_p_w_wk_ ...
,E_k_c_0_wk_ ...
,E_k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
);
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_(:)]);
E_k_c_2_wk_ = zeros(n_w_sum,1);
%%%%;
E_k_p_wk_ = sphere_val_A(parameter,k_p_r_max,n_w_sum,E_k_c_0_wk_,E_k_c_1_wk_,E_k_c_2_wk_);
E_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,E_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum ;
E_k_lim_ = prctile(abs(E_k_p_wk_),95,'all')*[-1,+1];
E_x_lim_ = prctile(abs(E_x_c_xx_),95,'all')*[-1,+1];
figure(1+nf);nf=nf+1;clf;figmed;p_row=1;p_col=3;np=0;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(E_k_p_wk_),E_k_lim_,colormap_80s);
title(sprintf('real(E_k_p_wk_) equi'),'Interpreter','none'); axis image; axisnotick;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(E_k_p_wk_),E_k_lim_,colormap_80s);
title(sprintf('imag(E_k_p_wk_) equi'),'Interpreter','none'); axis image; axisnotick;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(E_x_c_xx_),E_x_lim_,colormap_beach);
title(sprintf('real(E_x_c_xx_) equi'),'Interpreter','none'); axis image; axisnotick;
%subplot(p_row,p_col,1+np);np=np+1;
%imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,imag(E_x_c_xx_),E_x_lim_,colormap_beach);
%title(sprintf('imag(E_x_c_xx_) equi'),'Interpreter','none'); axis image; axisnotick;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGE',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024*1.5]); p_row=4;p_col=3;np=0;
for prow=0:p_row-1;
azimu_b = (pi*(prow+0.5))/max(1,p_row);
S_k_c_0_wk_ = E_k_c_0_wk_; S_k_c_1_wk_ = E_k_c_1_wk_; S_k_c_2_wk_ = E_k_c_2_wk_;
[S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_] = Rz(S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_,+0.5*pi/n_w_max);
[S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_] = Ry(S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_,+pi/2);
[S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_] = Rz(S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_,+azimu_b);
S_k_p_wk_ = sphere_val_A(parameter,k_p_r_max,n_w_sum,S_k_c_0_wk_,S_k_c_1_wk_,S_k_c_2_wk_);
S_x_c_xx_ = interp_k_p_to_x_c_xxnufft(n_x_c,diameter_x_c,n_x_c,diameter_x_c,n_k_p_r,k_p_r_,n_w_,S_k_p_wk_.*weight_2d_wk_*(2*pi)^2)*sqrt(n_x_c*n_x_c) * n_w_sum ;
S_k_lim_ = prctile(abs(S_k_p_wk_),95,'all')*[-1,+1];
S_x_lim_ = prctile(abs(S_x_c_xx_),95,'all')*[-1,+1];
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wk_),S_k_lim_,colormap_80s);
title(sprintf('real(S_k_p_wk_) azimu_b %3d',round(180*azimu_b/pi)),'Interpreter','none'); axis image; axisnotick;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,imag(S_k_p_wk_),S_k_lim_,colormap_80s);
title(sprintf('imag(S_k_p_wk_) azimu_b %3d',round(180*azimu_b/pi)),'Interpreter','none'); axis image; axisnotick;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,real(S_x_c_xx_),S_x_lim_,colormap_beach);
title(sprintf('real(S_x_c_xx_) azimu_b %3d',round(180*azimu_b/pi)),'Interpreter','none'); axis image; axisnotick;
%subplot(p_row,p_col,1+np);np=np+1;
%imagesc_c(n_x_c,x_c_0_,n_x_c,x_c_1_,imag(S_x_c_xx_),S_x_lim_,colormap_beach);
%title(sprintf('imag(S_x_c_xx_) azimu_b %3d',round(180*azimu_b/pi)),'Interpreter','none'); axis image; axisnotick;
end;% for prow=0:p_row-1;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGS',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%disp('returning'); return;
end;%if flag_disp;

tmp_vval = 92.5;
tmp_vlim_ = [0.5,99.5];
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
%%%%;
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGA',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;figsml;
fontsize_use = 24;
tmp_parameter = struct('type','parameter');
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),tmp_vval,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),tmp_vlim_,'all');
subplot(1,1,1); isosurface_f_x_u_1(tmp_parameter,a_x_c_xxx_); title('F(x)');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
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
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGABC',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
figure(1+nf);nf=nf+1;figmed;
fontsize_use = 24;
tmp_parameter = struct('type','parameter');
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),tmp_vval,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),tmp_vlim_,'all');
subplot(1,3,1); isosurface_f_x_u_1(tmp_parameter,a_x_c_xxx_); title('F(x) ori');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,2); isosurface_f_x_u_1(tmp_parameter,b_x_c_xxx_); title('F(x) rot+');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,3); isosurface_f_x_u_1(tmp_parameter,c_x_c_xxx_); title('F(x) rot-');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
%close(gcf);

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
val_A_all_ = zeros(n_k_all,1);
%%%%%%%%;
if ( ~isfield(parameter,'str_shape') );
return;
end;%if ( ~isfield(parameter,'str_shape') );
%%%%;
if (strcmp(parameter.str_shape,'beachball'));
val_symm_all_ = ...
1.5 + ...
+1.000*sin( 2.0*k_n_azimu_b_all_) ...
+0.500*sin( 4.0*k_n_azimu_b_all_) ...
+0.250*sin( 8.0*k_n_azimu_b_all_) ...
+0.125*sin(16.0*k_n_azimu_b_all_) ...
;
val_asym_all_ = ...
+1.000*cos(1.0*k_n_polar_a_all_) ...
+0.333*cos(3.0*k_n_polar_a_all_) ...
;
val_A_all_ = val_symm_all_ + i.*val_asym_all_.*val_symm_all_;
return;
end;%if (strcmp(parameter.str_shape,'beachball'));
%%%%%%%%;

function val_B_all_ = shell_val_B(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
gamma_z = +pi/6;
[k_c_0_rot_,k_c_1_rot_,k_c_2_rot_] = Rz(k_c_0_all_,k_c_1_all_,k_c_2_all_,gamma_z);
val_B_all_ = shell_val_A(parameter,k_p_r_max,n_k_all,k_c_0_rot_,k_c_1_rot_,k_c_2_rot_);

function val_C_all_ = shell_val_C(parameter,k_p_r_max,n_k_all,k_c_0_all_,k_c_1_all_,k_c_2_all_);
gamma_z = -pi/6;
[k_c_0_rot_,k_c_1_rot_,k_c_2_rot_] = Rz(k_c_0_all_,k_c_1_all_,k_c_2_all_,gamma_z);
val_C_all_ = shell_val_A(parameter,k_p_r_max,n_k_all,k_c_0_rot_,k_c_1_rot_,k_c_2_rot_);

function [k_c_0_rot_,k_c_1_rot_,k_c_2_rot_] = Rz(k_c_0_all_,k_c_1_all_,k_c_2_all_,gamma_z);
cz = cos(gamma_z); sz = sin(gamma_z);
Rz__ = [ ...
 +cz -sz 0 ...
;+sz +cz 0 ...
; 0   0  1 ...
];
tmp__ = [k_c_0_all_(:),k_c_1_all_(:),k_c_2_all_(:)]*transpose(Rz__);
k_c_0_rot_ = tmp__(:,1+0);
k_c_1_rot_ = tmp__(:,1+1);
k_c_2_rot_ = tmp__(:,1+2);

function [k_c_0_rot_,k_c_1_rot_,k_c_2_rot_] = Ry(k_c_0_all_,k_c_1_all_,k_c_2_all_,beta_a);
ca = cos(beta_a); sa = sin(beta_a);
Ry__ = [ ...
 +ca 0 +sa ...
; 0  1  0  ...
;-sa 0 +ca ...
];
tmp__ = [k_c_0_all_(:),k_c_1_all_(:),k_c_2_all_(:)]*transpose(Ry__);
k_c_0_rot_ = tmp__(:,1+0);
k_c_1_rot_ = tmp__(:,1+1);
k_c_2_rot_ = tmp__(:,1+2);

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

%}


