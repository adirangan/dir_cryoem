function test_heterogeneity_spurious_2();
%%%%%%%%;
% Sets up idealized single shell examples. ;
% Beach-ball example for niko and kexin. ;
% Equatorial projections are interchangable, ;
% but polar projection is unique. ;
%%%%%%%%;

str_thisfunction = 'test_heterogeneity_spurious_2';

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
k_p_r_max = (2.0)*48/(2*pi); k_eq_d = sqrt(1.0)*1.0/(2*pi);
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

if (flag_verbose> 0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;

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
subplot(2,3,1); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(a_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('real(F(k)) ori','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,2); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(b_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('real(F(k)) rot','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,3); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,real(c_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('real(F(k)) rot','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,4); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(a_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('imag(F(k)) ori','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,5); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(b_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('imag(F(k)) rot','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
subplot(2,3,6); imagesc_polar_a_azimu_b_0(k_p_polar_a_sub_,k_p_azimu_b_sub_,imag(c_k_p_sub_),a_k_p_sub_lim_,colormap_80s(),0); title('imag(F(k)) rot','Interpreter','none'); set(gca,'XTick',[],'YTick',[],'ZTick',[]);
dir_jpg = sprintf('%s/dir_jpg',pwd); if ~exist(dir_jpg,'dir'); mkdir(dir_jpg); end;
fname_fig_pre = sprintf('%s/test_heterogeneity_spurious_%s_FIGK',dir_jpg,parameter.infix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% writing %s',fname_fig_jpg));
print('-djpeg',fname_fig_jpg);
disp('returning'); return;
end;%if flag_disp;

flag_disp=1;
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
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),95,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),[ 1,99],'all');
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
tmp_parameter.vval_ = prctile(abs(a_x_c_xxx_),95,'all');
tmp_parameter.vlim_ = prctile(abs(a_x_c_xxx_),[ 1,99],'all');
subplot(1,3,1); isosurface_f_x_u_1(tmp_parameter,a_x_c_xxx_); title('F(x) ori');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,2); isosurface_f_x_u_1(tmp_parameter,b_x_c_xxx_); title('F(x) rot');
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(1,3,3); isosurface_f_x_u_1(tmp_parameter,c_x_c_xxx_); title('F(x) rot');
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


