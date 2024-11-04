%%%%%%%%;
% intended for use with test_slice_vs_volume_integral_trpv1_11.m ;
% This is substantially simpler than test_slice_vs_volume_integral_helper_eig_equa_band_3, ;
% and also simpler than test_slice_vs_volume_integral_helper_eig_equa_band_4. ;
% In this case we never limit ourselves to l_max, and do not project onto principal-modes. ;
%%%%%%%%;

[~,str_hostname] = system('hostname');
flag_128G = 0 ...
| ~isempty(strfind(str_hostname,'xcalibr8')) ...
| ~isempty(strfind(str_hostname,'crunchy')) ...
| ~isempty(strfind(str_hostname,'linserv')) ...
;

str_dir_jpg = sprintf('%s_jpg',dir_ssnll);
if ~exist(str_dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg)); mkdir(str_dir_jpg); end;
str_dir_jpg_stripped = sprintf('%s_jpg_stripped',dir_ssnll);
if ~exist(str_dir_jpg_stripped,'dir'); disp(sprintf(' %% mkdir %s',str_dir_jpg_stripped)); mkdir(str_dir_jpg_stripped); end;
str_infix = 'p_equa_band';

flag_calc = flag_128G;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% Now set up hi_k-quadrature on sphere. ;
%%%%%%%%;
hi_k_int = 1*k_int;
hi_k_eq_d_double = 1*k_eq_d_double;
hi_k_p_r_max = hi_k_int/(2*pi); hi_k_eq_d = hi_k_eq_d_double/(2*pi); str_T_vs_L = 'L';
flag_unif_vs_adap = 1; flag_tensor_vs_adap = 1;
[ ...
 n_hi_k_all ...
,n_hi_k_all_csum_ ...
,hi_k_p_r_all_ ...
,hi_k_p_azimu_b_all_ ...
,hi_k_p_polar_a_all_ ...
,weight_3d_hi_k_all_ ...
,weight_shell_hi_k_ ...
,n_hi_k_p_r ...
,hi_k_p_r_ ...
,weight_3d_hi_k_p_r_ ...
,hi_k_c_0_all_ ...
,hi_k_c_1_all_ ...
,hi_k_c_2_all_ ...
] = ...
sample_sphere_7( ...
 0*flag_verbose ...
,hi_k_p_r_max ...
,hi_k_eq_d ...
,str_T_vs_L ...
,flag_unif_vs_adap ...
,flag_tensor_vs_adap ...
) ; %<-- sum(weight_3d_hi_k_p_r_)*(4*pi) = (4/3)*pi*hi_k_p_r_max^3 --> sum(weight_3d_hi_k_p_r_) = (1/3)*hi_k_p_r_max^3 ;
%%%%%%%%;
nb=0; tmp_b=hi_k_p_azimu_b_all_(1+nb);
while hi_k_p_azimu_b_all_(1+nb+1)>tmp_b; nb=nb+1; tmp_b = hi_k_p_azimu_b_all_(1+nb); end;
n_hi_k_p_azimu_b = 1+nb;
n_hi_k_p_polar_a = n_hi_k_all/n_hi_k_p_azimu_b/n_hi_k_p_r;
hi_k_p_azimu_b_ = hi_k_p_azimu_b_all_(1+[0:n_hi_k_p_azimu_b-1]);
%%%%%%%%;
str_tolerance_pm = sprintf('nltInfpm%d',n_hi_k_p_r);
str_fname_nopath_prefix = sprintf('eig_i1_from_synth_%s_%s',str_tolerance_pm,str_infix);

flag_disp=0;
if flag_disp;
%%%%%%%%;
hi_k_upb = prctile(hi_k_p_r_all_,100);
tmp_index_ = efind(hi_k_p_r_all_==hi_k_upb);
figure(1+nf);nf=nf+1;clf;figsml;
parameter_grid = struct('type','parameter');
parameter_grid.k_max = 0.9975*hi_k_upb;
parameter_grid.k_mid = 0.95*hi_k_upb;
parameter_grid.flag_solid = 1;
plot_sphere_grid_0(parameter_grid);
hold on;
plot3( ...
 hi_k_c_0_all_(1+tmp_index_) ...
,hi_k_c_1_all_(1+tmp_index_) ...
,hi_k_c_2_all_(1+tmp_index_) ...
,'.' ...
);
hold off;
axis equal; axis vis3d;
%%%%%%%%;
end;%if flag_disp;

%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_hi_k_p_quad_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_base_(:).*xxx_u_weight_(:),-1,1e-12,n_hi_k_all,2*pi*hi_k_c_0_all_/eta,2*pi*hi_k_c_1_all_/eta,2*pi*hi_k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_hi_k_p_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/hi_k_p_r_max; tmp_t = tic;
a_x_u_rec0_ = xxnufft3d3(n_hi_k_all,2*pi*hi_k_c_0_all_*eta,2*pi*hi_k_c_1_all_*eta,2*pi*hi_k_c_2_all_*eta,a_hi_k_p_quad_.*(2*pi)^3.*weight_3d_hi_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_rec0_ time %0.2fs',tmp_t));
%%%%;
if (flag_verbose>0); disp(sprintf(' %% a_x_u_reco_ vs a_x_u_rec0_: %0.16f',fnorm(a_x_u_reco_-a_x_u_rec0_)/fnorm(a_x_u_reco_))); end;
%%%%;
dtau = 1e-1;
perturb_hi_k_p_azimu_b_ = hi_k_p_azimu_b_ + dtau*sin(2*hi_k_p_azimu_b_);
[~,a_hi_k_p_pert_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_,a_hi_k_p_quad_,n_hi_k_p_azimu_b,hi_k_p_azimu_b_);
[~,a_hi_k_p_reco_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,hi_k_p_azimu_b_,a_hi_k_p_pert_,n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_);
[~,a_hi_k_p_prec_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_,a_hi_k_p_reco_,n_hi_k_p_azimu_b,hi_k_p_azimu_b_);
disp(sprintf(' %% a_hi_k_p_quad_ vs a_hi_k_p_reco_: %0.16f',fnorm(a_hi_k_p_quad_-a_hi_k_p_reco_)/max(1e-12,fnorm(a_hi_k_p_quad_))));
disp(sprintf(' %% a_hi_k_p_pert_ vs a_hi_k_p_prec_: %0.16f',fnorm(a_hi_k_p_pert_-a_hi_k_p_prec_)/max(1e-12,fnorm(a_hi_k_p_pert_))));
%%%%;
eta = pi/hi_k_p_r_max; tmp_t = tic;
a_x_u_pert_ = xxnufft3d3(n_hi_k_all,2*pi*hi_k_c_0_all_*eta,2*pi*hi_k_c_1_all_*eta,2*pi*hi_k_c_2_all_*eta,a_hi_k_p_pert_.*(2*pi)^3.*weight_3d_hi_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_pert_ time %0.2fs',tmp_t));
%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make frames for perturbation movie. ;
%%%%%%%%;
interp_k = 1;
tmp_dtau_mag = 0.5;
dtau_ = tmp_dtau_mag*[-1:0.125:+1]; n_dtau = numel(dtau_); ddtau = mean(diff(dtau_));
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(a_x_u_reco_(:)),prct);
val_zoom = 2.0;
%%%%%%%%;
for ndtau=0:n_dtau-1;
%%%%%%%%;
dtau_mid = dtau_(1+ndtau);
dtau_pos = dtau_mid + 0.5*ddtau;
dtau_neg = dtau_mid - 0.5*ddtau;
%%%%;
for ntype=0:3-1;
if ntype==0; dtau=dtau_mid; end; if ntype==1; dtau=dtau_pos; end; if ntype==2; dtau=dtau_neg; end;
perturb_hi_k_p_azimu_b_ = hi_k_p_azimu_b_ + dtau*sin(2*hi_k_p_azimu_b_);
[~,a_hi_k_p_pert_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_,a_hi_k_p_quad_,n_hi_k_p_azimu_b,hi_k_p_azimu_b_);
eta = pi/hi_k_p_r_max; tmp_t = tic;
a_x_u_pert_ = xxnufft3d3(n_hi_k_all,2*pi*hi_k_c_0_all_*eta,2*pi*hi_k_c_1_all_*eta,2*pi*hi_k_c_2_all_*eta,a_hi_k_p_pert_.*(2*pi)^3.*weight_3d_hi_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_pert_ time %0.2fs',tmp_t));
if ntype==0; a_x_u_pmid_ = a_x_u_pert_; end; if ntype==1; a_x_u_ppos_ = a_x_u_pert_; end; if ntype==2; a_x_u_pneg_ = a_x_u_pert_; end;
end;%for ntype=0:3-1;
%%%%;
v_x_u_reco_ = a_x_u_ppos_ - a_x_u_pneg_;
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 16;
%%;
subplot(1,2,1);
isosurface_f_x_u_1( ...
 struct('vval_',[vval]) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),interp3(reshape(a_x_u_pmid_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dtau %+.04f',prct,dtau),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%;
subplot(1,2,2);
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),interp3(reshape(a_x_u_pmid_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),interp3(reshape(v_x_u_reco_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('prct %5.2f dtau %+.04f',prct,dtau),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%;
set(gcf,'Position',1+[0,0,1024*2,1.5*(768-128)]);
str_subplot = sprintf('p%.4d_dtau%d',100*prct,ndtau);
fname_fig_pre = sprintf('%s/%s_FIGK_%s',str_dir_jpg,str_fname_nopath_prefix,str_subplot);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGK_%s_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix,str_subplot);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;
end;%for ndtau=0:n_dtau-1;
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
subplot(1,2,[1]);
flim_ = 4.0*[0,2.0/(4*pi)];
flag_2d_vs_3d=0;
imagesc_polar_a_azimu_b_0( ...
 viewing_polar_a_S_ ... 
,viewing_azimu_b_S_ ... 
,viewing_polar_a_S_==pi/2 ...
,flim_ ... 
,colormap_beach ... 
,flag_2d_vs_3d ...
,k_p_r_max ...
);
title('$\mu(\tau)$','Interpreter','latex');
axisnotick3d; axis equal; axis vis3d;
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,2,[2]);
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(a_x_u_reco_(:)),prct);
val_zoom = 2.0;
isosurface_f_x_u_1( ...
 struct('vval_',[vval]) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),a_x_u_reco_) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGL',str_dir_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGL_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
close(gcf);
%%%%%%%%;  
end;%if flag_disp;


%%%%;
n_M_use = n_S;
euler_polar_a_M_use_ = viewing_polar_a_S_;
euler_azimu_b_M_use_ = viewing_azimu_b_S_;
tmp_dtau_euler_polar_a_M_use_ = zeros(n_M_use,1);
tmp_dtau_euler_azimu_b_M_use_ = sin(2*viewing_azimu_b_S_).*exp(-abs(viewing_polar_a_S_-pi/2).^2/max(1e-12,2*pi/128));
tmp_dtau_euler_gamma_z_M_use_ = zeros(n_M_use,1);
weight_imagecount_M_use_ = viewing_weight_S_;
tmp_dtau_euler_fnorm = sum(bsxfun(@times,[tmp_dtau_euler_polar_a_M_use_,tmp_dtau_euler_azimu_b_M_use_,tmp_dtau_euler_gamma_z_M_use_].^2,weight_imagecount_M_use_),'all');
tmp_dtau_euler_nrm_polar_a_M_ = tmp_dtau_euler_polar_a_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_azimu_b_M_ = tmp_dtau_euler_azimu_b_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
tmp_dtau_euler_nrm_gamma_z_M_ = tmp_dtau_euler_gamma_z_M_use_/max(1e-12,tmp_dtau_euler_fnorm);
factor_amplify = 2.0;
%%%%;
ddtau = 1e-2; dtau_mid = 0; dtau_pos = dtau_mid + 0.5*ddtau; dtau_neg = dtau_mid - 0.5*ddtau;
for ntype=0:3-1;
if ntype==0; dtau=dtau_mid; end; if ntype==1; dtau=dtau_pos; end; if ntype==2; dtau=dtau_neg; end;
perturb_hi_k_p_azimu_b_ = hi_k_p_azimu_b_ + dtau*sin(2*hi_k_p_azimu_b_);
[~,a_hi_k_p_pert_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_,a_hi_k_p_quad_,n_hi_k_p_azimu_b,hi_k_p_azimu_b_);
eta = pi/hi_k_p_r_max; tmp_t = tic;
a_x_u_pert_ = xxnufft3d3(n_hi_k_all,2*pi*hi_k_c_0_all_*eta,2*pi*hi_k_c_1_all_*eta,2*pi*hi_k_c_2_all_*eta,a_hi_k_p_pert_.*(2*pi)^3.*weight_3d_hi_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_pert_ time %0.2fs',tmp_t));
if ntype==0; a_x_u_pmid_ = a_x_u_pert_; end; if ntype==1; a_x_u_ppos_ = a_x_u_pert_; end; if ntype==2; a_x_u_pneg_ = a_x_u_pert_; end;
end;%for ntype=0:3-1;
%%%%;
v_x_u_reco_ = a_x_u_ppos_ - a_x_u_pneg_;
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
tmp_index_ = efind(abs(euler_polar_a_M_use_-pi/2)<1e-6);
tmp_n_M_use = numel(tmp_index_);
%%%%;
subplot(1,3,[1,2]);
hold on;
sphere_compass__0( ...
struct('type','sphere_compass','compass_r_base',1.0/(2*pi)/2*sqrt(0.5),'flag_2d_vs_3d',1,'flag_normalize',1) ...
,tmp_n_M_use ...
,euler_polar_a_M_use_(1+tmp_index_) ...
,euler_azimu_b_M_use_(1+tmp_index_) ...
,tmp_dtau_euler_nrm_polar_a_M_(1+tmp_index_) ...
,tmp_dtau_euler_nrm_azimu_b_M_(1+tmp_index_) ...
,tmp_dtau_euler_nrm_gamma_z_M_(1+tmp_index_) ...
);
hold off;
xlim([0,2*pi]); xlabel('azimuthal','Interpreter','latex');
set(gca,'XTick',0:pi/4:2*pi,'XTickLabel',[]);
ylim([0,1*pi]); ylabel('polar','Interpreter','latex');
set(gca,'YTick',0:pi/4:1*pi,'YTickLabel',[]);
grid on;
title('$\Delta\tau$: equatorial','Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,[3]);
prct = 98.5;
isosurface_f_x_u_1(struct('percent_threshold_',[100-prct,prct],'c_use__',colormap_pm),v_x_u_reco_);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$'),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
%%%%;
set(gcf,'Position',1+[0,0,1024*1.75,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGN',str_dir_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGN_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

%%%%;
ddtau = 0.25; dtau_mid = 0; dtau_pos = dtau_mid + 0.5*ddtau; dtau_neg = dtau_mid - 0.5*ddtau;
for ntype=0:3-1;
if ntype==0; dtau=dtau_mid; end; if ntype==1; dtau=dtau_pos; end; if ntype==2; dtau=dtau_neg; end;
perturb_hi_k_p_azimu_b_ = hi_k_p_azimu_b_ + dtau*sin(2*hi_k_p_azimu_b_);
[~,a_hi_k_p_pert_] = interp1_azimu_b_0([],n_hi_k_p_azimu_b,perturb_hi_k_p_azimu_b_,a_hi_k_p_quad_,n_hi_k_p_azimu_b,hi_k_p_azimu_b_);
eta = pi/hi_k_p_r_max; tmp_t = tic;
a_x_u_pert_ = xxnufft3d3(n_hi_k_all,2*pi*hi_k_c_0_all_*eta,2*pi*hi_k_c_1_all_*eta,2*pi*hi_k_c_2_all_*eta,a_hi_k_p_pert_.*(2*pi)^3.*weight_3d_hi_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_pert_ time %0.2fs',tmp_t));
if ntype==0; a_x_u_pmid_ = a_x_u_pert_; end; if ntype==1; a_x_u_ppos_ = a_x_u_pert_; end; if ntype==2; a_x_u_pneg_ = a_x_u_pert_; end;
end;%for ntype=0:3-1;
%%%%;
v_x_u_reco_ = 0.5*(a_x_u_ppos_ - a_x_u_pneg_);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
n_dvol = 2; dvol_mag = fnorm(v_x_u_reco_)/fnorm(a_x_u_pmid_);
%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% Make splash figure (i.e., summary). ;
%%%%%%%%;
interp_k = 1;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 16;
%%%%;
prct_ = [98.75]; prct = prct_(1+0);
vval = prctile(real(a_x_u_reco_(:)),prct);
vlim_g_ = prctile(real(v_x_u_reco_(:)),[0.25,99.75]);
val_zoom = 2.0;
for ndvol=0:n_dvol-1;
if ndvol==0; subplot(1,2,[1]); tmp_a_x_u_reco_ = a_x_u_pneg_; end;
if ndvol==1; subplot(1,2,[2]); tmp_a_x_u_reco_ = a_x_u_ppos_; end;
isosurface_f_x_u_2( ...
 struct('vval_',[vval],'vlim_g_',vlim_g_) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),interp3(reshape(tmp_a_x_u_reco_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
,a_x_u_from_a_x_u_zoom_0(struct('val_zoom',val_zoom),interp3(reshape(v_x_u_reco_,[n_x_u_pack,n_x_u_pack,n_x_u_pack]),interp_k,'spline')) ...
);
xlabel(''); ylabel(''); zlabel('');
title(sprintf('$\\Delta \\hat{F}$: $%+0.3f$',(-1)^(1+ndvol)*dvol_mag),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
end;%for ndvol=0:n_dvol-1;
%%%%;
set(gcf,'Position',1+[0,0,1024*1.0,512]);
%%%%;
fname_fig_pre = sprintf('%s/%s_FIGO',str_dir_jpg,str_fname_nopath_prefix);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_eps = sprintf('%s.eps',fname_fig_pre);
fname_fig_stripped_pre = sprintf('%s/%s_FIGO_stripped',str_dir_jpg_stripped,str_fname_nopath_prefix);
fname_fig_stripped_jpg = sprintf('%s.jpg',fname_fig_stripped_pre);
fname_fig_stripped_eps = sprintf('%s.eps',fname_fig_stripped_pre);
if flag_replot | ~exist(fname_fig_jpg,'file');
disp(sprintf(' %% writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
disp(sprintf(' %% writing %s',fname_fig_stripped_pre));
sgtitle('');
print('-djpeg',fname_fig_stripped_jpg);
end;%if flag_replot | ~exist(fname_fig_jpg,'file');
%%%%%%%%;  
close(gcf);
end;%if flag_disp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%if flag_calc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



