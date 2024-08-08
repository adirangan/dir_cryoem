%%%%%%%%;
% applying eig_ddssnll_lanczos_1 to a single shell. ;
%%%%%%%%;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

str_thisfunction = 'test_slice_vs_volume_integral_shell_10';
flag_recalc=0; flag_replot=1;
flag_verbose=1; flag_disp=1; nf=0;

k_int = 16;
if k_int==16;
k_eq_d_double = 0.50;
t_eq_d_double = 0.50;
n_w_int = 2;
KAPPA_flag_kernel_full = 1;
KAPPA_pole_north_double = 12*pi/24;
KAPPA_pole_south_double = 12*pi/24;
KAPPA_qref_k_eq_d_double = 0.5;
lanczos_n_iteration_max = 128;
end;%if k_int==16;
if k_int==48;
k_eq_d_double = 1.00;
t_eq_d_double = 1.00;
n_w_int = 1;
KAPPA_flag_kernel_full = 1;
KAPPA_pole_north_double = 12*pi/24;
KAPPA_pole_south_double = 12*pi/24;
KAPPA_qref_k_eq_d_double = 1.00;
lanczos_n_iteration_max = 128;
end;%if k_int==48;

%%%%%%%%;
fname_prefix = 'shell';
fname_prefix_xfix = sprintf('%s',fname_prefix);
dir_ssnll = sprintf('/%s/rangan/dir_cryoem/dir_%s/dir_ssnll_k%d',string_root,fname_prefix_xfix,k_int);
if (~exist(sprintf('%s_mat',dir_ssnll),'dir')); disp(sprintf(' %% mkdir %s_mat',dir_ssnll)); mkdir(sprintf('%s_mat',dir_ssnll)); end;
if (~exist(sprintf('%s_jpg',dir_ssnll),'dir')); disp(sprintf(' %% mkdir %s_jpg',dir_ssnll)); mkdir(sprintf('%s_jpg',dir_ssnll)); end;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

half_diameter_x_c = 1.0d0;
diameter_x_c = 2.0d0*half_diameter_x_c;
x_p_r_max = 1.0;
n_x_u = 256;
n_x_u_pack = 64;
n_pack = n_x_u/n_x_u_pack;
x_u_0_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_1_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
x_u_2_ = linspace(-x_p_r_max,+x_p_r_max,n_x_u_pack);
[x_u_0___,x_u_1___,x_u_2___] = ndgrid(x_u_0_,x_u_1_,x_u_2_); n_xxx_u = n_x_u_pack^3;
xxx_u_weight_ = (2*x_p_r_max/n_x_u_pack)^3;

%%%%%%%%;
% Now set up k-quadrature on shell. ;
%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi); str_T_vs_L = 'C';
flag_tensor_vs_adap = 0; %<-- This is set to match test_ssnll_from_a_k_Y_12 ;
[ ...
 n_k_all ...
,k_p_azimu_b_all_ ...
,k_p_polar_a_all_ ...
,weight_shell_k_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,n_polar_a ...
,polar_a_ ...
,n_azimu_b_ ...
] = ...
sample_shell_6( ...
 k_p_r_max ...
,k_eq_d ...
,str_T_vs_L ...
,flag_tensor_vs_adap ...
) ;
n_k_all_csum_ = [0;n_k_all];
weight_3d_k_all_ = weight_shell_k_/sum(weight_shell_k_) * (1/3)*k_p_r_max^3 * (4*pi);
weight_3d_k_p_r_ = (1/3)*k_p_r_max^3; %<-- sum(weight_3d_k_p_r_)*(4*pi) = (4/3)*pi*k_p_r_max^3 --> sum(weight_3d_k_p_r_) = (1/3)*k_p_r_max^3 ;
nk_p_r = 0; n_k_p_r = 1; k_p_r_ = k_p_r_max;
k_p_r_all_ = k_p_r_max*ones(n_k_all,1);
%%%%;
if (flag_disp>1);
figure(1+nf);nf=nf+1;clf;figsml;
subplot(p_row,p_col,1+nplot);
plot3(k_c_0_all_,k_c_1_all_,k_c_2_all_,'.');
axis equal; axis vis3d; axisnotick3d;
title(sprintf('nk_p_r %d/%d',nk_p_r,n_k_p_r),'Interpreter','none');
close(gcf);
end;%if flag_disp;
%%%%;

%%%%%%%%;
% Now set up polar-quadrature-weights on disk. ;
%%%%%%%%;
l_max_upb = round(2*pi*k_p_r_max);
l_max_max = min(l_max_upb,1+ceil(2*pi*k_p_r_(end)));
n_w_max = n_w_int*2*(l_max_max+1); n_w_0in_ = n_w_max*ones(n_k_p_r,1);
[ ...
 n_w_ ...
,weight_2d_k_p_r_ ...
,weight_2d_wk_ ...
,k_p_r_wk_ ...
,k_p_w_wk_ ...
,k_c_0_wk_ ...
,k_c_1_wk_ ...
] = ...
get_weight_2d_2( ...
 0*flag_verbose ...
,n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,-1 ...
,n_w_0in_ ...
); %<-- sum(weight_2d_k_p_r_) = pi*k_p_r_max^2/(4*pi^2) ;
n_w_sum = sum(n_w_); n_w_csum_ = cumsum([0;n_w_]);
%%%%%%%%;
% Now set up spherical-harmonics. ;
%%%%%%%%;
l_max_ = zeros(n_k_p_r,1);
for nk_p_r=0:n_k_p_r-1;
l_max_(1+nk_p_r) = max(0,min(l_max_upb,1+ceil(2*pi*k_p_r_(1+nk_p_r))));
end;%for nk_p_r=0:n_k_p_r-1;
n_lm_ = (l_max_+1).^2;
n_lm_max = max(n_lm_);
n_lm_sum = sum(n_lm_);
n_lm_csum_ = cumsum([0;n_lm_]);
l_max_max = max(l_max_);
m_max_ = -l_max_max : +l_max_max;
n_m_max = length(m_max_);
%%%%;
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

%%%%%%%%;
% Now define a_k_p_form_ ;
%%%%%%%%;
sigma_source = k_p_r_max/2.0;
n_source = 24;
rng(0);
n_3 = 3;
delta_source_3s__ = randn(n_3,n_source);
for nsource=0:n_source-1;
delta_source_3_ = delta_source_3s__(:,1+nsource);
delta_source_3_ = delta_source_3_/max(1e-12,fnorm(delta_source_3_));
delta_source_3s__(:,1+nsource) = k_p_r_max*delta_source_3_;
end;%for nsource=0:n_source-1;
a_k_p_form_ = zeros(n_k_all,1);
for nsource=0:n_source-1;
delta_source_3_ = delta_source_3s__(:,1+nsource);
tmp_p2_0_ = (k_c_0_all_ - delta_source_3_(1+0)).^2;
tmp_p2_1_ = (k_c_1_all_ - delta_source_3_(1+1)).^2;
tmp_p2_2_ = (k_c_2_all_ - delta_source_3_(1+2)).^2;
tmp_p2_ = tmp_p2_0_ + tmp_p2_1_ + tmp_p2_2_ ;
epos_all_ = exp(-tmp_p2_./(2*sigma_source^2));
tmp_n2_0_ = (k_c_0_all_ + delta_source_3_(1+0)).^2;
tmp_n2_1_ = (k_c_1_all_ + delta_source_3_(1+1)).^2;
tmp_n2_2_ = (k_c_2_all_ + delta_source_3_(1+2)).^2;
tmp_n2_ = tmp_n2_0_ + tmp_n2_1_ + tmp_n2_2_ ;
eneg_all_ = exp(-tmp_n2_./(2*sigma_source^2));
a_k_p_form_ = a_k_p_form_ + epos_all_ + i*epos_all_ + eneg_all_ - i*eneg_all_ ;
end;%for nsource=0:n_source-1;
if flag_disp;
figure(1+nf);nf=nf+1;clf;figmed;
alim_ = prctile(abs(a_k_p_form_),[  0,100],'all');
alim_ = mean(alim_) + 1.25*0.5*diff(alim_)*[-1,+1];
flag_2d_vs_3d = 0;
subplot(1,2,1);
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_ ... 
,k_p_azimu_b_all_ ... 
,real(a_k_p_form_) ... 
,alim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,1+0*k_p_r_max ...
);
if flag_2d_vs_3d==0; axisnotick3d; axis vis3d; end;
if flag_2d_vs_3d==1;
xlim([0,2*pi]); ylim([0,1*pi]);
xlabel('azimu_b','Interpreter','none');
ylabel('polar_a','Interpreter','none');
axisnotick;
end;%if flag_2d_vs_3d==1;
title('real(a_k_p_form_)','Interpreter','none');
subplot(1,2,2);
imagesc_polar_a_azimu_b_0( ...
 k_p_polar_a_all_ ... 
,k_p_azimu_b_all_ ... 
,imag(a_k_p_form_) ... 
,alim_ ... 
,colormap_80s ... 
,flag_2d_vs_3d ...
,1+0*k_p_r_max ...
);
if flag_2d_vs_3d==0; axisnotick3d; axis vis3d; end;
if flag_2d_vs_3d==1;
xlim([0,2*pi]); ylim([0,1*pi]);
xlabel('azimu_b','Interpreter','none');
ylabel('polar_a','Interpreter','none');
axisnotick;
end;%if flag_2d_vs_3d==1;
title('imag(a_k_p_form_)','Interpreter','none');
close(gcf);
end;%if flag_disp;
%%%%;
eta = pi/k_p_r_max; tmp_t = tic;
a_x_u_quad_ = xxnufft3d3(n_k_all,2*pi*k_c_0_all_*eta,2*pi*k_c_1_all_*eta,2*pi*k_c_2_all_*eta,a_k_p_form_.*(2*pi)^3.*weight_3d_k_all_,+1,1e-12,n_xxx_u,x_u_0___(:)/eta,x_u_1___(:)/eta,x_u_2___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_x_u_quad_ time %0.2fs',tmp_t));
%%%%;
eta = pi/x_p_r_max; tmp_t = tic;
a_k_p_reco_ = xxnufft3d3(n_xxx_u,x_u_0___(:)*eta,x_u_1___(:)*eta,x_u_2___(:)*eta,a_x_u_quad_(:).*xxx_u_weight_(:),-1,1e-12,n_k_all,2*pi*k_c_0_all_/eta,2*pi*k_c_1_all_/eta,2*pi*k_c_2_all_/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_p_reco_ time %0.2fs',tmp_t));
%%%%;
a_k_p_quad_ = a_k_p_form_;

%%%%%%%%;
% Now convert to a_k_Y_ ; 
%%%%%%%%;
tmp_t = tic;
if ~exist('Ylm_uklma___','var'); Ylm_uklma___ = []; end;
if ~exist('k_p_azimu_b_sub_uka__','var'); k_p_azimu_b_sub_uka__ = []; end;
if ~exist('k_p_polar_a_sub_uka__','var'); k_p_polar_a_sub_uka__ = []; end;
if ~exist('l_max_uk_','var'); l_max_uk_ = []; end;
if ~exist('index_nu_n_k_per_shell_from_nk_p_r_','var'); index_nu_n_k_per_shell_from_nk_p_r_ = []; end;
if ~exist('index_k_per_shell_uka__','var'); index_k_per_shell_uka__ = []; end;
[ ...
 a_k_Y_quad_yk_ ...
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
,a_k_p_quad_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_quad_yk_ time %0.2fs',tmp_t));
tmp_t = tic;
[ ...
 a_k_p_reco_ ...
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
,a_k_Y_quad_yk_ ...
,Ylm_uklma___ ...
,k_p_azimu_b_sub_uka__ ...
,k_p_polar_a_sub_uka__ ...
,l_max_uk_ ...
,index_nu_n_k_per_shell_from_nk_p_r_ ...
,index_k_per_shell_uka__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_p_reco_ time %0.2fs',tmp_t));
disp(sprintf(' %% convert_spharm_to_k_p_4: a_k_p_reco error: %0.16f',fnorm(a_k_p_quad_-a_k_p_reco_)/fnorm(a_k_p_quad_))); %<-- this should be 2-3 digits. ;
%%%%%%%%;
a_k_Y_quad_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_quad_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_quad_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_Y_quad_A',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
subplot(1,3,1); plot(Y_l_val_,log10(abs(a_k_Y_quad_yk_)),'.'); xlabel('Y_l_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,2); plot(Y_m_val_,log10(abs(a_k_Y_quad_yk_)),'.'); xlabel('Y_m_val_','Interpreter','none'); ylabel('log10(abs(a))');
subplot(1,3,3); plot(Y_k_val_,log10(abs(a_k_Y_quad_yk_)),'.'); xlabel('Y_k_val_','Interpreter','none'); ylabel('log10(abs(a))');
figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/a_k_Y_quad_yk_',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;
clf;
imagesc_Y(k_p_r_max,n_k_p_r,k_p_r_,l_max_,log10(abs(a_k_Y_quad_yk_)),[-10,0],colormap_beach());
xlabel('l_val','Interpreter','none');
ylabel('m_val','Interpreter','none');
zlabel('k_p_r','Interpreter','none');
axis vis3d; view([-65,+20]);
title('a_k_Y_quad_yk_','Interpreter','none');
%figbig;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% generate templates S_k_p_wk_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
tmp_t = tic();
template_k_eq_d = t_eq_d_double/k_p_r_max;
flag_tensor_vs_adap = 1; %<-- tensor grid. ;
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,viewing_k_c_0_S_ ...
,viewing_k_c_1_S_ ...
,viewing_k_c_2_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,template_k_eq_d ...
,str_T_vs_L ...
,flag_tensor_vs_adap ...
) ;
n_S = n_viewing_S;
if (flag_verbose>0); disp(sprintf(' %% n_S %d, n_viewing_polar_a %d, n_viewing_azimu_b [%d,..,%d]',n_S,n_viewing_polar_a,n_viewing_azimu_b_(1+0),n_viewing_azimu_b_(end))); end;
tmp_t = tic();
[ ...
 S_k_p_wkS__ ...
,n_w_max ...
,n_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_quad_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,n_viewing_polar_a ...
,viewing_polar_a_ ...
,n_viewing_azimu_b_ ...
);
S_k_p_wkS__ = reshape(S_k_p_wkS__,[n_w_sum,n_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% S_k_p_wkS__ (pm_template_2): %0.6fs',tmp_t)); end;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_p_wkS__',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = prctile(real(S_k_p_wkS__),[  0,100],'all');
Slim_ = mean(Slim_) + 1.25*0.5*diff(Slim_)*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(S_k_p_wkS__(:,1+nS)),Slim_,colormap_80s);
axis equal; axisnotick; title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('Sample S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% Now convert templates to S_k_q_wkS__. ;
%%%%%%%%;
S_k_q_wkS__ = zeros(n_w_sum,n_S);
for nS=0:n_S-1;
S_k_q_wkS__(:,1+nS) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,S_k_p_wkS__(:,1+nS));
end;%for nS=0:n_S-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/S_k_q_wkS__',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Slim_ = prctile(real(S_k_p_wkS__),[  0,100],'all');
Slim_ = mean(Slim_) + 1.25*0.5*diff(Slim_)*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nS = max(0,min(n_S-1,floor(n_S*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(S_k_q_wkS__(:,1+nS)),Slim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(S_k_q_wkS__(:,1+nS))/max(abs(S_k_q_wkS__(:)))),[-4,0],colormap_80s);
title(sprintf('nS %d',nS));
end;%for nl=0:15-1;
sgtitle(sprintf('S_k_p_wkS__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now set up images. ;
%%%%%%%%;
tmp_index_ = efind(abs(viewing_polar_a_S_-pi/2)<1e-12); %<-- equatorial viewing-angles. ;
n_M = numel(tmp_index_);
M_k_p_wkM__ = S_k_p_wkS__(:,1+tmp_index_);
euler_polar_a_M_ = viewing_polar_a_S_(1+tmp_index_);
euler_azimu_b_M_ = viewing_azimu_b_S_(1+tmp_index_);
euler_gamma_z_M_ = zeros(n_M,1);

%%%%%%%%;
% Now convert images to M_k_q_wkM__. ;
%%%%%%%%;
tmp_t = tic();
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_q_wkM__(:,1+nM) = interp_p_to_q(n_k_p_r,n_w_,n_w_sum,M_k_p_wkM__(:,1+nM));
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% M_k_q_wkM__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_k_q_wkM__',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
Mlim_ = prctile(real(S_k_p_wkS__),[  0,100],'all');
Mlim_ = mean(Mlim_) + 1.25*0.5*diff(Mlim_)*[-1,+1];
for nl=0:15-1;
subplot(3,5,1+nl); nM = max(0,min(n_M-1,floor(n_M*nl/(15-1))));
%imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,abs(M_k_q_wkM__(:,1+nM)),Mlim_,colormap_80s);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(abs(M_k_q_wkM__(:,1+nM))/max(abs(M_k_q_wkM__(:)))),[-4,0],colormap_80s);
title(sprintf('nM %d',nM));
end;%for nl=0:15-1;
sgtitle(sprintf('M_k_q_wkM__'),'Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/M_k_q_norm_',dir_ssnll);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;
subplot(2,2,1);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(S_k_q_norm_/max(S_k_q_norm_)),[-4,0],colormap_80s);
title('S_k_q_norm_','Interpreter','none');
subplot(2,2,2);
M_k_q_norm_ = sqrt(sum(abs(M_k_q_wkM__).^2,2));
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(M_k_q_norm_/max(M_k_q_norm_)),[-4,0],colormap_80s);
title('M_k_q_norm_','Interpreter','none');
subplot(2,2,3);
S_k_q_norm_ = sqrt(sum(abs(S_k_q_wkS__).^2,2)); tmp_eps = S_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,S_k_q_norm_-tmp_eps)/max(S_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('S_k_q_norm_ - tmp_eps','Interpreter','none');
subplot(2,2,4);
M_k_q_norm_ = sqrt(sum(abs(M_k_q_wkM__).^2,2)); tmp_eps = M_k_q_norm_(1+n_w_csum_(1+n_k_p_r-1)+n_w_(1+n_k_p_r-1)/2);
imagesc_q(n_k_p_r,k_p_r_,n_w_,n_w_sum,log10(max(0,M_k_q_norm_-tmp_eps)/max(M_k_q_norm_-tmp_eps)),[-4,0],colormap_80s);
title('M_k_q_norm_ - tmp_eps','Interpreter','none');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
 
%%%%%%%%;
hist2dab_k_eq_d = 0.5/k_p_r_max;
[ ...
 n_hist2dab_S ...
,hist2dab_azimu_b_S_ ...
,hist2dab_polar_a_S_ ...
,hist2dab_weight_S_ ...
,hist2dab_k_c_0_S_ ...
,hist2dab_k_c_1_S_ ...
,hist2dab_k_c_2_S_ ...
,n_hist2dab_polar_a ...
,hist2dab_polar_a_ ...
,n_hist2dab_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,hist2dab_k_eq_d ...
,'L' ... %<-- exclude pole. ;
,0 ... %<-- adaptive grid. ;
) ;
%%%%%%%%;
% generate refined templates hist2dab_S_k_p_wk_ on k_p_ grid with uniform n_w_. ;
%%%%%%%%;
template_k_eq_d = -1;
tmp_t = tic();
[ ...
 hist2dab_S_k_p_wkS__ ...
] = ...
pm_template_2( ...
 0*flag_verbose ...
,l_max ...
,n_k_p_r ...
,a_k_Y_quad_yk__ ...
,[] ...
,-1 ...
,n_w_max ...
,n_hist2dab_S ...
,hist2dab_azimu_b_S_ ...
,hist2dab_polar_a_S_ ...
,hist2dab_weight_S_ ...
,n_hist2dab_polar_a ...
,hist2dab_polar_a_ ...
,n_hist2dab_azimu_b_ ...
);
hist2dab_S_k_p_wkS__ = reshape(hist2dab_S_k_p_wkS__,[n_w_sum,n_hist2dab_S]);
tmp_t = toc(tmp_t); if (flag_verbose>0); disp(sprintf(' %% hist2dab_S_k_p_wkS__ (pm_template_2): %0.6fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% ')); end;
if (flag_verbose>0); disp(sprintf(' %% construct weight_3d_riesz')); end;
%%%%%%%%;
weight_3d_riesz_k_p_r_ = weight_3d_k_p_r_;
weight_3d_riesz_k_all_ = weight_3d_k_all_;
for nk_p_r=0:n_k_p_r-1;
k_p_r = k_p_r_(1+nk_p_r);
weight_3d_k_p_r = weight_3d_k_p_r_(1+nk_p_r);
weight_2d_k_p_r = weight_2d_k_p_r_(1+nk_p_r);
weight_3d_riesz_k_p_r_(1+nk_p_r) = weight_3d_k_p_r_(1+nk_p_r) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
tmp_index_ = n_k_all_csum_(1+nk_p_r):n_k_all_csum_(1+nk_p_r+1)-1;
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_k_all_(1+tmp_index_))/(weight_3d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_k_all_(1+tmp_index_))/(4*pi*weight_3d_k_p_r))); end;
weight_3d_riesz_k_all_(1+tmp_index_) = weight_3d_k_all_(1+tmp_index_) * weight_2d_k_p_r / max(1e-16,weight_3d_k_p_r);
if (flag_verbose> 1); disp(sprintf(' %% k_p_r %0.6f: sum(weight_3d_riesz_k_all_(1+tmp_index_))/(weight_2d_k_p_r * 4*pi): %0.16f',k_p_r,sum(weight_3d_riesz_k_all_(1+tmp_index_))/(4*pi*weight_2d_k_p_r))); end;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%%%%%;

%%%%%%%%;
% test_slice_vs_volume_integral_helper_eig_equa_band_3;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
