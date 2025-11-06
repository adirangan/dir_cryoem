%%%%%%%%;
% designed to test rotations. ;
%%%%%%%%;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

str_thisfunction = 'test_rotate_spharm_to_spharm_4';
flag_verbose = 1; rng(0);
flag_disp = 1; nf=0;
k_int = 48;
k_eq_d_double = 0.25;
dir_jpg = sprintf('/%s/rangan/dir_cryoem/dir_CryoBIFE_MD',string_root);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [entering %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

fname_map_png = sprintf('/%s/rangan/dir_cryoem/World_elevation_map.png',string_root);
map_ab__ = imread(fname_map_png);
map_ab__ = mean(map_ab__,3);
[n_a,n_b] = size(map_ab__);
n_10 = 10; %<-- downsample by 10. ;
assert(mod(n_a,n_10)==0);
map_ab__ = reshape(mean(reshape(map_ab__,[n_10,n_a*n_b/n_10]),1),[n_a/n_10,n_b]);
assert(mod(n_b,n_10)==0);
map_ba__ = transpose(reshape(mean(reshape(transpose(map_ab__),[n_10,n_b*n_a/n_10/n_10]),1),[n_b/n_10,n_a/n_10]));
map_ab__ = cast(map_ab__,'double'); 

[n_a,n_b] = size(map_ab__);
a_ = transpose(linspace(0,1*pi,n_a));
b_ = transpose(linspace(0,2*pi,n_b+1)); b_ = b_(1:n_b);
[a__,b__] = meshgrid(a_,b_); %<-- for later use in fixed.interp2. ;

%%%%%%%%;
k_p_r_max = k_int/(2*pi); k_eq_d = k_eq_d_double/(2*pi);
flag_uniform_over_polar_a = 1;
[ ...
 n_qk ...
,k_p_azimu_b_qk_ ...
,k_p_polar_a_qk_ ...
,weight_shell_k_ ...
,k_c_0_qk_ ...
,k_c_1_qk_ ...
,k_c_2_qk_ ...
] = ...
sample_shell_6( ...
 k_p_r_max ...
,k_eq_d ...
,'L' ...
,flag_uniform_over_polar_a ...
) ;
%%%%%%%%;
%map_qk_ = interp2(b_,a_,flipud(map_ab__),k_p_azimu_b_qk_,k_p_polar_a_qk_);
map_qk_ = interp2(b_,a_,map_ab__,k_p_azimu_b_qk_,k_p_polar_a_qk_);

%%%%%%%%;
n_k_p_r = 1;
n_qk_csum_ = [0;n_qk];
k_p_r_qk_ = k_p_r_max*ones(n_qk,1);
k_p_r_ = k_p_r_max;
weight_3d_k_p_r_ = 1;
weight_3d_k_p_qk_ = weight_shell_k_;
%%%%%%%%;
l_max_upb = k_int;
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
if (flag_verbose); disp(sprintf(' %% n_qk %d, l_max_max %d',n_qk,l_max_max)); end;
%%%%%%%%;

%%%%%%%%;
% set a_k_p_true_ to be map_qk_ and then project onto spherical-harmonic subspace. ;
%%%%%%%%;
a_k_p_true_ = map_qk_;
tmp_t = tic;
[a_k_Y_true_] = convert_k_p_to_spharm_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_true_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_true_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_p_reco_] = real(convert_spharm_to_k_p_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_true_ --> a_k_p_reco_ time %0.2fs',tmp_t));
tmp_t = tic;
[a_k_Y_reco_] = convert_k_p_to_spharm_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,a_k_p_reco_);
tmp_t = toc(tmp_t); disp(sprintf(' %% a_k_Y_reco_ time %0.2fs',tmp_t));
fnorm_disp(flag_verbose,'a_k_p_true_',a_k_p_true_,'a_k_p_reco_',a_k_p_reco_,' %<-- can be large');
fnorm_disp(flag_verbose,'a_k_Y_true_',a_k_Y_true_,'a_k_Y_reco_',a_k_Y_reco_,' %<-- should be small');
%%%%%%%%;

%%%%%%%%;
eulerbypi_d_ = [+0.2;-0.3;-0.4];
%%%%%%%%;
eulerbypi_b_ = [eulerbypi_d_(1+0);0;0];
b_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,eulerbypi_b_*pi);
tmp_t = tic;
[b_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,b_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% b_k_Y_true_ --> b_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
eulerbypi_c_ = [eulerbypi_d_(1+0);eulerbypi_d_(1+1);0];
c_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,eulerbypi_c_*pi);
tmp_t = tic;
[c_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,c_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% c_k_Y_true_ --> c_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
d_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,a_k_Y_true_,eulerbypi_d_*pi);
tmp_t = tic;
[d_k_p_quad_] = real(convert_spharm_to_k_p_1(flag_verbose,n_qk,n_qk_csum_,k_p_r_qk_,k_p_azimu_b_qk_,k_p_polar_a_qk_,weight_3d_k_p_qk_,weight_shell_k_,n_k_p_r,k_p_r_,weight_3d_k_p_r_,l_max_,d_k_Y_true_));
tmp_t = toc(tmp_t); disp(sprintf(' %% d_k_Y_true_ --> d_k_p_quad_ time %0.2fs',tmp_t));
%%%%%%%%;
e_k_Y_true_ = a_k_Y_true_;
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[eulerbypi_d_(1+0);0;0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;eulerbypi_d_(1+1);0]*pi);
e_k_Y_true_ = rotate_spharm_to_spharm_3(n_k_p_r,k_p_r_,l_max_,e_k_Y_true_,[0;0;eulerbypi_d_(1+2)]*pi);
disp(sprintf(' %% e_k_Y_true vs d_k_Y_true error: %0.16f',fnorm(e_k_Y_true_-d_k_Y_true_)/fnorm(e_k_Y_true_)));
%%%%%%%%;

%%%%%%%%;
% plot results. ;
%%%%%%%%;

str_step0 = sprintf('original $F$');
str_step1 = sprintf('$R_{z}(\\gamma=%+0.1f\\pi)F$',eulerbypi_d_(1+0));
str_step2 = sprintf('$R_{y}(\\beta=%+0.1f\\pi)R_{z}(\\gamma=%+0.1f\\pi)F$',eulerbypi_d_(1+1),eulerbypi_d_(1+0));
str_step3 = sprintf('$R_{z}(\\alpha=%+0.1f\\pi)R_{y}(\\beta=%+0.1f\\pi)R_{z}(\\gamma=%+0.1f\\pi)F$',eulerbypi_d_(1+2),eulerbypi_d_(1+1),eulerbypi_d_(1+0));
disp(sprintf('%s',str_step0));
disp(sprintf('%s',str_step1));
disp(sprintf('%s',str_step2));
disp(sprintf('%s',str_step3));

%%%%%%%%;
% define rotations in 3d. ;
%%%%%%%%;
Rz = @(azimu_b) ...
[ +cos(azimu_b) -sin(azimu_b) 0 ; ...
  +sin(azimu_b) +cos(azimu_b) 0 ; ...
   0             0            1 ; ...
] ;
%%%%%%%%;
Ry = @(polar_a) ...
[ +cos(polar_a) 0 +sin(polar_a) ; ...
   0            1  0            ; ...
  -sin(polar_a) 0 +cos(polar_a) ; ...
];
%%%%%%%%;
pole_0_ = [0;0;1];
pole_1_ = Rz(pi*eulerbypi_d_(1+0))*pole_0_;
pole_2_ = Ry(pi*eulerbypi_d_(1+1))*pole_1_;
pole_3_ = Rz(pi*eulerbypi_d_(1+2))*pole_2_;
r_pole = 1.125;

sealevel = 142.5 + 0.0; maplim_ = sealevel + 64*[-1,+1];
c_base__ = colormap(parula); close(gcf); 
c_81s__ = colormap_81s(256);
c_elev__ = 0.65*c_base__.^1.0 + 0.35*c_81s__;
% test with: ;
% imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(a_k_p_reco_),maplim_,c_elev__); xlim([0,2*pi]); ylim([0,1*pi]); set(gca,'ydir','reverse'); set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90); set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');

%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
flag_2d_vs_3d = 1;
fontsize_use = 12;
subplot(2,2,1);
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(a_k_p_reco_),maplim_,c_elev__);
xlim([0,2*pi]); ylim([0,1*pi]); set(gca,'ydir','reverse');
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
%title(sprintf('\\tau [0;0;0]\\pi'));
title(str_step0,'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
subplot(2,2,2);
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(b_k_p_quad_),maplim_,c_elev__);
xlim([0,2*pi]); ylim([0,1*pi]); set(gca,'ydir','reverse');
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_b_));
title(str_step1,'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
subplot(2,2,3);
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(c_k_p_quad_),maplim_,c_elev__);
xlim([0,2*pi]); ylim([0,1*pi]); set(gca,'ydir','reverse');
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_c_));
title(str_step2,'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
subplot(2,2,4);
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(d_k_p_quad_),maplim_,c_elev__);
xlim([0,2*pi]); ylim([0,1*pi]); set(gca,'ydir','reverse');
set(gca,'XTick',linspace(0,2*pi,21),'XTickLabel',linspace(0,2,21)); xlabel('azimu/\pi'); xtickangle(90);
set(gca,'YTick',linspace(0,1*pi,11),'YTickLabel',linspace(0,1,11)); ylabel('polar/\pi');
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_d_));
title(str_step3,'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
figbig;
fname_fig_pre = sprintf('%s/test_rotate_spharm_to_spharm_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
%%%%%%%%;

%%%%%%%%;
flag_2d_vs_3d = 0;
sphere_grid_k_max = 1.0 + 1.0/64;
fontsize_use = 24;
linewidth_use = 4*2;
markersize_use = 12*2;
%%%%;
figure(1+nf);nf=nf+1;clf;
subplot(1,1,1);%subplot(2,2,1);
hold on; plot_sphere_grid_0(struct('k_max',sphere_grid_k_max));
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(a_k_p_reco_),maplim_,c_elev__,flag_2d_vs_3d);
hold on;
plot3(r_pole*[0;pole_0_(1+0)],r_pole*[0;pole_0_(1+1)],r_pole*[0;pole_0_(1+2)],'k-','Linewidth',linewidth_use);
plot3(r_pole*[pole_0_(1+0)],r_pole*[pole_0_(1+1)],r_pole*[pole_0_(1+2)],'ko','MarkerSize',markersize_use,'MarkerFaceColor','r');
xlabel('x');ylabel('y');zlabel('z'); axis equal;% axis vis3d;
%title(sprintf('\\tau [0;0;0]\\pi'));
title(str_step0,'Interpreter','latex');
set(gca,'XTick',[-1,+1],'XTickLabel',{'-','+'});
set(gca,'YTick',[-1,+1],'YTickLabel',{'-','+'});
set(gca,'ZTick',[-1,+1],'ZTickLabel',{'-','+'});
set(gca,'FontSize',fontsize_use);
figbig;
fname_fig_pre = sprintf('%s/test_rotate_spharm_to_spharm_FIGB00',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
%%%%;
figure(1+nf);nf=nf+1;clf;
subplot(1,1,1);%subplot(2,2,2);
hold on; plot_sphere_grid_0(struct('k_max',sphere_grid_k_max));
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(b_k_p_quad_),maplim_,c_elev__,flag_2d_vs_3d);
hold on;
plot3(r_pole*[0;pole_1_(1+0)],r_pole*[0;pole_1_(1+1)],r_pole*[0;pole_1_(1+2)],'k-','Linewidth',linewidth_use);
plot3(r_pole*[pole_1_(1+0)],r_pole*[pole_1_(1+1)],r_pole*[pole_1_(1+2)],'ko','MarkerSize',markersize_use,'MarkerFaceColor','r');
xlabel('x');ylabel('y');zlabel('z'); axis equal;% axis vis3d;
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_b_));
title(str_step1,'Interpreter','latex');
set(gca,'XTick',[-1,+1],'XTickLabel',{'-','+'});
set(gca,'YTick',[-1,+1],'YTickLabel',{'-','+'});
set(gca,'ZTick',[-1,+1],'ZTickLabel',{'-','+'});
set(gca,'FontSize',fontsize_use);
figbig;
fname_fig_pre = sprintf('%s/test_rotate_spharm_to_spharm_FIGB01',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
%%%%;
figure(1+nf);nf=nf+1;clf;
subplot(1,1,1);%subplot(2,2,3);
hold on; plot_sphere_grid_0(struct('k_max',sphere_grid_k_max));
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(c_k_p_quad_),maplim_,c_elev__,flag_2d_vs_3d);
hold on;
plot3(r_pole*[0;pole_2_(1+0)],r_pole*[0;pole_2_(1+1)],r_pole*[0;pole_2_(1+2)],'k-','Linewidth',linewidth_use);
plot3(r_pole*[pole_2_(1+0)],r_pole*[pole_2_(1+1)],r_pole*[pole_2_(1+2)],'ko','MarkerSize',markersize_use,'MarkerFaceColor','r');
xlabel('x');ylabel('y');zlabel('z'); axis equal;% axis vis3d;
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_c_));
title(str_step2,'Interpreter','latex');
set(gca,'XTick',[-1,+1],'XTickLabel',{'-','+'});
set(gca,'YTick',[-1,+1],'YTickLabel',{'-','+'});
set(gca,'ZTick',[-1,+1],'ZTickLabel',{'-','+'});
set(gca,'FontSize',fontsize_use);
figbig;
fname_fig_pre = sprintf('%s/test_rotate_spharm_to_spharm_FIGB10',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
%%%%;
figure(1+nf);nf=nf+1;clf;
subplot(1,1,1);%subplot(2,2,4);
hold on; plot_sphere_grid_0(struct('k_max',sphere_grid_k_max));
imagesc_polar_a_azimu_b_0(k_p_polar_a_qk_,k_p_azimu_b_qk_,real(d_k_p_quad_),maplim_,c_elev__,flag_2d_vs_3d);
hold on;
plot3(r_pole*[0;pole_3_(1+0)],r_pole*[0;pole_3_(1+1)],r_pole*[0;pole_3_(1+2)],'k-','Linewidth',linewidth_use);
plot3(r_pole*[pole_3_(1+0)],r_pole*[pole_3_(1+1)],r_pole*[pole_3_(1+2)],'ko','MarkerSize',markersize_use,'MarkerFaceColor','r');
xlabel('x');ylabel('y');zlabel('z'); axis equal;% axis vis3d;
%title(sprintf('\\tau [%+0.2f;%+0.2f;%+0.2f]\\pi',euler_d_));
title(str_step3,'Interpreter','latex');
set(gca,'XTick',[-1,+1],'XTickLabel',{'-','+'});
set(gca,'YTick',[-1,+1],'YTickLabel',{'-','+'});
set(gca,'ZTick',[-1,+1],'ZTickLabel',{'-','+'});
set(gca,'FontSize',fontsize_use);
figbig;
fname_fig_pre = sprintf('%s/test_rotate_spharm_to_spharm_FIGB11',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
disp(sprintf(' %% Writing %s',fname_fig_pre));
print('-djpeg',fname_fig_jpg);
close(gcf);
%%%%%%%%;
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
if (flag_verbose>0); disp(sprintf(' %% [finished %s]',str_thisfunction)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

disp('returning'); return;
