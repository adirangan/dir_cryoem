%%%%%%%%;
% Make fig for Beamer_empm_5.tex. ;
%%%%%%%%;

%%%%%%%%;
% first recapitulates test_pm_trpv1_2. ;
%%%%%%%%;
test_slice_vs_volume_integral_trpv1_11;

flag_invert=0;
fname_emd = sprintf('%s/%s',dir_data_star,fname_nopath_volume);
a_x_c_load_ = cast(ReadMRC(fname_emd),'double');
n_x_c = numel(x_c_0_); x_c_2_ = x_c_1_;
[x_c_0___,x_c_1___,x_c_2___] = ndgrid(x_c_0_,x_c_1_,x_c_2_); n_xxx_c = n_x_c.^3; xxx_c_weight = (2*x_p_r_max/n_x_c)^3;
%%%%;
n_k_c = 1+128; n_kkk_c = n_k_c.^3;
k_c_0_3d_ = linspace(-1.5*k_p_r_max,+1.5*k_p_r_max,n_k_c);
k_c_1_3d_ = linspace(-1.5*k_p_r_max,+1.5*k_p_r_max,n_k_c);
k_c_2_3d_ = linspace(-1.5*k_p_r_max,+1.5*k_p_r_max,n_k_c);
[k_c_0_3d___,k_c_1_3d___,k_c_2_3d___] = ndgrid(k_c_0_3d_,k_c_1_3d_,k_c_2_3d_);
eta = pi/x_p_r_max; tmp_t = tic;
tmp_t = tic();
a_k_c_load_ = xxnufft3d3(n_xxx_c,x_c_0___(:)*eta,x_c_1___(:)*eta,x_c_2___(:)*eta,a_x_c_load_(:).*xxx_c_weight,-1,1e-12,n_kkk_c,2*pi*k_c_0_3d___(:)/eta,2*pi*k_c_1_3d___(:)/eta,2*pi*k_c_2_3d___(:)/eta)/sqrt(2*pi)/sqrt(2*pi)/sqrt(2*pi);
tmp_t = toc(tmp_t); disp(sprintf(' %% xxnufft3d3: a_k_c_load_ time %0.2fs',tmp_t));
a_k_c_load_ = reshape(a_k_c_load_,[n_k_c,n_k_c,n_k_c]);
%%%%;
a_x_c_zavg_ = repmat(mean(a_x_c_load_,3),[1,1,n_x_c]);
tmp_sigma_g = 4.0/(2*pi); tmp_g_3d_ = (1/sqrt(2*pi) / tmp_sigma_g) * exp(-(k_c_2_3d_.^2)/(2*tmp_sigma_g.^2)); tmp_g_3d_ = tmp_g_3d_/max(1e-12,max(tmp_g_3d_));
a_k_c_zavg_ = bsxfun(@times,a_k_c_load_,reshape(tmp_g_3d_,[1,1,n_k_c]));
%%%%;

%%%%;
k_c_0_2d_ = zeros(n_w_sum,1);
k_c_1_2d_ = zeros(n_w_sum,1);
na=0;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
n_w_csum = n_w_csum_(1+nk_p_r);
tmp_k_p_r = k_p_r_(1+nk_p_r);
for nw=0:n_w-1;
tmp_gamma = (2*pi*nw)/n_w;
tmp_k_c_0 = tmp_k_p_r * cos(tmp_gamma);
tmp_k_c_1 = tmp_k_p_r * sin(tmp_gamma);
k_c_0_2d_(1+na) = tmp_k_c_0;
k_c_1_2d_(1+na) = tmp_k_c_1;
na = na+1;
end;%for nw=0:n_w-1;
end;%for nk_p_r=0:n_k_p_r-1;
%%%%;

%%%%%%%%;
% Now generates figure. ;
%%%%%%%%;

n_M = 1024;
[ ...
 M_x_c___ ...
,index_nCTF_from_nM_ ...
,index_nM_from_nCTF_ ...
,Voltage_CTF_ ...
,DefocusU_CTF_ ...
,DefocusV_CTF_ ...
,DefocusAngle_CTF_ ...
,SphericalAberration_CTF_ ...
,AmplitudeContrast_CTF_ ...
] = ...
rlnImageName_from_star_1( ...
 dir_data_star ...
,fname_nopath_star ...
,n_M ...
);
if (fnorm(Voltage_CTF_)< 1e-3); disp(sprintf(' %% Warning, Voltage not set, setting Voltage to 300kV')); Voltage_CTF_ = 300*ones(n_M,1); end;
if flag_invert; M_x_c___ = -M_x_c___; end;
%%%%%%%%;
% Remove any edge artefacts, mean center and normalize each image. ;
%%%%%%%%;
disp(sprintf(' %% Removing edge-artefacts'));
n_M_ext_ = zeros(n_M,1);
for nM=0:n_M-1;
if (mod(nM,128)==0); disp(sprintf(' %% nM %d/%d',nM,n_M)); end;
n_pixel = 4; edge_tolerance = 0.5; n_edge_overshoot = 8; rseed = 0;
[M_x_c___(:,:,1+nM),n_M_ext_(1+nM)] = image_replace_edge_artefact_0(M_x_c___(:,:,1+nM),4,0.5,2,0);
end;%for nM=0:n_M-1;
disp(sprintf(' %% edge-artefacts detected in %d/%d images.',numel(find(n_M_ext_>0)),n_M));

%%%%%%%%;
% find index_nS_from_nM_. ;
% Also construct T_k_p__ (templates not shifted to align with images, but are rotated). ;
% And construct U_k_p__ (images shifted to align with templates, but not rotated). ;
%%%%%%%%;
n_S = n_viewing_S;
T_k_p__ = zeros(n_w_sum,n_M);
U_k_p__ = zeros(n_w_sum,n_M);
X_TU_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
viewing_k_c_0_S_ = sin(viewing_polar_a_S_).*cos(viewing_azimu_b_S_);
viewing_k_c_1_S_ = sin(viewing_polar_a_S_).*sin(viewing_azimu_b_S_);
viewing_k_c_2_S_ = cos(viewing_polar_a_S_);
for nM=0:n_M-1;
tmp_euler_polar_a = +euler_polar_a_empi_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_empi_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_empi_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_empi_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_empi_(1+nM);
M_k_p_ = M_k_p_wkM__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_S_,viewing_k_c_1_S_,viewing_k_c_2_S_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
index_nS_from_nM_(1+nM) = nS;
S_k_p_ = S_k_p_wkS__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
%T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p_wkC__(:,1+index_nCTF_from_nM_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
T_k_p__(:,1+nM) = T_k_p_;
U_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,M_k_p_,+tmp_image_delta_x,+tmp_image_delta_y);
%U_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,U_k_p_,-tmp_euler_gamma_z);
U_k_p__(:,1+nM) = U_k_p_;
tmp_TT = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,T_k_p_);
tmp_UU = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,U_k_p_,U_k_p_);
tmp_TU = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,T_k_p_,U_k_p_);
X_TU_(1+nM) = real(tmp_TU)/sqrt(tmp_TT*tmp_UU);
end;%for nM=0:n_M-1;
%%%%%%%%;

index_nS_from_nM_SM__ = sparse(1+index_nS_from_nM_,1:n_M,1,n_S,n_M);
[~,nS_target] = max(sum(index_nS_from_nM_SM__,2)); nS_target=nS_target-1;
index_nM_target_ = efind(index_nS_from_nM_==nS_target);
%%%%;
nMA=index_nM_target_(1+0); %<-- first image. ;
euler_polar_a_empi_nMA = euler_polar_a_empi_(1+nMA);
euler_azimu_b_empi_nMA = euler_azimu_b_empi_(1+nMA);
euler_gamma_z_empi_nMA = euler_gamma_z_empi_(1+nMA);
image_delta_x_empi_nMA = image_delta_x_empi_(1+nMA);
image_delta_y_empi_nMA = image_delta_y_empi_(1+nMA);
%%%%;
nMB=index_nM_target_(1+3); %<-- second image. ;
euler_polar_a_empi_nMB = euler_polar_a_empi_(1+nMB);
euler_azimu_b_empi_nMB = euler_azimu_b_empi_(1+nMB);
euler_gamma_z_empi_nMB = euler_gamma_z_empi_(1+nMB);
image_delta_x_empi_nMB = image_delta_x_empi_(1+nMB);
image_delta_y_empi_nMB = image_delta_y_empi_(1+nMB);
%%%%;
nMC=index_nM_target_(1+5); %<-- third image. ;
euler_polar_a_empi_nMC = euler_polar_a_empi_(1+nMC);
euler_azimu_b_empi_nMC = euler_azimu_b_empi_(1+nMC);
euler_gamma_z_empi_nMC = euler_gamma_z_empi_(1+nMC);
image_delta_x_empi_nMC = image_delta_x_empi_(1+nMC);
image_delta_y_empi_nMC = image_delta_y_empi_(1+nMC);
%%%%;
nMD=index_nM_target_(1+7); %<-- fourth image. ;
euler_polar_a_empi_nMD = euler_polar_a_empi_(1+nMD);
euler_azimu_b_empi_nMD = euler_azimu_b_empi_(1+nMD);
euler_gamma_z_empi_nMD = euler_gamma_z_empi_(1+nMD);
image_delta_x_empi_nMD = image_delta_x_empi_(1+nMD);
image_delta_y_empi_nMD = image_delta_y_empi_(1+nMD);
%%%%;

tmp_n_x = round(k_p_r_max*2*pi*sqrt(2));
tmp_index_ = n_x_u/2 + [-tmp_n_x:+tmp_n_x];
%%%%;
M_x_c_nMA_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMA);
S_x_c_nMA_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMA)));
S_x_c_nMA_ = S_x_c_nMA_(1+tmp_index_,1+tmp_index_);
U_k_p_nMA_ = U_k_p__(:,1+nMA);
T_k_p_nMA_ = T_k_p__(:,1+nMA);
%%%%;
M_x_c_nMB_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMB);
S_x_c_nMB_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMB)));
S_x_c_nMB_ = S_x_c_nMB_(1+tmp_index_,1+tmp_index_);
U_k_p_nMB_ = U_k_p__(:,1+nMB);
T_k_p_nMB_ = T_k_p__(:,1+nMB);
%%%%;
M_x_c_nMC_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMC);
S_x_c_nMC_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMC)));
S_x_c_nMC_ = S_x_c_nMC_(1+tmp_index_,1+tmp_index_);
U_k_p_nMC_ = U_k_p__(:,1+nMC);
T_k_p_nMC_ = T_k_p__(:,1+nMC);
%%%%;
M_x_c_nMD_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMD);
S_x_c_nMD_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMD)));
S_x_c_nMD_ = S_x_c_nMD_(1+tmp_index_,1+tmp_index_);
U_k_p_nMD_ = U_k_p__(:,1+nMD);
T_k_p_nMD_ = T_k_p__(:,1+nMD);
%%%%;

dir_fig = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig',string_root);
fname_fig = sprintf('%s/empm_fig_image_sample_FIGB_',dir_fig);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
fontsize_use = 18;
p_row = 2; p_col = 4; ns=0;
%%%%%%%%;
ns=0; subplot_{1+0} = subplot(p_row,p_col,1+ns);
M_x_c_lim_ = mean(M_x_c_nMA_,'all') + std(M_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMA_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
ns=1; subplot_{1+1} = subplot(p_row,p_col,1+ns);
M_x_c_lim_ = mean(M_x_c_nMB_,'all') + std(M_x_c_nMB_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMB_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
ns=4; subplot_{1+2} = subplot(p_row,p_col,1+ns);
M_x_c_lim_ = mean(M_x_c_nMC_,'all') + std(M_x_c_nMC_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMC_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMC,1+nMC),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
ns=5; subplot_{1+3} = subplot(p_row,p_col,1+ns);
M_x_c_lim_ = mean(M_x_c_nMD_,'all') + std(M_x_c_nMD_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMD_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMD,1+nMD),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
ns_=[2,3,6,7];
subplot_{1+4} = subplot(p_row,p_col,1+ns_);
tmp_index_ = 128 + [-64:+64]; isosurface_f_x_u_0(a_x_c_load_(1+tmp_index_,1+tmp_index_,1+tmp_index_),95); axisnotick;
title(sprintf('$F^{\\mbox{signal}}(\\vec{x})$'),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%%%%%;
colormap(subplot_{1+0},colormap_80s);
colormap(subplot_{1+1},colormap_80s);
colormap(subplot_{1+2},colormap_80s);
colormap(subplot_{1+3},colormap_80s);
colormap(subplot_{1+4},colormap_80s);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

fname_fig = sprintf('%s/empm_fig_SliceTheorem_2d_FIGB_',dir_fig);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,1024]);fig80s;
fontsize_use = 18;
colormap_80s_b__ = colormap_80s(); colormap_80s_b__(1,:) = [1,1,1];
p_row = 2; p_col = 2; ns=0;
%%%%%%%%;
ns=0; subplot_{1+ns} = subplot(p_row,p_col,1+ns);
S_x_c_lim_ = mean(S_x_c_nMA_,'all') + std(S_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(S_x_c_nMA_,S_x_c_lim_);
title(sprintf('$B_{%d} := A_{%d}^{\\mbox{signal}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
ns=1; subplot_{1+ns} = subplot(p_row,p_col,1+ns);
T_x_c_lim_ = 0*mean(real(T_k_p_nMA_),'all') + std(real(T_k_p_nMA_),1,'all')*2.5*[-1,+1];
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMA_),T_x_c_lim_,colormap_80s);
title(sprintf('${\\cal R} [ \\hat{B}_{%d}(\\vec{k}) := FFT\\circ B_{%d}(\\vec{x}) ]$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
ns=2; subplot_{1+ns} = subplot(p_row,p_col,1+ns);
tmp_S_x_c_nMA_ = min(S_x_c_lim_)*ones(size(S_x_c_nMA_));
tmp_index_ = floor(size(S_x_c_nMA_,1)/2) + [-4:+4];
tmp_S_x_c_nMA_(1+tmp_index_,:) = repmat(mean(S_x_c_nMA_,1),[numel(tmp_index_),1]);
imagesc(tmp_S_x_c_nMA_,S_x_c_lim_);
title(sprintf('$\\int dy \\ B_{%d}(\\vec{x}) $',1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
ns=3; subplot_{1+ns} = subplot(p_row,p_col,1+ns);
T_x_c_lim_ = 0*mean(real(T_k_p_nMA_),'all') + std(real(T_k_p_nMA_),1,'all')*2.5*[-1,+1];
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMA_),T_x_c_lim_,colormap_80s);
hold on;
tmp_k = max(k_p_r_); tmp_e = 0.25;
p=patch([-tmp_k;+tmp_k;+tmp_k;-tmp_k],[+tmp_e;+tmp_e;+tmp_k;+tmp_k],'w','EdgeColor','none');
p=patch([-tmp_k;+tmp_k;+tmp_k;-tmp_k],[-tmp_e;-tmp_e;-tmp_k;-tmp_k],'w','EdgeColor','none');
title(sprintf('${\\cal R} [ FFT \\circ [ \\int dy \\ B_{%d}(\\vec{x}) ] ] ]$',1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%%%%%;
colormap(subplot_{1+0},colormap_80s);
colormap(subplot_{1+1},colormap_80s);
colormap(subplot_{1+2},colormap_80s_b__);
colormap(subplot_{1+3},colormap_80s_b__);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));



