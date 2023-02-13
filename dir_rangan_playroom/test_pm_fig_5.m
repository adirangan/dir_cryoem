%%%%%%%%;
% makes figs for Principled_Marching_2?.tex. ;
%%%%%%%%;

%%%%%%%%;
% first recapitulates test_pm_trpv1_2. ;
%%%%%%%%;
test_pm_trpv1c_9b;
flag_invert=0;

%%%%%%%%;
% Now generates figures for paper. ;
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
n_S = n_viewing_all;
T_k_p__ = zeros(n_w_sum,n_M);
U_k_p__ = zeros(n_w_sum,n_M);
X_TU_ = zeros(n_M,1);
index_nS_from_nM_ = zeros(n_M,1);
viewing_k_c_0_all_ = sin(viewing_polar_a_all_).*cos(viewing_azimu_b_all_);
viewing_k_c_1_all_ = sin(viewing_polar_a_all_).*sin(viewing_azimu_b_all_);
viewing_k_c_2_all_ = cos(viewing_polar_a_all_);
for nM=0:n_M-1;
tmp_euler_polar_a = +euler_polar_a_true_(1+nM);
tmp_euler_azimu_b = +euler_azimu_b_true_(1+nM);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
tmp_image_delta_x = +1.0*image_delta_x_true_(1+nM);
tmp_image_delta_y = +1.0*image_delta_y_true_(1+nM);
M_k_p_ = M_k_p__(:,1+nM);
tmp_k_c_0 = sin(tmp_euler_polar_a)*cos(tmp_euler_azimu_b);
tmp_k_c_1 = sin(tmp_euler_polar_a)*sin(tmp_euler_azimu_b);
tmp_k_c_2 = cos(tmp_euler_polar_a);
nS = knnsearch([viewing_k_c_0_all_,viewing_k_c_1_all_,viewing_k_c_2_all_],[tmp_k_c_0,tmp_k_c_1,tmp_k_c_2]) - 1;
index_nS_from_nM_(1+nM) = nS;
S_k_p_ = S_k_p__(:,1+nS);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,S_k_p_,+tmp_euler_gamma_z);
%T_k_p_ = transf_p_to_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,T_k_p_,-tmp_image_delta_x,-tmp_image_delta_y);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_index_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
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
nMA=index_nM_target_(1+0); %<-- first image. ;
euler_polar_a_true_nMA = euler_polar_a_true_(1+nMA);
euler_azimu_b_true_nMA = euler_azimu_b_true_(1+nMA);
euler_gamma_z_true_nMA = euler_gamma_z_true_(1+nMA);
image_delta_x_true_nMA = image_delta_x_true_(1+nMA);
image_delta_y_true_nMA = image_delta_y_true_(1+nMA);
nSA = index_nS_from_nM_(1+nMA);
nMB=index_nM_target_(1+5); %<-- second image. ;
euler_polar_a_true_nMB = euler_polar_a_true_(1+nMB);
euler_azimu_b_true_nMB = euler_azimu_b_true_(1+nMB);
euler_gamma_z_true_nMB = euler_gamma_z_true_(1+nMB);
image_delta_x_true_nMB = image_delta_x_true_(1+nMB);
image_delta_y_true_nMB = image_delta_y_true_(1+nMB);

tmp_n_x = round(k_p_r_max*2*pi*sqrt(2));
tmp_index_ = n_x_u/2 + [-tmp_n_x:+tmp_n_x];
M_x_c_nMA_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMA);
M_x_c_nMB_ = M_x_c___(1+tmp_index_,1+tmp_index_,1+nMB);
S_x_c_nMA_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMA)));
S_x_c_nMA_ = S_x_c_nMA_(1+tmp_index_,1+tmp_index_);
S_x_c_nMB_ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p__(:,1+nMB)));
S_x_c_nMB_ = S_x_c_nMB_(1+tmp_index_,1+tmp_index_);
U_k_p_nMA_ = U_k_p__(:,1+nMA);
U_k_p_nMB_ = U_k_p__(:,1+nMB);
T_k_p_nMA_ = T_k_p__(:,1+nMA);
T_k_p_nMB_ = T_k_p__(:,1+nMB);

fname_fig = sprintf('%s_jpg/pm_fig_image_pair_FIGA_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
fontsize_use = 18;
p_row = 2; p_col = 4; ns=0;
%%%%%%%%;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
M_x_c_lim_ = mean(M_x_c_nMA_,'all') + std(M_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMA_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
S_x_c_lim_ = mean(S_x_c_nMA_,'all') + std(S_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(S_x_c_nMA_,S_x_c_lim_);
title(sprintf('$B_{%d} := A_{%d}^{\\mbox{signal}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(U_k_p_nMA_),[],colormap_80s);
title(sprintf('$\\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k}) + \\hat{A}_{%d}^{\\mbox{noise}}(\\vec{k})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMA_),[],colormap_80s);
title(sprintf('$\\hat{B}_{%d} := \\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
%%%%%%%%;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
M_x_c_lim_ = mean(M_x_c_nMB_,'all') + std(M_x_c_nMB_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMB_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
S_x_c_lim_ = mean(S_x_c_nMB_,'all') + std(S_x_c_nMB_,1,'all')*2.5*[-1,+1];
imagesc(S_x_c_nMB_,S_x_c_lim_);
title(sprintf('$B_{%d} := A_{%d}^{\\mbox{signal}}(\\vec{x})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(U_k_p_nMB_),[],colormap_80s);
title(sprintf('$\\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k}) + \\hat{A}_{%d}^{\\mbox{noise}}(\\vec{k})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMB_),[],colormap_80s);
title(sprintf('$\\hat{B}_{%d} := \\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
%%%%;
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Define pair of rings with high discriminability. ;
%%%%%%%%;
psi_ = linspace(0,2*pi,1+n_w_max); psi_ = psi_(1:end-1); psi_ = psi_(:);
T_k_p_wk__ = reshape(T_k_p_nMA_,[n_w_max,n_k_p_r]);
nk_p_r_0 = 19;
nk_p_r_1 = 27;
%%%%%%%%;
% Now find cost-matrix. ;
%%%%%%%%;
[ ...
 X_2d_S2_d0__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 2 ...
,[1;1] ...
,[1;1] ...
,n_w_max*[1;1] ...
,1 ...
,reshape(T_k_p_wk__(:,1+[nk_p_r_0,nk_p_r_1]),[2*n_w_max,1]) ...
);
[UX_2d_S2__,SX_2d_S2__,VX_2d_S2__] = svds(X_2d_S2_d0__,2); SX_2d_S2_ = diag(SX_2d_S2__);
UX_T_k_p_wn__ = T_k_p_wk__(:,1+[nk_p_r_0,nk_p_r_1]) * UX_2d_S2__;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_S2_FIGB_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024,768]);fontsize_use = 16;
%title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C^{\\mbox{signal}}(\\vec{u}) = %0.3f$ , $C^{\\mbox{signal}}(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
%title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C(\\vec{u}) = %0.3f$ , $C(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
title(sprintf('$k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C(\\vec{u}) = %0.3f$ , $C(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
%%%%;
hold on;
plot(psi_,real(T_k_p_wk__(:,1+nk_p_r_0)),'k:','LineWidth',4);
plot(psi_,real(T_k_p_wk__(:,1+nk_p_r_1)),'k-','LineWidth',3);
plot(psi_,real(UX_T_k_p_wn__(:,1+0)),'r-','LineWidth',3);
plot(psi_,real(UX_T_k_p_wn__(:,1+1)),'c-','LineWidth',3);
plot(psi_,zeros(size(psi_)),'k-','LineWidth',0.5);
hold off;
xlim([0,2*pi]); xlabel('\psi'); set(gca,'XTick',[0,2*pi],'XTickLabel',{'0','2\pi'});
set(gca,'YTick',[]);
ylim(8e-3*[-1,+1]);ylabel('real value');
legend({ ...
    ,'$\hat{A}^{signal}(k_{1},\psi)$' ...
    ,'$\hat{A}^{signal}(k_{2},\psi)$' ...
    ,'$\left[\vec{u}^{T}\hat{A}^{signal}\right](\psi)$' ...
    ,'$\left[\vec{v}^{T}\hat{A}^{signal}\right](\psi)$' ...
   },'Interpreter','latex','Location','North');
set(gca,'FontSize',fontsize_use);
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Define pair of rings with high value but low discriminability. ;
%%%%%%%%;
psi_ = linspace(0,2*pi,1+n_w_max); psi_ = psi_(1:end-1); psi_ = psi_(:);
T_k_p_wk__ = reshape(T_k_p_nMA_,[n_w_max,n_k_p_r]);
nk_p_r_0 = 10;
nk_p_r_1 = 32;
%%%%%%%%;
% Now find cost-matrix. ;
%%%%%%%%;
[ ...
 X_2d_S2_d0__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 2 ...
,[1;1] ...
,[1;1] ...
,n_w_max*[1;1] ...
,1 ...
,reshape(T_k_p_wk__(:,1+[nk_p_r_0,nk_p_r_1]),[2*n_w_max,1]) ...
);
[UX_2d_S2__,SX_2d_S2__,VX_2d_S2__] = svds(X_2d_S2_d0__,2); SX_2d_S2_ = diag(SX_2d_S2__);
UX_T_k_p_wn__ = T_k_p_wk__(:,1+[nk_p_r_0,nk_p_r_1]) * UX_2d_S2__;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_S2_FIGC_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024,768]);fontsize_use = 16;
%title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C^{\\mbox{signal}}(\\vec{u}) = %0.3f$ , $C^{\\mbox{signal}}(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
%title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C(\\vec{u}) = %0.3f$ , $C(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
title(sprintf('$k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C(\\vec{u}) = %0.3f$ , $C(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
%%%%;
hold on;
plot(psi_,real(T_k_p_wk__(:,1+nk_p_r_0)),'k:','LineWidth',4);
plot(psi_,real(T_k_p_wk__(:,1+nk_p_r_1)),'k-','LineWidth',3);
plot(psi_,real(UX_T_k_p_wn__(:,1+0)),'r-','LineWidth',3);
plot(psi_,real(UX_T_k_p_wn__(:,1+1)),'c-','LineWidth',3);
plot(psi_,zeros(size(psi_)),'k-','LineWidth',0.5);
hold off;
xlim([0,2*pi]); xlabel('\psi'); set(gca,'XTick',[0,2*pi],'XTickLabel',{'0','2\pi'});
set(gca,'YTick',[]);
ylim(8e-3*[-1,+1]);ylabel('real value');
legend({ ...
    ,'$\hat{A}^{signal}(k_{1},\psi)$' ...
    ,'$\hat{A}^{signal}(k_{2},\psi)$' ...
    ,'$\left[\vec{u}^{T}\hat{A}^{signal}\right](\psi)$' ...
    ,'$\left[\vec{v}^{T}\hat{A}^{signal}\right](\psi)$' ...
   },'Interpreter','latex','Location','SouthEast');
set(gca,'FontSize',fontsize_use);
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now find template cost-matrix. ;
%%%%%%%%;
[ ...
 X_2d_Sall_d0__ ...
,X_weight_r_ ...
] = ...
principled_marching_empirical_cost_matrix_0( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_w_ ...
,1 ...
,reshape(T_k_p_wk__,[n_w_sum,1]) ...
);
[UX_2d_Sall__,SX_2d_Sall__,VX_2d_Sall__] = svds(X_2d_Sall_d0__,n_k_p_r); SX_2d_Sall_ = diag(SX_2d_Sall__);
UX_T_k_p_wn__ = T_k_p_wk__ * diag(X_weight_r_) * UX_2d_Sall__;
%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_Sall_FIGD_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figsml; fontsize_use = 16; markersize_use = 8;
n_svd = 8;
subplot(1,1,1);
plot(SX_2d_Sall_,'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',markersize_use);
set(gca,'XTick',1:n_svd); xlim([0.5,n_svd+0.5]);
xlabel('rank');
ylabel('singular-value','Interpreter','latex');
grid on;
title(sprintf('$C(\\vec{u})$ for the first %d ranks',n_svd),'Interpreter','latex');
set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_Sall_FIGD2_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figsml;
markersize_use = 12;
n_svd = 8;
subplot(1,1,1);
plot(SX_2d_Sall_,'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',markersize_use);
set(gca,'XTick',1:n_svd); xlim([0.5,n_svd+0.5]);
xlabel('rank');
ylabel('singular-value','Interpreter','latex');
grid on;
title(sprintf('$C(\\vec{u})$ for the first %d ranks',n_svd),'Interpreter','latex');
set(gca,'FontSize',18);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));

%%%%%%%%;
% Now calculate image-template innerproduct-landscape for each rank. ;
% Be sure to point out that l2-error is not the goal here. ;
%%%%%%%%;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max);
UX_T_k_p_wn__ = T_k_p_wk__ * diag(X_weight_r_) * UX_2d_Sall__;
tmp_n_M = 2; n_svd = 8;
n_UX_rank = n_svd; %<-- just to check dimensions. ;
X_wMn___ = zeros(n_w_max,tmp_n_M,n_UX_rank);
%%%%;
for nUX_rank=0:n_UX_rank-1;
pm_n_UX_rank = 1+nUX_rank;
pm_n_w_ = n_w_max*ones(pm_n_UX_rank,1);
pm_n_w_sum = sum(pm_n_w_);
UX_T_k_p_ = reshape(UX_T_k_p_wn__(:,1:pm_n_UX_rank),[pm_n_w_sum,1]);
UX_T_k_q_ = ...
interp_p_to_q( ...
 pm_n_UX_rank ...
,pm_n_w_ ...
,pm_n_w_sum ...
,UX_T_k_p_ ...
);
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M,M_k_q__(:,1+[nMA,nMB]),pm_n_UX_rank,UX_2d_Sall__,X_weight_r_);
parameter = struct('type','parameter');
[ ...
 parameter ...
,X_wSM___ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,1 ...
,UX_T_k_q_ ...
,1 ...
,tmp_n_M ...
,svd_VUXM_lwnM____ ...
,[1,1] ...
);
X_wMn___(:,:,1+nUX_rank) = squeeze(X_wSM___);
%%%%;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;

gamma_ = linspace(0,2*pi,1+n_w_max); gamma_ = gamma_(1:end-1); gamma_ = gamma_(:);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_X_wSM_FIGE_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,1024,768]);
fontsize_use = 16;
intensity_gamma = 0.10;
intensity_scale_ = (1-intensity_gamma) + intensity_gamma*transpose(linspace(-1,1,64)).^2;
intensity_scale_ = intensity_scale_./max(1e-12,intensity_scale_(end));
c_beach__ = bsxfun(@times,colormap_beach(),intensity_scale_); n_c_beach = size(c_beach__,1);
for tmp_nM=[0,1];
nM=nMA; if (tmp_nM==1); nM=nMB; end;
tmp_str = sprintf('\\hat{A}_{%d}',1+nM);
subplot(1,2,1+tmp_nM);
hold on;
for nUX_rank=0:n_UX_rank-1;
nc_beach = max(0,min(n_c_beach-1,floor(n_c_beach*nUX_rank/n_UX_rank)));
plot(gamma_,X_wMn___(:,1+tmp_nM,1+nUX_rank),'-','Color',c_beach__(1+nc_beach,:),'LineWidth',2);
end;%for nUX_rank=0:n_UX_rank-1;
plot(gamma_,zeros(size(gamma_)),'k-','LineWidth',0.5);
hold off;
xlim([0,2*pi]); xlabel('\gamma'); set(gca,'XTick',[0,2*pi],'XTickLabel',{'0','2\pi'});
ylabel('$\cal{X}(\gamma)$','Interpreter','latex');
set(gca,'YTick',[]);
%title(sprintf('$\\langle %s , CTF\\odot\\left[ R_{\\gamma}\\circ\\hat{S}\\right] \\rangle$ for the first %d ranks',tmp_str,n_svd),'Interpreter','latex');
title(sprintf('$\\langle %s , CTF\\odot\\left[ R_{\\gamma}\\circ\\hat{S}\\right] \\rangle$',tmp_str),'Interpreter','latex');
str_legend_ = cell(n_svd,1);
for nsvd=0:n_svd-1;
str_legend_{1+nsvd} = sprintf('rank %.2d',1+nsvd);
end;%for nsvd=0:n_svd-1;
legend(str_legend_,'Location','North','FontSize',12);
set(gca,'FontSize',fontsize_use);
end;%for tmp_nM=[0,1];
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
% Now calculate image-image innerproduct-landscape for each rank. ;
%%%%%%%%;
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max);
UX_M_k_q_nMA_wq__ = reshape(M_k_q__(:,1+nMA),[n_w_max,n_k_p_r])* diag(X_weight_r_) * UX_2d_Sall__;
UX_M_k_q_nMB_wq__ = reshape(M_k_q__(:,1+nMB),[n_w_max,n_k_p_r])* diag(X_weight_r_) * UX_2d_Sall__;
tmp_n_M = 1; n_svd = 8;
n_UX_rank = n_svd; %<-- just to check dimensions. ;
X_wn__ = zeros(n_w_max,n_UX_rank);
%%%%;
for nUX_rank=0:n_UX_rank-1;
pm_n_UX_rank = 1+nUX_rank;
pm_n_w_ = n_w_max*ones(pm_n_UX_rank,1);
pm_n_w_sum = sum(pm_n_w_);
UX_M_k_q_nMA_wq_ = reshape(UX_M_k_q_nMA_wq__(:,1:pm_n_UX_rank),[pm_n_w_sum,1]);
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M,M_k_q__(:,1+nMB),pm_n_UX_rank,UX_2d_Sall__,X_weight_r_);
parameter = struct('type','parameter');
[ ...
 parameter ...
,X_wSM___ ...
,~ ...
,~ ...
,~ ...
,~ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,1 ...
,UX_M_k_q_nMA_wq_ ...
,1 ...
,tmp_n_M ...
,svd_VUXM_lwnM____ ...
,1 ...
);
X_wn__(:,1+nUX_rank) = squeeze(X_wSM___);
%%%%;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;

gamma_ = linspace(0,2*pi,1+n_w_max); gamma_ = gamma_(1:end-1); gamma_ = gamma_(:);
n_w_max_ups = n_w_max*4; gamma_ups_ = linspace(0,2*pi,1+n_w_max_ups); gamma_ups_ = gamma_ups_(1:end-1); gamma_ups_ = gamma_ups_(:);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_X_wSM_FIGF_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,512*2,768]);
fontsize_use = 16;
intensity_gamma = 0.10;
intensity_scale_ = (1-intensity_gamma) + intensity_gamma*transpose(linspace(-1,1,64)).^2;
intensity_scale_ = intensity_scale_./max(1e-12,intensity_scale_(end));
c_beach__ = bsxfun(@times,colormap_beach(),intensity_scale_); n_c_beach = size(c_beach__,1);
subplot(1,1,1);
hold on;
for nUX_rank=0:n_UX_rank-1;
nc_beach = max(0,min(n_c_beach-1,floor(n_c_beach*nUX_rank/n_UX_rank)));
%%%%;
%plot(gamma_,X_wn__(:,1+nUX_rank),'-','Color',c_beach__(1+nc_beach,:),'LineWidth',2);
n_w_max_ups = n_w_max*64; gamma_ups_ = linspace(0,2*pi,1+n_w_max_ups); gamma_ups_ = gamma_ups_(1:end-1); gamma_ups_ = gamma_ups_(:);
tmp_X_w_ = X_wn__(:,1+nUX_rank);
tmp_X_q_ = fft(tmp_X_w_)/n_w_max;
tmp_X_uq_ = zeros(n_w_max_ups,1);
tmp_n_w_max_2 = floor(min(n_w_max,n_w_max_ups)/2);
tmp_X_uq_(1:tmp_n_w_max_2) = tmp_X_q_(1:tmp_n_w_max_2);
tmp_X_uq_(end-tmp_n_w_max_2+1:end) = tmp_X_q_(end-tmp_n_w_max_2+1:end);
tmp_X_uw_ = ifft(tmp_X_uq_)*n_w_max_ups; tmp_X_uw_ = real(tmp_X_uw_);
plot(gamma_ups_,tmp_X_uw_,'-','Color',c_beach__(1+nc_beach,:),'LineWidth',2);
%%%%;
end;%for nUX_rank=0:n_UX_rank-1;
plot(gamma_ups_,zeros(size(gamma_ups_)),'k-','LineWidth',0.5);
hold off;
xlim([0,2*pi]); xlabel('\gamma'); set(gca,'XTick',[0,2*pi],'XTickLabel',{'0','2\pi'});
ylabel('$\cal{X}(\gamma)$','Interpreter','latex');
set(gca,'YTick',[]);
%title(sprintf('$\\langle \\hat{A}_{%d} , R_{\\gamma}\\circ\\hat{A}_{%d} \\rangle$ for the first %d ranks',1+nMA,1+nMB,n_svd),'Interpreter','latex');
title(sprintf('$\\langle \\hat{A}_{%d} , R_{\\gamma}\\circ\\hat{A}_{%d} \\rangle$',1+nMA,1+nMB),'Interpreter','latex');
str_legend_ = cell(n_svd,1);
for nsvd=0:n_svd-1;
str_legend_{1+nsvd} = sprintf('rank %.2d',1+nsvd);
end;%for nsvd=0:n_svd-1;
legend(str_legend_,'Location','NorthWest','Orientation','vertical','FontSize',12);
set(gca,'FontSize',fontsize_use);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

flag_check=0;
if flag_check;
%%%%%%%%;
% check to see if ampmh_X_wSM___8 and ampmh_X_single_cluster__10 do the same thing. ;
%%%%%%%%;
tmp8_X_n_ = zeros(n_UX_rank,1);
tmp9_X_n_ = zeros(n_UX_rank,1);
for nUX_rank=0:n_UX_rank-1;
pm_n_UX_rank = 1+nUX_rank;
pm_n_w_ = n_w_max*ones(pm_n_UX_rank,1);
pm_n_w_sum = sum(pm_n_w_);
UX_M_k_q_nMA_wq_ = reshape(UX_M_k_q_nMA_wq__(:,1:pm_n_UX_rank),[pm_n_w_sum,1]);
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M,M_k_q__(:,1+nMB),pm_n_UX_rank,UX_2d_Sall__,X_weight_r_);
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z=1;
[ ...
 parameter ...
,tmp8_X_SM__ ...
] = ...
ampmh_X_wSM___8( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,1 ...
,UX_M_k_q_nMA_wq_ ...
,1 ...
,tmp_n_M ...
,svd_VUXM_lwnM____ ...
,1 ...
);
tmp8_X_n_(1+nUX_rank) = squeeze(tmp8_X_SM__);
%%%%;
parameter = struct('type','parameter');
parameter.flag_optimize_over_gamma_z=1;
[ ...
 parameter ...
,tmp9_X_SM__ ...
] = ...
ampmh_X_single_cluster_SM__10( ...
 parameter ...
,FTK ...
,n_w_max ...
,pm_n_UX_rank ...
,1 ...
,UX_M_k_q_nMA_wq_ ...
,1 ...
,tmp_n_M ...
,svd_VUXM_lwnM____ ...
,1 ...
);
tmp9_X_n_(1+nUX_rank) = squeeze(tmp9_X_SM__);
%%%%;
end;%for nUX_rank=0:n_UX_rank-1;
disp(sprintf(' %% tmp8_X_n_ vs tmp9_X_n_: %0.16f',fnorm(tmp8_X_n_-tmp9_X_n_)/fnorm(tmp8_X_n_)));
end;%if flag_check;

%%%%%%%%;
% Now display some more images and templates. ;
%%%%%%%%;
n_nM_sub = 4;
nM_sub_ = max(0,min(n_M-1,floor(n_M*linspace(0,1,n_nM_sub))));
nM_sub_ = [0 , 221 , 750 , 900];
tmp_n_x = round(k_p_r_max*2*pi*sqrt(2));
tmp_index_ = n_x_u/2 + [-tmp_n_x:+tmp_n_x];
M_x_c_sub___ = zeros(numel(tmp_index_),numel(tmp_index_),n_nM_sub);
S_x_c_sub___ = zeros(numel(tmp_index_),numel(tmp_index_),n_nM_sub);
for nnM_sub=0:n_nM_sub-1;
nM = nM_sub_(1+nnM_sub);
M_x_c_sub___(:,:,1+nnM_sub) = M_x_c___(1+tmp_index_,1+tmp_index_,1+nM);
nS = index_nS_from_nM_(1+nM);
T_k_p_ = S_k_p__(:,1+nS);
tmp_euler_gamma_z = +euler_gamma_z_true_(1+nM);
T_k_p_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,T_k_p_,+tmp_euler_gamma_z);
T_k_p_ = sparse(1:n_w_sum,1:n_w_sum,CTF_k_p__(:,1+CTF_index_(1+nM)),n_w_sum,n_w_sum)*T_k_p_;
tmp_S_x_c__ = real(interp_k_p_to_x_c_xxnufft(n_x_u,diameter_x_c,n_x_u,diameter_x_c,n_k_p_r,k_p_r_,n_w_,T_k_p_));
tmp_S_x_c__ = tmp_S_x_c__(1+tmp_index_,1+tmp_index_);
S_x_c_sub___(:,:,1+nnM_sub) = tmp_S_x_c__;
end;%for nnM_sub=0:n_nM_sub-1;
%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_M_and_S_FIGG_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;clf;figbig;fig80s;
fontsize_use = 18;
p_row = 2; p_col = n_nM_sub; ns=0;
%%%%;
for nnM_sub=0:n_nM_sub-1;
nM = nM_sub_(1+nnM_sub);
nS = index_nS_from_nM_(1+nM);
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_M_x_c__ = M_x_c_sub___(:,:,1+nnM_sub);
M_x_c_lim_ = mean(tmp_M_x_c__,'all') + std(tmp_M_x_c__,1,'all')*2.5*[-1,+1];
imagesc(tmp_M_x_c__,M_x_c_lim_);
title(sprintf('$A_{%d}(\\vec{x})$',1+nM),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
end;%for nnM_sub=0:n_nM_sub-1;
%%%%;
for nnM_sub=0:n_nM_sub-1;
nM = nM_sub_(1+nnM_sub);
nS = index_nS_from_nM_(1+nM);
subplot(p_row,p_col,1+ns); ns=ns+1;
tmp_S_x_c__ = S_x_c_sub___(:,:,1+nnM_sub);
S_x_c_lim_ = mean(tmp_S_x_c__,'all') + std(tmp_S_x_c__,1,'all')*2.5*[-1,+1];
imagesc(tmp_S_x_c__,S_x_c_lim_);
title(sprintf('$A^{\\mbox{signal}}_{%d}(\\vec{x})$',1+nM),'Interpreter','Latex');
axis image; axisnotick; set(gca,'FontSize',fontsize_use);
end;%for nnM_sub=0:n_nM_sub-1;
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now perform stripped down innerproduct. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%;
% First cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
tmp_t = tic();
CTF_k_p_r_kC__ = CTF_k_p_r__;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.tolerance_cluster = 5e-2; %<-- value chosen for figure. ;
[ ...
 parameter ...
,index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% index_ncluster_from_nCTF_: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
n_cluster = 1+max(index_ncluster_from_nCTF_);
index_ncluster_from_nM_ = index_ncluster_from_nCTF_(1+index_nCTF_from_nM_);
index_nM_from_ncluster__ = cell(n_cluster,1);
n_index_nM_from_ncluster_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster__{1+ncluster} = efind(index_ncluster_from_nM_==ncluster);
n_index_nM_from_ncluster_(1+ncluster) = numel(index_nM_from_ncluster__{1+ncluster});
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% n_cluster: %0.3fs',tmp_t)); end;
%%%%%%%%;
% construct CTF of same size as images. ;
%%%%%%%%
CTF_k_p_wkC__ = zeros(n_w_sum,n_CTF);
for nCTF=0:n_CTF-1;
for nk_p_r=0:n_k_p_r-1;
n_w = n_w_(1+nk_p_r);
tmp_index_ = n_w_csum_(1+nk_p_r) + [0:n_w-1];
CTF_k_p_wkC__(1+tmp_index_,1+nCTF) = CTF_k_p_r_kC__(1+nk_p_r,1+nCTF);
end;%for nk_p_r=0:n_k_p_r-1;
end;%for nCTF=0:n_CTF-1;
%%%%%%%%;
% then calculate average CTFs for each cluster. ;
%%%%%%%%;
CTF_k_p_r_xavg_kc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
CTF_k_p_r_xavg_kc__(:,1+ncluster) = CTF_k_p_r_xavg_k_;
end;%for ncluster=0:n_cluster-1;
%%%%%%%%;
% prepare images. ;
%%%%%%%%;
tmp_t = tic();
M_k_p_wkM__ = M_k_p__;
M_k_q_wkM__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_wn_ = M_k_p_wkM__(:,1+nM);
M_k_q_wn_ = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_wn_ ...
);
M_k_q_wkM__(:,1+nM) = M_k_q_wn_;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q_wkM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
M_k_q_cwkM___ = cell(n_cluster,1);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
M_k_q_cwkM___{1+ncluster} = M_k_q_wkM__(:,1+index_nM_from_ncluster_);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q_cwkM___: %0.3fs',tmp_t)); end;
%%%%%%%%;
% prepare templates. ;
%%%%%%%%;
S_k_p_wkS__ = S_k_p__;
S_k_q_wkS__ = S_k_q__;
%%%%%%%%;
% Set up FTK. ;
%%%%%%%%;
delta_r_max = 0*delta_sigma; svd_eps = tolerance_master; n_delta_v_requested = 32;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
%%%%%%%%;

%%%%%%%%;
% Now calculate kernel for each cluster. ;
%%%%%%%%;
tot_t_kern_a = 0;
tmp_t = tic();
delta_sigma_base = 0;
a_k_Y_base_yk_ = a_k_Y_quad_;
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
CTF_k_p_r_xavg_kk__ = CTF_k_p_r_xavg_k_*transpose(CTF_k_p_r_xavg_k_);
[ ...
 X_2d_xavg_dx_kk__ ...
,X_2d_xavg_dx_weight_r_ ...
] = ...
principled_marching_cost_matrix_6( ...
 n_k_p_r ...
,k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_ ...
,[] ...
,[] ...
,a_k_Y_base_yk_ ...
,CTF_k_p_r_xavg_kk__ ...
,delta_sigma_base ...
);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
X_a_kkc___ = X_2d_xavg_dx_kkc___;
X_weight_rc__ = X_2d_xavg_dx_weight_rc__;
clear X_2d_xavg_dx_kkc__ X_2d_xavg_dx_weight_rc__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_a_kkc___: %0.3fs',tmp_t)); end;
tot_t_kern_a = tmp_t;
%%%%%%%%;
tmp_X_weight_r_ = sqrt(weight_2d_k_p_r_);
tot_t_kern_M = 0;
tmp_t0_kern_M = 0;
tmp_t1_kern_M = 0;
tmp_t2_kern_M = 0;
tmp_t3_kern_M = 0;
tmp_t = tic();
for ncluster=0:n_cluster-1;
%index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
%assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
tmp_n_M = n_index_nM_from_ncluster;
%tmp_CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
%%%%;
tmp_t0 = tic();
tmp_M_k_q_wkM___ = reshape(M_k_q_cwkM___{1+ncluster},[n_w_max,n_k_p_r,tmp_n_M]);
tmp_M_k_q_kwM___ = permute(tmp_M_k_q_wkM___,[2,1,3]);
tmp_t0 = toc(tmp_t0); tmp_t0_kern_M = tmp_t0_kern_M + tmp_t0;
%%%%;
tmp_t1 = tic();
tmp_X_00__ = zeros(n_k_p_r,n_k_p_r);
for tmp_nM=0:tmp_n_M-1;
tmp_M_k_q_kw__ = tmp_M_k_q_kwM___(:,:,1+tmp_nM);
tmp_X_00_single__ = conj(tmp_M_k_q_kw__)*transpose(tmp_M_k_q_kw__);
tmp_X_00__ = tmp_X_00__ + tmp_X_00_single__;
end;%for tmp_nM=0:tmp_n_M-1;
tmp_X_00__ = (2*pi)^2 * tmp_X_00__ / tmp_n_M ;
tmp_t1 = toc(tmp_t1); tmp_t1_kern_M = tmp_t1_kern_M + tmp_t1;
%%%%;
tmp_t2 = tic();
tmp_X_01__ = zeros(n_k_p_r,n_k_p_r);
tmp_M_k_q_k0M_ = sum(reshape(tmp_M_k_q_kwM___(:,1+0,:),[n_k_p_r,tmp_n_M]),2);
tmp_X_01__ = (2*pi)^2 * conj(tmp_M_k_q_k0M_) * transpose(tmp_M_k_q_k0M_) / tmp_n_M^2 ;
tmp_t2 = toc(tmp_t2); tmp_t2_kern_M = tmp_t2_kern_M + tmp_t2;
%%%%;
tmp_t3 = tic();
tmp_X__ = diag(tmp_X_weight_r_) * (2*real(tmp_X_00__) - 2*real(tmp_X_01__)) * diag(tmp_X_weight_r_) ;
tmp_t3 = toc(tmp_t3); tmp_t3_kern_M = tmp_t3_kern_M + tmp_t3;
%%%%;
X_M_kkc___(:,:,1+ncluster) = tmp_X__;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_M_kkc___: %0.3fs',tmp_t)); end;
tot_t_kern_M = tmp_t;
clear tmp_M* tmp_X*;
if (verbose>1); disp(sprintf(' %% X_M_kkc___: %0.3fs %0.3fs %0.3fs %0.3fs --> %0.3fs',tmp_t0_kern_M,tmp_t1_kern_M,tmp_t2_kern_M,tmp_t3_kern_M,tot_t_kern_M)); end;
%%%%%%%%;
S_k_q_wkS___ = reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]);
tmp_X_weight_r_ = sqrt(weight_2d_k_p_r_);
tot_t_kern_S = 0;
tot_t0_kern_S = 0;
tot_t1_kern_S = 0;
tot_t2_kern_S = 0;
tot_t3_kern_S = 0;
tmp_t = tic();
%%%%%%%%;
% Note that CTF is included in X_S_kkc___, but not X_S_kk__. ;
%%%%%%%%;
tmp_t1 = tic();
tmp_X_00__ = zeros(n_k_p_r,n_k_p_r);
for nS=0:n_S-1;
tmp_S_k_q_wk__ = S_k_q_wkS___(:,:,1+nS);
tmp_X_00_single__ = ctranspose(tmp_S_k_q_wk__)*tmp_S_k_q_wk__;
tmp_X_00__ = tmp_X_00__ + tmp_X_00_single__;
end;%for nS=0:n_S-1;
tmp_X_00__ = (2*pi)^2 * tmp_X_00__ / n_S ;
tmp_t1 = toc(tmp_t1); tot_t1_kern_S = tot_t1_kern_S + tmp_t1;
%%%%;
tmp_t2 = tic();
tmp_X_01__ = zeros(n_k_p_r,n_k_p_r);
tmp_S_k_q_k0S_ = sum(reshape(S_k_q_wkS___(1+0,:,:),[n_k_p_r,n_S]),2);
tmp_X_01__ = (2*pi)^2 * conj(tmp_S_k_q_k0S_) * transpose(tmp_S_k_q_k0S_) / n_S^2 ;
tmp_t2 = toc(tmp_t2); tot_t2_kern_S = tot_t2_kern_S + tmp_t2;
%%%%;
tmp_t3 = tic();
tmp_X__ = diag(tmp_X_weight_r_) * (2*real(tmp_X_00__) - 2*real(tmp_X_01__)) * diag(tmp_X_weight_r_) ;
tmp_t3 = toc(tmp_t3); tot_t3_kern_S = tot_t3_kern_S + tmp_t3;
%%%%;
X_S_kk__ = tmp_X__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_S_kk__: %0.3fs',tmp_t)); end;
clear tmp_S* tmp_X*;
%%%%%%%%;
tmp_t = tic;
for ncluster=0:n_cluster-1;
tmp_CTF_k_p_r_xavg_k_ = CTF_k_p_r_xavg_kc__(:,1+ncluster);
X_S_kkc___(:,:,1+ncluster) = diag(CTF_k_p_r_xavg_k_)*X_S_kk__*diag(CTF_k_p_r_xavg_k_);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_S_kkc___: %0.3fs',tmp_t)); end;
tot_t_kern_S = tmp_t;
if (verbose>1); disp(sprintf(' %% X_S_kkc___: %0.3fs %0.3fs %0.3fs %0.3fs --> %0.3fs',tot_t0_kern_S,tot_t1_kern_S,tot_t2_kern_S,tot_t3_kern_S,tot_t_kern_S)); end;
%%%%%%%%;
X_S_avg_kk__ = mean(X_S_kkc___(:,:,1+index_ncluster_from_nM_),3);
%%%%%%%%;
% Now use X_S_avg_kk__ to find single set of principal-modes for all clusters. ;
%%%%%%%%;
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
tmp_t = tic();
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_S_avg_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
UX_S_avg_kn__ = tmp_UX__;
SX_S_avg_k_ = diag(tmp_SX__);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_S_avg_kn__: %0.3fs',tmp_t)); end;
tot_t_UX_s_avg = tmp_t;
%%%%%%%%;
% Now calculate alignment, aggregating across clusters. ;
%%%%%%%%;
fname_mat = sprintf('%s_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
weight_2d_k_all_scaled_ = weight_2d_k_all_*4*pi^2;
CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_),2);
S_k_q_wkS___ = reshape(S_k_q_wkS__,[n_w_max,n_k_p_r,n_S]);
CTF_S_k_q_wkS___ = bsxfun(@times,S_k_q_wkS___,reshape(CTF_k_p_r_xavg_k_,[1,n_k_p_r,1]));
CTF_S_k_q_kwS___ = permute(CTF_S_k_q_wkS___,[2,1,3]);
M_k_q_wkM___ = reshape(M_k_q_wkM__,[n_w_max,n_k_p_r,n_M]);
M_k_q_kwM___ = permute(M_k_q_wkM___,[2,1,3]);
%%%%%%%%;
%set temporary array dimensions. ;
%%%%%%%%;
tmp_n_S = n_S;
tmp_n_M = 1;
tmp_n_S = min(tmp_n_S,n_S);
tmp_n_M = min(tmp_n_M,n_M);
%%%%%%%%;
tot_t_fill_ = zeros(n_UX_rank,1);
tot_t_mult_ = zeros(n_UX_rank,1);
tot_t_ifft_ = zeros(n_UX_rank,1);
tot_o_fill_ = zeros(n_UX_rank,1);
tot_o_mult_ = zeros(n_UX_rank,1);
tot_o_ifft_ = zeros(n_UX_rank,1);
X_wSMe___ = zeros(n_w_max,n_S,n_M);
X_wSM0___ = zeros(n_w_max,n_S,n_M);
l2er_X_X_n_ = zeros(n_UX_rank,1);
corr_X_X_n_ = zeros(n_UX_rank,1);
prct_X_X_n_ = zeros(n_UX_rank,1);
for nUX_rank=[n_UX_rank-1:-1:0];%for nUX_rank=n_UX_rank-1:-1:0
pm_n_UX_rank = 1+nUX_rank;
if (verbose); disp(sprintf(' %% pm_n_UX_rank %d/%d',pm_n_UX_rank,n_UX_rank)); end;
tmp_UX_weight_kn__ = diag(X_weight_r_)*UX_S_avg_kn__(:,1:pm_n_UX_rank);
tmp_SX_k_ = SX_S_avg_k_(1:pm_n_UX_rank);
%%%%%%%%;
%%%%;
n_M_batch = ceil(n_M/tmp_n_M); n_S_batch = ceil(n_S/tmp_n_S);
for nS_batch=0:n_S_batch-1;
tmp_index_nS_ = nS_batch*tmp_n_S + [0:tmp_n_S-1]; tmp_index_nS_ = intersect(tmp_index_nS_,[0:n_S-1]);
tmp_CTF_S_k_q_kwS___ = CTF_S_k_q_kwS___(:,:,1+tmp_index_nS_);
tmp_t = tic();
UX_weight_CTF_S_k_q_nwS___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(tmp_CTF_S_k_q_kwS___,[n_k_p_r,n_w_max*tmp_n_S]),[pm_n_UX_rank,n_w_max,tmp_n_S]);
UX_weight_conj_CTF_S_k_q_Snw___ = conj(permute(UX_weight_CTF_S_k_q_nwS___,[3,1,2]));
tmp_o = pm_n_UX_rank*n_k_p_r*n_w_max*(tmp_n_M + tmp_n_S);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% fill: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_fill_(1+nUX_rank) = tot_o_fill_(1+nUX_rank) + tmp_o;
tot_t_fill_(1+nUX_rank) = tot_t_fill_(1+nUX_rank) + tmp_t;
for nM_batch=0:n_M_batch-1;
tmp_index_nM_ = nM_batch*tmp_n_M + [0:tmp_n_M-1]; tmp_index_nM_ = intersect(tmp_index_nM_,[0:n_M-1]);
tmp_M_k_q_kwM___ = M_k_q_kwM___(:,:,1+tmp_index_nM_);
%%%%;
tmp_t = tic();
UX_weight_M_k_q_nwM___ = reshape(transpose(tmp_UX_weight_kn__)*reshape(tmp_M_k_q_kwM___,[n_k_p_r,n_w_max*tmp_n_M]),[pm_n_UX_rank,n_w_max,tmp_n_M]);
UX_weight_M_k_q_nMw___ = permute(UX_weight_M_k_q_nwM___,[1,3,2]);
tmp_o = pm_n_UX_rank*n_k_p_r*n_w_max*(tmp_n_M + tmp_n_S);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% fill: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_fill_(1+nUX_rank) = tot_o_fill_(1+nUX_rank) + tmp_o;
tot_t_fill_(1+nUX_rank) = tot_t_fill_(1+nUX_rank) + tmp_t;
%%%%;
tmp_t = tic();
conj_S_CTF_weight_weight_M_k_q_SMw___ = zeros(tmp_n_S,tmp_n_M,n_w_max);
for nw=0:n_w_max-1;
conj_S_CTF_weight_weight_M_k_q_SMw___(:,:,1+nw) = UX_weight_conj_CTF_S_k_q_Snw___(:,:,1+nw) * UX_weight_M_k_q_nMw___(:,:,1+nw);
end;%for nw=0:n_w_max-1;
tmp_o = n_w_max * tmp_n_S * tmp_n_M * pm_n_UX_rank;
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% mult: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_mult_(1+nUX_rank) = tot_o_mult_(1+nUX_rank) + tmp_o;
tot_t_mult_(1+nUX_rank) = tot_t_mult_(1+nUX_rank) + tmp_t;
%%%%;
tmp_t = tic();
ifft_conj_S_CTF_weight_weight_M_k_q_wSM___ = ifft(permute(conj_S_CTF_weight_weight_M_k_q_SMw___,[3,1,2]),[],1);
tmp_o = tmp_n_S * tmp_n_M * n_w_max * log(n_w_max);
tmp_t = toc(tmp_t); if (verbose>2); disp(sprintf(' %% ifft: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
tot_o_ifft_(1+nUX_rank) = tot_o_ifft_(1+nUX_rank) + tmp_o;
tot_t_ifft_(1+nUX_rank) = tot_t_ifft_(1+nUX_rank) + tmp_t;
X_wSM___(:,1+tmp_index_nS_,1+tmp_index_nM_) = real(ifft_conj_S_CTF_weight_weight_M_k_q_wSM___);
%%%%;
end;end;%for nM_batch=0:n_M_batch-1; for nS_batch=0:n_S_batch-1;
if (verbose>1); disp(sprintf(' %% tot_t_fill %0.3fs --> %0.2fGH',tot_t_fill_(1+nUX_rank),tot_o_fill_(1+nUX_rank)/tot_t_fill_(1+nUX_rank)/1e9)); end;
if (verbose>1); disp(sprintf(' %% tot_t_mult %0.3fs --> %0.2fGH',tot_t_mult_(1+nUX_rank),tot_o_mult_(1+nUX_rank)/tot_t_mult_(1+nUX_rank)/1e9)); end;
if (verbose>1); disp(sprintf(' %% tot_t_ifft %0.3fs --> %0.2fGH',tot_t_ifft_(1+nUX_rank),tot_o_ifft_(1+nUX_rank)/tot_t_ifft_(1+nUX_rank)/1e9)); end;
if (nUX_rank==n_UX_rank-1);
X_wSMe___ = X_wSM___;
end;%if (nUX_rank==n_UX_rank-1);
X_wSM0___ = X_wSM___;
l2er_X_X = fnorm(X_wSMe___ - X_wSM0___) / fnorm(X_wSMe___);
l2er_X_X_n_(1+nUX_rank) = l2er_X_X;
corr_X_X = corr(X_wSMe___(:),X_wSM0___(:));
corr_X_X_n_(1+nUX_rank) = corr_X_X;
[~,tmp_nw_SM_] = max(X_wSM0___,[],1); tmp_nw_SM_ = tmp_nw_SM_-1;
tmp_X_SM___ = zeros(1,n_S,n_M);
for nM=0:n_M-1; for nS=0:n_S-1;
tmp_nw = tmp_nw_SM_(1,1+nS,1+nM);
tmp_X = X_wSMe___(1+tmp_nw,1+nS,1+nM);
tmp_X_SM___(1,1+nS,1+nM) = tmp_X;
end;end;%for nM=0:n_M-1; for nS=0:n_S-1;
tmp_X_SM___ = repmat(tmp_X_SM___,[n_w_max,1,1]);
tmp_X_SM___ = X_wSMe___ > tmp_X_SM___;
tmp_p__ = sum(tmp_X_SM___,1)/n_w_max;
tmp_p = mean(tmp_p__,'all');
prct_X_X_n_(1+nUX_rank) = tmp_p;
%%%%%%%%;
tmp_t = tic();
gamma_z_ = 2*pi*[0:n_w_max-1]/n_w_max;
nS=3; nM=1;
X1_w_ = X_wSM___(:,1+nS,1+nM);
X0_w_ = zeros(n_w_max,1);
for nw=0:n_w_max-1;
gamma_z = 2*pi*nw/n_w_max;
S_k_p_wk_ = S_k_p_wkS__(:,1+nS); M_k_p_wk_ = M_k_p_wkM__(:,1+nM);
T_k_p_wk_ = reshape(reshape(S_k_p_wk_,[n_w_max,n_k_p_r])*diag(CTF_k_p_r_xavg_k_),[n_w_max*n_k_p_r,1]);
N_k_p_wk_ = rotate_p_to_p_fftw(n_k_p_r,n_w_,n_w_sum,M_k_p_wk_,-gamma_z);
X0_w_(1+nw) = real(sum(conj(T_k_p_wk_).*N_k_p_wk_.*weight_2d_k_all_scaled_));
end;%for nw=0:n_w_max-1;
tmp_o = n_w_max * n_k_p_r * n_w_max;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% test: %0.3fs --> %0.2fGHz',tmp_t,tmp_o/tmp_t/1e9)); end;
if (verbose); disp(sprintf(' %% X0_w_ vs X1_w_: %0.16f',fnorm(X0_w_-X1_w_)/fnorm(X0_w_))); end;
end;%for nUX_rank=0:n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'weight_2d_k_all_scaled_','CTF_k_p_r_xavg_k_' ...
     ,'prct_X_X_n_','corr_X_X_n_','l2er_X_X_n_' ...
     ,'tot_t_fill_','tot_t_mult_','tot_t_ifft_','tot_o_fill_','tot_o_mult_','tot_o_ifft_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_X_2d_Semp_d1_FIGH__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figmed;
fontsize_use = 12;
%%%%;
subplot(1,3,1);
plot(1:n_UX_rank,l2er_X_X_n_,'ko-','MarkerFaceColor','k');
xlim([0,n_k_p_r]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylabel('error'); grid on;
legend({'frob'}); set(gca,'FontSize',fontsize_use);
title('relative error');
%%%%;
subplot(1,3,2);
hold on;
plot(1:n_UX_rank,0+corr_X_X_n_,'ko-','MarkerFaceColor','r');
plot(1:n_UX_rank,1-prct_X_X_n_,'ko-','MarkerFaceColor','c');
hold off;
xlim([0,n_k_p_r]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'}); set(gca,'FontSize',fontsize_use);
title('correlation');
%%%%;
subplot(1,3,3);
tot_t_all_ = tot_t0_kern_S + tot_t1_kern_S + tot_t2_kern_S + tot_t3_kern_S + tot_t_UX_s_avg + tot_t_kern_S + tot_t_fill_ + tot_t_mult_ + tot_t_ifft_ ;
hold on;
plot(1:n_UX_rank,max(tot_t_mult_+tot_t_ifft_)./tot_t_all_,'ko-','MarkerFaceColor','r');
plot(1:n_UX_rank,max(tot_o_mult_)./tot_o_mult_,'ko-','MarkerFaceColor','c');
plot(1:n_UX_rank,ones(1,n_UX_rank),'k-','LineWidth',0.5);
hold off;
xlim([0,n_k_p_r]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0,16]); ylabel('factor'); set(gca,'YTick',0:1:16); grid on;
legend({'time','ops'}); set(gca,'FontSize',fontsize_use);
title('speedup');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_X_2d_Semp_d1_FIGHAB__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;figmed;
markersize_use = 8; fontsize_use = 10;
%%%%;
subplot(1,2,1);
plot(1:n_UX_rank,l2er_X_X_n_,'ko-','MarkerFaceColor','k','MarkerSize',markersize_use);
xlim([0,n_k_p_r]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylabel('error'); grid on;
legend({'frob'});
title('relative error');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,2,2);
hold on;
plot(1:n_UX_rank,0+corr_X_X_n_,'ko-','MarkerFaceColor','r','MarkerSize',markersize_use);
plot(1:n_UX_rank,1-prct_X_X_n_,'ko-','MarkerFaceColor','c','MarkerSize',markersize_use);
hold off;
xlim([0,n_k_p_r]); set(gca,'XTick',0:7:49); xlabel('rank $H$','Interpreter','latex');
ylim([0.5,1]); ylabel('value'); set(gca,'YTick',0.5:0.05:1.0); grid on;
legend({'corr','prct'});
title('correlation');
set(gca,'FontSize',fontsize_use);
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
%close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;



