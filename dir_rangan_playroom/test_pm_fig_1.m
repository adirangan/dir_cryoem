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
p_row = 2; p_col = 4; ns=0;
%%%%%%%%;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
M_x_c_lim_ = mean(M_x_c_nMA_,'all') + std(M_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMA_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
S_x_c_lim_ = mean(S_x_c_nMA_,'all') + std(S_x_c_nMA_,1,'all')*2.5*[-1,+1];
imagesc(S_x_c_nMA_,S_x_c_lim_);
title(sprintf('$B_{%d} := A_{%d}^{\\mbox{signal}}(\\vec{x})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(U_k_p_nMA_),[],colormap_80s);
title(sprintf('$\\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k}) + \\hat{A}_{%d}^{\\mbox{noise}}(\\vec{k})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMA_),[],colormap_80s);
title(sprintf('$\\hat{B}_{%d} := \\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k})$',1+nMA,1+nMA),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
%%%%%%%%;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
M_x_c_lim_ = mean(M_x_c_nMB_,'all') + std(M_x_c_nMB_,1,'all')*2.5*[-1,+1];
imagesc(M_x_c_nMB_,M_x_c_lim_);
title(sprintf('$A_{%d}^{\\mbox{signal}}(\\vec{x}) + A_{%d}^{\\mbox{noise}}(\\vec{x})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
S_x_c_lim_ = mean(S_x_c_nMB_,'all') + std(S_x_c_nMB_,1,'all')*2.5*[-1,+1];
imagesc(S_x_c_nMB_,S_x_c_lim_);
title(sprintf('$B_{%d} := A_{%d}^{\\mbox{signal}}(\\vec{x})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(U_k_p_nMB_),[],colormap_80s);
title(sprintf('$\\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k}) + \\hat{A}_{%d}^{\\mbox{noise}}(\\vec{k})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
subplot(p_row,p_col,1+ns); ns=ns+1;
imagesc_p(n_k_p_r,k_p_r_,n_w_,n_w_sum,real(T_k_p_nMB_),[],colormap_80s);
title(sprintf('$\\hat{B}_{%d} := \\hat{A}_{%d}^{\\mbox{signal}}(\\vec{k})$',1+nMB,1+nMB),'Interpreter','Latex');
axis image; axisnotick;
%%%%;
%%%%%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
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
figure(1+nf);nf=nf+1;figmed;
title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C^{\\mbox{signal}}(\\vec{u}) = %0.3f$ , $C^{\\mbox{signal}}(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
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
       },'Interpreter','latex');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
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
figure(1+nf);nf=nf+1;figmed;
title(sprintf('Two image-rings $k_{1}=%0.2f$ and $k_{2}=%0.2f$: $\\vec{u} = [ %0.2f , %0.2f ]^{T}$, $C^{\\mbox{signal}}(\\vec{u}) = %0.3f$ , $C^{\\mbox{signal}}(\\vec{v}) = %0.3f$',2*pi*k_p_r_(1+[nk_p_r_0,nk_p_r_1]),UX_2d_S2__(1+0,1),UX_2d_S2__(1+1,1),SX_2d_S2_(1+0),SX_2d_S2_(1+1)),'Interpreter','latex');
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
       },'Interpreter','latex');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
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
figure(1+nf);nf=nf+1;figsml;
n_svd = 8;
subplot(1,1,1);
plot(SX_2d_Sall_,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
set(gca,'XTick',1:n_svd); xlim([0.5,n_svd+0.5]);
xlabel('rank');
ylabel('singular-value','Interpreter','latex');
grid on;
title(sprintf('$C(\\vec{u})$ for the first %d ranks',n_svd),'Interpreter','latex');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
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
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
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
title(sprintf('$\\langle %s , CTF\\odot\\left[ R_{\\gamma}\\circ\\hat{S}\\right] \\rangle$ for the first %d ranks',tmp_str,n_svd),'Interpreter','latex');
str_legend_ = cell(n_svd,1);
for nsvd=0:n_svd-1;
str_legend_{1+nsvd} = sprintf('rank %.2d',1+nsvd);
end;%for nsvd=0:n_svd-1;
legend(str_legend_,'Location','NorthWest');
end;%for tmp_nM=[0,1];
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
%%%%%%%%;
fname_fig = sprintf('%s_jpg/pm_fig_image_pair_X_wSM_FIGF_',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,512,768]);
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1);
subplot(1,1,1);
hold on;
for nUX_rank=0:n_UX_rank-1;
nc_beach = max(0,min(n_c_beach-1,floor(n_c_beach*nUX_rank/n_UX_rank)));
plot(gamma_,X_wn__(:,1+nUX_rank),'-','Color',c_beach__(1+nc_beach,:),'LineWidth',2);
end;%for nUX_rank=0:n_UX_rank-1;
plot(gamma_,zeros(size(gamma_)),'k-','LineWidth',0.5);
hold off;
xlim([0,2*pi]); xlabel('\gamma'); set(gca,'XTick',[0,2*pi],'XTickLabel',{'0','2\pi'});
ylabel('$\cal{X}(\gamma)$','Interpreter','latex');
set(gca,'YTick',[]);
title(sprintf('$\\langle \\hat{A} , R_{\\gamma}\\circ\\hat{B} \\rangle$ for the first %d ranks',n_svd),'Interpreter','latex');
str_legend_ = cell(n_svd,1);
for nsvd=0:n_svd-1;
str_legend_{1+nsvd} = sprintf('rank %.2d',1+nsvd);
end;%for nsvd=0:n_svd-1;
legend(str_legend_,'Location','NorthWest');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
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
axis image; axisnotick;
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
axis image; axisnotick;
end;%for nnM_sub=0:n_nM_sub-1;
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

fname_mat = sprintf('%s_mat/pm_fig_X_2d_Semp_d1_FIGH__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
% Now calculate principal-modes across all images. ;
%%%%%%%%;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+index_nCTF_from_nM_(1:n_M)),n_k_p_r);
SCTF_c_ = diag(SCTF_c__);
n_CTF_rank = min(efind(SCTF_c_/max(SCTF_c_)<1e-2));
n_CTF_rank = 1;
[UCTF_kc__,SCTF_c__,VCTF_Mc__] = svds(CTF_k_p_r__(:,1+index_nCTF_from_nM_(1:n_M)),n_CTF_rank);
VSCTF_Mc__ = VCTF_Mc__*SCTF_c__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% VSCTF_Mc__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_ = M_k_p__(:,1+nM);
M_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;
%%%%%%%%;
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
[X_2d_Semp_d1__,X_2d_Semp_d1_weight_r_] = principled_marching_empirical_cost_matrix_0(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_M,T_k_p__);
[UX_2d_Semp_d1__,SX_2d_Semp_d1__,VX_2d_Semp_d1__] = svds(X_2d_Semp_d1__,n_UX_rank); SX_2d_Semp_d1_ = diag(SX_2d_Semp_d1__);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_2d_Semp_d1__: %0.3fs',tmp_t)); end;
%%%%%%%%;
UX__ = UX_2d_Semp_d1__;
X_weight_r_ = X_2d_Semp_d1_weight_r_;
%%%%%%%%;
pm_n_UX_rank = 16;
pm_n_k_p_r = pm_n_UX_rank;
pm_l_max_ = l_max_max*ones(pm_n_k_p_r,1);
pm_n_w_ = n_w_max*ones(pm_n_k_p_r,1);
pm_k_p_r_ = ones(pm_n_k_p_r,1);
pm_weight_k_p_r_ = ones(pm_n_k_p_r,1);
%%%%%%%%;
tmp_t = tic();
svd_VUXM_lwnM____ = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,n_M,M_k_q__,pm_n_UX_rank,UX__,X_weight_r_);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
tmp_t = tic();
UX_M_l2_dM__ = ampmh_UX_M_l2_dM__1(FTK,n_w_,n_M,pm_n_UX_rank,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
%%%%%%%%;
% Now form principal-images. ;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank,n_M,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
%%%%%%%%;
% use current euler-angles and displacements to solve for current model. ;
%%%%%%%%;
tmp_t = tic();
[ ...
 parameter ...
,a_UCTF_UX_Y_ync__ ... 
] = ...
a_UCTF_UX_Y_wrap_ync__0( ...
 parameter ...
,pm_n_k_p_r ...
,pm_l_max_ ...
,pm_n_w_ ...
,n_M ...
,reshape(UX_M_k_p_wnM___,[n_w_max*pm_n_k_p_r,n_M]) ...
,n_CTF_rank ...
,VSCTF_Mc__ ...
,euler_polar_a_true_ ...
,euler_azimu_b_true_...
,euler_gamma_z_true_...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% a_UCTF_UX_Y_ync__: %0.3fs',tmp_t)); end;
%%%%%%%%;
a_UCTF_UX_Y_ync__ = spharm__normalize_1(pm_n_k_p_r,pm_k_p_r_,pm_weight_k_p_r_,pm_l_max_,a_UCTF_UX_Y_ync__);
%%%%%%%%;
a_UCTF_UX_Y_ync___ = reshape(a_UCTF_UX_Y_ync__,[n_lm_max,pm_n_UX_rank,n_CTF_rank]);
%%%%%%%%;
verbose=2;
X_SMn___ = zeros(n_S,n_M,pm_n_UX_rank);
X_t_ = zeros(pm_n_UX_rank,1);
for pm_nUX_rank=0:pm_n_UX_rank-1;
%%%%%%%%;
% Use current principal-model to align principal-images. ;
% Groups principal-images by micrograph (i.e., inefficient if there are only a few images per micrograph). ;
% Calculates principal-templates associated with each micrograph. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
parameter = struct('type','parameter');
parameter.template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
parameter.flag_compress_S = 0;
tmp_a_UCTF_UX_Y_ync__ = reshape(a_UCTF_UX_Y_ync___(:,1:1+pm_nUX_rank,:),[n_lm_max*(1+pm_nUX_rank),n_CTF_rank]);
tmp_t = tic();
[ ...
 parameter ...
,n_S ...
,template_viewing_polar_a_all_ ...
,template_viewing_azimu_b_all_ ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_wrap_wrap_SM__8( ...
 parameter ...
,FTK ...
,n_w_max ...
,l_max_max ...
,1+pm_nUX_rank ...
,n_CTF_rank ...
,tmp_a_UCTF_UX_Y_ync__ ...
,n_M ...
,index_nCTF_from_nM_ ...
,VSCTF_Mc__ ...
,svd_VUXM_lwnM____(:,:,1:1+pm_nUX_rank,:) ...
,UX_M_l2_dM__ ...
,[] ...
,zeros(n_M,1) ...
,zeros(n_M,1) ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
X_SMn___(:,:,1+pm_nUX_rank) = X_SM__;
X_t_(1+pm_nUX_rank) = tmp_t;
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','n_UX_rank' ...
     ,'FTK','n_CTF_rank','SCTF_c_','VSCTF_Mc__','X_2d_Semp_d1__','X_2d_Semp_d1_weight_r_' ...
     ,'UX_2d_Semp_d1__','SX_2d_Semp_d1_' ...
     ,'pm_n_UX_rank' ...
     ,'a_UCTF_UX_Y_ync__' ...
     ,'X_SMn___','X_t_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_fig_X_2d_Semp_d1_FIGH__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));

figure(1+nf);nf=nf+1;figsml;

prctile_ = [5,15,50,85,95]; n_prctile = numel(prctile_);
cap_pn__ = zeros(n_prctile,pm_n_UX_rank);
X_SM_end__ = X_SMn___(:,:,end);
for nprctile=0:n_prctile-1;
cut_X_SM_end__ = zeros(size(X_SM_end__));
prc = prctile_(1+nprctile);
cut_X_SM_end__ = X_SM_end__ > repmat(prctile(X_SM_end__,prc),[n_S,1]);
for pm_nUX_rank=0:pm_n_UX_rank-1;
X_SM_sub__ = X_SMn___(:,:,1+pm_nUX_rank);
cut_X_SM_sub__ = zeros(size(X_SM_sub__));
cut_X_SM_sub__ = X_SM_sub__ > repmat(prctile(X_SM_sub__,prc),[n_S,1]);
cap_pn__(1+nprctile,1+pm_nUX_rank) = sum(cut_X_SM_end__.*cut_X_SM_sub__,'all');
end;%for pm_nUX_rank=0:pm_n_UX_rank-1;
end;%for nprctile=0:n_prctile-1;
cap_nrm_pn__ = cap_pn__./repmat(cap_pn__(:,end),[1,pm_n_UX_rank]);

subplot(1,3,1);
imagesc(X_SMn___(:,:,1+15)); colorbar;
subplot(1,3,2);
imagesc(X_SMn___(:,:,1+14)); colorbar;
subplot(1,3,3);
imagesc(X_SMn___(:,:,1+15) - X_SMn___(:,:,1+14)); colorbar;


%%%%;
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
% First cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
tmp_t = tic();
CTF_k_p_r_kC__ = CTF_k_p_r__;
parameter = struct('type','parameter');
parameter.tolerance_master = tolerance_master;
parameter.tolerance_cluster = 1e-3;
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
M_k_q__ = zeros(n_w_sum,n_M);
for nM=0:n_M-1;
M_k_p_ = M_k_p__(:,1+nM);
M_k_q__(:,1+nM) = ...
interp_p_to_q( ...
 n_k_p_r ...
,n_w_ ...
,n_w_sum ...
,M_k_p_ ...
);
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% M_k_q__: %0.3fs',tmp_t)); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% repeat using clusters. ;
% this time with idealized principal-modes. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
fname_mat = sprintf('%s_mat/pm_fig_X_2d_cluster_d0_FIGI__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_t = tic();
delta_sigma_base = 0;
a_k_Y_base_yk_ = a_k_Y_quad_;
X_2d_xavg_dx_kkc___ = zeros(n_k_p_r,n_k_p_r,n_cluster);
X_2d_xavg_dx_weight_rc__ = zeros(n_k_p_r,n_cluster);
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
n_index_nM_from_ncluster = n_index_nM_from_ncluster_(1+ncluster);
assert(n_index_nM_from_ncluster==numel(index_nM_from_ncluster_));
tmp_CTF_k_p_r_xavg_k_ = mean(CTF_k_p_r_kC__(:,1+index_nCTF_from_nM_(1+index_nM_from_ncluster_)),2);
tmp_CTF_k_p_r_xavg_kk__ = tmp_CTF_k_p_r_xavg_k_*transpose(tmp_CTF_k_p_r_xavg_k_);
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
,tmp_CTF_k_p_r_xavg_kk__ ...
,delta_sigma_base ...
);
X_2d_xavg_dx_kkc___(:,:,1+ncluster) = X_2d_xavg_dx_kk__;
X_2d_xavg_dx_weight_rc__(:,1+ncluster) = X_2d_xavg_dx_weight_r_;
clear X_2d_xavg_dx_kk__ X_2d_xavg_dx_weight_r_ ;
end;%for ncluster=0:n_cluster-1;
X_kkc___ = X_2d_xavg_dx_kkc___;
X_weight_rc__ = X_2d_xavg_dx_weight_rc__;
clear X_2d_xavg_dx_kkc__ X_2d_xavg_dx_weight_rc__;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_kkc___: %0.3fs',tmp_t)); end;
%%%%%%%%;
delta_r_max = 0*delta_sigma; svd_eps = tolerance_master; n_delta_v_requested = 32;
tmp_t = tic();
FTK = ampmh_FTK_1(n_k_p_r,k_p_r_,k_p_r_max,delta_r_max,svd_eps,n_delta_v_requested);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% FTK: %0.3fs',tmp_t)); end;
%%%%%%%%;
tolerance_pm_ = 0.1.^[0.5:0.5:18]; n_tolerance_pm = numel(tolerance_pm_);
X_SMp__ = zeros(n_S,n_M,n_tolerance_pm);
euler_gamma_z_SMp__ = zeros(n_S,n_M,n_tolerance_pm);
X_t_ = zeros(n_tolerance_pm,1);
X_t0_ = zeros(n_tolerance_pm,1);
X_t1_ = zeros(n_tolerance_pm,1);
X_t2_ = zeros(n_tolerance_pm,1);
X_t3_ = zeros(n_tolerance_pm,1);
pm_n_UX_rank_cp__ = zeros(n_cluster,n_tolerance_pm);
%%%%%%%%;
for ntolerance_pm=0:n_tolerance_pm-1;
tolerance_pm = tolerance_pm_(1+ntolerance_pm);
tmp_t = tic();
n_UX_rank = n_k_p_r-1; %<-- just to check dimensions. ;
UX_knc___ = zeros(n_k_p_r,n_UX_rank,n_cluster);
SX_kc__ = zeros(n_UX_rank,n_cluster);
pm_n_UX_rank_c_ = zeros(n_cluster,1);
for ncluster=0:n_cluster-1;
X_kk__ = X_kkc___(:,:,1+ncluster);
[tmp_UX__,tmp_SX__,tmp_VX__] = svds(X_kk__,n_UX_rank); tmp_SX_ = diag(tmp_SX__);
pm_n_UX_rank = max(find(tmp_SX_/max(tmp_SX_)> tolerance_pm));
if isempty(pm_n_UX_rank); pm_n_UX_rank = 1; end;
UX_knc___(:,:,1+ncluster) = tmp_UX__(:,1+[0:n_UX_rank-1]);
SX_kc__(:,1+ncluster) = tmp_SX_(1+[0:n_UX_rank-1]);
pm_n_UX_rank_c_(1+ncluster) = pm_n_UX_rank;
if (verbose>1); disp(sprintf(' %% ncluster %.2d/%.2d: pm_n_UX_rank %d/%d',ncluster,n_cluster,pm_n_UX_rank,n_UX_rank)); end;
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_knc___: %0.3fs',tmp_t)); end;
pm_n_UX_rank_cp__(:,1+ntolerance_pm) = pm_n_UX_rank_c_;
%%%%%%%%;
pm_n_UX_rank_max = max(pm_n_UX_rank_c_);
pm_n_lm_max_max = (1+l_max_max)^2;
M_k_q_wkM__ = M_k_q__;
svd_VUXM_lwnM____ = zeros(FTK.n_svd_l,n_w_max,pm_n_UX_rank_max,n_M);
%%%%%%%%;
tmp_t = tic();
tmp_M_index_ = 0:n_M-1; tmp_n_M = n_M;
for ncluster=0:n_cluster-1;
index_nM_from_ncluster_ = index_nM_from_ncluster__{1+ncluster};
pm_n_UX_rank = pm_n_UX_rank_c_(1+ncluster);
UX_kn__ = UX_knc___(:,1:pm_n_UX_rank,1+ncluster);
X_weight_r_ = X_weight_rc__(:,1+ncluster);
tmp_M_index_sub_ = intersect(tmp_M_index_,index_nM_from_ncluster_);
tmp_n_M_sub = numel(tmp_M_index_sub_);
if (tmp_n_M_sub> 0);
svd_VUXM_lwnM____(:,:,1:pm_n_UX_rank,1+tmp_M_index_sub_) = tpmh_VUXM_lwnM____3(FTK,n_k_p_r,n_w_,tmp_n_M_sub,M_k_q_wkM__(:,1+tmp_M_index_sub_),pm_n_UX_rank,UX_kn__,X_weight_r_);
end;%if (tmp_n_M_sub> 0);
end;%for ncluster=0:n_cluster-1;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_VUXM_lwnM____: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'svd_VUXM_lwnM____',tmp_t);
%%%%%%%%;
% Now calculate norms of the translated images. ;
%%%%%%%%;
UX_M_l2_dM__ = zeros(FTK.n_delta_v,n_M);
tmp_t = tic();
UX_M_l2_dM__(:,1+tmp_M_index_) = ampmh_UX_M_l2_dM__1(FTK,n_w_,tmp_n_M,pm_n_UX_rank_max,svd_VUXM_lwnM____(:,:,:,1+tmp_M_index_));
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_l2_dM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_l2_dM__1',tmp_t);
%%%%%%%%;
% Now, form principal-images;
%%%%%%%%;
tmp_t = tic();
[UX_M_k_q_wnM___,UX_M_k_p_wnM___] = ampmh_UX_M_k_p_wnM___0(FTK,n_w_,pm_n_UX_rank_max,n_M,svd_VUXM_lwnM____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% UX_M_k_q_wnM___: %0.6fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_UX_M_k_p_wnM___0',tmp_t);
%%%%%%%%;
a_k_Y_reco_yk_ = a_k_Y_quad_;
a_k_Y_reco_yk__ = zeros(n_lm_max,n_k_p_r);
for nk_p_r=0:n_k_p_r-1;
tmp_index_ = n_lm_csum_(1+nk_p_r) + (0:n_lm_(1+nk_p_r)-1);
a_k_Y_reco_yk__(1:n_lm_(1+nk_p_r),1+nk_p_r) = a_k_Y_reco_yk_(1+tmp_index_);
end;%for nk_p_r=0:n_k_p_r-1;
template_viewing_k_eq_d = 1.0/max(1e-12,k_p_r_max);
S_k_p_wkS__ = S_k_p__;
S_k_q_wkS__ = S_k_q__;
%%%%%%%%;
%%%%%%%%;
% storing all the principal-templates takes lots of memory. ;
%%%%%%%%;
CTF_UX_S_k_q_wnSc___ = [];
CTF_UX_S_l2_Sc__ = [];
%%%%%%%%;
% Use given volume to align principal-images. ;
% Groups principal-images by cluster. ;
% Calculates principal-templates associated with each cluster. ;
% Batches images into batches of size n_M_per_Mbatch (default 24). ;
% Batches templates into batches of size n_S_per_Sbatch (default 24). ;
% Only stores the optimal translation for each principal-image. ;
%%%%%%%%;
tmp_parameter = struct('type','parameter');
tmp_parameter.tolerance_cluster = parameter.tolerance_cluster;
tmp_t = tic();
[ ...
 tmp_parameter ...
,X_SM__ ...
,delta_x_SM__ ...
,delta_y_SM__ ...
,gamma_z_SM__ ...
,I_value_SM__ ...
] = ...
ampmh_X_cluster_wrap_SM__10( ...
 tmp_parameter ...
,FTK ...
,n_w_max ...
,n_k_p_r ...
,n_S ...
,S_k_q_wkS__ ...
,n_CTF ...
,CTF_k_p_r_kC__ ...
,index_nCTF_from_nM_ ...
,n_M ...
,n_cluster ...
,index_ncluster_from_nCTF_ ...
,pm_n_UX_rank_c_ ...
,UX_knc___ ...
,X_weight_rc__ ...
,svd_VUXM_lwnM____ ...
,UX_M_l2_dM__ ...
,CTF_k_p_r_xavg_kc__ ...
,CTF_UX_S_k_q_wnSc___ ...
,CTF_UX_S_l2_Sc__ ...
,index_ncluster_from_nM_ ...
,index_nM_from_ncluster__ ...
,n_index_nM_from_ncluster_ ...
);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% X_SM__: %0.3fs',tmp_t)); end;
parameter = parameter_timing_update(parameter,'ampmh_X_cluster_wrap_SM__10',tmp_t);
X_SMp___(:,:,1+ntolerance_pm) = X_SM__;
euler_gamma_z_SMp___(:,:,1+ntolerance_pm) = gamma_z_SM__;
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10'));
X_t_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_VUXM_nMwl____'));
X_t0_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_SVUXM_SMwl____'));
X_t1_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_SVUXM_lwSM____'));
X_t2_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
tmp_ij = find(strcmp(tmp_parameter.timing(:,1),'ampmh_X_single_cluster_SM__10: svd_USESVUXM_dwSM____'));
X_t3_(1+ntolerance_pm) = tmp_parameter.timing{tmp_ij,2};
end;%for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
save(fname_mat ...
     ,'delta_sigma','delta_sigma_base','n_UX_rank' ...
     ,'FTK' ...
     ,'tolerance_pm_','n_tolerance_pm' ...
     ,'pm_n_UX_rank_cp__' ...
     ,'X_SMp___','X_t_','X_t0_','X_t1_','X_t2_','X_t3_' ...
     ,'euler_gamma_z_SMp___' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

%%%%%%%%;
% Calculate error in euler_gamma_z. ;
% mean absolute error (in terms of angle). ;
% as well as percentile (over SM). ;
%%%%%%%%;
prctile_ = [90,95,99]; n_prctile = numel(prctile_);
error_euler_gamma_z_avg_p_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_prc_pp__ = zeros(n_tolerance_pm,n_prctile);
for ntolerance_pm=0:n_tolerance_pm-1;
tmp_ = abs(periodize(euler_gamma_z_SMp___(:,:,end) - euler_gamma_z_SMp___(:,:,1+ntolerance_pm),-pi,pi));
error_euler_gamma_z_avg_p_(1+ntolerance_pm) = mean(tmp_,'all');
error_euler_gamma_z_prc_pp__(1+ntolerance_pm,:) = prctile(tmp_,prctile_,'all');
clear tmp_;
end%;for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
% fraction of image-template-pairs within discretization error in gamma_z. ;
%%%%%%%%;
error_euler_gamma_z_f1_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_f2_ = zeros(n_tolerance_pm,1);
error_euler_gamma_z_f4_ = zeros(n_tolerance_pm,1);
dgamma = 2*pi/n_w_max;
for ntolerance_pm=0:n_tolerance_pm-1;
tmp_ = abs(periodize(euler_gamma_z_SMp___(:,:,end) - euler_gamma_z_SMp___(:,:,1+ntolerance_pm),-pi,pi));
error_euler_gamma_z_f1_(1+ntolerance_pm) = mean( tmp_> 1*dgamma + 1e-12 , 'all' );
error_euler_gamma_z_f2_(1+ntolerance_pm) = mean( tmp_> 2*dgamma + 1e-12 , 'all' );
error_euler_gamma_z_f4_(1+ntolerance_pm) = mean( tmp_> 4*dgamma + 1e-12 , 'all' );
clear tmp_;
end;%for ntolerance_pm=0:n_tolerance_pm-1;

figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,512,768]);
subplot(2,1,1);
hold on;
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_avg_p_,'k.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_prc_pp__,'.-')
hold off;
subplot(2,1,2);
hold on;
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f1_,'r.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f2_,'r.-');
plot(pm_n_UX_rank_avg_p_,error_euler_gamma_z_f4_,'r.-');
hold off;


pm_n_UX_rank_avg_p_ = mean(pm_n_UX_rank_cp__(1+index_ncluster_from_nM_,:),1);
tmp_ = reshape(X_SMp___,[n_S*n_M,n_tolerance_pm]);
corr_X_pp__ = corr(tmp_,tmp_);
corr_X_pe_ = corr_X_pp__(:,end);
%%%%%%%%;
prctile_ = [5,15,50,85,95]; n_prctile = numel(prctile_);
cap_pn__ = zeros(n_prctile,n_tolerance_pm);
X_SM_end__ = X_SMp___(:,:,end);
for nprctile=0:n_prctile-1;
cut_X_SM_end__ = zeros(size(X_SM_end__));
prc = prctile_(1+nprctile);
cut_X_SM_end__ = X_SM_end__ > repmat(prctile(X_SM_end__,prc),[n_S,1]);
for ntolerance_pm=0:n_tolerance_pm-1;
X_SM_sub__ = X_SMp___(:,:,1+ntolerance_pm);
cut_X_SM_sub__ = zeros(size(X_SM_sub__));
cut_X_SM_sub__ = X_SM_sub__ > repmat(prctile(X_SM_sub__,prc),[n_S,1]);
cap_pn__(1+nprctile,1+ntolerance_pm) = sum(cut_X_SM_end__.*cut_X_SM_sub__,'all');
end;%for ntolerance_pm=0:n_tolerance_pm-1;
end;%for nprctile=0:n_prctile-1;
cap_nrm_pn__ = cap_pn__./repmat(cap_pn__(:,end),[1,n_tolerance_pm]);
%%%%%%%%;

fname_mat = sprintf('%s_mat/pm_fig_X_2d_cluster_d0_timing_FIGI__.mat',dir_pm);
if (flag_recalc | ~exist(fname_mat,'file'));
disp(sprintf(' %% %s not found, creating',fname_mat));
%%%%%%%%;
tmp_t_mult_ = zeros(n_tolerance_pm,1);
tmp_t_ifft_ = zeros(n_tolerance_pm,1);
tmp_o_mult_ = zeros(n_tolerance_pm,1);
tmp_o_ifft_ = zeros(n_tolerance_pm,1);
for ntolerance_pm=0:n_tolerance_pm-1;%for ntolerance_pm=[0,5,n_tolerance_pm-1];
n_n = round(pm_n_UX_rank_avg_p_(1+ntolerance_pm));
tmp_n_w_max = n_w_max;
tmp_n_S = 1024; tmp_n_M = 1; 
%tmp_n_S = 1; tmp_n_M = 1024; 
tmp_n_S = min(tmp_n_S,n_S); tmp_n_M = min(tmp_n_M,n_M);
n_M_batch = ceil(n_M/tmp_n_M); n_S_batch = ceil(n_S/tmp_n_S);
tmp_t_mult = 0; tmp_t_ifft = 0;
tmp_o_mult = 0; tmp_o_ifft = 0;
for nS_batch=0:n_S_batch-1;for nM_batch=0:n_M_batch-1;
tmp_CTF_UX_S_k_q_Snw___ = randn(tmp_n_S,n_n,tmp_n_w_max) + i;
tmp_svd_VUXM_nMwl____ = randn(n_n,tmp_n_M,tmp_n_w_max) + i;
tmp_svd_SVUXM_SMwl____ = zeros(tmp_n_S,tmp_n_M,tmp_n_w_max);
tmp_t = tic();
for nw=0:tmp_n_w_max-1;
%tmp_svd_SVUXM_SMwl____(:,:,1+nw) = tmp_svd_SVUXM_SMwl____(:,:,1+nw) + tmp_CTF_UX_S_k_q_Snw___(:,:,1+nw)*tmp_svd_VUXM_nMwl____(:,:,1+nw);
tmp_svd_SVUXM_SMwl____(:,:,1+nw) = tmp_CTF_UX_S_k_q_Snw___(:,:,1+nw)*tmp_svd_VUXM_nMwl____(:,:,1+nw);
end;%for nw=0:tmp_n_w_max-1;
%tmp_svd_SVUXM_SMwl____ = multiprod(tmp_CTF_UX_S_k_q_Snw___,tmp_svd_VUXM_nMwl____);
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% tmp_svd_SVUXM_SMwl____: %0.6f',tmp_t)); end;
tmp_t_mult = tmp_t_mult + tmp_t;
tmp_o = tmp_n_w_max*tmp_n_S*tmp_n_M*n_n;
tmp_o_mult = tmp_o_mult + tmp_o;
tmp_t = tic();
tmp_svd_SVUXM_wSM___ = ifft(permute(tmp_svd_SVUXM_SMwl____,[3,1,2]),[],1)*tmp_n_w_max;
tmp_t = toc(tmp_t); if (verbose>1); disp(sprintf(' %% svd_SVUXM_wSM___: %0.6f',tmp_t)); end;
tmp_t_ifft = tmp_t_ifft + tmp_t;
tmp_o = tmp_n_w_max*log(tmp_n_w_max)*tmp_n_S*tmp_n_M;
tmp_o_ifft = tmp_o_ifft + tmp_o;
end;end;%for nS_batch=0:n_S_batch-1;for nM_batch=0:n_M_batch-1;
tmp_t_mult_(1+ntolerance_pm) = tmp_t_mult;
tmp_t_ifft_(1+ntolerance_pm) = tmp_t_ifft;
tmp_o_mult_(1+ntolerance_pm) = tmp_o_mult;
tmp_o_ifft_(1+ntolerance_pm) = tmp_o_ifft;
disp(sprintf(' %% ntolerance_pm %d n_n %d tmp_t_mult %0.5fs tmp_s_mult %0.5fs tmp_t_ifft %0.5fs tmp_s_ifft %0.5fs',ntolerance_pm,n_n,tmp_t_mult,tmp_o_mult/tmp_t_mult/1e9,tmp_t_ifft,tmp_o_ifft/tmp_t_ifft/1e9));
end;%for ntolerance_pm=0:n_tolerance_pm-1;
%%%%%%%%;
save(fname_mat ...
     ,'n_tolerance_pm','tolerance_pm_','n_w_max','n_S','n_M','tmp_n_w_max','tmp_n_S','tmp_n_M' ...
     ,'tmp_t_mult_','tmp_o_mult_','tmp_t_ifft_','tmp_o_ifft_' ...
     );
%%%%%%%%;
end;%if (~exist(fname_mat,'file'));
if ( exist(fname_mat,'file'));
disp(sprintf(' %% %s found, not creating',fname_mat));
load(fname_mat);
end;%if ( exist(fname_mat,'file'));
%%%%%%%%;

fname_fig = sprintf('%s_jpg/pm_fig_X_2d_cluster_d0_FIGI__',dir_pm);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1+nf);nf=nf+1;set(gcf,'Position',1+[0,0,512,768]);
p_row=3;p_col=1;ns=0;
tmp_ij=1:1:n_tolerance_pm; 
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
plot(pm_n_UX_rank_avg_p_(tmp_ij),-log10(tolerance_pm_(tmp_ij)),'ko-','MarkerFaceColor',0.85*[1,1,1]);
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylabel('$-\log_{10}(\mbox{tolerance})$','Interpreter','latex'); set(gca,'Ytick',0:2:20);
legend({'tolerance'},'Location','NorthWest');
grid on;
title('user-specified tolerance','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(pm_n_UX_rank_avg_p_(tmp_ij),corr_X_pe_(tmp_ij),'kd-','MarkerFaceColor',0.65*[1,1,1]);
plot(pm_n_UX_rank_avg_p_(tmp_ij),cap_nrm_pn__(end,tmp_ij),'ko-','MarkerFaceColor',[0.85,0.25,0.15]);
hold off;
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylabel('value','Interpreter','latex');ylim([0.5,1.0]);set(gca,'Ytick',0.5:0.1:1.0);
legend({'correlation','top 5%'},'Location','SouthEast');
grid on;
title('accuracy','Interpreter','latex');
%%%%;
subplot(p_row,p_col,1+ns);ns=ns+1;
hold on;
plot(pm_n_UX_rank_avg_p_(tmp_ij),max(tmp_t_mult_)./tmp_t_mult_(tmp_ij),'ks-','MarkerFaceColor',[0.15,0.35,0.85]);
plot(pm_n_UX_rank_avg_p_(tmp_ij),max(tmp_o_mult_)./tmp_o_mult_(tmp_ij),'k^-','MarkerFaceColor',[0.15,0.95,0.15]);
hold off;
xlim([0,n_k_p_r]); xlabel('average rank $H$','Interpreter','latex'); set(gca,'XTick',0:6:48);
ylim([0.0,8.0]);set(gca,'Ytick',0:1:8); ylabel('factor');
legend({'timing','operations'},'Location','NorthEast');
grid on;
title('efficiency','Interpreter','latex');
%%%%;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;

disp('returning'); return; 


