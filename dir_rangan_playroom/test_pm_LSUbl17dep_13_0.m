%%%%%%%%;
% First run test_pm_LSUbl17dep_13.m to load the data. ;
% Now we investigate heterogeneity. ;
%%%%%%%%;

verbose = 1;
nf=0;

%%%%%%%%;
% Measure correlation between volumes. ;
%%%%%%%%;
C_het_k_Y_quad__ = zeros(n_volume,n_volume);
C_het_k_p_quad__ = zeros(n_volume,n_volume);
C_het_x_u_pack__ = zeros(n_volume,n_volume);
for nvolume0=0:n_volume-1;
for nvolume1=nvolume0:n_volume-1;
tmp_t = tic();
[ ...
 tmp_C ...
] = ...
register_spharm_to_spharm_wigner_1( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,0 ...
,l_max_ ...
,a_het_k_Y_quad_hmlk__{1+nvolume0} ...
,a_het_k_Y_quad_hmlk__{1+nvolume1} ...
);
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% nvolume (%d,%d): register_spharm_to_spharm_wigner_1 time %0.2fs',nvolume0,nvolume1,tmp_t)); end;
C_het_k_Y_quad__(1+nvolume0,1+nvolume1) = tmp_C;
C_het_k_Y_quad__(1+nvolume1,1+nvolume0) = tmp_C;
tmp_C = real(corr(a_het_k_p_quad_hbak__{1+nvolume0},a_het_k_p_quad_hbak__{1+nvolume1}));
C_het_k_p_quad__(1+nvolume0,1+nvolume1) = tmp_C;
C_het_k_p_quad__(1+nvolume1,1+nvolume0) = tmp_C;
tmp_C = real(corr(a_het_x_u_pack_hxxx__{1+nvolume0}(:),a_het_x_u_pack_hxxx__{1+nvolume1}(:)));
C_het_x_u_pack__(1+nvolume0,1+nvolume1) = tmp_C;
C_het_x_u_pack__(1+nvolume1,1+nvolume0) = tmp_C;
end;%for nvolume1=nvolume0:n_volume-1;
end;%for nvolume0=0:n_volume-1;
figure(1+nf);nf=nf+1;clf;figmed;figbeach();
subplot(1,3,1); imagesc(C_het_k_Y_quad__,[0.5,1.0]); axis image; axisnotick; title('C_het_k_Y_quad__','Interpreter','none'); colorbar; 
subplot(1,3,2); imagesc(C_het_k_p_quad__,[0.5,1.0]); axis image; axisnotick; title('C_het_k_p_quad__','Interpreter','none'); colorbar; 
subplot(1,3,3); imagesc(C_het_x_u_pack__,[0.5,1.0]); axis image; axisnotick; title('C_het_x_u_quad__','Interpreter','none'); colorbar; 
%%%%%%%%;
% Note that volumes 0 and 5 are actually quite dissimilar from one another, ;
% while volumes 1-4 are quite similar. ;
%%%%%%%%;

flag_nrm = 0;
n_loading = 3;
n_loading_iteration = 16;
tmp_t = tic();
N_k_p_l2_ = zeros(n_M,1);
N_k_p_nrm__ = N_k_p__;
for nM=0:n_M-1;
N_k_p_ = N_k_p__(:,1+nM);
tmp_l2 = innerproduct_p_quad(n_k_p_r,k_p_r_,weight_2d_k_p_r_,n_w_,n_w_sum,N_k_p_,N_k_p_);
if flag_nrm; N_k_p_nrm__(:,1+nM) = N_k_p_/sqrt(tmp_l2); end;
N_k_p_l2_(1+nM) = tmp_l2;
end;%for nM=0:n_M-1;
tmp_t = toc(tmp_t); if (verbose); disp(sprintf(' %% N_k_p_nrm__ time %0.2fs',tmp_t)); end;
%%%%%%%%;
% test out get_loading_qbp_0 using the 'true' parameters. ;
%%%%%%%%;
for tmp_nvolume=0:n_volume-1;
%%%%%%%%;
tmp_euler_polar_a_true_ = transpose(euler_polar_a_true_hM__(1+tmp_nvolume,:));
tmp_euler_azimu_b_true_ = transpose(euler_azimu_b_true_hM__(1+tmp_nvolume,:));
tmp_euler_gamma_z_true_ = transpose(euler_gamma_z_true_hM__(1+tmp_nvolume,:));
tmp_image_delta_x_true_ = transpose(image_delta_x_true_hM__(1+tmp_nvolume,:));
tmp_image_delta_y_true_ = transpose(image_delta_y_true_hM__(1+tmp_nvolume,:));
tmp_fname_pos = sprintf('class_%c_from_ori','A'+tmp_nvolume);
tmp_fname_pre = sprintf('%s_mat/%s',dir_pm,tmp_fname_pos);
%%%%%%%%;
[tmp_flag_skip,tmp_fname_mat] = open_fname_tmp(tmp_fname_pre);
if ~tmp_flag_skip;
parameter = struct('type','parameter');
parameter.n_loading = n_loading;
parameter.n_loading_iteration = n_loading_iteration;
parameter.flag_loading_svd_vs_iterate = 0;
parameter.flag_SV_lsq_vs_dot = 1;
parameter.flag_U_lsq_vs_dot = 1;
[ ...
 parameter ...
,SV_loading_Ml__ ...
] = ...
get_loading_qbp_2( ...
 parameter ...
,n_k_p_r ...
,k_p_r_ ...
,weight_3d_k_p_r_ ...
,weight_2d_k_p_r_ ...
,l_max_max*ones(n_k_p_r,1) ...
,n_w_ ...
,n_M ...
,N_k_p_nrm__ ...
,index_nCTF_from_nM_ ...
,CTF_k_p__ ...
,tmp_euler_polar_a_true_ ...
,tmp_euler_azimu_b_true_ ...
,tmp_euler_gamma_z_true_ ...
,tmp_image_delta_x_true_ ...
,tmp_image_delta_y_true_ ...
);
save(tmp_fname_mat ...
     ,'n_loading','n_loading_iteration' ...
     ,'SV_loading_Ml__' ...
     );
close_fname_tmp(tmp_fname_pre);
end;%if ~tmp_flag_skip;
%%%%%%%%;
if  exist(tmp_fname_mat);
load(tmp_fname_mat);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/loading_het_%s_FIGA',dir_pm,tmp_fname_pos);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figsml;colormap('lines');
scatter3(SV_loading_Ml__(:,1),SV_loading_Ml__(:,2),SV_loading_Ml__(:,3),24,index_nvolume_from_nM_,'filled');
axis vis3d;
title('loading via true');
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/loading_het_%s_FIGB',dir_pm,tmp_fname_pos);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;figbeach();
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach__(end,:) = [1,1,1];
markersize_use = 3;
npc0 = 0;
npc1 = 1;
llim_0_ = prctile(SV_loading_Ml__(:,1+npc0),[2,98]);
llim_1_ = prctile(SV_loading_Ml__(:,1+npc1),[2,98]);
for nvolume=0:(1+n_volume)-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*nvolume/(n_volume))));
tmp_index_ = 0:n_M-1;
if (nvolume<n_volume); tmp_index_ = efind(index_nvolume_from_nM_==nvolume); end;
subplot(2,(1+n_volume),1+nvolume+0*(1+n_volume));
plot(SV_loading_Ml__(1+tmp_index_,1+npc0),SV_loading_Ml__(1+tmp_index_,1+npc1),'ko','MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc,:));
xlim(llim_0_); ylim(llim_1_); grid on ; axisnotick;
title(sprintf('%d',nvolume));
subplot(2,(1+n_volume),1+nvolume+1*(1+n_volume));
imagesc(log(1+hist2d_0(SV_loading_Ml__(1+tmp_index_,1+npc0),SV_loading_Ml__(1+tmp_index_,1+npc1),32,32,llim_0_,llim_1_)),[0,5]);
axis image; axisnotick;
title(sprintf('%d',nvolume));
end;%for nvolume=0:(1+n_volume)-1;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
fname_fig = sprintf('%s_jpg/loading_het_%s_FIGC',dir_pm,tmp_fname_pos);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figbig;figbeach();
c_beach__ = colormap_beach(); n_c_beach = size(c_beach__,1); c_beach__(end,:) = [1,1,1];
markersize_use = 12;
npc0 = 0;
npc1 = 1;
llim_0_ = prctile(SV_loading_Ml__(:,1+npc0),[2,98]);
llim_1_ = prctile(SV_loading_Ml__(:,1+npc1),[2,98]);
subplot(1,1,1);
hold on;
for nvolume=0:n_volume-1;
nc = max(0,min(n_c_beach-1,floor(n_c_beach*nvolume/(n_volume))));
tmp_index_ = 0:n_M-1;
if (nvolume<n_volume); tmp_index_ = efind(index_nvolume_from_nM_==nvolume); end;
plot(SV_loading_Ml__(1+tmp_index_,1+npc0),SV_loading_Ml__(1+tmp_index_,1+npc1),'ko','MarkerSize',markersize_use,'MarkerFaceColor',c_beach__(1+nc,:));
xlim(llim_0_); ylim(llim_1_); grid on ; axisnotick;
axis equal;
title(sprintf('%d',nvolume));
end;%for nvolume=0:(1+n_volume)-1;
hold off;
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
% check for correlation between loading and X_value. ;
%%%%%%%%;
SV_loading_l2_ = sum(SV_loading_Ml__.^2,2);
plot(log(SV_loading_l2_),image_X_value_true_,'.');
%%%%%%%%;
% Now construct and plot histogram of various molecular types. ;
%%%%%%%%;
npc0=0;
n_h = 48;
hlim_ = prctile(SV_loading_Ml__(:,1),[  0,100]); hlim_ = mean(hlim_) + 1.25*0.5*diff(hlim_)*[-1,+1];
hbin_ = linspace(min(hlim_),max(hlim_),n_h);
h_het__ = zeros(n_h,1+n_volume);
for nvolume=0:(1+n_volume)-1;
tmp_index_ = 0:n_M-1;
if (nvolume<n_volume); tmp_index_ = efind(index_nvolume_from_nM_==nvolume); end;
h_het__(:,1+nvolume) = hist(SV_loading_Ml__(1+tmp_index_,1+npc0),hbin_);
end;%for nvolume=0:(1+n_volume)-1;
h_all_ = h_het__(:,1+n_volume);
dh = mean(diff(hbin_));
h_denominator = sum(h_all_)*dh;
h_all_ = h_all_/h_denominator;
h_het__ = h_het__/h_denominator;
%%%%%%%%;
[ ...
 tmp_p_threshold ...
,tmp_param_out_ ...
,tmp_log_likelihood_ratio ...
,tmp_label_ ...
,tmp_label_auc ...
] = ...
so2g_mle_fminsearch(SV_loading_Ml__(:,1+npc0));
tmp_m_A = tmp_param_out_(1);
tmp_s_A = tmp_param_out_(2);
tmp_m_B = tmp_param_out_(3);
tmp_s_B = tmp_param_out_(4);
tmp_l = tmp_param_out_(5);
fp_p = @(x,m,s) 1/sqrt(2*pi)./s .* exp(-(x-m).^2./(2*s.^2));
fp_lp = @(x,m,s) log(fp_p(x,m,s));
fp_lA = @(l) 1./(1+exp(-l)) ; 
fp_lB = @(l) exp(-l)./(1+exp(-l)) ;
fp_rho = @(x,m_A,s_A,m_B,s_B,l) fp_lA(l).*fp_p(x,m_A,s_A) + fp_lB(l).*fp_p(x,m_B,s_B) ;
%%%%%%%%;
ylim_ = [0,max(h_all_)*1.125];
tmp_x_ = hbin_ - mean(diff(hbin_))*0.5;
%%%%%%%%;
index_AF_ = unionall({efind(index_nvolume_from_nM_==0),efind(index_nvolume_from_nM_==5)});
h_AF_ = sum(h_het__(:,1+[0,5]),2);
index_BE_ = unionall({efind(index_nvolume_from_nM_==1),efind(index_nvolume_from_nM_==2),efind(index_nvolume_from_nM_==3),efind(index_nvolume_from_nM_==4)});
h_BE_ = sum(h_het__(:,1+[1,2,3,4]),2);
auc_AF_vs_BE = auc_0(SV_loading_Ml__(1+index_BE_,1),SV_loading_Ml__(1+index_AF_,1));
auc_AF_vs_BE = max(auc_AF_vs_BE,1-auc_AF_vs_BE);
%%%%%%%%;
fname_fig = sprintf('%s_jpg/loading_het_%s_FIGD',dir_pm,tmp_fname_pos);
if (flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s not found, creating',fname_fig));
figure(1);clf;figsml;figbeach();
hold on;
stairs(tmp_x_,h_all_,'k-','LineWidth',1);
stairs(tmp_x_,h_AF_ ,'r-','LineWidth',3);
stairs(tmp_x_,h_BE_ ,'g-','LineWidth',3);
stairs(tmp_x_,h_all_,'k-');
%plot(tmp_x_,fp_rho(tmp_x_,tmp_m_A,tmp_s_A,tmp_m_B,tmp_s_B,tmp_l),'g-');
%plot(tmp_x_,fp_lA(tmp_l).*fp_p(tmp_x_,tmp_m_A,tmp_s_A),'r-');
%plot(tmp_x_,fp_lB(tmp_l).*fp_p(tmp_x_,tmp_m_B,tmp_s_B),'b-');
%plot(tmp_p_threshold*[1,1],[0,fp_rho(tmp_p_threshold,tmp_m_A,tmp_s_A,tmp_m_B,tmp_s_B,tmp_l)],'g-','LineWidth',3);
plot(tmp_p_threshold*[1,1],[0,max(ylim_)],'b-','LineWidth',3);
%stairs(tmp_x_,h_het__(:,1+0),'k-','LineWidth',1);
%stairs(tmp_x_,h_het__(:,1+1),'r-','LineWidth',1);
%stairs(tmp_x_,h_het__(:,1+2),'m-','LineWidth',1);
%stairs(tmp_x_,h_het__(:,1+3),'g-','LineWidth',1);
%stairs(tmp_x_,h_het__(:,1+4),'b-','LineWidth',1);
%stairs(tmp_x_,h_het__(:,1+5),'y-','LineWidth',1);
hold off;
title(sprintf('AUC=%0.3f',auc_AF_vs_BE),'Interpreter','none');
xlim(hlim_);ylim(ylim_);
disp(sprintf(' %% writing %s',fname_fig));
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close(gcf);
end;%if (~exist(sprintf('%s.jpg',fname_fig),'file'));
if ( exist(sprintf('%s.jpg',fname_fig),'file'));
disp(sprintf(' %% %s found, not creating',fname_fig));
end;%if ( exist(sprintf('%s.jpg',fname_fig),'file'));
%%%%%%%%;
end;%if  exist(tmp_fname_mat);
%%%%%%%%;
end;%for tmp_nvolume=0:n_volume-1;
