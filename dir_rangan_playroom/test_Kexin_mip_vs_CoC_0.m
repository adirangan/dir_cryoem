%%%%%%%%;
% testing dir_Kexin_mip_vs_CoC: ;
%{
result_search_5_2_radius_10.csv
result_search_5_2_radius_50.csv
result_search_5_3_radius_10.csv
result_search_5_3_radius_50.csv
result_search_9_2_radius_10.csv
result_search_9_2_radius_50.csv
result_search_9_3_radius_10.csv
result_search_9_3_radius_50.csv
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_2.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_3.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_2.mrc
simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_3.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_2.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_3.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_2.mrc
simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_3.mrc
  %}
%%%%%%%%;

%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;

flag_verbose=1; flag_disp=1; nf=0;
flag_replot = 1;

dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
dir_Kexin = sprintf('%s/dir_Kexin_mip_vs_CoC',dir_base);
dir_jpg = sprintf('%s/dir_jpg',dir_Kexin);
if ~exist(dir_jpg,'dir'); disp(sprintf(' %% mkdir %s',dir_jpg)); mkdir(dir_jpg); end;
dir_star = sprintf('%s/dir_Ground_Truth_star',dir_Kexin);

fname_mip_9_2_mrc = sprintf('%s/simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_2.mrc',dir_Kexin);
tmp_mip_9_2_ = ReadMRC_0(fname_mip_9_2_mrc);
fname_CoC_9_2_mrc = sprintf('%s/simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_2.mrc',dir_Kexin);
tmp_CoC_9_2_ = ReadMRC_0(fname_CoC_9_2_mrc);
fname_mip_9_3_mrc = sprintf('%s/simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_scaled_mip_9_3.mrc',dir_Kexin);
tmp_mip_9_3_ = ReadMRC_0(fname_mip_9_3_mrc);
fname_CoC_9_3_mrc = sprintf('%s/simulator_default_clathrin_exp2_t1000A_d5000A_1.06Apix_45e_montage_boxsize_512_coc_9_3.mrc',dir_Kexin);
tmp_CoC_9_3_ = ReadMRC_0(fname_CoC_9_3_mrc);
fname_mip_5_2_mrc = sprintf('%s/simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_2.mrc',dir_Kexin);
tmp_mip_5_2_ = ReadMRC_0(fname_mip_5_2_mrc);
fname_CoC_5_2_mrc = sprintf('%s/simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_2.mrc',dir_Kexin);
tmp_CoC_5_2_ = ReadMRC_0(fname_CoC_5_2_mrc);
fname_mip_5_3_mrc = sprintf('%s/simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_scaled_mip_5_3.mrc',dir_Kexin);
tmp_mip_5_3_ = ReadMRC_0(fname_mip_5_3_mrc);
fname_CoC_5_3_mrc = sprintf('%s/simulator_default_mature60S_exp2_t1000A_d5000A_1.06Apix_45e_box512_montage_coc_5_3.mrc',dir_Kexin);
tmp_CoC_5_3_ = ReadMRC_0(fname_CoC_5_3_mrc);
infix_mip_9_2 = 'clathrin_mip';
infix_CoC_9_2 = 'clathrin_CoC';
infix_mip_9_3 = 'clathrin_mip';
infix_CoC_9_3 = 'clathrin_CoC';
infix_mip_5_2 = 'mature60S_mip';
infix_CoC_5_2 = 'mature60S_CoC';
infix_mip_5_3 = 'mature60S_mip';
infix_CoC_5_3 = 'mature60S_CoC';

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% display the mrc files. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;fig81s;
fontsize_use = 12;
p_row = 2; p_col = 4; np=0;
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_mip_9_2_);axis image; axisnotick; title(infix_mip_9_2,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_CoC_9_2_);axis image; axisnotick; title(infix_CoC_9_2,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_mip_9_3_);axis image; axisnotick; title(infix_mip_9_3,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_CoC_9_3_);axis image; axisnotick; title(infix_CoC_9_3,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_mip_5_2_);axis image; axisnotick; title(infix_mip_5_2,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_CoC_5_2_);axis image; axisnotick; title(infix_CoC_5_2,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_mip_5_3_);axis image; axisnotick; title(infix_mip_5_3,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;cla;
imagesc(tmp_CoC_5_3_);axis image; axisnotick; title(infix_CoC_5_3,'Interpreter','none'); set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/test_Kexin_mip_vs_CoC_display_FIGA',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% readtable. ;
%%%%%%%%;
Table_9_3_10_rc__ = readtable(sprintf('%s/result_search_9_3_radius_10.csv',dir_Kexin));
Array_9_3_10_rc__ = table2array(Table_9_3_10_rc__);
Table_9_3_10_VariableNames_c_ = Table_9_3_10_rc__.Properties.VariableNames;
% try: sprintf('%s\n',Table_9_3_10_VariableNames_c_{:}) ;
nc=0;
nc_10_Psi = nc; nc=nc+1; %<-- euler angles best matched ;
nc_10_Theta = nc; nc=nc+1; %<-- euler angles best matched ;
nc_10_Phi = nc; nc=nc+1; %<-- euler angles best matched ;
nc_10_X = nc; nc=nc+1; %<-- coordinates ;
nc_10_Y = nc; nc=nc+1; %<-- coordinates ;
nc_10_Z = nc; nc=nc+1; %<-- defocus ;
nc_10_PixelSize = nc; nc=nc+1;
nc_10_Peak = nc; nc=nc+1; %<-- Peak: z-score ;
nc_10_ScalingFactor = nc; nc=nc+1; %<-- ScalingFactor: the ratio between z-score/CC_max ;
nc_10_CoC = nc; nc=nc+1; %<-- coc: correlation of correlations ;
nc_10_mip = nc; nc=nc+1; %<-- mip: the maximum CC ;
nc_10_pval = nc; nc=nc+1; %<-- pval: pval based on chi-2 dist of the null hypothesis ;
nc_10_log_pval = nc; nc=nc+1; %<-- log_pval: when pval=inf, I replace it with the min_pval - 1 s.t. it is easier show the p values using color map ;
nc_10_neg_logpval = nc; nc=nc+1;
n_c_10 = nc;
%%%%%%%%;
Table_9_3_50_rc__ = readtable(sprintf('%s/result_search_9_3_radius_50.csv',dir_Kexin));
Array_9_3_50_rc__ = table2array(Table_9_3_50_rc__);
Table_9_3_50_VariableNames_c_ = Table_9_3_50_rc__.Properties.VariableNames;
% try: sprintf('%s\n',Table_9_3_50_VariableNames_c_{:}) ;
nc=0;
nc_50_Psi = nc; nc=nc+1; %<-- euler angles best matched ;
nc_50_Theta = nc; nc=nc+1; %<-- euler angles best matched ;
nc_50_Phi = nc; nc=nc+1; %<-- euler angles best matched ;
nc_50_X = nc; nc=nc+1; %<-- coordinates ;
nc_50_Y = nc; nc=nc+1; %<-- coordinates ;
nc_50_Z = nc; nc=nc+1; %<-- defocus ;
nc_50_PixelSize = nc; nc=nc+1;
nc_50_Peak = nc; nc=nc+1; %<-- Peak: z-score ;
nc_50_ScalingFactor = nc; nc=nc+1; %<-- ScalingFactor: the ratio between z-score/CC_max ;
nc_50_CoC = nc; nc=nc+1; %<-- coc: correlation of correlations ;
nc_50_mip = nc; nc=nc+1; %<-- mip: the maximum CC ;
nc_50_labels = nc; nc=nc+1; %<-- labels: true or false based on translational error to grid center+pose similarity to ground truth ;
nc_50_grid_id = nc; nc=nc+1; %<-- grid_id: which particle grid this peak belongs to (so that we know the ground truth) ;
nc_50_pose_sim = nc; nc=nc+1; %<-- pose_sim: similarity between 2dtm-derived template to ground truth template ;
nc_50_d_xy = nc; nc=nc+1; %<-- d_xy: distance to center of grid ;
nc_50_CoC_scaled = nc; nc=nc+1; %<-- coc_scaled: just out of curiosity, where I apply the scaling of z-score to coc ;
nc_50_pval = nc; nc=nc+1; %<-- pval: pval based on chi-2 dist of the null hypothesis ;
nc_50_log_pval = nc; nc=nc+1; %<-- log_pval: when pval=inf, I replace it with the min_pval - 1 s.t. it is easier show the p values using color map ;
nc_50_neg_logpval = nc; nc=nc+1;
n_c_50 = nc;
%%%%%%%%;

lab_50_r_ = Array_9_3_50_rc__(:,1+nc_50_labels);
mip_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Peak);
CoC_50_r_ = Array_9_3_50_rc__(:,1+nc_50_CoC);
mip_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Peak);
CoC_10_r_ = Array_9_3_10_rc__(:,1+nc_10_CoC);

%%%%%%%%;
% check for errors. ;
%%%%%%%%;
tmp_mip_error = 0;
tmp_CoC_error = 0;
tmp_index_1_ = efind(lab_50_r_==1);
for nl=0:numel(tmp_index_1_)-1;
tmp_index_1 = tmp_index_1_(1+nl);
tmp_nx = Array_9_3_50_rc__(1+tmp_index_1,1+nc_50_X);
tmp_ny = Array_9_3_50_rc__(1+tmp_index_1,1+nc_50_Y);
tmp_mip_error = tmp_mip_error + abs(mip_50_r_(1+tmp_index_1) - tmp_mip_9_3_(1+tmp_nx,1+tmp_ny));
tmp_CoC_error = tmp_CoC_error + abs(CoC_50_r_(1+tmp_index_1) - tmp_CoC_9_3_(1+tmp_nx,1+tmp_ny));
end;%for nl=0:numel(tmp_index_1_)-1;
disp(sprintf(' %% tmp_mip_error: %0.16f',tmp_mip_error));
disp(sprintf(' %% tmp_CoC_error: %0.16f',tmp_CoC_error));

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Align the detected points to the mip mrc. ;
% Note that the coordinates are transposed and zero-based. ;
%%%%%%%%;
figure(1);clf;figbig;fig81s;fontsize_use = 12;
markersize_big = 12;
markersize_med = 8;
markersize_sml = 4;
linewidth_med = 1;
linewidth_sml = 0.5;
subplot(1,2,1);
hold on;
imagesc(tmp_mip_9_3_);
plot(1+Array_9_3_50_rc__(:,1+nc_50_Y),1+Array_9_3_50_rc__(:,1+nc_50_X),'go','MarkerSize',markersize_sml,'LineWidth',linewidth_med);
%plot(1+Array_9_3_10_rc__(:,1+nc_10_Y),1+Array_9_3_10_rc__(:,1+nc_10_X),'go','MarkerSize',markersize_sml,'LineWidth',linewidth_sml);
hold off;
axis image; axisnotick; title(sprintf('%s (full): 50: only',infix_mip_9_3),'Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(1,2,2);
hold on;
imagesc(tmp_mip_9_3_);
plot(1+Array_9_3_50_rc__(:,1+nc_50_Y),1+Array_9_3_50_rc__(:,1+nc_50_X),'go');
plot(1+Array_9_3_50_rc__(:,1+nc_50_Y),1+Array_9_3_50_rc__(:,1+nc_50_X),'go','MarkerSize',markersize_med,'LineWidth',linewidth_med);
plot(1+Array_9_3_10_rc__(:,1+nc_10_Y),1+Array_9_3_10_rc__(:,1+nc_10_X),'go','MarkerSize',markersize_sml,'LineWidth',linewidth_sml);
hold off;
axis(1024*2+[0,512-1,0,512-1]); axis square;
axisnotick; title(sprintf('%s (512x512): 50 big+sml; 10 sml only',infix_mip_9_3),'Interpreter','none'); set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/test_Kexin_mip_9_3_50_display_FIGB',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Show the z-score vs CoC scatterplot, colored by label. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
fontsize_use = 12;
hold on;
plot(mip_10_r_,CoC_10_r_,'x','Color',0.75*[1,1,1]);
tmp_index_0_ = efind(lab_50_r_==0);
plot(mip_50_r_(1+tmp_index_0_),CoC_50_r_(1+tmp_index_0_),'o','Color',0.35*[1,1,1]);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'o','Color',1.00*[1,0,1]);
hold off;
legend({'50 xx','10 fp','10 tp'},'Location','SouthEast');
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_10 + 9_3_50','Interpreter','none');
set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/test_Kexin_mip_vs_CoC_9_3_FIGC',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% plot Peak vs CoC vs d_xy. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
tmp_index_0_ = efind(lab_50_r_==0);
plot3(mip_50_r_(1+tmp_index_0_),CoC_50_r_(1+tmp_index_0_),Array_9_3_50_rc__(1+tmp_index_0_,1+nc_50_d_xy),'ko');
tmp_index_1_ = efind(lab_50_r_==1);
plot3(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_d_xy),'mo');
hold off;
axis vis3d;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% plot pose_sim vs d_xy. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
hold on;
tmp_index_0_ = efind(lab_50_r_==0);
plot(Array_9_3_50_rc__(1+tmp_index_0_,1+nc_50_pose_sim),Array_9_3_50_rc__(1+tmp_index_0_,1+nc_50_d_xy),'ko');
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_pose_sim),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_d_xy),'mo');
hold off;
end;%if flag_disp;

flag_disp=0;
if flag_disp;
%%%%%%%%;
% plot Peak vs CoC as well as net. ;
%%%%%%%%;
[n_x,n_y] = size(tmp_mip_9_3_);
tmp_r = 1; %<-- radius of square net. ;
figure(1+nf);nf=nf+1;clf;figbig;
fontsize_use = 12;
linewidth_sml = 1.5; alpha_use = 0.5;
c_fp_ = 0.65*[1,1,1]; c_fp2_ = 0.35*[1,1,1];
c_tp_ = [1,0.85,1]; c_tp2_ = [1,0,1];
hold on;
%%%%%%%%;
tmp_index_0_ = efind(lab_50_r_==0);
%plot(mip_50_r_(1+tmp_index_0_),CoC_50_r_(1+tmp_index_0_),'ko');
for nl=0:numel(tmp_index_0_)-1;
tmp_index_0 = tmp_index_0_(1+nl);
tmp_mip = mip_50_r_(1+tmp_index_0);
tmp_CoC = Array_9_3_50_rc__(1+tmp_index_0,1+nc_50_CoC);
%%%%;
if (tmp_mip> 6.0) & (tmp_CoC> 0.01);
tmp_nx = Array_9_3_50_rc__(1+tmp_index_0,1+nc_50_X);
tmp_ny = Array_9_3_50_rc__(1+tmp_index_0,1+nc_50_Y);
tmp_nx_ = max(0,min(n_x-1,tmp_nx + [-tmp_r:+tmp_r]));
tmp_ny_ = max(0,min(n_y-1,tmp_ny + [-tmp_r:+tmp_r]));
[tmp_nx__,tmp_ny__] = ndgrid(tmp_nx_,tmp_ny_);
tmp_index_ = tmp_nx__(:) + tmp_ny__(:)*n_x;
tmp_mip__ = reshape(tmp_mip_9_3_(tmp_index_),[1+2*tmp_r,1+2*tmp_r]);
tmp_CoC__ = reshape(tmp_CoC_9_3_(tmp_index_),[1+2*tmp_r,1+2*tmp_r]);
l=line(          tmp_mip__ ,          tmp_CoC__ ,'Color',[c_fp_,alpha_use],'LineWidth',linewidth_sml);
l=line(transpose(tmp_mip__),transpose(tmp_CoC__),'Color',[c_fp_,alpha_use],'LineWidth',linewidth_sml);
plot(tmp_mip,tmp_CoC,'o','Color',[c_fp_],'MarkerFaceColor',c_fp2_);
end;%if (tmp_mip> 6.0) & (tmp_CoC> 0.01);
%%%%;
end;%for nl=0:numel(tmp_index_0_)-1;
%%%%%%%%;
tmp_index_1_ = efind(lab_50_r_==1);
%plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'mo');
for nl=0:numel(tmp_index_1_)-1;
tmp_index_1 = tmp_index_1_(1+nl);
tmp_mip = mip_50_r_(1+tmp_index_1);
tmp_CoC = CoC_50_r_(1+tmp_index_1);
%%%%;
if (tmp_mip> 6.0) & (tmp_CoC> 0.01);
tmp_nx = Array_9_3_50_rc__(1+tmp_index_1,1+nc_50_X);
tmp_ny = Array_9_3_50_rc__(1+tmp_index_1,1+nc_50_Y);
tmp_nx_ = max(0,min(n_x-1,tmp_nx + [-tmp_r:+tmp_r]));
tmp_ny_ = max(0,min(n_y-1,tmp_ny + [-tmp_r:+tmp_r]));
[tmp_nx__,tmp_ny__] = ndgrid(tmp_nx_,tmp_ny_);
tmp_index_ = tmp_nx__(:) + tmp_ny__(:)*n_x;
tmp_mip__ = reshape(tmp_mip_9_3_(tmp_index_),[1+2*tmp_r,1+2*tmp_r]);
tmp_CoC__ = reshape(tmp_CoC_9_3_(tmp_index_),[1+2*tmp_r,1+2*tmp_r]);
l=line(          tmp_mip__ ,          tmp_CoC__ ,'Color',[c_tp_,alpha_use],'LineWidth',linewidth_sml);
l=line(transpose(tmp_mip__),transpose(tmp_CoC__),'Color',[c_tp_,alpha_use],'LineWidth',linewidth_sml);
plot(tmp_mip,tmp_CoC,'o','Color',[c_tp_],'MarkerFaceColor',c_tp2_);
end;%if (tmp_mip> 6.0) & (tmp_CoC> 0.01);
%%%%;
end;%for nl=0:numel(tmp_index_1_)-1;
%%%%%%%%;
hold off;
%%%%%%%%;
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_50 + 3x3 net','Interpreter','none');
set(gca,'FontSize',fontsize_use);
fname_fig_pre = sprintf('%s/test_Kexin_mip_vs_CoC_9_3_FIGD',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
end;%if flag_disp;

%%%%%%%%;
% Now load the various star-files. ;
%%%%%%%%;
prefix_star = 'simulated_clathrin_exposure2_t1000A_d5000A_1.06Apix_45e_';
posfix_star = '.mrc.star';
n_cell = 10;
psi_c_ = zeros(n_cell^2,1);
the_c_ = zeros(n_cell^2,1);
phi_c_ = zeros(n_cell^2,1);
for ncell=0:n_cell^2-1;
fname_star = sprintf('%s/%s%d%s',dir_star,prefix_star,1+ncell,posfix_star);
%%%%;
% columns from star-file: ;
% #    POS     PSI   THETA     PHI       SHX       SHY      DF1      DF2  ANGAST  PSHIFT     OCC      LogP      SIGMA   SCORE  CHANGE    PSIZE    VOLT      Cs    AmpC  BTILTX  BTILTY  ISHFTX  ISHFTY 2DCLS  TGRP    PaGRP  SUBSET  PREEXP  TOTEXP ;
%%%%;
fp=fopen(fname_star);
tmp_ = textscan(fp,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',35);
fclose(fp);
tmp_ = cell2mat(tmp_);
tmp_psi = tmp_(1+1); tmp_the = tmp_(1+2); tmp_phi = tmp_(1+3);
psi_c_(1+ncell) = tmp_psi; the_c_(1+ncell) = tmp_the; phi_c_(1+ncell) = tmp_phi;
clear fname_star tmp_ tmp_psi tmp_the tmp_phi ;
end;%for ncell=0:n_cell^2-1;
psi_c_ = psi_c_*pi/180.0; the_c_ = the_c_*pi/180.0; phi_c_ = phi_c_*pi/180.0;
psi_c__ = reshape(psi_c_,[n_cell,n_cell]);
the_c__ = reshape(the_c_,[n_cell,n_cell]);
phi_c__ = reshape(phi_c_,[n_cell,n_cell]); %<-- periodize. ;
phi_c__ = periodize(phi_c__,0,2*pi);
%%%%%%%%;
% Now define nx and ny for each cell. ;
%%%%%%%%;
ny_c__ = repmat(5120*[0.5:1:n_cell-0.5]/n_cell,[n_cell,1]);
nx_c__ = transpose(ny_c__);

flag_disp=0;
%%%%%%%%;
% Now visualize the psi_c_ vs Psi. ;
%%%%%%%%;
if flag_disp;
tmp_index_1_ = efind(lab_50_r_==1);
n_l = numel(tmp_index_1_);
tmp_nx_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_X);
tmp_ny_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Y);
tmp_mip_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak);
tmp_CoC_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC);
tmp_psi_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Psi);
tmp_the_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Theta);
tmp_phi_l_ = Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Phi);
tmp_psi_l_ = tmp_psi_l_*pi/180.0;
tmp_the_l_ = tmp_the_l_*pi/180.0;
tmp_phi_l_ = tmp_phi_l_*pi/180.0;
figure(1+nf);nf=nf+1;clf;figbig;
c_use__ = colormap('hsv'); n_c_use = size(c_use__,1);
markersize_big = 12;
%%%%;
for ntype=0:3-1;
if ntype==0; tmp_str = 'psi'; tmp_c__ = psi_c__; tmp_index = 0; end;
if ntype==1; tmp_str = 'the'; tmp_c__ = the_c__; tmp_index = 1; end;
if ntype==2; tmp_str = 'phi'; tmp_c__ = phi_c__; tmp_index = 2; end;
%%%%;
subplot(3,3,1+ntype+0*3);
hold on;
imagesc_c(n_cell,0:n_cell-1,n_cell,0:n_cell-1,tmp_c__,[0,2*pi],colormap('hsv'));
hold off;
axis image; axisnotick; set(gca,'Ydir','normal');
xlabel('row'); ylabel('col'); title(sprintf('%s_c__',tmp_str),'Interpreter','none');
%%%%;
subplot(3,3,1+ntype+1*3);
hold on;
for ncell=0:n_cell^2-1;
tmp_nx = nx_c__(1+ncell); tmp_ny = ny_c__(1+ncell);
tmp_psi = psi_c__(1+ncell); tmp_the = the_c__(1+ncell); tmp_phi = phi_c__(1+ncell);
tmp_eul_ = [tmp_psi,tmp_the,tmp_phi];
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(tmp_eul_(1+tmp_index)-0)/(2*pi))));
plot(tmp_nx,tmp_ny,'ko','MarkerSize',markersize_big,'MarkerFaceColor',c_use__(1+nc_use,:));
end;%for ncell=0:n_cell^2-1;
hold off;
xlim(5120*[0,1]);ylim(5120*[0,1]);
axis square; axisnotick;
xlabel('nx'); ylabel('ny'); title(sprintf('Cell %s',tmp_str),'Interpreter','none');
%%%%;
subplot(3,3,1+ntype+2*3);
hold on;
for nl=0:n_l-1;
tmp_nx = tmp_nx_l_(1+nl); tmp_ny = tmp_ny_l_(1+nl);
tmp_mip = tmp_mip_l_(1+nl); tmp_CoC = tmp_CoC_l_(1+nl);
tmp_psi = tmp_psi_l_(1+nl); tmp_the = tmp_the_l_(1+nl); tmp_phi = tmp_phi_l_(1+nl);
tmp_eul_ = [tmp_psi,tmp_the,tmp_phi];
nc_use = max(0,min(n_c_use-1,floor(n_c_use*(tmp_eul_(1+tmp_index)-0)/(2*pi))));
plot(tmp_nx,tmp_ny,'ko','MarkerSize',markersize_big,'MarkerFaceColor',c_use__(1+nc_use,:));
end;%for nl=0:n_l-1;
hold off;
xlim(5120*[0,1]);ylim(5120*[0,1]);
axis square; axisnotick;
xlabel('nx'); ylabel('ny'); title(sprintf('Array %s',tmp_str),'Interpreter','none');
%%%%;
end;%for ntype=0:3-1;
%%%%;
fname_fig_pre = sprintf('%s/test_Kexin_euler_9_3_FIGE',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

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
n_gamma_z = 1024; gamma_z_ = linspace(0,2*pi,n_gamma_z+1); gamma_z_ = transpose(gamma_z_(1:n_gamma_z));
ring_k_c_0_ = cos(gamma_z_);
ring_k_c_1_ = sin(gamma_z_);
ring_k_c_2_ = zeros(n_gamma_z,1);
ring_k_c_3z__ = transpose([ring_k_c_0_,ring_k_c_1_,ring_k_c_2_]);
%%%%%%%%;

%%%%%%%%;
% Now estimate the difference in euler-angle between the data and the ground_truth. ;
% For this difference we calculate the average squared-distance between points in the two templates. ;
% I am not sure which angle is azimu_b and which is gamma_z, ;
% so I take the minimum of both interpretations. ;
% Note that this notion of distance is only valid in the regime where the molecule is 'orientable' ;
% i.e., in which there are no symmetries. ;
%%%%%%%%;
n_l = size(Array_9_3_10_rc__,1);
tmp_nx_l_ = Array_9_3_10_rc__(:,1+nc_10_X);
tmp_ny_l_ = Array_9_3_10_rc__(:,1+nc_10_Y);
tmp_mip_l_ = mip_10_r_;
tmp_CoC_l_ = CoC_10_r_;
tmp_psi_l_ = Array_9_3_10_rc__(:,1+nc_10_Psi);
tmp_the_l_ = Array_9_3_10_rc__(:,1+nc_10_Theta);
tmp_phi_l_ = Array_9_3_10_rc__(:,1+nc_10_Phi);
tmp_psi_l_ = tmp_psi_l_*pi/180.0;
tmp_the_l_ = tmp_the_l_*pi/180.0;
tmp_phi_l_ = tmp_phi_l_*pi/180.0;
[tmp_ij_,dxy_10_r_] = knnsearch([nx_c__(:),ny_c__(:)],[tmp_nx_l_,tmp_ny_l_]); tmp_index_ = tmp_ij_ - 1;
tmp_psi_n_ = reshape(psi_c__(1+tmp_index_),[n_l,1]);
tmp_the_n_ = reshape(the_c__(1+tmp_index_),[n_l,1]);
tmp_phi_n_ = reshape(phi_c__(1+tmp_index_),[n_l,1]);
tmp_err_phi_the_psi_l_ = zeros(n_l,1);
tmp_err_psi_the_phi_l_ = zeros(n_l,1);
tmp_err_l_ = zeros(n_l,1);
for nl=0:n_l-1;
tmp_ring_est_k_c_3z__ = Rz(tmp_phi_l_(1+nl))*Ry(tmp_the_l_(1+nl))*Rz(tmp_psi_l_(1+nl))*ring_k_c_3z__;
tmp_ring_tru_k_c_3z__ = Rz(tmp_phi_n_(1+nl))*Ry(tmp_the_n_(1+nl))*Rz(tmp_psi_n_(1+nl))*ring_k_c_3z__;
tmp_ring_l2 = sum((tmp_ring_est_k_c_3z__ - tmp_ring_tru_k_c_3z__).^2,'all')*2*pi/max(1,n_gamma_z);
tmp_err_phi_the_psi_l_(1+nl) = sqrt(tmp_ring_l2);
tmp_ring_est_k_c_3z__ = Rz(tmp_psi_l_(1+nl))*Ry(tmp_the_l_(1+nl))*Rz(tmp_phi_l_(1+nl))*ring_k_c_3z__;
tmp_ring_tru_k_c_3z__ = Rz(tmp_psi_n_(1+nl))*Ry(tmp_the_n_(1+nl))*Rz(tmp_phi_n_(1+nl))*ring_k_c_3z__;
tmp_ring_l2 = sum((tmp_ring_est_k_c_3z__ - tmp_ring_tru_k_c_3z__).^2,'all')*2*pi/max(1,n_gamma_z);
tmp_err_psi_the_phi_l_(1+nl) = sqrt(tmp_ring_l2);
end;%for nl=0:n_l-1;
tmp_err_l_ = min(tmp_err_phi_the_psi_l_,tmp_err_psi_the_phi_l_);
err_phi_the_psi_10_r_ = tmp_err_phi_the_psi_l_;
err_psi_the_phi_10_r_ = tmp_err_psi_the_phi_l_;
err_10_r_ = tmp_err_l_;
%%%%%%%%;
n_l = size(Array_9_3_50_rc__,1);
tmp_nx_l_ = Array_9_3_50_rc__(:,1+nc_50_X);
tmp_ny_l_ = Array_9_3_50_rc__(:,1+nc_50_Y);
tmp_mip_l_ = mip_50_r_;
tmp_CoC_l_ = CoC_50_r_;
tmp_psi_l_ = Array_9_3_50_rc__(:,1+nc_50_Psi);
tmp_the_l_ = Array_9_3_50_rc__(:,1+nc_50_Theta);
tmp_phi_l_ = Array_9_3_50_rc__(:,1+nc_50_Phi);
tmp_psi_l_ = tmp_psi_l_*pi/180.0;
tmp_the_l_ = tmp_the_l_*pi/180.0;
tmp_phi_l_ = tmp_phi_l_*pi/180.0;
[tmp_ij_,dxy_50_r_] = knnsearch([nx_c__(:),ny_c__(:)],[tmp_nx_l_,tmp_ny_l_]); tmp_index_ = tmp_ij_ - 1;
tmp_psi_n_ = reshape(psi_c__(1+tmp_index_),[n_l,1]);
tmp_the_n_ = reshape(the_c__(1+tmp_index_),[n_l,1]);
tmp_phi_n_ = reshape(phi_c__(1+tmp_index_),[n_l,1]);
tmp_err_phi_the_psi_l_ = zeros(n_l,1);
tmp_err_psi_the_phi_l_ = zeros(n_l,1);
tmp_err_l_ = zeros(n_l,1);
for nl=0:n_l-1;
tmp_ring_est_k_c_3z__ = Rz(tmp_phi_l_(1+nl))*Ry(tmp_the_l_(1+nl))*Rz(tmp_psi_l_(1+nl))*ring_k_c_3z__;
tmp_ring_tru_k_c_3z__ = Rz(tmp_phi_n_(1+nl))*Ry(tmp_the_n_(1+nl))*Rz(tmp_psi_n_(1+nl))*ring_k_c_3z__;
tmp_ring_l2 = sum((tmp_ring_est_k_c_3z__ - tmp_ring_tru_k_c_3z__).^2,'all')*2*pi/max(1,n_gamma_z);
tmp_err_phi_the_psi_l_(1+nl) = sqrt(tmp_ring_l2);
tmp_ring_est_k_c_3z__ = Rz(tmp_psi_l_(1+nl))*Ry(tmp_the_l_(1+nl))*Rz(tmp_phi_l_(1+nl))*ring_k_c_3z__;
tmp_ring_tru_k_c_3z__ = Rz(tmp_psi_n_(1+nl))*Ry(tmp_the_n_(1+nl))*Rz(tmp_phi_n_(1+nl))*ring_k_c_3z__;
tmp_ring_l2 = sum((tmp_ring_est_k_c_3z__ - tmp_ring_tru_k_c_3z__).^2,'all')*2*pi/max(1,n_gamma_z);
tmp_err_psi_the_phi_l_(1+nl) = sqrt(tmp_ring_l2);
end;%for nl=0:n_l-1;
tmp_err_l_ = min(tmp_err_phi_the_psi_l_,tmp_err_psi_the_phi_l_);
err_phi_the_psi_50_r_ = tmp_err_phi_the_psi_l_;
err_psi_the_phi_50_r_ = tmp_err_psi_the_phi_l_;
err_50_r_ = tmp_err_l_;
%%%%%%%%;
clear tmp_nx_l_ tmp_ny_l_ tmp_mip_l_ tmp_CoC_l_ tmp_psi_l_ tmp_the_l_ tmp_phi_l_ tmp_psi_l_ tmp_the_l_ tmp_phi_l_ tmp_ij_ tmp_psi_n_ tmp_the_n_ tmp_phi_n_ tmp_err_phi_the_psi_l_ tmp_err_psi_the_phi_l_ tmp_err_l_ tmp_ring_est_k_c_3z__ tmp_ring_tru_k_c_3z__ tmp_ring_l2 ;
%%%%%%%%;

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Show the z-score vs CoC scatterplot, this time colored by distances. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,768*2]);
markersize_dot = 16;
markersize_cir = 10;
fontsize_use = 12;
%c_use__ = flipud(colormap_80s));
c_use__ = colormap('cool');
%%%%;
subplot(2,1,1);
hold on;
scatter(mip_50_r_,CoC_50_r_,markersize_dot,err_50_r_,'filled');
colormap(c_use__);
%tmp_index_0_ = efind(lab_50_r_==0);
%plot(mip_50_r_(1+tmp_index_0_),CoC_50_r_(1+tmp_index_0_),'o','Color',0.35*[1,1,1],'MarkerSize',markersize_cir);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([5,9]);ylim([-0.01,+0.03]); grid on;
colorbar;
%legend({'10 xx','10 fp','10 tp'},'Location','SouthEast');
legend({'10 xx','10 tp'},'Location','SouthEast');
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_50 euler-error','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(2,1,2);
hold on;
scatter(mip_50_r_,CoC_50_r_,markersize_dot,dxy_50_r_,'filled');
colormap(c_use__);
%tmp_index_0_ = efind(lab_50_r_==0);
%plot(mip_50_r_(1+tmp_index_0_),CoC_50_r_(1+tmp_index_0_),'o','Color',0.35*[1,1,1],'MarkerSize',markersize_cir);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([5,9]);ylim([-0.01,+0.03]); grid on;
colorbar;
%legend({'10 xx','10 fp','10 tp'},'Location','SouthEast');
legend({'10 xx','10 tp'},'Location','SouthEast');
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_50 xy-error','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
fname_fig_pre = sprintf('%s/test_Kexin_mip_vs_CoC_9_3_50_FIGF',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Show the z-score vs CoC scatterplot, this time colored by distances. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;set(gcf,'Position',1+[0,0,1024,768*2]);
markersize_dot = 16;
markersize_cir = 10;
fontsize_use = 12;
%c_use__ = flipud(colormap_80s));
c_use__ = colormap('cool');
%%%%;
subplot(2,1,1);
hold on;
scatter(mip_10_r_,CoC_10_r_,markersize_dot,err_10_r_,'filled');
colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([5,9]);ylim([-0.01,+0.03]); grid on;
colorbar;
legend({'10 xx','50 tp'},'Location','SouthEast');
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_10 euler-error','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(2,1,2);
hold on;
scatter(mip_10_r_,CoC_10_r_,markersize_dot,dxy_10_r_,'filled');
colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_Peak),Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([5,9]);ylim([-0.01,+0.03]); grid on;
colorbar;
legend({'10 xx','50 tp'},'Location','SouthEast');
xlabel('z-score (mip)'); ylabel('CoC');
title('9_3_10 xy-error','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
fname_fig_pre = sprintf('%s/test_Kexin_mip_vs_CoC_9_3_10_FIGF',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Show err_r_ vs pose_sim. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figsml;
markersize_dot = 8;
markersize_cir = 6;
fontsize_use = 12;
c_use__ = colormap('cool');
%%%%;
subplot(1,1,1);
hold on;
scatter(Array_9_3_50_rc__(:,1+nc_50_pose_sim),err_50_r_,markersize_dot,dxy_50_r_,'filled');
colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_pose_sim),err_50_r_(1+tmp_index_1_),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([-0.21,+1.01]); ylim([0,2*pi]);grid on;
colorbar;
legend({'50 xx','50 tp'},'Location','NorthEast');
xlabel('pose_sim','Interpreter','none'); ylabel('euler_error','Interpreter','none');
title('9_3_50 pose_sim vs euler_error (colored by dxy_50_r_)','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_Kexin_pose_sim_vs_err_9_3_10_FIGH',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
%plot3(mip_50_r_,CoC_50_r_,dxy_50_r_,'.'); axis vis3d; xlabel('z-score');ylabel('CoC');zlabel('dxy');
%plot3(mip_50_r_,CoC_50_r_,err_50_r_,'.'); axis vis3d; xlabel('z-score');ylabel('CoC');zlabel('euler-error');
%plot3(mip_10_r_,CoC_10_r_,err_10_r_,'.'); axis vis3d; xlabel('z-score');ylabel('CoC');zlabel('euler-error');
%scatter3(dxy_10_r_,err_10_r_,mip_10_r_,12,CoC_10_r_); axis vis3d; xlabel('dxy');ylabel('err');zlabel('Peak'); title('colored by CoC');
%scatter3(dxy_10_r_,err_10_r_,CoC_10_r_,12,mip_10_r_); axis vis3d; xlabel('dxy');ylabel('err');zlabel('CoC'); title('colored by Peak');

flag_disp=0;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
% Show the CoC vs dxy_r_ scatterplot, this time colored by err_r_. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figbig;
markersize_dot = 8;
markersize_cir = 6;
fontsize_use = 12;
c_use__ = colormap('cool');
%%%%;
subplot(2,1,1);
hold on;
scatter(CoC_10_r_,dxy_10_r_,markersize_dot,err_10_r_,'filled');
colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),dxy_50_r_(1+tmp_index_1_),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([-0.01,+0.03]); ylim([-10,375]);grid on;
colorbar;
legend({'10 xx','50 tp'},'Location','NorthEast');
xlabel('CoC'); ylabel('dxy');
title('9_3_10 CoC vs dxy (colored by err_10_r_)','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(2,1,2);
hold on;
scatter(CoC_50_r_,dxy_50_r_,markersize_dot,err_50_r_,'filled');
colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),dxy_50_r_(1+tmp_index_1_),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([-0.01,+0.03]); ylim([-10,375]);grid on;
colorbar;
legend({'50 xx','50 tp'},'Location','NorthEast');
xlabel('CoC'); ylabel('dxy');
title('9_3_50 CoC vs dxy (colored by err_50_r_)','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%%%%%;
fname_fig_pre = sprintf('%s/test_Kexin_CoC_vs_dxy_9_3_10_FIGG',dir_jpg);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg));
disp(sprintf(' %% Writing %s',fname_fig_pre));
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg));
close(gcf);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% quick test of erf. ;
% Note: CDF of gaussian = 0.5 + 0.5*erf(x/sqrt(2));
% also: CDFinv(z) = sqrt(2)*erfinv(2*z-1);
%%%%%%%%;
flag_disp=0;
if flag_disp;
tmp_x_ = linspace(-3,+3,1024);
tmp_e_ = 0.5*(1+erf(tmp_x_/sqrt(2)));
tmp_g_ = 1/sqrt(2*pi) * exp(-tmp_x_.^2/2);
tmp_de_ = (tmp_e_(2:end)-tmp_e_(1:end-1))/max(1e-12,mean(diff(tmp_x_)));
tmp_x2_ = 0.5*tmp_x_(2:end) + 0.5*tmp_x_(1:end-1);
tmp_y3_ = linspace(0,1,1024);
tmp_x3_ = erfinv(2*tmp_y3_-1)*sqrt(2);
figure(1+nf);nf=nf+1;clf;figmed;
subplot(1,2,1); plot(tmp_x2_,tmp_de_,'rx',tmp_x_,tmp_g_,'go'); legend({'derf/dx','g(x)'});
subplot(1,2,2); plot(tmp_x_,tmp_e_,'rx',tmp_x3_,tmp_y3_,'go'); legend({'erf(x)','erfinv'});
end;%if flag_disp;

%%%%%%%%;
% Now apply probit to both the Peak and the CoC. ;
% For now we apply probit globally (to bulk and outliers). ;
%%%%%%%%;
[~,rank_mip_50_r_] = sort(mip_50_r_,'ascend'); [~,rank_mip_50_r_] = sort(rank_mip_50_r_,'ascend'); rank_mip_50_r_ = (rank_mip_50_r_-0.5)/max(1,numel(mip_50_r_)); pro_mip_50_r_ = sqrt(2.0)*erfinv(2.0*rank_mip_50_r_ - 1.0);
[~,rank_CoC_50_r_] = sort(CoC_50_r_,'ascend'); [~,rank_CoC_50_r_] = sort(rank_CoC_50_r_,'ascend'); rank_CoC_50_r_ = (rank_CoC_50_r_-0.5)/max(1,numel(CoC_50_r_)); pro_CoC_50_r_ = sqrt(2.0)*erfinv(2.0*rank_CoC_50_r_ - 1.0);
[~,rank_mip_10_r_] = sort(mip_10_r_,'ascend'); [~,rank_mip_10_r_] = sort(rank_mip_10_r_,'ascend'); rank_mip_10_r_ = (rank_mip_10_r_-0.5)/max(1,numel(mip_10_r_)); pro_mip_10_r_ = sqrt(2.0)*erfinv(2.0*rank_mip_10_r_ - 1.0);
[~,rank_CoC_10_r_] = sort(CoC_10_r_,'ascend'); [~,rank_CoC_10_r_] = sort(rank_CoC_10_r_,'ascend'); rank_CoC_10_r_ = (rank_CoC_10_r_-0.5)/max(1,numel(CoC_10_r_)); pro_CoC_10_r_ = sqrt(2.0)*erfinv(2.0*rank_CoC_10_r_ - 1.0);
figure(1+nf);nf=nf+1;clf;figbig;
markersize_dot = 8;
markersize_cir = 6;
fontsize_use = 12;
c_use__ = colormap('cool');
%%%%;
subplot(2,2,1);
hold on;
scatter(mip_50_r_,CoC_50_r_,markersize_dot,err_50_r_,'filled'); colormap(c_use__);
tmp_index_1_ = efind(lab_50_r_==1);
plot(Array_9_3_50_rc__(1+tmp_index_1_,1+nc_50_CoC),dxy_50_r_(1+tmp_index_1_),'o','Color',c_use__(1,:),'MarkerSize',markersize_cir);
hold off;
xlim([-0.01,+0.03]); ylim([-10,375]);grid on;
colorbar;
legend({'10 xx','50 tp'},'Location','NorthEast');
xlabel('CoC'); ylabel('dxy');
title('9_3_10 CoC vs dxy (colored by err_10_r_)','Interpreter','none');
set(gca,'FontSize',fontsize_use);
%%%%;
subplot













