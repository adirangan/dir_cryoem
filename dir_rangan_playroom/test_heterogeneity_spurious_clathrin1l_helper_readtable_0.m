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

%%%%%%%%;
% determine labels. ;
%%%%%%%%;
n_50_r = size(Array_9_3_50_rc__,1);
lab_50_r_ = Array_9_3_50_rc__(:,1+nc_50_labels);
mip_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Peak);
CoC_50_r_ = Array_9_3_50_rc__(:,1+nc_50_CoC);
nx_50_r_ = Array_9_3_50_rc__(:,1+nc_50_X);
ny_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Y);
n_10_r = size(Array_9_3_10_rc__,1);
mip_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Peak);
CoC_10_r_ = Array_9_3_10_rc__(:,1+nc_10_CoC);
nx_10_r_ = Array_9_3_10_rc__(:,1+nc_10_X);
ny_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Y);
tmp_index_1_ = efind(lab_50_r_==1);
[tmp_ij_] = knnsearch([mip_10_r_,CoC_10_r_],[mip_50_r_(1+tmp_index_1_),CoC_50_r_(1+tmp_index_1_)]);
lab_10_r_ = zeros(n_10_r,1); lab_10_r_(tmp_ij_) = 1; clear tmp_ij_;
%%%%;
psi_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Psi)*pi/180;
the_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Theta)*pi/180;
phi_10_r_ = Array_9_3_10_rc__(:,1+nc_10_Phi)*pi/180;
%%%%;
psi_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Psi)*pi/180;
the_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Theta)*pi/180;
phi_50_r_ = Array_9_3_50_rc__(:,1+nc_50_Phi)*pi/180;
%%%%%%%%;

flag_disp=1;
if flag_disp;
%%%%%%%%;
% distribution of estimated viewing-angles. ;
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;figmed;
[ ...
 n_viewing_S ...
,viewing_azimu_b_S_ ...
,viewing_polar_a_S_ ...
,viewing_weight_S_ ...
,viewing_k_c_0_S_ ...
,viewing_k_c_1_S_ ...
,viewing_k_c_2_S_ ...
,viewing_n_polar_a ...
,viewing_polar_a_ ...
,viewing_n_azimu_b_ ...
] = ...
sample_shell_6( ...
 1.0 ...
,0.75/(2*pi) ...
) ;
%%%%%%%%;
subplot(1,2,1);
hlim_ = [];
c_use__ = colormap_80s;
flag_2d_vs_3d = 0;
flag_loghist_vs_hist = 0;
[ ...
 h_raw_10_ ...
 h_w3d_10_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,the_10_r_ ...
,phi_10_r_ ...
,hlim_ ...
,c_use__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
%%%%%%%%;
subplot(1,2,2);
hlim_ = [];
c_use__ = colormap_80s;
flag_2d_vs_3d = 1;
flag_loghist_vs_hist = 0;
[ ...
 h_raw_10_ ...
 h_w3d_10_ ...
] = ...
hist2d_polar_a_azimu_b_0( ...
 viewing_polar_a_S_ ...
,viewing_azimu_b_S_ ...
,viewing_weight_S_ ...
,the_10_r_ ...
,phi_10_r_ ...
,hlim_ ...
,c_use__ ...
,flag_2d_vs_3d ...
,flag_loghist_vs_hist ...
);
%%%%%%%%;
fname_fig_pre = sprintf('%s_jpg/test_heterogeneity_spurious_clathrin1l_readtable_FIGA_',dir_pm);
fname_fig_jpg = sprintf('%s.jpg',fname_fig_pre);
fname_fig_jpg_stripped = sprintf('%s_stripped.jpg',fname_fig_pre);
if (flag_replot | ~exist(fname_fig_jpg,'file'));
disp(sprintf(' %% %s not found, creating',fname_fig_jpg));
print('-djpeg',fname_fig_jpg_stripped);
sgtitle(fname_fig_pre,'Interpreter','none');
print('-djpeg',fname_fig_jpg);
end;%if (flag_replot | ~exist(fname_fig_jpg,'file'));
%%%%%%%%;
close(gcf);
end;%if flag_disp;
