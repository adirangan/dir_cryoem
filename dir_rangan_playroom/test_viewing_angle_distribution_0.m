%%%%%%%%;
platform = 'rusty';
if (exist('platform.type','file')); fp=fopen('platform.type'); platform = fscanf(fp,'%s'); fclose(fp); end;
if (strcmp(platform,'access1')); setup_access1; string_root = 'data'; end;
if (strcmp(platform,'OptiPlex')); setup_OptiPlex; string_root = 'home'; end;
if (strcmp(platform,'eval1')); setup_eval1; string_root = 'home'; end;
if (strcmp(platform,'rusty')); setup_rusty; string_root = 'mnt/home'; end;
%%%%%%%%;
dir_base = sprintf('/%s/rangan/dir_cryoem',string_root);
flag_replot = 0;
flag_replot_vol = 0;
flag_verbose = 1; nf=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
table_data__ = { ...
%'p28hRPT1_x0' , 'p28hRP' , 0.98 , 'emd_8674.map' , 'T1.star' ; ...
'ISWINCP_x0' , 'ISWINCP' , 1.07 , 'emd_9718.map' , 'ADPBeF.star' ; ...
'trpv1_x0' , 'trpv1' , 1.2156 , 'emd_5778.mrc' , 'tv1_relion_data.star' ; ...
'rib80s_x0' , 'rib80s' , 1.34 , 'emd_2660.mrc' , 'shiny_2sets.star' ; ...
'MlaFEDB_x0' , 'MlaFEDB' , 0.832 , 'emd_22116.map' , 'Empiar_10536_00_to_23.star' ; ...
'LetB1_x0' , 'LetB1' , 1.31 , 'emd_20993.map' , 'job_569_model_1.star' ; ...
'TMEM16F_x0' , 'TMEM16F' , 1.059 , 'emd_20244.map' , 'All_8192.star' ; ...
'LSUbl17dep_x0' , 'LSUbl17dep' , 1.31 , 'emd_8434.map' , 'Parameters_negated.star' ; ...
'ps1_x0' , 'precatalytic_spliceosome' , 1.699 , 'consensus_half1_class001.mrc' , 'consensus_data.star' ; ...
};
n_experiment = size(table_data__,1);
ampm__ = cell(n_experiment,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
tmp_ = clock;rng(tmp_(end));
for nexperiment=(randperm(n_experiment)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
na=0;
fname_prefix = table_data__{1+nexperiment,1+na}; na=na+1;
dir_nopath_data_star = table_data__{1+nexperiment,1+na}; na=na+1;
Pixel_Spacing = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_volume = table_data__{1+nexperiment,1+na}; na=na+1;
fname_nopath_star = table_data__{1+nexperiment,1+na}; na=na+1;
if (flag_verbose); 
disp(sprintf(' %% nexperiment %d/%d: %16s %16s %0.3f %16s %32s' ...
	     ,nexperiment,n_experiment ...
	     ,fname_prefix,dir_nopath_data_star,Pixel_Spacing,fname_nopath_volume,fname_nopath_star ...
	     ));
end;%if (flag_verbose); 
global_parameter = struct('type','parameter');
global_parameter.flag_replot=flag_replot;
global_parameter.flag_replot_vol=flag_replot_vol;
global_parameter.flag_invert=0;
global_parameter.flag_center_image=0;
global_parameter.tolerance_master = 1e-2;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_invert = 0; end;
if (strcmp(dir_nopath_data_star,'LSUbl17dep')); global_parameter.flag_center_image = 1; end;
if (strcmp(dir_nopath_data_star,'precatalytic_spliceosome')); global_parameter.flag_center_image = 1; end;
dir_pm = sprintf('%s/dir_%s/dir_pm',dir_base,fname_prefix);
if (flag_verbose); disp(sprintf(' %% fname_prefix: %s',fname_prefix)); end;
[~,ampm_] = ampm_fromdisk_3([],dir_pm);
ampm__{1+nexperiment} = ampm_;

tolerance_master = global_parameter.tolerance_master;
flag_replot = global_parameter.flag_replot;
flag_replot_vol = global_parameter.flag_replot_vol;
N_k_p_use__ = ampm_.M_k_p__; if (global_parameter.flag_center_image==1); N_k_p_use__ = ampm_.N_k_p__; end;
tmp_delta_r_max = 0.030;
tmp_delta_r_upb = 0.250;
n_M = ampm_.n_M;
n_iteration_orig = size(ampm_.corr_a_k_Y_i_,1);
n_iteration_seco = max(1,min(n_iteration_orig-2,6));

%%%%%%%%;
% sort. ;
%%%%%%%%;
[~,ij_X_best_ampm_srt_] = sort(ampm_.X_best_ampm_ia__(end,:),'descend'); index_X_best_ampm_srt_ = ij_X_best_ampm_srt_-1;

%%%%%%%%;
% Print out some of the volumes. ;
%%%%%%%%;
flag_disp = 1;
if flag_disp;
%%%%%%%%;
fname_vol_a_fig_pre = sprintf('%s_jpg/%s_X_2d_xcor_d0_a1t_vol_a_FIGE',ampm_.dir_pm,dir_nopath_data_star);
fname_vol_a_fig_jpg = sprintf('%s.jpg',fname_vol_a_fig_pre);
if (flag_replot_vol | ~exist(fname_vol_a_fig_jpg,'file'));
%%%%%%%%;
test_viewing_angle_distribution_vol_a_helper_0; %<-- construct a_x_u_ampm_. ;
%%%%%%%%;
if strfind(ampm_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(ampm_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(ampm_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(ampm_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
%percent_threshold_ = 85:1.25:97.5;
n_percent_threshold = numel(percent_threshold_);
%%%%;
n_x_u_pack = ampm_.n_x_u_pack;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
fontsize_use = 12;
if (n_percent_threshold> 1); p_row = tmp_n_a; p_col = n_percent_threshold; end;
if (n_percent_threshold==1); p_row = 4; p_col = ceil(tmp_n_a/p_row); end;
for npt=0:n_percent_threshold-1;
tmp_percent_threshold = percent_threshold_(1+npt); nnp=0;
for tmp_na=0:tmp_n_a-1;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(tmp_a_x_u_xxxa__(1+tmp_index_,1+tmp_na),tmp_percent_threshold);
%title(sprintf('Z %.3f',tmp_corr_a_(1+tmp_na))); 
xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
end;%for tmp_na=0:tmp_n_a-1;
%%%%;
end;%for npt=0:n_percent_threshold-1;
%%%%%%%%
set(gcf,'Position',1+[0,0,1024+512,1024]);
sgtitle(fname_vol_a_fig_pre,'Interpreter','none');
print('-djpeg',fname_vol_a_fig_jpg);
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_vol_stack',string_root);
fname_vol_a_fig_jpg_strip = sprintf('%s/%s_X_2d_xcor_d0_a1t_vol_a_FIGE_strip.jpg',tmp_dir,dir_nopath_data_star);
disp(sprintf(' %% writing %s',fname_vol_a_fig_jpg_strip));
print('-djpeg',fname_vol_a_fig_jpg_strip);
close gcf;
%%%%%%%%;
end;%if (flag_replot_vol | ~exist(fname_vol_a_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%;
% Choose the na. ;
%%%%%%%%;
na_ = index_X_best_ampm_srt_([1,2,round(numel(ij_X_best_ampm_srt_)/2),round(numel(ij_X_best_ampm_srt_)/2)+1]); %<-- take top 2 best as well as 2 near the median. ;
n_na = numel(na_);
%%%%%%%%%%%%%%%%;
for nna=0:n_na-1;
%%%%%%%%%%%%%%%%;
na = na_(1+nna);
str_fname_mat = ampm_.str_fname_mat_a_{1+na};
str_fname_mat_pre = str_fname_mat(1:strfind(str_fname_mat,'.mat')-1);
str_fname_mat_inf = str_fname_mat_pre(strfind(str_fname_mat_pre,'X_2d_xcor_d0'):end);
fname_prefix_str_fname_mat_inf = sprintf('%s_%s',fname_prefix,str_fname_mat_inf);
fname_pre = sprintf('%s_viewing_angle_distribution',str_fname_mat_pre);
[flag_skip,fname_mat] = open_fname_tmp(fname_pre);
if (flag_verbose); disp(sprintf(' %% fname_pre: %s <-- flag_skip %d',fname_pre,flag_skip)); end;
if ~flag_skip;

%%%%%%%%;
% First cluster the CTF based on tolerance_cluster. ;
%%%%%%%%;
parameter_CTF = struct('type','parameter');
parameter_CTF.tolerance_master = tolerance_master;
[ ...
 parameter_CTF ...
,ampm_.index_ncluster_from_nCTF_ ...
] = ...
knn_cluster_CTF_k_p_r_kC__0( ...
 parameter_CTF ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.weight_2d_k_p_r_ ...
,ampm_.n_CTF ...
,ampm_.CTF_k_p_r_kC__ ...
);
%%%%%%%%;
ampm_.n_cluster = 1+max(ampm_.index_ncluster_from_nCTF_);

%%%%%%%%;
% Now align M_k_p__ to a_k_Y_quad_. ;
%%%%%%%%;
parameter_align = struct('type','parameter');
parameter_align.fname_align_a_k_Y_pre = ''; %<-- empty for now. ;
parameter_align.delta_r_max = tmp_delta_r_max;
parameter_align.n_delta_v_requested = 24;
parameter_align.delta_r_upb = tmp_delta_r_upb;
parameter_align.n_iteration = n_iteration_seco;
[ ...
 parameter_align ...
,a_k_Y_reco_yki__ ...
,corr_a_k_Y_i_ ...
,euler_polar_a_Mi__ ...
,euler_azimu_b_Mi__ ...
,euler_gamma_z_Mi__ ...
,image_delta_x_acc_Mi__ ...
,image_delta_y_acc_Mi__ ...
,image_delta_x_upd_Mi__ ...
,image_delta_y_upd_Mi__ ...
,flag_image_delta_upd_Mi__ ...
,image_I_value_Mi__ ...
,image_X_value_Mi__ ...
,image_S_index_Mi__ ...
,n_S ...
,template_viewing_azimu_b_all_ ...
,template_viewing_polar_a_all_ ...
,X_SMi___ ...
,delta_x_SMi___ ...
,delta_y_SMi___ ...
,gamma_z_SMi___ ...
,I_value_SMi___ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter_align ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.k_p_r_max ...
,ampm_.weight_3d_k_p_r_ ...
,ampm_.weight_2d_k_p_r_ ...
,ampm_.n_w_ ...
,ampm_.n_M ...
,N_k_p_use__ ...
,ampm_.n_cluster ...
,ampm_.index_ncluster_from_nCTF_ ...
,ampm_.n_CTF ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_r_kC__ ...
,ampm_.l_max_ ...
,ampm_.a_k_Y_quad_ ...
,ampm_.euler_polar_a_Mi__(:,end-(parameter_align.n_iteration-1)) ...
,ampm_.euler_azimu_b_Mi__(:,end-(parameter_align.n_iteration-1)) ...
,ampm_.euler_gamma_z_Mi__(:,end-(parameter_align.n_iteration-1)) ...
,ampm_.image_delta_x_acc_Mi__(:,end-(parameter_align.n_iteration-1)) ...
,ampm_.image_delta_y_acc_Mi__(:,end-(parameter_align.n_iteration-1)) ...
);
if (flag_verbose); disp(sprintf(' %% corr(a_k_Y_reco_yki__(:,end),ampm_.a_k_Y_reco_yki__(:,end)): %0.16f',corr(a_k_Y_reco_yki__(:,end),ampm_.a_k_Y_reco_yki__(:,end)))); end;
flag_disp=0;
if flag_disp;
figure(1+nf);nf=nf+1;;clf;figbig;
p_row=2; p_col=ceil(n_iteration_orig/p_row); np=0;
for nl=0:n_iteration_orig-1;
subplot(p_row,p_col,1+np);np=np+1;
plot(ampm_.euler_azimu_b_Mi__(:,1+nl),ampm_.euler_polar_a_Mi__(:,1+nl),'.');
xlabel('azimu_b','Interpreter','none'); xlim([0,2*pi]);
ylabel('polar_a','Interpreter','none'); ylim([0,1*pi]);
grid on;
end;%for nl=0:n_iteration_orig-1;
%%%%;
figure(1+nf);nf=nf+1;;clf;figmed;
p_row=2; p_col=ceil(n_iteration_seco/p_row); np=0;
for nl=0:n_iteration_seco-1;
subplot(p_row,p_col,1+np);np=np+1;
plot(euler_azimu_b_Mi__(:,1+nl),euler_polar_a_Mi__(:,1+nl),'.');
xlabel('azimu_b','Interpreter','none'); xlim([0,2*pi]);
ylabel('polar_a','Interpreter','none'); ylim([0,1*pi]);
grid on;
end;%for nl=0:n_iteration_seco-1;
%%%%;
figure(1+nf);nf=nf+1;;clf;figbig;
p_row=2; p_col=ceil(n_iteration_orig/p_row); np=0;
for nl=0:n_iteration_orig-1;
subplot(p_row,p_col,1+np);np=np+1;
plot(ampm_.image_delta_x_acc_Mi__(:,1+nl),ampm_.image_delta_y_acc_Mi__(:,1+nl),'.');
xlabel('delta_x','Interpreter','none'); xlim(tmp_delta_r_upb*[-1,+1]);
ylabel('delta_y','Interpreter','none'); ylim(tmp_delta_r_upb*[-1,+1]);
grid on;
end;%for nl=0:n_iteration_orig-1;
%%%%;
figure(1+nf);nf=nf+1;;clf;figmed;
p_row=2; p_col=ceil(n_iteration_seco/p_row); np=0;
for nl=0:n_iteration_seco-1;
subplot(p_row,p_col,1+np);np=np+1;
plot(image_delta_x_acc_Mi__(:,1+nl),image_delta_y_acc_Mi__(:,1+nl),'.');
xlabel('delta_x','Interpreter','none'); xlim(tmp_delta_r_upb*[-1,+1]);
ylabel('delta_y','Interpreter','none'); ylim(tmp_delta_r_upb*[-1,+1]);
grid on;
end;%for nl=0:n_iteration_seco-1;
end;%if flag_disp;

%%%%%%%%;
% Now align M_k_p__ to the chosen a_k_Y_ampm_. ;
%%%%%%%%;
a_k_Y_ampm_ = ampm_.a_k_Y_ampm_yka__(:,1+na);
[ ... 
 a_k_Y_ampm_alig_ ...
,corr_a_k_Y_quad_vs_a_k_Y_ampm_alig ...
] = ...
spharm_register_and_rotate_2( ...
 ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.k_p_r_max ...
,ampm_.weight_3d_k_p_r_ ...
,ampm_.l_max_ ...
,ampm_.a_k_Y_quad_ ...
,ampm_.a_k_Y_ampm_yka__(:,1+na) ...
,0 ...
,ampm_.X_best_ampm_ia__(end,1+na) ...
,ampm_.X_best_flag_flip_ampm_ia__(end,1+na) ...
,ampm_.polar_a_best_ampm_ia__(end,1+na) ...
,ampm_.azimu_b_best_ampm_ia__(end,1+na) ...
,ampm_.gamma_z_best_ampm_ia__(end,1+na) ...
,ampm_.delta_best_ampm_dia__(1+0,end,1+na) ...
,ampm_.delta_best_ampm_dia__(1+1,end,1+na) ...
,ampm_.delta_best_ampm_dia__(1+2,end,1+na) ...
);
%%%%;
tmp_viewing_euler_gamma_z_ = ampm_.euler_gamma_z_ampm_Ma__(:,1+na);
if (ampm_.X_best_flag_flip_ampm_ia__(end,1+na));
[ ...
 tmp_viewing_euler_gamma_z_ ...
] = ...
rotate_viewing_angles_from_invert_spharm_0( ...
 flag_verbose ...
,tmp_viewing_euler_gamma_z_ ...
);
end;%if (ampm_.X_best_flag_flip_ampm_ia__(end,1+na));
%%%%;
[ ...
 ~ ...
,tmp_viewing_euler_polar_a_ ...
,tmp_viewing_euler_azimu_b_ ...
,tmp_viewing_euler_gamma_z_ ...
] = ...
rotate_viewing_angles_from_rotate_spharm_to_spharm_0( ...
 max(0,flag_verbose-1) ...
,ampm_.polar_a_best_ampm_ia__(end,1+na) ...
,ampm_.azimu_b_best_ampm_ia__(end,1+na) ...
,ampm_.gamma_z_best_ampm_ia__(end,1+na) ...
,ampm_.euler_polar_a_ampm_Ma__(:,1+na) ...
,ampm_.euler_azimu_b_ampm_Ma__(:,1+na) ...
,tmp_viewing_euler_gamma_z_ ...
);

parameter_align = struct('type','parameter');
parameter_align.fname_align_a_k_Y_pre = ''; %<-- empty for now. ;
parameter_align.delta_r_max = tmp_delta_r_max;
parameter_align.n_delta_v_requested = 24;
parameter_align.delta_r_upb = tmp_delta_r_upb;
parameter_align.n_iteration = n_iteration_seco;
[ ...
 parameter_align ...
,a_k_Y_ampm_alig_reco_yki__ ...
,corr_a_k_Y_ampm_alig_vs_a_k_Y_ampm_alig_reco_yk_i_ ...
,euler_polar_a_ampm_alig_Mi__ ...
,euler_azimu_b_ampm_alig_Mi__ ...
,euler_gamma_z_ampm_alig_Mi__ ...
,image_delta_x_acc_ampm_alig_Mi__ ...
,image_delta_y_acc_ampm_alig_Mi__ ...
,image_delta_x_upd_ampm_alig_Mi__ ...
,image_delta_y_upd_ampm_alig_Mi__ ...
,flag_image_delta_upd_ampm_alig_Mi__ ...
,image_I_value_ampm_alig_Mi__ ...
,image_X_value_ampm_alig_Mi__ ...
,image_S_index_ampm_alig_Mi__ ...
,~ ...
,~ ...
,~ ...
,X_ampm_alig_SMi___ ...
,delta_x_ampm_alig_SMi___ ...
,delta_y_ampm_alig_SMi___ ...
,gamma_z_ampm_alig_SMi___ ...
,I_value_ampm_alig_SMi___ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter_align ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.k_p_r_max ...
,ampm_.weight_3d_k_p_r_ ...
,ampm_.weight_2d_k_p_r_ ...
,ampm_.n_w_ ...
,ampm_.n_M ...
,N_k_p_use__ ...
,ampm_.n_cluster ...
,ampm_.index_ncluster_from_nCTF_ ...
,ampm_.n_CTF ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_r_kC__ ...
,ampm_.l_max_ ...
,a_k_Y_ampm_alig_ ...
,tmp_viewing_euler_polar_a_ ...
,tmp_viewing_euler_azimu_b_ ...
,tmp_viewing_euler_gamma_z_ ...
,ampm_.image_delta_x_ampm_Ma__(:,1+na) ...
,ampm_.image_delta_y_ampm_Ma__(:,1+na) ...
);
if (flag_verbose); disp(sprintf(' %% corr(a_k_Y_ampm_alig_,a_k_Y_ampm_alig_reco_yki__(:,end)): %0.16f',corr(a_k_Y_ampm_alig_,a_k_Y_ampm_alig_reco_yki__(:,end)))); end;

[~,corr_a_k_Y_quad_vs_a_k_Y_ampm_alig_reco] = ...
register_spharm_to_spharm_3( ...
 max(0,flag_verbose-1) ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.weight_3d_k_p_r_ ...
,ampm_.l_max_ ...
,ampm_.a_k_Y_quad_ ...
,a_k_Y_ampm_alig_reco_yki__(:,end) ...
);
if (flag_verbose); disp(sprintf(' %% ampm_.X_best_ampm_ia__(end,1+na):        %0.16f',ampm_.X_best_ampm_ia__(end,1+na))); end;
if (flag_verbose); disp(sprintf(' %% corr_a_k_Y_quad_vs_a_k_Y_ampm_alig:      %0.16f',real(corr_a_k_Y_quad_vs_a_k_Y_ampm_alig     ))); end;
if (flag_verbose); disp(sprintf(' %% corr_a_k_Y_quad_vs_a_k_Y_ampm_alig_reco: %0.16f',real(corr_a_k_Y_quad_vs_a_k_Y_ampm_alig_reco))); end;

parameter_align = struct('type','parameter');
parameter_align.fname_align_a_k_Y_pre = ''; %<-- empty for now. ;
parameter_align.delta_r_max = tmp_delta_r_max;
parameter_align.n_delta_v_requested = 24;
parameter_align.delta_r_upb = tmp_delta_r_upb;
parameter_align.n_iteration = n_iteration_seco;
[ ...
 parameter_align ...
,a_k_Y_ampm_reco_yki__ ...
,corr_a_k_Y_ampm_vs_a_k_Y_ampm_reco_yk_i_ ...
,euler_polar_a_ampm_Mi__ ...
,euler_azimu_b_ampm_Mi__ ...
,euler_gamma_z_ampm_Mi__ ...
,image_delta_x_acc_ampm_Mi__ ...
,image_delta_y_acc_ampm_Mi__ ...
,image_delta_x_upd_ampm_Mi__ ...
,image_delta_y_upd_ampm_Mi__ ...
,flag_image_delta_upd_ampm_Mi__ ...
,image_I_value_ampm_Mi__ ...
,image_X_value_ampm_Mi__ ...
,image_S_index_ampm_Mi__ ...
,~ ...
,~ ...
,~ ...
,X_ampm_SMi___ ...
,delta_x_ampm_SMi___ ...
,delta_y_ampm_SMi___ ...
,gamma_z_ampm_SMi___ ...
,I_value_ampm_SMi___ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 parameter_align ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.k_p_r_max ...
,ampm_.weight_3d_k_p_r_ ...
,ampm_.weight_2d_k_p_r_ ...
,ampm_.n_w_ ...
,ampm_.n_M ...
,N_k_p_use__ ...
,ampm_.n_cluster ...
,ampm_.index_ncluster_from_nCTF_ ...
,ampm_.n_CTF ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_r_kC__ ...
,ampm_.l_max_ ...
,a_k_Y_ampm_ ...
,ampm_.euler_polar_a_ampm_Ma__(:,1+na) ...
,ampm_.euler_azimu_b_ampm_Ma__(:,1+na) ...
,ampm_.euler_gamma_z_ampm_Ma__(:,1+na) ...
,ampm_.image_delta_x_ampm_Ma__(:,1+na) ...
,ampm_.image_delta_y_ampm_Ma__(:,1+na) ...
);
if (flag_verbose); disp(sprintf(' %% corr(a_k_Y_ampm_,a_k_Y_ampm_reco_yki__(:,end)): %0.16f',corr(a_k_Y_ampm_,a_k_Y_ampm_reco_yki__(:,end)))); end;

X_SM__ = X_SMi___(:,:,end);
X_ampm_alig_SM__ = X_ampm_alig_SMi___(:,:,end);
X_ampm_SM__ = X_ampm_SMi___(:,:,end);
%%%%;
P_SM__ = X_SM__;
for nM=0:n_M-1;
tmp_S_ = P_SM__(:,1+nM);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
P_SM__(:,1+nM) = tmp_ij_;
end;%for nM=0:n_M-1;
%%%%;
P_ampm_alig_SM__ = X_ampm_alig_SM__;
for nM=0:n_M-1;
tmp_S_ = P_ampm_alig_SM__(:,1+nM);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
P_ampm_alig_SM__(:,1+nM) = tmp_ij_;
end;%for nM=0:n_M-1;
%%%%;
P_ampm_SM__ = X_ampm_SM__;
for nM=0:n_M-1;
tmp_S_ = P_ampm_SM__(:,1+nM);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
P_ampm_SM__(:,1+nM) = tmp_ij_;
end;%for nM=0:n_M-1;
%%%%;
Q_SM__ = X_SM__;
for nS=0:n_S-1;
tmp_S_ = Q_SM__(1+nS,:);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
Q_SM__(1+nS,:) = tmp_ij_;
end;%for nS=0:n_S-1;
%%%%;
Q_ampm_alig_SM__ = X_ampm_alig_SM__;
for nS=0:n_S-1;
tmp_S_ = Q_ampm_alig_SM__(1+nS,:);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
Q_ampm_alig_SM__(1+nS,:) = tmp_ij_;
end;%for nS=0:n_S-1;
%%%%;
Q_ampm_SM__ = X_ampm_SM__;
for nS=0:n_S-1;
tmp_S_ = Q_ampm_SM__(1+nS,:);
[~,tmp_ij_] = sort(tmp_S_); [~,tmp_ij_] = sort(tmp_ij_);
Q_ampm_SM__(1+nS,:) = tmp_ij_;
end;%for nS=0:n_S-1;
%%%%;
CX_S_ = zeros(n_S,1);
CP_S_ = zeros(n_S,1);
CQ_S_ = zeros(n_S,1);
for nS=0:n_S-1;
CX_S_(1+nS) = corr(transpose(X_SM__(1+nS,:)),transpose(X_ampm_alig_SM__(1+nS,:)));
CP_S_(1+nS) = corr(transpose(P_SM__(1+nS,:)),transpose(P_ampm_alig_SM__(1+nS,:)));
CQ_S_(1+nS) = corr(transpose(Q_SM__(1+nS,:)),transpose(Q_ampm_alig_SM__(1+nS,:)));
end;%for nS=0:n_S-1;
CX_M_ = zeros(n_M,1);
CP_M_ = zeros(n_M,1);
CQ_M_ = zeros(n_M,1);
for nM=0:n_M-1;
CX_M_(1+nM) = corr((X_SM__(:,1+nM)),(X_ampm_alig_SM__(:,1+nM)));
CP_M_(1+nM) = corr((P_SM__(:,1+nM)),(P_ampm_alig_SM__(:,1+nM)));
CQ_M_(1+nM) = corr((Q_SM__(:,1+nM)),(Q_ampm_alig_SM__(:,1+nM)));
end;%for nM=0:n_M-1;

save(fname_mat ...
     ,'tmp_delta_r_max' ...
     ,'tmp_delta_r_upb' ...
     ,'n_M' ...
     ,'n_iteration_orig' ...
     ,'n_iteration_seco' ...
     ,'corr_a_k_Y_quad_vs_a_k_Y_ampm_alig' ...
     ,'corr_a_k_Y_ampm_alig_vs_a_k_Y_ampm_alig_reco_yk_i_' ...
     ,'corr_a_k_Y_ampm_vs_a_k_Y_ampm_reco_yk_i_' ...
     ,'corr_a_k_Y_quad_vs_a_k_Y_ampm_alig_reco' ...
     ,'n_S' ...
     ,'template_viewing_azimu_b_all_' ...
     ,'template_viewing_polar_a_all_' ...
     ,'X_SM__' ...
     ,'X_ampm_alig_SM__' ...
     ,'X_ampm_SM__' ...
     ,'P_SM__' ...
     ,'P_ampm_alig_SM__' ...
     ,'P_ampm_SM__' ...
     ,'Q_SM__' ...
     ,'Q_ampm_alig_SM__' ...
     ,'Q_ampm_SM__' ...
     ,'CX_S_' ...
     ,'CP_S_' ...
     ,'CQ_S_' ...
     ,'CX_M_' ...
     ,'CP_M_' ...
     ,'CQ_M_' ...
     );
%     ,'X_SMi___' ...
%     ,'delta_x_SMi___' ...
%     ,'delta_y_SMi___' ...
%     ,'gamma_z_SMi___' ...
%     ,'I_value_SMi___' ...
%     ,'X_ampm_alig_SMi___' ...
%     ,'delta_x_ampm_alig_SMi___' ...
%     ,'delta_y_ampm_alig_SMi___' ...
%     ,'gamma_z_ampm_alig_SMi___' ...
%     ,'I_value_ampm_alig_SMi___' ...
%     ,'X_ampm_SMi___' ...
%     ,'delta_x_ampm_SMi___' ...
%     ,'delta_y_ampm_SMi___' ...
%     ,'gamma_z_ampm_SMi___' ...
%     ,'I_value_ampm_SMi___' ...
close_fname_tmp(fname_pre);
end;%if ~flag_skip;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
% Now make figures. ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
if  exist(fname_mat);
%%%%%%%%%%%%%%%%;
tmp_ = load(fname_mat);

flag_disp=0;
if flag_disp;
fname_fig = sprintf('%s_FIGA',fname_pre);
if flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
figure(1+nf);nf=nf+1;; clf; figbig;
p_row = 6; p_col = 8; np=0;
tmp_n_M = min(tmp_.n_M,floor(p_row*p_col/2));
index_nM_ = 0:tmp_n_M-1;
Xlim_ = prctile(tmp_.X_SM__(:,1+index_nM_),[ 1,99],'all');
for nM=0:tmp_n_M-1;
subplot(p_row,p_col,1+np); np=np+1; if (np>=p_row*p_col); np=0; end;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,tmp_.X_SM__(:,1+nM));
axis image; axisnotick; title(sprintf('nM%.4d (true)',nM));
subplot(p_row,p_col,1+np); np=np+1; if (np>=p_row*p_col); np=0; end;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,tmp_.X_ampm_alig_SM__(:,1+nM));
axis image; axisnotick; title(sprintf('nM%.4d (ampm)',nM));
drawnow();
end;%for nM=0:tmp_n_M-1;
sgtitle(fname_fig,'Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
close gcf;
end;%if flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file');
end;%if flag_disp;

flag_disp=1;
if flag_disp;
fname_fig = sprintf('%s_FIGB',fname_pre);
if (flag_replot>=1) | ~exist(sprintf('%s.jpg',fname_fig),'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
figure(1+nf);nf=nf+1;;clf;set(gcf,'Position',1+[0,0,1024+512,512]);fig80s;
fontsize_use = 24;
n_h = 128; dh = 1/max(1,n_h); n_z = 1.0/max(1,n_h*n_h + numel(tmp_.X_SM__))/(dh*dh); %<-- now integrates to 1. ;
hlim_ = [-4,+4];
%%%%;
subplot(1,3,1);
imagesc(log(n_z) + log(1+hist2d_0(tmp_.X_SM__(:),tmp_.X_ampm_alig_SM__(:),n_h,n_h,[0,+1],[0,+1])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('corr (true)','Interpreter','none');
ylabel('corr (AMPM)','Interpreter','none');
title('Correlation','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,2);
imagesc(log(n_z) + log(1+hist2d_0(tmp_.P_SM__(:),tmp_.P_ampm_alig_SM__(:),n_h,n_h,[1,tmp_.n_S],[1,tmp_.n_S])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('rank (true)','Interpreter','none');
ylabel('rank (AMPM)','Interpreter','none');
title('template rank','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
subplot(1,3,3);
imagesc(log(n_z) + log(1+hist2d_0(tmp_.Q_SM__(:),tmp_.Q_ampm_alig_SM__(:),n_h,n_h,[1,tmp_.n_M],[1,tmp_.n_M])),hlim_);
set(gca,'Ydir','normal'); tmp_c_= colorbar; set(tmp_c_,'Ticks',[-4,+4]);
xlabel('rank (true)','Interpreter','none');
ylabel('rank (AMPM)','Interpreter','none');
title('image rank','Interpreter','none');
axis image; set(gca,'XTick',[1,n_h],'XTickLabel',[0,1],'YTick',[1,n_h],'YTickLabel',[0,1]);
set(gca,'FontSize',fontsize_use);
%%%%;
sgtitle(fname_fig,'Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
sgtitle('');
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_viewing_angle_distribution',string_root);
fname_fig_jpg_strip = sprintf('%s/%s_viewing_angle_distribution_FIGB_strip.jpg',tmp_dir,fname_prefix_str_fname_mat_inf);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close gcf;
end;%if flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file');
end;%if flag_disp;

flag_disp=1;
if flag_disp;
fname_fig = sprintf('%s_FIGC',fname_pre);
if (flag_replot>=1) | ~exist(sprintf('%s.jpg',fname_fig),'file');
if (flag_verbose); disp(sprintf(' %% %s not found, creating',fname_fig)); end;
%%%%;
figure(1+nf);nf=nf+1;;clf;set(gcf,'Position',1+[0,0,768,512]);fig80s; fontsize_use = 16;
p_row = 2; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.X_SM__,2));
xlabel('azimuthal','Interpreter','none');
ylabel('polar','Interpreter','none');
axis image; axisnotick; title('corr (true)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.X_ampm_alig_SM__,2));
xlabel('azimuthal','Interpreter','none');
ylabel('polar','Interpreter','none');
axis image; axisnotick; title('corr (AMPM)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
%subplot(p_row,p_col,1+np);np=np+1;
%imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.X_ampm_SM__,2));
%xlabel('azimuthal','Interpreter','none');
%ylabel('polar','Interpreter','none');
%axis image; axisnotick; title('corr (unal)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.P_SM__,2));
xlabel('azimuthal','Interpreter','none');
ylabel('polar','Interpreter','none');
axis image; axisnotick; title('template rank (true)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np);np=np+1;
imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.P_ampm_alig_SM__,2));
xlabel('azimuthal','Interpreter','none');
ylabel('polar','Interpreter','none');
axis image; axisnotick; title('template rank (AMPM)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
%subplot(p_row,p_col,1+np);np=np+1;
%imagesc_polar_a_azimu_b_0(tmp_.template_viewing_polar_a_all_,tmp_.template_viewing_azimu_b_all_,mean(tmp_.P_ampm_SM__,2));
%xlabel('azimuthal','Interpreter','none');
%ylabel('polar','Interpreter','none');
%axis image; axisnotick; title('template rank (unal)','Interpreter','none'); set(gca,'FontSize',fontsize_use);
%%%%;
sgtitle(fname_fig,'Interpreter','none');
print('-djpeg',sprintf('%s.jpg',fname_fig));
print('-depsc',sprintf('%s.eps',fname_fig));
sgtitle(sprintf('%s',fname_prefix(1:strfind(fname_prefix,'_')-1)),'Interpreter','none','FontSize',fontsize_use+8);
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_viewing_angle_X_S',string_root);
fname_fig_jpg_strip = sprintf('%s/%s_viewing_angle_distribution_FIGC_strip.jpg',tmp_dir,fname_prefix_str_fname_mat_inf);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close gcf;
end;%if flag_replot | ~exist(sprintf('%s.jpg',fname_fig),'file');
end;%if flag_disp;

flag_disp = 0;
if flag_disp;
%%%%%%%%;
fname_vol_pre = sprintf('%s_vol',str_fname_mat_pre);
fname_vol_fig_pre = sprintf('%s_FIGD',fname_vol_pre);
fname_vol_fig_jpg = sprintf('%s.jpg',fname_vol_fig_pre);
if (flag_replot_vol | ~exist(fname_vol_fig_jpg,'file'));
%%%%%%%%;
test_viewing_angle_distribution_vol_helper_0; %<-- construct a_x_u_ampm_. ;
%%%%%%%%;
if strfind(ampm_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(ampm_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(ampm_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(ampm_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
%percent_threshold_ = 85:1.25:97.5;
n_percent_threshold = numel(percent_threshold_);
%%%%;
n_x_u_pack = ampm_.n_x_u_pack;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
fontsize_use = 32;
if (n_percent_threshold> 1); p_row = 6; p_col = n_percent_threshold; end;
if (n_percent_threshold==1); p_row = 1; p_col = 6; end;
for npt=0:n_percent_threshold-1;
tmp_percent_threshold = percent_threshold_(1+npt); nnp=0;
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(ampm_.a_x_u_base_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('True %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('True')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(ampm_.a_x_u_0qbp_reco_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('OAM %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OAM')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_reco_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('OAMS %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OAMS')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
if exist('a_x_u_sher_from_quad_unst_alig_','var');
isosurface_f_x_u_0(a_x_u_sher_from_quad_unst_alig_(1+tmp_index_),tmp_percent_threshold);
end;%if exist('a_x_u_sher_from_quad_unst_alig_','var');
%title(sprintf('OBI %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OBI')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
if exist('a_x_u_sher_from_quad_stab_alig_','var');
isosurface_f_x_u_0(a_x_u_sher_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
end;%if exist('a_x_u_sher_from_quad_stab_alig_','var');
%title(sprintf('OBIS %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OBIS')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
%subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
%if exist('a_x_u_zero_from_quad_stab_alig_','var');
%isosurface_f_x_u_0(a_x_u_zero_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%if exist('a_x_u_zero_from_quad_stab_alig_','var');
%title(sprintf('zero %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%title(sprintf('zero')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
%subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
%if exist('a_x_u_0qbp_from_quad_stab_alig_','var');
%isosurface_f_x_u_0(a_x_u_0qbp_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%if exist('a_x_u_0qbp_from_quad_stab_alig_','var');
%title(sprintf('0qbp %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%title(sprintf('0qbp')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_ampm_alig_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('AMPM %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('AMPM')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
end;%for npt=0:n_percent_threshold-1;
set(gcf,'Position',1+[0,0,512*5,512+64]);
sgtitle(fname_vol_fig_pre,'Interpreter','none');
print('-djpeg',fname_vol_fig_jpg);
sgtitle(sprintf('%s',fname_prefix(1:strfind(fname_prefix,'_')-1)),'Interpreter','none','FontSize',fontsize_use+16);
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_vol',string_root);
fname_fig_jpg_strip = sprintf('%s/%s_vol_FIGD_strip.jpg',tmp_dir,fname_prefix_str_fname_mat_inf);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close gcf;
%%%%%%%%;
end;%if (flag_replot_vol | ~exist(fname_vol_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%;

flag_disp = 1;
if flag_disp;
%%%%%%%%;
fname_vol_pre = sprintf('%s_vol',str_fname_mat_pre);
fname_vol_fig_pre = sprintf('%s_FIGF',fname_vol_pre);
fname_vol_fig_jpg = sprintf('%s.jpg',fname_vol_fig_pre);
if (flag_replot_vol | ~exist(fname_vol_fig_jpg,'file'));
%%%%%%%%;
test_viewing_angle_distribution_vol_helper_0; %<-- construct a_x_u_ampm_. ;
%%%%%%%%;
if strfind(ampm_.dir_pm,'ISWINCP');    percent_threshold_ = 94.50; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'trpv1');      percent_threshold_ = 91.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'rib80s');     percent_threshold_ = 86.25; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'MlaFEDB');    percent_threshold_ = 95.00; tmp_nx = 8; end;
if strfind(ampm_.dir_pm,'LetB1');      percent_threshold_ = 91.75; tmp_nx = 14; end;
if strfind(ampm_.dir_pm,'TMEM16F');    percent_threshold_ = 94.50; tmp_nx = 16; end;
if strfind(ampm_.dir_pm,'LSUbl17dep'); percent_threshold_ = 86.25; tmp_nx = 15; end;
if strfind(ampm_.dir_pm,'ps1');        percent_threshold_ = 96.00; tmp_nx = 11; end;
%percent_threshold_ = 85:1.25:97.5;
n_percent_threshold = numel(percent_threshold_);
%%%%;
n_x_u_pack = ampm_.n_x_u_pack;
tmp_window_ = zeros(n_x_u_pack,n_x_u_pack,n_x_u_pack);
tmp_index_ = tmp_nx:n_x_u_pack-1-tmp_nx;
tmp_window_(1+tmp_index_,1+tmp_index_,1+tmp_index_)=1;
tmp_index_ = efind(tmp_window_);
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
fontsize_use = 32;
if (n_percent_threshold> 1); p_row = 5; p_col = n_percent_threshold; end;
if (n_percent_threshold==1); p_row = 1; p_col = 5; end;
for npt=0:n_percent_threshold-1;
tmp_percent_threshold = percent_threshold_(1+npt); nnp=0;
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(ampm_.a_x_u_base_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('True %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('True')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(ampm_.a_x_u_0qbp_reco_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('Oracle %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('Oracle')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_reco_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('OAM %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OAM')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
%subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
%if exist('a_x_u_sher_from_quad_unst_alig_','var');
%isosurface_f_x_u_0(a_x_u_sher_from_quad_unst_alig_(1+tmp_index_),tmp_percent_threshold);
%end;%if exist('a_x_u_sher_from_quad_unst_alig_','var');
%%title(sprintf('OBI %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%title(sprintf('OBI')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
if exist('a_x_u_sher_from_quad_stab_alig_','var');
isosurface_f_x_u_0(a_x_u_sher_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
end;%if exist('a_x_u_sher_from_quad_stab_alig_','var');
%title(sprintf('OBI %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('OBI')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
%subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
%if exist('a_x_u_zero_from_quad_stab_alig_','var');
%isosurface_f_x_u_0(a_x_u_zero_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%if exist('a_x_u_zero_from_quad_stab_alig_','var');
%title(sprintf('zero %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%title(sprintf('zero')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
%subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
%if exist('a_x_u_0qbp_from_quad_stab_alig_','var');
%isosurface_f_x_u_0(a_x_u_0qbp_from_quad_stab_alig_(1+tmp_index_),tmp_percent_threshold);
%if exist('a_x_u_0qbp_from_quad_stab_alig_','var');
%title(sprintf('0qbp %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%title(sprintf('0qbp')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
subplot(p_row,p_col,1+npt+nnp*n_percent_threshold);nnp=nnp+1;
isosurface_f_x_u_0(a_x_u_ampm_alig_(1+tmp_index_),tmp_percent_threshold);
%title(sprintf('EMPM %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
title(sprintf('EMPM')); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%;
end;%for npt=0:n_percent_threshold-1;
set(gcf,'Position',1+[0,0,512*5,512+64]);
sgtitle(fname_vol_fig_pre,'Interpreter','none');
print('-djpeg',fname_vol_fig_jpg);
sgtitle(sprintf('%s',fname_prefix(1:strfind(fname_prefix,'_')-1)),'Interpreter','none','FontSize',fontsize_use+16);
tmp_dir = sprintf('/%s/rangan/dir_cryoem/dir_ampm_manuscript/dir_ampm_fig_vol',string_root);
fname_fig_jpg_strip = sprintf('%s/%s_vol_FIGF_strip.jpg',tmp_dir,fname_prefix_str_fname_mat_inf);
disp(sprintf(' %% writing %s',fname_fig_jpg_strip));
print('-djpeg',sprintf('%s',fname_fig_jpg_strip));
close gcf;
%%%%%%%%;
end;%if (flag_replot_vol | ~exist(fname_vol_fig_jpg,'file'));
end;%if flag_disp;
%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%if  exist(fname_mat);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%for nna=0:n_na-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



