%%%%%%%%;
% based on test_viewing_angle_distribution_0.m ;
% trying ring-restriction applied to LetB1. ;
%%%%%%%%;

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
flag_replot_vol = 1;
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
%for nexperiment=(randperm(n_experiment)-1);
nexperiment = 4; %<-- LetB1. ;
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
n_iteration_seco = max(1,min(n_iteration_orig-2,6)); %<-- keep at 6?. ;

%%%%%%%%;
% sort. ;
%%%%%%%%;
[~,ij_X_best_ampm_srt_] = sort(ampm_.X_best_ampm_ia__(end,:),'descend'); index_X_best_ampm_srt_ = ij_X_best_ampm_srt_-1;

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
fname_pre = sprintf('%s_viewing_angle_ring',str_fname_mat_pre);
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
delta_x_SM__ = delta_x_SMi___(:,:,end);
delta_x_ampm_alig_SM__ = delta_x_ampm_alig_SMi___(:,:,end);
delta_x_ampm_SM__ = delta_x_ampm_SMi___(:,:,end);
delta_y_SM__ = delta_y_SMi___(:,:,end);
delta_y_ampm_alig_SM__ = delta_y_ampm_alig_SMi___(:,:,end);
delta_y_ampm_SM__ = delta_y_ampm_SMi___(:,:,end);
gamma_z_SM__ = gamma_z_SMi___(:,:,end);
gamma_z_ampm_alig_SM__ = gamma_z_ampm_alig_SMi___(:,:,end);
gamma_z_ampm_SM__ = gamma_z_ampm_SMi___(:,:,end);
I_value_SM__ = I_value_SMi___(:,:,end);
I_value_ampm_alig_SM__ = I_value_ampm_alig_SMi___(:,:,end);
I_value_ampm_SM__ = I_value_ampm_SMi___(:,:,end);
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
     ,'delta_x_SM__' ...
     ,'delta_x_ampm_alig_SM__' ...
     ,'delta_x_ampm_SM__' ...
     ,'delta_y_SM__' ...
     ,'delta_y_ampm_alig_SM__' ...
     ,'delta_y_ampm_SM__' ...
     ,'gamma_z_SM__' ...
     ,'gamma_z_ampm_alig_SM__' ...
     ,'gamma_z_ampm_SM__' ...
     ,'I_value_SM__' ...
     ,'I_value_ampm_alig_SM__' ...
     ,'I_value_ampm_SM__' ...
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

%%%%%%%%;
% indices. ;
%%%%%%%%;
n_k_all = ampm_.n_k_all;
n_k_all_csum_ = ampm_.n_k_all_csum_;
k_p_r_all_ = ampm_.k_p_r_all_;
k_p_azimu_b_all_ = ampm_.k_p_azimu_b_all_;
k_p_polar_a_all_ = ampm_.k_p_polar_a_all_;
weight_3d_k_all_ = ampm_.weight_3d_k_all_;
weight_shell_k_ = ampm_.weight_shell_k_;
n_k_p_r = ampm_.n_k_p_r;
k_p_r_ = ampm_.k_p_r_;
k_p_r_max = ampm_.k_p_r_max;
weight_3d_k_p_r_ = ampm_.weight_3d_k_p_r_;
l_max_ = ampm_.l_max_;
n_w_ = ampm_.n_w_;
n_w_max = ampm_.n_w_max;
%%%%%%%%;
n_k_all = ampm_.n_k_all;
weight_3d_k_all_ = ampm_.weight_3d_k_all_;
k_c_0_all_ = ampm_.k_c_0_all_;
k_c_1_all_ = ampm_.k_c_1_all_;
k_c_2_all_ = ampm_.k_c_2_all_;
k_p_r_max = ampm_.k_p_r_max;
n_xxx_u = ampm_.n_xxx_u;
xxx_u_weight_ = ampm_.xxx_u_weight_;
x_u_0___ = ampm_.x_u_0___;
x_u_1___ = ampm_.x_u_1___;
x_u_2___ = ampm_.x_u_2___;
x_p_r_max = ampm_.x_p_r_max;
%%%%%%%%;
qbp_eps = 1e-3;

%%%%;
tmp_ = load(fname_mat);
%%%%;

%%%%%%%%;
tmp_t = tic();
parameter_ampmh_ring = struct('type','parameter');
parameter_ampmh_ring.flag_euler_polar_a_restrict = 1;
parameter_ampmh_ring.euler_polar_a_restrict_band = pi/12;
[ ...
 parameter_ampmh_ring ...
,tmp_ring_euler_polar_a_ ...
,tmp_ring_euler_azimu_b_ ...
,tmp_ring_euler_gamma_z_ ...
,tmp_ring_image_delta_x_ ...
,tmp_ring_image_delta_y_ ...
,tmp_ring_image_I_value_ ...
,tmp_ring_image_X_value_ ...
,tmp_ring_image_S_index_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter_ampmh_ring ...
,n_w_max ...
,tmp_.n_S ...
,tmp_.template_viewing_polar_a_all_ ...
,tmp_.template_viewing_azimu_b_all_ ...
,tmp_.n_M ...
,tmp_.X_ampm_alig_SM__ ...
,tmp_.delta_x_ampm_alig_SM__ ...
,tmp_.delta_y_ampm_alig_SM__ ...
,tmp_.gamma_z_ampm_alig_SM__ ...
,tmp_.I_value_ampm_alig_SM__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% ampmh_MS_vs_SM_2: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
[ ...
 tmp_ring_a_k_Y_ampm_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.l_max_ ...
,ampm_.n_w_ ...
,ampm_.n_M ...
,N_k_p_use__ ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_wkC__ ...
,tmp_ring_euler_polar_a_ ...
,tmp_ring_euler_azimu_b_ ...
,tmp_ring_euler_gamma_z_ ...
,tmp_ring_image_delta_x_ ...
,tmp_ring_image_delta_y_ ...
,tmp_ring_image_I_value_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% qbp_6: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
[ ... 
 tmp_ring_a_k_Y_ampm_alig_ ...
,tmp_corr_ring_a_k_Y_ampm_alig_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,ampm_.a_k_Y_quad_ ...
,tmp_ring_a_k_Y_ampm_ ...
,0 ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% spharm_register_and_rotate_2: %0.6fs',tmp_t));
%%%%%%%%
if ~exist('Ylm_klma___','var'); Ylm_klma___ = []; end;
tmp_t = tic();
[ ...
 tmp_ring_a_k_p_ampm_alig_  ...
,Ylm_klma___ ...
] =  ...
convert_spharm_to_k_p_3( ...
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
,tmp_ring_a_k_Y_ampm_alig_ ...
,Ylm_klma___ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_3: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic(); 
[ ...
 tmp_ring_a_x_u_ampm_alig_ ...
] = ...
convert_k_p_to_x_c_1( ...
 0*flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,xxx_u_weight_ ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,x_p_r_max ...
,tmp_ring_a_k_p_ampm_alig_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_x_c_1: %0.6fs',tmp_t));
%%%%%%%%;

%%%%%%%%;
tmp_t = tic();
parameter_ampmh_unif = struct('type','parameter');
parameter_ampmh_unif.flag_euler_polar_a_restrict = 0;
parameter_ampmh_unif.euler_polar_a_restrict_band = pi/12;
[ ...
 parameter_ampmh_unif ...
,tmp_unif_euler_polar_a_ ...
,tmp_unif_euler_azimu_b_ ...
,tmp_unif_euler_gamma_z_ ...
,tmp_unif_image_delta_x_ ...
,tmp_unif_image_delta_y_ ...
,tmp_unif_image_I_value_ ...
,tmp_unif_image_X_value_ ...
,tmp_unif_image_S_index_ ...
] = ...
ampmh_MS_vs_SM_2( ...
 parameter_ampmh_unif ...
,n_w_max ...
,tmp_.n_S ...
,tmp_.template_viewing_polar_a_all_ ...
,tmp_.template_viewing_azimu_b_all_ ...
,tmp_.n_M ...
,tmp_.X_ampm_alig_SM__ ...
,tmp_.delta_x_ampm_alig_SM__ ...
,tmp_.delta_y_ampm_alig_SM__ ...
,tmp_.gamma_z_ampm_alig_SM__ ...
,tmp_.I_value_ampm_alig_SM__ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% ampmh_MS_vs_SM_2: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
[ ...
 tmp_unif_a_k_Y_ampm_ ...
] = ...
qbp_6( ...
 qbp_eps ...
,ampm_.n_k_p_r ...
,ampm_.k_p_r_ ...
,ampm_.l_max_ ...
,ampm_.n_w_ ...
,ampm_.n_M ...
,N_k_p_use__ ...
,ampm_.index_nCTF_from_nM_ ...
,ampm_.CTF_k_p_wkC__ ...
,tmp_unif_euler_polar_a_ ...
,tmp_unif_euler_azimu_b_ ...
,tmp_unif_euler_gamma_z_ ...
,tmp_unif_image_delta_x_ ...
,tmp_unif_image_delta_y_ ...
,tmp_unif_image_I_value_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% qbp_6: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic();
[ ... 
 tmp_unif_a_k_Y_ampm_alig_ ...
,tmp_corr_unif_a_k_Y_ampm_alig_ ...
] = ...
spharm_register_and_rotate_2( ...
 n_k_p_r ...
,k_p_r_ ...
,k_p_r_max ...
,weight_3d_k_p_r_ ...
,l_max_ ...
,ampm_.a_k_Y_quad_ ...
,tmp_unif_a_k_Y_ampm_ ...
,0 ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% spharm_register_and_rotate_2: %0.6fs',tmp_t));
%%%%%%%%
if ~exist('Ylm_klma___','var'); Ylm_klma___ = []; end;
tmp_t = tic();
[ ...
 tmp_unif_a_k_p_ampm_alig_  ...
,Ylm_klma___ ...
] =  ...
convert_spharm_to_k_p_3( ...
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
,tmp_unif_a_k_Y_ampm_alig_ ...
,Ylm_klma___ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_3: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic(); 
[ ...
 tmp_unif_a_x_u_ampm_alig_ ...
] = ...
convert_k_p_to_x_c_1( ...
 0*flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,xxx_u_weight_ ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,x_p_r_max ...
,tmp_unif_a_k_p_ampm_alig_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_x_c_1: %0.6fs',tmp_t));
%%%%%%%%;

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
fontsize_use = 12; tmp_percent_threshold = percent_threshold_(1);
p_row = 1; p_col = 2; np=0;
%%%%;
subplot(p_row,p_col,1+np); np=np+1;
isosurface_f_x_u_0(tmp_unif_a_x_u_ampm_alig_(1+tmp_index_),tmp_percent_threshold);
title(sprintf('Unif %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
subplot(p_row,p_col,1+np); np=np+1;
isosurface_f_x_u_0(tmp_ring_a_x_u_ampm_alig_(1+tmp_index_),tmp_percent_threshold);
title(sprintf('Ring %.3f',tmp_percent_threshold)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%%%%%;

%%%%%%%%;
% now iterate ring. ;
%%%%%%%%;
n_iteration_thir = 6;
tmp_ring_a_k_Y_ampm_pre_ = tmp_ring_a_k_Y_ampm_;
tmp_ring_euler_polar_a_pre_ = tmp_ring_euler_polar_a_;
tmp_ring_euler_azimu_b_pre_ = tmp_ring_euler_azimu_b_;
tmp_ring_euler_gamma_z_pre_ = tmp_ring_euler_gamma_z_;
tmp_ring_image_delta_x_pre_ = tmp_ring_image_delta_x_;
tmp_ring_image_delta_y_pre_ = tmp_ring_image_delta_y_;
tmp_ring_image_I_value_pre_ = tmp_ring_image_I_value_;
%%%%%%%%%%%%%%%%;
for niteration_thir=0:n_iteration_thir-1;
%%%%%%%%%%%%%%%%;
if (flag_verbose); disp(sprintf(' %% niteration_thir = %d/%d',niteration_thir,n_iteration_thir)); end;
%%%%%%%%;
% Now align M_k_p__ to a_k_Y_quad_. ;
%%%%%%%%;
tmp_parameter_align = struct('type','parameter');
tmp_parameter_align.fname_align_a_k_Y_pre = ''; %<-- empty for now. ;
tmp_parameter_align.delta_r_max = tmp_delta_r_max;
tmp_parameter_align.n_delta_v_requested = 24;
tmp_parameter_align.delta_r_upb = tmp_delta_r_upb;
tmp_parameter_align.n_iteration = 1;
tmp_parameter_align.flag_euler_polar_a_restrict = 1;
[ ...
 tmp_parameter_align ...
,tmp_ring_a_k_Y_ampm_pos_ ...
,tmp_corr_ring_a_k_Y_ampm_pos ...
,tmp_ring_euler_polar_a_pos_ ...
,tmp_ring_euler_azimu_b_pos_ ...
,tmp_ring_euler_gamma_z_pos_ ...
,tmp_ring_image_delta_x_acc_pos_ ...
,tmp_ring_image_delta_y_acc_pos_ ...
,tmp_ring_image_delta_x_upd_pos_ ...
,tmp_ring_image_delta_y_upd_pos_ ...
,tmp_ring_flag_image_delta_upd_pos_ ...
,tmp_ring_image_I_value_pos_ ...
,tmp_ring_image_X_value_pos_ ...
,~ ...
,~ ...
,~ ...
,~ ...
,tmp_ring_X_SM__ ...
,tmp_ring_delta_x_SM___ ...
,tmp_ring_delta_y_SM___ ...
,tmp_ring_gamma_z_SM___ ...
,tmp_ring_I_value_SM__ ...
] = ...
pm_align_M_k_p_to_a_k_Y_3( ...
 tmp_parameter_align ...
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
,tmp_ring_a_k_Y_ampm_pre_ ...
,tmp_ring_euler_polar_a_pre_ ...
,tmp_ring_euler_azimu_b_pre_ ...
,tmp_ring_euler_gamma_z_pre_ ...
,tmp_ring_image_delta_x_pre_ ...
,tmp_ring_image_delta_y_pre_ ...
);
%%%%%%%%;
tmp_ring_a_k_Y_ampm_pre_ = tmp_ring_a_k_Y_ampm_pos_;
tmp_ring_euler_polar_a_pre_ = tmp_ring_euler_polar_a_pos_;
tmp_ring_euler_azimu_b_pre_ = tmp_ring_euler_azimu_b_pos_;
tmp_ring_euler_gamma_z_pre_ = tmp_ring_euler_gamma_z_pos_;
tmp_ring_image_delta_x_pre_ = tmp_ring_image_delta_x_acc_pos_ + tmp_ring_image_delta_x_upd_pos_ ;
tmp_ring_image_delta_y_pre_ = tmp_ring_image_delta_y_acc_pos_ + tmp_ring_image_delta_y_upd_pos_ ;
%%%%%%%%;
if flag_disp;
%%%%%%%%;
if ~exist('Ylm_klma___','var'); Ylm_klma___ = []; end;
tmp_t = tic();
[ ...
 tmp_ring_a_k_p_ampm_pre_  ...
,Ylm_klma___ ...
] =  ...
convert_spharm_to_k_p_3( ...
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
,tmp_ring_a_k_Y_ampm_pre_ ...
,Ylm_klma___ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_spharm_to_k_p_3: %0.6fs',tmp_t));
%%%%%%%%;
tmp_t = tic(); 
[ ...
 tmp_ring_a_x_u_ampm_pre_ ...
] = ...
convert_k_p_to_x_c_1( ...
 0*flag_verbose ...
,n_k_all ...
,weight_3d_k_all_ ...
,k_c_0_all_ ...
,k_c_1_all_ ...
,k_c_2_all_ ...
,k_p_r_max ...
,n_xxx_u ...
,xxx_u_weight_ ...
,x_u_0___ ...
,x_u_1___ ...
,x_u_2___ ...
,x_p_r_max ...
,tmp_ring_a_k_p_ampm_pre_ ...
);
tmp_t = toc(tmp_t); disp(sprintf(' %% convert_k_p_to_x_c_1: %0.6fs',tmp_t));
%%%%%%%%;
figure(1+nf);nf=nf+1;clf;
fontsize_use = 12; tmp_percent_threshold = percent_threshold_(1);
subplot(1,1,1); isosurface_f_x_u_0(tmp_ring_a_x_u_ampm_pre_(1+tmp_index_),tmp_percent_threshold);
title(sprintf('Ring %.3f ni %d',tmp_percent_threshold,niteration_thir)); xlabel('');ylabel('');zlabel(''); set(gca,'FontSize',fontsize_use);
%%%%%%%%;
end;%if flag_disp;
%%%%%%%%;
%%%%%%%%%%%%%%%%;
end;%for niteration_thir=0:n_iteration_thir-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%if  exist(fname_mat);
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%;
end;%for nna=0:n_na-1;
%%%%%%%%%%%%%%%%;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
%end;%for nexperiment=0:n_experiment-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;



